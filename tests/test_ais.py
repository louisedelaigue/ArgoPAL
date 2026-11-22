# %%
"""Optional ship-noise check with AIS.

All AIS data here are synthetic (made in the tests), in the formats
argopal accepts: a table, decoded CSV (like the 2025 DYFAMED files) and
raw NMEA with tag blocks.
"""

import gzip
import warnings
from functools import reduce

import numpy as np
import pandas as pd
import pytest

import argopal


T0 = pd.Timestamp("2025-02-14 12:00")
LAT0, LON0 = 43.38, 7.83  # float position (DYFAMED)
KM_LAT = 1 / 111.2  # degrees of latitude per km


def _pal(times, lat=LAT0, lon=LON0):
    df = pd.DataFrame(
        {
            "datetime": pd.to_datetime(times),
            "latitude": lat,
            "longitude": lon,
            "f_8000Hz": 45.0,
            "spl_8000Hz": 46.0,
            "U_nystuen2015": 7.0,
        }
    )
    df.attrs = {"bands_hz": [8000]}
    return df


def _ais(rows):
    """rows: (minutes after T0, mmsi, km north of the float)."""
    return pd.DataFrame(
        {
            "datetime": [T0 + pd.Timedelta(minutes=m) for m, _, _ in rows],
            "mmsi": [s for _, s, _ in rows],
            "latitude": [LAT0 + km * KM_LAT for _, _, km in rows],
            "longitude": LON0,
        }
    )


# --- flag_ships -------------------------------------------------------------


def test_ship_within_radius_and_window_is_flagged():
    ais = _ais([(10, 111, 5.0), (-20, 222, 15.0), (0, 333, 30.0)])
    out = argopal.flag_ships(_pal([T0]), ais)
    assert out["ships_nearby"].iloc[0] == 2  # 333 is 30 km away
    assert bool(out["ship_flag"].iloc[0])


def test_ship_counted_once_per_mmsi():
    ais = _ais([(m, 111, 2.0) for m in range(0, 30, 2)])
    assert argopal.flag_ships(_pal([T0]), ais)["ships_nearby"].iloc[0] == 1


def test_ship_outside_window_not_counted():
    ais = _ais([(45, 111, 1.0), (-45, 222, 1.0), (0, 333, 50.0)])
    out = argopal.flag_ships(_pal([T0]), ais)
    assert out["ships_nearby"].iloc[0] == 0
    assert not bool(out["ship_flag"].iloc[0])
    wide = argopal.flag_ships(_pal([T0]), ais, window="1h")
    assert wide["ships_nearby"].iloc[0] == 2


def test_no_ais_coverage_is_unknown_not_zero():
    ais = _ais([(0, 111, 1.0)])
    far_in_time = T0 + pd.Timedelta("3h")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = argopal.flag_ships(_pal([T0, far_in_time]), ais)
    assert out["ships_nearby"].iloc[0] == 1
    assert np.isnan(out["ships_nearby"].iloc[1])
    assert out.attrs["ais"]["n_without_coverage"] == 1


def test_remove_blanks_acoustic_values_of_flagged_rows():
    ais = _ais([(0, 111, 1.0), (180, 222, 50.0)])
    pal = _pal([T0, T0 + pd.Timedelta("3h")])
    out = argopal.flag_ships(pal, ais, action="remove")
    assert (
        out[["f_8000Hz", "spl_8000Hz", "U_nystuen2015"]].iloc[0].isna().all()
    )
    assert (
        out[["f_8000Hz", "spl_8000Hz", "U_nystuen2015"]].iloc[1].notna().all()
    )
    with pytest.raises(ValueError, match="action"):
        argopal.flag_ships(pal, ais, action="drop")


def test_needs_positions():
    pal = _pal([T0]).drop(columns=["latitude", "longitude"])
    with pytest.raises(KeyError, match="attach_positions"):
        argopal.flag_ships(pal, _ais([(0, 111, 1.0)]))


# --- read_ais: table and CSV -----------------------------------------------------


def test_read_table():
    ais = argopal.read_ais(_ais([(5, 111, 1.0), (0, 222, 2.0)]))
    assert list(ais.columns) == ["datetime", "mmsi", "latitude", "longitude"]
    assert ais["datetime"].is_monotonic_increasing


def _write_2025_style_csv(path):
    """Like the decoded files of the 2025 DYFAMED study: x/y positions,
    UNIX receiver time, AIS 'timestamp' second, base stations (type 4)."""
    t = int(T0.timestamp())
    pd.DataFrame(
        {
            "id": [1, 3, 18, 4, 21],
            "mmsi": [111, 222, 333, 2275000, 992276001],
            "x": [7.83, 7.84, 7.85, 7.86, 7.87],
            "y": [43.38, 43.39, 43.40, 43.41, 43.42],
            "timestamp": [14, 15, 16, 60, 60],  # second of the minute
            "tagblock_timestamp": [t, t + 60, t + 120, t, t],
            "utc_hour": [np.nan, 12, np.nan, np.nan, np.nan],
        }
    ).to_csv(path, index=False)


def test_read_2025_style_csv(tmp_path):
    path = tmp_path / "ais_positions_20250215.csv"
    _write_2025_style_csv(path)
    ais = argopal.read_ais(str(path))
    # base station (4) and aid to navigation (21) dropped; time from the
    # receiver (tagblock_timestamp), not from the AIS 'timestamp' second
    assert sorted(ais["mmsi"]) == [111, 222, 333]
    assert ais["datetime"].iloc[0] == T0
    assert ais["longitude"].iloc[0] == pytest.approx(7.83)


@pytest.mark.parametrize(
    "columns",
    [
        {"BaseDateTime": "t", "LAT": "lat", "LON": "lon", "MMSI": "mmsi"},
        {
            "# Timestamp": "t",
            "Latitude": "lat",
            "Longitude": "lon",
            "MMSI": "mmsi",
        },
    ],
)
def test_read_common_exports(tmp_path, columns):
    # NOAA MarineCadastre and Danish Maritime Authority column names
    data = {"t": [str(T0)], "lat": [43.38], "lon": [7.83], "mmsi": [111]}
    df = pd.DataFrame({k: data[v] for k, v in columns.items()})
    path = tmp_path / "export.csv"
    df.to_csv(path, index=False)
    ais = argopal.read_ais(str(path))
    assert ais["datetime"].iloc[0] == T0 and ais["mmsi"].iloc[0] == 111


def test_csv_without_usable_columns(tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame({"time": [str(T0)], "lat": [43.0], "lon": [7.0]}).to_csv(
        path, index=False
    )
    with pytest.raises(ValueError, match="mmsi"):
        argopal.read_ais(str(path))


def test_csv_with_only_second_of_minute_is_refused(tmp_path):
    path = tmp_path / "seconds.csv"
    pd.DataFrame(
        {"timestamp": [14], "mmsi": [1], "lat": [43.0], "lon": [7.0]}
    ).to_csv(path, index=False)
    with pytest.raises(ValueError, match="UNIX seconds"):
        argopal.read_ais(str(path))


def test_bbox_and_folder(tmp_path):
    _write_2025_style_csv(tmp_path / "a.csv")
    ais = argopal.read_ais(str(tmp_path), bbox=(7.0, 7.845, 43.0, 44.0))
    assert sorted(ais["mmsi"]) == [111, 222]


# --- read_ais: raw NMEA ---------------------------------------------------------


def _tag_block(unix_time):
    body = f"c:{unix_time}"
    checksum = reduce(lambda a, b: a ^ b, body.encode())
    return f"\\{body}*{checksum:02X}\\"


def _write_nmea(path, compress=False):
    pyais = pytest.importorskip("pyais")
    from pyais.encode import encode_dict

    lines = []
    for k, (mmsi, km) in enumerate([(111, 1.0), (222, 3.0)]):
        msg = {"type": 1, "mmsi": mmsi, "lat": LAT0 + km * KM_LAT, "lon": LON0}
        for sentence in encode_dict(msg, radio_channel="A"):
            lines.append(_tag_block(int(T0.timestamp()) + 60 * k) + sentence)
    # a base-station report, which must not count as a ship
    base = {"type": 4, "mmsi": 2275000, "lat": LAT0, "lon": LON0}
    lines.append(_tag_block(int(T0.timestamp())) + encode_dict(base)[0])
    text = "\n".join(lines) + "\n"
    if compress:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return pyais


@pytest.mark.parametrize("compress", [False, True])
def test_read_nmea_with_tag_blocks(tmp_path, compress):
    path = tmp_path / ("messages.gz" if compress else "messages.nmea")
    _write_nmea(path, compress)
    ais = argopal.read_ais(str(path))
    assert sorted(ais["mmsi"]) == [111, 222]  # base station dropped
    assert ais["datetime"].min() == T0
    out = argopal.flag_ships(_pal([T0]), ais)
    assert out["ships_nearby"].iloc[0] == 2


def test_nmea_without_tag_blocks(tmp_path):
    pytest.importorskip("pyais")
    from pyais.encode import encode_dict

    path = tmp_path / "plain.nmea"
    msg = {"type": 1, "mmsi": 111, "lat": LAT0, "lon": LON0}
    path.write_text(encode_dict(msg)[0] + "\n")
    with pytest.raises(ValueError, match="tag-block"):
        argopal.read_ais(str(path))


# --- Pipeline ------------------------------------------------------------------------


def test_pipeline_with_ais(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = argopal.process_float(7902280, mode="raw")
    park = base["park"]
    # one ship 2 km from the float during the first hour of PAL data
    t = park["datetime"].min()
    ais = pd.DataFrame(
        {
            "datetime": [
                t + pd.Timedelta(minutes=m) for m in range(0, 60, 10)
            ],
            "mmsi": 111,
            "latitude": park["latitude"].iloc[0] + 2 * KM_LAT,
            "longitude": park["longitude"].iloc[0],
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = argopal.process_float(7902280, mode="raw", ais=ais)
    p = r["park"]
    assert {"ships_nearby", "ship_flag"} <= set(p.columns)
    assert p["ship_flag"].any()
    assert "ais" in r
    assert "ais" in argopal.plot_all(r)
    # without AIS: no ship columns, no figure
    assert "ships_nearby" not in park
    assert "ais" not in argopal.plot_all(base)


def test_coverage_comes_from_all_received_positions_not_the_box():
    # The only ship is 300 km away: reading with a box around the float
    # drops it, but the AIS feed was receiving, so "0 ships", not unknown.
    far = _ais([(0, 111, 300.0)])
    ais = argopal.read_ais(far, bbox=(7.5, 8.2, 43.0, 43.8))
    assert len(ais) == 0
    out = argopal.flag_ships(_pal([T0]), ais)
    assert out["ships_nearby"].iloc[0] == 0
    assert out.attrs["ais"]["n_without_coverage"] == 0


def test_several_files_coverage_is_merged(tmp_path):
    for name, shift in (("a.csv", 0), ("b.csv", 3600)):
        t = int(T0.timestamp()) + shift
        pd.DataFrame(
            {
                "mmsi": [111],
                "x": [7.83],
                "y": [43.38],
                "tagblock_timestamp": [t],
            }
        ).to_csv(tmp_path / name, index=False)
    ais = argopal.read_ais(str(tmp_path))
    assert len(ais) == 2
    assert len(ais.attrs["received_minutes"]) == 2
