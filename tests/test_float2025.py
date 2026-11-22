# %%
"""Second float: the 2025 sea trial near DYFAMED (Delaigue et al. 2026).

Not on the GDAC, so it is loaded from tests/data/float2025 by the
`float2025` helper.  It covers what 7902280 cannot: full 23-band spectra
(incl. 3.15 kHz), band levels not normalised on board, three park depths,
and real AIS along the track.
"""

import gzip
import math
import warnings

import numpy as np
import pandas as pd
import pytest

import argopal as ap

from . import float2025 as F


@pytest.fixture(scope="module")
def result():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return F.process(bands=(3150, 8000), ais=F.AIS_CSV)


# --- Spectra ------------------------------------------------------------------------


def test_full_23_band_spectra():
    pal = F.read_pal()
    assert pal.attrs["band_layout_counts"] == {
        "23": len(pal),
        "9": 0,
        "other": 0,
    }
    assert pal.attrs["band_order_score"] < -0.95
    assert ap.available_bands(pal) == list(ap.PAL_BANDS_HZ)


def test_3150_hz_is_kept_when_the_float_sends_it():
    with warnings.catch_warnings():
        warnings.simplefilter("error", ap.BandCoverageWarning)
        sel = ap.select_bands(F.read_pal(), bands=(3150, 8000))
    assert sel.attrs["bands_hz"] == [3150, 8000]
    assert sel.attrs["bands_dropped"] == []


def test_despike_on_real_spectra_matches_reference():
    pal = ap.select_bands(F.read_pal(), bands=(8000,))
    pal = pal[pal["cycle"] == 41].reset_index(drop=True)
    pal.attrs = {**pal.attrs, "bands_hz": [8000]}
    fast = ap.despike(pal)
    ref = []
    for t in pal["datetime"]:
        w = pal.loc[
            pal["datetime"].between(
                t - pd.Timedelta("1.5h"), t + pd.Timedelta("1.5h")
            ),
            "f_8000Hz",
        ]
        ref.append(w[w <= w.quantile(0.95)].mean())
    np.testing.assert_allclose(fast["f_8000Hz"], ref)


# --- Band normalisation off --------------------------------------------------------


def test_bandwidth_is_subtracted_for_band_levels(result):
    park = result["park"]
    assert park.attrs["spl"]["bandwidth_subtracted"] is True
    for band, paper in ((8000, -32.68), (3150, -28.66)):
        term = (
            park[f"spl_{band}Hz"]
            - park[f"f_{band}Hz"]
            - park[f"beta_{band}Hz"]
        )
        np.testing.assert_allclose(
            term, -10 * np.log10(ap.third_octave_bandwidth(band))
        )
        # paper code used 10 log10 of the base-2 nominal bandwidth: 0.03-0.05 dB apart
        assert term.median() == pytest.approx(paper, abs=0.06)


# --- Three park depths --------------------------------------------------------------


def test_each_park_depth_uses_its_own_profile(result):
    # Each profile starts at its cycle's park depth; it serves the PAL
    # measurements it reaches (deepest point + 100 dbar).  Deeper
    # excursions use the next profile that reaches them.
    park, prof = result["park"], F.read_profiles()
    assert result["ts_profiles"].attrs["cycles"] == [41, 42, 45]
    reach = prof.groupby("cycle")["PRES"].max() + 100
    for cycle in (41, 42, 45):
        rows = park[park["cycle"] == cycle]
        own = rows[rows["pressure"] <= reach[cycle]]
        assert len(own) > 0.8 * len(rows)  # cycle 41: 17 deeper rows
        assert (own["beta_cycle_after"] == cycle).all()
        deeper = rows[rows["pressure"] > reach[cycle]]
        assert (deeper["beta_cycle_after"] != cycle).all()
        assert rows["beta_8000Hz"].notna().all()


def test_beta_grows_with_park_depth(result):
    med = (
        result["park"]
        .groupby("cycle")[["beta_3150Hz", "beta_8000Hz"]]
        .median()
    )
    for band in ("beta_3150Hz", "beta_8000Hz"):
        assert med.loc[41, band] < med.loc[42, band] < med.loc[45, band]
    assert (med["beta_8000Hz"] > med["beta_3150Hz"]).all()


def test_beta_at_1000_m_matches_corrected_published_value(result):
    # After the unit fix, β at 1000 m on this float's Ligurian winter
    # profile is ≈ 0.96 dB at 8 kHz and ≈ 0.26 dB at 3.15 kHz
    # (Xiaogang; Delaigue), computed from a different (mean) profile and
    # with the original MATLAB sound speed.
    b = result["beta"]
    b45 = b[b["cycle"] == 45].set_index("depth")
    at = b45.index[np.argmin(np.abs(b45.index - 1000))]
    assert b45.loc[at, "beta_8000Hz"] == pytest.approx(0.96, abs=0.03)
    assert b45.loc[at, "beta_3150Hz"] == pytest.approx(0.26, abs=0.03)


# --- Park data, positions, wind -------------------------------------------------


def test_park_ctd_interpolated_on_every_cycle(result):
    park, samples = result["park"], result["park_samples"]
    for cycle, rows in park.groupby("cycle"):
        src = samples.loc[samples["cycle"] == cycle, "TEMP"]
        assert rows["TEMP"].notna().mean() > 0.8
        assert rows["TEMP"].dropna().between(src.min(), src.max()).all()


def test_positions_between_gps_fixes(result):
    park, fixes = result["park"], result["positions"]
    assert park["latitude"].notna().all()  # pattern-0 fix precedes cycle 41
    assert (
        park["latitude"]
        .between(fixes["latitude"].min(), fixes["latitude"].max())
        .all()
    )
    assert (
        park["longitude"]
        .between(fixes["longitude"].min(), fixes["longitude"].max())
        .all()
    )


def test_literature_winds(result):
    park = result["park"]
    assert set(park.attrs["wind_models"]["applied"]) == {
        "vagle1990",
        "nystuen2015",
        "pensieri2015",
    }
    for col in ("U_vagle1990", "U_nystuen2015", "U_pensieri2015"):
        assert park[col].dropna().between(0, 40).all()
    outside = ~park["spl_8000Hz"].between(30, 60, inclusive="left")
    assert park.loc[outside, "U_pensieri2015"].isna().all()


def test_figures(result):
    figs = ap.plot_all(result)
    assert {"ais", "wind", "beta", "profiles"} <= set(figs)


# --- Real AIS ----------------------------------------------------------------------


def test_ais_times_are_receiver_times_not_file_dates():
    # The positions come from files named ..0215, ..0216, ..0217, but each
    # daily file mostly holds the day before: times must follow the receiver.
    ais = ap.read_ais(F.AIS_CSV)
    assert ais["datetime"].min() >= pd.Timestamp("2025-02-14 12:30")
    assert ais["datetime"].max() <= pd.Timestamp("2025-02-16 02:15")
    assert (ais["datetime"].dt.day == 14).any()


def test_ais_csv_and_raw_nmea_agree():
    pytest.importorskip("pyais")
    # inside the area cut around the track (outside it, the slices keep
    # one position per minute only, for coverage)
    lat0, lon0, km = 43.36, 7.77, 24.0
    box = (
        lon0 - km / (111.2 * np.cos(np.radians(lat0))),
        lon0 + km / (111.2 * np.cos(np.radians(lat0))),
        lat0 - km / 111.2,
        lat0 + km / 111.2,
    )
    csv = ap.read_ais(F.AIS_CSV, bbox=box)
    nmea = ap.read_ais(F.AIS_NMEA, bbox=box)
    assert len(csv) > 10_000
    key = ["datetime", "mmsi"]
    a = csv.assign(
        latitude=csv.latitude.round(4), longitude=csv.longitude.round(4)
    )
    b = nmea.assign(
        latitude=nmea.latitude.round(4), longitude=nmea.longitude.round(4)
    )
    merged = a.merge(
        b, on=[*key, "latitude", "longitude"], how="outer", indicator=True
    )
    # every decoded CSV position is in the raw NMEA; the raw decoding finds
    # a handful more that the 2025 parser missed
    assert (merged["_merge"] != "left_only").all()
    assert (merged["_merge"] == "right_only").sum() <= 0.001 * len(merged)


def test_base_station_in_raw_nmea_is_not_a_ship():
    pyais = pytest.importorskip("pyais")
    from pyais.stream import IterMessages

    with gzip.open(F.AIS_NMEA, "rb") as f:
        lines = [line.rstrip(b"\r\n") for line in f if line.strip()]
    base = {
        m.decode().mmsi
        for m in IterMessages(lines)
        if m.decode().msg_type == 4
    }
    assert len(base) == 1
    assert not base & set(ap.read_ais(F.AIS_NMEA)["mmsi"])
    assert pyais  # used


def _brute_force_count(t, lat, lon, ais, radius_km=20, minutes=30):
    """Independent, slow count of ships near one measurement."""
    seen = set()
    for row in ais.itertuples():
        if abs((row.datetime - t).total_seconds()) > minutes * 60:
            continue
        p1, p2 = math.radians(lat), math.radians(row.latitude)
        dphi = p2 - p1
        dl = math.radians(row.longitude - lon)
        h = (
            math.sin(dphi / 2) ** 2
            + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
        )
        if 2 * 6371.0088 * math.asin(math.sqrt(h)) <= radius_km:
            seen.add(row.mmsi)
    return len(seen)


def test_ship_counts_match_brute_force(result):
    park, ais = result["park"], result["ais"]
    rows = park[park["cycle"] == 45].iloc[::40]
    assert len(rows) >= 10
    for r in rows.itertuples():
        expected = _brute_force_count(r.datetime, r.latitude, r.longitude, ais)
        assert r.ships_nearby == expected


def test_ais_coverage_only_where_ais_exists(result):
    park = result["park"]
    # AIS slice covers cycle 45 only: other cycles are unknown, not "no ship"
    assert park.loc[park["cycle"].isin([41, 42]), "ships_nearby"].isna().all()
    c45 = park.loc[park["cycle"] == 45, "ships_nearby"]
    assert c45.notna().all()
    assert 0 < park.loc[park["cycle"] == 45, "ship_flag"].mean() < 1


def test_ais_from_csv_or_nmea_gives_same_flags():
    pytest.importorskip("pyais")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = F.process(ais=F.AIS_CSV)["park"]
        b = F.process(ais=F.AIS_NMEA)["park"]
    np.testing.assert_array_equal(
        a["ships_nearby"].to_numpy(), b["ships_nearby"].to_numpy()
    )


def test_ais_remove_blanks_flagged_values():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        park = F.process(ais=F.AIS_CSV, ais_action="remove")["park"]
    flagged = park[park["ship_flag"]]
    assert len(flagged) > 0
    assert (
        flagged[["f_8000Hz", "spl_8000Hz", "U_nystuen2015"]].isna().all().all()
    )
    kept = park[~park["ship_flag"] & park["ships_nearby"].notna()]
    assert kept["spl_8000Hz"].notna().all()
