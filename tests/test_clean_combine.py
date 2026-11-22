# %%
"""PAL cleaning, park interpolation, export and the offline pipeline."""

import os
import warnings

import numpy as np
import pandas as pd
import pytest

import argopal
from argopal.clean import _nanpercentile_axis1


DATA = os.path.join(os.path.dirname(__file__), "data")
PAL = os.path.join(DATA, "7902280_Rtraj_aux_first200.nc")
META = os.path.join(DATA, "7902280_meta_aux.nc")
RTRAJ = os.path.join(DATA, "7902280_Rtraj_cycles1-2.nc")
SPROF = os.path.join(DATA, "7902280_Sprof_4profiles.nc")


def _pal_table(times, values, cycle=1):
    """Minimal 23-band PAL table with the same value in every band."""
    df = pd.DataFrame(
        {
            "cycle": cycle,
            "measurement_code": 290,
            "datetime": pd.to_datetime(times),
            "pressure": 1000.0,
            "spectrum_qc": "0",
        }
    )
    for f in argopal.PAL_BANDS_HZ:
        df[f"f_{f}Hz"] = np.asarray(values, dtype=float)
    df.attrs = {"bands_hz": list(argopal.PAL_BANDS_HZ), "wmo": 0}
    return df


# --- Despike ----------------------------------------------------------------


def test_fast_percentile_matches_numpy():
    rng = np.random.default_rng(0)
    w = rng.normal(size=(50, 30, 4))
    w[rng.random(w.shape) < 0.3] = np.nan
    w[0, :, 0] = np.nan  # an all-NaN window
    with np.errstate(all="ignore"), pytest.warns(RuntimeWarning):
        expected = np.nanpercentile(w, 95, axis=1, keepdims=True)
    np.testing.assert_allclose(_nanpercentile_axis1(w, 95), expected)


def test_despike_removes_a_spike():
    times = pd.date_range("2026-01-01", periods=37, freq="5min")
    values = np.full(37, 50.0)
    values[18] = 80.0  # a passing ship
    out = argopal.despike(_pal_table(times, values))
    # The spike is above the window's 95th percentile everywhere, so every
    # despiked value is the quiet background
    np.testing.assert_allclose(out["f_8000Hz"], 50.0)


def test_despike_matches_paper_definition():
    # mean of values <= 95th percentile within +/- 1.5 h, as in the paper
    rng = np.random.default_rng(1)
    times = pd.date_range("2026-01-01", periods=60, freq="5min")
    values = rng.normal(50, 2, 60)
    out = argopal.despike(_pal_table(times, values))
    s = pd.Series(values, index=times)
    i = 30
    w = s[times[i] - pd.Timedelta("1.5h") : times[i] + pd.Timedelta("1.5h")]
    expected = w[w <= w.quantile(0.95)].mean()
    assert out["f_8000Hz"].iloc[i] == pytest.approx(expected)


def test_despike_does_not_mix_cycles():
    t1 = pd.date_range("2026-01-01 00:00", periods=5, freq="5min")
    t2 = pd.date_range("2026-01-01 00:25", periods=5, freq="5min")
    a = _pal_table(t1, np.full(5, 40.0), cycle=1)
    b = _pal_table(t2, np.full(5, 60.0), cycle=2)
    both = pd.concat([a, b], ignore_index=True)
    both.attrs = a.attrs
    out = argopal.despike(both)
    assert (out.loc[out["cycle"] == 1, "f_8000Hz"] == 40.0).all()
    assert (out.loc[out["cycle"] == 2, "f_8000Hz"] == 60.0).all()


def test_clean_pal_drops_codes_and_bad_qc():
    times = pd.date_range("2026-01-01", periods=6, freq="5min")
    pal = _pal_table(times, np.full(6, 50.0))
    pal.loc[0, "measurement_code"] = 590
    pal.loc[1, "spectrum_qc"] = "4"
    out = argopal.clean_pal(pal)
    assert len(out) == 4
    assert out.attrs["clean"]["removed_measurement_code"] == 1
    assert out.attrs["clean"]["removed_bad_qc"] == 1


# --- Linear interpolation of park data --------------------------------------


def _park(times, values, cycle=1, qc="1", dm="R"):
    df = pd.DataFrame(
        {
            "cycle": cycle,
            "measurement_code": 290,
            "datetime": pd.to_datetime(times),
            "TEMP": values,
            "TEMP_qc": qc,
            "TEMP_data_mode": dm,
        }
    )
    df.attrs = {"variables": ["TEMP"], "mode": "raw"}
    return df


def test_attach_park_linear_in_time():
    park = _park(["2026-01-01 00:00", "2026-01-01 01:00"], [13.0, 14.0])
    pal = _pal_table(
        ["2026-01-01 00:00", "2026-01-01 00:15", "2026-01-01 01:00"],
        [50, 50, 50],
    )
    out = argopal.attach_park(pal, park)
    np.testing.assert_allclose(out["TEMP"], [13.0, 13.25, 14.0])
    assert list(out["TEMP_data_mode"]) == ["R"] * 3


def test_attach_park_no_extrapolation_or_long_gaps():
    park = _park(["2026-01-01 00:00", "2026-01-01 03:00"], [13.0, 14.0])
    pal = _pal_table(
        ["2025-12-31 23:55", "2026-01-01 01:00", "2026-01-01 03:05"],
        [50, 50, 50],
    )
    out = argopal.attach_park(pal, park, max_gap="2h")
    assert out["TEMP"].isna().all()  # before, inside a 3 h gap, after
    out = argopal.attach_park(pal, park, max_gap="4h")
    assert out["TEMP"].iloc[1] == pytest.approx(13 + 1 / 3)


def test_attach_park_never_bridges_cycles():
    park = pd.concat(
        [
            _park(["2026-01-01 00:00"], [13.0], cycle=1),
            _park(["2026-01-01 00:30"], [14.0], cycle=2),
        ],
        ignore_index=True,
    )
    park.attrs = {"variables": ["TEMP"]}
    pal = _pal_table(["2026-01-01 00:15"], [50], cycle=1)
    assert np.isnan(argopal.attach_park(pal, park)["TEMP"].iloc[0])


def test_attach_park_skips_bad_qc():
    park = _park(
        ["2026-01-01 00:00", "2026-01-01 00:30", "2026-01-01 01:00"],
        [13.0, 99.0, 14.0],
    )
    park.loc[1, "TEMP_qc"] = "4"
    pal = _pal_table(["2026-01-01 00:30"], [50])
    assert argopal.attach_park(pal, park)["TEMP"].iloc[0] == pytest.approx(
        13.5
    )


def test_attach_park_on_real_data():
    pal = argopal.clean_pal(argopal.read_pal(PAL))
    park = argopal.read_park(RTRAJ, mode="raw")
    out = argopal.attach_park(pal, park)
    # Park T is hourly, so nearly every PAL row gets a value
    assert out["TEMP"].notna().mean() > 0.9
    # Values stay within the range of the park samples they come from
    src = park["TEMP"].dropna()
    assert out["TEMP"].dropna().between(src.min(), src.max()).all()


# --- Offline pipeline and export --------------------------------------------


@pytest.mark.parametrize("mode", argopal.MODES)
def test_process_float_offline(offline, mode):
    if mode == "adjusted":
        # The fixture's adjusted-mode profiles are from June (cycles 38-39)
        # but its PAL rows from January: no β rather than June's β.
        with pytest.warns(UserWarning, match="usable profile"):
            result = argopal.process_float(7902280, mode=mode)
        assert result["park"]["beta_8000Hz"].isna().all()
        assert result["park"]["spl_8000Hz"].isna().all()
        return
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # e.g. β not covering PAL depths
        result = argopal.process_float(7902280, mode=mode)
    park, prof = result["park"], result["profiles"]
    assert len(park) == 200
    for col in ["f_8000Hz", "TEMP_data_mode", "latitude", "longitude",
                "depth", "beta_8000Hz", "spl_8000Hz"]:  # fmt: skip
        assert col in park
    # Every step ran on every row
    assert park[["latitude", "beta_8000Hz", "spl_8000Hz"]].notna().all().all()
    assert park.attrs["park"]["mode"] == mode
    assert prof.attrs["mode"] == mode
    assert result["ts_profiles"].attrs["mode"] == mode
    assert result["beta"].attrs["method"] == "interpolate"
    if mode == "delayed":
        assert park["TEMP"].isna().all()
        assert set(prof["TEMP_data_mode"]) <= {"D", ""}
        # Depth correction from the two delayed-mode profiles only
        assert result["ts_profiles"].attrs["cycles"] == [1, 2]
        assert set(park["beta_cycle_before"]) <= {1, 2}
    # 7902280 normalises bands on board: spl = band level + β
    np.testing.assert_allclose(
        park["spl_8000Hz"], park["f_8000Hz"] + park["beta_8000Hz"]
    )


def test_process_float_tables(offline):
    result = argopal.process_float(7902280)
    assert set(result) == {
        "park",
        "pal_raw",
        "park_samples",
        "profiles",
        "ts_profiles",
        "beta",
        "positions",
        "meta",
        "files",
    }
    mean = argopal.process_float(7902280, beta_method="mean")
    assert "mean_profile" in mean
    assert mean["beta"].attrs["method"] == "mean"
    with pytest.raises(ValueError, match="beta_method"):
        argopal.process_float(7902280, beta_method="nearest")


# --- Export at any stage ----------------------------------------------------


@pytest.mark.parametrize("fmt", ["csv", "netcdf"])
def test_export_any_single_table(tmp_path, fmt):
    pal = argopal.read_pal(PAL)
    stages = {
        "read": pal,
        "clean": argopal.clean_pal(pal),
        "profiles": argopal.read_profiles(SPROF),
    }
    for name, table in stages.items():
        ext = ".csv" if fmt == "csv" else ".nc"
        (fname,) = argopal.export(
            table, str(tmp_path / (name + ext)), format=fmt
        )
        back = (
            pd.read_csv(fname) if fmt == "csv" else argopal.read_export(fname)
        )
        assert isinstance(back, pd.DataFrame)
        assert len(back) == len(table)
        assert list(back.columns) == list(table.columns)


def test_export_csv_adds_extension(tmp_path):
    (fname,) = argopal.export(argopal.read_pal(PAL), str(tmp_path / "pal"))
    assert fname.endswith("pal.csv")


def test_export_dict_csv(offline, tmp_path):
    result = argopal.process_float(7902280, mode="raw")
    files = argopal.export(result, str(tmp_path / "wmo7902280"), format="csv")
    assert sorted(os.path.basename(f) for f in files) == [
        "wmo7902280_beta.csv",
        "wmo7902280_pal_raw.csv",
        "wmo7902280_park.csv",
        "wmo7902280_park_samples.csv",
        "wmo7902280_positions.csv",
        "wmo7902280_profiles.csv",
        "wmo7902280_ts_profiles.csv",
    ]
    back = pd.read_csv(str(tmp_path / "wmo7902280_park.csv"))
    assert len(back) == len(result["park"])
    np.testing.assert_allclose(
        back["spl_8000Hz"], result["park"]["spl_8000Hz"]
    )


def test_export_netcdf_roundtrip(offline, tmp_path):
    result = argopal.process_float(7902280, mode="raw")
    fname = str(tmp_path / "wmo7902280.nc")
    argopal.export(result, fname, format="netcdf")
    back = argopal.read_export(fname)
    assert set(back) == {
        "park",
        "pal_raw",
        "park_samples",
        "profiles",
        "ts_profiles",
        "beta",
        "positions",
    }
    for name in back:
        assert len(back[name]) == len(result[name])
    for col in ["f_8000Hz", "spl_8000Hz", "TEMP", "latitude"]:
        np.testing.assert_allclose(
            back["park"][col], result["park"][col], equal_nan=True
        )
    assert (
        back["park"]["TEMP_data_mode"] == result["park"]["TEMP_data_mode"]
    ).all()
    assert back["park"].attrs["wmo"] == 7902280


def test_export_partial_dict(offline, tmp_path):
    # Any subset of tables, e.g. only what has been computed so far
    pal = argopal.clean_pal(argopal.read_pal(PAL))
    files = argopal.export(
        {"pal": pal, "profiles": argopal.read_profiles(SPROF)},
        str(tmp_path / "early.nc"),
        format="netcdf",
    )
    assert set(argopal.read_export(files[0])) == {"pal", "profiles"}


def test_export_bad_input(tmp_path):
    with pytest.raises(ValueError, match="format"):
        argopal.export(
            argopal.read_pal(PAL), str(tmp_path / "x"), format="xlsx"
        )
    with pytest.raises(TypeError):
        argopal.export({"meta": {}}, str(tmp_path / "x"))


# --- Band selection -----------------------------------------------------------


def test_select_bands_default_is_8khz():
    pal = argopal.select_bands(argopal.read_pal(PAL))
    assert pal.attrs["bands_hz"] == [8000]
    assert [c for c in pal if c.startswith("f_")] == ["f_8000Hz"]


def test_select_several_bands_low_to_high():
    pal = argopal.select_bands(argopal.read_pal(PAL), bands=(8000, 3150))
    assert pal.attrs["bands_hz"] == [3150, 8000]


def test_select_unknown_band():
    with pytest.raises(ValueError, match="3000"):
        argopal.select_bands(argopal.read_pal(PAL), bands=3000)


def test_band_without_data_is_dropped():
    pal = argopal.read_pal(PAL)
    pal["f_3150Hz"] = np.nan  # as in 9-band spectra
    with pytest.warns(
        argopal.BandCoverageWarning, match="3150 Hz band is not used"
    ):
        out = argopal.select_bands(pal, bands=(3150, 8000))
    assert out.attrs["bands_hz"] == [8000]
    assert out.attrs["bands_dropped"] == [3150]
    assert "f_3150Hz" not in out
    assert 3150 not in argopal.available_bands(pal)
    with (
        pytest.warns(argopal.BandCoverageWarning),
        pytest.raises(ValueError, match="None of the requested"),
    ):
        argopal.select_bands(pal, bands=3150)


def test_process_float_bands(offline):
    r = argopal.process_float(7902280, mode="raw", bands=(3150, 8000))
    park = r["park"]
    for b in (3150, 8000):
        for prefix in ("f_", "beta_", "spl_"):
            assert f"{prefix}{b}Hz" in park
    assert "f_1000Hz" not in park  # only the selected bands are processed
    assert "f_1000Hz" in r["pal_raw"]  # but the raw table keeps them all
    assert set(r["beta"].filter(like="beta_").columns) == {
        "beta_3150Hz",
        "beta_8000Hz",
    }
    default = argopal.process_float(7902280, mode="raw")
    assert default["park"].attrs["bands_hz"] == [8000]
