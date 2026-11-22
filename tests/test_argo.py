# %%
"""Reading core Argo files with a chosen data mode."""

import os

import numpy as np
import pytest

import argopal as ap


DATA = os.path.join(os.path.dirname(__file__), "data")
SPROF = os.path.join(
    DATA, "7902280_Sprof_4profiles.nc"
)  # cycles 1,2 D; 38,39 A
RTRAJ = os.path.join(DATA, "7902280_Rtraj_cycles1-2.nc")


# --- select_mode on hand-made values ---------------------------------------

#            R      A      D     D (value removed by operator, QC 4)
RAW = [10.0, 20.0, 30.0, 40.0]
ADJ = [np.nan, 21.0, 31.0, np.nan]
RAW_QC = ["1", "1", "1", "1"]
ADJ_QC = ["", "1", "2", "4"]
DM = ["R", "A", "D", "D"]


def _select(mode):
    return ap.select_mode(RAW, ADJ, RAW_QC, ADJ_QC, DM, mode)


def test_default_mode_is_delayed():
    assert ap.DEFAULT_MODE == "delayed"


def test_delayed_uses_only_D_and_never_falls_back():
    v, qc, dm = _select("delayed")
    np.testing.assert_array_equal(v, [np.nan, np.nan, 31.0, np.nan])
    assert list(dm) == ["", "", "D", "D"]
    assert list(qc) == ["", "", "2", "4"]  # removed value keeps its flag


def test_adjusted_uses_only_A():
    # Real-time adjusted only: delayed-mode values are not mixed in
    v, _, dm = _select("adjusted")
    np.testing.assert_array_equal(v, [np.nan, 21.0, np.nan, np.nan])
    assert list(dm) == ["", "A", "", ""]


def test_raw_uses_raw_everywhere():
    v, qc, dm = _select("raw")
    np.testing.assert_array_equal(v, RAW)
    assert list(dm) == ["R"] * 4
    assert list(qc) == RAW_QC


def test_no_best_mode():
    assert ap.MODES == ("delayed", "adjusted", "raw")
    with pytest.raises(ValueError):
        _select("best")


def test_unknown_mode():
    with pytest.raises(ValueError, match="mode must be one of"):
        _select("DELAYED")


# --- Profiles ---------------------------------------------------------------


def test_profiles_default_is_delayed():
    prof = ap.read_profiles(SPROF)
    assert prof.attrs["mode"] == "delayed"
    # Only the two delayed-mode profiles survive
    assert sorted(prof["cycle"].unique()) == [1, 2]
    assert set(prof["TEMP_data_mode"]) <= {"D", ""}
    # Oxygen (A) and pH (R) have no delayed-mode values on this float
    assert prof["DOXY"].isna().all()
    assert prof["PH_IN_SITU_TOTAL"].isna().all()


def test_profiles_adjusted_keeps_only_A():
    prof = ap.read_profiles(SPROF, mode="adjusted")
    # T/S of cycles 1-2 are delayed mode, so only cycles 38-39 have T/S;
    # oxygen is in A mode in all four profiles
    assert sorted(prof.dropna(subset=["TEMP"])["cycle"].unique()) == [38, 39]
    assert set(prof["TEMP_data_mode"]) <= {"A", ""}
    assert prof["DOXY"].notna().any()
    assert prof["PH_IN_SITU_TOTAL"].isna().all()


def test_profiles_raw_has_everything():
    prof = ap.read_profiles(SPROF, mode="raw")
    for var in ["TEMP", "PSAL", "DOXY", "PH_IN_SITU_TOTAL"]:
        assert prof[var].notna().any()
        assert set(prof[var + "_data_mode"]) <= {"R", ""}


def test_profiles_columns_and_variable_selection():
    prof = ap.read_profiles(SPROF, mode="raw", variables=["TEMP", "NOPE"])
    assert prof.attrs["variables"] == ["TEMP"]
    for col in [
        "cycle",
        "direction",
        "datetime",
        "latitude",
        "longitude",
        "PRES",
        "TEMP",
        "TEMP_qc",
        "TEMP_data_mode",
    ]:
        assert col in prof.columns
    assert "PSAL" not in prof.columns


# --- Park drift (Rtraj) -----------------------------------------------------


def test_park_only_park_drift():
    park = ap.read_park(RTRAJ, mode="raw")
    assert set(park["measurement_code"]) == {290}
    assert park["datetime"].is_monotonic_increasing
    # T/S and O2 are measured at different times: rows hold one or other
    assert park["TEMP"].notna().any() and park["DOXY"].notna().any()


def test_park_default_delayed_is_empty_for_this_float():
    # 7902280's trajectory has no delayed-mode data at all
    park = ap.read_park(RTRAJ)
    for var in park.attrs["variables"]:
        assert park[var].isna().all()


# --- Inventory --------------------------------------------------------------


def test_available_variables():
    inv = ap.available_variables(
        {"7902280_Sprof.nc": SPROF, "7902280_Rtraj.nc": RTRAJ}
    )
    sprof = inv[inv["file"] == "7902280_Sprof.nc"].set_index("variable")
    assert sprof.loc["TEMP", "n_D"] == 2 and sprof.loc["TEMP", "n_A"] == 2
    assert sprof.loc["DOXY", "n_A"] == 4
    assert sprof.loc["PH_IN_SITU_TOTAL", "n_R"] == 4
    traj = inv[inv["file"] == "7902280_Rtraj.nc"]
    assert "DOXY" in set(traj["variable"])
    assert (traj["counted_per"] == "measurement").all()
