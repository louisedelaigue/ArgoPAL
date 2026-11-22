# %%
"""Depth correction β, absorption, spectral levels and positions.

Reference values: the MATLAB depth correction of Delaigue et al. (2026)
(beta_h.m), run on their Ligurian Sea mean profile.  The two files in
tests/data/matlab_ligurian_* are its input sound-speed/T/S profile and its
published β table; they are used only here, to check the Python port.
"""

import os

import numpy as np
import pandas as pd
import pytest

import argopal as ap
from argopal.acoustics import _rays_on_range_grid, ray_tables
from argopal.positions import _slerp


DATA = os.path.join(os.path.dirname(__file__), "data")
PROFILE = pd.read_csv(
    os.path.join(DATA, "matlab_ligurian_sound_speed_profile.csv")
)
PUBLISHED = pd.read_csv(
    os.path.join(DATA, "matlab_ligurian_beta_published.csv")
)
for _c in ("Beta_dB_f3150", "Beta_dB_f8000"):
    PUBLISHED[_c] = PUBLISHED[_c].map(
        lambda x: complex(str(x).replace("i", "j")).real
    )


RECEIVERS = PUBLISHED["Depth_m"].to_numpy(float)


def _beta(units, receivers=RECEIVERS):
    return ap.beta_profile(
        PROFILE["Depth_m"].to_numpy(float),
        PROFILE["SoundSpeed_mps"].to_numpy(float),
        PROFILE["Temperature_degC"].to_numpy(float),
        PROFILE["Salinity_PSU"].to_numpy(float),
        [3150, 8000],
        receiver_depths=receivers,
        absorption_units=units,
    )


# --- β: port of beta_h.m ----------------------------------------------------


def test_beta_reproduces_published_matlab_with_its_unit_bug():
    # absorption_units="Np" reproduces exp(-alpha L) of beta_h.m.  From
    # 400 m down MATLAB's values are real and match to 1e-3 dB; above,
    # MATLAB's path length went complex (it took the real part), so only
    # approximate agreement is expected there.
    b = _beta("Np")
    deep = PUBLISHED["Depth_m"] >= 400
    np.testing.assert_allclose(
        b.loc[deep, "beta_8000Hz"],
        PUBLISHED.loc[deep, "Beta_dB_f8000"],
        atol=1e-3,
    )
    np.testing.assert_allclose(
        b.loc[deep, "beta_3150Hz"],
        PUBLISHED.loc[deep, "Beta_dB_f3150"],
        atol=1e-3,
    )
    np.testing.assert_allclose(
        b["beta_8000Hz"], PUBLISHED["Beta_dB_f8000"], atol=0.01
    )


def test_beta_with_unit_fix():
    # α in dB/m -> loss 10**(-αL/10).  Values reported after the fix
    # (Xiaogang; Delaigue): ≈ 0.96 dB at 8 kHz, ≈ 0.26 dB at 3.15 kHz,
    # at 1000 m on this profile.
    b = _beta("dB").set_index("depth")
    assert b.loc[1000, "beta_8000Hz"] == pytest.approx(0.96, abs=0.01)
    assert b.loc[1000, "beta_3150Hz"] == pytest.approx(0.26, abs=0.01)


def test_beta_default_is_the_fixed_version():
    default = ap.beta_profile(
        PROFILE["Depth_m"].to_numpy(float),
        PROFILE["SoundSpeed_mps"].to_numpy(float),
        PROFILE["Temperature_degC"].to_numpy(float),
        PROFILE["Salinity_PSU"].to_numpy(float),
        [8000],
        receiver_depths=[1000.0],
    )
    assert default.attrs["absorption_units"] == "dB"
    assert default["beta_8000Hz"].iloc[0] == pytest.approx(0.96, abs=0.01)


def test_buggy_version_overestimates_by_about_4():
    fixed = _beta("dB").set_index("depth").loc[1000, "beta_8000Hz"]
    buggy = _beta("Np").set_index("depth").loc[1000, "beta_8000Hz"]
    assert 3.5 < buggy / fixed < 5


def test_beta_grows_with_depth_and_frequency():
    b = _beta("dB")
    assert (np.diff(b["beta_8000Hz"]) > 0).all()
    assert (b["beta_8000Hz"] >= b["beta_3150Hz"]).all()


def test_ray_inversion_matches_direct_root_finding():
    # beta_h.m solves r(θ) = r for each range with fzero; argopal inverts
    # r(θ) on a dense θ grid.  Check against bisection at 500 m.
    c = PROFILE["SoundSpeed_mps"].to_numpy(float)
    ih = 499  # receiver at 500 m (layers 1..500)
    theta, r_theta, l_theta = ray_tables(c, 1.0, [ih])
    ranges = np.array([1.0, 50.0, 500.0, 2000.0, 7999.0])
    th_fast, _ = _rays_on_range_grid(theta, r_theta[0], l_theta[0], ranges)

    def r_of(th):
        x = np.minimum(c[: ih + 1, None] / c[0] * np.cos(th), 0.999999)
        return (x / np.sqrt(1 - x * x)).sum(axis=0)

    lo = np.full(len(ranges), 1e-9)
    hi = np.full(len(ranges), np.pi / 2 - 1e-9)
    for _ in range(80):  # r falls as θ grows
        mid = (lo + hi) / 2
        too_far = r_of(mid) > ranges
        lo = np.where(too_far, mid, lo)
        hi = np.where(too_far, hi, mid)
    np.testing.assert_allclose(th_fast, (lo + hi) / 2, rtol=1e-4)


def test_beta_rejects_bad_grids():
    d = np.arange(1, 11.0)
    c = np.full(10, 1500.0)
    t = np.full(10, 13.0)
    s = np.full(10, 38.0)
    with pytest.raises(ValueError, match="evenly spaced"):
        ap.beta_profile(d + 0.5, c, t, s, [8000], [5.0])
    with pytest.raises(ValueError, match="missing"):
        ap.beta_profile(d, np.where(d > 5, np.nan, c), t, s, [8000], [5.0])
    with pytest.raises(ValueError, match="absorption_units"):
        ap.beta_profile(d, c, t, s, [8000], [5.0], absorption_units="neper")


# --- Absorption --------------------------------------------------------------


def test_absorption_increases_with_frequency():
    f = [1, 3.15, 8, 20]
    a = [ap.absorption_francois_garrison(x, 1000, 13.5, 38.6) for x in f]
    assert np.all(np.diff(a) > 0)


def test_absorption_order_of_magnitude():
    # Francois-Garrison at 8 kHz in warm Mediterranean water: ~0.5-1 dB/km
    a = ap.absorption_francois_garrison(8, 1000, 13.5, 38.6)
    assert 0.3e-3 < a < 1.5e-3


# --- Mean profile, depth correction, spectral level -------------------------


def _profiles(cycles, bottom=1000.0, top=2.0):
    rows = []
    for c in cycles:
        p = np.arange(top, bottom + 1, 2.0)
        rows.append(
            pd.DataFrame(
                {
                    "cycle": c,
                    "datetime": pd.Timestamp("2026-01-01")
                    + pd.Timedelta(days=5 * c),
                    "latitude": 43.5,
                    "longitude": 7.5,
                    "PRES": p,
                    "PRES_qc": "1",
                    "TEMP": 14 - 0.001 * p + 0.1 * c,
                    "TEMP_qc": "1",
                    "PSAL": 38.5 + 0.0001 * p,
                    "PSAL_qc": "1",
                }
            )
        )
    df = pd.concat(rows, ignore_index=True)
    df.attrs = {"mode": "delayed", "wmo": 0}
    return df


def test_mean_profile_averages_profiles():
    mp = ap.mean_profile(_profiles([1, 3]), max_pressure=1000)
    row = mp.set_index("pressure").loc[500]
    assert row["temperature"] == pytest.approx(14 - 0.5 + 0.2)
    assert mp.attrs["cycles"] == [1, 3]
    assert mp["pressure"].iloc[0] == 0  # extrapolated to the surface


def test_profile_held_constant_below_deepest_point():
    grid = ap.gridded_profiles(_profiles([1]), max_pressure=1020)
    deep = grid[grid["pressure"] > 1000]
    at_1000 = grid.set_index("pressure").loc[1000, "temperature"]
    np.testing.assert_allclose(deep["temperature"], at_1000)
    assert not deep["measured"].any()


def test_surface_held_constant_by_default():
    # Shallowest point 8 dbar: constant above it, no overshoot
    grid = ap.gridded_profiles(_profiles([1], top=8.0), max_pressure=1000)
    g = grid.set_index("pressure")
    assert g.loc[0, "temperature"] == g.loc[8, "temperature"]
    lin = ap.gridded_profiles(
        _profiles([1], top=8.0), max_pressure=1000, surface="linear"
    ).set_index("pressure")
    assert lin.loc[0, "temperature"] > lin.loc[8, "temperature"]


def test_profiles_skipped_with_reasons():
    prof = _profiles([1, 2, 3, 4])
    prof.loc[prof["cycle"] == 2, "TEMP_qc"] = "4"  # all bad
    prof = prof[~((prof["cycle"] == 3) & (prof["PRES"] < 20))]  # no surface
    prof = prof[~((prof["cycle"] == 4) & (prof["PRES"] > 80))]  # too short
    prof.attrs = {"mode": "delayed"}
    grid = ap.gridded_profiles(prof, max_pressure=1000)
    assert grid.attrs["cycles"] == [1]
    assert set(grid.attrs["skipped"]) == {2, 3, 4}
    with pytest.raises(ValueError, match="No profile"):
        ap.mean_profile(_profiles([1], bottom=60), max_pressure=1000)


def test_short_gap_at_bottom_is_accepted():
    # 975 dbar profile, PAL to 1028 dbar: held constant, not dropped
    grid = ap.gridded_profiles(_profiles([1], bottom=975), max_pressure=1028)
    assert grid.attrs["cycles"] == [1]


def _two_profile_betas():
    """β of two synthetic profiles 10 days apart: one warm, one cold."""
    prof = _profiles([1, 3])
    prof.loc[prof["cycle"] == 3, "TEMP"] += 8  # strong warming near surface
    grid = ap.gridded_profiles(prof, max_pressure=1000)
    return ap.compute_profile_betas(grid, [8000])


def test_profile_betas_one_per_profile():
    b = _two_profile_betas()
    assert sorted(b["cycle"].unique()) == [1, 3]
    assert {"datetime", "depth", "beta_8000Hz"} <= set(b.columns)
    assert b.attrs["method"] == "interpolate"


def test_depth_correct_interpolates_beta_in_time():
    b = _two_profile_betas()
    t1, t3 = sorted(b["datetime"].unique())
    times = [
        t1 - pd.Timedelta("1D"),
        t1,
        t1 + (t3 - t1) / 4,
        t3,
        t3 + pd.Timedelta("1D"),
    ]
    pal = pd.DataFrame(
        {"datetime": times, "pressure": 1000.0, "f_8000Hz": 45.0}
    )
    pal.attrs = {"bands_hz": [8000]}
    out = ap.depth_correct(pal, b)
    b1 = out["beta_8000Hz"].iloc[1]
    b3 = out["beta_8000Hz"].iloc[3]
    assert b1 != pytest.approx(b3)
    # before the first / after the last profile: that profile's β
    assert out["beta_8000Hz"].iloc[0] == pytest.approx(b1)
    assert out["beta_8000Hz"].iloc[4] == pytest.approx(b3)
    # a quarter of the way: linear in time
    assert out["beta_8000Hz"].iloc[2] == pytest.approx(0.75 * b1 + 0.25 * b3)
    assert list(out["beta_cycle_before"]) == [1, 1, 1, 3, 3]
    assert list(out["beta_cycle_after"]) == [1, 1, 3, 3, 3]


def test_depth_correct_no_beta_far_from_profiles():
    # e.g. delayed mode: profiles only for the first months of a float
    b = _two_profile_betas()
    t3 = b["datetime"].max()
    pal = pd.DataFrame(
        {
            "datetime": [t3 + pd.Timedelta("5D"), t3 + pd.Timedelta("60D")],
            "pressure": 1000.0,
            "f_8000Hz": 45.0,
        }
    )
    pal.attrs = {"bands_hz": [8000]}
    with pytest.warns(UserWarning, match="no usable profile"):
        out = ap.depth_correct(pal, b)
    assert np.isfinite(out["beta_8000Hz"].iloc[0])
    assert np.isnan(out["beta_8000Hz"].iloc[1])


def test_mean_beta_still_available():
    mp = ap.mean_profile(_profiles([1, 3]), max_pressure=1000)
    b = ap.compute_beta(mp, [8000])
    assert b.attrs["method"] == "mean"
    pal = pd.DataFrame(
        {
            "datetime": [pd.Timestamp("2026-01-08")],
            "pressure": [500.0],
            "f_8000Hz": [45.0],
        }
    )
    pal.attrs = {"bands_hz": [8000]}
    out = ap.depth_correct(pal, b)
    assert np.isfinite(out["beta_8000Hz"].iloc[0])
    assert "beta_cycle_before" not in out


def _pal(pressures):
    df = pd.DataFrame({"pressure": pressures, "f_8000Hz": 45.0})
    df.attrs = {"bands_hz": [8000], "band_normalization": True}
    return df


def test_depth_correct_interpolates_beta():
    beta = pd.DataFrame(
        {"depth": [10.0, 20.0, 1000.0], "beta_8000Hz": [0.1, 0.2, 1.0]}
    )
    beta.attrs = {"latitude": 43.5}
    out = ap.depth_correct(_pal([0.0, 504.0]), beta)
    # 0 dbar: linear extrapolation from the two shallowest β
    assert out["beta_8000Hz"].iloc[0] == pytest.approx(0.0, abs=1e-9)
    # ~500 m: between 20 m (0.2) and 1000 m (1.0)
    d = out["depth"].iloc[1]
    assert out["beta_8000Hz"].iloc[1] == pytest.approx(
        0.2 + 0.8 * (d - 20) / 980
    )


def test_depth_correct_warns_below_beta():
    beta = pd.DataFrame({"depth": [10.0, 20.0], "beta_8000Hz": [0.1, 0.2]})
    beta.attrs = {"latitude": 43.5}
    with pytest.warns(UserWarning, match="deeper than the deepest"):
        out = ap.depth_correct(_pal([500.0]), beta)
    assert np.isnan(out["beta_8000Hz"].iloc[0])


def test_third_octave_bandwidth():
    # 8 kHz band: exact centre 7943 Hz, width 1833 Hz -> 32.63 dB
    assert ap.third_octave_bandwidth(8000) == pytest.approx(1833, abs=1)
    assert 10 * np.log10(ap.third_octave_bandwidth(8000)) == pytest.approx(
        32.63, abs=0.01
    )


def test_spectral_level_with_and_without_band_normalisation():
    pal = _pal([1000.0])
    pal["beta_8000Hz"] = 1.0
    on = ap.spectral_level(pal)  # attrs: normalised on board
    assert on["spl_8000Hz"].iloc[0] == pytest.approx(46.0)
    off = ap.spectral_level(pal, band_normalization=False)
    assert off["spl_8000Hz"].iloc[0] == pytest.approx(46.0 - 32.63, abs=0.01)


def test_spectral_level_needs_normalisation_and_beta():
    pal = _pal([1000.0])
    pal["beta_8000Hz"] = 1.0
    pal.attrs.pop("band_normalization")
    with pytest.raises(ValueError, match="normalised"):
        ap.spectral_level(pal)
    with pytest.raises(KeyError, match="depth_correct"):
        ap.spectral_level(_pal([1000.0]), band_normalization=True)


# --- Positions ----------------------------------------------------------------


def test_slerp_midpoint_on_equator():
    lat, lon = _slerp(
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([90.0]),
        0.5,
    )
    assert lat[0] == pytest.approx(0, abs=1e-9)
    assert lon[0] == pytest.approx(45)


def test_attach_positions_interpolates_without_extrapolation():
    fixes = pd.DataFrame(
        {
            "datetime": pd.to_datetime(
                ["2026-01-01 00:00", "2026-01-01 10:00"]
            ),
            "latitude": [43.0, 44.0],
            "longitude": [7.0, 7.0],
        }
    )
    pal = pd.DataFrame(
        {
            "datetime": pd.to_datetime(
                ["2025-12-31 23:00", "2026-01-01 05:00", "2026-01-01 11:00"]
            )
        }
    )
    out = ap.attach_positions(pal, fixes)
    assert np.isnan(out["latitude"].iloc[0]) and np.isnan(
        out["latitude"].iloc[2]
    )
    assert out["latitude"].iloc[1] == pytest.approx(43.5, abs=1e-6)
    assert out["position_gap_h"].iloc[1] == pytest.approx(10)


def test_read_positions_includes_launch():
    pos = ap.read_positions(os.path.join(DATA, "7902280_Rtraj_cycles1-2.nc"))
    assert pos["datetime"].is_monotonic_increasing
    assert (pos["measurement_code"] == 0).sum() == 1  # launch
    assert set(pos["position_qc"]) <= {"1", "2"}


def test_shallow_profiles_serve_only_shallow_parks():
    # Like the 2025 float: each profile starts at that cycle's park depth.
    # Cycle 1 reaches 500 dbar, cycle 3 reaches 1000 dbar.
    prof = pd.concat(
        [
            _profiles([1], bottom=500),
            _profiles([3], bottom=1000),
        ],
        ignore_index=True,
    )
    prof.loc[prof["cycle"] == 3, "TEMP"] += 8  # warm: different β
    prof.attrs = {"mode": "raw"}
    grid = ap.gridded_profiles(prof, max_pressure=1000)
    assert grid.attrs["cycles"] == [1, 3]  # the 500 dbar profile is kept
    assert grid.groupby("cycle")["pressure"].max().to_dict() == {
        1: 600,
        3: 1000,
    }
    b = ap.compute_profile_betas(grid, [8000])
    t1, t3 = sorted(b["datetime"].unique())
    mid = t1 + (t3 - t1) / 2
    pal = pd.DataFrame(
        {"datetime": [mid, mid], "pressure": [450.0, 950.0], "f_8000Hz": 45.0}
    )
    pal.attrs = {"bands_hz": [8000]}
    out = ap.depth_correct(pal, b)
    # 450 dbar: between the two profiles; 950 dbar: only cycle 3 reaches it
    assert list(out["beta_cycle_before"]) == [1, 3]
    assert list(out["beta_cycle_after"]) == [3, 3]
    assert out["beta_weight_after"].iloc[0] == pytest.approx(0.5)
    assert out["beta_8000Hz"].notna().all()
