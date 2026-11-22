# %%
"""Band frequencies are assigned by argopal (the GDAC files do not store
them), so these tests make sure the assignment fails loudly instead of
silently giving wrong frequencies if a float sends a different layout."""

import os
import warnings

import numpy as np
import pytest
import xarray as xr

import argopal


DATA = os.path.join(os.path.dirname(__file__), "data")
TRAJ = os.path.join(DATA, "7902280_Rtraj_aux_first200.nc")
SPEC = "PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM"


def _write_spectrum(tmp_path, spectrum, dim):
    """Copy of the fixture with its spectrum replaced, saved to tmp_path."""
    keep = [
        "PLATFORM_NUMBER",
        "CYCLE_NUMBER",
        "CYCLE_NUMBER_ADJUSTED",
        "MEASUREMENT_CODE",
        "JULD",
        "JULD_ADJUSTED",
        "PRES",
        "PRES_ADJUSTED",
        SPEC + "_QC",
    ]
    ds = xr.open_dataset(TRAJ)[keep].assign(
        {SPEC: (("N_MEASUREMENT", dim), np.asarray(spectrum, dtype="f4"))}
    )
    ds.encoding.pop("unlimited_dims", None)  # N_HISTORY was dropped
    fname = str(tmp_path / "Rtraj_aux_modified.nc")
    ds.to_netcdf(fname)
    return fname


def _fixture_spectrum():
    ds = xr.open_dataset(TRAJ)
    return np.asarray(ds[SPEC].values, dtype=float)


# --- The assigned frequency list itself ------------------------------------


def _nine_band_spectrum():
    """The fixture's 23-band spectra reduced to the 9-band layout, as a
    9-band float would send them."""
    idx = [argopal.PAL_BANDS_HZ.index(f) for f in argopal.PAL_BANDS_9_HZ]
    return _fixture_spectrum()[:, idx]


def test_only_two_layouts():
    assert sorted(argopal.BAND_LAYOUTS) == [9, 23]
    assert argopal.BAND_LAYOUTS[23] == argopal.PAL_BANDS_HZ
    assert argopal.BAND_LAYOUTS[9] == argopal.PAL_BANDS_9_HZ


def test_23_band_list():
    assert argopal.PAL_BANDS_HZ == (
        63, 100, 125, 160, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
        2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000,
        25000,
    )  # fmt: skip


def test_9_band_list():
    assert argopal.PAL_BANDS_9_HZ == (
        63, 125, 400, 1000, 2000, 5000, 8000, 12500, 20000,
    )  # fmt: skip


def test_9_bands_are_a_subset_of_23():
    assert set(argopal.PAL_BANDS_9_HZ) <= set(argopal.PAL_BANDS_HZ)


@pytest.mark.parametrize("n", [23, 9])
def test_band_list_is_strictly_increasing(n):
    f = np.asarray(argopal.BAND_LAYOUTS[n])
    assert (np.diff(f) > 0).all()


def _third_octave_steps(bands):
    # Base-10 third-octave centres are 1000 * 10**(k/10) for integer k;
    # nominal values (63, 125, 12500 ...) are within 3% of the exact ones.
    f = np.asarray(bands, dtype=float)
    k = np.round(10 * np.log10(f / 1000))
    exact = 1000 * 10 ** (k / 10)
    assert (np.abs(f / exact - 1) < 0.03).all()
    return np.diff(k).astype(int).tolist()


def test_23_bands_are_standard_third_octave_centres():
    # One third-octave apart, except the known gap 160 -> 400 Hz
    # (200, 250, 315 Hz not sent): 4 steps.  63 -> 100 Hz is 2 steps
    # (80 Hz not sent).
    assert _third_octave_steps(argopal.PAL_BANDS_HZ) == [2, 1, 1, 4] + [1] * 18


def test_9_bands_are_standard_third_octave_centres():
    assert _third_octave_steps(argopal.PAL_BANDS_9_HZ) == [
        3,
        5,
        4,
        3,
        4,
        2,
        2,
        2,
    ]


# --- The number of values sent by the float --------------------------------


@pytest.mark.parametrize("n", [23, 9])
def test_bands_for_known_layouts(n):
    assert argopal.bands_for(n) == argopal.BAND_LAYOUTS[n]


@pytest.mark.parametrize("n", [0, 1, 8, 10, 22, 24, 31])
def test_bands_for_any_other_count_is_refused(n):
    with pytest.raises(ValueError, match=f"{n} values"):
        argopal.bands_for(n)


def test_read_pal_accepts_23_bands():
    pal = argopal.read_pal(TRAJ)
    assert pal.attrs["bands_hz"] == list(argopal.PAL_BANDS_HZ)
    assert "f_25000Hz" in pal.columns


def test_read_pal_accepts_9_bands(tmp_path):
    spectrum = _nine_band_spectrum()
    fname = _write_spectrum(tmp_path, spectrum, "N_VALUES9")
    with warnings.catch_warnings():
        warnings.simplefilter("error", argopal.BandOrderWarning)
        pal = argopal.read_pal(fname)
    assert pal.attrs["bands_hz"] == list(argopal.PAL_BANDS_9_HZ)
    bands = [c for c in pal.columns if c.startswith("f_")]
    assert bands == [f"f_{f}Hz" for f in argopal.PAL_BANDS_9_HZ]
    # Values land in the right columns: 8 kHz of the 9-band file is the
    # 8 kHz of the original 23-band spectrum
    full = argopal.read_pal(TRAJ)
    np.testing.assert_array_equal(pal["f_8000Hz"], full["f_8000Hz"])
    np.testing.assert_array_equal(pal["f_63Hz"], full["f_63Hz"])
    assert "f_100Hz" not in pal.columns


@pytest.mark.parametrize("n", [8, 10, 22, 24])
def test_read_pal_refuses_any_other_count(tmp_path, n):
    s = _fixture_spectrum()
    spectrum = np.hstack([s, s])[:, :n]
    fname = _write_spectrum(tmp_path, spectrum, f"N_VALUES{n}")
    with pytest.raises(ValueError, match=f"{n} values"):
        argopal.read_pal(fname)


def test_read_pal_refuses_dimension_name_mismatch(tmp_path):
    # 23 values, but the dimension claims 9: the file is not what we expect
    fname = _write_spectrum(tmp_path, _fixture_spectrum(), "N_VALUES9")
    with pytest.raises(ValueError, match="N_VALUES9"):
        argopal.read_pal(fname)


# --- The order of the values -----------------------------------------------


def test_real_data_falls_with_frequency():
    with warnings.catch_warnings():
        warnings.simplefilter("error", argopal.BandOrderWarning)
        pal = argopal.read_pal(TRAJ)
    assert pal.attrs["band_order_score"] < -0.95


@pytest.mark.parametrize("n", [23, 9])
def test_reversed_bands_are_flagged(tmp_path, n):
    s = _fixture_spectrum() if n == 23 else _nine_band_spectrum()
    fname = _write_spectrum(tmp_path, s[:, ::-1], f"N_VALUES{n}")
    with pytest.warns(argopal.BandOrderWarning, match="wrong order"):
        pal = argopal.read_pal(fname)
    assert pal.attrs["band_order_score"] > 0.9


@pytest.mark.parametrize("n", [23, 9])
@pytest.mark.parametrize("seed", range(5))
def test_shuffled_bands_are_flagged(tmp_path, n, seed):
    s = _fixture_spectrum() if n == 23 else _nine_band_spectrum()
    order = np.random.default_rng(seed).permutation(n)
    fname = _write_spectrum(tmp_path, s[:, order], f"N_VALUES{n}")
    with pytest.warns(argopal.BandOrderWarning):
        argopal.read_pal(fname)


def test_order_check_skips_too_few_rows():
    pal = argopal.read_pal(TRAJ).head(5)
    assert np.isnan(argopal.band_order_score(pal))
    with warnings.catch_warnings():
        warnings.simplefilter("error", argopal.BandOrderWarning)
        argopal.check_band_order(pal)


def test_select_all_bands_with_data():
    pal = argopal.read_pal(TRAJ)
    out = argopal.select_bands(pal, bands="all")
    assert out.attrs["bands_hz"] == argopal.available_bands(pal)
    assert out.attrs["bands_dropped"] == []
    with pytest.raises(ValueError, match='"all"'):
        argopal.select_bands(pal, bands="everything")
