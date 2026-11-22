# %%
"""Figures for each processing step, and the band-pattern check."""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import argopal as ap


DATA = os.path.join(os.path.dirname(__file__), "data")
PAL = os.path.join(DATA, "7902280_Rtraj_aux_first200.nc")
STEPS = [
    "spectra",
    "despike",
    "park",
    "positions",
    "profiles",
    "beta",
    "depth_correction",
    "wind",
]


@pytest.fixture
def result(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ap.process_float(7902280, mode="raw", verbose=False)


def test_plot_all_draws_every_step(result, tmp_path):
    figs = ap.plot_all(result)
    assert list(figs) == STEPS
    for name, fig in figs.items():
        fname = tmp_path / f"{name}.png"
        fig.savefig(fname)
        assert fname.stat().st_size > 10_000  # something was drawn
        plt.close(fig)


@pytest.mark.parametrize(
    "func, args",
    [
        (ap.plot_pal_spectra, lambda r: (r["park"],)),
        (ap.plot_despike, lambda r: (r["pal_raw"], r["park"])),
        (ap.plot_park, lambda r: (r["park"], r["park_samples"])),
        (ap.plot_positions, lambda r: (r["park"], r["positions"])),
        (ap.plot_profiles, lambda r: (r["ts_profiles"],)),
        (ap.plot_beta, lambda r: (r["beta"],)),
        (ap.plot_depth_correction, lambda r: (r["park"],)),
    ],
)
def test_each_plot_returns_figure_and_axes(result, func, args):
    fig, axes = func(*args(result))
    assert isinstance(fig, plt.Figure)
    axes = axes if isinstance(axes, (list, tuple)) else [axes]
    assert all(a.has_data() or not a.axison for a in axes)
    plt.close(fig)


def test_plots_work_on_partial_results(offline):
    # delayed mode on the fixture: no park values, β only for some rows
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ap.process_float(7902280, verbose=False)
    for fig in ap.plot_all(r).values():
        plt.close(fig)


def test_plot_mean_beta(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ap.process_float(7902280, mode="raw", beta_method="mean")
    fig, axes = ap.plot_beta(r["beta"])
    assert len(axes) == 1
    plt.close(fig)


def test_plot_unknown_band(result):
    with pytest.raises(KeyError, match="No 3000 Hz band"):
        ap.plot_despike(result["pal_raw"], result["park"], band=3000)


def test_despike_plot_picks_a_cycle(result):
    fig, (ax1, _) = ap.plot_despike(result["pal_raw"], result["park"], cycle=2)
    assert "cycle 2" in ax1.get_ylabel()
    plt.close(fig)


# --- Band pattern check -------------------------------------------------------


def _pal_with_pattern(filled_bands):
    pal = ap.read_pal(PAL)
    for f in ap.PAL_BANDS_HZ:
        if f not in filled_bands:
            pal[f"f_{f}Hz"] = np.nan
    return pal


def test_band_layout_full_spectra():
    assert ap.read_pal(PAL).attrs["band_layout_counts"] == {
        "23": 200,
        "9": 0,
        "other": 0,
    }


def test_nine_band_spectra_in_23_array_are_recognised():
    # As on 7902280: 9-band spectra stored with the 14 other bands empty
    pal = _pal_with_pattern(ap.PAL_BANDS_9_HZ)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ap.BandLayoutWarning)
        counts = ap.check_band_layout(pal)
    assert counts == {"23": 0, "9": 200, "other": 0}


def test_unexpected_band_pattern_warns():
    # e.g. a shifted 9-band set: 100 Hz instead of 63 Hz filled
    shifted = (100,) + ap.PAL_BANDS_9_HZ[1:]
    pal = _pal_with_pattern(shifted)
    with pytest.warns(ap.BandLayoutWarning, match="200 spectra"):
        counts = ap.check_band_layout(pal)
    assert counts["other"] == 200


def test_band_layout_survives_mixed_rows():
    pal = ap.read_pal(PAL)
    others = [f for f in ap.PAL_BANDS_HZ if f not in ap.PAL_BANDS_9_HZ]
    pal.loc[:49, [f"f_{f}Hz" for f in others]] = np.nan
    assert ap.band_layout_counts(pal) == {"23": 150, "9": 50, "other": 0}
    assert isinstance(pal, pd.DataFrame)


def test_plot_all_saves_at_300_dpi(result, tmp_path):
    from PIL import Image

    ap.plot_all(result, save_dir=str(tmp_path))
    with Image.open(tmp_path / "beta.png") as im:
        dpi = im.info["dpi"][0]
    assert round(dpi) == 300
    assert sorted(p.stem for p in tmp_path.glob("*.png")) == sorted(STEPS)
    plt.close("all")


def test_beta_time_series_at_1000_m_by_default(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ap.process_float(7902280, mode="raw")
    _, (_, ax2) = ap.plot_beta(r["beta"])
    assert "1000 m" in ax2.get_ylabel()
    _, (_, ax2) = ap.plot_beta(r["beta"], depth=500)
    assert "500 m" in ax2.get_ylabel()
    plt.close("all")


def test_despike_options_reach_the_pipeline(offline):
    r = ap.process_float(
        7902280, mode="raw", despike_window="1h", despike_percentile=90
    )
    assert r["park"].attrs["despike"] == "window=1h, percentile=90"
    default = ap.process_float(7902280, mode="raw")
    assert default["park"].attrs["despike"] == "window=3h, percentile=95"


def test_annotate_adds_title_and_notes(result):
    plain = ap.plot_profiles(result["ts_profiles"])[0]
    assert plain._suptitle is None
    noted = ap.plot_profiles(result["ts_profiles"], annotate=True)[0]
    assert "profiles used" in noted._suptitle.get_text()
    plt.close("all")


def test_panel_labels(result):
    fig, axes = ap.plot_beta(result["beta"])
    assert [t.get_text() for ax in axes for t in ax.texts] == ["a)", "b)"]
    plt.close(fig)
