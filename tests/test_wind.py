# %%
"""Literature wind models (Delaigue et al. 2026, Eqs. 2-4)."""

import warnings

import numpy as np
import pandas as pd
import pytest

import argopal


SPL = np.linspace(25, 70, 91)


# Reference implementations: the code of Delaigue et al. (2026)
# (github.com/.../WIND-FROM-FLOAT-ACOUSTICS, "Code for analysis and
# figures.ipynb", section 5), copied verbatim, to check the transcription.
def _paper_nystuen(SPL):
    return 0.0005 * SPL**3 - 0.0310 * SPL**2 + 0.4904 * SPL + 2.0871


def _paper_pensieri(SPL):
    SPL = np.asarray(SPL)
    result = np.full_like(SPL, np.nan, dtype=np.float64)
    mask1 = (SPL >= 30) & (SPL < 38)
    mask2 = (SPL >= 38) & (SPL < 60)
    result[mask1] = 0.1458 * SPL[mask1] - 3.146
    result[mask2] = 0.044642 * SPL[mask2] ** 2 - 3.2917 * SPL[mask2] + 63.016
    return result


def _paper_vagle(SPL):
    A, B, C = 7.38, -38.70, SPL - 21.69
    disc = B**2 - 4 * A * C
    result = np.full_like(SPL, np.nan, dtype=np.float64)
    mask = disc >= 0
    result[mask] = 10 ** ((-B - np.sqrt(disc[mask])) / (2 * A))
    return result


def test_nystuen_matches_paper_code():
    np.testing.assert_allclose(
        argopal.wind_nystuen2015(SPL), _paper_nystuen(SPL)
    )


def test_pensieri_matches_paper_code():
    np.testing.assert_allclose(
        argopal.wind_pensieri2015(SPL), _paper_pensieri(SPL), equal_nan=True
    )


def test_pensieri_regime_change_at_38_db():
    # "38 dB corresponding to a wind speed of 2.39 m/s" (Sect. 2.6)
    lo = argopal.wind_pensieri2015([37.999999])[0]
    hi = argopal.wind_pensieri2015([38.0])[0]
    assert lo == pytest.approx(2.39, abs=0.01)
    assert hi == pytest.approx(2.39, abs=0.01)
    assert np.isnan(argopal.wind_pensieri2015([29.9, 60.0])).all()
    assert np.isfinite(
        argopal.wind_pensieri2015([25.0], valid_range=None)
    ).all()


def test_vagle_matches_paper_code():
    np.testing.assert_allclose(
        argopal.wind_vagle1990(SPL), _paper_vagle(SPL), equal_nan=True
    )


def test_vagle_inverts_its_relation():
    u = np.array([1.0, 5.0, 10.0, 20.0])
    x = np.log10(u)
    spl = 21.69 + 38.70 * x - 7.38 * x**2
    np.testing.assert_allclose(argopal.wind_vagle1990(spl), u)
    assert np.isnan(argopal.wind_vagle1990([80.0])[0])  # no real root


def _pal(bands, spl=50.0, n=5):
    df = pd.DataFrame({f"spl_{b}Hz": np.full(n, spl) for b in bands})
    df.attrs = {"bands_hz": list(bands)}
    return df


def test_wind_literature_applies_available_models_only():
    out = argopal.wind_literature(_pal([8000]))
    for m in ("vagle1990", "nystuen2015", "pensieri2015"):
        assert f"U_{m}" in out
    info = out.attrs["wind_models"]
    assert set(info["applied"]) == {"vagle1990", "nystuen2015", "pensieri2015"}
    assert info["skipped"] == {}


def test_wind_literature_explicit_model_without_band_warns():
    with pytest.warns(UserWarning, match="nystuen2015"):
        out = argopal.wind_literature(_pal([5000]), models=["nystuen2015"])
    assert not [c for c in out if c.startswith("U_")]
    with pytest.raises(ValueError, match="Unknown"):
        argopal.wind_literature(_pal([8000]), models=["smith2020"])


def test_pipeline_adds_wind(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = argopal.process_float(7902280, mode="raw")
    park = r["park"]
    np.testing.assert_allclose(
        park["U_nystuen2015"], argopal.wind_nystuen2015(park["spl_8000Hz"])
    )
    assert set(argopal.LITERATURE_MODELS) == {
        "vagle1990",
        "nystuen2015",
        "pensieri2015",
    }
    figs = argopal.plot_all(r)
    assert "wind" in figs
