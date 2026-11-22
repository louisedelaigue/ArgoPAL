# %%
"""Uncertainty of the refitted wind models (wind_uncertainty)."""

import warnings

import numpy as np
import pandas as pd
import pytest

import argopal

from . import float2025 as F


TRUE = {"a3": 0.0004, "a2": -0.02, "a1": 0.3, "a0": 1.0}


def _float(n_cycles=20, per_cycle=60, noise=0.5, relative=False, seed=0):
    """Synthetic float: SPL in 35-60 dB, reference wind from Nystuen with
    known coefficients plus noise (fixed, or proportional to the wind)."""
    rng = np.random.default_rng(seed)
    n = n_cycles * per_cycle
    spl = rng.uniform(35, 60, n)
    u_true = argopal.wind_nystuen2015(spl, **TRUE)
    sd = noise * u_true if relative else noise
    pal = pd.DataFrame(
        {
            "cycle": np.repeat(np.arange(1, n_cycles + 1), per_cycle),
            "datetime": pd.date_range("2026-01-01", periods=n, freq="5min"),
            "spl_8000Hz": spl,
            "U_nystuen2015": argopal.wind_nystuen2015(spl),
            "U_buoy": u_true + rng.normal(0, 1, n) * sd,
            "buoy_distance_km": 5.0,
        }
    )
    pal.attrs = {
        "bands_hz": [8000],
        "wind_models": {
            "applied": {"nystuen2015": {"band_hz": 8000, "reference": "N"}}
        },
        "references": {"buoy": {"column": "U_buoy"}},
    }
    out, _ = argopal.fit_wind_models(pal)
    return out


def test_four_columns_and_bins():
    p, bins = argopal.wind_uncertainty(_float())
    for suffix in ("_cv", "_coef", "_ref", ""):
        assert f"U_nystuen2015_fit_sigma{suffix}" in p
    assert set(bins.columns) == {"model", "wind", "sigma_cv", "n", "lo", "hi"}
    assert (bins["n"] >= 30).all()


def test_empirical_sigma_recovers_the_noise():
    p, bins = argopal.wind_uncertainty(_float(noise=0.5))
    assert bins["sigma_cv"].to_numpy() == pytest.approx(0.5, abs=0.12)
    assert p["U_nystuen2015_fit_sigma_cv"].median() == pytest.approx(
        0.5, abs=0.1
    )


def test_empirical_sigma_grows_with_wind_when_noise_does():
    _, bins = argopal.wind_uncertainty(_float(noise=0.1, relative=True))
    assert bins["sigma_cv"].iloc[-1] > 2 * bins["sigma_cv"].iloc[0]


def test_coefficient_sigma_small_inside_large_outside_fitted_range():
    p, _ = argopal.wind_uncertainty(_float())
    inside = p["U_nystuen2015_fit_sigma_coef"].median()
    assert inside < 0.2
    # extrapolation: apply the bootstrap refits to a level far outside
    q = p.copy()
    q.loc[0, "spl_8000Hz"] = 75.0
    q.loc[0, "U_nystuen2015_fit"] = argopal.wind_nystuen2015(
        75.0, **q.attrs["wind_fits"]["coefficients"]["nystuen2015"]
    )
    q, _ = argopal.wind_uncertainty(q)
    assert q.loc[0, "U_nystuen2015_fit_sigma_coef"] > 5 * inside


def test_reference_error_is_removed_in_quadrature():
    p, _ = argopal.wind_uncertainty(_float(), reference_sigma=0.3)
    cv = p["U_nystuen2015_fit_sigma_cv"]
    coef = p["U_nystuen2015_fit_sigma_coef"]
    assert (p["U_nystuen2015_fit_sigma_ref"] == 0.3).all()
    np.testing.assert_allclose(
        p["U_nystuen2015_fit_sigma"], np.sqrt(cv**2 - 0.3**2 + coef**2)
    )
    # a reference error larger than the empirical one leaves only coef
    big, _ = argopal.wind_uncertainty(_float(), reference_sigma=5.0)
    np.testing.assert_allclose(
        big["U_nystuen2015_fit_sigma"], big["U_nystuen2015_fit_sigma_coef"]
    )


def test_needs_fit_cycles_and_enough_of_them():
    p = _float()
    p.attrs.pop("wind_fits")
    with pytest.raises(KeyError, match="fit_wind_models"):
        argopal.wind_uncertainty(p)
    with pytest.raises(KeyError, match="cycle"):
        argopal.wind_uncertainty(_float().drop(columns="cycle"))
    few = _float(n_cycles=2)
    with pytest.warns(UserWarning, match="at least 3"):
        q, _ = argopal.wind_uncertainty(few)
    assert q["U_nystuen2015_fit_sigma"].isna().all()
    with pytest.warns(UserWarning, match="unreliable"):
        argopal.wind_uncertainty(_float(n_cycles=5))


def test_2025_float_three_cycles():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = F.process()["park"]
        p = argopal.collocate_reference(
            p, F.DATA + "/dyfamed_wind.csv", name="dyfamed"
        )
        p, _ = argopal.fit_wind_models(p)
    with pytest.warns(UserWarning, match="unreliable"):
        q, bins = argopal.wind_uncertainty(p, reference_sigma=0.5)
    nys = bins[bins["model"] == "nystuen2015"]
    assert 0.5 < nys["sigma_cv"].min() and nys["sigma_cv"].max() < 4
    # three cycles cannot pin the coefficients down: a heavy tail of
    # large bootstrap spreads (several m/s) at the edges of the record
    assert q["U_nystuen2015_fit_sigma_coef"].quantile(0.95) > 2


def test_pipeline_adds_uncertainty(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = argopal.process_float(7902280, mode="raw")
        park = base["park"]
        t = pd.date_range(
            park["datetime"].min().floor("h"),
            park["datetime"].max().ceil("h"),
            freq="1h",
        )
        ref = pd.DataFrame(
            {
                "datetime": t,
                "wind_speed": np.linspace(3, 12, len(t)),
                "latitude": park["latitude"].iloc[0],
                "longitude": park["longitude"].iloc[0],
            }
        )
        r = argopal.process_float(
            7902280, mode="raw", reference={"buoy": ref}, n_bootstrap=10
        )
    assert "wind_uncertainty" in r
    assert "U_nystuen2015_fit_sigma" in r["park"]
