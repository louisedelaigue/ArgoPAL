# %%
"""Wind footprints from a reanalysis, footprint matching and selection,
and the per-reference diagnostics."""

import warnings

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import argopal

from .test_reference import _ten_minute_float


T0 = pd.Timestamp("2026-02-01")
TIMES = pd.date_range(T0, periods=300, freq="1h")
A, B = (43.38, 5.0), (43.38, 9.0)  # two buoys, 4° of longitude apart


def _reanalysis(names=("time", "latitude", "longitude"), uv=True, seed=0):
    """Wind that follows buoy A west of 7.5°E and buoy B east of it."""
    rng = np.random.default_rng(seed)
    lat = np.arange(40, 46.01, 0.25)
    lon = np.arange(4, 10.01, 0.25)
    sa = 8 + np.cumsum(rng.normal(0, 0.5, len(TIMES)))
    sb = 8 + np.cumsum(rng.normal(0, 0.5, len(TIMES)))
    west = (lon < 7.5)[None, None, :]
    field = np.where(west, sa[:, None, None], sb[:, None, None])
    field = field + rng.normal(0, 0.05, (len(TIMES), len(lat), len(lon)))
    t, y, x = names
    coords = {t: TIMES, y: lat, x: lon}
    if uv:  # as ERA5: components, wind from the south-west
        u = xr.DataArray(field / np.sqrt(2), coords, dims=names)
        return xr.Dataset({"u10": u, "v10": u})
    return xr.DataArray(field, coords, dims=names).to_dataset(name="si10")


def _buoys():
    rows = []
    for name, (la, lo) in (("A", A), ("B", B)):
        rows.append(
            pd.DataFrame(
                {
                    "reference_name": name,
                    "datetime": TIMES,
                    "latitude": la,
                    "longitude": lo,
                    "wind_speed": 5.0 if name == "A" else 9.0,
                    "height_m": 10.0,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def test_read_reanalysis_names_and_speed():
    ws = argopal.read_reanalysis(_reanalysis())
    assert ws.dims == ("time", "latitude", "longitude")
    assert ws.name == "wind_speed"
    ws2 = argopal.read_reanalysis(
        _reanalysis(names=("valid_time", "lat", "lon"), uv=False)
    )
    assert ws2.dims == ("time", "latitude", "longitude")
    np.testing.assert_allclose(ws.values, ws2.values)
    with pytest.raises(KeyError, match="wind speed"):
        argopal.read_reanalysis(xr.Dataset({"t2m": ws}))


def test_footprints():
    fp = argopal.wind_footprints(_reanalysis(), {"A": A, "B": B})
    assert fp.dims == ("reference", "latitude", "longitude")
    west = fp.sel(longitude=6.0, latitude=43.5)
    east = fp.sel(longitude=8.5, latitude=43.5)
    assert west.sel(reference="A") > 0.99 and east.sel(reference="B") > 0.99
    assert west.sel(reference="B") < 0.9 and east.sel(reference="A") < 0.9
    # from a reference table: each reference at its position
    fp2 = argopal.wind_footprints(_reanalysis(), _buoys())
    assert list(fp2["reference"].values) == ["A", "B"]
    # R at given positions; NaN outside the grid
    r = argopal.footprint_r_at(fp, [43.5, 43.5, 60.0], [6.0, 8.5, 6.0])
    assert r["A"].iloc[0] > 0.99 and r["B"].iloc[1] > 0.99
    assert r.iloc[2].isna().all()


def test_monthly_anomaly_footprints():
    times = pd.date_range("2000-01-01", periods=48, freq="MS")
    rng = np.random.default_rng(1)
    s = rng.normal(0, 1, len(times))
    field = np.broadcast_to(s[:, None, None], (48, 3, 3)).copy()
    da = xr.DataArray(
        field,
        {
            "time": times,
            "latitude": [42.0, 43.0, 44.0],
            "longitude": [5, 6, 7],
        },
        dims=("time", "latitude", "longitude"),
    )
    fp = argopal.wind_footprints(da, {"X": (43, 6)}, anomaly="monthly")
    assert float(fp.min()) > 0.99
    with pytest.raises(ValueError, match="anomaly"):
        argopal.wind_footprints(da, {"X": (43, 6)}, anomaly="daily")
    # hourly data over a few months: monthly anomalies are not defined
    with pytest.raises(ValueError, match="at least twice"):
        argopal.wind_footprints(_reanalysis(), {"A": A}, anomaly="monthly")


def _float_between():
    """A float at 7.2°E: closer to B, but in A's footprint."""
    pal = pd.DataFrame(
        {
            "datetime": TIMES[10:40] + pd.Timedelta("5min"),
            "latitude": 43.38,
            "longitude": 7.2,
            "cycle": np.repeat([1, 2, 3], 10),
            "spl_8000Hz": np.linspace(40, 55, 30),
        }
    )
    pal.attrs = {"bands_hz": [8000]}
    return argopal.wind_literature(pal)


def test_matching_by_footprint_instead_of_distance():
    pal = _float_between()
    closest = argopal.collocate_reference(pal, _buoys(), name="b")
    assert set(closest["b_reference_name"]) == {"B"}
    fp = argopal.wind_footprints(_reanalysis(), _buoys())
    by_fp = argopal.collocate_reference(pal, _buoys(), name="b", footprints=fp)
    assert set(by_fp["b_reference_name"]) == {"A"}
    assert (by_fp["b_footprint_r"] > 0.99).all()
    assert (by_fp["U_b"] == 5.0).all()
    info = by_fp.attrs["references"]["b"]
    assert info["matching"] == "footprint"
    # fitting by footprint: within it, however far
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, fits = argopal.fit_wind_models(
            by_fp, max_distance_km=None, min_footprint_r=0.9
        )
        _, none = argopal.fit_wind_models(by_fp, max_distance_km=40)
    assert fits["n"].iloc[0] > 0
    assert none["n"].iloc[0] == 0  # 180 km from A: nothing within 40 km
    with pytest.raises(KeyError, match="footprints"):
        argopal.fit_wind_models(closest, min_footprint_r=0.9)
    # a footprint missing for a reference
    with pytest.raises(ValueError, match="No footprint"):
        argopal.collocate_reference(
            pal, _buoys(), name="c", footprints=fp.sel(reference=["A"])
        )


def test_reference_map_with_footprints():
    pal = _float_between()
    fp = argopal.wind_footprints(_reanalysis(), _buoys())
    out = argopal.collocate_reference(pal, _buoys(), name="b", footprints=fp)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, _ = argopal.fit_wind_models(
            out, max_distance_km=None, min_footprint_r=0.9
        )
    _, ax = argopal.plot_reference_map(out, footprints=fp)
    labels = [ln.get_label() for ln in ax.lines]
    assert any("inside footprint" in x for x in labels)
    assert [t.get_text() for t in ax.texts if t.get_text() in "AB"] == ["A"]


def test_reference_diagnostics_and_plot():
    pal = _ten_minute_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted, fits = argopal.fit_wind_models(pal)
        d = argopal.reference_diagnostics(fitted)
    t = d["fits"].set_index("reference")
    assert list(t.index) == ["all", "buoy"]
    assert t.at["all", "n"] == fits["n"].iloc[1]
    # the pooled fit is the refit itself
    assert t.at["all", "rmse_pooled"] == pytest.approx(
        fits.set_index("model").at["nystuen2015", "rmse_fitted"]
    )
    assert np.isfinite(t.at["buoy", "rmse_own_cv"])
    assert list(d["transfer"].index) == ["pooled", "buoy"]
    assert {"spl", "published", "pooled", "buoy"} <= set(d["curves"])
    _, ax = argopal.plot_reference_fits(fitted, annotate=True)
    labels = [ln.get_label() for ln in ax.lines]
    assert "All references (refit)" in labels


def test_pipeline_with_footprints(offline):
    from .test_reference import _reference_for

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = argopal.process_float(7902280, mode="raw")
        park = base["park"]
        ref = _reference_for(park)
        t = pd.date_range(
            park["datetime"].min().floor("h"),
            park["datetime"].max().ceil("h"),
            freq="1h",
        )
        rng = np.random.default_rng(2)
        s = 8 + np.cumsum(rng.normal(0, 0.5, len(t)))
        lat, lon = np.arange(41, 46.01, 0.5), np.arange(4, 12.01, 0.5)
        field = s[:, None, None] + rng.normal(0, 0.1, (len(t), 11, 17))
        rea = xr.DataArray(
            field, {"time": t, "latitude": lat, "longitude": lon},
            dims=("time", "latitude", "longitude"),
        )  # fmt: skip
        r = argopal.process_float(
            7902280, mode="raw", reference=ref, footprints=rea,
            n_bootstrap=5,
        )  # fmt: skip
    assert r["footprints"].name == "footprint_r"
    assert r["park"].attrs["wind_fits"]["min_footprint_r"] == 0.9
    assert r["park"].attrs["wind_fits"]["max_distance_km"] is None
    assert r["wind_fits"]["n"].iloc[0] > 0
    assert "reference_map" in argopal.plot_all(r)
