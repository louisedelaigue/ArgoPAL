# %%
"""Reference winds, model refitting, skill and the second ship criterion."""

import os
import warnings

import numpy as np
import pandas as pd
import pytest

import argopal as ap

from . import float2025 as F


DATA25 = os.path.join(os.path.dirname(__file__), "data", "float2025")
DYFAMED = os.path.join(DATA25, "dyfamed_wind.csv")
T0 = pd.Timestamp("2025-02-14 12:00")


def _pal(times, lat=43.38, lon=7.83, spl=50.0):
    df = pd.DataFrame(
        {
            "datetime": pd.to_datetime(times),
            "latitude": lat,
            "longitude": lon,
            "spl_8000Hz": spl,
        }
    )
    df.attrs = {"bands_hz": [8000]}
    return df


# --- Reading and collocating ------------------------------------------------------


def test_read_reference_recognises_columns():
    df = pd.DataFrame(
        {"Time": [T0], "ff": [7.5], "Lat": [43.38], "Lon": [7.83]}
    )
    ref = ap.read_reference(df)
    assert list(ref.columns) == [
        "datetime",
        "wind_speed",
        "latitude",
        "longitude",
    ]
    assert ref["wind_speed"].iloc[0] == 7.5


def test_read_reference_fixed_station_and_errors():
    df = pd.DataFrame({"datetime": [T0], "wind_speed": [5.0]})
    ref = ap.read_reference(df, latitude=43.38, longitude=7.83)
    assert ref["latitude"].iloc[0] == 43.38
    with pytest.raises(ValueError, match="latitude"):
        ap.read_reference(df)


def test_collocate_table_nearest_in_time_within_tolerance():
    ref = pd.DataFrame(
        {
            "datetime": [T0, T0 + pd.Timedelta("1h")],
            "wind_speed": [5.0, 9.0],
            "latitude": 43.38,
            "longitude": 7.83,
        }
    )
    pal = _pal(
        [
            T0 + pd.Timedelta("10min"),
            T0 + pd.Timedelta("50min"),
            T0 + pd.Timedelta("3h"),
        ]
    )
    out = ap.collocate_reference(pal, ref, name="buoy")
    np.testing.assert_allclose(out["U_buoy"], [5.0, 9.0, np.nan])
    np.testing.assert_allclose(out["buoy_dt_min"], [-10.0, 10.0, np.nan])
    assert out["buoy_distance_km"].iloc[0] == pytest.approx(0, abs=1e-6)


def test_collocate_max_distance():
    ref = pd.DataFrame(
        {
            "datetime": [T0],
            "wind_speed": [5.0],
            "latitude": 43.38,
            "longitude": 7.83,
        }
    )
    far = _pal([T0], lat=43.38 + 50 / 111.2)  # 50 km north
    out = ap.collocate_reference(far, ref, name="buoy", max_distance_km=40)
    assert np.isnan(out["U_buoy"].iloc[0])
    assert out["buoy_distance_km"].iloc[0] == pytest.approx(50, rel=0.01)


def test_collocate_refuses_name_clash_and_needs_positions():
    ref = pd.DataFrame(
        {
            "datetime": [T0],
            "wind_speed": [5.0],
            "latitude": 43.38,
            "longitude": 7.83,
        }
    )
    pal = _pal([T0])
    pal.attrs["wind_models"] = {"applied": {"nystuen2015": {}}}
    with pytest.raises(ValueError, match="already exists"):
        ap.collocate_reference(pal, ref, name="nystuen2015")
    with pytest.raises(KeyError, match="attach_positions"):
        ap.collocate_reference(pal.drop(columns="latitude"), ref)


def _several_references():
    """Two buoys: A at the float (43.38N 7.83E), B 30 km north; A reports
    hourly 12:00-16:00, B every 30 min 12:00-15:00; both give an error."""
    t_a = pd.date_range(T0, periods=5, freq="1h")
    t_b = pd.date_range(T0, periods=7, freq="30min")
    a = pd.DataFrame(
        {
            "reference_name": "A",
            "datetime": t_a,
            "latitude": 43.38,
            "longitude": 7.83,
            "wind_speed": 5.0,
            "wind_speed_sigma": 0.5,
        }
    )
    b = pd.DataFrame(
        {
            "reference_name": "B",
            "datetime": t_b,
            "latitude": 43.38 + 30 / 111.2,
            "longitude": 7.83,
            "wind_speed": 9.0,
            "wind_speed_sigma": 1.0,
        }
    )
    return pd.concat([a, b]).sample(frac=1, random_state=0)  # mixed rows


def test_read_reference_optional_columns():
    ref = ap.read_reference(_several_references())
    assert list(ref.columns) == [
        "datetime",
        "wind_speed",
        "latitude",
        "longitude",
        "reference_name",
        "wind_speed_sigma",
    ]
    alias = _several_references().rename(
        columns={"reference_name": "station", "wind_speed_sigma": "sigma"}
    )
    assert set(ap.read_reference(alias)["reference_name"]) == {"A", "B"}


def test_collocate_several_references_closest_in_distance():
    pal = _pal(
        [
            T0 + pd.Timedelta("10min"),  # A and B in time: A, closer
            T0 + pd.Timedelta("3h50min"),  # A (4h00) and B (3h00) ok: A
            T0 + pd.Timedelta("4h30min"),  # B ended at 3h00, A within 1h
            T0 + pd.Timedelta("6h"),  # neither
        ]
    )
    out = ap.collocate_reference(pal, _several_references(), name="buoys")
    np.testing.assert_allclose(out["U_buoys"], [5.0, 5.0, 5.0, np.nan])
    assert list(out["buoys_reference_name"].iloc[:3]) == ["A"] * 3
    assert pd.isna(out["buoys_reference_name"].iloc[3])
    np.testing.assert_allclose(out["buoys_sigma"], [0.5, 0.5, 0.5, np.nan])
    assert np.isnan(out["buoys_distance_km"].iloc[3])
    info = out.attrs["references"]["buoys"]
    assert info["references"] == {"A": 3, "B": 0}

    # far from A, nearer B: B is used, within its own time window only
    north = _pal(
        [T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("4h20min")],
        lat=43.38 + 25 / 111.2,
    )
    out = ap.collocate_reference(north, _several_references(), name="buoys")
    assert list(out["buoys_reference_name"]) == ["B", "A"]
    np.testing.assert_allclose(out["U_buoys"], [9.0, 5.0])
    np.testing.assert_allclose(out["buoys_distance_km"], [5, 25], rtol=0.02)
    np.testing.assert_allclose(out["buoys_dt_min"], [-10, -20])


def test_collocate_several_references_radius():
    north = _pal([T0], lat=43.38 + 50 / 111.2)  # 20 km from B, 50 from A
    out = ap.collocate_reference(
        north, _several_references(), name="buoys", max_distance_km=10
    )
    assert np.isnan(out["U_buoys"].iloc[0])
    assert out["buoys_reference_name"].iloc[0] == "B"
    assert out["buoys_distance_km"].iloc[0] == pytest.approx(20, rel=0.02)


def test_uncertainty_uses_table_sigma():
    pal, _ = _synthetic_float()
    pal["cycle"] = np.repeat(np.arange(1, 11), len(pal) // 10)
    # measurements matched to buoy A (sigma 0.3) or B (sigma 0.4)
    half = len(pal) // 2
    pal["buoy_reference_name"] = ["A"] * half + ["B"] * (len(pal) - half)
    pal["buoy_sigma"] = [0.3] * half + [0.4] * (len(pal) - half)
    out, _ = ap.fit_wind_models(pal, max_distance_km=None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        q, _ = ap.wind_uncertainty(out, n_bootstrap=5)
        given, _ = ap.wind_uncertainty(out, n_bootstrap=5, reference_sigma=1)
    rms = np.sqrt((0.3**2 + 0.4**2) / 2)
    assert q["U_nystuen2015_fit_sigma_ref"].dropna().unique() == (
        pytest.approx(rms)
    )
    assert (given["U_nystuen2015_fit_sigma_ref"].dropna() == 1).all()
    skill = ap.wind_skill(out, max_distance_km=None, by_reference=True)
    assert list(skill["reference_name"].unique()) == ["all", "A", "B"]


# --- Metrics --------------------------------------------------------------------------


def test_wind_metrics():
    ref = np.array([2.0, 4.0, 6.0, 8.0, np.nan])
    m = ap.wind_metrics(ref + 1, ref)
    assert m["n"] == 4
    assert m["bias"] == pytest.approx(1.0)
    assert m["rmse"] == pytest.approx(1.0)
    assert m["slope"] == pytest.approx(1.0)
    assert m["intercept"] == pytest.approx(1.0)
    assert m["pearson_r2"] == pytest.approx(1.0)
    assert m["r2"] == pytest.approx(1 - 4 / 20)  # biased: R² < r²
    perfect = ap.wind_metrics(ref, ref)
    assert perfect["rmse"] == 0 and perfect["r2"] == 1


# --- Refitting ---------------------------------------------------------------------------


def _synthetic_float(n=400, seed=0, true=None):
    """SPL and a reference wind made with known Nystuen coefficients."""
    rng = np.random.default_rng(seed)
    spl = rng.uniform(35, 60, n)
    true = true or {"a3": 0.0004, "a2": -0.02, "a1": 0.3, "a0": 1.0}
    u = ap.wind_nystuen2015(spl, **true) + rng.normal(0, 0.2, n)
    pal = _pal(pd.date_range(T0, periods=n, freq="5min"), spl=spl)
    pal["U_nystuen2015"] = ap.wind_nystuen2015(spl)
    pal["U_buoy"] = u
    pal["buoy_distance_km"] = np.where(np.arange(n) < n // 2, 10.0, 60.0)
    pal.attrs["wind_models"] = {
        "applied": {"nystuen2015": {"band_hz": 8000, "reference": "Nystuen"}}
    }
    pal.attrs["references"] = {"buoy": {"column": "U_buoy"}}
    return pal, true


def test_fit_recovers_known_coefficients():
    pal, true = _synthetic_float()
    out, fits = ap.fit_wind_models(pal)
    got = fits.set_index("model").loc["nystuen2015", "fitted"]
    pred_true = ap.wind_nystuen2015(pal["spl_8000Hz"], **true)
    np.testing.assert_allclose(out["U_nystuen2015_fit"], pred_true, atol=0.1)
    assert set(got) == {"a3", "a2", "a1", "a0"}
    row = fits.iloc[0]
    assert row["rmse_fitted"] < row["rmse_published"]
    assert row["rmse_fitted"] == pytest.approx(0.2, abs=0.05)  # the noise


def test_fit_uses_only_nearby_measurements():
    pal, _ = _synthetic_float()
    _, fits = ap.fit_wind_models(pal, max_distance_km=40)
    assert fits["n"].iloc[0] == 200
    _, all_ = ap.fit_wind_models(pal, max_distance_km=None)
    assert all_["n"].iloc[0] == 400
    use = np.zeros(len(pal), bool)
    use[:50] = True
    _, sub = ap.fit_wind_models(pal, use=use)
    assert sub["n"].iloc[0] == 50
    # nothing close enough: the model keeps its row, unfitted
    with pytest.warns(UserWarning, match="too few"):
        out, none = ap.fit_wind_models(pal, max_distance_km=5)
    assert none["n"].iloc[0] == 0 and none["fitted"].iloc[0] is None
    assert np.isnan(none["rmse_fitted"].iloc[0])
    assert "U_nystuen2015_fit" not in out


def test_wind_skill_table():
    pal, _ = _synthetic_float()
    out, _ = ap.fit_wind_models(pal)
    skill = ap.wind_skill(out)
    assert list(skill["wind"]) == ["U_nystuen2015", "U_nystuen2015_fit"]
    assert skill["rmse"].iloc[1] < skill["rmse"].iloc[0]


# --- Second ship criterion -------------------------------------------------------------


def test_ship_deviation_needs_both_criteria():
    pal, _ = _synthetic_float(n=200)
    pal["buoy_distance_km"] = 10.0
    pal["U_nystuen2015_fit"] = pal["U_buoy"]
    pal.loc[[0, 1], "U_nystuen2015_fit"] += 5.0  # large error
    pal["ship_flag"] = False
    pal.loc[[0, 2], "ship_flag"] = True  # 0: ship + error; 2: ship only
    out = ap.flag_ship_deviation(pal, reference="buoy")
    # 1 has a large error but no ship: not contaminated
    assert list(np.flatnonzero(out["ship_contaminated"])) == [0]
    assert out.attrs["ais"]["n_contaminated"] == 1
    rm = ap.flag_ship_deviation(pal, reference="buoy", action="remove")
    assert np.isnan(rm.loc[0, "spl_8000Hz"]) and np.isnan(
        rm.loc[0, "U_nystuen2015_fit"]
    )
    assert rm.loc[0, "U_buoy"] == pal.loc[0, "U_buoy"]  # reference kept
    with pytest.raises(KeyError, match="flag_ships"):
        ap.flag_ship_deviation(pal.drop(columns="ship_flag"), reference="buoy")


# --- Real data: the 2025 float with DYFAMED -----------------------------------


@pytest.fixture(scope="module")
def float2025():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = F.process(ais=F.AIS_CSV)["park"]
        p = ap.collocate_reference(p, DYFAMED, name="dyfamed")
        p, fits = ap.fit_wind_models(p, reference="dyfamed")
    return p, fits


def test_2025_collocation(float2025):
    p, _ = float2025
    assert p["U_dyfamed"].notna().mean() > 0.85  # hourly buoy, gaps
    assert p["dyfamed_distance_km"].max() < 10


def test_2025_published_coefficients_need_refitting(float2025):
    # The paper's finding: published models are biased on a float
    _, fits = float2025
    f = fits.set_index("model")
    assert f.loc["nystuen2015", "bias_published"] > 1.0
    for model in ("vagle1990", "nystuen2015", "pensieri2015"):
        assert f.loc[model, "rmse_fitted"] < f.loc[model, "rmse_published"]
        assert f.loc[model, "r2_fitted"] > 0.9
        assert abs(f.loc[model, "bias_fitted"]) < 0.05


def test_2025_second_ship_criterion(float2025):
    p, _ = float2025
    out = ap.flag_ship_deviation(p, reference="dyfamed")
    assert out["ship_contaminated"].sum() > 0
    assert (out["ship_contaminated"] <= out["ship_flag"]).all()  # a subset
    assert 0 < out.attrs["ais"]["deviation_rmse"] < 2


def test_2025_figures(float2025):
    p, _ = float2025
    _, axes = ap.plot_wind(p)
    assert len(axes) == 2  # published, refitted
    _, axes = ap.plot_wind_fit(p)
    assert len(axes) == 3


# --- Pipeline -----------------------------------------------------------------------------


def _reference_for(park):
    """A buoy next to the float, with wind close to the published Nystuen."""
    t = pd.date_range(
        park["datetime"].min().floor("h"),
        park["datetime"].max().ceil("h"),
        freq="1h",
    )
    spl = (
        park.set_index("datetime")["spl_8000Hz"]
        .resample("1h")
        .mean()
        .reindex(t)
        .interpolate()
    )
    return pd.DataFrame(
        {
            "datetime": t,
            "wind_speed": 0.8 * ap.wind_nystuen2015(spl.to_numpy()),
            "latitude": park["latitude"].iloc[0],
            "longitude": park["longitude"].iloc[0],
        }
    )


def test_pipeline_with_reference(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = ap.process_float(7902280, mode="raw")
        ref = _reference_for(base["park"])
        r = ap.process_float(7902280, mode="raw", reference={"buoy": ref})
    p = r["park"]
    assert {"U_buoy", "U_nystuen2015_fit"} <= set(p.columns)
    assert r["wind_fits"].attrs["reference"] == "buoy"
    assert "wind_skill" in r and "wind_fit" in ap.plot_all(r)
    # single source: named "ref"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        single = ap.process_float(7902280, mode="raw", reference=ref)
    assert "U_ref" in single["park"]


def test_pipeline_with_several_references_in_one_file(offline, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = ap.process_float(7902280, mode="raw")
    near = _reference_for(base["park"]).assign(
        reference_name="near", wind_speed_sigma=0.5
    )
    far = near.assign(
        reference_name="far",
        latitude=near["latitude"] + 1.0,
        wind_speed=near["wind_speed"] + 5,
    )
    path = tmp_path / "references.csv"
    pd.concat([far, near]).to_csv(path, index=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ap.process_float(
            7902280, mode="raw", reference=str(path), n_bootstrap=5
        )
    p = r["park"]
    assert set(p["ref_reference_name"].dropna()) == {"near"}
    assert "ref_sigma" in p
    assert set(r["wind_skill"]["reference_name"]) == {"all", "near"}
    sig = p["U_nystuen2015_fit_sigma_ref"].dropna().unique()
    assert sig == pytest.approx(0.5)


def test_pipeline_ais_remove_with_reference_uses_both_criteria(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = ap.process_float(7902280, mode="raw")
    park = base["park"]
    ref = _reference_for(park)
    t = park["datetime"].min()
    ais = pd.DataFrame(
        {
            "datetime": [
                t + pd.Timedelta(minutes=m) for m in range(0, 600, 10)
            ],
            "mmsi": 111,
            "latitude": park["latitude"].iloc[0] + 0.02,
            "longitude": park["longitude"].iloc[0],
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ap.process_float(
            7902280,
            mode="raw",
            ais=ais,
            ais_action="remove",
            reference={"buoy": ref},
        )
    p = r["park"]
    flagged, contaminated = p["ship_flag"], p["ship_contaminated"]
    assert flagged.sum() > contaminated.sum() >= 0
    # only contaminated measurements lose their acoustic values
    kept = flagged & ~contaminated
    assert p.loc[kept, "spl_8000Hz"].notna().all()
    assert p.loc[contaminated, "spl_8000Hz"].isna().all()
