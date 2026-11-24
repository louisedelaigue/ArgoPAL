# %%
"""Reference winds, model refitting, skill and the second ship criterion."""

import os
import warnings

import numpy as np
import pandas as pd
import pytest

import argopal

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
    ref = argopal.read_reference(df)
    assert list(ref.columns) == [
        "datetime",
        "wind_speed",
        "latitude",
        "longitude",
    ]
    assert ref["wind_speed"].iloc[0] == 7.5


def test_read_reference_fixed_station_and_errors():
    df = pd.DataFrame({"datetime": [T0], "wind_speed": [5.0]})
    ref = argopal.read_reference(df, latitude=43.38, longitude=7.83)
    assert ref["latitude"].iloc[0] == 43.38
    with pytest.raises(ValueError, match="latitude"):
        argopal.read_reference(df)


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
    out = argopal.collocate_reference(pal, ref, name="buoy")
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
    out = argopal.collocate_reference(
        far, ref, name="buoy", max_distance_km=40
    )
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
        argopal.collocate_reference(pal, ref, name="nystuen2015")
    with pytest.raises(KeyError, match="attach_positions"):
        argopal.collocate_reference(pal.drop(columns="latitude"), ref)


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
    ref = argopal.read_reference(_several_references())
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
    assert set(argopal.read_reference(alias)["reference_name"]) == {"A", "B"}


def test_collocate_several_references_closest_in_distance():
    pal = _pal(
        [
            T0 + pd.Timedelta("10min"),  # A and B in time: A, closer
            T0 + pd.Timedelta("3h50min"),  # A (4h00) and B (3h00) ok: A
            T0 + pd.Timedelta("4h30min"),  # B ended at 3h00, A within 1h
            T0 + pd.Timedelta("6h"),  # neither
        ]
    )
    out = argopal.collocate_reference(pal, _several_references(), name="buoys")
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
    out = argopal.collocate_reference(
        north, _several_references(), name="buoys"
    )
    assert list(out["buoys_reference_name"]) == ["B", "A"]
    np.testing.assert_allclose(out["U_buoys"], [9.0, 5.0])
    np.testing.assert_allclose(out["buoys_distance_km"], [5, 25], rtol=0.02)
    np.testing.assert_allclose(out["buoys_dt_min"], [-10, -20])


def test_collocate_several_references_radius():
    north = _pal([T0], lat=43.38 + 50 / 111.2)  # 20 km from B, 50 from A
    out = argopal.collocate_reference(
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
    out, _ = argopal.fit_wind_models(pal, max_distance_km=None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        q, _ = argopal.wind_uncertainty(out, n_bootstrap=5)
        given, _ = argopal.wind_uncertainty(
            out, n_bootstrap=5, reference_sigma=1
        )
    # one calibration per buoy: each buoy's measurements get its own sigma
    by_buoy = q.groupby("buoy_calibration")["U_nystuen2015_fit_sigma_ref"]
    assert by_buoy.first().to_dict() == pytest.approx({"A": 0.3, "B": 0.4})
    assert (given["U_nystuen2015_fit_sigma_ref"].dropna() == 1).all()
    # pooled: one sigma, the RMS of both
    pooled, _ = argopal.fit_wind_models(
        pal, max_distance_km=None, calibration="pooled"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        qp, _ = argopal.wind_uncertainty(pooled, n_bootstrap=5)
    rms = np.sqrt((0.3**2 + 0.4**2) / 2)
    assert qp["U_nystuen2015_fit_sigma_ref"].dropna().unique() == (
        pytest.approx(rms)
    )
    skill = argopal.wind_skill(out, max_distance_km=None, by_reference=True)
    assert list(skill["reference_name"].unique()) == ["all", "A", "B"]


def _ten_minute_float(hours=48, seed=1):
    """A float measuring every 10 min, an hourly buoy at the float, and a
    wind_literature table, as in `process_float`."""
    rng = np.random.default_rng(seed)
    t = pd.date_range(T0, periods=hours * 6, freq="10min")
    spl = 45 + 5 * np.sin(np.arange(len(t)) / 20) + rng.normal(0, 1, len(t))
    pal = _pal(t, spl=spl)
    pal["cycle"] = np.repeat(np.arange(1, 9), len(t) // 8)
    pal = argopal.wind_literature(pal)
    tb = pd.date_range(T0, periods=hours + 1, freq="1h")
    buoy = pd.DataFrame(
        {
            "datetime": tb,
            "wind_speed": np.linspace(3, 12, len(tb)),
            "latitude": 43.38,
            "longitude": 7.83,
        }
    )
    return argopal.collocate_reference(pal, buoy, name="buoy")


def test_average_to_reference_one_point_per_reference_value():
    pal = _ten_minute_float()
    assert "buoy_time" in pal
    pts = argopal.average_to_reference(pal)
    # each hourly value collects the measurements within ±30 min of it
    assert pts["buoy_time"].is_unique
    assert pts["n_measurements"].sum() == pal["U_buoy"].notna().sum()
    assert pts["n_measurements"].max() <= 6
    # the window: from 30 min before (included) to 30 min after (excluded)
    t5 = pts["buoy_time"].iloc[5]
    dt = pal["datetime"] - t5
    one = pal[(dt >= pd.Timedelta("-30min")) & (dt < pd.Timedelta("30min"))]
    assert pts["n_measurements"].iloc[5] == len(one) == 6
    assert pts["spl_8000Hz"].iloc[5] == pytest.approx(one["spl_8000Hz"].mean())
    # the reference value itself is kept, not averaged
    assert (
        pts["U_buoy"].iloc[5]
        == pal.loc[pal["buoy_time"] == t5, "U_buoy"].iloc[0]
    )
    # model winds from the averaged level, not averaged winds
    np.testing.assert_allclose(
        pts["U_nystuen2015"], argopal.wind_nystuen2015(pts["spl_8000Hz"])
    )
    single = argopal.average_to_reference(pal, average=False)
    assert len(single) == pal["U_buoy"].notna().sum()
    assert (single["n_measurements"] == 1).all()


def test_average_window_width():
    pal = _ten_minute_float()
    hourly = argopal.average_to_reference(pal)
    narrow = argopal.average_to_reference(pal, window="20min")
    assert narrow["n_measurements"].max() == 2  # -10 and 0 min
    assert len(narrow) == len(hourly)
    # windows wider than the time step overlap: ±1.5 h, 18 measurements
    wide = argopal.average_to_reference(pal, window="3h")
    assert wide["n_measurements"].max() == 18
    assert wide["n_measurements"].sum() > hourly["n_measurements"].sum()
    # nothing within the radius: an empty table, not an error
    assert argopal.average_to_reference(pal, max_distance_km=-1).empty
    # "nearest": the measurements matched to each value (nearest in time)
    nearest = argopal.average_to_reference(pal, window="nearest")
    assert nearest["n_measurements"].sum() == pal["U_buoy"].notna().sum()


def test_window_carried_from_fit_to_skill_and_uncertainty():
    pal = _ten_minute_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, fits = argopal.fit_wind_models(pal, window="20min")
        q, _ = argopal.wind_uncertainty(out, n_bootstrap=5)
    assert out.attrs["wind_fits"]["window"] == "20min"
    pts20 = argopal.average_to_reference(pal, window="20min")
    n20 = len(pts20)
    assert pts20["n_measurements"].max() == 2
    assert fits["n"].iloc[0] == n20
    assert fits["n_measurements"].iloc[0] == pts20["n_measurements"].sum()
    skill = argopal.wind_skill(out)
    assert skill.attrs["window"] == "20min"
    assert (
        q.attrs["wind_uncertainty"]["models"]["nystuen2015"]["all"]["n_points"]
        == n20
    )


def test_average_keeps_references_apart():
    pal = _pal(
        [T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("20min")],
        lat=[43.38, 43.38 + 25 / 111.2],
    )
    out = argopal.collocate_reference(pal, _several_references(), name="buoys")
    pts = argopal.average_to_reference(out, max_distance_km=None)
    assert list(pts["buoys_reference_name"]) == ["A", "B"]


def test_fit_and_skill_on_averaged_points():
    pal = _ten_minute_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # covariance of a synthetic fit
        out, fits = argopal.fit_wind_models(pal)
    n_hours = len(argopal.average_to_reference(pal))
    row = fits.set_index("model").loc["nystuen2015"]
    assert row["n"] == n_hours
    assert row["n_measurements"] == pal["U_buoy"].notna().sum()
    assert out["U_nystuen2015_fit"].notna().sum() == len(pal)  # all applied
    skill = argopal.wind_skill(out).set_index("wind")
    assert skill.loc["U_nystuen2015_fit", "n"] == n_hours
    assert skill.loc["U_nystuen2015_fit", "rmse"] == pytest.approx(
        row["rmse_fitted"]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, single = argopal.fit_wind_models(pal, average=False)
    assert single["n"].iloc[0] == pal["U_buoy"].notna().sum()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        q, _ = argopal.wind_uncertainty(out, n_bootstrap=5)
    info = q.attrs["wind_uncertainty"]["models"]["nystuen2015"]["all"]
    assert info["n_points"] == n_hours
    assert q["U_nystuen2015_fit_sigma"].notna().sum() == len(pal)


# --- Metrics --------------------------------------------------------------------------


def test_wind_metrics():
    ref = np.array([2.0, 4.0, 6.0, 8.0, np.nan])
    m = argopal.wind_metrics(ref + 1, ref)
    assert m["n"] == 4
    assert m["bias"] == pytest.approx(1.0)
    assert m["rmse"] == pytest.approx(1.0)
    assert m["slope"] == pytest.approx(1.0)
    assert m["intercept"] == pytest.approx(1.0)
    assert m["pearson_r2"] == pytest.approx(1.0)
    assert m["r2"] == pytest.approx(1 - 4 / 20)  # biased: R² < r²
    perfect = argopal.wind_metrics(ref, ref)
    assert perfect["rmse"] == 0 and perfect["r2"] == 1


# --- Refitting ---------------------------------------------------------------------------


def _synthetic_float(n=400, seed=0, true=None):
    """SPL and a reference wind made with known Nystuen coefficients."""
    rng = np.random.default_rng(seed)
    spl = rng.uniform(35, 60, n)
    true = true or {"a3": 0.0004, "a2": -0.02, "a1": 0.3, "a0": 1.0}
    u = argopal.wind_nystuen2015(spl, **true) + rng.normal(0, 0.2, n)
    pal = _pal(pd.date_range(T0, periods=n, freq="5min"), spl=spl)
    pal["U_nystuen2015"] = argopal.wind_nystuen2015(spl)
    pal["U_buoy"] = u
    pal["buoy_distance_km"] = np.where(np.arange(n) < n // 2, 10.0, 60.0)
    pal.attrs["wind_models"] = {
        "applied": {"nystuen2015": {"band_hz": 8000, "reference": "Nystuen"}}
    }
    pal.attrs["references"] = {"buoy": {"column": "U_buoy"}}
    return pal, true


def test_fit_recovers_known_coefficients():
    pal, true = _synthetic_float()
    out, fits = argopal.fit_wind_models(pal)
    got = fits.set_index("model").loc["nystuen2015", "fitted"]
    pred_true = argopal.wind_nystuen2015(pal["spl_8000Hz"], **true)
    np.testing.assert_allclose(out["U_nystuen2015_fit"], pred_true, atol=0.1)
    assert set(got) == {"a3", "a2", "a1", "a0"}
    row = fits.iloc[0]
    assert row["rmse_fitted"] < row["rmse_published"]
    assert row["rmse_fitted"] == pytest.approx(0.2, abs=0.05)  # the noise


def test_fit_uses_only_nearby_measurements():
    pal, _ = _synthetic_float()
    _, fits = argopal.fit_wind_models(pal, max_distance_km=40)
    assert fits["n"].iloc[0] == 200
    _, all_ = argopal.fit_wind_models(pal, max_distance_km=None)
    assert all_["n"].iloc[0] == 400
    use = np.zeros(len(pal), bool)
    use[:50] = True
    _, sub = argopal.fit_wind_models(pal, use=use)
    assert sub["n"].iloc[0] == 50
    # nothing close enough: the model keeps its row, unfitted
    with pytest.warns(UserWarning, match="too few"):
        out, none = argopal.fit_wind_models(pal, max_distance_km=5)
    assert none["n"].iloc[0] == 0 and none["fitted"].iloc[0] is None
    assert np.isnan(none["rmse_fitted"].iloc[0])
    assert "U_nystuen2015_fit" not in out


def test_wind_skill_table():
    pal, _ = _synthetic_float()
    out, _ = argopal.fit_wind_models(pal)
    skill = argopal.wind_skill(out)
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
    out = argopal.flag_ship_deviation(pal, reference="buoy")
    # 1 has a large error but no ship: not contaminated
    assert list(np.flatnonzero(out["ship_contaminated"])) == [0]
    assert out.attrs["ais"]["n_contaminated"] == 1
    rm = argopal.flag_ship_deviation(pal, reference="buoy", action="remove")
    assert np.isnan(rm.loc[0, "spl_8000Hz"]) and np.isnan(
        rm.loc[0, "U_nystuen2015_fit"]
    )
    assert rm.loc[0, "U_buoy"] == pal.loc[0, "U_buoy"]  # reference kept
    with pytest.raises(KeyError, match="flag_ships"):
        argopal.flag_ship_deviation(
            pal.drop(columns="ship_flag"), reference="buoy"
        )


# --- Real data: the 2025 float with DYFAMED -----------------------------------


@pytest.fixture(scope="module")
def float2025():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = F.process(ais=F.AIS_CSV)["park"]
        p = argopal.collocate_reference(p, DYFAMED, name="dyfamed")
        p, fits = argopal.fit_wind_models(p, reference="dyfamed")
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
    out = argopal.flag_ship_deviation(p, reference="dyfamed")
    assert out["ship_contaminated"].sum() > 0
    assert (out["ship_contaminated"] <= out["ship_flag"]).all()  # a subset
    assert 0 < out.attrs["ais"]["deviation_rmse"] < 2


def test_2025_figures(float2025):
    p, _ = float2025
    _, axes = argopal.plot_wind(p)
    assert len(axes) == 2  # published, refitted
    _, axes = argopal.plot_wind_fit(p)
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
            "wind_speed": 0.8 * argopal.wind_nystuen2015(spl.to_numpy()),
            "latitude": park["latitude"].iloc[0],
            "longitude": park["longitude"].iloc[0],
        }
    )


def test_pipeline_with_reference(offline):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = argopal.process_float(7902280, mode="raw")
        ref = _reference_for(base["park"])
        r = argopal.process_float(7902280, mode="raw", reference={"buoy": ref})
    p = r["park"]
    assert {"U_buoy", "U_nystuen2015_fit"} <= set(p.columns)
    assert r["wind_fits"].attrs["reference"] == "buoy"
    assert "wind_skill" in r and "wind_fit" in argopal.plot_all(r)
    # single source: named "ref"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        single = argopal.process_float(7902280, mode="raw", reference=ref)
    assert "U_ref" in single["park"]
    assert single["wind_fits"].attrs["average"] is True
    n_avg = single["wind_fits"]["n"].iloc[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        every = argopal.process_float(
            7902280, mode="raw", reference=ref, reference_average=False
        )
    assert every["wind_fits"]["n"].iloc[0] > n_avg
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        short = argopal.process_float(
            7902280, mode="raw", reference=ref, reference_window="20min"
        )
    fits = short["wind_fits"]
    assert fits.attrs["window"] == "20min"
    assert (
        fits["n_measurements"].iloc[0]
        < single["wind_fits"]["n_measurements"].iloc[0]
    )


def test_pipeline_with_several_references_in_one_file(offline, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = argopal.process_float(7902280, mode="raw")
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
        r = argopal.process_float(
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
        base = argopal.process_float(7902280, mode="raw")
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
        r = argopal.process_float(
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


def test_reference_positions_and_map():
    pal = _pal(
        [T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("20min")],
        lat=[43.38, 43.38 + 25 / 111.2],
    )
    refs = _several_references()
    # buoy B also reported once at a second position, 10 km away
    extra = refs[refs["reference_name"] == "B"].iloc[:1].copy()
    extra["latitude"] += 10 / 111.2
    out = argopal.collocate_reference(
        pal, pd.concat([refs, extra]), name="buoys"
    )
    pos = out.attrs["references"]["buoys"]["positions"]
    assert set(pos) == {"A", "B"}
    assert pos["A"]["latitude"] == pytest.approx(43.38)
    assert pos["A"]["n_positions"] == 1 and pos["A"]["spread_km"] == 0
    # the most frequent position, and how far the other one is
    assert pos["B"]["latitude"] == pytest.approx(43.38 + 30 / 111.2, abs=1e-3)
    assert pos["B"]["n_positions"] == 2
    assert pos["B"]["spread_km"] == pytest.approx(10, rel=0.05)

    _, ax = argopal.plot_reference_map(out, radius_km=20, annotate=True)
    texts = [t.get_text() for t in ax.texts]
    assert "A" in texts and "B" in texts  # both used within 20 km
    # within 1 km only the measurement at A is used: B is not shown
    _, ax = argopal.plot_reference_map(out, radius_km=1)
    texts = [t.get_text() for t in ax.texts]
    assert "A" in texts and "B" not in texts
    with pytest.raises(KeyError, match="collocate_reference"):
        argopal.plot_reference_map(pal.assign(U_x=1.0))


def test_reference_timeseries():
    pal = _ten_minute_float()
    # before refitting: the reference and the published model only
    _, axes = argopal.plot_reference_timeseries(pal)
    assert len(axes) == 1
    assert [ln.get_label() for ln in axes[0].lines][:2] == [
        "BUOY (reference)",
        "Nystuen et al. (2015), published",
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fitted, _ = argopal.fit_wind_models(pal)
    _, axes = argopal.plot_reference_timeseries(fitted, annotate=True)
    labels = [ln.get_label() for ln in axes[0].lines]
    assert "Nystuen et al. (2015), refitted" in labels
    # ±1σ shading once the uncertainty is computed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with_sigma, _ = argopal.wind_uncertainty(fitted, n_bootstrap=5)
    _, axes = argopal.plot_reference_timeseries(with_sigma)
    legend = [t.get_text() for t in axes[0].get_legend().get_texts()]
    assert legend[-1] == "Nystuen et al. (2015), refitted ±1σ"
    _, axes = argopal.plot_reference_timeseries(with_sigma, uncertainty=False)
    assert not axes[0].collections
    # one point per hourly buoy value, as in the fit
    n_points = len(argopal.average_to_reference(fitted))
    ref_line = axes[0].lines[0]
    assert np.isfinite(ref_line.get_ydata()).sum() == n_points
    with pytest.raises(ValueError, match="within"):
        argopal.plot_reference_timeseries(fitted, radius_km=-1)


def test_reference_timeseries_one_panel_per_reference():
    pal = _pal(
        [T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("4h20min")],
        lat=[43.38 + 25 / 111.2, 43.38],
    )
    pal["U_nystuen2015"] = 5.0
    out = argopal.collocate_reference(pal, _several_references(), name="b")
    _, axes = argopal.plot_reference_timeseries(out)
    # B first (it is used first in time), then A
    assert [ax.get_title().split(":")[0] for ax in axes] == ["B", "A"]


def test_reference_plots_leave_out_measurements_without_float_wind():
    # two measurements near buoy A, the second without a float wind (e.g.
    # after the delayed-mode data end): only the first is used
    pal = _pal([T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("4h10min")])
    pal["spl_8000Hz"] = [50.0, np.nan]  # no depth correction: no level
    pal["U_nystuen2015"] = [5.0, np.nan]
    pal.attrs["wind_models"] = {
        "applied": {"nystuen2015": {"band_hz": 8000, "reference": "N"}}
    }
    out = argopal.collocate_reference(pal, _several_references(), name="b")
    assert out["b_latitude"].iloc[0] == pytest.approx(43.38)
    _, ax = argopal.plot_reference_map(out, radius_km=40)
    used = next(ln for ln in ax.lines if "Used for fitting" in ln.get_label())
    assert len(used.get_xdata()) == 1
    labels = [ln.get_label() for ln in ax.lines]
    assert any(x.startswith("Trajectory with float wind") for x in labels)
    assert any(x.startswith("Trajectory without float wind") for x in labels)
    _, axes = argopal.plot_reference_timeseries(out, radius_km=40)
    assert len(axes) == 1
    assert np.isfinite(axes[0].lines[0].get_ydata()).sum() == 1


def test_reference_map_greys_out_buoys_passed_without_float_wind():
    # near A with a float wind, then near B only without one
    pal = _pal(
        [T0 + pd.Timedelta("10min"), T0 + pd.Timedelta("30min")],
        lat=[43.38, 43.38 + 30 / 111.2],
    )
    pal["spl_8000Hz"] = [50.0, np.nan]
    pal["U_nystuen2015"] = [5.0, np.nan]
    pal.attrs["wind_models"] = {
        "applied": {"nystuen2015": {"band_hz": 8000, "reference": "N"}}
    }
    out = argopal.collocate_reference(pal, _several_references(), name="b")
    assert list(out["b_reference_name"]) == ["A", "B"]
    _, ax = argopal.plot_reference_map(out, radius_km=10, annotate=True)
    labels = [ln.get_label() for ln in ax.lines]
    assert "Reference used for fitting" in labels
    assert "Reference not used (no float wind)" in labels
    colours = {t.get_text(): t.get_color() for t in ax.texts}
    assert colours["B"] != colours["A"]  # B greyed out


# --- Reference wind at 10 m --------------------------------------------------------


def test_wind_to_10m():
    assert argopal.wind_to_10m(7.0, 10.0) == pytest.approx(7.0)
    assert argopal.wind_to_10m(7.0, 10.0, method="power") == pytest.approx(7)
    low = argopal.wind_to_10m(5.0, 4.0)
    assert 5.0 < low < 6.0  # stronger at 10 m than at 4 m
    # log profile: U10 / U(z) = ln(10 / z0) / ln(z / z0)
    z0 = 1.52e-4
    assert low == pytest.approx(5.0 * np.log(10 / z0) / np.log(4 / z0))
    assert argopal.wind_to_10m(5.0, 4.0, method="power") == pytest.approx(
        5.0 * (10 / 4) ** 0.11
    )
    with pytest.raises(ValueError, match="positive"):
        argopal.wind_to_10m(5.0, 0)
    with pytest.raises(ValueError, match="method"):
        argopal.wind_to_10m(5.0, 4.0, method="cubic")


def test_read_reference_converts_to_10m():
    refs = _several_references()  # A and B, no height
    with pytest.warns(argopal.ReferenceHeightWarning, match="A, B"):
        same = argopal.read_reference(refs)
    assert same.attrs["height"]["height_m"]["A"] == "unknown, assumed 10 m"
    # one height for all
    four = argopal.read_reference(refs, height_m=4)
    a = four[four["reference_name"] == "A"]["wind_speed"].iloc[0]
    assert a == pytest.approx(argopal.wind_to_10m(5.0, 4.0))
    # per reference; the one without is warned about
    with pytest.warns(argopal.ReferenceHeightWarning, match="B"):
        per = argopal.read_reference(refs, height_m={"A": 4})
    b = per[per["reference_name"] == "B"]["wind_speed"].iloc[0]
    assert b == pytest.approx(9.0)  # B unchanged, assumed 10 m
    # a height column (e.g. 10 for winds already at 10 m)
    col = argopal.read_reference(refs.assign(sensor_height=10.0))
    assert col["wind_speed"].tolist() == same["wind_speed"].tolist()
    assert col.attrs["height"]["height_m"]["B"] == [10.0]


def test_collocation_records_the_height():
    pal = _pal([T0 + pd.Timedelta("10min")])
    out = argopal.collocate_reference(
        pal, _several_references(), name="b", height_m=4
    )
    assert out["U_b"].iloc[0] == pytest.approx(argopal.wind_to_10m(5.0, 4.0))
    h = out.attrs["references"]["b"]["height"]
    assert h["method"] == "log" and h["height_m"]["A"] == [4.0]


def _two_buoy_float(seed=3):
    """A float 24 h at buoy A, then 24 h at buoy B (30 km north), then a
    few measurements 80 km north of A (50 km from B): too far for the
    refit.  Buoy B's wind is twice as sensitive to the level as A's."""
    rng = np.random.default_rng(seed)
    t = pd.date_range(T0, periods=48 * 6, freq="10min")
    spl = rng.uniform(40, 55, len(t))
    lat = np.where(np.arange(len(t)) < 24 * 6, 43.38, 43.38 + 30 / 111.2)
    far = np.arange(len(t)) >= len(t) - 6
    lat = np.where(far, 43.38 + 80 / 111.2, lat)
    pal = _pal(t, spl=spl)
    pal["latitude"] = lat
    pal["cycle"] = np.repeat(np.arange(1, 9), len(t) // 8)
    pal = argopal.wind_literature(pal)
    hourly = pal.set_index("datetime")["spl_8000Hz"].resample("1h").mean()
    tables = []
    for name, la, gain in (("A", 43.38, 0.4), ("B", 43.38 + 30 / 111.2, 0.8)):
        tables.append(
            pd.DataFrame(
                {
                    "reference_name": name,
                    "datetime": hourly.index,
                    "latitude": la,
                    "longitude": 7.83,
                    "wind_speed": 2 + gain * (hourly.to_numpy() - 40),
                    "height_m": 10.0,
                }
            )
        )
    out = argopal.collocate_reference(pal, pd.concat(tables), name="b")
    return out, far


def test_calibration_per_reference():
    pal, far = _two_buoy_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, fits = argopal.fit_wind_models(pal, models=["nystuen2015"])
    coefs = out.attrs["wind_fits"]["coefficients"]["nystuen2015"]
    assert set(coefs) == {"A", "B"}
    assert out.attrs["wind_fits"]["calibration"] == "per_reference"
    assert list(fits["reference_name"]) == ["A", "B"]
    # each buoy's points use its coefficients; the far ones (not used in
    # the refit) use the closest buoy's, B's
    cal = out["b_calibration"]
    used = out["b_used"]
    assert not used[far].any()
    assert (cal[used] == out.loc[used, "b_reference_name"]).all()
    assert (cal[far] == "B").all()
    assert fits["n_applied"].sum() == len(out)
    spl = out.loc[far, "spl_8000Hz"]
    np.testing.assert_allclose(
        out.loc[far, "U_nystuen2015_fit"],
        argopal.wind_nystuen2015(spl, **coefs["B"]),
    )
    # B's wind rises faster with the level than A's
    at = np.array([45.0, 52.0])
    rise_a = np.diff(argopal.wind_nystuen2015(at, **coefs["A"]))[0]
    rise_b = np.diff(argopal.wind_nystuen2015(at, **coefs["B"]))[0]
    assert rise_b > 1.5 * rise_a
    # pooled: one set for every measurement
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pooled, pfits = argopal.fit_wind_models(
            pal, models=["nystuen2015"], calibration="pooled"
        )
    assert list(pfits["reference_name"]) == ["all"]
    assert (pooled["b_calibration"] == "all").all()
    with pytest.raises(ValueError, match="calibration"):
        argopal.fit_wind_models(pal, calibration="nearest")


def test_uncertainty_per_reference():
    pal, _ = _two_buoy_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, _ = argopal.fit_wind_models(pal, models=["nystuen2015"])
        q, bins = argopal.wind_uncertainty(out, n_bootstrap=5)
    models = q.attrs["wind_uncertainty"]["models"]["nystuen2015"]
    assert set(models) == {"A", "B"}
    assert models["A"]["n_applied"] + models["B"]["n_applied"] == len(q)
    assert set(bins["reference_name"]) <= {"A", "B"}
    assert q["U_nystuen2015_fit_sigma"].notna().any()


def test_wind_timeline():
    pal, _ = _two_buoy_float()
    # without a refit: the published wind only
    _, axes = argopal.plot_wind_timeline(pal)
    assert len(axes) == 1  # two days: one row
    labels = [ln.get_label() for ln in axes[0].lines]
    assert "Nystuen et al. (2015), published" in labels
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, _ = argopal.fit_wind_models(pal, models=["nystuen2015"])
        out, _ = argopal.wind_uncertainty(out, n_bootstrap=5)
    _, axes = argopal.plot_wind_timeline(out, annotate=True)
    legend = [t.get_text() for t in axes[0].get_legend().get_texts()]
    assert "Nystuen et al. (2015), refitted" in legend
    assert any(t.startswith("A: used for refit") for t in legend)
    assert any(t.startswith("B: used for refit") for t in legend)
    # shaded refit periods and the coefficient strip
    assert len(axes[0].patches) >= 4
    # several rows for a longer record
    longer = out.copy()
    longer["datetime"] = longer["datetime"] + pd.to_timedelta(
        np.arange(len(longer)) * 6, unit="h"
    )
    _, axes = argopal.plot_wind_timeline(longer, months=1)
    assert len(axes) >= 2


def test_minimum_cycles_for_own_calibration():
    pal, _ = _two_buoy_float()  # 4 cycles at A, 4 at B
    # B in only 2 cycles: no coefficients of its own, A's are used
    fewer = pal.copy()
    at_b = fewer["b_reference_name"] == "B"
    fewer.loc[at_b, "cycle"] = np.where(fewer.loc[at_b, "cycle"] < 7, 5, 6)
    with pytest.warns(UserWarning, match=r"B \(2\)"):
        out, fits = argopal.fit_wind_models(fewer, models=["nystuen2015"])
    coefs = out.attrs["wind_fits"]["coefficients"]["nystuen2015"]
    assert set(coefs) == {"A"}
    row_b = fits.set_index("reference_name").loc["B"]
    assert row_b["fitted"] is None and row_b["n_cycles"] == 2
    assert np.isfinite(row_b["rmse_published"])  # still compared
    # B's measurements (and those closest to B) get no refitted wind
    near_b = out["b_closest"] == "B"
    assert near_b.any() and (out["b_closest"] == "A").any()
    assert out.loc[near_b, "b_calibration"].isna().all()
    assert out.loc[near_b, "U_nystuen2015_fit"].isna().all()
    near_a = out["b_closest"] == "A"
    assert (out.loc[near_a, "b_calibration"] == "A").all()
    assert out.loc[near_a & out["U_nystuen2015"].notna(),
                   "U_nystuen2015_fit"].notna().all()  # fmt: skip
    # the timeline marks where there is no calibration
    _, axes = argopal.plot_wind_timeline(out)
    legend = [t.get_text() for t in axes[0].figure.legends[0].get_texts()] if (
        axes[0].figure.legends) else [
        t.get_text() for t in axes[0].get_legend().get_texts()]  # fmt: skip
    assert any(x.startswith("No calibration") for x in legend)
    # no reference with enough cycles: one fit for all, with a warning
    with pytest.warns(UserWarning, match="all references together"):
        pooled, pfits = argopal.fit_wind_models(
            pal, models=["nystuen2015"], min_cycles=5
        )
    assert list(pfits["reference_name"]) == ["all"]
    assert pooled.attrs["wind_fits"]["calibration"] == "pooled"
    assert (pooled["b_calibration"] == "all").all()


def test_plot_limits_follow_the_wind_not_sigma():
    pal, _ = _two_buoy_float()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out, _ = argopal.fit_wind_models(pal, models=["nystuen2015"])
        out, _ = argopal.wind_uncertainty(out, n_bootstrap=5)
    out["U_nystuen2015_fit_sigma"] = 1000.0  # an absurd uncertainty
    for plot in (
        argopal.plot_wind,
        argopal.plot_wind_timeline,
        argopal.plot_reference_timeseries,
    ):
        _, axes = plot(out)
        for ax in np.atleast_1d(axes):
            assert ax.get_ylim()[1] < 100, plot.__name__
