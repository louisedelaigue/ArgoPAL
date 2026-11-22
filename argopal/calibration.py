"""Refit the wind models to a reference wind, and measure their skill.

As in Delaigue et al. (2026, Sect. 2.6): the published models are refitted
by nonlinear least squares to the float's depth-corrected 8 kHz level and
a collocated reference wind (e.g. a buoy), using only measurements within
`max_distance_km` (40 km) of the reference.  The refitted coefficients are
then applied to the whole record.
"""

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit

from .wind import COEFFICIENTS, LITERATURE_MODELS


def wind_metrics(estimate, reference):
    """Skill of a wind estimate against a reference.

    Parameters
    ----------
    estimate, reference : array_like
        Wind speeds (m s-1); pairs with a NaN are ignored.

    Returns
    -------
    dict
        n; bias (mean of estimate - reference); rmse (root mean square of
        estimate - reference); r2 (coefficient of determination,
        1 - SS_res / SS_tot, with the reference as truth); pearson_r2
        (squared correlation); slope and intercept of the least-squares
        line estimate = slope · reference + intercept.
    """
    e = np.asarray(estimate, dtype=float)
    r = np.asarray(reference, dtype=float)
    ok = np.isfinite(e) & np.isfinite(r)
    e, r = e[ok], r[ok]
    n = int(ok.sum())
    if n < 2:
        return {
            "n": n,
            "bias": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "pearson_r2": np.nan,
            "slope": np.nan,
            "intercept": np.nan,
        }
    diff = e - r
    ss_tot = np.sum((r - r.mean()) ** 2)
    slope, intercept = np.polyfit(r, e, 1)
    return {
        "n": n,
        "bias": float(diff.mean()),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "r2": float(1 - np.sum(diff**2) / ss_tot) if ss_tot > 0 else np.nan,
        "pearson_r2": float(np.corrcoef(e, r)[0, 1] ** 2),
        "slope": float(slope),
        "intercept": float(intercept),
    }


def _reference_column(pal, reference):
    refs = pal.attrs.get("references", {})
    if reference is None:
        if len(refs) != 1:
            raise ValueError(
                "Name the reference to use: "
                f"{', '.join(refs) or 'none collocated (collocate_reference)'}."
            )
        reference = next(iter(refs))
    col = f"U_{reference}"
    if col not in pal:
        raise KeyError(f"No {col!r}: run collocate_reference first.")
    return reference, col


def _within(pal, reference, max_distance_km):
    dist_col = f"{reference}_distance_km"
    if max_distance_km is None or dist_col not in pal:
        return np.ones(len(pal), dtype=bool)
    return (pal[dist_col] <= max_distance_km).to_numpy()


def average_to_reference(
    pal,
    reference=None,
    max_distance_km=40.0,
    use=None,
    average=True,
    window="1h",
):
    """Average the measurements around each reference value.

    The float measures more often than most references (e.g. every
    5 min against an hourly buoy wind).  For each reference value, the
    measurements within ±`window`/2 of its time (from the same reference,
    and within `max_distance_km` of it) are averaged: one point per
    reference value, on the reference's own time scale.  This is what
    `fit_wind_models`, `wind_skill`, `wind_uncertainty` and
    `plot_wind_fit` compare by default.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `collocate_reference`.
    reference : str, optional
        Name given to `collocate_reference`; optional if only one.
    max_distance_km : float or None, optional
        Only measurements this close to their reference, by default
        40 km.  None uses all.
    use : array of bool, optional
        Further restricts the measurements.
    average : bool, optional
        False returns the selected measurements without averaging.
    window : str or pd.Timedelta, optional
        Width of the averaging window, centred on each reference time:
        "1h" (default) averages the measurements from 30 min before
        (included) to 30 min after (excluded), so that hourly windows do
        not share measurements.  Windows wider than the reference's time step
        overlap, so a measurement then counts towards several reference
        values.  Only measurements collocated with the reference (with a
        reference value within `max_time` of `collocate_reference`, 1 h
        by default) are averaged.  "nearest" averages instead the
        measurements matched to each reference value (the one nearest in
        time).

    Returns
    -------
    pd.DataFrame
        One row per reference value with at least one measurement in its
        window: the reference value and its time, the reference name if
        any, the mean of every numeric measurement column (band and
        spectral levels, distance, ...), the mean time, the first cycle,
        and "n_measurements".  Levels in dB are averaged as they are, as
        in despiking, and the model wind speeds ("U_<model>",
        "U_<model>_fit") are computed from the averaged spectral level.
    """
    reference, ref_col = _reference_column(pal, reference)
    rows_ok = _within(pal, reference, max_distance_km)
    if use is not None:
        rows_ok = rows_ok & np.asarray(use, dtype=bool)
    rows_ok = rows_ok & pal[ref_col].notna().to_numpy()
    sub = pal.loc[rows_ok]
    time_col = f"{reference}_time"
    if not average or time_col not in sub:
        out = sub.copy()
        out["n_measurements"] = 1
        out.attrs = dict(pal.attrs)
        return out.reset_index(drop=True)
    if sub.empty:  # nothing close enough: no points
        out = sub.iloc[:0].copy()
        out["n_measurements"] = pd.Series(dtype="int64")
        out.attrs = dict(pal.attrs)
        return out
    name_col = f"{reference}_reference_name"
    sigma_col = f"{reference}_sigma"
    # values of the reference itself, not averaged over measurements
    ref_cols = [c for c in (ref_col, sigma_col) if c in sub]
    keys = [c for c in (name_col, time_col) if c in sub]
    numeric = [
        c
        for c in sub.select_dtypes("number").columns
        if c not in keys + ref_cols
    ]
    if window == "nearest":
        grouped = sub.groupby(keys, sort=True)
        agg = {c: "mean" for c in numeric + ref_cols}
        if "cycle" in agg:
            agg["cycle"] = "first"
        if "datetime" in sub:
            agg["datetime"] = "mean"
        out = grouped.agg(agg)
        out["n_measurements"] = grouped.size()
        out = out.reset_index()
    else:
        half = pd.Timedelta(window) / 2
        groups = sub.groupby(name_col) if name_col in sub else [(None, sub)]
        out = pd.concat(
            [
                _window_means(g, time_col, ref_cols, numeric, half).assign(
                    _name=name
                )
                for name, g in groups
            ],
            ignore_index=True,
        )
        name = out.pop("_name")
        if name_col in sub:
            out.insert(0, name_col, name)
        out = out.sort_values(keys).reset_index(drop=True)
    # model winds from the averaged level, as the refit uses them
    fitted = pal.attrs.get("wind_fits", {}).get("coefficients", {})
    for name, (func, band, _) in LITERATURE_MODELS.items():
        spl = f"spl_{band}Hz"
        if spl not in out:
            continue
        level = out[spl].to_numpy(float)
        if f"U_{name}" in out:
            out[f"U_{name}"] = func(level, **COEFFICIENTS[name])
        if f"U_{name}_fit" in out and name in fitted:
            out[f"U_{name}_fit"] = func(level, **fitted[name])
    out.attrs = dict(pal.attrs)
    return out


def _window_means(g, time_col, ref_cols, numeric, half):
    """Means of the measurements of one reference within ±`half` of each
    of its values (cumulative sums: fast, windows may overlap)."""
    g = g.sort_values("datetime")
    t = g["datetime"].to_numpy(dtype="datetime64[ns]")
    refs = g.drop_duplicates(time_col).sort_values(time_col)
    tr = refs[time_col].to_numpy(dtype="datetime64[ns]")
    lo = np.searchsorted(t, tr - half.to_timedelta64(), side="left")
    # [t - half, t + half): back-to-back windows share no measurement
    hi = np.searchsorted(t, tr + half.to_timedelta64(), side="left")
    vals = g[numeric].astype(float).to_numpy()
    finite = np.isfinite(vals)
    zero = np.zeros((1, vals.shape[1]))
    csum = np.vstack([zero, np.cumsum(np.where(finite, vals, 0.0), axis=0)])
    ccnt = np.vstack([zero, np.cumsum(finite, axis=0)])
    cnt = ccnt[hi] - ccnt[lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        means = (csum[hi] - csum[lo]) / cnt
    out = pd.DataFrame(means, columns=numeric)
    secs = (t - t[0]).astype("timedelta64[ns]").astype(float) / 1e9
    csec = np.concatenate([[0.0], np.cumsum(secs)])
    n = hi - lo
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_secs = (csec[hi] - csec[lo]) / n
    if "datetime" in g:
        out["datetime"] = t[0] + pd.to_timedelta(mean_secs, unit="s")
    if "cycle" in g:
        cyc = g["cycle"].to_numpy()
        out["cycle"] = np.where(n > 0, cyc[np.minimum(lo, len(cyc) - 1)], 0)
    out[time_col] = tr
    for c in ref_cols:
        out[c] = refs[c].to_numpy()
    out["n_measurements"] = n
    return out[n > 0]


def wind_skill(
    pal,
    reference=None,
    columns=None,
    max_distance_km=40.0,
    by_reference=False,
    average=None,
    window=None,
):
    """Skill of every wind estimate in a table against a reference.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table with wind columns ("U_<model>", "U_<model>_fit") and a
        collocated reference.
    reference : str, optional
        Name given to `collocate_reference`; optional if only one.
    columns : list of str, optional
        Wind columns to assess; by default every model column.
    max_distance_km : float or None, optional
        Only measurements this close to the reference, by default 40 km.
    by_reference : bool, optional
        If the reference table names several references, also give the
        skill against each of them separately (column "reference_name",
        "all" for the whole set).  By default False.
    average : bool, optional
        Compare one point per reference value (`average_to_reference`);
        by default as for `fit_wind_models` (True if not run).
    window : str, optional
        Averaging window around each reference time; by default as for
        `fit_wind_models` ("1h" if not run).

    Returns
    -------
    pd.DataFrame
        One row per wind column (and per reference): the output of
        `wind_metrics`.
    """
    reference, ref_col = _reference_column(pal, reference)
    if columns is None:
        names = list(LITERATURE_MODELS)
        columns = [
            c for n in names for c in (f"U_{n}", f"U_{n}_fit") if c in pal
        ]
    if average is None:
        average = pal.attrs.get("wind_fits", {}).get("average", True)
    if window is None:
        window = pal.attrs.get("wind_fits", {}).get("window", "1h")
    pts = average_to_reference(
        pal, reference, max_distance_km, average=average, window=window
    )
    name_col = f"{reference}_reference_name"
    subsets = [("all", np.ones(len(pts), dtype=bool))]
    if by_reference and name_col in pts:
        names = pts[name_col].dropna().unique()
        subsets += [
            (str(r), (pts[name_col] == r).to_numpy()) for r in sorted(names)
        ]
    rows = []
    for label, rows_ok in subsets:
        for c in columns:
            m = wind_metrics(pts.loc[rows_ok, c], pts.loc[rows_ok, ref_col])
            row = {"wind": c, "reference": reference}
            if len(subsets) > 1:
                row["reference_name"] = label
            rows.append({**row, **m})
    out = pd.DataFrame(rows)
    out.attrs = {
        "max_distance_km": max_distance_km,
        "average": average,
        "window": window,
    }
    return out


_FIT_METRICS = ("rmse", "r2", "bias")
_FIT_STATS = ("published", "fitted")


def fit_wind_models(
    pal,
    reference=None,
    models=None,
    max_distance_km=40.0,
    use=None,
    average=True,
    window="1h",
):
    """Refit the literature wind models to a collocated reference wind.

    Each model keeps its published form; its coefficients are refitted by
    nonlinear least squares (scipy.optimize.curve_fit), starting from the
    published values, on the measurements within `max_distance_km` of the
    reference (as in Delaigue et al. 2026).  By default the measurements
    within ±30 min of each reference value are first averaged
    (`average_to_reference`, `window`), so the fit compares one point per
    reference value.  The refitted model is then applied to every measurement as
    "U_<model>_fit".

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `wind_literature` and `collocate_reference`.
    reference : str, optional
        Name given to `collocate_reference`; optional if only one.
    models : list of str, optional
        Models to refit; by default every one that was applied.
    max_distance_km : float or None, optional
        Only measurements this close to the reference are used for the fit
        and the skill, by default 40 km.  None uses all.
    use : array of bool, optional
        Further restricts the measurements used for fitting (e.g. to leave
        some out for validation).
    average : bool, optional
        Average the measurements around each reference value before
        fitting, by default True.  False fits every measurement against
        the reference value nearest in time.
    window : str or pd.Timedelta, optional
        Width of the averaging window, centred on each reference time, by
        default "1h" (30 min before to 30 min after); "nearest" averages
        the measurements matched to each reference value.  See
        `average_to_reference`.

    Returns
    -------
    pal : pd.DataFrame
        A copy with "U_<model>_fit" columns; `attrs["wind_fits"]` holds the
        coefficients.
    fits : pd.DataFrame
        One row per model: n points used (reference values, or
        measurements if not averaged), n_measurements behind them,
        published and refitted coefficients, and the skill of both on the
        fitting points (rmse, r2, bias).
        A model that could not be refitted (too few measurements) keeps
        its row, with fitted None and NaN skill.
    """
    reference, ref_col = _reference_column(pal, reference)
    applied = pal.attrs.get("wind_models", {}).get("applied", {})
    if models is None:
        models = list(applied)
    pts = average_to_reference(
        pal,
        reference,
        max_distance_km,
        use=use,
        average=average,
        window=window,
    )
    out = pal.copy()
    rows, coefs = [], {}
    for name in models:
        func, band, _ = LITERATURE_MODELS[name]
        spl_col = f"spl_{band}Hz"
        if spl_col not in pal:
            warnings.warn(f"{name}: {band} Hz not processed, not refitted.")
            continue
        x = pts[spl_col].to_numpy(float)
        y = pts[ref_col].to_numpy(float)
        p0 = COEFFICIENTS[name]
        lit = func(x, **p0)
        sel = np.isfinite(x) & np.isfinite(y) & np.isfinite(lit)
        info = {
            "model": name,
            "reference": reference,
            "n": int(sel.sum()),
            "n_measurements": int(pts["n_measurements"].to_numpy()[sel].sum()),
            "published": p0,
            "fitted": None,
            **{f"{m}_{k}": np.nan for m in _FIT_METRICS for k in _FIT_STATS},
        }
        if sel.sum() <= len(p0):
            warnings.warn(
                f"{name}: too few collocations ({sel.sum()}) to fit."
            )
            rows.append(info)
            continue
        try:
            popt, _ = curve_fit(
                lambda s, *p, f=func: f(s, *p),
                x[sel],
                y[sel],
                p0=list(p0.values()),
                maxfev=20000,
            )
        except (RuntimeError, ValueError) as e:
            warnings.warn(f"{name}: fit failed ({e}).")
            rows.append(info)
            continue
        fitted = dict(zip(p0, map(float, popt)))
        out[f"U_{name}_fit"] = func(pal[spl_col].to_numpy(float), **fitted)
        before = wind_metrics(lit[sel], y[sel])
        after = wind_metrics(func(x[sel], **fitted), y[sel])
        info["fitted"] = fitted
        for m in _FIT_METRICS:
            info[f"{m}_published"] = before[m]
            info[f"{m}_fitted"] = after[m]
        rows.append(info)
        coefs[name] = fitted
    out.attrs = dict(pal.attrs)
    out.attrs["wind_fits"] = {
        "reference": reference,
        "max_distance_km": max_distance_km,
        "average": average,
        "window": window,
        "coefficients": coefs,
    }
    fits = pd.DataFrame(rows)
    fits.attrs = {
        "reference": reference,
        "max_distance_km": max_distance_km,
        "average": average,
        "window": window,
    }
    return out, fits


def flag_ship_deviation(
    pal,
    reference=None,
    wind="nystuen2015_fit",
    max_distance_km=40.0,
    action="flag",
):
    """Second ship-noise criterion: a ship nearby AND a wind error larger
    than the RMSE of uncontaminated measurements (Delaigue et al. 2026,
    Sect. 2.5).

    The RMSE is computed between `wind` and the reference on measurements
    without a ship nearby (ship_flag False) within `max_distance_km`.  A
    measurement is "ship_contaminated" when ship_flag is True and the
    absolute difference between its wind and the reference exceeds that
    RMSE.  Only measurements with a
    reference value can be judged.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `flag_ships`, with a collocated reference and wind.
    reference : str, optional
        Name given to `collocate_reference`; optional if only one.
    wind : str, optional
        Wind estimate to compare, without the `"U_"` prefix; by default the
        refitted Nystuen et al. (2015) model, as in the paper.
    max_distance_km : float or None, optional
        By default 40 km.
    action : {"flag", "remove"}, optional
        "remove" also blanks the acoustic values (f_*, spl_* and wind-model
        columns) of contaminated measurements; reference columns are kept.

    Returns
    -------
    pd.DataFrame
        A copy with a boolean "ship_contaminated" column; `attrs["ais"]`
        gains the RMSE used and the number of contaminated measurements.
    """
    if "ship_flag" not in pal:
        raise KeyError("No ship_flag: run flag_ships first.")
    if action not in ("flag", "remove"):
        raise ValueError("action must be 'flag' or 'remove'.")
    reference, ref_col = _reference_column(pal, reference)
    wcol = wind if wind.startswith("U_") else f"U_{wind}"
    if wcol not in pal:
        raise KeyError(f"No {wcol!r} column.")
    near = _within(pal, reference, max_distance_km)
    err = (pal[wcol] - pal[ref_col]).abs()
    clean = near & ~pal["ship_flag"].to_numpy(bool)
    rmse = wind_metrics(pal.loc[clean, wcol], pal.loc[clean, ref_col])["rmse"]
    contaminated = pal["ship_flag"].to_numpy(bool) & (err > rmse).to_numpy()
    out = pal.copy()
    out["ship_contaminated"] = contaminated
    if action == "remove" and contaminated.any():
        ref_cols = {
            v["column"] for v in pal.attrs.get("references", {}).values()
        }
        cols = [
            c
            for c in out
            if c.startswith(("f_", "spl_", "U_")) and c not in ref_cols
        ]
        out.loc[contaminated, cols] = np.nan
    out.attrs = dict(pal.attrs)
    ais = dict(out.attrs.get("ais", {}))
    ais.update(
        {
            "deviation_wind": wcol,
            "deviation_reference": reference,
            "deviation_rmse": rmse,
            "n_contaminated": int(contaminated.sum()),
            "deviation_action": action,
        }
    )
    out.attrs["ais"] = ais
    return out


# --- Uncertainty --------------------------------------------------------------


MIN_BOOTSTRAP_CYCLES = 10


def _refit(func, x, y, p0):
    with warnings.catch_warnings():
        # only the coefficients are used, not their covariance
        warnings.simplefilter("ignore", OptimizeWarning)
        popt, _ = curve_fit(
            lambda s, *p, f=func: f(s, *p), x, y, p0=p0, maxfev=20000
        )
    return popt


def _binned_rmse(wind, resid, bin_width, min_per_bin):
    """RMSE of residuals in wind-speed bins of `bin_width`, merging
    neighbouring bins until each holds at least `min_per_bin` values."""
    ok = np.isfinite(wind) & np.isfinite(resid)
    wind, resid = wind[ok], resid[ok]
    if len(wind) == 0:
        return pd.DataFrame(columns=["wind", "sigma_cv", "n", "lo", "hi"])
    edges = np.arange(
        np.floor(wind.min() / bin_width) * bin_width,
        wind.max() + bin_width,
        bin_width,
    )
    which = np.clip(np.digitize(wind, edges) - 1, 0, len(edges) - 2)
    groups, current = [], []
    for k in range(len(edges) - 1):
        current.append(k)
        if np.isin(which, current).sum() >= min_per_bin:
            groups.append(current)
            current = []
    if current:  # leftover high-wind bins join the last group
        if groups:
            groups[-1] += current
        else:
            groups.append(current)
    rows = []
    for g in groups:
        m = np.isin(which, g)
        if not m.any():
            continue
        rows.append(
            {
                "wind": float(wind[m].mean()),
                "sigma_cv": float(np.sqrt(np.mean(resid[m] ** 2))),
                "n": int(m.sum()),
                "lo": float(edges[g[0]]),
                "hi": float(edges[g[-1] + 1]),
            }
        )
    return pd.DataFrame(rows)


def _reference_sigma(pal, sel, reference_sigma, sigma_col):
    """Reference error: given, or RMS of the table's over `sel`."""
    if reference_sigma is not None:
        return float(reference_sigma)
    if sigma_col not in pal:
        return 0.0
    sig = pal[sigma_col].to_numpy(float)[sel]
    sig = sig[np.isfinite(sig)]
    return float(np.sqrt(np.mean(sig**2))) if len(sig) else 0.0


def wind_uncertainty(
    pal,
    models=None,
    reference_sigma=None,
    n_bootstrap=100,
    bin_width=2.0,
    min_per_bin=30,
    use=None,
    seed=0,
):
    """Uncertainty of the refitted wind models, per measurement.

    The terms are computed on the same points as the refit: by default one
    per reference value, the mean of the measurements matched to it
    (`average_to_reference`), so σ is the error of the wind on the
    reference's time scale (e.g. an hourly mean against an hourly buoy
    wind); single measurements scatter more.

    Three terms, each a column, and their combination:

    "U_<model>_fit_sigma_cv"
        Empirical, out-of-sample, wind-dependent error: the model is
        refitted on all cycles but one and predicts the cycle left out
        (leave-one-cycle-out; measurements minutes apart are not
        independent, so whole cycles are left out).  The RMSE of these
        residuals in wind-speed bins gives σ(U), applied to each
        measurement from its refitted wind.  It contains every error
        source together (model shape, depth correction, acoustic noise,
        ships, rain) but also the reference's own error.
    "U_<model>_fit_sigma_coef"
        Coefficient uncertainty: spread of the refitted wind over
        `n_bootstrap` refits on cycles resampled with replacement.  Small
        inside the wind range used for fitting; large where the refit is
        extrapolated.
    "U_<model>_fit_sigma_ref"
        The reference's own error (e.g. the buoy's wind accuracy, plus any
        float-buoy distance effect): `reference_sigma` if given; otherwise
        the root mean square of the reference table's wind_speed_sigma
        over the measurements used for fitting (the ones σ_cv is computed
        from); 0 if neither.
    "U_<model>_fit_sigma"
        Combined: sqrt(max(sigma_cv² - sigma_ref², 0) + sigma_coef²).  The
        reference error is removed from the empirical term, which contains
        it, rather than added.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `fit_wind_models` (needs "cycle").
    models : list of str, optional
        By default every refitted model.
    reference_sigma : float, optional
        Error of the reference wind (m s-1).  By default taken from the
        reference table's wind_speed_sigma column ("<reference>_sigma"),
        else 0 (the empirical term is then an upper bound).
    n_bootstrap : int, optional
        Bootstrap refits for the coefficient term, by default 100.
    bin_width : float, optional
        Width of the wind-speed bins (m s-1), by default 2.
    min_per_bin : int, optional
        Bins are merged until each holds this many residuals, by default 30.
    use : array of bool, optional
        As for `fit_wind_models`.
    seed : int, optional
        Random seed of the bootstrap, by default 0.

    Returns
    -------
    pal : pd.DataFrame
        A copy with the four sigma columns per model.
    bins : pd.DataFrame
        The wind-dependent empirical σ: model, wind (bin mean), sigma_cv,
        n, lo, hi.
    """
    fits = pal.attrs.get("wind_fits")
    if not fits or not fits.get("coefficients"):
        raise KeyError("No refitted models: run fit_wind_models first.")
    if "cycle" not in pal:
        raise KeyError("Needs a 'cycle' column for cross-validation.")
    reference = fits["reference"]
    ref_col = f"U_{reference}"
    if models is None:
        models = list(fits["coefficients"])
    # the points compared with the reference: one per reference value
    # (averaged) unless the fit was made on single measurements
    pts = average_to_reference(
        pal,
        reference,
        fits.get("max_distance_km"),
        use=use,
        average=fits.get("average", True),
        window=fits.get("window", "1h"),
    )
    sigma_col = f"{reference}_sigma"
    rng = np.random.default_rng(seed)
    cycles_all = pts["cycle"].to_numpy()
    y = pts[ref_col].to_numpy(float)
    out = pal.copy()
    tables, info = [], {}
    for name in models:
        func, band, _ = LITERATURE_MODELS[name]
        x_all = pal[f"spl_{band}Hz"].to_numpy(float)
        u_fit = pal[f"U_{name}_fit"].to_numpy(float)
        x = pts[f"spl_{band}Hz"].to_numpy(float)
        p0 = list(fits["coefficients"][name].values())
        sel = np.isfinite(x) & np.isfinite(y)
        cycles = pd.unique(cycles_all[sel])
        sigma_ref = _reference_sigma(pts, sel, reference_sigma, sigma_col)
        sigma_cv = np.full(len(pal), np.nan)
        sigma_coef = np.full(len(pal), np.nan)
        if len(cycles) < 3:
            warnings.warn(
                f"{name}: collocations in only {len(cycles)} cycle(s); "
                "at least 3 are needed for cross-validation and bootstrap."
            )
        else:
            # leave-one-cycle-out residuals
            resid = np.full(len(pts), np.nan)
            pred_cv = np.full(len(pts), np.nan)
            for c in cycles:
                test = sel & (cycles_all == c)
                train = sel & (cycles_all != c)
                try:
                    popt = _refit(func, x[train], y[train], p0)
                except (RuntimeError, ValueError):
                    continue
                pred_cv[test] = func(x[test], *popt)
                resid[test] = pred_cv[test] - y[test]
            bins = _binned_rmse(pred_cv, resid, bin_width, min_per_bin)
            if len(bins):
                sigma_cv = np.interp(
                    u_fit, bins["wind"].to_numpy(), bins["sigma_cv"].to_numpy()
                )
                sigma_cv[~np.isfinite(u_fit)] = np.nan
                tables.append(bins.assign(model=name))
            # bootstrap over whole cycles
            if len(cycles) < MIN_BOOTSTRAP_CYCLES:
                warnings.warn(
                    f"{name}: bootstrap over only {len(cycles)} cycles; the "
                    "coefficient uncertainty is unreliable (and large) with "
                    f"fewer than {MIN_BOOTSTRAP_CYCLES}."
                )
            rows_by_cycle = {
                c: np.flatnonzero(sel & (cycles_all == c)) for c in cycles
            }
            preds = []
            for _ in range(n_bootstrap):
                pick = rng.choice(cycles, size=len(cycles), replace=True)
                idx = np.concatenate([rows_by_cycle[c] for c in pick])
                try:
                    popt = _refit(func, x[idx], y[idx], p0)
                except (RuntimeError, ValueError):
                    continue
                preds.append(func(x_all, *popt))
            if len(preds) > 1:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    sigma_coef = np.nanstd(np.array(preds), axis=0, ddof=1)
        out[f"U_{name}_fit_sigma_cv"] = sigma_cv
        out[f"U_{name}_fit_sigma_coef"] = sigma_coef
        out[f"U_{name}_fit_sigma_ref"] = np.where(
            np.isfinite(u_fit), sigma_ref, np.nan
        )
        empirical = np.clip(sigma_cv**2 - sigma_ref**2, 0, None)
        out[f"U_{name}_fit_sigma"] = np.sqrt(empirical + sigma_coef**2)
        info[name] = {
            "sigma_ref": sigma_ref,
            "n_cycles": len(cycles),
            "n_points": int(sel.sum()),
            "n_bootstrap_ok": len(preds) if len(cycles) >= 3 else 0,
            "fitted_spl_range": (
                [float(np.nanmin(x[sel])), float(np.nanmax(x[sel]))]
                if sel.any()
                else None
            ),
        }
    out.attrs = dict(pal.attrs)
    out.attrs["wind_uncertainty"] = {
        "reference": reference,
        "reference_sigma": (
            "table"
            if reference_sigma is None and sigma_col in pal
            else float(reference_sigma or 0.0)
        ),
        "n_bootstrap": n_bootstrap,
        "bin_width": bin_width,
        "min_per_bin": min_per_bin,
        "models": info,
    }
    table = (
        pd.concat(tables, ignore_index=True)[
            ["model", "wind", "sigma_cv", "n", "lo", "hi"]
        ]
        if tables
        else pd.DataFrame(
            columns=["model", "wind", "sigma_cv", "n", "lo", "hi"]
        )
    )
    return out, table
