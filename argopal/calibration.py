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


def _within(pal, reference, max_distance_km, min_footprint_r=None):
    """Measurements close enough to their reference: within
    `max_distance_km` of it and/or inside its wind footprint (R at least
    `min_footprint_r`, see `wind_footprints`)."""
    ok = np.ones(len(pal), dtype=bool)
    dist_col = f"{reference}_distance_km"
    if max_distance_km is not None and dist_col in pal:
        ok &= (pal[dist_col] <= max_distance_km).to_numpy()
    if min_footprint_r is not None:
        r_col = f"{reference}_footprint_r"
        if r_col not in pal:
            raise KeyError(
                f"No {r_col!r}: collocate the reference with footprints= "
                "(see wind_footprints) to select by footprint."
            )
        ok &= (pal[r_col] >= min_footprint_r).to_numpy()
    return ok


def average_to_reference(
    pal,
    reference=None,
    max_distance_km=40.0,
    use=None,
    average=True,
    window="1h",
    min_footprint_r=None,
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
    min_footprint_r : float, optional
        Instead of (or as well as) the distance, keep the measurements
        inside the wind footprint of their reference: R at least this
        (needs a collocation with `footprints=`, see `wind_footprints`).
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
    rows_ok = _within(pal, reference, max_distance_km, min_footprint_r)
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
    ref_cols = [
        c
        for c in (
            ref_col, sigma_col, f"{reference}_latitude",
            f"{reference}_longitude", f"{reference}_calibration",
        )
        if c in sub
    ]  # fmt: skip
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
        if f"U_{name}_fit" in out and fitted.get(name):
            by_group = fitted[name]
            cal_col = f"{reference}_calibration"
            group = (
                out[cal_col].astype(object).to_numpy() if cal_col in out
                else np.full(len(out), "all", dtype=object)
            )  # fmt: skip
            u = np.full(len(out), np.nan)
            for g, c in by_group.items():
                rows = group == g
                u[rows] = func(level[rows], **c)
            out[f"U_{name}_fit"] = u
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
    min_footprint_r=None,
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
    min_footprint_r : float, optional
        Select by wind footprint (see `average_to_reference`); by default
        as for `fit_wind_models`.

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
    if min_footprint_r is None:
        min_footprint_r = pal.attrs.get("wind_fits", {}).get("min_footprint_r")
    pts = average_to_reference(
        pal,
        reference,
        max_distance_km,
        average=average,
        window=window,
        min_footprint_r=min_footprint_r,
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
    min_footprint_r=None,
    calibration="per_reference",
    min_cycles=4,
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
    min_footprint_r : float, optional
        Fit on the measurements inside the wind footprint of their
        reference (R at least this) rather than, or as well as, within
        `max_distance_km` (pass None for no distance limit).  Needs a
        collocation with `footprints=` (see `wind_footprints`).
    calibration : {"per_reference", "pooled"}, optional
        "per_reference" (default): one set of coefficients per reference,
        fitted on the measurements selected for it (within its radius or
        footprint).  Those measurements use their reference's
        coefficients; every other measurement uses those of the closest
        reference.  Where that reference has no coefficients (too few
        cycles), there is no refitted wind.  "pooled": one set of
        coefficients fitted to all references together, used everywhere.
        With a single reference, both are the same.
    min_cycles : int, optional
        With one calibration per reference, a reference needs points in at
        least this many float cycles to get its own coefficients, by
        default 4: points from fewer cycles are too few independent
        samples, and a refit on them (and its bootstrap uncertainty) is
        unreliable.  Its measurements, and those closest to it, then get no
        refitted wind.  If no reference has enough cycles, one set of
        coefficients is fitted to all references together (as "pooled"),
        with a warning.

    Returns
    -------
    pal : pd.DataFrame
        A copy with, for every measurement, "U_<model>_fit",
        "<reference>_used" (True for the measurements compared with a
        reference in the refit), "<reference>_closest" (the reference it
        belongs to: its own where used in the refit, else the closest) and
        "<reference>_calibration" (the reference whose coefficients it
        uses: the same if calibrated, else None and no refitted wind;
        "all" if pooled);
        `attrs["wind_fits"]["coefficients"]` holds them, model ->
        reference -> coefficients.
    fits : pd.DataFrame
        One row per model and reference ("all" if pooled): n points used
        (reference values, or measurements if not averaged),
        n_measurements behind them, n_applied (measurements using these
        coefficients), published and refitted coefficients, and the skill
        of both on the fitting points (rmse, r2, bias).  A reference whose
        model could not be refitted (too few points) keeps its row, with
        fitted None and NaN skill; the measurements belonging to it get no
        refitted wind.
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
        min_footprint_r=min_footprint_r,
    )
    out = pal.copy()
    # which measurements the refit uses: close enough to their reference
    # (distance or footprint), with a reference value and a float level
    used = _within(pal, reference, max_distance_km, min_footprint_r)
    used &= pal[ref_col].notna().to_numpy()
    levels = [
        f"spl_{LITERATURE_MODELS[m][1]}Hz"
        for m in models
        if f"spl_{LITERATURE_MODELS[m][1]}Hz" in pal
    ]
    if levels:
        used &= pal[levels].notna().any(axis=1).to_numpy()
    if use is not None:
        used &= np.asarray(use, dtype=bool)
    out[f"{reference}_used"] = used
    if calibration not in ("per_reference", "pooled"):
        raise ValueError(
            'calibration must be "per_reference" or "pooled"; '
            f"got {calibration!r}."
        )
    name_col = f"{reference}_reference_name"
    per_ref = calibration == "per_reference" and name_col in pts
    group = (
        pts[name_col].astype(str).to_numpy() if per_ref
        else np.full(len(pts), "all", dtype=object)
    )  # fmt: skip
    groups = sorted(set(group)) if len(group) else []
    # a reference needs enough cycles to get its own coefficients
    cyc = pts["cycle"].to_numpy() if "cycle" in pts else np.zeros(len(pts))
    has_ref = pts[ref_col].notna().to_numpy().copy()
    lv = [c for c in levels if c in pts]
    if lv:  # cycles where the float has a level, not just the reference
        has_ref &= pts[lv].notna().any(axis=1).to_numpy()
    n_cycles = {g: len(pd.unique(cyc[(group == g) & has_ref])) for g in groups}
    skipped = set()
    if per_ref:
        skipped = {g for g in groups if n_cycles[g] < min_cycles}
        if skipped and len(skipped) == len(groups):
            warnings.warn(
                f"No reference has float data in {min_cycles} cycles or "
                "more ("
                + ", ".join(f"{g}: {n_cycles[g]}" for g in groups)
                + "): one set of coefficients is fitted to all references "
                "together instead.  Lower min_cycles for one calibration "
                "per reference."
            )
            per_ref = False
            group = np.full(len(pts), "all", dtype=object)
            groups = ["all"]
            n_cycles = {"all": len(pd.unique(cyc[has_ref]))}
            skipped = set()
        elif skipped:
            warnings.warn(
                "Fewer than "
                f"{min_cycles} cycles for "
                + ", ".join(f"{g} ({n_cycles[g]})" for g in sorted(skipped))
                + ": no coefficients; the measurements matched or closest "
                "to them get no refitted wind.  Lower calibration_min_cycles "
                "(min_cycles) to calibrate them anyway."
            )
            # their measurements were not used for any refit
            if name_col in out:
                skip_rows = out[name_col].astype(str).isin(skipped).to_numpy()
                out[f"{reference}_used"] = used & ~skip_rows
    rows, coefs = [], {}
    spl_of = {}
    for name in models:
        func, band, _ = LITERATURE_MODELS[name]
        spl_col = f"spl_{band}Hz"
        if spl_col not in pal:
            warnings.warn(f"{name}: {band} Hz not processed, not refitted.")
            continue
        spl_of[name] = spl_col
        x_all = pts[spl_col].to_numpy(float)
        y_all = pts[ref_col].to_numpy(float)
        p0 = COEFFICIENTS[name]
        coefs[name] = {}
        for g in groups or ["all"]:
            in_g = group == g
            x, y = x_all[in_g], y_all[in_g]
            lit = func(x, **p0)
            sel = np.isfinite(x) & np.isfinite(y) & np.isfinite(lit)
            n_meas = pts["n_measurements"].to_numpy()[in_g][sel].sum()
            info = {
                "model": name,
                "reference": reference,
                "reference_name": g,
                "n": int(sel.sum()),
                "n_cycles": n_cycles.get(g, 0),
                "n_measurements": int(n_meas),
                "published": p0,
                "fitted": None,
                **{
                    f"{m}_{k}": np.nan
                    for m in _FIT_METRICS
                    for k in _FIT_STATS
                },
            }
            label = name if g == "all" else f"{name}, {g}"
            if g in skipped:  # too few cycles: compared, not fitted
                for k, v in wind_metrics(lit[sel], y[sel]).items():
                    if f"{k}_published" in info:
                        info[f"{k}_published"] = v
                rows.append(info)
                continue
            if sel.sum() <= len(p0):
                warnings.warn(
                    f"{label}: too few collocations ({sel.sum()}) to fit."
                )
                rows.append(info)
                continue
            try:
                popt = _refit(func, x[sel], y[sel], list(p0.values()))
            except (RuntimeError, ValueError) as e:
                warnings.warn(f"{label}: fit failed ({e}).")
                rows.append(info)
                continue
            fitted = dict(zip(p0, map(float, popt)))
            before = wind_metrics(lit[sel], y[sel])
            after = wind_metrics(func(x[sel], **fitted), y[sel])
            info["fitted"] = fitted
            for m in _FIT_METRICS:
                info[f"{m}_published"] = before[m]
                info[f"{m}_fitted"] = after[m]
            rows.append(info)
            coefs[name][g] = fitted

    # which coefficients each measurement uses
    fitted_any = sorted({g for c in coefs.values() for g in c})
    closest, cal = _calibration_of(out, reference, used, fitted_any, per_ref)
    if per_ref:
        out[f"{reference}_closest"] = closest
    out[f"{reference}_calibration"] = cal
    for name, by_group in coefs.items():
        if not by_group:
            continue
        func = LITERATURE_MODELS[name][0]
        spl = out[spl_of[name]].to_numpy(float)
        u = np.full(len(out), np.nan)
        for g, c in by_group.items():
            rows_g = (cal == g).to_numpy()
            u[rows_g] = func(spl[rows_g], **c)
        out[f"U_{name}_fit"] = u
    for row in rows:
        row["n_applied"] = int(
            (cal == row["reference_name"]).sum()
            if row["fitted"] is not None
            else 0
        )
    out.attrs = dict(pal.attrs)
    settings = {
        "reference": reference,
        "calibration": "per_reference" if per_ref else "pooled",
        "max_distance_km": max_distance_km,
        "average": average,
        "window": window,
        "min_footprint_r": min_footprint_r,
    }
    out.attrs["wind_fits"] = {**settings, "coefficients": coefs}
    fits = pd.DataFrame(rows)
    fits.attrs = dict(settings)
    return out, fits


def _calibration_of(pal, reference, used, groups, per_ref):
    """Whose coefficients each measurement uses.

    Returns (closest, calibration): the reference each measurement belongs
    to (its own reference where it was used in the refit, otherwise the
    closest reference of all), and that reference if it has coefficients,
    else None: no calibration, no refitted wind.  "all" when pooled.
    """
    closest = pd.Series(None, index=pal.index, dtype=object)
    if not per_ref:
        cal = pd.Series(
            "all" if groups else None, index=pal.index, dtype=object
        )
        return closest, cal
    name_col = f"{reference}_reference_name"
    positions = (
        pal.attrs.get("references", {}).get(reference, {}).get("positions", {})
    )
    known = [
        g for g, v in positions.items() if v.get("latitude") is not None
    ]  # every reference, calibrated or not
    if known:
        lat = pal["latitude"].to_numpy(float)
        lon = pal["longitude"].to_numpy(float)
        dist = np.stack(
            [
                _haversine(
                    lat, lon, positions[g]["latitude"],
                    positions[g]["longitude"],
                )
                for g in known
            ]
        )  # fmt: skip
        dist = np.where(np.isfinite(dist), dist, np.inf)
        ok = np.isfinite(dist).any(axis=0)
        nearest = np.asarray(known, dtype=object)[dist.argmin(axis=0)]
        closest[:] = np.where(ok, nearest, None)
    if name_col in pal:
        own = pal[name_col].astype(object)
        mine = used & own.notna().to_numpy()
        closest[mine] = own[mine]
    cal = closest.where(closest.isin(groups), None)
    return closest, cal


def _haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (
        np.sin((p2 - p1) / 2) ** 2
        + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    )
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))


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
    near = _within(
        pal,
        reference,
        max_distance_km,
        pal.attrs.get("wind_fits", {}).get("min_footprint_r"),
    )
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
    components="cv",
):
    """Uncertainty of the refitted wind models, per measurement.

    The terms are computed on the same points as the refit: by default one
    per reference value, the mean of the measurements matched to it
    (`average_to_reference`), so σ is the error of the wind on the
    reference's time scale (e.g. an hourly mean against an hourly buoy
    wind); single measurements scatter more.  With one set of coefficients
    per reference (`fit_wind_models(calibration="per_reference")`), each
    reference's terms come from its own points and apply to the
    measurements that use its coefficients.

    Three terms, each a column (only those asked for in `components`),
    and their combination:

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
        Combined from the terms computed: sigma_cv alone by default;
        sqrt(max(sigma_cv² - sigma_ref², 0) + sigma_coef²) with all three.
        The reference error is removed from the empirical term, which
        contains it, rather than added.

    The empirical term is usually enough: it is out-of-sample and already
    includes the effect of the coefficients being fitted on few cycles.
    The bootstrap term is unreliable (and can be huge where the refit is
    extrapolated) with fewer than about 10 cycles.

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
    components : str or list of str, optional
        The terms to compute and combine: "cv" (default), "coef", "ref",
        any combination (e.g. ["cv", "ref"]), or "all".  "ref" needs "cv"
        (it is removed from it).  Without "coef", no bootstrap is run.

    Returns
    -------
    pal : pd.DataFrame
        A copy with, per model, the sigma column of each term computed and
        the combined "U_<model>_fit_sigma".
    bins : pd.DataFrame
        The wind-dependent empirical σ: model, reference_name (the
        reference whose coefficients, "all" if pooled), wind (bin mean),
        sigma_cv, n, lo, hi.
    """
    fits = pal.attrs.get("wind_fits")
    if not fits or not fits.get("coefficients"):
        raise KeyError("No refitted models: run fit_wind_models first.")
    if "cycle" not in pal:
        raise KeyError("Needs a 'cycle' column for cross-validation.")
    comps = _components(components)
    do_coef = "coef" in comps
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
        min_footprint_r=fits.get("min_footprint_r"),
    )
    sigma_col = f"{reference}_sigma"
    rng = np.random.default_rng(seed)
    cycles_all = pts["cycle"].to_numpy()
    y = pts[ref_col].to_numpy(float)
    cal_col = f"{reference}_calibration"

    def groups_of(table):
        if cal_col in table:
            return table[cal_col].astype(object).to_numpy()
        return np.full(len(table), "all", dtype=object)

    pt_group, row_group = groups_of(pts), groups_of(pal)
    out = pal.copy()
    tables, info = [], {}
    for name in models:
        func, band, _ = LITERATURE_MODELS[name]
        x_all = pal[f"spl_{band}Hz"].to_numpy(float)
        u_fit = pal[f"U_{name}_fit"].to_numpy(float)
        x = pts[f"spl_{band}Hz"].to_numpy(float)
        sigma_cv = np.full(len(pal), np.nan)
        sigma_coef = np.full(len(pal), np.nan)
        sigma_ref_rows = np.full(len(pal), np.nan)
        info[name] = {}
        for g, coef in fits["coefficients"][name].items():
            label = name if g == "all" else f"{name}, {g}"
            p0 = list(coef.values())
            rows = row_group == g  # measurements using these coefficients
            sel = np.isfinite(x) & np.isfinite(y) & (pt_group == g)
            cycles = pd.unique(cycles_all[sel])
            sigma_ref = _reference_sigma(pts, sel, reference_sigma, sigma_col)
            sigma_ref_rows[rows] = sigma_ref
            preds = []
            if len(cycles) < 3:
                warnings.warn(
                    f"{label}: collocations in only {len(cycles)} cycle(s); "
                    "at least 3 are needed for cross-validation and "
                    "bootstrap."
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
                    sigma_cv[rows] = np.interp(
                        u_fit[rows],
                        bins["wind"].to_numpy(),
                        bins["sigma_cv"].to_numpy(),
                    )
                    tables.append(bins.assign(model=name, reference_name=g))
                # bootstrap over whole cycles
                if do_coef and len(cycles) < MIN_BOOTSTRAP_CYCLES:
                    warnings.warn(
                        f"{label}: bootstrap over only {len(cycles)} cycles; "
                        "the coefficient uncertainty is unreliable (and "
                        f"large) with fewer than {MIN_BOOTSTRAP_CYCLES}."
                    )
                rows_by_cycle = {
                    c: np.flatnonzero(sel & (cycles_all == c)) for c in cycles
                }
                for _ in range(n_bootstrap if do_coef else 0):
                    pick = rng.choice(cycles, size=len(cycles), replace=True)
                    idx = np.concatenate([rows_by_cycle[c] for c in pick])
                    try:
                        popt = _refit(func, x[idx], y[idx], p0)
                    except (RuntimeError, ValueError):
                        continue
                    preds.append(func(x_all[rows], *popt))
                if len(preds) > 1:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        sigma_coef[rows] = np.nanstd(
                            np.array(preds), axis=0, ddof=1
                        )
            info[name][g] = {
                "sigma_ref": sigma_ref,
                "n_cycles": len(cycles),
                "n_points": int(sel.sum()),
                "n_applied": int(rows.sum()),
                "n_bootstrap_ok": len(preds) if do_coef else None,
                "fitted_spl_range": (
                    [float(np.nanmin(x[sel])), float(np.nanmax(x[sel]))]
                    if sel.any()
                    else None
                ),
            }
        sigma_cv[~np.isfinite(u_fit)] = np.nan
        sigma_coef[~np.isfinite(u_fit)] = np.nan
        sigma_ref_rows[~np.isfinite(u_fit)] = np.nan
        terms = {"cv": sigma_cv, "coef": sigma_coef, "ref": sigma_ref_rows}
        for k in ("cv", "coef", "ref"):
            col = f"U_{name}_fit_sigma_{k}"
            if k in comps:
                out[col] = terms[k]
            elif col in out:  # from an earlier run with other components
                out = out.drop(columns=col)
        var = np.zeros(len(pal))
        if "cv" in comps:
            ref2 = sigma_ref_rows**2 if "ref" in comps else 0.0
            var = var + np.clip(sigma_cv**2 - ref2, 0, None)
        if do_coef:
            var = var + sigma_coef**2
        sig = np.sqrt(var)
        sig[~np.isfinite(u_fit)] = np.nan
        out[f"U_{name}_fit_sigma"] = sig
    out.attrs = dict(pal.attrs)
    out.attrs["wind_uncertainty"] = {
        "reference": reference,
        "reference_sigma": (
            "table"
            if reference_sigma is None and sigma_col in pal
            else float(reference_sigma or 0.0)
        ),
        "components": comps,
        "n_bootstrap": n_bootstrap if do_coef else 0,
        "bin_width": bin_width,
        "min_per_bin": min_per_bin,
        "models": info,
    }
    cols = ["model", "reference_name", "wind", "sigma_cv", "n", "lo", "hi"]
    table = (
        pd.concat(tables, ignore_index=True)[cols]
        if tables
        else pd.DataFrame(columns=cols)
    )
    return out, table


UNCERTAINTY_COMPONENTS = ("cv", "coef", "ref")


def _components(components):
    """The uncertainty terms asked for, in a fixed order."""
    if components == "all":
        return list(UNCERTAINTY_COMPONENTS)
    asked = [components] if isinstance(components, str) else list(components)
    bad = [c for c in asked if c not in UNCERTAINTY_COMPONENTS]
    if bad or not asked:
        raise ValueError(
            'components must be "cv", "coef", "ref", a list of them or '
            f'"all"; got {components!r}.'
        )
    if "ref" in asked and "cv" not in asked:
        raise ValueError(
            '"ref" is removed from the empirical term: it needs "cv" too.'
        )
    return [c for c in UNCERTAINTY_COMPONENTS if c in asked]


def reference_diagnostics(park_table, reference=None, model="nystuen2015"):
    """How well each reference agrees with the others when refitting.

    One wind model is fitted to each reference separately, and to all of
    them together (as `fit_wind_models` does), on the same points as the
    refit (same selection and averaging).  For each reference: the skill
    of the published, pooled and own coefficients, in-sample and with
    leave-one-cycle-out cross-validation; and how a reference's own
    coefficients do on the other references.  A reference that needs a
    very different curve, or whose coefficients transfer badly, is worth
    checking (anemometer height, distance, exposure, data quality).

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `collocate_reference` and `fit_wind_models`.
    reference : str, optional
        Name given to `collocate_reference`; by default the one used for
        fitting.
    model : str, optional
        Wind model, by default "nystuen2015".

    Returns
    -------
    dict
        "fits" : pd.DataFrame, one row per reference and one for "all":
            n points and cycles, median distance, the coefficients fitted
            to that reference alone, and RMSE / R² / bias (m s-1) of the
            published, pooled and own coefficients, plus the
            cross-validated RMSE of the pooled and own fits
            ("rmse_pooled_cv", "rmse_own_cv": each cycle predicted by a fit
            without it).
        "transfer" : pd.DataFrame, RMSE (m s-1) of the coefficients fitted
            to one reference (rows, and "pooled") on each reference
            (columns).
        "curves" : pd.DataFrame, the wind from each set of coefficients
            over the range of spectral levels.
    """
    fits = park_table.attrs.get("wind_fits", {})
    if reference is None:
        reference = fits.get("reference")
    reference, ref_col = _reference_column(park_table, reference)
    func, band, _ = LITERATURE_MODELS[model]
    spl_col = f"spl_{band}Hz"
    pts = average_to_reference(
        park_table,
        reference,
        fits.get("max_distance_km", 40.0),
        average=fits.get("average", True),
        window=fits.get("window", "1h"),
        min_footprint_r=fits.get("min_footprint_r"),
    )
    pts = pts.dropna(subset=[spl_col, ref_col])
    name_col = f"{reference}_reference_name"
    names = (
        pts[name_col].astype(str).to_numpy() if name_col in pts
        else np.full(len(pts), reference)
    )  # fmt: skip
    x = pts[spl_col].to_numpy(float)
    y = pts[ref_col].to_numpy(float)
    cycles = pts["cycle"].to_numpy() if "cycle" in pts else np.zeros(len(x))
    dist_col = f"{reference}_distance_km"
    dist = pts[dist_col].to_numpy(float) if dist_col in pts else None
    p0 = list(COEFFICIENTS[model].values())
    npar = len(p0)

    def refit(sel):
        if sel.sum() <= npar:
            return None
        try:
            return _refit(func, x[sel], y[sel], p0)
        except (RuntimeError, ValueError):
            return None

    every = np.ones(len(x), dtype=bool)
    pooled = refit(every)
    groups = ["all"] + sorted(set(names))
    own = {g: refit(every if g == "all" else names == g) for g in groups}

    def metrics(sel, coef):
        if coef is None or not sel.any():
            return {"rmse": np.nan, "r2": np.nan, "bias": np.nan}
        m = wind_metrics(func(x[sel], *coef), y[sel])
        return {k: m[k] for k in ("rmse", "r2", "bias")}

    def cv_rmse(sel, pool):
        """Leave-one-cycle-out: each cycle of `sel` predicted by a fit
        without it (on `sel` only, or on everything if `pool`)."""
        pred = np.full(len(x), np.nan)
        for c in pd.unique(cycles[sel]):
            test = sel & (cycles == c)
            train = (every if pool else sel) & (cycles != c)
            coef = refit(train)
            if coef is not None:
                pred[test] = func(x[test], *coef)
        ok = sel & np.isfinite(pred)
        if not ok.any():
            return np.nan
        return float(np.sqrt(np.mean((pred[ok] - y[ok]) ** 2)))

    rows = []
    for g in groups:
        sel = every if g == "all" else names == g
        row = {
            "reference": g,
            "n": int(sel.sum()),
            "n_cycles": len(pd.unique(cycles[sel])),
            "median_distance_km": (
                float(np.nanmedian(dist[sel])) if dist is not None else np.nan
            ),
            "coefficients": (
                dict(zip(COEFFICIENTS[model], map(float, own[g])))
                if own[g] is not None
                else None
            ),
        }
        for label, coef in (
            ("published", p0),
            ("pooled", pooled),
            ("own", own[g]),
        ):
            for k, v in metrics(sel, coef).items():
                row[f"{k}_{label}"] = v
        row["rmse_pooled_cv"] = cv_rmse(sel, pool=True)
        row["rmse_own_cv"] = cv_rmse(sel, pool=False)
        rows.append(row)
    table = pd.DataFrame(rows)
    sources = {"pooled": pooled, **{g: own[g] for g in groups[1:]}}
    transfer = pd.DataFrame(
        {
            g: {s: metrics(names == g, c)["rmse"] for s, c in sources.items()}
            for g in groups[1:]
        }
    )
    spl = np.linspace(np.nanmin(x), np.nanmax(x), 50) if len(x) else []
    curves = pd.DataFrame({"spl": spl, "published": func(spl, *p0)})
    for s, c in sources.items():
        if c is not None:
            curves[s] = func(spl, *c)
    meta = {"reference": reference, "model": model}
    for t in (table, transfer, curves):
        t.attrs = dict(meta)
    return {"fits": table, "transfer": transfer, "curves": curves}
