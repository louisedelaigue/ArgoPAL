"""Attach park-drift CTD/BGC measurements to PAL measurements."""

import numpy as np
import pandas as pd


def interpolate_in_time(
    t_target, t_source, values, max_gap, qc=None, data_mode=None
):
    """Linearly interpolate one series onto target times.

    A target time gets a value only if it lies between two source samples
    that are at most `max_gap` apart; otherwise NaN.  No extrapolation.
    QC and data mode are taken from the nearer of the two samples.

    Parameters
    ----------
    t_target, t_source : array of datetime64[ns]
        Target and source times; `t_source` must be sorted.
    values : array of float
        Source values (finite only).
    max_gap : pd.Timedelta
        Largest allowed time between the two bracketing samples.
    qc, data_mode : array of str, optional
        Source QC flags and data modes.

    Returns
    -------
    values, qc, data_mode : arrays matching `t_target`
    """
    t_target = np.asarray(t_target, dtype="datetime64[ns]").astype("i8")
    t_source = np.asarray(t_source, dtype="datetime64[ns]").astype("i8")
    n = len(t_target)
    out_v = np.full(n, np.nan)
    out_q = np.full(n, "", dtype=object)
    out_m = np.full(n, "", dtype=object)
    if len(t_source) == 0:
        return out_v, out_q.astype(str), out_m.astype(str)
    if len(t_source) == 1:
        exact = t_target == t_source[0]
        out_v[exact] = values[0]
        if qc is not None:
            out_q[exact] = qc[0]
        if data_mode is not None:
            out_m[exact] = data_mode[0]
        return out_v, out_q.astype(str), out_m.astype(str)

    # index of the sample at or after each target (right bracket)
    right = np.searchsorted(t_source, t_target, side="left")
    right = np.clip(right, 1, len(t_source) - 1)
    left = right - 1
    tl, tr = t_source[left], t_source[right]
    inside = (t_target >= tl) & (t_target <= tr)
    short = (tr - tl) <= pd.Timedelta(max_gap).value
    ok = inside & short

    span = (tr - tl).astype(float)
    w = np.where(span > 0, (t_target - tl) / np.where(span > 0, span, 1), 0)
    v = values[left] + w * (values[right] - values[left])
    out_v[ok] = v[ok]
    nearer = np.where(w <= 0.5, left, right)
    if qc is not None:
        out_q[ok] = np.asarray(qc, dtype=object)[nearer[ok]]
    if data_mode is not None:
        out_m[ok] = np.asarray(data_mode, dtype=object)[nearer[ok]]
    return out_v, out_q.astype(str), out_m.astype(str)


def attach_park(pal, park, max_gap="2h", bad_qc=("3", "4", "9")):
    """Add park CTD/BGC values to each PAL measurement.

    Each parameter of `park` (except PRES; the PAL rows have their own
    pressure) is interpolated linearly in time onto the PAL times, within
    the same cycle only, and only between two samples at most `max_gap`
    apart.  Samples with a bad QC flag are left out first.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table (from `read_pal` or `clean_pal`).
    park : pd.DataFrame
        Park measurements (from `read_park`).
    max_gap : str or pd.Timedelta, optional
        Largest allowed gap between bracketing samples, by default "2h"
        (park T/S/O2 are sampled about hourly).
    bad_qc : tuple of str, optional
        QC flags of samples not to interpolate from, by default
        ("3", "4", "9").

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with, for every park parameter, the columns
        <PARAM>, <PARAM>_qc and <PARAM>_data_mode (QC and data mode of the
        nearer of the two samples used).
    """
    max_gap = pd.Timedelta(max_gap)
    variables = [v for v in park.attrs.get("variables", []) if v in park]
    out = pal.copy()
    for var in variables:
        out[var] = np.nan
        out[var + "_qc"] = ""
        out[var + "_data_mode"] = ""

    for cycle, rows in out.groupby("cycle", sort=False):
        src_cycle = park[park["cycle"] == cycle]
        t_pal = rows["datetime"].to_numpy(dtype="datetime64[ns]")
        for var in variables:
            src = src_cycle[
                np.isfinite(src_cycle[var].to_numpy(dtype=float))
                & ~src_cycle[var + "_qc"].isin(bad_qc)
            ].sort_values("datetime")
            v, q, m = interpolate_in_time(
                t_pal,
                src["datetime"].to_numpy(dtype="datetime64[ns]"),
                src[var].to_numpy(dtype=float),
                max_gap,
                qc=src[var + "_qc"].to_numpy(),
                data_mode=src[var + "_data_mode"].to_numpy(),
            )
            out.loc[rows.index, var] = v
            out.loc[rows.index, var + "_qc"] = q
            out.loc[rows.index, var + "_data_mode"] = m

    out.attrs = dict(pal.attrs)
    out.attrs["park"] = {
        "mode": park.attrs.get("mode"),
        "variables": variables,
        "max_gap": str(max_gap),
        "bad_qc": list(bad_qc),
    }
    return out
