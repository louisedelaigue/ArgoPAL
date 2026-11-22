"""Clean PAL spectra: select measurements, drop bad QC, despike."""

import warnings

import numpy as np
import pandas as pd

from .constants import band_column


def band_columns(pal):
    """Names of the band columns of a PAL table, low to high frequency."""
    return [band_column(f) for f in pal.attrs["bands_hz"]]


DEFAULT_BANDS = (8000,)


class BandCoverageWarning(UserWarning):
    """A selected band is empty in most spectra."""


def available_bands(pal, min_coverage=0.5):
    """Bands that actually hold data in a PAL table.

    A band counts as available if it has a value in at least
    `min_coverage` of the spectra.  For example a float sending 9-band
    spectra has values only at 63, 125, 400 Hz, 1, 2, 5, 8, 12.5 and 20 kHz.

    Returns
    -------
    list of int
        Band centre frequencies (Hz), low to high.
    """
    if not len(pal):
        return []
    return [
        b
        for b in pal.attrs["bands_hz"]
        if pal[band_column(b)].notna().mean() >= min_coverage
    ]


def select_bands(pal, bands=DEFAULT_BANDS, min_coverage=0.5):
    """Keep only some frequency bands of a PAL table.

    Wind and rain retrievals (e.g. WOTAN) usually use 8 kHz, and sometimes
    3.15 kHz; processing only the bands needed is also faster.  Only bands
    that hold data are kept: a requested band that is empty in most spectra
    (see `available_bands`) is dropped, with a warning.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal` (all bands).
    bands : float or sequence of float, optional
        Band centre frequencies (Hz) to keep, by default (8000,).
    min_coverage : float, optional
        A band must have a value in at least this fraction of spectra to be
        kept, by default 0.5.

    Returns
    -------
    pd.DataFrame
        The table with only the kept band columns; `attrs["bands_hz"]`
        lists them and `attrs["bands_dropped"]` the requested bands that
        were dropped.

    Raises
    ------
    ValueError
        If a band is not one the PAL sends, or if none of the requested
        bands holds data.

    Warns
    -----
    BandCoverageWarning
        For each requested band dropped for lack of data, e.g. 3.15 kHz on
        a float sending 9-band spectra.
    """
    if np.isscalar(bands):
        bands = [bands]
    known = list(pal.attrs["bands_hz"])
    bands = [int(b) for b in bands]
    unknown = [b for b in bands if b not in known]
    if unknown:
        raise ValueError(
            f"Band(s) {unknown} Hz not sent by the PAL. Bands: "
            f"{', '.join(str(b) for b in known)} Hz."
        )
    with_data = available_bands(pal, min_coverage)
    dropped = [b for b in known if b in bands and b not in with_data]
    keep = [b for b in known if b in bands and b in with_data]
    for b in dropped:
        coverage = pal[band_column(b)].notna().mean()
        warnings.warn(
            f"WMO {pal.attrs.get('wmo')}: the {b} Hz band is not used: it "
            f"has values in only {coverage:.1%} of spectra. Bands with "
            f"data: {', '.join(str(x) for x in with_data)} Hz.",
            BandCoverageWarning,
            stacklevel=2,
        )
    if not keep:
        raise ValueError(
            f"None of the requested bands {bands} Hz holds data. Bands with "
            f"data: {', '.join(str(x) for x in with_data)} Hz."
        )
    drop = [band_column(b) for b in known if b not in keep]
    out = pal.drop(columns=drop)
    out.attrs = dict(pal.attrs)
    out.attrs["bands_hz"] = keep
    out.attrs["bands_dropped"] = dropped
    return out


def _window_trimmed_mean(levels, lo, hi, percentile, chunk_cells=5e6):
    """For each row i, mean of levels[lo[i]:hi[i]] after discarding values
    above that window's percentile, per column.

    Windows are gathered into a NaN-padded (rows, width, columns) block so
    the percentiles are computed in one call; rows are processed in chunks
    to keep the block under `chunk_cells` values.
    """
    n, n_cols = levels.shape
    result = np.full(levels.shape, np.nan)
    if n == 0:
        return result
    width = int((hi - lo).max())
    padded = np.vstack([levels, np.full((1, n_cols), np.nan)])
    step = max(1, int(chunk_cells // max(1, width * n_cols)))
    offsets = np.arange(width)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN windows
        for s in range(0, n, step):
            e = min(n, s + step)
            idx = lo[s:e, None] + offsets[None, :]
            idx = np.where(idx < hi[s:e, None], idx, n)  # n -> NaN row
            w = padded[idx]  # (rows, width, columns)
            thr = _nanpercentile_axis1(w, percentile)
            result[s:e] = np.nanmean(np.where(w <= thr, w, np.nan), axis=1)
    return result


def _nanpercentile_axis1(w, percentile):
    """np.nanpercentile(w, percentile, axis=1, keepdims=True) with linear
    interpolation (same as pandas' quantile), but much faster: sorting puts
    NaNs last, so the percentile is read from the first k valid values."""
    s = np.sort(w, axis=1)
    k = np.isfinite(w).sum(axis=1, keepdims=True)
    pos = (k - 1).clip(min=0) * (percentile / 100)
    below = np.floor(pos).astype(int)
    above = np.minimum(below + 1, (k - 1).clip(min=0))
    v_below = np.take_along_axis(s, below, axis=1)
    v_above = np.take_along_axis(s, above, axis=1)
    thr = v_below + (v_above - v_below) * (pos - below)
    return np.where(k > 0, thr, np.nan)


def despike(
    pal,
    window="3h",
    percentile=95,
    by="cycle",
):
    """Despike PAL band levels with a moving window (Delaigue et al. 2026).

    For each measurement at time t, take all measurements of the same
    cycle within t +/- window/2.  In each band, values above the window's
    `percentile`-th percentile are discarded and the rest are averaged; the
    average replaces the value at t.  This removes short loud events
    (passing ships, float noise) while keeping the wind signal.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal`.
    window : str or pd.Timedelta, optional
        Full window width, by default "3h" (i.e. +/- 1.5 h).
    percentile : float, optional
        Percentile above which values are discarded, by default 95.
    by : str, optional
        Column grouping measurements that may share a window, by default
        "cycle" (windows never span two cycles).

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with despiked band columns.
    """
    half = pd.Timedelta(window) / 2
    cols = band_columns(pal)
    out = pal.copy()
    for _, g in pal.groupby(by, sort=False):
        g = g.sort_values("datetime")
        t = g["datetime"].to_numpy(dtype="datetime64[ns]")
        levels = g[cols].to_numpy(dtype=float)
        lo = np.searchsorted(t, t - half.to_timedelta64(), side="left")
        hi = np.searchsorted(t, t + half.to_timedelta64(), side="right")
        out.loc[g.index, cols] = _window_trimmed_mean(
            levels, lo, hi, percentile
        )
    out.attrs = dict(pal.attrs)
    out.attrs["despike"] = f"window={window}, percentile={percentile}"
    return out


def clean_pal(
    pal,
    measurement_codes=(290,),
    bad_qc=("3", "4", "9"),
    despike_window="3h",
    despike_percentile=95,
    despike_spectra=True,
):
    """Clean a PAL table before further processing.

    Steps:
      1. keep only the given measurement codes (default 290 = park drift);
      2. drop spectra whose QC flag is bad (default 3, 4 = bad, 9 =
         missing);
      3. despike every band (see `despike`), unless `despike_spectra` is
         False.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal`.
    measurement_codes : tuple of int, optional
        Measurement codes to keep, by default (290,).
    bad_qc : tuple of str, optional
        Spectrum QC flags to drop, by default ("3", "4", "9").
    despike_window : str, optional
        Full despiking window, by default "3h".
    despike_percentile : float, optional
        Despiking percentile, by default 95.
    despike_spectra : bool, optional
        Whether to despike, by default True.

    Returns
    -------
    pd.DataFrame
        The cleaned table; `attrs` records what was done and how many rows
        each step removed.
    """
    attrs = dict(pal.attrs)
    n0 = len(pal)
    out = pal[pal["measurement_code"].isin(measurement_codes)]
    n1 = len(out)
    if "spectrum_qc" in out:
        out = out[~out["spectrum_qc"].isin(bad_qc)]
    n2 = len(out)
    out = out.reset_index(drop=True)
    out.attrs = attrs
    if despike_spectra:
        out = despike(
            out, window=despike_window, percentile=despike_percentile
        )
    out.attrs["clean"] = {
        "measurement_codes": list(measurement_codes),
        "bad_qc": list(bad_qc),
        "removed_measurement_code": n0 - n1,
        "removed_bad_qc": n1 - n2,
    }
    return out
