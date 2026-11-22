"""Sanity checks on PAL spectra."""

import warnings

import numpy as np


class BandLayoutWarning(UserWarning):
    """Some spectra have bands filled in an unexpected pattern."""


class BandOrderWarning(UserWarning):
    """The spectrum does not fall with frequency as ambient noise should."""


# Spearman correlation between band frequency and median band level above
# which the band order is suspect.  Correctly ordered ambient-noise spectra
# fall with frequency (7902280: about -0.99); reversed order gives about
# +0.99 and a random order about 0.
BAND_ORDER_THRESHOLD = -0.8

# Fewer measurements than this are too few for a stable median spectrum.
BAND_ORDER_MIN_ROWS = 10


def _ranks(x):
    """Ranks (0..n-1) of x, ties broken by position."""
    return np.argsort(np.argsort(x, kind="stable"), kind="stable")


def band_order_score(pal):
    """Spearman correlation between band frequency and median band level.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal` (needs `attrs["bands_hz"]` and the band
        columns).

    Returns
    -------
    float
        Between -1 and 1; close to -1 when the levels fall steadily with
        frequency.  NaN if there are too few measurements or bands.
    """
    from .constants import band_column

    bands = pal.attrs["bands_hz"]
    levels = pal[[band_column(f) for f in bands]].to_numpy(dtype=float)
    if np.isfinite(levels).any(axis=1).sum() < BAND_ORDER_MIN_ROWS:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN bands
        median = np.nanmedian(levels, axis=0)
    ok = np.isfinite(median)
    if ok.sum() < 3:
        return np.nan
    r_f = _ranks(np.asarray(bands, dtype=float)[ok])
    r_l = _ranks(median[ok])
    return float(np.corrcoef(r_f, r_l)[0, 1])


def check_band_order(pal, threshold=BAND_ORDER_THRESHOLD):
    """Warn if the spectrum does not fall with frequency.

    Ambient sea noise at park depth is loudest in the lowest bands and
    quietest in the highest ones.  If the assigned band frequencies were
    in the wrong order (e.g. a firmware change sends the bands high to
    low), the median spectrum would no longer fall with frequency.  This
    catches reversed or scrambled orders; it cannot catch two neighbouring
    bands swapped, whose levels are too similar.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal`.
    threshold : float, optional
        Warn if the Spearman correlation between band frequency and median
        level is above this value, by default -0.8.

    Returns
    -------
    float
        The correlation (see `band_order_score`).

    Warns
    -----
    BandOrderWarning
        If the correlation is above `threshold`.
    """
    score = band_order_score(pal)
    if np.isfinite(score) and score > threshold:
        warnings.warn(
            f"WMO {pal.attrs.get('wmo')}: the median PAL spectrum does not "
            f"fall with frequency (Spearman r = {score:.2f}, expected "
            f"< {threshold}). The band frequencies may be assigned in the "
            "wrong order.",
            BandOrderWarning,
            stacklevel=2,
        )
    return score


def band_layout_counts(pal):
    """Count spectra by which bands hold a value.

    A float can send full 23-band spectra or reduced 9-band ones; in the
    GDAC aux files the 9-band spectra are stored in the 23-value array with
    the other 14 bands empty.  Each spectrum is classified as "23" (all
    bands), "9" (exactly the 9-band subset) or "other".

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal`.

    Returns
    -------
    dict
        {"23": n, "9": n, "other": n}.
    """
    from .constants import PAL_BANDS_9_HZ, band_column

    bands = list(pal.attrs["bands_hz"])
    filled = np.isfinite(pal[[band_column(f) for f in bands]].to_numpy(float))
    full = filled.all(axis=1)
    nine_mask = np.isin(bands, PAL_BANDS_9_HZ)
    if len(bands) == len(PAL_BANDS_9_HZ):
        nine = full
        full = np.zeros_like(full)
    else:
        nine = (filled == nine_mask).all(axis=1)
    return {
        "23": int(full.sum()),
        "9": int(nine.sum()),
        "other": int((~full & ~nine).sum()),
    }


def check_band_layout(pal):
    """Warn if some spectra are neither full nor the 9-band subset.

    Such spectra may mean bands were lost or shifted; their values are
    kept, but should be looked at.

    Returns
    -------
    dict
        As `band_layout_counts`.

    Warns
    -----
    BandLayoutWarning
        If any spectrum is classified "other".
    """
    counts = band_layout_counts(pal)
    if counts["other"]:
        warnings.warn(
            f"WMO {pal.attrs.get('wmo')}: {counts['other']} spectra have "
            "bands filled in a pattern that is neither all bands nor the "
            "9-band subset (63, 125, 400 Hz, 1, 2, 5, 8, 12.5, 20 kHz).",
            BandLayoutWarning,
            stacklevel=2,
        )
    return counts
