"""Wind speed from ambient noise: published empirical models.

One function per model, with its published coefficients, as applied in
Delaigue et al. (2026, Ocean Science, Sect. 2.6, Table 2, Eqs. 2-4).  All
take the depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹) and
return the 10 m wind speed U (m s-1).

| model                 | band  | form                              |
|-----------------------|-------|-----------------------------------|
| Vagle et al. 1990     | 8 kHz | log-quadratic                     |
| Nystuen et al. 2015   | 8 kHz | cubic                             |
| Pensieri et al. 2015  | 8 kHz | linear (< 38 dB) / quadratic      |

`wind_literature` applies them to a PAL table after `spectral_level`, when
the 8 kHz band was processed.
"""

import warnings

import numpy as np


def wind_vagle1990(spl_8khz, a=7.38, b=-38.70, c=21.69):
    """Vagle et al. (1990): log-quadratic model at 8 kHz.

    Equation 2 of Delaigue et al. (2026):

    SPL = 21.69 + 38.70 x - 7.38 x² with x = log10(U), solved for the
    lower root:

        U = 10 ** ((38.70 - sqrt(38.70² - 4 · 7.38 · (SPL - 21.69)))
                   / (2 · 7.38))

    NaN where the root does not exist (SPL > ≈ 72.4 dB).

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    a, b, c : float, optional
        Coefficients of SPL = c - b x - a x²; by default the published
        values (7.38, -38.70, 21.69).  Refitted values come from
        `fit_wind_models`.
    """
    spl = np.asarray(spl_8khz, dtype=float)
    disc = b**2 - 4 * a * (spl - c)
    with np.errstate(invalid="ignore"):
        root = np.sqrt(np.where(disc >= 0, disc, np.nan))
    return 10 ** ((-b - root) / (2 * a))


def wind_nystuen2015(spl_8khz, a3=0.0005, a2=-0.0310, a1=0.4904, a0=2.0871):
    """Nystuen et al. (2015): cubic model at 8 kHz.

    Equation 3 of Delaigue et al. (2026):

        U = a3 SPL³ + a2 SPL² + a1 SPL + a0
          = 0.0005 SPL³ - 0.0310 SPL² + 0.4904 SPL + 2.0871   (published)

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    a3, a2, a1, a0 : float, optional
        Coefficients; by default the published values.
    """
    spl = np.asarray(spl_8khz, dtype=float)
    return a3 * spl**3 + a2 * spl**2 + a1 * spl + a0


def wind_pensieri2015(
    spl_8khz,
    a1=0.1458,
    b1=-3.146,
    a2=0.044642,
    b2=-3.2917,
    c2=63.016,
    valid_range=(30.0, 60.0),
):
    """Pensieri et al. (2015): two-regime model at 8 kHz.

    Equation 4 of Delaigue et al. (2026):

        U = a1 SPL + b1                 30 <= SPL < 38 dB  (low regime)
        U = a2 SPL² + b2 SPL + c2       38 <= SPL < 60 dB  (high regime)

    The regimes meet at 38 dB (2.39 m s-1), as in the paper (Eq. 4).  The
    outer limits, 30 and 60 dB, come from the code of Delaigue et al.
    (2026) that produced its results; the paper's Eq. 4 leaves the regimes
    open-ended (`valid_range=None`).

    with published coefficients 0.1458, -3.146 and 0.044642, -3.2917,
    63.016.

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    a1, b1, a2, b2, c2 : float, optional
        Coefficients; by default the published values.
    valid_range : (float, float) or None, optional
        Lower limit of the low regime and upper limit of the high regime
        (dB), NaN outside; by default (30, 60).  None leaves both regimes
        open-ended, as written in the paper.
    """
    spl = np.asarray(spl_8khz, dtype=float)
    u = np.where(
        spl < 38,
        a1 * spl + b1,
        a2 * spl**2 + b2 * spl + c2,
    )
    if valid_range is not None:
        lo, hi = valid_range
        u = np.where((spl >= lo) & (spl < hi), u, np.nan)
    return np.where(np.isfinite(spl), u, np.nan)


# Published coefficients, in the order of each function's arguments
COEFFICIENTS = {
    "vagle1990": {"a": 7.38, "b": -38.70, "c": 21.69},
    "nystuen2015": {"a3": 0.0005, "a2": -0.0310, "a1": 0.4904, "a0": 2.0871},
    "pensieri2015": {
        "a1": 0.1458,
        "b1": -3.146,
        "a2": 0.044642,
        "b2": -3.2917,
        "c2": 63.016,
    },
}

# name -> (function, band in Hz, reference)
LITERATURE_MODELS = {
    "vagle1990": (wind_vagle1990, 8000, "Vagle et al. (1990)"),
    "nystuen2015": (wind_nystuen2015, 8000, "Nystuen et al. (2015)"),
    "pensieri2015": (wind_pensieri2015, 8000, "Pensieri et al. (2015)"),
}


def wind_literature(pal, models=None):
    """Wind speed from the literature models, where their band is available.

    Adds a column "U_<model>" (m s-1) per model, e.g. "U_nystuen2015".

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `spectral_level` (needs "spl_<f>Hz").
    models : list of str, optional
        Models to apply (keys of `LITERATURE_MODELS`); by default all whose
        band is available, silently skipping the others.  If models are
        named explicitly, those that cannot run raise a warning.

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with the wind columns; `attrs["wind_models"]` lists
        what was applied (band, reference) and what was skipped.
    """
    explicit = models is not None
    if models is None:
        models = list(LITERATURE_MODELS)
    unknown = [m for m in models if m not in LITERATURE_MODELS]
    if unknown:
        raise ValueError(
            f"Unknown model(s) {unknown}; available: "
            f"{', '.join(LITERATURE_MODELS)}."
        )
    out = pal.copy()
    applied, skipped = {}, {}
    for name in models:
        func, band, ref = LITERATURE_MODELS[name]
        col = f"spl_{band}Hz"
        if col not in pal:
            skipped[name] = f"{band} Hz band not processed or without data"
            continue
        out[f"U_{name}"] = func(pal[col].to_numpy(float))
        applied[name] = {"band_hz": band, "reference": ref}
    if explicit and skipped:
        warnings.warn(
            "Not applied: "
            + "; ".join(f"{m} ({why})" for m, why in skipped.items()),
            stacklevel=2,
        )
    out.attrs = dict(pal.attrs)
    out.attrs["wind_models"] = {"applied": applied, "skipped": skipped}
    return out
