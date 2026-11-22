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


def wind_vagle1990(spl_8khz):
    """Vagle et al. (1990): log-quadratic model at 8 kHz (Eq. 2).

    SPL = 21.69 + 38.70 x - 7.38 x² with x = log10(U), solved for the
    lower root:

        U = 10 ** ((38.70 - sqrt(38.70² - 4 · 7.38 · (SPL - 21.69)))
                   / (2 · 7.38))

    NaN where the root does not exist (SPL > ≈ 72.4 dB).

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    """
    a, b, c = 7.38, -38.70, 21.69
    spl = np.asarray(spl_8khz, dtype=float)
    disc = b**2 - 4 * a * (spl - c)
    with np.errstate(invalid="ignore"):
        root = np.sqrt(np.where(disc >= 0, disc, np.nan))
    return 10 ** ((-b - root) / (2 * a))


def wind_nystuen2015(spl_8khz):
    """Nystuen et al. (2015): cubic model at 8 kHz (Eq. 3).

        U = 0.0005 SPL³ - 0.0310 SPL² + 0.4904 SPL + 2.0871

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    """
    spl = np.asarray(spl_8khz, dtype=float)
    return 0.0005 * spl**3 - 0.0310 * spl**2 + 0.4904 * spl + 2.0871


def wind_pensieri2015(spl_8khz, valid_range=(30.0, 60.0)):
    """Pensieri et al. (2015): two-regime model at 8 kHz (Eq. 4).

        U = 0.1458 SPL - 3.146                        SPL < 38 dB
        U = 0.044642 SPL² - 3.2917 SPL + 63.016       SPL >= 38 dB

    Parameters
    ----------
    spl_8khz : array_like
        Depth-corrected spectral level at 8 kHz (dB re 1 µPa² Hz⁻¹).
    valid_range : (float, float) or None, optional
        SPL range (dB) where the model is applied, NaN outside; by default
        30 <= SPL < 60 dB, as in the code of Delaigue et al. (2026).  None
        applies it everywhere.
    """
    spl = np.asarray(spl_8khz, dtype=float)
    u = np.where(
        spl < 38,
        0.1458 * spl - 3.146,
        0.044642 * spl**2 - 3.2917 * spl + 63.016,
    )
    if valid_range is not None:
        lo, hi = valid_range
        u = np.where((spl >= lo) & (spl < hi), u, np.nan)
    return np.where(np.isfinite(spl), u, np.nan)


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
