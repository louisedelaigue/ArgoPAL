"""Fixed values used across argopal."""

# Argo GDAC (Ifremer mirror).  Core files live under /dac, auxiliary files
# (including the PAL passive acoustic spectra) under /aux.
GDAC_URL = "https://data-argo.ifremer.fr"

# Data assembly centres, in the order they are searched for a WMO number.
# Coriolis first because it processes the NKE floats that carry the PAL.
DACS = (
    "coriolis",
    "aoml",
    "bodc",
    "csio",
    "csiro",
    "incois",
    "jma",
    "kiost",
    "kma",
    "meds",
    "nmdis",
)

# Nominal centre frequencies (Hz) of the 23 one-third-octave bands in
# PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM (dimension N_VALUES23).
# The GDAC files do not store these, so they are assigned here.  They are
# the standard base-10 third-octave nominal frequencies, low to high, with
# a gap between 160 and 400 Hz (no 200, 250, 315 Hz bands).
PAL_BANDS_HZ = (
    63,
    100,
    125,
    160,
    400,
    500,
    630,
    800,
    1000,
    1250,
    1600,
    2000,
    2500,
    3150,
    4000,
    5000,
    6300,
    8000,
    10000,
    12500,
    16000,
    20000,
    25000,
)

# Nominal centre frequencies (Hz) of the reduced 9-band PAL transmission
# (dimension N_VALUES9): a subset of the 23 bands above, low to high.
PAL_BANDS_9_HZ = (
    63,
    125,
    400,
    1000,
    2000,
    5000,
    8000,
    12500,
    20000,
)

# The only band layouts the PAL sends, keyed by the number of values per
# spectrum: 23 bands, or the 9-band subset.  A float sending any other
# number of values is refused by `read_pal` rather than given wrong
# frequencies.
BAND_LAYOUTS = {
    23: PAL_BANDS_HZ,
    9: PAL_BANDS_9_HZ,
}


def bands_for(n_values):
    """Band centre frequencies (Hz) for a spectrum of `n_values` values.

    Raises
    ------
    ValueError
        If no band layout is known for that number of values.
    """
    if n_values not in BAND_LAYOUTS:
        raise ValueError(
            f"The spectrum has {n_values} values per measurement, but the "
            "PAL only sends "
            f"{' or '.join(str(n) for n in sorted(BAND_LAYOUTS))} bands. "
            "Refusing to guess the band frequencies."
        )
    return BAND_LAYOUTS[n_values]


SPECTRUM_VAR = "PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM"
CLASSIFICATION_VAR = "PASSIVE_ACOUSTICS_BAND_VARIANCE_CLASSIFICATION"


def band_column(frequency_hz):
    """Column name used for one band, e.g. 8000 -> "f_8000Hz"."""
    return f"f_{int(frequency_hz)}Hz"
