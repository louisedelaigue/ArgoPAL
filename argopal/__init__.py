"""
ArgoPAL
=======

Retrieve wind speed from passive acoustic (PAL) measurements on BGC-Argo
floats, starting from a WMO number.

Code examples assume that the import convention has been followed:

  >>> import argopal

The following subfunctions are then available:

Full pipeline
-------------
    process_float
    export
    read_export

Download from the Argo GDAC
---------------------------
    find_dac
    list_files
    download_aux
    download_core
    fetch_pal

Read data
---------
    read_pal
    read_pal_meta
    read_park
    read_profiles
    read_positions
    available_variables

Process (in pipeline order)
---------------------------
    available_bands
    select_bands
    clean_pal
    despike
    attach_park
    attach_positions
    gridded_profiles
    compute_profile_betas
    depth_correct
    mean_profile, compute_beta (single-profile β, as in the paper)
    spectral_level

Ship noise (optional, AIS)
--------------------------
    read_ais
    flag_ships

Wind speed (literature coefficients)
------------------------------------
    wind_literature
    wind_vagle1990
    wind_nystuen2015
    wind_pensieri2015

Reference wind and calibration
------------------------------
    read_reference
    collocate_reference
    fit_wind_models
    wind_skill
    wind_metrics
    wind_uncertainty
    flag_ship_deviation

Acoustics
---------
    beta_profile
    absorption_francois_garrison
    third_octave_bandwidth

Figures (one per step)
----------------------
    plot_pal_spectra
    plot_despike
    plot_park
    plot_positions
    plot_profiles
    plot_beta
    plot_depth_correction
    plot_ais
    plot_wind
    plot_wind_fit
    plot_all
    save_figure

Check data
----------
    check_band_order
    band_order_score
    check_band_layout
    band_layout_counts
"""

from .acoustics import absorption_francois_garrison, beta_profile
from .ais import flag_ships, read_ais, read_ais_csv, read_ais_nmea, track_bbox
from .argo import (
    DEFAULT_MODE,
    MODES,
    available_variables,
    read_park,
    read_profiles,
    select_mode,
)
from .calibration import (
    average_to_reference,
    fit_wind_models,
    flag_ship_deviation,
    wind_metrics,
    wind_skill,
    wind_uncertainty,
)
from .checks import (
    BandLayoutWarning,
    BandOrderWarning,
    band_layout_counts,
    band_order_score,
    check_band_layout,
    check_band_order,
)
from .clean import (
    DEFAULT_BANDS,
    BandCoverageWarning,
    available_bands,
    clean_pal,
    despike,
    select_bands,
)
from .combine import attach_park
from .constants import (
    BAND_LAYOUTS,
    PAL_BANDS_9_HZ,
    PAL_BANDS_HZ,
    bands_for,
)
from .correction import (
    compute_beta,
    compute_profile_betas,
    depth_correct,
    gridded_profiles,
    mean_profile,
    spectral_level,
    third_octave_bandwidth,
)
from .export import export, read_export
from .fetch import fetch_pal
from .gdac import download_aux, download_core, find_dac, list_files
from .meta import __version__, hello
from .pipeline import process_float
from .plot import (
    map_axes,
    plot_ais,
    plot_all,
    plot_beta,
    plot_depth_correction,
    plot_despike,
    plot_pal_spectra,
    plot_park,
    plot_positions,
    plot_profiles,
    plot_reference_map,
    plot_reference_timeseries,
    plot_wind,
    plot_wind_fit,
    save_figure,
)
from .positions import attach_positions, read_positions
from .read import read_pal, read_pal_meta
from .reference import (
    ReferenceHeightWarning,
    collocate_reference,
    read_reference,
    wind_to_10m,
)
from .search import find_pal_floats, read_profile_index
from .wind import (
    COEFFICIENTS,
    LITERATURE_MODELS,
    wind_literature,
    wind_nystuen2015,
    wind_pensieri2015,
    wind_vagle1990,
)


__all__ = [
    "BAND_LAYOUTS",
    "COEFFICIENTS",
    "DEFAULT_BANDS",
    "DEFAULT_MODE",
    "LITERATURE_MODELS",
    "MODES",
    "PAL_BANDS_9_HZ",
    "PAL_BANDS_HZ",
    "BandCoverageWarning",
    "BandLayoutWarning",
    "BandOrderWarning",
    "ReferenceHeightWarning",
    "__version__",
    "absorption_francois_garrison",
    "attach_park",
    "attach_positions",
    "available_bands",
    "available_variables",
    "average_to_reference",
    "band_layout_counts",
    "band_order_score",
    "bands_for",
    "beta_profile",
    "check_band_layout",
    "check_band_order",
    "clean_pal",
    "collocate_reference",
    "compute_beta",
    "compute_profile_betas",
    "depth_correct",
    "despike",
    "download_aux",
    "download_core",
    "export",
    "fetch_pal",
    "find_dac",
    "find_pal_floats",
    "fit_wind_models",
    "flag_ship_deviation",
    "flag_ships",
    "gridded_profiles",
    "hello",
    "list_files",
    "map_axes",
    "mean_profile",
    "plot_ais",
    "plot_all",
    "plot_beta",
    "plot_depth_correction",
    "plot_despike",
    "plot_pal_spectra",
    "plot_park",
    "plot_positions",
    "plot_profiles",
    "plot_reference_map",
    "plot_reference_timeseries",
    "plot_wind",
    "plot_wind_fit",
    "process_float",
    "read_ais",
    "read_ais_csv",
    "read_ais_nmea",
    "read_export",
    "read_pal",
    "read_pal_meta",
    "read_park",
    "read_positions",
    "read_profile_index",
    "read_profiles",
    "read_reference",
    "save_figure",
    "select_bands",
    "select_mode",
    "spectral_level",
    "third_octave_bandwidth",
    "track_bbox",
    "wind_literature",
    "wind_metrics",
    "wind_nystuen2015",
    "wind_pensieri2015",
    "wind_skill",
    "wind_to_10m",
    "wind_uncertainty",
    "wind_vagle1990",
]
