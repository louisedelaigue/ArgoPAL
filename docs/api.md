# API reference

Every public function, grouped by processing step.  All are available from
the top-level package:

```python
import argopal as ap
```

## Whole pipeline

Run every step from a WMO number.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: process_float
```

## Argo GDAC

Find, list and download a float's files from the Argo GDAC.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: fetch_pal
.. autofunction:: find_dac
.. autofunction:: list_files
.. autofunction:: download_aux
.. autofunction:: download_core
```

## Reading

Read PAL spectra, the hydrophone configuration and the core Argo files, in the chosen data mode.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: read_pal
.. autofunction:: read_pal_meta
.. autofunction:: available_variables
.. autofunction:: select_mode
.. autofunction:: read_park
.. autofunction:: read_profiles
.. autofunction:: read_positions
```

## Bands

Frequency bands: which hold data, which to keep, and checks of the band layout and order.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: available_bands
.. autofunction:: select_bands
.. autofunction:: bands_for
.. autofunction:: check_band_layout
.. autofunction:: band_layout_counts
.. autofunction:: check_band_order
.. autofunction:: band_order_score
```

## Cleaning and combining

Park-drift selection, QC, despiking, and park data and positions attached to each measurement.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: clean_pal
.. autofunction:: despike
.. autofunction:: attach_park
.. autofunction:: attach_positions
```

## Depth correction

T/S profiles, the depth correction β, and spectral levels.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: gridded_profiles
.. autofunction:: mean_profile
.. autofunction:: compute_profile_betas
.. autofunction:: compute_beta
.. autofunction:: depth_correct
.. autofunction:: spectral_level
.. autofunction:: third_octave_bandwidth
.. autofunction:: beta_profile
.. autofunction:: absorption_francois_garrison
```

## Wind speed

Wind speed from the published 8 kHz models.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: wind_literature
.. autofunction:: wind_vagle1990
.. autofunction:: wind_nystuen2015
.. autofunction:: wind_pensieri2015
```

## Reference winds and calibration

Collocation with reference winds, refitting of the models, skill and uncertainty.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: read_reference
.. autofunction:: collocate_reference
.. autofunction:: average_to_reference
.. autofunction:: fit_wind_models
.. autofunction:: wind_skill
.. autofunction:: wind_metrics
.. autofunction:: wind_uncertainty
```

## Ship noise (AIS)

Read AIS vessel positions and flag measurements with a ship nearby.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: read_ais
.. autofunction:: read_ais_csv
.. autofunction:: read_ais_nmea
.. autofunction:: track_bbox
.. autofunction:: flag_ships
.. autofunction:: flag_ship_deviation
```

## Figures

One figure per processing step; each returns `(fig, axes)`.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: plot_all
.. autofunction:: plot_pal_spectra
.. autofunction:: plot_despike
.. autofunction:: plot_park
.. autofunction:: plot_positions
.. autofunction:: plot_profiles
.. autofunction:: plot_beta
.. autofunction:: plot_depth_correction
.. autofunction:: plot_wind
.. autofunction:: plot_wind_fit
.. autofunction:: plot_ais
.. autofunction:: save_figure
```

## Export

Write any table, or all of them, to CSV or NetCDF, and read them back.

```{eval-rst}
.. currentmodule:: argopal

.. autofunction:: export
.. autofunction:: read_export
```

## Constants

| name | value |
|---|---|
| `ap.PAL_BANDS_HZ` | centre frequencies (Hz) of the 23 PAL bands |
| `ap.PAL_BANDS_9_HZ` | centre frequencies (Hz) of the 9-band subset |
| `ap.BAND_LAYOUTS` | number of values per spectrum -> band frequencies |
| `ap.DEFAULT_BANDS` | bands processed by default (8 kHz) |
| `ap.MODES` | data modes: `"delayed"`, `"adjusted"`, `"raw"` |
| `ap.DEFAULT_MODE` | `"delayed"` |
| `ap.COEFFICIENTS` | published coefficients of each wind model |
| `ap.LITERATURE_MODELS` | wind model name -> (function, band, reference) |

## Warnings

```{eval-rst}
.. autoclass:: argopal.BandLayoutWarning
.. autoclass:: argopal.BandOrderWarning
.. autoclass:: argopal.BandCoverageWarning
```
