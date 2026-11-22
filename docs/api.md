# API reference

All functions are available from the top-level package:

```python
import argopal as ap
```

## Whole pipeline

```{eval-rst}
.. autofunction:: argopal.process_float
```

## Argo GDAC

```{eval-rst}
.. autofunction:: argopal.fetch_pal
.. autofunction:: argopal.find_dac
.. autofunction:: argopal.list_files
.. autofunction:: argopal.download_aux
.. autofunction:: argopal.download_core
```

## Reading

```{eval-rst}
.. autofunction:: argopal.read_pal
.. autofunction:: argopal.read_pal_meta
.. autofunction:: argopal.available_variables
.. autofunction:: argopal.select_mode
.. autofunction:: argopal.read_park
.. autofunction:: argopal.read_profiles
.. autofunction:: argopal.read_positions
```

## Bands

```{eval-rst}
.. autofunction:: argopal.available_bands
.. autofunction:: argopal.select_bands
.. autofunction:: argopal.bands_for
.. autofunction:: argopal.check_band_layout
.. autofunction:: argopal.band_layout_counts
.. autofunction:: argopal.check_band_order
.. autofunction:: argopal.band_order_score
```

## Cleaning and combining

```{eval-rst}
.. autofunction:: argopal.clean_pal
.. autofunction:: argopal.despike
.. autofunction:: argopal.attach_park
.. autofunction:: argopal.attach_positions
```

## Depth correction

```{eval-rst}
.. autofunction:: argopal.gridded_profiles
.. autofunction:: argopal.mean_profile
.. autofunction:: argopal.compute_profile_betas
.. autofunction:: argopal.compute_beta
.. autofunction:: argopal.depth_correct
.. autofunction:: argopal.spectral_level
.. autofunction:: argopal.third_octave_bandwidth
.. autofunction:: argopal.beta_profile
.. autofunction:: argopal.absorption_francois_garrison
```

## Wind speed

```{eval-rst}
.. autofunction:: argopal.wind_literature
.. autofunction:: argopal.wind_vagle1990
.. autofunction:: argopal.wind_nystuen2015
.. autofunction:: argopal.wind_pensieri2015
```

## Reference winds and calibration

```{eval-rst}
.. autofunction:: argopal.read_reference
.. autofunction:: argopal.collocate_reference
.. autofunction:: argopal.fit_wind_models
.. autofunction:: argopal.wind_skill
.. autofunction:: argopal.wind_metrics
.. autofunction:: argopal.wind_uncertainty
```

## Ship noise (AIS)

```{eval-rst}
.. autofunction:: argopal.read_ais
.. autofunction:: argopal.read_ais_csv
.. autofunction:: argopal.read_ais_nmea
.. autofunction:: argopal.track_bbox
.. autofunction:: argopal.flag_ships
.. autofunction:: argopal.flag_ship_deviation
```

## Figures

```{eval-rst}
.. autofunction:: argopal.plot_all
.. autofunction:: argopal.plot_pal_spectra
.. autofunction:: argopal.plot_despike
.. autofunction:: argopal.plot_park
.. autofunction:: argopal.plot_positions
.. autofunction:: argopal.plot_profiles
.. autofunction:: argopal.plot_beta
.. autofunction:: argopal.plot_depth_correction
.. autofunction:: argopal.plot_wind
.. autofunction:: argopal.plot_wind_fit
.. autofunction:: argopal.plot_ais
.. autofunction:: argopal.save_figure
```

## Export

```{eval-rst}
.. autofunction:: argopal.export
.. autofunction:: argopal.read_export
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
