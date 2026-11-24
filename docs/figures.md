# Figures

Every step has a figure showing what it did (before → after, which data
were used).  Each function returns `(fig, axes)`:

```python
figs = argopal.plot_all(result, save_dir="figures")  # all of them, saved at 300 dpi
argopal.save_figure(figs["beta"], "beta.png")        # one figure, 300 dpi
argopal.plot_all(result, annotate=True)              # + titles and notes (mode, what was skipped)

argopal.plot_pal_spectra(result["park"])  # band levels + median spectrum
argopal.plot_despike(result["pal_raw"], result["park"], band=8000, cycle=3)
argopal.plot_park(result["park"], result["park_samples"], cycle=3)
argopal.plot_positions(result["park"], result["positions"])
argopal.plot_profiles(result["ts_profiles"])  # T/S profiles used for β
argopal.plot_beta(result["beta"], band=8000)  # β vs depth, β through time
argopal.plot_depth_correction(result["park"], beta=result["beta"])
argopal.plot_wind(result["park"])  # wind speed, each model
```

| step | figure |
|---|---|
| read / clean | `plot_pal_spectra`, `plot_despike` |
| park data | `plot_park` |
| positions | `plot_positions` |
| profiles | `plot_profiles` |
| β | `plot_beta` |
| depth correction, spectral level | `plot_depth_correction` |
| reference winds | `plot_reference_map`, `plot_reference_timeseries`, `plot_reference_fits`, `plot_wind_fit`, `plot_wind_timeline` |
| wind speed | `plot_wind` |

## Examples

Float 7902280 with the default settings (delayed mode, 8 kHz).

### PAL spectra

```{image} img/spectra.png
:alt: Band levels and median spectrum
```

### Despiking

```{image} img/despike.png
:alt: 8 kHz level before and after despiking
```

### Park data

```{image} img/park.png
:alt: Park temperature, salinity and oxygen
```

### Positions

```{image} img/positions.png
:alt: Float track and PAL measurement positions
```

### Profiles used for the depth correction

```{image} img/profiles.png
:alt: Temperature and salinity profiles
```

### β

```{image} img/beta.png
:alt: β against depth, and β at 1000 m through time
```

### Depth correction

```{image} img/depth_correction.png
:alt: Band level and depth-corrected spectral level
```

### Wind speed

```{image} img/wind.png
:alt: Wind speed from the published models
```
