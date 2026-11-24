# Figures

Each step has a figure showing what it did.  `plot_all` makes all of
them, and each function returns `(fig, axes)`:

```python
figs = argopal.plot_all(result, save_dir="figures")  # saved at 300 dpi
argopal.plot_wind(result["park"])                    # one figure
argopal.save_figure(fig, "wind.png")                 # 300 dpi
```

`annotate=True` adds titles and notes (data mode, what was skipped).

| step | function |
|---|---|
| read and clean | `plot_pal_spectra`, `plot_despike` |
| park data | `plot_park` |
| positions | `plot_positions` |
| depth correction | `plot_profiles`, `plot_beta`, `plot_depth_correction` |
| wind speed | `plot_wind` |
| reference winds | `plot_reference_map`, `plot_reference_timeseries`, `plot_wind_fit`, `plot_wind_timeline`, `plot_reference_fits` |
| ship noise | `plot_ais` |

## Examples

Float 7902280, default settings (delayed mode, 8 kHz).  The reference
wind figures are in the
[wind speed tutorial](tutorials/argopal_tutorial.ipynb).

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} PAL spectra: `plot_pal_spectra`
```{image} img/spectra.png
:alt: Band levels and median spectrum
```
:::

:::{grid-item-card} Despiking: `plot_despike`
```{image} img/despike.png
:alt: 8 kHz level before and after despiking
```
:::

:::{grid-item-card} Park data: `plot_park`
```{image} img/park.png
:alt: Park temperature, salinity and oxygen
```
:::

:::{grid-item-card} Positions: `plot_positions`
```{image} img/positions.png
:alt: Float track and PAL measurement positions
```
:::

:::{grid-item-card} Profiles used for β: `plot_profiles`
```{image} img/profiles.png
:alt: Temperature and salinity profiles
```
:::

:::{grid-item-card} β: `plot_beta`
```{image} img/beta.png
:alt: β against depth, and β at 1000 m through time
```
:::

:::{grid-item-card} Depth correction: `plot_depth_correction`
```{image} img/depth_correction.png
:alt: Band level and depth-corrected spectral level
```
:::

:::{grid-item-card} Wind speed: `plot_wind`
```{image} img/wind.png
:alt: Wind speed from the published models
```
:::

::::
