# Quick-start guide

All examples assume:

```python
import argopal
```

Each section shows the main choices; [Options](options.md) lists all of
them, and [How it works](method.md) explains them.

## Process a float

`process_float` does everything, from the GDAC download to wind speed:

```python
result = argopal.process_float(7902280)

result["park"]       # one row per acoustic measurement: band levels,
                     # depth correction, spectral level, wind speed,
                     # position and park CTD/BGC data
result["profiles"]   # the float's profiles
```

Files are downloaded into `argopal_data/<wmo>/`, and downloaded again only
when the GDAC copy is newer.  See [Outputs](outputs.md) for every table
and column.

## Choose the data

The Argo data mode applies to every CTD/BGC variable, including the
temperature and salinity used for the depth correction.  There is no
fallback from one mode to another.

::::{tab-set}
:::{tab-item} Delayed mode (default)
```python
result = argopal.process_float(7902280, mode="delayed")
```
Quality-controlled data.  Recent profiles are usually not in delayed
mode yet, so the wind stops where delayed-mode data stop.
:::
:::{tab-item} Raw
```python
result = argopal.process_float(7902280, mode="raw")
```
The whole record, with unadjusted data.
:::
:::{tab-item} Adjusted
```python
result = argopal.process_float(7902280, mode="adjusted")
```
Real-time adjusted data only.
:::
::::

The 8 kHz band is processed by default, as the wind models need it.  Other
bands can be added, e.g. `bands=(3150, 8000)`.

## Calibrate against reference winds

The published wind models are biased on floats.  With reference winds
(e.g. moored buoys), ArgoPAL refits them:

```python
result = argopal.process_float(7902280, reference="reference_winds.csv")

result["park"]["U_nystuen2015_fit"]        # refitted wind
result["park"]["U_nystuen2015_fit_sigma"]  # its uncertainty (1σ)
result["wind_fits"]                        # coefficients and skill
```

**The reference table.**  One CSV holds every reference, one row per
value:

| reference_name | datetime | latitude | longitude | wind_speed | height_m |
|---|---|---|---|---|---|
| DYFAMED | 2026-02-10 12:00 | 43.360 | 7.842 | 8.1 | 3.5 |
| LION | 2026-02-10 12:00 | 42.106 | 4.687 | 12.4 | 3.5 |

The wind models are calibrated for 10 m winds: `height_m`, the height of
the wind values, is used to convert them.  Without it, the winds are taken
as 10 m winds, with a warning.  `wind_speed_sigma`, the error of each
value, is optional.

**Which float data calibrate each reference.**

::::{tab-set}
:::{tab-item} By distance (default)
```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    fit_max_distance_km=40,   # default
)
```
The float data within 40 km of a reference calibrate it.
:::
:::{tab-item} By wind footprint
```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    footprints="era5_u10_v10.nc",  # reanalysis 10 m wind
    footprint_min_r=0.9,           # default
)
```
The float data where the reanalysis wind follows the wind at a reference
(correlation R ≥ 0.9) calibrate it, however far.  Useful near the coast,
where the wind can change over a few tens of kilometres.  The reanalysis
(e.g. ERA5 `u10`, `v10`) must cover the float's area and period.
:::
::::

**One calibration per reference, or one for all.**

::::{tab-set}
:::{tab-item} Per reference (default)
```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    calibration="per_reference",  # default
    calibration_min_cycles=4,     # default
)
```
Each reference gets its own coefficients, used for the measurements
closest to it.  A reference needs float data in at least 4 cycles; the
measurements closest to a reference with fewer get no refitted wind.
:::
:::{tab-item} Pooled
```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    calibration="pooled",
)
```
One set of coefficients, fitted on all references together and used
everywhere.
:::
::::

The float measures more often than most references: the measurements
within ±30 min of each reference value are averaged before fitting
(`reference_window="1h"`; any width, e.g. `"3h"`).

## Choose the uncertainty

::::{tab-set}
:::{tab-item} Cross-validated (default)
```python
result = argopal.process_float(
    7902280, reference="reference_winds.csv", uncertainty="cv"
)
```
The error of the refitted wind when predicting a cycle left out of the
fit, by wind speed.  It includes every source of error.
:::
:::{tab-item} All terms
```python
result = argopal.process_float(
    7902280, reference="reference_winds.csv", uncertainty="all"
)
```
Adds the spread of the coefficients (bootstrap over cycles; needs about
10 cycles) and removes the reference's own error from the
cross-validated term.  Any combination works, e.g. `["cv", "ref"]`.
:::
::::

## Flag ship noise

With AIS vessel positions (a table, decoded CSV or raw NMEA files):

::::{tab-set}
:::{tab-item} Flag (default)
```python
result = argopal.process_float(7902280, ais="AIS/messages-*.gz")
result["park"]["ship_flag"]  # a ship within 20 km and ±30 min
```
:::
:::{tab-item} Remove
```python
result = argopal.process_float(
    7902280, ais="AIS/messages-*.gz", ais_action="remove"
)
```
Blanks the flagged measurements; with a reference too, only those whose
wind also deviates from the reference.
:::
::::

## Figures

```python
figs = argopal.plot_all(result, save_dir="figures")  # one figure per step
argopal.plot_wind_timeline(result["park"])            # the whole record
argopal.plot_reference_map(result["park"])            # float and references
```

See [Figures](figures.md).

## Export

::::{tab-set}
:::{tab-item} NetCDF
```python
argopal.export(result, "wmo7902280.nc", format="netcdf")
result = argopal.read_export("wmo7902280.nc")  # read back
```
One file, one group per table, with the processing record.
:::
:::{tab-item} CSV
```python
argopal.export(result, "wmo7902280", format="csv")
```
One file per table.
:::
::::

## Find floats with PAL data

```python
floats = argopal.find_pal_floats()  # every PAL float on the GDAC
```

## Run the steps one by one

Each step of `process_float` is also a function, e.g.:

```python
pal, meta = argopal.fetch_pal(7902280)    # PAL spectra and hydrophone
pal = argopal.clean_pal(pal)              # park drift, QC, despiking
```

See [How it works](method.md) for the steps and the
[API reference](api.md) for every function.
