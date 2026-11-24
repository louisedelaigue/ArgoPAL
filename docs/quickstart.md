# Quick-start guide

All examples assume:

```python
import argopal
```

## Process a float

`process_float` does everything, from the GDAC download to wind speed:

```python
result = argopal.process_float(7902280)

result["park"]       # one row per acoustic measurement: band levels,
                     # depth correction, spectral level, wind speed,
                     # position and park CTD/BGC data
result["profiles"]   # the float's profiles
```

Files are downloaded into `argopal_data/<wmo>/` and downloaded again only
when the GDAC copy is newer.  See [Outputs](outputs.md) for every table
and column.

## Choose the data and bands

```python
result = argopal.process_float(
    7902280,
    mode="delayed",        # or "adjusted", "raw"
    bands=(3150, 8000),    # Hz; 8000 by default
)
```

`mode` applies to every CTD/BGC variable, including the T/S used for the
depth correction, with no fallback from one mode to another.  Recent
profiles are often not in delayed mode yet: `mode="raw"` covers the whole
record.

## Calibrate against reference winds

The published wind models are biased on floats.  With reference winds
(e.g. buoys), ArgoPAL refits them:

```python
result = argopal.process_float(7902280, reference="reference_winds.csv")

result["park"]["U_nystuen2015_fit"]        # refitted wind
result["park"]["U_nystuen2015_fit_sigma"]  # its uncertainty (1σ)
result["wind_fits"]                        # coefficients and skill
```

One CSV holds every reference, one row per value:

| reference_name | datetime | latitude | longitude | wind_speed | height_m |
|---|---|---|---|---|---|
| DYFAMED | 2026-02-10 12:00 | 43.360 | 7.842 | 8.1 | 3.5 |
| LION | 2026-02-10 12:00 | 42.106 | 4.687 | 12.4 | 3.5 |

`height_m` (the anemometer height) is needed to convert the winds to
10 m; `wind_speed_sigma` (error of each value) is optional.  Each
reference is calibrated on the float data within 40 km of it, and its
coefficients are used for the measurements closest to it.  See
[Options](options.md#reference-winds-and-calibration) and
[How it works](method.md#calibration-against-reference-winds).

## Uncertainty

By default the uncertainty is the cross-validated error, `sigma_cv`.  The
coefficient and reference terms can be added:

```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    uncertainty="all",     # or "cv" (default), ["cv", "coef"], ...
)
```

## Select by wind footprint

Instead of a fixed radius around each buoy, with a reanalysis (e.g. ERA5
10 m wind over the float's area and period):

```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    footprints="era5_u10_v10.nc",
)
```

## Flag ship noise

With AIS vessel positions:

```python
result = argopal.process_float(7902280, ais="AIS/messages-*.gz")
result["park"]["ship_flag"]  # True with a ship within 20 km and ±30 min
```

`ais_action="remove"` blanks the flagged measurements instead.

## Figures

```python
figs = argopal.plot_all(result, save_dir="figures")  # one figure per step
argopal.plot_wind_timeline(result["park"])            # the whole record
```

See [Figures](figures.md).

## Export

```python
argopal.export(result, "wmo7902280.nc", format="netcdf")  # one file
argopal.export(result, "wmo7902280", format="csv")        # one CSV per table
result = argopal.read_export("wmo7902280.nc")             # read back
```

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
