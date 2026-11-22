# Outputs

`process_float` returns a dict of tables (pandas DataFrames) and two dicts.
Which tables are present depends on the options:

```python
result = ap.process_float(7902280, reference="reference_winds.csv")
list(result)
# ['park', 'pal_raw', 'park_samples', 'profiles', 'ts_profiles', 'beta',
#  'positions', 'wind_fits', 'wind_skill', 'wind_uncertainty', 'meta', 'files']
```

| key | one row per | present |
|---|---|---|
| [`park`](#park-the-main-table) | PAL measurement at park depth | always |
| [`pal_raw`](#pal-raw) | PAL measurement, as read | always |
| [`park_samples`](#park-samples) | CTD/BGC sample during the park drift | always |
| [`profiles`](#profiles) | profile level | always |
| [`ts_profiles`](#ts-profiles) | level of the gridded T/S profiles used for β | always |
| [`mean_profile`](#mean-profile) | level of the mean T/S profile | `beta_method="mean"` |
| [`beta`](#beta) | depth (and profile) | always |
| [`positions`](#positions) | surface position fix | always |
| [`ais`](#ais) | AIS vessel position | `ais=` given |
| [`wind_fits`](#wind-fits) | wind model | `reference=` given |
| [`wind_skill`](#wind-skill) | wind estimate (and reference) | `reference=` given |
| [`wind_uncertainty`](#wind-uncertainty) | wind model and wind-speed bin | `reference=` given |
| [`meta`](#meta-and-files) | — (dict) | always |
| [`files`](#meta-and-files) | — (dict) | always |

In column names, `<f>` is a band centre frequency in Hz (e.g. `f_8000Hz`),
`<model>` a wind model (`vagle1990`, `nystuen2015`, `pensieri2015`) and
`<name>` the name of the reference (`ref` by default, or the key of the
`reference=` dict).  Times are UTC.  Argo variables keep their Argo names
and units (e.g. `TEMP` in °C, `PSAL` on the practical salinity scale,
`DOXY` in µmol kg⁻¹, `PRES` in dbar).

## park: the main table

One row per PAL measurement during the park drift (measurement code 290),
after cleaning.  This is the table to use for wind speed.

**Measurement**

| column | unit | description |
|---|---|---|
| `cycle` | | float cycle number |
| `measurement_code` | | Argo measurement code (290: park drift) |
| `datetime` | UTC | time of the measurement |
| `pressure` | dbar | pressure at the measurement |
| `spectrum_qc` | | Argo QC flag of the spectrum |
| `variance_class_1` … `variance_class_4` | | PAL classification of the variability of four bands, as sent by the float (NaN if not sent) |

**Acoustics**

| column | unit | description |
|---|---|---|
| `f_<f>Hz` | dB re 1 µPa | band level, as sent by the float, after QC and despiking; only the bands chosen with `bands=` |
| `depth` | m | depth of the measurement, from pressure (TEOS-10) |
| `beta_<f>Hz` | dB | depth correction β: how much the surface-generated sound is attenuated on its way down to the hydrophone |
| `beta_cycle_before`, `beta_cycle_after` | | the two profiles β was interpolated between (per-profile β only) |
| `beta_weight_after` | 0–1 | time weight of the later profile |
| `spl_<f>Hz` | dB re 1 µPa² Hz⁻¹ | depth-corrected spectral level: band level + β, minus the bandwidth term if the float did not normalise on board |

**Park CTD/BGC data** (one set per variable; which variables depends on
the float, or on `variables=`)

| column | description |
|---|---|
| `<PARAM>` | value interpolated linearly in time from the park samples, within the cycle, between samples at most 2 h apart (`max_gap`) |
| `<PARAM>_qc` | Argo QC flag of the nearer of the two samples used |
| `<PARAM>_data_mode` | data mode (R, A or D) of the nearer sample |

Values come from the chosen data mode only (`mode=`): with
`mode="delayed"`, a variable without delayed-mode values is NaN.

**Position**

| column | unit | description |
|---|---|---|
| `latitude`, `longitude` | °N, °E | position interpolated along the great circle between the surface fixes before and after the measurement |
| `position_gap_h` | h | time between those two fixes: the longer, the less certain the position |

**Wind speed** (published coefficients)

| column | unit | description |
|---|---|---|
| `U_vagle1990`, `U_nystuen2015`, `U_pensieri2015` | m s⁻¹ | wind speed from `spl_8000Hz`; NaN where the model is not defined |

**Reference wind** (with `reference=`)

| column | unit | description |
|---|---|---|
| `U_<name>` | m s⁻¹ | reference wind matched to the measurement: from the closest reference with a value within 1 h |
| `<name>_distance_km` | km | distance between the float and that reference |
| `<name>_dt_min` | min | reference time minus measurement time |
| `<name>_time` | UTC | time of the reference value used |
| `<name>_reference_name` | | reference used (if the table names its references) |
| `<name>_sigma` | m s⁻¹ | error of the reference value (if the table gives `wind_speed_sigma`) |
| `U_<model>_fit` | m s⁻¹ | wind speed from the refitted model, for every measurement |
| `U_<model>_fit_sigma_cv` | m s⁻¹ | empirical error (leave-one-cycle-out), depends on wind speed |
| `U_<model>_fit_sigma_coef` | m s⁻¹ | uncertainty from the fitted coefficients (bootstrap over cycles) |
| `U_<model>_fit_sigma_ref` | m s⁻¹ | the reference's own error |
| `U_<model>_fit_sigma` | m s⁻¹ | combined: √(max(σ_cv² − σ_ref², 0) + σ_coef²) |

By default the models are fitted on one point per reference value (the
average of the measurements matched to it), so σ is the error of the wind
on the reference's time scale, e.g. an hourly mean.  See
[Wind speed](wind.md).

**Ship noise** (with `ais=`)

| column | description |
|---|---|
| `ships_nearby` | number of vessels within 20 km and ±30 min; NaN where the AIS feed received nothing |
| `ship_flag` | True if at least one vessel was nearby |
| `ship_contaminated` | with a reference too: a ship nearby **and** a refitted wind error larger than the RMSE |

## pal_raw

The PAL table as read from the `Rtraj_aux` file, before any cleaning:
every measurement in the file, whatever its measurement code and QC flag,
with every band (`f_63Hz` … `f_25000Hz`; bands the float did not send are
NaN), not despiked, and with `cycle`, `measurement_code`, `datetime`,
`pressure`, `spectrum_qc` and `variance_class_*`.  Used to compare before
and after cleaning (`plot_despike`).

## park_samples

The park-drift samples of the core `Rtraj` file (measurement code 290),
before they are interpolated onto the PAL times: `cycle`,
`measurement_code`, `datetime`, and `PRES`, `TEMP`, `PSAL`, `DOXY`, … each
with `_qc` and `_data_mode`, in the chosen data mode.

## profiles

The synthetic profiles (`Sprof` file), one row per level: `cycle`,
`direction` (A ascending, D descending), `datetime`, `latitude`,
`longitude`, `position_qc`, and every variable with `_qc` and
`_data_mode`, in the chosen data mode.

## ts_profiles

The temperature and salinity profiles used for the depth correction, on a
1 dbar grid down to the park depth: `cycle`, `datetime`, `latitude`,
`longitude`, `pressure` (dbar), `temperature` (°C), `salinity`, and
`measured` (False where the value was extended above the shallowest or
below the deepest measured level).  Profiles that could not be used are
listed in `attrs["skipped"]`.

## mean_profile

With `beta_method="mean"` only: the mean T/S profile over the deployment:
`pressure` (dbar), `temperature` (°C), `salinity`, their standard
deviations across profiles (`temperature_std`, `salinity_std`), and
`n_measured`, the number of profiles measured at each level.

## beta

The depth correction at each depth: `depth` (m), `beta_<f>Hz` (dB) and
`alpha_<f>Hz`, the sound absorption coefficient at that depth (dB m⁻¹,
Francois & Garrison 1982).  With per-profile β (the default) there is one
set per profile, identified by `cycle`, `datetime`, `latitude` and
`longitude`; with `beta_method="mean"`, a single set from the mean
profile.

## positions

The float's surface positions: launch position and GPS fixes, with
`datetime`, `latitude`, `longitude`, `position_qc`, `cycle` and
`measurement_code`.

## ais

The AIS vessel positions read (within a box around the float's track):
`datetime`, `mmsi`, `latitude`, `longitude`.

## wind_fits

One row per wind model.

| column | description |
|---|---|
| `model` | wind model |
| `reference` | reference used |
| `n` | points used for fitting: reference values (averaged) or measurements |
| `n_measurements` | measurements behind those points |
| `published`, `fitted` | coefficients (dicts); `fitted` is None if the model could not be refitted |
| `rmse_published`, `rmse_fitted` | RMSE against the reference (m s⁻¹) |
| `r2_published`, `r2_fitted` | coefficient of determination, 1 − SS_res/SS_tot |
| `bias_published`, `bias_fitted` | mean of model − reference (m s⁻¹) |

## wind_skill

One row per wind estimate (`U_<model>` and `U_<model>_fit`), for all
references together and, with several references, for each one
(`reference_name`): `n`, `bias`, `rmse`, `r2`, `pearson_r2` (squared
correlation), and `slope` and `intercept` of the least-squares line
estimate = slope · reference + intercept.

## wind_uncertainty

The empirical error by wind-speed bin, for each refitted model: `model`,
`wind` (mean predicted wind in the bin, m s⁻¹), `sigma_cv` (RMSE of the
leave-one-cycle-out residuals, m s⁻¹), `n` (points in the bin), and the
bin edges `lo` and `hi`.  `U_<model>_fit_sigma_cv` in `park` is
interpolated from this table.

## meta and files

`meta` holds the hydrophone and its configuration from the `meta_aux`
file: `wmo`, `sensor`, `sensor_model`, `sensor_maker`,
`sensor_serial_no`, `band_normalization` (whether the float normalises
band levels by bandwidth), and the `CONFIG_Pal*` parameters at launch
(`launch_config`), fixed (`static_config`) and per mission
(`mission_config`).

`files` maps each downloaded file name to its local path.

## Processing metadata

Each table keeps a record of how it was made in its `attrs`, e.g.:

```python
park = result["park"]
park.attrs["wmo"], park.attrs["bands_hz"], park.attrs["bands_dropped"]
park.attrs["despike"]            # window and percentile
park.attrs["clean"]              # measurements removed, and why
park.attrs["spl"]                # whether the bandwidth term was subtracted
park.attrs["depth_correction"]   # method, profiles used, extrapolation limit
park.attrs["wind_fits"]          # reference, radius, averaging, coefficients
park.attrs["wind_uncertainty"]   # σ_ref, bootstrap, cycles and points per model
```

## Exported files

`ap.export` writes any table, or all of them, at any stage.

**CSV** (`format="csv"`, the default):

```python
ap.export(result, "wmo7902280")         # wmo7902280_park.csv, wmo7902280_beta.csv, ...
ap.export(result["park"], "park.csv")   # one table
```

- one file per table, `<path>_<table>.csv`, with the columns above;
- `meta` and `files` are not written, and neither are the `attrs`: use
  NetCDF to keep the processing record.

**NetCDF** (`format="netcdf"`):

```python
ap.export(result, "wmo7902280.nc", format="netcdf")
tables = ap.read_export("wmo7902280.nc")   # back to DataFrames, with attrs
```

- one file, with one group per table (`park`, `beta`, …) and the
  dimension `N_OBS` (one per row);
- global attributes: `title`, `argopal_version`, `wmo` and `tables` (the
  group names);
- each group's attributes are the table's `attrs` (nested values stored
  as JSON text);
- text columns are stored as character arrays, and missing integers as
  NaN.

The groups can also be opened without ArgoPAL, e.g. with xarray:

```python
import xarray as xr
park = xr.open_dataset("wmo7902280.nc", group="park")
```
