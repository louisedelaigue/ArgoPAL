# Options

Every option of `process_float`, by topic, with its default.  Only `wmo`
is required.  The same options exist on the functions that run each step
on its own (see the [API reference](api.md)).

```python
result = argopal.process_float(7902280, mode="raw", reference="buoys.csv")
```

## Data

- `wmo`: WMO number of the float.
- `mode`: Argo data mode of every CTD/BGC variable, including the T/S used
  for the depth correction.  Default **`"delayed"`**.
  - `"delayed"`: delayed-mode adjusted values only;
  - `"adjusted"`: real-time adjusted values only;
  - `"raw"`: raw values, for every profile.

  No fallback from one mode to another: values missing in the chosen mode
  are NaN.
- `bands`: frequency bands to process, in Hz.  Default **`(8000,)`**;
  e.g. `(3150, 8000)`.  The wind models need 8000.
- `variables`: CTD/BGC variables to keep, e.g. `["TEMP", "PSAL", "DOXY"]`.
  Default: all that the float measured.
- `data_dir`: download folder.  Default **`"argopal_data"`**.
- `update`: download again when the GDAC copy is newer.  Default **`True`**.
- `dac`, `gdac`: data centre and GDAC address.  Default: found
  automatically, Ifremer GDAC.
- `verbose`: print progress.  Default **`True`**.

## Cleaning

- `despike`: despike the PAL spectra.  Default **`True`**.
- `despike_window`: full width of the despiking window.  Default
  **`"3h"`** (±1.5 h).
- `despike_percentile`: in each window, values above this percentile are
  discarded before averaging.  Default **`95`**.
- `bad_qc`: QC flags of park CTD/BGC samples not to use.  Default
  **`("3", "4", "9")`**.  With `mode="raw"`, remove `"3"` to keep raw
  oxygen, which Argo flags 3.
- `max_gap`: largest gap between park samples to interpolate across.
  Default **`"2h"`**.

## Depth correction

- `beta_method`: how the depth correction β is computed.  Default
  **`"interpolate"`**.
  - `"interpolate"`: β from every profile, interpolated in time to each
    measurement, so it follows the season;
  - `"mean"`: one β from the deployment-mean profile, as in Delaigue et
    al. (2026).
- `band_normalization`: whether the float divided band levels by their
  bandwidth.  Default: read from the float's `meta_aux` file.

## Wind models

- `wind_models`: published models to apply.  Default: all three,
  **`["vagle1990", "nystuen2015", "pensieri2015"]`**.

## Reference winds and calibration

Only used with `reference`.

- `reference`: reference winds, as a CSV file or table (see the
  [quick-start guide](quickstart.md#calibrate-against-reference-winds)).
  A dict names it, e.g. `{"buoys": "buoys.csv"}` gives columns `U_buoys`,
  `buoys_*`; the default name is `ref`.  Default **`None`** (no
  calibration).
- `reference_height_m`: height of the reference winds (m), if the table
  has no `height_m` column: one value, or one per reference, e.g.
  `{"DYFAMED": 3.5}`.  The winds are converted to 10 m.  Default: taken as
  10 m, with a warning.
- `reference_max_time`: largest time difference to match a measurement
  with a reference value.  Default **`"1h"`**.
- `fit_max_distance_km`: the float must be this close to a reference for
  its data to be used to calibrate it (km).  Default **`40`**; `None`
  uses all.
- `reference_window`: the measurements in this window around each
  reference value are averaged before fitting.  Default **`"1h"`** (±30
  min).
- `reference_average`: average around each reference value at all.
  Default **`True`**; `False` compares every measurement.
- `calibration`: default **`"per_reference"`**.
  - `"per_reference"`: one set of coefficients per reference, fitted on
    the data near it and used for the measurements closest to it;
  - `"pooled"`: one set of coefficients for all references, used
    everywhere.
- `calibration_min_cycles`: float cycles with data near a reference that
  it needs to get coefficients.  Default **`4`**.  The measurements
  closest to a reference without coefficients get no refitted wind.
  3 is the lowest that still gives the cross-validated uncertainty (one
  cycle left out, fitted on the others): worth trying when a reference
  with a good fit is left out, e.g. one the float passed only briefly.
- `fit_reference`: with several reference sources (a dict), the one to
  calibrate against.  Default: the first.

## Wind footprints

Select the data used to calibrate each reference by wind footprint
instead of distance.

- `footprints`: reanalysis 10 m wind (file(s) or xarray dataset, e.g.
  ERA5 `u10`/`v10`) over the float's area and period, or footprints from
  `wind_footprints`.  Default **`None`** (select by distance).
- `footprint_min_r`: correlation defining a footprint.  Default
  **`0.9`**.
- `footprint_anomaly`: default **`None`**.
  - `None`: correlate the wind as given, e.g. hourly over the deployment;
  - `"monthly"`: correlate monthly anomalies, for a reanalysis over
    several years.

## Uncertainty

Only used with `reference`.

- `uncertainty`: terms of the refitted wind's uncertainty.  Default
  **`"cv"`**.
  - `"cv"`: cross-validated error (leave one cycle out), by wind speed;
  - `"coef"`: spread of the coefficients over a bootstrap on cycles
    (unreliable with fewer than about 10 cycles);
  - `"ref"`: the reference's own error, removed from `"cv"`;
  - a list, e.g. `["cv", "ref"]`, or `"all"`.
- `reference_sigma`: error of the reference winds (m s⁻¹), for `"ref"`.
  Default: the table's `wind_speed_sigma`, else 0.
- `n_bootstrap`: bootstrap refits, for `"coef"`.  Default **`100`**.

## Ship noise

Only used with `ais`.

- `ais`: AIS vessel positions: a table, decoded CSV files or raw NMEA
  files (see [How it works](method.md#ship-noise)).  Default **`None`**.
- `ais_radius_km`: distance within which a ship counts.  Default
  **`20`**.
- `ais_window`: time window either side of each measurement.  Default
  **`"30min"`**.
- `ais_action`: default **`"flag"`**.
  - `"flag"`: add `ships_nearby` and `ship_flag`;
  - `"remove"`: blank the flagged measurements (with a reference, only
    those whose wind also deviates from the reference by more than the
    RMSE).
