# Wind speed

## Published models

Step 8 applies the published models at 8 kHz, as in Delaigue et al. (2026,
Sect. 2.6, Eqs. 2-4), to the depth-corrected spectral level (dB re 1 µPa²
Hz⁻¹); one column per model:

| function | column | form |
|---|---|---|
| `wind_vagle1990` | `U_vagle1990` | log-quadratic |
| `wind_nystuen2015` | `U_nystuen2015` | cubic |
| `wind_pensieri2015` | `U_pensieri2015` | two regimes: linear for 30 ≤ SPL < 38 dB, quadratic for 38 ≤ SPL < 60 dB |

They need the 8 kHz band (processed by default);
`process_float(..., wind_models=[...])` picks a subset.  Each function can
also be called on its own, e.g. `argopal.wind_nystuen2015(park["spl_8000Hz"])`.

## Reference winds, refitting and uncertainty

ArgoPAL does not download reference winds: a table of buoy winds at hand
is collocated with each PAL measurement, then the models are refitted to
it, as in Delaigue et al. (2026, Sect. 2.6).  One CSV can hold every
reference:

| reference_name | datetime | latitude | longitude | wind_speed | wind_speed_sigma |
|---|---|---|---|---|---|
| DYFAMED | 2026-02-10 12:00:00 | 43.360 | 7.842 | 8.1 | 0.5 |
| LION | 2026-02-10 12:00:00 | 42.106 | 4.687 | 12.4 | 0.5 |

```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    fit_max_distance_km=40,  # radius for fitting, skill and uncertainty
)
result["wind_fits"]         # published vs refitted coefficients, RMSE, R², bias
result["wind_skill"]        # skill against all references, and each one
result["wind_uncertainty"]  # cross-validated σ by wind-speed bin
```

Columns are recognised by name: time (`datetime`, `time`, `date`, …), wind
speed (`wind_speed`, `wind`, `ff`, `u10`, …) and position (`latitude`/`lat`,
`longitude`/`lon`, or `latitude=`/`longitude=` for one fixed buoy).
Optional: `reference_name` (also `station`, `buoy_name`, `platform`),
needed when the table holds several references; `wind_speed_sigma` (also
`wind_sigma`, `sigma`), the error of each value; and `height_m` (also
`sensor_height`, `anemometer_height`), the height of the wind values.
Moving references (ships) are fine: positions are read per value.

**Reference winds must be 10 m winds.**  The wind models are calibrated
against wind at 10 m, so a reference measured lower (e.g. a buoy
anemometer at 4 m) must be converted first, or the refit would build the
height difference into the coefficients.  Give the height of the
`wind_speed` values, as a `height_m` column or as `reference_height_m=`
(one value, or one per reference: `{"DYFAMED": 4, "LION": 4}`), and
ArgoPAL converts them to 10 m with a neutral logarithmic profile
(`argopal.wind_to_10m`; `height_method="power"` for the power law of Hsu
et al. 1994).  Without a height, the winds are taken as 10 m winds and a
`ReferenceHeightWarning` names the references concerned.  For winds
already at 10 m, give `height_m` = 10.

Matching: for each reference, its value nearest in time within
`reference_max_time` (1 h) is taken; of these, the one from the reference
**closest to the float** is used (the nearer in time if two are equally
close).  The reference used, its distance and time difference are kept
(`ref_reference_name`, `ref_distance_km`, `ref_dt_min`, and `ref_sigma`
from the table).  Fitting, skill and uncertainty then use only the
measurements within `fit_max_distance_km` (40 km by default; None for all)
of the reference they were matched to, pooled across references.

A dict names the reference (`reference={"buoys": "reference_winds.csv"}`
gives `U_buoys`) or collocates several sources separately.

**Averaging around each reference value.**  The float usually measures
more often than the reference (e.g. every 5 min against an hourly buoy).
By default, for each reference value, the measurements from 30 min before
(included) to 30 min after (excluded) its time are averaged
(`average_to_reference`): one point per reference value, with the mean
spectral level of those measurements and the model winds computed from
it.  The models are fitted and assessed on these points, then applied to
every measurement; σ is then the error on the reference's time scale
(here an hourly mean).

| option | effect |
|---|---|
| `reference_window="1h"` (default) | ±30 min around each reference time |
| `reference_window="30min"`, `"3h"`, … | any width; windows wider than the reference's time step overlap |
| `reference_window="nearest"` | the measurements whose nearest reference value it is |
| `reference_average=False` | no averaging: every measurement against the reference value nearest in time |

`fit_wind_models`, `wind_skill`, `wind_uncertainty` and `plot_wind_fit`
take the same `window=` and `average=` options (by default those of the
fit).

**One calibration per reference** (the default, `calibration=
"per_reference"`): the models are refitted to each reference separately,
on the measurements selected for it (within its radius or inside its
footprint).  Those measurements use their reference's coefficients;
every other measurement, used in no refit, uses the coefficients of the
closest reference (`ref_closest`).  `ref_calibration` says whose
coefficients each measurement uses, and `ref_used` which measurements the
refit used.  A reference needs float data in at least 4 cycles
(`calibration_min_cycles`) to get its own coefficients: hourly points
within one cycle are not independent, and a refit on 2–3 cycles (and its
bootstrap uncertainty) is unreliable.  A reference with fewer cycles is
still compared with the float (`reference_diagnostics`), but it gets no
coefficients, and neither its measurements nor those closest to it get a
refitted wind (`U_<model>_fit` is NaN, `ref_calibration` empty): the
refitted wind is only given where a calibration applies, the
measurements used for a reference's refit and those closest to that
reference.  This happens, for example, when the float has no
delayed-mode data yet near a reference.  If no reference has enough
cycles, one set of coefficients is fitted to all together, with a
warning.  `calibration="pooled"` fits one set of coefficients to all
references together and uses it everywhere.

```python
argopal.plot_wind_timeline(result["park"])  # the whole record, 2 months per row
```

shows the refitted wind (±1σ) through the deployment, the periods used
for each reference's refit shaded in its colour, and along the top of
each row whose coefficients apply.

Steps: `collocate_reference`, `fit_wind_models` (nonlinear least squares
from the published coefficients; adds `U_<model>_fit`), `wind_skill` /
`wind_metrics` (n, bias, RMSE, R² = 1 − SS_res/SS_tot, Pearson r², slope,
intercept; `by_reference=True` adds one set of rows per reference), and
`wind_uncertainty`, which adds four columns per refitted model:

| column | what it is |
|---|---|
| `U_<model>_fit_sigma_cv` | empirical error, out of sample and wind-dependent: the model is refitted on all cycles but one and predicts the cycle left out; the RMSE of these residuals in 2 m s⁻¹ wind bins gives σ(U).  Includes every error source, and the buoy's own error |
| `U_<model>_fit_sigma_coef` | coefficient uncertainty: spread over 100 refits on cycles resampled with replacement.  Small inside the fitted wind range, large where the refit is extrapolated (unreliable with fewer than 10 cycles) |
| `U_<model>_fit_sigma_ref` | the reference's own error: `reference_sigma` if given, else the root mean square of `wind_speed_sigma` over the measurements used for fitting, else 0 |
| `U_<model>_fit_sigma` | combined: √(max(σ_cv² − σ_ref², 0) + σ_coef²); the buoy error is removed from the empirical term, which contains it |

With AIS as well, the second ship-noise criterion of the paper is applied
(`flag_ship_deviation`): a measurement is `ship_contaminated` when a ship
is nearby **and** the refitted Nystuen wind differs from the buoy by more
than the RMSE of uncontaminated measurements; `ais_action="remove"` then
blanks only those.  Figures: `argopal.plot_wind` (published and refitted, ±1σ,
with the buoy) and `argopal.plot_wind_fit` (float vs buoy, coloured by
distance).

## Do the references agree? Diagnostics

The refit gives one set of coefficients for all references together.  If
the references disagree about how the float's sound level relates to the
wind (anemometer height, distance, exposure, data quality), that single
curve is a compromise that fits none of them well.  To check:

```python
d = argopal.reference_diagnostics(result["park"])
d["fits"]       # per reference: own fit, pooled fit, published; in-sample
                # and leave-one-cycle-out RMSE
d["transfer"]   # RMSE of each reference's coefficients on each reference
argopal.plot_reference_fits(result["park"])  # each reference's curve
```

A reference whose own curve is far from the others, or whose
coefficients do badly on the other references, is worth checking.  Its
own fit always looks better in-sample; the cross-validated RMSE
(`rmse_own_cv`, `rmse_pooled_cv`) says whether it really is.

## Selecting by wind footprint instead of distance

Near the coast, the wind 40 km from a buoy can differ from the wind at the
buoy (fetch, orography), and a buoy further away can be a better match.
Instead of a fixed radius, the measurements compared with each reference
can be chosen by **wind footprint**: where a reanalysis shows the wind
varying as at the reference.

```python
result = argopal.process_float(
    7902280,
    reference="reference_winds.csv",
    footprints="era5_u10_v10.nc",   # reanalysis over the float's area and period
    footprint_min_r=0.9,            # default
)
argopal.plot_reference_map(result["park"], footprints=result["footprints"])
```

On the map, each reference's footprint is shaded in its colour, grid cell
by grid cell: the cells where that reference's R is the highest and at
least `footprint_min_r`, i.e. where a float would be matched to it.  With
a distance selection, a circle of that radius is drawn instead.

For each reference, the correlation R between the reanalysis wind at the
reference's position and at every grid cell is computed
(`wind_footprints`); the reference's own wind is not used.  By default
the wind is correlated as given: with hourly data over the deployment,
the footprint is where the wind followed the reference's during the
deployment.  `footprint_anomaly="monthly"` correlates monthly anomalies
instead (monthly means minus each calendar month's mean), for a
reanalysis over several years: a climatological footprint.  Each float measurement is
matched to the reference whose footprint has the highest R at its
position (instead of the closest one), and fitting, skill and uncertainty
use the measurements where that R is at least `footprint_min_r`, however
far from the reference; `fit_max_distance_km` is then not used.  The
reanalysis is a NetCDF file or xarray dataset with the 10 m wind ("u10"
and "v10", as in ERA5, or a speed such as "si10"); ArgoPAL does not
download it.
