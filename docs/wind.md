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
`longitude`/`lon`, or `latitude=`/`longitude=` for one fixed buoy); it
should be 10 m wind.  Optional: `reference_name` (also `station`,
`buoy_name`, `platform`), needed when the table holds several references,
and `wind_speed_sigma` (also `wind_sigma`, `sigma`), the error of each
value.  Moving references (ships) are fine: positions are read per value.

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
