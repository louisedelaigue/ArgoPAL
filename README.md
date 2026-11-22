<p align="center">
  <img src="docs/img/argopal_logo.png" alt="ArgoPAL" width="480">
</p>

[![Tests](https://github.com/louisedelaigue/ArgoPAL/actions/workflows/run-tests.yml/badge.svg)](https://github.com/louisedelaigue/ArgoPAL/actions/workflows/run-tests.yml)
[![Coverage](.misc/coverage.svg)](.misc/coverage.txt)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Wind from passive acoustics on BGC-Argo floats.

> [!NOTE]
> ArgoPAL is under active development.  Functions, options and results
> may change between versions, and the processing has so far been tested
> on two floats only.  Check the
> [status and limitations](#status-and-limitations) before relying on its
> output, and report problems on the
> [issue tracker](https://github.com/louisedelaigue/ArgoPAL/issues).

Documentation: https://argopal.readthedocs.io

ArgoPAL starts from a float's WMO number, downloads its passive acoustic
(PAL) spectra from the Argo GDAC, and processes them into surface wind
speed.

## Install

From GitHub:

    pip install git+https://github.com/louisedelaigue/ArgoPAL.git
    pip install "argopal[ais] @ git+https://github.com/louisedelaigue/ArgoPAL.git"  # with raw NMEA AIS

From a local copy, for development:

    pip install -e .

## Import

```python
import argopal as ap
```

## Tutorial

[`examples/argopal_tutorial.ipynb`](examples/argopal_tutorial.ipynb) walks
through every step on float 7902280: processing, options, the pipeline
function by function, figures, published wind models, calibration against
three Météo-France buoys with uncertainty, AIS and export.  It downloads
the float from the GDAC; the buoy winds it uses are in `examples/data/`.
Rebuild it with `python examples/make_tutorial.py`, then run it from
`examples/`.

## Status and limitations

Version 0.1.0b1 is a beta: the processing chain is complete and tested,
but function names and options may still change before 1.0.  Changes are
listed in [`docs/versions.md`](docs/versions.md).

- **Calibration needs a reference.**  The published wind models were not
  made for floats and are biased on them (e.g. +1 to +2 m s⁻¹ for Nystuen on
  7902280).  Refitting needs reference winds (buoys) that the float
  passes close to; ArgoPAL does not download them.  With fewer than 10
  cycles near a reference, the coefficient uncertainty is unreliable, and
  a refit is not to be trusted far outside the winds seen near the
  reference.
- **Delayed mode lags.**  The default mode uses delayed-mode profiles
  only; recent profiles are not in delayed mode yet, so the depth
  correction, and the wind, stop where delayed-mode data stop.
  `mode="raw"` covers the whole record with unadjusted data.
- **8 kHz only for wind.**  The three wind models are 8 kHz models.  Other
  bands are processed and depth-corrected but not turned into wind.
- **Band layouts.**  Only 23-band spectra and the 9-band subset are
  recognised; any other layout raises an error.
- **Ship noise needs AIS supplied by the user.**  Without AIS, ship noise
  is reduced by despiking only.
- **Rain is not flagged.**  Rain raises the 8 kHz level and the wind
  estimate with it.
- **Tested floats.**  Tested on two NKE PROVOR floats with PAL
  hydrophones in the Ligurian Sea; floats elsewhere or with other
  configurations may need adjustments.

## Get the acoustic data of a float

```python
pal, meta = ap.fetch_pal(7902280)
```

This downloads the float's auxiliary files from
`https://data-argo.ifremer.fr/aux/<dac>/<wmo>/` into
`argopal_data/<wmo>/` (re-downloading only files that changed on the
server) and returns:

- `pal`: a DataFrame with one row per acoustic measurement, with columns
  `cycle`, `measurement_code`, `datetime`, `pressure` and the
  third-octave band levels `f_63Hz` … `f_25000Hz` (dB re 1 µPa; 23 bands,
  or 9 for floats sending the reduced spectrum);
- `meta`: the hydrophone details and the float's PAL configuration
  (`CONFIG_Pal*` parameters).

The steps can also be run separately:

```python
ap.find_dac(7902280)                 # "coriolis"
paths = ap.download_aux(7902280)     # aux files, incl. Rtraj_aux.nc
pal = ap.read_pal(paths["7902280_Rtraj_aux.nc"])
meta = ap.read_pal_meta(paths["7902280_meta_aux.nc"])
ap.download_core(7902280)            # meta, Rtraj, Sprof (CTD, BGC, GPS)
```

## Process a float: PAL + all other float data

```python
ap.available_variables("argopal_data/7902280")  # what the float measured, per data mode
result = ap.process_float(7902280)              # 8 kHz, mode="delayed" by default
result = ap.process_float(7902280, bands=(3150, 8000))   # several bands
result["park"]          # one row per PAL measurement: bands, β, spectral levels,
                        # position, park T/S/O2/pH
result["profiles"]      # one row per profile level (Sprof)
result["ts_profiles"]   # the T/S profiles used for the depth correction
result["beta"]          # β of every band at each depth, per profile
result["positions"]     # launch + surface GPS fixes
```

The main table is `result["park"]`, one row per acoustic measurement:

| columns | content |
|---|---|
| `cycle`, `datetime`, `pressure` | when and where the float measured |
| `f_<f>Hz` | band level as sent, after QC and despiking (dB re 1 µPa) |
| `TEMP`, `PSAL`, `DOXY`, … with `_qc`, `_data_mode` | park CTD/BGC data at the time of the measurement |
| `latitude`, `longitude`, `position_gap_h` | position between surface fixes |
| `depth`, `beta_<f>Hz` | depth (m) and depth correction β (dB) |
| `spl_<f>Hz` | depth-corrected spectral level (dB re 1 µPa² Hz⁻¹) |
| `U_<model>` | wind speed from each published model (m s⁻¹) |
| `U_<model>_fit`, `U_<model>_fit_sigma*` | with a reference: refitted wind and its uncertainty |
| `ships_nearby`, `ship_flag` | with AIS: vessels nearby |

Every table and column, with units, is described in
[Outputs](https://argopal.readthedocs.io/en/latest/outputs.html).

`process_float` runs these steps; each is also a function on its own:

| step | function(s) | adds |
|---|---|---|
| 1 | `download_aux`, `download_core` | GDAC files |
| 2 | `read_pal`, `clean_pal` | park-drift spectra, bad QC dropped, despiked (±1.5 h, > 95th percentile discarded) |
| 3 | `read_park`, `attach_park` | park CTD/BGC, linear in time within a cycle (no gaps > 2 h) |
| 4 | `read_positions`, `attach_positions` | `latitude`, `longitude` between surface fixes |
| 5 | `read_profiles`, `gridded_profiles` | every usable T/S profile in the chosen data mode, on a 1 dbar grid |
| 6 | `compute_profile_betas`, `depth_correct` | `depth`, `beta_<f>Hz` for every band, from the profiles before and after each measurement |
| 7 | `spectral_level` | `spl_<f>Hz`: depth-corrected spectral level, dB re 1 µPa²/Hz |

### Bands

`bands=` selects the frequency bands to process (default 8 kHz; for wind
and rain, e.g. WOTAN, usually 8 and/or 3.15 kHz).  Only those bands are
despiked, depth-corrected and exported; `result["pal_raw"]` keeps every
band.  A warning is raised if a band is empty in most spectra: floats
sending 9-band spectra (63, 125, 400 Hz, 1, 2, 5, 8, 12.5, 20 kHz) have no
3.15 kHz values.

### Despiking

`despike_window="3h"` and `despike_percentile=95` (defaults, as in
Delaigue et al. 2026): within ±1.5 h of each measurement, values above the
95th percentile are discarded and the rest averaged.  Both can be changed,
or despiking turned off with `despike=False`.

### Data mode

`mode=` applies to every CTD/BGC variable, including the T/S used for the
depth correction.  Each mode takes values from one source only, with no
fallback:

| mode | values used |
|---|---|
| `"delayed"` (default) | delayed-mode (D) adjusted values; NaN otherwise |
| `"adjusted"` | real-time adjusted (A) values; NaN otherwise |
| `"raw"` | raw values (kept by Argo for every data mode) |

Every variable `X` comes with `X_qc` and `X_data_mode` columns.

### Export at any stage

`ap.export` writes one table or any set of tables, from any step:

```python
ap.export(result, "wmo7902280", format="csv")         # one CSV per table
ap.export(result, "wmo7902280.nc", format="netcdf")   # one file, one group per table
ap.export(ap.clean_pal(pal), "pal_clean.csv")         # a single intermediate table
ap.read_export("wmo7902280.nc")                       # back to DataFrames
```

CSV gives one file per table, without the processing record; NetCDF keeps
everything in one file, one group per table, with each table's processing
record as attributes.  See [Outputs](https://argopal.readthedocs.io/en/latest/outputs.html).

### Depth correction

β is computed for **every profile** and interpolated linearly in time
between the profile before and the profile after each PAL measurement, so
it follows the seasonal change of the sound-speed profile (for 7902280,
β at 1000 m and 8 kHz goes from 0.94 dB in February to 1.42 dB in
August).  More than 6 days before the first or after the last usable
profile, β is NaN rather than borrowed from another season.
`process_float(..., beta_method="mean")` uses one deployment-mean profile
instead, as in Delaigue et al. (2026).

T/S profiles are held constant above their shallowest point (mixed
layer) and below their deepest point; `gridded_profiles(...,
surface="linear")` reproduces the paper's linear extrapolation to the
surface, which can overshoot (e.g. +5 °C when the two shallowest points
are close together in a strong gradient).

β follows Vagle et al. (1990) as implemented by Cauchy et al. (2018), ported
from MATLAB, with the absorption applied in the right units:
Francois–Garrison α is in dB/m, so the loss along a ray is 10^(−αL/10).
The original code used exp(−αL), which overestimates β about 4.3 times (at
1000 m, 8 kHz: 0.96 dB, not 4.0 dB). The tests check both: with the old
convention the port reproduces the published table, and with the fix it
gives 0.96 dB.

## Ship noise (optional): AIS

Ship noise can dominate the 1-10 kHz band up to ~20 km away.  Given AIS
vessel positions, `process_float` flags each PAL measurement with a ship
within 20 km of the float and ±30 min (Delaigue et al. 2026, Sect. 2.5):

```python
result = ap.process_float(7902280, ais="AIS/messages-*.gz")           # flag only
result = ap.process_float(7902280, ais=my_table, ais_action="remove")  # blank flagged values
```

Adds `ships_nearby` (vessels counted once per MMSI; NaN when the AIS feed
received nothing at all in the window) and `ship_flag`.  `ais_radius_km`
and `ais_window` change the criteria.  Figure: `ap.plot_ais(result["park"])`.

**What AIS data is needed.** Vessel positions covering the float's area
and dates, as any of:

| format | what ArgoPAL needs |
|---|---|
| pandas DataFrame | columns `datetime` (UTC), `mmsi`, `latitude`, `longitude` |
| CSV of decoded positions (`.csv`, `.csv.gz`) | a time column (`datetime`, `time`, `BaseDateTime`, `# Timestamp`, `tagblock_timestamp` = UNIX s, …), `mmsi`, latitude (`latitude`/`lat`/`y`), longitude (`longitude`/`lon`/`x`). Covers NOAA MarineCadastre, Danish Maritime Authority and the decoded files of the 2025 DYFAMED study. |
| raw NMEA 0183 (`!AIVDM`, plain or `.gz`) | NMEA 4.0 tag blocks with the receiver time, e.g. `\c:1739500274,s:2954*14\!AIVDM,...`; needs `pip install argopal[ais]` (pyais) |

Only vessel position reports (message types 1, 2, 3, 18, 19, 27) count:
base stations and aids to navigation are dropped when the CSV has a
message-type column.  Times must be the time each position was received:
do not date positions from file names, since daily AIS files often span
two calendar days.

## Wind speed: literature coefficients

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
also be called on its own, e.g. `ap.wind_nystuen2015(park["spl_8000Hz"])`.

## Reference wind, refitting and uncertainty

ArgoPAL does not download reference winds: a table of buoy winds at hand
is collocated with each PAL measurement, then the models are refitted to
it, as in Delaigue et al. (2026, Sect. 2.6).  One CSV can hold every
reference:

| reference_name | datetime | latitude | longitude | wind_speed | wind_speed_sigma |
|---|---|---|---|---|---|
| DYFAMED | 2026-02-10 12:00:00 | 43.360 | 7.842 | 8.1 | 0.5 |
| LION | 2026-02-10 12:00:00 | 42.106 | 4.687 | 12.4 | 0.5 |

```python
result = ap.process_float(
    7902280,
    reference="reference_winds.csv",
    fit_max_distance_km=40,       # radius for fitting, skill and uncertainty
)
result["wind_fits"]          # published vs refitted coefficients, RMSE, R², bias
result["wind_skill"]         # skill against all references, and each one
result["wind_uncertainty"]   # cross-validated σ by wind-speed bin
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

The float usually measures more often than the reference (e.g. every
10 min against an hourly buoy), so several measurements share one
reference value.  By default they are averaged first
(`average_to_reference`): one point per reference value, the mean spectral
level of its measurements, with the model winds computed from it.  The
models are fitted and assessed on these points, then applied to every
measurement; σ is then the error on the reference's time scale (e.g. an
hourly mean).  `reference_average=False` uses every measurement against
the reference value nearest in time.  A dict
names the reference (`reference={"buoys": "reference_winds.csv"}` gives
`U_buoys`) or collocates several sources separately.

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
blanks only those.  Figures: `ap.plot_wind` (published and refitted, ±1σ,
with the buoy) and `ap.plot_wind_fit` (float vs buoy, coloured by
distance).

## Figures: follow each step

Every step has a figure showing what it did (before → after, which data
were used).  Each function returns `(fig, axes)`:

```python
figs = ap.plot_all(result, save_dir="figures")  # all of them, saved at 300 dpi
ap.save_figure(figs["beta"], "beta.png")        # one figure, 300 dpi
ap.plot_all(result, annotate=True)              # + titles and notes (mode, what was skipped)

ap.plot_pal_spectra(result["park"])             # band levels + median spectrum
ap.plot_despike(result["pal_raw"], result["park"], band=8000, cycle=3)
ap.plot_park(result["park"], result["park_samples"], cycle=3)
ap.plot_positions(result["park"], result["positions"])
ap.plot_profiles(result["ts_profiles"])          # T/S profiles used for β
ap.plot_beta(result["beta"], band=8000)          # β vs depth, β through time
ap.plot_depth_correction(result["park"], beta=result["beta"])
ap.plot_wind(result["park"])                   # wind speed, each model
```

| step | figure |
|---|---|
| read / clean | `plot_pal_spectra`, `plot_despike` |
| park data | `plot_park` |
| positions | `plot_positions` |
| profiles | `plot_profiles` |
| β | `plot_beta` |
| depth correction, spectral level | `plot_depth_correction` |
| wind speed | `plot_wind` |

Examples for float 7902280 (default settings: delayed mode, 8 kHz) are in `docs/img/`.

## Find current version

```python
ap.hello()
```

## Tests

    pytest                            # offline tests
    pytest -m network                 # tests that download from the GDAC
    sh .misc/update_coverage.sh       # all tests with coverage; updates the badge

The coverage badge and report (`.misc/coverage.svg`, `.misc/coverage.txt`)
are made locally with pytest-cov and genbadge, from all tests including
the GDAC ones; rerun the script before committing changes to the code.

The tests run on two floats:

- **7902280** (GDAC): small slices of its aux and core files in
  `tests/data/` (rebuild: `tests/make_fixtures.py`);
- **the 2025 sea-trial float** of Delaigue et al. (2026), which is not on the
  GDAC: its spectra, CTD, GPS fixes and real AIS around its track are in
  `tests/data/float2025/` (see the README there; rebuild:
  `tests/make_fixtures_2025.py`).  It covers full 23-band spectra
  (3.15 kHz), band levels not normalised on board, three park depths
  (500, 740, 1020 m) and real AIS (decoded CSV and raw NMEA).

The β port is also checked against the published MATLAB depth correction
(`tests/data/matlab_ligurian_*`).
