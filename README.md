# argopal

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Wind from passive acoustics on BGC-Argo floats.

argopal starts from a float's WMO number, downloads its passive acoustic
(PAL) spectra from the Argo GDAC, and will process them into surface wind
speed.

## Install

    pip install -e .

## Import

```python
import argopal as ap
```

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

| format | what argopal needs |
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
| `wind_pensieri2015` | `U_pensieri2015` | linear < 38 dB, quadratic above (applied for 30 ≤ SPL < 60 dB) |

They need the 8 kHz band (processed by default);
`process_float(..., wind_models=[...])` picks a subset.  Each function can
also be called on its own, e.g. `ap.wind_nystuen2015(park["spl_8000Hz"])`.

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

    pytest             # offline tests
    pytest -m network  # tests that download from the GDAC

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
