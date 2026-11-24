"""Build examples/argopal_tutorial.ipynb (run, then execute the notebook).

python examples/make_tutorial.py
"""

import os

import nbformat as nbf


HERE = os.path.dirname(os.path.abspath(__file__))
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip()))


md("""
# From PAL data to wind speed

**ArgoPAL** turns the passive acoustic (PAL) measurements of BGC-Argo floats
into surface wind speed, starting from a float's WMO number and using only
data from the Argo GDAC.

This notebook goes through the full pipeline, from the acoustic spectra to
depth-corrected and calibrated wind speed, on float **7902280** (Ligurian
Sea, 2026).  To get only the acoustic data (every band, with the float's
CTD/BGC data, no further processing), see
[Getting the PAL data](pal_data_tutorial.ipynb).

1. one-call processing
2. what the float measured
3. processing options
4. the pipeline step by step
5. figures to follow each step
6. wind speed from published models
7. calibration against reference buoys, with uncertainty
8. ship noise with AIS (optional)
9. export

Installation:

```
pip install git+https://github.com/louisedelaigue/ArgoPAL.git
pip install "argopal[ais] @ git+https://github.com/louisedelaigue/ArgoPAL.git"  # raw NMEA AIS
```
""")

code("""
import warnings

import argopal

warnings.simplefilter("ignore")  # ArgoPAL warnings are discussed in the text
argopal.hello()
""")

md("""
## 1. One call

`process_float` downloads the float's files from the GDAC (acoustic spectra
from the *aux* folder; CTD, oxygen, pH and positions from the core files),
then runs every processing step.  Defaults: delayed-mode data only, 8 kHz
band, β per profile, despiking over ±1.5 h at the 95th percentile.
""")

code("""
result = argopal.process_float(7902280, verbose=False)
list(result)
""")

md("""
`result["park"]` holds one row per acoustic measurement: band level, depth
correction β, depth-corrected spectral level, position, park T/S/O₂ and wind
speed from each published model.
""")

code("""
park = result["park"]
cols = ["cycle", "datetime", "latitude", "longitude", "pressure",
        "f_8000Hz", "beta_8000Hz", "spl_8000Hz", "U_nystuen2015"]
park[cols].dropna().head()
""")

md("""
In delayed mode, every value comes from delayed-mode Argo data only, with no
fallback.  For a young float, most profiles are not in delayed mode yet, so
the depth correction (and everything after it) is available only for the
first months:
""")

code("""
print(f"{park['spl_8000Hz'].notna().mean():.0%} of measurements have a "
      f"depth-corrected level in {park.attrs['park']['mode']!r} mode")
cycles = result["ts_profiles"].attrs["cycles"]
print(f"profiles used for β: cycles {cycles[0]} to {cycles[-1]}")
""")

md("""
## 2. What the float measured

The variables available in each file, with how many values are in delayed
(D), adjusted (A) and real-time (R) mode:
""")

code("""
argopal.available_variables("argopal_data/7902280")
""")

md("""
The hydrophone and its configuration come from the `meta_aux` file.  This
float normalises band levels by bandwidth on board, so they are already
spectral levels:
""")

code("""
meta = result["meta"]
print(meta["sensor_model"], meta["sensor_serial_no"])
print("band normalisation on board:", meta["band_normalization"])
print("spectra by band layout:", result["pal_raw"].attrs["band_layout_counts"])
""")

md("""
Almost all spectra are 9-band spectra (63, 125, 400 Hz, 1, 2, 5, 8, 12.5,
20 kHz), stored in the 23-band array with the other bands empty:
""")

code("""
fig, axes = argopal.plot_pal_spectra(result["pal_raw"])
""")

md("""
## 3. Options

- `mode`: `"delayed"` (default), `"adjusted"` or `"raw"`; one source each,
  no fallback.
- `bands`: frequency bands to process; default 8 kHz.  Bands without data
  are dropped with a warning (here 3.15 kHz, absent from 9-band spectra).
- `despike_window`, `despike_percentile`; `beta_method="mean"` for one β
  from the deployment-mean profile.
""")

code("""
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    raw = argopal.process_float(7902280, mode="raw", bands=(3150, 8000), verbose=False)
print([str(w.message)[:110] for w in caught if "band" in str(w.message)])
print("bands processed:", raw["park"].attrs["bands_hz"])
print(f"depth-corrected in raw mode: {raw['park']['spl_8000Hz'].notna().mean():.0%}")
""")

md("""
## 4. Step by step

Each step of `process_float` is a function of its own; tables can be
inspected, changed or exported between steps.
""")

code("""
files = argopal.download_aux(7902280, verbose=False) | argopal.download_core(7902280, verbose=False)

meta = argopal.read_pal_meta(files["7902280_meta_aux.nc"])
pal = argopal.read_pal(files["7902280_Rtraj_aux.nc"], meta=meta)           # spectra
pal = argopal.select_bands(pal, bands=8000)                                # 8 kHz
pal = argopal.clean_pal(pal)                                               # park drift, QC, despike

park_samples = argopal.read_park(files["7902280_Rtraj.nc"], mode="raw")    # park CTD/BGC
pal = argopal.attach_park(pal, park_samples)                               # linear in time
pal = argopal.attach_positions(pal, argopal.read_positions(files["7902280_Rtraj.nc"]))

profiles = argopal.read_profiles(files["7902280_Sprof.nc"], mode="raw")
ts = argopal.gridded_profiles(profiles, max_pressure=pal["pressure"].max() + 10)
beta = argopal.compute_profile_betas(ts, pal.attrs["bands_hz"])           # β per profile
pal = argopal.depth_correct(pal, beta)                                     # β in time and depth
pal = argopal.spectral_level(pal)                                          # dB re 1 µPa² Hz⁻¹
pal = argopal.wind_literature(pal)                                         # published models

pal[["datetime", "pressure", "f_8000Hz", "beta_8000Hz", "spl_8000Hz",
     "U_vagle1990", "U_nystuen2015", "U_pensieri2015"]].dropna().head()
""")

md("""
## 5. Figures

One figure per step, in the style of Delaigue et al. (2026).
`argopal.plot_all(result, save_dir="figures")` saves all of them at 300 dpi;
here a few, on the raw-mode run.
""")

code("""
fig, axes = argopal.plot_despike(raw["pal_raw"], raw["park"], cycle=10)
""")

code("""
fig, axes = argopal.plot_profiles(raw["ts_profiles"])
""")

md("""
β follows the season: the summer thermocline changes the sound-speed
profile and so the propagation of surface noise down to the hydrophone.
""")

code("""
fig, axes = argopal.plot_beta(raw["beta"])
""")

code("""
fig, axes = argopal.plot_depth_correction(raw["park"], beta=raw["beta"])
""")

md("""
## 6. Wind speed from published models

Vagle et al. (1990), Nystuen et al. (2015) and Pensieri et al. (2015), all at
8 kHz, with their published coefficients (`U_<model>` columns).  Each is also
a function: e.g. `argopal.wind_nystuen2015(spl)`.
""")

code("""
fig, axes = argopal.plot_wind(raw["park"])
""")

md("""
## 7. Calibration against reference buoys

Published coefficients are not tuned to floats.  With a reference wind, the
models are refitted (nonlinear least squares, measurements within 40 km of
a reference), and each refitted wind gets an uncertainty.

ArgoPAL does not download reference winds: one table holds them all, here
three Météo-France buoys (DYFAMED, LION, PACA), January–June 2026.  Columns
`reference_name`, `datetime`, `latitude`, `longitude`, `wind_speed`, and
optionally `wind_speed_sigma`, the error of each value, and `height_m`, the
height of the wind measurement.

The wind models are calibrated against **10 m wind**, so the reference
must be 10 m wind too.  These buoys measure wind at about 3.5 m
(`height_m` = 3.5); ArgoPAL converts it to 10 m when reading the table
(`argopal.wind_to_10m`, a factor of about 1.10 here).  The float wind,
from models calibrated against 10 m wind, is then 10 m wind as well.
""")

code("""
import pandas as pd

buoys = pd.read_csv("data/reference_winds_2026.csv")
buoys.groupby("reference_name")[["latitude", "longitude", "height_m"]].mean().round(2)
""")

md("""
Each measurement is matched to the closest buoy with a value within 1 h.

The float measures every 5 min and the buoys report hourly.  So for each
buoy value, the measurements from **30 min before to 30 min after** its
time are averaged first (`reference_window="1h"`, the default; any width
can be given, e.g. `"30min"` or `"3h"`).  The models are fitted and
assessed on these hourly points, then applied to every measurement.

`fit_max_distance_km` sets the radius used for fitting, skill and
uncertainty (40 km by default).  The buoys give no error per value, so one
error for all of them is passed as `reference_sigma`.
""")

code("""
cal = argopal.process_float(
    7902280,
    mode="raw",
    reference={"buoys": "data/reference_winds_2026.csv"},
    fit_max_distance_km=40,
    reference_sigma=0.5,   # the buoys' wind error (m s⁻¹)
    verbose=False,
)
cal["park"].attrs["references"]["buoys"]["references"]
""")

md("""
The float passed within 40 km of DYFAMED in February 2026 (cycles 9–12)
and of PACA in April 2026 (cycles 17–22); while the LION record lasts
(to June), the float is more than 100 km from it.  The refit pools the
measurements near both buoys.
""")

code("""
p = cal["park"]
p.groupby("buoys_reference_name")["buoys_distance_km"].describe()[
    ["count", "min", "50%", "max"]
].round(1)
""")

md("""
On a map: the trajectory, the measurements within 40 km of the buoy they
were matched to (black, those used for the refit), and the buoys used,
each with a 40 km circle (LION, never used, is not shown).  With
`annotate=True`, the figure also notes buoys reported at several
positions: PACA appears at two positions about 25 km apart in the
Météo-France file, and is shown at the more frequent one.
""")

code("""
fig, ax = argopal.plot_reference_map(cal["park"], annotate=True)
""")

code("""
cal["wind_fits"][["model", "reference_name", "n", "n_applied",
                  "rmse_published", "rmse_fitted", "r2_fitted",
                  "bias_published"]].round(2)
""")

md("""
**One calibration per buoy.**  Each buoy's coefficients are fitted on the
float data near it, and used there (`n` points).  Every other measurement
uses the coefficients of the closest buoy (`n_applied` in all), and
`buoys_calibration` says whose.  `calibration="pooled"` would fit one set
of coefficients to both buoys together.

The whole record, two months per row: the refitted wind with ±1σ, the
periods used for each buoy's refit shaded in its colour, and along the
top whose coefficients apply:
""")

code("""
fig, axes = argopal.plot_wind_timeline(cal["park"])
""")

md("""
Skill of every wind estimate against the buoys within 40 km, all together
and per buoy:
""")

code("""
cal["wind_skill"].round(2)
""")

md("""
Uncertainty, per measurement and per refitted model:

- `_sigma_cv`: empirical error, cross-validated by leaving one cycle out,
  binned by wind speed;
- `_sigma_coef`: coefficient uncertainty, bootstrap over cycles;
- `_sigma_ref`: the buoys' own error (`reference_sigma`, or the table's
  `wind_speed_sigma`);
- `_sigma`: combined, √(σ_cv² − σ_ref² + σ_coef²).

The empirical error by wind-speed bin:
""")

code("""
cal["wind_uncertainty"].round(2)
""")

md("""
Ten cycles pass within 40 km of a buoy, the minimum ArgoPAL asks for the
bootstrap; with fewer, it warns that the coefficient term is unreliable.
The refit should still not be trusted far outside the winds seen near the
buoys.
""")

code("""
cols = [c for c in cal["park"] if c.startswith("U_nystuen2015_fit_sigma")]
cal["park"][cols].quantile([0.5, 0.9]).round(2)
""")

md("""
Refitted wind against the buoys, coloured by distance (as in Fig. 5 of
Delaigue et al. 2026).  By default the figure shows every measurement within
80 km, so its statistics include measurements further away than those used
for fitting:
""")

code("""
fig, axes = argopal.plot_wind_fit(cal["park"])
""")

md("""
Restricted to the 40 km used for fitting, the statistics match `wind_fits`:
""")

code("""
fig, axes = argopal.plot_wind_fit(cal["park"], max_distance_km=40)
""")

md("""
The same points through time: the buoy wind, and Nystuen et al. (2015)
with its published and refitted coefficients.  Each point belongs to the
buoy the float was matched to, the closest one with a value within the
hour, so one row per buoy:
""")

code("""
fig, axes = argopal.plot_reference_timeseries(cal["park"], annotate=True)
""")

code("""
fig, axes = argopal.plot_wind(cal["park"])
""")

md("""
### More with references

The same steps can be run one by one on any processed table, here the
raw-mode `park` table from section 3, which has no reference yet.

**Radius.** `fit_max_distance_km` in `process_float`, or `max_distance_km`
in `fit_wind_models`, `wind_skill` and `wind_uncertainty`, sets how close a
measurement must be to its reference to be used; `None` uses all.  To
compare settings with one number each, these examples fit one set of
coefficients for both buoys (`calibration="pooled"`).
""")

code("""
park = argopal.collocate_reference(raw["park"], buoys, name="buoys")
for radius in (10, 25, 40, 80):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, fits = argopal.fit_wind_models(
            park, max_distance_km=radius, calibration="pooled"
        )
    row = fits.set_index("model").loc["nystuen2015"]
    print(f"{radius:>3} km: n = {row['n']:>6}, RMSE = {row['rmse_fitted']:.2f}, "
          f"R² = {row['r2_fitted']:.2f}")
""")

md("""
No measurement is within 10 km of a buoy, so nothing is refitted there
(the model keeps its row, with n = 0).  The RMSE grows with the radius, as
buoy winds further from the float represent its wind less well; R² also
depends on the range of winds covered, which is narrower within 25 km.
""")

md("""
**Which reference each measurement used.**  The closest buoy with a value
within `max_time` (1 h by default), with its distance and time difference:
""")

code("""
park[["datetime", "U_buoys", "buoys_reference_name", "buoys_distance_km",
      "buoys_dt_min"]].dropna().iloc[::5000]
""")

md("""
**Averaging window.**  `window` in `fit_wind_models` (or
`reference_window` in `process_float`) sets how much float data is
averaged around each buoy value: `"1h"` (default) is 30 min before to
30 min after.  `"nearest"` averages the measurements whose nearest buoy
value it is, and `average=False` uses every measurement on its own.
""")

code("""
for window in ("20min", "1h", "2h"):
    _, fits = argopal.fit_wind_models(
        park, max_distance_km=40, window=window, calibration="pooled"
    )
    row = fits.set_index("model").loc["nystuen2015"]
    print(f"{window:>6}: {row['n']} buoy values, "
          f"{row['n_measurements'] / row['n']:.1f} measurements each, "
          f"RMSE = {row['rmse_fitted']:.2f}")
_, fits = argopal.fit_wind_models(
    park, max_distance_km=40, average=False, calibration="pooled"
)
row = fits.set_index("model").loc["nystuen2015"]
print(f"  none: {row['n']} measurements, RMSE = {row['rmse_fitted']:.2f}")
""")

md("""
**Collocation time limit and distance.**  `max_time` sets how far in
time the nearest buoy value may be;
`max_distance_km` in `collocate_reference` blanks reference values further
away than that (by default they are kept, so that fitting and figures can
choose).
""")

code("""
strict = argopal.collocate_reference(
    raw["park"], buoys, name="buoys", max_time="30min", max_distance_km=40
)
print(f"{strict['U_buoys'].notna().sum()} measurements with a buoy wind "
      f"within 30 min and 40 km, of {len(strict)}")
""")

md("""
**Skill per reference.**  `by_reference=True` adds one set of rows per
buoy to the skill table:
""")

code("""
fitted, _ = argopal.fit_wind_models(park, max_distance_km=40)
skill = argopal.wind_skill(fitted, max_distance_km=40, by_reference=True)
skill[skill["wind"] == "U_nystuen2015_fit"].round(2)
""")

md("""
**Error of each reference value.**  A `wind_speed_sigma` column in the
table is carried to each measurement (`<name>_sigma`) and used as σ_ref by
`wind_uncertainty`, unless `reference_sigma` is given.  The values below
are only for illustration: these buoys do not give an error per value.
""")

code("""
illustrative = buoys.assign(
    wind_speed_sigma=buoys["reference_name"].map(
        {"DYFAMED": 0.5, "LION": 0.5, "PACA": 0.8}
    )
)
p = argopal.collocate_reference(raw["park"], illustrative, name="buoys")
p, _ = argopal.fit_wind_models(p, models=["nystuen2015"])
p, _ = argopal.wind_uncertainty(p, n_bootstrap=20)
by_buoy = p.attrs["wind_uncertainty"]["models"]["nystuen2015"]
print("σ_ref from the table:",
      {buoy: round(v["sigma_ref"], 2) for buoy, v in by_buoy.items()})
""")

md("""
**Do the buoys agree?**  `reference_diagnostics` fits the model to each
buoy separately and to all together.  A buoy's own fit always looks
better on its own data; the cross-validated RMSE (each cycle predicted
by a fit without it) says whether it really is, and the transfer table
how one buoy's coefficients do on the other:
""")

code("""
d = argopal.reference_diagnostics(cal["park"])
d["fits"][["reference", "n", "n_cycles", "rmse_published", "rmse_pooled",
           "rmse_own", "rmse_pooled_cv", "rmse_own_cv"]].round(2)
""")

code("""
d["transfer"].round(2)
""")

code("""
fig, ax = argopal.plot_reference_fits(cal["park"])
""")

md("""
Here DYFAMED's curve is close to the published one and PACA's lower: the
two buoys disagree about the float-wind relationship, and the pooled
refit lies between them.

**Wind footprints instead of a distance.**  Near the coast, the wind a
few tens of kilometres from a buoy can differ from the buoy's.  With a
reanalysis of the 10 m wind over the float's area and period (e.g. ERA5
`u10`/`v10`), each measurement can be matched to the buoy whose wind
footprint (where the reanalysis wind varies as at the buoy, R ≥ 0.9)
covers it, instead of the closest buoy within 40 km.  No reanalysis comes
with these examples; with one:

```python
result = argopal.process_float(
    7902280,
    reference="data/reference_winds_2026.csv",
    footprints="era5_u10_v10_2026.nc",
    footprint_min_r=0.9,
)
argopal.plot_reference_map(result["park"], footprints=result["footprints"])
```
""")

md("""
**Other table layouts.**

- Column names are recognised from common aliases (`time`, `ff`, `u10`,
  `lat`, `lon`, `station`, `buoy_name`, `sigma`, ...); see
  `argopal.read_reference`.
- A single fixed buoy whose table has no position columns:
  `argopal.collocate_reference(park, "buoy.csv", latitude=43.36, longitude=7.84)`.
- Moving references (ships, drifters): positions are read for every value.
- A dict collocates several sources separately, each under its own name:
  `reference={"dyfamed": "dyfamed.csv", "paca": "paca.csv"}`, with
  `fit_reference="paca"` to choose the one used for refitting.
""")

md("""
## 8. Ship noise with AIS (optional)

Ship noise can dominate the 8 kHz band up to ~20 km away.  Given AIS vessel
positions for the float's area and period, measurements with a vessel within
20 km and ±30 min are flagged (`ships_nearby`, `ship_flag`).  With a
reference buoy as well, measurements whose refitted wind also deviates from
the buoy by more than the RMSE are marked `ship_contaminated`
(`ais_action="remove"` blanks only those).

No AIS data cover 7902280's 2026 track here, so this cell is not run.  AIS
can be a table (`datetime`, `mmsi`, `latitude`, `longitude`), decoded CSV
files, or raw NMEA files with tag-block receiver times (needs
`pip install argopal[ais]`):

```python
res = argopal.process_float(
    7902280,
    ais="AIS/messages-*.gz",          # or a DataFrame / CSV files
    reference="reference_winds.csv",
    ais_action="remove",
)
fig, axes = argopal.plot_ais(res["park"])
```
""")

md("""
## 9. Export

Any table, or any set of tables, at any stage, as CSV or NetCDF:
""")

code("""
argopal.export(cal, "output/wmo7902280", format="csv")
""")

code("""
files = argopal.export(cal, "output/wmo7902280.nc", format="netcdf")
back = argopal.read_export(files[0])
list(back), back["park"].shape
""")

code("""
argopal.export(pal, "output/pal_step_by_step.csv")   # a single intermediate table
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    },
    "language_info": {"name": "python"},
}
nbf.write(nb, os.path.join(HERE, "argopal_tutorial.ipynb"))
print("written", len(cells), "cells")
