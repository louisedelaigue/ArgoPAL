# Processing a float

```python
# what the float measured, per data mode
argopal.available_variables("argopal_data/7902280")

result = argopal.process_float(7902280)  # 8 kHz, mode="delayed" by default
result = argopal.process_float(7902280, bands=(3150, 8000))  # several bands

result["park"]         # one row per PAL measurement: bands, β, spectral
                       # levels, position, park T/S/O2/pH
result["profiles"]     # one row per profile level (Sprof)
result["ts_profiles"]  # the T/S profiles used for the depth correction
result["beta"]         # β of every band at each depth, per profile
result["positions"]    # launch + surface GPS fixes
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

## Bands

`bands=` selects the frequency bands to process (default 8 kHz; for wind
and rain, e.g. WOTAN, usually 8 and/or 3.15 kHz).  Only those bands are
despiked, depth-corrected and exported; `result["pal_raw"]` keeps every
band.  A warning is raised if a band is empty in most spectra: floats
sending 9-band spectra (63, 125, 400 Hz, 1, 2, 5, 8, 12.5, 20 kHz) have no
3.15 kHz values.

## Despiking

`despike_window="3h"` and `despike_percentile=95` (defaults, as in
Delaigue et al. 2026): within ±1.5 h of each measurement, values above the
95th percentile are discarded and the rest averaged.  Both can be changed,
or despiking turned off with `despike=False`.

## Data mode

`mode=` applies to every CTD/BGC variable, including the T/S used for the
depth correction.  Each mode takes values from one source only, with no
fallback:

| mode | values used |
|---|---|
| `"delayed"` (default) | delayed-mode (D) adjusted values; NaN otherwise |
| `"adjusted"` | real-time adjusted (A) values; NaN otherwise |
| `"raw"` | raw values (kept by Argo for every data mode) |

Every variable `X` comes with `X_qc` and `X_data_mode` columns.

## Export at any stage

`argopal.export` writes one table or any set of tables, from any step:

```python
argopal.export(result, "wmo7902280", format="csv")        # one CSV per table
argopal.export(result, "wmo7902280.nc", format="netcdf")  # one file, one group per table
argopal.export(argopal.clean_pal(pal), "pal_clean.csv")   # a single intermediate table
argopal.read_export("wmo7902280.nc")                      # back to DataFrames
```

CSV gives one file per table, without the processing record; NetCDF keeps
everything in one file, one group per table, with each table's processing
record as attributes.  See [Outputs](https://argopal.readthedocs.io/en/latest/outputs.html).

## Depth correction

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
