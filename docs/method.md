# How it works

`process_float` runs the steps below in order; each is also a function on
its own.

| step | function(s) | result |
|---|---|---|
| 1. Download | `download_aux`, `download_core` | GDAC files |
| 2. Read and clean | `read_pal`, `clean_pal` | park-drift spectra, QC applied, despiked |
| 3. Park data | `read_park`, `attach_park` | CTD/BGC at each measurement |
| 4. Positions | `read_positions`, `attach_positions` | position at each measurement |
| 5. Profiles | `read_profiles`, `gridded_profiles` | T/S profiles for the depth correction |
| 6. Depth correction | `compute_profile_betas`, `depth_correct` | β at each measurement |
| 7. Spectral level | `spectral_level` | `spl_<f>Hz`, dB re 1 µPa² Hz⁻¹ |
| 8. Wind speed | `wind_literature` | `U_<model>`, published models |
| 9. Calibration | `collocate_reference`, `fit_wind_models`, `wind_skill`, `wind_uncertainty` | `U_<model>_fit` and its σ |
| optional | `read_ais`, `flag_ships` | ship noise flags |

## The acoustic data

Floats fitted with a PAL (Passive Acoustic Listener) hydrophone send
third-octave band levels (dB re 1 µPa), decoded into the GDAC's
**auxiliary** folder, `aux/<dac>/<wmo>/`.  The band frequencies are not
in the files: ArgoPAL assigns them from the number of values per
spectrum, and refuses any layout it does not know.

:::{dropdown} Files and band layouts
| file | contents |
|---|---|
| `<wmo>_Rtraj_aux.nc` | band levels per measurement (`PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM`), with time, pressure, cycle and measurement code |
| `<wmo>_meta_aux.nc` | hydrophone model and serial number, `CONFIG_Pal*` parameters |

| values | bands (Hz) | constant |
|---|---|---|
| 23 | 63, 100, 125, 160, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000, 25000 | `argopal.PAL_BANDS_HZ` |
| 9 | 63, 125, 400, 1000, 2000, 5000, 8000, 12500, 20000 | `argopal.PAL_BANDS_9_HZ` |

Each read checks that the median spectrum falls with frequency, as
ambient noise does, and raises a `BandOrderWarning` if not.

`CONFIG_PalBandNormalization_LOGICAL` says whether the float already
divided each band level by its bandwidth; `spectral_level` subtracts the
bandwidth term only if it did not.
:::

No GDAC index records which floats have a PAL: `find_pal_floats` looks
for PAL spectra in every float's aux files (a couple of minutes; the list
is saved and reused for a week).

## Cleaning

Only park-drift measurements are kept (measurement code 290), and
spectra with bad QC are dropped.  Despiking removes transient noise:
within ±1.5 h of each measurement, values above the 95th percentile are
discarded and the rest averaged (Delaigue et al. 2026).

## Park data and positions

The park CTD/BGC samples are interpolated linearly in time to each
acoustic measurement, within a cycle and across gaps of at most 2 h.
Positions are interpolated between surface GPS fixes; `position_gap_h`
gives the time between those fixes.

## Depth correction

Wind-generated sound is produced at the surface and attenuated on its
way down to the hydrophone.  β (dB) is that attenuation, computed from
the sound-speed profile with the ray method of Vagle et al. (1990), as
implemented by Cauchy et al. (2018), and Francois–Garrison absorption.

By default β is computed for **every profile** and interpolated in time
between the profiles before and after each measurement, so it follows
the seasonal change of the sound-speed profile (for 7902280, at 1000 m and
8 kHz: 0.94 dB in February, 1.42 dB in August).  More than 6 days from
any usable profile, β is NaN.  `beta_method="mean"` uses one
deployment-mean profile, as in Delaigue et al. (2026).

:::{dropdown} Details: profile ends and absorption units
T/S profiles are held constant above their shallowest point and below
their deepest point.  `gridded_profiles(..., surface="linear")`
reproduces the paper's linear extrapolation to the surface, which can
overshoot in a strong gradient.

The Francois–Garrison absorption α is in dB m⁻¹, so the loss along a ray
is 10^(−αL/10).  The original MATLAB code used exp(−αL), which
overestimates β about 4.3 times (1000 m, 8 kHz: 4.0 dB instead of
0.96 dB).  The tests check both: with the old convention the port
reproduces the published table.
:::

The depth-corrected spectral level, `spl_<f>Hz` (dB re 1 µPa² Hz⁻¹), is
the band level plus β, minus the bandwidth term if the float did not
normalise.

## Wind models

Three published models give wind speed from the 8 kHz spectral level
(Delaigue et al. 2026, Eqs. 2–4):

| model | column | form |
|---|---|---|
| Vagle et al. (1990) | `U_vagle1990` | log-quadratic |
| Nystuen et al. (2015) | `U_nystuen2015` | cubic |
| Pensieri et al. (2015) | `U_pensieri2015` | linear 30–38 dB, quadratic 38–60 dB |

Each is also a function, e.g. `argopal.wind_nystuen2015(spl)`.  These
models were not made for floats and are biased on them: calibrating
against reference winds corrects this.

:::{note}
Acoustic wind has two limits that calibration does not remove.  Below
about 2–3 m s⁻¹ few waves break, other noise sets the level, and the
acoustic wind levels off: low values mean "below about 3 m s⁻¹".  Above
about 15 m s⁻¹ bubble clouds absorb part of the sound, and floats seldom
see such winds near a reference, so the refitted wind tends to be low.
:::

## Calibration against reference winds

ArgoPAL does not download reference winds: a table of buoy (or ship)
winds is given.

**Matching.**  Each measurement is matched to the closest reference with
a value within 1 h (`reference_max_time`).

**10 m winds.**  The models are calibrated for 10 m winds, so reference
winds measured lower are converted to 10 m first, with a neutral
logarithmic profile (`wind_to_10m`), from `height_m` or
`reference_height_m`.

**Averaging.**  The float measures more often than most references
(e.g. every 5 min against an hourly buoy).  The measurements from 30 min
before to 30 min after each reference value are averaged
(`reference_window="1h"`), and the models are fitted on these points.

**One calibration per reference.**  The models are refitted to each
reference separately (nonlinear least squares, from the published
coefficients), on the float data within 40 km of it
(`fit_max_distance_km`) or inside its footprint.  Each measurement then
uses the coefficients of the reference closest to it.  A reference needs
float data in at least 4 cycles (`calibration_min_cycles`; 3 is the
lowest that still allows cross-validation) to be calibrated; the
measurements closest to a reference without coefficients get **no
refitted wind**, e.g. where the float has no delayed-mode data near it
yet.  Such a reference is not counted in `<name>_used` and is greyed out
on the map.  `calibration="pooled"` fits one set of coefficients to all
references together instead.

New columns tell what was used: `<name>_used` (the measurements used for
a refit), `<name>_closest` (the closest reference) and
`<name>_calibration` (whose coefficients apply).
`argopal.plot_wind_timeline` shows all three through the deployment.

:::{dropdown} Reference table columns
Columns are recognised by name:

| content | names |
|---|---|
| time (UTC) | `datetime`, `time`, `date`, … |
| wind speed (m s⁻¹) | `wind_speed`, `wind`, `ff`, `u10`, … |
| position | `latitude`/`lat`, `longitude`/`lon` (or `latitude=`/`longitude=` for one fixed buoy) |
| reference name (several references) | `reference_name`, `station`, `buoy_name`, `platform` |
| height of the wind (m) | `height_m`, `sensor_height`, `anemometer_height` |
| error of each value (m s⁻¹) | `wind_speed_sigma`, `wind_sigma`, `sigma` |

Positions are read per value, so moving references (ships) work.
:::

## Uncertainty

`wind_uncertainty` gives the uncertainty (1σ) of the refitted wind of each
measurement.  By default it is the **cross-validated error**,
`sigma_cv`: the model is refitted on all cycles but one and predicts the
cycle left out, and the RMSE of these residuals, binned by wind speed,
gives σ as a function of wind.  It is out of sample, so it includes every
source of error, including the uncertainty of coefficients fitted on few
cycles.

Two terms are optional (`uncertainty=`):

| term | what it is |
|---|---|
| `"coef"` | spread of the refitted wind over refits on cycles drawn at random with replacement (bootstrap).  Large where the refit is extrapolated, and unreliable with fewer than about 10 cycles |
| `"ref"` | the reference's own error, from `reference_sigma` or the table's `wind_speed_sigma`.  `sigma_cv` contains it, so it is removed: √(σ_cv² − σ_ref²) |

With all three, σ = √(max(σ_cv² − σ_ref², 0) + σ_coef²).

## Do the references agree?

`reference_diagnostics` fits the model to each reference separately and
to all together, and gives how each reference's coefficients do on the
others, in-sample and cross-validated.  A reference whose curve or
coefficients differ from the others is worth checking (anemometer
height, distance, exposure).  `plot_reference_fits` shows each curve.

## Wind footprints

Near the coast, the wind 40 km from a buoy can differ from the wind at
the buoy, and a buoy further away can be a better match.  With a
reanalysis (`footprints=`), each reference's **footprint** is where the
reanalysis wind correlates with the wind at the reference with R of at
least 0.9 (`wind_footprints`; the reference's own data are not used).
Each measurement is then matched to the reference with the highest R at
its position, and calibration uses the measurements inside the
footprint, however far.  `plot_reference_map` shades each footprint.

By default the wind is correlated as given (e.g. hourly over the
deployment).  `footprint_anomaly="monthly"` correlates monthly anomalies
of a reanalysis over several years: a climatological footprint.

## Ship noise

Ship noise can dominate the 1–10 kHz band up to about 20 km away.  With
AIS vessel positions, each measurement with a vessel within 20 km and
±30 min is flagged (`ships_nearby`, `ship_flag`; Delaigue et al. 2026).
With a reference too, a measurement is `ship_contaminated` only if its
wind also deviates from the reference by more than the RMSE.

:::{dropdown} AIS formats
| format | needs |
|---|---|
| table | `datetime` (UTC), `mmsi`, `latitude`, `longitude` |
| decoded CSV (`.csv`, `.csv.gz`) | a time column (`datetime`, `BaseDateTime`, `# Timestamp`, `tagblock_timestamp`, …), `mmsi`, latitude and longitude |
| raw NMEA (`!AIVDM`, plain or `.gz`) | tag blocks with the receiver time; needs the `ais` extra (pyais) |

Only vessel position reports count (message types 1, 2, 3, 18, 19, 27).
Times must be when each position was received, not taken from file
names: daily AIS files often span two days.
:::
