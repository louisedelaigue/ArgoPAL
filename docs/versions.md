# Version history

## Unreleased

- Reference winds: the measurements around each reference value are
  averaged before the models are fitted and assessed
  (`average_to_reference`; `reference_average=True` and
  `reference_window="1h"`, i.e. ±30 min, by default; any window width can
  be chosen).  Comparing every 5-min measurement with an hourly value
  gave vertical stripes in the fit figure and treated measurements that
  share one reference value as independent.
- `collocate_reference` adds `<name>_time`, the time of the reference
  value used.
- `select_bands(pal, bands="all")` keeps every band that holds data.
- Variable names from the core files are plain strings.
- `find_pal_floats`: every float on the GDAC with PAL data (found by
  reading the start of each `Rtraj_aux` file; saved and reused for a
  week), with where and when it measured; `read_profile_index`: the
  GDAC's synthetic-profile index.
- Tested on the three other PAL floats on the GDAC (2904023, 2904024,
  7902388), including full 23-band spectra.
- `plot_reference_map`: the float's trajectory and the references used to
  refit the wind models on a map, with the fitting radius and the
  measurements used (also in `plot_all` when there is a reference).
  `collocate_reference` records each reference's position (the most
  frequent one, and how far the others are).
- `reference_diagnostics` and `plot_reference_fits`: the wind model
  fitted to each reference separately and to all together, in-sample and
  cross-validated, and how each reference's coefficients do on the others.
- Wind footprints: with a reanalysis (`footprints=`), each measurement is
  matched to the reference whose wind footprint (correlation of the
  reanalysis wind with the wind at the reference, `wind_footprints`) has
  the highest R at its position, and the refit uses the measurements with
  R at least `footprint_min_r` (0.9) instead of those within 40 km;
  `plot_reference_map` outlines the footprints.
- Reference winds are converted to 10 m, the height the wind models are
  calibrated for: `wind_to_10m` (neutral logarithmic profile, or power
  law), from a `height_m` column or `reference_height_m=` (one value or
  one per reference).  Without a height, a `ReferenceHeightWarning` says
  the winds are assumed to be 10 m already.
- The example buoy winds (`examples/data/reference_winds_2026.csv`) give
  their anemometer height (3.5 m) and are converted to 10 m.
- `plot_reference_timeseries`: the reference wind and the float wind
  (published and refitted model) against time, for the points used in the
  refit, one row per reference, with ±1σ around the refitted wind (also
  in `plot_all` after a refit).
- `plot_reference_map` and `plot_reference_timeseries` show only the
  measurements used in the refit: within the radius and with a float wind
  (not, e.g., after the delayed-mode data end).  Circles are drawn at the
  positions of the reference values used (several for a buoy reported at
  several positions).  On the map, the trajectory is black where
  the float has a wind and grey where it has none (e.g. no delayed-mode
  profiles yet), and references passed only without a float wind are
  greyed out.
- `collocate_reference` adds `<name>_latitude` and `<name>_longitude`, the
  position of the reference value used.
- Maps use cartopy, with land and labelled gridlines (`map_axes`, used by
  `plot_positions` and the tutorials); cartopy is now a dependency.
- Tutorials: "Finding PAL floats and processing several" and
  "Understanding the depth correction".
- Tutorials: a second, shorter notebook, "Getting the PAL data" (every
  band in a park table with the CTD/BGC data, and a profiles table, with
  no further processing); the first is now "From PAL data to wind speed".
- Documentation: every function on one API page, grouped by processing
  step, and an Outputs page describing every table and column.

## 0.1.0b1 (2026-11-22)

First beta release.  The interface may still change before 1.0.

### Data from the Argo GDAC

- `download_aux`, `download_core`, `find_dac`, `list_files`: download a
  float's aux (PAL spectra, hydrophone metadata) and core (trajectory,
  synthetic profiles) files by WMO number; files are only downloaded again
  when the GDAC copy is newer.
- `read_pal`, `read_pal_meta`: PAL third-octave spectra and hydrophone
  configuration.  Band frequencies are assigned from the number of values
  (23 bands, or the 9-band subset); any other layout raises an error, and
  the band order is checked.
- `available_variables`: every variable in the core files, with its number
  of real-time, adjusted and delayed-mode values.
- Data mode: `"delayed"` (default), `"adjusted"` or `"raw"`, with no
  fallback from one to another.

### Processing

- `select_bands`: one or several bands (8 kHz by default); bands without
  data in a float's spectra are dropped with a warning.
- `clean_pal`: park-drift measurements only, QC flags, despiking (3 h
  window, 95th percentile by default, both configurable).
- `attach_park`, `attach_positions`: park T, S, O₂, pH and position
  interpolated in time to each acoustic measurement.
- Depth correction: β from each profile's sound-speed structure (ray
  method, Francois–Garrison absorption, TEOS-10), interpolated in time
  between the profiles around each measurement; `beta_method="mean"` for
  one β from the deployment-mean profile.  Absorption is applied as
  10^(−αL/10), with α in dB m⁻¹.
- `spectral_level`: bandwidth normalisation, unless the float already
  normalised on board.
- `process_float`: the whole chain from a WMO number.

### Wind speed

- Published 8 kHz models: Vagle et al. (1990), Nystuen et al. (2015),
  Pensieri et al. (2015, two regimes: linear 30–38 dB, quadratic
  38–60 dB).
- Reference winds (`read_reference`, `collocate_reference`): one table can
  hold several references (`reference_name`); each measurement is matched
  to the closest reference with a value within 1 h.  Optional
  `wind_speed_sigma` per value.
- `fit_wind_models`: refits each model to the references within a chosen
  radius (40 km by default).
- `wind_skill`, `wind_metrics`: bias, RMSE, R², Pearson r², slope,
  intercept; overall and per reference.
- `wind_uncertainty`: four columns per refitted model — leave-one-cycle-out
  empirical σ by wind speed, cycle-bootstrap coefficient σ, reference σ,
  and their combination.

### Ship noise (optional)

- `read_ais` (table, decoded CSV or raw NMEA with tag blocks; NMEA needs
  `pip install argopal[ais]`), `flag_ships`: measurements with a vessel
  within 20 km and ±30 min.
- `flag_ship_deviation`: with a reference wind as well, only measurements
  that also deviate from the reference by more than the RMSE are marked
  contaminated.

### Figures and export

- One figure per step (`plot_pal_spectra`, `plot_despike`, `plot_park`,
  `plot_positions`, `plot_profiles`, `plot_beta`, `plot_depth_correction`,
  `plot_wind`, `plot_wind_fit`, `plot_ais`), or all at once with
  `plot_all`; saved at 300 dpi.
- `export`, `read_export`: any table, or all of them, to CSV or NetCDF at
  any stage.

### Documentation and tests

- Tutorial notebook (`examples/argopal_tutorial.ipynb`) on float 7902280
  with three Météo-France buoys.
- Tests on two floats: 7902280 (GDAC) and the 2025 sea-trial float of
  Delaigue et al. (2026), with real AIS; the depth correction is checked
  against the published MATLAB code.
