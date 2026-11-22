"""The argopal processing pipeline, from a WMO number."""

from .ais import flag_ships, read_ais, track_bbox
from .argo import DEFAULT_MODE, read_park, read_profiles
from .calibration import (
    fit_wind_models,
    flag_ship_deviation,
    wind_skill,
    wind_uncertainty,
)
from .clean import DEFAULT_BANDS, clean_pal, select_bands
from .combine import attach_park
from .constants import GDAC_URL
from .correction import (
    BETA_METHODS,
    compute_beta,
    compute_profile_betas,
    depth_correct,
    gridded_profiles,
    mean_profile,
    spectral_level,
)
from .gdac import download_aux, download_core, find_dac
from .positions import attach_positions, read_positions
from .read import read_pal, read_pal_meta
from .reference import collocate_reference
from .wind import wind_literature


def process_float(
    wmo,
    data_dir="argopal_data",
    mode=DEFAULT_MODE,
    variables=None,
    bands=DEFAULT_BANDS,
    despike=True,
    despike_window="3h",
    despike_percentile=95,
    max_gap="2h",
    bad_qc=("3", "4", "9"),
    band_normalization=None,
    beta_method="interpolate",
    wind_models=None,
    ais=None,
    ais_radius_km=20.0,
    ais_window="30min",
    ais_action="flag",
    reference=None,
    reference_max_time="1h",
    fit_reference=None,
    fit_max_distance_km=40.0,
    reference_sigma=None,
    n_bootstrap=100,
    dac=None,
    update=True,
    gdac=GDAC_URL,
    verbose=True,
):
    """Download a float and run every processing step.

    Steps (each is also a standalone function):
      1. download the aux (PAL) and core (meta, Rtraj, Sprof) files;
      2. `read_pal` + `clean_pal`: park-drift spectra, bad QC dropped,
         despiked;
      3. `read_park` + `attach_park`: park CTD/BGC interpolated in time
         onto each PAL measurement;
      4. `read_positions` + `attach_positions`: float position at each PAL
         measurement, between surface GPS fixes;
         (optional) `read_ais` + `flag_ships`: ships within 20 km and
         ±30 min of each PAL measurement, if `ais` is given;
      5. `read_profiles` + `gridded_profiles`: every usable T/S profile
         in `mode`, on a regular grid;
      6. `compute_profile_betas` + `depth_correct`: depth correction β of
         every band for each profile (Francois-Garrison absorption in
         dB/m), interpolated in time between the profiles before and
         after each PAL measurement (or, with beta_method="mean",
         `mean_profile` + `compute_beta`: one β from the mean profile);
      7. `spectral_level`: depth-corrected spectral levels
         spl_<f>Hz (dB re 1 µPa²/Hz);
      8. `wind_literature`: wind speed U_<model> from the published
         models at 8 kHz (Vagle 1990, Nystuen 2015, Pensieri 2015);
      9. (if `reference` is given) `collocate_reference`, then
         `fit_wind_models` (refitted models U_<model>_fit, within 40 km of
         the reference), `wind_skill` and `wind_uncertainty` (σ per
         measurement: cross-validated, coefficient, reference, combined); with `ais` too,
         `flag_ship_deviation` (ship nearby and wind error > RMSE).

    Any of the returned tables can be written with `export`, as can the
    output of any single step.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    data_dir : str, optional
        Local folder for downloads, by default "argopal_data".
    mode : str, optional
        Argo data mode for all CTD/BGC variables, including the T/S used
        for the depth correction: "delayed" (default), "adjusted" or
        "raw".  Each mode uses one source only, with no fallback; see
        `argopal.argo.select_mode`.
    variables : list of str, optional
        CTD/BGC parameters to keep (e.g. ["TEMP", "PSAL", "DOXY"]); by
        default all that the float measured.  TEMP and PSAL are always
        read for the depth correction.
    bands : float or sequence of float, optional
        Frequency bands (Hz) to process, by default (8000,).  Several are
        allowed, e.g. (3150, 8000).  `pal_raw` keeps every band.
    despike : bool, optional
        Despike the PAL spectra, by default True.
    despike_window : str, optional
        Full width of the despiking window, by default "3h" (±1.5 h).
    despike_percentile : float, optional
        In each window, values above this percentile are discarded before
        averaging, by default 95.
    max_gap : str, optional
        Largest gap between park samples to interpolate across, by
        default "2h".
    bad_qc : tuple of str, optional
        QC flags of park samples not to interpolate from, by default
        ("3", "4", "9").  Raw (unadjusted) DOXY is flagged "3" by Argo
        convention, so with mode="raw" park oxygen stays empty unless "3"
        is removed from this list.
    band_normalization : bool, optional
        Whether the float normalised band levels by bandwidth; by default
        read from the meta_aux file.
    beta_method : {"interpolate", "mean"}, optional
        "interpolate" (default): β from each profile, interpolated linearly
        in time onto each PAL measurement, so it follows seasonal changes
        in the sound-speed profile.  "mean": one β from the deployment-mean
        profile, as in Delaigue et al. (2026).
    wind_models : list of str, optional
        Literature wind models to apply (see `argopal.wind`); by default
        all of them (they need the 8 kHz band).
    ais : DataFrame, str or list of str, optional
        AIS vessel positions for the ship-noise check (see `read_ais` for
        the accepted formats: a table, decoded CSV files, or raw NMEA files
        with tag-block times).  By default no check.
    ais_radius_km, ais_window : optional
        Distance and half time window for a ship to count, by default
        20 km and "30min" (Delaigue et al. 2026).
    ais_action : {"flag", "remove"}, optional
        "flag" (default) adds ships_nearby / ship_flag columns.  "remove"
        blanks the acoustic values of measurements with a ship nearby; with
        a `reference` too, only of those whose refitted Nystuen wind also
        deviates from the reference by more than the RMSE (both criteria,
        as in Delaigue et al. 2026).
    reference : DataFrame, str or dict, optional
        Reference wind to collocate (see `collocate_reference`): one table
        or CSV with reference_name, datetime, latitude, longitude,
        wind_speed and, optionally, wind_speed_sigma.  It can hold several
        references (buoys): each measurement is matched to the closest one
        with a value within `reference_max_time`.  A single source is named
        "ref"; a dict, name -> source, names it (e.g. {"buoys":
        "reference_winds.csv"}) or collocates several sources separately.
        argopal does not download reference winds.
    reference_max_time : str, optional
        Largest time difference for collocation, by default "1h".
    fit_reference : str, optional
        Reference used to refit the models; by default the first one.
    fit_max_distance_km : float, optional
        Measurements used for fitting, skill and uncertainty must be this
        close to the reference they are matched to, by default 40 km;
        None uses all.
    reference_sigma : float, optional
        Error of the reference wind (m s-1), removed from the empirical
        uncertainty; by default the table's wind_speed_sigma, else 0.  See
        `wind_uncertainty`.
    n_bootstrap : int, optional
        Bootstrap refits for the coefficient uncertainty, by default 100.
    dac, update, gdac, verbose
        As for `download_aux`.

    Returns
    -------
    dict
        "park" : one row per PAL measurement: band levels, β, spectral
            levels, wind speeds, position and park CTD/BGC;
        "pal_raw" : the PAL spectra as read (all bands), before cleaning;
        "park_samples" : the park-drift CTD/BGC samples, before
            interpolation onto the PAL measurements;
        "profiles" : one row per profile level;
        "ts_profiles" : the T/S profiles used for β, on a regular grid;
        "mean_profile" : (beta_method="mean" only) the mean T/S profile;
        "beta" : β of every band at each depth (per profile, with cycle
            and datetime columns, for beta_method="interpolate");
        "positions" : the surface fixes;
        "ais" : (if `ais` was given) the AIS positions used;
        "wind_fits" : (if `reference`) published and refitted
            coefficients of each model, with their skill;
        "wind_skill" : (if `reference`) bias, RMSE, R² of every wind
            estimate against the fitting reference;
        "wind_uncertainty" : (if `reference`) the cross-validated,
            wind-dependent σ of each refitted model, by wind-speed bin;
        "meta" : PAL sensor and configuration (`read_pal_meta`);
        "files" : file name -> local path of everything downloaded.
    """
    if beta_method not in BETA_METHODS:
        raise ValueError(f"beta_method must be one of {BETA_METHODS}.")
    if dac is None:
        dac = find_dac(wmo, gdac=gdac)
    kw = {
        "data_dir": data_dir,
        "dac": dac,
        "update": update,
        "gdac": gdac,
        "verbose": verbose,
    }
    files = {**download_aux(wmo, **kw), **download_core(wmo, **kw)}
    wmo = int(wmo)

    meta_path = files.get(f"{wmo}_meta_aux.nc")
    meta = read_pal_meta(meta_path) if meta_path else None
    pal_raw = read_pal(files[f"{wmo}_Rtraj_aux.nc"], meta=meta)
    pal = clean_pal(
        select_bands(pal_raw, bands),
        despike_spectra=despike,
        despike_window=despike_window,
        despike_percentile=despike_percentile,
    )

    park = read_park(files[f"{wmo}_Rtraj.nc"], mode=mode, variables=variables)
    pal = attach_park(pal, park, max_gap=max_gap, bad_qc=bad_qc)

    positions = read_positions(files[f"{wmo}_Rtraj.nc"])
    pal = attach_positions(pal, positions)

    ais_positions = None
    if ais is not None:
        ais_positions = read_ais(
            ais, bbox=track_bbox(pal, ais_radius_km), verbose=verbose
        )
        pal = flag_ships(
            pal,
            ais_positions,
            radius_km=ais_radius_km,
            window=ais_window,
            # with a reference, removal waits for the deviation criterion
            action="flag" if reference is not None else ais_action,
        )

    if variables is not None:
        variables = list(dict.fromkeys(["TEMP", "PSAL", *variables]))
    profiles = read_profiles(
        files[f"{wmo}_Sprof.nc"], mode=mode, variables=variables
    )
    # 10 dbar margin so β covers the deepest PAL measurement
    max_pressure = pal["pressure"].max() + 10
    bands = pal.attrs["bands_hz"]
    ts_profiles = gridded_profiles(profiles, max_pressure)
    extra = {}
    if beta_method == "interpolate":
        beta = compute_profile_betas(ts_profiles, bands)
    else:
        extra["mean_profile"] = mean_profile(profiles, max_pressure)
        beta = compute_beta(extra["mean_profile"], bands)
    pal = depth_correct(pal, beta)
    pal = spectral_level(pal, band_normalization=band_normalization)
    pal = wind_literature(pal, models=wind_models)

    if reference is not None:
        refs = reference if isinstance(reference, dict) else {"ref": reference}
        for name, source in refs.items():
            opts = {}
            if isinstance(source, dict):
                opts = {k: v for k, v in source.items() if k != "source"}
                source = source["source"]
            pal = collocate_reference(
                pal, source, name=name, max_time=reference_max_time, **opts
            )
        fit_ref = fit_reference or next(iter(refs))
        pal, extra["wind_fits"] = fit_wind_models(
            pal, reference=fit_ref, max_distance_km=fit_max_distance_km
        )
        extra["wind_skill"] = wind_skill(
            pal,
            reference=fit_ref,
            max_distance_km=fit_max_distance_km,
            by_reference=True,
        )
        pal, extra["wind_uncertainty"] = wind_uncertainty(
            pal, reference_sigma=reference_sigma, n_bootstrap=n_bootstrap
        )
        if ais is not None and "U_nystuen2015_fit" in pal:
            pal = flag_ship_deviation(
                pal,
                reference=fit_ref,
                max_distance_km=fit_max_distance_km,
                action=ais_action,
            )

    return {
        "park": pal,
        "pal_raw": pal_raw,
        "park_samples": park,
        "profiles": profiles,
        "ts_profiles": ts_profiles,
        **extra,
        "beta": beta,
        "positions": positions,
        **({"ais": ais_positions} if ais_positions is not None else {}),
        "meta": meta,
        "files": files,
    }
