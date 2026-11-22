"""Depth correction and spectral levels of PAL measurements."""

import warnings

import gsw
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .acoustics import beta_profile
from .clean import band_columns


GOOD_QC = ("1", "2", "5", "8")
SURFACE_METHODS = ("constant", "linear")
BETA_METHODS = ("interpolate", "mean")


# --- T/S profiles on a regular grid -----------------------------------------


def _good_levels(g, good_qc):
    """Pressure-sorted good T/S of one profile (duplicates averaged)."""
    ok = (
        np.isfinite(g[["PRES", "TEMP", "PSAL"]].to_numpy(float)).all(1)
        & g["PRES_qc"].isin(good_qc)
        & g["TEMP_qc"].isin(good_qc)
        & g["PSAL_qc"].isin(good_qc)
    )
    g = g[ok].groupby("PRES", as_index=False)[["TEMP", "PSAL"]].mean()
    return (
        g["PRES"].to_numpy(float),
        g["TEMP"].to_numpy(float),
        g["PSAL"].to_numpy(float),
    )


def _profile_on_grid(p, t, s, grid, surface):
    """T/S of one profile on `grid` (dbar).

    PCHIP inside the measured range.  Above the shallowest point: held
    constant (surface="constant", a mixed layer) or extrapolated linearly
    from the two shallowest points (surface="linear", as in Delaigue et al.
    2026; this can overshoot when those points are close together in a
    strong gradient).  Below the deepest point (down to the end of
    `grid`): held constant.
    """
    ti = PchipInterpolator(p, t, extrapolate=False)(grid)
    si = PchipInterpolator(p, s, extrapolate=False)(grid)
    above, below = grid < p[0], grid > p[-1]
    if surface == "linear":
        ti[above] = t[0] + (t[1] - t[0]) / (p[1] - p[0]) * (grid[above] - p[0])
        si[above] = s[0] + (s[1] - s[0]) / (p[1] - p[0]) * (grid[above] - p[0])
    else:
        ti[above], si[above] = t[0], s[0]
    ti[below], si[below] = t[-1], s[-1]
    return ti, si


def gridded_profiles(
    profiles,
    max_pressure,
    max_surface_pressure=10.0,
    max_bottom_gap=100.0,
    min_bottom_pressure=100.0,
    surface="constant",
    good_qc=GOOD_QC,
    dp=1.0,
):
    """Put every usable T/S profile on a regular pressure grid.

    A profile is used if its good T/S (QC in `good_qc`) start within
    `max_surface_pressure` of the surface and reach at least
    `min_bottom_pressure`.  Each profile is gridded from the surface down to
    its own deepest point plus `max_bottom_gap` (at most `max_pressure`):
    profiles from a shallow park depth are kept, and only serve PAL
    measurements as deep as they reach (see `depth_correct`).  Above its
    shallowest point T/S are held constant (or extrapolated linearly,
    `surface="linear"`); below its deepest point, within `max_bottom_gap`,
    they are held constant.

    Parameters
    ----------
    profiles : pd.DataFrame
        Output of `read_profiles`; its data mode is the one used.
    max_pressure : float
        Deepest pressure needed (dbar), e.g. the deepest PAL measurement.
    max_surface_pressure : float, optional
        By default 10 dbar.
    max_bottom_gap : float, optional
        How far below its deepest point a profile is extended (T/S held
        constant), by default 100 dbar.
    min_bottom_pressure : float, optional
        Shallower profiles are skipped, by default 100 dbar.
    surface : {"constant", "linear"}, optional
        T/S above the shallowest point; by default "constant".
    good_qc : tuple of str, optional
        QC flags accepted for PRES, TEMP and PSAL.
    dp : float, optional
        Grid step (dbar), by default 1.

    Returns
    -------
    pd.DataFrame
        Long table: cycle, datetime, latitude, longitude, pressure,
        temperature, salinity, measured (False where held/extrapolated).
        Profiles end at different pressures.
        `attrs["skipped"]` gives the reason for each unused profile.

    Raises
    ------
    ValueError
        If no profile is usable.
    """
    if surface not in SURFACE_METHODS:
        raise ValueError(f"surface must be one of {SURFACE_METHODS}.")
    grid = np.arange(0, max_pressure + dp / 2, dp)
    tables, skipped = [], {}
    for cycle, g in profiles.groupby("cycle"):
        p, t, s = _good_levels(g, good_qc)
        if len(p) < 3:
            skipped[int(cycle)] = "fewer than 3 good T/S levels"
            continue
        if p[0] > max_surface_pressure:
            skipped[int(cycle)] = f"shallowest good level {p[0]:.0f} dbar"
            continue
        if p[-1] < min_bottom_pressure:
            skipped[int(cycle)] = f"deepest good level {p[-1]:.0f} dbar"
            continue
        own = grid[grid <= p[-1] + max_bottom_gap]
        ti, si = _profile_on_grid(p, t, s, own, surface)
        tables.append(
            pd.DataFrame(
                {
                    "cycle": int(cycle),
                    "datetime": g["datetime"].iloc[0],
                    "latitude": g["latitude"].iloc[0],
                    "longitude": g["longitude"].iloc[0],
                    "pressure": own,
                    "temperature": ti,
                    "salinity": si,
                    "measured": (own >= p[0]) & (own <= p[-1]),
                }
            )
        )
    if not tables:
        raise ValueError(
            f"No profile in data mode {profiles.attrs.get('mode')!r} has "
            f"good T/S from <= {max_surface_pressure:.0f} dbar down to "
            f">= {min_bottom_pressure:.0f} dbar, so no depth correction can "
            "be computed."
        )
    out = pd.concat(tables, ignore_index=True)
    out.attrs = {
        "wmo": profiles.attrs.get("wmo"),
        "mode": profiles.attrs.get("mode"),
        "cycles": sorted(out["cycle"].unique().tolist()),
        "skipped": skipped,
        "surface": surface,
    }
    return out


def mean_profile(
    profiles,
    max_pressure,
    max_surface_pressure=10.0,
    max_bottom_gap=100.0,
    min_bottom_pressure=100.0,
    surface="constant",
    good_qc=GOOD_QC,
    dp=1.0,
):
    """Mean T/S profile of a float (all usable profiles averaged).

    Profiles are first put on a regular grid with `gridded_profiles` (same
    parameters), then averaged level by level over the profiles reaching
    each level.

    Returns
    -------
    pd.DataFrame
        Columns pressure, temperature, salinity, temperature_std,
        salinity_std, n_measured (profiles with a measurement at each
        level, as opposed to held or extrapolated values).  `attrs` holds
        the cycles used, their mean position and the data mode.
    """
    grid = gridded_profiles(
        profiles,
        max_pressure,
        max_surface_pressure=max_surface_pressure,
        max_bottom_gap=max_bottom_gap,
        min_bottom_pressure=min_bottom_pressure,
        surface=surface,
        good_qc=good_qc,
        dp=dp,
    )
    by = grid.groupby("pressure")
    out = pd.DataFrame(
        {
            "temperature": by["temperature"].mean(),
            "salinity": by["salinity"].mean(),
            "temperature_std": by["temperature"].std(ddof=0),
            "salinity_std": by["salinity"].std(ddof=0),
            "n_measured": by["measured"].sum(),
        }
    ).reset_index()
    first = grid.groupby("cycle").first()
    out.attrs = {
        **grid.attrs,
        "n_profiles": len(grid.attrs["cycles"]),
        "latitude": float(first["latitude"].mean()),
        "longitude": float(first["longitude"].mean()),
    }
    return out


# --- β from T/S -------------------------------------------------------------


def _beta_from_ts(
    pressure, temperature, salinity, lat, lon, frequencies_hz, dz, step
):
    """β of each frequency from one T/S profile on a pressure grid."""
    z = -gsw.z_from_p(pressure, lat)
    depth = np.arange(dz, np.ceil(z.max()) + dz / 2, dz)
    t = np.interp(depth, z, temperature)
    s = np.interp(depth, z, salinity)
    p_grid = gsw.p_from_z(-depth, lat)
    sa = gsw.SA_from_SP(s, p_grid, lon, lat)
    ct = gsw.CT_from_t(sa, t, p_grid)
    c = gsw.sound_speed(sa, ct, p_grid)
    receivers = np.arange(step, depth[-1] + dz / 2, step)
    receivers = receivers[receivers <= depth[-1]]
    if receivers[-1] < depth[-1]:  # always reach the bottom of the grid
        receivers = np.append(receivers, depth[-1])
    return beta_profile(depth, c, t, s, frequencies_hz, receivers)


def compute_beta(mean_prof, frequencies_hz, dz=1.0, receiver_step=10.0):
    """Depth correction β for each band, from one (mean) T/S profile.

    Pressure is converted to depth (TEOS-10, at the profile's mean
    latitude); sound speed uses TEOS-10 (Absolute Salinity and Conservative
    Temperature from the measured practical salinity and in-situ
    temperature); absorption uses Francois & Garrison (1982) in dB/m.

    Parameters
    ----------
    mean_prof : pd.DataFrame
        Output of `mean_profile`.
    frequencies_hz : list of float
        Band centre frequencies, e.g. `pal.attrs["bands_hz"]`.
    dz : float, optional
        Layer thickness for the ray computation (m), by default 1.
    receiver_step : float, optional
        Spacing of the depths at which β is computed (m), by default 10.

    Returns
    -------
    pd.DataFrame
        From `beta_profile`: "depth", "beta_<f>Hz", "alpha_<f>Hz".
    """
    lat = mean_prof.attrs["latitude"]
    lon = mean_prof.attrs["longitude"]
    beta = _beta_from_ts(
        mean_prof["pressure"].to_numpy(float),
        mean_prof["temperature"].to_numpy(float),
        mean_prof["salinity"].to_numpy(float),
        lat,
        lon,
        frequencies_hz,
        dz,
        receiver_step,
    )
    beta.attrs.update(
        {
            "method": "mean",
            "wmo": mean_prof.attrs.get("wmo"),
            "mode": mean_prof.attrs.get("mode"),
            "profile_cycles": mean_prof.attrs.get("cycles"),
            "latitude": lat,
            "longitude": lon,
        }
    )
    return beta


def compute_profile_betas(
    ts_profiles, frequencies_hz, dz=1.0, receiver_step=10.0
):
    """Depth correction β for each band, from every profile separately.

    Parameters
    ----------
    ts_profiles : pd.DataFrame
        Output of `gridded_profiles`.
    frequencies_hz, dz, receiver_step
        As for `compute_beta`.

    Returns
    -------
    pd.DataFrame
        Long table: cycle, datetime (of the profile), latitude, longitude,
        depth, "beta_<f>Hz" and "alpha_<f>Hz".
    """
    tables = []
    for cycle, g in ts_profiles.groupby("cycle"):
        lat, lon = g["latitude"].iloc[0], g["longitude"].iloc[0]
        b = _beta_from_ts(
            g["pressure"].to_numpy(float),
            g["temperature"].to_numpy(float),
            g["salinity"].to_numpy(float),
            lat,
            lon,
            frequencies_hz,
            dz,
            receiver_step,
        )
        b.insert(0, "longitude", lon)
        b.insert(0, "latitude", lat)
        b.insert(0, "datetime", g["datetime"].iloc[0])
        b.insert(0, "cycle", int(cycle))
        tables.append(b)
    out = pd.concat(tables, ignore_index=True)
    out.attrs = {
        "method": "interpolate",
        "wmo": ts_profiles.attrs.get("wmo"),
        "mode": ts_profiles.attrs.get("mode"),
        "profile_cycles": ts_profiles.attrs.get("cycles"),
        "absorption_units": tables[0].attrs.get("absorption_units"),
        "latitude": float(out["latitude"].mean()),
    }
    return out


# --- Applying β to PAL measurements ------------------------------------------


def _beta_at(depth_table, values, d):
    """β of one profile at depths d: linear in depth, extrapolated linearly
    to the surface, NaN below the deepest β depth."""
    val = np.interp(d, depth_table, values, right=np.nan)
    above = d < depth_table[0]
    slope = (values[1] - values[0]) / (depth_table[1] - depth_table[0])
    val[above] = values[0] + slope * (d[above] - depth_table[0])
    return val


def _covering_neighbours(t, d, tp, zmax):
    """For each measurement (time t, depth d): index of the latest profile
    at or before t and of the earliest at or after t, among the profiles
    whose β reaches depth d (-1 if none)."""
    covers = d[:, None] <= zmax[None, :] + 1e-6
    before = covers & (tp[None, :] <= t[:, None])
    after = covers & (tp[None, :] >= t[:, None])
    k = len(tp)
    left = np.where(before.any(1), k - 1 - np.argmax(before[:, ::-1], 1), -1)
    right = np.where(after.any(1), np.argmax(after, 1), -1)
    return left, right


def depth_correct(pal, beta, max_extrapolation="6D"):
    """Add the depth correction β of each band to a PAL table.

    Each measurement's pressure is converted to depth (TEOS-10).  β is
    interpolated linearly in depth.

    With per-profile β (from `compute_profile_betas`), β is then
    interpolated linearly in time between the nearest profile before and
    the nearest profile after each measurement, among the profiles that
    reach its depth (on floats changing park depth, shallow profiles only
    serve shallow parks).  With a profile on one side only, its β is used if
    it is within `max_extrapolation`; otherwise β is NaN (the sound-speed
    profile may have changed with the season).  With a single β (from
    `compute_beta`), it applies to all measurements.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table (from `read_pal`, `clean_pal`, `attach_park`, ...).
    beta : pd.DataFrame
        Output of `compute_profile_betas` or `compute_beta`.
    max_extrapolation : str or pd.Timedelta, optional
        Per-profile β only: how far from the only profile available on one
        side its β may be used, by default "6D" (about one cycle).

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "depth" and "beta_<f>Hz" columns (dB); with
        per-profile β also "beta_cycle_before", "beta_cycle_after" (the
        profiles used) and "beta_weight_after" (time weight of the later
        one).
    """
    lat = pal.get("latitude")
    if lat is None or lat.isna().all():
        lat = beta.attrs["latitude"]
    else:
        lat = lat.fillna(beta.attrs["latitude"]).to_numpy(float)
    out = pal.copy()
    out["depth"] = -gsw.z_from_p(out["pressure"].to_numpy(float), lat)
    d = out["depth"].to_numpy(float)
    bands = pal.attrs["bands_hz"]

    if "cycle" in beta:
        groups = sorted(
            beta.groupby("cycle"), key=lambda x: x[1]["datetime"].iloc[0]
        )
        cycles = np.array([c for c, _ in groups])
        tp = np.array(
            [g["datetime"].iloc[0] for _, g in groups], dtype="datetime64[ns]"
        ).astype("i8")
        zmax = np.array([g["depth"].max() for _, g in groups])
        t = out["datetime"].to_numpy(dtype="datetime64[ns]").astype("i8")
        left, right = _covering_neighbours(t, d, tp, zmax)
        limit = pd.Timedelta(max_extrapolation).value
        only_left = (left >= 0) & (right < 0)
        only_right = (right >= 0) & (left < 0)
        left = np.where(only_left & (t - tp[left] > limit), -1, left)
        right = np.where(only_right & (tp[right] - t > limit), -1, right)
        # one side only: use it on both sides
        left = np.where((left < 0) & (right >= 0) & only_right, right, left)
        right = np.where((right < 0) & (left >= 0) & only_left, left, right)
        ok = (left >= 0) & (right >= 0)
        li, ri = np.where(ok, left, 0), np.where(ok, right, 0)
        span = (tp[ri] - tp[li]).astype(float)
        w = np.where(span > 0, (t - tp[li]) / np.where(span > 0, span, 1), 0.0)
        w = np.where(ok, w, np.nan)
        for f in bands:
            col = f"beta_{int(f)}Hz"
            b_left = np.full(len(d), np.nan)
            b_right = np.full(len(d), np.nan)
            for k, (_, g) in enumerate(groups):
                zb = g["depth"].to_numpy(float)
                bv = g[col].to_numpy(float)
                for side, idx in ((b_left, li), (b_right, ri)):
                    rows = ok & (idx == k)
                    if rows.any():
                        side[rows] = _beta_at(zb, bv, d[rows])
            out[col] = (1 - w) * b_left + w * b_right
        out["beta_cycle_before"] = pd.array(
            np.where(ok, cycles[li], 0), dtype="Int64"
        )
        out["beta_cycle_after"] = pd.array(
            np.where(ok, cycles[ri], 0), dtype="Int64"
        )
        out.loc[~ok, ["beta_cycle_before", "beta_cycle_after"]] = pd.NA
        out["beta_weight_after"] = w
        n_none = int((~ok).sum())
        if n_none:
            warnings.warn(
                f"{n_none} PAL measurements have no usable profile (data "
                f"mode {beta.attrs.get('mode')!r}) reaching their depth "
                f"before and after them, nor one within {max_extrapolation}"
                "; their β is NaN.",
                stacklevel=2,
            )
    else:
        zb = beta["depth"].to_numpy(float)
        for f in bands:
            col = f"beta_{int(f)}Hz"
            out[col] = _beta_at(zb, beta[col].to_numpy(float), d)
        n_deep = int((d > zb[-1]).sum())
        if n_deep:
            warnings.warn(
                f"{n_deep} PAL measurements are deeper than the deepest β "
                f"({zb[-1]:.0f} m); their β is NaN.",
                stacklevel=2,
            )
    out.attrs = dict(pal.attrs)
    out.attrs["depth_correction"] = {
        "method": beta.attrs.get("method"),
        "max_extrapolation": str(max_extrapolation),
        "absorption_units": beta.attrs.get("absorption_units"),
        "profile_cycles": beta.attrs.get("profile_cycles"),
        "mode": beta.attrs.get("mode"),
    }
    return out


def third_octave_bandwidth(frequency_hz):
    """Bandwidth (Hz) of a base-10 one-third-octave band (IEC 61260).

    The nominal frequency (e.g. 8000) is mapped to its exact centre
    1000 * 10**(k/10) (7943.3 Hz); the band spans exact centre times
    10**(±1/20).
    """
    f = np.asarray(frequency_hz, dtype=float)
    k = np.round(10 * np.log10(f / 1000))
    exact = 1000 * 10 ** (k / 10)
    return exact * (10 ** (1 / 20) - 10 ** (-1 / 20))


def spectral_level(pal, band_normalization=None):
    """Depth-corrected spectral levels (dB re 1 µPa²/Hz) of every band.

    spl_<f>Hz = f_<f>Hz + beta_<f>Hz [- 10 log10(bandwidth)]

    The bandwidth term is subtracted only if the float did not already
    normalise its band levels by bandwidth before sending them.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `depth_correct`.
    band_normalization : bool, optional
        Whether the float already normalised band levels by bandwidth
        (CONFIG_PalBandNormalization_LOGICAL).  By default taken from
        `pal.attrs["band_normalization"]` (set when the meta_aux file is
        read); an error is raised if it is unknown.

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "spl_<f>Hz" columns.
    """
    if band_normalization is None:
        band_normalization = pal.attrs.get("band_normalization")
    if band_normalization is None:
        raise ValueError(
            "Unknown whether the float normalised band levels by bandwidth "
            "(CONFIG_PalBandNormalization_LOGICAL). Read the meta_aux file "
            "with read_pal(..., meta=read_pal_meta(...)) or pass "
            "band_normalization=True/False."
        )
    out = pal.copy()
    for f, col in zip(pal.attrs["bands_hz"], band_columns(pal)):
        beta_col = f"beta_{int(f)}Hz"
        if beta_col not in pal:
            raise KeyError(f"{beta_col} missing: run depth_correct first.")
        level = pal[col] + pal[beta_col]
        if not band_normalization:
            level = level - 10 * np.log10(third_octave_bandwidth(f))
        out[f"spl_{int(f)}Hz"] = level
    out.attrs = dict(pal.attrs)
    out.attrs["spl"] = {
        "units": "dB re 1 uPa^2/Hz",
        "band_normalization_on_float": bool(band_normalization),
        "bandwidth_subtracted": not band_normalization,
    }
    return out
