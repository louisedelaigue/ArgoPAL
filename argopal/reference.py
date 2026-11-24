"""Reference winds: collocate buoy winds with PAL measurements.

ArgoPAL does not download reference winds.  The reference is one table
(DataFrame or CSV file) of measured wind speed at known times and
positions, from one or several references (buoys, ships, stations).
Columns, recognised by name (case-insensitive):

  time:            datetime, time, date, time_utc, validity_time
  wind speed:      wind_speed, wind, windspeed, u, u10, ff, ws  (m s-1)
  position:        latitude / lat, longitude / lon / long; for a single
                   fixed station they can instead be given as numbers.
  name (optional): reference_name, station, buoy_name, platform
                   Needed when the table holds several references.
  error (optional): wind_speed_sigma, wind_sigma, sigma  (m s-1)
                   Error of each reference value, used by
                   `wind_uncertainty`.
  height (optional): height_m, sensor_height, anemometer_height  (m)
                   Height of the wind_speed values above the sea surface.

The wind models are calibrated against wind at 10 m, so the reference must
be 10 m wind.  With a height (a column, or `height_m=`), the wind is
converted to 10 m (`wind_to_10m`); without one, it is assumed to be 10 m
already, with a `ReferenceHeightWarning`.

Each PAL measurement is matched to the closest reference (in distance)
among those with a value within `max_time` of it; for each reference, its
value nearest in time is taken.  The reference should be the wind at 10 m,
like the models.
"""

import warnings

import numpy as np
import pandas as pd


class ReferenceHeightWarning(UserWarning):
    """The height of a reference wind is unknown: assumed to be 10 m."""


TABLE_ALIASES = {
    "datetime": ("datetime", "time", "date", "time_utc", "validity_time"),
    "wind_speed": ("wind_speed", "wind", "windspeed", "u", "u10", "ff", "ws"),
    "latitude": ("latitude", "lat"),
    "longitude": ("longitude", "lon", "long"),
    "reference_name": ("reference_name", "station", "buoy_name", "platform"),
    "wind_speed_sigma": ("wind_speed_sigma", "wind_sigma", "sigma"),
    "height_m": ("height_m", "sensor_height", "anemometer_height"),
}
OPTIONAL = ("reference_name", "wind_speed_sigma", "height_m")
Z0_OPEN_SEA = 1.52e-4  # m, roughness length of the open sea (neutral)
POWER_LAW_EXPONENT = 0.11  # Hsu et al. (1994), over the sea


def wind_to_10m(wind_speed, height_m, method="log", z0=Z0_OPEN_SEA):
    """Wind speed at 10 m from wind measured at another height.

    For a neutral atmosphere over the sea:

    - "log" (default): logarithmic profile,
      U10 = U(z) ln(10 / z0) / ln(z / z0), with roughness length `z0`
      (1.52e-4 m by default, typical of the open sea);
    - "power": power law, U10 = U(z) (10 / z) ** 0.11 (Hsu et al. 1994).

    Parameters
    ----------
    wind_speed : array_like
        Wind speed (m s-1) at `height_m`.
    height_m : float or array_like
        Height of the measurement above the sea surface (m).
    method : {"log", "power"}, optional
        Profile, by default "log".
    z0 : float, optional
        Roughness length (m) for "log", by default 1.52e-4.

    Returns
    -------
    np.ndarray
        Wind speed at 10 m (m s-1).
    """
    u = np.asarray(wind_speed, dtype=float)
    z = np.asarray(height_m, dtype=float)
    if np.any(z <= 0):
        raise ValueError("height_m must be positive.")
    if method == "log":
        return u * np.log(10.0 / z0) / np.log(z / z0)
    if method == "power":
        return u * (10.0 / z) ** POWER_LAW_EXPONENT
    raise ValueError(f'method must be "log" or "power"; got {method!r}.')


def _haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (
        np.sin((p2 - p1) / 2) ** 2
        + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    )
    return 2 * r * np.arcsin(np.sqrt(a))


def read_reference(
    source, latitude=None, longitude=None, height_m=None, height_method="log"
):
    """Read a reference wind table (see module docstring), with the wind at
    10 m.

    Parameters
    ----------
    source : pd.DataFrame or str
        A table or a CSV file.
    latitude, longitude : float, optional
        Position of a fixed station, if the table has no position columns.
    height_m : float or dict, optional
        Height (m) of the wind_speed values, if the table has no height
        column: one value, or reference_name -> height.  10 for winds
        already at 10 m.
    height_method : {"log", "power"}, optional
        Profile used to convert to 10 m (see `wind_to_10m`), by default
        "log".

    Returns
    -------
    pd.DataFrame
        Columns datetime (UTC), wind_speed (at 10 m), latitude, longitude,
        and reference_name / wind_speed_sigma / height_m when known; sorted
        by time, without missing winds.  `attrs["height"]` records the
        height of each reference and how the wind was converted.

    Warns
    -----
    ReferenceHeightWarning
        If the height of some values is unknown: they are assumed to be
        10 m wind already.
    """
    df = source if isinstance(source, pd.DataFrame) else pd.read_csv(source)
    lower = {c.strip().lower(): c for c in df.columns}
    found = {}
    for target, aliases in TABLE_ALIASES.items():
        for a in aliases:
            if a in lower:
                found[target] = lower[a]
                break
    if latitude is not None:
        found.pop("latitude", None)
    if longitude is not None:
        found.pop("longitude", None)
    needed = ["datetime", "wind_speed"]
    if latitude is None:
        needed.append("latitude")
    if longitude is None:
        needed.append("longitude")
    missing = [t for t in needed if t not in found]
    if missing:
        raise ValueError(
            f"Reference wind table: no column for {', '.join(missing)}. "
            "Accepted names: "
            + "; ".join(f"{t}: {', '.join(TABLE_ALIASES[t])}" for t in missing)
            + ". For a fixed station, pass latitude= and longitude=."
        )
    out = df[list(found.values())].rename(
        columns={v: k for k, v in found.items()}
    )
    if latitude is not None:
        out["latitude"] = float(latitude)
    if longitude is not None:
        out["longitude"] = float(longitude)
    t = pd.to_datetime(out["datetime"], errors="coerce", utc=True)
    out["datetime"] = t.dt.tz_localize(None).astype("datetime64[ns]")
    for c in ("wind_speed", "latitude", "longitude"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if "wind_speed_sigma" in out:
        out["wind_speed_sigma"] = pd.to_numeric(
            out["wind_speed_sigma"], errors="coerce"
        )
    if "reference_name" in out:
        out["reference_name"] = out["reference_name"].astype(str).str.strip()
    out = _to_10m(out, height_m, height_method)
    out = out.dropna(subset=["datetime", "wind_speed"])
    cols = ["datetime", "wind_speed", "latitude", "longitude"]
    attrs = out.attrs
    out = out[cols + [c for c in OPTIONAL if c in out]]
    out = out.sort_values("datetime").reset_index(drop=True)
    out.attrs = attrs
    return out


def _to_10m(out, height_m, method):
    """Convert wind_speed to 10 m where the height is known; warn where it
    is not."""
    if height_m is not None:
        if isinstance(height_m, dict):
            if "reference_name" not in out:
                raise ValueError(
                    "height_m as a dict needs a reference_name column."
                )
            h = out["reference_name"].map(
                {str(k): float(v) for k, v in height_m.items()}
            )
        else:
            h = pd.Series(float(height_m), index=out.index)
        if "height_m" in out:  # a column fills what is not given
            h = h.fillna(pd.to_numeric(out["height_m"], errors="coerce"))
        out["height_m"] = h
    elif "height_m" in out:
        out["height_m"] = pd.to_numeric(out["height_m"], errors="coerce")
    known = out["height_m"].notna() if "height_m" in out else None
    if known is not None and known.any():
        out.loc[known, "wind_speed"] = wind_to_10m(
            out.loc[known, "wind_speed"], out.loc[known, "height_m"], method
        )
    names = (
        out["reference_name"] if "reference_name" in out
        else pd.Series("reference", index=out.index)
    )  # fmt: skip
    unknown = sorted(
        names[~known].unique() if known is not None else names.unique()
    )
    if unknown:
        warnings.warn(
            "Reference wind height unknown for "
            + ", ".join(str(u) for u in unknown)
            + ": assumed to be 10 m.  The wind models are calibrated "
            "against 10 m wind; give the height of the wind_speed values "
            "(a height_m column, or height_m=) to convert them.",
            ReferenceHeightWarning,
            stacklevel=3,
        )
    heights = {}
    for name in names.unique():
        h = out.loc[names == name, "height_m"] if known is not None else None
        h = h.dropna().unique() if h is not None else []
        heights[str(name)] = (
            [float(x) for x in h] if len(h) else "unknown, assumed 10 m"
        )
    out.attrs["height"] = {"method": method, "height_m": heights}
    return out


def _position_summary(g):
    pos = g[["latitude", "longitude"]].round(3).dropna()
    if pos.empty:
        return {"latitude": None, "longitude": None, "n_positions": 0,
                "spread_km": None}  # fmt: skip
    counts = pos.value_counts()
    lat, lon = counts.index[0]
    spread = _haversine_km(lat, lon, pos["latitude"], pos["longitude"])
    return {
        "latitude": float(lat),
        "longitude": float(lon),
        "n_positions": len(counts),
        "spread_km": round(float(np.max(spread)), 1),
    }


def collocate_reference(
    pal,
    reference,
    name="ref",
    max_time="1h",
    max_distance_km=None,
    latitude=None,
    longitude=None,
    height_m=None,
    height_method="log",
    footprints=None,
):
    """Add a reference wind (at 10 m) to each PAL measurement.

    For each reference in the table (see the module docstring), its value
    nearest in time is taken, if within `max_time`; of these, the one from
    the reference closest to the float is kept (the nearer in time if two
    are equally close).  The distance, time difference and, for several
    references, the name of the reference used are recorded.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table with datetime, and latitude / longitude from
        `attach_positions`.
    reference : pd.DataFrame or str
        A wind table or CSV file (see `read_reference`), with one or
        several references.
    name : str, optional
        Short name used in the new columns, e.g. "dyfamed" or "buoys"; by
        default "ref".
    max_time : str or pd.Timedelta, optional
        Largest time difference allowed, by default "1h".
    max_distance_km : float, optional
        If given, reference values further than this from the float are
        set to NaN.  By default no limit: the distance is kept so that
        fitting and validation can choose (`fit_wind_models` uses 40 km).
    latitude, longitude : float, optional
        Position of a fixed station whose table has no position columns.
    height_m : float or dict, optional
        Height (m) of the reference winds if the table has no height
        column (one value, or reference_name -> height); they are
        converted to 10 m.  See `read_reference`.
    height_method : {"log", "power"}, optional
        Profile used for the conversion (`wind_to_10m`), by default "log".
    footprints : xr.DataArray, optional
        Wind footprints of the references (`wind_footprints`).  Each
        measurement is then matched to the reference whose footprint has
        the highest R at the float's position (the closest if equal),
        rather than to the closest reference, and that R is recorded
        ("<name>_footprint_r") for selecting by footprint when fitting.

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "U_<name>" (m s-1, at 10 m), "<name>_distance_km",
        "<name>_dt_min" (reference time minus measurement time, minutes),
        "<name>_time", "<name>_latitude", "<name>_longitude" (time and
        position of the reference value used),
        "<name>_reference_name" if the table names its references, and
        "<name>_sigma" (m s-1) if it gives their error.  `attrs
        ["references"][<name>]["positions"]` records where each reference
        is: its most frequent position, how many distinct positions it
        has and how far the furthest is ("spread_km").
    """
    col = f"U_{name}"
    models = pal.attrs.get("wind_models", {}).get("applied", {})
    if name in models or col in pal:
        raise ValueError(
            f"Column {col!r} already exists; choose another name."
        )
    if "latitude" not in pal:
        raise KeyError("No float positions: run attach_positions first.")
    tol_min = pd.Timedelta(max_time).total_seconds() / 60
    out = pal.copy()
    t = out["datetime"].to_numpy(dtype="datetime64[ns]")
    lat = out["latitude"].to_numpy(float)
    lon = out["longitude"].to_numpy(float)
    n = len(out)

    ref = read_reference(
        reference,
        latitude=latitude,
        longitude=longitude,
        height_m=height_m,
        height_method=height_method,
    )
    named = "reference_name" in ref
    groups = (
        list(ref.groupby("reference_name", sort=True)) if named
        else [(None, ref)]
    )  # fmt: skip
    has_sigma = "wind_speed_sigma" in ref

    best_dist = np.full(n, np.inf)
    best_dt = np.full(n, np.nan)
    best_t = np.full(n, np.datetime64("NaT", "ns"))
    best_lat = np.full(n, np.nan)
    best_lon = np.full(n, np.nan)
    best_r = np.full(n, -np.inf)
    if footprints is not None:
        from .footprint import footprint_r_at

        fp_r = footprint_r_at(footprints, lat, lon)
        if not named and len(fp_r.columns) == 1:  # one unnamed reference
            fp_r.columns = [None]
        missing = [str(r) for r, _ in groups if r not in fp_r.columns]
        if missing:
            raise ValueError(
                "No footprint for " + ", ".join(missing) + ": compute "
                "wind_footprints with the same references."
            )
    best_u = np.full(n, np.nan)
    best_sigma = np.full(n, np.nan)
    best_name = np.full(n, None, dtype=object)
    for ref_name, g in groups:
        tr = g["datetime"].to_numpy(dtype="datetime64[ns]")
        if len(tr) > 1:
            it = np.clip(np.searchsorted(tr, t), 1, len(tr) - 1)
            it = np.where(
                np.abs(tr[it - 1] - t) <= np.abs(tr[it] - t), it - 1, it
            )
        else:
            it = np.zeros(n, dtype=int)
        dt = (tr[it] - t).astype("timedelta64[s]").astype(float) / 60
        rlat = g["latitude"].to_numpy(float)[it]
        rlon = g["longitude"].to_numpy(float)[it]
        dist = _haversine_km(lat, lon, rlat, rlon)
        ok = np.abs(dt) <= tol_min
        if footprints is not None:
            r = np.nan_to_num(fp_r[ref_name].to_numpy(float), nan=-np.inf)
            closer = ok & ((r > best_r) | ((r == best_r) & (dist < best_dist)))
            best_r = np.where(closer, r, best_r)
        else:
            closer = ok & (
                (dist < best_dist)
                | ((dist == best_dist) & (np.abs(dt) < np.abs(best_dt)))
            )
        best_dist = np.where(closer, dist, best_dist)
        best_dt = np.where(closer, dt, best_dt)
        best_t = np.where(closer, tr[it], best_t)
        best_lat = np.where(closer, rlat, best_lat)
        best_lon = np.where(closer, rlon, best_lon)
        best_u = np.where(closer, g["wind_speed"].to_numpy(float)[it], best_u)
        if has_sigma:
            sig = g["wind_speed_sigma"].to_numpy(float)[it]
            best_sigma = np.where(closer, sig, best_sigma)
        best_name[closer] = ref_name

    best_dist[~np.isfinite(best_dist)] = np.nan
    if max_distance_km is not None:
        far = ~(best_dist <= max_distance_km)
        best_u[far] = np.nan
        best_sigma[far] = np.nan
    out[col] = best_u
    out[f"{name}_distance_km"] = best_dist
    out[f"{name}_dt_min"] = best_dt
    out[f"{name}_time"] = best_t
    out[f"{name}_latitude"] = best_lat
    out[f"{name}_longitude"] = best_lon
    if footprints is not None:
        out[f"{name}_footprint_r"] = np.where(
            np.isfinite(best_r), best_r, np.nan
        )
    if named:
        out[f"{name}_reference_name"] = best_name
    if has_sigma:
        out[f"{name}_sigma"] = best_sigma
    out.attrs = dict(pal.attrs)
    refs = dict(out.attrs.get("references", {}))
    matched = np.isfinite(best_u)
    refs[name] = {
        "column": col,
        "max_time": str(max_time),
        "max_distance_km": max_distance_km,
        "n_matched": int(matched.sum()),
        "height": ref.attrs.get("height"),
        "matching": "footprint" if footprints is not None else "closest",
    }
    if footprints is not None:
        refs[name]["footprints"] = dict(footprints.attrs)
    if named:
        counts = pd.Series(best_name[matched]).value_counts()
        refs[name]["references"] = {
            str(r): int(counts.get(r, 0)) for r, _ in groups
        }
    # where each reference is (for maps): its most frequent position, and
    # how far its other positions are from it (a ship, or a buoy reported
    # at several positions)
    refs[name]["positions"] = {
        str(name if r is None else r): _position_summary(g) for r, g in groups
    }
    out.attrs["references"] = refs
    return out
