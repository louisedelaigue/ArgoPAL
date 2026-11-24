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

Each PAL measurement is matched to the closest reference (in distance)
among those with a value within `max_time` of it; for each reference, its
value nearest in time is taken.  The reference should be the wind at 10 m,
like the models.
"""

import numpy as np
import pandas as pd


TABLE_ALIASES = {
    "datetime": ("datetime", "time", "date", "time_utc", "validity_time"),
    "wind_speed": ("wind_speed", "wind", "windspeed", "u", "u10", "ff", "ws"),
    "latitude": ("latitude", "lat"),
    "longitude": ("longitude", "lon", "long"),
    "reference_name": ("reference_name", "station", "buoy_name", "platform"),
    "wind_speed_sigma": ("wind_speed_sigma", "wind_sigma", "sigma"),
}
OPTIONAL = ("reference_name", "wind_speed_sigma")


def _haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (
        np.sin((p2 - p1) / 2) ** 2
        + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    )
    return 2 * r * np.arcsin(np.sqrt(a))


def read_reference(source, latitude=None, longitude=None):
    """Read a reference wind table (see module docstring).

    Parameters
    ----------
    source : pd.DataFrame or str
        A table or a CSV file.
    latitude, longitude : float, optional
        Position of a fixed station, if the table has no position columns.

    Returns
    -------
    pd.DataFrame
        Columns datetime (UTC), wind_speed, latitude, longitude, and
        reference_name / wind_speed_sigma when the table has them; sorted
        by time, without missing winds.
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
    out = out.dropna(subset=["datetime", "wind_speed"])
    cols = ["datetime", "wind_speed", "latitude", "longitude"]
    out = out[cols + [c for c in OPTIONAL if c in out]]
    return out.sort_values("datetime").reset_index(drop=True)


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
):
    """Add a reference wind to each PAL measurement.

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

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "U_<name>" (m s-1), "<name>_distance_km",
        "<name>_dt_min" (reference time minus measurement time, minutes),
        "<name>_time" (time of the reference value used),
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

    ref = read_reference(reference, latitude=latitude, longitude=longitude)
    named = "reference_name" in ref
    groups = (
        list(ref.groupby("reference_name", sort=True)) if named
        else [(None, ref)]
    )  # fmt: skip
    has_sigma = "wind_speed_sigma" in ref

    best_dist = np.full(n, np.inf)
    best_dt = np.full(n, np.nan)
    best_t = np.full(n, np.datetime64("NaT", "ns"))
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
        dist = _haversine_km(
            lat, lon, g["latitude"].to_numpy(float)[it],
            g["longitude"].to_numpy(float)[it],
        )  # fmt: skip
        ok = np.abs(dt) <= tol_min
        closer = ok & (
            (dist < best_dist)
            | ((dist == best_dist) & (np.abs(dt) < np.abs(best_dt)))
        )
        best_dist = np.where(closer, dist, best_dist)
        best_dt = np.where(closer, dt, best_dt)
        best_t = np.where(closer, tr[it], best_t)
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
    }
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
