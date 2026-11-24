"""Wind footprints of reference stations, from a reanalysis.

Near the coast, the wind at a buoy can be a poor proxy for the wind a few
tens of kilometres away (fetch, orography), and a buoy further away can
be a better one.  Instead of a fixed distance, a reference can be
attached to the float where the wind varies in the same way as at the
reference: its footprint.

For each reference, the correlation R between the reanalysis wind at the
reference's position and the wind at every grid cell is computed over the
period given.  The footprint is where R is at least a threshold (0.9 by
default).  Each float measurement is then compared with the reference
whose footprint has the highest R at the float's position.

ArgoPAL does not download reanalyses: give a file or dataset (e.g. ERA5
10 m wind, "u10" and "v10") covering the float's area and period.
"""

import numpy as np
import pandas as pd
import xarray as xr


TIME_NAMES = ("time", "valid_time", "date", "datetime")
LAT_NAMES = ("latitude", "lat")
LON_NAMES = ("longitude", "lon")
SPEED_NAMES = ("si10", "ws10", "wind_speed", "ws", "wspd", "u10_speed")


def _first(names, available, what):
    for n in names:
        if n in available:
            return n
    raise KeyError(
        f"No {what} in the reanalysis; looked for {', '.join(names)}."
    )


def read_reanalysis(source, variable=None):
    """Wind speed (m s-1) from a reanalysis, on (time, latitude, longitude).

    Parameters
    ----------
    source : str, list of str, xr.Dataset or xr.DataArray
        NetCDF file(s) or dataset.  The wind speed is computed from "u10"
        and "v10" (ERA5), or read from a speed variable ("si10", "ws10",
        "wind_speed", ...), or is the DataArray itself.
    variable : str, optional
        Name of the wind-speed variable, if not found automatically.

    Returns
    -------
    xr.DataArray
        "wind_speed" on (time, latitude, longitude), longitudes in
        [-180, 180), sorted.
    """
    if isinstance(source, (str, list, tuple)):
        files = [source] if isinstance(source, str) else list(source)
        parts = [xr.open_dataset(f) for f in files]
        tname = _first(TIME_NAMES, parts[0].dims, "time dimension")
        ds = xr.concat(parts, dim=tname) if len(parts) > 1 else parts[0]
    else:
        ds = source
    if isinstance(ds, xr.DataArray):
        ws = ds
    elif variable is not None:
        ws = ds[variable]
    elif "u10" in ds and "v10" in ds:
        ws = np.hypot(ds["u10"], ds["v10"])
    else:
        ws = ds[_first(SPEED_NAMES, ds.data_vars, "wind speed variable")]
    dims = set(ws.dims)
    tname = _first(TIME_NAMES, dims, "time dimension")
    yname = _first(LAT_NAMES, dims, "latitude dimension")
    xname = _first(LON_NAMES, dims, "longitude dimension")
    extra = [d for d in ws.dims if d not in (tname, yname, xname)]
    if extra:  # e.g. a single pressure level or ensemble member
        ws = ws.isel({d: 0 for d in extra}, drop=True)
    ws = ws.rename({tname: "time", yname: "latitude", xname: "longitude"})
    ws = ws.assign_coords(
        longitude=((ws["longitude"] + 180) % 360) - 180
    ).sortby(["time", "latitude", "longitude"])
    ws = ws.transpose("time", "latitude", "longitude").astype(float)
    ws.name = "wind_speed"
    return ws


def _reference_positions(references):
    """name -> (lat, lon), from a dict or a reference table."""
    if isinstance(references, dict) and all(
        isinstance(v, (tuple, list)) and len(v) == 2
        for v in references.values()
    ):
        return {
            str(k): (float(v[0]), float(v[1])) for k, v in references.items()
        }
    import warnings

    from .reference import _position_summary, read_reference

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # heights do not matter here
        table = read_reference(references)
    groups = (
        table.groupby("reference_name") if "reference_name" in table
        else [("reference", table)]
    )  # fmt: skip
    out = {}
    for name, g in groups:
        pos = _position_summary(g)
        if pos["latitude"] is not None:
            out[str(name)] = (pos["latitude"], pos["longitude"])
    return out


def wind_footprints(reanalysis, references, anomaly=None, variable=None):
    """Correlation of the reanalysis wind with the wind at each reference.

    Parameters
    ----------
    reanalysis : str, list, xr.Dataset or xr.DataArray
        Reanalysis wind covering the float's area and period (see
        `read_reanalysis`).
    references : dict or pd.DataFrame or str
        name -> (latitude, longitude), or a reference wind table / CSV
        (each reference at its most frequent position).
    anomaly : {None, "monthly"}, optional
        None (default): correlate the wind as given (e.g. hourly) over the
        period given, typically the deployment: where the wind followed
        the reference's during that period.  "monthly": correlate monthly
        anomalies (monthly means minus each calendar month's mean) over
        several years: where the wind varies with the reference's month to
        month, a climatological footprint.  Needs each calendar month at
        least twice.  Only the reanalysis is used, at each reference's
        position, not the reference's own wind.
    variable : str, optional
        Wind-speed variable, if not found automatically.

    Returns
    -------
    xr.DataArray
        "footprint_r" on (reference, latitude, longitude): R between the
        wind at each reference's grid cell and at every grid cell.
        `attrs` record the positions, period and anomaly used.
    """
    ws = read_reanalysis(reanalysis, variable=variable)
    if anomaly == "monthly":
        monthly = ws.resample(time="1MS").mean()
        per_month = pd.Series(
            pd.to_datetime(monthly["time"].values).month
        ).value_counts()
        if per_month.min() < 2:
            raise ValueError(
                'anomaly="monthly" needs each calendar month at least twice '
                "(several years of data): with "
                f"{len(monthly['time'])} months, the monthly anomalies are "
                "not defined.  For a reanalysis over the deployment only, "
                "use anomaly=None (the wind as given, e.g. hourly)."
            )
        ws = monthly.groupby("time.month") - monthly.groupby(
            "time.month"
        ).mean("time")
        ws = ws.drop_vars("month", errors="ignore")
    elif anomaly is not None:
        raise ValueError(
            f'anomaly must be None or "monthly"; got {anomaly!r}.'
        )
    positions = _reference_positions(references)
    maps = []
    for name, (lat, lon) in positions.items():
        at = ws.sel(latitude=lat, longitude=lon, method="nearest")
        maps.append(xr.corr(at, ws, dim="time").expand_dims(reference=[name]))
    r = xr.concat(maps, dim="reference")
    r.name = "footprint_r"
    times = pd.to_datetime(ws["time"].values)
    r.attrs = {
        "positions": {k: list(v) for k, v in positions.items()},
        "period": [str(times.min()), str(times.max())],
        "anomaly": str(anomaly),
    }
    return r


def footprint_r_at(footprints, latitude, longitude):
    """R of each reference's footprint at given positions (nearest grid
    cell; NaN outside the grid).

    Returns
    -------
    pd.DataFrame
        One column per reference, one row per position.
    """
    lat = np.asarray(latitude, dtype=float)
    lon = ((np.asarray(longitude, dtype=float) + 180) % 360) - 180
    glat = footprints["latitude"].to_numpy()
    glon = footprints["longitude"].to_numpy()
    out = {}
    iy = _nearest(glat, lat)
    ix = _nearest(glon, lon)
    ok = (iy >= 0) & (ix >= 0)
    for name in footprints["reference"].to_numpy():
        field = footprints.sel(reference=name).to_numpy()
        vals = np.full(lat.shape, np.nan)
        vals[ok] = field[iy[ok], ix[ok]]
        out[str(name)] = vals
    return pd.DataFrame(out)


def _nearest(grid, values):
    """Index of the nearest grid value, -1 beyond half a step outside."""
    order = np.argsort(grid)
    g = grid[order]
    step = np.median(np.diff(g)) if len(g) > 1 else np.inf
    i = np.clip(np.searchsorted(g, values), 1, max(len(g) - 1, 1))
    left = np.abs(values - g[i - 1]) <= np.abs(
        values - g[np.minimum(i, len(g) - 1)]
    )
    j = np.where(left, i - 1, np.minimum(i, len(g) - 1))
    far = ~np.isfinite(values) | (np.abs(values - g[j]) > step / 2 + 1e-9)
    return np.where(far, -1, order[j])
