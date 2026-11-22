"""Float positions: GPS fixes from the core Rtraj file, interpolated."""

import numpy as np
import pandas as pd
import xarray as xr

from .read import _strings


def read_positions(fname_or_ds, good_qc=("1", "2")):
    """Read the float's surface positions (launch and GPS fixes).

    Parameters
    ----------
    fname_or_ds : str or xr.Dataset
        Path to `<wmo>_Rtraj.nc`, or the opened dataset.
    good_qc : tuple of str, optional
        POSITION_QC flags to keep, by default ("1", "2").

    Returns
    -------
    pd.DataFrame
        Columns datetime, latitude, longitude, position_qc, cycle,
        measurement_code; sorted by time.
    """
    ds = (
        fname_or_ds
        if isinstance(fname_or_ds, xr.Dataset)
        else xr.open_dataset(fname_or_ds)
    )
    lat = np.asarray(ds["LATITUDE"].values, dtype=float)
    lon = np.asarray(ds["LONGITUDE"].values, dtype=float)
    # Times: adjusted where available, else raw (e.g. the launch position)
    t_raw = ds["JULD"].values
    t = ds["JULD_ADJUSTED"].values if "JULD_ADJUSTED" in ds else t_raw
    t = np.where(pd.isna(t), t_raw, t)
    qc = _strings(ds["POSITION_QC"])
    keep = np.isfinite(lat) & np.isfinite(lon) & ~pd.isna(t)
    keep &= np.isin(qc, good_qc)
    out = pd.DataFrame(
        {
            "datetime": pd.to_datetime(t[keep]),
            "latitude": lat[keep],
            "longitude": lon[keep],
            "position_qc": qc[keep],
            "cycle": np.asarray(ds["CYCLE_NUMBER"].values, float)[keep],
            "measurement_code": np.asarray(
                ds["MEASUREMENT_CODE"].values, float
            )[keep],
        }
    )
    return out.sort_values("datetime").reset_index(drop=True)


def _slerp(lat1, lon1, lat2, lon2, w):
    """Point a fraction w along the great circle from 1 to 2 (degrees)."""
    p1 = np.radians(np.stack([lat1, lon1]))
    p2 = np.radians(np.stack([lat2, lon2]))

    def xyz(p):
        return np.stack(
            [
                np.cos(p[0]) * np.cos(p[1]),
                np.cos(p[0]) * np.sin(p[1]),
                np.sin(p[0]),
            ]
        )

    a, b = xyz(p1), xyz(p2)
    omega = np.arccos(np.clip((a * b).sum(axis=0), -1, 1))
    small = omega < 1e-12
    so = np.where(small, 1, np.sin(omega))
    ka = np.where(small, 1 - w, np.sin((1 - w) * omega) / so)
    kb = np.where(small, w, np.sin(w * omega) / so)
    v = ka * a + kb * b
    lat = np.degrees(np.arctan2(v[2], np.hypot(v[0], v[1])))
    lon = np.degrees(np.arctan2(v[1], v[0]))
    return lat, lon


def attach_positions(pal, positions):
    """Add the float's position at each PAL measurement.

    Positions are interpolated along the great circle between the two
    surface fixes around each measurement, linearly in time.  The float
    drifts at depth between fixes, so this is an estimate; no
    extrapolation before the first or after the last fix.

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table.
    positions : pd.DataFrame
        Output of `read_positions`.

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "latitude" and "longitude" columns, and
        "position_gap_h": hours between the two fixes used (a measure of
        how uncertain the position is).
    """
    out = pal.copy()
    t = out["datetime"].to_numpy(dtype="datetime64[ns]").astype("i8")
    tf = positions["datetime"].to_numpy(dtype="datetime64[ns]").astype("i8")
    lat = np.full(len(t), np.nan)
    lon = np.full(len(t), np.nan)
    gap = np.full(len(t), np.nan)
    if len(tf) >= 2:
        right = np.clip(np.searchsorted(tf, t, side="left"), 1, len(tf) - 1)
        left = right - 1
        inside = (t >= tf[left]) & (t <= tf[right])
        span = (tf[right] - tf[left]).astype(float)
        w = np.where(span > 0, (t - tf[left]) / np.where(span > 0, span, 1), 0)
        la, lo = _slerp(
            positions["latitude"].to_numpy()[left],
            positions["longitude"].to_numpy()[left],
            positions["latitude"].to_numpy()[right],
            positions["longitude"].to_numpy()[right],
            w,
        )
        lat[inside], lon[inside] = la[inside], lo[inside]
        gap[inside] = span[inside] / 3.6e12
    out["latitude"] = lat
    out["longitude"] = lon
    out["position_gap_h"] = gap
    out.attrs = dict(pal.attrs)
    return out
