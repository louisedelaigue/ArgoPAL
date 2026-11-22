"""Optional ship-noise check with AIS vessel positions.

Ship noise can dominate the ambient sound field in the 1-10 kHz band up
to ~20 km away.  Following Delaigue et al. (2026, Sect. 2.5), a PAL
measurement is flagged when an AIS-reported vessel is within `radius_km`
of the float's (interpolated) position and within `window` of the
measurement time.

AIS input
---------
`read_ais` accepts, and turns into one table with columns
``datetime`` (UTC), ``mmsi``, ``latitude``, ``longitude``:

1. a pandas DataFrame with those four columns;
2. CSV files of decoded AIS positions (one row per position report).
   Columns are recognised by name (case-insensitive), e.g.
   time: datetime, time, timestamp_utc, basedatetime, "# timestamp",
         tagblock_timestamp (UNIX seconds), unix_time;
   latitude: latitude, lat, y;  longitude: longitude, lon, long, x;
   vessel: mmsi.
   This covers the decoded files of the 2025 DYFAMED study
   (``ais_positions_YYYYMMDD.csv``), NOAA MarineCadastre
   (BaseDateTime, LAT, LON, MMSI) and Danish Maritime Authority exports
   (# Timestamp, Latitude, Longitude, MMSI).
3. raw NMEA 0183 files (``!AIVDM`` sentences, plain or ``.gz``) with
   NMEA 4.0 tag blocks giving the receiver time (``\\c:1739500274,...\\``),
   as delivered by most shore-station networks.  Decoding needs the
   optional ``pyais`` package: ``pip install argopal[ais]``.

Only vessel position reports (message types 1, 2, 3, 18, 19, 27) are
used: in CSV files with a message-type column (msg_type, message_type,
message_id, id), other rows such as base stations (type 4) and aids to
navigation (type 21) are dropped.

Time: the time of each position must be the time it was received (UTC).
Do not rebuild dates from file names: daily AIS files often span two
calendar days.
"""

import glob
import gzip
import os
import warnings

import numpy as np
import pandas as pd


POSITION_TYPES = (1, 2, 3, 18, 19, 27)

COLUMN_ALIASES = {
    "datetime": (
        "datetime",
        "time",
        "time_utc",
        "timestamp_utc",
        "basedatetime",
        "# timestamp",
        "tagblock_timestamp",
        "unix_time",
        "epoch",
        "timestamp",
    ),
    "mmsi": ("mmsi",),
    "latitude": ("latitude", "lat", "y"),
    "longitude": ("longitude", "lon", "long", "x"),
}
# A message-type column, if present, is used to keep vessel position
# reports only (not base stations, type 4, or aids to navigation, type 21).
MESSAGE_TYPE_ALIASES = (
    "msg_type",
    "message_type",
    "message_id",
    "msgtype",
    "id",
)

# Numeric time columns must be UNIX seconds; anything below this (e.g. the
# AIS "timestamp" field, the UTC second 0-59) is rejected.
_MIN_UNIX_TIME = 1e8  # 1973


def _empty():
    return pd.DataFrame(
        {
            "datetime": pd.Series(dtype="datetime64[ns]"),
            "mmsi": pd.Series(dtype="int64"),
            "latitude": pd.Series(dtype=float),
            "longitude": pd.Series(dtype=float),
        }
    )


def _minutes(times):
    """Sorted unique UTC minutes (int64, minutes since 1970) of `times`."""
    t = pd.to_datetime(pd.Series(times)).dropna().astype("datetime64[ns]")
    return np.unique(t.to_numpy().astype("i8") // 60_000_000_000)


def _clean(df, bbox=None):
    """Valid positions only, inside bbox (west, east, south, north).

    `attrs["received_minutes"]` records every minute in which the AIS feed
    received a position anywhere (before the bbox cut), to tell "no ship
    nearby" from "no AIS coverage"."""
    df = df.dropna(subset=["datetime", "mmsi", "latitude", "longitude"])
    df = df[
        df["latitude"].between(-90, 90) & df["longitude"].between(-180, 180)
    ]
    received = _minutes(df["datetime"])
    if bbox is not None:
        w, e, s, n = bbox
        df = df[df["longitude"].between(w, e) & df["latitude"].between(s, n)]
    df = df.astype({"mmsi": "int64"})
    df["datetime"] = pd.to_datetime(df["datetime"]).astype("datetime64[ns]")
    df = df[["datetime", "mmsi", "latitude", "longitude"]].copy()
    df.attrs = {"received_minutes": received.tolist()}
    return df


def _is_nmea(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        first = f.readline().lstrip()
    return first[:1] in (b"\\", b"!", b"$")


def _find_columns(columns):
    lower = {c.strip().lower(): c for c in columns}
    found = {}
    for target, aliases in COLUMN_ALIASES.items():
        for a in aliases:
            if a in lower:
                found[target] = lower[a]
                break
    missing = [t for t in COLUMN_ALIASES if t not in found]
    return found, missing


def read_ais_csv(path, bbox=None):
    """Read decoded AIS positions from a CSV file (see module docstring).

    Parameters
    ----------
    path : str
        CSV file (may be compressed, e.g. ``.csv.gz``).
    bbox : (west, east, south, north), optional
        Keep only positions inside this box (degrees).

    Returns
    -------
    pd.DataFrame
        Columns datetime, mmsi, latitude, longitude.
    """
    header = pd.read_csv(path, nrows=0).columns
    found, missing = _find_columns(header)
    lower = {c.strip().lower(): c for c in header}
    type_col = next(
        (lower[a] for a in MESSAGE_TYPE_ALIASES if a in lower), None
    )
    if missing:
        raise ValueError(
            f"{os.path.basename(path)}: no column for {', '.join(missing)}. "
            "Accepted names (case-insensitive): "
            + "; ".join(
                f"{t}: {', '.join(COLUMN_ALIASES[t])}" for t in missing
            )
        )
    cols = list(found.values()) + ([type_col] if type_col else [])
    df = pd.read_csv(path, usecols=cols, low_memory=False)
    if type_col is not None:
        types = pd.to_numeric(df[type_col], errors="coerce")
        if types.dropna().between(1, 27).all():  # AIS message types
            df = df[types.isin(POSITION_TYPES)]
        df = df.drop(columns=type_col)
    df = df.rename(columns={v: k for k, v in found.items()})
    t = df["datetime"]
    if pd.api.types.is_numeric_dtype(t):
        if t.max() < _MIN_UNIX_TIME:
            raise ValueError(
                f"{os.path.basename(path)}: time column {found['datetime']!r} "
                "holds small numbers, not UNIX seconds (the AIS 'timestamp' "
                "field is only the second of the minute). Provide a full "
                "date/time or UNIX-seconds column, e.g. tagblock_timestamp."
            )
        df["datetime"] = pd.to_datetime(t, unit="s", errors="coerce")
    else:
        df["datetime"] = pd.to_datetime(t, errors="coerce", utc=True)
        df["datetime"] = df["datetime"].dt.tz_localize(None)
    for c in ("latitude", "longitude", "mmsi"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return _clean(df, bbox)


def read_ais_nmea(path, bbox=None):
    """Decode AIS position reports from a raw NMEA file with tag blocks.

    Each sentence must carry an NMEA 4.0 tag block with the receiver time
    (``c:`` field, UNIX seconds), e.g.
    ``\\c:1739500274,s:2954*14\\!AIVDM,1,1,,A,14hEMK0vP0PPaDpHt?1f49RL0H8A,0*3D``.
    Needs the optional ``pyais`` package (``pip install argopal[ais]``).

    Parameters
    ----------
    path : str
        NMEA file, plain text or ``.gz``.
    bbox : (west, east, south, north), optional
        Keep only positions inside this box (degrees).

    Returns
    -------
    pd.DataFrame
        Columns datetime, mmsi, latitude, longitude.
    """
    try:
        from pyais.exceptions import AISBaseException
        from pyais.stream import IterMessages
    except ImportError as e:
        raise ImportError(
            "Reading raw NMEA AIS files needs the pyais package: "
            "pip install argopal[ais]  (or: pip install pyais). "
            "Alternatively, give decoded positions as a CSV or DataFrame."
        ) from e
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        lines = [line.rstrip(b"\r\n") for line in f if line.strip()]
    rows, no_time, undecodable = [], 0, 0
    for msg in IterMessages(lines):
        try:
            d = msg.decode()
        except (AISBaseException, ValueError):  # corrupt or unsupported
            undecodable += 1
            continue
        if d.msg_type not in POSITION_TYPES:
            continue
        tb = msg.tag_block
        t = None
        if tb is not None:
            tb.init()
            t = tb.receiver_timestamp
        if not t:
            no_time += 1
            continue
        rows.append((int(t), d.mmsi, d.lat, d.lon))
    if no_time and not rows:
        raise ValueError(
            f"{os.path.basename(path)}: no tag-block receiver times (c:) "
            "found; argopal needs them to time AIS positions."
        )
    if not rows:
        return _empty()
    if undecodable > 0.5 * (len(rows) + undecodable):
        warnings.warn(
            f"{os.path.basename(path)}: {undecodable} AIS sentences could "
            "not be decoded.",
            stacklevel=2,
        )
    df = pd.DataFrame(rows, columns=["unix", "mmsi", "latitude", "longitude"])
    df["datetime"] = pd.to_datetime(df.pop("unix"), unit="s")
    return _clean(df, bbox)


def read_ais(source, bbox=None, verbose=False):
    """Read AIS vessel positions from any supported source.

    Parameters
    ----------
    source : pd.DataFrame, str or list of str
        A table (columns datetime, mmsi, latitude, longitude), a file, a
        folder, a glob pattern (e.g. "AIS/messages-*.gz") or a list of
        files.  CSV and raw NMEA files are told apart by their content.
    bbox : (west, east, south, north), optional
        Keep only positions inside this box (degrees); reading many days of
        AIS is much lighter with it.
    verbose : bool, optional
        Print one line per file.

    Returns
    -------
    pd.DataFrame
        Columns datetime (UTC), mmsi, latitude, longitude; sorted by time.
    """
    if isinstance(source, pd.DataFrame):
        found, missing = _find_columns(source.columns)
        if missing:
            raise ValueError(
                f"AIS table has no column for {', '.join(missing)}; it needs "
                "datetime, mmsi, latitude, longitude."
            )
        df = source.rename(columns={v: k for k, v in found.items()})
        out = _clean(df, bbox)
        received = out.attrs["received_minutes"]
        out = out.sort_values("datetime").reset_index(drop=True)
        out.attrs = {"received_minutes": received}
        return out
    if isinstance(source, (str, os.PathLike)):
        source = str(source)
        if os.path.isdir(source):
            files = sorted(
                os.path.join(source, f)
                for f in os.listdir(source)
                if not f.startswith(".")
            )
        else:
            files = sorted(glob.glob(source)) or [source]
    else:
        files = [str(f) for f in source]
    tables = []
    for path in files:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        reader = read_ais_nmea if _is_nmea(path) else read_ais_csv
        df = reader(path, bbox=bbox)
        if verbose:
            print(f"{os.path.basename(path)}: {len(df):,} positions")
        tables.append(df)
    if not tables:
        return _empty()
    received = (
        np.unique(
            np.concatenate(
                [
                    np.asarray(t.attrs.get("received_minutes", []), "int64")
                    for t in tables
                ]
            )
        )
        .astype("int64")
        .tolist()
    )
    for t in tables:  # pandas compares attrs when concatenating
        t.attrs = {}
    out = pd.concat(tables, ignore_index=True)
    out = out.sort_values("datetime").reset_index(drop=True)
    out.attrs = {"received_minutes": received}
    return out


def track_bbox(pal, margin_km=20.0):
    """Box (west, east, south, north) around the float's positions, widened
    by `margin_km`, to read only the AIS data that matter."""
    lat = pal["latitude"].to_numpy(float)
    lon = pal["longitude"].to_numpy(float)
    if not np.isfinite(lat).any():
        raise ValueError("No float positions: run attach_positions first.")
    dlat = margin_km / 111.0
    mid = np.nanmean(lat)
    dlon = margin_km / (111.0 * np.cos(np.radians(mid)))
    return (
        np.nanmin(lon) - dlon,
        np.nanmax(lon) + dlon,
        np.nanmin(lat) - dlat,
        np.nanmax(lat) + dlat,
    )


def _haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def flag_ships(pal, ais, radius_km=20.0, window="30min", action="flag"):
    """Flag PAL measurements with a ship nearby (AIS).

    A measurement is flagged when at least one vessel (counted once per
    MMSI) reported a position within `radius_km` of the float's position
    and within ±`window` of the measurement time.  Where no AIS message at
    all exists within ±`window` (no AIS coverage), the result is unknown
    (NaN), not "no ship".

    Parameters
    ----------
    pal : pd.DataFrame
        PAL table after `attach_positions` (needs latitude, longitude).
    ais : pd.DataFrame
        Output of `read_ais`.
    radius_km : float, optional
        By default 20 km.
    window : str or pd.Timedelta, optional
        Half-width of the time window, by default "30min".
    action : {"flag", "remove"}, optional
        "flag" (default) only adds the columns below.  "remove" also sets
        the acoustic values (f_*, spl_*, U_* columns) of flagged
        measurements to NaN.

    Returns
    -------
    pd.DataFrame
        A copy of `pal` with "ships_nearby" (number of vessels; NaN without
        AIS coverage or float position) and "ship_flag" (True / False; False
        also where unknown, see "ships_nearby").
    """
    if action not in ("flag", "remove"):
        raise ValueError("action must be 'flag' or 'remove'.")
    half = pd.Timedelta(window).value
    out = pal.copy()
    t_ais = ais["datetime"].to_numpy(dtype="datetime64[ns]").astype("i8")
    order = np.argsort(t_ais, kind="stable")
    t_ais = t_ais[order]
    lat_a = ais["latitude"].to_numpy(float)[order]
    lon_a = ais["longitude"].to_numpy(float)[order]
    mmsi = ais["mmsi"].to_numpy()[order]

    t = out["datetime"].to_numpy(dtype="datetime64[ns]").astype("i8")
    lat = out["latitude"].to_numpy(float) if "latitude" in out else None
    if lat is None:
        raise KeyError("No float positions: run attach_positions first.")
    lon = out["longitude"].to_numpy(float)
    lo = np.searchsorted(t_ais, t - half, side="left")
    hi = np.searchsorted(t_ais, t + half, side="right")
    # AIS coverage: was anything received (anywhere) within the window?
    received = np.asarray(ais.attrs.get("received_minutes", []), "int64")
    if not len(received):
        received = _minutes(ais["datetime"])
    half_min = half // 60_000_000_000
    t_min = t // 60_000_000_000
    covered = np.searchsorted(
        received, t_min + half_min, side="right"
    ) > np.searchsorted(received, t_min - half_min, side="left")
    count = np.full(len(out), np.nan)
    count[covered & np.isfinite(lat) & np.isfinite(lon)] = 0
    for i in np.flatnonzero((hi > lo) & np.isfinite(lat) & np.isfinite(lon)):
        sl = slice(lo[i], hi[i])
        near = _haversine_km(lat[i], lon[i], lat_a[sl], lon_a[sl]) <= radius_km
        count[i] = len(np.unique(mmsi[sl][near])) if near.any() else 0
    out["ships_nearby"] = count
    out["ship_flag"] = count > 0
    n_flag = int(out["ship_flag"].sum())
    n_unknown = int(np.isnan(count).sum())
    if action == "remove" and n_flag:
        cols = [c for c in out if c.startswith(("f_", "spl_", "U_"))]
        out.loc[out["ship_flag"], cols] = np.nan
    if n_unknown == len(out):
        warnings.warn(
            "No AIS message overlaps the PAL measurements in time: check "
            "the AIS files' dates.",
            stacklevel=2,
        )
    out.attrs = dict(pal.attrs)
    out.attrs["ais"] = {
        "radius_km": float(radius_km),
        "window": str(window),
        "action": action,
        "n_positions": len(ais),
        "n_vessels": int(ais["mmsi"].nunique()),
        "n_flagged": n_flag,
        "n_without_coverage": n_unknown,
    }
    return out
