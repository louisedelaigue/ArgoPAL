"""Read core Argo files (Sprof, Rtraj) with a chosen data mode."""

import os

import numpy as np
import pandas as pd
import xarray as xr

from .read import _adjusted_or_raw, _strings, _wmo


# Data modes a user can ask for.  Each mode takes values from one source
# only and never falls back to another; "delayed" is the default.
MODES = ("delayed", "adjusted", "raw")
DEFAULT_MODE = "delayed"


def _check_mode(mode):
    if mode not in MODES:
        raise ValueError(
            f"mode must be one of {', '.join(MODES)}; got {mode!r}."
        )


def select_mode(raw, adjusted, raw_qc, adjusted_qc, data_mode, mode):
    """Pick values according to the requested data mode.

    All inputs are arrays of the same shape.  `data_mode` holds the Argo
    data mode letter of each value ("R", "A", "D", or "" if unknown).

    Each mode takes values from one source only; there is no fallback.

    Modes
    -----
    "delayed"
        PARAM_ADJUSTED where the data mode is "D"; NaN everywhere else.
    "adjusted"
        PARAM_ADJUSTED where the data mode is "A" (real-time adjusted);
        NaN everywhere else, including delayed-mode values.
    "raw"
        PARAM, the unadjusted values, which Argo keeps for every data mode.

    Where the data centre removed an adjusted value (PARAM_ADJUSTED empty,
    QC "4"), the result is NaN; the raw value is never used instead.

    Returns
    -------
    values : np.ndarray of float
    qc : np.ndarray of str
        The QC flag that goes with the chosen value.  Where an adjusted
        value was removed by the data centre, the value is NaN but its
        adjusted QC flag (e.g. "4") is kept, so the reason stays visible.
        "" where the mode does not apply.
    used_mode : np.ndarray of str
        Data mode of the chosen value: "D", "A" (adjusted) or "R" (raw);
        "" where the mode does not apply.
    """
    _check_mode(mode)
    raw = np.asarray(raw, dtype=float)
    adjusted = np.asarray(adjusted, dtype=float)
    data_mode = np.asarray(data_mode, dtype=str)
    if mode == "delayed":
        take_adj = data_mode == "D"
        take_raw = np.zeros(raw.shape, dtype=bool)
    elif mode == "adjusted":
        take_adj = data_mode == "A"
        take_raw = np.zeros(raw.shape, dtype=bool)
    else:  # raw
        take_adj = np.zeros(raw.shape, dtype=bool)
        take_raw = np.ones(raw.shape, dtype=bool)

    values = np.full(raw.shape, np.nan)
    values[take_adj] = adjusted[take_adj]
    values[take_raw] = raw[take_raw]
    raw_qc = np.asarray(raw_qc, dtype=object)
    adjusted_qc = np.asarray(adjusted_qc, dtype=object)
    qc = np.full(raw.shape, "", dtype=object)
    qc[take_adj] = adjusted_qc[take_adj]
    qc[take_raw] = raw_qc[take_raw]
    used = np.full(raw.shape, "", dtype=object)
    used[take_adj] = data_mode[take_adj]
    used[take_raw] = "R"
    # Nothing to report where neither the value nor a QC flag exists
    empty = ~np.isfinite(values) & (qc == "")
    used[empty] = ""
    return values, qc.astype(str), used.astype(str)


def _qc(ds, name, shape):
    if name in ds:
        return _strings(ds[name])
    return np.full(shape, "", dtype=str)


def _pick(ds, var, data_mode, mode):
    raw = np.asarray(ds[var].values, dtype=float)
    adj_name = var + "_ADJUSTED"
    if adj_name in ds:
        adjusted = np.asarray(ds[adj_name].values, dtype=float)
    else:
        adjusted = np.full(raw.shape, np.nan)
    return select_mode(
        raw,
        adjusted,
        _qc(ds, var + "_QC", raw.shape),
        _qc(ds, adj_name + "_QC", raw.shape),
        data_mode,
        mode,
    )


def _open(fname_or_ds):
    if isinstance(fname_or_ds, xr.Dataset):
        return fname_or_ds
    return xr.open_dataset(fname_or_ds)


# --- Profiles (Sprof) ------------------------------------------------------


def _profile_parameters(ds):
    """Per-profile parameter names and their data modes (N_PROF, N_PARAM)."""
    names = _strings(ds["STATION_PARAMETERS"])
    modes = _strings(ds["PARAMETER_DATA_MODE"])
    return names, modes


def _profile_data_mode(names, modes, var):
    """Data mode of `var` in each profile ("" if not measured)."""
    out = np.full(names.shape[0], "", dtype=object)
    for i in range(names.shape[0]):
        hits = np.flatnonzero(names[i] == var)
        if hits.size:
            out[i] = modes[i, hits[0]]
    return out.astype(str)


def profile_variables(fname_or_ds):
    """Parameters measured in the profiles of an Sprof file (incl. PRES)."""
    ds = _open(fname_or_ds)
    names, _ = _profile_parameters(ds)
    found = [n for n in dict.fromkeys(names.ravel()) if n]
    return [v for v in found if v in ds]


def read_profiles(fname_or_ds, mode=DEFAULT_MODE, variables=None):
    """Read an Sprof file into a table, one row per profile level.

    Parameters
    ----------
    fname_or_ds : str or xr.Dataset
        Path to `<wmo>_Sprof.nc`, or the opened dataset.
    mode : str, optional
        Data mode to use: "delayed" (default), "adjusted" or "raw";
        see `select_mode`.  Pressure uses the same mode, and levels without
        a pressure in that mode are dropped.
    variables : list of str, optional
        Parameters to keep, e.g. ["TEMP", "PSAL", "DOXY"].  By default all
        parameters in the file.  Names not in the file are ignored.

    Returns
    -------
    pd.DataFrame
        Columns: cycle, direction, datetime, latitude, longitude,
        position_qc, then PRES and each parameter, each followed by
        <PARAM>_qc and <PARAM>_data_mode.
    """
    _check_mode(mode)
    ds = _open(fname_or_ds)
    names, modes = _profile_parameters(ds)
    available = [v for v in profile_variables(ds) if v != "PRES"]
    if variables is None:
        variables = available
    else:
        variables = [v for v in variables if v in available]

    n_prof, n_lev = ds.sizes["N_PROF"], ds.sizes["N_LEVELS"]
    cycle, _ = _adjusted_or_raw(ds, "CYCLE_NUMBER")
    meta = {
        "cycle": np.asarray(cycle, dtype=float),
        "direction": _strings(ds["DIRECTION"]),
        "datetime": pd.to_datetime(ds["JULD"].values),
        "latitude": np.asarray(ds["LATITUDE"].values, dtype=float),
        "longitude": np.asarray(ds["LONGITUDE"].values, dtype=float),
        "position_qc": _strings(ds["POSITION_QC"]),
    }
    df = pd.DataFrame(
        {k: np.repeat(np.asarray(v), n_lev) for k, v in meta.items()}
    )

    for var in ["PRES"] + variables:
        dm = _profile_data_mode(names, modes, var)
        dm = np.repeat(dm[:, None], n_lev, axis=1)
        values, qc, used = _pick(ds, var, dm, mode)
        df[var] = values.ravel()
        df[var + "_qc"] = qc.ravel()
        df[var + "_data_mode"] = used.ravel()

    keep = np.isfinite(df["PRES"].to_numpy())
    if variables:
        keep &= np.isfinite(df[variables].to_numpy()).any(axis=1)
    df = df.loc[keep].reset_index(drop=True)
    df["cycle"] = df["cycle"].astype("Int64")
    df.attrs = {
        "wmo": _wmo(ds),
        "mode": mode,
        "variables": variables,
        "n_profiles": int(n_prof),
    }
    return df


# --- Park drift measurements (Rtraj) ---------------------------------------


def trajectory_variables(fname_or_ds):
    """Parameters in an Rtraj file that hold at least one value."""
    ds = _open(fname_or_ds)
    names = [n for n in _strings(ds["TRAJECTORY_PARAMETERS"]) if n]
    out = []
    for n in names:
        if n in ds and np.isfinite(np.asarray(ds[n].values, float)).any():
            out.append(n)
    return out


def read_park(
    fname_or_ds,
    mode=DEFAULT_MODE,
    variables=None,
    measurement_codes=(290,),
):
    """Read the park-drift measurements of an Rtraj file.

    Parameters
    ----------
    fname_or_ds : str or xr.Dataset
        Path to `<wmo>_Rtraj.nc`, or the opened dataset.
    mode : str, optional
        Data mode: "delayed" (default), "adjusted" or "raw".
    variables : list of str, optional
        Parameters to keep; by default all parameters with data.  PRES is
        always included.
    measurement_codes : tuple of int, optional
        Argo measurement codes to keep, by default (290,) = park drift.

    Returns
    -------
    pd.DataFrame
        Columns: cycle, measurement_code, datetime, then PRES and each
        parameter, each followed by <PARAM>_qc and <PARAM>_data_mode.  One
        row per measurement; different sensors are often measured at
        different times, so most rows hold only some parameters.
    """
    _check_mode(mode)
    ds = _open(fname_or_ds)
    params = list(_strings(ds["TRAJECTORY_PARAMETERS"]))
    available = [v for v in trajectory_variables(ds) if v != "PRES"]
    if variables is None:
        variables = available
    else:
        variables = [v for v in variables if v in available]

    mc = np.asarray(ds["MEASUREMENT_CODE"].values, dtype=float)
    rows = np.isin(mc, measurement_codes)
    datetime, _ = _adjusted_or_raw(ds, "JULD")
    cycle, _ = _adjusted_or_raw(ds, "CYCLE_NUMBER")
    df = pd.DataFrame(
        {
            "cycle": np.asarray(cycle, dtype=float)[rows],
            "measurement_code": mc[rows],
            "datetime": pd.to_datetime(np.asarray(datetime)[rows]),
        }
    )
    all_modes = _strings(ds["TRAJECTORY_PARAMETER_DATA_MODE"])
    for var in ["PRES"] + variables:
        dm = all_modes[:, params.index(var)]
        values, qc, used = _pick(ds, var, dm, mode)
        df[var] = values[rows]
        df[var + "_qc"] = qc[rows]
        df[var + "_data_mode"] = used[rows]
    df["cycle"] = df["cycle"].astype("Int64")
    df["measurement_code"] = df["measurement_code"].astype("Int64")
    df = df.sort_values("datetime").reset_index(drop=True)
    df.attrs = {"wmo": _wmo(ds), "mode": mode, "variables": variables}
    return df


# --- Inventory -------------------------------------------------------------


def _mode_counts(letters):
    s = pd.Series(letters)
    s = s[s != ""]
    return {m: int((s == m).sum()) for m in ("D", "A", "R")}


def available_variables(paths):
    """List the variables a float has, with how much data in each mode.

    Parameters
    ----------
    paths : dict or str
        File name -> local path (as returned by `download_core` /
        `download_aux`), or a folder holding a float's files.

    Returns
    -------
    pd.DataFrame
        One row per (file, variable), with columns file, variable,
        counted_per ("profile" or "measurement"), n (with data),
        n_D, n_A, n_R (how many are in delayed / adjusted / real-time
        mode) and n_adjusted_values.
    """
    if isinstance(paths, str):
        folder = paths
        paths = {
            f: os.path.join(folder, f)
            for f in sorted(os.listdir(folder))
            if f.endswith(".nc")
        }
    rows = []
    for name, path in sorted(paths.items()):
        if name.endswith("_Sprof.nc"):
            ds = _open(path)
            names, modes = _profile_parameters(ds)
            for var in profile_variables(ds):
                dm = _profile_data_mode(names, modes, var)
                has = np.isfinite(np.asarray(ds[var].values, float)).any(
                    axis=1
                )
                adj = var + "_ADJUSTED"
                n_adj = (
                    int(np.isfinite(np.asarray(ds[adj].values, float)).sum())
                    if adj in ds
                    else 0
                )
                rows.append(
                    {
                        "file": name,
                        "variable": var,
                        "counted_per": "profile",
                        "n": int(has.sum()),
                        **{
                            f"n_{k}": v
                            for k, v in _mode_counts(dm[has]).items()
                        },
                        "n_adjusted_values": n_adj,
                    }
                )
        elif name.endswith("_Rtraj.nc"):
            ds = _open(path)
            params = list(_strings(ds["TRAJECTORY_PARAMETERS"]))
            all_modes = _strings(ds["TRAJECTORY_PARAMETER_DATA_MODE"])
            for var in trajectory_variables(ds):
                has = np.isfinite(np.asarray(ds[var].values, float))
                dm = all_modes[:, params.index(var)]
                adj = var + "_ADJUSTED"
                n_adj = (
                    int(np.isfinite(np.asarray(ds[adj].values, float)).sum())
                    if adj in ds
                    else 0
                )
                rows.append(
                    {
                        "file": name,
                        "variable": var,
                        "counted_per": "measurement",
                        "n": int(has.sum()),
                        **{
                            f"n_{k}": v
                            for k, v in _mode_counts(dm[has]).items()
                        },
                        "n_adjusted_values": n_adj,
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "file",
            "variable",
            "counted_per",
            "n",
            "n_D",
            "n_A",
            "n_R",
            "n_adjusted_values",
        ],
    )
