"""Write ArgoPAL results to CSV or NetCDF."""

import json
import os

import numpy as np
import pandas as pd
import xarray as xr

from .meta import __version__


def _attr(value):
    """NetCDF attributes must be str or numbers; encode the rest as JSON."""
    if isinstance(value, (bool, np.bool_)):  # before int: bool is an int
        return int(value)
    if isinstance(value, (str, int, float, np.integer, np.floating)):
        return value
    return json.dumps(value, default=str)


def _to_dataset(df, dim):
    df = df.copy()
    for col in df.columns:
        if isinstance(df[col].dtype, pd.Int64Dtype):
            if df[col].isna().any():
                df[col] = df[col].astype(float)
            else:
                df[col] = df[col].astype("int64")
        elif df[col].dtype == object:
            df[col] = df[col].fillna("").astype(str)
    ds = xr.Dataset.from_dataframe(df.reset_index(drop=True))
    ds = ds.rename({"index": dim}).drop_vars(dim)
    ds.attrs = {k: _attr(v) for k, v in df.attrs.items() if v is not None}
    return ds


def _tables(obj):
    """DataFrame or dict -> {name: DataFrame}; non-table values skipped."""
    if isinstance(obj, pd.DataFrame):
        return {"data": obj}
    if isinstance(obj, dict):
        tables = {k: v for k, v in obj.items() if isinstance(v, pd.DataFrame)}
        if tables:
            return tables
    raise TypeError(
        "export needs a DataFrame or a dict of DataFrames (e.g. the output "
        "of process_float)."
    )


def export(obj, path, format="csv"):
    """Write ArgoPAL tables to CSV or NetCDF, at any stage of processing.

    Works on the output of any ArgoPAL step: a single table (e.g. from
    `read_pal`, `clean_pal`, `read_profiles`, `depth_correct`) or a dict of
    tables (e.g. from `process_float`; entries that are not tables, such as
    "meta", are skipped).

    Parameters
    ----------
    obj : pd.DataFrame or dict
        One table, or name -> table.
    path : str
        For "csv": one table is written to `path` (".csv" added if
        missing); a dict is written to `<path>_<name>.csv` per table.
        For "netcdf": the file name; one table goes in the root group, a
        dict gets one group per table.
    format : str, optional
        "csv" (default) or "netcdf".

    Returns
    -------
    list of str
        The files written.
    """
    if format not in ("csv", "netcdf"):
        raise ValueError(f"format must be 'csv' or 'netcdf'; got {format!r}.")
    single = isinstance(obj, pd.DataFrame)
    tables = _tables(obj)
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)

    if format == "csv":
        if single:
            fname = path if path.endswith(".csv") else path + ".csv"
            obj.to_csv(fname, index=False)
            return [fname]
        base = path.removesuffix(".csv")
        written = []
        for name, df in tables.items():
            fname = f"{base}_{name}.csv"
            df.to_csv(fname, index=False)
            written.append(fname)
        return written

    first = next(iter(tables.values()))
    root_attrs = {
        "title": "ArgoPAL: BGC-Argo float data with PAL acoustics",
        "argopal_version": __version__,
        "wmo": _attr(first.attrs.get("wmo", "")),
        "tables": ",".join(tables) if not single else "",
    }
    if single:
        ds = _to_dataset(obj, "N_OBS")
        ds.attrs = {**root_attrs, **ds.attrs}
        ds.to_netcdf(path, mode="w", engine="netcdf4", encoding=_encoding(ds))
        return [path]
    xr.Dataset(attrs=root_attrs).to_netcdf(path, mode="w", engine="netcdf4")
    for name, df in tables.items():
        ds = _to_dataset(df, "N_OBS")
        ds.to_netcdf(
            path,
            mode="a",
            group=name,
            engine="netcdf4",
            encoding=_encoding(ds),
        )
    return [path]


def _encoding(ds):
    encoding = {}
    for v in ds.data_vars:
        enc = {"zlib": True, "complevel": 4}
        if ds[v].dtype.kind in "OU":  # text: store as char array
            enc["dtype"] = "S1"
        encoding[v] = enc
    return encoding


def _to_frame(ds):
    df = ds.to_dataframe().reset_index(drop=True)
    df.attrs = dict(ds.attrs)
    return df


def read_export(path):
    """Read a NetCDF file written by `export` back into DataFrames.

    Returns
    -------
    pd.DataFrame or dict
        One table if a single table was exported; otherwise name -> table.
    """
    with xr.open_dataset(path) as root:
        names = root.attrs.get("tables", "")
        if not names:
            return _to_frame(root.load())
    out = {}
    for name in names.split(","):
        with xr.open_dataset(path, group=name) as ds:
            out[name] = _to_frame(ds.load())
    return out
