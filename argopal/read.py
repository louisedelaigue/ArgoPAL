"""Read PAL passive acoustic data from Argo auxiliary NetCDF files."""

import os

import numpy as np
import pandas as pd
import xarray as xr

from .checks import check_band_layout, check_band_order
from .constants import (
    CLASSIFICATION_VAR,
    SPECTRUM_VAR,
    band_column,
    bands_for,
)


def _strings(da):
    """Argo char variables (bytes, padded) -> numpy array of clean str."""
    values = np.asarray(da.values)
    out = []
    for v in values.ravel():
        if isinstance(v, bytes):
            v = v.decode("utf-8", errors="replace")
        elif v is None or (isinstance(v, float) and np.isnan(v)):
            v = ""  # fill value of a padded char array
        out.append(str(v).strip())
    return np.array(out).reshape(values.shape)


def _wmo(ds):
    """WMO number from PLATFORM_NUMBER, stored either as one string or as
    an array of single characters, depending on how the file was decoded."""
    chars = _strings(ds["PLATFORM_NUMBER"]).ravel()
    if all(len(c) <= 1 for c in chars):
        return int("".join(chars))
    return int(chars[0])


def _number(value):
    """Config values come as str or float; return int/float where possible."""
    if isinstance(value, (bytes, str)):
        s = value.decode() if isinstance(value, bytes) else value
        s = s.strip()
        try:
            f = float(s)
        except ValueError:
            return s
    else:
        f = float(value)
    if np.isnan(f):
        return np.nan
    return int(f) if f.is_integer() else f


def _open(fname_or_ds):
    if isinstance(fname_or_ds, xr.Dataset):
        return fname_or_ds, None
    return xr.open_dataset(fname_or_ds), str(fname_or_ds)


def _adjusted_or_raw(ds, name):
    """Use NAME_ADJUSTED where it is finite, else NAME.  Return (values,
    source), where source is "adjusted", "raw" or "mixed"."""
    raw = ds[name].values
    adj_name = name + "_ADJUSTED"
    if adj_name not in ds:
        return raw, "raw"
    adj = ds[adj_name].values
    has_adj = ~pd.isna(adj)
    if has_adj.all():
        return adj, "adjusted"
    if not has_adj.any():
        return raw, "raw"
    return np.where(has_adj, adj, raw), "mixed"


def read_pal_meta(fname_or_ds):
    """Read the PAL sensor description and configuration from a meta_aux file.

    Parameters
    ----------
    fname_or_ds : str or xr.Dataset
        Path to `<wmo>_meta_aux.nc`, or the opened dataset.

    Returns
    -------
    dict
        With keys:

        - "wmo" : int
        - "sensor", "sensor_model", "sensor_maker", "sensor_serial_no" :
          str
        - "launch_config" : dict of `CONFIG_Pal*` name -> value at launch
        - "static_config" : dict of `CONFIG_Pal*` name -> value (fixed)
        - "mission_config" : pd.DataFrame of `CONFIG_Pal*` values that
          change between missions, indexed by mission number
        - "band_normalization" : bool or None, from
          CONFIG_PalBandNormalization_LOGICAL (whether the float
          normalises band levels by bandwidth before transmission)
    """
    ds, _ = _open(fname_or_ds)

    def config(names, values, pal_only=True):
        out = {}
        for n, v in zip(_strings(ds[names]), np.asarray(ds[values].values)):
            if n and (n.startswith("CONFIG_Pal") or not pal_only):
                out[n] = _number(v)
        return out

    hydrophone = 0
    sensors = list(_strings(ds["SENSOR"]))
    if "HYDROPHONE" in sensors:
        hydrophone = sensors.index("HYDROPHONE")

    def sensor(var):
        return str(_strings(ds[var])[hydrophone]) if var in ds else ""

    mission_config = pd.DataFrame()
    if "CONFIG_PARAMETER_NAME" in ds:
        names = list(_strings(ds["CONFIG_PARAMETER_NAME"]))
        values = np.asarray(ds["CONFIG_PARAMETER_VALUE"].values, dtype=float)
        mission_config = pd.DataFrame(
            values,
            columns=names,
            index=pd.Index(
                np.asarray(ds["CONFIG_MISSION_NUMBER"].values).astype(int),
                name="mission",
            ),
        )
        mission_config = mission_config[
            [c for c in names if c.startswith("CONFIG_Pal")]
        ]

    launch = config(
        "LAUNCH_CONFIG_PARAMETER_NAME", "LAUNCH_CONFIG_PARAMETER_VALUE"
    )
    norm = launch.get("CONFIG_PalBandNormalization_LOGICAL")
    return {
        "wmo": _wmo(ds),
        "sensor": sensor("SENSOR"),
        "sensor_model": sensor("SENSOR_MODEL"),
        "sensor_maker": sensor("SENSOR_MAKER"),
        "sensor_serial_no": sensor("SENSOR_SERIAL_NO"),
        "launch_config": launch,
        "static_config": config(
            "STATIC_CONFIG_PARAMETER_NAME", "STATIC_CONFIG_PARAMETER_VALUE"
        ),
        "mission_config": mission_config,
        "band_normalization": None if norm is None else bool(norm),
    }


def read_pal(fname_or_ds, meta=None):
    """Read PAL third-octave spectra from an Rtraj_aux file.

    One row per acoustic measurement, one column per frequency band.

    Parameters
    ----------
    fname_or_ds : str or xr.Dataset
        Path to `<wmo>_Rtraj_aux.nc`, or the opened dataset.
    meta : dict, optional
        Output of `read_pal_meta`; if given, its "band_normalization" and
        sensor details are copied into the DataFrame's `attrs`.

    Returns
    -------
    pd.DataFrame
        Columns:

        - cycle : int, float cycle number (adjusted if available)
        - measurement_code : int, Argo reference table 15 (290 = park
          drift)
        - datetime : datetime64, UTC (JULD_ADJUSTED if available)
        - pressure : float, dbar (PRES_ADJUSTED if available)
        - f_63Hz ... f_25000Hz : float, band level in dB re 1 uPa; 23
          columns, or 9 (f_63Hz ... f_20000Hz) for a 9-band float
        - spectrum_qc : str, Argo QC flag of the spectrum
        - variance_class_1 ... variance_class_4 : float, PAL
          classification of the standard deviation of 4 bands (NaN if
          not sent)

        `attrs` holds the WMO number, band frequencies, which time /
        pressure / cycle source was used, and the source file.

    Raises
    ------
    KeyError
        If the file has no PAL spectrum variable.
    ValueError
        If the number of values per spectrum is not 23 or 9 (the only PAL
        layouts, see `argopal.constants.BAND_LAYOUTS`), or does not match
        the number in the dimension name (e.g. N_VALUES23).

    Warns
    -----
    BandOrderWarning
        If the median spectrum does not fall with frequency, which would
        mean the band frequencies are assigned in the wrong order (see
        `check_band_order`).
    """
    ds, fname = _open(fname_or_ds)
    if SPECTRUM_VAR not in ds:
        raise KeyError(
            "{} not in {}; this float may have no PAL data.".format(
                SPECTRUM_VAR, fname or "dataset"
            )
        )
    spectrum = np.asarray(ds[SPECTRUM_VAR].values, dtype=float)
    n_values = spectrum.shape[1]
    # The values dimension is named after its size (N_VALUES23); a mismatch
    # means the file is not laid out the way ArgoPAL expects.
    values_dim = ds[SPECTRUM_VAR].dims[1]
    suffix = values_dim.removeprefix("N_VALUES")
    if suffix.isdigit() and int(suffix) != n_values:
        raise ValueError(
            f"{SPECTRUM_VAR} has {n_values} values per measurement "
            f"but its dimension is named {values_dim}."
        )
    bands = bands_for(n_values)
    keep = np.isfinite(spectrum).any(axis=1)

    datetime, time_source = _adjusted_or_raw(ds, "JULD")
    pressure, pres_source = _adjusted_or_raw(ds, "PRES")
    cycle, cycle_source = _adjusted_or_raw(ds, "CYCLE_NUMBER")

    df = pd.DataFrame(
        {
            "cycle": cycle,
            "measurement_code": ds["MEASUREMENT_CODE"].values,
            "datetime": pd.to_datetime(datetime),
            "pressure": np.asarray(pressure, dtype=float),
        }
    )
    for i, f in enumerate(bands):
        df[band_column(f)] = spectrum[:, i]
    qc_var = SPECTRUM_VAR + "_QC"
    if qc_var in ds:
        df["spectrum_qc"] = _strings(ds[qc_var])
    if CLASSIFICATION_VAR in ds:
        classes = np.asarray(ds[CLASSIFICATION_VAR].values, dtype=float)
        for j in range(classes.shape[1]):
            df[f"variance_class_{j + 1}"] = classes[:, j]

    df = df.loc[keep].reset_index(drop=True)
    df["cycle"] = df["cycle"].astype("Int64")
    df["measurement_code"] = df["measurement_code"].astype("Int64")

    df.attrs = {
        "wmo": _wmo(ds),
        "bands_hz": list(bands),
        "units": "dB re 1 uPa",
        "time_source": time_source,
        "pressure_source": pres_source,
        "cycle_source": cycle_source,
        "source_file": os.path.basename(fname) if fname else None,
    }
    if meta is not None:
        for k in ("band_normalization", "sensor_model", "sensor_serial_no"):
            df.attrs[k] = meta.get(k)
    df.attrs["band_order_score"] = check_band_order(df)
    df.attrs["band_layout_counts"] = check_band_layout(df)
    return df
