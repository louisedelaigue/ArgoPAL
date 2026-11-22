"""Load the 2025 sea-trial float (tests/data/float2025) as argopal tables.

The float is not on the GDAC, so argopal's readers cannot fetch it; this
helper builds the same tables they would (`read_pal`, `read_park`,
`read_profiles`, `read_positions`) from the copied test files.
"""

import os

import numpy as np
import pandas as pd

import argopal


DATA = os.path.join(os.path.dirname(__file__), "data", "float2025")
AIS_CSV = os.path.join(DATA, "ais_positions.csv.gz")
AIS_NMEA = os.path.join(DATA, "ais_messages.nmea.gz")
PARK_DEPTHS = {41: 500, 42: 740, 45: 1020}  # approximate, dbar


def read_pal():
    """Like `argopal.read_pal`: one row per park-drift spectrum, 23 bands."""
    df = pd.read_csv(
        os.path.join(DATA, "pal.csv.gz"), parse_dates=["datetime"]
    )
    df.insert(1, "measurement_code", 290)
    df["spectrum_qc"] = "0"
    df["cycle"] = df["cycle"].astype("Int64")
    df["measurement_code"] = df["measurement_code"].astype("Int64")
    df.attrs = {
        "wmo": None,
        "bands_hz": list(argopal.PAL_BANDS_HZ),
        "units": "dB re 1 uPa",
        "band_normalization": False,  # the 2025 firmware sent band levels
        "source_file": "float2025/pal.csv.gz",
    }
    df.attrs["band_order_score"] = argopal.check_band_order(df)
    df.attrs["band_layout_counts"] = argopal.check_band_layout(df)
    return df


def _ctd():
    return pd.read_csv(
        os.path.join(DATA, "ctd.csv.gz"), parse_dates=["datetime"]
    )


def read_park():
    """Like `argopal.read_park(mode="raw")`: park-drift CTD samples."""
    c = _ctd()
    c = c[c["phase"] == "PAR"]
    df = pd.DataFrame(
        {
            "cycle": c["cycle"].astype("Int64").to_numpy(),
            "measurement_code": 290,
            "datetime": c["datetime"].to_numpy(),
            "PRES": c["pressure"].to_numpy(),
            "PRES_qc": "1",
            "PRES_data_mode": "R",
            "TEMP": c["temperature"].to_numpy(),
            "TEMP_qc": "1",
            "TEMP_data_mode": "R",
            "PSAL": c["salinity"].to_numpy(),
            "PSAL_qc": "1",
            "PSAL_data_mode": "R",
        }
    )
    df["measurement_code"] = df["measurement_code"].astype("Int64")
    df = df.sort_values("datetime").reset_index(drop=True)
    df.attrs = {"wmo": None, "mode": "raw", "variables": ["TEMP", "PSAL"]}
    return df


def read_positions():
    """Like `argopal.read_positions`: surface GPS fixes."""
    p = pd.read_csv(
        os.path.join(DATA, "positions.csv"), parse_dates=["datetime"]
    )
    p["position_qc"] = "1"
    p["measurement_code"] = 703.0
    return p[
        [
            "datetime",
            "latitude",
            "longitude",
            "position_qc",
            "cycle",
            "measurement_code",
        ]
    ].reset_index(drop=True)


def read_profiles():
    """Like `argopal.read_profiles(mode="raw")`: ascent T/S profiles, each at
    the GPS fix that ended it."""
    c = _ctd()
    c = c[c["phase"] == "ASC"].dropna(subset=["pressure"])
    fixes = read_positions()
    end = (
        fixes[fixes["cycle"].isin(c["cycle"].unique())].groupby("cycle").last()
    )
    df = pd.DataFrame(
        {
            "cycle": c["cycle"].astype("Int64").to_numpy(),
            "direction": "A",
            "datetime": c["cycle"].map(end["datetime"]).to_numpy(),
            "latitude": c["cycle"].map(end["latitude"]).to_numpy(),
            "longitude": c["cycle"].map(end["longitude"]).to_numpy(),
            "position_qc": "1",
            "PRES": c["pressure"].to_numpy(),
            "PRES_qc": "1",
            "PRES_data_mode": "R",
            "TEMP": c["temperature"].to_numpy(),
            "TEMP_qc": "1",
            "TEMP_data_mode": "R",
            "PSAL": c["salinity"].to_numpy(),
            "PSAL_qc": "1",
            "PSAL_data_mode": "R",
        }
    )
    df = df[np.isfinite(df[["TEMP", "PSAL"]].to_numpy(float)).all(axis=1)]
    df = df.sort_values(["cycle", "PRES"]).reset_index(drop=True)
    df.attrs = {"wmo": None, "mode": "raw", "variables": ["TEMP", "PSAL"]}
    return df


def process(bands=(8000,), ais=None, ais_action="flag"):
    """Run argopal's steps 2-8 on the 2025 float, as `process_float` would
    after downloading."""
    pal_raw = read_pal()
    pal = argopal.clean_pal(argopal.select_bands(pal_raw, bands))
    park = read_park()
    pal = argopal.attach_park(pal, park)
    positions = read_positions()
    pal = argopal.attach_positions(pal, positions)
    ais_positions = None
    if ais is not None:
        ais_positions = argopal.read_ais(ais, bbox=argopal.track_bbox(pal, 20))
        pal = argopal.flag_ships(pal, ais_positions, action=ais_action)
    profiles = read_profiles()
    ts = argopal.gridded_profiles(
        profiles, max_pressure=pal["pressure"].max() + 10
    )
    beta = argopal.compute_profile_betas(ts, pal.attrs["bands_hz"])
    pal = argopal.depth_correct(pal, beta)
    pal = argopal.spectral_level(pal)
    pal = argopal.wind_literature(pal)
    out = {
        "park": pal,
        "pal_raw": pal_raw,
        "park_samples": park,
        "profiles": profiles,
        "ts_profiles": ts,
        "beta": beta,
        "positions": positions,
    }
    if ais_positions is not None:
        out["ais"] = ais_positions
    return out
