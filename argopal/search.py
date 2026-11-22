"""Find the floats that carry a PAL, and where and when they measured.

No GDAC index records which floats have passive acoustic data, so
`find_pal_floats` looks for it in the aux folders: for each float with a
`<wmo>_Rtraj_aux.nc` file, the first part of the file (where NetCDF-4
keeps the variable names) is read and searched for the PAL spectrum
variable.  A full scan takes a couple of minutes; the result is saved and
reused.  Where and when each float measured comes from the GDAC's
synthetic-profile index (`read_profile_index`).
"""

import concurrent.futures
import os
import re
import time
import urllib.error

import pandas as pd

from .constants import GDAC_URL
from .gdac import _request, download_file


PAL_VARIABLE = b"PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM"
HEAD_BYTES = 262144  # the variable names sit in the first ~30 kB
PROFILE_INDEX = "argo_synthetic-profile_index.txt.gz"


def _get(url, nbytes=None):
    headers = {"Range": f"bytes=0-{nbytes - 1}"} if nbytes else None
    with _request(url, headers=headers) as r:
        return r.read()


def _folders(url):
    """Sub-folder names in a GDAC directory listing."""
    return re.findall(r'href="([\w.-]+)/"', _get(url).decode())


def _has_pal(gdac, dac, wmo):
    """True / False if the float has PAL spectra, None if not checked."""
    folder = f"{gdac}/aux/{dac}/{wmo}/"
    try:
        listing = _get(folder).decode()
        name = f"{wmo}_Rtraj_aux.nc"
        if f'"{name}"' not in listing:
            return False
        return PAL_VARIABLE in _get(folder + name, HEAD_BYTES)
    except (urllib.error.URLError, OSError):
        return None


def find_pal_floats(
    data_dir="argopal_data",
    refresh=False,
    max_age="7D",
    summary=True,
    workers=8,
    gdac=GDAC_URL,
    verbose=True,
):
    """List the floats on the Argo GDAC that have PAL acoustic data.

    Every float folder under `<gdac>/aux/` is checked (a couple of minutes
    for the whole GDAC, with `workers` parallel requests).  The list is
    saved to `<data_dir>/pal_floats.csv` and reused while it is younger
    than `max_age`.

    Parameters
    ----------
    data_dir : str, optional
        Folder for the saved list and the profile index, by default
        "argopal_data".
    refresh : bool, optional
        Scan the GDAC again even if a recent list is saved, by default
        False.
    max_age : str or pd.Timedelta, optional
        Age after which a saved list is scanned again, by default "7D".
    summary : bool, optional
        Add where and when each float measured, from the synthetic-profile
        index (`read_profile_index`), by default True.
    workers : int, optional
        Parallel requests to the GDAC, by default 8.
    gdac : str, optional
        GDAC root URL.
    verbose : bool, optional
        Print progress.

    Returns
    -------
    pd.DataFrame
        One row per float: "wmo", "dac", and with `summary` also
        "n_profiles", "first_date", "last_date" and the first and last
        positions ("first_latitude", "first_longitude", "last_latitude",
        "last_longitude").  `attrs["scanned"]` is when the GDAC was
        scanned, `attrs["n_checked"]` how many floats were checked.
    """
    path = os.path.join(data_dir, "pal_floats.csv")
    floats = None
    if os.path.isfile(path) and not refresh:
        age = time.time() - os.path.getmtime(path)
        if age < pd.Timedelta(max_age).total_seconds():
            floats = pd.read_csv(path, dtype={"wmo": int, "dac": str})
            scanned = pd.Timestamp(os.path.getmtime(path), unit="s")
            floats.attrs = {"scanned": str(scanned.floor("s"))}
            if verbose:
                print(f"PAL floats from {path} (scanned {scanned:%Y-%m-%d}).")
    if floats is None:
        floats = _scan(gdac, workers, verbose)
        os.makedirs(data_dir, exist_ok=True)
        floats[["wmo", "dac"]].to_csv(path, index=False)
    if summary:
        index = read_profile_index(data_dir=data_dir, gdac=gdac)
        floats = _summarise(floats, index)
    return floats


def _scan(gdac, workers, verbose):
    t0 = time.time()
    folders = [
        (dac, wmo)
        for dac in _folders(f"{gdac}/aux/")
        for wmo in _folders(f"{gdac}/aux/{dac}/")
        if wmo.isdigit()
    ]
    if verbose:
        print(f"Checking {len(folders)} floats in {gdac}/aux/ ...")
    with concurrent.futures.ThreadPoolExecutor(workers) as ex:
        found = list(ex.map(lambda f: _has_pal(gdac, *f), folders))
    failed = sum(f is None for f in found)
    rows = [
        {"wmo": int(wmo), "dac": dac}
        for (dac, wmo), pal in zip(folders, found)
        if pal
    ]
    out = pd.DataFrame(rows, columns=["wmo", "dac"])
    out = out.sort_values("wmo").reset_index(drop=True)
    out.attrs = {
        "scanned": str(pd.Timestamp.now().floor("s")),
        "n_checked": len(folders),
        "n_failed": failed,
    }
    if verbose:
        print(
            f"{len(out)} floats with PAL data ({time.time() - t0:.0f} s"
            + (f"; {failed} folders could not be read" if failed else "")
            + ")."
        )
    return out


def read_profile_index(
    wmos=None, data_dir="argopal_data", update=True, gdac=GDAC_URL
):
    """The GDAC's synthetic-profile index: where and when each BGC profile
    was measured.

    The index (`argo_synthetic-profile_index.txt.gz`, ~8 MB) is
    downloaded to `data_dir`, and again only when the GDAC copy is newer.

    Parameters
    ----------
    wmos : int or list of int, optional
        Keep only these floats.
    data_dir : str, optional
        Where to save the index, by default "argopal_data".
    update : bool, optional
        Check the GDAC for a newer index, by default True.
    gdac : str, optional
        GDAC root URL.

    Returns
    -------
    pd.DataFrame
        One row per profile: "wmo", "cycle", "direction" ("A" ascending,
        "D" descending), "dac", "datetime" (UTC), "latitude",
        "longitude", "ocean", "profiler_type", "file".
    """
    fname = os.path.join(data_dir, PROFILE_INDEX)
    download_file(f"{gdac}/{PROFILE_INDEX}", fname, update=update)
    df = pd.read_csv(
        fname,
        comment="#",
        usecols=["file", "date", "latitude", "longitude", "ocean",
                 "profiler_type"],
        dtype={"file": str, "date": str, "ocean": str},
    )  # fmt: skip
    parts = df["file"].str.extract(
        r"^(?P<dac>[^/]+)/(?P<wmo>\d+)/profiles/S[RD]?\d+_(?P<cycle>\d+)"
        r"(?P<direction>D?)\.nc$"
    )
    df["dac"] = parts["dac"]
    df["wmo"] = pd.to_numeric(parts["wmo"], errors="coerce").astype("Int64")
    df["cycle"] = pd.to_numeric(parts["cycle"], errors="coerce").astype(
        "Int64"
    )
    df["direction"] = parts["direction"].replace("", "A")
    df["datetime"] = pd.to_datetime(
        df["date"], format="%Y%m%d%H%M%S", errors="coerce"
    )
    df = df.drop(columns="date")
    if wmos is not None:
        wmos = [wmos] if pd.api.types.is_scalar(wmos) else list(wmos)
        df = df[df["wmo"].isin([int(w) for w in wmos])]
    cols = ["wmo", "cycle", "direction", "dac", "datetime", "latitude",
            "longitude", "ocean", "profiler_type", "file"]  # fmt: skip
    return df[cols].sort_values(["wmo", "datetime"]).reset_index(drop=True)


def _summarise(floats, index):
    idx = index[index["wmo"].isin(floats["wmo"])].dropna(subset=["datetime"])
    g = idx.sort_values("datetime").groupby("wmo")
    summ = pd.DataFrame(
        {
            "n_profiles": g.size(),
            "first_date": g["datetime"].first(),
            "last_date": g["datetime"].last(),
            "first_latitude": g["latitude"].first(),
            "first_longitude": g["longitude"].first(),
            "last_latitude": g["latitude"].last(),
            "last_longitude": g["longitude"].last(),
        }
    )
    out = floats.merge(summ, left_on="wmo", right_index=True, how="left")
    out.attrs = dict(floats.attrs)
    return out
