"""Download float files from the Argo GDAC by WMO number."""

import os
import re
import shutil
import urllib.error
import urllib.request
from email.utils import parsedate_to_datetime

from .constants import DACS, GDAC_URL
from .meta import __version__


_HEADERS = {"User-Agent": f"argopal/{__version__}"}
_RE_NC_LINK = re.compile(r'href="([^"?/]+\.nc)"')
_TIMEOUT = 60


def _request(url, method="GET", headers=None):
    req = urllib.request.Request(
        url, headers={**_HEADERS, **(headers or {})}, method=method
    )
    return urllib.request.urlopen(req, timeout=_TIMEOUT)


def _url_exists(url):
    try:
        with _request(url, method="HEAD") as r:
            return r.status == 200
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return False
        raise


def float_url(wmo, dac, aux=False, gdac=GDAC_URL):
    """URL of the GDAC folder holding one float's files.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    dac : str
        Data assembly centre, e.g. "coriolis".
    aux : bool, optional
        If True, point to the auxiliary folder (/aux), which holds the PAL
        acoustic data; otherwise the core folder (/dac).  By default False.
    gdac : str, optional
        Root URL of the GDAC, by default the Ifremer mirror.

    Returns
    -------
    str
        The folder URL, ending with "/".
    """
    root = "aux" if aux else "dac"
    return "{}/{}/{}/{}/".format(gdac.rstrip("/"), root, dac, int(wmo))


def find_dac(wmo, gdac=GDAC_URL):
    """Find which data assembly centre serves a float.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    gdac : str, optional
        Root URL of the GDAC, by default the Ifremer mirror.

    Returns
    -------
    str
        The DAC name, e.g. "coriolis".

    Raises
    ------
    FileNotFoundError
        If no DAC has a folder for this WMO number.
    """
    for dac in DACS:
        if _url_exists(float_url(wmo, dac, gdac=gdac)):
            return dac
    raise FileNotFoundError(f"WMO {wmo} not found on {gdac} under any DAC.")


def list_files(wmo, dac=None, aux=False, gdac=GDAC_URL):
    """List the NetCDF files in one float's GDAC folder.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    dac : str, optional
        Data assembly centre; looked up with `find_dac` if not given.
    aux : bool, optional
        List the auxiliary folder instead of the core one, by default False.
    gdac : str, optional
        Root URL of the GDAC, by default the Ifremer mirror.

    Returns
    -------
    list of str
        File names (not full URLs), sorted.  Empty if the folder does not
        exist (e.g. a float without auxiliary data).
    """
    if dac is None:
        dac = find_dac(wmo, gdac=gdac)
    try:
        with _request(float_url(wmo, dac, aux=aux, gdac=gdac)) as r:
            html = r.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise
    return sorted(set(_RE_NC_LINK.findall(html)))


def _remote_info(url):
    with _request(url, method="HEAD") as r:
        size = r.headers.get("Content-Length")
        modified = r.headers.get("Last-Modified")
    size = int(size) if size is not None else None
    if modified is not None:
        modified = parsedate_to_datetime(modified).timestamp()
    return size, modified


def download_file(url, fname, update=True):
    """Download one file, skipping it if the local copy is up to date.

    The file is written to a temporary ".part" file first and renamed when
    complete, so an interrupted download never leaves a truncated file.
    The local modification time is set to the server's Last-Modified time,
    which is what the up-to-date check compares against.

    Parameters
    ----------
    url : str
        URL of the file.
    fname : str
        Local path to save it to.
    update : bool, optional
        If True (default), re-download when the server copy differs in size
        or is newer than the local copy.  If False, never re-download an
        existing file.

    Returns
    -------
    bool
        True if the file was downloaded, False if the local copy was kept.
    """
    if os.path.isfile(fname):
        if not update:
            return False
        size, modified = _remote_info(url)
        same_size = size is None or size == os.path.getsize(fname)
        not_newer = modified is None or modified <= os.path.getmtime(fname)
        if same_size and not_newer:
            return False
    os.makedirs(os.path.dirname(os.path.abspath(fname)), exist_ok=True)
    part = fname + ".part"
    with _request(url) as r, open(part, "wb") as f:
        shutil.copyfileobj(r, f)
        modified = r.headers.get("Last-Modified")
    if modified is not None:
        t = parsedate_to_datetime(modified).timestamp()
        os.utime(part, (t, t))
    os.replace(part, fname)
    return True


def _download_folder(wmo, data_dir, dac, aux, files, update, gdac, verbose):
    if dac is None:
        dac = find_dac(wmo, gdac=gdac)
    available = list_files(wmo, dac=dac, aux=aux, gdac=gdac)
    if files is not None:
        wanted = {f"{int(wmo)}_{f}.nc" for f in files}
        missing = wanted - set(available)
        if missing:
            raise FileNotFoundError(
                "Not on the GDAC for WMO {}: {}".format(
                    wmo, ", ".join(sorted(missing))
                )
            )
        available = sorted(wanted)
    folder = os.path.join(data_dir, str(int(wmo)))
    base = float_url(wmo, dac, aux=aux, gdac=gdac)
    paths = {}
    for name in available:
        fname = os.path.join(folder, name)
        got = download_file(base + name, fname, update=update)
        if verbose:
            print("{} {}".format("downloaded" if got else "up to date", name))
        paths[name] = fname
    return paths


def download_aux(
    wmo,
    data_dir="argopal_data",
    dac=None,
    update=True,
    gdac=GDAC_URL,
    verbose=True,
):
    """Download all auxiliary files of a float, including the PAL acoustics.

    Files are saved as `data_dir/<wmo>/<wmo>_<kind>_aux.nc`, e.g.
    `argopal_data/7902280/7902280_Rtraj_aux.nc`.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    data_dir : str, optional
        Local folder to save into, by default "argopal_data".
    dac : str, optional
        Data assembly centre; looked up with `find_dac` if not given.
    update : bool, optional
        Re-download files that changed on the server, by default True.
    gdac : str, optional
        Root URL of the GDAC, by default the Ifremer mirror.
    verbose : bool, optional
        Print one line per file, by default True.

    Returns
    -------
    dict
        File name -> local path, for every file in the aux folder.

    Raises
    ------
    FileNotFoundError
        If the float has no auxiliary folder.
    """
    paths = _download_folder(
        wmo, data_dir, dac, True, None, update, gdac, verbose
    )
    if not paths:
        raise FileNotFoundError(
            f"WMO {wmo} has no auxiliary (aux) files on the GDAC."
        )
    return paths


def download_core(
    wmo,
    data_dir="argopal_data",
    files=("meta", "Rtraj", "Sprof"),
    dac=None,
    update=True,
    gdac=GDAC_URL,
    verbose=True,
):
    """Download core Argo files of a float (CTD, BGC, positions, metadata).

    Parameters
    ----------
    wmo : int or str
        WMO number of the float.
    data_dir : str, optional
        Local folder to save into, by default "argopal_data".
    files : tuple of str or None, optional
        Which `<wmo>_<kind>.nc` files to get, by default
        ("meta", "Rtraj", "Sprof").  None downloads every file in the folder
        (not the per-cycle profiles/ subfolder).
    dac, update, gdac, verbose
        As for `download_aux`.

    Returns
    -------
    dict
        File name -> local path.
    """
    return _download_folder(
        wmo, data_dir, dac, False, files, update, gdac, verbose
    )
