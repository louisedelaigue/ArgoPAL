"""One-call download-and-read of a float's acoustic data."""

from .constants import GDAC_URL
from .gdac import download_aux
from .read import read_pal, read_pal_meta


def fetch_pal(
    wmo,
    data_dir="argopal_data",
    dac=None,
    update=True,
    gdac=GDAC_URL,
    verbose=True,
):
    """Download (or update) a float's aux files and read its PAL spectra.

    Parameters
    ----------
    wmo : int or str
        WMO number of the float, e.g. 7902280.
    data_dir : str, optional
        Local folder for the downloads, by default "argopal_data".
    dac : str, optional
        Data assembly centre; looked up automatically if not given.
    update : bool, optional
        Re-download files that changed on the server, by default True.
    gdac : str, optional
        Root URL of the GDAC, by default the Ifremer mirror.
    verbose : bool, optional
        Print one line per file, by default True.

    Returns
    -------
    pal : pd.DataFrame
        The spectra, as returned by `read_pal`.
    meta : dict
        The PAL sensor and configuration, as returned by `read_pal_meta`.
    """
    paths = download_aux(
        wmo,
        data_dir=data_dir,
        dac=dac,
        update=update,
        gdac=gdac,
        verbose=verbose,
    )
    wmo = int(wmo)
    traj = paths.get(f"{wmo}_Rtraj_aux.nc")
    if traj is None:
        raise FileNotFoundError(
            "WMO {} has aux files but no Rtraj_aux.nc: {}".format(
                wmo, ", ".join(paths)
            )
        )
    meta_path = paths.get(f"{wmo}_meta_aux.nc")
    meta = read_pal_meta(meta_path) if meta_path else None
    return read_pal(traj, meta=meta), meta
