"""Rebuild the small test files in tests/data from full GDAC downloads.

The fixtures are slices of real files of float 7902280.  They are cut with
the netCDF4 library so every variable keeps its exact on-disk type and
layout (char arrays, fill values, compression), i.e. they read exactly like
files straight from the GDAC.

Usage (needs the full files, e.g. from `ap.process_float(7902280)`):

    python tests/make_fixtures.py argopal_data/7902280
"""

import os
import sys

import netCDF4
import numpy as np


HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def subset(src, dst, dim, index, comment):
    """Copy `src` to `dst`, keeping only `index` along dimension `dim`."""
    index = np.asarray(index)
    with netCDF4.Dataset(src) as a, netCDF4.Dataset(dst, "w") as b:
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        b.setncattr("comment_argopal", comment)
        for name, d in a.dimensions.items():
            size = len(index) if name == dim else len(d)
            b.createDimension(name, None if d.isunlimited() else size)
        for name, v in a.variables.items():
            fill = (
                v.getncattr("_FillValue")
                if "_FillValue" in v.ncattrs()
                else None
            )
            w = b.createVariable(
                name, v.datatype, v.dimensions, zlib=True, fill_value=fill
            )
            w.setncatts(
                {k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"}
            )
            v.set_auto_maskandscale(False)
            w.set_auto_maskandscale(False)
            data = v[:]
            if dim in v.dimensions:
                data = np.take(data, index, axis=v.dimensions.index(dim))
            w[:] = data


def char_rows(var):
    """Char variable (..., STRINGn) -> array of stripped str."""
    raw = np.asarray(var[:])
    return np.char.strip(netCDF4.chartostring(raw).astype(str))


def main(folder):
    wmo = os.path.basename(os.path.normpath(folder))

    # PAL: first 200 spectra (cycle 1 and the start of cycle 2)
    subset(
        os.path.join(folder, f"{wmo}_Rtraj_aux.nc"),
        os.path.join(HERE, f"{wmo}_Rtraj_aux_first200.nc"),
        "N_MEASUREMENT",
        np.arange(200),
        "First 200 N_MEASUREMENT of the GDAC file, for argopal tests",
    )

    # Core trajectory: all measurements of cycles 1-2, plus the launch
    # position (cycle -1) and the cycle-0 GPS fixes before cycle 1
    src = os.path.join(folder, f"{wmo}_Rtraj.nc")
    with netCDF4.Dataset(src) as a:
        cycle = np.asarray(a["CYCLE_NUMBER"][:])
    subset(
        src,
        os.path.join(HERE, f"{wmo}_Rtraj_cycles1-2.nc"),
        "N_MEASUREMENT",
        np.flatnonzero(np.isin(cycle, [-1, 0, 1, 2])),
        "Measurements of cycles -1 to 2 of the GDAC file, for argopal tests",
    )

    # Profiles: the first two delayed-mode and first two adjusted-mode
    src = os.path.join(folder, f"{wmo}_Sprof.nc")
    with netCDF4.Dataset(src) as a:
        names = char_rows(a["STATION_PARAMETERS"])
        modes = np.char.strip(
            np.asarray(a["PARAMETER_DATA_MODE"][:]).astype(str)
        )
    i_temp = list(names[0]).index("TEMP")
    d = np.flatnonzero(modes[:, i_temp] == "D")[:2]
    m = np.flatnonzero(modes[:, i_temp] == "A")[:2]
    keep = np.sort(np.concatenate([d, m]))
    subset(
        src,
        os.path.join(HERE, f"{wmo}_Sprof_4profiles.nc"),
        "N_PROF",
        keep,
        f"Profiles {keep.tolist()} (2 delayed-mode, 2 adjusted-mode) of "
        "the GDAC file, for argopal tests",
    )


if __name__ == "__main__":
    main(sys.argv[1])
