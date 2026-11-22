"""Build examples/pal_data_tutorial.ipynb (run, then execute the notebook).

    python examples/make_pal_data_tutorial.py
"""

import os

import nbformat as nbf


HERE = os.path.dirname(os.path.abspath(__file__))
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip()))


md("""
# Getting the PAL data

This notebook gets the acoustic (PAL) data of a BGC-Argo float into two
tables, with no further processing:

- **park**: one row per acoustic measurement during the park drift, with
  every frequency band the float sent, and the float's CTD/BGC data and
  position at that time;
- **profiles**: the float's CTD/BGC profiles.

No depth correction, no wind speed, and no despiking: the band levels are
the ones the float sent, which is what is needed to look at other sounds
than wind, e.g. ships or marine mammals.  The
[ArgoPAL tutorial](argopal_tutorial.ipynb) goes on to wind speed.

Float **7902280** (Ligurian Sea, 2026).
""")

code("""
import warnings

import matplotlib.pyplot as plt

import argopal

warnings.simplefilter("ignore")
""")

md("""
## 1. Download the files

Two sets of files from the Argo GDAC: the **aux** files, which hold the
acoustic spectra and the hydrophone configuration, and the **core** files
(trajectory, profiles, metadata).  They are saved in `argopal_data/7902280/`
and only downloaded again when the GDAC copy is newer.
""")

code("""
files = argopal.download_aux(7902280, verbose=False)
files |= argopal.download_core(7902280, verbose=False)
sorted(files)
""")

md("""
## 2. Read the acoustic spectra

`read_pal` gives one row per acoustic measurement: `cycle`,
`measurement_code`, `datetime`, `pressure`, and one column per band,
`f_63Hz` … `f_25000Hz`.
""")

code("""
meta = argopal.read_pal_meta(files["7902280_meta_aux.nc"])
pal = argopal.read_pal(files["7902280_Rtraj_aux.nc"], meta=meta)
pal.head()
""")

md("""
The PAL sends either 23 bands or a 9-band subset; this float sends mostly
9 bands, so the other columns are nearly empty.  Whether the band levels
were divided by the bandwidth on board is in the hydrophone
configuration: here they were, so the values are spectral levels
(dB re 1 µPa² Hz⁻¹) at the hydrophone's depth.
""")

code("""
print("spectra by number of bands:", pal.attrs["band_layout_counts"])
print("bands with data (Hz):", argopal.available_bands(pal))
print("normalised by bandwidth on board:", meta["band_normalization"])
print("measurement codes:", sorted(pal["measurement_code"].unique().tolist()))
""")

md("""
This float records acoustics during the park drift only (measurement code
290).  A float that also records during descent or ascent would show other
codes here; `clean_pal(..., measurement_codes=...)` keeps them.

## 3. Keep the bands with data, and clean

`select_bands(pal, bands="all")` keeps every band that holds data;
`bands=(1000, 8000)` would keep only those.  `clean_pal` keeps the park
drift and drops spectra with a bad QC flag.

Despiking is **off** here (`despike_spectra=False`): it removes short,
loud events, which is right for wind but removes ships passing or
marine-mammal sounds.
""")

code("""
pal = argopal.select_bands(pal, bands="all")
park = argopal.clean_pal(pal, despike_spectra=False)
print(len(park), "measurements;", park.attrs["clean"])
""")

md("""
One cycle, every band: the measurements are 5 min apart, and the short
peaks are what despiking would remove.
""")

code("""
one = park[park["cycle"] == 10]
fig, ax = plt.subplots(figsize=(12, 4.5), layout="constrained")
for f in park.attrs["bands_hz"]:
    ax.plot(one["datetime"], one[f"f_{f}Hz"], lw=0.8, label=f"{f:g} Hz")
ax.set_ylabel("Band level (dB re 1 µPa)")
ax.legend(ncol=5, fontsize=8, frameon=False)
ax.grid(alpha=0.3)
""")

md("""
## 4. Add the CTD/BGC data and positions

The park-drift CTD/BGC samples are interpolated in time onto each acoustic
measurement (within a cycle, between samples at most 2 h apart), and the
position between the surface fixes before and after.

`mode` chooses which Argo values are used, with no mixing: `"delayed"`
(quality-controlled, the default), `"adjusted"` or `"raw"`.  This float is
recent, so its park data are not in delayed mode yet:
""")

code("""
park_delayed = argopal.read_park(files["7902280_Rtraj.nc"], mode="delayed")
park_raw = argopal.read_park(files["7902280_Rtraj.nc"], mode="raw")
print("park temperature, delayed mode:", park_delayed["TEMP"].notna().sum(), "values")
print("park temperature, raw:", park_raw["TEMP"].notna().sum(), "values")
""")

md("""
So raw values here (each column `X` comes with `X_qc` and `X_data_mode`,
to filter on later):
""")

code("""
park = argopal.attach_park(park, park_raw)
park = argopal.attach_positions(
    park, argopal.read_positions(files["7902280_Rtraj.nc"])
)
park[["datetime", "pressure", "f_1000Hz", "f_8000Hz", "TEMP", "PSAL",
      "DOXY", "latitude", "longitude"]].head()
""")

md("""
## 5. The profiles

One row per profile level, every variable with its QC flag and data mode,
in the same data mode:
""")

code("""
profiles = argopal.read_profiles(files["7902280_Sprof.nc"], mode="raw")
print(profiles["cycle"].nunique(), "profiles")
profiles.head()
""")

md("""
## 6. A first look

Daily median of each band over the deployment, and the median spectrum:
""")

code("""
fig, axes = argopal.plot_pal_spectra(park)
""")

md("""
## 7. Save

Both tables as CSV files, or in one NetCDF file (one group per table, with
the processing record as attributes):
""")

code("""
tables = {"park": park, "profiles": profiles}
argopal.export(tables, "output/wmo7902280", format="csv")
""")

code("""
argopal.export(tables, "output/wmo7902280.nc", format="netcdf")
""")

md("""
Every column is described on the
[Outputs](https://argopal.readthedocs.io/en/latest/outputs.html) page.  To
go on to the depth correction and wind speed, see the
[ArgoPAL tutorial](argopal_tutorial.ipynb), or `argopal.process_float`.
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    },
    "language_info": {"name": "python"},
}
nbf.write(nb, os.path.join(HERE, "pal_data_tutorial.ipynb"))
print("written", len(cells), "cells")
