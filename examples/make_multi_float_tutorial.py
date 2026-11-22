"""Build examples/multi_float_tutorial.ipynb (run, then execute the notebook).

    python examples/make_multi_float_tutorial.py
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
# Finding PAL floats and processing several

Which floats carry a passive acoustic listener (PAL), where and when did
they measure, and what wind did they see?  This notebook:

1. finds every float on the Argo GDAC with PAL data;
2. maps where they measured;
3. selects floats by region and period;
4. processes each of them;
5. compares them.
""")

code("""
import warnings

import matplotlib.pyplot as plt
import pandas as pd

import argopal

warnings.simplefilter("ignore")
""")

md("""
## 1. Find the floats

No GDAC index says which floats have acoustic data, so
`find_pal_floats` checks every float folder in the GDAC's `aux/`
directory for the PAL spectrum variable.  This takes a couple of minutes
the first time; the list is then saved in `argopal_data/pal_floats.csv`
and reused for a week (`refresh=True` scans again).

Where and when each float measured comes from the GDAC's index of BGC
profiles.
""")

code("""
floats = argopal.find_pal_floats()
floats
""")

md("""
## 2. Where they measured

Every profile of these floats, from the same index:
""")

code("""
index = argopal.read_profile_index(floats["wmo"])

fig, ax = plt.subplots(figsize=(11, 4.5), layout="constrained")
for wmo, g in index.groupby("wmo"):
    ax.plot(g["longitude"], g["latitude"], ".-", ms=3, lw=0.5, label=wmo)
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_aspect("equal")
ax.grid(alpha=0.3)
_ = ax.legend(title="WMO", fontsize=8)
""")

md("""
## 3. Select by region and period

For example the floats with profiles in the western Mediterranean in
2026.  Any condition on the index works (position, date, cycle, ...).
""")

code("""
region = index["longitude"].between(-6, 17) & index["latitude"].between(35, 45)
period = index["datetime"].between("2026-01-01", "2026-12-31")
selected = sorted(int(w) for w in index.loc[region & period, "wmo"].unique())
selected
""")

md("""
## 4. Process each float

`process_float` on each, in a loop that carries on if one float fails
(e.g. a file being updated on the GDAC), with the WMO number added to
each table.

These floats are recent, so most of their data are not in delayed mode
yet: `mode="raw"` uses the real-time values (see the
[tutorial on the full pipeline](argopal_tutorial.ipynb) for the modes).
""")

code("""
results, failed = {}, {}
for wmo in selected:
    try:
        results[wmo] = argopal.process_float(wmo, mode="raw", verbose=False)
    except Exception as e:
        failed[wmo] = f"{type(e).__name__}: {e}"
print("processed:", list(results), "failed:", failed)
""")

md("""
All the park tables in one, with a `wmo` column:
""")

code("""
park = pd.concat(
    [r["park"].assign(wmo=wmo) for wmo, r in results.items()],
    ignore_index=True,
)
park.groupby("wmo").agg(
    measurements=("datetime", "size"),
    first=("datetime", "min"),
    last=("datetime", "max"),
    park_depth_m=("depth", "median"),
).round(0)
""")

md("""
The floats do not all send the same spectra: some send all 23 bands,
others mostly the 9-band subset (`pal_raw` keeps every band):
""")

code("""
pd.DataFrame(
    {wmo: r["pal_raw"].attrs["band_layout_counts"] for wmo, r in results.items()}
).T
""")

md("""
## 5. Compare the floats

Daily median wind speed from the published Nystuen et al. (2015) model.
Without a reference buoy, these winds are not calibrated: the published
models are biased on floats (see the
[full-pipeline tutorial](argopal_tutorial.ipynb), section 7), so this
compares the floats with each other rather than giving absolute winds.
""")

code("""
def daily_median(park):
    return (
        park.set_index("datetime")
        .groupby("wmo")["U_nystuen2015"]
        .resample("1D")
        .median()
        .unstack(0)
    )


def plot_daily(daily, ylim=None):
    fig, ax = plt.subplots(figsize=(12, 4), layout="constrained")
    for wmo in daily:
        ax.plot(daily.index, daily[wmo], lw=1, label=wmo)
    ax.set_ylabel("Wind speed (m s$^{-1}$), daily median")
    ax.set_ylim(ylim)
    ax.grid(alpha=0.3)
    ax.legend(title="WMO", fontsize=8)


plot_daily(daily_median(park))
""")

md("""
Some days reach 80–120 m s⁻¹, which is not wind.  They come from a few
cycles where the float recorded loud, **flat** spectra: wind noise falls
with frequency, these do not.  Float 2904023, cycle 17 against the cycle
before:
""")

code("""
raw = results[2904023]["pal_raw"]
bands = [c for c in raw if c.startswith("f_") and raw[c].notna().mean() > 0.5]
freqs = [int(c[2:-2]) for c in bands]

fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
for cycle in (16, 17):
    ax.plot(freqs, raw.loc[raw["cycle"] == cycle, bands].median(), ".-",
            label=f"cycle {cycle}")
ax.set_xscale("log")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Median band level (dB re 1 µPa)")
ax.grid(alpha=0.3)
_ = ax.legend()
""")

md("""
ArgoPAL does not flag such spectra yet (despiking removes short events,
not a whole cycle).  For the comparison below, measurements above
30 m s⁻¹ are left out:
""")

code("""
clean = park[park["U_nystuen2015"] <= 30]
print(f"left out: {len(park) - len(clean)} of {len(park)} measurements")
plot_daily(daily_median(clean))
""")

md("""
Wind along each track (daily medians at the float's daily position):
""")

code("""
track = (
    clean.set_index("datetime")
    .groupby("wmo")[["longitude", "latitude", "U_nystuen2015"]]
    .resample("1D")
    .median()
    .dropna()
    .reset_index()
)

fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
sc = ax.scatter(
    track["longitude"], track["latitude"], c=track["U_nystuen2015"],
    cmap="viridis", s=14, vmin=0, vmax=15,
)
for wmo, g in track.groupby("wmo"):
    ax.annotate(wmo, (g["longitude"].iloc[0], g["latitude"].iloc[0]),
                fontsize=8, xytext=(4, 4), textcoords="offset points")
fig.colorbar(sc, ax=ax, label="Wind speed (m s$^{-1}$), daily median")
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_aspect("equal")
ax.grid(alpha=0.3)
""")

md("""
Monthly statistics per float:
""")

code("""
(
    clean.assign(month=clean["datetime"].dt.to_period("M"))
    .groupby(["wmo", "month"])["U_nystuen2015"]
    .median()
    .unstack(0)
    .round(1)
    .tail(6)
)
""")

md("""
## 6. Save

One NetCDF file per float, or the combined park table as one CSV file:
""")

code("""
for wmo, r in results.items():
    argopal.export(r, f"output/wmo{wmo}.nc", format="netcdf")
argopal.export(park, "output/park_all_floats.csv")
""")

md("""
With reference winds for the region (a table with `reference_name`,
`datetime`, `latitude`, `longitude`, `wind_speed`), the same loop
calibrates every float: `argopal.process_float(wmo, reference=...)`.
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
nbf.write(nb, os.path.join(HERE, "multi_float_tutorial.ipynb"))
print("written", len(cells), "cells")
