"""Figures to follow each processing step.

One function per step; each draws one or two panels and returns
``(fig, axes)``, so figures can be shown, adjusted or saved:

    fig, ax = argopal.plot_despike(result["pal_raw"], result["park"])
    argopal.save_figure(fig, "despike.png")        # 300 dpi by default

    argopal.plot_all(result, save_dir="figures")   # every step, saved at 300 dpi

Style follows Delaigue et al. (2026, Ocean Science): white background,
framed axes, light grid, black lines (grey = before a step, dashed = a
second series), panel labels a), b), framed legends.  Many profiles at once
are coloured by date (viridis).  Pass ``annotate=True`` to add a title and
notes (float, data mode, what was skipped and why).
"""

import os
from contextlib import contextmanager

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from .clean import band_columns


DPI = 300  # resolution used when saving figures

AFTER = "black"  # data after a step
BEFORE = "#9e9e9e"  # data before a step
TIME_CMAP = "viridis"  # many profiles at once, early -> late
LEVEL_CMAP = "viridis"  # band level, quiet -> loud

_RC = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "axes.labelcolor": "black",
    "axes.labelsize": 11,
    "axes.grid": True,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "grid.color": "#b0b0b0",
    "grid.alpha": 0.3,
    "grid.linewidth": 0.8,
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "legend.frameon": True,
    "legend.framealpha": 1.0,
    "legend.edgecolor": "#cccccc",
    "lines.linewidth": 1.5,
    "font.size": 10,
}


def save_figure(fig, path, dpi=DPI):
    """Save a figure (300 dpi by default) and return the path."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


@contextmanager
def _style():
    # Only the style keys are saved and restored: plt.rc_context restores
    # every rcParam, which stops Jupyter's inline backend from showing
    # figures made after the first one.
    saved = {k: mpl.rcParams[k] for k in _RC}
    mpl.rcParams.update(_RC)
    try:
        yield
    finally:
        mpl.rcParams.update(saved)


def _panel_labels(axes):
    """a), b), ... above the top-left corner of each panel."""
    for k, ax in enumerate(axes):
        ax.text(
            0, 1.02, f"{chr(97 + k)})", transform=ax.transAxes,
            fontsize=12, va="bottom", ha="left",
        )  # fmt: skip


def _date_axis(ax):
    loc = mdates.AutoDateLocator(minticks=3, maxticks=7)
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))


def _annotate(fig, title=None, note=None):
    """Optional title above and note below the figure."""
    if title:
        fig.suptitle(title, x=0.01, ha="left", fontsize=11)
    if note:
        fig.get_layout_engine().set(rect=(0, 0.06, 1, 0.94))
        fig.text(
            0.01, 0.01, note, fontsize=9, ha="left", va="bottom", wrap=True
        )


def _hours(td):
    try:
        return f"{pd.Timedelta(td).total_seconds() / 3600:g} h"
    except (TypeError, ValueError):
        return str(td)


def _days(td):
    try:
        return f"{pd.Timedelta(td).total_seconds() / 86400:g} days"
    except (TypeError, ValueError):
        return "range"


def _khz(band):
    return f"{band / 1000:g} kHz"


def _default_band(df, band):
    """`band`, or the 8 kHz band if processed, else the first band."""
    if band is not None:
        return band
    bands = df.attrs.get("bands_hz") or [8000]
    return 8000 if 8000 in bands else bands[0]


def _band(pal, band):
    col = f"f_{int(band)}Hz"
    if col not in pal:
        raise KeyError(
            f"No {band} Hz band; available: "
            f"{', '.join(str(f) for f in pal.attrs.get('bands_hz', []))}."
        )
    return col


def _pick_cycle(df, cycle):
    if cycle is not None:
        return cycle
    return int(df["cycle"].value_counts().index[0])  # most measurements


def _wmo(df):
    w = df.attrs.get("wmo")
    return f"WMO {w}: " if w else ""


def _depth_axis(ax, bottom):
    ax.set_ylim(bottom, 0)


def _time_colors(times):
    t = mdates.date2num(pd.to_datetime(times))
    norm = Normalize(t.min(), t.max() if t.max() > t.min() else t.min() + 1)
    cmap = plt.get_cmap(TIME_CMAP)
    return norm, [cmap(norm(x)) for x in t]


def _time_colorbar(fig, axes, norm):
    sm = plt.cm.ScalarMappable(norm=norm, cmap=TIME_CMAP)
    cb = fig.colorbar(sm, ax=axes, pad=0.02, aspect=30)
    loc = mdates.AutoDateLocator(minticks=3, maxticks=6)
    cb.ax.yaxis.set_major_locator(loc)
    cb.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    cb.set_label("Profile date")
    return cb


# --- Step 2: read_pal / clean_pal ------------------------------------------


def plot_pal_spectra(pal, annotate=False):
    """What the float sent: (a) daily median band levels over the
    deployment; (b) the median spectrum, which should fall with frequency.

    Parameters
    ----------
    pal : pd.DataFrame
        Output of `read_pal` or `clean_pal`.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, (ax_time, ax_spectrum)
    """
    cols = band_columns(pal)
    bands = np.asarray(pal.attrs["bands_hz"], dtype=float)
    units = pal.attrs.get("units", "dB").replace("uPa", "µPa")
    df = pal.sort_values("datetime")
    daily = df.set_index("datetime")[cols].resample("1D").median()
    with _style():
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(14, 5), width_ratios=[2, 1], layout="constrained"
        )
        edges = mdates.date2num(
            daily.index.append(daily.index[-1:] + pd.Timedelta("1D"))
        )
        mesh = ax1.pcolormesh(
            edges,
            np.arange(len(bands) + 1),
            daily.to_numpy(float).T,
            cmap=LEVEL_CMAP,
            shading="flat",
            rasterized=True,
        )
        ax1.set_yticks(np.arange(len(bands)) + 0.5)
        ax1.set_yticklabels([f"{b / 1000:g}" for b in bands])
        ax1.set_ylabel("Band centre frequency (kHz)")
        ax1.grid(False)
        _date_axis(ax1)
        cb = fig.colorbar(mesh, ax=ax1, pad=0.02)
        cb.set_label(f"Band level ({units})")

        with np.errstate(all="ignore"):
            q = np.nanpercentile(
                df[cols].to_numpy(float), [10, 50, 90], axis=0
            )
        ax2.fill_between(
            bands, q[0], q[2], color=BEFORE, alpha=0.3, lw=0,
            label="10–90th percentile",
        )  # fmt: skip
        ax2.plot(bands, q[1], color=AFTER, marker="o", ms=4, label="Median")
        ax2.set_xscale("log")
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel(f"Band level ({units})")
        ax2.legend(loc="upper right")
        _panel_labels([ax1, ax2])
        if annotate:
            layout = pal.attrs.get("band_layout_counts") or {}
            score = pal.attrs.get("band_order_score", np.nan)
            _annotate(
                fig,
                f"{_wmo(pal)}{len(pal):,} PAL measurements",
                f"Spectra with all bands: {layout.get('23', 0):,}; 9-band "
                f"spectra (other bands empty): {layout.get('9', 0):,}; other "
                f"patterns: {layout.get('other', 0)}. Band-order check: "
                f"Spearman r = {score:.2f} (expected < −0.8).",
            )
    return fig, (ax1, ax2)


def plot_despike(raw, clean, band=None, cycle=None, annotate=False):
    """Before (grey) and after (black) despiking, for one band: (a) one
    cycle in detail; (b) the whole deployment as the daily 95th
    percentile, where spikes show.

    Parameters
    ----------
    raw : pd.DataFrame
        PAL table before cleaning (`read_pal`; `result["pal_raw"]`).
    clean : pd.DataFrame
        PAL table after `clean_pal` (or any later step).
    band : float, optional
        Band centre frequency (Hz); by default 8 kHz if processed, else
        the first processed band.
    cycle : int, optional
        Cycle shown in detail; by default the one with most measurements.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, (ax_cycle, ax_deployment)
    """
    band = _default_band(clean, band)
    col = _band(clean, band)
    cycle = _pick_cycle(clean, cycle)
    with _style():
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(12, 8), layout="constrained"
        )
        r = raw[raw["cycle"] == cycle].sort_values("datetime")
        c = clean[clean["cycle"] == cycle].sort_values("datetime")
        ax1.plot(r["datetime"], r[col], color=BEFORE, lw=1, label="Before")
        ax1.plot(
            c["datetime"], c[col], color=AFTER, lw=1.5, label="After despiking"
        )
        ax1.set_ylabel(f"{_khz(band)} level, cycle {cycle} (dB)")
        ax1.legend(loc="upper right")
        _date_axis(ax1)

        p95_r = raw.set_index("datetime")[col].resample("1D").quantile(0.95)
        p95_c = clean.set_index("datetime")[col].resample("1D").quantile(0.95)
        ax2.plot(p95_r.index, p95_r, color=BEFORE, lw=1, label="Before")
        ax2.plot(
            p95_c.index, p95_c, color=AFTER, lw=1.5, label="After despiking"
        )
        ax2.set_ylabel(f"{_khz(band)} level, daily 95th percentile (dB)")
        ax2.legend(loc="upper right")
        _date_axis(ax2)
        _panel_labels([ax1, ax2])
        if annotate:
            info = clean.attrs.get("clean", {})
            _annotate(
                fig,
                f"{_wmo(clean)}despiking "
                f"({clean.attrs.get('despike', 'not applied')})",
                f"Rows removed: {info.get('removed_measurement_code', 0)} "
                f"(measurement code), {info.get('removed_bad_qc', 0)} (bad "
                f"QC); {len(raw):,} → {len(clean):,} measurements.",
            )
    return fig, (ax1, ax2)


# --- Step 3: attach_park -------------------------------------------------------


def plot_park(
    park_table, park_samples, variables=None, cycle=None, annotate=False
):
    """Park CTD/BGC samples (open circles) and the values interpolated onto
    the PAL measurements (line), for one cycle; one panel per variable.

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `attach_park` (`result["park"]`).
    park_samples : pd.DataFrame
        Output of `read_park` (`result["park_samples"]`).
    variables : list of str, optional
        Variables to show; by default up to 3 core variables (TEMP, PSAL,
        DOXY, ...) with data in this cycle.
    cycle : int, optional
        Cycle shown; by default the one with most measurements.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, axes
    """
    cycle = _pick_cycle(park_table, cycle)
    p = park_table[park_table["cycle"] == cycle].sort_values("datetime")
    s = park_samples[park_samples["cycle"] == cycle]
    info = park_table.attrs.get("park", {})
    mode = info.get("mode")
    if variables is None:
        candidates = info.get("variables", [])
        core = ["TEMP", "PSAL", "DOXY", "PH_IN_SITU_TOTAL", "NITRATE", "CHLA"]
        candidates = [v for v in core if v in candidates] + [
            v for v in candidates if v not in core
        ]
        variables = [v for v in candidates if s[v].notna().any()][:3]
    with _style():
        if not variables:
            fig, ax = plt.subplots(figsize=(12, 2.5), layout="constrained")
            ax.text(
                0.5, 0.5,
                f"No park CTD/BGC values in data mode {mode!r} for cycle "
                f"{cycle}.",
                ha="center", va="center", transform=ax.transAxes,
            )  # fmt: skip
            ax.set_axis_off()
            return fig, [ax]
        fig, axes = plt.subplots(
            len(variables), 1, figsize=(12, 2.8 * len(variables)),
            sharex=True, layout="constrained", squeeze=False,
        )  # fmt: skip
        axes = list(axes[:, 0])
        for ax, var in zip(axes, variables):
            ax.plot(
                p["datetime"], p[var], color=AFTER, lw=1.5,
                label="Interpolated onto PAL measurements",
            )  # fmt: skip
            sv = s.dropna(subset=[var])
            ax.plot(
                sv["datetime"], sv[var], ls="none", marker="o", ms=5,
                mfc="white", mec=AFTER, mew=1, label="Park samples",
            )  # fmt: skip
            ax.set_ylabel(var)
        axes[0].legend(loc="upper right")
        _date_axis(axes[-1])
        _panel_labels(axes)
        if annotate:
            bad = info.get("bad_qc", [])
            skipped = [
                v
                for v in variables
                if s[v].notna().any() and p[v].isna().all()
            ]
            note = None
            if skipped:
                note = (
                    f"Not interpolated: {', '.join(skipped)} (samples flagged "
                    f"QC {'/'.join(bad)}; raw DOXY is flagged 3 by Argo "
                    "convention)."
                )
            _annotate(
                fig,
                f"{_wmo(park_table)}park data, cycle {cycle} (mode {mode!r}, "
                f"linear in time, gaps ≤ {_hours(info.get('max_gap'))})",
                note,
            )
    return fig, axes


# --- Step 4: attach_positions ----------------------------------------------------


_GRID_STEPS = (0.5, 1, 2, 5, 10, 20, 30, 60)


def _lon_label(lon):
    if np.isclose(lon, 0):
        return "0°"
    return f"{abs(lon):g}°W" if lon < 0 else f"{lon:g}°E"


def map_axes(
    longitude=None,
    latitude=None,
    ax=None,
    global_map=None,
    pad=0.15,
    figsize=None,
):
    """Map axes in ArgoPAL's style: land from Natural Earth (1:50m, light
    grey), dashed labelled gridlines (°W / °E, °N / °S).

    A global map (Robinson projection, gridlines every 60° and 30°) when
    the positions span more than 120° of longitude or 60° of latitude;
    otherwise a regional map (plate carrée) around them, with gridlines
    at a spacing that suits the extent.  Plot on it with
    `transform=cartopy.crs.PlateCarree()`.

    Parameters
    ----------
    longitude, latitude : array_like, optional
        Positions to show (°E, °N); without them, a global map.
    ax : cartopy GeoAxes, optional
        Axes to draw on (made with the projection of your choice);
        by default a new figure.
    global_map : bool, optional
        Force a global (True) or regional (False) map.
    pad : float, optional
        Margin around the positions, as a fraction of their extent, by
        default 0.15 (at least 0.5°).
    figsize : tuple, optional
        Figure size, by default (11, 6) global, (8, 7) regional.

    Returns
    -------
    fig, ax
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import LatitudeFormatter

    lon = np.asarray(longitude if longitude is not None else [], float)
    lat = np.asarray(latitude if latitude is not None else [], float)
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    if global_map is None:
        global_map = len(lon) == 0 or np.ptp(lon) > 120 or np.ptp(lat) > 60
    if ax is None:
        proj = ccrs.Robinson() if global_map else ccrs.PlateCarree()
        size = figsize or ((11, 6) if global_map else (8, 7))
        fig = plt.figure(figsize=size, layout="constrained")
        ax = fig.add_subplot(1, 1, 1, projection=proj)
    else:
        fig = ax.figure
    if global_map:
        ax.set_global()
        xstep, ystep = 60, 30
    else:
        span = max(np.ptp(lon), np.ptp(lat), 1.0)
        m = max(pad * span, 0.5)
        x0, x1 = lon.min() - m, lon.max() + m
        y0, y1 = max(lat.min() - m, -90), min(lat.max() + m, 90)
        ax.set_extent([x0, x1, y0, y1], crs=ccrs.PlateCarree())
        step = next(
            (g for g in _GRID_STEPS if max(x1 - x0, y1 - y0) / g <= 6),
            60,
        )
        xstep = ystep = step
    ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical", "land", "50m", edgecolor="none",
            facecolor="lightgrey",
        ),
        zorder=0,
    )  # fmt: skip
    gl = ax.gridlines(
        draw_labels=True, alpha=0.5, linestyle="--", x_inline=False,
        y_inline=False,
    )  # fmt: skip
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180.1, xstep))
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 90.1, ystep))
    gl.xformatter = mticker.FuncFormatter(lambda x, pos: _lon_label(x))
    gl.yformatter = LatitudeFormatter()
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 8}
    return fig, ax


def plot_positions(park_table, positions, annotate=False):
    """Float track on a map: surface fixes (open circles) and the
    positions interpolated onto the PAL measurements (line).

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `attach_positions`.
    positions : pd.DataFrame
        Output of `read_positions` (`result["positions"]`).
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, ax
    """
    import cartopy.crs as ccrs

    pc = ccrs.PlateCarree()
    with _style():
        p = park_table.dropna(subset=["latitude", "longitude"])
        fig, ax = map_axes(
            np.concatenate([p["longitude"], positions["longitude"]]),
            np.concatenate([p["latitude"], positions["latitude"]]),
        )
        ax.plot(
            p["longitude"], p["latitude"], color=AFTER, lw=1.5,
            label="Interpolated at PAL measurements", transform=pc,
            zorder=3,
        )  # fmt: skip
        ax.plot(
            positions["longitude"], positions["latitude"], ls="none",
            marker="o", ms=5, mfc="white", mec=AFTER, mew=1,
            label="Surface fixes", transform=pc, zorder=4,
        )  # fmt: skip
        first, last = positions.iloc[0], positions.iloc[-1]
        for row, text, dy in ((first, "START", 6), (last, "LATEST", -14)):
            ax.annotate(
                text, (row["longitude"], row["latitude"]), xytext=(6, dy),
                xycoords=pc._as_mpl_transform(ax),
                textcoords="offset points", fontsize=9, fontweight="bold",
                zorder=5,
            )  # fmt: skip
        ax.legend(loc="best")
        if annotate:
            gap = park_table.get("position_gap_h")
            note = None
            if gap is not None and gap.notna().any():
                note = (
                    f"{len(positions)} fixes; time between the two fixes "
                    f"around each PAL measurement: median {gap.median():.0f} "
                    f"h, max {gap.max():.0f} h."
                )
            _annotate(fig, f"{_wmo(park_table)}float track", note)
    return fig, ax


def _has_float_wind(park_table, p):
    """True where the float has a wind (or spectral level) to compare."""
    applied = park_table.attrs.get("wind_models", {}).get("applied", {})
    winds = [f"U_{m}" for m in applied if f"U_{m}" in p]
    if not winds:
        winds = [c for c in p if c.startswith("spl_")]
    if not winds:
        return pd.Series(True, index=p.index)
    return p[winds].notna().any(axis=1)


def _near(park_table, reference, radius_km, min_r=None):
    """Measurements close enough to the reference value matched to them:
    within `radius_km`, or inside its footprint (R >= `min_r`) if given."""
    p = park_table.dropna(subset=["latitude", "longitude"])
    if min_r is not None:
        r = p.get(f"{reference}_footprint_r")
        if r is None:
            return p.iloc[:0]
        ok = r >= min_r
    else:
        dist = p.get(f"{reference}_distance_km")
        if dist is None:
            return p.iloc[:0]
        ok = dist <= radius_km
    return p[ok & p[f"U_{reference}"].notna()]


def _used_measurements(park_table, reference, radius_km, min_r=None):
    """Measurements compared with the reference when refitting: close
    enough to the reference value matched to them (`_near`), and with a
    float wind (or spectral level) to compare."""
    near = _near(park_table, reference, radius_km, min_r)
    return near[_has_float_wind(park_table, near)]


def _station_positions(used, reference, name, fallback, merge_km=2.0):
    """Positions of a reference for the values actually used: one if
    they are all within `merge_km`, else one per group of nearby
    positions (a buoy reported at several positions)."""
    lat_c, lon_c = f"{reference}_latitude", f"{reference}_longitude"
    if lat_c not in used or used.empty:
        return [(fallback["latitude"], fallback["longitude"])]
    pos = used[[lat_c, lon_c]].dropna().round(3)
    counts = pos.value_counts()
    groups = []  # [lat, lon, n]
    for (la, lo), n in counts.items():
        for g in groups:
            if _km(la, lo, g[0], g[1]) <= merge_km:
                g[2] += n
                break
        else:
            groups.append([la, lo, n])
    groups.sort(key=lambda g: -g[2])
    return [(g[0], g[1]) for g in groups]


def _km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (
        np.sin((p2 - p1) / 2) ** 2
        + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    )
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))


def plot_reference_map(
    park_table,
    reference=None,
    radius_km=None,
    footprints=None,
    annotate=False,
):
    """Map of the float's trajectory and of the references used to refit
    the wind models.

    The trajectory, black where the float has a wind and grey where it
    has none (e.g. no delayed-mode profiles yet, so no depth
    correction), with the measurements used to refit the wind models
    (red): within `radius_km` of the reference value matched to them,
    and with a float wind.  Each reference used is drawn (red triangle) at
    the positions of the values used, with a circle of that radius around
    each (several for a buoy reported at several positions).  References
    the float came within `radius_km` of only without a float wind are
    greyed out; those it never came close to are not shown.

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `attach_positions` and `collocate_reference`
        (`result["park"]` with a `reference`).
    reference : str, optional
        Name given to `collocate_reference`; by default the one used for
        fitting, or the only one.
    radius_km : float, optional
        Radius of the circles and of the measurements shown as used; by
        default the one used by `fit_wind_models`, else 40 km.  Not used
        when the refit selected by footprint.
    footprints : xr.DataArray, optional
        Wind footprints (`wind_footprints`, or `result["footprints"]`).
        When the refit selected by footprint, each reference's footprint
        (R at least the threshold used) is outlined instead of a circle.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, ax
    """
    import cartopy.crs as ccrs
    from cartopy.geodesic import Geodesic

    reference = reference or _fit_reference(park_table)
    info = park_table.attrs.get("references", {}).get(reference, {})
    stations = {
        k: v for k, v in info.get("positions", {}).items()
        if v.get("latitude") is not None
    }  # fmt: skip
    if not stations:
        raise KeyError(
            "No reference positions: run collocate_reference first."
        )
    if radius_km is None:
        radius_km = park_table.attrs.get("wind_fits", {}).get(
            "max_distance_km"
        )
        radius_km = 40.0 if radius_km is None else radius_km
    pc = ccrs.PlateCarree()
    p = park_table.dropna(subset=["latitude", "longitude"])
    wind = _has_float_wind(park_table, p)
    min_r = park_table.attrs.get("wind_fits", {}).get("min_footprint_r")
    used = _used_measurements(park_table, reference, radius_km, min_r)
    near = _near(park_table, reference, radius_km, min_r)
    within = (
        f"inside footprint, R ≥ {min_r:g}" if min_r is not None
        else f"within {radius_km:g} km"
    )  # fmt: skip
    name_col = f"{reference}_reference_name"

    def of(table, name):
        if name_col in table:
            return table[table[name_col].astype(str) == name]
        return table  # a single reference

    # references used for fitting, and references the float came close to
    # only without a float wind (e.g. before the delayed-mode data exist)
    shown, unused = {}, {}
    for name, v in stations.items():
        mine = of(used, name)
        if len(mine):
            shown[name] = _station_positions(mine, reference, name, v)
        elif len(of(near, name)):
            unused[name] = _station_positions(
                of(near, name), reference, name, v
            )
    every = list(shown.values()) + list(unused.values())
    slat = [la for pos in every for la, _ in pos]
    slon = [lo for pos in every for _, lo in pos]
    mode = park_table.attrs.get("depth_correction", {}).get("mode")
    with _style():
        fig, ax = map_axes(
            np.concatenate([p["longitude"], slon]),
            np.concatenate([p["latitude"], slat]),
        )
        # trajectory: black with a float wind, grey without
        lon, lat = p["longitude"].to_numpy(), p["latitude"].to_numpy()
        has = wind.to_numpy()
        ax.plot(
            np.where(has, lon, np.nan), np.where(has, lat, np.nan),
            color=AFTER, lw=1.2, transform=pc, zorder=2,
            label="Trajectory with float wind"
            + (f" ({mode} mode)" if mode else ""),
        )  # fmt: skip
        if (~has).any():
            # draw the grey part with its neighbours so that it joins up
            edge = ~has | np.r_[False, ~has[:-1]] | np.r_[~has[1:], False]
            ax.plot(
                np.where(edge, lon, np.nan), np.where(edge, lat, np.nan),
                color=BEFORE, lw=1.2, transform=pc, zorder=1.9,
                label="Trajectory without float wind"
                + (f" (no {mode}-mode profiles)" if mode else ""),
            )  # fmt: skip
        ax.plot(
            used["longitude"], used["latitude"], ls="none", marker=".",
            ms=3, color="tab:red", transform=pc, zorder=3,
            label=f"Used for fitting ({within})",
        )  # fmt: skip
        geod = Geodesic()
        labelled = set()
        for group, colour, label in (
            (shown, "tab:red", "Reference used for fitting"),
            (unused, BEFORE, "Reference not used (no float wind)"),
        ):
            for name, positions in group.items():
                if min_r is not None and (
                    footprints is not None
                    and name in footprints["reference"].to_numpy()
                ):
                    f = footprints.sel(reference=name)
                    ax.contour(
                        f["longitude"], f["latitude"], f, levels=[min_r],
                        colors=[colour], linestyles="--",
                        linewidths=0.9, transform=pc, zorder=3,
                    )  # fmt: skip
                for k, (la, lo) in enumerate(positions):
                    if min_r is None:
                        ring = np.asarray(geod.circle(lo, la, radius_km * 1e3))
                        ax.plot(
                            np.append(ring[:, 0], ring[0, 0]),
                            np.append(ring[:, 1], ring[0, 1]),
                            color=colour, ls="--" if k == 0 else ":",
                            lw=0.8, transform=pc, zorder=3,
                        )  # fmt: skip
                    ax.plot(
                        lo, la, ls="none", marker="^",
                        ms=9 if k == 0 else 6, color=colour, mec="white",
                        mew=0.8, transform=pc, zorder=5,
                        label=None if label in labelled else label,
                    )  # fmt: skip
                    labelled.add(label)
                la, lo = positions[0]
                ax.annotate(
                    name, (lo, la), xytext=(7, 4),
                    xycoords=pc._as_mpl_transform(ax),
                    textcoords="offset points", fontsize=9,
                    fontweight="bold", color=AFTER if colour != BEFORE
                    else BEFORE, zorder=6,
                )  # fmt: skip
        if len(p):
            ax.annotate(
                "START", (p["longitude"].iloc[0], p["latitude"].iloc[0]),
                xytext=(6, -12), xycoords=pc._as_mpl_transform(ax),
                textcoords="offset points", fontsize=8, zorder=6,
            )  # fmt: skip
        ax.legend(loc="best", fontsize=8)
        if annotate:
            several = [k for k, v in shown.items() if len(v) > 1]
            note = (
                f"{len(used)} of {len(p)} measurements used for fitting, "
                f"{within.replace('inside footprint', 'inside the footprint')} "
                "of " + (", ".join(shown) if shown else "no reference") + "."
            )
            if unused:
                note += (
                    " Passed without a float wind, not used: "
                    + ", ".join(unused)
                    + "."
                )
            if several:
                note += (
                    " Reported at several positions (one circle each): "
                    + ", ".join(several)
                    + "."
                )
            _annotate(fig, f"{_wmo(park_table)}float and references", note)
    return fig, ax


# --- Step 5: gridded_profiles ------------------------------------------------------


def plot_profiles(ts_profiles, annotate=False):
    """The T/S profiles used for the depth correction, coloured by date.

    Solid: measured; dotted: held constant above the shallowest or below
    the deepest measurement.

    Parameters
    ----------
    ts_profiles : pd.DataFrame
        Output of `gridded_profiles` (`result["ts_profiles"]`).
    annotate : bool, optional
        Add a title and notes (incl. skipped profiles), by default False.

    Returns
    -------
    fig, (ax_temperature, ax_salinity)
    """
    first = ts_profiles.groupby("cycle")["datetime"].first().sort_values()
    norm, colors = _time_colors(first.values)
    with _style():
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(12, 6), sharey=True, layout="constrained"
        )
        for (cycle, _), color in zip(first.items(), colors):
            g = ts_profiles[ts_profiles["cycle"] == cycle]
            for ax, var in ((ax1, "temperature"), (ax2, "salinity")):
                ax.plot(
                    g[var].where(g["measured"]), g["pressure"], color=color,
                    lw=1,
                )  # fmt: skip
                ax.plot(
                    g[var].where(~g["measured"]), g["pressure"], color=color,
                    lw=1, ls=":",
                )  # fmt: skip
        _depth_axis(ax1, ts_profiles["pressure"].max())
        ax1.set_ylabel("Pressure (dbar)")
        ax1.set_xlabel("Temperature (°C)")
        ax2.set_xlabel("Salinity")
        _time_colorbar(fig, [ax1, ax2], norm)
        _panel_labels([ax1, ax2])
        if annotate:
            cycles = ts_profiles.attrs.get("cycles", list(first.index))
            skipped = ts_profiles.attrs.get("skipped") or {}
            note = (
                "Solid: measured. Dotted: held constant beyond the measured "
                "range."
            )
            if skipped:
                note += " Skipped: " + "; ".join(
                    f"cycle {c}: {why}" for c, why in skipped.items()
                )
            _annotate(
                fig,
                f"{_wmo(ts_profiles)}{len(cycles)} profiles used for the "
                f"depth correction (mode {ts_profiles.attrs.get('mode')!r})",
                note,
            )
    return fig, (ax1, ax2)


# --- Step 6: compute_profile_betas / compute_beta --------------------------------


def plot_beta(beta, band=None, depth=None, annotate=False):
    """The depth correction β: (a) against depth for every profile,
    coloured by date; (b) at one depth (1000 m) through the deployment.
    For a single (mean) β: β against depth, 3.15 kHz solid and `band`
    dashed, as in Delaigue et al. (2026, Fig. 2b).

    Parameters
    ----------
    beta : pd.DataFrame
        Output of `compute_profile_betas` or `compute_beta`
        (`result["beta"]`).
    band : float, optional
        Band centre frequency (Hz); by default 8 kHz if processed, else
        the first processed band.
    depth : float, optional
        Depth (m) for panel (b); by default 1000 m (or the deepest β depth
        common to all profiles, if shallower).
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, axes
        Two panels for per-profile β, one for a single (mean) β.
    """
    if band is None:
        cols = [c for c in beta if c.startswith("beta_")]
        band = 8000 if "beta_8000Hz" in cols else int(cols[0][5:-2])
    col = f"beta_{int(band)}Hz"
    if col not in beta:
        raise KeyError(f"No β for {band} Hz in this table.")
    with _style():
        if "cycle" not in beta:
            fig, ax = plt.subplots(figsize=(6, 6), layout="constrained")
            if "beta_3150Hz" in beta and band != 3150:
                ax.plot(
                    beta["beta_3150Hz"], beta["depth"], color=AFTER,
                    label="3.15 kHz",
                )  # fmt: skip
            ax.plot(
                beta[col],
                beta["depth"],
                color=AFTER,
                ls="--",
                label=_khz(band),
            )
            _depth_axis(ax, beta["depth"].max())
            ax.set_xlim(left=0)
            ax.set_xlabel("β (dB)")
            ax.set_ylabel("Depth (m)")
            ax.legend(loc="upper right")
            if annotate:
                _annotate(fig, f"{_wmo(beta)}β from the mean profile")
            return fig, [ax]

        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(14, 6), width_ratios=[1, 1.5], layout="constrained"
        )
        fig.get_layout_engine().set(wspace=0.08)
        first = beta.groupby("cycle")["datetime"].first().sort_values()
        norm, colors = _time_colors(first.values)
        for (cycle, _), color in zip(first.items(), colors):
            g = beta[beta["cycle"] == cycle]
            ax1.plot(g[col], g["depth"], color=color, lw=1)
        _depth_axis(ax1, beta["depth"].max())
        ax1.set_xlim(left=0)
        ax1.set_xlabel(f"β at {_khz(band)} (dB)")
        ax1.set_ylabel("Depth (m)")
        _time_colorbar(fig, ax1, norm)

        if depth is None:
            deepest = beta.groupby("cycle")["depth"].max().min()
            depth = 1000.0 if deepest >= 1000 else deepest
        at = (
            beta.assign(_d=(beta["depth"] - depth).abs())
            .sort_values("_d")
            .groupby("cycle")
            .first()
            .sort_values("datetime")
        )
        ax2.plot(
            at["datetime"], at[col], color=AFTER, marker="o", ms=4,
            mfc="white", mew=1,
        )  # fmt: skip
        ax2.set_ylabel(f"β at {_khz(band)}, {depth:.0f} m (dB)")
        _date_axis(ax2)
        _panel_labels([ax1, ax2])
        if annotate:
            _annotate(
                fig,
                f"{_wmo(beta)}depth correction β per profile (absorption in "
                f"{beta.attrs.get('absorption_units', 'dB')}/m, mode "
                f"{beta.attrs.get('mode')!r})",
            )
    return fig, [ax1, ax2]


# --- Steps 6-7: depth_correct / spectral_level ---------------------------------------


def plot_depth_correction(park_table, band=None, beta=None, annotate=False):
    """Depth correction: (a) β applied to each PAL measurement, with the
    dates of the profiles used (ticks) and the periods without β (grey); (b) the band
    level at park depth (grey) and the corrected spectral level (black).

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `depth_correct` (and ideally `spectral_level`).
    band : float, optional
        Band centre frequency (Hz); by default 8 kHz if processed, else
        the first processed band.
    beta : pd.DataFrame, optional
        Output of `compute_profile_betas`, to mark the profile dates.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, (ax_beta, ax_level)
    """
    band = _default_band(park_table, band)
    col = _band(park_table, band)
    bcol = f"beta_{int(band)}Hz"
    scol = f"spl_{int(band)}Hz"
    df = park_table.sort_values("datetime")
    dc = df.attrs.get("depth_correction", {})
    with _style():
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(12, 8), sharex=True, height_ratios=[1, 1.4],
            layout="constrained",
        )  # fmt: skip
        ax1.plot(
            df["datetime"], df[bcol], color=AFTER, lw=1.5, label="β applied"
        )
        if beta is not None and "cycle" in beta:
            times = beta.groupby("cycle")["datetime"].first()
            y0 = np.nanmin(df[bcol]) if df[bcol].notna().any() else 0
            ax1.plot(
                times, np.full(len(times), y0), ls="none", marker="|", ms=10,
                mew=1.2, color=AFTER, label="Profiles used",
            )  # fmt: skip
        missing = df[bcol].isna()
        if missing.any():
            ax1.fill_between(
                df["datetime"], 0, 1, where=missing,
                transform=ax1.get_xaxis_transform(), color=BEFORE,
                alpha=0.25, lw=0,
                label=(
                    f"No β: no {dc.get('mode', '')}-mode profile within "
                    f"{_days(dc.get('max_extrapolation'))}"
                ),
            )  # fmt: skip
        ax1.set_ylabel(f"β at {_khz(band)} (dB)")
        ax1.legend(loc="upper left")

        series = df.set_index("datetime")
        ax2.plot(
            series[col].resample("6h").median(), color=BEFORE, lw=1,
            label="Band level at park depth",
        )  # fmt: skip
        units = "dB"
        if scol in df:
            ax2.plot(
                series[scol].resample("6h").median(), color=AFTER, lw=1.2,
                label="Corrected level at the surface",
            )  # fmt: skip
            units = df.attrs.get("spl", {}).get("units", "dB")
        units = units.replace("uPa^2/Hz", "µPa² Hz⁻¹").replace("uPa", "µPa")
        ax2.set_ylabel(f"{_khz(band)} level, 6-h median ({units})")
        ax2.legend(loc="upper right")
        _date_axis(ax2)
        _panel_labels([ax1, ax2])
        if annotate:
            spl = df.attrs.get("spl", {})
            note = None
            if spl:
                on = spl.get("band_normalization_on_float")
                sub = spl.get("bandwidth_subtracted")
                note = (
                    f"Band normalisation on the float: {'yes' if on else 'no'}"
                    f"; bandwidth term "
                    f"{'subtracted' if sub else 'not subtracted'}."
                )
            _annotate(
                fig,
                f"{_wmo(df)}depth correction ({dc.get('method', '?')} β, "
                f"mode {dc.get('mode')!r})",
                note,
            )
    return fig, (ax1, ax2)


_WIND_ORDER = ["vagle1990", "nystuen2015", "pensieri2015"]
_WIND_SHADES = dict(zip(_WIND_ORDER, (0.1, 0.5, 0.7)))  # as in Fig. 4


def _reference_label(park_table, reference):
    """Label of a reference: "Reference" for a set of references or the
    default name, else its name in capitals."""
    info = park_table.attrs.get("references", {}).get(reference, {})
    if reference == "ref" or len(info.get("references", {})) > 1:
        return "Reference"
    return reference.upper()


def _fit_reference(park_table):
    ref = park_table.attrs.get("wind_fits", {}).get("reference")
    if ref is None:
        refs = park_table.attrs.get("references", {})
        ref = next(iter(refs), None)
    return ref


def plot_wind(park_table, resample="6h", reference=None, annotate=False):
    """Wind speed from each model through the deployment, as in figure 4
    of Delaigue et al. (2026): (a) published coefficients; (b) refitted
    coefficients, if `fit_wind_models` was run, shaded ±1σ if
    `wind_uncertainty` was run.  The reference wind, if collocated, is
    drawn in black.

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `wind_literature` (`result["park"]`).
    resample : str or None, optional
        Median over this period before plotting, by default "6h"; None
        plots every measurement.
    reference : str, optional
        Reference to draw; by default the one used for fitting.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, axes
    """
    info = park_table.attrs.get("wind_models", {})
    applied = info.get("applied", {})
    if not applied:
        raise KeyError("No wind columns: run wind_literature first.")
    names = [m for m in _WIND_ORDER if m in applied] + [
        m for m in applied if m not in _WIND_ORDER
    ]
    fitted = [m for m in names if f"U_{m}_fit" in park_table]
    reference = reference or _fit_reference(park_table)
    cmap = plt.get_cmap("viridis")
    series = park_table.sort_values("datetime").set_index("datetime")

    def _series(col):
        u = series[col]
        return u.resample(resample).median() if resample else u

    panels = [("", "published")] + ([("_fit", "refitted")] if fitted else [])
    with _style():
        fig, axes = plt.subplots(
            len(panels), 1, figsize=(12, 4.5 * len(panels)), sharex=True,
            layout="constrained", squeeze=False,
        )  # fmt: skip
        axes = list(axes[:, 0])
        for ax, (suffix, label) in zip(axes, panels):
            if reference and f"U_{reference}" in series:
                u = _series(f"U_{reference}")
                ax.plot(
                    u.index, u, color=AFTER, lw=1.2,
                    label=_reference_label(park_table, reference),
                )  # fmt: skip
            for name in names if not suffix else fitted:
                color = cmap(_WIND_SHADES.get(name, 0.3))
                u = _series(f"U_{name}{suffix}")
                ax.plot(
                    u.index, u, color=color, lw=1.2,
                    label=applied[name]["reference"],
                )  # fmt: skip
                sig_col = f"U_{name}_fit_sigma"
                if suffix and sig_col in series:
                    sig = _series(sig_col).reindex(u.index)
                    ax.fill_between(
                        u.index, u - sig, u + sig, color=color, alpha=0.15,
                        lw=0,
                    )  # fmt: skip
            ax.set_ylim(bottom=0)
            ax.set_ylabel(
                f"Wind speed, {label} (m s$^{{-1}}$)"
                + (
                    f", {resample.replace('h', '-h')} median"
                    if resample
                    else ""
                )
            )
            ax.legend(loc="upper right", ncols=2)
        _date_axis(axes[-1])
        if len(axes) > 1:
            _panel_labels(axes)
        if annotate:
            skipped = info.get("skipped", {})
            note = (
                "Not applied: "
                + "; ".join(f"{m} ({why})" for m, why in skipped.items())
                if skipped
                else None
            )
            _annotate(fig, f"{_wmo(park_table)}wind speed", note)
    return fig, axes


def plot_wind_fit(
    park_table,
    reference=None,
    max_distance_km=80.0,
    average=None,
    window=None,
    annotate=False,
):
    """Refitted float wind against the reference wind, one panel per model,
    coloured by distance to the reference (as in Delaigue et al. 2026,
    Fig. 5), with the 1:1 line and the regression slope, intercept and R².

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `fit_wind_models`.
    reference : str, optional
        Reference to compare with; by default the one used for fitting.
    max_distance_km : float, optional
        Colour scale limit and largest distance shown, by default 80 km.
    average : bool, optional
        One point per reference value, the mean of the measurements
        around it (`average_to_reference`); by default as for
        `fit_wind_models`.
    window : str, optional
        Averaging window around each reference time; by default as for
        `fit_wind_models`.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, axes
    """
    from .calibration import average_to_reference, wind_metrics

    reference = reference or _fit_reference(park_table)
    ref_col = f"U_{reference}"
    fitted = [m for m in _WIND_ORDER if f"U_{m}_fit" in park_table]
    if not fitted or ref_col not in park_table:
        raise KeyError("Run collocate_reference and fit_wind_models first.")
    if average is None:
        average = park_table.attrs.get("wind_fits", {}).get("average", True)
    if window is None:
        window = park_table.attrs.get("wind_fits", {}).get("window", "1h")
    dist = park_table.get(f"{reference}_distance_km")
    min_r = park_table.attrs.get("wind_fits", {}).get("min_footprint_r")
    df = average_to_reference(
        park_table,
        reference,
        max_distance_km if dist is not None and min_r is None else None,
        average=average,
        window=window,
        min_footprint_r=min_r,
    )
    refs = {
        m: park_table.attrs["wind_models"]["applied"][m]["reference"]
        for m in fitted
    }
    lim = float(
        np.nanmax(df[[ref_col] + [f"U_{m}_fit" for m in fitted]].to_numpy())
    )
    lim = np.ceil(lim / 5) * 5
    label = _reference_label(park_table, reference)
    with _style():
        fig, axes = plt.subplots(
            1, len(fitted), figsize=(5 * len(fitted) + 1, 5),
            layout="constrained", squeeze=False,
        )  # fmt: skip
        axes = list(axes[0])
        sc = None
        for ax, name in zip(axes, fitted):
            y = df[f"U_{name}_fit"]
            c = df[f"{reference}_distance_km"] if dist is not None else None
            sc = ax.scatter(
                df[ref_col], y, c=c, cmap="viridis", vmin=0,
                vmax=max_distance_km, s=12, alpha=0.7,
            )  # fmt: skip
            ax.plot([0, lim], [0, lim], color=AFTER, ls="--", lw=1)
            m = wind_metrics(y, df[ref_col])
            ax.text(
                0.04, 0.96,
                f"Slope: {m['slope']:.2f}\nIntercept: {m['intercept']:.2f}"
                f"\nR²: {m['r2']:.2f}\nRMSE: {m['rmse']:.2f}",
                transform=ax.transAxes, va="top",
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "#cccccc"},
            )  # fmt: skip
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_aspect("equal")
            ax.set_xlabel(f"{label} wind speed (m s$^{{-1}}$)")
            ax.set_title(refs[name], fontsize=10)
        axes[0].set_ylabel("Float wind speed, refitted (m s$^{-1}$)")
        _panel_labels(axes)
        if dist is not None and sc is not None:
            cb = fig.colorbar(sc, ax=axes, pad=0.02, aspect=30)
            cb.set_label(
                f"Distance from {label.replace('Reference', 'reference')} (km)"
            )
        if annotate:
            fits = park_table.attrs.get("wind_fits", {})
            _annotate(
                fig,
                f"{_wmo(park_table)}refitted models vs {reference}",
                f"Fitted on measurements within "
                f"{fits.get('max_distance_km')} km of the reference"
                + (
                    f", averaged over {fits.get('window', '1h')} around "
                    "each reference value."
                    if fits.get("average", True)
                    else "."
                ),
            )
    return fig, axes


def _break_gaps(t, y, gap):
    """Insert NaN where consecutive times are more than `gap` apart, so
    that lines are not drawn across gaps."""
    t = pd.to_datetime(pd.Series(t)).reset_index(drop=True)
    y = pd.Series(np.asarray(y, float))
    cut = np.flatnonzero(t.diff() > gap)
    if not len(cut):
        return t, y
    tt = list(t)
    yy = list(y)
    for k in cut[::-1]:
        tt.insert(k, t[k - 1] + (t[k] - t[k - 1]) / 2)
        yy.insert(k, np.nan)
    return pd.Series(tt), pd.Series(yy)


def plot_reference_timeseries(
    park_table,
    reference=None,
    model="nystuen2015",
    radius_km=None,
    average=None,
    window=None,
    gap="3h",
    uncertainty=True,
    annotate=False,
):
    """Wind against time where the float is near a reference: the
    reference wind, and the float wind from one model with its published
    and its refitted coefficients.

    The points are those compared with the reference when fitting: the
    measurements within `radius_km` of the reference they were matched to
    (the closest one with a value in time) and with a float wind (none,
    e.g., after the delayed-mode data end), averaged around each reference
    value (`average_to_reference`).  One row per reference, each over its
    own period, as the float is near each at different times.

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `collocate_reference` (and `fit_wind_models` for
        the refitted wind).
    reference : str, optional
        Name given to `collocate_reference`; by default the one used for
        fitting, or the only one.
    model : str, optional
        Wind model, by default "nystuen2015".
    radius_km : float, optional
        By default the one used by `fit_wind_models`, else 40 km.
    average, window : optional
        Averaging around each reference value; by default as for
        `fit_wind_models` (see `average_to_reference`).
    gap : str, optional
        Lines are broken where points are further apart in time, by
        default "3h".
    uncertainty : bool, optional
        Shade ±1σ around the refitted wind (the combined uncertainty,
        "U_<model>_fit_sigma", see `wind_uncertainty`), if computed; by
        default True.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, axes
    """
    from .calibration import average_to_reference, wind_metrics

    reference = reference or _fit_reference(park_table)
    ref_col = f"U_{reference}"
    if ref_col not in park_table:
        raise KeyError("No reference wind: run collocate_reference first.")
    pub_col, fit_col = f"U_{model}", f"U_{model}_fit"
    if pub_col not in park_table:
        raise KeyError(f"No {pub_col!r}: run wind_literature first.")
    fits = park_table.attrs.get("wind_fits", {})
    if radius_km is None:
        radius_km = fits.get("max_distance_km")
        radius_km = 40.0 if radius_km is None else radius_km
    if average is None:
        average = fits.get("average", True)
    if window is None:
        window = fits.get("window", "1h")
    min_r = fits.get("min_footprint_r")
    pts = average_to_reference(
        park_table,
        reference,
        None if min_r is not None else radius_km,
        average=average,
        window=window,
        min_footprint_r=min_r,
    )
    # only points with a float wind: those compared in the refit (e.g. not
    # after the delayed-mode data end)
    pts = pts.dropna(subset=[ref_col, "datetime", pub_col])
    if pts.empty:
        raise ValueError(
            f"No measurement within {radius_km:g} km of the reference "
            "with a float wind."
        )
    name_col = f"{reference}_reference_name"
    groups = (
        list(pts.groupby(name_col, sort=False)) if name_col in pts
        else [(_reference_label(park_table, reference), pts)]
    )  # fmt: skip
    groups.sort(key=lambda g: g[1]["datetime"].min())
    applied = park_table.attrs.get("wind_models", {}).get("applied", {})
    label = applied.get(model, {}).get("reference", model)
    colour = plt.get_cmap("viridis")(_WIND_SHADES.get(model, 0.3))
    gap = pd.Timedelta(gap)
    with _style():
        fig, axes = plt.subplots(
            len(groups), 1, figsize=(12, 3.2 * len(groups) + 0.5),
            sharey=True, layout="constrained", squeeze=False,
        )  # fmt: skip
        axes = list(axes[:, 0])
        for ax, (name, g) in zip(axes, groups):
            g = g.sort_values("datetime")
            series = [
                (ref_col, AFTER, "-", f"{name} (reference)"),
                (pub_col, BEFORE, "--", f"{label}, published"),
            ]
            if fit_col in g:
                series.append((fit_col, colour, "-", f"{label}, refitted"))
            sig_col = f"{fit_col}_sigma"
            if uncertainty and fit_col in g and sig_col in g:
                t, lo = _break_gaps(
                    g["datetime"], g[fit_col] - g[sig_col], gap
                )
                _, hi = _break_gaps(
                    g["datetime"], g[fit_col] + g[sig_col], gap
                )
                ax.fill_between(
                    t, lo.clip(lower=0), hi, color=colour, alpha=0.2,
                    lw=0, label=f"{label}, refitted ±1σ",
                )  # fmt: skip
            for col, c, ls, lab in series:
                t, y = _break_gaps(g["datetime"], g[col], gap)
                ax.plot(t, y, color=c, ls=ls, lw=1.2, label=lab)
            d = g[f"{reference}_distance_km"]
            ax.set_title(
                f"{name}: {len(g)} points, {d.min():.0f}–{d.max():.0f} km",
                fontsize=10,
            )
            locator = mdates.AutoDateLocator(maxticks=10)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
            ax.set_ylabel("Wind speed (m s$^{-1}$)")
            handles, labels = ax.get_legend_handles_labels()
            order = sorted(
                range(len(labels)), key=lambda i: "±1σ" in labels[i]
            )  # the shading last
            ax.legend(
                [handles[i] for i in order], [labels[i] for i in order],
                loc="upper left", fontsize=8, ncol=4,
            )  # fmt: skip
        cols = [c for c in (ref_col, pub_col, fit_col) if c in pts]
        top = np.nanmax(pts[cols].to_numpy(float))
        axes[0].set_ylim(0, np.ceil(top * 1.05))
        _panel_labels(axes)
        if annotate:
            parts = []
            for col, what in ((pub_col, "published"), (fit_col, "refitted")):
                if col in pts:
                    m = wind_metrics(pts[col], pts[ref_col])
                    parts.append(
                        f"{what}: RMSE {m['rmse']:.2f}, bias {m['bias']:+.2f}"
                    )
            how = (
                f"averaged over {window} around each reference value"
                if average
                else "single measurements"
            )
            _annotate(
                fig,
                f"{_wmo(park_table)}{label} against the reference",
                f"Within {radius_km:g} km, {how}; "
                + "; ".join(parts) + " (m s$^{-1}$).",
            )  # fmt: skip
    return fig, axes


def plot_reference_fits(
    park_table, reference=None, model="nystuen2015", annotate=False
):
    """Each reference's own fit, against the fit to all of them.

    The float's spectral level against the reference wind (one colour per
    reference, the points used in the refit), with the wind model fitted
    to each reference alone, to all of them (black), and with its
    published coefficients (grey, dashed).  References needing clearly
    different curves disagree about the float-wind relationship (see
    `reference_diagnostics`).

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `collocate_reference` and `fit_wind_models`.
    reference : str, optional
        By default the one used for fitting.
    model : str, optional
        Wind model, by default "nystuen2015".
    annotate : bool, optional
        Add a title and a note with each reference's RMSE.

    Returns
    -------
    fig, ax
    """
    from .calibration import average_to_reference, reference_diagnostics
    from .wind import LITERATURE_MODELS

    diag = reference_diagnostics(park_table, reference, model)
    reference = diag["fits"].attrs["reference"]
    _, band, _ = LITERATURE_MODELS[model]
    fits = park_table.attrs.get("wind_fits", {})
    pts = average_to_reference(
        park_table,
        reference,
        fits.get("max_distance_km", 40.0),
        average=fits.get("average", True),
        window=fits.get("window", "1h"),
        min_footprint_r=fits.get("min_footprint_r"),
    ).dropna(subset=[f"spl_{band}Hz", f"U_{reference}"])
    name_col = f"{reference}_reference_name"
    curves = diag["curves"]
    names = [c for c in curves if c not in ("spl", "published", "pooled")]
    colours = plt.get_cmap("tab10")
    applied = park_table.attrs.get("wind_models", {}).get("applied", {})
    label = applied.get(model, {}).get("reference", model)
    with _style():
        fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
        for i, name in enumerate(names):
            g = (
                pts[pts[name_col].astype(str) == name]
                if name_col in pts
                else pts
            )
            c = colours(i % 10)
            ax.plot(
                g[f"spl_{band}Hz"], g[f"U_{reference}"], ls="none",
                marker=".", ms=4, alpha=0.4, color=c,
            )  # fmt: skip
            ax.plot(
                curves["spl"], curves[name], color=c, lw=2,
                label=f"{name} ({len(g)} points)",
            )  # fmt: skip
        ax.plot(
            curves["spl"], curves["pooled"], color=AFTER, lw=2.5,
            label="All references (refit)",
        )  # fmt: skip
        ax.plot(
            curves["spl"], curves["published"], color=BEFORE, lw=1.5,
            ls="--", label=f"{label}, published",
        )  # fmt: skip
        ax.set_xlabel(
            f"Spectral level at {_khz(band)} (dB re 1 µPa$^2$ Hz$^{{-1}}$)"
        )
        ax.set_ylabel("Reference wind speed (m s$^{-1}$)")
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper left", fontsize=8)
        if annotate:
            t = diag["fits"].set_index("reference")
            note = "; ".join(
                f"{n}: RMSE own {t.at[n, 'rmse_own']:.2f}, all "
                f"{t.at[n, 'rmse_pooled']:.2f} (cross-validated "
                f"{t.at[n, 'rmse_own_cv']:.2f}, {t.at[n, 'rmse_pooled_cv']:.2f})"
                for n in names
            )
            _annotate(
                fig, f"{_wmo(park_table)}{label} fitted per reference",
                note + " m s$^{-1}$.",
            )  # fmt: skip
    return fig, ax


def plot_ais(park_table, band=None, annotate=False):
    """Ship-noise check: (a) the band level, with the measurements that had
    a ship nearby; (b) the number of ships nearby through the deployment.

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `flag_ships` (`result["park"]`).
    band : float, optional
        Band centre frequency (Hz); by default 8 kHz if processed.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, (ax_level, ax_ships)
    """
    if "ships_nearby" not in park_table:
        raise KeyError("No ship columns: run flag_ships first.")
    band = _default_band(park_table, band)
    info = park_table.attrs.get("ais", {})
    df = park_table.sort_values("datetime")
    # band level as sent (flagged values may have been removed)
    col = _band(df, band)
    flagged = df[df["ship_flag"]]
    with _style():
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(12, 7), sharex=True, height_ratios=[1.4, 1],
            layout="constrained",
        )  # fmt: skip
        ax1.plot(df["datetime"], df[col], color=BEFORE, lw=0.6, label="All")
        ax1.plot(
            flagged["datetime"], flagged[col], ls="none", marker="o", ms=2,
            color=AFTER, label=(
                f"Ship within {info.get('radius_km', 20):g} km, "
                f"±{_hours(info.get('window', '30min')).replace(' h', '')} h"
                if pd.Timedelta(info.get("window", "30min")) >= pd.Timedelta("1h")
                else f"Ship within {info.get('radius_km', 20):g} km, "
                f"±{pd.Timedelta(info.get('window', '30min')).total_seconds() / 60:g} min"
            ),
        )  # fmt: skip
        ax1.set_ylabel(f"{_khz(band)} level (dB)")
        ax1.legend(loc="upper right", markerscale=3)
        ax2.plot(
            df["datetime"], df["ships_nearby"], color=AFTER, lw=1,
            drawstyle="steps-mid",
        )  # fmt: skip
        unknown = df["ships_nearby"].isna()
        if unknown.any():
            ax2.fill_between(
                df["datetime"], 0, 1, where=unknown,
                transform=ax2.get_xaxis_transform(), color=BEFORE,
                alpha=0.25, lw=0, label="No AIS coverage",
            )  # fmt: skip
            ax2.legend(loc="upper right")
        ax2.set_ylabel("Ships nearby")
        ax2.set_ylim(bottom=0)
        ax2.yaxis.get_major_locator().set_params(integer=True)
        _date_axis(ax2)
        _panel_labels([ax1, ax2])
        if annotate:
            n = len(df)
            _annotate(
                fig,
                f"{_wmo(df)}ship-noise check (AIS, action "
                f"{info.get('action')!r})",
                f"Flagged: {info.get('n_flagged', 0):,} of {n:,} measurements; "
                f"without AIS coverage: {info.get('n_without_coverage', 0):,}; "
                f"{info.get('n_vessels', 0):,} vessels in "
                f"{info.get('n_positions', 0):,} AIS positions.",
            )
    return fig, (ax1, ax2)


def plot_all(
    result, band=None, cycle=None, save_dir=None, dpi=DPI, annotate=False
):
    """Every step's figure for a `process_float` result.

    Parameters
    ----------
    result : dict
        Output of `process_float`.
    band : float, optional
        Band shown where one band is plotted; by default 8 kHz if it was
        processed, else the first processed band.
    cycle : int, optional
        Cycle shown in detail; by default the one with most measurements.
    save_dir : str, optional
        If given, every figure is saved there as `<step>.png`.
    dpi : int, optional
        Resolution of saved figures, by default 300.
    annotate : bool, optional
        Add titles and notes to every figure, by default False.

    Returns
    -------
    dict
        Step name -> figure.
    """
    a = {"annotate": annotate}
    figs = {
        "spectra": plot_pal_spectra(result["pal_raw"], **a)[0],
        "despike": plot_despike(
            result["pal_raw"], result["park"], band=band, cycle=cycle, **a
        )[0],
        "park": plot_park(
            result["park"], result["park_samples"], cycle=cycle, **a
        )[0],
        "positions": plot_positions(result["park"], result["positions"], **a)[
            0
        ],
        "profiles": plot_profiles(result["ts_profiles"], **a)[0],
        "beta": plot_beta(result["beta"], band=band, **a)[0],
        "depth_correction": plot_depth_correction(
            result["park"], band=band, beta=result["beta"], **a
        )[0],
    }
    if "ships_nearby" in result["park"]:
        figs["ais"] = plot_ais(result["park"], band=band, **a)[0]
    if result["park"].attrs.get("wind_models", {}).get("applied"):
        figs["wind"] = plot_wind(result["park"], **a)[0]
    if result["park"].attrs.get("references"):
        figs["reference_map"] = plot_reference_map(
            result["park"], footprints=result.get("footprints"), **a
        )[0]
    if any(c.endswith("_fit") for c in result["park"]):
        figs["wind_fit"] = plot_wind_fit(result["park"], **a)[0]
        figs["reference_timeseries"] = plot_reference_timeseries(
            result["park"], **a
        )[0]
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        for name, fig in figs.items():
            save_figure(fig, os.path.join(save_dir, f"{name}.png"), dpi=dpi)
    return figs
