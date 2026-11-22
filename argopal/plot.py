"""Figures to follow each processing step.

One function per step; each draws one or two panels and returns
``(fig, axes)``, so figures can be shown, adjusted or saved:

    fig, ax = ap.plot_despike(result["pal_raw"], result["park"])
    ap.save_figure(fig, "despike.png")        # 300 dpi by default

    ap.plot_all(result, save_dir="figures")   # every step, saved at 300 dpi

Style follows Delaigue et al. (2026, Ocean Science): white background,
framed axes, light grid, black lines (grey = before a step, dashed = a
second series), panel labels a), b), framed legends.  Many profiles at once
are coloured by date (viridis).  Pass ``annotate=True`` to add a title and
notes (float, data mode, what was skipped and why).
"""

import os
from contextlib import contextmanager

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
    with plt.rc_context(_RC):
        yield


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


def plot_positions(park_table, positions, annotate=False):
    """Float track: surface fixes (open circles) and the positions
    interpolated onto the PAL measurements (line).

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
    with _style():
        fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
        p = park_table.dropna(subset=["latitude", "longitude"])
        ax.plot(
            p["longitude"], p["latitude"], color=AFTER, lw=1.5,
            label="Interpolated at PAL measurements",
        )  # fmt: skip
        ax.plot(
            positions["longitude"], positions["latitude"], ls="none",
            marker="o", ms=5, mfc="white", mec=AFTER, mew=1,
            label="Surface fixes",
        )  # fmt: skip
        first, last = positions.iloc[0], positions.iloc[-1]
        for row, text, dy in ((first, "START", 6), (last, "LATEST", -14)):
            ax.annotate(
                text, (row["longitude"], row["latitude"]), xytext=(6, dy),
                textcoords="offset points", fontsize=9, fontweight="bold",
            )  # fmt: skip
        lat0 = np.nanmean(positions["latitude"])
        ax.set_aspect(1 / np.cos(np.radians(lat0)))
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
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
    """(a) β applied to each PAL measurement, with the dates of the
    profiles used (ticks) and the periods without β (grey); (b) the band
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


def plot_wind(park_table, resample="6h", annotate=False):
    """Wind speed from each literature model, through the deployment
    (as in Delaigue et al. 2026, Fig. 4a).

    Parameters
    ----------
    park_table : pd.DataFrame
        PAL table after `wind_literature` (`result["park"]`).
    resample : str or None, optional
        Median over this period before plotting, by default "6h"; None
        plots every measurement.
    annotate : bool, optional
        Add a title and notes, by default False.

    Returns
    -------
    fig, ax
    """
    info = park_table.attrs.get("wind_models", {})
    applied = info.get("applied", {})
    if not applied:
        raise KeyError("No wind columns: run wind_literature first.")
    order = ["vagle1990", "nystuen2015", "pensieri2015"]
    shades = dict(zip(order, (0.1, 0.5, 0.7)))  # as in Fig. 4
    cmap = plt.get_cmap("viridis")
    series = park_table.sort_values("datetime").set_index("datetime")
    with _style():
        fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
        for name in [m for m in order if m in applied] + [
            m for m in applied if m not in order
        ]:
            u = series[f"U_{name}"]
            if resample:
                u = u.resample(resample).median()
            ax.plot(
                u.index, u, color=cmap(shades.get(name, 0.3)), lw=1.2,
                label=applied[name]["reference"],
            )  # fmt: skip
        ax.set_ylim(bottom=0)
        ax.set_ylabel(
            "Wind speed (m s$^{-1}$)"
            + (f", {resample.replace('h', '-h')} median" if resample else "")
        )
        ax.legend(loc="upper right", ncols=2)
        _date_axis(ax)
        if annotate:
            skipped = info.get("skipped", {})
            note = (
                "Not applied: "
                + "; ".join(f"{m} ({why})" for m, why in skipped.items())
                if skipped
                else None
            )
            _annotate(
                fig,
                f"{_wmo(park_table)}wind speed, literature coefficients",
                note,
            )
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
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        for name, fig in figs.items():
            save_figure(fig, os.path.join(save_dir, f"{name}.png"), dpi=dpi)
    return figs
