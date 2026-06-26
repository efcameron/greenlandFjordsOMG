"""
ts_diagram.py
-------------
Self-contained module for plotting Temperature-Salinity (TS) diagrams
for Greenland fjord CTD profiles.

Public API
----------
plot_ts(ax, profiles, fjord_row=None, ...)
    Plot a full TS diagram onto an existing Axes, including:
      - Density contours (sigma0)
      - Water mass labels and AWm shading
      - Gade (subglacial melt) lines
      - Freshwater mixing lines
      - Depth markers at standard depths
      - Sill and grounding line depth markers

Example
-------
    from ts_diagram import plot_ts

    fig, ax = plt.subplots(figsize=(8, 7))

    plot_ts(
        ax,
        profiles,           # dict with keys 'fjord_mouth', 'mid_fjord', 'near_glacier'
        fjord_row=meta_row, # pandas Series with lat, lon, sill_depth_m, grounding_line_depth_m
    )

    plt.savefig("ts_diagram.png", dpi=300, bbox_inches="tight")

Profile dict format
-------------------
Each key maps to a dict with 'depth', 'temp', 'sal' as 1-D numpy arrays, e.g.:

    profiles = {
        "fjord_mouth":  {"depth": ..., "temp": ..., "sal": ...},
        "mid_fjord":    {"depth": ..., "temp": ..., "sal": ...},
        "near_glacier": {"depth": ..., "temp": ..., "sal": ...},
    }
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import gsw
from shapely.geometry import Polygon
import pandas as pd

# ── Colour / label lookup ──────────────────────────────────────────────────────
_COLORS = {"fjord_mouth": "black", "mid_fjord": "orange", "near_glacier": "red"}
_LABELS = {"fjord_mouth": "Mouth", "mid_fjord": "Mid-fjord", "near_glacier": "Near-glacier"}

# Standard depth markers and their label offsets (in points)
_DEFAULT_MARKER_DEPTHS = (30, 130, 200, 300, 400)
_LABEL_OFFSETS = {
    30:  (0,   -20),
    130: (0,   -20),
    200: (-30,   7),
    300: (-30,  -4),
    400: (-30,  -4),
}

# ── Axis limits ────────────────────────────────────────────────────────────────
_XMIN, _XMAX = 31, 35.5 # was 30, 35
_YMIN, _YMAX = -2, 5 # was -1, 4.5


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

def _melt_temperature(S_melt: float) -> float:
    """
    Effective meltwater temperature (°C) following the three-equation approach.
    """
    L, c_w, c_i = 335_000, 3974, 2009
    T_i_K    = -15 + 273.15
    S_safe   = max(S_melt, 0.001)
    Theta_f_K = gsw.CT_freezing(S_safe, p=0, saturation_fraction=0) + 273.15
    return (Theta_f_K - L / c_w - (c_i / c_w) * (Theta_f_K - T_i_K)) - 273.15


def _plot_gade_lines(ax, S_farfield_vals, Theta_farfield: float = 10.0):
    """Overplot a family of Gade (subglacial melt) lines."""
    Theta_melt = _melt_temperature(0.0)
    for S_ff in S_farfield_vals:
        S = np.linspace(0.0, S_ff, 100)
        m = (Theta_farfield - Theta_melt) / (S_ff - 0.0)
        ax.plot(S, Theta_melt + m * S, "--", color="orange", alpha=0.5, lw=1)


def _plot_freshwater_mixing_lines(
    ax,
    p_grounding: float,
    seawater_salinity: float = 36,
    seawater_temps: tuple = (1, 3, 5, 7, 9, 11),
):
    """Freshwater mixing lines from subglacial runoff to ambient seawater."""
    theta_sub = gsw.CT_freezing(0, p=p_grounding, saturation_fraction=0)
    S_fw = np.linspace(seawater_salinity, 0, 100)
    for theta_ff in seawater_temps:
        m = (theta_ff - theta_sub) / (seawater_salinity - 0)
        ax.plot(S_fw, theta_sub + m * S_fw, "--", color="blue", alpha=0.2, lw=1)

def _add_depth_markers(
    ax,
    depths, salinity, temperature,
    color: str,
    marker_depths: tuple = _DEFAULT_MARKER_DEPTHS,
):
    """Plot circle markers at standard depths with offset depth labels."""
    if any(v is None for v in (depths, salinity, temperature)):
        return

    for d in marker_depths:
        idx = np.argmin(np.abs(depths - d))
        s, t = salinity[idx], temperature[idx]

        if np.isfinite(s) and np.isfinite(t):
            ax.scatter(s, t, color=color, edgecolors="black",
                       marker="o", s=80, zorder=10)
            dx, dy = _LABEL_OFFSETS.get(d, (6, 6))
            ax.annotate(
                f"{d} m", (s, t),
                textcoords="offset points", xytext=(dx, dy),
                ha="center", fontsize=16, color=color, zorder=11,
            )


def _add_geometry_markers(ax, depths, salinity, temperature, fjord_row):
    """
    Mark sill depth and grounding line depth on the TS diagram using
    the nearest binned data point.
    """
    if any(v is None for v in (depths, salinity, temperature)) or fjord_row is None:
        return

    for key, colour, label in [
        ("sill_depth_m",           "grey", "Sill depth"),
        ("grounding_line_depth_m", "blue", "Grounding line"),
    ]:
        if key in fjord_row and not np.isnan(fjord_row[key]):
            idx = np.argmin(np.abs(depths - fjord_row[key]))
            s, t = salinity[idx], temperature[idx]
            if np.isfinite(s) and np.isfinite(t):
                ax.scatter(s, t, color=colour, marker="o",
                           s=80, zorder=9, label=label)
                
def _plot_water_mass_markers(ax, water_mass_csv, fjord): # remove if we dont like the water mass markers
    """
    Plot literature water mass positions as black crosses with labels.
    Filters to the relevant fjord ('klu' or 'qim') from the TEOS-10 CSV.
    """
    df = pd.read_csv(water_mass_csv)
    df = df[df["fjord"] == fjord]

    for _, row in df.iterrows():
        ax.scatter(
            row["absolute_salinity_g_kg"],
            row["conservative_temperature"],
            marker="x", color="black", s=100, linewidths=1.5, zorder=12,
        )
        ax.text(
            row["absolute_salinity_g_kg"],
            row["conservative_temperature"] + 0.1,
            row["water_mass"],
            ha="center", va="bottom", fontsize=14, color="black", zorder=13,
        )


# ══════════════════════════════════════════════════════════════════════════════
# Public API
# ══════════════════════════════════════════════════════════════════════════════

def plot_ts(
    ax: plt.Axes,
    profiles: dict,
    fjord_row=None,
    *,
    marker_depths: tuple = _DEFAULT_MARKER_DEPTHS,
    gade_farfield_temp: float = 10.0,
    seawater_salinity: float = 36.0,
    seawater_temps: tuple = (1, 3, 5, 7, 9, 11),
    label_fontsize: int = 20,
    fjord=None, 
    water_mass_csv=None,
    xmin=_XMIN, xmax=_XMAX,        
    ymin=_YMIN, ymax=_YMAX,
  
) -> plt.Axes:
    """
    Plot a full TS diagram onto *ax*.

    Parameters
    ----------
    ax : matplotlib Axes
        Axes to draw on.
    profiles : dict
        Keyed by 'fjord_mouth', 'mid_fjord', 'near_glacier'. Each value is a
        dict with 'depth', 'temp', 'sal' arrays (binned, 1-D numpy arrays).
    fjord_row : pandas Series or dict, optional
        Metadata row containing at minimum:
          - 'latitude', 'longitude'
          - 'grounding_line_depth_m'
          - 'sill_depth_m'
        Used to draw mixing lines and geometry markers.
    marker_depths : tuple of int
        Depths (m) at which to plot labelled circle markers (default: 30, 130,
        200, 300, 400).
    gade_farfield_temp : float
        Far-field temperature used for Gade line family (default: 10 °C).
    seawater_salinity : float
        Far-field salinity for freshwater mixing lines (default: 36 g/kg).
    seawater_temps : tuple of float
        Far-field temperatures for the freshwater mixing line family.
    label_fontsize : int
        Font size for axis labels (default: 20).

    Returns
    -------
    ax : matplotlib Axes
        The same Axes, with the TS diagram drawn on it.
    """
    # ── Density contours ──────────────────────────────────────────────────────
    S_grid, T_grid = np.meshgrid(
        np.linspace(xmin, xmax, 100),
        np.linspace(ymin, ymax, 100),
    )
    D = gsw.sigma0(S_grid, T_grid)

    CS = ax.contour(
        S_grid, T_grid, D,
        levels=np.arange(np.floor(D.min()), np.ceil(D.max()) + 0.5, 0.5),
        colors="black", linewidths=0.5, zorder=0,
    )
    ax.clabel(CS, inline=True, fontsize=12, fmt="%.1f")


    if water_mass_csv is not None and fjord is not None:
        _plot_water_mass_markers(ax, water_mass_csv, fjord)

    # ── Profile scatter ───────────────────────────────────────────────────────
    draw_order = ["fjord_mouth", "mid_fjord", "near_glacier"]
    markers    = {"fjord_mouth": "ko", "mid_fjord": "o", "near_glacier": "ro"}

    for key in draw_order:
        p = profiles.get(key, {})
        s, t = p.get("sal"), p.get("temp")
        if s is not None:
            kw = dict(color="orange") if key == "mid_fjord" else {}
            ax.plot(s, t, markers[key], ms=5, label=_LABELS[key], **kw)

    # ── Depth markers on near-glacier profile ─────────────────────────────────
    ng = profiles.get("near_glacier", {})
    _add_depth_markers(
        ax,
        ng.get("depth"), ng.get("sal"), ng.get("temp"),
        color="red",
        marker_depths=marker_depths,
    )

    # ── Sill and grounding line markers ───────────────────────────────────────
    _add_geometry_markers(
        ax,
        ng.get("depth"), ng.get("sal"), ng.get("temp"),
        fjord_row,
    )

    # ── Gade lines ────────────────────────────────────────────────────────────
    _plot_gade_lines(ax, np.arange(26, 40, 1), Theta_farfield=gade_farfield_temp)

    # ── Freshwater mixing lines ───────────────────────────────────────────────
    if fjord_row is not None:
        lat    = fjord_row["latitude"]
        p_grnd = gsw.p_from_z(-fjord_row["grounding_line_depth_m"], lat)
        _plot_freshwater_mixing_lines(
            ax, p_grnd,
            seawater_salinity=seawater_salinity,
            seawater_temps=seawater_temps,
        )

    # ── Axes styling ──────────────────────────────────────────────────────────
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Absolute Salinity (g/kg)",      fontsize=label_fontsize)
    ax.set_ylabel("Conservative Temperature (°C)", fontsize=label_fontsize)
    ax.grid(False)
    ax.legend(loc="lower left")

    return ax