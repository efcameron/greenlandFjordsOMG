"""
fjord_case_studies.py
---------------------
Produces two separate 2×2 figures, one per fjord case study (QIM and KLU).

Layout per figure
-----------------
    [0,0] Temperature vs depth   |  [0,1] Salinity vs depth
    [1,0] TS diagram             |  [1,1] Bathymetry map

Dependencies
------------
Requires bathymetry_map.py and ts_diagram.py to be on your Python path.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4

from config_loader import load_paths
from bathymetry_map import load_bedmachine, plot_bathymetry, DEFAULTS
from ts_diagram import plot_ts

# ── Paths ──────────────────────────────────────────────────────────────────────
paths        = load_paths()
csv_dir      = paths["csv_dir"]
nc_dir       = paths["nc_dir"]
results_root = paths["results_dir"]

output_dir = results_root / "caseStudies_qimKlu"
output_dir.mkdir(parents=True, exist_ok=True)

qim_csv = csv_dir / "qim_fjord_data.csv"
klu_csv = csv_dir / "klu_fjord_data.csv"
water_mass_csv = paths["csv_dir"] / "rysgaard_water_masses_teos10.csv"

# ── Plot styling ───────────────────────────────────────────────────────────────
plt.rcParams["xtick.labelsize"] = 20
plt.rcParams["ytick.labelsize"] = 20
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["axes.titlesize"]  = 25
plt.rcParams["axes.labelsize"]  = 20

# ── Colour / label lookup ──────────────────────────────────────────────────────
COLORS = {"fjord_mouth": "black", "mid_fjord": "orange", "near_glacier": "red"}
LABELS = {"fjord_mouth": "Mouth", "mid_fjord": "Mid-fjord", "near_glacier": "Near-glacier"}


# ══════════════════════════════════════════════════════════════════════════════
# 1.  DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def isolate_descent_phase(depth, values):
    """Return only the downcast portion of a CTD profile."""
    deepest_idx = np.argmax(depth)
    return depth[: deepest_idx + 1], values[: deepest_idx + 1]


def load_nc_profile(nc_path):
    """Load depth, temperature, and salinity from a GSW NetCDF CTD file."""
    with netCDF4.Dataset(nc_path) as ds:
        depth = np.squeeze(ds.variables["depth"][:])
        temp  = np.squeeze(ds.variables["conservative_temperature"][:])
        sal   = np.squeeze(ds.variables["absolute_salinity"][:])

    depth = np.ma.masked_invalid(depth).filled(np.nan)
    temp  = np.ma.masked_invalid(temp).filled(np.nan)
    sal   = np.ma.masked_invalid(sal).filled(np.nan)
    return depth, temp, sal


def bin_profile(depth, values, bin_edges):
    """
    Bin profile data into fixed depth intervals (mean per bin).
    Requires >= 3 valid points per bin; otherwise returns NaN.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    binned = np.full(len(bin_centers), np.nan)

    for i in range(len(bin_edges) - 1):
        mask = (
            (depth >= bin_edges[i]) &
            (depth <  bin_edges[i + 1]) &
            np.isfinite(values)
        )
        if mask.sum() >= 3:
            binned[i] = np.nanmean(values[mask])

    return bin_centers, binned


def load_fjord_2016(csv_path, nc_dir):
    """
    Load, QC, and bin all 2016 profiles for one fjord CSV.

    Returns
    -------
    profiles : dict keyed by type ('fjord_mouth', 'mid_fjord', 'near_glacier'),
               each containing 'depth', 'temp', 'sal', and 'row' (metadata).
    """
    VALID_TYPES = {"fjord_mouth", "mid_fjord", "near_glacier"}

    df = pd.read_csv(csv_path)
    df["type"] = df["type"].str.strip().str.lower()
    df = df[(df["year"] == 2016) & df["type"].isin(VALID_TYPES)]

    profiles  = {}
    bin_edges = np.arange(1, 1001, 10)

    for _, row in df.iterrows():
        nc_path = os.path.join(nc_dir, row["GSW_filename"])
        if not os.path.exists(nc_path):
            print(f"Missing file: {row['GSW_filename']}")
            continue

        depth, temp, sal = load_nc_profile(nc_path)

        depth, temp = isolate_descent_phase(depth, temp)
        _,     sal  = isolate_descent_phase(depth, sal)

        sal[(sal < 30) | (sal > 40)] = np.nan
        temp[np.isnan(sal)]          = np.nan

        depth_bins, temp_binned = bin_profile(depth, temp, bin_edges)
        _,           sal_binned = bin_profile(depth, sal,  bin_edges)

        profiles[row["type"]] = {
            "depth": depth_bins,
            "temp":  temp_binned,
            "sal":   sal_binned,
            "row":   row,
        }

    return profiles


# ══════════════════════════════════════════════════════════════════════════════
# 2.  PROFILE PLOT HELPERS  (temperature and salinity vs depth)
# ══════════════════════════════════════════════════════════════════════════════

def _add_depth_reference_lines(ax, fjord_row):
    """
    Draw horizontal dashed lines for sill depth and grounding line depth,
    with text labels, if those depths fall within the visible y-range.
    """
    if fjord_row is None:
        return

    ymin, ymax = ax.get_ylim()

    refs = [
        ("sill_depth_m",           "grey", "Sill depth"),
        ("grounding_line_depth_m", "blue", "Grounding line"),
    ]
    for key, colour, label in refs:
        if key in fjord_row and not np.isnan(fjord_row[key]):
            depth_val = fjord_row[key]
            if ymin <= depth_val <= ymax:
                ax.axhline(depth_val, color=colour, linestyle="--", linewidth=1.5)
                ax.text(
                    0.01, depth_val + 10, label,
                    transform=ax.get_yaxis_transform(),
                    fontsize=16, color=colour, ha="left", va="top",
                )


def plot_temp_depth(ax, profiles, fjord_row=None):
    """Temperature vs depth panel with optional sill / grounding-line markers."""
    for key in ["fjord_mouth", "mid_fjord", "near_glacier"]:
        if key in profiles:
            p = profiles[key]
            ax.plot(p["temp"], p["depth"], color=COLORS[key], label=LABELS[key])

    # Add sill and grounding line depths
    ymin, ymax = ax.get_ylim()
    if fjord_row is not None:

        if 'sill_depth_m' in fjord_row and not np.isnan(fjord_row['sill_depth_m']):
            sill = fjord_row['sill_depth_m']
            
            if ymin <= sill <= ymax:   # only plot if visible  
                ax.axhline(sill, color='grey', linestyle='--', linewidth=1.5)

                ax.text(
                    0.01, sill+10,
                    'Sill depth',
                    transform=ax.get_yaxis_transform(),
                    fontsize=16,
                    color='grey',
                    ha='left',
                    va='top'
                )

        if 'grounding_line_depth_m' in fjord_row and not np.isnan(fjord_row['grounding_line_depth_m']):
            gl = fjord_row['grounding_line_depth_m']
            
            if ymin <= gl <= ymax:   # only plot if visible
                ax.axhline(gl, color='blue', linestyle='--', linewidth=1.5)

                ax.text(
                    0.01, gl+10,
                    'Grounding line',
                    transform=ax.get_yaxis_transform(),
                    color='blue',
                    fontsize=16,
                    ha='left',
                    va='top'
                )

    ax.set_xlim(-1, 4.5)
    ax.set_ylim(1000, 0)
    ax.set_xlabel("Conservative Temperature (°C)")
    ax.set_ylabel("Depth (m)")
    _add_depth_reference_lines(ax, fjord_row)


def plot_sal_depth(ax, profiles, fjord_row=None):
    """Salinity vs depth panel with optional sill / grounding-line markers."""
    for key in ["fjord_mouth", "mid_fjord", "near_glacier"]:
        if key in profiles:
            p = profiles[key]
            ax.plot(p["sal"], p["depth"], color=COLORS[key], label=LABELS[key])

    # Add sill and grounding line depths
    ymin, ymax = ax.get_ylim()
    if fjord_row is not None:

        if 'sill_depth_m' in fjord_row and not np.isnan(fjord_row['sill_depth_m']):
            sill = fjord_row['sill_depth_m']
            
            if ymin <= sill <= ymax:   # only plot if visible  
                ax.axhline(sill, color='grey', linestyle='--', linewidth=1.5)

                ax.text(
                    0.01, sill+10,
                    'Sill depth',
                    transform=ax.get_yaxis_transform(),
                    fontsize=16,
                    color='grey',
                    ha='left',
                    va='top'
                )

        if 'grounding_line_depth_m' in fjord_row and not np.isnan(fjord_row['grounding_line_depth_m']):
            gl = fjord_row['grounding_line_depth_m']
            
            if ymin <= gl <= ymax:   # only plot if visible
                ax.axhline(gl, color='blue', linestyle='--', linewidth=1.5)

                ax.text(
                    0.01, gl+10,
                    'Grounding line',
                    transform=ax.get_yaxis_transform(),
                    color='blue',
                    fontsize=16,
                    ha='left',
                    va='top'
                )

    ax.set_xlim(30, 35)
    ax.set_ylim(1000, 0)
    ax.set_xlabel("Absolute Salinity (g/kg)")
    ax.set_ylabel("Depth (m)")
    _add_depth_reference_lines(ax, fjord_row)


# ══════════════════════════════════════════════════════════════════════════════
# 3.  PER-FJORD FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def _locations_from_profiles(profiles):
    ordered = [
        ("near_glacier", "Near-glacier"),
        ("mid_fjord",    "Mid-fjord"),
        ("fjord_mouth",  "Mouth"),
    ]
    return {
        display_name: {"lon": profiles[ptype]["row"]["longitude"],
                       "lat": profiles[ptype]["row"]["latitude"]}
        for ptype, display_name in ordered
        if ptype in profiles
    }


def make_fjord_figure(
    fjord_name,
    fjord_id,
    profiles,
    csv_path,
    bed_cache,
    output_dir,
    ts_xmin=31, ts_xmax=35.5,         # ← hardcoded defaults, not _XMIN/_XMAX
    ts_ymin=-2, ts_ymax=5,    
    panel_labels=("a", "b", "c", "d"),
):
    """
    Build and save a 2×2 figure for one fjord.

    Layout
    ------
    [0,0] Temperature vs depth  |  [0,1] Salinity vs depth
    [1,0] TS diagram            |  [1,1] Bathymetry map

    Parameters
    ----------
    fjord_name : str
        Used as the figure title and in the output filename.
    profiles : dict
        Output of load_fjord_2016().
    csv_path : Path
        CSV file for this fjord (used to pull sill / grounding-line metadata).
    bed_cache : BedMachineCache
        Pre-loaded BedMachine object from load_bedmachine().
    output_dir : Path
        Directory to save the figure.
    panel_labels : tuple of str
        Four single-character labels for panels (a–d), (e–h), etc.
    """
    # ── Metadata ──────────────────────────────────────────────────────────────
    df = pd.read_csv(csv_path)
    df["type"] = df["type"].str.strip().str.lower()
    fjord_meta = df[df["type"] == "near_glacier"].iloc[0]

    # ── Figure and axes ───────────────────────────────────────────────────────
    fig, axes = plt.subplots(
        2, 2, figsize=(16, 14),
        gridspec_kw={"width_ratios": [1, 1]},
    )
    fig.subplots_adjust(hspace=0.30, wspace=0.30,
                        top=0.90, bottom=0.08, left=0.10, right=0.95)
    fig.suptitle(fjord_name, fontsize=24, fontweight="bold", y=0.96)

    # ── [0,0] Temperature vs depth ────────────────────────────────────────────
    plot_temp_depth(axes[0, 0], profiles, fjord_row=fjord_meta)

    # ── [0,1] Salinity vs depth ───────────────────────────────────────────────
    plot_sal_depth(axes[0, 1], profiles, fjord_row=fjord_meta)

    # ── [1,0] TS diagram (from ts_diagram module) ─────────────────────────────
    plot_ts(axes[1, 0], profiles, fjord_row=fjord_meta,
        fjord=fjord_id, water_mass_csv=water_mass_csv, xmin=ts_xmin, xmax=ts_xmax,
        ymin=ts_ymin, ymax=ts_ymax,
        )

    # ── [1,1] Bathymetry map (from bathymetry_map module) ────────────────────
    plot_bathymetry(
        fig,
        axes[1, 1],
        bed_cache,
        locations=_locations_from_profiles(profiles),
        shapefile_defs=DEFAULTS["shapefile_defs"],
        panel_label=f"({panel_labels[3]})", inset_offset=0.0185,  # fraction of figure width; tweak if needed # 0.021 ilu, 0.025 klu, 0.015 kng, 0.012 kng with 25% size,
        coarsen = 4,
    )

    # ── Panel labels for the three matplotlib axes ────────────────────────────
    for ax, lbl in zip([axes[0, 0], axes[0, 1], axes[1, 0]], panel_labels[:3]):
        ax.text(-0.18, 1.08, f"({lbl})", transform=ax.transAxes,
                fontsize=18, fontweight="bold", va="top", ha="left")

    # ── Legends ───────────────────────────────────────────────────────────────
    # Profile panels: show profile-type labels only (exclude sill / GL entries
    # which are added by _add_depth_reference_lines as text, not legend handles).
    for ax in [axes[0, 0], axes[0, 1]]:
        handles, labels = ax.get_legend_handles_labels()
        filtered = [
            (h, l) for h, l in zip(handles, labels)
            if "sill" not in l.lower() and "grounding" not in l.lower()
        ]
        if filtered:
            ax.legend(*zip(*filtered), loc="lower left")

    # ── Save ──────────────────────────────────────────────────────────────────
    safe_name = fjord_name.replace(" ", "_").replace("(", "").replace(")", "")
    out_path  = output_dir / f"{safe_name}_2x2.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved → {out_path}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# 4.  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print("Loading BedMachine …")
    bed_cache = load_bedmachine(DEFAULTS["bedmachine_path"])

    print("Loading QIM profiles …")
    qim_profiles = load_fjord_2016(qim_csv, nc_dir)

    print("Loading KLU profiles …")
    klu_profiles = load_fjord_2016(klu_csv, nc_dir)

    print("Plotting QIM figure …")
    make_fjord_figure(
        fjord_name   = "Qarasaap Imaa (2016)",
        profiles     = qim_profiles,
        csv_path     = qim_csv,
        bed_cache    = bed_cache,
        fjord_id     = "qim",        
        output_dir   = output_dir,
        panel_labels = ("a", "b", "c", "d"),
    )

    print("Plotting KLU figure …")
    make_fjord_figure(
        fjord_name   = "Kangerluluk (2016)",
        ts_xmin      = 30,
        profiles     = klu_profiles,
        csv_path     = klu_csv,
        bed_cache    = bed_cache,
        fjord_id     = "klu",
        output_dir   = output_dir,
        panel_labels = ("a", "b", "c", "d"),
    )

    print("\nAll done.")


if __name__ == "__main__":
    main()