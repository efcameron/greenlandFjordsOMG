"""
bathymetry_map.py
-----------------
Self-contained module for plotting fjord bathymetry maps with a Greenland inset.

Public API
----------
load_bedmachine(path)  →  BedMachineCache
    Load a BedMachine NetCDF file once; reuse the returned object across many plots.

plot_bathymetry(fig, ax, cache, locations, shapefile_defs, ...)  →  GeoAxes
    Replace a placeholder Axes with a polar-stereo bathymetry map.

Example
-------
    from bathymetry_map import load_bedmachine, plot_bathymetry, DEFAULTS

    cache = load_bedmachine(DEFAULTS["bedmachine_path"])

    fig, axes = plt.subplots(2, 2, figsize=(16, 20))

    ax_map = plot_bathymetry(
        fig,
        axes[1, 1],
        cache,
        locations={
            "Mouth":       {"lon": -45.0, "lat": 70.0},
            "Near-glacier":{"lon": -44.5, "lat": 70.3},
            "Mid-fjord":   {"lon": -44.7, "lat": 70.15},
        },
        shapefile_defs=DEFAULTS["shapefile_defs"],
        panel_label="(d)",
    )

    plt.savefig("output.png", dpi=150, bbox_inches="tight")
"""

from __future__ import annotations

import numpy as np
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pyproj import CRS, Transformer
from shapely.geometry import box

# ── Colour / marker style for known location types ────────────────────────────
# Keys are matched case-insensitively against the 'locations' dict keys.
_LOCATION_STYLES: dict[str, dict] = {
    "mouth":        {"color": "black",  "label": "Mouth"},
    "near-glacier": {"color": "red",    "label": "Near glacier"},
    "near_glacier": {"color": "red",    "label": "Near glacier"},
    "mid-fjord":    {"color": "orange", "label": "Mid-fjord"},
    "mid_fjord":    {"color": "orange", "label": "Mid-fjord"},
}
_FALLBACK_COLOR = "cyan"

# ── Default hardcoded paths (override by passing your own to plot_bathymetry) ─
_BASE = (
    r"C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD"
    r"\First_Year\Projects\Fjord_Shelf_Observations\Working_data"
)

DEFAULTS: dict = {
    "bedmachine_path": (
        r"C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD"
        r"\First_Year\Projects\Fjord_Shelf_Observations\Working_data"
        r"\Bathymetry\BedMachineGreenland-v5.nc"
    ),
    # Each entry: (path, label, colour, zorder, fill)
    "shapefile_defs": [
        (
            _BASE + r"\GEEDiT_terminus_lines\SE\GEEDiT_termini_SE_all.gpkg",
            "Terminus Line", "blue", 4, False,
        ),
        (
            _BASE + r"\Gerrish_2020_Greenland_coast\Gerrish_2020_Greenland_coast.gpkg",
            "Land", "grey", 2, True,
        ),
        (
            _BASE + r"\PROMICE2022IceMask\02-PROMICE-2022-IceMask-polygon.gpkg",
            "Ice Sheet", "white", 3, True,
        ),
    ],
}


# ══════════════════════════════════════════════════════════════════════════════
# BedMachine cache
# ══════════════════════════════════════════════════════════════════════════════

class BedMachineCache:
    """
    Loads BedMachine once and serves coarsened spatial subsets on demand.

    Parameters
    ----------
    path : str
        Path to BedMachineGreenland NetCDF file.
    """

    def __init__(self, path: str):
        self.ds = xr.open_dataset(
            path,
            decode_coords="all",
            drop_variables=["errbed", "source"],
        )
        self.bed = self.ds["bed"]

    def get_subset(
        self,
        x_min: float, x_max: float,
        y_min: float, y_max: float,
        coarsen: int = 10,
    ) -> xr.DataArray:
        """
        Return a coarsened tile of the bed grid.

        Parameters
        ----------
        x_min, x_max, y_min, y_max : float
            Bounding box in EPSG:3413 (polar stereo) metres.
        coarsen : int
            Spatial coarsening factor (higher = faster but lower resolution).
        """
        return (
            self.bed
            .sel(x=slice(x_min, x_max), y=slice(y_max, y_min))
            .coarsen(x=coarsen, y=coarsen, boundary="trim")
            .mean()
            .compute()
        )


def load_bedmachine(path: str) -> BedMachineCache:
    """
    Load BedMachine from *path* and return a reusable cache object.

    Intended to be called once at the top of a script and passed around.
    """
    return BedMachineCache(path)


# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

_PS3413 = CRS.from_epsg(3413)
_WGS84  = CRS.from_epsg(4326)


def _latlon_to_ps(lon: float, lat: float) -> tuple[float, float]:
    """Convert a single lon/lat pair to EPSG:3413 (x, y) in metres."""
    tf = Transformer.from_crs(_WGS84, _PS3413, always_xy=True)
    return tf.transform(lon, lat)


def _style_for_key(key: str) -> dict:
    """Look up marker style by location key (case-insensitive)."""
    return _LOCATION_STYLES.get(key.lower(), {"color": _FALLBACK_COLOR, "label": key})


def _draw_shapefiles(ax, shapefile_defs, extent_box, simplify_m: float = 500):
    """Clip, simplify, and render each shapefile onto *ax*."""
    for path, label, colour, zorder, fill in shapefile_defs:
        gdf = gpd.read_file(path).to_crs(epsg=3413)
        gdf = gdf[gdf.intersects(extent_box)]
        if gdf.empty:
            continue
        gdf["geometry"] = gdf["geometry"].simplify(simplify_m, preserve_topology=True)
        ax.add_geometries(
            gdf.geometry,
            crs=ccrs.NorthPolarStereo(),
            facecolor=colour if fill else "none",
            edgecolor=colour,
            linewidth=1,
            zorder=zorder,
            label=label,
            rasterized=True,
        )


def _add_scale_bar(ax, x_min, x_max, y_min, length_m: float = 20_000):
    """Draw a simple scale bar in the lower-left corner of *ax*."""
    s0 = (
        x_min + 0.05 * (x_max - x_min),
        y_min + 0.05 * (x_max - x_min),   # use same scale for a square offset
    )
    ax.add_artist(Line2D(
        [s0[0], s0[0] + length_m], [s0[1], s0[1]],
        color="black", linewidth=3,
        transform=ccrs.NorthPolarStereo(), zorder=20,
    ))
    ax.text(
        s0[0] + length_m / 2, s0[1] + 0.02,
        f"{int(length_m / 1000)} km",
        ha="center", va="bottom", fontsize=14, color="black",
        transform=ccrs.NorthPolarStereo(), zorder=20,
    )


def _add_greenland_inset(fig, ax_main, x_ng, y_ng, shapefile_defs, inset_offset=0.0185):
    """
    Add a small Greenland overview inset to the upper-left of *ax_main*.

    Only the coastline and ice-mask shapefiles are drawn (identified by
    'coast' or 'icemask' in their path strings, case-insensitive).
    """
    
    bpos   = ax_main.get_position()
    iw, ih = 0.30 * bpos.width, 0.30 * bpos.height
    ax_i   = fig.add_axes(
        [bpos.x0 - inset_offset, bpos.y1 - ih, iw, ih],
        projection=ccrs.NorthPolarStereo(),
    )

    # Determine Greenland extent from the coastline layer
    coast_paths = [p for p, *_ in shapefile_defs if "coast" in p.lower()]
    if coast_paths:
        gdf_land = gpd.read_file(coast_paths[0]).to_crs(epsg=3413)
        minx, miny, maxx, maxy = gdf_land.total_bounds
        pad = 200e3
        ax_i.set_extent(
            [minx - pad, maxx + pad, miny - pad, maxy + pad],
            crs=ccrs.NorthPolarStereo(),
        )

    for path, label, colour, zorder, fill in shapefile_defs:
        if "coast" in path.lower() or "icemask" in path.lower():
            gdf = gpd.read_file(path).to_crs(epsg=3413)
            gdf["geometry"] = gdf["geometry"].simplify(500, preserve_topology=True)
            ax_i.add_geometries(
                gdf.geometry,
                crs=ccrs.NorthPolarStereo(),
                facecolor=colour if fill else "none",
                edgecolor=colour,
                linewidth=0.5,
                zorder=zorder,
                rasterized=True,
            )

    ax_i.scatter(
        x_ng, y_ng, color="red", marker="*", s=80, # star on inset map
        transform=ccrs.NorthPolarStereo(), zorder=10,
    )

    for spine in ax_i.spines.values():
        spine.set_visible(True)
        spine.set_color("grey")
        spine.set_linewidth(0.5)

    return ax_i


# ══════════════════════════════════════════════════════════════════════════════
# Public plotting function
# ══════════════════════════════════════════════════════════════════════════════

def plot_bathymetry(
    fig: plt.Figure,
    ax: plt.Axes,
    cache: BedMachineCache,
    locations: dict[str, dict],
    shapefile_defs: list[tuple] | None = None,
    *,
    buffer_m: float = 50_000,
    coarsen: int = 10,
    scale_bar_m: float = 20_000,
    panel_label: str = "",
    inset: bool = True,
    inset_key: str = "near-glacier",
    cbar_label: str = "Depth (m)",
    simplify_m: float = 500,
    label_fontsize: int = 20,
    inset_offset=0.0185
) -> plt.Axes:
    """
    Replace placeholder *ax* with a polar-stereo bathymetry map.

    Parameters
    ----------
    fig : matplotlib Figure
        The parent figure.
    ax : matplotlib Axes
        Placeholder axes whose position is reused for the map.
    cache : BedMachineCache
        Pre-loaded BedMachine object (from ``load_bedmachine()``).
    locations : dict
        Mapping of location name → {"lon": float, "lat": float}.
        Known names ("Mouth", "Near-glacier" / "Near_glacier",
        "Mid-fjord" / "Mid_fjord") get predefined colours; any other
        name gets a cyan marker.
        Example::

            {
                "Mouth":        {"lon": -45.0, "lat": 70.0},
                "Near-glacier": {"lon": -44.5, "lat": 70.3},
                "Mid-fjord":    {"lon": -44.7, "lat": 70.15},
            }

    shapefile_defs : list of (path, label, colour, zorder, fill) tuples
        If None, falls back to ``DEFAULTS["shapefile_defs"]``.
    buffer_m : float
        Padding around the location bounding box in metres (default 50 km).
    coarsen : int
        BedMachine coarsening factor (default 10).
    scale_bar_m : float
        Scale bar length in metres (default 20 km).
    panel_label : str
        Panel label drawn in the upper-left corner, e.g. ``"(d)"``.
    inset : bool
        Whether to draw the Greenland overview inset (default True).
    inset_key : str
        Which location key to mark on the inset (default "near-glacier").
    cbar_label : str
        Colourbar label (default "Depth (m)").
    simplify_m : float
        Shapefile geometry simplification tolerance in metres (default 500).

    Returns
    -------
    ax_map : GeoAxes
        The new polar-stereo axes (replaces the placeholder *ax*).
    """
    if shapefile_defs is None:
        shapefile_defs = DEFAULTS["shapefile_defs"]

    # ── Convert all locations to polar stereo ─────────────────────────────────
    ps_coords: dict[str, tuple[float, float]] = {
        key: _latlon_to_ps(loc["lon"], loc["lat"])
        for key, loc in locations.items()
    }

    xs = [c[0] for c in ps_coords.values()]
    ys = [c[1] for c in ps_coords.values()]
    x_min, x_max = min(xs) - buffer_m, max(xs) + buffer_m
    y_min, y_max = min(ys) - buffer_m, max(ys) + buffer_m

    # ── Fetch BedMachine tile ─────────────────────────────────────────────────
    bed_subset = cache.get_subset(x_min, x_max, y_min, y_max, coarsen=coarsen)

    # ── Replace placeholder axes with a GeoAxes ───────────────────────────────
    pos   = ax.get_position()
    fig.delaxes(ax)
    ax_map = fig.add_axes(pos.bounds, projection=ccrs.NorthPolarStereo())

    if panel_label:
        ax_map.text(
            -0.15, 0.98, panel_label,
            transform=ax_map.transAxes,
            fontsize=label_fontsize, fontweight="bold", va="top", ha="left",
        )

    # ── Bathymetry pcolormesh ─────────────────────────────────────────────────
    bathy = ax_map.pcolormesh(
        bed_subset.x, bed_subset.y, bed_subset,
        cmap="terrain", shading="auto",
        transform=ccrs.NorthPolarStereo(),
        zorder=0,
        rasterized=True,
    )
    cbar = fig.colorbar(bathy, ax=ax_map, orientation="vertical",
                        fraction=0.02, pad=0.04)
    cbar.set_label(cbar_label, fontsize=label_fontsize)

    ax_map.set_extent([x_min, x_max, y_min, y_max], crs=ccrs.NorthPolarStereo())

    # ── Scale bar ─────────────────────────────────────────────────────────────
    _add_scale_bar(ax_map, x_min, x_max, y_min, length_m=scale_bar_m)

    # ── CTD scatter markers ───────────────────────────────────────────────────
    for key, (x, y) in ps_coords.items():
        style = _style_for_key(key)
        ax_map.scatter(
            x, y,
            color=style["color"], marker="*", s=350, # stars on main map
            transform=ccrs.NorthPolarStereo(),
            zorder=10, label=style["label"],
        )

    ax_map.legend(loc="lower right")

    # ── Shapefiles ────────────────────────────────────────────────────────────
    extent_box = box(x_min, y_min, x_max, y_max)
    _draw_shapefiles(ax_map, shapefile_defs, extent_box, simplify_m=simplify_m)

    # ── Greenland inset ───────────────────────────────────────────────────────
    if inset:
        inset_key_lower = inset_key.lower()
        inset_xy = next(
            (xy for k, xy in ps_coords.items() if k.lower() == inset_key_lower),
            next(iter(ps_coords.values())),   # fall back to first location
        )
        _add_greenland_inset(fig, ax_map, *inset_xy, shapefile_defs, inset_offset=inset_offset)

    return ax_map