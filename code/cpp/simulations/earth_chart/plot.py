from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from matplotlib.ticker import FuncFormatter
plt.style.use("seaborn-v0_8-paper")
plt.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{siunitx}",
        "font.size": 22,  # base size for TeX-rendered text
        "axes.titlesize": 22,
        "axes.labelsize": 22,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 20,
    }
)

##############################################################################
# PARAMETERS
INPUT_FILE = Path("results/earth_results.npz")
OUTPUT_FILE = Path("results/earth_map.pdf")
MARKERS = [
    (53.550556, 9.993682),
    (-6.200000, 106.826944),
    (35.689444, 139.691667),
    (-35.279722, 149.128998),
]
FIGSIZE = (14, 7)
DPI = 300
COLORMAP = "viridis"
VALUE_NAME = "x"      # or "t"
##############################################################################


def load_data(filename):
    data = np.load(filename)
    lat = data["lat"]
    lon = data["lon"]
    values = data[VALUE_NAME]
    return lat, lon, values


def create_figure():
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes([0.05, 0.1, 0.55, 0.6], projection=ccrs.Robinson())
    ax.set_global()
    return fig, ax


def draw_map(ax):
    ax.add_feature(cfeature.LAND, zorder=0)
    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3)
    gl = ax.gridlines(
        draw_labels=False,
        linewidth=0.3,
        alpha=0.5,
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    gl.xlabel_style = {"size": 18}
    gl.ylabel_style = {"size": 18}
    gl.xlocator = cticker.LongitudeLocator()
    gl.ylocator = cticker.LatitudeLocator()
    gl.xformatter = FuncFormatter(lambda x, pos: rf"${x:.0f}^\circ$")
    gl.yformatter = FuncFormatter(lambda y, pos: rf"${y:.0f}^\circ$")


def draw_markers(ax, markers):
    styles = ["s", "o", "^", "v"]
    for (lat, lon), m in zip(markers, styles):
        ax.scatter(
            lon,
            lat,
            transform=ccrs.PlateCarree(),
            marker=m,
            s=100,
            edgecolor="black",
            linewidth=1.0,
            zorder=10,
            color="red",
        )

def draw_data(ax, lat, lon, values):
    mesh = ax.pcolormesh(
        lon,
        lat,
        values,
        transform=ccrs.PlateCarree(),
        shading="auto",
        cmap=COLORMAP,
    )
    return mesh


def main():
    lat, lon, values = load_data(INPUT_FILE)
    fig, ax = create_figure()
    draw_map(ax)
    mesh = draw_data(
        ax,
        lon=lon,
        lat=lat,
        values=values,
    )
    draw_markers(ax, MARKERS)
    cbar = plt.colorbar(
        mesh,
        ax=ax,
        orientation="vertical",
        pad=0.03,
        shrink=0.75,
    )
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"${100*x:.1f}$"))
    cbar.set_label(r"$x|_{t=0.5 \, \unit{s}} \, [\unit{cm}]$")
    cbar.set_ticks([-0.08, -0.04, 0.0, 0.08, 0.04])
    plt.savefig(
        OUTPUT_FILE,
        dpi=DPI,
        bbox_inches="tight",
    )
    plt.close(fig)
    print(f"Saved {OUTPUT_FILE}")


if __name__ == "__main__":
    main()