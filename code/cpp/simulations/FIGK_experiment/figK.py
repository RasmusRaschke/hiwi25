import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
import simutils as su

plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{siunitx} \usepackage{bm}",
    "font.size": 20,
    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
})

R = 0.005  # sphere radius

def colored_line(ax, x, y, c, cmap, norm, lw=2.0, alpha=1.0, zorder=2):
    x = np.asarray(x)
    y = np.asarray(y)
    c = np.asarray(c)

    pts = np.column_stack([x, y]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    c_seg = 0.5 * (c[:-1] + c[1:])

    lc = LineCollection(
        segs, cmap=cmap, norm=norm,
        linewidth=lw, alpha=alpha, zorder=zorder
    )
    lc.set_array(c_seg)
    ax.add_collection(lc)
    return lc

datasets = su.extract()
names = ["hamburg", "jakata", "tokyo"]

# Marker style per pair
markers = {
    "hamburg": "o",   # filled circle
    "jakata": "s",    # filled square
    "tokyo": "^",     # filled triangle
}

all_O = []
all_mu = []
all_xy = []

for name in names:
    d = datasets[name]

    x = np.asarray(d.x)
    y = np.asarray(d.y)

    mu_x = np.asarray(d.mu_x)
    mu_y = np.asarray(d.mu_y)
    mu_z = np.asarray(d.mu_z)

    mu_norm = np.sqrt(mu_x**2 + mu_y**2 + mu_z**2)
    mu_norm_safe = np.where(mu_norm == 0, np.nan, mu_norm)

    point_x = x + R * mu_x / mu_norm_safe
    point_y = y + R * mu_y / mu_norm_safe

    if hasattr(d, "O_x"):
        Ox, Oy, Oz = d.O_x, d.O_y, d.O_z
    else:
        Ox, Oy, Oz = d.Ox, d.Oy, d.Oz

    O_norm = np.sqrt(Ox**2 + Oy**2 + Oz**2)

    # Keep only 1.0-scale fluctuations for the mu-coloring
    mu_norm_plot = np.round(mu_norm, 0)

    all_O.append(O_norm[np.isfinite(O_norm)])
    all_mu.append(mu_norm_plot[np.isfinite(mu_norm_plot)])

    all_xy.append(np.column_stack([x * 100, y * 100]))
    all_xy.append(np.column_stack([point_x * 100, point_y * 100]))

all_O = np.concatenate(all_O)
all_mu = np.concatenate(all_mu)

norm_O = Normalize(vmin=np.nanmin(all_O), vmax=np.nanmax(all_O))
norm_mu = Normalize(vmin=np.nanmin(all_mu), vmax=np.nanmax(all_mu))

cmap_O = plt.cm.viridis
cmap_mu = plt.cm.plasma

fig, ax = plt.subplots(figsize=(7, 5))
fig.subplots_adjust(right=0.83)

for name in names:
    d = datasets[name]

    x = np.asarray(d.x)
    y = np.asarray(d.y)

    mu_x = np.asarray(d.mu_x)
    mu_y = np.asarray(d.mu_y)
    mu_z = np.asarray(d.mu_z)

    mu_norm = np.sqrt(mu_x**2 + mu_y**2 + mu_z**2)
    mu_norm_safe = np.where(mu_norm == 0, np.nan, mu_norm)
    mu_norm_plot = np.round(mu_norm, 0)

    point_x = x + R * mu_x / mu_norm_safe
    point_y = y + R * mu_y / mu_norm_safe

    if hasattr(d, "O_x"):
        Ox, Oy, Oz = d.O_x, d.O_y, d.O_z
    else:
        Ox, Oy, Oz = d.Ox, d.Oy, d.Oz

    O_norm = np.sqrt(Ox**2 + Oy**2 + Oz**2)

    # COM trajectory
    colored_line(ax, x * 100, y * 100, O_norm, cmap_O, norm_O, lw=2.0)

    # mu-tip trajectory
    colored_line(ax, point_x * 100, point_y * 100, mu_norm_plot, cmap_mu, norm_mu, lw=2.0)

    m = markers[name]

    # Start/end points for COM trajectory
    ax.scatter(
        x[0] * 100, y[0] * 100,
        marker=m, s=80, facecolors="white", edgecolors="black", linewidths=1.2, zorder=5
    )
    ax.scatter(
        x[-1] * 100, y[-1] * 100,
        marker=m, s=80, facecolors="black", edgecolors="black", linewidths=1.2, zorder=5
    )

    # Start/end points for mu-tip trajectory
    ax.scatter(
        point_x[0] * 100, point_y[0] * 100,
        marker=m, s=80, facecolors="white", edgecolors="black", linewidths=1.2, zorder=5
    )
    ax.scatter(
        point_x[-1] * 100, point_y[-1] * 100,
        marker=m, s=80, facecolors="black", edgecolors="black", linewidths=1.2, zorder=5
    )

from matplotlib.lines import Line2D

shape_legend = [
    Line2D([0], [0], marker='o', linestyle='None',
           markerfacecolor='black', markeredgecolor='black',
           markersize=10, label='Hamburg'),
    Line2D([0], [0], marker='s', linestyle='None',
           markerfacecolor='black', markeredgecolor='black',
           markersize=10, label='Jakata'),
    Line2D([0], [0], marker='^', linestyle='None',
           markerfacecolor='black', markeredgecolor='black',
           markersize=10, label='Tokyo'),
]

ax.legend(handles=shape_legend, loc="best", frameon=True)
all_xy = np.concatenate(all_xy, axis=0)
xmin, ymin = np.nanmin(all_xy, axis=0)
xmax, ymax = np.nanmax(all_xy, axis=0)

pad_x = 0.05 * (xmax - xmin)
pad_y = 0.05 * (ymax - ymin)

ax.set_xlim(xmax + pad_x, xmin - pad_x)
ax.set_ylim(ymin - pad_y, ymax + pad_y)

ax.set_xlabel(r"$x \, [\unit{\centi\metre}]$")
ax.set_ylabel(r"$y \, [\unit{\centi\metre}]$")
ax.grid(True, alpha=0.3)

sm_O = ScalarMappable(norm=norm_O, cmap=cmap_O)
sm_O.set_array([])
sm_mu = ScalarMappable(norm=norm_mu, cmap=cmap_mu)
sm_mu.set_array([])

cax1 = fig.add_axes([0.86, 0.56, 0.025, 0.30])
cbar1 = fig.colorbar(sm_O, cax=cax1)
cbar1.set_label(r"$\| \bm{\Omega} \| \, [\unit{\second^{-1}}]$")
from matplotlib.ticker import MaxNLocator
cax2 = fig.add_axes([0.86, 0.14, 0.025, 0.30])
cbar2 = fig.colorbar(sm_mu, cax=cax2)
cbar2.set_label(r"$\|\bm{\mu}\| \, [\unit{\ampere\second\squared}]$")
cbar1.set_ticks([0, 50, 100, 150, 200])
cbar2.set_ticks([0.9, 1.0, 1.1])

plt.savefig("FIGK.pdf", dpi=300, bbox_inches="tight")