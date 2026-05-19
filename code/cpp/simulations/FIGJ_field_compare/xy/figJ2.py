import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import simutils as su
import matplotlib.ticker as mticker
from fractions import Fraction
import re
from collections import defaultdict
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{siunitx}",
    "font.size": 20,          # base size for TeX-rendered text
    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
})

g = 9.80665
B = 5.0e-5
R = 0.005
M = 0.004

data = su.extract_intervals("mag_moment_data", n_intervals=10)
keys = ["c0_5", "c1_0", "c2_0"]
titles = [r"$B = 0.5 \, \unit{\tesla}$", r"$B = 1 \, \unit{\tesla}$", r"$B = 2 \, \unit{\tesla}$"]
all_azimuth = np.concatenate([data[k]["angle"] for k in keys])
norm = Normalize(vmin=all_azimuth.min(), vmax=all_azimuth.max())
cmap = plt.cm.plasma
fig, axs = plt.subplots(1, 3, figsize=(15, 5), layout="constrained")
panel_labels = ["(a)", "(b)", "(c)"]
for key in keys:
    ang = np.array(data[key]["angle"])
    print(key, ang.min(), ang.max(), np.sum(ang > 0), np.sum(ang < 0))
for ax, key, title, panel in zip(axs, keys, titles, panel_labels):
    d = data[key]
    for azimuth, x, y in zip(d["angle"], d["x"], d["y"]):
        ax.plot(x * 100, y * 100, color=cmap(norm(azimuth)), lw=1.5)
    ax.set_title(title)
    ax.set_xlabel(r"$x \, [\unit{\centi\metre}]$")
    ax.set_ylabel(r"$y \, [\unit{\centi\metre}]$")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(15))
    ax.set_ylim(-70.0, 2.0)
    ax.set_xlim(-30.0, 30.0)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_inverted(True)
    ax.text(
        0.03, 0.97, panel,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
    )

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs, location="right", pad=0.02)

# New displayed parameter values
tick_pos = np.linspace(-np.pi /2, np.pi /2, 5)   # same colorbar span, but 5 labels
tick_lab = [
    r"$-\frac{\pi}{2}$",
    r"$-\frac{\pi}{4}$",
    r"$0$",
    r"$\frac{\pi}{4}$",
    r"$\frac{\pi}{2}$",
]

cbar.set_ticks(tick_pos)
cbar.set_ticklabels(tick_lab)
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r"$\overline{\theta} \, [\unit{rad}]$")
plt.savefig("FIGJ2_field_traj.pdf", dpi=300)