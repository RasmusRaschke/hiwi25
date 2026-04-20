import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import simutils as su
import matplotlib.ticker as mticker
import re
from collections import defaultdict
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "font.size": 16,        # general default
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
})

g = 9.80665
B = 5.0e-5
R = 0.005
M = 0.004

data = su.extract_intervals("mag_moment_data", n_intervals=10)
keys = ["c0_5", "c1_0", "c2_0"]
titles = [r"$\mu = 0.5 \, \text{Am}^2$", r"$\mu = 1 \, \text{Am}^2$", r"$\mu = 2 \, \text{Am}^2$"]
all_azimuth = np.concatenate([data[k]["angle"] for k in keys])
norm = Normalize(vmin=all_azimuth.min(), vmax=all_azimuth.max())
cmap = plt.cm.plasma
fig, axs = plt.subplots(1, 3, figsize=(15, 5), layout="constrained")

for ax, key, title in zip(axs, keys, titles):
    d = data[key]
    for azimuth, x, y in zip(d["angle"], d["x"], d["y"]):
        ax.plot(x * 100, y * 100, color=cmap(norm(azimuth)), lw=1.5)
    ax.set_title(title)
    ax.set_xlabel(r"$x \, [\text{cm}]$")
    ax.set_ylabel(r"$y \, [\text{cm}]$")
    ax.grid(True, alpha=0.3)

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs, location="right", pad=0.02)
ticks = np.linspace(0, np.pi, 6)
cbar.set_ticks(ticks)
def pi_fraction_formatter(x, pos):
    n = int(round(x / np.pi * 5))
    if n == 0:
        return r"$0$"
    if n == 5:
        return r"$\pi$"
    return rf"$\frac{{{n}\pi}}{{5}}$"
cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(pi_fraction_formatter))
cbar.set_label(r"$\varphi \, [\text{rad}]$")
plt.savefig("FIGF2_mag_mom_traj.pdf", dpi=300)