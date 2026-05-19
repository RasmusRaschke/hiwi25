import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import simutils as su
import matplotlib.ticker as mticker
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

datasets = su.extract()
t = datasets["data"].t
x = datasets["data"].x
y = datasets["data"].y
z = datasets["data"].z
mu_x = datasets["data"].mu_x
mu_y = datasets["data"].mu_y
mu_z = datasets["data"].mu_z
mu_norm = np.sqrt(mu_x**2 + mu_y**2 + mu_z**2)
point_x = x + R * mu_x / mu_norm
point_y = y + R * mu_y / mu_norm

datasets = su.extract()
fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
ax.plot(x * 100, y * 100, lw=1.5, color="black")
ax.set_xlabel(r"$x \, [\unit{\centi\metre}]$")
ax.set_ylabel(r"$y \, [\unit{\centi\metre}]$")
ax.grid(True, alpha=0.3)
ax.xaxis.set_inverted(True)
dx_pt = 100 * (point_x - x)
dy_pt = 100 * (point_y - y)
axins = inset_axes(ax, width="30%", height="30%", loc="lower right", borderpad=2)
axins.plot(dx_pt, dy_pt, lw=1.2, color="tab:orange")
axins.set_xlim(-0.01,0.11)
axins.set_ylim(-0.55,0.55)
axins.grid(True, alpha=0.3)
"""
    ax.text(
        0.03, 0.97, panel,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
    )
"""
plt.savefig("FIGK1_ramp_traj.pdf", dpi=300)