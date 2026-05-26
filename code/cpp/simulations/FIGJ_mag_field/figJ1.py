import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as tck
import simutils as su
from mpl_toolkits.mplot3d import axes3d
from fractions import Fraction
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

mus = [0.5]
"""
azimuth_data = {}
x_data = {}
for mu in mus:
    folder = f"mag_moment_data/c{mu:.1f}".replace(".", "_")
    azimuth, x = su.batch_extract_azimuth(folder)
    azimuth_data[mu] = azimuth
    x_data[mu] = x
    out_name = f"mag_moment_azimuth_c{mu:.1f}".replace(".", "_") + ".npz"
    np.savez(out_name, azimuth=azimuth, x=x)
"""
azimuth_sorted = {}
x_sorted = {}

for mu in mus:
    fname = f"mag_moment_azimuth_c{mu:.1f}".replace(".", "_") + ".npz"
    d = np.load(fname)

    azimuth = d["azimuth"]
    x = d["x"]

    order = np.argsort(azimuth)
    azimuth_sorted[mu] = azimuth[order]
    x_sorted[mu] = x[order]

#NOTE: this formatter is technically a pullback along a reparametrization of the input, 
# so yes, it does not make much sense out of context
def pi_formatter(x, pos):
    frac = Fraction(x / np.pi).limit_denominator(12)
    n, d = frac.numerator, frac.denominator
    if n == 0:
        return r"$0$"
    elif d == 1:
        if n == 1:
            return r"$\frac{\pi}{2}$"
        else:
            return rf"${n}\pi$"
    else:
        if d == 2:
            if x < 0:
                return r"$-\frac{\pi}{2}$"
            else:
                return r"$\frac{\pi}{2}$"
        elif d == 4:
            if x < 0:
                return r"$- \frac{\pi}{4}$"
            else:
                return r"$\frac{\pi}{4}$"
        elif n == 1:
            return rf"$\frac{{\pi}}{{{d}}}$"
        else:
            return rf"$\frac{{{n}\pi}}{{{d}}}$"
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
ax1.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
mu = 0.5
ax1.plot(
    azimuth_sorted[mu],
    x_sorted[mu],
    lw=2.5,
    color="Black",
)
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.set_xlabel(r"$\gamma \, [\text{rad}]$")
ax1.set_ylabel(r"$x \, [\text{m}]$")
ax1.set_xlim(-np.pi / 2,np.pi / 2)
ax1.set_ylim(-0.8, 0.25)
ax1.grid(True, alpha=0.3)
ax1.text(
        0.03, 0.97, "(a)",
        transform=ax1.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
)
data = su.extract_intervals("mag_moment_data", n_intervals=10)
keys = ["c0_5"]
titles = [r"$\mu_0 = 0.5 \, \text{Am}^2$"]
all_azimuth = np.concatenate([data[k]["angle"] for k in keys])
norm = Normalize(vmin=all_azimuth.min(), vmax=all_azimuth.max())
cmap = plt.cm.plasma
panel_labels = ["(b)"]
for key, title, panel in zip(keys, titles, panel_labels):
    d = data[key]
    for azimuth, x, y in zip(d["angle"], d["x"], d["y"]):
        ax2.plot(x, y, color=cmap(norm(azimuth)), lw=1.5)
    ax2.set_title(title)
    ax2.set_xlabel(r"$x \, [\text{m}]$")
    ax2.set_ylabel(r"$y \, [\text{m}]$")
    ax2.grid(True, alpha=0.3)
    ax2.text(
        0.03, 0.97, panel,
        transform=ax2.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
    )

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig1.colorbar(sm, ax=ax2, location="right", pad=0.02)

# New displayed parameter values
tick_pos = np.linspace(0, np.pi, 5)   # same colorbar span, but 5 labels
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
cbar.set_label(r"$\beta \, [\text{rad}]$")
plt.savefig("FIGJ_field.pdf", dpi=300)
