import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import simutils as su
import matplotlib.ticker as tck
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
mus = np.round(np.arange(0.0, 2.5 + 0.1, 0.1), 1)
"""

# All folders from c0_0 to c2_5 in 0.1 steps
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


def pi_formatter(x, pos=None):
    frac = Fraction(x / np.pi).limit_denominator(12)
    n, d = frac.numerator, frac.denominator

    if n == 0:
        return r"$0$"

    sign = "-" if n < 0 else ""
    n = abs(n)

    if d == 1:
        if n == 1:
            return rf"${sign}\pi$"
        return rf"${sign}{n}\pi$"

    if n == 1:
        return rf"${sign}\frac{{\pi}}{{{d}}}$"
    return rf"${sign}\frac{{{n}\pi}}{{{d}}}$"
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
ax1.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
ax2.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
plot_mus = np.round(np.arange(0.0, 2.5 + 0.1, 0.1), 1)
norm = Normalize(vmin=min(plot_mus), vmax=max(plot_mus))
plot_mus_1 = np.round(np.arange(0.0, 0.9 + 0.1, 0.1), 1)
plot_mus_2 = np.round(np.arange(1.0, 2.5 + 0.1, 0.1), 1)
cmap = plt.cm.coolwarm
for i, mu in enumerate(plot_mus_1):
    ax1.plot(
        azimuth_sorted[mu],
        x_sorted[mu] * 100,
        lw=2.5,
        color=cmap(norm(mu)),
    )
for i, mu in enumerate(plot_mus_2):
    ax2.plot(
        azimuth_sorted[mu],
        x_sorted[mu] * 100,
        lw=2.5,
        color=cmap(norm(mu)),
    )
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax2)
cbar.set_label(r"$\mu_0\,[\unit{\ampere\metre\squared}]$")
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.set_xlabel(r"$\overline{\theta} \, [\unit{rad}]$")
ax1.set_ylabel(r"$x|_{t=1 \, \unit{s}} \, [\unit{\centi\metre}]$")
ax1.grid(True, alpha=0.3)
ax2.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax2.set_xlabel(r"$\overline{\theta} \, [\unit{rad}]$")
ax2.set_ylabel(r"$x|_{t=1 \, \unit{s}} \, [\unit{\centi\metre}]$")
ax2.grid(True, alpha=0.3)
ax1.set_xlim(-np.pi /2, np.pi /2)
ax2.set_xlim(-np.pi /2, np.pi /2)
ax1.text(
        0.03, 0.97, "(a)",
        transform=ax1.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
)
ax2.text(
        0.03, 0.97, "(b)",
        transform=ax2.transAxes,
        ha="left", va="top",
        fontsize=16, fontweight="bold"
)
#ax1.legend()
plt.savefig("FIGJ1_field_azimuth.pdf", dpi=300)
#"""
