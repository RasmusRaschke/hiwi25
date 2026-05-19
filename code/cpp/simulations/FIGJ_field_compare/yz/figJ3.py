import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import simutils as su
import matplotlib.ticker as tck
import re
from fractions import Fraction
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
"""
polar, y = su.batch_extract_polar("mag_moment_data/c1_0")
np.savez("mag_moment_polar.npz", polar=polar, y=y)
"""
data = np.load("mag_moment_polar.npz")
angle = data["polar"]
y = data["y"]

datasets = su.extract()
t = datasets["sample1"].t
y1 = datasets["sample1"].y
y2 = datasets["sample2"].y
y3 = datasets["sample3"].y

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
def pi_formatter(x, pos):
    frac = Fraction(x / np.pi).limit_denominator(12)

    sign = "-" if frac < 0 else ""
    frac = abs(frac)

    n, d = frac.numerator, frac.denominator

    if n == 0:
        return "0"

    # integer multiples of pi
    if d == 1:
        if n == 1:
            return rf"${sign}\pi$"
        else:
            return rf"${sign}{n}\pi$"

    # fractional multiples of pi
    if n == 1:
        return rf"${sign}\frac{{\pi}}{{{d}}}$"
    else:
        return rf"${sign}\frac{{{n}\pi}}{{{d}}}$"
ax1.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
order = np.argsort(angle)
ax1.plot(angle[order], y[order]*100, lw=3.5, color="black")
x_points = np.array([1.178, -0.141, -1.178])
idx = np.abs(angle[:, None] - x_points).argmin(axis=0)
y_points = y[idx]
ax1.plot(x_points[0], y_points[0] * 100, "o", color="blue", ms=10)
ax1.plot(x_points[1], y_points[1] * 100, "o", color="orange", ms=10)
ax1.plot(x_points[2], y_points[2] * 100, "o", color="red", ms=10)
ax1.set_xlabel(r"$\overline{\psi} \, [\unit{rad}]$")
ax1.set_ylabel(r"$y|_{t=1 \, \unit{s}} \, [\unit{\centi\metre}]$")
ax1.set_xlim(-np.pi,np.pi)
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.grid(True, alpha=0.3)
ax2.plot(t, y1*100, label=r"$\overline{\psi} = 1.18 \, \unit{rad}$", lw=3.5, color="blue")
ax2.plot(t, y2*100, label=r"$\overline{\psi} = -0.14 \, \unit{rad}$", lw=3.5, color="orange")
ax2.plot(t, y3*100, label=r"$\overline{\psi} = -1.18 \, \unit{rad}$", lw=3.5, color="red")
ax2.set_xlabel(r"$t \, [\unit{\second}]$")
ax2.set_ylabel(r"$y \, [\unit{\centi\metre}]$")
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlim(0,1)
ax2.set_ylim(-10,1.5)
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
plt.savefig("FIGJ3_field_yz.pdf", dpi=300)
