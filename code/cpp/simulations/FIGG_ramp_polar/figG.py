from fractions import Fraction

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np
import simutils as su
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import axes3d

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


g = 9.80665
B = 5.0e-5
R = 0.005
M = 0.004

# polar, y = su.batch_extract_polar("mag_moment_data1/c1_0")
#polar2, y2 = su.batch_extract_polar2("mag_moment_data2/c1_0")
#np.savez("mag_moment_polar2.npz", polar=polar2, y=y2)
data = np.load("mag_moment_polar1.npz")
data2 = np.load("mag_moment_polar2.npz")
angle = data["polar"]
y = data["y"]
angle2 = data2["polar"]
z = data2["y"]

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
ax1.plot(angle[order], y[order] * 100, lw=3.5, color="black")
order = np.argsort(angle2)
ax1.plot(angle2[order], z[order] * 100, lw=1.5, color="gray", ls="--")
x_points = np.array([1.178, 0.118, -1.178])
idx = np.abs(angle[:, None] - x_points).argmin(axis=0)
y_points = y[idx]
ax1.plot(x_points[0], y_points[0] * 100, "o", color="blue", ms=10)
ax1.plot(x_points[1], y_points[1] * 100, "o", color="orange", ms=10)
ax1.plot(x_points[2], y_points[2] * 100, "o", color="red", ms=10)
ax1.set_xlabel(r"$\psi \, [\unit{rad}]$")
ax1.set_ylabel(r"$y|_{t=1 \, \unit{s}} \, [\unit{\centi\metre}]$")
ax1.set_xlim(-np.pi, np.pi)
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.grid(True, alpha=0.3)
ax2.plot(t, y1 * 100, label=r"$\psi = 1.2 \, \unit{rad}$", lw=3.5, color="blue")
ax2.plot(t, y2 * 100, label=r"$\psi = 0.1 \, \unit{rad}$", lw=3.5, color="orange")
ax2.plot(t, y3 * 100, label=r"$\psi = -1.2 \, \unit{rad}$", lw=3.5, color="red")
ax2.set_xlabel(r"$t \, [\unit{\second}]$")
ax2.set_ylabel(r"$y \, [\unit{\centi\metre}]$")
ax2.grid(True, alpha=0.3)
ax2.legend(loc="lower left")
ax2.set_xlim(0, 1)
ax2.set_ylim(-17, 1.80)
ax1.text(
    0.03,
    0.97,
    "(a)",
    transform=ax1.transAxes,
    ha="left",
    va="top",
    fontsize=16,
    fontweight="bold",
)
ax2.text(
    0.03,
    0.97,
    "(b)",
    transform=ax2.transAxes,
    ha="left",
    va="top",
    fontsize=16,
    fontweight="bold",
)
plt.savefig("FIGG_polar.pdf", dpi=300)
