import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as tck
from fractions import Fraction
import simutils as su
from mpl_toolkits.mplot3d import axes3d
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
    n, d = frac.numerator, frac.denominator
    if n == 0:
        return "0"
    elif d == 1:
        if n == 1:
            return r"$\pi$"
        if n == -1:
            return r"$-\pi$"
        else:
            return rf"${n}\pi$"
    else:
        if n == 1:
            return rf"$\frac{{\pi}}{{{d}}}$"
        else:
            return rf"$\frac{{{n}\pi}}{{{d}}}$"
ax1.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
order = np.argsort(angle)
ax1.plot(angle[order], y[order]*100, lw=3.5, color="black")
x_points = np.array([1.178, 0.118, -1.178])
idx = np.abs(angle[:, None] - x_points).argmin(axis=0)
y_points = y[idx]
ax1.plot(x_points[0], y_points[0] * 100, "o", color="blue", ms=10)
ax1.plot(x_points[1], y_points[1] * 100, "o", color="orange", ms=10)
ax1.plot(x_points[2], y_points[2] * 100, "o", color="red", ms=10)
ax1.set_xlabel(r"$\vartheta \, [rad]$")
ax1.set_ylabel(r"$y \, [\text{cm}]$")
ax1.set_xlim(-np.pi,np.pi)
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.grid(True)
ax2.plot(t, y1*100, label=r"$\vartheta = 1.178 \, \text{rad}$", lw=3.5, color="blue")
ax2.plot(t, y2*100, label=r"$\vartheta = 0.118 \, \text{rad}$", lw=3.5, color="orange")
ax2.plot(t, y3*100, label=r"$\vartheta = -1.178 \, \text{rad}$", lw=3.5, color="red")
ax2.set_xlabel(r"$t \, [\text{s}]$")
ax2.set_ylabel(r"$y \, [\text{cm}]$")
ax2.grid(True)
ax2.legend()
ax2.set_xlim(0,1)
ax2.set_ylim(-10,1)
plt.savefig("FIGG_polar.pdf", dpi=300)
#"""
