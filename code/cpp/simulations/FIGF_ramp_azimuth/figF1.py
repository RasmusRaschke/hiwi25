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
"""
azimuth0, x0 = su.batch_extract_azimuth("mag_moment_data/c0_0")
azimuth1, x1 = su.batch_extract_azimuth("mag_moment_data/c0_5")
azimuth2, x2 = su.batch_extract_azimuth("mag_moment_data/c1_0")
azimuth3, x3 = su.batch_extract_azimuth("mag_moment_data/c1_5")
azimuth4, x4 = su.batch_extract_azimuth("mag_moment_data/c2_0")
np.savez("mag_moment_azimuth_0.npz", azimuth=azimuth0, x=x0)
np.savez("mag_moment_azimuth_1.npz", azimuth=azimuth1, x=x1)
np.savez("mag_moment_azimuth_2.npz", azimuth=azimuth2, x=x2)
np.savez("mag_moment_azimuth_3.npz", azimuth=azimuth3, x=x3)
np.savez("mag_moment_azimuth_4.npz", azimuth=azimuth4, x=x4)
"""

azimuth_sorted = []
x_sorted = []
for i in [0,1,2,3,4]:
    d = np.load(f"mag_moment_azimuth_{i}.npz")
    azimuth = d["azimuth"]
    x = d["x"]
    order = np.argsort(azimuth)
    azimuth_sorted.append(azimuth[order])
    x_sorted.append(x[order])

colors = plt.cm.viridis(np.linspace(0,1,5))
fig1, ax1 = plt.subplots(figsize=(7, 5), layout="constrained")
def pi_formatter(x, pos):
    frac = Fraction(x / np.pi).limit_denominator(12)
    n, d = frac.numerator, frac.denominator
    if n == 0:
        return "0"
    elif d == 1:
        if n == 1:
            return r"$\pi$"
        else:
            return rf"${n}\pi$"
    else:
        if n == 1:
            return rf"$\frac{{\pi}}{{{d}}}$"
        else:
            return rf"$\frac{{{n}\pi}}{{{d}}}$"
ax1.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
ax1.plot(azimuth_sorted[0], x_sorted[0] * 100, lw=2.5, color=colors[0], label=r"$\mu = 0.0 \, [\text{Am}^2]$")
ax1.plot(azimuth_sorted[1], x_sorted[1] * 100, lw=2.5, color=colors[1], label=r"$\mu = 0.5 \, [\text{Am}^2]$")
ax1.plot(azimuth_sorted[2], x_sorted[2] * 100, lw=2.5, color=colors[2], label=r"$\mu = 1.0 \, [\text{Am}^2]$")
ax1.plot(azimuth_sorted[3], x_sorted[3] * 100, lw=2.5, color=colors[3], label=r"$\mu = 1.5 \, [\text{Am}^2]$")
ax1.plot(azimuth_sorted[4], x_sorted[4] * 100, lw=2.5, color=colors[4], label=r"$\mu = 2.0 \, [\text{Am}^2]$")
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.set_xlabel(r"$\varphi \, [\text{rad}]$")
ax1.set_ylabel(r"$x \, [\text{cm}]$")
ax1.set_xlim(0,np.pi)
ax1.grid(True)
ax1.legend()
plt.savefig("FIGF1_mag_azimuth.pdf", dpi=300)
#"""