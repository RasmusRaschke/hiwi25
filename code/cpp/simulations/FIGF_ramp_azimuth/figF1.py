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

#NOTE: this formatter is technically a pullback along a reparametrization of the input, 
# so yes, it does not make much sense out of context
def pi_formatter(x, pos):
    frac = Fraction(x / np.pi).limit_denominator(12)
    n, d = frac.numerator, frac.denominator
    if n == 0:
        return r"$-\frac{\pi}{2}$"
    elif d == 1:
        if n == 1:
            return r"$\frac{\pi}{2}$"
        else:
            return rf"${n}\pi$"
    else:
        if d == 2:
            return r"$0$"
        elif d == 4:
            if n == 1:
                return r"$- \frac{\pi}{4}$"
            else:
                return r"$\frac{\pi}{4}$"
        elif n == 1:
            return rf"$\frac{{\pi}}{{{d}}}$"
        else:
            return rf"$\frac{{{n}\pi}}{{{d}}}$"
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
cbar.set_label(r"$\mu\,[\mathrm{Am}^2]$")
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.set_xlabel(r"$\varphi \, [\text{rad}]$")
ax1.set_ylabel(r"$x \, [\text{cm}]$")
ax1.set_xlim(0,np.pi)
ax1.grid(True, alpha=0.3)
ax2.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax2.set_xlabel(r"$\varphi \, [\text{rad}]$")
ax2.set_ylabel(r"$x \, [\text{cm}]$")
ax2.set_xlim(0,np.pi)
ax2.grid(True, alpha=0.3)
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
plt.savefig("FIGF1_mag_azimuth.pdf", dpi=300)
"""
fig2, ax2 = plt.subplots(figsize=(7, 5), layout="constrained")
ax2.xaxis.set_major_formatter(tck.FuncFormatter(pi_formatter))
plot_mus = np.round(np.arange(1.0, 2.5 + 0.1, 0.1), 1)
colors = plt.cm.plasma(np.linspace(0, 1, len(plot_mus)))

for i, mu in enumerate(plot_mus):
    ax1.plot(
        azimuth_sorted[mu],
        x_sorted[mu] * 100,
        lw=2.5,
        color=colors[i],
        label=rf"$\mu = {mu:.1f} \, [\text{{Am}}^2]$"
    )
ax1.xaxis.set_major_locator(tck.MultipleLocator(np.pi / 4))
ax1.set_xlabel(r"$\varphi \, [\text{rad}]$")
ax1.set_ylabel(r"$x \, [\text{cm}]$")
ax1.set_xlim(0,np.pi)
ax1.grid(True)
#ax1.legend()
plt.savefig("FIGF1_mag_azimuth.pdf", dpi=300)
"""