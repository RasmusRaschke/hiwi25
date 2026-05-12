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

datasets = su.extract()
t = datasets["azi_pi_5"].t
mu1 = np.array([datasets["azi_pi_5"].mu_x, datasets["azi_pi_5"].mu_y, datasets["azi_pi_5"].mu_z])
pos1 = np.array([datasets["azi_pi_5"].x, datasets["azi_pi_5"].y, datasets["azi_pi_5"].z])
mu2 = np.array([datasets["azi_pi_2"].mu_x, datasets["azi_pi_2"].mu_y, datasets["azi_pi_2"].mu_z])
pos2 = np.array([datasets["azi_pi_2"].x, datasets["azi_pi_2"].y, datasets["azi_pi_2"].z])
mus = [mu1, mu2]
mu_perp = [np.sqrt(m[0]**2 + m[1]**2) for m in mus]
mu_par  = [abs(m[2]) for m in mus]

def centered_projection(mu, pos, t, tmax, frac):
    mask = t <= tmax
    mu_x = np.asarray(mu[0])[mask]
    mu_y = np.asarray(mu[1])[mask]
    x = np.asarray(pos[0])[mask]
    y = np.asarray(pos[1])[mask]
    mu_span = max(np.ptp(mu_x), np.ptp(mu_y))
    pos_span = max(np.ptp(x), np.ptp(y))
    if pos_span == 0:
        scale = 1.0
    else:
        scale = frac * mu_span / pos_span
    x_plot = scale * x
    y_plot = scale * y
    return mu_x, mu_y, x_plot, y_plot

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
ax1.plot(t, mu_perp[0] * 1000, lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\|\mu_\perp \| \, [\text{mAm}^2]$")
ax1.set_xlim(0,1.0)
ax1.set_ylim(40,105)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax1.grid(True)
ax2.plot(t, mu_par[0] * 1000, lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\| \mu_\| \| \, [\text{mAm}^2]$")
ax2.set_xlim(0,1.0)
ax2.set_ylim(-5,95)
ax2.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax2.grid(True)
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
plt.savefig("FIGH1_magnet_rolling.pdf", dpi=300)

fig2, ax3 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x, mu_y, x, y = centered_projection(mu2, pos2, t, 3.0, 0.75)
ax3.plot(mu_x, mu_y)
ax3.plot(x, y)
idx = np.arange(0, len(mu_x), 300)
for i in idx:
    ax3.plot([x[i], mu_x[i]], [y[i], mu_y[i]], "k-", alpha=0.25, linewidth=0.8, color="black")
ax3.scatter(x[0], y[0], marker="o", color="orange")
ax3.set_xlabel(r"$\mu_x - c x \, [\text{arb. u.}]$")
ax3.set_ylabel(r"$\mu_y - c y \, [\text{arb. u.}]$")
ax3.grid(True, alpha=0.3)
plt.savefig("FIGH3_projection_looping.pdf", dpi=300)

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
ax1.plot(t, mu_perp[1], lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\| \mu_\perp \| \, [\text{Am}^2]$")
ax1.set_xlim(0,1.0)
ax1.set_ylim(0.0,1.6)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax1.grid(True)
ax2.plot(t, mu_par[1], lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\|\mu_\| \| \, [\text{Am}^2]$")
ax2.set_xlim(0,1.0)
ax2.set_ylim(-0.05,1.65)
ax2.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax2.grid(True)
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
plt.savefig("FIGH4_magnet_looping.pdf", dpi=300)

fig2, ax3 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x, mu_y, x, y = centered_projection(mu1, pos1, t, 1.0, 0.75)
ax3.plot(mu_x, mu_y)
ax3.plot(x, y)
#ax3.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] - y[:-1], angles="xy", scale_units="xy", scale=1, alpha=0.4)
ax3.scatter(x[0], y[0], marker="o", color="orange")
idx = np.arange(0, len(mu_x), 300)
for i in idx:
    ax3.plot([x[i], mu_x[i]], [y[i], mu_y[i]], "k-", alpha=0.25, linewidth=0.8, color="black")
ax3.set_xlabel(r"$\mu_x - c x \, [\text{arb. u.}]$")
ax3.set_ylabel(r"$\mu_y - c y \, [\text{arb. u.}]$")
ax3.grid(True, alpha=0.3)
plt.savefig("FIGH2_projection_rolling.pdf", dpi=300)




fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
ax1.plot(t, mu_perp[0], "--", lw=3.5 ,color="blue", label=r"$\|\mu_\perp \| \, [\text{Am}^2]$")
ax1.plot(t, mu_par[0], "--", lw=3.5, color="red", label=r"$\|\mu_\| \| \, [\text{Am}^2]$")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\mu \, [\text{Am}^2]$")
ax1.set_xlim(0,0.5)
ax1.legend()
ax1.set_ylim(-0.005,0.12)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.1))
ax1.grid(True)
ax2.plot(t, mu_perp[1], "--", lw=3.5, color="blue", label=r"$\|\mu_\perp \| \, [\text{Am}^2]$")
ax2.plot(t, mu_par[1], "--", lw=3.5, color="red", label=r"$\|\mu_\| \| \, [\text{Am}^2]$")
ax2.legend()
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\mu \, [\text{Am}^2]$")
ax2.set_xlim(0,0.5)
ax2.set_ylim(-0.01,1.75)
ax2.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax2.grid(True)
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
plt.savefig("FIGH5_overlay.pdf", dpi=300)