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
mu_par  = [m[2] for m in mus]

def centered_projection(mu, pos, scale=1.0, stretch=1.0):
    mu_x = np.asarray(mu[0])
    mu_y = np.asarray(mu[1])
    x = np.asarray(pos[0])
    y = np.asarray(pos[1])
    return stretch * mu_x - scale * x, stretch * mu_y - scale * y

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
ax1.plot(t, mu_perp[0] * 1000, lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\mu_\perp \, [\text{mAm}^2]$")
ax1.set_xlim(0,1.0)
ax1.set_ylim(40,105)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax1.grid(True)
ax2.plot(t, mu_par[0] * 1000, lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\mu_\| \, [\text{mAm}^2]$")
ax2.set_xlim(0,1.0)
ax2.set_ylim(-90,95)
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

mask = t <= 1.0
fig2, ax3 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x_c, mu_y_c = centered_projection(mu2, pos2, scale=1.0, stretch=1.0)
ax3.plot(mu_x_c[mask], mu_y_c[mask])
#ax3.plot(pos1[0]*0.001, pos1[1]*0.001)
ax3.scatter(mu_x_c[0], mu_y_c[0], marker="o")
ax3.set_xlabel(r"$\mu_x - x_c \, [\mathrm{Am}^2]$")
ax3.set_ylabel(r"$\mu_y - y_c \, [\mathrm{Am}^2]$")
ax3.grid(True, alpha=0.3)
plt.savefig("FIGH3_projection_looping.pdf", dpi=300)

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), layout="constrained")
ax1.plot(t, mu_perp[1], lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\mu_\perp \, [\text{Am}^2]$")
ax1.set_xlim(0,1.0)
ax1.set_ylim(0.0,1.6)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.25))
ax1.grid(True)
ax2.plot(t, mu_par[1], lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\mu_\| \, [\text{Am}^2]$")
ax2.set_xlim(0,1.0)
ax2.set_ylim(-0.4,1.65)
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

mask = t <= 1.0
fig2, ax3 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x_c, mu_y_c = centered_projection(mu1, pos1, scale=1.0, stretch=1.0)
ax3.plot(mu_x_c[mask], mu_y_c[mask])
#ax3.plot(pos1[0]*0.001, pos1[1]*0.001)
ax3.scatter(mu_x_c[0], mu_y_c[0], marker="o")
ax3.set_xlabel(r"$\mu_x - x_c \, [\mathrm{Am}^2]$")
ax3.set_ylabel(r"$\mu_y - y_c \, [\mathrm{Am}^2]$")
ax3.grid(True, alpha=0.3)
plt.savefig("FIGH2_projection_rolling.pdf", dpi=300)

