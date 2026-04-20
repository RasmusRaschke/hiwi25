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
t = datasets["azi_0"].t
mu1 = np.array([datasets["azi_0"].mu_x, datasets["azi_0"].mu_y, datasets["azi_0"].mu_z])
mu2 = np.array([datasets["azi_pi_5"].mu_x, datasets["azi_pi_5"].mu_y, datasets["azi_pi_5"].mu_z])
mu3 = np.array([datasets["azi_4pi_5"].mu_x, datasets["azi_4pi_5"].mu_y, datasets["azi_4pi_5"].mu_z])
mu4 = np.array([datasets["azi_pi"].mu_x, datasets["azi_pi"].mu_y, datasets["azi_pi"].mu_z])
mus = [mu1, mu2, mu3, mu4]
mu_perp = [np.sqrt(m[0]**2 + m[1]**2) for m in mus]
mu_par  = [m[2] for m in mus]

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
ax1.plot(t, mu_perp[1], lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\mu_\perp \, [\text{Am}^2]$")
ax1.set_xlim(0,0.5)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.1))
ax1.grid(True)
ax2.plot(t, mu_par[1], lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\mu_\| \, [\text{Am}^2]$")
ax2.set_xlim(0,0.5)
ax2.xaxis.set_major_locator(tck.MultipleLocator(0.1))
ax2.grid(True)
plt.savefig("FIGH1_magnet_rolling.pdf", dpi=300)

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
ax1.plot(t, mu_perp[0], lw=3.5, color="black")
ax1.set_xlabel(r"$t \, [s]$")
ax1.set_ylabel(r"$\mu_\perp \, [\text{Am}^2]$")
ax1.set_xlim(0,1)
ax1.xaxis.set_major_locator(tck.MultipleLocator(0.2))
ax1.grid(True)
ax2.plot(t, mu_par[0], lw=3.5, color="black")
ax2.set_xlabel(r"$t \, [s]$")
ax2.set_ylabel(r"$\mu_\| \, [\text{Am}^2]$")
ax2.set_xlim(0,1)
ax2.xaxis.set_major_locator(tck.MultipleLocator(0.2))
ax2.grid(True)
plt.savefig("FIGH2_magnet_oscillating.pdf", dpi=300)

mask = t <= 1
fig2, ax3 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x = np.asarray(mu2[0])
mu_y = np.asarray(mu2[1])
ax3.plot(mu_x[mask], mu_y[mask])
ax3.scatter(mu_x[0], mu_y[0], marker="o")   # start point
#ax3.scatter(mu_x[-1] * 100, mu_y[-1] * 100, marker="x") # end point
ax3.set_xlabel(r"$\mu_x \, [\text{Am}^2]$")
ax3.set_ylabel(r"$\mu_y \, [\text{Am}^2]$")
#ax3.xaxis.set_major_locator(tck.MultipleLocator(0.2))
ax3.grid(True, alpha=0.3)
plt.savefig("FIGH3_projection_rolling.pdf", dpi=300)

mask = t <= 3
fig3, ax4 = plt.subplots(figsize=(7, 5), layout="constrained")
mu_x = np.asarray(mu1[0])
mu_y = np.asarray(mu1[1])
ax4.plot(mu_x[mask], mu_y[mask])
ax4.scatter(mu_x[0], mu_y[0], marker="o")   # start point
#ax4.scatter(mu_x[-1], mu_y[-1], marker="x") # end point
ax4.set_xlabel(r"$\mu_x \, [\text{Am}^2]$")
ax4.set_ylabel(r"$\mu_y \, [\text{Am}^2]$")
#ax3.xaxis.set_major_locator(tck.MultipleLocator(0.2))
ax4.grid(True, alpha=0.3)
plt.savefig("FIGH4_projection_oscillating.pdf", dpi=300)