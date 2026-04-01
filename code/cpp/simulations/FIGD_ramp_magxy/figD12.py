import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
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
mx, my, z = su.batch_extract("mag_moment_data")
np.savez("mag_moment_grid.npz", mx=mx, my=my, z=z)
"""
data = np.load("mag_moment_grid.npz")
mx = data["mx"]
my = data["my"]
z = data["z"]
mx_sort = np.sort(np.unique(mx))
my_sort = np.sort(np.unique(my))
Z = np.full((len(my_sort), len(mx_sort)), np.nan)
mx_idx = {v: i for i,v in enumerate(mx_sort)}
my_idx = {v: i for i,v in enumerate(my_sort)}

for mxi, myi, zi in zip(mx, my, z):
    Z[my_idx[myi], mx_idx[mxi]] = zi
cmap = plt.cm.plasma
X, Y = np.meshgrid(mx_sort, my_sort)
fig1 = plt.figure(figsize=(9,7), constrained_layout=True)
ax1 = fig1.add_subplot(111, projection="3d")
ax1.plot_surface(X, Y, Z*100, cmap=cmap, edgecolor="royalblue", lw=0.5, rstride=1, cstride=1, alpha=0.4)
ax1.contour(X, Y, Z*100, zdir='x', offset=-1.5, cmap=cmap)
ax1.contour(X, Y, Z*100, zdir='y', offset=1.5, cmap=cmap)
#ax1.contour(X, Y, Z*100, zdir='z', offset=-40, cmap="coolwarm")
ax1.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), zlim=(-40, 40))
ax1.set_xlabel(r"$\mu_x \, [\text{Am}^2]$", labelpad=14)
ax1.set_ylabel(r"$\mu_y \, [\text{Am}^2]$", labelpad=14)
ax1.set_zlabel(r"$x \, [\text{cm}]$", labelpad=14)
#plt.tight_layout()
ax1.tick_params(axis="both", which="major", pad=7)
plt.savefig("FIGD1_mag_mom_3d.pdf", dpi=300)


fig2, (ax2, ax3) = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
vmin = min(mx_sort.min(), my_sort.min())
vmax = max(mx_sort.max(), my_sort.max())
norm = Normalize(vmin=vmin, vmax=vmax)
for target in my_sort:
    mask = np.isclose(my, target)
    mx_line = mx[mask]
    z_line = z[mask]
    order = np.argsort(mx_line)
    ax2.plot(mx_line[order], z_line[order]*100, color=cmap(norm(target)), lw=1.5)
ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
ax2.set_xlabel(r"$\mu_x \, [\text{Am}^2]$")
ax2.set_ylabel(r"$x \, [\text{cm}]$")
ax2.set_xlim(-1,1)
ax2.grid(True, alpha=0.3)

for target in mx_sort:
    mask = np.isclose(mx, target)
    my_line = my[mask]
    z_line = z[mask]
    order = np.argsort(my_line)
    ax3.plot(my_line[order], z_line[order]*100, color=cmap(norm(target)), lw=1.5)
ax3.xaxis.set_major_locator(plt.MaxNLocator(5))
ax3.set_xlabel(r"$\mu_y \, [\text{Am}^2]$")
ax3.set_xlim(-1,1)
#ax3.set_ylabel(r"$x \, [\text{cm}]$")
ax3.grid(True, alpha=0.3)
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig2.colorbar(sm, ax=[ax2, ax3], location="right", pad=0.02)
cbar.set_label(r"$\mu_{x,y}^{\text{slice}} \, [\text{Am}^2]$")
cbar.ax.yaxis.set_major_locator(plt.MaxNLocator(4))
plt.savefig("FIGD3_mag_mom_2d.pdf", dpi=300)

'''

fig1, ax1 = plt.subplots(layout='constrained')
ax1.plot(t, y_ho_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax1.plot(t, y_ho*1000, 'o-', markevery=400, color='black', ms=7, label=r"numerical")
ax1.grid(True)
ax1.set_xlabel(r"$t \, [\text{s}]$")
ax1.set_ylabel(r"$y \, [\text{mm}]$")
ax1.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax1.set_xlim([0,3])
ax1.set_ylim([-2.1,0.1])
ax1.legend(loc='lower center', ncols=2)
plt.savefig("FIGB1_harm_osz_vgl.pdf", dpi=100)


fig2, ax2 = plt.subplots(layout='constrained')
#ax2.plot(t, x_deg5, color='darkred', label=r"$\tilde{x}(t)$")
ax2.plot(t5, y_deg5, color='grey', label=r"$\varphi = 5^\circ$", lw=3.5)
#ax2.plot(t, x_deg10, color='darkred')
ax2.plot(t5, y_deg10, color='blue', label=r"$\varphi = 10^\circ$", lw=3.5)
#ax2.plot(t, x_deg15, color='darkred')
ax2.plot(t5, y_deg15, color='red', label=r"$\varphi = 15^\circ$", lw=3.5)
#ax2.plot(t, x, 'o', markevery=2000, color='tomato', label=r"$x(t)$")
ax2.plot(t5, y5, 'o', markevery=2000, color='black', ms=7, label="numerical")
ax2.plot(t5, y10, 'o', markevery=2000, color='black', ms=7)
ax2.plot(t5, y15, 'o', markevery=2000, color='black', ms=7)
ax2.grid(True)
ax2.set_xlabel(r"$t \, [\text{s}]$")
ax2.set_ylabel(r"$y \, [\text{m}]$")
ax2.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax2.set_xlim([0,10])
ax2.set_ylim([-60,1])
ax2.legend()
plt.savefig("FIGA1_no_mag_traj.pdf", dpi=150)

fig3, ax3 = plt.subplots(layout='constrained')
ax3.plot(t10, E10, color='black', label=r"no friction", lw=3.5)
ax3.plot(t10, E10c, color='blue', label=r"coulomb", lw=3.5)
ax3.plot(t10, E10v, color='red', label=r"viscous", lw=3.5)
ax3.plot(t10, E10a, color='magenta', label=r"aerodynamic", lw=3.5)
ax3.grid(True)
ax3.set_xlabel(r"$t \, [\text{s}]$")
ax3.set_ylabel(r"$E \, [\text{J}]$")
ax3.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax3.set_xlim([0,10])
#ax3.set_ylim([-60,1])
ax3.legend()
plt.savefig("FIGA2_no_mag_fric_energ.pdf", dpi=150)

fig4, ax4 = plt.subplots(layout='constrained')
ax4.plot(t10, y10, color='black', label=r"no friction", lw=3.5)
ax4.plot(t10, y10c, color='blue', label=r"coulomb", lw=3.5)
ax4.plot(t10, y10v, color='red', label=r"viscous", lw=3.5)
ax4.plot(t10, y10a, color='magenta', label=r"aerodynamic", lw=3.5)
ax4.grid(True)
ax4.set_xlabel(r"$t \, [\text{s}]$")
ax4.set_ylabel(r"$y \, [\text{m}]$")
ax4.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax4.set_xlim([0,10])
ax4.set_ylim([-60,1])
ax4.legend()
plt.savefig("FIGA3_no_mag_fric_traj.pdf", dpi=150)
'''
