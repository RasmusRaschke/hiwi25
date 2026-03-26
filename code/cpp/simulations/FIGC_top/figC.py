import numpy as np 
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks
import simutils as su

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

t = datasets["test"].t
x_test, y_test, z_test = magnetic_top(0,0,0,[0.05,0.05,1.0],np.deg2rad(1000),t)

fig1, ax1 = plt.subplots(layout='constrained')
ax1.plot(x_test*1000, y_test*1000, lw=3, color='red')
ax1.plot(datasets["test"].x *1000, datasets["test"].y *1000, "-o", markevery=1000, color='black', lw=1.0)
ax1.grid(True)
ax1.set_xlabel(r"$x \, [\text{mm}]$")
ax1.set_ylabel(r"$y \, [\text{mm}]$")
ax1.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
plt.savefig("test.pdf", dpi=100)
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
