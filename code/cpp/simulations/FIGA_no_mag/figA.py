import numpy as np 
import matplotlib.pyplot as plt 
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

datasets = su.extract()
t = datasets["data5"].t

x_deg5, y_deg5, z_deg5 = su.incline(0,0,0,0,0,0,5,t,g)
x_deg10, y_deg10, z_deg10 = su.incline(0,0,0,0,0,0,10,t,g)
x_deg15, y_deg15, z_deg15 = su.incline(0,0,0,0,0,0,15,t,g)

fig2, ax2 = plt.subplots(layout='constrained')
ax2.plot(t, y_deg5, color='grey', label=r"$\varphi = 5^\circ$", lw=3.5)
ax2.plot(t, y_deg10, color='blue', label=r"$\varphi = 10^\circ$", lw=3.5)
ax2.plot(t, y_deg15, color='red', label=r"$\varphi = 15^\circ$", lw=3.5)
ax2.plot(t, datasets["data5"].y, 'o', markevery=2000, color='black', ms=7, label="numerical")
ax2.plot(t, datasets["data10"].y, 'o', markevery=2000, color='black', ms=7)
ax2.plot(t, datasets["data15"].y, 'o', markevery=2000, color='black', ms=7)
ax2.grid(True)
ax2.set_xlabel(r"$t \, [\text{s}]$")
ax2.set_ylabel(r"$y \, [\text{m}]$")
ax2.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax2.set_xlim([0,10])
ax2.set_ylim([-60,1])
ax2.legend()
plt.savefig("FIGA1_no_mag_traj.pdf", dpi=150)

fig3, ax3 = plt.subplots(layout='constrained')
ax3.plot(t, datasets["data10"].E, color='black', label=r"no friction", lw=3.5)
ax3.plot(t, datasets["data10c"].E, color='blue', label=r"Coulomb", lw=3.5)
ax3.plot(t, datasets["data10v"].E, color='red', label=r"viscous", lw=3.5)
ax3.plot(t, datasets["data10a"].E, color='magenta', label=r"aerodynamic", lw=3.5)
ax3.grid(True)
ax3.set_xlabel(r"$t \, [\text{s}]$")
ax3.set_ylabel(r"$E \, [\text{J}]$")
ax3.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax3.set_xlim([0,10])
#ax3.set_ylim([-60,1])
ax3.legend()
plt.savefig("FIGA2_no_mag_fric_energ.pdf", dpi=150)

fig4, ax4 = plt.subplots(layout='constrained')
ax4.plot(t, datasets["data10"].y, color='black', label=r"no friction", lw=3.5)
ax4.plot(t, datasets["data10c"].y, color='blue', label=r"Coulomb", lw=3.5)
ax4.plot(t, datasets["data10v"].y, color='red', label=r"viscous", lw=3.5)
ax4.plot(t, datasets["data10a"].y, color='magenta', label=r"aerodynamic", lw=3.5)
ax4.grid(True)
ax4.set_xlabel(r"$t \, [\text{s}]$")
ax4.set_ylabel(r"$y \, [\text{m}]$")
ax4.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax4.set_xlim([0,10])
ax4.set_ylim([-60,1])
ax4.legend()
plt.savefig("FIGA3_no_mag_fric_traj.pdf", dpi=150)
