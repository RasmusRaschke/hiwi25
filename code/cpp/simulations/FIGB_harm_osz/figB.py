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
mu_init = [0,0.174,0.985]
beta0 = 0.174533

datasets = su.extract()

x_ho_vgl, y_ho_vgl, z_ho_vgl = su.get_pos(1,0,0,0,0,0,0,mu_init,1,100.0, t, B, M, R)
x_uho_sth_vgl, y_uho_sth_vgl, z_uho_sth_vgl = harm_osz(2,0,0,0,0,0,0, mu_init,0.002, 100.0, t, B, M, R)
x_dry_sth_vgl, y_dry_sth_vgl, z_dry_sth_vgl = harm_osz(5,0,0,0,0,0,0, mu_init,0.002, 100.0, t, B, M, R)
x_uho_vgl, y_uho_vgl, z_uho_vgl = harm_osz(2,0,0,0,0,0,0, mu_init,54.0, 0.01, t, B, M, R)
ev1_uho_upper, ev1_uho_lower = envelope(1, 0.0, beta0, np.linalg.norm(mu_init), 54.0, 0.01, t, B, M, R)
x_cho_vgl, y_cho_vgl, z_cho_vgl = harm_osz(3,0,0,0,0,0,0, mu_init,540.0, 0.01, t, B, M, R)
x_oho_vgl, y_oho_vgl, z_oho_vgl = harm_osz(4,0,0,0,0,0,0, mu_init,5400.0, 0.01, t, B, M, R)
x_dry_vgl, y_dry_vgl, z_dry_vgl = harm_osz(5,0,0,0,0,0,0, mu_init,0.002, 1, t, B, M, R)

fig1, ax1 = plt.subplots(layout='constrained')
ax1.plot(t, y_ho_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax1.plot(t, datasets["harm_osz"].y*1000, 'o-', markevery=400, color='black', ms=7, label=r"numerical")
ax1.grid(True)
ax1.set_xlabel(r"$t \, [\text{s}]$")
ax1.set_ylabel(r"$y \, [\text{mm}]$")
ax1.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax1.set_xlim([0,3])
ax1.set_ylim([-2.1,0.1])
ax1.legend(loc='lower center', ncols=2)
plt.savefig("FIGB1_harm_osz_vgl.pdf", dpi=100)

fig2, ax2 = plt.subplots(layout='constrained')
ax2.plot(t, y_uho_sth_vgl*1000, color='red', label=r"underdamped", lw=3.5)
ax2.plot(t, y_dry_sth_vgl*1000, color='blue', label=r"dry", lw=3.5)
ax2.plot(t, datasets["smooth_harm_osz"].y*1000, 'o-', markevery=500, color='black', ms=7, lw=2, label=r"numerical")
ax2.grid(True)
ax2.set_xlabel(r"$t \, [\text{s}]$")
ax2.set_ylabel(r"$y \, [\text{mm}]$")
ax2.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax2.set_xlim([0,3])
ax2.set_ylim([-2.1,0.1])
ax2.legend(loc='lower center', ncols=2)
plt.savefig("FIGB2_smooth_harm_osz_vgl.pdf", dpi=100)

fig3, ax3 = plt.subplots(layout='constrained')
ax3.plot(t, y_uho_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax3.plot(t, ev1_uho_upper*1000, '--', color='red', label=r"envelope", lw=2)
ax3.plot(t, datasets["ud_harm_osz"].y*1000, 'o-', markevery=500, color='black', ms=7,lw=2, label=r"numerical")
ax3.grid(True)
ax3.set_xlabel(r"$t \, [\text{s}]$")
ax3.set_ylabel(r"$y \, [\text{mm}]$")
ax3.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax3.set_xlim([0,3])
ax3.set_ylim([-1.75,0.1])
ax3.legend(loc='lower center', ncols=2)
plt.savefig("FIGB3_udamp_harm_osz_vgl.pdf", dpi=100)

fig4, ax4 = plt.subplots(layout='constrained')
ax4.plot(t, y_cho_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax4.plot(t, datasets["cd_harm_osz"].y*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
ax4.grid(True)
ax4.set_xlabel(r"$t \, [\text{s}]$")
ax4.set_ylabel(r"$y \, [\text{mm}]$")
ax4.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax4.set_xlim([0,3])
ax4.set_ylim([-1.1,0.1])
ax4.legend(loc='lower center', ncols=2)
plt.savefig("FIGB4_cdamp_harm_osz_vgl.pdf", dpi=100)

fig5, ax5 = plt.subplots(layout='constrained')
ax5.plot(t, y_oho_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax5.plot(t, datasets["od_harm_osz"].y*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
ax5.grid(True)
ax5.set_xlabel(r"$t \, [\text{s}]$")
ax5.set_ylabel(r"$y \, [\text{mm}]$")
ax5.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax5.set_xlim([0,3])
ax5.set_ylim([-1,0.1])
ax5.legend(loc='lower left', ncols=2)
plt.savefig("FIGB5_odamp_harm_osz_vgl.pdf", dpi=100)

fig6, ax6 = plt.subplots(layout='constrained')
ax6.plot(t, y_dry_vgl*1000, color='red', label=r"analytical", lw=3.5)
ax6.plot(t, datasets["dry_harm_osz"].y*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
ax6.grid(True)
ax6.set_xlabel(r"$t \, [\text{s}]$")
ax6.set_ylabel(r"$y \, [\text{mm}]$")
ax6.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax6.set_xlim([0,3])
ax6.set_ylim([-2,0.1])
ax6.legend(loc='lower center', ncols=2)
plt.savefig("FIGB6_dry_harm_osz_vgl.pdf", dpi=100)

fig7, ax7 = plt.subplots(layout="constrained")
ax7.plot(t, datasets["smooth_harm_osz"].E*1e6, color='red', label=r"smoothed", lw=3.5)
ax7.plot(t, datasets["ud_harm_osz"].E*1e6, color='blue', label=r"underdamped", lw=3.5)
ax7.plot(t, datasets["cd_harm_osz"].E*1e6, color='aqua', label=r"critial", lw=3.5)
ax7.plot(t, datasets["od_harm_osz"].E*1e6, color='darkblue', label=r"overdamped", lw=3.5)
ax7.plot(t, datasets["dry_harm_osz"].E*1e6, color='deeppink', label=r"dry", lw=3.5)
ax7.set_xlabel(r"$t \, [\text{s}]$")
ax7.set_ylabel(r"$E \, [\text{µJ}]$")
ax7.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax7.grid(True)
ax7.set_xlim([0,3])
ax7.legend(loc='upper right', ncols=2)
plt.savefig("FIGB7_energies.pdf", dpi=100)