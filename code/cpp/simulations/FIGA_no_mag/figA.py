import numpy as np 
import matplotlib.pyplot as plt 
import simutils as su

plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{siunitx} \usepackage{bm}",
    "font.size": 22,
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
})

g = 9.80665

datasets = su.extract()
t = datasets["data5"].t

x_deg5, y_deg5, z_deg5 = su.incline(0,0,0,0,0,0,5,t,g)
x_deg10, y_deg10, z_deg10 = su.incline(0,0,0,0,0,0,10,t,g)
x_deg15, y_deg15, z_deg15 = su.incline(0,0,0,0,0,0,15,t,g)


def print_vector_error_stats(name, x_ref, y_ref, z_ref, x_test, y_test, z_test):
    dx = np.asarray(x_ref) - np.asarray(x_test)
    dy = np.asarray(y_ref) - np.asarray(y_test)
    dz = np.asarray(z_ref) - np.asarray(z_test)

    err = np.sqrt(dx**2 + dy**2 + dz**2)

    print(f"\n{name}")
    print("-" * (len(name) + 1))
    print(f"N            = {err.size}")
    print(f"Mean         = {np.mean(err):.6e}")
    print(f"Median       = {np.median(err):.6e}")
    print(f"Std. Dev.    = {np.std(err, ddof=1):.6e}")
    print(f"RMSE         = {np.sqrt(np.mean(err**2)):.6e}")
    print(f"Min          = {np.min(err):.6e}")
    print(f"Max          = {np.max(err):.6e}")
    print(f"95th pct.    = {np.percentile(err, 95):.6e}")
"""
print_vector_error_stats(
    "5° Inclination",
    datasets["data5"].x, datasets["data5"].y, datasets["data5"].z,
    x_deg5, y_deg5, z_deg5
)

print_vector_error_stats(
    "10° Inclination",
    datasets["data10"].x, datasets["data10"].y, datasets["data10"].z,
    x_deg10, y_deg10, z_deg10
)

print_vector_error_stats(
    "15° Inclination",
    datasets["data15"].x, datasets["data15"].y, datasets["data15"].z,
    x_deg15, y_deg15, z_deg15
)
"""
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
ax1.tick_params(which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax1.plot(t, y_deg5, color='grey', label=r"$\varphi = 5^\circ$", lw=3.5)
ax1.plot(t, y_deg10, color='blue', label=r"$\varphi = 10^\circ$", lw=3.5)
ax1.plot(t, y_deg15, color='red', label=r"$\varphi = 15^\circ$", lw=3.5)
ax1.plot(t, datasets["data5"].y, 'o', markevery=2000, color='black', ms=7, label=r"$\textrm{numerical}$")
ax1.plot(t, datasets["data10"].y, 'o', markevery=2000, color='black', ms=7)
ax1.plot(t, datasets["data15"].y, 'o', markevery=2000, color='black', ms=7)
ax1.set_xlabel(r"$t \, [\unit{s}]$")
ax1.set_ylabel(r"$y \, [\unit{m}]$")
ax1.set_xlim([0,10])
ax1.set_ylim([-60,7.5])
ax1.grid(True, alpha=0.3)
ax1.legend(loc="lower left")
ax1.text(
    0.03,
    0.97,
    "(a)",
    transform=ax1.transAxes,
    ha="left",
    va="top",
    fontsize=16,
    fontweight="bold",
)
ax2.tick_params(which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax2.set_xlim([0,10])
ax2.set_ylim([-60,7.5])
ax2.plot(t, datasets["data10"].y, color='black', label=r"$\textrm{no friction}$", lw=3.5)
ax2.plot(t, datasets["data10c"].y, color='blue', label=r"$\textrm{dry}$", lw=3.5)
ax2.plot(t, datasets["data10cs"].y, '--', color='magenta', label=r"$\textrm{smoothed}$", lw=2.5)
ax2.plot(t, datasets["data10v"].y, color='orange', label=r"$\textrm{viscous}$", lw=3.5)
ax2.plot(t, datasets["data10a"].y, color='red', label=r"$\textrm{aerodynamic}$", lw=3.5)
ax2.set_xlabel(r"$t \, [\unit{s}]$")
ax2.set_ylabel(r"$y \, [\unit{m}]$")
ax2.grid(True, alpha=0.3)
ax2.legend(loc="lower left")
ax2.text(
    0.03,
    0.97,
    "(b)",
    transform=ax2.transAxes,
    ha="left",
    va="top",
    fontsize=16,
    fontweight="bold",
)
plt.savefig("FIGA_apdx.pdf", dpi=300)
