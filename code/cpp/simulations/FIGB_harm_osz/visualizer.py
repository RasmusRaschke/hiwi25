import numpy as np 
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks
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
'''
omg_test = np.sqrt(5 * np.linalg.norm(mu_init) * B / (7 * M * R**2))
nu_test = 14 * omg_test / (5 * g * k)
print(nu_test)
'''
def beta_ana(system, nu, mu_norm, beta0, k, t):
    if system == 1:
        print("System is harmonic oscillator")
        omega = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        return beta0 * np.cos(omega * t)
    elif system == 2:
        print("System is underdamped harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g * k) / (14 * omega0) 
        print(omega0)
        print(zeta)
        if zeta > 1:
            print("Warning: omega_0 will be imaginary.")
        omegad = omega0 * np.sqrt(1.0 - zeta**2)
        print(np.exp(-zeta * omega0 * t))
        return beta0 * np.exp(-zeta * omega0 * t) * (np.cos(omegad * t) + (zeta * omega0 / omegad) * np.sin(omegad * t))
    elif system == 3:
        print("System is critically damped hamonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        return beta0 * (1 + omega0 * t) * np.exp(- omega0 * t)
    elif system == 4:
        print("System is overdamped harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g * k) / (14 * omega0) 
        if zeta < 1:
            print("Warning: lambda will be imaginary.")
        lambda1 = -zeta * omega0 + omega0 * np.sqrt(zeta**2 - 1)
        lambda2 = -zeta * omega0 - omega0 * np.sqrt(zeta**2 - 1)
        return (beta0 / (lambda1 - lambda2)) * (-lambda2 * np.exp(lambda1 * t) + lambda1 * np.exp(lambda2 * t))
    elif system == 5:
        print("System is dry harmonic oscillator.")
        omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
        zeta = (5 * nu * g * k) / (14 * omega0)
        T_half = np.pi / omega0
        stick = 2.0 * zeta / omega0
        t_arr = np.atleast_1d(t)
        beta = np.zeros_like(t_arr, dtype=float)
        n = np.floor(t_arr / T_half).astype(int)
        tau = t_arr - n * T_half
        A_n = beta0 - 2.0 * stick * n
        mask_stick = A_n <= stick
        beta[mask_stick] = 0.0
        mask = ~mask_stick
        sign = (-1)**n[mask]
        beta[mask] = sign * (stick + (A_n[mask] - stick) * np.cos(omega0 * tau[mask]))
        return beta if np.ndim(t) else beta[0]
    else:
        raise ValueError("Invalid System!")
            

def harm_osz(system, x0, y0, z0, vx0, vy0, vz0, mu, nu, k, t):
    mu_norm = np.linalg.norm(mu)
    beta = beta_ana(system, nu, mu_norm, beta0, k, t)
    x = vx0 * t + x0
    y = y0 + R * (beta - beta0)
    z = vz0 * t + z0
    return x, y, z

def envelope(system, y0, beta0, mu_norm, nu, k, t):
    omega0 = np.sqrt((5 * mu_norm * B) / (7 * R**2 * M))
    if system == 1:
        zeta = (5 * nu * g * k) / (14 * omega0) 
        A = beta0 / (np.sqrt(1-zeta**2)) * np.exp(- zeta * omega0 * t)
        upper = y0 + R * (A - beta0)
        lower = y0 + R * (-A - beta0)
        return upper, lower
    elif system == 2:
        zeta = (5 * nu * g * k) / (14 * omega0) 
        A = np.maximum(beta0 - (4.0 * zeta / np.pi) * t, 0.0)
        upper = y0 + R * (A - beta0)
        lower = y0 + R * (-A - beta0)
        return upper, lower
    else:
        raise ValueError("Invalid System!")

#Comapre analytical solutions without magnetic field

####################No Friction#############################
data = np.genfromtxt('harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
t = cols['t']; x_ho = cols['x']; y_ho = cols['y']; z_ho = cols['z']
Omx_ho = cols['Ox']; Omy_ho = cols['Oy']; Omz_ho = cols['Oz']
T_trans_ho = cols['T_trans']; T_rot_ho = cols['T_rot']; U_gr_ho = cols['U_gr']; U_em_ho = cols['U_em']; E_ho = cols['E']


######################Friction##########################

data = np.genfromtxt('smooth_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_smooth = cols['x']; y_smooth = cols['y']; z_smooth = cols['z']
Omx_smooth = cols['Ox']; Omy_smooth = cols['Oy']; Omz_smooth = cols['Oz']
T_trans_smooth = cols['T_trans']; T_rot_smooth = cols['T_rot']; U_gr_smooth = cols['U_gr']; U_em_smooth = cols['U_em']; E_smooth = cols['E']

data = np.genfromtxt('ud_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_uho = cols['x']; y_uho = cols['y']; z_uho = cols['z']
Omx_uho = cols['Ox']; Omy_uho = cols['Oy']; Omz_uho = cols['Oz']
T_trans_uho = cols['T_trans']; T_rot_uho = cols['T_rot']; U_gr_uho = cols['U_gr']; U_em_uho = cols['U_em']; E_uho = cols['E']
peaks, props = find_peaks(y_uho, distance=50, prominence=1e-8)
t_peak = t[peaks]
y_peak = y_uho[peaks]

data = np.genfromtxt('cd_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_cho = cols['x']; y_cho = cols['y']; z_cho = cols['z']
Omx_cho = cols['Ox']; Omy_cho = cols['Oy']; Omz_cho = cols['Oz']
T_trans_cho = cols['T_trans']; T_rot_cho = cols['T_rot']; U_gr_cho = cols['U_gr']; U_em_cho = cols['U_em']; E_cho = cols['E']

data = np.genfromtxt('od_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_oho = cols['x']; y_oho = cols['y']; z_oho = cols['z']
Omx_oho = cols['Ox']; Omy_oho = cols['Oy']; Omz_oho = cols['Oz']
T_trans_oho = cols['T_trans']; T_rot_oho = cols['T_rot']; U_gr_oho = cols['U_gr']; U_em_oho = cols['U_em']; E_oho = cols['E']

data = np.genfromtxt('dry_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_dry = cols['x']; y_dry = cols['y']; z_dry = cols['z']
Omx_dry = cols['Ox']; Omy_dry = cols['Oy']; Omz_dry = cols['Oz']
T_trans_dry = cols['T_trans']; T_rot_dry = cols['T_rot']; U_gr_dry = cols['U_gr']; U_em_dry = cols['U_em']; E_dry = cols['E']

data = np.genfromtxt('air_harm_osz.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x_aho = cols['x']; y_aho = cols['y']; z_aho = cols['z']
Omx_aho = cols['Ox']; Omy_aho = cols['Oy']; Omz_aho = cols['Oz']
T_trans_aho = cols['T_trans']; T_rot_aho = cols['T_rot']; U_gr_aho = cols['U_gr']; U_em_aho = cols['U_em']; E_aho = cols['E']

x_ho_vgl, y_ho_vgl, z_ho_vgl = harm_osz(1,0,0,0,0,0,0,mu_init,1,100.0, t)

x_uho_sth_vgl, y_uho_sth_vgl, z_uho_sth_vgl = harm_osz(2,0,0,0,0,0,0, mu_init,0.002, 100.0, t)
x_dry_sth_vgl, y_dry_sth_vgl, z_dry_sth_vgl = harm_osz(5,0,0,0,0,0,0, mu_init,0.002, 100.0, t)

x_uho_vgl, y_uho_vgl, z_uho_vgl = harm_osz(2,0,0,0,0,0,0, mu_init,54.0, 0.01, t)
ev1_uho_upper, ev1_uho_lower = envelope(1, 0.0, beta0, np.linalg.norm(mu_init), 54.0, 0.01, t)
x_cho_vgl, y_cho_vgl, z_cho_vgl = harm_osz(3,0,0,0,0,0,0, mu_init,540.0, 0.01, t)
x_oho_vgl, y_oho_vgl, z_oho_vgl = harm_osz(4,0,0,0,0,0,0, mu_init,5400.0, 0.01, t)

x_dry_vgl, y_dry_vgl, z_dry_vgl = harm_osz(5,0,0,0,0,0,0, mu_init,0.002, 1 / R, t) #technically, k is 1 as it does not enter the dry model, but the dry model has an additional 1/R, so this compensates for that

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
ax2.plot(t, y_uho_sth_vgl*1000, color='red', label=r"underdamped", lw=3.5)
ax2.plot(t, y_dry_sth_vgl*1000, color='blue', label=r"dry", lw=3.5)
ax2.plot(t, y_smooth*1000, 'o-', markevery=500, color='black', ms=7, lw=2, label=r"numerical")
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
ax3.plot(t, y_uho*1000, 'o-', markevery=500, color='black', ms=7,lw=2, label=r"numerical")
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
ax4.plot(t, y_cho*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
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
ax5.plot(t, y_oho*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
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
ax6.plot(t, y_dry*1000, 'o-', markevery=1000, color='black', ms=7,lw=2, label=r"numerical")
ax6.grid(True)
ax6.set_xlabel(r"$t \, [\text{s}]$")
ax6.set_ylabel(r"$y \, [\text{mm}]$")
ax6.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax6.set_xlim([0,3])
ax6.set_ylim([-2,0.1])
ax6.legend(loc='lower center', ncols=2)
plt.savefig("FIGB6_dry_harm_osz_vgl.pdf", dpi=100)



'''
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
