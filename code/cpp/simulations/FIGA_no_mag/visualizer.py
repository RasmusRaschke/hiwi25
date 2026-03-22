import numpy as np 
import matplotlib.pyplot as plt 
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

def test_pos(x0, y0, z0, vx0, vy0, vz0, p, t):
    phi = p * np.pi / 180.0
    x = vx0 * t + x0
    y = (-5*g / 14) * np.sin(phi) * np.cos(phi) * t**2 + vy0 * t + y0
    z = (5*g / 14) * (1-np.cos(phi)) * t**2 + vz0 * t + z0
    return x, y, z


def test_ang(Ox0, Oy0, Oz0, p, t):
    phi = p * np.pi / 180.0
    Ox = - (5*g / 7) * np.sin(phi) * t + Ox0
    Oy = Oy0
    Oz = Oz0
    return Ox, Oy, Oz

#Comapre analytical solutions without magnetic field
data = np.genfromtxt('data5.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
t5 = cols['t']; x5 = cols['x']; y5 = cols['y']; z5 = cols['z']
Omx5 = cols['Ox']; Omy5 = cols['Oy']; Omz5 = cols['Oz']
T_trans5 = cols['T_trans']; T_rot5 = cols['T_rot']; U_gr5 = cols['U_gr']; U_em5 = cols['U_em']; E5 = cols['E']

data = np.genfromtxt('data10.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
t10 = cols['t']; x10 = cols['x']; y10 = cols['y']; z10 = cols['z']
Omx10 = cols['Ox']; Omy10 = cols['Oy']; Omz10 = cols['Oz']
T_trans10 = cols['T_trans']; T_rot10 = cols['T_rot']; U_gr10 = cols['U_gr']; U_em10 = cols['U_em']; E10 = cols['E']

data = np.genfromtxt('data15.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
t15 = cols['t']; x15 = cols['x']; y15 = cols['y']; z15 = cols['z']
Omx15 = cols['Ox']; Omy15 = cols['Oy']; Omz15 = cols['Oz']
T_trans15 = cols['T_trans']; T_rot15 = cols['T_rot']; U_gr15 = cols['U_gr']; U_em15 = cols['U_em']; E15 = cols['E']

###############################friction##########################
data = np.genfromtxt('data10c.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x10c = cols['x']; y10c = cols['y']; z10c = cols['z']
Omx10c = cols['Ox']; Omy10c = cols['Oy']; Omz10c = cols['Oz']
T_trans10c = cols['T_trans']; T_rot10c = cols['T_rot']; U_gr10c = cols['U_gr']; U_em10c = cols['U_em']; E10c = cols['E']

data = np.genfromtxt('data10v.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x10v = cols['x']; y10v = cols['y']; z10v = cols['z']
Omx10v = cols['Ox']; Omy10v = cols['Oy']; Omz10v = cols['Oz']
T_trans10v = cols['T_trans']; T_rot10v = cols['T_rot']; U_gr10v = cols['U_gr']; U_em10v = cols['U_em']; E10v = cols['E']

data = np.genfromtxt('data10a.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T
x10a = cols['x']; y10a = cols['y']; z10a = cols['z']
Omx10a = cols['Ox']; Omy10a = cols['Oy']; Omz10a = cols['Oz']
T_trans10a = cols['T_trans']; T_rot10a = cols['T_rot']; U_gr10a = cols['U_gr']; U_em10a = cols['U_em']; E10a = cols['E']

x_deg5, y_deg5, z_deg5 = test_pos(0,0,0,0,0,0,5,t5)
x_deg10, y_deg10, z_deg10 = test_pos(0,0,0,0,0,0,10,t10)
x_deg15, y_deg15, z_deg15 = test_pos(0,0,0,0,0,0,15,t15)

'''
fig1, ax1 = plt.subplots(layout='constrained')
delta_x = np.absolute(np.subtract(x10, x_deg10)) * 1000
delta_y = np.absolute(np.subtract(y10, y_deg10)) * 1000
ax1.plot(t5, delta_x, color='blue', label=r"$\Delta x$")
ax1.plot(t5, delta_y, color='red', label=r"$\Delta y$")
ax1.grid(True)
ax1.set_xlabel(r"$t \, [\text{s}]$")
ax1.set_ylabel(r"$\Delta r \, [\text{mm}]$")
ax1.tick_params(which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True, labelright=False, labeltop=False)
ax1.set_xlim([0,10])
#ax1.set_ylim([-60,1])
ax1.legend()
plt.savefig("FIG1_no_mag_diff.pdf", dpi=100)
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
