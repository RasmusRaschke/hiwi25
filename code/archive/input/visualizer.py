import numpy as np 
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description='Simulate the motion of a rolling ball with magnetic moment.')
parser.add_argument('--input', type=str, default='input', help='Path to the input file containing system parameters.')
args = parser.parse_args()

#Load input parameters
with open(args.input, 'r') as f:
    lines = f.readlines()
    M = float(lines[2].split()[0])
    R = float(lines[3].split()[0])
    g = float(lines[4].split()[0])
    phi = float(lines[5].split()[0])* np.pi / 180.0
    mag_mom_x = float(lines[6].split()[0])
    mag_mom_y = float(lines[7].split()[0])
    mag_mom_z = float(lines[8].split()[0])
    Bx = float(lines[9].split()[0])
    By = float(lines[10].split()[0])
    Bz = float(lines[11].split()[0])
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    mag_mom = np.sqrt(mag_mom_x**2 + mag_mom_y**2 + mag_mom_z**2)

TEST = 2 if phi == 0 else 1


def test_pos(x0, y0, z0, vx0, vy0, vz0, t):
    x = vx0 * t + x0
    y = (-5*g / 14) * np.sin(phi) * np.cos(phi) * t**2 + vy0 * t + y0
    z = (5*g / 14) * (1-np.cos(phi)) * t**2 + vz0 * t + z0
    return x, y, z

def get_pos(t):
    I = 2/5*M*R**2
    a = g*np.sin(phi)/(1+I/(M*R**2))
    y = 0.5*a*t**2
    return -y

def test_ang(Ox0, Oy0, Oz0, t):
    Ox = - (5*g / 7) * np.sin(phi) * t + Ox0
    Oy = Oy0
    Oz = Oz0
    return Ox, Oy, Oz

def test_osc_pos(t):
    I = 2/5 * M * R**2
    omega = np.sqrt(1/I * mag_mom*B)
    alpha_0 = np.arccos(mag_mom_z/mag_mom)
    return 2*R*(alpha_0*(np.sin(omega*t/2)**2))



data = np.genfromtxt('data.csv', delimiter=",", names=True, dtype=None, encoding=None)
headers = data.dtype.names
names = [n.strip() for n in headers]
cols = { clean: data[orig] for clean, orig in zip(names, headers)}
data_np = np.vstack([cols[n] for n in names]).T

t = cols['t']
x = cols['x'] 
y = cols['y'] 
z = cols['z']
Omx = cols['Ox'] 
Omy = cols['Oy'] 
Omz = cols['Oz']
qw = cols['q0'] 
qx = cols['q1'] 
qy = cols['q2'] 
qz = cols['q3']
vx = cols['vx'] 
vy = cols['vy'] 
vz = cols['vz']
T_trans = cols['T_trans'] 
T_rot = cols['T_rot']
U_gr = cols['U_gr'] 
U_em = cols['U_em'] 
E = cols['E']
qnorm = cols['q_norm']

fig, axs = plt.subplots(3,2, figsize=(10,10), constrained_layout=True)
axs[0,0].plot(t, x, label=r'$x \, [\text{m}]$') 
axs[0,0].plot(t, y, label=r'$y \, [\text{m}]$') 
if TEST == 1:
    xt = get_pos(t)
    axs[0,0].plot(t, xt, label=r"$\tilde{x}$", linestyle='dashed')
    # axs[0,0].plot(t, yt, label=r"$\tilde{y}$", linestyle='dashed')
elif TEST == 2:
    xt = test_osc_pos(t)
    axs[0,0].plot(t, xt, label=r"$\tilde{x}$", linestyle='dashed')

axs[0,0].set_title('Centre of Mass Motion')
axs[1,0].plot(t, vx, label=r'$v_x \, [\text{ms}^{-1}]$')
axs[1,0].plot(t, vy, label=r'$v_y \, [\text{ms}^{-1}]$') 
axs[1,0].set_title('Centre of Mass Velocities')
axs[2,0].plot(t, Omx, label=r'$\Omega_x \, [\text{s}^{-1}]$') 
axs[2,0].plot(t, Omy, label=r'$\Omega_y \, [\text{s}^{-1}]$') 
axs[2,0].plot(t, Omz, label=r'$\Omega_z \, [\text{s}^{-1}]$') 
axs[2,0].set_title(r'Angular Velocities')

axs[0,1].plot(x, y)
axs[0,1].set_title('Trajectory')
axs[1,1].plot(t, T_trans+T_rot, label=r'$T$')
axs[1,1].plot(t, U_gr+U_em, label=r'$U$') 
axs[1,1].set_title('Energy Transfer')
axs[2,1].plot(t, E)
axs[2,1].set_title('Total Energy')

for ax in axs.flat:
    ax.grid(True); ax.legend(fontsize='small')

plt.show()