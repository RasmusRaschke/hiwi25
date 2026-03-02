import numpy as np 
import matplotlib.pyplot as plt 

TEST = 0
g = 9.80665
phi = 20.0 * np.pi / 180.0

def test_pos(x0, y0, z0, vx0, vy0, vz0, t):
    x = vx0 * t + x0
    y = (-5*g / 14) * np.sin(phi) * np.cos(phi) * t**2 + vy0 * t + y0
    z = (5*g / 14) * (1-np.cos(phi)) * t**2 + vz0 * t + z0
    return x, y, z


def test_ang(Ox0, Oy0, Oz0, t):
    Ox = - (5*g / 7) * np.sin(phi) * t + Ox0
    Oy = Oy0
    Oz = Oz0
    return Ox, Oy, Oz


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
if TEST:
    xt, yt, zt = test_pos(0,0,0,0,0,0,t)
    axs[0,0].plot(t, xt, label=r"$\tilde{x}$")
    axs[0,0].plot(t, yt, label=r"$\tilde{y}$")
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