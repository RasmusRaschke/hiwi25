import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

import odes as ode

###Parameters###
R = 0.635  # radius
M = 10  # mass
m = 1090.0  # dipole moment
alpha = -1  # slope
g = 981.0  # gravitational acceleration
B1, B2, B3 = .0, .1, 0.0  # earth magnetic field
A1, A2 = 0.0, 0.0  # constants of integration
t0, t1, dt = .0, 1.0, .00001 #time in seconds


###Quaternion Normalizer, input quaternion as numpy array###
def norm(quat):
    qnorm = np.linalg.norm(quat)
    if qnorm == 0:
        raise ValueError("Zero quaternion not on S3!")
    return quat / qnorm


###Initial Data with normalized quaternion###
pos0 = np.array([0.0, 0.0])
quat0 = np.array([1.0, 0.0, 0.0, 0.0])
quat0 = norm(quat0)
vel0 = np.array([0.0, 0.01, 0.0, 0.0])
state0 = np.concatenate((pos0, quat0, vel0))

    
def equations(t, X, Y, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2):
    f = np.zeros(10)
    f[0] = ode.dX(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    f[1] = ode.dY(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    f[2] = Vq0
    f[3] = Vq1
    f[4] = Vq2
    f[5] = Vq3
    f[6] = ode.ddQ0(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    f[7] = ode.ddQ1(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    f[8] = ode.ddQ2(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    f[9] = ode.ddQ3(t, Q0, Q1, Q2, Q3, Vq0, Vq1, Vq2, Vq3, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
    return f



def rk4(state0, t0, t1, dt, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2):
    time = np.arange(t0, t1+dt, dt)
    state = np.zeros((len(time), len(state0)))
    state[0] = state0
    for i in range(1, len(time)):
        print("Step:", i)
        t = time[i-1]
        f = state[i-1]
        k1 = equations(t, *f, radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
        k2 = equations(t + .5*dt, *(f + .5*dt*k1), radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
        k3 = equations(t + .5*dt, *(f + .5*dt*k2), radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
        k4 = equations(t + dt, *(f + dt*k3), radius, mass, dipole, slope, grav, field1, field2, field3, const1, const2)
        f_step = f + (dt / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        q = f_step[2:6]
        q /= np.linalg.norm(q)
        f_step[2:6] = q
        '''
        v = f_step[6:10]
        v -= np.dot(q, v) * q
        f_step[6:10] = v
        '''
        state[i] = f_step
    return time, state

time, motion = rk4(state0, t0, t1, dt, R, M, m, alpha, g, B1, B2, B3, A1, A2)

limit = 0
threshold = 10000000000
mask = np.all(np.abs(motion)<threshold, axis=1)
t_reduced = time[mask]
motion_reduced = motion[mask]

if limit:
    time = t_reduced
    motion = motion_reduced
fig, axs = plt.subplots(3, 2, sharex=True)
axs[0,0].plot(time, motion[:, 0], label=r"$x(t) \, [\text{s}]$")
axs[0,0].plot(time, motion[:, 1], label=r"$y(t) \, [\text{s}]$")
axs[1,0].plot(time, motion[:, 2], label=r"$q_0(t) \, [\text{s}]$")
axs[1,0].plot(time, motion[:, 3], label=r"$q_1(t) \, [\text{s}]$")
axs[1,0].plot(time, motion[:, 4], label=r"$q_2(t) \, [\text{s}]$")
axs[1,0].plot(time, motion[:, 5], label=r"$q_3(t) \, [\text{s}]$")
axs[2,0].plot(time, motion[:, 6], label=r"$\dot{q}_0(t) \, [\text{cm}^{-1}\text{s}]$")
axs[2,0].plot(time, motion[:, 7], label=r"$\dot{q}_1(t) \, [\text{cm}^{-1}\text{s}]$")
axs[2,0].plot(time, motion[:, 8], label=r"$\dot{q}_2(t) \, [\text{cm}^{-1}\text{s}]$")
axs[2,0].plot(time, motion[:, 9], label=r"$\dot{q}_3(t) \, [\text{cm}^{-1}\text{s}]$")
axs[0,0].set_ylim(-10,10)
axs[1,0].set_ylim(0, 1)
axs[2,0].set_ylim(-10, 10)
axs[0,1].plot(motion[:, 1], motion[:, 0])
axs[0,1].set_xlabel(r"$y (t) \, [\text{cm}]$")
axs[0,1].set_ylabel(r"$x (t) \, [\text{cm}]$")
for ax in axs.flat:
    ax.legend()
plt.show()

def quat_to_vec(q):
    q0, q1, q2, q3 = q
    x = 1 - 2*(q2**2 + q3**2)
    y = 2*(q1*q2 + q0*q3)
    z = 2*(q1*q3 - q0*q2)
    vec = np.array([x, y, z])
    vec /= np.linalg.norm(vec) + 1e-16
    return vec

x = motion[:, 0]
y = motion[:, 1]
z = np.tan(np.radians(alpha)) * x
quats = motion[:, 2:6]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Xsurf = np.linspace(np.min(x), np.max(x), 10)
Ysurf = np.linspace(np.min(y), np.max(y), 10)
Xsurf, Ysurf = np.meshgrid(Xsurf, Ysurf)
Zsurf = np.tan(np.radians(alpha)) * Xsurf
ax.plot_surface(Xsurf, Ysurf, Zsurf, alpha=0.2, color='gray')
dot, = ax.plot([], [], [], 'ro', markersize=5)
arrow = ax.quiver(0,0,0,0,0,0,length=1,color='blue')
ax.set_xlim(np.min(x)-1, np.max(x)+1)
ax.set_ylim(np.min(y)-1, np.max(y)+1)
ax.set_zlim(-1, 1)
ax.set_xlabel(r"x")
ax.set_ylabel(r"y")
ax.set_zlabel(r"z")
def update(i):
    dot.set_data([x[i]], [y[i]])
    dot.set_3d_properties([z[i]])
    for coll in ax.collections[:]:
        coll.remove()
    vec = quat_to_vec(quats[i])
    ax.quiver(x[i], y[i], z[i], vec[0], vec[1], vec[2], length=0.5, color='blue')
    return dot,
ani = FuncAnimation(fig, update, frames=len(x), interval=10, blit=False)
plt.show()

