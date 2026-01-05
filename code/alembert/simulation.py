'''
This code simulates a magnetic ball, modeled as an ideal dipole, rolling without slipping on an incline in a homogenous, constant magnetic field.
It solves the Lagrange-d'Alembert equations with a Lie group integrator.
'''
import numpy as np
import matplotlib.pyplot as plt
import os

#####---Constants---#####
M      = 8.15                           # mass of the ball
m      = 1090.0                         # magnetic dipole modulus
R      = 0.635                          # radius of the ball
alpha  = .2                           # slope of the incline 
g      = 981                            # gravitational acceleration
B      = np.array([0.0, 0.25, 0.1])    # magnetic field B=(B1,B2,B3)
mu     = np.array([m,0,0])              # magnetic dipole moment

#####---Initial Values---#####
x0, y0 = 0.0, 0.0                       # initial position
vx0, vy0 = 0.0, 0.0                     # initial linear velocity
s0 = 0.0                                # initial spin
q0 = np.array([1.0,0.0,0.0,0.0])        # initial orientation

#####---Integrator---#####
t0, t1 = 0.0, 10.0                      # start and end time
dt = .002                               # timestep
rk4 = True                             # Runge-Kutta 4 or simple differential?
save = True                            # Save plots?
output_dir = os.getcwd()                

#####---Quaternion Algebra---##### TODO: maybe make a class for quaternions?
def q_mult(q1, q2):
    # multipy two quaternions
    re_q1, im_q1 = q1[0], q1[1:4]
    re_q2, im_q2 = q2[0], q2[1:4]
    return np.array([re_q1 * re_q2 - np.dot(im_q1, im_q2), *(re_q1 * im_q2 + re_q2 * im_q1 + np.cross(im_q1, im_q2))])

def q_norm(q):
    # returns norm of quaternion q
    return np.sqrt(np.dot(q,q))

def q_exp(w):
    # returns exponential Lie map of w in so(3)
    angle = np.linalg.norm(w)
    if angle < 1e-12:
        return np.array([1.0,0.0,0.0,0.0])
    axis = w / angle
    return np.array([np.cos(angle), *(np.sin(angle) * axis)])

def q_rot(q):
    # do chart transition to rotation matrix
    q0, q1, q2, q3 = q
    return np.array([[1-2*(q2*q2+q3*q3), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1*q1+q2*q2)]])

#####---Constraints---#####
J = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, alpha]]) #Jacobian of center of mass
n = np.array([0.0, -alpha, 1.0])
n /= np.linalg.norm(n) # unit surface normal
I = (2 / 5) * M * R**2 # inertia tensor (for sphere)
W = np.zeros((3,2))
for i in range(2):
    W[:, i] = np.cross(n, J[:, i]) / R # auxilliary cross-product jacobian operator
M_eff = M * (J.T @ J) + I * (W.T @ W) # effective mass matrix
M_eff_inv = np.linalg.inv(M_eff)
Q_gr = -M * g * np.array([0.0, alpha])

#####---Energies---#####
def get_energies(state):
    x, y, vx, vy, s = state[:5]
    q = state[5:9]
    d_xi = np.array([vx, vy])
    d_r = J @ d_xi # center of mass velocity
    omega = W @ d_xi + s * n # angular velocity
    mu_lab = q_rot(q) @ mu # lab frame dipole (since B is in the lab frame)
    
    T_lin = 0.5 * M * np.dot(d_r, d_r)
    T_rot = 0.5 * I * np.dot(omega, omega)
    U_gr = M * g * alpha * y
    U_em = -np.dot(mu_lab, B)

    return(T_lin + T_rot + U_gr + U_em, T_lin, T_rot, U_gr, U_em)

#####---Integrator---#####
def integrate():
    steps = int(np.ceil((t1 - t0) / dt))
    x, y, vx, vy, s = x0, y0, vx0, vy0, s0
    q = q0 / q_norm(q0) # enforce unit quaternion

    motion, times, energies = [], [], []
    def get_constraint(vx, vy, s, q):
        omega = W @ np.array([vx, vy]) + s*n
        torque = np.cross(q_rot(q) @ mu, B)
        d_spin = np.dot(torque, n) / I
        Q_em = J.T @ (np.cross(torque, n) / R)
        spin_coupling = I * d_spin * (W.T @ n)
        Q = M_eff_inv @ (Q_gr + Q_em - spin_coupling)
        return Q[0], Q[1], d_spin, omega
    if rk4:
        for i in range(steps+1):
            motion.append([x, y, vx, vy, s, *q])
            times.append(t0 + i * dt)
            energies.append(get_energies(np.array([x, y, vx, vy, s, *q])))

            # Runge-Kutta 4
            Qx1, Qy1, d_spin1, omega1 = get_constraint(vx, vy, s, q)
            k1 = np.array([vx, vy, Qx1, Qy1, d_spin1])
            Qx2, Qy2, d_spin2, omega2 = get_constraint(vx + 0.5*dt*Qx1, vy + 0.5*dt*Qy1, s + 0.5*dt*d_spin1, q)
            k2 = np.array([vx + 0.5*dt*Qx1, vy + 0.5*dt*Qy1, Qx2, Qy2, d_spin2])
            Qx3, Qy3, d_spin3, omega3 = get_constraint(vx + 0.5*dt*Qx2, vy + 0.5*dt*Qy2, s + 0.5*dt*d_spin2, q)
            k3 = np.array([vx + 0.5*dt*Qx2, vy + 0.5*dt*Qy2, Qx3, Qy3, d_spin3])
            Qx4, Qy4, d_spin4, omega4 = get_constraint(vx + dt*Qx3, vy + dt*Qy3, s + dt*d_spin3, q)
            k4 = np.array([vx + dt*Qx3, vy + dt*Qy3, Qx4, Qy4, d_spin4])
            update = (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
            x += update[0]
            y += update[1]
            vx += update[2]
            vy += update[3]
            s += update[4]

            # Lie integrator
            omega_midpoint = W @ np.array([vx, vy]) + s * n
            dq = q_exp(0.5 * dt * omega_midpoint)
            q = q_mult(q, dq)
            q /= q_norm(q)
    else:
        for i in range(steps+1):
            motion.append([x, y, vx, vy, s, *q])
            times.append(t0 + i * dt)
            energies.append(get_energies(np.array([x, y, vx, vy, s, *q])))
            Q0, Q1, d_spin, omega = get_constraint(vx,vy,s,q)
            vx += dt * Q0
            vy += dt * Q1
            x += dt * vx
            y += dt * vy
            s += dt * d_spin
            omega = W @ np.array([vx,vy]) + s * n
            dq = q_exp(0.5 * dt * omega)
            q = q_mult(q, dq)
            q /= q_norm(q)


    return np.array(times), np.array(motion), np.array(energies)


#####---Run---#####

times, motion, energies = integrate()
x, y = motion[:,0], motion[:,1]
vx, vy = motion[:,2], motion[:, 3]
s = motion[:, 4]
T_lin, T_rot = energies[:, 1], energies[:, 2]
T = T_lin + T_rot
U_gr, U_em = energies[:, 3], energies[:, 4]
U = U_gr + U_em
E = energies[:, 0]

fig, axs = plt.subplots(3, 2, sharex=False, layout="constrained")
axs[0,0].plot(times, x, label=r"$x(t) \, [\text{cm}]$")
axs[0,0].plot(times, y, label=r"$y(t) \, [\text{cm}]$")
axs[0,0].set_title('Linear Motion')
axs[1,0].plot(times, vx, label=r"$v_x(t) \, [\text{cm}\text{s}^-1]$")
axs[1,0].plot(times, vy, label=r"$v_y(t) \, [\text{cm}\text{s}^-1]$")
axs[1,0].set_title('Linear Velocities')
axs[2,0].plot(times, s, label=r"$s(t) \, [\text{s}^-1]$")
axs[2,0].set_title('Spin')
axs[0,1].plot(x, y)
axs[0,1].set_title('Trajectory')
axs[1,1].plot(T, U)
axs[1,1].set_title('Energy Phase Space')
axs[2,1].plot(times, E, label=r"$E(t) \, [\text{erg}]$")
axs[2,1].set_title("Total Energy")
for ax in axs.flat:
    ax.legend()
    ax.grid(True)
if save:
    plt.savefig(os.path.join(output_dir, 'data.png'), dpi=300)
plt.show()

'''
axs[0,0].set_ylim(-10,10)
axs[1,0].set_ylim(0, 1)
axs[2,0].set_ylim(-10, 10)
axs[0,1].plot(motion[:, 1], motion[:, 0])
axs[0,1].set_xlabel(r"$y (t) \, [\text{cm}]$")
axs[0,1].set_ylabel(r"$x (t) \, [\text{cm}]$")
for ax in axs.flat:
    ax.legend()
plt.show()
'''