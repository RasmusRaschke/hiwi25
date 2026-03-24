import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Simulate the motion of a rolling ball with magnetic moment.')
parser.add_argument('--input', type=str, default='input', help='Path to the input file containing system parameters.')
args = parser.parse_args()

with open(args.input, 'r') as f:
    lines = f.readlines()
    t_max = float(lines[1].split()[0])
    M = float(lines[2].split()[0])
    R = float(lines[3].split()[0])
    g = float(lines[4].split()[0])
    phi = float(lines[5].split()[0])*np.pi/180 # convert to radians
    mag_mom_x = float(lines[6].split()[0])
    mag_mom_y = float(lines[7].split()[0])
    mag_mom_z = float(lines[8].split()[0])
    Bx = float(lines[9].split()[0])
    By = float(lines[10].split()[0])
    Bz = float(lines[11].split()[0])
    start_spin = float(lines[12].split()[0]) # in degrees per second
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    mag_mom = np.sqrt(mag_mom_x**2 + mag_mom_y**2 + mag_mom_z**2)

def get_n(phi):
    return np.array([np.sin(phi), 0, np.cos(phi)])

def get_e_x(phi):
    return np.array([np.cos(phi), 0, -np.sin(phi)])

def get_e_y(phi):
    return np.array([0.0, 1.0, 0.0])

def magnetic_moment(alpha,beta,mu):
    return mu*np.array([np.sin(alpha)*np.cos(beta),np.sin(alpha)*np.sin(beta),np.cos(alpha)])

def f(t,x,system_para,epsilon= 0.01):
    R                   = system_para[0]
    phi                 = system_para[1]
    m                   = system_para[2]
    g                   = system_para[3]
    eta                 = system_para[4]
    t_max               = system_para[5]
    mu                  = system_para[6] # ergs/gauss
    B_x                 = system_para[7] # gauss
    B_y                 = system_para[8] # gauss
    B_z                 = system_para[9] # gauss

    I = 2/5*m*R**2

    x_coor,y_corr,alpha,beta,L_x,L_y,L_z = x
    L = np.array([L_x,L_y,L_z])
    B = np.array([B_x,B_y,B_z])

    grav_force = -m*g*np.array([0,0,1])
    n = get_n(phi)
    e_x = get_e_x(phi)
    e_y = get_e_y(phi)

    n_grav_force = grav_force-np.dot(grav_force,n)*n

    grav_torque = -R*np.cross(grav_force,n)

    mag_torque =  - np.cross(magnetic_moment(alpha,beta,mu),B)

    omega = L/(I+m*R**2)

    normal_gravity = -np.dot(grav_force,n)


    total_torque = grav_torque + mag_torque


    dx = np.dot(np.cross(R*n,omega),e_x)
    dy = np.dot(np.cross(R*n,omega),e_y)


    dalpha = omega[1]*np.cos(beta)-omega[0]*np.sin(beta)
    dbeta = omega[2] - (1/(np.tan(alpha)))*(omega[0]*np.cos(beta)-omega[1]*np.sin(beta))

    return np.array([dx,dy,dalpha,dbeta,total_torque[0],total_torque[1],total_torque[2]])

def plot_trajectory_with_angles(sol,system_para, plot_analytical_solution=True, start_alpha=None):
    R, phi, m, g, eta, t_max, mu, Bx, By, Bz = system_para
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    I = 2/5*m*R**2
    a = g*np.sin(phi)/(1+I/(m*R**2))
    axs[0].plot(sol.t,sol.y[0],linewidth=10, label="Numerics")
    if plot_analytical_solution:
        omega = np.sqrt(1/(I+m*R**2)*mu*Bz)
        if start_alpha is None:
            axs[0].plot(sol.t,-0.5*a*sol.t**2,linewidth=3,linestyle="--",label="Analytical solution")
        else:
            axs[0].plot(sol.t, 2*R*((np.pi-start_alpha)*(np.sin(omega*sol.t/2)**2)),linewidth=3,linestyle="--",label="Analytical solution")
    axs[0].set_xlabel("t in s")
    axs[0].set_ylabel("x in cm")
    axs[0].legend()

    a = np.arccos(np.cos(sol.y[2]))
    scatter = axs[1].scatter(sol.y[1],sol.y[0],c=a,cmap="hsv",vmin=0,vmax=np.pi)
    max_v = np.max([np.abs(sol.y[0][0]-sol.y[0][-1]),np.abs(sol.y[1][0]-sol.y[1][-1])])
    axs[1].set_xlim(-max_v,max_v)
    axs[1].set_ylim(-max_v,max_v)
    axs[1].set_xlabel("y in cm")
    axs[1].set_ylabel("x in cm")
    axs[1].set_aspect('equal', adjustable='box')
    cbar = fig.colorbar(scatter, ax=axs[1], ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
    cbar.set_label('$\\alpha$')
    cbar.ax.set_yticklabels(["$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])

    scatter = axs[2].scatter(sol.y[1],sol.y[0],c=sol.y[3],cmap="hsv",vmin=0,vmax=2*np.pi)
    max_v = np.max([np.abs(sol.y[0][0]-sol.y[0][-1]),np.abs(sol.y[1][0]-sol.y[1][-1])])
    axs[2].set_xlim(-max_v,max_v)
    axs[2].set_ylim(-max_v,max_v)
    axs[2].set_xlabel("y in cm")
    axs[2].set_ylabel("x in cm")
    axs[2].set_aspect('equal', adjustable='box')
    cbar = fig.colorbar(scatter, ax=axs[2], ticks=[0,2*np.pi/4,2*np.pi/2,2*3*np.pi/4,2*np.pi])
    cbar.set_label('$\\beta$')
    cbar.ax.set_yticklabels(["$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":


    M = M*1000 # convert to grams
    R = R*100 # convert to cm
    g = g*100 # convert to cm/s^2
    mu = mag_mom*1e3 # convert to ergs/gauss
    alpha = np.pi-np.arccos(mag_mom_z/mag_mom)
    if mag_mom != 0.0:
        alpha = np.pi - np.arccos(mag_mom_z/mag_mom)
    else:
        alpha = 0.0
    beta = np.arctan2(mag_mom_y,mag_mom_x) - np.pi
    Bx = Bx*1e4 # convert to gauss
    By = By*1e4 # convert to gauss
    Bz = Bz*1e4 # convert to gauss
    B = np.sqrt(Bx**2 + By**2 + Bz**2)

    #                       R     phi  m     g         eta      t_max    mu      Bx   By   Bz
    system_para = np.array([R, phi,    M,    g,        0.0,      t_max,   mu,  Bx,   By,   Bz])
    t_max = system_para[5]
    n_eval = 1000
    # Add a spin to the ball (Angular velocity vector perpendicular to the inclined plane) 
    n = get_n(phi)
    spin_omega = start_spin/180.0*np.pi # convert to radians per second
    I = 2/5*M*R**2
    L = I/(I+M*R**2)*spin_omega*n

    x_0 = np.array([0.0,0.0,alpha,beta,L[0],L[1],L[2]])
    sol = solve_ivp(f,[0,t_max],x_0,args=(system_para,),t_eval=np.linspace(0,t_max,n_eval),method='Radau',rtol=1e-5,atol=1e-5)

    plot_trajectory_with_angles(sol,system_para,start_alpha=alpha if phi == 0 else None)
