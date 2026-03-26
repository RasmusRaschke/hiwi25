import numpy as np 
import matplotlib.pyplot as plt
import simutils as su

plt.style.use('seaborn-v0_8-paper')

g = 9.80665

datasets = su.extract() 

t = datasets["data"].t
fig, axs = plt.subplots(3,2, figsize=(10,10), constrained_layout=True)
axs[0,0].plot(t, datasets["data"].x, label=r'$x \, [\text{m}]$') 
axs[0,0].plot(t, datasets["data"].y, label=r'$y \, [\text{m}]$') 
axs[0,0].plot(t, datasets["data"].z, label=r'$z \, [\text{m}]$') 
axs[0,0].set_title('Centre of Mass Motion')
axs[1,0].plot(t, datasets["data"].vx, label=r'$v_x \, [\text{ms}^{-1}]$')
axs[1,0].plot(t, datasets["data"].vy, label=r'$v_y \, [\text{ms}^{-1}]$') 
axs[1,0].plot(t, datasets["data"].vz, label=r'$v_z \, [\text{ms}^{-1}]$') 
axs[1,0].set_title('Centre of Mass Velocities')
axs[2,0].plot(t, datasets["data"].Ox, label=r'$\Omega_x \, [\text{s}^{-1}]$') 
axs[2,0].plot(t, datasets["data"].Oy, label=r'$\Omega_y \, [\text{s}^{-1}]$') 
axs[2,0].plot(t, datasets["data"].Oz, label=r'$\Omega_z \, [\text{s}^{-1}]$') 
axs[2,0].set_title(r'Angular Velocities')
axs[0,1].plot(datasets["data"].x, datasets["data"].y)
axs[0,1].set_title('Trajectory')
axs[1,1].plot(t, datasets["data"].T_trans + datasets["data"].T_rot, label=r'$T$')
axs[1,1].plot(t, datasets["data"].U_gr + datasets["data"].U_em, label=r'$U$') 
axs[1,1].set_title('Energy Transfer')
axs[2,1].plot(t, datasets["data"].E)
axs[2,1].set_title('Total Energy')
for ax in axs.flat:
    ax.grid(True); ax.legend(fontsize='small')
plt.show()