import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import simutils as su
from mpl_toolkits.mplot3d import axes3d
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
"""
mz, y = su.batch_extract_mz("mag_moment_data")
np.savez("mag_moment_grid.npz", mz=mz, y=y)
"""
data = np.load("mag_moment_grid.npz")
mz = data["mz"]
y = data["y"]

datasets = su.extract()
t = datasets["sample1"].t
y1 = datasets["sample1"].y
y2 = datasets["sample2"].y
y3 = datasets["sample3"].y
y4 = datasets["sample4"].y
y5 = datasets["sample5"].y

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
order = np.argsort(mz)
ax1.plot(mz[order], y[order]*100, lw=3.5, color="black")
ax1.set_xlabel(r"$\mu_z \, [\text{Am}^2]$")
ax1.set_ylabel(r"$y \, [\text{cm}]$")
ax1.set_xlim(0,1)
ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
ax1.grid(True)
ax2.plot(t, y1*100, label=r"$\mu_z = 0.0$", lw=3.5, color="blue")
ax2.plot(t, y2*100, label=r"$\mu_z = 0.95$", lw=3.5, color="orange")
#ax2.plot(t, y3*100, label=r"$\mu_z = 0.96$", lw=2.5) both already oscillatory
#ax2.plot(t, y4*100, label=r"$\mu_z = 0.97$", lw=2.5)
ax2.plot(t, y5*100, label=r"$\mu_z = 1.0$", lw=3.5, color="red")
ax2.set_xlabel(r"$t \, [\text{s}]$")
ax2.set_ylabel(r"$y \, [\text{cm}]$")
ax2.grid(True)
ax2.legend()
ax2.set_xlim(0,1)
ax2.set_ylim(-10,1)
plt.savefig("FIGE_mag_mom_z.pdf", dpi=300)
#"""
