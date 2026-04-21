import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import simutils as su
import matplotlib.ticker as mticker
import re
from collections import defaultdict
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

data = su.extract()
t = data["data_1"].t
mask = t <= 1
x = [data["data_1"].x, data["data_10"].x, data["data_100"].x, data["data_1000"].x, data["data_0"].x]
y = [data["data_1"].y, data["data_10"].y, data["data_100"].y, data["data_1000"].y, data["data_0"].y]
z = [data["data_1"].z, data["data_10"].z, data["data_100"].z, data["data_1000"].z, data["data_0"].z]
titles = [r"$\phi = \frac{\pi}{2} \, [\text{rad}]$", r"$\phi = \frac{\pi}{20} \, [\text{rad}]$", r"$\phi = \frac{\pi}{200} \, [\text{rad}]$", r"$\phi = \frac{\pi}{2000} \, [\text{rad}]$", r"$\phi = 0 \, [\text{rad}]$"]
cmap = plt.cm.plasma
colors = cmap(np.linspace(0,1,5))
fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
#ax.plot(x[0], y[0], color=colors[0], lw=1.5, label=titles[0]) 
#ax.plot(x[1], y[1], color=colors[1], lw=1.5, label=titles[1]) 
ax.plot(x[2], y[2], color=colors[2], lw=1.5, label=titles[2]) 
ax.plot(x[3], y[3], color=colors[3], lw=1.5, label=titles[3]) 
ax.plot(x[4], y[4], color=colors[4], lw=1.5, label=titles[4]) 
ax.set_xlabel(r"$x \, [\text{m}]$")
ax.set_ylabel(r"$y \, [\text{m}]$")
ax.legend()
ax.grid(True, alpha=0.3)
plt.savefig("FIGI_comp_dev.pdf", dpi=300)