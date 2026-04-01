import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import simutils as su
from mpl_toolkits.mplot3d import axes3d
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

datasets = su.extract("mag_moment_data")
slices_y = [-1.0, -0.3, 0.3, 1.0]
fname_re = re.compile(r"data_mx_(?P<mx>.+)_my_(?P<my>.+)$")
selected = {my: [] for my in slices_y}
all_mx = []
for name, ds in datasets.items():
    m = fname_re.fullmatch(name)
    if m is None:
        continue
    mx = su.decode(m.group("mx"))
    my = su.decode(m.group("my"))
    if my in selected:
        selected[my].append((mx, ds.x, ds.y))
        all_mx.append(mx)
for my in selected:
    selected[my].sort(key=lambda t: t[0])

all_mx = np.array(all_mx)
norm = Normalize(vmin=all_mx.min(), vmax=all_mx.max())
cmap = plt.cm.plasma

fig3, axs = plt.subplots(2, 2, figsize=(14, 14), layout="constrained")
axs = np.atleast_1d(axs).ravel()
for ax, my in zip(axs, slices_y):
    for mx, x, y in selected[my]:
        ax.plot(x*100, y*100, color=cmap(norm(mx)), lw=1.2)
    ax.set_title(fr"fixed $\mu_y = {my:g}$")
    ax.set_xlabel(r"$x \, [\text{cm}]$")
    ax.set_ylabel(r"$y \, [\text{cm}]$")
    ax.set_box_aspect(1)
    ax.grid(True, alpha=0.3)

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig3.colorbar(sm, ax=axs, location="right", pad=0.02)
cbar.set_label(r"$\mu_x \,[\text{Am}^2]$")
plt.savefig("FIGD3_mag_mom_traj.pdf", dpi=300)