"""Self-convergence test for magsphere

"""

from pathlib import Path
import argparse
import csv
import subprocess
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
EXE = Path('./magsphere.out').resolve()
INPUT = Path('input.base')
H = np.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001])
H_REF = 1e-6
SAVE_DT = 1e-3   # output cadence assumed in the patched C++ code
ALPHAS = [5.0, 10.0, 15.0]

# Add/remove magnetic cases here.
# Each entry is: name: (mu-vector, B-vector, marker)
FIELD_CASES = {
    'B=2T ez, mu=5e-5 ex': ((5e-5, 0.0, 0.0), (0.0, 0.0, 2.0), 'o'),
    'B=5e-5T ex, mu=2 ez': ((0.0, 0.0, 2.0), (5e-5, 0.0, 0.0), 's'),
    # example low-field variants you can activate if wanted:
    # 'B=0.02T ez, mu=5e-5 ex': ((5e-5, 0.0, 0.0), (0.0, 0.0, 0.02), '^'),
    # 'B=5e-7T ex, mu=2 ez': ((0.0, 0.0, 2.0), (5e-7, 0.0, 0.0), 'd'),
}

# color per alpha; extend if you add more alpha values
ALPHA_COLORS = {
    5.0: 'tab:gray',
    10.0: 'tab:blue',
    15.0: 'tab:red',
}

T_END = 1.0
OUT = Path(f'convergence_results_{T_END:g}s')
OUT.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# STYLE
# -----------------------------------------------------------------------------
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{siunitx} \usepackage{bm}",
    "font.size": 24,
    "axes.titlesize": 24,
    "axes.labelsize": 24,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 24,
})

# -----------------------------------------------------------------------------
# COMPARE
# -----------------------------------------------------------------------------
base = [s for line in INPUT.read_text().splitlines()
        if (s := line.split('#')[0].strip())]


def simulate(alpha_deg, case_name, h, mu, B):
    """Run the compiled solver and load:
      t, x, y, z, Ox, Oy, Oz, q0, q1, q2, q3
    """
    p = base.copy()
    p[0] = str(h)
    p[1] = str(T_END)
    p[5] = str(alpha_deg)
    p[6:9] = map(str, mu)
    p[9:12] = map(str, B)

    folder = OUT / f'alpha_{alpha_deg:g}' / case_name.replace(' ', '_').replace(',', '')
    folder.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(EXE)], input='\n'.join(p) + '\n', text=True,
                   cwd=folder, stdout=subprocess.DEVNULL, check=True)

    target = folder / f'dt_{h:.9g}.csv'
    (folder / 'data.csv').replace(target)
    return np.loadtxt(target, delimiter=',', skiprows=1, usecols=range(11), ndmin=2)


def so3_error(q, qref):
    """Geodesic distance on SO(3). q and -q represent the same rotation."""
    q = q / np.linalg.norm(q, axis=1)[:, None]
    qref = qref / np.linalg.norm(qref, axis=1)[:, None]
    inner = np.abs(np.sum(q * qref, axis=1))
    inner = np.clip(inner, 0.0, 1.0)
    return 2.0 * np.arccos(inner)


def trajectory_errors(a, ref):
    """Compare a coarse run against the reference at matching saved times."""
    # reference is sampled every SAVE_DT starting at t = SAVE_DT
    idx = np.rint(a[:, 0] / SAVE_DT).astype(int) - 1
    ref_rows = ref[idx]

    dr = np.linalg.norm(a[:, 1:4] - ref_rows[:, 1:4], axis=1)
    domega = np.linalg.norm(a[:, 4:7] - ref_rows[:, 4:7], axis=1)
    dquat = so3_error(a[:, 7:11], ref_rows[:, 7:11])
    return dr, domega, dquat

raw = {}
results = {}
rows = []

for alpha in ALPHAS:
    raw[alpha] = {}
    results[alpha] = {}
    for case_name, (mu, B, marker) in FIELD_CASES.items():
        raw[alpha][case_name] = {}
        for h in [*H, H_REF]:
            raw[alpha][case_name][h] = simulate(alpha, case_name, h, mu, B)

        ref = raw[alpha][case_name][H_REF]
        arr = np.zeros((3, len(H)))
        for j, h in enumerate(H):
            dr, domega, dquat = trajectory_errors(raw[alpha][case_name][h], ref)
            arr[:, j] = [dr[-1], domega[-1], dquat[-1]]
            rows.append([alpha, case_name, h, dr[-1], domega[-1], dquat[-1]])
        results[alpha][case_name] = arr

with (OUT / 'summary.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['alpha_deg', 'case', 'dt',
                'position_final_m', 'omega_final_rad_s', 'orientation_final_rad'])
    w.writerows(rows)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
ylabels = [
    r'$\|\Delta_t \mathbf{r}\|\ [\mathrm{m}]$',
    r'$\|\Delta_t \boldsymbol{\Omega}\|\ [\unit{\radian\per\second}]$',
    r'$\|\Delta_t \boldsymbol{q}\|\ [\mathrm{rad}]$'
]
letters = ['(a)', '(b)', '(c)']

# Keep the plotting area at 16 x 5 inches. Reserve only 1.45 inches
# above it for the two-row legend.
fig = plt.figure(figsize=(16, 6.45), layout='constrained')
legend_area, plot_area = fig.subfigures(2, 1, height_ratios=[1.45, 5], hspace=0)
axes = plot_area.subplots(1, 3)

# one O(dt^2) guide per panel, anchored to the first available finest-step datum
for k, ax in enumerate(axes):
    anchor = None
    for alpha in ALPHAS:
        for case_name in FIELD_CASES:
            y = results[alpha][case_name][k]
            if np.all(y > 0):
                anchor = y[-1] / (H[-1] ** 2 - H_REF ** 2)
                break
        if anchor is not None:
            break
    if anchor is not None:
        ax.loglog(H, anchor * (H ** 2 - H_REF ** 2), 'k--', lw=1.8,
                  label=r'$\Delta t^2$')

    for alpha in ALPHAS:
        color = ALPHA_COLORS.get(alpha, None)
        for case_name, (_, _, marker) in FIELD_CASES.items():
            y = results[alpha][case_name][k]
            ax.loglog(H, y, marker=marker, linestyle='-', color=color,
                      lw=2.0, ms=6,
                      label=rf'$\alpha={alpha:g}^\circ$: {case_name}')

    ax.set_xlabel(r'$\Delta t\ [\mathrm{s}]$')
    ax.set_ylabel(ylabels[k])
    ax.set_xticks([1e-5, 1e-3, 1e-1])
    ax.grid(True, which='both', alpha=0.3)
    ax.text(0.01, 0.98, letters[k], transform=ax.transAxes,
            ha='left', va='top', fontsize=20)

# Single boxed legend with two rows. Matplotlib fills columns first, so
# interleave the entries: [top-left, bottom-left, top-middle, bottom-middle, ...].
case_markers = [case[2] for case in FIELD_CASES.values()]
upper_handles = [
    Line2D([], [], color='black', linestyle='--', lw=1.8),
    Line2D([], [], color='black', marker=case_markers[0], linestyle='None', ms=10),
    Line2D([], [], color='black', marker=case_markers[1], linestyle='None', ms=10),
]
upper_labels = [
    r'$\Delta t^2$',
    r'$B_z=2\,\mathrm{T},\ \mu_x=5\cdot10^{-5}\,\mathrm{A\,m^2}$',
    r'$B_x=5\cdot10^{-5}\,\mathrm{T},\ \mu_z=2\,\mathrm{A\,m^2}$',
]
lower_handles = [Line2D([], [], color=ALPHA_COLORS[a], lw=2) for a in ALPHAS]
lower_labels = [rf'$\varphi={a:g}^\circ$' for a in ALPHAS]

handles = [item for pair in zip(upper_handles, lower_handles) for item in pair]
labels = [item for pair in zip(upper_labels, lower_labels) for item in pair]
legend_area.legend(handles, labels, ncol=3, loc='center',
                   frameon=True, fancybox=True, framealpha=1,
                   edgecolor='black', fontsize=20, borderpad=0.3, labelspacing=0.35,
                   columnspacing=1.2, handletextpad=0.5)

fig.savefig(OUT / 'convergence.pdf')
fig.savefig(OUT / 'convergence.png')
plt.close(fig)