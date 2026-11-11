"""Effect of tanh smoothing k on the rolling-sphere trajectory.
"""

from pathlib import Path
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
DT = 1e-5
T_END = 1.0
K0 = 100.0  # nominal, physically chosen smoothing parameter [s/m]
K_VALUES = np.array([10, 20, 30, 50, 75, 100, 125, 150, 170, 180, 200, 500, 750, 1000], dtype=float)
ALPHAS = [5.0, 10.0, 15.0]

# name: (magnetic moment, magnetic field, marker)
FIELD_CASES = {
    'B=2T ez, mu=5e-5 ex': ((5e-5, 0.0, 0.0), (0.0, 0.0, 2.0), 'o'),
    'B=5e-5T ex, mu=2 ez': ((0.0, 0.0, 2.0), (5e-5, 0.0, 0.0), 's'),
    # 'B=0.02T ez, mu=5e-5 ex': ((5e-5, 0, 0), (0, 0, 0.02), '^'),
    # 'B=5e-7T ex, mu=2 ez': ((0, 0, 2.0), (5e-7, 0, 0), 'd'),
}
ALPHA_COLORS = {
    5.0: 'tab:gray',
    10.0: 'tab:blue',
    15.0: 'tab:red',
}

OUT = Path(f'k_sensitivity_results_{T_END:g}s')
OUT.mkdir(exist_ok=True)

# Exact style block from figA.py.
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

base = [s for line in INPUT.read_text().splitlines()
        if (s := line.split('#')[0].strip())]


def simulate(alpha, case_name, k, mu, B):
    p = base.copy()
    p[0], p[1], p[5] = str(DT), str(T_END), str(alpha)
    p[6:9] = map(str, mu)
    p[9:12] = map(str, B)
    p[13] = '1'              # smoothed Coulomb friction
    p[15] = str(k)           # tanh steepness

    folder = OUT / f'alpha_{alpha:g}' / case_name.replace(' ', '_').replace(',', '')
    folder.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(EXE)], input='\n'.join(p) + '\n', text=True,
                   cwd=folder, stdout=subprocess.DEVNULL, check=True)

    target = folder / f'k_{k:g}.csv'
    (folder / 'data.csv').replace(target)
    return np.loadtxt(target, delimiter=',', skiprows=1,
                      usecols=(0, 1, 2, 3), ndmin=2)

results = {}
rows = []

for alpha in ALPHAS:
    results[alpha] = {}
    for case_name, (mu, B, marker) in FIELD_CASES.items():
        ref = simulate(alpha, case_name, K0, mu, B)
        deviations = []

        for k in K_VALUES:
            data = ref if k == K0 else simulate(alpha, case_name, k, mu, B)
            displacement = np.linalg.norm(data[:, 1:4] - ref[:, 1:4], axis=1)
            maximum_mm = 1000.0 * displacement.max()
            deviations.append(maximum_mm)
            rows.append([alpha, case_name, k, maximum_mm])

        results[alpha][case_name] = np.array(deviations)

with (OUT / 'summary.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['alpha_deg', 'case', 'k_s_per_m', 'max_separation_from_k100_mm'])
    w.writerows(rows)

fig, ax = plt.subplots(figsize=(16, 6.45), layout='constrained')

for alpha in ALPHAS:
    for case_name, (_, _, marker) in FIELD_CASES.items():
        ax.semilogx(K_VALUES, results[alpha][case_name],
                    '-', color=ALPHA_COLORS[alpha], marker=marker,
                    lw=2.5, ms=8)

ax.set_xlabel(r'$k \, [\unit{\second\per\metre}]$')
ax.set_ylabel(r'$\max_t\|\Delta_k \mathbf{r}\|\, [\unit{\milli\metre}]$')
ax.set_xticks([10, 100, 1000])
ax.grid(True, which='both', alpha=0.3)

# Single, boxed two-row legend as in the previous convergence plot.
fields = list(FIELD_CASES.values())
upper_handles = [
    Line2D([], [], color='black', marker=fields[0][2], linestyle='None', ms=10),
    Line2D([], [], color='black', marker=fields[1][2], linestyle='None', ms=10),
    Line2D([], [], linestyle='None'),
]
upper_labels = [
    r'$B_z=2\,\mathrm{T},\ \mu_x=5\cdot10^{-5}\,\mathrm{A\,m^2}$',
    r'$B_x=5\cdot10^{-5}\,\mathrm{T},\ \mu_z=2\,\mathrm{A\,m^2}$',
    '',
]
lower_handles = [Line2D([], [], color=ALPHA_COLORS[a], lw=2.5) for a in ALPHAS]
lower_labels = [rf'$\varphi={a:g}^\circ$' for a in ALPHAS]

handles = [item for pair in zip(upper_handles, lower_handles) for item in pair]
labels = [item for pair in zip(upper_labels, lower_labels) for item in pair]
ax.legend(handles, labels, ncol=3, loc='upper right',
                   frameon=True, fancybox=True, framealpha=1,
                   edgecolor='black', borderpad=0.3, labelspacing=0.35,
                   columnspacing=1.2, handletextpad=0.5)

fig.savefig(OUT / 'k_sensitivity.pdf')
fig.savefig(OUT / 'k_sensitivity.png', dpi=300)
plt.close(fig)
