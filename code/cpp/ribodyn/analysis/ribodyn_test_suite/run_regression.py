#!/usr/bin/env python3
"""Run and validate the ribodyn regression suite.

Usage:
    python3 run_regression.py /absolute/or/relative/path/to/solver
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


def load_csv(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    if data.dtype.names is None:
        raise RuntimeError(f"No CSV header found in {path}")
    return {name.strip(): np.asarray(data[name]) for name in data.dtype.names}


def maximum_relative_drift(values: np.ndarray) -> float:
    scale = max(abs(float(values[0])), 1.0)
    return float(np.max(np.abs(values - values[0])) / scale)


def vector(data, names):
    return np.column_stack([data[name] for name in names])


def close(value, expected, atol, message):
    error = float(np.max(np.abs(np.asarray(value) - np.asarray(expected))))
    if error > atol:
        raise AssertionError(f"{message}: error {error:.3e} > {atol:.3e}")


def common_checks(data):
    required = (
        "t", "x", "y", "z", "vx", "vy", "vz",
        "qw", "qx", "qy", "qz", "Ox", "Oy", "Oz",
        "T_trans", "T_rot", "E_total",
        "constraint_residual", "quaternion_norm",
    )
    missing = [name for name in required if name not in data]
    if missing:
        raise AssertionError(f"Missing output columns: {missing}")

    for name, values in data.items():
        if not np.all(np.isfinite(values)):
            raise AssertionError(f"Column {name} contains non-finite values")

    q_error = float(np.max(np.abs(data["quaternion_norm"] - 1.0)))
    if q_error > 1e-10:
        raise AssertionError(f"Quaternion norm error {q_error:.3e}")

    residual = float(np.max(np.abs(data["constraint_residual"])))
    if residual > 1e-9:
        raise AssertionError(f"Constraint residual {residual:.3e}")


def compare(a, b, names, atol, label):
    for name in names:
        close(a[name], b[name], atol, f"{label}: column {name}")


def validate(results):
    # 01: free particle
    d = results["01_free_particle"]
    tf = d["t"][-1]
    close([d["x"][-1], d["y"][-1], d["z"][-1]],
          np.array([1.0, -2.0, 0.5]) * tf, 2e-10,
          "free particle final position")
    compare(d, {**d, "vx": np.full_like(d["vx"], 1.0),
                      "vy": np.full_like(d["vy"], -2.0),
                      "vz": np.full_like(d["vz"], 0.5)},
            ("vx", "vy", "vz"), 2e-12, "free particle velocity")

    # 02: torque-free asymmetric top
    d = results["02_free_asymmetric_top"]
    if maximum_relative_drift(d["E_total"]) > 2e-7:
        raise AssertionError("asymmetric-top energy drift too large")
    Rmom = np.column_stack((d["Ox"], 2.0*d["Oy"], 3.0*d["Oz"]))
    body_L2 = np.sum(Rmom * Rmom, axis=1)
    if maximum_relative_drift(body_L2) > 5e-7:
        raise AssertionError("asymmetric-top angular-momentum magnitude drift too large")

    # 03/04: uniform gravity
    for key in ("03_uniform_gravity_dalembert", "04_uniform_gravity_lagrangian"):
        d = results[key]
        tf = d["t"][-1]
        close(d["z"][-1], 10.0 - 0.5*9.81*tf*tf, 2e-8, f"{key} final z")
        close(d["vz"][-1], -9.81*tf, 2e-8, f"{key} final vz")
        if maximum_relative_drift(d["E_total"]) > 2e-10:
            raise AssertionError(f"{key} total-energy drift too large")
    compare(results["03_uniform_gravity_dalembert"],
            results["04_uniform_gravity_lagrangian"],
            ("x","y","z","vx","vy","vz","E_total"), 2e-10,
            "uniform gravity mode agreement")

    # 05/06: rolling
    for key in ("05_rolling_gravity_dalembert", "06_rolling_gravity_lagrangian"):
        d = results[key]
        close(d["z"], 0.1, 2e-10, f"{key} height")
        close(d["vx"], 1.0, 2e-9, f"{key} vx")
        close(d["Oy"], 10.0, 2e-8, f"{key} Oy")
        if np.max(np.abs(d["constraint_residual"])) > 1e-10:
            raise AssertionError(f"{key} rolling residual too large")
    compare(results["05_rolling_gravity_dalembert"],
            results["06_rolling_gravity_lagrangian"],
            ("x","y","z","vx","vy","vz","Ox","Oy","Oz"), 5e-9,
            "rolling mode agreement")

    # 07/08: circular central orbit
    for key in ("07_central_orbit_dalembert", "08_central_orbit_lagrangian"):
        d = results[key]
        radius = np.sqrt(d["x"]**2 + d["y"]**2 + d["z"]**2)
        if np.max(np.abs(radius - 1.0)) > 3e-8:
            raise AssertionError(f"{key} orbital-radius error too large")
        if maximum_relative_drift(d["E_total"]) > 5e-10:
            raise AssertionError(f"{key} energy drift too large")
        close([d["x"][-1], d["y"][-1]], [1.0, 0.0], 2e-6,
              f"{key} one-orbit return")
    compare(results["07_central_orbit_dalembert"],
            results["08_central_orbit_lagrangian"],
            ("x","y","z","vx","vy","vz"), 2e-9,
            "central-gravity mode agreement")

    # 09/10: uniform electric field
    for key in ("09_uniform_electric_dalembert", "10_uniform_electric_lagrangian"):
        d = results[key]
        tf = d["t"][-1]
        close(d["x"][-1], tf*tf, 2e-9, f"{key} final x")
        close(d["vx"][-1], 2.0*tf, 2e-9, f"{key} final vx")
        if maximum_relative_drift(d["E_total"]) > 2e-10:
            raise AssertionError(f"{key} energy drift too large")
    compare(results["09_uniform_electric_dalembert"],
            results["10_uniform_electric_lagrangian"],
            ("x","y","z","vx","vy","vz","E_total"), 2e-9,
            "electric-field mode agreement")

    # 11/12: cyclotron motion, omega_c = qB/m = 1
    for key in ("11_uniform_magnetic_dalembert", "12_uniform_magnetic_lagrangian"):
        d = results[key]
        speed = np.sqrt(d["vx"]**2 + d["vy"]**2 + d["vz"]**2)
        if np.max(np.abs(speed - 1.0)) > 3e-9:
            raise AssertionError(f"{key} speed is not conserved")
        close([d["x"][-1], d["y"][-1]], [0.0, 0.0], 3e-6,
              f"{key} one-cyclotron-period return")
        if maximum_relative_drift(d["E_total"]) > 2e-10:
            raise AssertionError(f"{key} energy drift too large")
    compare(results["11_uniform_magnetic_dalembert"],
            results["12_uniform_magnetic_lagrangian"],
            ("x","y","z","vx","vy","vz"), 3e-8,
            "magnetic-field mode agreement")

    # 13: permanent dipole in uniform B
    d = results["13_magnetic_dipole"]
    if maximum_relative_drift(d["E_total"]) > 3e-6:
        raise AssertionError("magnetic-dipole total-energy drift too large")
    mu_body = vector(d, ("mu_body_x", "mu_body_y", "mu_body_z"))
    close(mu_body, np.array([1.0, 0.0, 0.0]), 2e-12,
          "body magnetic moment")
    mu_world = vector(d, ("mu_world_x", "mu_world_y", "mu_world_z"))
    mu_norm = np.linalg.norm(mu_world, axis=1)
    close(mu_norm, 1.0, 2e-10, "world magnetic-moment norm")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("solver", type=Path, help="Path to the solver executable")
    parser.add_argument("--inputs", type=Path,
                        default=Path(__file__).resolve().parent / "inputs")
    parser.add_argument("--results", type=Path,
                        default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()

    solver = args.solver.expanduser().resolve()
    if not solver.is_file():
        raise SystemExit(f"Solver executable not found: {solver}")

    inputs = args.inputs.resolve()
    results_dir = args.results.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    failures = []

    for input_file in sorted(inputs.glob("*.in")):
        key = input_file.stem
        case_dir = results_dir / key
        if case_dir.exists():
            shutil.rmtree(case_dir)
        case_dir.mkdir(parents=True)

        print(f"[RUN ] {key}")
        process = subprocess.run(
            [str(solver), str(input_file.resolve())],
            cwd=case_dir,
            text=True,
            capture_output=True,
        )
        if process.returncode != 0:
            failures.append((key, process.stderr or process.stdout))
            print(f"[FAIL] {key}: executable returned {process.returncode}")
            continue

        output = case_dir / "output.csv"
        if not output.is_file():
            failures.append((key, "output.csv was not produced"))
            print(f"[FAIL] {key}: no output.csv")
            continue

        try:
            data = load_csv(output)
            common_checks(data)
            results[key] = data
            print(f"[ OK ] {key}")
        except Exception as exc:
            failures.append((key, str(exc)))
            print(f"[FAIL] {key}: {exc}")

    if not failures:
        try:
            validate(results)
        except Exception as exc:
            failures.append(("cross-case validation", str(exc)))

    if failures:
        print("\nRegression suite failed:")
        for key, error in failures:
            print(f"  - {key}: {error.strip()}")
        raise SystemExit(1)

    print(f"\nAll {len(results)} cases passed.")
    print(f"Results: {results_dir}")


if __name__ == "__main__":
    main()
