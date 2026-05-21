import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import multiprocessing as mp

import numpy as np


_FILENAME_RE = re.compile(
    r"data_azimuth_(?P<psi>[-+]?\d+(?:\.\d+)?)_polar_(?P<phi>[-+]?\d+(?:\.\d+)?)\.csv"
)


def _replace_leading_number(line, value):
    return re.sub(
        r"^\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        f"{value:.16g}",
        line,
        count=1,
    )


def _num_tag(x, digits=6):
    """
    Convert numbers into filetags.
    Examples:
        0.0   -> 0_0
        1.25  -> 1_25
        -0.5  -> m0_5
    """
    s = f"{float(x):.{digits}f}".rstrip("0").rstrip(".")
    if "." not in s:
        s += ".0"
    if s == "-0.0":
        s = "0.0"
    return s.replace("-", "m").replace(".", "_")


def _as_values(x):
    if np.isscalar(x):
        return [float(x)]
    return [float(v) for v in x]


def make_input_file(template, c, theta, phi, out_path):
    lines = template.read_text().splitlines(True)
    mx = c * np.sin(theta) * np.sin(phi)
    my = c * np.sin(theta) * np.cos(phi)
    mz = c * np.cos(theta)
    Bx = 0.0
    By = 0.0
    Bz = 5e-5
    lines[6] = _replace_leading_number(lines[6], mx)
    lines[7] = _replace_leading_number(lines[7], my)
    lines[8] = _replace_leading_number(lines[8], mz)
    lines[9] = _replace_leading_number(lines[9], Bx)
    lines[10] = _replace_leading_number(lines[10], By)
    lines[11] = _replace_leading_number(lines[11], Bz)
    out_path.write_text("".join(lines))


def run_one_case(
    exe,
    template,
    c,
    theta,
    phi,
    output_dir,
    input_name="input.base",
    angle_digits=6,
    c_digits=6,
):
    exe = Path(exe)
    template = Path(template)
    output_dir = Path(output_dir)

    c_tag = _num_tag(c, digits=c_digits)
    theta_tag = _num_tag(theta, digits=angle_digits)
    phi_tag = _num_tag(phi, digits=angle_digits)

    c_dir = output_dir / f"c{c_tag}"
    c_dir.mkdir(parents=True, exist_ok=True)

    final_csv = c_dir / f"data_azimuth_{theta_tag}_polar_{phi_tag}.csv"

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        input_path = tmpdir / input_name

        make_input_file(template=template, c=c, theta=theta, phi=phi, out_path=input_path)

        with input_path.open("rb") as f:
            subprocess.run([str(exe)], cwd=tmpdir, stdin=f, check=True)

        csvs = list(tmpdir.glob("*.csv"))
        if len(csvs) != 1:
            raise RuntimeError(f"Expected exactly one CSV output, found {len(csvs)} in {tmpdir}")

        shutil.move(str(csvs[0]), str(final_csv))

    return final_csv


def _run_one_case_star(args):
    return run_one_case(*args)


def sweep_parameters(
    exe="./magsphere.out",
    template="./input.base",
    c_values=1.0,
    theta_values=0.0,
    phi_values=0.0,
    output_dir="./results",
    angle_digits=6,
    c_digits=6,
    max_cores=12,
):
    exe = Path(exe).resolve()
    template = Path(template).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    c_values = _as_values(c_values)
    theta_values = _as_values(theta_values)
    phi_values = _as_values(phi_values)
    tasks = [
        (exe, template, c, theta, phi, output_dir, "input.base", angle_digits, c_digits)
        for c in c_values
        for theta in theta_values
        for phi in phi_values
    ]
    outputs = []
    with ProcessPoolExecutor(max_workers=max_cores) as pool:
        for csv_path in pool.map(_run_one_case_star, tasks):
            outputs.append(csv_path)
    return outputs


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    sweep_parameters(
        exe="./magsphere.out",
        template="./input.base",
        c_values=[0.0],
        theta_values=[0.0],
        phi_values=[0.0],
        output_dir="./results",
        angle_digits=6,
        c_digits=2,
        max_cores=11,
    )