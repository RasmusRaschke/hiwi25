import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from os import PathLike
import numpy as np
import multiprocessing as mp
mp.set_start_method("spawn", force=True)

_FILENAME_RE = re.compile(
    r"data_c_(?P<c>[-+]?\d+(?:\.\d+)?)_azimuth_(?P<azimuth>[-+]?\d+(?:\.\d+)?)\.csv"
)


def _replace_leading_number(line, value):
    return re.sub(
        r"^\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        f"{value:.16g}",
        line,
        count=1,
    )

"""
def make_input_file(template, theta, c, out_path):
    lines = template.read_text().splitlines(True)
    mx = c * np.sin(theta)
    my = c * np.cos(theta)
    mz = 0.0
    lines[6] = _replace_leading_number(lines[6], mx)
    lines[7] = _replace_leading_number(lines[7], my)
    lines[8] = _replace_leading_number(lines[8], mz)
    out_path.write_text("".join(lines))
"""
#for changing outer magnetic field
def make_input_file(template, theta, c, out_path):
    lines = template.read_text().splitlines(True)
    mx = 0.5 * np.sin(np.pi / 200)
    my = 0.0
    mz = 0.5 * np.cos(theta)
    Bx = 0.0
    By = c * np.sin(theta)
    Bz = c * np.cos(theta)
    lines[6] = _replace_leading_number(lines[6], mx)
    lines[7] = _replace_leading_number(lines[7], my)
    lines[8] = _replace_leading_number(lines[8], mz)
    lines[9] = _replace_leading_number(lines[9], Bx)
    lines[10] = _replace_leading_number(lines[10], By)
    lines[11] = _replace_leading_number(lines[11], Bz)
    out_path.write_text("".join(lines))


def run_one_case(exe, template, theta, c, output_dir, input_name="input.base", angle_digits=6, c_digits=6):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    c_tag = f"{c:.{c_digits}f}"
    theta_tag = f"{theta:.{angle_digits}f}"
    final_csv = output_dir / f"data_c_{c_tag}_azimuth_{theta_tag}.csv"
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        input_path = tmpdir / input_name
        make_input_file(template=template, theta=theta, c=c, out_path=input_path)
        with input_path.open("rb") as f:
            subprocess.run([str(exe)], cwd=tmpdir, stdin=f, check=True)
        csvs = list(tmpdir.glob("*.csv"))
        shutil.move(str(csvs[0]), str(final_csv))

    return final_csv

def _run_one_case_star(args):
    return run_one_case(*args)

def sweep_theta_and_c(
    exe="./magsphere.out",
    template="./input.base",
    c_values = (1.0,),
    n_steps=100,
    output_dir="./results",
    angle_digits=6,
    c_digits=6,
    max_cores=12
):
    exe = Path(exe).resolve()
    template = Path(template).resolve()
    output_dir = Path(output_dir).resolve()
    thetas = np.linspace(-np.pi, np.pi, n_steps + 1)
    #thetas = np.array([np.pi / 2, np.pi / 20, np.pi / 200, np.pi/2000, 0.0])
    #thetas = np.array([0.0])
    print(thetas)
    tasks = [
        (exe, template, float(theta), float(c), output_dir, "input.base", angle_digits)
        for c in c_values
        for theta in thetas
    ]
    outputs: list[Path] = []
    with ProcessPoolExecutor(max_workers=max_cores) as pool:
        for csv_path in pool.map(_run_one_case_star, tasks):
            outputs.append(csv_path)

    return outputs

if __name__ == "__main__":
    sweep_theta_and_c(
        exe = "./magsphere.out",
        template = "./input.base",
        c_values = [5.0e-5],
        n_steps = 400,
        output_dir = "./results",
        angle_digits = 6,
        c_digits = 2,
        max_cores=11
    )