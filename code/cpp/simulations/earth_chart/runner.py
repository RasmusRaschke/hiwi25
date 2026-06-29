import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import multiprocessing as mp

import numpy as np


_FILENAME_RE = re.compile(
    r"data_lat_(?P<lat>[-+]?\d+(?:\.\d+)?)_lon_(?P<lon>[-+]?\d+(?:\.\d+)?)\.csv"
)


def _replace_leading_number(line, value):
    return re.sub(
        r"^\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        f"{value:.16g}",
        line,
        count=1,
    )


def _num_tag(x, digits=6):
    s = f"{float(x):.{digits}f}".rstrip("0").rstrip(".")
    if "." not in s:
        s += ".0"
    if s == "-0.0":
        s = "0.0"
    return s.replace("-", "m").replace(".", "_")


def _as_values(x):
    if np.isscalar(x):
        return [x]
    return list(x)


def make_input_file(
    template,
    c,
    Bx,
    By,
    Bz,
    out_path,
):
    """
    Parameters
    ----------
    Bx,By,Bz : float
        Magnetic field components in Tesla.
    c : float
        Magnitude of magnetic moment.

    The magnetic moment is chosen parallel to B.
    """

    lines = template.read_text().splitlines(True)

    B = np.array([Bx, By, Bz], dtype=float)

    Bnorm = np.linalg.norm(B)
    if Bnorm == 0:
        raise ValueError("Magnetic field magnitude is zero.")

    mu = c * B / Bnorm

    mx, my, mz = mu

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
    lat,
    lon,
    Bx,
    By,
    Bz,
    output_dir,
    input_name="input.base",
    coord_digits=4,
    c_digits=4,
):
    exe = Path(exe)
    template = Path(template)
    output_dir = Path(output_dir)

    c_tag = _num_tag(c, digits=c_digits)
    lat_tag = _num_tag(lat, digits=coord_digits)
    lon_tag = _num_tag(lon, digits=coord_digits)

    c_dir = output_dir / f"c{c_tag}"
    c_dir.mkdir(parents=True, exist_ok=True)

    final_csv = c_dir / f"data_lat_{lat_tag}_lon_{lon_tag}.csv"

    with tempfile.TemporaryDirectory() as tmp:

        tmpdir = Path(tmp)

        input_path = tmpdir / input_name

        make_input_file(
            template=template,
            c=c,
            Bx=Bx,
            By=By,
            Bz=Bz,
            out_path=input_path,
        )

        with input_path.open("rb") as f:
            subprocess.run(
                [str(exe)],
                cwd=tmpdir,
                stdin=f,
                check=True,
            )

        csvs = list(tmpdir.glob("*.csv"))

        if len(csvs) != 1:
            raise RuntimeError(
                f"Expected exactly one CSV output, found {len(csvs)}."
            )

        shutil.move(str(csvs[0]), str(final_csv))

    return final_csv


def _run_one_case_star(args):
    return run_one_case(*args)


def sweep_parameters(
    exe="./magsphere.out",
    template="./input.base",
    c_values=1.0,
    cases=None,
    output_dir="./results",
    coord_digits=4,
    c_digits=4,
    max_cores=12,
):
    """
    Parameters
    ----------
    cases : iterable

        Each element must be

        (lat, lon, Bx, By, Bz)

        where B is already given in Tesla.
    """

    exe = Path(exe).resolve()
    template = Path(template).resolve()
    output_dir = Path(output_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    c_values = _as_values(c_values)

    if cases is None:
        raise ValueError("cases must be supplied.")

    tasks = []

    for c in c_values:
        for lat, lon, Bx, By, Bz in cases:

            tasks.append(
                (
                    exe,
                    template,
                    c,
                    lat,
                    lon,
                    Bx,
                    By,
                    Bz,
                    output_dir,
                    "input.base",
                    coord_digits,
                    c_digits,
                )
            )

    outputs = []

    with ProcessPoolExecutor(max_workers=max_cores) as pool:

        for csv in pool.map(_run_one_case_star, tasks):
            outputs.append(csv)

    return outputs


if __name__ == "__main__":

    mp.set_start_method("spawn", force=True)

    # Simple test field (50 µT northward)

    cases = [
        (
            0.0,
            0.0,
            50e-6,
            0.0,
            0.0,
        )
    ]

    sweep_parameters(
        exe="./magsphere.out",
        template="./input.base",
        c_values=[1.0],
        cases=cases,
        output_dir="./results",
        max_cores=1,
    )