### Written by Rasmus Raschke, Universität Hamburg, 2026 ###
## Contains a utility to quickly extract all data from any given csv file ##

import numpy as np
from pathlib import Path
from typing import Union
from types import SimpleNamespace
import re
from os import PathLike

import numpy as np


_FILENAME_RE = re.compile(r"data_mx_(?P<mx>.+)_my_(?P<my>.+)\.csv$")
_MZ_FILENAME_RE = re.compile(r"data_mz_(?P<mz>\d+)\.csv")
_AZIMUTH_FILENAME_RE = re.compile(r"data_azimuth_(?P<azimuth>[-+]?\d+(?:\.\d+)?)\.csv")
_AZIMUTH_FILENAME_RE2 = re.compile(r"data_azimuth_(?P<sign>[mp]?)(?P<int>\d+)_(?P<frac>\d+)_?\.csv")
PathLike = Union[str, Path]
_NEW_C_DIR_RE = re.compile(r"^c(?P<c>[-+m]?\d+(?:_\d+)?)$") 
_NEW_FILE_RE = re.compile(
    r"data_azimuth_(?P<azimuth>[-+m]?\d+(?:_\d+)?)_polar_(?P<polar>[-+m]?\d+(?:_\d+)?)\.csv$"
)

def load(path: PathLike) -> dict[str, np.ndarray]:
    """Extracts data with headers from a csv file

    Parameters
    ----------
    path : str
        Absolute or relative path to csv

    Returns
    -------
    cols : dict[str, np.ndarray]
        Contains csv columns with headers as {headers: colum data list} 
    """

    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    headers = data.dtype.names
    names = [h.strip() for h in headers]
    cols = {clean: data[orig] for clean, orig in zip(names, headers)}
    return cols

def extract(path: PathLike = ".") -> dict[str, SimpleNamespace]:
    """Build a dictionary with data from an arbitrary amount of csv files in a directory callable by the SimpleNamespace utility

    Parameters
    ----------
    path : str
        Either directory path or path to a single csv
    
    Returns
    -------
    datasets : dict[str, SimpleNamespace]
        Dictionary of extracted csv data accessible like datasets["file name"].column
    """

    folder = Path(path)
    datasets = {}

    if folder.is_dir():
        csv_files = folder.glob("*.csv")
    else:
        csv_files = folder.parent.glob(folder.name)

    for csv in csv_files:
        key = csv.stem
        cols = load(csv)
        datasets[key] = SimpleNamespace(**cols)
    return datasets

def decode(tag: str) -> float:
    """Decode filename tags to floats

    Parameters
    ----------
    tag : str 
        Tag to convert

    Returns
    -------
    float
        Tag as float
    """
    return float(tag.replace("m", "-").replace("p", "."))

def decode2(tag: str) -> float:
    """Decode filename tags to floats when p occures elsewhere in tag

    Parameters
    ----------
    tag : str 
        Tag to convert

    Returns
    -------
    float
        Tag as float
    """
    return float(tag.replace("m", "-").replace("_", "."))

def batch_extract(path: PathLike = ".") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract x values and mx, my tag from a big amount of data

    Parameters
    ----------
    path : PathLike
        Path to the data files
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        mx, my and last x value
    """
    folder = Path(path)
    mx_list: list[float] = []
    my_list: list[float] = []
    x_list: list[float] = []

    for csv in folder.glob("data_mx_*_my_*.csv"):
        match = _FILENAME_RE.fullmatch(csv.name)
        if match is None:
            continue
        
        mx = decode(match.group("mx"))
        my = decode(match.group("my"))
        cols = load(csv)
        x = np.asarray(cols["x"])[-1]
        mx_list.append(mx)
        my_list.append(my)
        x_list.append(float(x))
    return np.array(mx_list), np.array(my_list), np.array(x_list)

def batch_extract_mz(path: PathLike = ".") -> tuple[np.ndarray, np.ndarray]:
    """Extract y values and mz tag from a big amount of data

    Parameters
    ----------
    path : PathLike
        Path to the data files
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        mz and last y value
    """
    folder = Path(path)
    mz_list: list[float] = []
    y_list: list[float] = []

    for csv in folder.glob("data_mz_*.csv"):
        match = _MZ_FILENAME_RE.fullmatch(csv.name)
        if match is None:
            continue

        mz = int(match.group("mz")) / 100.0
        cols = load(csv)
        y = np.asarray(cols["y"])[-1]

        mz_list.append(mz)
        y_list.append(float(y))

    return np.array(mz_list), np.array(y_list)

def batch_extract_azimuth(path: PathLike = ".") -> tuple[np.ndarray, np.ndarray]:
    """Extract x values and azimuth tag from a big amount of data

    Parameters
    ----------
    path : PathLike
        Path to the data files
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        azimuth angle and last x value
    """
    folder = Path(path)
    azimuth_list: list[float] = []
    x_list: list[float] = []

    for csv in folder.glob("data_azimuth_*.csv"):
        match = _AZIMUTH_FILENAME_RE.fullmatch(csv.name)
        if match is None:
            continue

        azimuth = float(match.group("azimuth"))
        cols = load(csv)
        x = np.asarray(cols["x"])[-1]

        azimuth_list.append(azimuth)
        x_list.append(float(x))

    return np.array(azimuth_list), np.array(x_list)

def batch_extract_polar(path: PathLike = ".") -> tuple[np.ndarray, np.ndarray]:
    """Extract y values and polar tag from a big amount of data

    Parameters
    ----------
    path : PathLike
        Path to the data files
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        polar angle and last x value
    """
    folder = Path(path)
    azimuth_list: list[float] = []
    y_list: list[float] = []

    for csv in folder.glob("data_azimuth_*.csv"):
        match = _AZIMUTH_FILENAME_RE.fullmatch(csv.name)
        if match is None:
            continue

        azimuth = float(match.group("azimuth"))
        cols = load(csv)
        y = np.asarray(cols["y"])[-1]

        azimuth_list.append(azimuth)
        y_list.append(float(y))

    return np.array(azimuth_list), np.array(y_list)

def batch_extract_polar2(path: PathLike = ".") -> tuple[np.ndarray, np.ndarray]:
    """Extract y values and polar tag from a big amount of data

    Parameters
    ----------
    path : PathLike
        Path to the data files
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        polar angle and last x value
    """
    folder = Path(path)
    azimuth_list: list[float] = []
    y_list: list[float] = []

    for csv in folder.glob("data_azimuth_*.csv"):
        match = _AZIMUTH_FILENAME_RE2.fullmatch(csv.name)
        if match is None:
            continue

        sign = -1 if match.group("sign") == "m" else 1
        azimuth = sign * float(f"{match.group('int')}.{match.group('frac')}")
        cols = load(csv)
        y = np.asarray(cols["y"])[-1]

        azimuth_list.append(azimuth)
        y_list.append(float(y))

    return np.array(azimuth_list), np.array(y_list)


def extract_intervals(path: PathLike = "mag_moment_data", n_intervals: int = 10):
    """
    Extract t, x, y arrays for angles on a grid from 0 to pi.

    Parameters
    ----------
    path : PathLike
        Path to the data files

    Returns
    -------
    dict[str, dict[str, list[np.ndarray] | np.ndarray]]
        One entry per subfolder, e.g. "c0_0".
    """
    folder = Path(path)
    target_angles = np.linspace(-np.pi / 2, np.pi / 2, n_intervals + 1)

    result = {}

    for subfolder in sorted(p for p in folder.iterdir() if p.is_dir()):
        files: list[tuple[float, Path]] = []

        for csv in subfolder.glob("data_azimuth_*.csv"):
            match = _AZIMUTH_FILENAME_RE.fullmatch(csv.name)
            if match is None:
                continue
            angle = float(match.group("azimuth"))
            files.append((angle, csv))

        if not files:
            continue

        files.sort(key=lambda item: item[0])
        file_angles = np.array([a for a, _ in files], dtype=float)

        chosen_angles = []
        t_list = []
        x_list = []
        y_list = []
        for target in target_angles:
            idx = int(np.argmin(np.abs(file_angles - target)))
            angle, csv = files[idx]

            cols = load(csv)  # your existing CSV loader
            t = np.asarray(cols["t"])
            x = np.asarray(cols["x"])
            y = np.asarray(cols["y"])

            chosen_angles.append(angle)
            t_list.append(t)
            x_list.append(x)
            y_list.append(y)

        result[subfolder.name] = {
            "angle": np.array(chosen_angles),
            "t": t_list,
            "x": x_list,
            "y": y_list,
        }

    return result

def batch_extract2(
    path: PathLike = "results",
    out_path: PathLike | None = None,
    value_columns: tuple[str, ...] = ("x", "y"),
    reduce: str = "last",
):
    """
    Extract data from the new folder/file layout and optionally save it as an NPZ.

    Expected layout
    ---------------
    results/
        c0_0/
            data_azimuth_0_0_polar_1_0.csv
            data_azimuth_0_5_polar_1_0.csv
        c1_0/
            ...

    Parameters
    ----------
    path:
        Root directory containing the c-subfolders.
    out_path:
        If given, write an .npz file there.
    value_columns:
        CSV columns to extract.
    reduce:
        How to reduce each column if it is an array:
        - "last"  -> take the last entry
        - "all"   -> store the full array
        - "first" -> take the first entry

    Returns
    -------
    dict[str, np.ndarray]
        Arrays for c, azimuth, polar, and requested value columns.
    """
    path = Path(path)

    c_list: list[float] = []
    azimuth_list: list[float] = []
    polar_list: list[float] = []
    values: dict[str, list] = {col: [] for col in value_columns}

    def _reduce_array(arr: np.ndarray):
        if reduce == "last":
            return arr[-1]
        if reduce == "first":
            return arr[0]
        if reduce == "all":
            return arr
        raise ValueError("reduce must be 'last', 'first', or 'all'")

    for csv in path.rglob("data_azimuth_*_polar_*.csv"):
        if not csv.is_file():
            continue

        c_match = _NEW_C_DIR_RE.fullmatch(csv.parent.name)
        f_match = _NEW_FILE_RE.fullmatch(csv.name)
        if c_match is None or f_match is None:
            continue

        c = decode2(c_match.group("c"))
        azimuth = decode2(f_match.group("azimuth"))
        polar = decode2(f_match.group("polar"))

        data = np.genfromtxt(csv, delimiter=",", names=True, dtype=None, encoding=None)
        if data.dtype.names is None:
            continue

        c_list.append(c)
        azimuth_list.append(azimuth)
        polar_list.append(polar)

        for col in value_columns:
            if col not in data.dtype.names:
                raise KeyError(f"Column '{col}' not found in {csv}")
            values[col].append(_reduce_array(np.asarray(data[col])))

    result = {
        "c": np.asarray(c_list, dtype=float),
        "azimuth": np.asarray(azimuth_list, dtype=float),
        "polar": np.asarray(polar_list, dtype=float),
    }

    for col in value_columns:
        result[col] = np.asarray(values[col], dtype=object if reduce == "all" else float)

    if out_path is not None:
        np.savez(Path(out_path), **result)

    return result