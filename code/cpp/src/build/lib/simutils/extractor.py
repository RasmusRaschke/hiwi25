### Written by Rasmus Raschke, Universität Hamburg, 2026 ###
## Contains a utility to quickly extract all data from any given csv file ##

import numpy as np
from pathlib import Path
from typing import Union
from types import SimpleNamespace
import re
from functools import lru_cache

_FILENAME_RE = re.compile(r"data_mx_(?P<mx>.+)_my_(?P<my>.+)\.csv$")
PathLike = Union[str, Path]

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

@lru_cache(maxsize=None)
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