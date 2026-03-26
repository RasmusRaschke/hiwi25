### Written by Rasmus Raschke, Universität Hamburg, 2026 ###
## Contains a utility to quickly extract all data from any given csv file ##

import numpy as np
from pathlib import Path
from typing import Union
from types import SimpleNamespace

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