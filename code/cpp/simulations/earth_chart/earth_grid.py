from pathlib import Path
import numpy as np
import pandas as pd
from wmm import wmm_calc
import runner


###########################################################################
# PARAMETERS
YEAR = 2026
MONTH = 6
DAY = 29
ALTITUDE_M = 0.0
DLAT = 5.0          # latitude spacing [deg]
DLON = 5.0          # longitude spacing [deg]
C_VALUES = [1.0]
OUTPUT_DIR = Path("results")
OUTFILE = OUTPUT_DIR / "earth_results.npz"
MAX_CORES = 12
###########################################################################


def create_grid(dlat, dlon):
    lats = np.arange(-90.0, 90.0 + dlat, dlat)
    lons = np.arange(-180.0, 180.0 + dlon, dlon)
    Lon, Lat = np.meshgrid(lons, lats)
    return Lat, Lon


def create_wmm():
    model = wmm_calc()
    model.setup_time(YEAR, MONTH, DAY)
    return model


def field_at(model, lat, lon):
    model.setup_env(
        lat=lat,
        lon=lon,
        alt=ALTITUDE_M,
        unit="m",
        msl=True,
    )
    try:
        B = model.get_all()
        bx = float(np.asarray(B["x"]).squeeze())
        by = float(np.asarray(B["y"]).squeeze())
        bz = -float(np.asarray(B["z"]).squeeze())
    except Exception:
        bx = float(np.asarray(model.get_Bx()).squeeze())
        by = float(np.asarray(model.get_By()).squeeze())
        bz = -float(np.asarray(model.get_Bz()).squeeze())
    # convert nT -> Tesla
    bx *= 1e-9
    by *= 1e-9
    bz *= 1e-9

    return bx, by, bz


def build_cases(model, Lat, Lon):
    cases = []
    nlat, nlon = Lat.shape
    for i in range(nlat):
        for j in range(nlon):
            lat = float(Lat[i, j])
            lon = float(Lon[i, j])
            bx, by, bz = field_at(model, lat, lon)
            cases.append(
                (
                    lat,
                    lon,
                    bx,
                    by,
                    bz,
                )
            )

    return cases


def read_last_values(csv_file):
    df = pd.read_csv(csv_file)
    x_last = float(df["x"].iloc[-1])
    t_last = float(df["t"].iloc[-1])
    return x_last, t_last


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("Creating latitude-longitude grid...")
    Lat, Lon = create_grid(DLAT, DLON)
    print("Initialising WMM...")
    model = create_wmm()
    print("Computing magnetic field...")
    cases = build_cases(model, Lat, Lon)
    print(f"{len(cases)} grid points")
    print("Running simulations...")
    csv_files = runner.sweep_parameters(
        exe="./magsphere.out",
        template="./input.base",
        c_values=C_VALUES,
        cases=cases,
        output_dir=OUTPUT_DIR,
        max_cores=MAX_CORES,
    )
    nlat, nlon = Lat.shape

    x_last = np.empty((nlat, nlon))
    t_last = np.empty((nlat, nlon))

    Bx = np.empty((nlat, nlon))
    By = np.empty((nlat, nlon))
    Bz = np.empty((nlat, nlon))

    k = 0

    for i in range(nlat):
        for j in range(nlon):
            bx, by, bz = cases[k][2:]
            Bx[i, j] = bx
            By[i, j] = by
            Bz[i, j] = bz
            x, t = read_last_values(csv_files[k])
            x_last[i, j] = x
            t_last[i, j] = t
            k += 1

    np.savez(
        OUTFILE,
        lat=Lat,
        lon=Lon,
        Bx=Bx,
        By=By,
        Bz=Bz,
        x=x_last,
        t=t_last,
    )
    print("Saved")
    print(OUTFILE)


if __name__ == "__main__":
    main()