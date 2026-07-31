# ribodyn regression suite

This suite exercises the currently exposed executable features:

1. Free translation.
2. Torque-free asymmetric rigid-body rotation.
3. Uniform gravity in d'Alembert and Lagrangian modes.
4. Rolling constraint under gravity in both modes.
5. Central gravity and a circular orbit in both modes.
6. Uniform electric field in both modes.
7. Uniform magnetic field / cyclotron motion in both modes.
8. Permanent magnetic dipole torque in a uniform magnetic field.
9. Quaternion normalization, energy diagnostics, and constraint residuals.

## Run

From the suite directory:

```bash
python3 run_regression.py /path/to/ribodyn/build/solver
```

For the project path used in the examples:

```bash
python3 run_regression.py   /home/lyserg/Documents/hiwi25/code/cpp/ribodyn/build/solver
```

Each case runs in its own directory under `results/`, because the executable
always writes `output.csv` to its current working directory.

## Plot a result

Copy `overview_full.py` and `extractor.py` into this directory, then run:

```bash
python3 overview_full.py results/13_magnetic_dipole/output.csv
```

## Important coverage gap

The current executable recognizes only `emType uniform`. Therefore this suite
does not exercise:

- finite-difference magnetic-field Jacobians,
- translational force on a dipole in a nonuniform magnetic field,
- explicitly time-dependent vector potentials / induced electric fields.

Add a nonuniform concrete EM field and a corresponding `main.cpp` factory entry
before treating those branches as tested.
