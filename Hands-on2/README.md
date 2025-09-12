# Hands-on Tutorial 2: Generating EPR File

## Introduction

The first and crucial step in the Perturbo workflow is generating the `epr.h5` file. This file contains the electron-phonon (e-ph) matrix elements on coarse Brillouin zone grids for both electrons (**k**-points) and phonons (**q**-points), represented in the Wannier function basis. The `qe2pert.x` interface program reads DFPT results from Quantum ESPRESSO and the unitary transformation matrices (*U* matrices) from Wannier90 to create this HDF5 format file.

## Workflow for Polar Material: GaAs

### Step 1: Ground-state DFT Calculation

```bash
pw.x -in scf.in > scf.out
```

This performs a self-consistent field (SCF) calculation to obtain the ground-state electronic structure. The calculation determines the equilibrium charge density and Kohn-Sham wavefunctions, which form the foundation for all subsequent calculations.

### Step 2: Phonon Calculation with DFPT

```bash
ph.x -in ph.in > ph.out
```

This calculates phonon frequencies and eigenvectors using Density Functional Perturbation Theory. For GaAs (a polar material), this step computes both the short-range interatomic force constants and identifies the long-range polar contributions from the Born effective charges and dielectric tensor.

### Step 3: Non-self-consistent Calculation on Dense k-grid

```bash
pw.x -in nscf.in > nscf.out
```

This non-self-consistent calculation computes electronic states on a denser k-point grid needed for Wannierization. The dense grid ensures accurate interpolation of electronic properties when constructing Wannier functions.

### Step 4: Wannierization Process

```bash
wannier90-3.1.0/wannier90.x -pp gaas.win
pw2wannier90.x -inp pw2wan.in > pw2wan.out
wannier90-3.1.0/wannier90.x gaas.win
```

- **First command**: Pre-processes to determine the required overlaps between Bloch states
- **Second command**: Extracts wavefunctions from QE and computes overlap matrices M(k,b) and A(k,n,m)
- **Third command**: Performs the actual Wannierization, creating maximally localized Wannier functions and the transformation matrices needed for interpolation

### Step 5: Generate the EPR File

```bash
mpirun -n 4 qe2pert.x -npools 4 -in qe2pert.in > qe2pert.out
```

This crucial step reads:
- Electronic wavefunctions and energies from QE
- Phonon perturbation potentials from DFPT
- Wannier transformation matrices from Wannier90

And produces the `<prefix>_epr.h5` file containing e-ph matrix elements in the Wannier basis. 

⚠️ **Important:** The number of MPI tasks must equal the number of pools or the calculation will fail.

### Step 6: Run Perturbo Calculations

Before running any calculations with `perturbo.x`, link the EPR file to your working directory:

```bash
# Link the epr file to your working directory
ln -sf /path/to/<prefix>_epr.h5 .

# Then run perturbo
mpirun -n 4 perturbo.x -npools 4 -in pert.in > pert.out
```

With the `epr.h5` file ready, you can now perform various Perturbo calculations (band interpolation, transport, dynamics, etc.) by specifying different `calc_mode` options in the input file.

## Special Note for Polar Materials

For polar semiconductors like GaAs, Perturbo automatically handles the separation of e-ph interactions into:

- **Long-range (Fröhlich) part:** Arising from macroscopic electric fields due to polar optical phonons
- **Short-range part:** Local atomic displacements

By default, both contributions are computed and included in calculations like `imsigma` and `imsigma_spin`, ensuring accurate treatment of carrier scattering in polar materials.
