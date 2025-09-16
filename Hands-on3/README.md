# Hands-on Tutorial 2: Generating EPR File

## Introduction

The first and crucial step in the Perturbo workflow is generating the `epr.h5` file. This file contains the electron-phonon (e-ph) matrix elements on coarse Brillouin zone grids for both electrons (**k**-points) and phonons (**q**-points), represented in the Wannier function basis. The `qe2pert.x` interface program reads DFPT results from Quantum ESPRESSO and the unitary transformation matrices (*U* matrices) from Wannier90 to create this HDF5 format file.

To start, let us activate the dockerized Perturbo with:

```bash
docker run -v /path/to/your/Hands-on3:/home/user/run/Hands-on3 --user 500 -it --rm --name perturbo perturbo/perturbo:gcc_openmp_3.0

# Use 4 OMP threads:
export OMP_NUM_THREADS=4 
```
Here I will be using the serial with OpenMP Docker version of Perturbo, but you are welcome to use the others.

## Workflow for Polar Material: GaAs

### Step 1: Ground-state DFT Calculation

```bash
# First, locate the /pw-ph-wann/scf/ directory
cd /Hands-on3/GaAs_polar/pw-ph-wann/scf/

pw.x -in scf.in > scf.out
```

This performs a self-consistent field (SCF) calculation to obtain the equilibrium charge density for DFPT.

### Step 2: Phonon Calculation with DFPT

```bash
# Now we go to ../ph/
cd /Hands-on3/GaAs_polar/pw-ph-wann/ph/

ph.x -in ph.in > ph.out

# Collect the dfpt data to save/
sh ph-collect.sh
```

In this step, we calculate the phonon frequencies and eigenvectors using Density Functional Perturbation Theory.

⚠️ **Note:** If you are using the parallel version, you may need to modify `ph-collect.sh` by changing all the file extension `*.dvscf` to `*.dvscf1`.

### Step 3: Non-self-consistent Calculation on Dense k-grid

```bash
# Next, we go to the ../nscf/ directory to compute the WFs for Wannier90
cd /Hands-on3/GaAs_polar/pw-ph-wann/nscf/

pw.x -in nscf.in > nscf.out
```

This non-self-consistent calculation computes electronic states on a denser k-point grid needed for Wannierization. The dense grid ensures accurate interpolation of electronic properties when constructing Wannier functions.

⚠️ **Note:**
- The full BZ grid must be provided to the input files. In this examples, we used the Wannier90's `kmesh.pl` utilities to generate the grid.
- In the following Wannierization step, you also need to use the exact same grid.

### Step 4: Wannierization Process

```bash
wannier90-3.1.0/wannier90.x -pp gaas.win

pw2wannier90.x -inp pw2wan.in > pw2wan.out

wannier90-3.1.0/wannier90.x gaas.win
```

- **First line**: Pre-processes to determine the required overlaps between Bloch states
- **Second line**: Extracts wavefunctions from QE and computes overlap matrices M(k,b) and A(k,n,m)
- **Third line**: Performs the actual Wannierization, creating maximally localized Wannier functions and the transformation matrices needed for interpolation

### Step 5: Generate the EPR File
Run `qe2pert.x`:

```bash
cd ./qe2pert
qe2pert.x -in qe2pert.in > qe2pert.out
```
or with in parallel with `mpirun`:

```bash
mpirun -n 4 qe2pert.x -npools 4 -in qe2pert.in > qe2pert.out
```

This very crucial step reads:
- Electronic wavefunctions and energies from QE
- Phonon perturbation potentials from DFPT
- Wannier transformation matrices from Wannier90

And produces the `<prefix>_epr.h5` file containing e-ph matrix elements in the Wannier basis. 

⚠️ **Important:** If you run with `mpirun`, the number of MPI tasks must equal the number of pools or the calculation will fail.

### Step 6: Run Perturbo Calculations

Before running any calculations with `perturbo.x`, make sure link or copy the EPR file to your working directory:

```bash
# Link the epr file to your working directory
ln -sf /path/to/<prefix>_epr.h5 .

# Then run perturbo
perturbo.x -in pert.in > pert.out
```
or in parallel:

```bash
# npools must be equal to ntasks
mpirun -n 4 perturbo.x -npools 4 -in pert.in > pert.out
```

With the `epr.h5` file ready, you can now perform various Perturbo calculations (band interpolation, transport, dynamics, etc.) by specifying different `calc_mode` options in the input file.

## Special Note for Polar Materials

For polar semiconductors like GaAs, Perturbo automatically handles the separation of e-ph interactions into:

- **Long-range (Fröhlich) part:** Arising from macroscopic electric fields due to polar optical phonons
- **Short-range part:** Local atomic displacements

By default, both contributions are computed and included in calculations like `imsigma` and `imsigma_spin`, separating these 2 parts could help with reducing the computation cost.

## Part 2: Electron-Phonon with Quadrupole Corrections
While standard electron-phonon calculations capture the monopole and dipole contributions to the interaction, certain materials require higher-order multipole corrections for accurate results. Quadrupole corrections become particularly important in materials with strong covalent bonding or when studying properties that are sensitive to the detailed shape of the electron-phonon coupling. These corrections account for the spatial variation of the electric field gradient produced by phonons, which can significantly affect carrier scattering rates in some semiconductors.

## Workflow for Quadrupole Corrections

### Step 1: Perform Linear Response Calculations with ABINIT

In addition to `QE` and `Wannier90`, you need to calculate the quadrupole tensors using ABINIT's linear response capabilities. This involves computing third-order energy derivatives with respect to two electronic perturbations and one phonon perturbation.

For detailed instructions on setting up these calculations, please visit the [ABINIT Long-Wave Quadrupole Tutorial](https://docs.abinit.org/tutorial/lw_quad/).

The ABINIT calculation will produce output files containing the quadrupole tensor elements. In this example, we will use the file named `tlw_3.abo` as the output for demo. This file contains the matrix elements needed to correct the electron-phonon coupling for quadrupole effects.

⚠️ **Note:** You open the file to read, but make sure not to edit it or error may occur.

### Step 2: Extract Quadrupole Tensors from ABINIT Output
Pick either method A or B below that apply to you.

### Version A: If you only have `hdf5` (which is true if you are running Perturbocker):

```bash
cd ./qtensor/qtensor-utils/

# Compile the "qtensor_hdf5.f90" version using the h5fc wrapper:
/opt/hdf5/bin/h5fc -o qtensor qtensor_hdf5.f90

# Run qtensor in the same directory as your ABINIT output
cd ../

./qtensor-utils/qtensor

# Enter the file name "tlw_3.abo"
```

### Version B: If you know you have `h5fortran` library installed:
Navigate to the `qtensor` folder and modify the `Makefile` to match your system configuration.

```makefile
# Change to gfortran/mpif90 depending on your HDF5
FC = mpifort

# Set the path to your HDF5 installation
# This should point to the root directory containing 'include' and 'lib' subdirectories
H5_ROOT = /path/to/root/hdf5
```

Once you've configured the Makefile, compile and run the `qtensor` extraction tool:

```bash
# Compile the qtensor utility (if not already done)
make

# Run qtensor in the same directory as your ABINIT output
cd ../

./qtensor

# Enter the file name (tlw_3.abo)
```

Note that he `qtensor` utility must be executed in the same directory as your ABINIT output file (`tlw_3.abo`). It will parse the ABINIT output and extract the quadrupole tensor elements.

After this, you should obtain the file `qtensor.h5`. Now rename this file and copy to where you would run `qe2pert.x`.

```bash
# Rename the file to match your calculation prefix
mv qtensor.h5 <prefix>_qtensor.h5

# Copy to the directory where you'll run qe2pert.x
cp <prefix>_qtensor.h5 /path/to/qe2pert/directory/
```

### Step 4: Generating EPR file
Now, run `qe2pert.x` normally as above. The code will automatically detect this file and apply quadrupole correction to the eph calculations when generating the `epr.h5` file.

```bash
qe2pert.x -in qe2pert.in > qe2pert.out
```

You can verify that quadrupole corrections were applied by checking the `qe2pert.out` file for messages indicating that the quadrupole tensor file was successfully read and processed.

Done. You are now ready to perform subsequent `PERTURBO` calculations.

The quadrupole corrections are now embedded in the electron-phonon matrix elements and will automatically be included in all calculations of scattering rates, transport properties, and carrier dynamics.

## Conclusion

The inclusion of quadrupole corrections can lead to:
- More accurate carrier mobilities, especially at low temperatures where the detailed structure of the electron-phonon coupling matters most
- Better agreement with experimental data for materials with strong covalent character
- Improved description of hot carrier relaxation in ultrafast dynamics simulations

These corrections are particularly important for materials like silicon, germanium, and III-V semiconductors where the quadrupole contribution can be comparable to the standard electron-phonon coupling in certain scattering channels.
