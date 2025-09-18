# Hands-on Tutorial 2: Generating EPR File

## Introduction
In general, the typical Perturbo workflow can be visualized as follow:
<img width="1000" height="800" alt="image" src="https://github.com/user-attachments/assets/6f349399-d2fa-4648-a147-93d91679d2f2" />

After running QE and Wannier90, the next crucial step is generating the `epr.h5` file. This file contains the electron-phonon (e-ph) matrix elements on coarse Brillouin zone grids for both electrons (**k**-points) and phonons (**q**-points), represented in the Wannier function basis. The `qe2pert.x` interface program reads DFPT results from Quantum ESPRESSO and the unitary transformation matrices (*U* matrices) from Wannier90 to create this HDF5 format file for us.

To start, let us activate the dockerized Perturbo with:

```bash
docker run -v /path/to/your/Hands-on3:/home/user/run/Hands-on3 --user 500 -it --rm --name perturbo perturbo/perturbo:gcc_openmp_3.0

# Use 4 OMP threads:
export OMP_NUM_THREADS=4 
```
For this tutorial, I will be using the serial with OpenMP Docker version of Perturbo, but you are welcome to use other versions.

## Workflow for Polar Material: GaAs

### Step 1: Ground-state DFT Calculation

```bash
# First, locate the /pw-ph-wann/scf/ directory
cd Hands-on3/GaAs_polar/pw-ph-wann/scf/

# 'tee' command is to print outputs to screen for demo purpose
pw.x -in scf.in | tee scf.out
cp -r tmp ../ph/
```

This performs a self-consistent field (SCF) calculation to obtain the equilibrium charge density for DFPT.

### Step 2: Phonon Calculation with DFPT

```bash
# Now switch to "Hands-on3/GaAs_polar/pw-ph-wann/ph/"
cd ../ph/
ph.x -in ph.in | tee ph.out

# Collect the dfpt data to save/
chmod +x ph-collect.sh
./ph-collect.sh
```

In this step, we calculate the phonon frequencies and eigenvectors using Density Functional Perturbation Theory.

⚠️ **Note:** If you are using the parallel version, you may need to modify `ph-collect.sh` by changing all the file extension `*.dvscf` to `*.dvscf1`.

### Step 3: Non-self-consistent Calculation on Dense k-grid

```bash
# Next, we go to "Hands-on3/GaAs_polar/pw-ph-wann/nscf/" directory
cd ../nscf/
cp -r ../scf/tmp .
pw.x -in nscf.in | tee nscf.out
```

This non-self-consistent calculation computes electronic states on a denser k-point grid needed for Wannierization. The dense grid ensures accurate interpolation of electronic properties when constructing Wannier functions.

⚠️ **Note:**
- The full BZ grid must be provided to the input files. In this examples, we used the Wannier90's `kmesh.pl` utilities to generate the grid.
```
K_POINTS crystal
64
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  0.00000000  0.00000000  0.50000000  1.562500e-02
  0.00000000  0.00000000  0.75000000  1.562500e-02
  0.00000000  0.25000000  0.00000000  1.562500e-02
...
```
- In the following Wannierization step, you also need to use the exact same grid or Wannier90 will complain.

### Step 4: Wannierization
First, we link the previous `nscf/tmp/<prefix>.save` directory to the `wann`:
```bash
cd ../wann/

mkdir tmp
cd tmp && ln -s ../../nscf/tmp/gaas.save
cd ../
```
Now we are ready to run Wannier90:

```bash
wannier90.x -pp gaas.win
pw2wannier90.x -inp pw2wan.in | tee pw2wan.out
wannier90.x gaas.win
```

- **First line**: Pre-processes to determine the required overlaps between Bloch states
- **Second line**: Extracts wavefunctions from QE and computes overlap matrices M(k,b) and A(k,n,m)
- **Third line**: Performs the actual Wannierization, creating maximally localized Wannier functions and the transformation matrices needed for interpolation

### Step 5: Generate the EPR File
In `qe2pert` directory, we link or copy the following QE and Wannier90 outputs:
```bash
cd ../../qe2pert/

ln -s ../pw-ph-wann/wann/gaas_centres.xyz
ln -s ../pw-ph-wann/wann/gaas_u.mat
ln -s ../pw-ph-wann/wann/gaas_u_dis.mat

# Also link the nscf WFs:
mkdir tmp
cd tmp/ && ln -s ../../pw-ph-wann/nscf/tmp/gaas.save
cd ../
```

Run `qe2pert.x`:

```bash
qe2pert.x -in qe2pert.in | tee qe2pert.out
```
or with in parallel with `mpirun`:

```bash
mpirun -n 4 qe2pert.x -npools 4 -in qe2pert.in | tee qe2pert.out
```

This very crucial step reads:
- Electronic wavefunctions and energies from QE
- Phonon perturbation potentials from DFPT
- Wannier transformation matrices from Wannier90

And produces the `<prefix>_epr.h5` file containing e-ph matrix elements in the Wannier basis.

There is also a little `python3` code `access_epr.py` that print out some information contained in the `epr.h5` file for you to peek in to.

⚠️ **Important:** If you run with `mpirun`, the number of MPI tasks must equal the number of pools or the calculation will fail.

### Step 6: Run Perturbo Calculations

Before running any calculations with `perturbo.x`, make sure link or copy the EPR file to your working directory:

```bash
# Link the epr file to your working directory
ln -s /path/to/<prefix>_epr.h5 .

# Then run perturbo
perturbo.x -in pert.in | tee pert.out
```
or in parallel:

```bash
# npools must be equal to ntasks
mpirun -n 4 perturbo.x -npools 4 -in pert.in | tee pert.out
```

With the `epr.h5` file, you can perform various Perturbo calculations (band interpolation, transport, dynamics, etc.) by specifying different `calc_mode` options in the input file.

As a demo, let us interpolate the el-ph matrix on a dense grid:
```bash
cd GaAs_polar/perturbo/pert-ephmat/
perturbo.x -in pert.in | tee pert.out
python3 plot_ephmat.py
```
Which would result in the following plot:
<img width="800" height="600" alt="image" src="https://github.com/user-attachments/assets/c32af3f5-8bad-4990-93bb-91121d3c3534" />

I also included other sample calculation modes in the same directory for you to explore later.

## Special Note for Polar Materials

For polar semiconductors like GaAs, Perturbo automatically handles the separation of e-ph interactions into:

- **Long-range (Fröhlich) part:** Arising from macroscopic electric fields due to polar optical phonons
- **Short-range part:** Local atomic displacements

By default, both contributions are computed and included in calculations like `imsigma` and `imsigma_spin`, separating these 2 parts could help with reducing the computation cost.

## Part 2: Electron-Phonon with Quadrupole Corrections
While standard electron-phonon calculations capture the monopole and dipole contributions to the interaction, certain materials require higher-order multipole corrections for accurate results. Quadrupole corrections become particularly important in materials with strong covalent bonding or when studying properties that are sensitive to the detailed shape of the electron-phonon coupling. These corrections account for the spatial variation of the electric field gradient produced by phonons, which can significantly affect carrier scattering rates in some semiconductors.

## Workflow for Quadrupole Corrections
In this section, we will demonstrate how to add the quadrupole correction to the e-ph coupling matrices from first-principles.
This is particularly necessary to piezoelectric (PE) e-ph interaction and long-range scattering mechanism due to acoustic phonons in noncentrosymmetric polar materials.

The quadrupole correction is computed as:
<img width="1152" height="106" alt="image" src="https://github.com/user-attachments/assets/b1019be9-4241-46c4-8221-7fc522a65686" />

where $Q_{\kappa,\alpha \beta \gamma}$ is the dynamical quadrupole tensor.

For more information, see:
* [Detailed theories by Miquel Royo and Massimiliano Stengel](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.041027)

and also our publications:
* [Theory of Quadrupole Correction to electron-phonon interactions (Si and PbTiO3)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.136602)
* [Applications to GaN](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.125203)


### Step 1: Perform Linear Response Calculations with ABINIT

In addition to `QE` and `Wannier90`, you need to calculate the quadrupole tensors using ABINIT's linear response capabilities. This involves computing third-order energy derivatives with respect to two electronic perturbations and one phonon perturbation.

For detailed instructions on setting up these calculations, please visit the [ABINIT Long-Wave Quadrupole Tutorial](https://docs.abinit.org/tutorial/lw_quad/).

The ABINIT calculation will produce output files containing the quadrupole tensor elements. In this example, we will use the file named `tlw_3.abo` as the output for demo. This file contains the matrix elements needed to correct the electron-phonon coupling for quadrupole effects.

⚠️ **Note:** You open the file to read, but make sure not to edit it or error may occur.

### Step 2: Extract Quadrupole Tensors from ABINIT Output
Pick either method A or B below that apply to you.

### Version A: If you only have `hdf5` (which is true if you are running Perturbocker):

```bash
cd Si_with_quadrupole/qtensor/qtensor-utils/

# Compile the "qtensor_hdf5.f90" version using the h5fc wrapper:
/opt/hdf5/bin/h5fc -o qtensor qtensor_hdf5.f90

# Run qtensor in the same directory as your ABINIT output
cd ../

./qtensor-utils/qtensor

# Enter the file name "tlw_3.abo"
```

### Version B: If you know you have `h5fortran` library installed:
Navigate to the `qtensor` folder and modify the `Makefile` to match your system configuration.
```bash
cd Si_with_quadrupole/qtensor/qtensor-utils/
```
Open `Makefile`:
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

./qtensor-utils/qtensor

# Enter the file name (tlw_3.abo)
```

Note that he `qtensor` utility must be executed in the same directory as your ABINIT output file (`tlw_3.abo`). It will parse the ABINIT output and extract the quadrupole tensor elements.

After this, you should obtain the file `qtensor.h5`. Now rename this file and copy to where you would run `qe2pert.x`.

```bash
# Rename the file to match your calculation prefix
mv qtensor.h5 si_qtensor.h5

# Copy to the directory where you'll run qe2pert.x
cp si_qtensor.h5 ../qe2pert/
```

### Step 4: Generating EPR file
Now, perform `QE` and `Wannier90`, then run `qe2pert.x` normally as above. The code will automatically detect `<prefix>_qtensor.h5` file and apply quadrupole correction to the eph calculations when generating the `epr.h5` file.

```bash
qe2pert.x -in qe2pert.in | tee qe2pert.out
```

You can verify that quadrupole corrections has been applied by checking the `qe2pert.out` file for messages indicating that the quadrupole tensor file was successfully read and processed.

**Done!** :tada::tada::tada:

Now, you are all set to perform subsequent `PERTURBO` calculations.

The quadrupole corrections are now embedded in the el-ph matrix elements (the `epr.h5` file) and will automatically be included in all calculations of scattering rates, transport properties, and carrier dynamics.

Below is the result comparing eph with and without QC from [our GaN paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.125203).
We could see it is particularly important for the acoustic modes near the **q** = Γ point.

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/5b47764b-9636-45ea-975e-3f6b4ee2b3e3" />

