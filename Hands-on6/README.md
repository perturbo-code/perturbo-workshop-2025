# Hands-on - Spin Relaxation Times

## Introduction

This hands-on session will show how to compute spin relaxation times using PERTURBO. We will be running calculations on diamond with fully-relativistic pseudopotentials.

The materials used in this tutorial can be found at this box [link](https://caltech.box.com/s/ue5pcts8vmbhp5ko1arxdantxm2luagu) or in the Hands-on6 folder of the [github repository](https://github.com/perturbo-code/perturbo-workshop-2025/tree/main/Hands-on6).

## Setting the Docker environment

Similar to the previous Hands-on tutorials, we will be using the [docker environment](https://github.com/perturbo-code/perturbo-workshop-2025/tree/main/Hands-on1) for PERTURBO.

```bash
docker run -v "$(PWD):/home/user/run/workshop2025" -h perturbocker -it --rm --name perturbo perturbo/perturbo:gcc_openmp_3.0
cd Hands-on6
```

Then set OpenMP commands (if used)
```bash
export OMP_NUM_THREADS=4
```

## Option 1: Generate the epr.h5 from inputs
### SCF/NSCF/Wannier90 workflow

The epr.h5 file must be generated with fully-relativistic pseudopotentials, spin-orbit coupling, and noncollinear magnetism. There are additional spinor-related flags that need to be specified in the scf, nscf, and Wannier90 calculations.

The workflow for generating the epr.h5 file is the same as shown in earlier hands-on tutorials. The inputs include the required spinor-related variables and are provided at this box [link](https://caltech.box.com/s/5n47c18v0xetp2xltixvt5fvsdpon66n). 

Some comments about the workflow:
* Fully-relativistic pseudopotentials must be used
* In the SCF step, `noncolin` and `lspinorb` flags must be set to `.true.`
* In the NSCF step, `noncolin` and `lspinorb` flags must be set to `.true.`
* In the phonon step, there are no additional spinor-related variables specified.
* In the Wannier90 step, `spinors` flag must be `.true.` (in the **diam.win** input file) and `write_spn` flag must be `.true.` (in the **pw2wan.in** input file)

Note, the Wannier90 calculation flag `write_spn = .true.` will generate the **diam.spn** file. This file is required in the qe2pert step.

### Running qe2pert

Once the scf, nscf, Wannier90, and phonon calculations are complete, we run `qe2pert.x` to generate the **diam_epr.h5** file. The inputs can be found at this [link](https://caltech.box.com/s/clbiusws5uox5ltyjxzrd0bfacojk5pj).

Here is the input file (**qe2pert.in**):

```
&qe2pert
 prefix = 'diam'
 outdir = '.tmp'
 phdir = '../pw-ph-wann/phonon/save'
 nk1 = 4, nk2=4, nk3=4
 dft_band_min = 1
 dft_band_max = 17
 num_wann = 16
 lwannier = .true.
 lspinmat = .true.
/
```

Note, the flag `lspinmat` must be set to `.true.` to perform the spin calculations.

Before running `qe2pert.x`, ensure that the files **diam.spn**, **diam_centres.xyz**, **diam_u.mat**, and, when present, **diam_u_dis.mat** are in the run directory. The **diam.spn** file is generated from the Wannier90 calculation flag `write_spn`. Create a directory called `tmp` and inside it soft link the NSCF output directory **diam.save** and **diam.xml** (located in `pw-ph-wann/nscf/tmp`). 

```bash
# Go to the qe2pert directory
cd qe2pert

# Link Wannier matrices/centers if not already present
ln -sf ../pw-ph-wann/wann/diam.spn
ln -sf ../pw-ph-wann/wann/diam_centers.xyz
ln -sf ../pw-ph-wann/wann/diam_u.mat
ln -sf ../pw-ph-wann/wann/diam_u_dis.mat

# Link nscf files
mkdir tmp 
cd tmp 
ln -sf ../../pw-ph-wann/nscf/tmp/diam.xml
ln -sf ../../pw-ph-wann/nscf/tmp/diam.save/
```

Set OpenMP commands (if used) and run `qe2pert.x`:

```bash
export OMP_NUM_THREADS=4
qe2pert.x -i qe2pert.in | tee qe2pert.out
```

For faster performance with `qe2pert.x` using MPI parallelization, try:
```bash
mpirun -n 2 qe2pert.x -npools 2 -i qe2pert.in | tee qe2pert.out
```

Now, we can perform calculations using the **diam_epr.h5** file.

## Option 2: Download the existing epr file

Instead of running the above steps, one can download the **diam_epr.h5** file from this [link](https://caltech.box.com/s/20dn5dbcl8koqiypxrc0mqifgezkijza).

Download and move the **diam_epr.h5** file into the directory `/diamond/perturbo/`. We will create soft links to this file throughout the tutorial.

## E-ph spin-flip matrix elements

The inputs can be found in the `pert-ephmat-spin` directory:

```
cd perturbo-workshop-2025/Hands-on6/diamond/perturbo/pert-ephmat-spin
```

We will now compute the absolute values of the spin-flip e-ph matrix elements, summed over the number of electronic bands, given two lists of **k**- and **q**-points. In a typical scenario, one computes the e-ph matrix elements for a chosen **k**-point as a function of **q**-point. In this example, we are computing the e-ph spin-flip matrix elements summed over the bands from 9 to 10.

Here is the input file (**pert.in**):

```
&perturbo
 prefix = 'diam'
 calc_mode = 'ephmat_spin'
 fklist = 'diam_band.kpt'
 fqlist = 'diam_band.qpt'

 band_min = 9
 band_max = 10

 phfreq_cutoff = 3          !meV
/
```

Before running `perturbo.x`, ensure that these files (**diam_epr.h5**, **diam_band.kpt**, and **diam_band.qpt**) are in the run directory. Create a soft link to the **diam_epr.h5** file to the current directory:

```bash
ln -sf ../diam_epr.h5
```

Now, we can run `perturbo.x`:

```bash
perturbo.x -i pert.in | tee pert.out 
```

The calculation outputs 2 files - **diam.ephmat_flip** and **diam_ephmat_spin.yml**
- **diam.ephmat_flip** contains the spin-flip e-ph matrix elements in meV for the list of **k**- and **q**-points.
- **diam_ephmat_spin.yml** contains the inputs and outputs of the calculation in a YAML format for easier postprocessing.

We can use PerturboPy next to export the data from the YAML file to Python for postprocessing. We will use the plotting script **plot.py** to plot the phonon dispersion overlaid with a color map of the spin-flip e-ph matrix elements.

![diam-ephmat-spin.png](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on6/diamond/perturbo/pert-ephmat-spin/References/diam-ephmat-spin.png)

## Spin texture

The inputs can be found in the spins directory:

```
cd perturbo-workshop-2025/Hands-on6/diamond/perturbo/pert-spins
```

We will now perform a `calc_mode = 'spins'` calculation to compute the $\langle n \vert \sigma_z \vert n \rangle$ values for each **k**-point specified. 

Here is the input file(**pert.in**):

```
&perturbo
 prefix = 'diam'
 calc_mode = 'spins'
 fklist = 'diam_band.kpt'
/
```

Similar to the `'bands'` calculation, the **k**-point list (here **diam_band.kpt**) must be included in the run directory. Before running `perturbo.x`, remember to soft link the **diam_epr.h5** in the run directory.

```bash
ln -sf ../diam_epr.h5
```

Run `perturbo.x`:

```bash
perturbo.x -i pert.in | tee pert.out 
```

The calculation outputs 2 files - **diam.spins** and **diam_spins.yml**
- **diam.spins** contains the interpolated band structure and the corresponding $\langle n \vert \sigma_z \vert n \rangle$ values.
- **diam_spins.yml** contains the inputs and outputs of the calculation in a YAML format for easier postprocessing.

We can use PerturboPy next to export the data from the YAML file to Python for postprocessing. We will use the plotting script **plot.py** to plot the band structure overlaid with a color map of the $\langle n \vert \sigma_z \vert n \rangle$ values.

![diam-spins.png](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on6/diamond/perturbo/pert-spins/References/diam-spins.png)


Note, that the color bar can be plotted as a linear or logarithmic scale. Here are plotting the spin texture values on a linear scale.

## Spin Relaxation Times

### Step 1 - Setup

The inputs can be found in the setup directory:

```
cd perturbo-workshop-2025/Hands-on6/diamond/perturbo/pert-setup
```

To calculate the spin relaxation times, we first need to perform the `calc_mode = 'setup'` calculation.

Here is the input file (**pert.in**):

```
&perturbo
 prefix      = 'diam'
 calc_mode   = 'setup'

 boltz_kdim(1) = 80   !kgrid along three crystal axes
 boltz_kdim(2) = 80
 boltz_kdim(3) = 80

 boltz_emin = 17.3  !Energy window
 boltz_emax = 17.7  !Energy window
 band_min = 9      !Band index
 band_max = 10

 find_efermi = .true.
 ftemper  = 'diam.temper' !Name of .temper file
/
```

The **diam.temper** file contains information about temperature, chemical potential, and carrier concentration.

Before running `perturbo.x`, ensure that these files (**diam_epr.h5**, **diam_band.kpt**, and **diam.temper**) are in the run directory. We can create a soft link the **diam_epr.h5** file location to the current directory:

```bash
ln -sf ../diam_epr.h5
```

Now, we can run `perturbo.x`:

```bash
perturbo.x -i pert.in | tee pert.out 
```

The `'setup'` calculation takes the information about the grid, energy window, temperature and fermi level, and outputs 4 files - **diam_tet.h5**, **diam_tet.kpt**, **diam.doping**, **diam.dos** and **diam_setup.yml**
- **diam_tet.h5** and **diam_tet.kpt** contain information about the grid points after applying the specified energy window cutoffs and symmetry operations
- **diam.doping** contains information on the carrier concentrations from the Fermi level
- **diam.dos** contains the density of states as a function of carrier energies
- **diam_setup.yml** is a YAML file with the inputs and outputs from the calculation

We will be using the **diam_tet.kpt** and **diam.temper** files for the next step.

### Step 2 - Imsigma Spin

The inputs can be found in the imsigma-spin directory:

```
cd perturbo-workshop-2025/Hands-on6/diamond/perturbo/pert-imsigma-spin
```

Now, we will compute the imaginary part of the e-ph self-energy from the spin-flip process for states in a range of bands and with crystal momenta **k**.

Here is the input file (**pert.in**):

```
&perturbo
 prefix = 'diam'
 calc_mode = 'imsigma_spin'

 fklist = 'diam_tet.kpt'     !kpt file from setup calculation
 ftemper = 'diam.temper'

 band_min = 9
 band_max = 10

 phfreq_cutoff = 3           ! meV !Phonon frequency cutoff
 delta_smear = 20            ! meV      !Smearing value

 sampling = 'uniform'        !Type of q-point sampling
 nsamples = 10000          !Number of q-points
 /
```

We will use the **diam.temper** and **diam_tet.kpt** files from the previous `'setup'` calculation. Copy these files to the current run directory:

```bash
cp ../pert-setup/diam.temper .
cp ../pert-setup/diam_tet.kpt .
```

Make sure to create a soft link from the **diam_epr.h5** file location to the current directory:

```bash
ln -sf ../diam_epr.h5
```

Now, run `perturbo.x`:

```bash
export OMP_NUM_THREADS=4
perturbo.x -i pert.in | tee pert.out
```

The `'imsigma_spin'` calculation is the most expensive calculation in this tutorial. When executed serially, the calculation can take around 5-10 minutes. It's best to set OpenMP threads to the maximum number of cores in your machine. Note that this is not a completely converged calculation, as we have scaled down the parameters to reduce the computational cost.

If you use MPI parallelization, the calculation will be much faster.

```bash
mpirun -n 8 perturbo.x -npools 8 -i pert.in | tee pert.out
```

The calculation outputs 5 files - **diam.imsigma**, **diam.imsigma_mode**, **diam.imsigma_flip**, **diam.imsigma_flip_mode**, and **diam_imsigma_spin.yml**
- **diam.imsigma** contains $\mathrm{Im}\Sigma$ in meV as a function of carrier energies.
- **diam.imsigma_mode** contains mode resolved values of $\mathrm{Im}\Sigma$ in meV as a function of carrier energies.
- **diam.imsigma_flip** contains $\mathrm{Im}\Sigma_{flip}$ in meV as a function of carrier energies.
- **diam.imsigma_flip_mode** contains mode resolved values of $\mathrm{Im}\Sigma_{flip}$ in meV as a function of carrier energies.
- **diam_imsigma_spin.yml** contains the inputs and outputs from the calculation in a YAML format for easier postprocessing.

We can use PerturboPy next to export the data from the YAML file to Python for postprocessing. Use the plotting script **plot.py** to plot the values of $\mathrm{Im}\Sigma_{flip}$ vs energy.

![diam-imsigma-flip.png](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on6/diamond/perturbo/pert-imsigma-spin/References/diam-imsigma-flip.png)

We will need the **diam.imsigma_flip** file for the next and final step.

### Step 3 - Spin Lifetime

The inputs can be found in the spinlifetime directory

```
cd perturbo-workshop-2025/Hands-on6/diamond/perturbo/pert-spinlifetime
```

We are now ready to compute the spin relaxation times.

Here is the input file (**pert.in**):

```
&perturbo
 prefix = 'diam'
 calc_mode = 'spinlifetime'
  
 boltz_kdim(1) = 80
 boltz_kdim(2) = 80
 boltz_kdim(3) = 80

 ftemper  = 'diam.temper'

 band_min = 9
 band_max = 10
/
```

We will use the **diam.temper** and **diam_tet.h5** files from the `'setup'` calculation, and the **diam.imsigma_flip** file from the previous `'imsigma_spin'` calculation. Copy these files to the current directory:

```bash
cp ../pert-setup/diam.temper .
cp ../pert-setup/diam_tet.h5 .
cp ../pert-imsigma-spin/diam.imsigma_flip .
```

Make sure to create a soft link the **diam_epr.h5** file location to the current directory:

```bash
ln -sf ../diam_epr.h5
```

Now, run `perturbo.x`:

```bash
perturbo.x -i pert.in | tee pert.out
```

The calculation outputs 2 files - **diam.spin** and **diam_spinlifetime.yml**
- **diam.spin** contains the thermally averaged spin-flip scattering rates in nanoseconds for each temperature specified
- **diam_spinlifetime.yml** contains the inputs and outputs of the `'spinlifetime'` calculation in a YAML file
