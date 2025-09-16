# Hands-on Session 5 Ultrafast Carrier Dynamics

The following tutorial demonstrates the feature to simulate ultrafast carrier dynamics with PERTURBO by simulating a pump-probe spectroscopy experiment in GaAs.

## Description of the example problem ##

We compute the transient absorption signal during a pump-probe type experiment in GaAs. 

For a typical pump-probe experiment, shortly a few tens of fs after the excitation, electrons can be described by a hot Fermi-Dirac distribution with a temperature of a few thousand degrees. But in this example, we will explicitly model how electron and hole dynamics contribute to the transient absorption signal.

From $t=0$, we simulate the pump pulse by computing the additional excited electron and hole populations at each time snapshot during the time of the pump. Then we time-evolve the dynamics of the electron and hole populations separately, where the excited populations are added during the pulse, then cooled in a phonon bath of 300 K. Phonon populations in this simulation are set to be time-independent. The carrier concentration is determined by the user-defined pump magnitude.

## Objectives ##
- Prepare input files and run PERTURBO for `calc_mode = 'setup'`, `calc_mode = 'dynamics-run'`, and `calc_mode = 'dynamics-pp'`.
- Analyze and visualize simulation results using PerturboPy.
- Explore advanced features such as using a initial population generated from a pump pulse.

## Required files ##

The materials used in this tutorial can be found [here](https://github.com/perturbo-code/perturbo-workshop-2025/tree/main/Hands-on5), Hands-on5 of this workshop. The folder contains scripts `generate_pump_pulse.py`, `compute_trans_abs.py`, and `plot_trans_abs.py`, and `plot_bands.py`, `gaas_band.kpt` and `pert.in` in `bands` subdirectory. To follow through the example with reference results, one can download the box folder [here](https://caltech.box.com/s/p42xq9w7p33z56k71kdfkv14xy6r8lo3).

`gaas_epr.h5` file used in this tutorial can be found [here](https://caltech.app.box.com/s/va0568719rn0y28j1s2mfk09ihw3buuj/folder/298661392249), with example input files from [Ultrafast Spectroscopy](https://perturbo-code.github.io/mydoc_spectroscopy.html).

Download and move `gaas_epr.h5` into the directory `Hands-on5`.


## Setting the environment 

Set up [docker environment](https://github.com/perturbo-code/perturbo-workshop-2025/tree/main/Hands-on1) for PERTURBO.

```bash
docker run -v "$(PWD):/home/user/run/Hands-on5" -h perturbocker -it --rm --name perturbo perturbo/perturbo:gcc_openmp_3.0
cd Hands-on5
```

Then set OpenMP commands (if used)
```bash
export OMP_NUM_THREADS=4
```

## Workflow


### Setting up momentum grids with setup

We first need to run [bands](https://perturbo-code.github.io/mydoc_interpolation.html#electronic-bandscalc_mode--bands) calculation to determine the appropriate energy windows for electrons and holes determined by the electronic structure of GaAs.

```bash
cd bands
```

First, link `gaas_epr.h5` file.

```bash
ln -sf ../gaas_epr.h5
```

The file `gaas_band.kpt` contains a high-symmetry path to visualize the band structure:
```
6
  0.500 0.500 0.500  50 !L
  0.000 0.000 0.000  50 !G
  0.500 0.000 0.500  20 !X
  0.500 0.250 0.750  20 !W
  0.375 0.375 0.750  50 !K
  0.000 0.000 0.000  1  !G
```

Run `bands` calculation with the following input file `pert.in` 
```
&perturbo
 prefix      = 'gaas'
 calc_mode   = 'bands'
 fklist   = 'gaas_band.kpt'
/
```

```bash
perturbo.x -i pert.in | tee pert.out
```

For faster performance with `perturbo.x` using MPI parallelization, try

```bash
mpirun -n 4 perturbo.x -npools 4 -i pert.in | tee pert.out
```

The number of `-npools` has to be the same as number of MPI tasks.

We use `PerturboPy` package to visualize the electronic and hole bands stored in `gaas_bands.yml` by running the following Python script `plot_bands.py`. 

```python
import numpy as np
import perturbopy.postproc as ppy
import matplotlib.pyplot as plt
plt.rcParams.update(ppy.plot_tools.plotparams)

gaas_bands = ppy.Bands.from_yaml('gaas_bands.yml')

fig, ax  = plt.subplots(figsize=(10, 8))
labels_dict = {'L': [0.5, 0.5, 0.5], 'G': [0.0, 0.0, 0.0], 'X': [0.5, 0.0, 0.5], 'W': [0.5, 0.25, 0.75], 'K': [0.375, 0.375, 0.75], 'G': [0.0, 0.0, 0.0]}
gaas_bands.kpt.add_labels(labels_dict)
gaas_bands.plot_bands(ax, show_kpoint_labels=True)

bandgap_size, point = gaas_bands.direct_bandgap(4,5)
CBM = min(gaas_bands.bands[5])
VBM = max(gaas_bands.bands[4])
print(f'CBM = {CBM} eV')
print(f'VBM = {VBM} eV')
print(f'direct band gap is {bandgap_size} eV')

emin_holes, emax_holes = 4.9, 5.3 # eV
emin_electrons, emax_electrons = 5.9, 6.7 # eV

xmin, xmax = ax.get_xlim()
ax.fill_between([xmin, xmax], [emin_holes, emin_holes],
                [emax_holes, emax_holes], color='blue', alpha=0.2)
ax.fill_between([xmin, xmax], [emin_electrons, emin_electrons],
                [emax_electrons, emax_electrons], color='red', alpha=0.2)
ax.set_xlim(xmin, xmax)

plt.savefig('bands.png', dpi=200)
```

For dynamics calculation, we need to select an energy range to include all carriers within the range of the excitation. After running 
```bash
python3 plot_bands.py
```

![bands](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on5/images/bands.png)

We select the bottom ~0.7 eV of the conduction bands for electrons, and the top ~0.4 eV of valence bands for holes based on the phonon bath temperature (300 K) and the pump pulse energy (1.5 eV) in this example.

Now we proceed to `calc_mode = 'setup'` to create k-grids for both electrons and holes.

```bash
cd ../
mkdir cdyna-elec
cd cdyna-elec
ln -sf ../gaas_epr.h5
```

First we create the `gaas.temper` file for a given carrier concentration and phonon bath temperature.

```bash
1         
300.0     0.0     1.0E+13
```

Create a new input file `setup.in` from the [interative workflow](https://perturbo-code.github.io/mydoc_interactive_workflow.html) for `setup`, copy the information into the input file `setup.in`, and edit the following input parameters with our desired `band_min`, `band_max`, `boltz_emin` and `boltz_emax`. 

```bash=
&perturbo
! ***Mandatory parameters***
 calc_mode           = 'setup'
 prefix              = 'gaas'
 boltz_kdim(1)       = 60
 boltz_kdim(2)       = 60
 boltz_kdim(3)       = 60
 ftemper             = 'gaas.temper'
 

! ***Optional parameters***
 boltz_qdim(1)      = 20
 boltz_qdim(2)      = 20
 boltz_qdim(3)      = 20
 hole               = .false.
 find_efermi        = .true.
 band_min           = 5
 band_max           = 5
 boltz_emin         = 5.9         ! eV
 boltz_emax         = 6.7         ! eV
/
```

Set `find_efermi = .true.` to determine the chemical potential. Although, in `dynamics-run`, the only used variables in `prefix.temper` are the line number 1 and the phonon bath temperature.

If parallelization is not available for running this example, use a smaller grid size, such as `boltz_kdim(1) = 24` and `boltz_qdim(1) = 12`. The workflow will run, although the results will not be physical.

Run `setup` calculation
```bash
perturbo.x -i setup.in | tee setup.out
```

Now we generate `gaas_tet.kpt` and `gaas_tet.h5`, which store the information of k-grid. `gaas_tet.h5` will be read for the following dynamics calculations.

Similarly, we set up the hole grids with the following commands

```bash
cd ../
mkdir cdyna-hole
cd cdyna-hole
ln -sf ../gaas_epr.h5
```

We similarly create the following `gaas.temper`
```
1         
300.0     0.0     1.0E+13
```

In `setup.in`, we modify the band and energy windows as well as `hole = .true.`: 

```
&perturbo
! ***Mandatory parameters***
 calc_mode           = 'setup'
 prefix              = 'gaas'
 boltz_kdim(1)       = 60
 boltz_kdim(2)       = 60
 boltz_kdim(3)       = 60
 ftemper             = 'gaas.temper'

! ***Optional parameters***
 boltz_qdim(1)      = 20
 boltz_qdim(2)      = 20
 boltz_qdim(3)      = 20
 hole               = .true.
 find_efermi        = .true.
 band_min           = 2
 band_max           = 4
 boltz_emin         = 4.9         ! eV
 boltz_emax         = 5.3         ! eV
/
```

Note: 
  - We set the phonon q-grid `boltz_qdim` to be less than `boltz_kdim` for a faster calculation. One should confirm that the resulting dynamics converge with respect to these grid size settings. 
  - The grids for holes need to be the same as those for electrons to use an initial population generated by a pump pulse using PerturboPy.

Run `setup` calculation
```bash
perturbo.x -i setup.in | tee setup.out
```

### Set initial carrier populations

#### Initial populations by a pump pulse (this example)

In this example, we specify initial populations with a pump pulse. We first run a mock dynamics with a single time step, to compute the electron-phonon matrix elements and create a sample `prefix_cdyna.h5` and `prefix_dyna_init.yml`, which will be used for PerturboPy to generate the initial carrier distributions due to the pulse. 

In the `cdyna-elec` directory, we create a new input file `pulse.in` from the [interative workflow](https://perturbo-code.github.io/mydoc_interactive_workflow.html) for `dynamics-run`, and edit the input parameters.

```
&perturbo
! ***Mandatory parameters***
 calc_mode           = 'dynamics-run'
 prefix              = 'gaas'
 boltz_kdim(1)       = 60
 boltz_kdim(2)       = 60
 boltz_kdim(3)       = 60
 ftemper             = 'gaas.temper'
 boltz_nstep         = 1
 output_nstep        = 1
 time_step           = 1.0             ! fs
 boltz_init_dist     = 'gaussian'
 scat_impl           = 'std'

! ***Optional parameters***
 solver             = 'rk4'
 boltz_init_e0      = -9999.0         ! eV
 boltz_init_smear   = 20.0            ! meV
 boltz_init_ampl    = 0.0             ! arbitrary
 delta_smear        = 8.0             ! meV
 phfreq_cutoff      = 1.0             ! meV
 boltz_qdim(1)      = 20
 boltz_qdim(2)      = 20
 boltz_qdim(3)      = 20
 hole               = .false.
 band_min           = 5
 band_max           = 5
 boltz_emin         = 5.9             ! eV
 boltz_emax         = 6.7             ! eV
 load_scatter_eph   = .false.
 tmp_dir            = './tmp'
 pump_pulse         = .false.
 pump_pulse_fname   = 'pump_pulse.h5'
 yaml_fname = 'gaas_dyna_init.yml'
/
```

Note:
  - The energy and band window should be the same as that in `calc_mode = 'setup'` calculation done earlier.
  - We set the output YAML file name as `yaml_fname = 'gaas_dyna_init.yml'` to distinguish from the default name for `calc_mode = 'dynamics-run'`.
  - We set `boltz_init_ampl = 0` to disable the default distribution. 
  - For this mock run, we also set `pump_pulse = .false.`
  - Since we are only using `gaas_cdyna.h5` for pump pulse input generation, we set `boltz_nstep = 1`.

We run the mock dynamics, which is a relatively expensive calculation. Parallelization is strongly recommended. Alternatively, one can also include the `tmp` folder under `Refs/cdyna-elec` and `Refs/cdyna-hole` from [here](https://caltech.box.com/s/p42xq9w7p33z56k71kdfkv14xy6r8lo3), use 1 MPI task, and set `load_scatter_eph = .true.` to skip part of the expensive calculation.

```bash
perturbo.x -i pulse.in | tee pulse.out
```

We can look at result `gaas_cdyna.h5`
```bash
h5ls gaas_cdyna.h5
```
And the yaml output file

```bash
vim gaas_dyna_init.yml
```

Similarly, we perform the initial dynamics for holes in the `cdyna-hole` directory with `pulse.in` slightly changed from that for electrons:

```
&perturbo
! ***Mandatory parameters***
 calc_mode           = 'dynamics-run'
 prefix              = 'gaas'
 boltz_kdim(1)       = 60
 boltz_kdim(2)       = 60
 boltz_kdim(3)       = 60
 ftemper             = 'gaas.temper'
 boltz_nstep         = 1
 output_nstep        = 1
 time_step           = 1.0             ! fs
 boltz_init_dist     = 'gaussian'
 scat_impl           = 'std'

! ***Optional parameters***
 solver             = 'rk4'
 boltz_init_e0      = -9999.0         ! eV
 boltz_init_smear   = 20.0            ! meV
 boltz_init_ampl    = 0.0             ! arbitrary
 delta_smear        = 8.0             ! meV
 phfreq_cutoff      = 1.0             ! meV
 boltz_qdim(1)      = 20
 boltz_qdim(2)      = 20
 boltz_qdim(3)      = 20
 hole               = .true.
 band_min           = 2
 band_max           = 4
 boltz_emin         = 4.9         ! eV
 boltz_emax         = 5.3         ! eV
 load_scatter_eph   = .false.
 tmp_dir            = './tmp'
 pump_pulse         = .false.
 pump_pulse_fname   = 'pump_pulse.h5'
 yaml_fname = 'gaas_dyna_init.yml'
/
```

```bash
perturbo.x -i pulse.in | tee pulse.out
```

### Generate pump pulse

Now we are ready to generate the pump pulse `hdf5` file with [PerturboPy](https://perturbopy.readthedocs.io/en/latest/index.html). First, we make sure that we have `gaas_dyna_init.yml` and `gaas_cdyna.h5` under both `cdyna-elec` and `cdyna-hole` directories.

Then we create the following script `generate_pump_pulse.py`

```python
import perturbopy.postproc as ppy

prefix = 'gaas'
elec_folder = 'cdyna-elec' 
hole_folder = 'cdyna-hole' 

# Define electron paths and DynaRun object
cdyna_elec_path = f'{elec_folder}/{prefix}_cdyna.h5'
elec_tet_path = f'{elec_folder}/{prefix}_tet.h5'
elec_yaml_path = f'{elec_folder}/{prefix}_dyna_init.yml'
elec_dyna_run = ppy.DynaRun.from_hdf5_yaml(cdyna_elec_path, elec_tet_path, elec_yaml_path)

# Define hole paths and DynaRun object
cdyna_hole_path = f'{hole_folder}/{prefix}_cdyna.h5'
hole_tet_path = f'{hole_folder}/{prefix}_tet.h5'
hole_yaml_path = f'{hole_folder}/{prefix}_dyna_init.yml'
hole_dyna_run = ppy.DynaRun.from_hdf5_yaml(cdyna_hole_path, hole_tet_path, hole_yaml_path)

# Print DynaRun info
print(elec_dyna_run)
print(hole_dyna_run)

# Generate the pump pulse HDF5 files in the
# elec_pump_pulse_path and hole_pump_pulse_path folders
pump_energy = 1.5 		# energy of pump in eV
elec_pump_pulse_path = f'{elec_folder}/pump_pulse_elec_Epump_{pump_energy:05.2f}.h5'
hole_pump_pulse_path = f'{hole_folder}/pump_pulse_hole_Epump_{pump_energy:05.2f}.h5'
pump_time_step = 1.0 		# pump pulse time step in fs
pump_factor = 1.0 		# factor multiplied to generate desired carrier concentration
pump_spectral_width_fwhm = 0.09 # full width half maximum of spectral width in eV
pump_duration_fwhm = 20.0 	# full width half maximum of pump in fs
pump_time_window = 50.0 	# duration of pulse in fs
finite_width = True 		# specify a gaussian pulse; if false use step function
tr_dipoles = None 		# transient dipole squared in .npy file
ppy.utils.spectra_generate_pulse.setup_pump_pulse(
    elec_pump_pulse_path,
    hole_pump_pulse_path,
    elec_dyna_run,
    hole_dyna_run,
    pump_energy=pump_energy,
    pump_time_step=pump_time_step,
    pump_duration_fwhm=pump_duration_fwhm,
    pump_spectral_width_fwhm=pump_spectral_width_fwhm,
    pump_time_window=pump_time_window,
    pump_factor=pump_factor,
    finite_width=finite_width,
    animate=True, 
    plot_scale=1e3,
    cnum_check=True,
    tr_dipoles_sqr=tr_dipoles)

elec_dyna_run.close_hdf5_files()
hole_dyna_run.close_hdf5_files()
```

and run

```bash
python3 generate_pump_pulse.py
```

Note:
  - One can incorporate transient dipole amplitudes by computing them using software such as [YAMBO](https://www.yambo-code.eu). PertuboPy can read the values stored in an `.npy` file with dimensions (number of electron bands, number of hole bands, number of k points). With `tr_dipoles_sqr = None`, a constant value `1.0` will be used.
  - `cnum_check = True` enables carrier concentration check (see the section below).

The script above will generate `cdyna-elec/pump_pulse_elec_Epump_1.50.h5` and `cdyna-hole/pump_pulse_hole_Epump_1.50.h5`.

Run 
```bash
h5ls cdyna-elec/pump_pulse_elec_Epump_1.50.h5
```
to see the list of data stored in the pump pulse file.

```
energy_profile           Dataset {200, 2}
hole                     Dataset {SCALAR}
num_bands                Dataset {SCALAR}
num_kpoints              Dataset {SCALAR}
num_steps                Dataset {SCALAR}
optional_params          Dataset {10}
pump_duration_fwhm       Dataset {SCALAR}
pump_energy              Dataset {SCALAR}
pump_factor              Dataset {SCALAR}
pump_pulse_snaps         Group
pump_spectral_width_fwhm Dataset {SCALAR}
pump_time_step           Dataset {SCALAR}
time_profile             Dataset {101, 2}
time_window              Dataset {SCALAR}
```

### Computing excited carrier concentration

In the previous section, we generated the pump pulse with `cnum_check = True`, which creates a `cnum_check` subdirectory in `cdyan-elec` and `cdyna-hole` directory with `gaas_cdyna.h5`. This `gaas_cdyna.h5` file contains the cumulated carrier population after the pump pulse.

To compute the carrier concentration, one can use the PERTURBO postprocessing feature `calc_mode = 'dynamics-pp'`.

First link and copy necessary files
```bash
cd cdyna-elec/cnum_check
ln -sf ../../gaas_epr.h5
ln -sf ../gaas_tet.h5
ln -sf ../gaas.temper
vim pp.in
```

`pp.in` contains the following for electrons, 
```
&perturbo
 calc_mode           = 'dynamics-pp'
 prefix              = 'gaas'
 boltz_kdim(1)       = 60
 boltz_kdim(2)       = 60
 boltz_kdim(3)       = 60
 ftemper             = 'gaas.temper'
 boltz_qdim(1)      = 20
 boltz_qdim(2)      = 20
 boltz_qdim(3)      = 20
 hole               = .false.
 band_min           = 5
 band_max           = 5
 boltz_emin         = 5.9             ! eV
 boltz_emax         = 6.7             ! eV
 boltz_de 	    = 1 	      ! meV
/
```

`boltz_de` specifies the grid size in meV for tetrahedron integration over carrier energy.

Run
```bash
perturbo.x -i pp.in | tee pp.out
```
The carrier concentration per unit cell is contained in
```
vim gaas_cdyna.dat
```
Combined with the unit cell volume stored in `gaas_epr.h5/basic_data/volume`, one can compute the carrier concentration.

To match a desired carrier concentration, simply scale the `pump_factor` in `generate_pump_pulse.py`, which is directly proportional to carrier concentration, and regenerate the pump pulse.

### Simulate dynamics with a pump pulse signal

Now we are ready to simulate the dynamics with a pump pulse separately for electrons and holes. We first create an input file based on `pulse.in`.

```bash
cd cdyna-elec
cp pulse.in dyna.in
```

In `dyna.in`, change the following lines
```
 boltz_nstep         = 1000
 output_nstep        = 20
 load_scatter_eph   = .true.
 pump_pulse         = .true.
 pump_pulse_fname   = 'pump_pulse_elec_Epump_01.50.h5'
 yaml_fname = 'gaas_dynamics-run.yml'
```

Note:
  - `time_step` should be the same size as `pump_time_step` in `generate_pump_pulse.py` for first the `pump_time_window` / `time_step` + 1 steps. If one hopes to use a different time step size, first run with the same `time_step` beyond the pulse duration, and restart the calculation with a different `time_step` and `boltz_init_dist = 'restart'`.  
  - `yaml_fname` should changed to prevent overwriting `gaas_dyna_init.yml`.

Run dynamics calculation again with the new input file, where `prefix_cdyna.h5` file will now be over written with the real electron dynamics. This will be a relatively expensive calculation as well.

```bash
perturbo.x -i dyna.in | tee dyna.out
```

Now perform the same calculations for holes.

```bash
cd cdyna-hole
cp pulse.in dyna.in
```

In `dyna.in`, similarly modify the following inputs
```
 boltz_nstep         = 1000
 output_nstep        = 20
 load_scatter_eph   = .true.
 pump_pulse         = .true.
 pump_pulse_fname   = 'pump_pulse_hole_Epump_01.50.h5'
 yaml_fname = 'gaas_dynamics-run.yml'
```

```bash
perturbo.x -i dyna.in | tee dyna.out
```

### Result analysis with Perturbopy

To analyze the population with PerturboPy, see the [PERTURBO 2023 workshop tutorial](https://github.com/perturbo-code/perturbo-workshop-2023/tree/main/Tutorial4). In this example, we will focus on computing the transient absorption signal.

Create the following python script `compute_trans_abs.py` in the parent directory, which is also contained in this repository.

```python
import perturbopy.postproc as ppy

prefix = 'gaas'
elec_folder = 'cdyna-elec'
hole_folder = 'cdyna-hole'

# Define electron paths and DynaRun object
cdyna_elec_path = f'{elec_folder}/{prefix}_cdyna.h5'
elec_tet_path = f'{elec_folder}/{prefix}_tet.h5'
elec_yaml_path = f'{elec_folder}/{prefix}_dynamics-run.yml'
elec_dyna_run = ppy.DynaRun.from_hdf5_yaml(cdyna_elec_path, elec_tet_path, elec_yaml_path)

# Define hole paths and DynaRun object
cdyna_hole_path = f'{hole_folder}/{prefix}_cdyna.h5'
hole_tet_path = f'{hole_folder}/{prefix}_tet.h5'
hole_yaml_path = f'{hole_folder}/{prefix}_dynamics-run.yml'
hole_dyna_run = ppy.DynaRun.from_hdf5_yaml(cdyna_hole_path, hole_tet_path, hole_yaml_path)

# Print the PumpPulse info
ppy.utils.spectra_trans_abs.compute_trans_abs(elec_dyna_run,
    hole_dyna_run,
    de_grid=0.02,
    energy_grid_max=None,
    eta=0.02,
    save_npy=True,
    tr_dipoles_sqr=None)

elec_dyna_run.close_hdf5_files()
hole_dyna_run.close_hdf5_files()
```

The following command will generate and save the transient absorption signal A(E,t) with dimension (number of energy grid points, number of time snapshots).

```bash
python3 compute_trans_abs.py
```

To visualize the transient absorption signals, one can use a simple plot script `plot_trans_abs.py` utilizing PerturboPy functions

```python
import numpy as np
import perturbopy.postproc as ppy
import matplotlib.pyplot as plt

pump_energy = '1.5000'

# Load the numpy binary files
dA_elec = np.load(f'trans_abs_dA_elec_Epump_{pump_energy}.npy') # electron signal
dA_hole = np.load(f'trans_abs_dA_hole_Epump_{pump_energy}.npy') # hole signal
time_grid = np.load(f'trans_abs_T_Epump_{pump_energy}.npy')
energy_grid = np.load(f'trans_abs_E_Epump_{pump_energy}.npy')

# Plot the total transient absorption
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ppy.utils.spectra_plots.plot_trans_abs_map(ax, time_grid, energy_grid, dA_elec + dA_hole, vmin=-2.5, vmax=0.0)
ax.set_title(f'Transient absorption in GaAs, Epump={pump_energy} eV')
ax.set_ylim(ymax=2.0)
plt.savefig('transient_absorption_total.png', dpi=400)
```

```bash
python3 plot_trans_abs.py
```

![Transient Absorption](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on5/images/transient_absorption_total.png)

The grid sizes used in this example does not lead to a converged result. The following plot shows a converged result using `boltz_kdim(1) = boltz_kdim(2) = boltz_kdim(3) = 180` and `boltz_qdim(1) = boltz_qdim(2) = boltz_qdim(3) = 90` with the same `delta_smear` and `phfreq_cutoff` as in the example.

![Converged_absorption](https://github.com/perturbo-code/perturbo-workshop-2025/blob/main/Hands-on5/images/converged_trans_abs.png)

---
>#### Initial populations at t = 0 as a defined distribution (no pump pulse)
>With `calc_mode = 'dynamics-run'`, one can specify a Fermi-Dirac distribution, a Gaussian or Lorentzian excitation at initial time t = 0 with the input parameters `boltz_init_dist`, `boltz_init_e0`, `boltz_init_smear` and `boltz_init_ampl`. 
>
>```
>&perturbo
> calc_mode           = 'dynamics-run'
> ...
> boltz_init_dist     = 'fermi'
> boltz_init_e0      = -4.2281057129   ! eV
> boltz_init_smear   = 146.57          ! meV
> boltz_init_ampl    = 1.0             !
> ...
>```
>
>`boltz_init_ampl` is a new parameter in the PERTURBO v3.0 release that specifies a factor multiplied by initial distribution. No initial carrier distributions are generated with  `boltz_init_dist` if `boltz_init_ampl = 0`.
>
>For details, see PERTURBO website [tutorial](https://perturbo-code.github.io/mydoc_dynamics.html#zero-field-ultrafast-dynamicscalc_mode--dynamics-run) or past [workshop](https://github.com/perturbo-code/perturbo-workshop-2023/tree/main/Tutorial4) on ultrafast dynamics.

---

### Additional Resources
To explore more features of ultrafast dynamics, refer to PERTURBO tutorials on [Ultrafast dynamics](https://perturbo-code.github.io/mydoc_dynamics.html#zero-field-ultrafast-dynamicscalc_mode--dynamics-run) and [Ultrafast spectroscopy](https://perturbo-code.github.io/mydoc_spectroscopy.html).