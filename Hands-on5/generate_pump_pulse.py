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
tr_dipoles=None 		# transient dipole squared in .npy file
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
