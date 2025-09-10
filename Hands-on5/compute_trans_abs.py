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
