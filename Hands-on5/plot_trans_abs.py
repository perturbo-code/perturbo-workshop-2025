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
