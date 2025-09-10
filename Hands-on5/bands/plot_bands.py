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

emin_holes, emax_holes = 4.8, 5.3 # eV
emin_electrons, emax_electrons = 5.9, 6.7 # eV

xmin, xmax = ax.get_xlim()
ax.fill_between([xmin, xmax], [emin_holes, emin_holes],
                [emax_holes, emax_holes], color='blue', alpha=0.2)
ax.fill_between([xmin, xmax], [emin_electrons, emin_electrons],
                [emax_electrons, emax_electrons], color='red', alpha=0.2)
ax.set_xlim(xmin, xmax)

plt.savefig('bands.png', dpi=200)
