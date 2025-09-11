import perturbopy.postproc as ppy
import matplotlib.pyplot as plt

# Example using the imsigma calculation mode
diam_ephmat_spin = ppy.EphmatSpin.from_yaml('diam_ephmat_spin.yml')

plt.rcParams.update(ppy.plot_tools.plotparams)
fig, ax = plt.subplots()

diam_ephmat_spin.qpt.add_labels(ppy.lattice.points_fcc)

diam_ephmat_spin.plot_ephmat(ax)

plt.show()

plt.savefig('diam-ephmat-spin.png')


