import perturbopy.postproc as ppy
import matplotlib.pyplot as plt

# Example using the imsigma calculation mode
diam_spins = ppy.Spins.from_yaml('diam_spins.yml')

plt.rcParams.update(ppy.plot_tools.plotparams)
fig, ax = plt.subplots()

diam_spins.kpt.add_labels(ppy.lattice.points_fcc)

diam_spins.plot_spins(ax, log = False)

plt.show()

plt.savefig('diam-spins.png')


