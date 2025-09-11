import perturbopy.postproc as ppy
import matplotlib.pyplot as plt

# Example using the imsigma calculation mode
diam_imsigma_flip = ppy.ImsigmaSpin.from_yaml('diam_imsigma_spin.yml')


plt.rcParams.update(ppy.plot_tools.plotparams)
fig = plt.figure()
ax = fig.add_subplot(111)

#Bands 1 and 2
ax.scatter(diam_imsigma_flip.bands[1][:],diam_imsigma_flip.imsigma[1][1][:],color='b')
ax.scatter(diam_imsigma_flip.bands[2][:],diam_imsigma_flip.imsigma[1][2][:],color='b')

ax.set_xlim(17.4,18.0)
ax.set_xlabel('Energy(eV)')
ax.set_ylabel('ImSigma(meV)')
plt.savefig('diam-imsigma-flip.png')


