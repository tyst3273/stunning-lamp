import numpy as np
import matplotlib.pyplot as plt


freq = np.loadtxt('freqs')

nb = freq.shape[1]
nq = freq.shape[0]

for bb in range(nb):

    plt.plot(freq[:,bb],marker='o',ms=2,mfc='w',mec='k',mew=0.5,ls='')

plt.show()













