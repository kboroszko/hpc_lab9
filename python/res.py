import numpy as np
from matplotlib import pyplot as plt
#%%
seq = np.genfromtxt('results/seq/all.dat', delimiter='\t')
noopt2  = np.genfromtxt('results/noopt/par2.err', delimiter='\t')
noopt5  = np.genfromtxt('results/noopt/par5.err', delimiter='\t')
noopt12 = np.genfromtxt('results/noopt/par12.err', delimiter='\t')
noopt24 = np.genfromtxt('results/noopt/par24.err', delimiter='\t')
noopt48 = np.genfromtxt('results/noopt/par48.err', delimiter='\t')
#%%
plt.plot(seq[:,0], seq[:, 1], label='seq')
plt.plot(noopt2[:, 0], noopt2[:, 2], label='noopt2')
plt.legend()
plt.show()
