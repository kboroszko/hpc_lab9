import numpy as np
from matplotlib import pyplot as plt
#%%
seq = np.genfromtxt('data/seq/all.dat', delimiter='\t')

noopt  = [np.genfromtxt('data/noopt/par{}.err'.format(i), delimiter='\t') for i in [2,5,12,24,48] ]
allred  = [np.genfromtxt('data/allred/par{}.err'.format(i), delimiter='\t') for i in [2,5,12,24,48] ]
asyn  = [np.genfromtxt('data/async/par{}.err'.format(i), delimiter='\t') for i in [2,5,12,24,48] ]

#%%
plt.plot(seq[:,0], seq[:, 1], label='seq')
# nodes = [2,5,12,24,48]
# for i in range(len(noopt)):
#     arr = noopt[i]
#     lab = "noopt{}".format(nodes[i])
#     plt.plot(arr[:, 0], arr[:, 2], label=lab)
plt.plot(noopt[0][:, 0], noopt[0][:, 2], label='n=2')
plt.plot(noopt[1][:, 0], noopt[1][:, 2], label='n=5')
plt.plot(noopt[3][:, 0], noopt[3][:, 2], label='n=24')
plt.legend()
plt.show()

#%%

n1000 = np.genfromtxt('data/noopt/n1000.dat', delimiter='\t')

#%%
plt.figure()
plt.plot(n1000[:,0], n1000[:, 1])
plt.plot()

#%%
plt.plot(seq[:,0], seq[:, 1], label='seq')
# plt.plot(noopt[0][:, 0], noopt[0][:, 2], label='noopt n=2')
# plt.plot(asyn[0][:, 0], asyn[0][:, 2], label='async n=2')
plt.plot(allred[0][:, 0], allred[0][:, 2], label='allred n=2')
plt.legend()
plt.show()