import matplotlib.pyplot as plt
import numpy as np
import pickle

legends = [
    'N = [10, 7, 4]',
    'beta = [1, 2, 0.4]',
    'alpha = [1, 2, 0.4]',
    'Identical',
]

paths = [
    'data/data_tau_111_111_100704.bin',
    'data/data_tau_111_1204_101010.bin',
    'data/data_tau_1204_111_101010.bin',
    'data/data_tau_111_111_101010.bin',
]

cm = plt.get_cmap('gist_rainbow')
colors = [cm(float(i) / float(len(paths))) for i in range(len(paths))]


plt.figure()
for i, p in enumerate(paths):
    min_epsilons = pickle.load(open(p, 'r'))
    values = np.array(sorted(min_epsilons.iteritems()))
    plt.plot(values[:, 0], values[:, 1], linewidth=2, color=colors[i], label=legends[i])
plt.xlabel('tau')
plt.ylabel('Minimum acheivable \epsilon')
plt.legend(loc=2)
plt.show()
