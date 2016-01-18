import math
import matplotlib.pyplot as plt
import numpy as np
import itertools
import pickle
import matplotlib.colors as colors

import seaborn as sns
from scipy.stats import kendalltau

run = 3
save_epsilon = 'data/epsilon_data_run' + str(run) + '.p'
save_params = 'data/params_run' + str(run) + '.p'

data = pickle.load(open(save_epsilon, 'r'))
params = pickle.load(open(save_params, 'r'))
# fix this
z1 = params['z1']
z2 = params['z2']
z4 = params['z4']
s_range = params['s_range']
r_range = params['r_range']
sim = params['sim']

if sim=='SPECIES':
    par_range = s_range
elif sim=='RATES':
    par_range = r_range


name = 'Purples'
name = 'RdPu'

cmap=plt.get_cmap(name)

fig = plt.figure(figsize=(5,4))

dx = (np.max(par_range) - np.min(par_range)) / (2. * float(len(par_range) - 1))
plt.imshow(data,cmap=cmap,interpolation='nearest', origin='lower',
    extent=[np.min(par_range) - dx, np.max(par_range) + dx, np.min(par_range) - dx,
    np.max(par_range) + dx])
#plt.xlabel('%s_2' % args.sweep_type)
#plt.ylabel('%s_3' % args.sweep_type)
plt.colorbar()
plt.show()


fig = plt.figure(figsize=(5,4))
y, x = np.meshgrid(par_range, par_range)                
pcm = plt.pcolor(x, y, data,
                   norm=colors.PowerNorm(gamma=2.5),
                   cmap=cmap)   
plt.colorbar(pcm)
plt.tight_layout()
plt.axis([x.min(), x.max(), y.min(), y.max()])      

plt.show()



