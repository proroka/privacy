import math
import matplotlib.pyplot as plt
import numpy as np
import itertools
import pickle
import matplotlib.colors as colors


run = 5

save_epsilon = 'data/epsilon_data_run' + str(run) + '.p'
save_params = 'data/params_run' + str(run) + '.p'

data = pickle.load(open(save_epsilon, 'r'))
params = pickle.load(open(save_params, 'r'))

# fix this
if run!=2:
    z1 = params['z1']
    z2 = params['z2']
    r_range = params['r_range']
    sim = params['sim']
z4 = params['z4']
s_range = params['s_range']


if run==2:
    par_range = s_range
elif sim=='SPECIES':
    par_range = s_range
elif sim=='RATES':
    par_range = r_range


name = 'Purples'
name = 'YlGnBu'
name = 'RdPu'


cmap=plt.get_cmap(name)

fig = plt.figure(figsize=(5,4))

dx = (np.max(par_range) - np.min(par_range)) / (2. * float(len(par_range) - 1))
if run==2:
    plt.imshow(data,cmap=cmap,interpolation='nearest', origin='lower',
    extent=[np.min(par_range) - dx, np.max(par_range) + dx, np.min(par_range) - dx,
    np.max(par_range) + dx],norm=colors.PowerNorm(gamma=2.5), clim=(2.,11.))
else:
    plt.imshow(data,cmap=cmap,interpolation='nearest', origin='lower')#,
#        extent=[np.min(par_range) - dx, np.max(par_range) + dx, np.min(par_range) - dx,
#        np.max(par_range) + dx])

plt.tight_layout()
plt.colorbar()
plt.show()
