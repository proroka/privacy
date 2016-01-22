# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:44:31 2016

@author: amanda
"""

from solve_aggregation import *


# default values
# species
z1 = 120
z2 = 130
z4 = 100
# rates
a1 = 1.0
a2 = 1.0
b1 = 1.0
b2 = 1.0



p = get_stationary_distr(z1, z2, z4, a1, a2, b1, b2)
p_obs = BuildObservableDistr(p)

fig = plot_distr(p) 
plt.show() 
#plot_distr(p_obs, ('a', 'b', 'c')), plt.show() 


figname = 'plots/agg_distribution.eps'
fig.savefig(figname) 