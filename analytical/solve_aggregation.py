from sympy.solvers import solve
from sympy import Symbol
import math
import matplotlib.pyplot as plt
import numpy as np
import itertools
import pickle


a1 = 1.0
a2 = 1.0
b1 = 1.0
b2 = 1.0

_EPSILON = 1e-5

#-----------------------------------------------------------
# Define functions

# def get_stationary_c(z1, z2, z4, a1, a2, b1, b2):
#     #x3 = Symbol('x3', real=True, nonnegative=True)
#     x3 = Symbol('x3')
#     solutions = solve(
#         [
#             a1*(z1-x3-z4+(b2*z4 / (b1*x3 + b2)))*(z2-x3-z4+(b2*z4 / (b1*x3 + b2))) - a2*x3
#         ],
#         dict=True)
#
#     print solutions
#
#     cnt = 0
#     for s in solutions:
#         _x3 = float(s[x3])
#         x4 = b2 * z4 / (b1 * _x3 + b2)
#         x5 = z4 - x4
#         x1 = z1 - _x3 - x5
#         x2 = z2 - _x3 - x5
#
#         print x1
#         print x2
#         print x4
#         print x5
#
#         if (x1>=0. and x2>=0. and x4>=0. and x5>=0.):
#             cnt += 1
#             sol = [x1, x2, _x3, x4, x5]
#
#     #print sol
#     if cnt == 1:
#         return sol
#     return None

def get_stationary_c(z1, z2, z4, a1, a2, b1, b2):
    x3 = Symbol('x3')
    solutions = solve(
        [
            a1*(z1-x3-z4+(b2*z4 / (b1*x3 + b2)))*(z2-x3-z4+(b2*z4 / (b1*x3 + b2))) - a2*x3
        ],
        dict=True)

    cnt = 0
    for s in solutions:
        #print s
        _x3 = s[x3]
        x3_value = None
        re, im = _x3.as_real_imag()
        if float(im) < 1e-6:
            x3_value = float(re)
    
        if x3_value is None or x3_value < 0.:
            continue

        x4 = b2 * z4 / (b1 * x3_value + b2)
        x5 = z4 - x4
        x1 = z1 - x3_value - x5
        x2 = z2 - x3_value - x5

        if (x1>=0. and x2>=0. and x4>=0. and x5>=0.):
            cnt += 1
            sol = [x1, x2, x3_value, x4, x5]

    if cnt == 1:
        return sol
    print 'Oops. Number of solutions:', cnt
    return None
        
        
def get_stationary_distr(z1, z2, z4, a1, a2, b1, b2):
    sol = get_stationary_c(float(z1), float(z2), float(z4), a1, a2, b1, b2)    
    if sol is None: return None
    
    prob = dict()
    
    # compute all possible values for x_i according to rule:
    # x3 == [0, min(z1, z2)]
    # x4 == [0, z4]
    sum_p = 0
    for x3 in range(min(z1,z2) + 1):
        for x4 in range(z4 + 1):
            x5 = z4 - x4
            x1 = z1 - x3 - x5
            x2 = z2 - x3 - x5
            # need to check that all xi are positive, else skip this loop
            if x1<0 or x2<0 or x5<0: continue
            currx = [x1, x2, x3, x4, x5]
            
            # compute formula for stationary distribution f(c_i, x_i):
            p = 1
            for ci, xi in zip(sol, currx):
                #p *= math.pow(ci,xi) / math.factorial(xi) 
                for i in range(1, xi + 1):
                    p *= ci / float(i)
                    
            # return dict with key(x_i values) value(eval pi(x_i values))
            prob[(x1, x2, x3, x4, x5)] = p
            sum_p += p
    
    for k in prob:
        prob[k] /= sum_p
                   
    return prob 
 
 

def plot_distr(prob, species_name=('a', 'b', 'c', 'd', 'e')):
    def _GetColors(n):
       cm = plt.get_cmap('gist_rainbow')
       return [cm(float(i) / float(n)) for i in range(n)]

    min_population = float('inf')
    max_population = -1.
    for k in prob:
       for s in k:
           max_population = max(s, max_population)
           min_population = min(s, min_population)
    colors = _GetColors(len(species_name))
    fig = plt.figure(figsize=(20,5))
    for i in xrange(len(species_name)):
        plt.subplot(1,len(species_name),i+1)
        values = np.zeros((max_population - min_population + 1, ))
        for k, v in prob.iteritems():
            values[k[i] - min_population] += v
        plt.bar(np.arange(min_population - 0.5, max_population - 0.4, 1.), values, width=1.,
               label=species_name[i], color=colors[i], alpha=0.3)
        plt.legend(loc='upper right', shadow=False, fontsize='x-large')
        plt.xlabel('Population')
        plt.ylabel('Probability')
    plt.tight_layout()
    #return fig, ax
 
# Compares two distributions computed by BuildDistribution().
def CompareDistributions(A, B, prune=_EPSILON):
   all_values = []
   for k, v in A.iteritems():
       if k in B:
           all_values.append(np.abs(np.log(A[k] + prune) - np.log(B[k] + prune)))
       elif k not in B:
           all_values.append(np.abs(np.log(A[k] + prune) - np.log(prune)))
       # Ignore when one is larger and the other one is smaller or when both are smaller.
   for k, v in B.iteritems():
       if k not in A:
           all_values.append(np.abs(np.log(prune) - np.log(B[k] + prune)))
   return max(all_values)
    
def BuildObservableDistr(p):
    new_p = {}
    for k, v in p.iteritems():
        new_k = (k[0]+k[1]+k[3], k[2]+k[4])
        if new_k not in new_p:
            new_p[new_k] = 0.
        new_p[new_k] += v
    return new_p
    
def MaxLeakage(pop, p_obs):
    # number of epsilon values
    epsilons = []
    # loop for all adjacent databases for a given population
    for change_ind in itertools.permutations(range(3), 2): # first: +1, second: -1
        bt = np.array([0, 0, 0])
        bt[change_ind[0]] = 1
        bt[change_ind[1]] = -1
        # adjacent population
        pop_adj = pop + bt
        
        p_adj = get_stationary_distr(pop_adj[0], pop_adj[1], pop_adj[2], a1, a2, b1, b2)
        if p_adj is None: continue
        p_adj_obs = BuildObservableDistr(p_adj)
        
        max_diff = CompareDistributions(p_obs, p_adj_obs, prune=_EPSILON)
        epsilons.append(max_diff)

        print 'Adj', pop_adj
        print max_diff
        
    if not epsilons: return np.nan  
    return max(epsilons)
       
              
#-----------------------------------------------------------
# Main
plot_on = False
    
#p = get_stationary_distr(z1, z2, z4, a1, a2, b1, b2)
#p_obs = BuildObservableDistr(p)
#if plot_on: plot_distr(p), plt.show() 
#if plot_on: plot_distr(p_obs, ('a', 'b', 'c')), plt.show() 


run = 4
sim = 'RATES' #sim = {'RATES', 'SPECIES'}

z1 = 220
z2 = 220
z4 = 200

s_range = range(100,301,10)
r_range = np.arange(0.1, 2.0, 0.1)
#r_range = np.arange(0.1, 0.3, 0.1)


if sim=='RATES':
    num_pars = len(r_range)
elif sim=='SPECIES':
    num_pars = len(s_range)

epsilon = np.zeros((num_pars, num_pars))
params = {'z1': z1,'z2': z2, 'z4': z4, 's_range': s_range, 'r_range': r_range, 'sim': sim}

save_params = 'data/params_run' + str(run) + '.p'
save_epsilon = 'data/epsilon_data_run' + str(run) + '.p'


# loop for varying parameters
for i in range(num_pars):
    if sim=='RATES':
        a1 = r_range[i]
        a2 = r_range[i]
        range_pars = range(num_pars) 
    elif sim=='SPECIES':
        z1 = s_range[i]
        range_pars = range(i,num_pars) 
 
    for j in range_pars:
        if sim=='RATES':
            b1 = r_range[j]
            b2 = r_range[j]
        elif sim=='SPECIES':
            z2 = s_range[j] 
        
        pop = np.array([z1, z2, z4])
        if pop[0]<1 or pop[1]<1 or pop[2]<1: print "Species cannot be 0"
        if sim=='SPECIES':
            print "Solving population (", i*num_pars+j, '/', num_pars**2,'): ', pop
        elif sim=='RATES':
            print "Solving for rates (", i*num_pars+j, '/', num_pars**2,'): ', a1, a2, b1, b2
        p = get_stationary_distr(pop[0], pop[1], pop[2], a1, a2, b1, b2)
        if p is None: 
            epsilon[i,j] = np.nan
            if sim=='SPECIES':
                epsilon[j,i] = np.nan
            continue
            
        p_obs = BuildObservableDistr(p)
        # builds all adjacent populations and gets max leakage
        epsilon[i,j] = MaxLeakage(pop, p_obs)
        if sim=='SPECIES':
            epsilon[j,i] = epsilon[i,j]
        
pickle.dump(epsilon, open(save_epsilon, 'w'))
pickle.dump(params, open(save_params, 'w'))
print epsilon














