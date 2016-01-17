from sympy.solvers import solve
from sympy import Symbol
import math
import matplotlib.pyplot as plt
import numpy as np
import itertools

a1 = 1.0
a2 = 1.0
b1 = 1.0
b2 = 1.0

z1 = 5
z2 = 6
z4 = 8

_EPSILON = 1e-5

#-----------------------------------------------------------
# Define functions

def get_stationary_c(z1, z2, z4, a1, a2, b1, b2):
    x3 = Symbol('x3', real=True, nonnegative=True)
    solutions = solve(
        [
            a1*(z1-x3-z4+(b2*z4 / (b1*x3 + b2)))*(z2-x3-z4+(b2*z4 / (b1*x3 + b2))) - a2*x3
        ],
        dict=True)

    #print solutions
    
    cnt = 0
    for s in solutions: 
        _x3 = float(s[x3])
        x4 = b2 * z4 / (b1 * _x3 + b2)
        x5 = z4 - x4
        x1 = z1 - _x3 - x5
        x2 = z2 - _x3 - x5
        
        if (x1>=0. and x2>=0. and x4>=0. and x5>=0.):
            cnt += 1
            sol = [x1, x2, _x3, x4, x5]

    if cnt == 1:        
        return sol
    return None
        
        
def get_stationary_distr(z1, z2, z4, a1, a2, b1, b2):
    sol = get_stationary_c(float(z1), float(z2), float(z4), a1, a2, b1, b2)    
    #print sol
    
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
                
                #p *= ci^xi / xi! * exp(-ci):
                p *= math.pow(ci,xi) / math.factorial(xi) #* math.exp(-ci))
                
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
        if k in B and (v > prune or B[k] > prune):
            all_values.append(np.abs(np.log(A[k]) - np.log(B[k])))
        elif k not in B and v > prune:
            all_values.append(float('inf'))
        # Ignore when one is larger and the other one is smaller or when both are smaller.
    for k, v in B.iteritems():
        if k not in A and v > prune:
            all_values.append(float('inf'))
    return max(all_values)


def choose_iter(elements, length):
    for i in xrange(len(elements)):
        if length == 1:
            yield (elements[i],)
        else:
            for next in choose_iter(elements[i+1:len(elements)], length-1):
                yield (elements[i],) + next
def choose(l, k):
    return list(choose_iter(l, k))
    
#-----------------------------------------------------------
# Main
plot_on = False
    
p = get_stationary_distr(z1, z2, z4, a1, a2, b1, b2)

if plot_on: plot_distr(p), plt.show()

# loop for varying populations
# loop for all adjacent databases for a given population
ind = [0, 1, 2] # z1, z2, z4
base = np.array([0, 0, 0])
pop = np.array([z1, z2, z4])
# no species can be 0
# CHECK

# number of epsilon values: number of permutations
epsilons = []#np.zeros(num_species * (num_species-1))

for change_ind in itertools.permutations(ind, 2): # first: +1, second: -1
    bt = base.copy()
    bt[change_ind[0]] += 1
    bt[change_ind[1]] -= 1
    # adjacent population
    pop_adj = pop + bt
    # check if any negative
    print pop_adj
    
    p = get_stationary_distr(z1, z2, z4, a1, a2, b1, b2)
    p_adj = get_stationary_distr(pop_adj[0], pop_adj[1], pop_adj[2], a1, a2, b1, b2)
    
    max_diff = CompareDistributions(p, p_adj, prune=_EPSILON)
    max_diff = CompareDistributions(p, p, prune=_EPSILON)
    print max_diff
    
    epsilons.append(max_diff)
    
print epsilons













