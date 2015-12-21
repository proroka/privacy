# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:53:42 2015
@author: amandaprorok

"""


import numpy as np
import scipy as sp
import pylab as pl
import matplotlib.pyplot as plt
import sys
import time
import pickle


#----------------------------------------------------------------------
# Define classes and utility functions

class Species:
    def __init__(self, N, M, prob):
        self.num_robots = N
        self.num_states = M
        self.prob = prob

        
        self.robots = np.zeros((self.num_robots, self.num_states))
        
        # all robots start in state 1 (search)
        self.robots[:,0] = np.ones((self.num_robots,))


# input: a row from transition probability matrix
def pick_transition(p):
    rand = np.random.rand(1)
    v = 0
    for i in range(np.size(p)):
        v += p[i]
        if rand <= v:
            return i
    # Should not happen (unless probabilities do not sum to 1 exactly).
    return np.size(p) - 1



#----------------------------------------------------------------------
# Create robot community

num_species = 2
num_robots = [10, 20];
num_states = 2;
# set the transitioning probabilities
prob = [];
prob_s0 = np.array([[0.8, 0.2], [0.01, 0.99]]); # species 0
prob_s1 = np.array([[0.8, 0.2],[0.1, 0.9]]); # species 1

prob.append(prob_s0)
prob.append(prob_s1)

# Create species
species = []
for s in range(num_species):
    st = Species(num_robots[s], num_states, prob[s])
    species.append(st)
    
 
tmax = 100
evol = np.zeros((num_species, tmax))
   
# Run simulation
for t in range(tmax):
    
    # run robots
    for s in range(num_species):
        p = species[s].prob
        robots = species[s].robots
        
        # iterate through robots
        for r in range(species[s].num_robots):
            # get current state
            curr = np.where(robots[r,:])[0][0]
            new = pick_transition(p[curr,:])
            if curr != new: # move robot
                    robots[r,curr] = 0
                    robots[r,new] = 1
        state = 0                
        evol[s,t] = np.mean(species[s].robots[:,state])
        
    # update system state


fig = plt.figure()
plt.plot(evol[0,:],'blue')
plt.plot(evol[1,:],'red')

plt.show()    
    
