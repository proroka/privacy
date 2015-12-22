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
# Define class

class Species:
    def __init__(self, N, M, prob):
        self.num_robots = N
        self.num_states = M
        self.prob = prob

        
        self.robots = np.zeros((self.num_robots, self.num_states))
        
        # all robots start in state 1 (search)
        self.robots[:,0] = np.ones((self.num_robots,))


    def get_state(self, state):
        return np.mean(self.robots[:,state])


    def update_states(self):
        # iterate through robots
        for r in range(self.num_robots):
            # get current state
            curr = np.where(self.robots[r,:])[0][0] # where() returns list of arrays
            new = pick_transition(self.prob[curr,:])
            if curr != new: # move robot
                    self.robots[r,curr] = 0
                    self.robots[r,new] = 1
                    
    def transition_to_search(self, task, num): 
        ind = np.where(self.robots[:,task])[0]
        for i in range(num):
            rand = np.random.randint(0,np.size(ind))
            self.robots[ind[rand],task] = 0
            self.robots[ind[rand],0] = 1
            
                       
 
                   
#----------------------------------------------------------------------
# Define utility functions
                   
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


def build_probabilities(num_states, p00, p_wait):
    p = np.zeros((num_states, num_states))
    # first row: search to tasks
    p0x = (1. - p00) / (num_states - 1.)
    p[0, 0] = p00
    p[0,1:] = p0x
    # tasks
    p[1:,1:] = p[1:,1:] + np.eye(num_states -1) * p_wait
    p_fail = 1 - p_wait
    p[1:,0] = p_fail
    
    return p


#----------------------------------------------------------------------
# Create robot community

verbose = False

num_species = 2
num_robots = [20, 10];
# states: 0 is search, 1-M are tasks
num_tasks = 30;
num_states = num_tasks + 1;

# set the transitioning probabilities
prob = [];

p00 = 0.6
p_wait = 0.8
prob_s0 = build_probabilities(num_states, p00, p_wait)
p00 = 0.7
p_wait = 0.9
prob_s1  = build_probabilities(num_states, p00, p_wait)

prob.append(prob_s0)
prob.append(prob_s1)

# Create list of species
species = []
for s in range(num_species):
    st = Species(num_robots[s], num_states, prob[s])
    species.append(st)

# Check that correctly set   
if verbose: 
    for m in range(0,num_states):
        n = []
        for s in range(num_species):
            nt = np.sum(species[s].robots[:,m])
            print 'State: ', m, 'Species: ', s, 'Sum: ', nt
            n.append(nt)
    print '***'
 
tmax = 100
evol = np.zeros((num_species, tmax, num_states))

   
# Run simulation
for t in range(tmax):
    
    # run robots
    for s in range(num_species):        
        species[s].update_states()
        for a in range(num_states):               
            evol[s,t,a] = species[s].get_state(a)  
                   
    # collaboration dynamics: move collaborative robots back to search
    for m in range(num_states):
        at_task = [] # for given task, all sum of robots per species
        for s in range(num_species):
            temp = np.sum(species[s].robots[:,m])
            if verbose:
                print 'T:',t, 'State: ', m, 'Species: ', s, 'Sum: ', temp
            at_task.append(temp)
        #print at_task    
        # get collaborative tasks (1 of each species needed)
        num_collab = int(min(at_task))
        #print num_collab
        if (num_collab > 0 and m!=0):
            if verbose:
                print 'At state: ', m, 'moving ', num_collab, 'robots back to search'
            # get number of robots per species to be removed
            #new_at_task = at_task - min(at_task)
            # for each species, update states (put back into search)
            for s in range(num_species):
                species[s].transition_to_search(m, num_collab)
            

    if verbose: print '***'



#----------------------------------------------------------------------    
plot_on = True

if plot_on:
    fig = plt.figure(figsize=(6,4))
    st = 0
    plt.plot(evol[0,:,st],'blue')
    plt.plot(evol[1,:,st],'red')
    plt.ylim([0, 1])
    plt.show()    

    fig = plt.figure(figsize=(6,4))
    st = 1
    plt.plot(evol[0,:,st],'blue')
    plt.plot(evol[1,:,st],'red')
    plt.ylim([0, 1])
    plt.show()        
