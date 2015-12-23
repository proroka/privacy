# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:53:42 2015
@author: amandaprorok

"""

# Standard modules
import numpy as np
#import scipy as sp
#import pylab as pl
import matplotlib.pyplot as plt
#import sys
#import time
#import pickle

# My modules
from utility_functions import *

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
            if (curr != new): # move robot back to search
                if new == 0:
                    self.robots[r,curr] = 0
                    self.robots[r,new] = 1
                else: # check if any other robots present at task
                    if (np.sum(self.robots[:,new]) < 1):
                        self.robots[r,curr] = 0
                        self.robots[r,new] = 1
                    
    # move numerous robots back to search                
    def transition_to_search_n(self, task, num): 
        ind = np.where(self.robots[:,task])[0]
        for i in range(num):
            rand = np.random.randint(0,np.size(ind))
            self.robots[ind[rand],task] = 0
            self.robots[ind[rand],0] = 1
            
    # move only robot back to search        
    def transition_to_search(self, task): 
        ind = np.where(self.robots[:,task])[0][0]
        self.robots[ind,task] = 0
        self.robots[ind,0] = 1                  
 
                   

#----------------------------------------------------------------------
# Create robot community

verbose = False

# Global settings
tmax = 100
num_species = 2
num_robots = [20, 10];
# states: 0 is search, 1-M are tasks
num_tasks = 30;
num_states = num_tasks + 1;

# robots are either searching or waiting to collaborate
num_wait_states = np.power(2, num_species) - 2 # power set without empty and complete sets
num_searching = np.zeros((tmax, num_species))
num_waiting = np.zeros((tmax, num_wait_states))


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
        # get collaborative tasks (1 of each species needed)
        num_collab = int(min(at_task))
        #print num_collab
        if (num_collab > 0 and m!=0):
            if verbose:
                print 'At state: ', m, 'moving ', num_collab, 'robot back to search'
            # get number of robots per species to be removed
            #new_at_task = at_task - min(at_task)
            # for each species, update states (put back into search)
            for s in range(num_species):
                species[s].transition_to_search(m)
            
    if verbose: print '***'

    # summarize 
#    for s in range(num_species):
#        num_searching[t,s] = np.sum(species[s].robots[:,0])
#    for i in range(num_wait_states):
#        num_waiting[t,i]
    




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
