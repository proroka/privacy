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
from class_species import *

                  

#----------------------------------------------------------------------
# Create robot community

verbose = False

# Global settings
tmax = 100
num_species = 3
num_robots = [10, 10, 10];

# set the transitioning probabilities
prob = [];
ps = np.array([[0.6, 0.8], [0.7, 0.9], [0.5, 0.5]]) # p00 and p_wait for each species
for s in range(num_species):
    temp = build_probabilities(num_states, ps[s,0], ps[s,1])
    prob.append(temp)

# states: 0 is search, 1-M are tasks
num_tasks = 30;
num_states = num_tasks + 1;

# robots are either searching or waiting to collaborate
num_searching = np.zeros((tmax, num_species))

# generate wait states
wait_states = list(subsets(range(num_species)))
# remove empty set and full set
wait_states = wait_states[1:-1]
num_wait_states = len(wait_states)
num_waiting = np.zeros((tmax, num_wait_states))


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

    # summarize number of robots in search
    for s in range(num_species):
        num_searching[t,s] = np.sum(species[s].robots[:,0])
        
    # summarize number of robots waiting to collaborate    
    for m in range(1, num_states):
        s_list = []
        # get list of speices waiting at this task
        for s in range(num_species):
            if (np.sum(species[s].robots[:,m])) > 0:
                s_list.append(s)
        # find wait-state and increment
        if s_list in wait_states:
            ind = wait_states.index(s_list)
            num_waiting[t,ind] =  num_waiting[t,ind] + 1 




#----------------------------------------------------------------------   
# Plots
 
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

    fig = plt.figure(figsize=(6,4))
    plt.plot(num_waiting[:,0],'blue')
    plt.plot(num_waiting[:,1],'red')
    plt.plot(num_waiting[:,3],'green')
    plt.plot(num_waiting[:,4],'cyan')
    #plt.ylim([0, 1])
    plt.show()    