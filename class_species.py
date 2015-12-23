# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 11:45:11 2015

@author: amandaprorok
"""

import numpy as np
import matplotlib.pyplot as plt

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
        return np.sum(self.robots[:,state])


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