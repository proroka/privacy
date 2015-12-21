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
# Define classes

class Species:
    def __init__(self, N, M, prob):
        self.num_robots = N
        self.num_states = M
        self.prob = prob

        
        self.robots = np.zeros((self.num_robots, self.num_states))
        # all robots start in state 1 (search)
        self.robots[:,0] = np.ones((self.num_robots,))



#----------------------------------------------------------------------
# Create robot community

num_species = 2
num_robots = [10, 20];
num_states = 2;
# set the transitioning probabilities
prob = [];
prob_s0 = np.array([[0.8, 0.2], [0.01, 0.99]]); # species 0
prob_s1 = np.array([[0.8, 0.2],[ 0.1, 0.9]]); # species 1

prob.append(prob_s0)
prob.append(prob_s1)

species = []
for s in range(num_species):
    st = Species(num_robots[s], num_states, prob[s])
    species.append(st)
