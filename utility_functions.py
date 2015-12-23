# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 11:35:51 2015
@author: amandaprorok

"""

import numpy as np
import matplotlib.pyplot as plt
import itertools

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


def subsets(iterable):
    "Generate the subsets of elements in the iterable, in order by length."
    items = list(iterable)
    for k in xrange(len(items) + 1):
        for subset in itertools.combinations(items, k):
            yield subset