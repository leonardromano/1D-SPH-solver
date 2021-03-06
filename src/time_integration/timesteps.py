#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 12:22:49 2020

@author: leonard
"""
from math import ceil
from numpy import log2
from sys import exit

from Constants import Dt
from Parameters import NTimebins


def assign_timestep_classes(timeBins, currentLevel):
    """
    For all active particles check the timestep criterion and populate 
    bins with shorter timesteps if necessary
    """
    lowestPopulatedBin = currentLevel
    #create empty instances of higher time bins
    for i in range(currentLevel+1, NTimebins):
        timeBins[i] = list()
    #run over all active particles
    for particle in timeBins[currentLevel]:
        #make sure each particles courant timestep is up to date
        particle.update_timestep_criterion()
        if particle.timestepCriterion < Dt/2**(currentLevel): 
            #need to reduce timestep
            #determine the new timebin
            k = ceil(log2(Dt/particle.timestepCriterion))
            particle.timeBin = k
            if k <= NTimebins: #make sure the new timebin exists
                #need to add particle to all bins with timestep >= necessary timestep
                for i in range(currentLevel+1, k):
                    timeBins[i].append(particle)
                #update lowest populated bin
                if k>lowestPopulatedBin:
                    lowestPopulatedBin = k
            else:
                print("Too small timestep detected! Reduce dt or increase Ntimesteps!")
                exit()
    return lowestPopulatedBin
        
def get_active_time_bin(localTimeStep):
    "determine the current active time bin"
    for level in range(NTimebins):
        if (localTimeStep*2**level)%2 == 1:
            return level