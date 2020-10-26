#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 12:22:49 2020

@author: leonard
"""
from Constants import Dt, Ntimebins
from numpy import log2
from sys import exit

def get_min_timestep(particles):
    minTimestep = Dt
    for particle in particles:
        particle.updateCourantTimestep()
        if particle.courantTimestep<minTimestep:
            minTimestep = particle.courantTimestep
    return minTimestep

def assign_timestep_classes(timeBins, currentLevel):
    """
    For all active particles check the timestep criterion and populate 
    bins with shorter timesteps if necessary
    """
    lowestPopulatedBin = currentLevel
    #create empty instances of higher time bins
    for i in range(currentLevel+1, Ntimebins):
        timeBins[i] = list()
    #run over all active particles
    for particle in timeBins[currentLevel]:
        #make sure each particles courant timestep is up to date
        particle.update_courant_timestep()
        if particle.courantTimestep<Dt/2**(currentLevel): #need to reduce timestep
            k = int(log2(Dt/particle.courantTimestep)) #determine the new timebin
            if k<=Ntimebins: #make sure the new timebin exists
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
    for level in range(Ntimebins):
        if (localTimeStep*2**level)%2 == 1:
            return level