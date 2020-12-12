#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:25:36 2020

@author: leonard
"""

from src.time_integration.Drift import Kick, Drift
from src.forces.forceCalculation import force_step
import src.init.Initialization as init
from src.sph.density import density
from Constants import NumberOfTimesteps, Dt, StepsBetweenSnapshots, FinalTime
from src.writing.writing import write_data, write_header
from src.time_integration.timesteps import assign_timestep_classes, \
    get_active_time_bin
from src.tree.tree import ngbtree

def main():
    "This function initializes the simulation and contains the main loop"
    write_header()
    #initialize particles and perform first force calculation
    particles = init.initialize_particles()
    timeBins = init.initialise_time_bins(particles)
    NgbTree = ngbtree(particles)
    init.sph_quantities(particles, NgbTree)
    force_step(timeBins[0])
    #assign particles to time bins
    lowestPopulatedTimeBin = assign_timestep_classes(timeBins, 0)
    #loop over all timesteps
    for i in range(NumberOfTimesteps):
        #Write a snapshot if required
        if i%StepsBetweenSnapshots == 0:
            print("writing Snapshot no. %d" %(i//StepsBetweenSnapshots))
            write_data(particles, i//StepsBetweenSnapshots, i*Dt)
        #j is the local 'clock'
        j=0
        activeTimeBin = 0
        #cycle integration step over all timebins
        while j<1:
            #first kick all active particles by the current active timestep
            Kick(timeBins[activeTimeBin])
            #then synchronize all particles with the next synchronization point
            Drift(timeBins[0], Dt/2**lowestPopulatedTimeBin)
            
            #update the local clock to the next synchronization point and 
            #determine which particles are active next
            j += 2**(-lowestPopulatedTimeBin)
            activeTimeBin = get_active_time_bin(j)
            
            #at the new synchronisation point update sph quantities and forces
            NgbTree = ngbtree(timeBins[0])
            density(timeBins[0], NgbTree, True)
            force_step(timeBins[activeTimeBin])
            
            #kick the active particles by the current active timestep
            Kick(timeBins[activeTimeBin])
            
            #determine the new smallest timestep
            lowestPopulatedTimeBin = assign_timestep_classes(timeBins, activeTimeBin)
            
            #repeat the previous steps
    #after we're done we want to write the final results in a snapshot
    write_data(particles, NumberOfTimesteps//StepsBetweenSnapshots, FinalTime)