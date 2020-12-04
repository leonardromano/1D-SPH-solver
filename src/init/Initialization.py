#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:52:07 2020

@author: leonard
"""
from src.data.Particle_Class import Particle
from src.data.int_conversion import convert_to_int_position, \
    find_minimum_offset_left, find_stepsize_right
import numpy as np
from Constants import NumberOfGridCells, NumberOfParticles, Ntimebins, \
    ICSpecifier, ICfile, K, BinWidth, FacIntToCoord
from src.sph.density import density
from src.sph.neighbors import find_neighbors
from sys import exit


def initialize_particles():
    "creates instances of particle objects and writes them in a list"
    if ICSpecifier == "from_file":
        particles = initialize_from_file()
    else:
        initialPositions = get_initial_positions()
        initialVelocities = get_initial_velocities()
        initialEntropies = get_initial_entropies()
        particles = list()
        for i in range(NumberOfParticles):
            particles.append(Particle(initialPositions[i], initialVelocities[i], \
                                     initialEntropies[i], i))
    return particles

def initialise_time_bins(particles):
    """
    creates instance of list of timebins and fills the first entry 
    with list of particles
    """
    timeTree = dict()
    for i in range(Ntimebins):
        timeTree[i] = list()
    timeTree[0] = particles
    return timeTree
        
def initialize_grid(particles):
    "creates instance of grid which subdivides space into bins which contain particles"
    grid = list()
    for i in range(NumberOfGridCells):
        cell = list()
        for particle in particles:
            if particle.currentBin == i:
                cell.append(particle)
        grid.append(cell)
    return grid

def get_initial_positions():
    "returns the initial positions according to the specified initial conditions"
    initialPositions = np.zeros(NumberOfParticles)
    if ICSpecifier == "shocktube":   #initialize positions for shocktube ICs
        Nl = int(NumberOfParticles*0.8)
        dx0l = int(find_minimum_offset_left(Nl))
        dxl = int((2**31 - 2*dx0l)//(Nl-1))
        Nr = int(NumberOfParticles*0.2)
        dx0r = int(dx0l + 0.8*NumberOfParticles*dxl)
        dxr = int(find_stepsize_right(Nr, dx0r))
        for i in range(NumberOfParticles):
            if i <= 0.8*NumberOfParticles-1:
                initialPositions[i] += dx0l + i*dxl
            else:
                initialPositions[i] += dx0r + int(i-0.8*NumberOfParticles)*dxr
    else:   #if the wanted ICs aren't defined stop
        print("The specified initial conditions cannot be found.")
        exit()
    return initialPositions

def get_initial_velocities():
    "returns the initial velocities"
    initialVelocities = np.zeros(NumberOfParticles)
    if ICSpecifier == "shocktube":
        return initialVelocities
    else:
        print("The specified initial conditions cannot be found.")
        exit()

def get_initial_entropies():
    "returns the initial entropy functions"
    initialEntropies = np.zeros(NumberOfParticles)
    if ICSpecifier == "shocktube":
        for i in range(NumberOfParticles):
            if i<=0.8*NumberOfParticles-1:
                initialEntropies[i] += 1.352
            else:
                initialEntropies[i] += 1.656
    else:
        print("The specified initial conditions cannot be found.")
        exit()
    return initialEntropies

def initialize_from_file():
    "initializes the particles from data provided by a file"
    file = np.loadtxt(ICfile)
    if (NumberOfParticles != file.shape[0]):
        print("NumberOfParticles = %d, expected: %d" \
              %(file.shape[0], NumberOfParticles))
        exit()
    particles = list()
    if file.size:
        for line in file:
            particles.append(Particle(convert_to_int_position(line[1]), \
                                      line[2], line[3], line[0]))
    else:
        print("Initial condition file cannot be found.")
        exit()
    return particles
    
def update_sph_quantities(particles, ahead = False):
    "This function initializes the sph properties of the particles"
    grid = initialize_grid(particles)
    #ahead=False means we need to provide an initial guess
    if not ahead:
        initial_guess_hsml(grid, particles)
    #compute density and hsml and set thermodynamic quantities
    density(grid, particles, ahead)
    #Now we just need to reassign the neighbors with correct densities
    find_neighbors(grid, particles, ahead)
    
def initial_guess_hsml(grid, particles):
    "computes an initial guess for the smoothing lengths"
    for particle in particles:
        Nest = len(grid[particle.currentBin])
        l=0
        while(3*K>Nest):
            l += 1
            if (l<=particle.currentBin<NumberOfGridCells-l):
                #add number of particles in both neighboring bins
                Nest += len(grid[particle.currentBin - l])
                Nest += len(grid[particle.currentBin + l])
            elif  (l<=particle.currentBin):
                #add number of particles in the next lower bin twice
                Nest += 2*len(grid[particle.currentBin - l])
            elif (particle.currentBin<NumberOfGridCells-l):
                #add number of particles in the next higher bin twice
                Nest += 2*len(grid[particle.currentBin + l])
            else:
                #already enclosing all bins something wrent wong!
                print("Want to add more Neighbors than there are particles! Abort!")
                exit()
        particle.smoothingLength = K/Nest*(1+2*l)*BinWidth*FacIntToCoord
    