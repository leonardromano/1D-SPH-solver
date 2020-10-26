#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:52:07 2020

@author: leonard
"""
from src.data.Particle_Class import Particle
import numpy as np
from Constants import NumberOfGridCells, NumberOfParticles, Ntimebins, \
    ICSpecifier, ICfile, Mp
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
        for i in range(NumberOfParticles):
            if i <= 0.8*NumberOfParticles-1:
                dx = 1/(1 + 0.8*NumberOfParticles)
                initialPositions[i] += (1/2 + i)*dx - 1
            else:
                dx = 1/(1+0.2*NumberOfParticles)
                initialPositions[i] += (1/2 + (i-0.8*NumberOfParticles))*dx
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
            particles.append(Particle(line[1], line[2], line[3], line[0]))
    else:
        print("Initial condition file cannot be found.")
        exit()
    return particles
    
    