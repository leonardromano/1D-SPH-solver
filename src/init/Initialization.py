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
from Constants import NumberOfParticles, Ntimebins, TotalSideLength, \
    ICSpecifier, ICfile, K, FacIntToCoord, BITS_FOR_POSITIONS, Noise_Var_Fac
from src.sph.density import density
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
    timeBins = dict()
    for i in range(Ntimebins):
        timeBins[i] = list()
    timeBins[0] = particles
    return timeBins

def get_initial_positions():
    "returns the initial positions according to the specified initial conditions"
    initialPositions = np.zeros(NumberOfParticles, dtype = int)
    if ICSpecifier == "shocktube":   #initialize positions for shocktube ICs
        Nl = int(NumberOfParticles*0.8)
        dx0l = int(find_minimum_offset_left(Nl))
        dxl = int((2**(BITS_FOR_POSITIONS-1) - 2*dx0l)//(Nl-1))
        Nr = int(NumberOfParticles*0.2)
        dx0r = int(dx0l + Nl*dxl)
        dxr = int(find_stepsize_right(Nr, dx0r))
        for i in range(NumberOfParticles):
            if i <= 0.8*NumberOfParticles-1:
                initialPositions[i] += dx0l + i*dxl
                #perturb with random displacement
                initialPositions[i] += int(np.random.normal(0, dxl * Noise_Var_Fac))
            else:
                initialPositions[i] += dx0r + int(i-0.8*NumberOfParticles)*dxr
                #perturb with random displacement
                initialPositions[i] += int(np.random.normal(0, dxr * Noise_Var_Fac))
    elif ICSpecifier == "homogeneous":
        dx = int(2**BITS_FOR_POSITIONS // (NumberOfParticles+1))
        for i in range(NumberOfParticles):
            initialPositions[i] += (i+1)*dx
            #perturb with random displacement
            initialPositions[i] += int(np.random.normal(0, dx * Noise_Var_Fac))
    else:   #if the wanted ICs aren't defined stop
        print("The specified initial conditions cannot be found.")
        exit()
    return initialPositions

def get_initial_velocities():
    "returns the initial velocities"
    initialVelocities = np.zeros(NumberOfParticles)
    if ICSpecifier == "shocktube":
        initialVelocities += np.random.normal(0, Noise_Var_Fac, NumberOfParticles)
    elif ICSpecifier == "homogeneous":
        initialVelocities += np.random.normal(0, Noise_Var_Fac, NumberOfParticles)
    else:
        print("The specified initial conditions cannot be found.")
        exit()
    return initialVelocities

def get_initial_entropies():
    "returns the initial entropy functions"
    initialEntropies = np.zeros(NumberOfParticles)
    if ICSpecifier == "shocktube":
        for i in range(NumberOfParticles):
            if i<=0.8*NumberOfParticles-1:
                initialEntropies[i] += 1.352
            else:
                initialEntropies[i] += 1.656
    elif ICSpecifier == "homogeneous":
        for i in range(NumberOfParticles):
            initialEntropies[i] += 1.
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
    
def sph_quantities(particles, NgbTree):
    "This function initializes the sph properties of the particles"
    initial_guess_hsml(particles, NgbTree)
    #compute density, smoothing length and thermodynamic quantities and 
    #finds neighbors
    density(particles, NgbTree)
    
def initial_guess_hsml(particles, NgbTree):
    "computes an initial guess for the smoothing lengths"
    for i in range(NgbTree.MaxPart):
        no = NgbTree.Father[i]
        while(2 * K * particles[i].mass > NgbTree.get_nodep(no).Mass):
            p = NgbTree.get_nodep(no).father
            if p<0:
                break
            no = p
        length = 0
        if(NgbTree.get_nodep(no).level > 0):
            length = 2**(BITS_FOR_POSITIONS - NgbTree.get_nodep(no).level) * FacIntToCoord
        else:
            length = TotalSideLength

        particles[i].smoothingLength = K * particles[i].mass * length / \
                                       NgbTree.get_nodep(no).Mass /2
    