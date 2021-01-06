#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:52:07 2020

@author: leonard
"""
import h5py
import numpy as np
from sys import exit

from Constants import FacIntToCoord, BITS_FOR_POSITIONS, NORM_COEFF, Mp
from Parameters import NPart, NTimebins, BoxSize, ICSpecifier, ICfile, \
    DESNNGBS, Noise_Var_Fac, DIM

from src.data.int_conversion import convert_to_int_position, \
    find_minimum_offset_left, find_stepsize_right
from src.data.Particle_Class import Particle
from src.sph.density import density



def initialize_particles():
    "creates instances of particle objects and writes them in a list"
    if ICSpecifier == "from_file":
        particles = initialize_from_file()
    else:
        initialPositions = get_initial_positions()
        initialVelocities = get_initial_velocities()
        initialEntropies = get_initial_entropies()
        particles = list()
        for i in range(NPart):
            particles.append(Particle(initialPositions[i], initialVelocities[i], \
                                     initialEntropies[i], i))
    return particles

def initialise_time_bins(particles):
    """
    creates instance of list of timebins and fills the first entry 
    with list of particles
    """
    timeBins = dict()
    for i in range(NTimebins):
        timeBins[i] = list()
    timeBins[0] = particles
    return timeBins

def get_initial_positions():
    "returns the initial positions according to the specified initial conditions"
    initialPositions = np.zeros((NPart, DIM), dtype = int)
    if ICSpecifier == "shocktube" and DIM == 1:   #initialize positions for shocktube ICs
        Nl = int(NPart*0.8)
        dx0l = int(find_minimum_offset_left(Nl))
        dxl = int(((1 << (BITS_FOR_POSITIONS-1)) - 2*dx0l)//(Nl-1))
        Nr = int(NPart*0.2)
        dx0r = int(dx0l + Nl*dxl)
        dxr = int(find_stepsize_right(Nr, dx0r))
        for i in range(NPart):
            if i <= 0.8 * NPart-1:
                initialPositions[i] += dx0l + i*dxl
                #perturb with random displacement
                initialPositions[i] += int(np.random.normal(0, dxl * Noise_Var_Fac))
            else:
                initialPositions[i] += dx0r + int(i-0.8*NPart)*dxr
                #perturb with random displacement
                initialPositions[i] += int(np.random.normal(0, dxr * Noise_Var_Fac))
    elif ICSpecifier == "homogeneous" and DIM == 1:
        dx = int(2**BITS_FOR_POSITIONS // (NPart+1))
        for i in range(NPart):
            initialPositions[i] += (i+1)*dx
            #perturb with random displacement
            initialPositions[i] += int(np.random.normal(0, dx * Noise_Var_Fac))      
    else:   #if the wanted ICs aren't defined stop
        print("The specified initial conditions cannot be found.")
        exit()
    return initialPositions

def get_initial_velocities():
    "returns the initial velocities"
    initialVelocities = np.zeros((NPart, DIM))
    if ICSpecifier == "shocktube" and DIM == 1:
        initialVelocities += np.random.normal(0, Noise_Var_Fac, (NPart, DIM))
    elif ICSpecifier == "homogeneous" and DIM == 1:
        initialVelocities += np.random.normal(0, Noise_Var_Fac, (NPart, DIM))
    else:
        print("The specified initial conditions cannot be found.")
        exit()
    return initialVelocities

def get_initial_entropies():
    "returns the initial entropy functions"
    initialEntropies = np.zeros(NPart)
    if ICSpecifier == "shocktube":
        for i in range(NPart):
            if i<=0.8*NPart-1:
                initialEntropies[i] += 1.352
            else:
                initialEntropies[i] += 1.656
    elif ICSpecifier == "homogeneous":
        for i in range(NPart):
            initialEntropies[i] += 1.
    else:
        print("The specified initial conditions cannot be found.")
        exit()
    return initialEntropies

def initialize_from_file():
    "initializes the particles from data provided by a file"
    file = h5py.File(ICfile, "r")
    header = file["Header"].attrs
    if (NPart != header["NumPart"]):
        print("NumberOfParticles = %d, expected: %d" \
              %(header["NumPart"], NPart))
        exit()
    positions  = np.asarray(file["PartData/Coordinates"])
    velocities = np.asarray(file["PartData/Velocity"]) 
    entropies  = np.asarray(file["PartData/Entropy"])
    particles = list()
    for i in range(header["NumPart"]):
        particles.append(Particle(convert_to_int_position(positions[i]), \
                                  velocities[i], entropies[i], i))
    file.close()
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
        while(2 * DESNNGBS * Mp > NgbTree.get_nodep(no).Mass):
            p = NgbTree.get_nodep(no).father
            if p < 0:
                break
            no = p
        length = 0
        if(NgbTree.get_nodep(no).level > 0):
            length = (1 << (BITS_FOR_POSITIONS - NgbTree.get_nodep(no).level)) * \
                     FacIntToCoord.max()
        else:
            length = BoxSize.max()

        particles[i].smoothingLength =  length * \
            (DESNNGBS * Mp / NgbTree.get_nodep(no).Mass / NORM_COEFF)**(1/DIM)