#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:28:41 2020

@author: leonard
"""
from Constants import Viscosity, Order, AdiabaticIndex, ViscositySoftening, \
    FacIntToCoord
from src.sph.Kernel import kernel
from numpy import heaviside, sign

def force_step(particles):
    "Updates the particles force and rate of entropy change"
    for particle in particles:
        #initialize results
        particle.acceleration = 0
        particle.entropyChange = 0
        for [neighbor, intDistance, vij] in particle.neighbors:
            #we will need these multiple times
            dist = intDistance*FacIntToCoord
            dkern_i = sign(dist)*kernel(dist, particle.smoothingLength, Order, True)
            dkern_j = sign(dist)*kernel(dist, neighbor.smoothingLength, Order, True)
            visc = viscosity_tensor(particle, neighbor, dist, vij)
            
            particle.acceleration -= neighbor.mass * visc * \
                                               (dkern_i + dkern_j)/2
            particle.entropyChange += neighbor.mass * visc * vij * \
                                      (dkern_i + dkern_j)/2
            
            dkern_i *= particle.dhsmlDensityFactor*particle.pressure/particle.density**2
            dkern_j *= neighbor.dhsmlDensityFactor*neighbor.pressure/neighbor.density**2
            particle.acceleration -= neighbor.mass*(dkern_i +dkern_j)
        particle.entropyChange *= (AdiabaticIndex - 1)/2/particle.density**(AdiabaticIndex - 1)

def viscosity_tensor(particle1, particle2, distance, vij):
    "Returns the viscosity tensor for two particles"
    hij = (particle1.smoothingLength + particle2.smoothingLength)/2
    rhoij = (particle1.density + particle2.density)/2
    cij = (particle1.soundspeed + particle2.soundspeed)/2
    muij = hij*vij*distance/(distance**2 + ViscositySoftening*hij**2)
    return heaviside(-distance*vij, 0)* Viscosity*muij*(2*muij - cij)/rhoij

        