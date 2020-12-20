#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:28:41 2020

@author: leonard
"""
import numpy as np
from Constants import Viscosity, AdiabaticIndex, ViscositySoftening, \
    FacIntToCoord, DIM
from src.sph.Kernel import kernel
from numpy.linalg import norm

def force_step(particles):
    "Updates the particles force and rate of entropy change"
    for particle in particles:
        #initialize results
        particle.acceleration = np.zeros(DIM, dtype=float)
        particle.entropyChange = 0
        for [neighbor, intDistance, vij] in particle.neighbors:
            if neighbor.index != particle.index:
                #we will need these multiple times
                dist = intDistance*FacIntToCoord
                dkern_i = dist/norm(dist) * kernel(norm(dist)/particle.smoothingLength, \
                                                   particle.smoothingLength, \
                                                   True)
                dkern_j = dist/norm(dist) * kernel(norm(dist)/neighbor.smoothingLength, \
                                                   neighbor.smoothingLength, \
                                                   True)
                visc = viscosity_tensor(particle, neighbor, dist, vij)
                #add viscosity contribution
                particle.acceleration -= neighbor.mass * visc * \
                                               (dkern_i + dkern_j)/2
                particle.entropyChange += 0.5 * neighbor.mass * visc * \
                                          np.dot(vij, dkern_i + dkern_j)
                #add pressure force contribution
                dkern_i *= particle.dhsmlDensityFactor*\
                           particle.pressure/particle.density**2
                dkern_j *= neighbor.dhsmlDensityFactor*\
                           neighbor.pressure/neighbor.density**2
                particle.acceleration -= neighbor.mass*(dkern_i +dkern_j)
        particle.entropyChange *= (AdiabaticIndex - 1)/2/particle.density**(AdiabaticIndex - 1)
        
def viscosity_tensor(particle1, particle2, distance, vij):
    "Returns the viscosity tensor for two particles"
    if np.dot(distance, vij) < 0:
        hij = (particle1.smoothingLength + particle2.smoothingLength)/2
        rhoij = (particle1.density + particle2.density)/2
        cij = (particle1.soundspeed + particle2.soundspeed)/2
        muij = hij*np.dot(vij, distance)/\
               (norm(distance)**2 + ViscositySoftening*hij**2)
        return Viscosity*muij*(2*muij - cij)/rhoij
    else:
        return 0

        