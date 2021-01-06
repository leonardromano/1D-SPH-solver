#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:28:41 2020

@author: leonard
"""
from numpy import zeros, dot
from numpy.linalg import norm

from Constants import FacIntToCoord, Mp
from Parameters import Viscosity, AdiabaticIndex, ViscositySoftening, DIM, \
    ExternalForce, GravAcceleration, Floor, GravAxis
from src.data.int_conversion import convert_to_phys_position
from src.sph.Kernel import kernel


def force_step(particles):
    "Updates the particles force and rate of entropy change"
    for particle in particles:
        #initialize results
        particle.acceleration = zeros(DIM, dtype=float)
        particle.entropyChange = 0
        for [neighbor, intDistance, vij] in particle.neighbors:
            if neighbor.index != particle.index:
                #we will need these multiple times
                dist = convert_to_phys_position(intDistance)
                dkern_i = dist/norm(dist) * kernel(norm(dist)/particle.hsml, \
                                                   particle.hsml, True)
                dkern_j = dist/norm(dist) * kernel(norm(dist)/neighbor.hsml, \
                                                   neighbor.hsml, True)
                visc = viscosity_tensor(particle, neighbor, dist, vij)
                #add viscosity contribution
                particle.acceleration -= visc * (dkern_i + dkern_j)/2
                particle.entropyChange += 0.5 * Mp * visc * \
                                          dot(vij, dkern_i + dkern_j)
                #add pressure force contribution
                dkern_i *= particle.dhsmldrho*\
                           particle.pressure/particle.density**2
                dkern_j *= neighbor.dhsmldrho*\
                           neighbor.pressure/neighbor.density**2
                particle.acceleration -= (dkern_i +dkern_j)
        particle.entropyChange *= (AdiabaticIndex - 1)/2/particle.density**(AdiabaticIndex - 1)
        particle.acceleration  *= Mp 
        if ExternalForce:
            if Floor and particle.position[GravAxis] * \
                FacIntToCoord[GravAxis] < particle.hsml:
                   #particle touching the ground
                   continue
            #particle is not touching the ground: apply gravity
            particle.acceleration[GravAxis] -= GravAcceleration
        
def viscosity_tensor(particle1, particle2, rij, vij):
    "Returns the viscosity tensor for two particles"
    if dot(rij, vij) < 0:
        hij = (particle1.hsml + particle2.hsml)/2
        rhoij = (particle1.density + particle2.density)/2
        cij = (particle1.soundspeed + particle2.soundspeed)/2
        muij = hij * dot(vij, rij)/(norm(rij)**2 + ViscositySoftening*hij**2)
        return Viscosity*muij*(2*muij - cij)/rhoij
    else:
        return 0

        