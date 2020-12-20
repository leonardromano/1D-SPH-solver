#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:06:22 2020

@author: leonard
"""
from Constants import FacIntToCoord, Dt, BITS_FOR_POSITIONS, DIM, BoundaryPeriodic

def Drift(particles, timestep):
    "drifts all particle's positions to the next synchronisation point"
    for particle in particles:
        for i in range(DIM):
            particle.position[i] += int((particle.velocity[i]*timestep)//FacIntToCoord)
            #Now enforce boundary conditions
            if BoundaryPeriodic[i]:
                #periodic boundaries
                if particle.position[i] < 0:
                   particle.position[i] += (1 << BITS_FOR_POSITIONS)
                elif particle.position[i] > (1 << BITS_FOR_POSITIONS):
                    particle.position[i] -= (1 << BITS_FOR_POSITIONS)
            else:
                #reflective boundaries
                if particle.position[i] < 0:
                   particle.position[i] -= 2*particle.position[i]
                   particle.velocity[i] -= 2*particle.velocity[i]
                elif particle.position[i] > (1 << BITS_FOR_POSITIONS):
                    particle.position[i] += int(2*((1 << BITS_FOR_POSITIONS) \
                                                   - particle.position[i]))
                    particle.velocity[i] -= 2*particle.velocity[i]

def Kick(particles):
    "drifts all particles velocities and entropy functions to the next sync point"
    for particle in particles:
        tb = particle.timeBin
        #update velocity
        particle.velocity_ahead = particle.velocity + particle.acceleration*\
                                  Dt/(1 << tb)
        particle.velocity += particle.acceleration*Dt/(1 << (1+tb))
        #update entropy
        particle.entropy_ahead = particle.entropy + particle.entropyChange*\
            Dt/(1 << tb)
        particle.entropy += particle.entropyChange*Dt/(1 << (1+tb))
    

