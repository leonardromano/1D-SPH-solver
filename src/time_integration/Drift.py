#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:06:22 2020

@author: leonard
"""
from Constants import FacIntToCoord, Dt, BITS_FOR_POSITIONS

def Drift(particles, timestep):
    "drifts all particle's positions to the next synchronisation point"
    for particle in particles:
        particle.position += int((particle.velocity*timestep)//FacIntToCoord)
        if particle.position < 0:
            particle.position -= 2*particle.position
            particle.velocity -= 2*particle.velocity
        elif particle.position > int(2**BITS_FOR_POSITIONS):
            particle.position += int(2*(2**BITS_FOR_POSITIONS - particle.position))
            particle.velocity -= 2*particle.velocity

def Kick(particles):
    "drifts all particles velocities and entropy functions to the next sync point"
    for particle in particles:
        tb = particle.timeBin
        #update velocity
        particle.velocity_ahead = particle.velocity + particle.acceleration*\
                                  Dt/2**tb
        particle.velocity += particle.acceleration*Dt/2**(1+tb)
        #update entropy
        particle.entropy_ahead = particle.entropy + particle.entropyChange*\
            Dt/2**tb
        particle.entropy += particle.entropyChange*Dt/2**(1+tb)
    

