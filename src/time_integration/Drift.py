#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:06:22 2020

@author: leonard
"""
from Constants import LeftBoundary, RightBoundary, Dt

def Drift(particles, timestep):
    "drifts all particle's positions to the next synchronisation point"
    for particle in particles:
        particle.position += particle.velocity*timestep
        if particle.position <LeftBoundary:
            particle.position += 2*(LeftBoundary - particle.position)
            particle.velocity -= 2*particle.velocity
        elif particle.position >RightBoundary:
            particle.position +=2*(RightBoundary - particle.position)
            particle.velocity -= 2*particle.velocity
        particle.update_current_bin()

def Kick(particles):
    "drifts all particles velocities and entropy functions to the next sync point"
    for particle in particles:
        particle.velocity_ahead = particle.velocity + \
            (particle.acceleration_pressure + particle.acceleration_viscosity)*\
            Dt/2**(particle.timeBin)
        particle.velocity += (particle.acceleration_pressure + \
                              particle.acceleration_viscosity)*\
                             Dt/2**(1+particle.timeBin)
        particle.entropy_ahead = particle.entropy + particle.entropyChange*\
            Dt/2**(particle.timeBin)
        particle.entropy += particle.entropyChange*Dt/2**(1+particle.timeBin)
    

