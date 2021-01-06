#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:06:22 2020

@author: leonard
"""
from Constants import Dt, BITS_FOR_POSITIONS, MAX_INT
from Parameters import DIM, Periodic

from src.data.int_conversion import convert_to_int_positions as phys2int

def Drift(particles, timestep):
    "drifts all particle's positions to the next synchronisation point"
    for particle in particles:
        particle.position += phys2int(particle.velocity * timestep)
        keep_inside_box(particle)

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
    

def keep_inside_box(particle):
    "Makes sure the particle is within the domain"
    for axis in range(DIM):
        if Periodic[axis]:
            while particle.position[axis] < 0:
                particle.position[axis] += (1 << BITS_FOR_POSITIONS)
            while particle.position[axis] > MAX_INT:
                particle.position[axis] -= (1 << BITS_FOR_POSITIONS)
        else:
            while particle.position[axis] < 0 or particle.position[axis] > MAX_INT:
                if particle.position[axis] < 0:
                    particle.position[axis] *= -1
                    particle.velocity[axis] *= -1
                else:
                    particle.position[axis] += 2 * (MAX_INT - particle.position[axis])
                    particle.velocity[axis] *= -1
