#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:28:41 2020

@author: leonard
"""
from Constants import Viscosity, Order, AdiabaticIndex, ViscositySoftening
from src.sph.Kernel import kernel
from numpy import heaviside

def force_step(particles):
    "Calls the functions for force and entropy change calculations"
    calculate_forces(particles)
    calculate_entropy_change(particles)

def calculate_forces(particles):
    "Computes the forces for all particles by summation over all neighbors"
    for particle in particles:
        particle.acceleration_pressure = 0
        particle.acceleration_viscosity = 0
        for [neighbor,distance, vij] in particle.neighbors:
            particle.acceleration_viscosity -= neighbor.mass * \
                viscosity_tensor(particle, neighbor, distance, vij)* \
                (kernel(distance, particle.smoothingLength, Order, True) + \
                 kernel(-distance, neighbor.smoothingLength, Order, True))/2
            particle.acceleration_pressure -= neighbor.mass*\
                (particle.pressure/particle.density**2 + \
                 neighbor.pressure/neighbor.density**2)*\
                kernel(distance,particle.smoothingLength,Order,True)

def calculate_entropy_change(particles):
    "Computes the rate of entropy change for all particles"
    for particle in particles:
        particle.entropyChange = 0
        for [neighbor,distance, vij] in particle.neighbors:
            particle.entropyChange += (AdiabaticIndex - 1)/2/\
                particle.density**(AdiabaticIndex - 1)*neighbor.mass * \
                viscosity_tensor(particle, neighbor, distance, vij)*vij* \
                (kernel(distance, particle.smoothingLength, Order, True) + \
                 kernel(-distance, neighbor.smoothingLength, Order, True))/2

def viscosity_tensor(particle1, particle2, distance, vij):
    "returns the viscosity tensor for two particles"
    hij = (particle1.smoothingLength + particle2.smoothingLength)/2
    rhoij = (particle1.density + particle2.density)/2
    cij = (particle1.soundspeed + particle2.soundspeed)/2
    muij = hij*vij*distance/(distance**2 + ViscositySoftening*hij**2)
    return heaviside(-distance*vij, 0)* Viscosity*muij*(2*muij - cij)/rhoij

        