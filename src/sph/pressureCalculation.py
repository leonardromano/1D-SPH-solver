#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:15:09 2020

@author: leonard
"""

from src.sph.neighbors import find_neighbors
from src.init.Initialization import initialize_grid
from src.sph.smoothingLengths import set_smoothing_lengths

from Constants import Order
from src.sph.Kernel import kernel

def calculate_density(particles, ahead = False):
    """
    for each particle adds up each neighbors contribution to the density 
    and then updates pressure and soundspeed
    """
    for particle in particles:
        particle.density = 0
        for [neighbor, distance, vij] in particle.neighbors:
            particle.density += neighbor.mass*\
                kernel(distance, particle.smoothingLength, Order)
        particle.update_pressure(ahead)
        particle.update_soundspeed()

def compute_sph_quantities(particles, ahead = False):
    "Fills the grid with particles and uses it to compute sph quantities"
    grid = initialize_grid(particles)
    set_smoothing_lengths(grid, particles)
    find_neighbors(grid, particles, ahead)
    calculate_density(particles, ahead)