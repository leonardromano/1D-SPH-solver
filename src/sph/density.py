#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:19:30 2020

@author: leonard
"""
from Constants import K, DK, Mp, Order, FacIntToCoord
from src.sph.neighbors import find_neighbors_density
from src.sph.Kernel import kernel 
from sys import exit

def density(grid, particles, ahead = False):
    "For each particle compute density, smoothing length and thermodynamic variables"
    for particle in particles:
        numberOfNeighbors = 0
        lowerBound = 0
        upperBound = 0
        while(abs(numberOfNeighbors-K)>DK):
            #clear all density data
            particle.density = 0
            particle.dhsmlDensityFactor = 0
            particle.neighbors = list()
            #find list of neighbors used for density calculation
            find_neighbors_density(grid, particle)
            #compute density and dh/drho
            evaluate_kernel(particle)
            #now update numberOfNeighbors
            numberOfNeighbors = 2*particle.smoothingLength*particle.density/Mp
            if (abs(numberOfNeighbors-K)>DK):
                #check whether we're done
                if lowerBound>0 and upperBound >0:
                    if upperBound-lowerBound < 10e-3 * lowerBound:
                        continue
                lowerBound, upperBound = update_bounds(lowerBound, upperBound, \
                                                       numberOfNeighbors, \
                                                       particle.smoothingLength)
                update_smoothing_length(lowerBound, upperBound, \
                                        numberOfNeighbors, particle)
        #after the density has been computed we can now set pressure and soundspeed
        set_thermodynamic_variables(particle, ahead)
        

def update_bounds(lowerBound, upperBound, numberOfNeighbors, h):
    "Update the bounds for the smoothing length in the bisection algorithm"
    if numberOfNeighbors < K-DK:
        lowerBound = max(lowerBound, h)
    else:
        if upperBound != 0:
            if h < upperBound:
                upperBound = h
        else:
            upperBound = h
    return lowerBound, upperBound

def update_smoothing_length(lowerBound, upperBound, numberOfNeighbors, particle):
    "Perform the bisection part of the bisection algorithm"
    if lowerBound >0 and upperBound >0:
        particle.smoothingLength = ((lowerBound**3 + upperBound**3)/2)**(1/3)
    else:
        if upperBound == 0 and lowerBound == 0:
            print("Upper and Lower bounds not updated!")
            exit()
            
        if upperBound == 0 and lowerBound>0:
            if abs(numberOfNeighbors - K) < 0.5*K:
                fac = 1-(numberOfNeighbors - K)/numberOfNeighbors*\
                    particle.dhsmlDensityFactor
                if fac<1.26:
                    particle.smoothingLength *=fac
                else:
                    particle.smoothingLength *=1.26
            else:
                particle.smoothingLength *=1.26
                    
        if upperBound > 0 and lowerBound  == 0:
            if abs(numberOfNeighbors - K) < 0.5*K:
                fac = 1-(numberOfNeighbors - K)/numberOfNeighbors*\
                    particle.dhsmlDensityFactor
                if fac > 1/1.26:
                    particle.smoothingLength *=fac
                else:
                    particle.smoothingLength /=1.26
            else:
                particle.smoothingLength /=1.26      

def set_thermodynamic_variables(particle, ahead = False):
    "Update the particle's pressure and sound speed"
    particle.update_pressure(ahead)
    particle.update_soundspeed()
        
def evaluate_kernel(particle):
    "Perform the neighbor sum to compute density and SPH correction factor"
    for [neighbor, intDistance] in particle.neighbors:
        #precompute the kernel and its derivative
        dist = abs(intDistance*FacIntToCoord)
        kern = kernel(dist, particle.smoothingLength, Order)
        dkern = kernel(dist, particle.smoothingLength, Order, True)
        #add neighbors contribution to density and dhsml
        particle.density += neighbor.mass*kern
        particle.dhsmlDensityFactor -= neighbor.mass*dist*dkern
    #now if we have more than one neighbor update the dhsml factor
    if particle.dhsmlDensityFactor > 0:
        particle.dhsmlDensityFactor = particle.density/particle.dhsmlDensityFactor
        if particle.dhsmlDensityFactor > 10:
            print(print("Number of neighbors = %d, index = %d"\
              %(len(particle.neighbors), particle.index)))
    else:
        print("Number of neighbors = %d, index = %d"\
              %(len(particle.neighbors), particle.index))
        particle.dhsmlDensityFactor = 1