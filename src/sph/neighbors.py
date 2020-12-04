#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:56:47 2020

@author: leonard
"""
from math import ceil
from Constants import NumberOfGridCells, BinWidth, FacIntToCoord
from src.data.Particle_Class import Particle

def find_neighbors(grid, particles, ahead = False):
    "For each particle creates a list of neighboring particles"
    #first need to empty all information about neighbors
    clear_neighbors(particles)
    #Run over all particles
    for i in range(len(particles)):
        intHsml = ceil(particles[i].smoothingLength/FacIntToCoord)
        l = ceil(intHsml/BinWidth)
        #Run over all bins
        for j in range(max(particles[i].currentBin-l, 0), \
                       min(particles[i].currentBin+1+l, NumberOfGridCells)):
            #In each bin run over all particles
            for citizen in grid[j]:
                #In order to save some effort when a particle is another
                #particles neighbor and they are both within each others
                #smoothing radius add each other to their list of neighbors.
                #Naturally this may happen only once.
                if citizen.index > i:
                    intDistance = particles[i].position - citizen.position
                    if abs(intDistance) < intHsml:
                        if ahead == False:
                            vij = particles[i].velocity - citizen.velocity
                        else:
                             vij = particles[i].velocity_ahead - \
                                 citizen.velocity_ahead
                        particles[i].neighbors.append([citizen, intDistance, vij])
                        if abs(intDistance) < ceil(citizen.smoothingLength/FacIntToCoord):
                            citizen.neighbors.append([particles[i], -intDistance, -vij])
                elif citizen.index < i:
                    #A particle with large smoothing length can have a neighbor
                    #with short smoothing length who does not have the particle
                    #as a neighbor
                    if citizen.smoothingLength < particles[i].smoothingLength:
                        intDistance = particles[i].position - citizen.position
                        if int(citizen.smoothingLength//FacIntToCoord) < \
                            abs(intDistance) < intHsml:
                            if ahead == False:
                                vij = particles[i].velocity - citizen.velocity
                            else:
                                vij = particles[i].velocity_ahead - \
                                    citizen.velocity_ahead
                            particles[i].neighbors.append([citizen, intDistance, vij])
        #Now add all ghost particles to list of neighbors of particles that
        #are close to the boundaries
        intDistanceFromWall = get_distance_from_wall(particles[i])
        if intDistanceFromWall<intHsml:
            minimumIntDistanceFromWall = \
                get_minimum_distance_from_wall(particles[i].neighbors)
            listOfNeighbors = copy_list(particles[i].neighbors)
            for [neighbor, intDistance, vij] in listOfNeighbors:
                neighborIntDistanceFromWall = get_distance_from_wall(neighbor)
                if neighborIntDistanceFromWall > intDistanceFromWall and \
                    intDistanceFromWall < abs(intDistance) + minimumIntDistanceFromWall:
                    particles[i].neighbors.\
                        append(mirror_particle(particles[i], neighbor, ahead))

def find_neighbors_density(grid, particle):
    "Creates a list of neighbors of particle"
    intHsml = ceil(particle.smoothingLength/FacIntToCoord)
    l = ceil(intHsml/BinWidth)
    #Run over all near bins
    for j in range(max(particle.currentBin-l, 0), \
                   min(particle.currentBin+1+l, NumberOfGridCells)):
        #In each bin run over all particles
        for citizen in grid[j]:
            intDistance = particle.position - citizen.position
            if abs(intDistance) < intHsml:
                particle.neighbors.append([citizen, intDistance])
    #Now add all ghost particles to list of neighbors if the particle is
    #close to the boundary
    intDistanceFromWall = get_distance_from_wall(particle)
    if intDistanceFromWall < intHsml:
       minimumIntDistanceFromWall = get_minimum_distance_from_wall(particle.neighbors)
       listOfNeighbors = copy_list(particle.neighbors)
       for [neighbor, intDistance] in listOfNeighbors:
           neighborIntDistanceFromWall = get_distance_from_wall(neighbor)
           if neighborIntDistanceFromWall > intDistanceFromWall and \
              intDistanceFromWall < abs(intDistance) + minimumIntDistanceFromWall:
                  particle.neighbors.append(mirror_particle_density(particle, neighbor))            

def clear_neighbors(particles):
    "Resets the list of neighbors for all particles"
    for particle in particles:
        particle.neighbors = list()
        
def mirror_particle(particle, neighbor, ahead):
    """
    creates an instance of a ghost particle mirroring a neighbor 
    of a particle close to the wall
    """
    if ahead:
        velocity = -neighbor.velocity_ahead
        vij = particle.velocity_ahead - velocity
    else:
        velocity = -neighbor.velocity
        vij = particle.velocity - velocity
    mirrorParticle = Particle(int(2*particle.position - neighbor.position), \
                              velocity, neighbor.entropy, neighbor.index)
    mirrorParticle.entropy = neighbor.entropy
    mirrorParticle.entropy_ahead = neighbor.entropy_ahead
    mirrorParticle.smoothingLength = neighbor.smoothingLength
    mirrorParticle.density = neighbor.density
    mirrorParticle.dhsmlDensityFactor = neighbor.dhsmlDensityFactor
    mirrorParticle.update_pressure(ahead)
    mirrorParticle.update_soundspeed()
    return [mirrorParticle, particle.position - mirrorParticle.position, vij]

def mirror_particle_density(particle, neighbor):
    """
    creates an instance of a ghost particle mirroring a neighbor 
    of a particle close to the wall
    """
    mirrorParticle = Particle(int(2*particle.position - neighbor.position), \
                              0, neighbor.entropy, neighbor.index)
    return [mirrorParticle, particle.position - mirrorParticle.position]

def get_distance_from_wall(particle):
    "returns the distance from the nearest boundary"
    return int(min(particle.position, 2**32 - particle.position))

def get_minimum_distance_from_wall(particles):
    "Determines the minimum distance from the wall among a list of particles"
    minimumDistanceFromWall = 2**32
    for particle in particles:
        dist = get_distance_from_wall(particle[0])
        if dist<minimumDistanceFromWall:
            minimumDistanceFromWall = dist
    return minimumDistanceFromWall

def copy_list(oldList):
    "creates a hard copy of a list"
    newList = list()
    for item in oldList:
        newList.append(item)
    return newList