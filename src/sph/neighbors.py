#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:56:47 2020

@author: leonard
"""
from Constants import NumberOfGridCells, BinWidth, K, Mp
from src.data.Particle_Class import Particle
from src.sph.smoothingLengths import get_distance_from_wall

def find_neighbors(grid, particles, ahead = False):
    "For each particle creates a list of neighboring particles"
    clear_neighbors(particles)
    #Run over all particles
    for i in range(len(particles)):
        l = 1 + int((particles[i].smoothingLength)//BinWidth)
        #Run over all bins
        for j in range(max(particles[i].currentBin-l, 0), \
                       min(particles[i].currentBin+1+l, NumberOfGridCells)):
            #In each bin run over all particles
            for citizen in grid[j]:
                #In order to save some effort when a particle is another
                #particles neighbor and they are both within each others
                #smoothing radius add each other to their list of neighbors.
                #Naturally this may happen only once.
                if citizen.index>=i:
                    distance = particles[i].position - citizen.position
                    if abs(distance) < particles[i].smoothingLength:
                        if ahead == False:
                            vij = particles[i].velocity - citizen.velocity
                        else:
                             vij = particles[i].velocity_ahead - \
                                 citizen.velocity_ahead
                        particles[i].neighbors.append([citizen, distance, vij])
                        if abs(distance)<citizen.smoothingLength and \
                            citizen.index !=i:
                            citizen.neighbors.append([particles[i], -distance, -vij])
                else:
                    #A particle with large smoothing length can have a neighbor
                    #with short smoothing length who does not have the particle
                    #as a neighbor
                    if citizen.smoothingLength<=particles[i].smoothingLength:
                        distance = particles[i].position - citizen.position
                        if citizen.smoothingLength < abs(distance) < \
                            particles[i].smoothingLength:
                            if ahead == False:
                                vij = particles[i].velocity - citizen.velocity
                            else:
                                vij = particles[i].velocity_ahead - \
                                    citizen.velocity_ahead
                            particles[i].neighbors.append([citizen, distance, vij])
        #Now add all ghost particles to list of neighbors of particles that
        #are close to the boundaries
        distanceFromWall = get_distance_from_wall(particles[i])
        if distanceFromWall<particles[i].smoothingLength:
            minimumDistanceFromWall = \
                get_minimum_distance_from_wall(particles[i].neighbors)
            listOfNeighbors = copy_list(particles[i].neighbors)
            for [neighbor, distance, vij] in listOfNeighbors:
                neighborDistanceFromWall = get_distance_from_wall(neighbor)
                if neighborDistanceFromWall>distanceFromWall and \
                    distanceFromWall<abs(distance)+minimumDistanceFromWall:
                    particles[i].neighbors.\
                        append(mirror_particle(particles[i], neighbor, ahead))
            
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
    mirrorParticle = Particle(2*particle.position - neighbor.position, \
                              velocity, neighbor.entropy, neighbor.index)
    mirrorParticle.density = K*Mp/2*neighbor.smoothingLength
    mirrorParticle.update_pressure(ahead)
    mirrorParticle.update_soundspeed()
    return [mirrorParticle, particle.position - mirrorParticle.position, vij]

def get_minimum_distance_from_wall(particles):
    "Determines the minimum distance from the wall among a list of particles"
    minimumDistanceFromWall = get_distance_from_wall(particles[0][0])
    for [particle, distance, vij] in particles:
        if get_distance_from_wall(particle)<minimumDistanceFromWall:
            minimumDistanceFromWall = get_distance_from_wall(particle)
    return minimumDistanceFromWall

def copy_list(oldList):
    "creates a hard copy of a list"
    newList = list()
    for item in oldList:
        newList.append(item)
    return newList