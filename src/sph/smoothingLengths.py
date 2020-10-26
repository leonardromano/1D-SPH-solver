#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 19:13:02 2020

@author: leonard
"""
from Constants import Dh, K, Order, NumberOfGridCells, Hmax, BinWidth, \
    LeftBoundary, RightBoundary
from src.sph.Kernel import kernel 

def set_smoothing_lengths(grid, particles):
    "For each particle find the smoothing length by minimizing the sph functional"
    for particle in particles:
        particle.smoothingLength = minimize_sph_functional(grid, particle)

def minimize_sph_functional(grid, particle):
    "Iteratively minimizes the sph functional"
    i=1
    if sph_functional(grid, particle, Dh)<=0:
        return Dh
    else:
        res=1
        while res>0 and i*Dh<Hmax:
            res = sph_functional(grid, particle, i*Dh)
            i+=1
        if res == 0:
            return i*Dh
        else:
            return (i-1)*Dh

def get_distance_from_wall(particle):
    "returns the distance from the nearest boundary"
    return min(abs(particle.position - LeftBoundary), \
               abs(particle.position - RightBoundary))

def sph_functional(grid, particle, hsml):
    "Returns the value of the sph functional including contributions from ghosts"
    y = K/2
    l = 1 + int(hsml//BinWidth)
    neighbors = list()
    distanceFromWall = get_distance_from_wall(particle)
    minimumDistanceFromWall = distanceFromWall
    for j in range(max(particle.currentBin-l, 0), \
                       min(particle.currentBin+l+1, NumberOfGridCells)):
            for citizen in grid[j]:
                distance = abs(particle.position - citizen.position)
                if distance <= hsml:
                    neighbors.append([citizen, distance])
                    y -= hsml*kernel(distance, hsml, Order)
                    if get_distance_from_wall(citizen)<minimumDistanceFromWall:
                        minimumDistanceFromWall = get_distance_from_wall(citizen)
    #Now add the contribution from all neighbors who are further away 
    #than the neighbor who is closest to the wall
    for [neighbor, distance] in neighbors:
        if distanceFromWall-distance<minimumDistanceFromWall:
            y -= hsml*kernel(distance, hsml, Order)
    return y