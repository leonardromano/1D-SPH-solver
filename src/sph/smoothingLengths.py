#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 19:13:02 2020

@author: leonard
"""
from Constants import Dh, K, Order, NumberOfGridCells, Hmax, BinWidth, \
    LeftBoundary, RightBoundary, DesRes
from src.sph.Kernel import kernel 
from numpy import sign

def set_smoothing_lengths(grid, particles):
    "For each particle find the smoothing length by minimizing the sph functional"
    for particle in particles:
        particle.smoothingLength = minimize_sph_functional(grid, particle)

def minimize_sph_functional(grid, particle):
    "Iteratively minimizes the sph functional"
    h = Dh
    if abs(sph_functional(grid, particle, h))<DesRes:
        return h
    else:
        i=0
        res=1
        resprev = 1
        hprev = h
        while abs(res)>DesRes and h<Hmax:
            #equal sign need to increment or decrement
            if sign(res) == sign(resprev):
                hprev = h
                h += sign(res)*Dh/2**i
                if h <= 0:
                    h = abs(hprev)/2
                resprev = res
                res = sph_functional(grid, particle, h)
            else: #different sign need to reduce steps
                htemp = h
                h = (hprev + h)/2
                hprev = htemp
                resprev = res
                res = sph_functional(grid, particle, h)
                i+=1
        if h>Hmax:
            return Hmax
        else:
            return h

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
        if distanceFromWall-distance+10**(-10)<minimumDistanceFromWall:
            y -= hsml*kernel(distance, hsml, Order)
    return y