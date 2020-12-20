#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:35:37 2020

@author: leonard
"""
from Constants import FacIntToCoord, BITS_FOR_POSITIONS, DIM, \
    BoundaryPeriodic
from numpy import zeros 

def convert_to_int_position(position):
    """
    converts a coordinate with double accuracy to an integer position
    The left boundary corresponds to 0 and the right boundary to 2^32
    """
    intpos = zeros(DIM, dtype = int)
    for i in range(DIM):
        intpos[i] += int(position[i]/FacIntToCoord)
    return intpos

def convert_to_phys_position(intpos):
    "converts an integer coordinate vector to a physical position vector"
    pos = zeros(DIM, dtype = float)
    for i in range(DIM):
        pos[i] += intpos[i] * FacIntToCoord
    return pos

def find_minimum_offset_left(N):
    "finds the minimum offset i for which one has equal spacing for a region divided into 2^31 positions"
    for i in range(N):
        if ((1 << (BITS_FOR_POSITIONS-1)) - 2*i)%(N-1) == 0:
            return i+(1 << (BITS_FOR_POSITIONS-1))//N

def find_stepsize_right(N, dx_left):
    for i in range(N):
        if((1 << BITS_FOR_POSITIONS) - dx_left - i)%(N-1) == 0:
            return ((1 << BITS_FOR_POSITIONS) - dx_left - i)//(N-1)
        
        
def get_distance_vector(x, y):
    "returns the minimum distance vector taking into account periodic boundaries"
    dx = zeros(DIM, dtype = int)
    for i in range(DIM):
        if BoundaryPeriodic[i]:
            if x[i] <= y[i]:
                if abs(x[i]-y[i]) < (1 << BITS_FOR_POSITIONS) + x[i] - y[i]:
                    dx[i] += x[i] - y[i]
                else:
                    dx[i] += (1 << BITS_FOR_POSITIONS) + x[i] - y[i]
            else:
                if x[i] - y[i] < abs(x[i] - y[i] - (1 << BITS_FOR_POSITIONS)):
                    dx[i] += x[i] - y[i]
                else:
                    dx[i] += x[i] - y[i] - (1 << BITS_FOR_POSITIONS) 
        else:
            dx[i] += x[i] - y[i]
    return dx