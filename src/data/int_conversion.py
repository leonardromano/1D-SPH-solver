#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:35:37 2020

@author: leonard
"""
from Constants import LeftBoundary, FacIntToCoord

def convert_to_int_position(position):
    """
    converts a coordinate with double accuracy to an integer position
    The left boundary corresponds to 0 and the right boundary to 2^32
    """
    return int((position - LeftBoundary)//FacIntToCoord)

def find_minimum_offset_left(N):
    "finds the minimum offset i for which one has equal spacing for a region divided into 2^31 positions"
    for i in range(N):
        if (2**(31) - 2*i)%(N-1) == 0:
            return i

def find_stepsize_right(N, dx_left):
    for i in range(N):
        if(2**32 - dx_left - i)%(N-1) == 0:
            return (2**32 - dx_left - i)//(N-1)
        