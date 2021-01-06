#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:47:22 2020

@author: leonard
"""
from sys import exit

from Constants import NORM_FAC
from Parameters import DIM, Kernel, DESNNGB

"""
This file contains the implementation for the choice and evaluation of 
smoothing kernels. The implementation is general enough so that different 
implementations of the kernels can be used by the code simply by specifying
the name of the kernel function by 'Kernel'.
"""
if Kernel == "cubic":     
    def kernel(x, h, derivative = False):
        if derivative:
            if x < 0.5:
                return -6 * NORM_FAC * (2 * x - 3 * x**2) * h**(-(DIM+1))
            elif x < 1:
                return -6 * NORM_FAC * (1-x)**2 * h**(-(DIM+1))
            else:
                return 0
        else:
            if x < 0.5:
                return NORM_FAC * (1-6*x**2 + 6*x**3) * h**(-DIM)
            elif x < 1:
                return 2 * NORM_FAC * (1-x)**3 * h**(-DIM)
            else:
                return 0
elif Kernel == "Wendland_C2":
    if DIM == 1:
        def kernel(x, h, derivative = False):
            if derivative:
                if x < 1:
                    return -12 * NORM_FAC * (1-x)**2 * x / h**2
                else:
                    return 0
            else:
                if x < 1:
                    return NORM_FAC * (1-x)**3 * (1 + 3 * x) /h
                else:
                    return 0
    else:
        def kernel(x, h, derivative = False):
            if derivative:
                if x < 1:
                    return -20 * NORM_FAC * (1-x)**3 * x * h**(-(DIM+1))
                else:
                    return 0
            else:
                if x < 1:
                    return NORM_FAC * (1-x)**4 * (1 + 4 * x) * h**(-DIM)
                else:
                    return 0 
else:
    print("Kernel function not defined!")
    exit()

def bias_correction(h):
    if DIM == 3 and Kernel == "Wendland_C2":
        return -0.0294 * (DESNNGB * 0.01)**(-0.977) * NORM_FAC/h**3
    else:
        return 0.0
