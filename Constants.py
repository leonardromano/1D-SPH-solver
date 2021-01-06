#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPH simulation code by Leonard Romano
"""
from Parameters import DIM, TotalMass, NPart, BoxSize, \
    InitialTime, FinalTime, NTimesteps, Kernel
from numpy.random import seed
from src.data.math_utility import factorial
from numpy import pi
from sys import exit

"""
This file stores some constant global parameters.
Specify the parameters of your simulation in this file 
before running the simulation with 'run'.
"""

#Integer coordinate utility
BITS_FOR_POSITIONS = 32
MAX_INT = (1 << BITS_FOR_POSITIONS) - 1

#tree utility
TREE_NUM_BEFORE_NODESPLIT = 3

#General setup
Mpart = TotalMass/NPart               #Particle Mass

#domain boundaries fixed to be [0, TotalSideLength]^DIM
FacIntToCoord = BoxSize/(1 << BITS_FOR_POSITIONS) #Conversion factor for int32 coords to double coords

#Time parameters
Dt = (FinalTime - InitialTime)/NTimesteps#Maximum Timestep

#initial conditions                                
seed(69420)                                     #Seed for random number generator

#######################################
#Geometrical constants and normalization constants
#Volume of DIM-dimensional unit ball
if DIM % 2 == 0:
    NORM_COEFF =pi**(DIM//2)/factorial(DIM//2)
    print("Set up NORM_COEFF = %g"%(NORM_COEFF))
else:
    NORM_COEFF = 2 * factorial((DIM - 1)//2) * (4 *pi)**((DIM-1)//2) / factorial(DIM)
    print("Set up NORM_COEFF = %g"%(NORM_COEFF))
    
#Kernel normalisation constant
if Kernel == "cubic":
    if DIM == 1:
        NORM_FAC = 4/3
        print("Set up NORM_FAC = %g"%NORM_FAC)
    elif DIM % 2 == 0:
        NORM_FAC = factorial(DIM + 3) * factorial(DIM//2 - 1) * pi**(-DIM//2)/12\
            /factorial(DIM-1)/(2 - 1/(1 << DIM))
        print("Set up NORM_FAC = %g"%NORM_FAC)
    else:
        NORM_FAC = factorial(DIM + 3) * pi**(-(DIM-1)//2) / 3/(DIM-1)/\
              factorial((DIM-1)//2 - 1)/((1 << (DIM+1)) - 1)
        print("Set up NORM_FAC = %g"%NORM_FAC)
elif Kernel == "Wendland_C2":
    if DIM == 1:
        NORM_FAC = 5/4
        print("Set up NORM_FAC = %g"%NORM_FAC)
    elif DIM % 2 == 0:
        NORM_FAC = factorial(DIM + 5) * DIM * factorial(DIM//2 - 1) * \
                   pi**(-DIM//2) / factorial(DIM + 1) / 240
        print("Set up NORM_FAC = %g"%NORM_FAC)
    else:
        NORM_FAC = factorial(DIM + 5) * (4 * pi)**(-(DIM-1)//2) / 120 / \
                   factorial((DIM-1)//2 - 1) / (DIM**2 - 1)
        print("Set up NORM_FAC = %g"%NORM_FAC)
else:
    print("Kernel function not defined!")
    exit()