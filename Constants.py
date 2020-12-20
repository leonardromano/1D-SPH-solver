#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPH simulation code by Leonard Romano
"""
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
DIM       = 1                                   #number of spatial dimensions
BoundaryPeriodic = [False]                      #list of DIM bools stating whether ith boundary is periodic or not
TotalMass = 1.                                  #Sum of particle masses
NumberOfParticles = 1440                        #Particle refinement
Mp = TotalMass/NumberOfParticles                #Particle Mass

#domain boundaries fixed to be [0, TotalSideLength]^DIM
TotalSideLength = 2.                                  #Total size of spatial domain
FacIntToCoord = TotalSideLength/(1 << BITS_FOR_POSITIONS) #Conversion factor for int32 coords to double coords

#Time parameters
InitialTime = 0                                 #Starting Time
FinalTime = 0.35                                #Final Time
NumberOfTimesteps = 128                         #Time refinement
Dt = (FinalTime - InitialTime)/NumberOfTimesteps#Maximum Timestep
Ntimebins = 16                                  #Maximum depth of timestep hierarchy

#initial conditions                                
seed(69420)                                     #Seed for random number generator
Noise_Var_Fac = 1e-3                            #Noise variance (relative to particle spacing) for initial condition noise
"""
'ICSpecifier' specifies which ICs should be initialised. If ICs are 
loaded from a file should be "from_file". This then ignores NumberOfParticles.
Dh should then however be set to a sensible value. Also initial and final time
should be set accordingly. The format of the IC file should be just as the
output produced by this code. Note however that smoothing length, density, 
pressure as well as number of neighbors aren't used for the ICs and thus can 
be ignored. If one reads in ICs like this one should make sure that all physical
parameters like total mass, adiabatic index etc. are consistent.
"""
ICSpecifier = "shocktube"                     #currently "shocktube" or "from_file"
ICfile      = "../ICs/initial_condition.txt"    #/path/to/file.txt
Output      = "../Results"                      #path to outputdirectory

#SPH parameters
AdiabaticIndex = 1.4                            #adiabatic Index gamma P~rho^gamma
Kernel = "cubic"                                #SPH Kernel function new Kernel options can be defined in Kernel.py
Viscosity = 1.                                  #Viscosity Parameter \alpha
ViscositySoftening = 0.01                       #Softening parameter to prevent blow up of viscosity force
K = 55                                          #Desired number of neighbors
DK = 1                                          #Limit how much the actual number of neighbors may deviate from the desired value
CourantParameter = 0.3                          #Courant parameter
TimestepLimiter =  0.01                         #Limiter for kinematic timestep

#Output Parameter
StepsBetweenSnapshots = 8                       #Number of timesteps between Snapshots


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