#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPH simulation code by Leonard Romano
"""
from numpy.random import seed

"""
This file stores some constant global parameters.
Specify the parameters of your simulation in this file 
before running the simulation with 'run'.
"""

#Properties of Domain, general setup
TotalMass = 1.                                  #Sum of particle masses
LeftBoundary = -1.                              #Left spatial domain boundary
RightBoundary = 1.                              #Right spatial domain boundary  
TotalSideLength = (RightBoundary-LeftBoundary)  #Total size of spatial domain
InitialTime = 0                                 #Starting Time
FinalTime = 1.                                  #Final Time
NumberOfGridCells = 32                          #Number of grid cells for neighbor search
BinWidth = TotalSideLength/NumberOfGridCells    #Spatial extent of grid cells
NumberOfTimesteps = 128                         #Time refinement
NumberOfParticles = 256                         #Particle refinement
Mp = TotalMass/NumberOfParticles                #Particle Mass
Dt = (FinalTime - InitialTime)/NumberOfTimesteps#Maximum Timestep
Ntimebins = 16                                  #Maximum depth of timestep hierarchy

#initial conditions                                
seed(69420)                                     #Seed for random number generator
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
ICSpecifier = "from_file"                       
ICfile      = "../ICs/sph_16.txt"               #/path/to/file.txt

#SPH parameters
AdiabaticIndex = 1.4                            #adiabatic Index gamma P~rho^gamma
Order = "cubic"                                 #SPH Kernel function new Kernel options can be defined in Kernel.py
Viscosity = 1                                   #Viscosity Parameter \alpha
ViscositySoftening = 0.01                       #Softening parameter to prevent blow up of viscosity force
K = 4.                                          #Mean number of neighbors
Dh = K/16/(1+0.8*NumberOfParticles)             #Stepsize for smoothing length ~smaller than smallest particle spacing in IC
Hmin = Dh                                       #Minimum smoothing length
Hmax = TotalSideLength                          #Maximum smoothing length
CourantParameter = 0.3                          #Courant parameter

#Output Parameter
StepsBetweenSnapshots = 8                       #Number of timesteps between Snapshots