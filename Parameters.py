#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:09:40 2021

@author: leonard
"""
"""
SPH simulation code by Leonard Romano
"""

from numpy import array

"""
This file stores runtime parameters.
Specify the parameters of your simulation in this file 
before running the simulation with 'run'.
"""

#General setup
DIM               = 2              #Number of spatial dimensions
NPart             = 40**2          #Particle refinement
TotalMass         = 0.75           #Sum of particle masses


#domain boundaries fixed to be [0, TotalSideLength]^DIM
BoxSize          = array([2., 1.]) #Size of Box in physical units
Periodic = [True, False]           #Axis periodic?

#Time parameters
InitialTime = 0                    #Starting Time
FinalTime   = 1                    #Final Time
NTimesteps  = 128                  #Time refinement
NTimebins   = 16                   #Maximum depth of timestep hierarchy

#initial conditions                                
Noise_Var_Fac = 1e-3               #Noise variance (relative to particle spacing) for initial condition noise
"""
'ICSpecifier' specifies which ICs should be initialised. If ICs are 
loaded from a file should be "from_file". 'NumberOfParticles' should then be 
set to the number of particles in the file. Initial and final time should be 
set accordingly. The IC file is expected to be of the hdf5-type with groups 
"Header" containing an attribute "NumberOfParticles" and "PartData" containing 
the datasets "Coordinates", "Velocities" and "Entropy". Note however that 
smoothing length, density, pressure as well as number of neighbors aren't 
used for the ICs and thus can be ignored. If one reads in ICs like this one 
should make sure that all physical parameters like total mass, 
adiabatic index etc. are consistent.
"""
ICSpecifier = "from_file"                       #currently "shocktube" or "from_file"
ICfile      = "../ICs/initial_condition.txt"    #/path/to/file.txt
Output      = "../Results"                      #path to outputdirectory

#SPH parameters
AdiabaticIndex     = 1.4                        #adiabatic Index gamma P~rho^gamma
Kernel             = "cubic"                    #SPH Kernel function new Kernel options can be defined in Kernel.py
Viscosity          = 1.                         #Viscosity Parameter \alpha
ViscositySoftening = 0.01                       #Softening parameter to prevent blow up of viscosity force
DESNNGBS = 55                                   #Desired number of neighbors
NNGBSDEV = 1                                    #Limit how much the actual number of neighbors may deviate from the desired value
CourantParameter = 0.3                          #Courant parameter
TimestepLimiter =  0.01                         #Limiter for kinematic timestep

#gravity parameter
ExternalForce = True                            #Enable external forces true/false
Floor = True                                    #If true particles touching the floor don't experience gravity
GravAxis = 1                                    #Axis along which gravitational field lies
GravAcceleration = 0.5                          #magnitude of external force

#Output Parameter
OutputFrequency = 8                       #Number of timesteps between Snapshots