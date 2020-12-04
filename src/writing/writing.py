#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:16:08 2020

@author: leonard
"""
from Constants import Mp, AdiabaticIndex, NumberOfParticles, FinalTime, \
    Viscosity, ViscositySoftening, K, CourantParameter, LeftBoundary, \
    RightBoundary, Output, FacIntToCoord
from sys import exit

def write_data(particles, label, Timestamp):
    "writes all the particle data in a text file. Eventually will be replaced by hdf5."
    f = open("%s/sph_%d.txt"%(Output,label), "w")
    if f:
        for particle in particles:
            position = (particle.position - 2**31)*FacIntToCoord
            f.write("%d %g %g %g %g %g %g %d\n" \
                    %(particle.index, position, particle.velocity, \
                      particle.entropy, particle.density, particle.pressure, \
                      particle.smoothingLength, len(particle.neighbors)))
        f.close()
    else:
        print("Results/sph_%d.txt could not be opened."%(label))
        exit()

def write_header():
    "writes some header information in a text file."
    f = open("%s/sph_header.txt"%(Output), "w")
    if f:
        f.write("Number of particles:%d\nFinal time:%g\nMass unit:%g\n\
Adiabatic index:%g\nViscosity:%g\nViscosity Softening:%g\n\
Average number of Neighbors:%g\nCourant parameter:%g\n\
Left boundary:%g\nRight boundary:%g" \
                %(NumberOfParticles, FinalTime, Mp, AdiabaticIndex, \
                  Viscosity, ViscositySoftening, K, CourantParameter, \
                  LeftBoundary, RightBoundary))
        f.close()
    else:
        print("Results/sph_header.txt could not be opened.")
        exit()