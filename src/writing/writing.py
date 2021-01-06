#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:16:08 2020

@author: leonard
"""
import h5py
from numpy import zeros

from Constants import Mpart
from Parameters import AdiabaticIndex, NPart, FinalTime, Viscosity, \
    ViscositySoftening, DESNNGBS, CourantParameter, BoxSize, Output, DIM

from src.data.int_conversion import convert_to_phys_position


def write_data(particles, label, Time):
    "writes all the particle data in a hdf5 file."
    f = h5py.File("%s/sph_%d.hdf5"%(Output,label), "w")
    #first dump all the header info
    header = f.create_group("Header")
    h_att = header.attrs
    h_att.create("NumPart", NPart)
    h_att.create("FinalTime", FinalTime)
    h_att.create("Mass", Mpart)
    h_att.create("AdiabaticIndex", AdiabaticIndex)
    h_att.create("Viscosity", Viscosity)
    h_att.create("ViscositySoftening", ViscositySoftening)
    h_att.create("Time", Time)
    h_att.create("NumNeighbors", DESNNGBS)
    h_att.create("CourantParameter", CourantParameter)
    h_att.create("BoxSize", BoxSize)
    h_att.create("NDIM", DIM)
    #now make the data sets for the particle data
    IDs        = zeros((NPart), dtype = int)
    positions  = zeros((NPart, DIM), dtype = float)
    velocities = zeros((NPart, DIM), dtype = float)
    entropies  = zeros((NPart), dtype = float)
    densities  = zeros((NPart), dtype = float)
    pressures  = zeros((NPart), dtype = float)
    hsml       = zeros((NPart), dtype = float)
    NN         = zeros((NPart), dtype = int)
    for i in range(NPart):
        IDs[i]        += particles[i].index
        positions[i]  += convert_to_phys_position(particles[i].position)
        velocities[i] += particles[i].velocity
        entropies[i]  += particles[i].entropy
        densities[i]  += particles[i].density
        pressures[i]  += particles[i].pressure
        hsml[i]       += particles[i].hsml
        NN[i]         += len(particles[i].neighbors)
    f.create_dataset("PartData/IDs", data = IDs,  dtype = "u4")
    f.create_dataset("PartData/Coordinates", data = positions, dtype = "f4")
    f.create_dataset("PartData/Velocity", data = velocities, dtype = "f4")
    f.create_dataset("PartData/Entropy", data = entropies, dtype = "f4")
    f.create_dataset("PartData/Density", data = densities, dtype = "f4")
    f.create_dataset("PartData/Pressure", data = pressures, dtype = "f4")
    f.create_dataset("PartData/SmoothingLength", data = hsml, dtype = "f4")
    f.create_dataset("PartData/NumberOfNeighbors", data = NN, dtype = "f4")
    f.close()