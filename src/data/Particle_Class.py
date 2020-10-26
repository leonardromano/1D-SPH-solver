#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:15:30 2020

@author: leonard
"""
from Constants import LeftBoundary, BinWidth, Mp, Hmin, AdiabaticIndex, \
    CourantParameter, Dt
from numpy import sqrt

class Particle:
    "A class to store all the properties of a single particle"
    def __init__(self, position, velocity, entropy, index):
        
        self.position = position
        self.velocity = velocity
        self.velocity_ahead = 0
        self.acceleration_pressure = 0
        self.acceleration_viscosity = 0
        self.entropy = entropy
        self.entropy_ahead = 0
        self.entropyChange = 0
        self.mass = Mp
        self.index = index
        self.smoothingLength = Hmin
        self.density = 0
        self.pressure = 0
        self.soundspeed = 0
        self.courantTimestep = Dt
        self.neighbors = list()
        self.currentBin = self.get_current_bin()
    
    def get_current_bin(self):
        return int((self.position-LeftBoundary)//BinWidth)
    
    def update_current_bin(self):
        self.currentBin = self.get_current_bin()
        
    def update_pressure(self, ahead):
        if ahead == False:
            self.pressure = self.entropy*self.density**(AdiabaticIndex)
        else:
            self.pressure = self.entropy_ahead*self.density**(AdiabaticIndex)
        
    def update_soundspeed(self):
        self.soundspeed = sqrt(AdiabaticIndex*self.pressure/self.density)
    
    def update_courant_timestep(self):
        cmax = 0
        for neighbor in self.neighbors:
            if neighbor[0].soundspeed>cmax:
                cmax = neighbor[0].soundspeed
        self.courantTimestep =  CourantParameter*2*\
            self.smoothingLength/(self.soundspeed + cmax)
    