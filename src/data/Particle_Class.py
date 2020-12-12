#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:15:30 2020

@author: leonard
"""
from Constants import Mp, AdiabaticIndex, CourantParameter, Dt, \
                      TimestepLimiter, FacIntToCoord
from numpy import sqrt

class Particle:
    "A class to store all the properties of a single particle"
    def __init__(self, position, velocity, entropy, index):
        #dynamic quantities
        self.position = position
        self.velocity = velocity
        self.velocity_ahead = 0
        self.acceleration = 0
        self.entropy = entropy
        self.entropy_ahead = 0
        self.entropyChange = 0
        #constants
        self.mass = Mp
        self.index = index
        #SPH quantities
        self.smoothingLength = 0
        self.dhsmlDensityFactor = 0
        self.density = 0
        self.pressure = 0
        self.soundspeed = 0
        #access variables
        self.timestepCriterion = Dt
        self.neighbors = list()
        self.timeBin = 0
        
    def update_pressure(self, ahead):
        if ahead == False:
            self.pressure = self.entropy*self.density**(AdiabaticIndex)
        else:
            self.pressure = self.entropy_ahead*self.density**(AdiabaticIndex)
        
    def update_soundspeed(self):
        if self.density > 0 and self.pressure > 0:
            self.soundspeed = sqrt(AdiabaticIndex*self.pressure/self.density)
        else:
            print("pressure = %g\ndensity = %g\nacceleration = \
%g\nvelocity = %g\nposition = %d"%(self.pressure, self.density, \
    self.acceleration, self.velocity, self.position * FacIntToCoord))
    
    def update_timestep_criterion(self):
        cmax = 0
        for neighbor in self.neighbors:
            if neighbor[0].soundspeed > cmax:
                cmax = neighbor[0].soundspeed
        courantTimestep =  CourantParameter*2*\
            self.smoothingLength/(self.soundspeed + cmax)
        if abs(self.acceleration) >0:
            kinematicTimestep = sqrt(2*TimestepLimiter*self.smoothingLength/\
                                 abs(self.acceleration))
        else:
            kinematicTimestep = courantTimestep
        self.timestepCriterion = min(kinematicTimestep, courantTimestep)
    