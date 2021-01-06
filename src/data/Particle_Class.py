#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 20:15:30 2020

@author: leonard
"""
from numpy import sqrt, zeros
from numpy.linalg import norm

from Constants import Dt
from Parameters import AdiabaticIndex, CourantParameter, TimestepLimiter, DIM                     
                      
from src.data.int_conversion import convert_to_phys_position


class Particle:
    "A class to store all the properties of a single particle"
    def __init__(self, position, velocity, entropy, index):
        #dynamic quantities
        self.position = position
        self.velocity = velocity
        self.velocity_ahead = zeros(DIM, dtype = float)
        self.acceleration = zeros(DIM, dtype = float)
        self.entropy = entropy
        self.entropy_ahead = 0
        self.entropyChange = 0
        #constants
        self.index = index
        #SPH quantities
        self.hsml = 0
        self.dhsmldrho = 0
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
            output_string = "index = %d\npressure = %g\ndensity = %g\nacceleration = ("\
                             %(self.index, self.pressure, self.density)
            for a in self.acceleration:
                output_string += " %g "%(a)
            output_string += ")\nvelocity = ("
            for v in self.velocity:
                output_string += " %g "%(v)
            output_string += ")\nposition = ("
            for p in convert_to_phys_position(self.position):
                output_string += " %g "%(p)
            output_string += ")\n"
            print(output_string)
    
    def update_timestep_criterion(self):
        cmax = 0
        for neighbor in self.neighbors:
            if neighbor[0].soundspeed > cmax:
                cmax = neighbor[0].soundspeed
        courantTimestep =  CourantParameter * 2 *\
            self.smoothingLength/(self.soundspeed + cmax)
        if norm(self.acceleration) > 0:
            kinematicTimestep = sqrt(2 * TimestepLimiter * self.smoothingLength/\
                                 norm(self.acceleration))
        else:
            kinematicTimestep = courantTimestep
        self.timestepCriterion = min(kinematicTimestep, courantTimestep)
    