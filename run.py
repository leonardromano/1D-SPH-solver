#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPH simulation code by Leonard Romano
"""
from src.main import main
from time import time
from Constants import NumberOfParticles, NumberOfTimesteps

"""
This file is used to start the simulation. After specifying all the parameters
in the file 'Constants' simply run the code in this file to start the simulation.

"""

t0 = time()
main()
t1 = time()
print("Needed %01.3f s (%01.3f min) in order to simulate %d timesteps \
      for %d Particles." %(t1-t0, (t1-t0)/60, NumberOfTimesteps, \
          NumberOfParticles))