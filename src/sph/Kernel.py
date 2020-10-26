#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:47:22 2020

@author: leonard
"""
import numpy as np
from sys import exit
"""
This file contains the implementation for the choice and evaluation of 
smoothing kernels. The implementation is general enough so that different 
implementations of the kernels can be used by the code simply by specifying
the name of the kernel function by 'order' and the name of the derivative by
'del_order'. The function kernel then parses this expression and evaluates
the expression.
"""

def cubic(r, h):
    "cubic Kernel"
    x = r/h
    return (4/3/h)*(np.heaviside(0.5-x, 1)*(1-6*x**2 + 6*x**3) + \
                    np.heaviside(1-x, 0)*np.heaviside(x-0.5,0)*2*(1-x)**3)

def del_cubic(r,h):
    "derivative of cubic Kernel"
    x = r/h
    return -8/h * (np.heaviside(0.5-x, 1)*(2*x-3*x**2) + \
                    np.heaviside(1-x, 0)*np.heaviside(x-0.5,0)*(1-x)**2)
        

def kernel(r, h, order, derivative = False):
    "Returns the Kernel function of order 'order' or its derivative"
    if not derivative:
        try:
            return eval("%s(%g,%g)" %(order, abs(r), h))
        except NameError:
            print("%s(r, h) is not defined"%(order))
            exit()
    else:
        try:
            return np.sign(r)*eval("del_%s(%g,%g)"%(order, abs(r), h))
        except NameError:
            print("del_%s(r, h) is not defined"%(order))
            exit()