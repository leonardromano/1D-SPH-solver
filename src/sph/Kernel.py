#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:47:22 2020

@author: leonard
"""
from sys import exit
"""
This file contains the implementation for the choice and evaluation of 
smoothing kernels. The implementation is general enough so that different 
implementations of the kernels can be used by the code simply by specifying
the name of the kernel function by 'order' and the name of the derivative by
'del_order'. The function kernel then parses this expression and evaluates
the expression.
"""

def cubic(x):
    "cubic Kernel"
    if x < 0.5:
        return 4/3 *(1-6*x**2 + 6*x**3) 
    elif x < 1:
        return 8/3 * (1-x)**3
    else:
        return 0

def del_cubic(x):
    "derivative of cubic Kernel"
    if x < 0.5:
        return -8 * (2 * x - 3 * x**2)
    elif x < 1:
        return -8 * (1-x)**2
    else:
        return 0
        
def Wendland_C2(x):
    if x < 1:
        return 5/4 * (1-x)**3 * (1 + 3 * x)
    else:
        return 0

def del_Wendland_C2(x):
    if x < 1:
        return -15 * (1-x)**2 * x
    else:
        return 0
        

def kernel(r, h, order, derivative = False):
    "Returns the Kernel function of order 'order' or its derivative"
    if not derivative:
        try:
            return eval("%s(%g)" %(order, abs(r)/h))/h
        except NameError:
            print("%s(r, h) is not defined"%(order))
            exit()
    else:
        try:
            return eval("del_%s(%g)"%(order, abs(r)/h))/h**2
        except NameError:
            print("del_%s(r, h) is not defined"%(order))
            exit()