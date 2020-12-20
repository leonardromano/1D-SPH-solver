#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 15:21:28 2020

@author: leonard
"""

def factorial(n):
    result = 1
    for i in range(2, n+1):
        result *= i
    return result