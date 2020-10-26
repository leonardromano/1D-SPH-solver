# 1D-SPH-solver

A simulation can be started by executing run.py with python.
All the simulation parameters and initial conditions can be specified in Constants.py.

Currently there are two options to load initial conditions:
1) from a file that has the following pattern: ID position velocity entropy
2) from implemented initial conditions.
Currently only sod shocktube initial conditions are implemented but it is straightforward to 
implement other initial conditions by simply adding new options in Initialization.py.

Currently only a cubic B-spline kernel function is implemented. However the implementation of other
kernel functions is straightforward by simply defining the function and its derivative in Kernel.py
and preserving the pattern "name", "del_name" for the kernel's and its derivative's function name. 
The code parses the name of the function making choosing different kernel functions very easy.
