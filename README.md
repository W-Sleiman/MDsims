# MDsims
Molecular dynamics simulation code using velocity-verlet algorithm to simulate a polymer with one end grafted to an attractive surface.

Code sructure: 

makefile = makefile
defs.h = definitions/initiation of all global variables
vv_md.c = main body of code/executable
measure.c, io.c, integrate.c = for functions and methods
gasdev.c = RNG
in = input parameters

Compiler: icc

