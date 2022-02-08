# MDsims
Molecular dynamics simulation code using velocity-verlet algorithm to simulate a polymer with one end grafted to an attractive surface.

## Code sructure: 
```
makefile = makefile
defs.h = definitions/initiation of all global variables
vv_md.c = main body of code/executable
measure.c, io.c = for functions and methods
integrate.c = Langevin Dynamics velocity-Verlet methods
gasdev.c = RNG
in = input parameters

Compiler: icc
```
## How to build:

Follow these steps:
```
make 
./vv
```
