# rfd-solver

The project is aimed at developing a solver for simulating radiation-dominant electron-positron plasma dynamics in the limit of particles moving along radiation-free direction (RFD).  

# How the code works

To simply run the code for one of the predefine experiments,

1. Use the run.sh bash script which
2. Compiles to code,
3. Creates directories for data and figures depending on arguments specified in
the args variable in the script. The currently available experiments are
-langm,
-slab,
-Jz,
-Ez,
-wave,
-gauss
4. These can all be ran in the same go using -all, but this takes a lot of time
and hard drive space (>=8Gb).

The code consists of a main file which creates a child class of IC_struct
containing all the initial conditions for the simulation. A instance of
a solver class is then initialized using this IC_struct.

The solver class contains all the methods necessary for simulating the plasma
time evolution.

Depending on the variable use_RFD, the simulation is iterated using either the
boris stepper or RFD stepper via the methods Iterate_boris() or Iterate_RFD()
respectively.

These methods in turn call class methods which propagate particles, interpolate
grid quantities and propagate fields.
