#!/bin/bash

make

echo ./program.bin
./program.bin

echo Plotting!
#python3 ./plot_files/plot_RFD_debug.py
#python3 ./plot_files/plot_propagation_debug.py
#python3 ./plot_files/plot_FDTD_debug.py
#python3 ./plot_files/plot_solver_debug.py
#python3 ./plot_files/plot_interpol_debug.py
python3 ./plot_files/plot_boris_debug.py
python3 ./plot_files/plot_1dslice_debug.py

echo Plotting done!
