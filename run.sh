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
#python3 ./plot_files/plot_boris_debug.py
#python3 ./plot_files/plot_1dslice_debug.py

#python3 ./v2_plot_files/plot_density.py
#python3 ./v2_plot_files/plot_Ey.py
python3 ./v2_plot_files/plot_Ex.py
python3 ./v2_plot_files/plot_vel_ic.py
python3 ./v2_plot_files/plot_scatter.py
#python3 ./v2_plot_files/plot_Jy.py
python3 ./v2_plot_files/plot_Jx.py
#python3 ./v2_plot_files/plot_interpol_density.py
python3 ./v2_plot_files/plot_momentum.py
python3 ./v2_plot_files/plot_freq_est.py

echo Plotting done!
