#!/bin/bash


make

# To run a specific experiment, set args -Experiment_name
# To run all, set args -all, but note this takes time and about 10Gb storage
args="-wave"

# Dir is only relevant when args is not -all
dir="wave"

echo 'Running program using args: '
echo $args

all="-all"
if [ "$args" = "$all" ]; then
    mkdir ./data/gauss
    mkdir ./data/wave
    mkdir ./data/Jz
    mkdir ./data/Ez
    mkdir ./data/slab
    mkdir ./data/langm

    mkdir ./figures/gauss
    mkdir ./figures/wave
    mkdir ./figures/Jz
    mkdir ./figures/Ez
    mkdir ./figures/slab
    mkdir ./figures/langm

    py_args='gauss wave Jz Ez slab langm'
else
    mkdir ./data/$dir
    mkdir ./figures/$dir
    py_args=$dir
fi

echo ./program.bin $args
./program.bin $args

echo Plotting!
python3 Final_plotting_scripts/main_plot.py $py_args

echo Plotting done!
