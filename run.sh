#!/bin/bash


make

# To run a specific experiment, set args -Experiment_name
# To run all, set args -all, but note this takes time and about 10Gb storage

# To only plot, comment out the line ./program.bin
experiment="-pulse"

# Model: -Boris or -RFD 
model="-RFD"
args=$experiment" "$model

# Dir is only relevant when args is not -all
dir="pulse-RFD"

echo 'Running program using args: '
echo $args

all="-all"
if [ "$experiment" = "$all" ]; then
    mkdir ./data/gauss$model
    mkdir ./data/wave$model
    mkdir ./data/Jz$model
    mkdir ./data/Ez$model
    mkdir ./data/slab$model
    mkdir ./data/langm$model
    mkdir ./data/pulse$model

    mkdir ./figures/gauss$model
    mkdir ./figures/wave$model
    mkdir ./figures/Jz$model
    mkdir ./figures/Ez$model
    mkdir ./figures/slab$model
    mkdir ./figures/langm$model
    mkdir ./figures/pulse$model

    py_args="gauss"$model" wave"$model" Jz"$model" Ez"$model" slab"$model" pulse"$model
else
    mkdir ./data/$dir
    mkdir ./figures/$dir
    py_args=$dir
fi

echo ./program.bin $args
./program.bin $args

echo Plotting! $py_args
python3 Final_plotting_scripts/main_plot.py $py_args

echo Plotting done!
