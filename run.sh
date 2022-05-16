#!/bin/bash


make

#echo 'write_param.py \n'
#python3 v3_plot_files/write_param.py

echo 'Using Experiment object to init!'

args="-all"

all="-all"
if [ "$args"="$all" ]; then
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
fi

echo ./program.bin $args
./program.bin $args

echo Plotting!
python3 v3_plot_files/main_plot.py wave gauss Jz Ez slab langm
#./plot.sh

echo Plotting done!
