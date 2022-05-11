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

mkdir ./figures/gauss
mkdir ./figures/wave
mkdir ./figures/Jz
fi

echo ./program.bin $args
./program.bin $args

echo Plotting!
python3 v3_plot_files/main_plot.py wave gauss Jz
#./plot.sh

echo Plotting done!
