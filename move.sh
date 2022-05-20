#!/bin/bash

now=$(date +"%m-%d-%X")
filename='figures'$now

mkdir $filename
cp ./figures/* $filename
cp ./src/Experiment* $filename

echo 'Figures moved to new dir '$filename
