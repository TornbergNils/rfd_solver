#!/bin/bash

now=$(date +"%m-%d-%X")
filename='figures'$now

mkdir $filename
mv ./figures/* $filename

echo Files moved to new dir!
