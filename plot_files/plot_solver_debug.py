import numpy as np
import matplotlib.pyplot as plt
import csv

with open('config.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('config_new.csv', mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[0]:rows[1] for rows in reader}

print( mydict ) 
