import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl

mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'Spectral'

with open('config.csv', mode='r') as infile:
    reader = csv.reader(infile)
    mydict = {rows[0]:rows[1] for rows in reader}

tmax = float(mydict['tmax'])
n_tsteps = int(mydict['n_tsteps'])
dt = float(mydict['dt'])
save_rate = int(mydict['save_rate'])
nx = int(mydict['nx'])
ny = int(mydict['ny'])
n_particles = int(mydict['n_particles'])
delta_x = float( mydict["delta_x"])
delta_y = float( mydict["delta_y"])

n_frames = int(  n_tsteps / save_rate ) + 1
extent = ( 0, nx*delta_x, 0, ny*delta_y)

# Load and format data
J_z = np.fromfile( "./data/J_z", dtype="double", count=-1 ) 

J_z = np.reshape(J_z, ( n_frames, ny, nx ) )
time = np.linspace(0, tmax, n_frames)

fig1, ax1 = plt.subplots()

c_cgs = 2.99792458 * 1e10
abstatampere_to_ampere = 10/c_cgs

J_z_sum = np.mean( J_z[:,int(nx*3/8):int(nx*5/8),
                         int(ny*3/8):int(ny*5/8) ],
                          axis=(1,2) )

conductor_area = (np.pi*(nx*delta_x/4)**2)
J_z_sum = J_z_sum*conductor_area*abstatampere_to_ampere

plt.plot(time, J_z_sum )

plt.savefig("./figures/Jz_max_evolution.png" )
plt.close()
print("Jz tot plot done!" )
