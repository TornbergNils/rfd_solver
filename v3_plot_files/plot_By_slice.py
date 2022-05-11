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
EMB_y = np.fromfile( "./data/EMB_y", dtype="double", count=-1 ) 

EMB_y = np.reshape(EMB_y, ( n_frames, ny, nx ) )

power_current = EMB_y[0,:,:]

fig1, ax1 = plt.subplots()

xrange = np.linspace(0, nx*delta_x, nx )
im1, = ax1.plot(xrange, power_current[ int(ny/2),:] )
im2, = ax1.plot(xrange, np.mean(power_current[:, :], axis=0 ) )

def update_B_slice(frame):

    power_current = EMB_y[frame,:,:]
    im1.set_data(xrange, power_current[ int(ny/2),:] )
    im2.set_data(xrange, np.mean(power_current[:, :], axis=0 ) )
    v_max = np.max( power_current )
    v_min = np.min( power_current )

    ax1.set_ylim( [v_min, v_max] )





ani = anim.FuncAnimation(fig1, update_B_slice,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/By_slice_evolution.mp4", fps=5 )
plt.close()
print("By slice plot done!" )
