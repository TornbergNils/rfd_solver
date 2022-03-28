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
EME_y = np.fromfile( "./data/EME_y", dtype="double", count=-1 )

EME_y = np.reshape(EME_y, ( n_frames, ny, nx ) )


fig1, ax1 = plt.subplots()
im1 = ax1.imshow( EME_y[0,:,:], extent=extent, vmax=np.max(EME_y) )

cbar1 = fig1.colorbar(im1, ax = ax1)

def update_density(frame):
    
    im1.set_data( EME_y[frame,:,:] )
    im1.set_clim( np.min(EME_y), np.max(EME_y) )

ani = anim.FuncAnimation(fig1, update_density,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/Ey_evolution.mp4", fps=5 )
plt.close()
print("Ey plot done!" )