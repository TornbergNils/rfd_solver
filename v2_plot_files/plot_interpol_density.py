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
rho_q = np.fromfile( "./data/rho_q", dtype="double", count=-1 ) 
rho_q = np.reshape(rho_q, ( n_frames, nx, ny ) )

fig1, ax1 = plt.subplots()
im1 = ax1.imshow( rho_q[0,:,:], extent=extent, vmax=np.max(rho_q) )

cbar1 = fig1.colorbar(im1, ax = ax1)

def update_density(frame):
    im1.set_data( rho_q[frame,:,:] )

ani = anim.FuncAnimation(fig1, update_density,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/interpol_density_evolution.mp4", fps=5 )
plt.close()
print("Density interpol plot done!" )