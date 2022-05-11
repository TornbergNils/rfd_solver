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
electron_pos = np.fromfile( "./data/particle_electron", dtype="double", count=-1 ) 
electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )

x_data = electron_pos[0,0::3]
y_data = electron_pos[0,1::3]

density,edges1,edges2 = np.histogram2d(x_data, y_data, [nx, ny] )

fig1, ax1 = plt.subplots()
im1 = ax1.imshow( np.transpose(density), extent=extent, vmax=np.max(density), aspect='auto' )

cbar1 = fig1.colorbar(im1, ax = ax1)

def update_density(frame):
    x_data = electron_pos[frame,0::3]
    y_data = electron_pos[frame,1::3]

    density,edges1,edges2 = np.histogram2d(x_data, y_data, [nx, ny] )
    
    im1.set_data( density )
    im1.set_clim( 0, np.max(density) )

ani = anim.FuncAnimation(fig1, update_density,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/histo_density_evolution.mp4", fps=5 )
plt.close()
print("Density histo plot done!" )

def update_frac_density(frame):
    x_data = electron_pos[frame,0::3]
    y_data = electron_pos[frame,1::3]

    density,edges1,edges2 = np.histogram2d(x_data, y_data, [nx, ny] )
    density = density / np.mean( density )

    im1.set_data( np.transpose(density) )
    im1.set_clim( 0, np.max(density) )

ani = anim.FuncAnimation(fig1, update_frac_density,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/histo_frac_density_evolution.mp4", fps=5 )
plt.close()
print("Density frac plot done!" )