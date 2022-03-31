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
EME_x = np.fromfile( "./data/EME_x", dtype="double", count=-1 )

EME_x = np.reshape(EME_x, ( n_frames, ny, nx ) )


fig1, ax1 = plt.subplots()
im1 = ax1.imshow( EME_x[0,:,:], extent=extent, vmax=np.max(EME_x) )

cbar1 = fig1.colorbar(im1, ax = ax1)

def update_density(frame):
    
    im1.set_data( EME_x[frame,:,:] )
    im1.set_clim( np.min(EME_x), np.max(EME_x) )

ani = anim.FuncAnimation(fig1, update_density,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/Ex_evolution.mp4", fps=5 )
plt.close()
print("Ex plot done!" )

fig1, ax1 = plt.subplots()

power_current = EME_x[0,:,:]
xrange = np.linspace(0, nx*delta_x, len(power_current[int(ny/2), :]) )
im21, = ax1.plot( xrange, power_current[int(ny/2),:] )
im22, = ax1.plot( xrange, np.mean(power_current[:, :], axis=0 ) )
#ax1.hlines( 5.35*1e4, 0, nx*delta_x )


def update_current_slice(frame):

    power_current = EME_x[frame,:,:]
    im21.set_data(xrange, power_current[ int(ny/2),:] )
    im22.set_data(xrange, np.mean(power_current[:, :], axis=0 ) )
    v_max = np.max( power_current )
    v_min = np.min( power_current )
    
    ax1.set_ylim( [v_min, v_max] )




ani = anim.FuncAnimation(fig1, update_current_slice,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/Ex_slice_evolution.mp4", fps=5 )
plt.close()
print("Ex slice plot done!" )