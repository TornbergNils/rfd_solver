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
J_x = np.fromfile( "./data/J_x", dtype="double", count=-1 ) 

J_x = np.reshape(J_x, ( n_frames, ny, nx ) )

fig1, ax1 = plt.subplots()
im1 = ax1.imshow( J_x[0,:,:], extent=extent )

cbar1 = fig1.colorbar(im1, ax = ax1)
im1.set_clim( np.min(J_x), np.max(J_x) )

def update_Jx(frame):
    im1.set_data( J_x[frame,:,:] )

ani = anim.FuncAnimation(fig1, update_Jx,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/Jx_everywhere.mp4", fps=5 )
plt.close()
print("Jx density plot done!" )

fig1, ax1 = plt.subplots()

power_current = J_x[0,:,:]
xrange = np.linspace(0, nx*delta_x, nx )
im21, = ax1.plot( xrange, power_current[int(ny/2),:] )
im22, = ax1.plot( xrange, np.mean(power_current[:, :], axis=0 ) )

Ne = 6.39*1e18
ampl = 0.1
q_e_cgs = -4.8*1e-10
wavelen = 1e-4
c_cgs = 2.99792*1e10
v_therm = 0.047
PI=3.14

def func( x ):
    return ampl*c_cgs*v_therm*q_e_cgs*Ne*np.sin(2*PI/wavelen*x)

im23, = ax1.plot( xrange, func( xrange ) )



def update_current_slice(frame):

    power_current = J_x[frame,:,:]
    im21.set_data(xrange, power_current[ int(ny/2),:] )
    im22.set_data(xrange, np.mean(power_current[:, :], axis=0 ) )
    v_max = np.max( power_current )
    v_min = np.min( power_current )

    #ax1.set_ylim( [v_min, v_max] )




ani = anim.FuncAnimation(fig1, update_current_slice,
        frames=range(1, n_frames ), repeat=False, blit=False)


ani.save("./figures/Jx_slice_evolution.mp4", fps=5 )
plt.close()
print("Jx slice plot done!" )