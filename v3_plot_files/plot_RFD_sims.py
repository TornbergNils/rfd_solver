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


print( mydict ) 


EME_x = np.fromfile( "./data/EME_x", dtype="double", count=-1 )
EME_y = np.fromfile( "./data/EME_y", dtype="double", count=-1 )
EME_z = np.fromfile( "./data/EME_z", dtype="double", count=-1 )

EMB_x = np.fromfile( "./data/EMB_x", dtype="double", count=-1 )
EMB_y = np.fromfile( "./data/EMB_y", dtype="double", count=-1 )
EMB_z = np.fromfile( "./data/EMB_z", dtype="double", count=-1 )

#RFD_x = np.fromfile( "./data/RFD_x", dtype="double", count=-1 ) 
#RFD_y = np.fromfile( "./data/RFD_y", dtype="double", count=-1 ) 
#RFD_z = np.fromfile( "./data/RFD_z", dtype="double", count=-1 ) 

electron_pos = np.fromfile( "./data/particle_electron", dtype="double", count=-1 ) 
positron_pos = np.fromfile( "./data/particle_positron", dtype="double", count=-1 ) 

J_x = np.fromfile( "./data/J_x", dtype="double", count=-1 ) 
J_y = np.fromfile( "./data/J_y", dtype="double", count=-1 ) 
J_z = np.fromfile( "./data/J_z", dtype="double", count=-1 ) 

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
print( n_frames ) 

EME_x = np.reshape(EME_x, ( n_frames, ny, nx ) )
EME_y = np.reshape(EME_y, ( n_frames, ny, nx ) )
EME_z = np.reshape(EME_z, ( n_frames, ny, nx ) )


EMB_x = np.reshape(EMB_x, ( n_frames, ny, nx ) )
EMB_y = np.reshape(EMB_y, ( n_frames, ny, nx ) )
EMB_z = np.reshape(EMB_z, ( n_frames, ny, nx ) )

#RFD_x = np.reshape(RFD_x, (  n_frames, n_particles  ) )
#RFD_y = np.reshape(RFD_y, (  n_frames, n_particles  ) )
#RFD_z = np.reshape(RFD_z, (  n_frames, n_particles  ) )

electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )
positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
                                          
J_x = np.reshape(J_x, ( n_frames, ny, nx ) )
J_y = np.reshape(J_y, ( n_frames, ny, nx ) )
J_z = np.reshape(J_z, ( n_frames, ny, nx ) )

rho_q = np.fromfile( "./data/rho_q", dtype="double", count=-1 ) 
rho_q = np.reshape(rho_q, ( n_frames, ny, nx ) )

##########################################
# Plotting setup
# First init entire figure
fig, (( ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

x_data = electron_pos[0, 0::3]
y_data = electron_pos[0, 1::3]
x_posit_data = positron_pos[0, 0::3]
y_posit_data = positron_pos[0, 1::3]
#RFDx_data = RFD_x[0, :]
#RFDy_data = RFD_y[0, :]
Ez_grid = EME_z[0, :, :]

extent = ( 0, nx*delta_x, 0, ny*delta_y)


# Plot density using histogram
density = rho_q[0,:,:]

im1 = ax1.imshow( density, extent=extent, vmax=np.max(density*0.8), aspect='auto' )
title1 = ax1.text(0.5, 0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
        transform=ax1.transAxes, ha='center' )

cbar1 = fig.colorbar(im1, ax = ax1)

# Initialize 2nd plot, positron histogram density
xbins = np.linspace(0, nx*delta_x, nx)
ybins = np.linspace(0, ny*delta_y, ny) 

p_density,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )

im2 = ax2.imshow( np.transpose(p_density), extent=extent, vmax=np.max(density*0.8), aspect='auto')
cbar2 = fig.colorbar(im2, ax = ax2)

# Initialize 3rd plot, E^2 + B^2
power_grid = ( np.square( EME_x[0,:,:] ) + np.square( EME_y[0,:,:] )
    + np.square( EME_z[0,:,:] ) + np.square( EMB_x[0,:,:] )
    + np.square( EMB_y[0,:,:] ) + np.square( EMB_z[0,:,:] ) )

ax3.text( 1.0, 1.08, "E^2 + B^2", c='w')
im3 = ax3.imshow( power_grid, extent=extent, aspect='auto' )

cbar3 = fig.colorbar(im3, ax = ax3)

# Initialize 4th plot, Ex plot

    
power_current = EME_x[0,:,:]
im4 = ax4.imshow( power_current, extent=extent, aspect='auto' )    
im4.set_data( power_current )
cbar4 = fig.colorbar(im4, ax = ax4)

vmax_j=1e-16
vmin_j=0
vmax_p=1e-16
vmin_p=0

###############################################################################

def init_anim():

    vmax=1
    vmin=1
    return


def update(frame):    
    global vmax_j, vmin_j, vmax_p, vmin_p

    x_posit_data = positron_pos[frame, 0::3]
    y_posit_data = positron_pos[frame, 1::3]

    # Update plot 1, interpolated charge density
    titlestring = str(frame * dt * save_rate)
    title1.set_text( "t = " + titlestring )
    density = rho_q[frame,:,:]
    im1.set_data( density )
    
    # Update plot 2, positron histo density
    p_density,temp2,temp3 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )
    im2.set_data( np.transpose(p_density) )
    im2.set_clim( np.min(p_density), np.max(p_density) )
    
    # Update plot3, E^2 + B^2
    power_grid = ( np.square( EME_x[frame,:,:] )
        + np.square( EME_y[frame,:,:] )
        + np.square( EME_z[frame,:,:] )
        + np.square( EMB_x[frame,:,:] )
        + np.square( EMB_y[frame,:,:] )
        + np.square( EMB_z[frame,:,:] ) )

    im3.set_data( power_grid )

    temp = np.max( power_grid )
    if( temp > vmax_p ):
        vmax_p = temp

    temp = np.min( power_grid )
    if( temp < vmin_p ):
        vmin_p = temp
    
    im3.set_clim(vmin_p, vmax_p )

    # Update plot4, Ex plot
    power_current = EME_x[frame,:,:]
    im4.set_data( power_current )
    
    temp = np.max( power_current )
    if( temp > vmax_j ):
        vmax_j = temp

    temp = np.min( power_current )
    if( temp < vmin_j ):
        vmin_j = temp
    
    im4.set_clim(vmin_j, vmax_j )
    
    progress = str( (frame / len( EME_x[:,0,0] ) ) * 100 ) + "%"
    print( progress )
    return


ani = anim.FuncAnimation(fig, update,
        frames=range(1, len( EME_x[:,0,0] ) ), 
        init_func=init_anim, repeat=False, blit=False)

myWriter=anim.FFMpegWriter( fps=5 )
ani.save("./figures/RFD_sims_evolution.mp4", writer=myWriter)
plt.close()


plt.imshow( EME_x[0,:,:] )
fname = "EME_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EME_y[0,:,:] )
fname = "EME_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EME_z[0,:,:] )
fname = "EME_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_x[0,:,:] )
fname = "EMB_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_y[0,:,:] )
fname = "EMB_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_z[0,:,:] )
fname = "EMB_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

"""
plt.imshow( RFD_x[:,0] )
fname = "RFD_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( RFD_y[:,0] )
fname = "RFD_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( RFD_z[:,0] )
fname = "RFD_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()
"""
