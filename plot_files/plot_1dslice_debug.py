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

tmax = int(mydict['tmax'])
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

##########################################
# Plotting setup
# First init entire figure
fig, (( ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

# initialize 1st plot, scatter on Ez
#filename = Ez_data_files[0]
#Ez_data = np.fromfile( filename, dtype="double", count=-1 )
#Ez_grid = np.reshape( Ez_data, ( ny, nx ) )
#
#filename = RFDx_data_files[0]
#RFDx_data = np.fromfile( filename, dtype="double", count=-1 )
#
#filename = RFDy_data_files[0]
#RFDy_data = np.fromfile( filename, dtype="double", count=-1 )
#
#
#filenamex = xpos_data_files[0]
#filenamey = ypos_data_files[0]
x_data = electron_pos[0, 0::3]
y_data = electron_pos[0, 1::3]
x_posit_data = electron_pos[0, 0::3]
y_posit_data = electron_pos[0, 1::3]
#RFDx_data = RFD_x[0, :]
#RFDy_data = RFD_y[0, :]
Ez_grid = EME_z[0, :, :]

extent = ( 0, nx*delta_x, 0, ny*delta_y)


# Plot density using histogram
density,edges1,edges2 = np.histogram2d(x_data, y_data, [nx, ny] )

#im1 = ax1.imshow( np.transpose(density), extent=extent, vmax=np.max(density*0.8) )
#xrange = np.linspace(0, ny*delta_y, len(density[int(nx/2), :] ) )
#im1, = ax1.plot(xrange, density[int(nx/2), :])
scatter1 = ax1.scatter( x_data[::50], y_data[::50], c='b', s=4 )

#scatter_pos = ax1.scatter( x_posit_data, y_posit_data, c='r', s=4 )
#quiver1 = ax1.quiver( x_data, y_data, RFDx_data, RFDy_data )
title1 = ax1.text(0.5, 0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
        transform=ax1.transAxes, ha='center' )

#cbar1 = fig.colorbar(im1, ax = ax1)

# Initialize 2nd plot, 3d scatterplot
#z_data = electron_pos[0, 2::3]
#
#ax2.remove()
#ax2 = fig.add_subplot(222, projection='3d' )
#scatter2 = ax2.scatter( x_data, y_data, z_data )
#ax2.set_xlabel("x")
#ax2.set_ylabel("y")
#ax2.set_zlabel("z")

p_density,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [nx, ny] )


power_current = J_y[0,:,:]
xrange = np.linspace(0, ny*delta_y, len(power_current[:, int(nx/2)] ) )
im21, = ax2.plot( xrange, power_current[:, int(nx/2)] )
im22, = ax2.plot( xrange, np.mean(power_current[:, :], axis=1 ) )
#im2 = ax2.imshow( np.transpose(p_density), extent=extent, vmax=np.max(density*0.8) )
#cbar2 = fig.colorbar(im2, ax = ax2)

# Initialize 3rd plot, E^2 + B^2
power_grid = ( np.square( EME_x[0,:,:] ) + np.square( EME_y[0,:,:] )
    + np.square( EME_z[0,:,:] ) + np.square( EMB_x[0,:,:] )
    + np.square( EMB_y[0,:,:] ) + np.square( EMB_z[0,:,:] ) )

ax3.text( 1.0, 1.08, "E^2 + B^2", c='w')
im3 = ax3.imshow( power_grid, extent=extent )
#scatter3 = ax3.scatter( x_data, y_data, c='b', s=3)

cbar3 = fig.colorbar(im3, ax = ax3)

# Initialize 4th plot, J plot

#x_for_J = np.linspace( 0, nx * delta_x, nx )
#y_for_J =  np.linspace( 0, ny * delta_y, ny )
#quiver4 = ax4.quiver( x_for_J, y_for_J, J_x[0,:], J_y[0,:] )
    
power_current = J_y[0,:,:]

#np.sqrt(( np.square( J_x[0,:,:] )
#    + np.square( J_y[0,:,:] )
#    + np.square( J_z[0,:,:] ) ))

power_current[0,0] = 50.0

xrange = np.linspace(0, ny*delta_y, len(power_current[:, int(nx/2)] ) )

#im4, = ax4.plot( xrange, power_current[:, int(nx/2)] )
im4 = ax4.imshow( power_current, extent=extent )    
#im4.set_data( power_current[int(nx/2), :] )

cbar4 = fig.colorbar(im4, ax = ax4)
vmax_j=0.1
vmin_j=0
vmax_p=0.1
vmin_p=0
#ax4.remove()
#ax4 = fig.add_subplot(224, projection='3d' )

#filename = Ex_data_files[0]
#Ex_data = np.fromfile( filename, dtype="double", count=-1 )
#Ex_grid = np.reshape( Ex_data, ( ny, nx ) )
#
#filename = Ey_data_files[0]
#Ey_data = np.fromfile( filename, dtype="double", count=-1 )
#Ey_grid = np.reshape( Ey_data, ( ny, nx ) )
#
#x,y = np.meshgrid( np.arange( -12, 12, 24/nx ), 
#        np.arange( -12, 12, 24/ny ) )
#quiver4 = ax4.quiver( x[::10], y[::10], Ex_grid[::10], Ey_grid[::10]  )

###############################################################################

def init_anim():

    vmax=1
    vmin=1
    return


def update(frame):    
    global vmax_j, vmin_j, vmax_p, vmin_p

    x_data = electron_pos[frame,0::3]
    y_data = electron_pos[frame,1::3]
    z_data = electron_pos[frame,2::3]
    x_posit_data = positron_pos[frame, 0::3]
    y_posit_data = positron_pos[frame, 1::3]
    # Note that currently RFD is only accesible for e-
    #RFDx_data = RFD_x[frame,:]
    #RFDy_data = RFD_y[frame,:]
    Ez_grid = EME_z[frame,:,:]

    # Update plot 1, scatter on Ez background
    titlestring = str(frame * dt * save_rate)
    title1.set_text( "t = " + titlestring )
    #density,edges1,edges2 = np.histogram2d(x_data, y_data, [nx, ny] )
    
    #im1.set_data( xrange, density[int(nx/2), :] )
    
    
    #temp = np.max( density )
    #if( temp > vmax ):
    #    vmax = temp

    #np.min( density )
    #if( temp < vmin ):
    #    vmin = temp
    
    #im1.set_clim(vmin, vmax )


    scatter1.set_offsets( np.c_[ x_data[::50], y_data[::50]] )
    #scatter_pos.set_offsets( np.c_[ x_posit_data, y_posit_data] )
    #quiver1.set_offsets( np.c_[ x_data, y_data] )
    #quiver1.set_UVC( RFDx_data, RFDy_data )


    # Update plot 2, 3d scatter plot
    #scatter2._offsets3d = ( x_data, y_data, z_data )
    
    #p_density,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [nx, ny] )
    #im2.set_data( np.transpose(p_density) )
    
    
    power_current = J_y[frame,:,:]
    im21.set_data(xrange, power_current[:, int(nx/2)] )
    im22.set_data(xrange, np.mean(power_current[:, :]) )
    
    # Update plot3, E^2 + B^2
    power_grid = ( np.square( EME_x[frame,:,:] )
        + np.square( EME_y[frame,:,:] )
        + np.square( EME_z[frame,:,:] )
        + np.square( EMB_x[frame,:,:] )
        + np.square( EMB_y[frame,:,:] )
        + np.square( EMB_z[frame,:,:] ) )


    #im3.set_data( power_grid )

    #temp = np.mean( power_grid )
    #vmax = temp * 1.5
    #vmin = np.min( power_grid )
    
    temp = np.max( power_grid )
    if( temp > vmax_p ):
        vmax_p = temp

    temp = np.min( power_grid )
    if( temp < vmin_p ):
        vmin_p = temp
    
    im3.set_clim(vmin_p, vmax_p )
    #ax3.set_ylim(vmin_p, vmax_p )
    #scatter3.set_offsets( np.c_[ x_data, y_data] )

    # Update plot4, quiver Ex Ey
    #filename = Ex_data_files[frame]
    #Ex_data = np.fromfile( filename, dtype="double", count=-1 )
    #Ex_grid = np.reshape( Ex_data, ( ny, nx ) )

    #filename = Ey_data_files[frame]
    #Ey_data = np.fromfile( filename, dtype="double", count=-1 )
    #Ey_grid = np.reshape( Ey_data, ( ny, nx ) )
    
    power_current = J_y[frame,:,:]
    #power_current = np.sqrt(( np.square( J_x[frame,:,:] )
    #    + np.square( J_y[frame,:,:] )
    #    + np.square( J_z[frame,:,:] ) ))
    
    #print( np.max( power_current))
    print( np.sum( J_y[frame,:,:] ))
    
    im4.set_data( power_current )
    #im4.set_data(xrange, power_current[:, int(nx/2)] )
    
    #temp = np.mean( power_grid )
    #vmax = temp * 1.5
    #vmin = np.min( power_grid )
    
    temp = np.max( power_current )
    if( temp > vmax_j ):
        vmax_j = temp

    temp = np.min( power_current )
    if( temp < vmin_j ):
        vmin_j = temp
    
    im4.set_clim(vmin_j, vmax_j )
    ax2.set_ylim( [vmin_j, vmax_j] )
    #quiver4.set_UVC( J_x[frame,:], J_y[frame,:] )
    
    progress = str( (frame / len( EME_x[:,0,0] ) ) * 100 ) + "%"
    print( progress )
    return


ani = anim.FuncAnimation(fig, update,
        frames=range(1, len( EME_x[:,0,0] ) ), 
        init_func=init_anim, repeat=False, blit=False)


ani.save("./figures/E_field_evolution.mp4", fps=5 )
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
