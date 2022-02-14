import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.mplot3d import Axes3D
import glob
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

Ex_data_files = [I for I in glob.glob("./data/time/E_data*_x.dat")]
Ey_data_files = [I for I in glob.glob("./data/time/E_data*_y.dat")]
Ez_data_files = [I for I in glob.glob("./data/time/E_data*_z.dat")]

Bx_data_files = [I for I in glob.glob("./data/time/B_data*_x.dat")]
By_data_files = [I for I in glob.glob("./data/time/B_data*_y.dat")]
Bz_data_files = [I for I in glob.glob("./data/time/B_data*_z.dat")]

RFDx_data_files = [I for I in glob.glob("./data/time/RFD_x*")]
RFDy_data_files = [I for I in glob.glob("./data/time/RFD_y*")]
RFDz_data_files = [I for I in glob.glob("./data/time/RFD_z*")]

xpos_data_files = [I for I in glob.glob("./data/time/pos_x*")]
ypos_data_files = [I for I in glob.glob("./data/time/pos_y*")]
zpos_data_files = [I for I in glob.glob("./data/time/pos_z*")]

power_files = [I for I in glob.glob("./data/time/power_dense*")]

Ex_data_files = natural_sort( Ex_data_files )
Ey_data_files = natural_sort( Ey_data_files )
Ez_data_files = natural_sort( Ez_data_files )

RFDx_data_files = natural_sort( RFDx_data_files )
RFDy_data_files = natural_sort( RFDy_data_files )
RFDz_data_files = natural_sort( RFDz_data_files )

xpos_data_files = natural_sort( xpos_data_files )
ypos_data_files = natural_sort( ypos_data_files )
zpos_data_files = natural_sort( zpos_data_files )

power_files = natural_sort( power_files )

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]


##########################################
# Plotting setup
# First init entire figure
fig, (( ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

# initialize 1st plot, scatter on Ez
filename = Ez_data_files[0]
Ez_data = np.fromfile( filename, dtype="double", count=-1 )
Ez_grid = np.reshape( Ez_data, ( ny, nx ) )

filename = RFDx_data_files[0]
RFDx_data = np.fromfile( filename, dtype="double", count=-1 )

filename = RFDy_data_files[0]
RFDy_data = np.fromfile( filename, dtype="double", count=-1 )


filenamex = xpos_data_files[0]
filenamey = ypos_data_files[0]
x_data = np.fromfile(filenamex, dtype="double", count=-1 )
y_data = np.fromfile(filenamey, dtype="double", count=-1 )
im1 = ax1.imshow( Ez_grid, extent=(-12.0,12.0,-12.0,12.0) )
scatter1 = ax1.scatter( x_data, y_data, c='r', s=4 )
quiver1 = ax1.quiver( x_data, y_data, RFDx_data, RFDy_data )
title1 = ax1.text(0.5, 0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
        transform=ax1.transAxes, ha='center' )

cbar1 = fig.colorbar(im1, ax = ax1)

# Initialize 2nd plot, 3d scatterplot
filenamez = zpos_data_files[0]
z_data = np.fromfile( filenamez, dtype="double", count=-1 )

ax2.remove()
ax2 = fig.add_subplot(222, projection='3d' )
scatter2 = ax2.scatter( x_data, y_data, z_data )
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")


# Initialize 3rd plot, E^2 + B^2
filename = power_files[0]
power_data = np.fromfile( filename, dtype="double", count=-1 )
power_grid = np.reshape( power_data, ( ny, nx ) )

ax3.text( 1.0, 1.08, "E^2 + B^2", c='w')
im3 = ax3.imshow( power_grid, extent=(-12.0,12.0,-12.0,12.0) )
scatter3 = ax3.scatter( x_data, y_data, c='r', s=3)

cbar3 = fig.colorbar(im3, ax = ax3)

# Initialize 4th plot, 3d E field plot

#ax4.remove()
#ax4 = fig.add_subplot(224, projection='3d' )

filename = Ex_data_files[0]
Ex_data = np.fromfile( filename, dtype="double", count=-1 )
Ex_grid = np.reshape( Ex_data, ( ny, nx ) )

filename = Ey_data_files[0]
Ey_data = np.fromfile( filename, dtype="double", count=-1 )
Ey_grid = np.reshape( Ey_data, ( ny, nx ) )

x,y = np.meshgrid( np.arange( -12, 12, 24/nx ), 
                   np.arange( -12, 12, 24/ny ) )
#quiver4 = ax4.quiver( x[::10], y[::10], Ex_grid[::10], Ey_grid[::10]  )

###############################################################################

def init_anim():
    pass

def update(frame):    
    # Update plot 1, scatter on Ez background
    filename = Ez_data_files[frame]
    titlestring = str(re.findall( r"[-+]?\d*\.?\d+|\d+", filename ))
    title1.set_text( "t = " + titlestring )
    Ez_data = np.fromfile(filename, dtype="double", count=-1 )
    filename = xpos_data_files[frame]
    x_data = np.fromfile(filename, dtype="double", count=-1 )
    filename = ypos_data_files[frame]
    y_data = np.fromfile(filename, dtype="double", count=-1 )
    
    Ez_grid = np.reshape( Ez_data, ( ny, nx ) )

    filename = RFDx_data_files[frame]
    RFDx_data = np.fromfile( filename, dtype="double", count=-1 )
    
    filename = RFDy_data_files[frame]
    RFDy_data = np.fromfile( filename, dtype="double", count=-1 )

    im1.set_data( Ez_grid )
    scatter1.set_offsets( np.c_[ x_data, y_data] )
    quiver1.set_offsets( np.c_[ x_data, y_data] )
    quiver1.set_UVC( RFDx_data, RFDy_data )
    

    # Update plot 2, 3d scatter plot
    filename = zpos_data_files[frame]
    z_data = np.fromfile(filename, dtype="double", count=-1 )
    scatter2._offsets3d = ( x_data, y_data, z_data )

    # Update plot3, E^2 + B^2
    filename = power_files[frame]
    power_data = np.fromfile( filename, dtype="double", count=-1 )
    power_grid = np.reshape( power_data, ( ny, nx ) )
    
    im3.set_data( power_grid )
    scatter3.set_offsets( np.c_[ x_data, y_data] )

    # Update plot4, quiver Ex Ey
    filename = Ex_data_files[frame]
    Ex_data = np.fromfile( filename, dtype="double", count=-1 )
    Ex_grid = np.reshape( Ex_data, ( ny, nx ) )
    
    filename = Ey_data_files[frame]
    Ey_data = np.fromfile( filename, dtype="double", count=-1 )
    Ey_grid = np.reshape( Ey_data, ( ny, nx ) )

    #quiver4.set_UVC( Ex_grid[::10], Ey_grid[::10] )
    
    return


ani = anim.FuncAnimation(fig, update,
        frames=range(1, len( Ez_data_files ) ), 
        init_func=init_anim, repeat=False, blit=False)


ani.save("./figures/E_field_evolution.mp4", fps=5 )
plt.close()
