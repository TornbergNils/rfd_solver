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

#Bx_data_files = [I for I in glob.glob("./data/time/B_data*_x.dat")]
#By_data_files = [I for I in glob.glob("./data/time/B_data*_y.dat")]
#Bz_data_files = [I for I in glob.glob("./data/time/B_data*_z.dat")]

Ex_data_files = natural_sort( Ex_data_files )
Ey_data_files = natural_sort( Ey_data_files )
Ez_data_files = natural_sort( Ez_data_files )

E_databefore_x = np.fromfile( "./data/E_databefore_x.dat",
        dtype="double", count=-1 )
E_databefore_y = np.fromfile( "./data/E_databefore_y.dat",
        dtype="double", count=-1 )
E_databefore_z = np.fromfile( "./data/E_databefore_z.dat",
        dtype="double", count=-1 )

B_databefore_x = np.fromfile( "./data/B_databefore_x.dat",
        dtype="double", count=-1 )
B_databefore_y = np.fromfile( "./data/B_databefore_y.dat",
        dtype="double", count=-1 )
B_databefore_z = np.fromfile( "./data/B_databefore_z.dat",
        dtype="double", count=-1 )

E_dataafter_x = np.fromfile( "./data/E_dataafter_x.dat",
        dtype="double", count=-1 )
E_dataafter_y = np.fromfile( "./data/E_dataafter_y.dat",
        dtype="double", count=-1 )
E_dataafter_z = np.fromfile( "./data/E_dataafter_z.dat",
        dtype="double", count=-1 )

B_dataafter_x = np.fromfile( "./data/B_dataafter_x.dat",
        dtype="double", count=-1 )
B_dataafter_y = np.fromfile( "./data/B_dataafter_y.dat",
        dtype="double", count=-1 )
B_dataafter_z = np.fromfile( "./data/B_dataafter_z.dat",
        dtype="double", count=-1 )

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )


nx = settings[0]
ny = settings[1]

E_grid_x_before = np.reshape( E_databefore_x, ( ny, nx ) )
E_grid_y_before = np.reshape( E_databefore_y, ( ny, nx ) )
E_grid_z_before = np.reshape( E_databefore_z, ( ny, nx ) )

B_grid_x_before = np.reshape( B_databefore_x, ( ny, nx ) )
B_grid_y_before = np.reshape( B_databefore_y, ( ny, nx ) )
B_grid_z_before = np.reshape( B_databefore_z, ( ny, nx ) )

E_grid_x_after = np.reshape( E_dataafter_x, ( ny, nx ) )
E_grid_y_after = np.reshape( E_dataafter_y, ( ny, nx ) )
E_grid_z_after = np.reshape( E_dataafter_z, ( ny, nx ) )
                                        
B_grid_x_after = np.reshape( B_dataafter_x, ( ny, nx ) )
B_grid_y_after = np.reshape( B_dataafter_y, ( ny, nx ) )
B_grid_z_after = np.reshape( B_dataafter_z, ( ny, nx ) )

np.savetxt( "./data/E_grid_z_before.csv", E_grid_z_before, fmt='%f.3' )
np.savetxt( "./data/E_grid_z_after.csv", E_grid_z_after, fmt='%f.3' )

gridlist = [ E_grid_x_before,
            E_grid_y_before,
            E_grid_z_before,
            B_grid_x_before,
            B_grid_y_before,
            B_grid_z_before,
            E_grid_x_after,
            E_grid_y_after,
            E_grid_z_after,
            B_grid_x_after,
            B_grid_y_after,
            B_grid_z_after ]

gridnames = [ "E_grid_x_before",
              "E_grid_y_before",
              "E_grid_z_before",
              "B_grid_x_before",
              "B_grid_y_before",
              "B_grid_z_before",
              "E_grid_x_after",
              "E_grid_y_after",
              "E_grid_z_after",
              "B_grid_x_after",
              "B_grid_y_after",
              "B_grid_z_after" ]


for grid, name in zip( gridlist, gridnames ):
    fig, ax = plt.subplots()
    im = ax.pcolormesh(  grid,
            shading='auto', cmap='coolwarm', vmin=-1.0, vmax=1.0 )
    fig.colorbar(im, ax=ax)
    
    
    plt.savefig( "./figures/" + name + ".png" )
    plt.close()

##########################################
# Plotting setup
# First init entire figure
fig, ax1 = plt.subplots()

# initialize 1st plot, scatter on Ez
filename = Ez_data_files[0]
Ez_data = np.fromfile( filename, dtype="double", count=-1 )
Ez_grid = np.reshape( Ez_data, ( ny, nx ) )

im = ax.pcolormesh(  Ez_grid,
            shading='auto', cmap='coolwarm', vmin=-1.0, vmax=1.0 )
fig.colorbar(im, ax=ax)

im1 = ax1.imshow( Ez_grid )
title1 = ax1.text(0.5, 0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
        transform=ax1.transAxes, ha='center' )


def init_anim():
    pass

def update(frame):    
    # Update plot 1, scatter on Ez background
    filename = Ez_data_files[frame]
    titlestring = str(re.findall( r"[-+]?\d*\.?\d+|\d+", filename ))
    title1.set_text( "t = " + titlestring )
    Ez_data = np.fromfile(filename, dtype="double", count=-1 )
    Ez_grid = np.reshape( Ez_data, ( ny, nx ) )
    im1.set_data( Ez_grid )
    print( "{:.2f}".format(float(frame) * 100 / float( len( Ez_data_files ) )), "%" )
    
ani = anim.FuncAnimation(fig, update,
        frames=range(1, len( Ez_data_files ) ), 
        init_func=init_anim, repeat=False, blit=False)


ani.save("./figures/E_field_evolution.mp4", fps=5 )
plt.close()
