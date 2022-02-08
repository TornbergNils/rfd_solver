import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import glob
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


E_data_files = [I for I in glob.glob("./data/E_data*_z.dat")]
xpos_data_files = [I for I in glob.glob("./data/pos_x*")]
ypos_data_files = [I for I in glob.glob("./data/pos_y*")]
E_data_files = natural_sort( E_data_files )
xpos_data_files = natural_sort( xpos_data_files )
ypos_data_files = natural_sort( ypos_data_files )
print( E_data_files )

#RFD_data_x = np.fromfile( "./data/RFD_x.dat", dtype="double", count=-1 ) 
#RFD_data_y = np.fromfile( "./data/RFD_y.dat", dtype="double", count=-1 ) 
#RFD_data_z = np.fromfile( "./data/RFD_z.dat", dtype="double", count=-1 ) 
#
#u_data = np.fromfile( "./data/u.dat", dtype="double", count=-1 ) 
#magnitudes = np.fromfile( "./data/magnitudes.dat", dtype="double", count=-1 ) 

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]

#RFD_grid_x = np.reshape( RFD_data_x, ( ny, nx ) )
#RFD_grid_y = np.reshape( RFD_data_y, ( ny, nx ) )
#RFD_grid_z = np.reshape( RFD_data_z, ( ny, nx ) )
#
#u_grid = np.reshape( u_data, ( ny, nx ) )
#
#magnitudes_grid = np.reshape( magnitudes, (ny, nx ) )


#np.savetxt( './data/RFD_z.csv', RFD_grid_z, fmt='%.2f', delimiter=',')
#np.savetxt( './data/RFD_x.csv', RFD_grid_x, fmt='%.2f', delimiter=',')
#np.savetxt( './data/E_data_x.csv', E_grid_x, fmt='%.2f', delimiter=',')
#np.savetxt( './data/E_data_y.csv', E_grid_y, fmt='%.2f', delimiter=',')
#np.savetxt( './data/magnitudes.csv', magnitudes_grid, fmt='%.2f', delimiter=',')


color_pos_x = np.linspace( -3.0, 3.0, nx+1 )
color_pos_y = np.linspace( -3.0, 3.0, ny+1 )

arrow_pos_x = np.linspace( -3.0, 3.0, nx )
arrow_pos_y = np.linspace( -3.0, 3.0, ny )


##########################################
# Plotting setup
fig, ax = plt.subplots()

filename = E_data_files[0]
E_data = np.fromfile( filename, dtype="double", count=-1 )
E_grid = np.reshape( E_data, ( ny, nx ) )

filename = xpos_data_files[0]
filename = ypos_data_files[0]
x_data = np.fromfile(filename, dtype="double", count=-1 )
y_data = np.fromfile(filename, dtype="double", count=-1 )

im = ax.imshow( E_grid )
scatter = ax.scatter( x_data, y_data, c='r' )

cbar = fig.colorbar(im, ax = ax)

def init_anim():
    pass

def update(frame):    
    filename = E_data_files[frame]
    E_data = np.fromfile(filename, dtype="double", count=-1 )
    filename = xpos_data_files[frame]
    x_data = np.fromfile(filename, dtype="double", count=-1 )
    filename = ypos_data_files[frame]
    y_data = np.fromfile(filename, dtype="double", count=-1 )
    
    E_grid = np.reshape( E_data, ( ny, nx ) )

    im = ax.imshow( E_grid, extent=[-6,6,-6,6] )
    scatter.set_offsets( np.c_[ x_data, y_data] )

ani = anim.FuncAnimation(fig, update,
        frames=range(1, len( E_data_files ) ), 
        init_func=init_anim, repeat=False, blit=False)


ani.save("./figures/E_field_evolution.mp4", fps=5 )
cbar.remove()
plt.close()
