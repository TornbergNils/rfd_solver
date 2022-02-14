import numpy as np
import matplotlib.pyplot as plt

E_data_x = np.fromfile( "./data/E_data_x.dat", dtype="double", count=-1 )
E_data_y = np.fromfile( "./data/E_data_y.dat", dtype="double", count=-1 )
E_data_z = np.fromfile( "./data/E_data_z.dat", dtype="double", count=-1 )

B_data_x = np.fromfile( "./data/B_data_x.dat", dtype="double", count=-1 )
B_data_y = np.fromfile( "./data/B_data_y.dat", dtype="double", count=-1 )
B_data_z = np.fromfile( "./data/B_data_z.dat", dtype="double", count=-1 )

RFD_data_x = np.fromfile( "./data/RFD_x.dat", dtype="double", count=-1 ) 
RFD_data_y = np.fromfile( "./data/RFD_y.dat", dtype="double", count=-1 ) 
RFD_data_z = np.fromfile( "./data/RFD_z.dat", dtype="double", count=-1 ) 

magnitudes = np.fromfile( "./data/magnitudes.dat", dtype="double", count=-1 ) 

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]

E_grid_x = np.reshape( E_data_x, ( ny, nx ) )
E_grid_y = np.reshape( E_data_y, ( ny, nx ) )
E_grid_z = np.reshape( E_data_z, ( ny, nx ) )

B_grid_x = np.reshape( B_data_x, ( ny, nx ) )
B_grid_y = np.reshape( B_data_y, ( ny, nx ) )
B_grid_z = np.reshape( B_data_z, ( ny, nx ) )

RFD_grid_x = np.reshape( RFD_data_x, ( ny, nx ) )
RFD_grid_y = np.reshape( RFD_data_y, ( ny, nx ) )
RFD_grid_z = np.reshape( RFD_data_z, ( ny, nx ) )


magnitudes_grid = np.reshape( magnitudes, (ny, nx ) )

print( magnitudes_grid ) 

np.savetxt( './data/RFD_z.csv', RFD_grid_z, fmt='%.2f', delimiter=',')
np.savetxt( './data/RFD_x.csv', RFD_grid_x, fmt='%.2f', delimiter=',')
np.savetxt( './data/E_data_x.csv', E_grid_x, fmt='%.2f', delimiter=',')
np.savetxt( './data/E_data_y.csv', E_grid_y, fmt='%.2f', delimiter=',')
np.savetxt( './data/magnitudes.csv', magnitudes_grid, fmt='%.2f', delimiter=',')


color_pos_x = np.linspace( -3.0, 3.0, nx+1 )
color_pos_y = np.linspace( -3.0, 3.0, ny+1 )

arrow_pos_x = np.linspace( -3.0, 3.0, nx )
arrow_pos_y = np.linspace( -3.0, 3.0, ny )

fig, ax = plt.subplots()



im = ax.pcolormesh( color_pos_x, color_pos_y, RFD_grid_z,
        shading='auto', cmap='coolwarm' )
fig.colorbar(im, ax=ax)
# Create quiver-plot part of RFD image
step1 = 1
q = ax.quiver( arrow_pos_x[::step1], arrow_pos_y[::step1],
        RFD_grid_x[::step1,::step1], RFD_grid_y[::step1,::step1], scale=30 )


plt.savefig( "./figures/RFD_grid_quiver.png" )
plt.close()


fig, ax = plt.subplots()

im = ax.pcolormesh(magnitudes_grid, shading='auto', cmap='coolwarm' )
fig.colorbar(im, ax=ax)
plt.savefig( "./figures/magnitudes.png" )
plt.close()

def save_quiver_plot( filename, grid_x, grid_y, step=1 ):
    # Create heatmap part of plot
    fig, ax = plt.subplots()
    q = ax.quiver( arrow_pos_x[::step], arrow_pos_y[::step],
        grid_x[::step,::step], grid_y[::step,::step], scale=30 )

    plt.savefig( filename )
    plt.close()
dirname = "./figures/"
save_quiver_plot( dirname + "E_xz.png", E_grid_x, E_grid_y, step1  )
save_quiver_plot( dirname + "B_xz.png", B_grid_x, B_grid_y, step1  )
