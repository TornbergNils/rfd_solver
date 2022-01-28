import numpy as np
import matplotlib.pyplot as plt

E_data_x = np.fromfile( "./data/E_data_x.dat", dtype="double", count=-1 )
E_data_y = np.fromfile( "./data/E_data_y.dat", dtype="double", count=-1 )
E_data_z = np.fromfile( "./data/E_data_z.dat", dtype="double", count=-1 )

RFD_data_x = np.fromfile( "./data/RFD_x.dat", dtype="double", count=-1 ) 
RFD_data_y = np.fromfile( "./data/RFD_y.dat", dtype="double", count=-1 ) 
RFD_data_z = np.fromfile( "./data/RFD_z.dat", dtype="double", count=-1 ) 

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]

E_grid_x = np.reshape( E_data_x, ( ny, nx ) )
E_grid_y = np.reshape( E_data_y, ( ny, nx ) )
E_grid_z = np.reshape( E_data_z, ( ny, nx ) )

RFD_grid_x = np.reshape( RFD_data_x, ( ny, nx ) )
RFD_grid_y = np.reshape( RFD_data_y, ( ny, nx ) )
RFD_grid_z = np.reshape( RFD_data_z, ( ny, nx ) )

np.savetxt( './data/RFD_z.csv', RFD_grid_z, fmt='%.2f', delimiter=',')
np.savetxt( './data/E_data_x.csv', E_grid_x, fmt='%.2f', delimiter=',')
np.savetxt( './data/E_data_y.csv', E_grid_y, fmt='%.2f', delimiter=',')

# Create heatmap part of plot
arrow_pos_x = np.linspace( -3.0, 3.0, nx )
arrow_pos_y = np.linspace( -3.0, 3.0, ny )

fig, ax = plt.subplots()

ax.pcolormesh( arrow_pos_x, arrow_pos_y, RFD_grid_z, shading='auto', cmap='seismic_r' )

# Create quiver-plot part of RFD image
q = ax.quiver( arrow_pos_x[::3], arrow_pos_y[::3],
        RFD_grid_x[::3,::3], RFD_grid_y[::3,::3], scale=30 )


plt.savefig( "./figures/RFD_grid_quiver.png" )
