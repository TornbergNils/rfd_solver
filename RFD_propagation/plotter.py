import numpy as np
import matplotlib.pyplot as plt

E_data0_x = np.fromfile( "./data/E_data0_x.dat", dtype="double", count=-1 )
E_data1_x = np.fromfile( "./data/E_data1_x.dat", dtype="double", count=-1 )
E_data2_x = np.fromfile( "./data/E_data2_x.dat", dtype="double", count=-1 )
E_data3_x = np.fromfile( "./data/E_data3_x.dat", dtype="double", count=-1 )
E_data4_x = np.fromfile( "./data/E_data4_x.dat", dtype="double", count=-1 )
E_data5_x = np.fromfile( "./data/E_data5_x.dat", dtype="double", count=-1 )

E_data_y = np.fromfile( "./data/E_data0_y.dat", dtype="double", count=-1 )
E_data_z = np.fromfile( "./data/E_data0_z.dat", dtype="double", count=-1 )

B_data_x = np.fromfile( "./data/B_data0_x.dat", dtype="double", count=-1 )
B_data_y = np.fromfile( "./data/B_data0_y.dat", dtype="double", count=-1 )
B_data_z = np.fromfile( "./data/B_data0_z.dat", dtype="double", count=-1 )


#RFD_data_x = np.fromfile( "./data/RFD_x.dat", dtype="double", count=-1 ) 
#RFD_data_y = np.fromfile( "./data/RFD_y.dat", dtype="double", count=-1 ) 
#RFD_data_z = np.fromfile( "./data/RFD_z.dat", dtype="double", count=-1 ) 
#
#u_data = np.fromfile( "./data/u.dat", dtype="double", count=-1 ) 
#magnitudes = np.fromfile( "./data/magnitudes.dat", dtype="double", count=-1 ) 

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]

E_grid0_x = np.reshape( E_data0_x, ( ny, nx ) )
E_grid1_x =  np.reshape( E_data1_x, ( ny, nx ) )
E_data2_x =  np.reshape( E_data2_x, ( ny, nx ) )
E_data3_x =  np.reshape( E_data3_x, ( ny, nx ) )
E_data4_x =  np.reshape( E_data4_x, ( ny, nx ) )
E_data5_x =  np.reshape( E_data5_x, ( ny, nx ) )


E_grid_y = np.reshape( E_data_y, ( ny, nx ) )
E_grid_z = np.reshape( E_data_z, ( ny, nx ) )

B_grid_x = np.reshape( B_data_x, ( ny, nx ) )
B_grid_y = np.reshape( B_data_y, ( ny, nx ) )
B_grid_z = np.reshape( B_data_z, ( ny, nx ) )

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


myGrids = [E_grid0_x,E_grid1_x,E_data2_x,E_data3_x, E_data4_x, E_data5_x  ]

for grid, time in zip( myGrids, range(0, 6 ) ):

    fig, ax = plt.subplots()
    im = ax.pcolormesh( color_pos_x, color_pos_y, grid,
            shading='auto', cmap='coolwarm' )
    fig.colorbar(im, ax=ax)
    plt.savefig( "./figures/EM_field" + str(time) + ".png" )
    plt.close()
