import numpy as np
import matplotlib.pyplot as plt

E_data = np.fromfile( "./data/E_data_x.dat", dtype="double", count=-1 )
B_data = np.fromfile( "./data/B_data_x.dat", dtype="double", count=-1 )

settings = np.loadtxt( "./data/header.csv", delimiter=',', dtype="int64" )

nx = settings[0]
ny = settings[1]

my_E_grid = np.reshape( E_data, ( ny, nx ) )
my_B_grid = np.reshape( B_data, ( ny, nx ) )

plt.imshow( my_E_grid )
plt.savefig( "./figures/E_grid.png" )
plt.close()

print( my_E_grid )
print( my_B_grid )
