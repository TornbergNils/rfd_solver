import numpy as np
import matplotlib.pyplot as plt

E_data = np.fromfile( "./data/E_data_x.dat", dtype="double", count=-1 )
B_data = np.fromfile( "./data/B_data_x.dat", dtype="double", count=-1 )

my_E_grid = np.reshape( E_data, ( 4, 4 ) )
my_B_grid = np.reshape( B_data, ( 4, 4 ) )

print( my_E_grid )
print( my_B_grid )
