import numpy as np
import matplotlib.pyplot as plt
import csv

with open('config.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('config_new.csv', mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[0]:rows[1] for rows in reader}


print( mydict ) 


EME_x = np.fromfile( "./data/EME_x", dtype="double", count=-1 )
EME_y = np.fromfile( "./data/EME_y", dtype="double", count=-1 )
EME_z = np.fromfile( "./data/EME_z", dtype="double", count=-1 )
                             
EMB_x = np.fromfile( "./data/EMB_x", dtype="double", count=-1 )
EMB_y = np.fromfile( "./data/EMB_y", dtype="double", count=-1 )
EMB_z = np.fromfile( "./data/EMB_z", dtype="double", count=-1 )

RFD_x = np.fromfile( "./data/RFD_x", dtype="double", count=-1 ) 
RFD_y = np.fromfile( "./data/RFD_y", dtype="double", count=-1 ) 
RFD_z = np.fromfile( "./data/RFD_z", dtype="double", count=-1 ) 

tmax = int(mydict['tmax'])
nx = int(mydict['nx'])
ny = int(mydict['ny'])

EME_x = np.reshape(EME_x, ( ny, nx ) )
EME_y = np.reshape(EME_y, ( ny, nx ) )
EME_z = np.reshape(EME_z, ( ny, nx ) )

EMB_x = np.reshape(EMB_x, ( ny, nx ) )
EMB_y = np.reshape(EMB_y, ( ny, nx ) )
EMB_z = np.reshape(EMB_z, ( ny, nx ) )

RFD_x = np.reshape(RFD_x, ( ny, nx ) )
RFD_y = np.reshape(RFD_y, ( ny, nx ) )
RFD_z = np.reshape(RFD_z, ( ny, nx ) )

plt.imshow( EME_x )
fname = "EME_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EME_y )
fname = "EME_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EME_z )
fname = "EME_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_x )
fname = "EMB_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_y )
fname = "EMB_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( EMB_z )
fname = "EMB_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( RFD_x )
fname = "RFD_x"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( RFD_y )
fname = "RFD_y"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.imshow( RFD_z )
fname = "RFD_z"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()
