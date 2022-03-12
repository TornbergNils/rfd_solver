import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv

# Load all data and read parameters from config.csv
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

RFD_x = np.fromfile( "./data/RFD_x", dtype="double", count=-1 ) 
RFD_y = np.fromfile( "./data/RFD_y", dtype="double", count=-1 ) 
RFD_z = np.fromfile( "./data/RFD_z", dtype="double", count=-1 ) 

electron_pos = np.fromfile( "./data/particle_electron", dtype="double", count=-1 ) 
positron_pos = np.fromfile( "./data/particle_positron", dtype="double", count=-1 ) 

tmax = int(mydict['tmax'])
n_tsteps = int(mydict['n_tsteps'])
dt = float(mydict['dt'])
save_rate = int(mydict['save_rate'])
nx = int(mydict['nx'])
ny = int(mydict['ny'])
n_particles = int(mydict['n_particles'])

n_frames = int(  n_tsteps / save_rate ) + 1

EME_x = np.reshape(EME_x, ( n_frames, ny, nx ) )
EME_y = np.reshape(EME_y, ( n_frames, ny, nx ) )
EME_z = np.reshape(EME_z, ( n_frames, ny, nx ) )

EMB_x = np.reshape(EMB_x, ( n_frames, ny, nx ) )
EMB_y = np.reshape(EMB_y, ( n_frames, ny, nx ) )
EMB_z = np.reshape(EMB_z, ( n_frames, ny, nx ) )

RFD_x = np.reshape(RFD_x, (  n_frames, n_particles  ) )
RFD_y = np.reshape(RFD_y, (  n_frames, n_particles  ) )
RFD_z = np.reshape(RFD_z, (  n_frames, n_particles  ) )

electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )
positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
                                          
# Plot particle trajectory
p1x_traj = electron_pos[:,500]
p1y_traj = electron_pos[:,500]

plt.plot( p1x_traj, p1y_traj )


fname = "e_traj_cartesian"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

radius = np.sqrt( electron_pos[:,0]**2 + electron_pos[:,1]**2 )
p1r_traj = radius - np.mean( radius )
p1phi_traj = np.abs(np.arctan( electron_pos[:,1] / electron_pos[:,0] ))

distance = np.sqrt( ( p1x_traj[0] - p1x_traj[1:-1] )**2 
        + ( p1y_traj[0] -p1y_traj[1:-1] )**2 )

gyroradius = np.max( distance ) / 2.0
print( "Gyroradius" )
print( gyroradius )

print( len( distance ) )
index_where_rotation_completes = 2
print( np.min( distance ) )
for ix in range( 100 , len( distance ) ):
    if( distance[ix] < 10**-2 ):
        index_where_rotation_completes = ix
        break

print( index_where_rotation_completes )

period = index_where_rotation_completes * save_rate * dt
ang_freq = 2*np.pi / period

print( "period and ang_freq:" )
print( period )
print( ang_freq )

plt.plot( p1r_traj, p1phi_traj )
fname = "e_traj_polar"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

time = np.linspace( 0, tmax, n_frames )
plt.plot( time[1:-1], distance )

fname = "distance_vs_time"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.plot( time, p1x_traj )

fname = "x_vs_time"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()
