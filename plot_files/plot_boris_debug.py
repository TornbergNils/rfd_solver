import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import Fit_sine

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

#RFD_x = np.fromfile( "./data/RFD_x", dtype="double", count=-1 ) 
#RFD_y = np.fromfile( "./data/RFD_y", dtype="double", count=-1 ) 
#RFD_z = np.fromfile( "./data/RFD_z", dtype="double", count=-1 ) 

electron_pos = np.fromfile( "./data/particle_electron", dtype="double", count=-1 ) 
positron_pos = np.fromfile( "./data/particle_positron", dtype="double", count=-1 ) 

tmax = float(mydict['tmax'])
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

#RFD_x = np.reshape(RFD_x, (  n_frames, n_particles  ) )
#RFD_y = np.reshape(RFD_y, (  n_frames, n_particles  ) )
#RFD_z = np.reshape(RFD_z, (  n_frames, n_particles  ) )

electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )
positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
                                          
# Plot particle trajectory
p1x_traj = electron_pos[:,0]
p2x_traj = electron_pos[:,3]
p3x_traj = electron_pos[:,6]
p4x_traj = electron_pos[:,9]

p1y_traj = electron_pos[:,1]
p2y_traj = electron_pos[:,4]
p3y_traj = electron_pos[:,7]
p4y_traj = electron_pos[:,10]

p5x_traj = electron_pos[:,12]
p6x_traj = electron_pos[:,15]
p7x_traj = electron_pos[:,18]
p8x_traj = electron_pos[:,21]

p5y_traj = electron_pos[:,13]
p6y_traj = electron_pos[:,16]
p7y_traj = electron_pos[:,19]
p8y_traj = electron_pos[:,22]



plt.plot( p1x_traj, p1y_traj )
plt.plot( p2x_traj, p2y_traj )
plt.plot( p3x_traj, p3y_traj )
plt.plot( p4x_traj, p4y_traj )
fname = "e_traj_cartesian1"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.plot( p5x_traj, p5y_traj )
plt.plot( p6x_traj, p6y_traj )
plt.plot( p7x_traj, p7y_traj )
plt.plot( p8x_traj, p8y_traj )
fname = "e_traj_cartesian2"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()


posi1x_traj = positron_pos[:,0]
posi2x_traj = positron_pos[:,3]
posi3x_traj = positron_pos[:,6]
posi4x_traj = positron_pos[:,9]

posi1y_traj = positron_pos[:,1]
posi2y_traj = positron_pos[:,4]
posi3y_traj = positron_pos[:,7]
posi4y_traj = positron_pos[:,10]

plt.plot( posi1x_traj, posi1y_traj )
plt.plot( posi2x_traj, posi2y_traj )
plt.plot( posi3x_traj, posi3y_traj )
plt.plot( posi4x_traj, posi4y_traj )
fname = "p_traj_cartesian2"
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

plt.plot( time, p1x_traj - p1x_traj[0] )
plt.plot( time, p2x_traj - p2x_traj[0] )
plt.plot( time, p3x_traj - p3x_traj[0] )
plt.plot( time, p4x_traj - p4x_traj[0] )
plt.legend()

fname = "x_vs_time"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

plt.plot( time, p1y_traj - p1y_traj[0] )
plt.plot( time, p2y_traj - p2y_traj[0] )
plt.plot( time, p3y_traj - p3y_traj[0] )
plt.plot( time, p4y_traj - p4y_traj[0] )

fname = "y_vs_time"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()

# that is, we take every sample_rate:th particle
x_pos1 = electron_pos[:,0] - np.mean( electron_pos[:,0])
x_pos2 = electron_pos[:,3] - np.mean( electron_pos[:,3])
x_pos3 = electron_pos[:,6] - np.mean( electron_pos[:,6])
x_pos4 = electron_pos[:,9] - np.mean( electron_pos[:,9])
x_pos5 = electron_pos[:,12] - np.mean( electron_pos[:,12])
#all_x_starts = electron_pos[0,0::3*sample_rate]

my_ps_es1 = np.square(np.abs( np.fft.rfft( x_pos1, 30 ) ))
my_ps_es2 = np.square(np.abs( np.fft.rfft( x_pos2, 30 ) ))
my_ps_es3 = np.square(np.abs( np.fft.rfft( x_pos3, 30 ) ))
my_ps_es4 = np.square(np.abs( np.fft.rfft( x_pos4, 30 ) ))
my_ps_es5 = np.square(np.abs( np.fft.rfft( x_pos5, 30 ) ))
powerspectrum = my_ps_es1 + my_ps_es2 + my_ps_es3 + my_ps_es4 + my_ps_es5 
time_step = dt * save_rate 
freqs = np.fft.rfftfreq( 30, time_step )
idx = np.argsort(freqs)

frq = freqs[idx]
pwr = powerspectrum[idx]

frq = frq[0:100]
pwr = pwr[0:100]


#plt.plot( frq[idx], pwr[idx] )
#plt.plot( freqs[idx], my_ps_es1[idx] )
#plt.plot( freqs[idx], my_ps_es2[idx] )
#plt.plot( freqs[idx], my_ps_es3[idx] )
#plt.plot( freqs[idx], my_ps_es4[idx] )
plt.plot( freqs[idx], powerspectrum[idx] )

fname = "powerspectrum"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()




results = Fit_sine.fit_sin( time[0:20], electron_pos[0:20,1] )
print( "Best guess for frequency is: ")
bestguess_omega = results["omega"]
print( bestguess_omega )


A=results["amp"]
w = results["omega"]
p = results["phase"]
offsetc  = results["offset"] 
f =results["freq"]
period = results["period"]
fitfunc =results["fitfunc"] 

funkvals = fitfunc( time[0:20] )

plt.plot(time[0:20], funkvals )
plt.plot(time[0:20], electron_pos[0:20,1] )
plt.legend(["guess", "real"])

fname = "Best guess"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()