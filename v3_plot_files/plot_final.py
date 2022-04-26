import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl


mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'Spectral'


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

J_x = np.fromfile( "./data/J_x", dtype="double", count=-1 ) 
J_y = np.fromfile( "./data/J_y", dtype="double", count=-1 ) 
J_z = np.fromfile( "./data/J_z", dtype="double", count=-1 ) 

tmax = float(mydict['tmax'])
n_tsteps = int(mydict['n_tsteps'])
dt = float(mydict['dt'])
save_rate = int(mydict['save_rate'])
nx = int(mydict['nx'])
ny = int(mydict['ny'])
n_particles = int(mydict['n_particles'])
delta_x = float( mydict["delta_x"])
delta_y = float( mydict["delta_y"])

n_frames = int(  n_tsteps / save_rate ) + 1
print( n_frames ) 

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
                                          
J_x = np.reshape(J_x, ( n_frames, ny, nx ) )
J_y = np.reshape(J_y, ( n_frames, ny, nx ) )
J_z = np.reshape(J_z, ( n_frames, ny, nx ) )

rho_q = np.fromfile( "./data/rho_q", dtype="double", count=-1 ) 
rho_q = np.reshape(rho_q, ( n_frames, ny, nx ) )

## Settings
frame_start = 0
frame_stop = n_frames - 1
change_in_time = dt*save_rate*(frame_stop - frame_stop)
extent = ( 0, nx*delta_x, 0, ny*delta_y)

time = np.linspace(0, tmax, n_frames)

## Plot charge density

fig_q, ax_q = plt.subplots()

density = rho_q[frame_start,:,:]
im1 = ax_q.imshow( density, extent=extent, vmax=np.max(density), aspect='auto' )

fig_q.savefig("figures/charge_density_before.png")


density = rho_q[frame_stop,:,:]
im1 = ax_q.imshow( density, extent=extent, vmax=np.max(density), aspect='auto' )

fig_q.savefig("figures/charge_density_after.png")

## Plot positron density
fig_den, ax_den = plt.subplots()

xbins = np.linspace(0, nx*delta_x, nx)
ybins = np.linspace(0, ny*delta_y, ny) 

x_posit_data = positron_pos[0, 0::3]
y_posit_data = positron_pos[0, 1::3]

p_density_before,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )

x_posit_data = positron_pos[frame_stop, 0::3]
y_posit_data = positron_pos[frame_stop, 1::3]

p_density_after,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )

vmax=np.max( ( np.max(p_density_before*0.8),
               np.max(p_density_after*0.8) ))

im2 = ax_den.imshow( np.transpose(p_density_before), extent=extent, vmax=vmax, aspect='auto')
cbar2 = fig_den.colorbar(im2, ax = ax_den)


fig_den.savefig("figures/positron_density_before.png")


im2 = ax_den.imshow( np.transpose(p_density_after), extent=extent, vmax=vmax, aspect='auto')
fig_den.savefig("figures/positron_density_after.png")

## Plot 3 random electron and positron trajectories

fig_traj, ax_traj = plt.subplots()

# Electrons
elec1_trajx = electron_pos[0:20, 0]
elec1_trajy = electron_pos[0:20, 1]

elec2_trajx = electron_pos[0:20, 3]
elec2_trajy = electron_pos[0:20, 4]

elec3_trajx = electron_pos[0:20, 6]
elec3_trajy = electron_pos[0:20, 7]

# Positrons
posi1_trajx = positron_pos[0:20, 0]
posi1_trajy = positron_pos[0:20, 1]

posi2_trajx = positron_pos[0:20, 3]
posi2_trajy = positron_pos[0:20, 4]

posi3_trajx = positron_pos[0:20, 6]
posi3_trajy = positron_pos[0:20, 7]

ax_traj.plot( posi1_trajx,
     posi1_trajy, '--r' )
ax_traj.plot( posi2_trajx,
     posi2_trajy, '--r' )
ax_traj.plot( posi3_trajx,
     posi3_trajy, '--r' )

ax_traj.plot( elec1_trajx,
     elec1_trajy, 'b' )
ax_traj.plot( elec2_trajx,
     elec2_trajy, 'b' )
ax_traj.plot( elec3_trajx,
     elec3_trajy, 'b' )

#
#ax_traj.plot( [posi1_trajx, posi2_trajx, posi3_trajx],
#     [posi1_trajy, posi2_trajy, posi3_trajy], '--r' )
#
#ax_traj.plot( [elec1_trajx, elec2_trajx, elec3_trajx],
#     [elec1_trajy, elec2_trajy, elec3_trajy], 'b' )

fig_traj.savefig("figures/trajectories.png")
ax_traj.clear()

elec1_trajx = electron_pos[:, 0]
elec2_trajx = electron_pos[:, 3]
elec3_trajx = electron_pos[:, 6]
ax_traj.plot(time, elec1_trajx, 'b')
ax_traj.plot(time, elec2_trajx, 'b')
ax_traj.plot(time, elec3_trajx, 'b')
fig_traj.savefig("figures/elec_x_vs_time.png")
ax_traj.clear()

elec1_trajy = electron_pos[:, 1]
elec2_trajy = electron_pos[:, 4]
elec3_trajy = electron_pos[:, 7]
elec4_trajy = electron_pos[:, 10]
elec5_trajy = electron_pos[:, 13]
ax_traj.plot(time, elec1_trajy, 'b')
ax_traj.plot(time, elec2_trajy, 'b')
ax_traj.plot(time, elec3_trajy, 'b')
ax_traj.plot(time, elec4_trajy, 'b')
ax_traj.plot(time, elec5_trajy, 'b')
fig_traj.savefig("figures/elec_y_vs_time.png")


## Plot E, B field snapshots

def plot_EM( EM, frame, name ):    
     im = plt.imshow( EM[frame,:,:] )
     fname = name
     cbar = plt.colorbar( im )
     plt.savefig( "./figures/" + fname + ".png" )
     plt.close()

plot_EM(EME_x, frame_start, "EME_x")
plot_EM(EME_y, frame_start, "EME_y")
plot_EM(EME_z, frame_start, "EME_z")

plot_EM(EMB_x, frame_start, "EMB_x")
plot_EM(EMB_y, frame_start, "EMB_y")
plot_EM(EMB_z, frame_start, "EMB_z")

plot_EM(EME_x, frame_stop, "EME_x_after")
plot_EM(EME_y, frame_stop, "EME_y_after")
plot_EM(EME_z, frame_stop, "EME_z_after")

plot_EM(EMB_x, frame_stop, "EMB_x_after")
plot_EM(EMB_y, frame_stop, "EMB_y_after")
plot_EM(EMB_z, frame_stop, "EMB_z_after")

def plot_EM_history_left_middle(EM, time, name ):
     history = EM[:,int((nx*1)/4),int(ny/2)]
     plt.plot(time, history)
     fname = name
     plt.savefig( "./figures/" + fname + ".png" )
     plt.close()

plot_EM_history_left_middle(EME_x, time, "history_EME_x")
plot_EM_history_left_middle(EME_y, time, "history_EME_y")
plot_EM_history_left_middle(EME_z, time, "history_EME_z")
                                                
plot_EM_history_left_middle(EMB_x, time, "history_EMB_x")
plot_EM_history_left_middle(EMB_y, time, "history_EMB_y")
plot_EM_history_left_middle(EMB_z, time, "history_EMB_z")