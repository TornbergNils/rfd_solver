import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches
import csv
import matplotlib as mpl
import fit_sine

# Plotting settings for matplotlib 
mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'Spectral'
mpl.rcParams["axes.formatter.limits"] = [-2,2]
## Examples
## To load EM:
# EME_x = np.fromfile( "./data/EME_x", dtype="double", count=-1 )
# EME_x = np.reshape( EME_x, ( n_frames, ny, nx ) )


def qdense_at_frame(mydict, frame, fname, rho, extent, vmin, vmax):
    fig_q, ax_q = plt.subplots()
    plt.tight_layout()
    density = rho[frame,:,:]
    im1 = ax_q.imshow( density, extent=extent,
            vmax=vmax,vmin=vmin, aspect='auto' )
    
    cbar = fig_q.colorbar(im1, ax = ax_q)

    save_rate = int(mydict['save_rate'])
    dt = float(mydict['dt'])
    time = frame*dt*save_rate

    filename = "figures/" + fname + "/qdense_f" + str(frame) + "t" + str(time) + ".png"
    fig_q.savefig(filename)


def plot_qdensity( mydict, n_frames, n_plots, data_dir, fname, nx, ny ):
    
    nx = int(mydict['nx'])
    ny = int(mydict['ny'])
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)

    rho = np.fromfile( data_dir + fname + "/rho_q", dtype="double", count=-1 ) 
    rho = np.reshape(rho, ( n_frames, ny, nx ) )
    vmax = np.max(rho )
    vmin = np.min(rho )


    for frame in range(0, n_frames, int(n_frames/n_plots) ):
        qdense_at_frame(mydict, frame, fname, rho, extent, vmin, vmax )

def particle_dense_at_frame(frame, fname, extent,
    position, mydict ):

    nx = int(mydict['nx'])
    ny = int(mydict['ny'])
    n_tsteps = int(mydict['n_tsteps'])
    dt = float(mydict['dt'])
    save_rate = int(mydict['save_rate'])
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    weight = float( mydict["weight"])
    
    n_frames = int(  n_tsteps / save_rate ) + 1
    fig_den, ax_den = plt.subplots()
    x_posit_data = position[frame, 0::3]
    y_posit_data = position[frame, 1::3]
    
    xbins = np.linspace(0, nx*delta_x, nx)
    ybins = np.linspace(0, ny*delta_y, ny) 

    p_density, edges1, edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )
    
    # TODO: Read weight into program and calculate actual density
    p_density = p_density * weight / ( delta_x * delta_y)
    vmax = np.max( p_density )
    #vmin = np.min( p_density )

    im2 = ax_den.imshow( np.transpose(p_density), extent=extent, vmax=vmax, aspect='auto')
    plt.tight_layout()
    cbar2 = fig_den.colorbar(im2, ax = ax_den)
    
    save_rate = int(mydict['save_rate'])
    dt = float(mydict['dt'])
    time = frame*dt*save_rate
    
    filename = "figures/" + fname + "/pdense_f"+ str(frame) + "t" + str(time) + ".png"
    fig_den.savefig(filename)


def plot_positron_density( mydict, n_plots, data_dir, fname ):

    nx = int(mydict['nx'])
    ny = int(mydict['ny'])
    n_tsteps = int(mydict['n_tsteps'])
    dt = float(mydict['dt'])
    n_particles = int(mydict['n_particles'])
    save_rate = int(mydict['save_rate'])
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)
    
    n_frames = int(  n_tsteps / save_rate ) + 1

    positron_pos = np.fromfile( data_dir + fname + "/particle_positron", dtype="double", count=-1 ) 
    positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
    position = positron_pos

    for frame in range(0, n_frames, int(n_frames/n_plots) ):
        particle_dense_at_frame(frame, fname, extent,
            position, mydict )


def plot_mean_density_v_time( mydict, data_dir, fname ):


    n_particles = int(mydict['n_particles'])
    nx = int(mydict['nx'])
    ny = int(mydict['ny'])
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    tmax = float(mydict['tmax'])
    xbins = np.linspace(0, nx*delta_x, nx)
    ybins = np.linspace(0, ny*delta_y, ny) 

    n_frames = int(  n_tsteps / save_rate ) + 1

    fig_dense_vs_time, ax_dense_vs_time = plt.subplots()
    plt.tight_layout()
    
    positron_pos = np.fromfile( data_dir + fname + "/particle_positron", dtype="double", count=-1 ) 
    positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
    position = positron_pos

    p_dense_at_left_middle = []
    for frame in range(0, n_frames):
         x_posit_data = position[frame, 0::3]
         y_posit_data = position[frame, 1::3]
         p_density_after,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )
         mean_dense_along_line = np.mean(p_density_after[int(nx/4), :])
         p_dense_at_left_middle.append(mean_dense_along_line)
    
    time = np.linspace(0, tmax, n_frames)
    ax_dense_vs_time.plot( time, np.array(p_dense_at_left_middle) )
    figure_name = "figures/"  + fname +  "/positron_density_v_time.png"
    fig_dense_vs_time.savefig(figure_name)

def plot_slice_density_along_middle( mydict, data_dir, fname, n_plots ):


    n_particles = int(mydict['n_particles'])
    nx = int(mydict['nx'])
    ny = int(mydict['ny'])
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    tmax = float(mydict['tmax'])
    dt = float(mydict['dt'])

    xbins = np.linspace(0, nx*delta_x, nx)
    ybins = np.linspace(0, ny*delta_y, ny) 
    n_frames = int(  n_tsteps / save_rate ) + 1

    positron_pos = np.fromfile( data_dir + fname + "/particle_positron", dtype="double", count=-1 ) 
    positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
    position = positron_pos

    for frame in range(0, n_frames, int(n_frames/n_plots) ):
        fig_dense_vs_time, ax_dense_vs_time = plt.subplots()
        plt.tight_layout()
        x_posit_data = position[frame, 0::3]
        y_posit_data = position[frame, 1::3]
        p_density_after,edges1,edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )
        mean_dense_along_line = np.mean( p_density_after[:, :], axis=1 )
    
        ax_dense_vs_time.plot( mean_dense_along_line )
        
        time = frame*dt*save_rate
        figure_name = "figures/"  + fname +  "/positron_density_slice_f"+ str(frame) + "t" + str(time) + ".png"
        fig_dense_vs_time.savefig(figure_name)
        plt.close()

class grid_mp4():
    def __init__(self, grid, extent, fname):
        self.fig, self.ax = plt.subplots()
        plt.xlabel("x, (cm)")
        plt.ylabel("y, (cm)")
        plt.title( fname[-4:-1] )
        self.grid = grid
        self.extent = extent
    
    def ani_init(self ):
        self.im = self.ax.imshow( self.grid[0,:,:], extent=self.extent, vmax=np.max(self.grid), aspect='auto' )
        self.cbar = self.fig.colorbar(self.im, ax = self.ax)
        self.im.set_clim( np.min(self.grid), np.max(self.grid) )
        plt.tight_layout()

    def ani_update(self, frame):
        self.im.set_data( self.grid[frame,:,:] )
        self.im.set_clim( np.min(self.grid[frame,:,:]), np.max(self.grid[frame,:,:]) )
    
    def animate(self, n_frames):
        self.ani = anim.FuncAnimation( self.fig, self.ani_update, init_func=self.ani_init,
            frames=range(1, n_frames ), repeat=False, blit=False)
    
    def ani_save(self, fname, grid_quantity):
        self.ani.save("./figures/"+ fname + grid_quantity +"_evolution.mp4", fps=5 )



def plot_grid_mp4( mydict, data_dir, fname, grid_quantity ):

    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)
    
    grid = np.fromfile( data_dir + fname + grid_quantity,
        dtype="double", count=-1 )
    grid = np.reshape(grid, ( n_frames, ny, nx ) )

    my_movie = grid_mp4(grid, extent, fname)
    my_movie.animate( n_frames )
    my_movie.ani_save(fname, grid_quantity)

    plt.close()


def load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, grid_quantity):
    grid = np.fromfile( data_dir + fname + grid_quantity,
        dtype="double", count=-1 )
    grid = np.reshape(grid, ( n_frames, ny, nx ) )
    return grid


def plot_energy_density_mp4( mydict, data_dir, fname ):

    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)
    
    energy_density = (
    load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_x")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_y")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_z")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_x")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_y")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_z")**2 )


    my_movie = grid_mp4(energy_density, extent, fname + "EB2")
    my_movie.animate( n_frames )
    my_movie.ani_save(fname, "/energy_density")

    plt.close()

def plot_trajectories( mydict, data_dir, fname, n_trajs, traj_len_fraction ):
    
    fig_traj, ax_traj = plt.subplots()

    n_particles = int(mydict['n_particles'])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)
    
    positron_pos = np.fromfile( data_dir + fname + "/particle_positron", dtype="double", count=-1 ) 
    positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
    
    electron_pos = np.fromfile( data_dir + fname + "/particle_electron", dtype="double", count=-1 ) 
    electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )

    
    i_particle = (np.random.randint(0, n_particles-1, size=n_trajs))*3

    traj_stop_frame = int(traj_len_fraction*n_frames)
    elec_x_trajs = np.transpose( electron_pos[0:traj_stop_frame:2, i_particle]   )
    elec_y_trajs = np.transpose( electron_pos[0:traj_stop_frame:2, i_particle+1] )

    posi_x_trajs = np.transpose( positron_pos[0:traj_stop_frame:2, i_particle]   )
    posi_y_trajs = np.transpose( positron_pos[0:traj_stop_frame:2, i_particle+1] )


    for trajx, trajy in zip( elec_x_trajs, elec_y_trajs ):
         label1 = ax_traj.plot( trajx,
              trajy, '.b' )
    
    for trajx, trajy in zip( posi_x_trajs, posi_y_trajs ):
         label2 = ax_traj.plot( trajx,
              trajy, '.r' )
    
    plt.title("Electron and positron trajectories")
    red_patch = mpatches.Patch(color='red', label='Positrons')
    blue_patch = mpatches.Patch(color='blue', label='Electrons')
    plt.legend(handles=[red_patch,blue_patch])

    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    figure_name = "figures/"  + fname +  "/particle_trajectories" + ".png"
    fig_traj.savefig(figure_name)

def plot_velocities( mydict, data_dir, fname, n_trajs, traj_len_fraction ):
    
    fig_vel, ax_vel = plt.subplots()

    n_particles = int(mydict['n_particles'])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    tmax = float(mydict['tmax'])
    
    electron_vel = np.fromfile( data_dir + fname + "/particle_e_velocities", dtype="double", count=-1 ) 
    electron_vel = np.reshape(electron_vel, ( n_frames, 3*n_particles ) )

    i_particle = (np.random.randint(0, n_particles-1, size=n_trajs))*3

    elec_x_vels = np.transpose( electron_vel[:, i_particle]   )
    elec_y_vels = np.transpose( electron_vel[:, i_particle+1] )


    time = np.linspace(0, tmax, n_frames)
    for vel_x in elec_x_vels:
         ax_vel.plot( time,
              vel_x, '.b' )
    
    plt.title( "Electron velocity along x-axis vs time")
    plt.xlabel( "time (s)" )
    plt.ylabel( "x-velocity (cm/s)" )
    figure_name = "figures/"  + fname +  "/electron_x_vel_v_time" + ".png"
    fig_vel.savefig(figure_name)
    plt.close()
    
    fig_vel, ax_vel = plt.subplots()
    for vel_y in elec_y_vels:
         ax_vel.plot( time,
              vel_y, '.r' )
    
    plt.title( "Electron velocity along y-axis vs time")
    plt.xlabel( "time (s)" )
    plt.ylabel( "y-velocity (cm/s)" )
    figure_name = "figures/"  + fname +  "/electron_y_vel_v_time" + ".png"
    fig_vel.savefig(figure_name)
    plt.close()

def plot_grid_snapshot(mydict, data_dir, fname, EM, frame, grid_quantity ):    
    figt, axt = plt.subplots()
    im = axt.imshow( EM[frame,:,:] )
    cbar = figt.colorbar( im )
    
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1
    
    save_rate = int(mydict['save_rate'])
    dt = float(mydict['dt'])
    time = frame*dt*save_rate
    
    plt.title( grid_quantity[-3:-1] + " at " + str( frame * save_rate ) + " of " + str( n_frames ) )
    plt.xlabel( "x (cm)" )
    plt.xlabel( "y (cm)" )
    figt.savefig( "./figures/" + fname + grid_quantity + "_f" + str(frame) + "t" + str(time) + ".png" )
    axt.clear()

def plot_grid_evolution(mydict, data_dir, fname, n_plots, grid_quantity ):
    
    n_particles = int(mydict['n_particles'])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)

    EM = load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, grid_quantity)
    for frame in range(0, n_frames, int(n_frames/n_plots)):
        plot_grid_snapshot(mydict, data_dir, fname, EM, frame, grid_quantity )


def plot_pos_v_time( mydict, data_dir, fname, n_trajs, traj_len_fraction ):
    
    fig_pos, ax_pos = plt.subplots()

    n_particles = int(mydict['n_particles'])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    tmax = float(mydict['tmax'])
    
    electron_pos = np.fromfile( data_dir + fname + "/particle_electron", dtype="double", count=-1 ) 
    electron_pos = np.reshape(electron_pos, ( n_frames, 3*n_particles ) )

    i_particle = (np.random.randint(0, n_particles-1, size=n_trajs))*3

    elec_x_pos = np.transpose( electron_pos[:, i_particle]   )
    elec_y_pos = np.transpose( electron_pos[:, i_particle+1] )


    time = np.linspace(0, tmax, n_frames)
    for pos_x in elec_x_pos:
         ax_pos.plot( time,
              pos_x, '.b' )
    
    figure_name = "figures/"  + fname +  "/electron_x_pos_v_time" + ".png"
    
    plt.title("Electron x position vs time")
    plt.ylabel("Position x (cm)")
    plt.xlabel("Time since sim. start (s)")
    fig_pos.savefig(figure_name)
    plt.close()
    
    fig_pos, ax_pos = plt.subplots()
    for pos_y in elec_y_pos:
         ax_pos.plot( time,
              pos_y, '.r' )
    
    plt.title("Electron y position vs time")
    plt.ylabel("Position y (cm)")
    plt.xlabel("Time since sim. start (s)")
    figure_name = "figures/"  + fname +  "/electron_y_pos_v_time" + ".png"
    fig_pos.savefig(figure_name)
    plt.close()

def plot_grid_snapshot(mydict, data_dir, fname, EM, frame, grid_quantity ):    
    figt, axt = plt.subplots()
    
    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)
    
    im = axt.imshow( EM[frame,:,:], extent=extent )
    cbar = figt.colorbar( im )
    
    save_rate = int(mydict['save_rate'])
    dt = float(mydict['dt'])
    time = frame*dt*save_rate
    
    figt.savefig( "./figures/" + fname + grid_quantity + "_f" + str(frame) + "t" + str(time) + ".png" )
    axt.clear()

class multi_movie():
    def __init__(self, grid1, grid2, grid3, grid4, extent ):
        self.fig, ((self.ax1, self.ax2),
                   (self.ax3, self.ax4 )) = plt.subplots(2, 2)
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid3 = grid3
        self.grid4 = grid4
        self.extent = extent
    
    def ani_init( self ):
        self.im1 = self.ax1.imshow( self.grid1[0,:,:], extent=self.extent,
                                   vmax=np.max(self.grid1), aspect='auto')
        self.cbar1 = self.fig.colorbar(self.im1, ax = self.ax1)
        plt.tight_layout()
        
        self.im2 = self.ax2.imshow( self.grid2[0,:,:], extent=self.extent,
                                   vmax=np.max(self.grid2), aspect='auto')
        self.cbar2 = self.fig.colorbar(self.im2, ax = self.ax2)
        self.im2.set_clim( np.min(self.grid2), np.max(self.grid2) )
        plt.tight_layout()
        
        self.im3 = self.ax3.imshow( self.grid3[0,:,:], extent=self.extent,
                                   vmax=np.max(self.grid3), aspect='auto')
        self.cbar3 = self.fig.colorbar(self.im3, ax = self.ax3)
        self.im3.set_clim( np.min(self.grid3), np.max(self.grid3) )
        plt.tight_layout()
        
        self.im4 = self.ax4.imshow( self.grid4[0,:,:], extent=self.extent,
                                   vmax=np.max(self.grid4), aspect='auto')
        self.cbar4 = self.fig.colorbar(self.im4, ax = self.ax4)
        self.im4.set_clim( np.min(self.grid4), np.max(self.grid4) )
        plt.tight_layout()
    
    def ani_update(self, frame):
        self.im1.set_data( self.grid1[frame,:,:] )
        self.im1.set_clim( np.min(self.grid1[frame,:,:]), np.max(self.grid1[frame,:,:]) )
        self.im2.set_data( self.grid2[frame,:,:] )
        self.im2.set_clim( np.min(self.grid2[frame,:,:]), np.max(self.grid2[frame,:,:]) )
        self.im3.set_data( self.grid3[frame,:,:] )
        self.im3.set_clim( np.min(self.grid3[frame,:,:]), np.max(self.grid3[frame,:,:]) )
        self.im4.set_data( self.grid4[frame,:,:] )
        self.im4.set_clim( np.min(self.grid4[frame,:,:]), np.max(self.grid4[frame,:,:]) )
    
    def animate(self, n_frames):
        self.ani = anim.FuncAnimation( self.fig, self.ani_update, init_func=self.ani_init,
            frames=range(1, n_frames ), repeat=False, blit=False)
    
    def ani_save(self, fname, grid_quantity):
        self.ani.save("./figures/"+ fname + grid_quantity + ".mp4", fps=5 )
        

def plot_multi_movie_mp4( mydict, data_dir, fname ):

    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1
    weight = float( mydict["weight"])
    n_particles = int(mydict['n_particles'])

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    extent = ( 0, nx*delta_x, 0, ny*delta_y)

    rho = np.fromfile( data_dir + fname + "/rho_q", dtype="double", count=-1 ) 
    rho = np.reshape(rho, ( n_frames, ny, nx ) )

    positron_pos = np.fromfile( data_dir + fname + "/particle_positron", dtype="double", count=-1 ) 
    positron_pos = np.reshape(positron_pos, ( n_frames, 3*n_particles ) )
    
    
    xbins = np.linspace(0, nx*delta_x, nx)
    ybins = np.linspace(0, ny*delta_y, ny) 
    
    temp_list = []
    for ix in range( n_frames ):
        x_posit_data = positron_pos[ix, 0::3]
        y_posit_data = positron_pos[ix, 1::3]

        p_density, edges1, edges2 = np.histogram2d(x_posit_data, y_posit_data, [xbins, ybins] )
        temp_list.append( np.transpose(p_density))
        # Calculate density in cgs units
        p_density = p_density * weight / ( delta_x * delta_y)
    
    p_density_grid = np.array( temp_list )



    energy_density = (
    load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_x")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_y")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EMB_z")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_x")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_y")**2
    + load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/EME_z")**2 )

    J_z = load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, "/J_z" )


    my_movie = multi_movie( rho, p_density_grid, energy_density, J_z, extent )
    my_movie.animate( n_frames )
    my_movie.ani_save(fname, "/multimovie")

    plt.close()

def plot_EM_slice( mydict, data_dir, fname, grid_quantity, n_plots ):

    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    nx = int(mydict['nx'])
    ny = int(mydict['ny']) 
    delta_x = float( mydict["delta_x"])
    delta_y = float( mydict["delta_y"])
    dt = float( mydict["dt"])

    EM = load_and_reshape_EM( data_dir, fname, nx, ny, n_frames, grid_quantity )
    xvalues = np.linspace( 0, nx*delta_x, nx )

    for frame in range(0, n_frames, int(n_frames/n_plots) ):
        fig_slice, ax_slice = plt.subplots()
        ax_slice.plot( xvalues, EM[frame,int(ny/2), : ] )
        ax_slice.plot( xvalues, np.mean( EM[frame,:,:], axis=0 ) )
        plt.ylabel( fname[-4:-1] )
        plt.xlabel("x (cm)")
        time = frame*dt*save_rate
        figname = ("./figures/" + fname + grid_quantity 
            + "_slice_f" + str(frame) + "t" + str(time) + ".png")
        fig_slice.savefig( figname )


# class scatter_mp4():
#     def __init__(self, grid, extent, fname):
#         self.fig, self.ax = plt.subplots()
#         plt.xlabel("x, (cm)")
#         plt.ylabel("y, (cm)")
#         plt.title( fname[-4:-1] )
#         self.grid = grid
#         self.extent = extent
#     
#     def ani_init(self ):
#         self.im = self.ax.imshow( self.grid[0,:,:], extent=self.extent, vmax=np.max(self.grid), aspect='auto' )
#         self.cbar = self.fig.colorbar(self.im, ax = self.ax)
#         self.im.set_clim( np.min(self.grid), np.max(self.grid) )
#         plt.tight_layout()
# 
#     def ani_update(self, frame):
#         self.im.set_data( self.grid[frame,:,:] )
#         self.im.set_clim( np.min(self.grid[frame,:,:]), np.max(self.grid[frame,:,:]) )
#     
#     def animate(self, n_frames):
#         self.ani = anim.FuncAnimation( self.fig, self.ani_update, init_func=self.ani_init,
#             frames=range(1, n_frames ), repeat=False, blit=False)
#     
#     def ani_save(self, fname, grid_quantity):
#         self.ani.save("./figures/"+ fname + grid_quantity +"_evolution.mp4", fps=5 )


def fit_e_momentum( mydict, data_dir, fname ):

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

    time = np.linspace(0, tmax, n_frames)

    electron_vel = np.fromfile( data_dir + fname + "/particle_e_velocities", dtype="double", count=-1 ) 
    electron_vel = np.reshape(electron_vel, ( n_frames, 3*n_particles ) )

    tot_momentum = np.sum( electron_vel[:,0::3], axis=1 ) / np.sum( electron_vel[0,0::3] )

    results1 = fit_sine.fit_sin( time, tot_momentum )


    bestguess_omega1 = results1["omega"]

    fitfunc1 =results1["fitfunc"] 
    funkvals1 = fitfunc1( time )

    plt.plot( time, tot_momentum )
    plt.plot(time, funkvals1 )
    best_guess_legend = "Fitted sine \omega = "  + "{:2.2e}".format(bestguess_omega1)
    plt.legen( ["Data", best_guess_legend] )
    plt.savefig( "./figures/" + fname + "/Plasma_freq" + ".png" )
    plt.close()


def plot_velocities( mydict, data_dir, fname, n_trajs, traj_len_fraction ):
    
    fig_vel, ax_vel = plt.subplots()

    n_particles = int(mydict['n_particles'])
    n_tsteps = int(mydict['n_tsteps'])
    save_rate = int(mydict['save_rate'])
    n_frames = int(  n_tsteps / save_rate ) + 1

    tmax = float(mydict['tmax'])
    
    electron_vel = np.fromfile( data_dir + fname + "/particle_e_velocities", dtype="double", count=-1 ) 
    electron_vel = np.reshape(electron_vel, ( n_frames, 3*n_particles ) )

    i_particle = (np.random.randint(0, n_particles-1, size=n_trajs))*3

    elec_x_vels = np.transpose( electron_vel[:, i_particle]   )
    elec_y_vels = np.transpose( electron_vel[:, i_particle+1] )


    time = np.linspace(0, tmax, n_frames)
    for vel_x in elec_x_vels:
         ax_vel.plot( time,
              vel_x, '.b' )
    
    plt.title( "Electron velocity along x-axis vs time")
    plt.xlabel( "time (s)" )
    plt.ylabel( "x-velocity (cm/s)" )
    figure_name = "figures/"  + fname +  "/electron_x_vel_v_time" + ".png"
    fig_vel.savefig(figure_name)
    plt.close()
    
    fig_vel, ax_vel = plt.subplots()
    for vel_y in elec_y_vels:
         ax_vel.plot( time,
              vel_y, '.r' )
    
    plt.title( "Electron velocity along y-axis vs time")
    plt.xlabel( "time (s)" )
    plt.ylabel( "y-velocity (cm/s)" )
    figure_name = "figures/"  + fname +  "/electron_y_vel_v_time" + ".png"
    fig_vel.savefig(figure_name)
    plt.close()