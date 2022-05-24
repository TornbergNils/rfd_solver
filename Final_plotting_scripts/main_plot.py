import numpy as np
import sys
import gc
import csv
import matplotlib.pyplot as plt
import funcs_for_plotting


## Plotting settings

data_dir = "./data/"
print("Args are:")
print(sys.argv) 
if( len(sys.argv) > 1):
    experiments = sys.argv[1:]
else:
    data_dir = "./data"
    experiments = [""]
# Number of equally spaced snapshots of the
# particle and charge densities to plot
n_plots = 5
# What fraction of the entire particle trajectories to plot
traj_len_fraction = 0.5
# Number of particle trajectories to plot
n_trajs = 8


## For syncing colorbars across rho
#rho_max_list = list()
#rho_min_list = list()
#rho_im_list = list()

def plot_everything( experiment ):
    with open( data_dir + experiment + '/config.csv', mode='r') as infile:
        reader = csv.reader(infile)
        mydict = {rows[0]:rows[1] for rows in reader}
    
    
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

    #funcs_for_plotting.plot_qdensity( mydict, n_frames, n_plots, data_dir, 
    #    experiment,  nx, ny )
    #plt.close('all')

    #funcs_for_plotting.plot_positron_density( mydict, n_plots, data_dir,
    #    experiment )
    #plt.close('all')
    
    #funcs_for_plotting.plot_mean_density_v_time( mydict, data_dir, experiment)
    #plt.close('all')

    #funcs_for_plotting.plot_slice_density_along_middle( mydict, data_dir,
    #    experiment, n_plots)
    #plt.close('all')
    
    ## Call garbage collector to avoid too much memory use
    #gc.collect()

    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EME_x") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EME_y") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EME_z") 
    #plt.close('all')
    ## Call garbage collector to avoid too much memory use
    #gc.collect()
    
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EMB_x") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EMB_y") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/EMB_z") 
    #plt.close('all')
    ## Call garbage collector to avoid too much memory use
    #gc.collect()

    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/J_x") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/J_y") 
    #plt.close('all')
    #funcs_for_plotting.plot_grid_mp4(mydict, data_dir, experiment, "/J_z") 
    #plt.close('all')
    ## Call garbage collector to avoid too much memory use
    #gc.collect()
    
    #funcs_for_plotting.plot_energy_density_mp4(mydict, data_dir, experiment ) 
    
    #funcs_for_plotting.plot_trajectories( mydict, data_dir, experiment, n_trajs, traj_len_fraction )
    #funcs_for_plotting.plot_velocities( mydict, data_dir, experiment, n_trajs, traj_len_fraction )
    #plt.close('all')
    #gc.collect()

    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EME_x" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EME_y" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EME_z" )
    #plt.close('all')
    #gc.collect()
    
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EMB_x" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EMB_y" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/EMB_z" )
    
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/J_x" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/J_y" )
    #funcs_for_plotting.plot_grid_evolution( mydict, data_dir, experiment, n_plots, "/J_z" )

    #funcs_for_plotting.plot_pos_v_time( mydict, data_dir, experiment, n_trajs, traj_len_fraction )
    #plt.close('all')
    #gc.collect()
    
    funcs_for_plotting.plot_multi_movie_mp4(mydict, data_dir, experiment ) 
    
    funcs_for_plotting.plot_EM_slice(mydict, data_dir, experiment, "/EMB_x", n_plots )
    funcs_for_plotting.fit_e_momentum( mydict, data_dir, experiment )



for experiment in experiments:
    plot_everything( experiment )
