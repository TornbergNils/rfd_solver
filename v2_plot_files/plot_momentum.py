import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl
import scipy.stats as stats
import math
import Fit_sine


mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.cmap'] = 'Spectral'

with open('config.csv', mode='r') as infile:
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

time = np.linspace(0, tmax, n_frames-1)

vel_ic = np.fromfile( "./data/e_momenta", dtype="double", count=-1 )

vel_ic = np.reshape( vel_ic, ((n_frames-1), n_particles * 3 ))

tot_momentum = np.sum( vel_ic[:,0::3], axis=1 ) / np.sum( vel_ic[0,0::3] )

print( "initial avg momentum: ", np.mean( vel_ic[0,0::3]) )
print( "final avg momentum: ", np.mean( vel_ic[n_frames-2,0::3]) )


results1 = Fit_sine.fit_sin( time, tot_momentum )


print( "Best guess for p frequency is: ")
bestguess_omega1 = results1["omega"]

fitfunc1 =results1["fitfunc"] 
funkvals1 = fitfunc1( time )
print( "{:2.2e}".format(bestguess_omega1) )

plt.plot( time, tot_momentum )
plt.plot(time, funkvals1 )
plt.savefig("./figures/tot_x_momentum.png" )
plt.close()
