import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl
import scipy.stats as stats
import math


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

vel_ic = np.fromfile( "./data/initial_velocities", dtype="double", count=-1 )

x_vels = vel_ic[0::3]
y_vels = vel_ic[1::3]

print( y_vels )


c_cgs = 2.99792458 * 1e10
v_thermal = 0.05

mu = 0
variance = v_thermal * c_cgs
sigma = math.sqrt(variance)
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
plt.hist( y_vels, bins=50, density=True )
#plt.plot(x, stats.norm.pdf(x, mu, sigma))
plt.legend( ["y_vels", "ref_pdf" ] )
plt.savefig("./figures/y_vels_vs_gaussian.png" )
plt.close()

plt.hist( x_vels, bins=50, density=True )
#plt.plot(x, stats.norm.pdf(x, mu, sigma))
plt.legend( ["y_vels", "ref_pdf" ] )
plt.savefig("./figures/x_vels_vs_gaussian.png" )
plt.close()