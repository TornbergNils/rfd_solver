import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl
import Fit_sine
import GetFFt

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
extent = ( 0, nx*delta_x, 0, ny*delta_y)

# Load and format data
EME_x = np.fromfile( "./data/EME_x", dtype="double", count=-1 )

EME_x = np.reshape(EME_x, ( n_frames, ny, nx ) )

EME_x_at_pt1=EME_x[:, int(ny/4), int(nx/4) ]
EME_x_at_pt2=EME_x[:, int(2*ny/4), int(nx/4) ]
EME_x_at_pt3=EME_x[:, int(3*ny/4), int(2*nx/4) ]
EME_x_at_pt4=EME_x[:, int(2*ny/4), int(3*nx/4) ]


time = np.linspace(0, tmax, n_frames)

results1 = Fit_sine.fit_sin( time, EME_x_at_pt1 )
results2 = Fit_sine.fit_sin( time, EME_x_at_pt2 )
results3 = Fit_sine.fit_sin( time, EME_x_at_pt3 )
results4 = Fit_sine.fit_sin( time, EME_x_at_pt4 )
print( "Best guess for E frequency is: ")

bestguess_omega1 = results1["omega"]
bestguess_omega2 = results1["omega"]
bestguess_omega3 = results1["omega"]
bestguess_omega4 = results1["omega"]

print( "{:2.2e}".format(bestguess_omega1) )
print( "{:2.2e}".format(bestguess_omega2) )
print( "{:2.2e}".format(bestguess_omega3) )
print( "{:2.2e}".format(bestguess_omega4) )


fig1, ax1 = plt.subplots()


fitfunc1 =results1["fitfunc"] 
fitfunc2 =results2["fitfunc"] 
fitfunc3 =results3["fitfunc"] 
fitfunc4 =results4["fitfunc"] 

funkvals1 = fitfunc1( time )

plt.plot(time, funkvals1 )
plt.plot(time, EME_x_at_pt1 )
plt.legend(["guess", "real"])

fname = "Best guessess"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()


fft_results = GetFFt.get_fft( time, EME_x_at_pt1 )

plt.plot( fft_results["ff"], fft_results["Fyy"]  )
fname = "EMEx_fft"
plt.savefig( "./figures/" + fname + ".png" )
plt.close()