import numpy as np

file = open( "./parameters.txt", 'w')
file.close()

file = open( "./parameters.txt", 'a')
## Constants
PI = np.pi
k_boltz_ev = (8.617 * 1e-5)
k_boltz_erg = 1.38064 * 1e-16

c_cgs = 2.99792458 * 1e10
m_e_cgs = 9.1093819*1e-28
q_e_cgs = 4.80320425e-10



######## Input param
set_wave_ic = 1
positrons_enabled = 1
use_RFD = 0
n_species = 2 if positrons_enabled==1 else 1


max_density = 0.005 # units of crit density 

plasma_wavelen = 2e-4 # cm
plasma_wavenum = 2*PI/plasma_wavelen
plasma_period = (plasma_wavelen/c_cgs) / np.sqrt(max_density)
Wave_amplitude = 0.1

nx = 108 # matrixSize_x
ny = 108
n_particles = 100000
weight = 4000

x_min = -1e-4 # cm
x_max = 1e-4 # cm

dx = (x_max - x_min ) / nx
dy = (x_max - x_min ) / nx

n_tsteps = 100
save_rate = 1
dt = dx / (2*c_cgs)
tmax = n_tsteps * dt


Te = 1e9 # eV
T_kelvin = Te / k_boltz_ev # kelvin
print( "T_kelvin = ", "{:2.2e}".format(T_kelvin) )
thermal_velocity = np.sqrt( k_boltz_erg * T_kelvin / m_e_cgs )
print( "thermal_velocity = ", "{:2.2e}".format(thermal_velocity) )

# NOTE: momentum would be T*kb*3/2 but we are in 2 dimensions
electron_momentum = np.sqrt( T_kelvin * k_boltz_erg  / ( m_e_cgs ) )
print( "electron_momentum / m_s = ", "{:2.2e}".format(electron_momentum) )

######## Derived param 

W0 = 2 * PI * c_cgs / plasma_wavelen # in HZ ??
print("W0 = ", "{:2.2e}".format(W0) )

crit_density = ( W0**2 ) * m_e_cgs / ( 4 * PI * q_e_cgs**2 ) # n / cm^3
#print("crit_density = ", "{:2.2e}".format(crit_density ) )

E0 = m_e_cgs * c_cgs * W0 / ( - q_e_cgs ) # CGSE
#print("E0 = ", "{:2.2e}".format(E0 ) )

N_cr = crit_density
print("crit density = ", "{:2.2e}".format(crit_density ) )
Ne = n_species * max_density * crit_density # In cm^3
print("Desired density = ", "{:2.2e}".format(Ne ) )
actual_density = n_species * n_particles * weight / ( nx * ny * dx * dy )
print("Actual density = ", "{:2.2e}".format(actual_density ) )

z_velocity = 0.1*c_cgs
abstatampere_to_ampere = 10/c_cgs
expected_total_current = n_particles*weight*q_e_cgs*z_velocity * abstatampere_to_ampere
print("expected_total_current = ", "{:2.2e}".format(expected_total_current ) )

#Emax = 2*plasma_wavelen*Wave_amplitude*Ne*(-q_e_cgs)
#print("Emax = ", "{:2.2e}".format(Emax ) )
#pmax = (Emax*(-q_e_cgs) * plasma_period/(2*PI)) / ( m_e_cgs * c_cgs)

# particle_concen = Ne*(1 + Wave_amplitude*np.cos( 2 * PI * X / plasma_wavelen))
# particle pxc = pmax * cos( 2 * PI * X/ plasma_wavelen )

L_debye = np.sqrt( Te * 1.6e-12 /( 4 * PI * q_e_cgs**2*Ne))
L_debye_actual = np.sqrt( Te * 1.6e-12 /( 4 * PI * q_e_cgs**2*actual_density))
print("Actual debye length = ", "{:2.2e}".format(L_debye_actual ) )
print("dx = ", "{:2.2e}".format(dx ) )
print("dy = ", "{:2.2e}".format(dy ) )

#print("Te * 1.6e-12 = ", Te * 1.6e-12  )
print("Te converted to vel", "{:2.2e}".format(np.sqrt(Te *c_cgs**2/0.511e6))  )


#est_EM_wave_Freq = 

Wp = np.sqrt(4*PI*q_e_cgs**2*Ne/m_e_cgs)
actual_plasma_freq = np.sqrt(4*PI*q_e_cgs**2*actual_density/m_e_cgs)
print("actual_plasma_freq", "{:2.2e}".format(actual_plasma_freq)  )
dispersion_rel_freq = actual_plasma_freq*np.sqrt( 1 + 3*L_debye_actual**2 * plasma_wavenum**2)
print("dispersion_rel_freq", "{:2.2e}".format(dispersion_rel_freq)  )
print("1/dt ", "{:2.2e}".format( 1 / dt) )

EM_wave_freq = np.sqrt( c_cgs**2 * (2*PI/plasma_wavelen)**2 + Wp**2 )
print("EM_wave_freq", "{:2.2e}".format(EM_wave_freq)  )
EM_crit_density = m_e_cgs / (q_e_cgs**2)*EM_wave_freq**2
print("EM_wave_crit_density", "{:2.2e}".format(EM_crit_density)  )

est_E_from_v_deviation = 4*PI*q_e_cgs*Ne*Wave_amplitude*thermal_velocity / actual_plasma_freq
print("my estimated Emax", "{:2.2e}".format(est_E_from_v_deviation)  )


## Wave ic settings

wave1_amplitude =  1e14 #est_E_from_v_deviation
wave1_wavevect = 9*plasma_wavenum
wave1_freq = c_cgs*wave1_wavevect

wave2_amplitude = 0.0
wave2_wavevect = 0.0 #PI
wave2_freq = c_cgs*wave2_wavevect
Ex_raw = 0.0; #10000.0 # est_E_from_v_deviation
Ex_wavevect = plasma_wavenum # plasma_wavenum

file.write( "set_wave_ic " + str(set_wave_ic) + "\n" )
file.write( "use_RFD " + str(use_RFD) + "\n" )
file.write( "wave1_amplitude " + str(wave1_amplitude) + "\n" )
file.write( "wave1_wavevect " + str(wave1_wavevect) + "\n" )
file.write( "wave2_amplitude " + str(wave2_amplitude) + "\n" )
file.write( "wave2_wavevect " + str(wave2_wavevect) + "\n" )
file.write( "Ex_raw " + str(Ex_raw) + "\n" )
file.write( "Ex_wavevect " + str(Ex_wavevect) + "\n" )

file.write( "nx " + str(nx) + "\n" )
file.write( "ny " + str(ny) + "\n" )
file.write( "dx " + str(dx) + "\n" )
file.write( "dy " + str(dy) + "\n" )
file.write( "weight " + str(weight) + "\n" )
file.write( "n_particles " + str(n_particles) + "\n" )
file.write( "n_species " + str(n_species)   + "\n" )

file.write( "n_tsteps " + str(n_tsteps) + "\n" )
file.write( "dt " + str(dt) + "\n" )
file.write( "save_rate " + str(save_rate) + "\n" )
file.write( "tmax " + str(tmax) + "\n" )

file.write( "c " +   str(c_cgs)   + "\n" )
file.write( "m_e " + str(m_e_cgs) + "\n" )
file.write( "q_e " + str(q_e_cgs) + "\n" )

file.write( "v_thermal " + str(thermal_velocity) + "\n" )
file.write( "electron_momentum " + str(electron_momentum) + "\n" )
file.write( "plasma_wavenum " +   str(plasma_wavenum)   + "\n" )
file.write( "plasma_wavelen " +   str(plasma_wavelen)   + "\n" )
file.write( "plasma_amplitude " + str(Wave_amplitude)   + "\n" )

file.write( "plasma_freq " + "{:2.2e}".format(actual_plasma_freq)   + "\n" )
file.write( "density " + "{:2.2e}".format(actual_density)   + "\n" )


file.close()

print()
print("########################################")
print()