import numpy as np

## Constants
PI = np.pi
k_boltz_ev = (8.617 * 1e-5)
k_boltz_erg = 1.38064 * 1e-16
c_cgs = 2.99792458 * 1e10
m_e_cgs = 9.1093819*1e-28
q_e_cgs = 4.80320425e-10

######## Input param

max_density = 0.005 # units of crit density 

plasma_wavelen = 1e-4 # cm
plasma_period = (plasma_wavelen/c_cgs) / np.sqrt(max_density)
Wave_amplitude = 0.1

nx = 1024  # matrixSize_x
x_min = -2e-4 # cm
x_max = 2e-4 # cm

steps_per_plasma_period = 2**6 
Te = 1000 # eV
T_kelvin = Te / k_boltz_ev # kelvin
print( "T_kelvin = ", "{:2.2e}".format(T_kelvin) )
thermal_velocity = np.sqrt( k_boltz_erg * T_kelvin / m_e_cgs )
print( "thermal_velocity = ", "{:2.2e}".format(thermal_velocity) )




######## Derived param 

W0 = 2 * PI * c_cgs / plasma_wavelen # in HZ ??
print("W0 = ", "{:2.2e}".format(W0) )

crit_density = ( W0**2 ) * m_e_cgs / ( 4 * PI * q_e_cgs**2 ) # n / cm^3
print("crit_density = ", "{:2.2e}".format(crit_density ) )

E0 = m_e_cgs * c_cgs * W0 / ( - q_e_cgs ) # CGSE
print("E0 = ", "{:2.2e}".format(E0 ) )

N_cr = crit_density
Ne = max_density * crit_density # In cm^3
print("Ne = ", "{:2.2e}".format(Ne ) )

Emax = 2*plasma_wavelen*Wave_amplitude*Ne*(-q_e_cgs)
print("Emax = ", "{:2.2e}".format(Emax ) )
pmax = (Emax*(-q_e_cgs) * plasma_period/(2*PI)) / ( m_e_cgs * c_cgs)

# particle_concen = Ne*(1 + Wave_amplitude*np.cos( 2 * PI * X / plasma_wavelen))
# particle pxc = pmax * cos( 2 * PI * X/ plasma_wavelen )

dx = (x_max - x_min ) / nx
print("dx = ", "{:2.2e}".format(dx ) )
L_debye = np.sqrt( Te * 1.6e-12 /( 4 * PI * q_e_cgs**2*Ne))
print("L_debye = ", "{:2.2e}".format(L_debye ) )

print("Te * 1.6e-12 = ", Te * 1.6e-12  )
print("Te converted to vel", "{:2.2e}".format(np.sqrt(Te *c_cgs**2/0.511e6))  )


#est_EM_wave_Freq = 

Wp = np.sqrt(4*PI*q_e_cgs**2*Ne/m_e_cgs)
print("Wp", "{:2.2e}".format(Wp)  )

EM_wave_freq = np.sqrt( c_cgs**2 * (2*PI/plasma_wavelen)**2 + Wp**2 )
print("EM_wave_freq", "{:2.2e}".format(EM_wave_freq)  )

est_E_from_v_deviation = 4*PI*q_e_cgs*Ne*0.1*0.047*c_cgs / Wp
print("my estimated Emax", "{:2.2e}".format(est_E_from_v_deviation)  )