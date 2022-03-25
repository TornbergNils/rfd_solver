import numpy as np

nx = 64
ny = 64

delta_x = 0.02
delta_y = 0.02

x_len = nx * delta_x
y_len = ny * delta_y
n_particles = 50000

print("particle per cell ", n_particles / (nx*ny) )

PI = 3.14159265358979

# cgs units
c_cgs = 2.99792458 * 1e10
q_e_cgs = 1.70269007*1e-9
m_e_cgs = 9.1093819*1e-28
boltzmann_cgs = 1.380658*1e-16
v_thermal_cgs = 0.13 * c_cgs



c = 1.0
v_thermal = 0.05 * c
q_e = 0.30282212088
m_e = 511 # in KeV

# Conversion factors
len_KeV_to_cm = 1.9732705*1e-8
sq_len_to_cm2 = len_KeV_to_cm * len_KeV_to_cm
nat_KeV_to_Hz = 1.519*1e18



density_nat = n_particles / (nx * ny * delta_x * delta_y)
density_cgs = n_particles / (nx * ny * delta_x * delta_y * sq_len_to_cm2 )
print( "density nat:  {:.2e}".format(density_nat) + " \n")
print( "density cgs:  {:.2e}".format(density_cgs) + " \n")

plasma_freq_cgs = np.sqrt(  density_cgs * q_e_cgs * q_e_cgs / ( m_e_cgs ) )
plasma_freq_nat = np.sqrt( density_nat * q_e * q_e / m_e  )
plasma_freq_converted = 2*PI * plasma_freq_cgs / nat_KeV_to_Hz

print( "plasma freq cgs: {:.2e}".format(plasma_freq_cgs))
print( "plasma freq cnv: {:.2e}".format(plasma_freq_converted))
print( "plasma freq nat: {:.2e}".format(plasma_freq_nat))
print( "nat d converted: {:.2f}".format( plasma_freq_nat / plasma_freq_converted ))

temperature_kelvin = v_thermal_cgs**2 * m_e_cgs / boltzmann_cgs
debye_length = v_thermal_cgs / plasma_freq_cgs

print( "Electron temperature: {:.2e}".format( temperature_kelvin ))
print( "Debye length: {:.2e}".format( debye_length ))

"""
printf("plasma freq nat 2 %2.2e \n", plasma_freq_nat2)
debye_length = v_thermal / (stdsqrt(2) * plasma_freq_nat )


debye_length_cgs = (m_e_cgs * (v_thermal * c_cgs )* (v_thermal * c_cgs) / 2  ) / ( density_cgs * q_e_cgs * q_e_cgs)

printf("debye length cgs = %2.2e \n", debye_length_cgs )
debye_test = debye_length_cgs / len_KeV_to_cm
debye_test2 = stdsqrt( ( v_thermal * m_e ) / ( density * q_e * q_e ) )
print( "debye_test = %2.2e \n", debye_test )
print( "debye_test2 = %2.2e \n", debye_test2 )
nat_temp = v_thermal*v_thermal*m_e
k_boltz_Kev = 8.617333262*1e-8
kelvin_temp = nat_temp / k_boltz_Kev
"""