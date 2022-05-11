#ifndef IC_WAVE
#define IC_WAVE

#include <vector>
#include <random>
#include <map>
#include <chrono>
#include "classes.hpp"
#include "propagation.hpp"
#include "generate_IC.hpp"


/*
    Struct containing initial conditions for simulation. 

    Use by creating object with zero initialized vectors and
    then run the member functions to get the boundary conditions
    as defined in the member functions.
*/
class Experiment_wave : public IC_struct
{
public:
    // Physical constants
    

    Experiment_wave() : IC_struct(
    
    100000,   // n_particles
    216,     // nx         
    54,       // ny         
    4000,    // weight     
    0,       // use_RFD    
    280,    // n_tsteps  
    2,      // save_rate 
          
    2e-4,    // plasma_wavelen
         
    -1e-4,   // x_min         
    1e-4,    // x_max         
                   
    1e9,     // Te, in eV            
                 
    0.0,     // wave1_A       
    0.0,     // wave1_k       
    0.0,     // wave2_A       
    0.0,     // wave2_k       
    0.0,     // Ex_A          
    0.0     // Ex_k          
       ) {

        
        wave1_A = 1e14;
        wave1_k = 2*plasma_wavenum;
        if( use_RFD==1 ) {
            dt = dt*1/25;
            n_tsteps = 7000;
            save_rate = 50;
        }

        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc
        print_primitives();
        print_derived_quantities();
    
        Generate_electron_positions();
        Generate_positron_positions();
        Generate_electron_velocities();
        Generate_positron_velocities();
        Set_EM_field();
    


    }

    void Generate_electron_positions()
    {

        double x_len = nx * delta_x;
        double y_len = ny * delta_y;

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            
            e_pos_ic[ip] = global_random() * x_len; ///3 + x_len/3;
            e_pos_ic[ip + 1] = global_random() * y_len;
            e_pos_ic[ip + 2] = 0.0;
            
        }
    }

    void Generate_positron_positions()
    {
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;
        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            p_pos_ic[ip] = global_random() * x_len; ///3 + x_len/3;
            p_pos_ic[ip + 1] = global_random() * y_len;
            p_pos_ic[ip + 2] = 0.0;
            
        }
    }

    void Generate_electron_velocities()
    {

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            std::vector<double> vel_and_gamma = Get_relativistic_vel_and_gamma(electron_momentum, PI, c);

            
            e_vel_ic[ip] = vel_and_gamma[0];
            e_vel_ic[ip + 1] = vel_and_gamma[1];
            e_vel_ic[ip + 2] = 0.0;

            double vel_squared_e = e_vel_ic[ip] * e_vel_ic[ip] 
                + e_vel_ic[ip + 1] * e_vel_ic[ip + 1] 
                + e_vel_ic[ip + 2] * e_vel_ic[ip + 2];
            
            //double gamma_e = 1.0 / std::sqrt(1.0 - vel_squared_e / (c * c));
            e_gamma_ic[ip / 3] = vel_and_gamma[3];

            if (vel_squared_e > c * c)
            {
                printf("e particle speed exceeds c! c = %lf \n", c);
            }
        }
    }

    void Generate_positron_velocities()
    {

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            // vel and gamma contains vx, vy, vz, gamma
            std::vector<double> vel_and_gamma = Get_relativistic_vel_and_gamma(electron_momentum, PI, c);
            

            p_vel_ic[ip] = vel_and_gamma[0];
            p_vel_ic[ip + 1] = vel_and_gamma[1]; //-0.5 - 0.1 * std::sin( positron_pos[ip+1] / 100 );
            p_vel_ic[ip + 2] = 0.0;

            
            double vel_squared_p = p_vel_ic[ip] * p_vel_ic[ip] 
            + p_vel_ic[ip + 1] * p_vel_ic[ip + 1] 
            + p_vel_ic[ip + 2] * p_vel_ic[ip + 2];
            /*
            double gamma_p = 1.0 / std::sqrt(1.0 - vel_squared_p / (c * c));
            */

            p_gamma_ic[ip / 3] = vel_and_gamma[3];
            if (vel_squared_p > c * c)
            {
                printf("p particle speed exceeds c! \n");
            }
        }
    }

    void Set_EM_field()
    {

        printf("Ex_A %lf \n", Ex_A);
        printf("delta_x %2.2e \n", delta_x);

        std::vector<double>::size_type num_waves = 2;
        std::vector<std::vector<double>> wave_config_init{num_waves,
                                                          std::vector<double>(4)};

        wave_config_init[0][0] = wave1_A; // amplitude
        wave_config_init[0][1] = wave1_k; // wavevect = ang_freq, ok for c=1 or t=0
        wave_config_init[0][2] = 0.0;     // prop angle
        wave_config_init[0][3] = 0.0;     // phase

        wave_config_init[1][0] = 0.0; // wave2_A;
        wave_config_init[1][1] = wave2_k;
        wave_config_init[1][2] = 0.0;
        wave_config_init[1][3] = 0.0;
        EM_wave_config config(wave_config_init);

        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * delta_x;
                double y = iy * delta_y;

                EM_ic.E_x[ix + iy * nx] += Get_EM_wave_component(0, config, x, y, 0);
                EM_ic.E_y[ix + iy * nx] += Get_EM_wave_component(1, config, x, y, 0);
                EM_ic.E_z[ix + iy * nx] += Get_EM_wave_component(2, config, x, y, 0);

                EM_ic.B_x[ix + iy * nx] += 0.0; //Get_EM_wave_component(3, config, x, y, 0);
                EM_ic.B_y[ix + iy * nx] += Get_EM_wave_component(4, config, x, y, 0);
                EM_ic.B_z[ix + iy * nx] += Get_EM_wave_component(5, config, x, y, 0);
            }
        }
        
        
    }
};

#endif //IC_WAVE