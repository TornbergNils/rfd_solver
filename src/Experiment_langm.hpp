#ifndef IC_LANGM
#define IC_LANGM

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
class Experiment_langm : public IC_struct
{
public:
    // Physical constants
    

    Experiment_langm() : IC_struct(
    
    50000,   // n_particles
    512,     // nx         
    2,       // ny         
    8000,    // weight     
    0,       // use_RFD    
    20000,    // n_tsteps  
    200,      // save_rate 
          
    2e-4,    // plasma_wavelen
         
    -1e-4,   // x_min         
    1e-4,    // x_max         
                   
    1e3,     // Te, in eV            
                 
    0.0,     // wave1_A       
    0.0,     // wave1_k       
    0.0,     // wave2_A       
    0.0,     // wave2_k       
    0.0,     // Ex_A          
    0.0     // Ex_k          
       ) {

        Generate_electron_positions();
        Generate_positron_positions();
        Generate_electron_velocities();
        Generate_positron_velocities();
        Set_EM_field();

        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc
        print_primitives();
    
    


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

            
            double v1 = 0.1 * electron_momentum * std::sin(wavevector * e_pos_ic[ip]);
            e_vel_ic[ip] = vel_and_gamma[0] + v1;
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
            
            double v1 = -0.1 * electron_momentum * std::sin(wavevector * p_pos_ic[ip]);

            p_vel_ic[ip] = vel_and_gamma[0] + v1;
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
                // printf( "(%.2lf, %.2lf)", x, y );

                EM_ic.E_x[ix + iy * nx] = 0.0;
                EM_ic.E_y[ix + iy * nx] = Ex_A * std::sin(Ex_k * x);
                EM_ic.E_z[ix + iy * nx] = 0.0; //4e4; // 2000; // Ex_A*std::cos( Ex_k * x );

                EM_ic.B_x[ix + iy * nx] = 0;
                EM_ic.B_y[ix + iy * nx] = 0.0; // Ex_A*std::sin( Ex_k * x );
                EM_ic.B_z[ix + iy * nx] = Ex_A * std::sin(Ex_k * x);

                /*
                EM_IC.E_x[ix + iy * nx] += std::copysign(1.0,x-nx*delta_x/2)
                  * Gaussian( x-nx*delta_x/2, y-ny*delta_y/2, wave1_A*0.2, 10*delta_x );

                EM_IC.E_y[ix + iy * nx] += std::copysign(1.0,y-ny*delta_y/2)
                  * Gaussian( x-nx*delta_x/2, y-ny*delta_y/2, wave1_A*0.2, 10*delta_x );
                */

                EM_ic.E_x[ix + iy * nx] += Get_EM_wave_component(0, config, x, y, 0);
                EM_ic.E_y[ix + iy * nx] += Get_EM_wave_component(1, config, x, y, 0);
                EM_ic.E_z[ix + iy * nx] += Get_EM_wave_component(2, config, x, y, 0);
                // printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

                EM_ic.B_x[ix + iy * nx] += Get_EM_wave_component(3, config, x, y, 0);
                EM_ic.B_y[ix + iy * nx] += Get_EM_wave_component(4, config, x, y, 0);
                EM_ic.B_z[ix + iy * nx] += Get_EM_wave_component(5, config, x, y, 0);
            }
        }

        
        wave_config_init[0][0] = wave2_A; // wave2_A;
        wave_config_init[0][1] = wave1_k;
        wave_config_init[0][2] = -PI;
        wave_config_init[0][3] = 0.0;
        EM_wave_config config2(wave_config_init);

        for (int iy = 0; iy < ny; iy++)
        {
          for (int ix = nx*2/3; ix < nx; ix++)
          {
            double x = ix * delta_x;
            double y = iy * delta_y;
            EM_ic.E_x[ix + iy * nx] += Get_EM_wave_component(0, config2, x, y, 0);
            EM_ic.E_y[ix + iy * nx] += Get_EM_wave_component(1, config2, x, y, 0);
            EM_ic.E_z[ix + iy * nx] += Get_EM_wave_component(2, config2, x, y, 0);
            //printf( "%lf, ", EM_ic.E_z[ix + iy * nx] );

            EM_ic.B_x[ix + iy * nx] += Get_EM_wave_component(3, config2, x, y, 0);
            EM_ic.B_y[ix + iy * nx] += Get_EM_wave_component(4, config2, x, y, 0);
            EM_ic.B_z[ix + iy * nx] += Get_EM_wave_component(5, config2, x, y, 0);

          // printf( "\n" );
        }
        }
        
        
    }
};

#endif //IC_LANGM