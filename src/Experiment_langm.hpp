#ifndef IC_LANGM
#define IC_LANGM

#include <vector>
#include <random>
#include <map>
#include <chrono>
#include "EM.hpp"
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
    

    Experiment_langm( int model ) : IC_struct(
    
    100000,   // n_particles
    512,     // nx         
    12,       // ny         
    5000000,    // weight     
    model,       // use_RFD    
    1000,    // n_tsteps  
    10,      // save_rate 
          
    1e-2,    // plasma_wavelen
         
    -2e-2,   // x_min         
    2e-2,    // x_max         
                   
    1e4,     // Te, in eV            
                 
    0.0,     // wave1_A       
    0.0,     // wave1_k       
    0.0,     // wave2_A       
    0.0,     // wave2_k       
    0.0,     // Ex_A          
    0.0     // Ex_k          
       ) {
        
        wave1_k = (2*PI)/plasma_wavelen;

        Generate_electron_positions();
        Generate_positron_positions();
        Generate_electron_velocities();
        Generate_positron_velocities();
        Set_EM_field();

        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc
        double num_megabytes = (n_tsteps / save_rate * (nx * ny * 8.0 * 8.0 + n_particles * 2.0 * 12.0 * 8.0)) / 1e6;
        printf("Simulation will require %lf megabytes of harddrive space! \n", num_megabytes);
        print_primitives();
        print_derived_quantities();
    
    


    }

    double sine_density_perturb( double x, double y ) {
        
        double uniform_frac = 0.6;
        double sine_frac = 0.05;
        int n_cells = nx * ny;
        double macrop_per_cell = n_particles / (n_cells);

        double macrop = macrop_per_cell;
        macrop = macrop + macrop_per_cell * sine_frac * std::sin( x * wave1_k );

        return macrop;
    }

    void smooth_start( std::vector<double> &positions ) {
        
        int assigned_particles = 0;
        for( int ix = 0; ix < nx; ix++ ) {
            for( int iy = 0; iy < ny; iy++ ) {
                double x = ix * delta_x + delta_x/2;
                double y = iy * delta_y + delta_y/2;

                double macrop = sine_density_perturb(x, y);
                std::cout << macrop << "\n";
                int p_cell = std::floor(macrop) + ( global_random() < (macrop - std::floor(macrop)));

                bool unassigned_part = assigned_particles < n_particles - p_cell;
                if( !unassigned_part ) {
                    p_cell = n_particles - assigned_particles - 1;
                }
                for( int iparticle = 0; ( iparticle < p_cell ); iparticle++ ) {
                    positions[assigned_particles * 3] =  global_random() * delta_x - delta_x/2 + x;
                    positions[assigned_particles * 3 + 1] = global_random() * delta_y - delta_y/2 + y;
                    assigned_particles = assigned_particles + 1;
                }
            }
        }
        std::cout << "Assigned fraction: " << assigned_particles/ (double) n_particles << "\n";
    }

    void Generate_electron_positions()
    {
        smooth_start( e_pos_ic );
    }  

    void Generate_positron_positions()
    {
        smooth_start( p_pos_ic );
    }

    void Generate_electron_velocities()
    {

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            std::vector<double> vel_and_gamma = Get_relativistic_vel_and_gamma(electron_momentum, PI, c);

            
            double v1 = 0.1 * electron_momentum * std::sin(wavevector * e_pos_ic[ip]);
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
            
            double v1 = -0.1 * electron_momentum * std::sin(wavevector * p_pos_ic[ip]);

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
          for (int ix = 0; ix < nx; ix++)
          {
            double x = ix * delta_x;
            double y = iy * delta_y;
            double rho_q = n_particles * (double) weight / ( nx*ny*delta_x*delta_x);
            EM_ic.E_x[ix + iy * nx] = 8.0 * PI * q_e_cgs/wave1_k * 0.05  * rho_q * std::cos( wave1_k * x);

            //EM_ic.E_x[ix + iy * nx] += Get_EM_wave_component(0, config2, x, y, 0);
            //EM_ic.E_y[ix + iy * nx] += Get_EM_wave_component(1, config2, x, y, 0);
            //EM_ic.E_z[ix + iy * nx] += Get_EM_wave_component(2, config2, x, y, 0);
            ////printf( "%lf, ", EM_ic.E_z[ix + iy * nx] );

            //EM_ic.B_x[ix + iy * nx] += Get_EM_wave_component(3, config2, x, y, 0);
            //EM_ic.B_y[ix + iy * nx] += Get_EM_wave_component(4, config2, x, y, 0);
            //EM_ic.B_z[ix + iy * nx] += Get_EM_wave_component(5, config2, x, y, 0);

          // printf( "\n" );
        }
        }
        
        
    }
};

#endif //IC_LANGM