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
    100000,    // weight     
    model,       // use_RFD    
    10000,    // n_tsteps  
    100,      // save_rate 
          
    2e-2,    // plasma_wavelen
         
    -2e-2,   // x_min         
    2e-2,    // x_max         
                   
    1e3,     // Te, in eV            
                 
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
        
        if( use_RFD==1 ) {
            dt = dt*1/25;
            n_tsteps = 2000;
            save_rate = 50;
        }
        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc
        double num_megabytes = (n_tsteps / save_rate * (nx * ny * 8.0 * 8.0 + n_particles * 2.0 * 12.0 * 8.0)) / 1e6;
        printf("Simulation will require %lf megabytes of harddrive space! \n", num_megabytes);
        print_primitives();
        print_derived_quantities();
    
    


    }

    double sine_density_perturb( double x, double y ) {
        
        double uniform_frac = 0.6;
        double sine_frac = 0.1;
        int n_cells = nx * ny;
        double macrop_per_cell = n_particles / (double ) n_cells;

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
                int p_cell = std::round(std::floor(macrop) + ( global_random() < (macrop - std::floor(macrop))));
                std::cout << p_cell << "\n";

                bool unassigned_part = assigned_particles < n_particles - p_cell;
                for( int iparticle = 0; ( iparticle < p_cell ) && unassigned_part; iparticle++ ) {
                    positions[assigned_particles * 3] =  global_random() * delta_x - delta_x/2 + x;
                    positions[assigned_particles * 3 + 1] = global_random() * delta_y - delta_y/2 + y;
                    assigned_particles = assigned_particles + 1;
                }
            }
        }
        while( assigned_particles < n_particles ) {
            positions[assigned_particles * 3] =  global_random() * delta_x*nx;
            positions[assigned_particles * 3 + 1] = global_random() * delta_y*ny;
            assigned_particles = assigned_particles + 1;
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

        for (int iy = 0; iy < ny; iy++)
        {
          for (int ix = 0; ix < nx; ix++)
          {
            double x = ix * delta_x;
            double y = iy * delta_y;
            double rho_q = n_particles * (double) weight / ( nx*ny*delta_x*delta_x);
            EM_ic.E_x[ix + iy * nx] = 8.0 * PI * q_e_cgs/wave1_k * 0.1  * rho_q * std::cos( wave1_k * x);

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