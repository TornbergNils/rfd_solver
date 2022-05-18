#ifndef IC_EZ_PINCH
#define IC_EZ_PINCH

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
class Experiment_Ez_pinch : public IC_struct
{
public:
    // Physical constants
    

    Experiment_Ez_pinch( int model ) : IC_struct(
    
    100000,   // n_particles
    128,     // nx         
    128,       // ny         
    4000,    // weight     
    model,       // use_RFD    
    560,    // n_tsteps  
    4,      // save_rate 
          
    2e-4,    // plasma_wavelen
         
    -1e-4,   // x_min         
    1e-4,    // x_max         
                   
    1e4,     // Te, in eV            
                 
    0.0,     // wave1_A       
    0.0,     // wave1_k       
    0.0,     // wave2_A       
    0.0,     // wave2_k       
    1e3,     // Ex_A          
    0.0     // Ex_k          
       ) {

        
        if( use_RFD==1 ) {
            dt = dt*1/25;
            n_tsteps = 2000;
            save_rate = 20;
        }

        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc to text file
        double num_megabytes = (n_tsteps / save_rate * (nx * ny * 8.0 * 8.0 + n_particles * 2.0 * 12.0 * 8.0)) / 1e6;
        printf("Simulation will require %lf megabytes of harddrive space! \n", num_megabytes);
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
            
            
            double radius = x_len * 3 / 8 * std::sqrt(global_random());
            double theta = global_random() * 2 * PI;

            e_pos_ic[ip] = x_len / 2 + radius * std::cos(theta); /// 2 + x_len/2;
            e_pos_ic[ip + 1] = y_len / 2 + radius * std::sin(theta);
            e_pos_ic[ip + 2] = 0.0;
            
        }
    }

    void Generate_positron_positions()
    {
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;
        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            double radius = x_len * 3 / 8 * std::sqrt(global_random());
            double theta = global_random() * 2 * PI;

            p_pos_ic[ip] = x_len / 2 + radius * std::cos(theta); /// 2 + x_len/2;
            p_pos_ic[ip + 1] = y_len / 2 + radius * std::sin(theta);
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

        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * delta_x;
                double y = iy * delta_y;

                EM_ic.E_x[ix + iy * nx] += 0.000;
                EM_ic.E_z[ix + iy * nx] += Ex_A;
            }
        }
        
        
    }
};

#endif //IC_EZ_PINCH