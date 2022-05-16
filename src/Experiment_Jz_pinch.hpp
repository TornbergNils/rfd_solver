#ifndef IC_JZ_PINCH
#define IC_JZ_PINCH

#include <vector>
#include <random>
#include <map>
#include <chrono>
#include <cmath>
#include "classes.hpp"
#include "propagation.hpp"
#include "generate_IC.hpp"


/*
    Struct containing initial conditions for simulation. 

    Use by creating object with zero initialized vectors and
    then run the member functions to get the boundary conditions
    as defined in the member functions.
*/
class Experiment_Jz_pinch : public IC_struct
{
public:
    double p_z_global;
    double elec_dist_radius;
    

    Experiment_Jz_pinch() : IC_struct(
    
    200000,   // n_particles
    250,     // nx         
    250,       // ny         
    1,    // weight     
    1,       // use_RFD    
    2000,    // n_tsteps  
    25,      // save_rate 
          
    2e-4,    // plasma_wavelen
         
    -5e-1,   // x_min         
    5e-1,    // x_max         
                   
    1e3,     // Te, in eV            
                 
    0.0,     // wave1_A       
    0.0,     // wave1_k       
    0.0,     // wave2_A       
    0.0,     // wave2_k       
    0.0,     // Ex_A          
    0.0     // Ex_k          
       ) {
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;
        p_z_global = 0.5 * c;
        elec_dist_radius = x_len * 1.5 / 8.0;
        
        // E_B_param governs ratio |E|/|B|, changing dynamics for RFD case
        double E_B_param = 0.02;

        double B_max = 4.0 * n_particles * weight
            * p_z_global * q_e_cgs / (c * elec_dist_radius );
        Ex_A = B_max * E_B_param;
        

        double Iz = 2.0 * n_particles * weight * p_z_global * q_e_cgs;
        double Iz_ampere = Iz * 10.0/c;
        double mu_0_SI = 1.25663706212*1e-6;
        double k_boltz_SI = 1.380649* 1e-23;
        // Bennett conditon
        double Bennet_dense = Iz_ampere*Iz_ampere*mu_0_SI/( 8.0*PI*k_boltz_SI*T_kelvin );
        // Convert line density per m to per cm
        Bennet_dense = Bennet_dense / 1e2;
        // Convert line density to density
        Bennet_dense = Bennet_dense / ( elec_dist_radius * elec_dist_radius * PI );
        
        // Default density calculation is incorrect, this yields
        // the correct value
        double actual_density = 2.0 * n_particles * weight 
        / ( nx * ny * delta_x * delta_y );
        double disk_area = elec_dist_radius*elec_dist_radius*PI;
        double disk_fraction = disk_area / (x_len * y_len);
        actual_density = actual_density / disk_fraction;

        printf( "Iz, ampere = %2.2e \n", Iz_ampere);
        printf( "Bennet density, cgs = %2.2e \n", Bennet_dense);
        printf( "Density at sim start, cgs = %2.2e \n", actual_density);
        if( use_RFD==1 ) {
            dt = dt*1/25;
            n_tsteps = 5000;
            save_rate = 100;
            Ex_A = 0.1;
        }

        // TODO: Print all interesting variables and quantities such
        // as debye length, density etc
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
    
    // Note we overwrite this function as we require a velocity
    // in the z-direction in this test case
    std::vector<double> Get_relativistic_vel_and_gamma( double thermal_momentum,
     const double PI, const double c ) {
        
        double p_x = Get_maxwellian(thermal_momentum, PI );
        double p_y = Get_maxwellian(thermal_momentum, PI );
        double p_z = p_z_global;

        double total_momentum = std::sqrt( p_x * p_x + p_y * p_y + p_z * p_z );
        double gamma = std::sqrt( 1 + total_momentum*total_momentum / ( c * c ));
        //printf("p_x: %2.2e p_y: %2.2e gamma: %2.2e \n", p_x, p_y, gamma );
        std::vector<double> vel_and_gamma{ p_x / gamma, p_y / gamma, p_z / gamma, gamma };

        return vel_and_gamma;
    }

    void Generate_electron_positions()
    {

        double x_len = nx * delta_x;
        double y_len = ny * delta_y;

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            
            double radius = elec_dist_radius * std::sqrt(global_random());
            double theta = global_random() * 2.0 * PI;

            e_pos_ic[ip] = x_len / 2.0 + radius * std::cos(theta); /// 2 + x_len/2;
            e_pos_ic[ip + 1] = y_len / 2.0 + radius * std::sin(theta);
            e_pos_ic[ip + 2] = 0.0;
            
        }
    }

    void Generate_positron_positions()
    {
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;
        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            double radius = elec_dist_radius * std::sqrt(global_random());
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
            e_vel_ic[ip + 2] = -vel_and_gamma[2];

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
            p_vel_ic[ip + 2] = vel_and_gamma[2];

            
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

    double Calculate_B_theta(double x, double y ) {
        
        double x_len = nx * delta_x; 
        double radius = std::sqrt( x*x + y*y );
        
        double B_theta;
        if( radius < elec_dist_radius ) {
            return B_theta = 4.0 * n_particles * weight 
            * p_z_global * q_e_cgs * radius
            / ( c * elec_dist_radius * elec_dist_radius );
        } else {
            double disp_pow = std::pow( (radius-elec_dist_radius) /(40*delta_x), 2);
            //printf( "diff: %2.2e \n", (radius-elec_dist_radius)/(17.0*delta_x) );
            return B_theta = 4.0 * n_particles * weight
            * p_z_global * q_e_cgs / (c * radius ) 
            * std::exp(-(disp_pow) ); // Shape factor
        }

    }

    void Set_EM_field()
    {

        double x_len = nx * delta_x; 
        double y_len = ny * delta_y; 
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * delta_x - x_len / 2.0;
                double y = iy * delta_y - y_len / 2.0;

                double theta = std::atan2( y+ 0.9e-9 ,x+1.1e-9 );
                double beta  = PI/2 - theta + PI;

                double B_theta = Calculate_B_theta(x, y);
                

                double B_x = std::cos(beta) * B_theta;
                double B_y = -std::sin(beta) * B_theta;



                EM_ic.E_x[ix + iy * nx] += 0.0;
                EM_ic.E_z[ix + iy * nx] += Ex_A;
                
                EM_ic.B_x[ix + iy * nx] += B_x;
                EM_ic.B_y[ix + iy * nx] += B_y;
            }
        }
        
        
    }
};

#endif //IC_JZ_PINCH
