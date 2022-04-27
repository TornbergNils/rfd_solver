#ifndef IC_GEN_H
#define IC_GEN_H

#include <vector>
#include <random>
#include <map>


/*
    Class containing initial conditions for simulation. 

    Use by creating a parameters.txt via the python script and
    then loading them into a object of this type. If you want
    to run a specific experiment, use one of the experiment classes
    instead, which have experiment parameters hard-coded.
*/
class IC_struct
{
public:
    double k_boltz_erg = 1.38064 * 1e-16;
    double k_boltz_ev = (8.617 * 1e-5);
    double c = 2.99792458 * 1e10;
    double m_e_cgs = 9.1093819*1e-28;
    double q_e_cgs = 4.80320425e-10;

    int n_particles;
    int nx;
    int ny;
    int weight;
    double m_e;
    double q_e;
    
    bool use_RFD;
    int n_tsteps;
    int save_rate;

    double plasma_wavelen;
    double plasma_wavenum;
    double wavevector;

    double x_min;
    double x_max;

    // Note same step length for x, y
    double delta_x;
    double delta_y;
    double dt;
    double tmax;

    // Physical properties

    double Te; // eV
    double T_kelvin;
    double electron_momentum;
        
    double wave1_A;
    double wave1_k;
    double wave2_A;
    double wave2_k;
    double Ex_A;
    double Ex_k;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::vector<double> e_pos_ic;
    std::vector<double> p_pos_ic;

    std::vector<double> e_vel_ic;
    std::vector<double> p_vel_ic;

    std::vector<double> e_gamma_ic;
    std::vector<double> p_gamma_ic;

    EM_field_matrix EM_ic;

    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;

    const double PI = 3.14159265358979;

    IC_struct(
        int n_particles, 
        int nx,
        int ny,
        int weight,
        bool use_RFD,
        int n_tsteps,
        int save_rate,

        double plasma_wavelen,

        double x_min,
        double x_max,
        // Physical properties

        double Te, // eV

        double wave1_A,
        double wave1_k,
        double wave2_A,
        double wave2_k,
        double Ex_A,
        double Ex_k
        ) 
        : n_particles(n_particles),
        nx(nx), ny(ny), 
        weight(weight), use_RFD(use_RFD),
        n_tsteps(n_tsteps), save_rate(save_rate),
        plasma_wavelen(plasma_wavelen),
        x_min(x_min), x_max(x_max), Te(Te),
        wave1_A(wave1_A), wave1_k(wave1_k),
        wave2_A(wave2_A), wave2_k(wave2_k),
        Ex_A(Ex_A), Ex_k(Ex_k),
        e_pos_ic(n_particles*3),
        p_pos_ic(n_particles*3),
        e_vel_ic(n_particles*3),
        p_vel_ic(n_particles*3),
        e_gamma_ic(n_particles),
        p_gamma_ic(n_particles),
        EM_ic(nx, ny),
        generator(seed),
        distribution(0.0, 1.0 ) {

        m_e = m_e_cgs * weight;
        q_e = q_e_cgs * weight;
    
        plasma_wavenum = 2*PI/plasma_wavelen;
        wavevector = plasma_wavenum;
    
        delta_x = (x_max - x_min ) / nx;
        delta_y = (x_max - x_min ) / nx;
        dt = delta_x / (2*c);
        tmax = n_tsteps * dt;
    
        T_kelvin = Te / k_boltz_ev;
        electron_momentum = std::sqrt( T_kelvin * k_boltz_erg  / ( m_e_cgs ) );

         }

    double global_random() { return distribution(generator); }
    
    double Get_maxwellian( double stdev, const double PI)
    {
        // Set particle positions and velocities
        double U1 = global_random();
        double U2 = global_random();
        // printf( "U1, U2 = %lf, %lf, \n", U1, U2  );
        double maxw1 = stdev * std::sqrt(-2 * std::log(U1)) * std::cos(2 * PI * U2);
        // double maxw2 = v_thermal * std::sqrt( -2*std::log(U1))*std::sin(2*PI*U2);
        return maxw1;
    }

    //double Velocity_from_gamma( double gamma, double c ) {
    //    return std::sqrt(1.0-1.0/(gamma*gamma))*c;
    //}

    std::vector<double> Get_relativistic_vel_and_gamma( double thermal_momentum,
     const double PI, const double c ) {
        
        double p_x = Get_maxwellian(thermal_momentum, PI );
        double p_y = Get_maxwellian(thermal_momentum, PI );
        //double p_z = Get_maxwellian();

        double total_momentum = std::sqrt( p_x * p_x + p_y * p_y );
        double gamma = std::sqrt( 1 + total_momentum*total_momentum / ( c * c ));
        //printf("p_x: %2.2e p_y: %2.2e gamma: %2.2e \n", p_x, p_y, gamma );
        std::vector<double> vel_and_gamma{ p_x / gamma, p_y / gamma, 0, gamma };

        return vel_and_gamma;
    }

    void Generate_electron_positions(
        std::map<std::string, double> &ic_param,
        std::vector<double> &electron_pos)
    {

        int nx = std::lrint(ic_param["nx"]);
        int ny = std::lrint(ic_param["ny"]);
        double delta_x = ic_param["dx"];
        double delta_y = ic_param["dy"];
        int n_particles = ic_param["n_particles"];
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {

            
            // Place particles in a circle
            // double radius = x_len / 8 * std::sqrt( dist(gen) );
            // double theta = dist(gen) * 2 * PI;
            
            //electron_pos[ip] = x_len / 2 + radius * std::cos(theta);     // global_random() * x_len/2 + x_len/4; // /2 + x_len/2;
            //electron_pos[ip + 1] = y_len / 2 + radius * std::sin(theta); // global_random() * y_len/2 + y_len/4;
            //electron_pos[ip + 2] = 0.0;
            
            electron_pos[ip] =     global_random() * x_len/3 + x_len/3; ///3 + x_len/3;
            electron_pos[ip + 1] = global_random() * y_len;
            electron_pos[ip + 2] = 0.0;
            
        }
    }

    void Generate_positron_positions(
        std::map<std::string, double> &ic_param,
        std::vector<double> &positron_pos)
    {
        double electron_momentum = ic_param["electron_momentum"];
        int nx = std::lrint(ic_param["nx"]);
        int ny = std::lrint(ic_param["ny"]);
        double delta_x = ic_param["dx"];
        double delta_y = ic_param["dy"];
        int n_particles = ic_param["n_particles"];
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            double radius = x_len / 8 * std::sqrt(global_random());
            double theta = global_random() * 2 * PI;

            // positron_pos[ip] = x_len / 2 + radius * std::cos(theta); /// 2 + x_len/2;
            // positron_pos[ip + 1] = y_len / 2 + radius * std::sin(theta);
            // positron_pos[ip + 2] = 0.0;
            
            positron_pos[ip] =     global_random() * x_len/3 + x_len/3;
            positron_pos[ip + 1] = global_random() * y_len;
            positron_pos[ip + 2] = 0.0;
            
        }
    }

    void Generate_electron_velocities(
        std::map<std::string, double> &ic_param,
        std::vector<double> &electron_pos,
        std::vector<double> &electron_vel,
        std::vector<double> &electron_gamma)
    {

        double electron_momentum = ic_param["electron_momentum"];
        double wavevector = ic_param["Ex_wavevect"];
        double v_thermal = ic_param["v_thermal"];
        int n_particles = ic_param["n_particles"];
        double c = ic_param["c"];

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            std::vector<double> vel_and_gamma = Get_relativistic_vel_and_gamma(electron_momentum, PI, c);

            /*
            double v0 = Get_maxwellian( v_thermal, PI);
            double v1 = 0.0 * v_thermal * std::sin(wavevector * electron_pos[ip]);
            // double v2 = 0.05 * c * std::sin( 2 * PI *electron_pos[ip] / x_len );
            */

            double v1 = 0.05 * electron_momentum * std::sin(wavevector * electron_pos[ip]);
            electron_vel[ip] = vel_and_gamma[0];
            electron_vel[ip + 1] = vel_and_gamma[1];
            electron_vel[ip + 2] = 0.0;

            double vel_squared_e = electron_vel[ip] * electron_vel[ip] 
                + electron_vel[ip + 1] * electron_vel[ip + 1] 
                + electron_vel[ip + 2] * electron_vel[ip + 2];
            
            //double gamma_e = 1.0 / std::sqrt(1.0 - vel_squared_e / (c * c));
            electron_gamma[ip / 3] = vel_and_gamma[3];

            if (vel_squared_e > c * c)
            {
                printf("e particle speed exceeds c! c = %lf \n", c);
            }
        }
    }

    void Generate_positron_velocities(
        std::map<std::string, double> &ic_param,
        std::vector<double> &positron_pos,
        std::vector<double> &positron_vel,
        std::vector<double> &positron_gamma)
    {

        double electron_momentum = ic_param["electron_momentum"];
        double wavevector = ic_param["Ex_wavevect"];
        double v_thermal = ic_param["v_thermal"];
        int n_particles = ic_param["n_particles"];
        double c = ic_param["c"];

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            
            // vel and gamma contains vx, vy, vz, gamma
            std::vector<double> vel_and_gamma = Get_relativistic_vel_and_gamma(electron_momentum, PI, c);
           
           /* 
            double v0 = Get_maxwellian( v_thermal, PI);
            v0 = Get_maxwellian( v_thermal, PI);
            */
            double v1 = -0.05 * electron_momentum * std::sin(wavevector * positron_pos[ip]);
            positron_vel[ip] = vel_and_gamma[0];
            positron_vel[ip + 1] = vel_and_gamma[1]; //-0.5 - 0.1 * std::sin( positron_pos[ip+1] / 100 );
            positron_vel[ip + 2] = 0.0;

            
            double vel_squared_p = positron_vel[ip] * positron_vel[ip] 
            + positron_vel[ip + 1] * positron_vel[ip + 1] 
            + positron_vel[ip + 2] * positron_vel[ip + 2];
            /*
            double gamma_p = 1.0 / std::sqrt(1.0 - vel_squared_p / (c * c));
            */

            positron_gamma[ip / 3] = vel_and_gamma[3];
            if (vel_squared_p > c * c)
            {
                printf("p particle speed exceeds c! \n");
            }
        }
    }

    double Gaussian(double x, double y, double ampli, double stdev)
    {
        double variable = std::exp(-(x * x + y * y) / (stdev * stdev));
        return ampli * variable;
    }

    void Set_EM_field(
        EM_field_matrix &EM_IC,
        std::map<std::string, double> &ic_param)
    {

        int nx = ic_param["nx"];
        int ny = ic_param["ny"];
        double delta_x = ic_param["dx"];
        double delta_y = ic_param["dy"];
        double wave1_A = ic_param["wave1_amplitude"];
        double wave1_k = ic_param["wave1_wavevect"];
        double wave2_A = ic_param["wave2_amplitude"];
        double wave2_k = ic_param["wave2_wavevect"];
        double Ex_A = ic_param["Ex_raw"];
        double Ex_k = ic_param["Ex_wavevect"];
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
            for (int ix = 0; ix < nx/3; ix++)
            {
                double x = ix * delta_x;
                double y = iy * delta_y;
                // printf( "(%.2lf, %.2lf)", x, y );

                EM_IC.E_x[ix + iy * nx] = 0.0;
                EM_IC.E_y[ix + iy * nx] = Ex_A * std::sin(Ex_k * x);
                EM_IC.E_z[ix + iy * nx] = 0.0; //4e4; // 2000; // Ex_A*std::cos( Ex_k * x );

                EM_IC.B_x[ix + iy * nx] = 0;
                EM_IC.B_y[ix + iy * nx] = 0.0; // Ex_A*std::sin( Ex_k * x );
                EM_IC.B_z[ix + iy * nx] = Ex_A * std::sin(Ex_k * x);

                /*
                EM_IC.E_x[ix + iy * nx] += std::copysign(1.0,x-nx*delta_x/2)
                  * Gaussian( x-nx*delta_x/2, y-ny*delta_y/2, wave1_A*0.2, 10*delta_x );

                EM_IC.E_y[ix + iy * nx] += std::copysign(1.0,y-ny*delta_y/2)
                  * Gaussian( x-nx*delta_x/2, y-ny*delta_y/2, wave1_A*0.2, 10*delta_x );
                */

                EM_IC.E_x[ix + iy * nx] += Get_EM_wave_component(0, config, x, y, 0);
                EM_IC.E_y[ix + iy * nx] += Get_EM_wave_component(1, config, x, y, 0);
                EM_IC.E_z[ix + iy * nx] += Get_EM_wave_component(2, config, x, y, 0);
                // printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

                EM_IC.B_x[ix + iy * nx] += Get_EM_wave_component(3, config, x, y, 0);
                EM_IC.B_y[ix + iy * nx] += Get_EM_wave_component(4, config, x, y, 0);
                EM_IC.B_z[ix + iy * nx] += Get_EM_wave_component(5, config, x, y, 0);
            }
        }

        
        wave_config_init[0][0] = wave1_A; // wave2_A;
        wave_config_init[0][1] = -wave1_k;
        wave_config_init[0][2] = -PI;
        wave_config_init[0][3] = 0.0;
        EM_wave_config config2(wave_config_init);

        for (int iy = 0; iy < ny; iy++)
        {
          for (int ix = (nx*2)/3; ix < nx; ix++)
          {
            double x = ix * delta_x;
            double y = iy * delta_y;
            EM_IC.E_x[ix + iy * nx] += Get_EM_wave_component(0, config2, x, y, 0);
            EM_IC.E_y[ix + iy * nx] += Get_EM_wave_component(1, config2, x, y, 0);
            EM_IC.E_z[ix + iy * nx] += Get_EM_wave_component(2, config2, x, y, 0);
            //printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

            //EM_IC.B_x[ix + iy * nx] += Get_EM_wave_component(3, config2, x, y, 0);
            EM_IC.B_y[ix + iy * nx] += Get_EM_wave_component(4, config2, x, y, 0);
            EM_IC.B_z[ix + iy * nx] += Get_EM_wave_component(5, config2, x, y, 0);

          // printf( "\n" );
        }
        }
        
        
    }
        
    void print_primitives() {
        std::cout             << n_particles << "\n" <<
        "nx             =" <<  nx  << "\n" <<
        "ny             =" <<  ny  << "\n" <<
        "weight         =" <<  weight  << "\n" <<
        "m_e            =" <<  m_e  << "\n" <<
        "q_e            =" <<  q_e  << "\n" <<
    
        "use_RFD        =" <<  use_RFD  << "\n" <<
        "n_tsteps       =" <<  n_tsteps  << "\n" <<
        "save_rate      =" <<  save_rate  << "\n" <<
    
        "plasma_wavelen =" <<  plasma_wavelen  << "\n" <<
        "plasma_wavenum =" <<  plasma_wavenum  << "\n" <<
        "wavevector     =" <<  wavevector  << "\n" <<
    
        "x_min          =" <<  x_min  << "\n" <<
        "x_max          =" <<  x_max  << "\n" <<
    
        "delta_x        =" <<  delta_x  << "\n" <<
        "delta_y        =" <<  delta_y  << "\n" <<
        "dt             =" <<  dt  << "\n" <<
        "tmax           =" <<  tmax  << "\n" <<
        "Te             =" <<  Te  << "\n" <<
        "T_kelvin       =" <<  T_kelvin  << "\n" <<
        "electron_p     =" <<  electron_momentum  << "\n" <<
    
        "wave1_A        =" <<  wave1_A  << "\n" <<
        "wave1_k        =" <<  wave1_k  << "\n" <<
        "wave2_A        =" <<  wave2_A  << "\n" <<
        "wave2_k        =" <<  wave2_k  << "\n" <<
        "Ex_A           =" <<  Ex_A  << "\n" <<
        "Ex_k           =" <<  Ex_k  << "\n";
       }
};

#endif //IC_GEN_H