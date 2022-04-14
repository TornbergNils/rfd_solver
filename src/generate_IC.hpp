#include <vector>
#include <random>
#include <map>

struct IC_struct
{
    std::vector<double> e_pos_ic;
    std::vector<double> p_pos_ic;

    std::vector<double> e_vel_ic;
    std::vector<double> p_vel_ic;

    std::vector<double> e_gamma_ic;
    std::vector<double> p_gamma_ic;

    EM_field_matrix EM_ic;

    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;


    double Get_maxwellian_vel( double v_thermal, const double PI)
    {
        // Set particle positions and velocities
        double U1 = distribution(generator);
        double U2 = distribution(generator);
        // printf( "U1, U2 = %lf, %lf, \n", U1, U2  );
        double maxw1 = v_thermal * std::sqrt(-2 * std::log(U1)) * std::cos(2 * PI * U2);
        // double maxw2 = v_thermal * std::sqrt( -2*std::log(U1))*std::sin(2*PI*U2);
        return maxw1;
    }

    void Generate_electron_positions(
        std::map<std::string, double> &ic_param,
        std::vector<double> &electron_pos)
    {

        auto global_random = std::bind(distribution, generator);
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
            /*
            double radius = x_len / 4 * std::sqrt(global_random());
            double theta = global_random() * 2 * PI;
            electron_pos[ip] = x_len / 2 + radius * std::cos(theta);     // global_random() * x_len/2 + x_len/4; // /2 + x_len/2;
            electron_pos[ip + 1] = y_len / 2 + radius * std::sin(theta); // global_random() * y_len/2 + y_len/4;
            electron_pos[ip + 2] = 0.0;
            */
            electron_pos[ip] = global_random() * x_len;
            electron_pos[ip + 1] = global_random() * y_len;
            electron_pos[ip + 2] = 0.0;
        }
    }

    void Generate_positron_positions(
        std::map<std::string, double> &ic_param,
        std::vector<double> &positron_pos)
    {
        auto global_random = std::bind(distribution, generator);
        int nx = std::lrint(ic_param["nx"]);
        int ny = std::lrint(ic_param["ny"]);
        double delta_x = ic_param["dx"];
        double delta_y = ic_param["dy"];
        int n_particles = ic_param["n_particles"];
        double x_len = nx * delta_x;
        double y_len = ny * delta_y;

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {
            /*
            double radius = x_len / 4 * std::sqrt(global_random());
            double theta = global_random() * 2 * PI;

            positron_pos[ip] = x_len / 2 + radius * std::cos(theta); /// 2 + x_len/2;
            positron_pos[ip + 1] = y_len / 2 + radius * std::sin(theta);
            positron_pos[ip + 2] = 0.0;
            */
            positron_pos[ip] = global_random() * x_len;
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

        auto global_random = std::bind(distribution, generator);
        double wavevector = ic_param["Ex_wavevect"];
        double v_thermal = ic_param["v_thermal"];
        int n_particles = ic_param["n_particles"];
        double c = ic_param["c"];

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {

            double v0 = Get_maxwellian_vel( v_thermal, PI);
            double v1 = 0.1 * v_thermal * std::sin(wavevector * electron_pos[ip]);
            // double v2 = 0.05 * c * std::sin( 2 * PI *electron_pos[ip] / x_len );
            electron_vel[ip] = v0 + v1;
            electron_vel[ip + 1] = 0.0;
            electron_vel[ip + 2] = 0.0;

            double vel_squared_e = electron_vel[ip] * electron_vel[ip] + electron_vel[ip + 1] * electron_vel[ip + 1] + electron_vel[ip + 2] * electron_vel[ip + 2];
            double gamma_e = 1.0 / std::sqrt(1.0 - vel_squared_e / (c * c));
            electron_gamma[ip / 3] = gamma_e;

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

        auto global_random = std::bind(distribution, generator);
        double wavevector = ic_param["Ex_wavevect"];
        double v_thermal = ic_param["v_thermal"];
        int n_particles = ic_param["n_particles"];
        double c = ic_param["c"];

        for (int ip = 0; ip < n_particles * 3; ip += 3)
        {

            double v0 = Get_maxwellian_vel( v_thermal, PI);
            double v1 = -0.1 * v_thermal * std::sin(wavevector * positron_pos[ip]);
            positron_vel[ip] = v0 + v1;
            positron_vel[ip + 1] = 0.0; //-0.5 - 0.1 * std::sin( positron_pos[ip+1] / 100 );
            positron_vel[ip + 2] = 0.0;

            double vel_squared_p = positron_vel[ip] * positron_vel[ip] + positron_vel[ip + 1] * positron_vel[ip + 1] + positron_vel[ip + 2] * positron_vel[ip + 2];
            double gamma_p = 1.0 / std::sqrt(1.0 - vel_squared_p / (c * c));
            positron_gamma[ip / 3] = gamma_p;
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
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * delta_x;
                double y = iy * delta_y;
                // printf( "(%.2lf, %.2lf)", x, y );

                EM_IC.E_x[ix + iy * nx] = 0.0;
                EM_IC.E_y[ix + iy * nx] = Ex_A * std::sin(Ex_k * x);
                EM_IC.E_z[ix + iy * nx] = 0.0; // 2000; // Ex_A*std::cos( Ex_k * x );

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

        /*
        wave_config_init[0][0] = wave1_A; // wave2_A;
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
            EM_IC.E_x[ix + iy * nx] += Get_EM_wave_component(0, config2, x, y, 0);
            EM_IC.E_y[ix + iy * nx] += Get_EM_wave_component(1, config2, x, y, 0);
            EM_IC.E_z[ix + iy * nx] += Get_EM_wave_component(2, config2, x, y, 0);
            //printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

            EM_IC.B_x[ix + iy * nx] += Get_EM_wave_component(3, config2, x, y, 0);
            EM_IC.B_y[ix + iy * nx] += Get_EM_wave_component(4, config2, x, y, 0);
            EM_IC.B_z[ix + iy * nx] += Get_EM_wave_component(5, config2, x, y, 0);

          // printf( "\n" );
        }
        }
        */
    }
};
