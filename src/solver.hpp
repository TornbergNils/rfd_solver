
//const double q_e_cgs = 6.03587913*1e-9;
//const double m_e_cgs = 9.1094*1e-28;
//const double PI = 3.14159265358979;
//const double c_SI = 2.99792458 * 1e8;
//// ???
//const double Me = 9.1093819*1e-31;
//const double Kb = 1.380650*1e-23;
//const double Me_by_Kb = 6.59789976*1e-8;

class Solver
{
  // Given parameters
  int nx;
  int ny;
  int n_particles;
  int n_tsteps;
  int save_rate;
  double dt;
  std::vector<double>::size_type n_elec;
  std::vector<double>::size_type n_posi;

  double delta_x;
  double delta_y;
  

  // Physical constants
  double c;
  double v_thermal;
  double q_e;
  double m_e;
  
  // Derived parameters
  double tmax;


  // Data
  std::vector<double> electron_pos;
  std::vector<double> positron_pos;
  std::vector<double> electron_vel;
  std::vector<double> positron_vel;

  std::vector<double> Jx;
  std::vector<double> Jy;
  std::vector<double> Jz;
  
  std::vector<double> rho_q;

  EM_field_matrix EM;
  RFD_matrix RFD;

public:
  Solver(int nx, int ny, int n_particles, double dt, int n_tsteps, int save_rate,
         double delta_x, double delta_y, std::map< std::string, double > ic_param )
      : nx(nx), ny(ny), n_particles(n_particles), n_tsteps(n_tsteps), 
        save_rate(save_rate), dt(dt), n_elec(n_particles),
        n_posi(n_particles), delta_x(delta_x), delta_y(delta_y),
        electron_pos(n_particles * 3), positron_pos(n_particles * 3),
        electron_vel(n_particles * 3), positron_vel(n_particles * 3),
        Jx(nx * ny), Jy(nx * ny), Jz(nx * ny), rho_q(nx * ny),
        EM(nx, ny), RFD(EM, 1) {

          v_thermal = ic_param["v_thermal"];
          c = ic_param["c"];
          q_e = ic_param["q_e"];
          m_e = ic_param["m_e"];
          printf("Constructiong! param: %lf, %lf, %lf, %lf", v_thermal, c, q_e, m_e );

        }

  void Interpolate_charge_boris()
  {
    // For each electron, round down position to get grid position
    // J is co-located with Ez
    double sign = -1 * q_e; // * q_e_cgs;
    
    for (long unsigned ip = 0; ip < n_elec * 3; ip += 3)
    {
      double x = electron_pos[ip];
      double y = electron_pos[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;

      rho_q[Get_index(ix, iy)] += sign * w00;
      rho_q[Get_index(ix+1, iy)] += sign * w10;
      rho_q[Get_index(ix+1, iy+1)] += sign * w11;
      rho_q[Get_index(ix, iy+1)] += sign * w01;
      
    }
  }

  void Interpolate_current_component(std::vector<double> &current, double vel, int ix,
                                       int iy, double w00, double w10, double w11, double w01, double gamma, double sign )
  {
      //printf("half current at pt: %2.2e \n", sign * vel / ( gamma * 2 )  );
      current[Get_index(ix, iy)] += sign * w00 * vel / ( gamma * 2 );
      current[Get_index(ix+1, iy)] += sign * w10 * vel / ( gamma * 2 );
      current[Get_index(ix+1, iy+1)] += sign * w11 * vel / ( gamma * 2 );
      current[Get_index(ix, iy+1)] += sign * w01 * vel / ( gamma * 2 );
  }
  void Interpolate_half_current_boris()
  {
    // For each electron, round down position to get grid position
    // J is co-located with Ez
    double sign = -1* q_e / (delta_x * delta_y); // * q_e_cgs;
    
    for (long unsigned ip = 0; ip < n_elec * 3; ip += 3)
    {
      double x = electron_pos[ip];
      double y = electron_pos[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;

      double vel_squared = electron_vel[ip] * electron_vel[ip] 
      + electron_vel[ip+1] * electron_vel[ip+1] 
      + electron_vel[ip+2] * electron_vel[ip+2];

      double gamma = std::sqrt(1.0 + vel_squared/(c*c));
      
      Interpolate_current_component( Jx, electron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jy, electron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jz, electron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
    }
    
    // repeat for positrons
    /*
    sign = 1 * q_e; // * q_e_cgs;
    for (long unsigned ip = 0; ip < n_posi * 3; ip += 3)
    {
      double x = positron_pos[ip];
      double y = positron_pos[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;



      double vel_squared = positron_vel[ip] * positron_vel[ip] 
      + positron_vel[ip+1] * positron_vel[ip+1] 
      + positron_vel[ip+2] * positron_vel[ip+2];
      double gamma = std::sqrt(1.0 + vel_squared/(c*c));

      Interpolate_current_component( Jx, positron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jy, positron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jz, positron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
    }
  */    
  }
  void Interpolate_half_current_RFD_elec() {

    // For each electron, round down position to get grid position
    // J is co-located with Ez
    double sign = -1 * q_e;
    
    for (long unsigned ip = 0; ip < n_elec * 3; ip += 3)
    {
      double x = electron_pos[ip];
      double y = electron_pos[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;

      Interpolate_current_component( Jx, RFD.RFD_x[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
      Interpolate_current_component( Jy, RFD.RFD_y[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
      Interpolate_current_component( Jz, RFD.RFD_z[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
    }
  }
  void Interpolate_half_current_RFD_posi() {

    // For each positron, round down position to get grid position
    // J is co-located with Ez
    double sign = 1 * q_e;
    
    for (long unsigned ip = 0; ip < n_elec * 3; ip += 3)
    {
      double x = electron_pos[ip];
      double y = electron_pos[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;

      Interpolate_current_component( Jx, RFD.RFD_x[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
      Interpolate_current_component( Jy, RFD.RFD_y[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
      Interpolate_current_component( Jz, RFD.RFD_z[ip/3], ix,
                                       iy, w00, w10, w11, w01, 1, sign );
    }
  }

  double Get_maxwellian_vel( std::mt19937 &gen, std::uniform_real_distribution<double> &dist,
    double v_thermal, const double PI ) {
      double U1 = dist(gen);
      double U2 = dist(gen);
      //printf( "U1, U2 = %lf, %lf, \n", U1, U2  );
      double maxw1 = v_thermal * std::sqrt( -2*std::log(U1))*std::cos(2*PI*U2);
      //double maxw2 = v_thermal * std::sqrt( -2*std::log(U1))*std::sin(2*PI*U2);
      return maxw1;
  }

  void Initialize(EM_field_matrix EM_IC)
  {
    printf("Initializing! param: %lf, %lf, %lf, %lf", v_thermal, c, q_e, m_e );
    // Setup rng
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //seed = 3;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    auto random = std::bind(distribution, generator);

    tmax = dt * n_tsteps;
    //dt = (double)tmax / (double)n_tsteps;
    double x_len = nx * delta_x;
    double y_len = ny * delta_y;

    // set EM to equal to EM_IC, since EM is 0-initialized by constructor
    EM = EM_IC;
    RFD.Update(EM, 1);
    
    // Constants for calculating system properties and debugging
    double gamma_e;
    
    // Set particle positions and velocities
    for (int ip = 0; ip < n_particles * 3; ip += 3)
    {

      // Place particles uniformly
      electron_pos[ip] = random() * x_len; 
      electron_pos[ip + 1] = random() * y_len;
      electron_pos[ip + 2] = 0.0;

      positron_pos[ip] = random() * x_len;
      positron_pos[ip + 1] = random() * y_len;
      positron_pos[ip + 2] = 0.0;

      double v0 = Get_maxwellian_vel( generator, distribution, v_thermal, PI );
      double v1 = 0.1 * v_thermal * std::sin( 4 * PI *electron_pos[ip] / x_len );
      //double v2 = 0.05 * c * std::sin( 2 * PI *electron_pos[ip] / x_len );
      electron_vel[ip] = v0 + v1;
      electron_vel[ip + 1] = 0.0;
      electron_vel[ip + 2] = 0.0;
      

      double vel_squared_e = electron_vel[ip] * electron_vel[ip]
      + electron_vel[ip+1] * electron_vel[ip+1]
      + electron_vel[ip+2] * electron_vel[ip+2];
      gamma_e = 1.0 / std::sqrt(1.0 - vel_squared_e / (c*c) );
      
      if( vel_squared_e > c*c ) {printf( "e particle speed exceeds c! \n" );}
      electron_vel[ip] *= gamma_e;
      electron_vel[ip + 1] *= gamma_e; 
      electron_vel[ip + 2] *= gamma_e;

      positron_vel[ip] = 0.0;
      positron_vel[ip + 1] = 0.0; //-0.5 - 0.1 * std::sin( positron_pos[ip+1] / 100 );
      positron_vel[ip + 2] = 0.0;
      
      double vel_squared_p = positron_vel[ip] * positron_vel[ip]
      + positron_vel[ip+1] * positron_vel[ip+1]
      + positron_vel[ip+2] * positron_vel[ip+2];
      double gamma_p = 1.0 / std::sqrt(1.0 - vel_squared_p/(c*c) );
      if( vel_squared_p > c*c ) {printf( "p particle speed exceeds c! \n" );}
      
      positron_vel[ip] *= gamma_p;
      positron_vel[ip + 1] *= gamma_p; 
      positron_vel[ip + 2] *= gamma_p;


    }
    Write_vector_to_binary( std::string("./data/initial_velocities"), electron_vel, 0 );

    // Wonky fix to show approx current on first iter
    Interpolate_half_current_boris();
    Interpolate_half_current_boris();
    Interpolate_charge_boris();
    //Interpolate_half_current_RFD_elec();
    //Interpolate_half_current_RFD_elec();
    //Interpolate_half_current_RFD_posi();
    //Interpolate_half_current_RFD_posi();
  }
  void Save_parameters_to_text(std::string filename, int form)
  {
    if (form == 0)
    {
      std::ofstream filestream(filename);
      filestream << "nx, " << nx << "\n"
                 << "ny, " << ny << "\n"
                 << "n_particles, " << n_particles << "\n"
                 << "tmax, " << tmax << "\n"
                 << "n_tsteps, " << n_tsteps << "\n"
                 << "save_rate, " << save_rate << "\n"
                 << "n_elec, " << n_elec << "\n"
                 << "n_posi, " << n_posi << "\n"
                 << "delta_x, " << delta_x << "\n"
                 << "delta_y, " << delta_y << "\n"
                 << "dt, " << dt << "\n";
      filestream.close();
    }
    else if (form == 1)
    {
      std::ofstream filestream(filename);
      filestream << "nx, "
                 << "ny, "
                 << "n_particles, "
                 << "tmax, "
                 << "n_tsteps, "
                 << "save_rate, "
                 << "n_elec, "
                 << "n_posi, "
                 << "delta_x, "
                 << "delta_y, "
                 << "dt \n"
                 << nx << ", " << ny << ", " << n_particles << ", " << tmax
                 << ", " << n_tsteps << ", " << save_rate << ", " << n_elec
                 << ", " << n_posi << ", " << delta_x << ", " << delta_y << ", "
                 << dt << "\n";
    }
  }
  int Write_vector_to_binary(std::string filename, std::vector<double> vect,
                             bool append)
  {
    std::ofstream filestream;
    if (append == true)
    {
      filestream.open(filename,
                      std::ios::out | std::ios::app | std::ios::binary);
    }
    else
    {
      filestream.open(filename,
                      std::ios::out | std::ios::trunc | std::ios::binary);
    }

    filestream.write((char *)&vect[0], vect.size() * sizeof(double));
    filestream.close();

    return 0;
  }

  int Propagate_particles(std::vector<double> &particles,
                          const RFD_matrix &RFD_ap, const double dt)
  {
    int ip_max = particles.size();
    if (RFD_ap.RFD_x.size() != ip_max / 3)
    {
      printf(" RFD vector and particle vector size mismatch!");
    }
    for (int ip = 0, irf = 0; ip < ip_max; ip += 3, irf++)
    {
      particles[ip] += RFD_ap.RFD_x[irf] * c * dt;
      particles[ip + 1] += RFD_ap.RFD_y[irf] * c * dt;
      particles[ip + 2] += RFD_ap.RFD_z[irf] * c * dt;
      
      // Periodic Bc for particles
      particles[ip] =      std::fmod( particles[ip] + nx * delta_x, nx * delta_x );
      particles[ip + 1] =  std::fmod( particles[ip+1] + ny * delta_y, ny * delta_y );
    }
    return 0;
  }

  std::vector<double> Get_cross_product(std::vector<double> u,
                                        std::vector<double> v)
  {
    double x = u[1] * v[2] - u[2] * v[1];
    double y = -1 * (u[0] * v[2] - u[2] * v[0]);
    double z = u[0] * v[1] - u[1] * v[0];
    std::vector<double> cprod{x, y, z};
    return cprod;
  }

  int Boris_velocity(std::vector<double> &particles,
                                std::vector<double> &vel,
                                const EM_field_matrix &EM_ap,
                                const int sign)
  {
    // Currently natural units w scale by eV
    // double q_div_by_m = sign * q_e / m_e; //q_e_cgs / m_e_cgs;
    // const double prop_factor = q_div_by_m / 2.0 * dt;
    int ip_max = particles.size();
    if (EM_ap.E_x.size() != ip_max / 3)
    {
      printf(" EM vector and particle vector size mismatch!");
    }

    for (int ip = 0, iem = 0; ip < ip_max; ip += 3, iem++)
    {

      // Important! note u = gamma * v, gamma is the lorentz factor
      std::vector<double> u{vel[ip],
                            vel[ip + 1],
                            vel[ip + 2]};
      std::vector<double> v(3);

      // Add half impulse from E-field to get u- in v
      v[0] = u[0] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_x[iem];  //prop_factor * EM_ap.E_x[iem];
      v[1] = u[1] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_y[iem];  //prop_factor * EM_ap.E_y[iem];
      v[2] = u[2] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_z[iem];  //prop_factor * EM_ap.E_z[iem];

      double u_minus_squared = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      // Recalculate gamma at different time
      double gamma = std::sqrt(1.0 + u_minus_squared/(c*c));
      if( gamma > 1000 ) {
        printf( "gamma1: %lf \n", gamma );
        }
      //  Get t vector and put it in u
      u[0] = (sign * q_e * dt ) / ( 2 * m_e * gamma * c ) * EM_ap.B_x[iem];
      u[1] = (sign * q_e * dt ) / ( 2 * m_e * gamma * c ) * EM_ap.B_y[iem];
      u[2] = (sign * q_e * dt ) / ( 2 * m_e * gamma * c ) * EM_ap.B_z[iem];

      // Get s = 2t/(1+t^2)
      double t_squared = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      std::vector<double> s{2.0 * u[0] / (1.0 + t_squared),
                            2.0 * u[1] / (1.0 + t_squared),
                            2.0 * u[2] / (1.0 + t_squared)};

      // replace u by cross product of u- and t
      u = Get_cross_product(v, u);
      // perform first part of boris rotation putting u prime in u
      u[0] = v[0] + u[0];
      u[1] = v[1] + u[1];
      u[2] = v[2] + u[2];

      // u prime cross s
      s = Get_cross_product(u, s);

      // perform second part of boris rotation putting u+ in u
      u[0] = v[0] + s[0];
      u[1] = v[1] + s[1];
      u[2] = v[2] + s[2];

      // Add remaining half impulse from E-field to get u at next timestep
      // from u+
      u[0] = u[0] + (sign * q_e * dt ) / ( 2 * m_e) * EM_ap.E_x[iem];
      u[1] = u[1] + (sign * q_e * dt ) / ( 2 * m_e) * EM_ap.E_y[iem];
      u[2] = u[2] + (sign * q_e * dt ) / ( 2 * m_e) * EM_ap.E_z[iem];

      // double u_now_squared = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      // gamma = std::sqrt(1.0 + u_now_squared/(c*c));
      //printf(", gamma2: %lf \n", gamma);
      // Finally update velocity vector
      vel[ip] = u[0];
      vel[ip + 1] = u[1];
      vel[ip + 2] = u[2];
      }


    return 0;
  }
  void Boris_move_particles( std::vector<double> &pos, std::vector<double> &vel ) {
    
    int ip_max = electron_pos.size();
    for (int ip = 0, iem = 0; ip < ip_max; ip += 3, iem++)
    {
      // note +=
      // Note also particles and vel known at offset times!
      
      double vel_squared = vel[ip] * vel[ip] 
      + vel[ip+1] * vel[ip+1] 
      + vel[ip+2] * vel[ip+2];

      double gamma = std::sqrt(1.0 + vel_squared/(c*c));
      
      pos[ip] += vel[ip] * dt / gamma;
      pos[ip + 1] += vel[ip + 1] * dt / gamma;
      pos[ip + 2] += vel[ip + 2] * dt / gamma;

      // Periodic Bc for particles
      pos[ip] =      std::fmod( pos[ip] + nx * delta_x, nx * delta_x );
      pos[ip + 1] =  std::fmod( pos[ip+1] + ny * delta_y, ny * delta_y );
    }
  }
  inline int Get_index(int ix, int iy)
  {
    return ((ix + nx) % nx) + ((iy + ny) % ny) * nx;
  }


  double Interpolate_field_component(const std::vector<double> &field, int ix,
                                     int iy, double x, double y)
  {
    if (x > delta_x || y > delta_y || x < 0.0 || y < 0.0 )
    {
      printf("Error! Interpolation coordinate sizes exceed bounds!");
      printf("x,y %lf, %lf, ix, iy %d, %d, \n", x, y, ix, iy);
    }
    double cell_size = delta_x * delta_y;

    double w00 = (delta_x - x) * (delta_y - y) / cell_size;
    double w10 = x * (delta_y - y) / cell_size;
    double w11 = x * y / cell_size;
    double w01 = (delta_x - x) * y / cell_size;

    return field[Get_index(ix, iy)] * w00 + field[Get_index(ix + 1, iy)] * w10 +
           field[Get_index(ix + 1, iy + 1)] * w11 +
           field[Get_index(ix, iy + 1)] * w01;
  }

  EM_field_matrix
  Interpolate_EM_at_particles(const std::vector<double> &particle_pos)
  {
    int num_particles = particle_pos.size() / 3;
    EM_field_matrix interpolated_field(num_particles, 1);

    // For each paticle, round down position to get nearest lower left cornet
    // indices
    for (int ip = 0; ip < num_particles * 3; ip += 3)
    {
      double xp = particle_pos[ip];
      double yp = particle_pos[ip + 1];
      int Ez_ix = std::floor(xp / delta_x);
      int Ez_iy = std::floor(yp / delta_y);

      int ExHy_ix = std::floor((xp - delta_x / 2) / delta_x);
      int ExHy_iy = std::floor(yp / delta_y);

      int EyHx_ix = std::floor(xp / delta_x);
      int EyHx_iy = std::floor((yp  - delta_y / 2) / delta_y);

      int Hz_ix = std::floor((xp  - delta_x / 2) / delta_x);
      int Hz_iy = std::floor((yp  - delta_y / 2) / delta_y);
      // double z = particle_pos[ip+2];
      // int iz = std::floor( z );

      // need new x,y relative to lower left corner of box
      // Ex and Hy are displaced by delta_x / 2 in x-direction
      double x = xp - ExHy_ix * delta_x - delta_x / 2;
      double y = yp - ExHy_iy * delta_y;
      
      //double temp = Interpolate_field_component(EM.E_x, ExHy_ix, ExHy_iy, x, y);
      //if( std::isnan( temp )) { printf("nan at particle %d", ip/3 ); }
      
      interpolated_field.E_x[ip / 3] = Interpolate_field_component(EM.E_x, ExHy_ix, ExHy_iy, x, y);


      // Ey and Hx are displaced by delta_y / 2 in y-direction
      x = xp - EyHx_ix * delta_x;
      y = yp - EyHx_iy * delta_y - delta_y / 2;
      interpolated_field.E_y[ip / 3] =
          Interpolate_field_component(EM.E_y, EyHx_ix, EyHx_iy, x, y);

      // Ez is not displaced
      x = xp - Ez_ix * delta_x;
      y = yp - Ez_iy * delta_y;
      interpolated_field.E_z[ip / 3] =
          Interpolate_field_component(EM.E_z, Ez_ix, Ez_iy, x, y);

      // Ey and Hx are displaced by delta_y / 2 in y-direction
      x = xp - EyHx_ix * delta_x;
      y = yp - EyHx_iy * delta_y - delta_y / 2;
      interpolated_field.B_x[ip / 3] =
          Interpolate_field_component(EM.B_x, EyHx_ix, EyHx_iy, x, y);

      // Ex and Hy are displaced by delta_x / 2 in x-direction
      x = xp - ExHy_ix * delta_x - delta_x / 2;
      y = yp - ExHy_iy * delta_y;
      interpolated_field.B_y[ip / 3] =
          Interpolate_field_component(EM.B_y, ExHy_ix, ExHy_iy, x, y);

      // Hz is displaced by delta/2 in along both axes
      x = xp - Hz_ix * delta_x - delta_x / 2;
      y = yp - Hz_iy * delta_y - delta_y / 2;
      interpolated_field.B_z[ip / 3] =
          Interpolate_field_component(EM.B_z, Hz_ix, Hz_iy, x, y);
    }

    return interpolated_field;
  }

  RFD_matrix Calculate_RFD_at_particles(const EM_field_matrix &EM, int sign)
  {
    RFD_matrix RFD_temp(EM, sign);
    return RFD_temp;
  }
  void FDTD()
  {

    // Start with H parts of modes, "1/2 timesteps"
    for (int iy = 0; iy < ny; iy++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        int index = Get_index(ix, iy);

        EM.B_z[index] =
            EM.B_z[index] +
            c*dt  *
                ((EM.E_x[Get_index(ix, iy + 1)] - EM.E_x[Get_index(ix, iy)]) / delta_y -
                 (EM.E_y[Get_index(ix + 1, iy)] - EM.E_y[Get_index(ix, iy)]) / delta_x);

        EM.B_x[index] = EM.B_x[index] - c*dt *
                                            (EM.E_z[Get_index(ix, iy + 1)] -
                                             EM.E_z[Get_index(ix, iy)]) / delta_y;

        EM.B_y[index] = EM.B_y[index] + c*dt *
                                            (EM.E_z[Get_index(ix + 1, iy)] -
                                             EM.E_z[Get_index(ix, iy)]) / delta_x;
      }
    }
    // Update E, remaining "1/2 timestep"
    for (int iy = 0; iy < ny; iy++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        int index = Get_index(ix, iy);

        EM.E_z[index] =
            EM.E_z[index] +
            c * dt *
                ((EM.B_y[Get_index(ix, iy)] - EM.B_y[Get_index(ix - 1, iy)]) / delta_x -
                 (EM.B_x[Get_index(ix, iy)] - EM.B_x[Get_index(ix, iy - 1)]) / delta_y) - 4*PI*dt * Jz[index];

        EM.E_x[index] = EM.E_x[index] + c * dt  *
                                            (EM.B_z[Get_index(ix, iy)] -
                                             EM.B_z[Get_index(ix, iy - 1)]) / delta_y - 4*PI*dt * ( Jx[Get_index(ix, iy)] + Jx[Get_index(ix, iy)] ) / 2;

        EM.E_y[index] = EM.E_y[index] - c * dt *
                                            (EM.B_z[Get_index(ix, iy)] -
                                             EM.B_z[Get_index(ix - 1, iy)]) / delta_x - 4*PI*dt * ( Jy[Get_index(ix, iy)] + Jy[Get_index(ix, iy+1)] ) / 2;
      }
    }
    
  }
  void Test_nan() {
    int ix = 0;
    for( const auto& x: RFD.RFD_x ) {
      ix++;
      if( std::isnan(x) ) { printf( "nan in RFDx: %lf, particle %d", x, ix ); std::exit(0); }
    }
    for( const auto& x: RFD.RFD_y ) {
      if( std::isnan(x) ) { printf( "nan in RFDy: %lf", x ); std::exit(0); }
    }
    for( const auto& x: RFD.RFD_z ) {
      if( std::isnan(x) ) { printf( "nan in RFDz: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.E_x ) {
      if( std::isnan(x) ) { printf( "nan in E_x: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.E_y ) {
      if( std::isnan(x) ) { printf( "nan in E_y: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.E_z ) {
      if( std::isnan(x) ) { printf( "nan in E_z: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.B_x ) {
      if( std::isnan(x) ) { printf( "nan in B_x: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.B_y ) {
      if( std::isnan(x) ) { printf( "nan in B_y: %lf", x ); std::exit(0); }
    }
    for( const auto& x: EM.B_z ) {
      if( std::isnan(x) ) { printf( "nan in B_z: %lf", x ); std::exit(0); }
    }
    for( const auto& x: electron_pos ) {
      if( std::isnan(x) ) { printf( "nan in e_pos: %lf", x ); std::exit(0); }
    }
    for( const auto& x: positron_pos ) {
      if( std::isnan(x) ) { printf( "nan in p_pos: %lf", x ); std::exit(0); }
      
    }

  }

  void Reset_current() {
    for( int ix = 0; ix < nx * ny; ix++ ) {
      Jx[ix] = 0.0;
      Jy[ix] = 0.0;
      Jz[ix] = 0.0;
    }
  }
  void Reset_charge() {
    for( auto &x: rho_q ) {
      x = 0.0;
    }
  }

  void Iterate_boris()
  {
    Reset_current();
    // 1st half of current, no effect on EM-fields yet
    Interpolate_half_current_boris();
    
    Boris_move_particles( electron_pos, electron_vel );
    Interpolate_half_current_boris();
    // Interpolate 2nd half of contribution to current
    
    FDTD();

    Reset_charge();
    Interpolate_charge_boris();


    //EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(positron_pos);
    //Boris_velocity(positron_pos, positron_vel, EM_at_particles, 1);

    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    Boris_velocity(electron_pos, electron_vel, EM_at_particles, -1);

    //Boris_move_particles( positron_pos, positron_vel );
    
  }
  void Iterate_RFD()
  {
    Reset_current();

    printf("Iterating!\n");
    // 1st half of current, no effect on EM-fields yet
    Interpolate_half_current_RFD_posi();
    Interpolate_half_current_RFD_elec();


    EM_field_matrix EM_at_positrons = Interpolate_EM_at_particles(positron_pos);
    EM_field_matrix EM_at_electrons = Interpolate_EM_at_particles(electron_pos);
    Test_nan();
    
    // Move positrons
    RFD = Calculate_RFD_at_particles(EM_at_positrons, 1);
    Propagate_particles(positron_pos, RFD, dt);
    
    // Move electrons
    RFD = Calculate_RFD_at_particles(EM_at_electrons, -1);
    Propagate_particles(electron_pos, RFD, dt);
    
    // Get rest of current, yielding a current that uses the average shape
    // factor of the particles
    Interpolate_half_current_RFD_posi();
    Interpolate_half_current_RFD_elec();

    // Using currents evaluate fields
    FDTD();

  }

  void Save_current_state(std::string EM_filename,
                          std::string particle_filename,
                          std::string RFD_filename,
                          std::string current_filename,
                          std::string charge_filename ) {
    bool append = false;
    EM.Save(EM_filename, append);
    // printf( "\n\n" );
    // for( int ix = 0; ix < nx * ny; ix++ ) {
    //  printf( "%lf, ", EM.E_z[ix] );
    //}
    // bad fix

    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    //EM_at_particles.Save("data/interpol", append);
    //RFD = Calculate_RFD_at_particles(EM_at_particles, 1);

    //RFD.Save(RFD_filename, append);
    Write_vector_to_binary(particle_filename + "_electron", electron_pos,
                           append);
    Write_vector_to_binary(particle_filename + "_positron", positron_pos,
                           append);
    Write_vector_to_binary(current_filename + "_x", Jx,
                           append);
    Write_vector_to_binary(current_filename + "_y", Jy,
                           append);
    Write_vector_to_binary(current_filename + "_z", Jz,
                           append);
    Write_vector_to_binary(charge_filename, rho_q,
                           append);
  }

  void Append_current_state(std::string EM_filename,
                            std::string particle_filename,
                            std::string RFD_filename,
                            std::string current_filename,
                            std::string charge_filename )
  {
    bool append = true;
    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    EM_at_particles.Save("data/interpol", append);

    EM.Save(EM_filename, append);
    RFD.Save(RFD_filename, append);
    Write_vector_to_binary(particle_filename + "_electron", electron_pos,
                           append);
    Write_vector_to_binary(particle_filename + "_positron", positron_pos,
                           append);
    Write_vector_to_binary(current_filename + "_x", Jx,
                           append);
    Write_vector_to_binary(current_filename + "_y", Jy,
                           append);
    Write_vector_to_binary(current_filename + "_z", Jz,
                           append);
    Write_vector_to_binary(charge_filename, rho_q,
                           append);
      
    Write_vector_to_binary( std::string("./data/e_momenta"), electron_vel, 1 );
  }
};
