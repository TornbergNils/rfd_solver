
class Solver
{
  // Given parameters
  int nx;
  int ny;
  int n_particles;
  int tmax;
  int n_tsteps;
  int save_rate;
  std::vector<double>::size_type n_elec;
  std::vector<double>::size_type n_posi;

  double delta_x;
  double delta_y;

  // Derived parameters
  double dt;

  // Data
  std::vector<double> electron_pos;
  std::vector<double> positron_pos;
  std::vector<double> electron_vel;
  std::vector<double> positron_vel;

  std::vector<double> Jx;
  std::vector<double> Jy;
  std::vector<double> Jz;

  EM_field_matrix EM;
  RFD_matrix RFD;

public:
  Solver(int nx, int ny, int n_particles, int tmax, int n_tsteps, int save_rate,
         double delta_x, double delta_y)
      : nx(nx), ny(ny), n_particles(n_particles), tmax(tmax),
        n_tsteps(n_tsteps), save_rate(save_rate), n_elec(n_particles),
        n_posi(n_particles), delta_x(delta_x), delta_y(delta_y),
        electron_pos(n_particles * 3), positron_pos(n_particles * 3),
        electron_vel(n_particles * 3), positron_vel(n_particles * 3),
        Jx(nx * ny), Jy(nx * ny), Jz(nx * ny),
        EM(nx, ny), RFD(EM, 1) {}

  void Interpolate_current_component(std::vector<double> &current, double vel, int ix,
                                       int iy, double w00, double w10, double w11, double w01, double gamma, int sign )
  {
      current[Get_index(ix, iy)] += sign * w00 * vel / ( gamma * 2 );
      current[Get_index(ix+1, iy)] += sign * w10 * vel / ( gamma * 2 );
      current[Get_index(ix+1, iy+1)] += sign * w11 * vel / ( gamma * 2 );
      current[Get_index(ix, iy+1)] += sign * w01 * vel / ( gamma * 2 );
  }
  void Interpolate_half_current()
  {
    // For each electron, round down position to get grid position
    // J is co-located with Ez
    int sign = -1;
    
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

      double gamma = std::sqrt(1.0 + vel_squared);
      
      Interpolate_current_component( Jx, electron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jy, electron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jz, electron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
    }
    
    // repeat for positrons
    /*
    sign = 1;
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
      double gamma = std::sqrt(1.0 + vel_squared);

      Interpolate_current_component( Jx, positron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jy, positron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
      Interpolate_current_component( Jz, positron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, gamma, sign );
    }
  */    
  }
  void Initialize(EM_field_matrix EM_IC)
  {
    // Setup rng
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution;
    auto random = std::bind(distribution, generator);

    dt = (double)tmax / (double)n_tsteps;
    double x_len = nx * delta_x;
    double y_len = ny * delta_y;

    // set EM to equal to EM_IC, since EM is 0-initialized by constructor
    EM = EM_IC;
    RFD.Update(EM, 1);
    
    // Constants for calculating system properties and debugging
    const double PI = 3.14159265358979;
    const double Me = 9.1093819*1e-31;
    const double Kb = 1.380650*1e-23;
    const double Me_by_Kb = 6.59789976*1e-8;
    const double c_SI = 2.99792458 * 1e8;
    double electron_temp = 0.0;
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

      electron_vel[ip] = 0.0;
      electron_vel[ip + 1] =  0.5 + 0.1 * std::sin( electron_pos[ip+1] / 100 );
      electron_vel[ip + 2] = 0.0;
      
      double vel_squared_e = electron_vel[ip] * electron_vel[ip]
      + electron_vel[ip+1] * electron_vel[ip+1]
      + electron_vel[ip+2] * electron_vel[ip+2];
      gamma_e = 1.0 / std::sqrt(1.0 - vel_squared_e);
      
      if( vel_squared_e > 1.0 ) {printf( "e particle speed exceeds c! \n" );}
      electron_vel[ip] *= gamma_e;
      electron_vel[ip + 1] *= gamma_e; 
      electron_vel[ip + 2] *= gamma_e;

      positron_vel[ip] = 0.0;
      positron_vel[ip + 1] = -0.5 - 0.1 * std::sin( positron_pos[ip+1] / 100 );
      positron_vel[ip + 2] = 0.0;
      
      double vel_squared_p = positron_vel[ip] * positron_vel[ip]
      + positron_vel[ip+1] * positron_vel[ip+1]
      + positron_vel[ip+2] * positron_vel[ip+2];
      double gamma_p = 1.0 / std::sqrt(1.0 - vel_squared_p);
      if( vel_squared_p > 1.0 ) {printf( "p particle speed exceeds c! \n" );}
      
      positron_vel[ip] *= gamma_p;
      positron_vel[ip + 1] *= gamma_p; 
      positron_vel[ip + 2] *= gamma_p;

      electron_temp += (gamma_e -1 ) * (c_SI * c_SI ) * ( Me_by_Kb ) / (3*n_particles);

    }
    double Debye_length = 0;
    
    
    printf("Ratio: %lf \n", ( 0.14 * c_SI * c_SI ) * ( Me / Kb ) );
    printf("Electron temp: %.5e \n", electron_temp);
    // Wonky fix to show approx current on first iter
    Interpolate_half_current();
    Interpolate_half_current();
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
      particles[ip] += RFD_ap.RFD_x[irf] * dt;
      particles[ip + 1] += RFD_ap.RFD_y[irf] * dt;
      particles[ip + 2] += RFD_ap.RFD_z[irf] * dt;
      
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

  int Boris_propagate_particles(std::vector<double> &particles,
                                std::vector<double> &vel,
                                const EM_field_matrix &EM_ap,
                                const int sign)
  {

    // set unit s.t -q / m_e = 1
    double q_div_by_m = sign * 1.0;
    const double prop_factor = q_div_by_m / 2.0 * dt;
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
      v[0] = u[0] + prop_factor * EM_ap.E_x[iem];
      v[1] = u[1] + prop_factor * EM_ap.E_y[iem];
      v[2] = u[2] + prop_factor * EM_ap.E_z[iem];

      double u_minus_squared = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      // Recalculate gamma at different time
      double gamma = std::sqrt(1.0 + u_minus_squared);
      // printf( "gamma1: %lf \n", gamma );
      //  Get t vector and put it in u
      u[0] = prop_factor / gamma * EM_ap.B_x[iem];
      u[1] = prop_factor / gamma * EM_ap.B_y[iem];
      u[2] = prop_factor / gamma * EM_ap.B_z[iem];

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
      u[0] = u[0] + prop_factor * EM_ap.E_x[iem];
      u[1] = u[1] + prop_factor * EM_ap.E_y[iem];
      u[2] = u[2] + prop_factor * EM_ap.E_z[iem];

      double u_now_squared = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      gamma = std::sqrt(1.0 + u_now_squared);
      //printf(", gamma2: %lf \n", gamma);
      // Finally update vectors
      vel[ip] = u[0];
      vel[ip + 1] = u[1];
      vel[ip + 2] = u[2];

      // note +=
      // Note also particles and vel known at offset times!
      particles[ip] += u[0] * dt / gamma;
      particles[ip + 1] += u[1] * dt / gamma;
      particles[ip + 2] += u[2] * dt / gamma;

      // Periodic Bc for particles
      particles[ip] =      std::fmod( particles[ip] + nx * delta_x, nx * delta_x );
      particles[ip + 1] =  std::fmod( particles[ip+1] + ny * delta_y, ny * delta_y );
      }

    /*
    int ip_max = particles.size();
    if (RFD_ap.RFD_x.size() != ip_max / 3) {
      printf(" RFD vector and particle vector size mismatch!");
    }
    for (int ip = 0, irf = 0; ip < ip_max; ip += 3, irf++) {
      particles[ip] += RFD_ap.RFD_x[irf] * dt;
      particles[ip + 1] += RFD_ap.RFD_y[irf] * dt;
      particles[ip + 2] += RFD_ap.RFD_z[irf] * dt;
    }
    */
    return 0;
  }

  inline int Get_index(int ix, int iy)
  {
    return ((ix + nx) % nx) + ((iy + ny) % ny) * nx;
  }


  double Interpolate_field_component(const std::vector<double> &field, int ix,
                                     int iy, double x, double y)
  {
    if (x > 1.0 || y > 1.0)
    {
      printf("Error! Interpolation coordinate sizes exceed 1!");
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
  };

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

      int ExHy_ix = std::floor(xp / delta_x - delta_x / 2);
      int ExHy_iy = std::floor(yp / delta_y);

      int EyHx_ix = std::floor(xp / delta_x);
      int EyHx_iy = std::floor(yp / delta_y - delta_y / 2);

      int Hz_ix = std::floor(xp / delta_x - delta_x / 2);
      int Hz_iy = std::floor(yp / delta_y - delta_y / 2);
      // double z = particle_pos[ip+2];
      // int iz = std::floor( z );

      // need new x,y relative to lower left corner of box
      // Ex and Hy are displaced by delta_x / 2 in x-direction
      double x = xp - ExHy_ix * delta_x - delta_x / 2;
      double y = yp - ExHy_iy * delta_y;

      interpolated_field.E_x[ip / 3] =
          Interpolate_field_component(EM.E_x, ExHy_ix, ExHy_iy, x, y);

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
            dt / delta_x *
                ((EM.E_x[Get_index(ix, iy + 1)] - EM.E_x[Get_index(ix, iy)]) -
                 (EM.E_y[Get_index(ix + 1, iy)] - EM.E_y[Get_index(ix, iy)]));

        EM.B_x[index] = EM.B_x[index] - dt / delta_x *
                                            (EM.E_z[Get_index(ix, iy + 1)] -
                                             EM.E_z[Get_index(ix, iy)]);

        EM.B_y[index] = EM.B_y[index] + dt / delta_x *
                                            (EM.E_z[Get_index(ix + 1, iy)] -
                                             EM.E_z[Get_index(ix, iy)]);
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
            dt / delta_x *
                ((EM.B_y[Get_index(ix, iy)] - EM.B_y[Get_index(ix - 1, iy)]) -
                 (EM.B_x[Get_index(ix, iy)] - EM.B_x[Get_index(ix, iy - 1)])) - dt * Jz[index];

        EM.E_x[index] = EM.E_x[index] + dt / delta_x *
                                            (EM.B_z[Get_index(ix, iy)] -
                                             EM.B_z[Get_index(ix, iy - 1)]) - dt * ( Jx[Get_index(ix, iy)] + Jx[Get_index(ix+1, iy)] ) / 2;

        EM.E_y[index] = EM.E_y[index] - dt / delta_x *
                                            (EM.B_z[Get_index(ix, iy)] -
                                             EM.B_z[Get_index(ix - 1, iy)]) - dt * ( Jy[Get_index(ix, iy)] + Jy[Get_index(ix, iy+1)] ) / 2;
      }
    }
  }

  void Reset_current() {
    for( int ix = 0; ix < nx * ny; ix++ ) {
      Jx[ix] = 0.0;
      Jy[ix] = 0.0;
      Jz[ix] = 0.0;
    }
  }

  void Iterate()
  {
    Reset_current();
    Interpolate_half_current();
    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(positron_pos);
    //RFD = Calculate_RFD_at_particles(EM_at_particles, 1);
    //Propagate_particles(positron_pos, RFD, dt);
    Boris_propagate_particles(positron_pos, positron_vel, EM_at_particles, 1);

    EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    //RFD = Calculate_RFD_at_particles(EM_at_particles, -1);
    //Propagate_particles(electron_pos, RFD, dt);
    Boris_propagate_particles(electron_pos, electron_vel, EM_at_particles, -1);

    // Interpolate 2nd half of contribution to current
    Interpolate_half_current();
    FDTD();
  }

  void Save_current_state(std::string EM_filename,
                          std::string particle_filename,
                          std::string RFD_filename,
                          std::string current_filename)
  {
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
  }

  void Append_current_state(std::string EM_filename,
                            std::string particle_filename,
                            std::string RFD_filename,
                            std::string current_filename)
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
  }
};
