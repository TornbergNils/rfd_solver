#ifndef SOLVER_H
#define SOLVER_H

#include "generate_IC.hpp"
#include "Experiment_slab.hpp"
#include "EM.hpp"
/*
  Main solver class for simulations using the particle in cell
  method with either the Boris pusher or the RFD for moving particles.

  To iterate the state of the solver call either the Iterate_boris or the Iterate_RFD
  command to advance the state of the simulation by one timestep.

  Initialize must be called before running any simulation in order to load 
  the initial conditions and calculate derived parameters.


*/

class Solver
{
  // Given parameters
  int nx;
  int ny;
  int weight;
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
  double q_e;
  double m_e;
  
  // Derived parameters
  double tmax;


  // Data
  std::vector<double> electron_pos;
  std::vector<double> positron_pos;
  std::vector<double> electron_vel;
  std::vector<double> positron_vel;
  std::vector<double> electron_gamma;
  std::vector<double> positron_gamma;

  std::vector<double> Jx;
  std::vector<double> Jy;
  std::vector<double> Jz;
  
  std::vector<double> rho_q;
  std::vector<double> RFD_frac;

  EM_field_matrix EM;
  RFD_matrix RFD;

public:
  // Constructor. Initializes most parameters via initializer list and some using
  // a map.
  Solver(int nx, int ny, int weight, int n_particles, double dt, int n_tsteps, int save_rate,
         double delta_x, double delta_y, std::map< std::string, double > ic_param )
      : nx(nx), ny(ny), n_particles(n_particles), n_tsteps(n_tsteps), 
        save_rate(save_rate), dt(dt), n_elec(n_particles),
        n_posi(n_particles), delta_x(delta_x), delta_y(delta_y),
        electron_pos(n_particles * 3), positron_pos(n_particles * 3),
        electron_vel(n_particles * 3), positron_vel(n_particles * 3),
        electron_gamma(n_particles), positron_gamma(n_particles),
        Jx(nx * ny), Jy(nx * ny), Jz(nx * ny), rho_q(nx * ny), RFD_frac(nx * ny),
        EM(nx, ny), RFD(EM, 1) {

          c = ic_param["c"];
          q_e = ic_param["q_e"] * ic_param["weight"];
          m_e = ic_param["m_e"] * ic_param["weight"];

        }
  Solver( const IC_struct& IC ) : nx( IC.nx ), ny( IC.ny ), weight( IC.weight),
        n_particles( IC.n_particles ), n_tsteps( IC.n_tsteps ), 
        save_rate( IC.save_rate ), dt( IC.dt ),
        n_elec( IC.n_particles ), n_posi( IC.n_particles ),
        delta_x( IC.delta_x ), delta_y( IC.delta_y ),
        electron_pos( IC.e_pos_ic ),
        positron_pos( IC.p_pos_ic ),
        electron_vel( IC.e_vel_ic ),
        positron_vel( IC.p_vel_ic ),
        electron_gamma( IC.e_gamma_ic ),
        positron_gamma( IC.p_gamma_ic ),
        Jx( IC.nx * IC.ny ), // zero.initialized
        Jy( IC.nx * IC.ny ),
        Jz( IC.nx * IC.ny ),
        rho_q( IC.nx * IC.ny ),
        EM( IC.EM_ic ), // Copy constructor
        RFD( EM, 1)  {
          c = IC.c;
          q_e = IC.q_e;
          m_e = IC.m_e;
        }

  // Interpolates the charge density using the particle in cell method
  // Counts both contribution from electrons and positrons

  // TODO: add parameter that regulates inclusion of positrons/electrons
  void Interpolate_charge_boris()
  {
    // For each electron, round down position to get grid position
    // J is co-located with Ez
    double sign = -1 * q_e / ( delta_x * delta_y ); // * q_e_cgs;
    
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
    sign = 1 * q_e / ( delta_x * delta_y ); // * q_e_cgs;
    
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

      rho_q[Get_index(ix, iy)] += sign * w00;
      rho_q[Get_index(ix+1, iy)] += sign * w10;
      rho_q[Get_index(ix+1, iy+1)] += sign * w11;
      rho_q[Get_index(ix, iy+1)] += sign * w01;
      
    }
  }

  // Interpolates one of Jx/Jy/Jz using arguments
  void Interpolate_current_component(std::vector<double> &current, double vel, int ix,
                                       int iy, double w00, double w10, double w11, double w01, double sign )
  {
      //printf("half current at pt: %2.2e \n", sign * vel / ( gamma * 2 )  );
      current[Get_index(ix, iy)] += sign * w00 * vel / 2;
      current[Get_index(ix+1, iy)] += sign * w10 * vel / 2;
      current[Get_index(ix+1, iy+1)] += sign * w11 * vel / 2;
      current[Get_index(ix, iy+1)] += sign * w01 * vel / 2;
  }
  
  /* Interpolates the current density using the particle in cell method
   Counts both contribution from electrons and positrons

   In order to use the average of the weights before and after the 
   particles are moved, only half the contribution is taken and
   the function needs to be used before and after particle movement.

   TODO: add parameter that regulates inclusion of positrons/electrons
  */
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

      Interpolate_current_component( Jx, electron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, sign );
      Interpolate_current_component( Jy, electron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, sign );
      Interpolate_current_component( Jz, electron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, sign );
    }
    
    // repeat for positrons
    sign = 1 * q_e / (delta_x * delta_y ); // * q_e_cgs;
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

      Interpolate_current_component( Jx, positron_vel[ip], ix,
                                       iy, w00, w10, w11, w01, sign );
      Interpolate_current_component( Jy, positron_vel[ip+1], ix,
                                       iy, w00, w10, w11, w01, sign );
      Interpolate_current_component( Jz, positron_vel[ip+2], ix,
                                       iy, w00, w10, w11, w01, sign );
    }
      
  }
  /* Interpolates the current density using the particle in cell method
   
   Takes the contribution to total current from the species whose
   positions are given as argument. 

   Sign should be +1 for positrons and -1 for electrons.

   The function needs to be used before and after particle movement
   to use the average of the particle shape functions

  */
  void Interpolate_half_current_RFD(std::vector<double> position, int sign ) {

    // For each particle, round down position to get grid position
    // J is co-located with Ez
    double rho_per_cell = sign * q_e / ( delta_x * delta_y );
    
    for (long unsigned ip = 0; ip < n_elec * 3; ip += 3)
    {
      double x = position[ip];
      double y = position[ip + 1];
      int ix = std::floor(x / delta_x);
      int iy = std::floor(y / delta_y);
      x = x - ix * delta_x;
      y = y - iy * delta_y;

      double cell_size = delta_x * delta_y;

      double w00 = (delta_x - x) * (delta_y - y) / cell_size;
      double w10 = x * (delta_y - y) / cell_size;
      double w11 = x * y / cell_size;
      double w01 = (delta_x - x) * y / cell_size;

      Interpolate_current_component( Jx, RFD.RFD_x[ip/3] * c, ix,
                                       iy, w00, w10, w11, w01, rho_per_cell );
      Interpolate_current_component( Jy, RFD.RFD_y[ip/3] * c, ix,
                                       iy, w00, w10, w11, w01, rho_per_cell );
      Interpolate_current_component( Jz, RFD.RFD_z[ip/3] * c, ix,
                                       iy, w00, w10, w11, w01, rho_per_cell );
    }
  }

  // Transfer initial conditions from IC struct to solver
  // Set up any actions to be taken before iterating.
  void Initialize( const IC_struct& IC  )
  {
    tmax = dt * n_tsteps;
    
    EM = IC.EM_ic;
    electron_pos = IC.e_pos_ic;
    positron_pos = IC.p_pos_ic;

    electron_vel = IC.e_vel_ic;
    positron_vel = IC.p_vel_ic;

    electron_gamma = IC.e_gamma_ic;
    positron_gamma = IC.p_gamma_ic;

    EM_field_matrix mat( n_elec, 1);
    RFD_matrix temp( mat, 1);
    RFD = temp;
    
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
                 << "weight, " << weight << "\n"
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
                 << "weight, "
                 << "n_particles, "
                 << "tmax, "
                 << "n_tsteps, "
                 << "save_rate, "
                 << "n_elec, "
                 << "n_posi, "
                 << "delta_x, "
                 << "delta_y, "
                 << "dt \n"
                 << nx << ", " << ny << ", " << weight << ", " << n_particles << ", " << tmax
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

// Propagate particles along RFD at velocity c
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
  
  int Propagate_hybrid(std::vector<double> &particles, std::vector<double> &velocities,
                          const RFD_matrix &RFD_ap, const double dt, const double RFD_frac)
  {
    int ip_max = particles.size();
    if (RFD_ap.RFD_x.size() != ip_max / 3)
    {
      printf(" RFD vector and particle vector size mismatch!");
    }
    for (int ip = 0, irf = 0; ip < ip_max; ip += 3, irf++)
    {
      particles[ip] +=     (RFD_frac * RFD_ap.RFD_x[irf] * c 
        + (1.0 - RFD_frac) * velocities[ip]) * dt;

      particles[ip + 1] += (RFD_frac * RFD_ap.RFD_y[irf] * c 
        + (1.0 - RFD_frac) * velocities[ip + 1] ) * dt;

      particles[ip + 2] += (RFD_frac * RFD_ap.RFD_z[irf] * c
        + (1.0 - RFD_frac) * velocities[ip + 2] ) * dt;
      
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

// Calculate particle velocity and gamma factor according to Boris algorithm
  int Boris_velocity(std::vector<double> &particles,
                                std::vector<double> &vel,
                                std::vector<double> &gamma_vect,
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
      std::vector<double> u{vel[ip] * gamma_vect[iem],
                            vel[ip + 1] * gamma_vect[iem],
                            vel[ip + 2] * gamma_vect[iem]};
      std::vector<double> v(3);

      // Add half impulse from E-field to get u- in v
      v[0] = u[0] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_x[iem];  //prop_factor * EM_ap.E_x[iem];
      v[1] = u[1] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_y[iem];  //prop_factor * EM_ap.E_y[iem];
      v[2] = u[2] + ( sign * q_e * dt ) / ( 2 * m_e ) * EM_ap.E_z[iem];  //prop_factor * EM_ap.E_z[iem];

      double u_minus_squared = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      // Recalculate gamma at different time
      double gamma = std::sqrt(1.0 + u_minus_squared/(c*c));
      
      //if( gamma > 1000 ) {
      //  printf( "gamma1: %lf \n", gamma );
      //  }
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

      double u_now_squared = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
      gamma = std::sqrt(1.0 + u_now_squared/(c*c));
      //printf(", gamma2: %lf \n", gamma);
      // Finally update velocity vector
      vel[ip] = u[0] / gamma;
      vel[ip + 1] = u[1] / gamma;
      vel[ip + 2] = u[2] / gamma;
       
      gamma_vect[iem] = gamma;
      }


    return 0;
  }
// Moves particles and implements periodic boundary conditions
  void Boris_move_particles( std::vector<double> &pos, std::vector<double> &vel ) {
    
    int ip_max = electron_pos.size();
    for (int ip = 0; ip < ip_max; ip += 3 )
    {

      pos[ip] += vel[ip] * dt;
      pos[ip + 1] += vel[ip + 1] * dt;
      pos[ip + 2] += vel[ip + 2] * dt;

      // Periodic Bc for particles
      pos[ip] =      std::fmod( pos[ip] + nx * delta_x, nx * delta_x );
      pos[ip + 1] =  std::fmod( pos[ip+1] + ny * delta_y, ny * delta_y );
    }
  }

// Implements periodic boundary conditions for grid properties via indexing
  inline int Get_index(int ix, int iy)
  {
    return ((ix + nx) % nx) + ((iy + ny) % ny) * nx;
  }

// Interpolates any field component using particle in cell method
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

// Returns a object containing all interpolated 
// field components at all particle positions
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

// Returns a object containing the RFD vector at every particle position
  RFD_matrix Calculate_RFD_at_particles(const EM_field_matrix &EM, int sign)
  {
    RFD_matrix RFD_temp(EM, sign);
    return RFD_temp;
  }

// Propagates the field at the grid points using the FDTD
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
                                             EM.B_z[Get_index(ix, iy - 1)]) / delta_y - 4*PI*dt * ( Jx[Get_index(ix, iy)] + Jx[Get_index(ix+1, iy)] ) / 2;

        EM.E_y[index] = EM.E_y[index] - c * dt *
                                            (EM.B_z[Get_index(ix, iy)] -
                                             EM.B_z[Get_index(ix - 1, iy)]) / delta_x - 4*PI*dt * ( Jy[Get_index(ix, iy)] + Jy[Get_index(ix, iy+1)] ) / 2;
      }
    }
    
  }
  
  void Calculate_RFD_fraction( double omega ) {
    double ss_field = 4.413*1e13;
    double compton_freq = 7.7807*1e20;
    double om_by_compton = omega*omega/(compton_freq*compton_freq);

    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        double index = Get_index(ix, iy );
        double Ex = EM.E_x[index];
        double Ey = EM.E_y[index];
        double Ez = EM.E_z[index];

        double EM_abs = std::sqrt( Ex*Ex + Ey*Ey + Ez*Ez );
        RFD_frac[index] = ( EM_abs / ss_field ) * om_by_compton > 3e8 ? 1 : 0;
      }
    }
    
  }

// Iterates through all
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

// Sets current to 0
  void Reset_current() {
    for( int ix = 0; ix < nx * ny; ix++ ) {
      Jx[ix] = 0.0;
      Jy[ix] = 0.0;
      Jz[ix] = 0.0;
    }
  }

// Sets charge to 0
  void Reset_charge() {
    for( auto &x: rho_q ) {
      x = 0.0;
    }
  }

// Advances the system in time by dt using PIC with the Boris pusher
  void Iterate_boris()
  {
    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(positron_pos);
    Boris_velocity(positron_pos, positron_vel, positron_gamma, EM_at_particles, 1);

    EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    Boris_velocity(electron_pos, electron_vel, electron_gamma, EM_at_particles, -1);
    
    Reset_current();
    // 1st half of current, no effect on EM-fields yet
    Interpolate_half_current_boris();
    
    Boris_move_particles( electron_pos, electron_vel );
    Boris_move_particles( positron_pos, positron_vel );
    Interpolate_half_current_boris();
    // Interpolate 2nd half of contribution to current

    Reset_charge();
    Interpolate_charge_boris();

    FDTD();
    
  }
/*
  Advances the system in time by dt using PIC with RFD pusher
*/
  void Iterate_RFD()
  {
    Reset_current();
    Test_nan();

    // Move positrons and interpolate positron current
    RFD_move_particles_and_interpolJ( positron_pos, 1 );
    
    // Move electrons and interpolate positron current
    RFD_move_particles_and_interpolJ( electron_pos, -1 );
    
    // Update interpolated charge density, not used but
    // useful diagnostic
    Reset_charge();
    Interpolate_charge_boris();

    // Using currents update fields using FDTD scheme
    FDTD();

  }

  void Iterate_hybrid() {
    
    Reset_current();
    Reset_charge();
    double omega = 0.0;
    
    // Calculate velocity using both boris and RFD
    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(positron_pos);
    Boris_velocity(positron_pos, positron_vel, positron_gamma, EM_at_particles, 1);
    RFD = Calculate_RFD_at_particles(EM_at_particles, 1);
    
    Interpolate_half_current_boris();
    Propagate_hybrid( positron_pos, positron_vel, RFD, dt, RFD_frac );
    Interpolate_half_current_boris();
    
    // Do it again for electrons
    EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    Boris_velocity(electron_pos, electron_vel, electron_gamma, EM_at_particles, 1);
    RFD = Calculate_RFD_at_particles(EM_at_particles, 1);


    Interpolate_half_current_boris();
    Propagate_hybrid( electron_pos, electron_vel, RFD, dt, RFD_frac );
    Interpolate_half_current_boris();

    // New positron velocity is then (1-RFD_fraction) * boris_vel 
    // + RFD_fraction * RFD*c
    

    Interpolate_charge_boris();

    // Using currents update fields using FDTD scheme
    FDTD();

  }

  void RFD_move_particles_and_interpolJ(std::vector<double> &positions,
                                        int sign ) {
    
    EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(positions);
    RFD = Calculate_RFD_at_particles(EM_at_particles, sign);
    
    // NOTE: as interpolate uses current RFD solver var, it
    // is VERY important functions are executed in the correct order
    Interpolate_half_current_RFD(positions, sign);
    Propagate_particles(positions, RFD, dt);
    Interpolate_half_current_RFD(positions, sign);

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

    //EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(electron_pos);
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
    Write_vector_to_binary( particle_filename + "_e_velocities", electron_vel, append );
  }

  void Append_current_state(std::string EM_filename,
                            std::string particle_filename,
                            std::string RFD_filename,
                            std::string current_filename,
                            std::string charge_filename )
  {
    bool append = true;
    //EM_field_matrix EM_at_particles = Interpolate_EM_at_particles(electron_pos);
    //EM_at_particles.Save("data/interpol", append);

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
    Write_vector_to_binary( particle_filename + "_e_velocities", electron_vel, append );
  }
};

#endif // SOLVER_H
