#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "classes.hpp"
#include "RFD.hpp"
#include "propagation.hpp"
#include "FDTD.hpp"
#include "solver.hpp"

const double PI = 3.14159265358979;

int run_debug_RFD_function();
int run_debug_particle_propagation();
int run_debug_FDTD();
int run_debug_solver();

int main() {

  /*
  for( int iy = 0; iy < 10; iy++ ) {
    for( int ix = 0; ix < 10; ix++ ) {
      int idx = Get_index( ix, iy + 1, 10, 10 );
      printf("%-3d, ", idx );
    }
    printf("\n");
  }
  */
  // run_debug_RFD_function();
  // run_debug_particle_propagation();
  // run_debug_FDTD();
  run_debug_solver();

  return 0;
}
  

int run_debug_solver() {

  int nx = 82;
  int ny = 82;
  std::vector<double>::size_type n_particles = 1;
  int save_rate = 1;

  const double tmax = 21.0;
  const double n_tsteps = 5000;

  double delta_x = 1.0;
  double delta_y = 1.0;

  std::string EM_filename("./data/EM");
  std::string RFD_filename("./data/RFD");
  std::string particle_filename("./data/particle");

  EM_field_matrix EM_IC(nx, ny); 
  std::vector<double>::size_type num_waves = 2;
  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};
  
  wave_config_init[0][0] = 1.0; // amplitude
  wave_config_init[0][1] = 0.3; // ang_freq
  wave_config_init[0][2] = PI / 2.0; // prop angle
  wave_config_init[0][3] = 0.0; // phase

  wave_config_init[1][0] = 1.0;
  wave_config_init[1][1] = 0.3;
  wave_config_init[1][2] = 0.0;
  wave_config_init[1][3] = 0.0;
  EM_wave_config config(wave_config_init);

  EM_field_matrix all_modes(nx, ny);
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      double x = ix * delta_x;
      double y = iy * delta_y;
      // printf( "(%.2lf, %.2lf)", x, y );

      
      EM_IC.E_x[ix + iy * nx] = 0;
      EM_IC.E_y[ix + iy * nx] = 0;
      EM_IC.E_z[ix + iy * nx] = 0;

      EM_IC.B_x[ix + iy * nx] = 0;
      EM_IC.B_y[ix + iy * nx] = 0;
      EM_IC.B_z[ix + iy * nx] = 1.0;
      /*
      EM_IC.E_x[ix + iy * nx] = Gaussian(x, y);
      EM_IC.E_y[ix + iy * nx] = Gaussian(x, y);
      EM_IC.E_z[ix + iy * nx] = Gaussian(x, y);

      EM_IC.B_x[ix + iy * nx] = Gaussian(x, y);
      EM_IC.B_y[ix + iy * nx] = Gaussian(x, y);
      EM_IC.B_z[ix + iy * nx] = Gaussian(x, y);
      */
      /*
      EM_IC.E_x[ix + iy * nx] = Get_EM_wave_component(0, config, x, y, 0);
      EM_IC.E_y[ix + iy * nx] = Get_EM_wave_component(1, config, x, y, 0);
      EM_IC.E_z[ix + iy * nx] = Get_EM_wave_component(2, config, x, y, 0);
      //printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

      EM_IC.B_x[ix + iy * nx] = Get_EM_wave_component(3, config, x, y, 0);
      EM_IC.B_y[ix + iy * nx] = Get_EM_wave_component(4, config, x, y, 0);
      EM_IC.B_z[ix + iy * nx] = Get_EM_wave_component(5, config, x, y, 0);
      */
    }
    // printf( "\n" );
  }


  Solver mySolver(nx, ny, n_particles, tmax, n_tsteps, save_rate, delta_x,
                  delta_y);
  std::string filename("./config.csv");

  mySolver.Initialize( EM_IC );
  mySolver.Save_parameters_to_text(filename, 0);
  mySolver.Save_current_state(EM_filename, particle_filename, RFD_filename);
  

  for( int tx = 0; tx < n_tsteps; tx++ ) {
    //printf( "tx = %d \n", tx );
    if( tx % save_rate == 0 && tx != 0 ) {
      mySolver.Append_current_state( EM_filename, particle_filename, RFD_filename );
    }
    mySolver.Iterate();
  }

  return 0;
}

int write_vector_to_binary(std::string filename, std::vector<double> vect) {

  std::ofstream filestream(filename, std::ios::out | std::ios::binary);

  filestream.write((char *)&vect[0], vect.size() * sizeof(double));
  filestream.close();

  return 0;
}

// Takes EM matrix object and writes content to file,
// nx columns, ny rows
// adapted from old code, TODO: Refresh memory on how it works
int Write_EM_to_binary(std::string filename_E, std::string filename_B,
                       EM_field_matrix EM_field) {

  // write all components of EM_field to file
  write_vector_to_binary(filename_E + "_x.dat", EM_field.E_x);
  write_vector_to_binary(filename_E + "_y.dat", EM_field.E_y);
  write_vector_to_binary(filename_E + "_z.dat", EM_field.E_z);

  write_vector_to_binary(filename_B + "_x.dat", EM_field.B_x);
  write_vector_to_binary(filename_B + "_y.dat", EM_field.B_y);
  write_vector_to_binary(filename_B + "_z.dat", EM_field.B_z);

  return 0;
}

double Gaussian(double x, double y) { return 1.0 * std::exp(-x * x - y * y); }

double Point_dist(double x, double y) {
  if (x * x + y * y < 1.0)
    return 0.0;
  else
    return 1.0 / (x * x + y * y);
}

int run_debug_FDTD() {
  int nx = 101;
  int ny = nx;
  int save_rate = 10;
  int tmax = 1000;
  std::vector<double>::size_type num_waves = 2;

  double x_max = 3.0;
  double delta_x = (2 * x_max) / (nx - 1);
  double delta_t = 0.005000;
  printf("delta_x, delta_t: %lf, %lf \n", delta_x, delta_t);
  double t = 0;

  double max_ampl = 1.0;
  // double max_ang_freq = 1.0;
  // double max_phi = 2.0 * 3.14;
  // double max_phase = 3.0;

  std::ofstream header_filestream("./data/header.csv");
  std::ofstream point_filestream("./data/time_ev_at_pt.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  // Setup rng
  // unsigned seed =
  // std::chrono::system_clock::now().time_since_epoch().count();
  // std::default_random_engine generator(seed);
  // std::uniform_real_distribution<double> distribution;
  // auto random = std::bind(distribution, generator);

  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};

  /*
  std::ofstream config_filestream("./config.csv");
  config_filestream
      << "Wave, amplitude, ang_freq, propagation angle from +x, phase \n";
  for (size_t ix = 0; ix < num_waves; ix++) {
    wave_config_init[ix][0] = max_ampl;
    wave_config_init[ix][1] = 2.0;
    wave_config_init[ix][2] = 0.78;
    wave_config_init[ix][3] = 0.0;
    config_filestream << ix << ", " << wave_config_init[ix][0] << ", "
                      << wave_config_init[ix][1] << ", "
                      << wave_config_init[ix][2] << ", "
                      << wave_config_init[ix][3] << "\n";
  }
  config_filestream.close();
  */
  wave_config_init[0][0] = max_ampl;
  wave_config_init[0][1] = 2.0;
  wave_config_init[0][2] = -PI / 2.0;
  wave_config_init[0][3] = 0.0;
  wave_config_init[1][0] = max_ampl;
  wave_config_init[1][1] = 2.0;
  wave_config_init[1][2] = 0.0;
  wave_config_init[1][3] = 0.0;
  EM_wave_config config(wave_config_init);
  // FDTD_mode TMz( nx, ny );
  // FDTD_mode TEz( nx, ny );

  // TMz mode is ( Ez, Hz, Hy )
  // TEz mode is ( Hz, Ex, Ey )

  EM_field_matrix all_modes(nx, ny);
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      double x = 2.0 * ix / ((double)nx - 1) * x_max - x_max;
      double y = 2.0 * iy / ((double)ny - 1) * x_max - x_max;
      // printf( "(%.2lf, %.2lf)", x, y );

      /*
      all_modes.E_x[ix + iy * nx] = 1.0;
      all_modes.E_y[ix + iy * nx] = 1.0;
      all_modes.E_z[ix + iy * nx] = 1.0;

      all_modes.B_x[ix + iy * nx] = 1.0;
      all_modes.B_y[ix + iy * nx] = 1.0;
      all_modes.B_z[ix + iy * nx] = 1.0;

      all_modes.E_x[ix + iy * nx] = Gaussian(x, y);
      all_modes.E_y[ix + iy * nx] = Gaussian(x, y);
      all_modes.E_z[ix + iy * nx] = Gaussian(x, y);

      all_modes.B_x[ix + iy * nx] = Gaussian(x, y);
      all_modes.B_y[ix + iy * nx] = Gaussian(x, y);
      all_modes.B_z[ix + iy * nx] = Gaussian(x, y);
      */

      all_modes.E_x[ix + iy * nx] = Get_EM_wave_component(0, config, x, y, t);
      all_modes.E_y[ix + iy * nx] = Get_EM_wave_component(1, config, x, y, t);
      all_modes.E_z[ix + iy * nx] = Get_EM_wave_component(2, config, x, y, t);

      all_modes.B_x[ix + iy * nx] = Get_EM_wave_component(3, config, x, y, t);
      all_modes.B_y[ix + iy * nx] = Get_EM_wave_component(4, config, x, y, t);
      all_modes.B_z[ix + iy * nx] = Get_EM_wave_component(5, config, x, y, t);
    }
    // printf( "\n" );
  }
  const std::string E_filename("./data/E_data");
  const std::string B_filename("./data/B_data");

  const std::string E_history("./data/time/E_data");
  const std::string B_history("./data/time/B_data");

  // all_modes.E_z[Get_index( nx/2, ny/2, nx, ny )] = std::cos( PI/10.0 * t );

  Write_EM_to_binary(E_filename + "before", B_filename + "before", all_modes);
  for (int tx = 0; tx < tmax; tx++) {
    // all_modes.E_z[Get_index( nx/2, ny/2, nx, ny )] = std::cos( PI/10.0 * t );
    Update_mode(all_modes, delta_t, delta_x);
    t = t + delta_t;
    if (tx % save_rate == 0) {
      Write_EM_to_binary(E_history + std::to_string(t),
                         B_history + std::to_string(t), all_modes);
      point_filestream << all_modes.E_z[nx / 2 + nx * ny / 2] << ", " << t
                       << "\n";
    }
  }
  Write_EM_to_binary(E_filename + "after", B_filename + "after", all_modes);

  return 0;
}

int run_debug_particle_propagation() {

  int nx = 64;
  int ny = 64;
  std::vector<double>::size_type n_particles = 100;
  int pp = 1;
  int sign_plus = 1;

  std::vector<double>::size_type num_waves = 4;
  double max_ampl = 1.0;
  double max_ang_freq = 1.0;
  double max_phi = 2.0 * 3.14;
  double max_phase = 3.0;

  // double c = 1.0;

  double x_max = 12.0;
  double y_max = 12.0;
  double z_max = 2.0;
  const double tmax = 20.0;
  const double dt = 0.01;

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("./data/header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  const std::string E_filename("./data/time/E_data");
  const std::string B_filename("./data/time/B_data");
  const std::string xpos_filename("./data/time/pos_x");
  const std::string ypos_filename("./data/time/pos_y");
  const std::string zpos_filename("./data/time/pos_z");
  const std::string filename_RFD("./data/time/RFD");

  std::vector<std::vector<double>> particles{n_particles,
                                             std::vector<double>(3)};
  // Setup rng
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution;
  auto random = std::bind(distribution, generator);

  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};

  std::ofstream config_filestream("./config.csv");
  config_filestream
      << "Wave nr, amplitude, ang_freq, propagation angle from +x, phase \n";
  for (size_t ix = 0; ix < num_waves; ix++) {
    wave_config_init[ix][0] = random() * max_ampl;
    wave_config_init[ix][1] = random() * max_ang_freq;
    wave_config_init[ix][2] = random() * max_phi;
    wave_config_init[ix][3] = random() * max_phase;
    config_filestream << ix << ", " << wave_config_init[ix][0] << ", "
                      << wave_config_init[ix][1] << ", "
                      << wave_config_init[ix][2] << ", "
                      << wave_config_init[ix][3] << "\n";
  }
  config_filestream.close();
  EM_wave_config config(wave_config_init);

  for (size_t ix = 0; ix < n_particles; ix++) {
    particles[ix][0] = 2.0 * random() * x_max - x_max;
    particles[ix][1] = 2.0 * random() * y_max - y_max;
    particles[ix][2] = 2.0 * random() * z_max - z_max;
  }

  EM_field_matrix field_at_particles(n_particles, pp);
  EM_field_matrix EM_field(nx, ny);

  double t = 0.0;
  // Main simulation loop
  for (int tx = 0; t < tmax; t = t + dt, tx++) {

    // For visualizing the fields
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {

        std::vector<double> coords = {
            2.0 * iy / (double)(ny - 1) * x_max - x_max,
            2.0 * ix / (double)(nx - 1) * y_max - y_max, 0};
        // printf( "%lf, %lf, %lf \n", coords[0] - coords[1],
        // coords[1], coords[2] );

        EM_field.E_x[ix * ny + iy] =
            Get_EM_wave_component(0, config, coords[0], coords[1], t);
        EM_field.E_y[ix * ny + iy] =
            Get_EM_wave_component(1, config, coords[0], coords[1], t);
        EM_field.E_z[ix * ny + iy] =
            Get_EM_wave_component(2, config, coords[0], coords[1], t);

        EM_field.B_x[ix * ny + iy] =
            Get_EM_wave_component(3, config, coords[0], coords[1], t);
        EM_field.B_y[ix * ny + iy] =
            Get_EM_wave_component(4, config, coords[0], coords[1], t);
        EM_field.B_z[ix * ny + iy] =
            Get_EM_wave_component(5, config, coords[0], coords[1], t);
      }
    }
    Get_EM_field_from_positions(n_particles, particles, field_at_particles, t,
                                config);
    RFD_matrix RFD_at_particles(field_at_particles, sign_plus);

    Propagate_particles(particles, RFD_at_particles, dt);

    if (tx % 10 == 0) {
      // Temporary solution for storing particle positions
      std::vector<double> vect(n_particles);
      for (size_t ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][0];
      }
      write_vector_to_binary(xpos_filename + std::to_string(t), vect);
      for (size_t ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][1];
      }
      write_vector_to_binary(ypos_filename + std::to_string(t), vect);
      for (size_t ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][2];
      }
      write_vector_to_binary(zpos_filename + std::to_string(t), vect);
      Write_EM_to_binary(E_filename + std::to_string(t),
                         B_filename + std::to_string(t), EM_field);
      std::string pdense_filename("./data/time/power_dense");
      std::vector<double> p_dense(nx * ny);
      for (int kx = 0; kx < nx * ny; kx++) {
        double E_squared = EM_field.E_x[kx] * EM_field.E_x[kx] +
                           EM_field.E_y[kx] * EM_field.E_y[kx] +
                           EM_field.E_z[kx] * EM_field.E_z[kx];

        double B_squared = EM_field.B_x[kx] * EM_field.B_x[kx] +
                           EM_field.B_y[kx] * EM_field.B_y[kx] +
                           EM_field.B_z[kx] * EM_field.B_z[kx];

        p_dense[kx] = E_squared + B_squared;
      }

      write_vector_to_binary(pdense_filename + std::to_string(t), p_dense);

      write_vector_to_binary(filename_RFD + "_x" + std::to_string(t),
                             RFD_at_particles.RFD_x);
      write_vector_to_binary(filename_RFD + "_y" + std::to_string(t),
                             RFD_at_particles.RFD_y);
      write_vector_to_binary(filename_RFD + "_z" + std::to_string(t),
                             RFD_at_particles.RFD_z);
    }
  }

  return 0;
}

int run_debug_RFD_function() {
  int nx = 37;
  int ny = 37;
  double EB_max = 3.0;
  int sign = 1;
  EM_field_matrix EM_field(nx, ny);

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("./data/header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  std::string E_filename("./data/E_data");
  std::string B_filename("./data/B_data");

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      EM_field.E_x[ix * ny + iy] =
          2.0 * (iy / (double)(ny - 1)) * EB_max - EB_max;
      EM_field.E_y[ix * ny + iy] =
          2.0 * (ix / (double)(nx - 1)) * EB_max - EB_max;
      EM_field.B_x[ix * ny + iy] = 1.0;
    }
  }

  Write_EM_to_binary(E_filename, B_filename, EM_field);

  RFD_matrix RFD(EM_field, sign);

  std::vector<double> magnitudes(nx * ny);

  for (int ix = 0; ix < nx * ny; ix++) {
    magnitudes[ix] = std::sqrt(RFD.RFD_x[ix] * RFD.RFD_x[ix] +
                               RFD.RFD_y[ix] * RFD.RFD_y[ix] +
                               RFD.RFD_z[ix] * RFD.RFD_z[ix]);
  }

  std::string filename_RFD = "./data/RFD";
  write_vector_to_binary(filename_RFD + "_x.dat", RFD.RFD_x);
  write_vector_to_binary(filename_RFD + "_y.dat", RFD.RFD_y);
  write_vector_to_binary(filename_RFD + "_z.dat", RFD.RFD_z);
  write_vector_to_binary("./data/magnitudes.dat", magnitudes);

  return 0;
}
