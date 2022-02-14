#include "classes.hpp"
#include "RFD.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int run_debug_RFD_function();
int run_debug_particle_propagation();

int main() {

  run_debug_RFD_function();
  //run_debug_particle_propagation();

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
                       EM_field_matrix EM_field ) {

  // write all components of EM_field to file
  write_vector_to_binary(filename_E + "_x.dat", EM_field.E_x);
  write_vector_to_binary(filename_E + "_y.dat", EM_field.E_y);
  write_vector_to_binary(filename_E + "_z.dat", EM_field.E_z);

  write_vector_to_binary(filename_B + "_x.dat", EM_field.B_x);
  write_vector_to_binary(filename_B + "_y.dat", EM_field.B_y);
  write_vector_to_binary(filename_B + "_z.dat", EM_field.B_z);

  return 0;
}


int Propagate_particles(std::vector<std::vector<double>> &particles,
                        const RFD_matrix &RFD, const double dt) {

  for (int ip = 0; ip < particles.size(); ip++) {
    particles[ip][0] += RFD.RFD_x[ip] * dt;
    particles[ip][1] += RFD.RFD_y[ip] * dt;
    particles[ip][2] += RFD.RFD_z[ip] * dt;
    //printf("at particle %d, RFD: %lf, %lf, %lf \n", ip, RFD.RFD_x[ip],
    //      RFD.RFD_y[ip], RFD.RFD_z[ip]);
  }
  return 0;
}

class EM_wave_config {
public:
  int num_waves;
  std::vector<std::vector<double>> wave_configs;
  EM_wave_config(std::vector<std::vector<double>> init_vect) {
    if (init_vect[0].size() == 4) {
      num_waves = init_vect.size();
      wave_configs = init_vect;
    } else {
      printf("Invalid init vector!");
    }
  }
};

double Get_EM_wave_component(int dim, EM_wave_config config, double x, double y,
                             double t) {
  double component = 0.0;
  for (int ix = 0; ix < config.num_waves; ix++) {
    double magnitude = config.wave_configs[ix][0];
    double ang_freq = config.wave_configs[ix][1];
    double phi = config.wave_configs[ix][2];
    double alpha = config.wave_configs[ix][3];
    // Ex
    if (dim == 0) {
      component += 0.0;
      // Ey
    } else if (dim == 1) {
      component += 0.0;
      // Ez
    } else if (dim == 2) {
      component +=
          magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
      // Bx
    } else if (dim == 3) {
      component +=
          std::sin(phi) * magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
      // By
    } else if (dim == 4) {
      component +=
          -1.0 * std::cos(phi) * magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
      // Bz
    } else if (dim == 5) {
      component += 0.0;
    } else {
      printf("Invalid dimension!!");
    }
  }
  return component;
}

int Get_EM_field_from_positions(const int n_particles,
                                std::vector<std::vector<double>> particles,
                                EM_field_matrix &EM_field, double t,
                                EM_wave_config config) {

  for (int ix = 0; ix < n_particles; ix++) {
    EM_field.E_x[ix] =
        Get_EM_wave_component(0, config, particles[ix][0], particles[ix][1], t);
    EM_field.E_y[ix] =
        Get_EM_wave_component(1, config, particles[ix][0], particles[ix][1], t);
    EM_field.E_z[ix] =
        Get_EM_wave_component(2, config, particles[ix][0], particles[ix][1], t);

    EM_field.B_x[ix] =
        Get_EM_wave_component(3, config, particles[ix][0], particles[ix][1], t);
    EM_field.B_y[ix] =
        Get_EM_wave_component(4, config, particles[ix][0], particles[ix][1], t);
    EM_field.B_z[ix] =
        Get_EM_wave_component(5, config, particles[ix][0], particles[ix][1], t);
  }
  return 0;
}

int run_debug_particle_propagation() {

  int nx = 64;
  int ny = 64;
  int n_particles = 100;
  int pp = 1;
  int sign_plus = 1;

  int num_waves = 4;
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
  for (int ix = 0; ix < num_waves; ix++) {
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

  for (int ix = 0; ix < n_particles; ix++) {
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
      for (int ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][0];
      }
      write_vector_to_binary(xpos_filename + std::to_string(t), vect);
      for (int ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][1];
      }
      write_vector_to_binary(ypos_filename + std::to_string(t), vect);
      for (int ix = 0; ix < n_particles; ix++) {
        vect[ix] = particles[ix][2];
      }
      write_vector_to_binary(zpos_filename + std::to_string(t), vect);
      Write_EM_to_binary(E_filename + std::to_string(t),
                         B_filename + std::to_string(t), EM_field );
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

  Write_EM_to_binary(E_filename, B_filename, EM_field );

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
