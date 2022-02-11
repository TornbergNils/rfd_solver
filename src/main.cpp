#include "classes.hpp"
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

  // run_debug_RFD_function();
  run_debug_particle_propagation();

  return 0;
}

// Dimension indicates what part of cross prod, 0,1,2 => x,y,z
double Get_cross_product(const int ix, EM_field_matrix EM_field,
                         const int dimension) {
  double cross_product;

  // x-component of cross prod
  if (dimension == 0) {

    cross_product = (EM_field.E_y[ix] * EM_field.B_z[ix] -
                     EM_field.E_z[ix] * EM_field.B_y[ix]);
    // y-component
  } else if (dimension == 1) {
    cross_product = (EM_field.E_x[ix] * EM_field.B_z[ix] -
                     EM_field.E_z[ix] * EM_field.B_x[ix]);
    // z-component
  } else if (dimension == 2) {
    cross_product = (EM_field.E_x[ix] * EM_field.B_y[ix] -
                     EM_field.E_y[ix] * EM_field.B_x[ix]);
  } else {
    printf("Invalid dimension for cross product!");
  }

  return cross_product;
}

double get_w(const double E_cross_B_squared, const double E_squared,
             const double B_squared) {

  double E_sq_plus_B_sq = E_squared + B_squared;
  double w = 4.0 * E_cross_B_squared / (E_sq_plus_B_sq * E_sq_plus_B_sq);
  return w;
}

double get_u(const double w, const double E_squared, const double B_squared,
             const double eps) {

  double factor1 = 2.0 * B_squared / (E_squared + B_squared);
  double factor2;
  if (w > 1.0 - eps) {
    factor2 = 1.0 / w;
  } else if (w > eps) {
    factor2 = (1.0 - std::sqrt(1 - w)) / w;
  } else {
    factor2 = 1.0 / 2.0 + w / 8;
  }
  double u = factor1 * factor2;
  return u;
}

double get_RFD_component(const double u, const double w,
                         const double E_cross_B_component,
                         const double B_component, const double E_component,
                         const double E_squared, const double B_squared,
                         const double E_dot_B, const double E_cross_B_squared,
                         const double sign) {

  const double term1 = std::sqrt(u - u * u) * E_cross_B_component;

  const double term2 = (1 - u) * std::sqrt(B_squared) * E_component +
                       u * E_dot_B * B_component / std::sqrt(B_squared);

  const double factor1 = term1 + sign * term2;

  const double factor2 =
      1.0 / std::sqrt(E_squared * B_squared - u * E_cross_B_squared);

  const double RFD_component = factor1 * factor2;
  return RFD_component;
}

class RFD_matrix {
public:
  int nx;
  int ny;
  std::vector<double> RFD_x;
  std::vector<double> RFD_y;
  std::vector<double> RFD_z;
  int sign;

  RFD_matrix(EM_field_matrix EM_field, int sign_init) {
    nx = EM_field.nx;
    ny = EM_field.ny;
    sign = sign_init;
    std::vector<double> temp_RFD_x(nx * ny);
    std::vector<double> temp_RFD_y(nx * ny);
    std::vector<double> temp_RFD_z(nx * ny);

    const double eps = 10e-15;

    for (int ix = 0; ix < nx * ny; ix++) {
      const double E_cross_B_x = Get_cross_product(ix, EM_field, 0);
      const double E_cross_B_y = Get_cross_product(ix, EM_field, 1);
      const double E_cross_B_z = Get_cross_product(ix, EM_field, 2);

      const double E_squared = EM_field.E_x[ix] * EM_field.E_x[ix] +
                               EM_field.E_y[ix] * EM_field.E_y[ix] +
                               EM_field.E_z[ix] * EM_field.E_z[ix];

      const double B_squared = EM_field.B_x[ix] * EM_field.B_x[ix] +
                               EM_field.B_y[ix] * EM_field.B_y[ix] +
                               EM_field.B_z[ix] * EM_field.B_z[ix];

      const double E_cross_B_squared = E_cross_B_x * E_cross_B_x +
                                       E_cross_B_y * E_cross_B_y +
                                       E_cross_B_z * E_cross_B_z;

      const double E_dot_B = EM_field.E_x[ix] * EM_field.B_x[ix] +
                             EM_field.E_y[ix] * EM_field.B_y[ix] +
                             EM_field.E_z[ix] * EM_field.B_z[ix];

      // This case, E=0 needs more careful investigation
      if (E_squared < eps) {
        // use copysign as signum

        temp_RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_x[ix];
        temp_RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_y[ix];
        temp_RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.B_z[ix];

        // This case, B=0 needs more careful investigation
      } else if (B_squared < eps) {

        temp_RFD_x[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_x[ix];
        temp_RFD_y[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_y[ix];
        temp_RFD_z[ix] = sign * std::copysign(1.0, E_dot_B) * EM_field.E_z[ix];

        // Case 2b, E_dot_B = 0
      } else if ((E_dot_B * E_dot_B) < eps && E_squared <= B_squared + eps) {
        // printf( "%lf\n", E_dot_B );
        // Case 2bi, E_dot_B = 0 and |E| = |B|
        if (std::abs(std::sqrt(E_squared) - std::sqrt(B_squared)) < eps) {

          double denominator = std::sqrt(E_cross_B_squared);

          temp_RFD_x[ix] = E_cross_B_x / denominator;
          temp_RFD_y[ix] = E_cross_B_y / denominator;
          temp_RFD_z[ix] = E_cross_B_z / denominator;
          // Case 2bii, E_dot_B = 0 and E^2 < B^2
        } else {

          const double w = get_w(E_cross_B_squared, E_squared, B_squared);
          const double u = get_u(w, E_squared, B_squared, eps);

          double term1x = E_cross_B_x;
          double term1y = E_cross_B_y;
          double term1z = E_cross_B_z;

          double term2x;
          double term2y;
          double term2z;
          if ( u > 1.0 - eps) {
            term2x = 0.0;
            term2y = 0.0;
            term2z = 0.0;
          } else {
            term2x = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_x[ix];
            term2y = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_y[ix];
            term2z = sign * std::sqrt(1.0 - u) * std::sqrt(B_squared) *
                     EM_field.E_z[ix];
          }
          double term3x = sign * EM_field.B_x[ix];
          double term3y = sign * EM_field.B_y[ix];
          double term3z = sign * EM_field.B_z[ix];

          double numeratorx = term1x + term2x + term3x;
          double numeratory = term1y + term2y + term3y;
          double numeratorz = term1z + term2z + term3z;
          double denominator =
              std::sqrt(numeratorx * numeratorx + numeratory * numeratory +
                        numeratorz * numeratorz);

          temp_RFD_x[ix] = numeratorx / denominator;
          temp_RFD_y[ix] = numeratory / denominator;
          temp_RFD_z[ix] = numeratorz / denominator;
        }
      } else {
        const double w = get_w(E_cross_B_squared, E_squared, B_squared);

        const double u = get_u(w, E_squared, B_squared, eps);

        temp_RFD_x[ix] = get_RFD_component(
            u, w, E_cross_B_x, EM_field.B_x[ix], EM_field.E_x[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);

        temp_RFD_y[ix] = get_RFD_component(
            u, w, E_cross_B_y, EM_field.B_y[ix], EM_field.E_y[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);

        temp_RFD_z[ix] = get_RFD_component(
            u, w, E_cross_B_z, EM_field.B_z[ix], EM_field.E_z[ix], E_squared,
            B_squared, E_dot_B, E_cross_B_squared, sign);
      }
    }
    RFD_x = temp_RFD_x;
    RFD_y = temp_RFD_y;
    RFD_z = temp_RFD_z;
  }
};

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
                       EM_field_matrix EM_field, int nx, int ny) {

  // write all components of EM_field to file
  write_vector_to_binary(filename_E + "_x.dat", EM_field.E_x);
  write_vector_to_binary(filename_E + "_y.dat", EM_field.E_y);
  write_vector_to_binary(filename_E + "_z.dat", EM_field.E_z);

  write_vector_to_binary(filename_B + "_x.dat", EM_field.B_x);
  write_vector_to_binary(filename_B + "_y.dat", EM_field.B_y);
  write_vector_to_binary(filename_B + "_z.dat", EM_field.B_z);

  return 0;
}

const double amplitude1 = 0.33;
const double amplitude2 = 0.4;
const double amplitude3 = 0.84;
const double amplitude4 = 0.2;
const double alpha1 = 2.11;
const double alpha2 = 0.54;
const double alpha3 = 1.33;
const double alpha4 = 0.25;
const double omega1 = 0.7;
const double omega2 = 0.45;
const double omega3 = 1.2;
const double omega4 = 0.4;

double Ex_function(const std::vector<double> &particle, double t) {

  double wave3 = amplitude3 * std::cos(omega2 * (-particle[1] - t) + alpha3);

  return wave3;
}
double Ey_function(const std::vector<double> &particle, double t) {
  double wave1 = amplitude1 * std::cos(omega1 * (particle[0] - t) + alpha1);

  return wave1;
}
double Ez_function(const std::vector<double> &particle, double t) {

  double wave2 = amplitude2 * std::cos(omega2 * (particle[0] - t) + alpha2);
  double wave4 = amplitude4 * std::cos(omega1 * (-particle[1] - t) + alpha4);

  return wave2 + wave4;
}
double Bx_function(const std::vector<double> &particle, double t) {

  double wave4 = amplitude4 * std::cos(omega1 * (-particle[1] - t) + alpha4);
  return wave4;
}
double By_function(const std::vector<double> &particle, double t) {

  double wave2 = -amplitude2 * std::cos(omega2 * (particle[0] - t) + alpha2);

  return wave2;
}
double Bz_function(const std::vector<double> &particle, double t) {

  double wave1 = amplitude1 * std::cos(omega1 * (particle[0] - t) + alpha1);
  double wave3 = -amplitude3 * std::cos(omega2 * (-particle[1] - t) + alpha3);

  return wave1 + wave3;
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
                         B_filename + std::to_string(t), EM_field, nx, ny);
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

  Write_EM_to_binary(E_filename, B_filename, EM_field, nx, ny);

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
