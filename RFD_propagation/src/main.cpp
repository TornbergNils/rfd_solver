#include "classes.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <random>


// Dimension indicates what part of cross prod, 0,1,2 => x,y,z
double Get_cross_product(const int ix, EM_field_matrix EM_field,
                         const int dimension) {
  const int nx = EM_field.nx;
  const int ny = EM_field.ny;

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
  if (w > eps) {
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
      1.0 / std::sqrt(E_squared * B_squared + u * E_cross_B_squared);

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

  RFD_matrix(EM_field_matrix EM_field) {
    nx = EM_field.nx;
    ny = EM_field.ny;
    std::vector<double> temp_RFD_x(nx * ny);
    std::vector<double> temp_RFD_y(nx * ny);
    std::vector<double> temp_RFD_z(nx * ny);

    // Note that this fixes the sign of RFD to +!!
    const int sign = 1;
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
        temp_RFD_x[ix] = std::copysign(1.0, E_dot_B) * EM_field.B_x[ix];
        temp_RFD_y[ix] = std::copysign(1.0, E_dot_B) * EM_field.B_y[ix];
        temp_RFD_z[ix] = std::copysign(1.0, E_dot_B) * EM_field.B_z[ix];

        // This case, B=0 needs more careful investigation
      } else if (B_squared < eps) {
        temp_RFD_x[ix] = std::copysign(1.0, E_dot_B) * EM_field.E_x[ix];
        temp_RFD_y[ix] = std::copysign(1.0, E_dot_B) * EM_field.E_y[ix];
        temp_RFD_z[ix] = std::copysign(1.0, E_dot_B) * EM_field.E_z[ix];

        // Case 2b, E_dot_B = 0
      } else if ((E_dot_B * E_dot_B) < eps && E_squared <= B_squared) {
        // printf( "%lf\n", E_dot_B );
        // Case 2bi, E_dot_B = 0 and |E| = |B|
        if (std::abs(std::sqrt(E_squared) - std::sqrt(B_squared)) < eps) {

          temp_RFD_x[ix] = E_cross_B_x;
          temp_RFD_y[ix] = E_cross_B_y;
          temp_RFD_z[ix] = E_cross_B_z;
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
          if (std::abs(u - 1.0) < eps) {
            term2x = 0.0;
            term2y = 0.0;
            term2z = 0.0;
          } else {
            term2x =
                std::sqrt(1.0 - u) * std::sqrt(B_squared) * EM_field.E_x[ix];
            term2y =
                std::sqrt(1.0 - u) * std::sqrt(B_squared) * EM_field.E_y[ix];
            term2z =
                std::sqrt(1.0 - u) * std::sqrt(B_squared) * EM_field.E_z[ix];
          }
          double term3x = EM_field.B_x[ix];
          double term3y = EM_field.B_y[ix];
          double term3z = EM_field.B_z[ix];

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

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("./data/header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  // write all components of EM_field to file
  write_vector_to_binary(filename_E + "_x.dat", EM_field.E_x);
  write_vector_to_binary(filename_E + "_y.dat", EM_field.E_y);
  write_vector_to_binary(filename_E + "_z.dat", EM_field.E_z);

  write_vector_to_binary(filename_B + "_x.dat", EM_field.B_x);
  write_vector_to_binary(filename_B + "_y.dat", EM_field.B_y);
  write_vector_to_binary(filename_B + "_z.dat", EM_field.B_z);

  return 0;
}

const double amplitude1 = 0.63;
const double amplitude2 = 0.4;
const double amplitude3 = 0.54;
const double amplitude4 = 0.7;
const double alpha1 = 0.11;
const double alpha2 = 0.54;
const double alpha3 = 0.73;
const double alpha4 = 0.25;
const double k1 = 1.0;
const double k2 = 0.5;
const double k3 = 1.2;
const double k4 = 0.4;
const double omega1 = 1.0;
const double omega2 = 0.5;
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

  double wave2 = amplitude2 * std::cos(omega2 * (-particle[0] - t) + alpha2);
  double wave4 = amplitude4 * std::cos(omega1 * (-particle[1] - t) + alpha4);

  return wave2 + wave4;
}
double Bx_function(const std::vector<double> &particle, double t) {

  double wave4 = amplitude4 * std::cos(omega1 * (-particle[1] - t) + alpha4);
  return wave4;
}
double By_function(const std::vector<double> &particle, double t) {

  double wave2 = -amplitude2 * std::cos(omega2 * (-particle[0] - t) + alpha2);

  return wave2;
}
double Bz_function(const std::vector<double> &particle, double t) {

  double wave1 = amplitude1 * std::cos(omega1 * (particle[0] - t) + alpha1);
  double wave3 = -amplitude3 * std::cos(omega2 * (-particle[1] - t) + alpha3);

  return wave1 + wave3;
}

int Get_EM_field_from_positions(const int n_particles,
                                std::vector<std::vector<double>> particles,
                                EM_field_matrix &EM_field, double t) {

  for (int ix = 0; ix < n_particles; ix++) {
    EM_field.E_x[ix] = Ex_function(particles[ix], t);
    EM_field.E_y[ix] = Ey_function(particles[ix], t);
    EM_field.E_z[ix] = Ez_function(particles[ix], t);

    EM_field.B_x[ix] = Bx_function(particles[ix], t);
    EM_field.B_y[ix] = By_function(particles[ix], t);
    EM_field.B_z[ix] = Bz_function(particles[ix], t);
  }
  return 0;
}

int Propagate_particles(std::vector<std::vector<double>> &particles,
                        const RFD_matrix &RFD) {

  for (int ip = 0; ip < particles.size(); ip++) {
    particles[ip][0] += RFD.RFD_x[ip];
    particles[ip][1] += RFD.RFD_y[ip];
    particles[ip][2] += RFD.RFD_z[ip];
  }
  return 0;
}

int main() {

  int nx = 24;
  int ny = 24;
  int n_particles = 25;
  int pp = 1;
  // double c = 1.0;

  double x_max = 6.0;
  double y_max = 6.0;
  double z_max = 2.0;

  std::vector<std::vector<double>> particles{n_particles,
                                             std::vector<double>(3)};
  // Setup rng
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution;
  auto random = std::bind( distribution, generator );

  for( int ix = 0; ix < n_particles; ix++ ) {
    particles[ix][0] = 2.0 * random() * x_max - x_max;
    particles[ix][1] = 2.0 * random() * y_max - y_max;
    particles[ix][2] = 2.0 * random() * z_max - z_max;
  }

  EM_field_matrix field_at_particles(n_particles, pp);
  EM_field_matrix EM_field(nx, ny);

  // Main simulation loop
  for (double t = 0.0; t < 10; t = t + 0.2) {

    // For visualizing the fields
    for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {

        std::vector<double> coords = {
            2.0 * iy / (double)(ny - 1) * x_max - x_max,
            2.0 * ix / (double)(nx - 1) * y_max - y_max, 0};
        // printf( "%lf, %lf, %lf \n", coords[0] - coords[1],
        // coords[1], coords[2] );

        EM_field.E_x[ix * ny + iy] = Ex_function(coords, t);
        EM_field.E_y[ix * ny + iy] = Ey_function(coords, t);
        EM_field.E_z[ix * ny + iy] = Ez_function(coords, t);

        EM_field.B_x[ix * ny + iy] = Bx_function(coords, t);
        EM_field.B_y[ix * ny + iy] = By_function(coords, t);
        EM_field.B_z[ix * ny + iy] = Bz_function(coords, t);
      }
    }
    Get_EM_field_from_positions(n_particles, particles, field_at_particles, t);
    RFD_matrix RFD_at_particles(field_at_particles);

    Propagate_particles(particles, RFD_at_particles);
    std::string xpos_filename("./data/pos_x");
    std::string ypos_filename("./data/pos_y");
    std::string zpos_filename("./data/pos_z");

    // Temporary solution for storing particle positions
    std::vector<double> vect(n_particles);
    for( int ix = 0; ix < n_particles; ix++ ) {
      vect[ix] = particles[ix][0];
    }
    write_vector_to_binary(xpos_filename + std::to_string(t),
                            vect);
    for( int ix = 0; ix < n_particles; ix++ ) {
      vect[ix] = particles[ix][1];
    }
    write_vector_to_binary(ypos_filename + std::to_string(t),
                            vect);
    for( int ix = 0; ix < n_particles; ix++ ) {
      vect[ix] = particles[ix][2];
    }
    write_vector_to_binary(zpos_filename + std::to_string(t),
                            vect);

    std::string E_filename("./data/E_data");
    std::string B_filename("./data/B_data");

    Write_EM_to_binary(E_filename + std::to_string(t),
                       B_filename + std::to_string(t), EM_field, nx, ny);
  }

  return 0;
}
