#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

class EM_field_matrix {
public:
  int nx;
  int ny;
  std::vector<double> E_x;
  std::vector<double> E_y;
  std::vector<double> E_z;

  std::vector<double> B_x;
  std::vector<double> B_y;
  std::vector<double> B_z;

  // Construct with 6 nx * ny zero_init vectors
  EM_field_matrix(int init_nx, int init_ny)
      : nx(init_nx), ny(init_ny), E_x(nx * ny, 0.0), E_y(nx * ny, 0.0),
        E_z(nx * ny, 0.0), B_x(nx * ny, 0.0), B_y(nx * ny, 0.0),
        B_z(nx * ny, 0.0) {}
};

// Dimension indicates what part of cross prod, 0,1,2 => x,y,z
std::vector<double> Get_cross_product(EM_field_matrix EM_field, int dimension) {
  const int nx = EM_field.nx;
  const int ny = EM_field.ny;

  std::vector<double> cross_product(nx * ny);

  // x-component of cross prod
  if (dimension == 0) {

    for (int ix = 0; ix < nx * ny; ix++) {
      cross_product[ix] =
          (EM_field.E_y[ix] * EM_field.B_z[ix] - EM_field.E_z[ix] * EM_field.B_y[ix]);
    }
    // y-component
  } else if (dimension == 1) {

    for (int ix = 0; ix < nx * ny; ix++) {
      cross_product[ix] =
          (EM_field.E_x[ix] * EM_field.B_z[ix] - EM_field.E_z[ix] * EM_field.B_x[ix]);
    }
    // z-component
  } else if (dimension == 2) {

    for (int ix = 0; ix < nx * ny; ix++) {
      cross_product[ix] =
          (EM_field.E_x[ix] * EM_field.B_y[ix] - EM_field.E_y[ix] * EM_field.B_x[ix]);
    }
  } else {
    printf("Invalid dimension for cross product!");
  }

  return cross_product;
}

std::vector<double> get_w(const EM_field_matrix EM_field) {
  const int nx = EM_field.nx;
  const int ny = EM_field.ny;

  std::vector<double> w(nx * ny);

  std::vector<double> E_cross_B_x = Get_cross_product(EM_field, 0);
  std::vector<double> E_cross_B_y = Get_cross_product(EM_field, 1);
  std::vector<double> E_cross_B_z = Get_cross_product(EM_field, 2);

  for (int ix = 0; ix < nx * ny; ix++) {
    double E_cross_B_squared = E_cross_B_x[ix] * E_cross_B_x[ix] +
                               E_cross_B_y[ix] * E_cross_B_y[ix] +
                               E_cross_B_z[ix] * E_cross_B_z[ix];

    double E_squared = EM_field.E_x[ix] * EM_field.E_x[ix] +
                       EM_field.E_y[ix] * EM_field.E_y[ix] +
                       EM_field.E_z[ix] * EM_field.E_z[ix];

    double B_squared = EM_field.B_x[ix] * EM_field.B_x[ix] +
                       EM_field.B_y[ix] * EM_field.B_y[ix] +
                       EM_field.B_z[ix] * EM_field.B_z[ix];

    double E_sq_plus_B_sq = E_squared + B_squared;
    w[ix] = 4.0 * E_cross_B_squared / (E_sq_plus_B_sq * E_sq_plus_B_sq);
  }
  return w;
}

std::vector<double> get_u(const EM_field_matrix EM_field,
                          const std::vector<double> w) {

  const int nx = EM_field.nx;
  const int ny = EM_field.ny;

  std::vector<double> u( nx * ny );
  //std::vector<double> w = get_w( EM_field );
  
  for (int ix = 0; ix < nx * ny; ix++) {

    double E_squared = EM_field.E_x[ix] * EM_field.E_x[ix] +
                       EM_field.E_y[ix] * EM_field.E_y[ix] +
                       EM_field.E_z[ix] * EM_field.E_z[ix];

    double B_squared = EM_field.B_x[ix] * EM_field.B_x[ix] +
                       EM_field.B_y[ix] * EM_field.B_y[ix] +
                       EM_field.B_z[ix] * EM_field.B_z[ix];

    double factor1 = 2.0 * B_squared / (E_squared + B_squared );
    double factor2 = ( 1 - sqrt( 1 - w[ix] ) ) / w[ix];

    u[ix] =  factor1 * factor2;
  }
  
  return u;
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
  }

};

    int
    write_vector_to_binary(std::string filename, std::vector<double> vect) {

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

int main() {

  int nx = 8;
  int ny = 4;
  EM_field_matrix EM_field(nx, ny);

  std::string E_filename("./data/E_data");
  std::string B_filename("./data/B_data");

  for (int ix = 0; ix < nx * ny; ix++) {
    EM_field.E_x[ix] = ix;
    EM_field.B_x[ix] = -ix;
  }

  for (int ix = 0; ix < nx * ny; ix++) {
    std::cout << EM_field.E_x[ix] << ", ";
    std::cout << EM_field.B_x[ix] << "\n";
  }

  Write_EM_to_binary(E_filename, B_filename, EM_field, nx, ny);

  return 0;
}
