#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


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

double get_u(const double w, const double E_squared, const double B_squared, const double eps) {

  double factor1 = 2.0 * B_squared / (E_squared + B_squared);
  double factor2;
  if( w > eps ) {
    factor2 = (1.0 - std::sqrt(1 - w)) / w;
  }
  else {
    factor2 = 1.0/2.0;
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
    const double eps = 10e-14;
    
    printf("%.16lf \n", eps );

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
      if( E_squared < eps ) {
        temp_RFD_x[ix] = 1.0;
        temp_RFD_y[ix] = 0.0;
        temp_RFD_z[ix] = 0.0;
        printf( "%lf\n", temp_RFD_x[ix] );
        printf( "%lf\n", temp_RFD_y[ix] );
        printf( "%lf\n", temp_RFD_z[ix] );

      }
      else if( ( E_dot_B * E_dot_B ) < eps && E_squared <= B_squared ) {
        // printf( "%lf\n", E_dot_B );
        // Case 2bi, E_dot_B = 0 and |E| = |B|
        if( std::abs( std::sqrt( E_squared ) - std::sqrt(B_squared) ) < eps ) {
          printf( "at 0: %lf,",
              std::abs( std::sqrt( E_squared ) - std::sqrt(B_squared)));

          temp_RFD_x[ix] = E_cross_B_x;
          temp_RFD_y[ix] = E_cross_B_y;
          temp_RFD_z[ix] = E_cross_B_z;
        }  
        else {

        const double w = get_w(E_cross_B_squared, E_squared, B_squared);

        const double u = get_u(w, E_squared, B_squared, eps);
        
        temp_RFD_x[ix] = ( std::sqrt( u ) * E_cross_B_x 
          + std::sqrt( 1.0 - u ) * std::sqrt( B_squared ) * EM_field.E_x[ix]
          + u * E_dot_B * EM_field.B_x[ix] ) 
          / ( std::sqrt( E_squared ) * std::sqrt( B_squared ) );

        temp_RFD_y[ix] = ( std::sqrt( u ) * E_cross_B_y 
          + std::sqrt( 1.0 - u ) * std::sqrt( B_squared ) * EM_field.E_y[ix]
          + u * E_dot_B * EM_field.B_y[ix] ) 
          / ( std::sqrt( E_squared ) * std::sqrt( B_squared ) );
        
        temp_RFD_z[ix] = ( std::sqrt( u ) * E_cross_B_z 
          + std::sqrt( 1.0 - u ) * std::sqrt( B_squared ) * EM_field.E_z[ix]
          + u * E_dot_B * EM_field.B_z[ix] ) 
          / ( std::sqrt( E_squared ) * std::sqrt( B_squared ) );
        }
      }
      else {
        const double w = get_w(E_cross_B_squared, E_squared, B_squared);

        const double u = get_u(w, E_squared, B_squared, eps);
        
        //printf( "w = %lf, u = %lf, \n", w, u );
      
        temp_RFD_x[ix] = get_RFD_component(u, w, E_cross_B_x,
            EM_field.B_x[ix], EM_field.E_x[ix], E_squared, B_squared,
            E_dot_B, E_cross_B_squared, sign);
      
        temp_RFD_y[ix] = get_RFD_component(u, w, E_cross_B_y,
            EM_field.B_y[ix], EM_field.E_y[ix], E_squared, B_squared,
            E_dot_B, E_cross_B_squared, sign);

        temp_RFD_z[ix] = get_RFD_component(u, w, E_cross_B_z,
            EM_field.B_z[ix], EM_field.E_z[ix], E_squared, B_squared,
            E_dot_B, E_cross_B_squared, sign);
      }
      // DEBUG 
      /*
      bool checknan = std::isnan( temp_RFD_x[ix] ) ||
          std::isnan( temp_RFD_y[ix] ) || std::isnan( temp_RFD_z[ix] );
      if( checknan ) {
        printf( "NaN detected at %i", ix );
      }
      */

      

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

int main() {

  int nx = 25;
  int ny = 25;
  double EB_max = 3.0;
  EM_field_matrix EM_field(nx, ny);

  std::string E_filename("./data/E_data");
  std::string B_filename("./data/B_data");

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      EM_field.E_x[ix*ny + iy ] = 2.0 * ( iy/(double) (ny - 1 ) ) * EB_max - EB_max;
      EM_field.E_y[ix*ny + iy ] = 2.0 * ( ix/(double) (nx - 1 ) ) * EB_max - EB_max;
      EM_field.B_x[ix*ny + iy ] = 1.0;
    }
  }

  Write_EM_to_binary(E_filename, B_filename, EM_field, nx, ny);

  RFD_matrix RFD(EM_field);

  std::vector<double> magnitudes( nx * ny );

  for( int ix = 0; ix < nx * ny; ix++ ) {
    magnitudes[ix] = std::sqrt( RFD.RFD_x[ix] * RFD.RFD_x[ix] 
                         + RFD.RFD_y[ix] * RFD.RFD_y[ix] 
                         + RFD.RFD_z[ix] * RFD.RFD_z[ix] );
  }

  // NORMALIZE
  for( int ix = 0; ix < nx * ny; ix++ ) {
    RFD.RFD_x[ix] /= magnitudes[ix];
    RFD.RFD_y[ix] /= magnitudes[ix];
    RFD.RFD_z[ix] /= magnitudes[ix];
      
  }
  
  for( int ix = 0; ix < nx * ny; ix++ ) {
    magnitudes[ix] = std::sqrt( RFD.RFD_x[ix] * RFD.RFD_x[ix] 
                         + RFD.RFD_y[ix] * RFD.RFD_y[ix] 
                         + RFD.RFD_z[ix] * RFD.RFD_z[ix] );
    bool checknan = std::isnan( RFD.RFD_x[ix] ) ||
      std::isnan( RFD.RFD_y[ix] ) || std::isnan( RFD.RFD_z[ix] );
    if( checknan ) {
      printf( "NaN detected at %i", ix );
    }
  }
  
  std::string filename_RFD = "./data/RFD";
  write_vector_to_binary(filename_RFD + "_x.dat", RFD.RFD_x);
  write_vector_to_binary(filename_RFD + "_y.dat", RFD.RFD_y);
  write_vector_to_binary(filename_RFD + "_z.dat", RFD.RFD_z);
  write_vector_to_binary("./data/magnitudes.dat", magnitudes);

  return 0;
}
