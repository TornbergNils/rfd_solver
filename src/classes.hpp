#include <vector>
#include <algorithm>
#include <functional>

class Particle {
public:
  int charge;
  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
};

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
  
  int Elementwise_add( EM_field_matrix & EM_field  ) {
    if( nx != EM_field.nx || ny != EM_field.ny ) {
      printf( "Cannot add fields of different grid sizes!" );
      return -1;
    }
    for( int ix = 0; ix < nx * ny; ix++ ) {
      E_x[ix] += EM_field.E_x[ix];
      E_y[ix] += EM_field.E_x[ix];
      E_z[ix] += EM_field.E_x[ix];
      
      B_x[ix] += EM_field.E_x[ix];
      B_y[ix] += EM_field.E_x[ix];
      B_z[ix] += EM_field.E_x[ix];
    }
    
  }
  
};
