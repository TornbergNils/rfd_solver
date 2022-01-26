#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class EM_field_matrix {
public:
  int nx;
  int ny;
  std::vector<double> E_x;
  std::vector<double> B_x;
  std::vector<double> E_y;
  std::vector<double> B_y;
  std::vector<double> E_z;
  std::vector<double> B_z;

  // Construct with 6 nx * ny zero_init vectors
  EM_field_matrix(int init_nx, int init_ny)
      : nx(init_nx), ny(init_ny), E_x(nx * ny, 0.0), B_x(nx * ny, 0.0),
        E_y(nx * ny, 0.0), B_y(nx * ny, 0.0), E_z(nx * ny, 0.0),
        B_z(nx * ny, 0.0) {}
};

int write_vector_to_binary( std::string filename, std::vector<double> vect ) {
  
  std::ofstream filestream(filename, std::ios::out | std::ios::binary);

  filestream.write((char *)&vect[0],
                     vect.size() * sizeof(double));
  filestream.close();
  
  return 0;
}

// Takes EM matrix object and writes content to file,
// nx columns, ny rows
// adapted from old code, TODO: Refresh memory on how it works
int Write_EM_to_binary( std::string filename_E, std::string filename_B,
                       EM_field_matrix EM_field, int nx, int ny) {

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();
  
  // write all components of EM_field to file
  write_vector_to_binary( filename_E + "_x.dat", EM_field.E_x );
  write_vector_to_binary( filename_E + "_y.dat", EM_field.E_y );
  write_vector_to_binary( filename_E + "_z.dat", EM_field.E_z );
  
  write_vector_to_binary( filename_B + "_x.dat", EM_field.B_x );
  write_vector_to_binary( filename_B + "_y.dat", EM_field.B_y );
  write_vector_to_binary( filename_B + "_z.dat", EM_field.B_z );

  std::ofstream E_filestream(filename_E, std::ios::out | std::ios::binary);
  std::ofstream B_filestream(filename_B, std::ios::out | std::ios::binary);
 
  return 0;
}

int main() {

  int nx = 4;
  int ny = 4;
  EM_field_matrix EM_field(nx, ny);
  
  std::string E_filename( "./data/E_data" );
  std::string B_filename( "./data/B_data" );

  for (int ix = 0; ix < nx * ny; ix++) {
    EM_field.E_x[ix] = ix;
    EM_field.B_x[ix] = -ix;
  }

  for (int ix = 0; ix < nx * ny; ix++) {
    std::cout << EM_field.E_x[ix] << ", ";
    std::cout << EM_field.B_x[ix] << "\n";
  }

  Write_EM_to_binary(E_filename, B_filename, EM_field, nx,
                     ny);

  return 0;
}
