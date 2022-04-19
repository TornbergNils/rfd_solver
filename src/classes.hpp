#ifndef CLASSES_H
#define CLASSES_H

#include <algorithm>
#include <functional>
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

  void Save(std::string filename, bool append) {
    std::ofstream savestream_Ex;
    std::ofstream savestream_Ey;
    std::ofstream savestream_Ez;
    std::ofstream savestream_Bx;
    std::ofstream savestream_By;
    std::ofstream savestream_Bz;
    if (append == true) {
      savestream_Ex.open(filename + "E_x",
                         std::ios::out | std::ios::app | std::ios::binary);
      savestream_Ey.open(filename + "E_y",
                         std::ios::out | std::ios::app | std::ios::binary);
      savestream_Ez.open(filename + "E_z",
                         std::ios::out | std::ios::app | std::ios::binary);
      savestream_Bx.open(filename + "B_x",
                         std::ios::out | std::ios::app | std::ios::binary);
      savestream_By.open(filename + "B_y",
                         std::ios::out | std::ios::app | std::ios::binary);
      savestream_Bz.open(filename + "B_z",
                         std::ios::out | std::ios::app | std::ios::binary);
    } else {
      savestream_Ex.open(filename + "E_x",
                         std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_Ey.open(filename + "E_y",
                         std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_Ez.open(filename + "E_z",
                         std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_Bx.open(filename + "B_x",
                         std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_By.open(filename + "B_y",
                         std::ios::out | std::ios::trunc | std::ios::binary);
      savestream_Bz.open(filename + "B_z",
                         std::ios::out | std::ios::trunc | std::ios::binary);
    }

    savestream_Ex.write((char *)&E_x[0], E_x.size() * sizeof(double));
    savestream_Ey.write((char *)&E_y[0], E_y.size() * sizeof(double));
    savestream_Ez.write((char *)&E_z[0], E_z.size() * sizeof(double));
    savestream_Bx.write((char *)&B_x[0], B_x.size() * sizeof(double));
    savestream_By.write((char *)&B_y[0], B_y.size() * sizeof(double));
    savestream_Bz.write((char *)&B_z[0], B_z.size() * sizeof(double));
    savestream_Ex.close();
    savestream_Ey.close();
    savestream_Ez.close();
    savestream_Bx.close();
    savestream_By.close();
    savestream_Bz.close();
  }

  int Elementwise_add(EM_field_matrix &EM_field) {
    if (nx != EM_field.nx || ny != EM_field.ny) {
      printf("Cannot add fields of different grid sizes!");
      return -1;
    }
    for (int ix = 0; ix < nx * ny; ix++) {
      E_x[ix] += EM_field.E_x[ix];
      E_y[ix] += EM_field.E_y[ix];
      E_z[ix] += EM_field.E_z[ix];

      B_x[ix] += EM_field.B_x[ix];
      B_y[ix] += EM_field.B_y[ix];
      B_z[ix] += EM_field.B_z[ix];
    }
    return 0;
  }
};


#endif // CLASSES_H