#ifndef FDTD_H
#define FDTD_H


/*
  Test file for the FDTD algorithm before inclusion in solver.hpp
*/

inline int Get_index(int ix, int iy, int nx, int ny) {
  int idx = ((ix + nx) % nx) + ((iy + ny) % ny) * nx;
  //printf( "%d \n", idx );
  return idx;
}

int Update_mode(EM_field_matrix &all_modes, const double delta_t,
                const double delta_x) {

  const int nx = all_modes.nx;
  const int ny = all_modes.ny;

  // Start with H parts of modes, "1/2 timesteps"
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      int index = Get_index(ix, iy, nx, ny);

      all_modes.B_z[index] = all_modes.B_z[index] + delta_t / delta_x *
                              ((all_modes.E_x[Get_index(ix, iy + 1, nx, ny)] -
                                all_modes.E_x[Get_index(ix, iy, nx, ny)]) -
                               (all_modes.E_y[Get_index(ix + 1, iy, nx, ny)] -
                                all_modes.E_y[Get_index(ix, iy, nx, ny)]));

      all_modes.B_x[index] = all_modes.B_x[index] - delta_t / delta_x *
                              (all_modes.E_z[Get_index(ix, iy + 1, nx, ny)] -
                               all_modes.E_z[Get_index(ix, iy, nx, ny)]);

      all_modes.B_y[index] = all_modes.B_y[index] + delta_t / delta_x *
                              (all_modes.E_z[Get_index(ix + 1, iy, nx, ny)] -
                               all_modes.E_z[Get_index(ix, iy, nx, ny)]);
    }
  }
  // Update E, remaining "1/2 timestep"
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      int index = Get_index(ix, iy, nx, ny);

      all_modes.E_z[index] =
          all_modes.E_z[index] +
          delta_t / delta_x *
              ((all_modes.B_y[Get_index(ix, iy, nx, ny)] -
                all_modes.B_y[Get_index(ix - 1, iy, nx, ny)]) -
               (all_modes.B_x[Get_index(ix, iy, nx, ny)] -
                all_modes.B_x[Get_index(ix, iy - 1, nx, ny)]));

      all_modes.E_x[index] = all_modes.E_x[index] +
                             delta_t / delta_x *
                                 (all_modes.B_z[Get_index(ix, iy, nx, ny)] -
                                  all_modes.B_z[Get_index(ix, iy - 1, nx, ny)]);

      all_modes.E_y[index] = all_modes.E_y[index] -
                             delta_t / delta_x *
                                 (all_modes.B_z[Get_index(ix, iy, nx, ny)] -
                                  all_modes.B_z[Get_index(ix - 1, iy, nx, ny)]);

      // printf("%d, %d, %lf, %lf \n", ix, iy,
      //    all_modes.E_z[index], all_modes.E_y[index] );
    }
  }
  /*
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < 20; ix++) {
      
      all_modes.E_x[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_x[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
      
      all_modes.E_y[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_y[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
      
      all_modes.E_z[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_z[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
      
      all_modes.B_x[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_x[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
      
      all_modes.B_y[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_y[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
      
      all_modes.B_z[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_z[Get_index(nx-ix, iy, nx, ny)] *= 0.5;
    }
  }
  for (int iy = 0; iy < 20; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      
      all_modes.E_x[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_x[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
      
      all_modes.E_y[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_y[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
      
      all_modes.E_z[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.E_z[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
      
      all_modes.B_x[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_x[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
      
      all_modes.B_y[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_y[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
      
      all_modes.B_z[Get_index(ix, iy, nx, ny)] *= 0.5;
      all_modes.B_z[Get_index(ix, ny-iy, nx, ny)] *= 0.5;
    }
  }
  */
  return 0;
}

#endif // FDTD_H