
static inline int Get_index(int ix, int iy, int nx, int ny) {
  if (ix < 0) {
    ix = nx - 1;
  }
  if (iy < 0) {
    iy = ny - 1;
  }
  return (ix % nx) + (iy % ny) * nx;
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
  return 0;
}
