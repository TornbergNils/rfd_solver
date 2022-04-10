
int Propagate_particles(std::vector<std::vector<double>> &particles,
                        const RFD_matrix &RFD, const double dt) {
  int ip_max = particles.size();
  for (int ip = 0; ip < ip_max; ip++) {
    particles[ip][0] += RFD.RFD_x[ip] * dt;
    particles[ip][1] += RFD.RFD_y[ip] * dt;
    particles[ip][2] += RFD.RFD_z[ip] * dt;
    // printf("at particle %d, RFD: %lf, %lf, %lf \n", ip, RFD.RFD_x[ip],
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
      component +=
          magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
      // Ez
    } else if (dim == 2) {
      component += 0.0;
      // Bx
    } else if (dim == 3) {
      component +=
          std::sin(phi) * magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
      // By
    } else if (dim == 4) {
      component += 0.0;
      // Bz
    } else if (dim == 5) {
      component +=
          1.0 * std::cos(phi) * magnitude *
          std::cos(ang_freq * (x * std::cos(phi) + y * std::sin(phi) - t) +
                   alpha);
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
