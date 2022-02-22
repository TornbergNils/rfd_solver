
class Solver {
  // Given parameters
  int nx;
  int ny;
  int n_particles;
  int tmax;
  int n_tsteps;
  int save_rate;
  std::vector<double>::size_type n_elec;
  std::vector<double>::size_type n_posi;

  double delta_x;
  double delta_y;

  // Derived parameters
  double dt;

  // Data
  std::vector<double> electron_pos;
  std::vector<double> positron_pos;
  EM_field_matrix EM;
  RFD_matrix RFD;

public:
  Solver(int nx, int ny, int n_particles, int tmax, int n_tsteps, int save_rate)
      : nx(nx), ny(ny), n_particles(n_particles), tmax(tmax),
        n_tsteps(n_tsteps), save_rate(save_rate), n_elec(n_particles),
        n_posi(n_particles), electron_pos(n_particles * 3),
        positron_pos(n_particles * 3), EM(nx, ny), RFD(EM, 1) {}

  void Initialize(EM_field_matrix EM_IC) {
    // Setup rng
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution;
    auto random = std::bind(distribution, generator);

    dt = (double)tmax / (double)n_tsteps;
    double x_len = nx * delta_x;
    double y_len = ny * delta_y;

    // set EM to equal to EM_IC, since EM is 0-initialized by constructor
    EM.Elementwise_add(EM_IC);
    RFD.Update(EM, 1);

    // Set particle positions to randomly be in x-y plane
    for (int ix = 0; ix < n_particles * 3; ix += 3) {
      electron_pos[ix] = random() * x_len / 2 - x_len / 2;
      electron_pos[ix + 1] = random() * y_len / 2 - y_len / 2;
      electron_pos[ix + 2] = 0.0;

      positron_pos[ix] = random() * x_len / 2 - x_len / 2;
      positron_pos[ix + 1] = random() * y_len / 2 - y_len / 2;
      positron_pos[ix + 2] = 0.0;
    }
  }

  void Save_parameters_to_text(std::string filename, int form) {
    if (form == 0) {
      std::ofstream filestream(filename);
      filestream << "nx= " << nx << "\n"
                 << "ny= " << ny << "\n"
                 << "n_particles= " << n_particles << "\n"
                 << "tmax= " << tmax << "\n"
                 << "n_tsteps= " << n_tsteps << "\n"
                 << "save_rate= " << save_rate << "\n"
                 << "n_elec= " << n_elec << "\n"
                 << "n_posi= " << n_posi << "\n"
                 << "delta_x= " << delta_x << "\n"
                 << "delta_y= " << delta_y << "\n"
                 << "dt= " << dt << "\n";
      filestream.close();
    } else if (form == 1) {
      std::ofstream filestream(filename);
      filestream << "nx, "
                 << "ny, "
                 << "n_particles, "
                 << "tmax, "
                 << "n_tsteps, "
                 << "save_rate, "
                 << "n_elec, "
                 << "n_posi, "
                 << "delta_x, "
                 << "delta_y, "
                 << "dt \n"
                 << nx << ", " << ny << ", " << n_particles << ", "
                 << tmax << ", " << n_tsteps << ", " << save_rate << ", "
                 << n_elec << ", " << n_posi << ", " << delta_x << ", "
                 << delta_y << ", " << dt << "\n";
    }
  }
  void Save_current_state(std::string EM_filename,
                          std::string particle_filename,
                          std::string RFD_filename) {}

  void Append_current_state(std::string EM_filename,
                            std::string particle_filename,
                            std::string RFD_filename) {}
};
