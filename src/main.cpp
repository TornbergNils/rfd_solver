#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <sstream>


const double PI = 3.14159265358979;

#include "classes.hpp"
#include "RFD.hpp"
#include "propagation.hpp"
#include "FDTD.hpp"
#include "solver.hpp"


int run_debug_RFD_function();
int run_debug_particle_propagation();
int run_debug_FDTD();
int run_debug_solver();

// electron_temp += (gamma_e -1 ) * (c_SI * c_SI ) * ( Me_by_Kb ) / (3*n_particles);
void run_temp_test();

int main()
{

  /*
  for( int iy = 0; iy < 10; iy++ ) {
    for( int ix = 0; ix < 10; ix++ ) {
      int idx = Get_index( ix, iy + 1, 10, 10 );
      printf("%-3d, ", idx );
    }
    printf("\n");
  }
  */
  // run_debug_RFD_function();
  // run_debug_particle_propagation();
  // run_debug_FDTD();
  run_debug_solver();
  //run_temp_test();
  return 0;
}

void Set_EM_field( EM_field_matrix&, std::map<std::string, double>& );

std::map<std::string, double> Load_parameters( std::string param_file ) {
  std::ifstream param_stream;
  param_stream.open( param_file );
  std::string line;
  std::string key;
  double value;
  std::map<std::string, double> ic_param;
  while( param_stream.good() ) {
    std::getline( param_stream, line );
    std::istringstream string_stream( line );
    string_stream >> key >> value;
    ic_param[key] = value;
  }
  param_stream.close();
  //for( const auto& elem : ic_param ){
  //  std::cout << elem.first << " = " << elem.second << "\n";
  //}
  return ic_param;
}

void run_temp_test() {
  std::string param_file( "./parameters.txt" );
  std::ifstream param_stream;
  param_stream.open( param_file );
  std::string line;
  std::string key;
  double value;
  std::map<std::string, double> ic_param;
  while( param_stream.good() ) {
    std::getline( param_stream, line );
    std::istringstream string_stream( line );
    string_stream >> key >> value;
    ic_param[key] = value;
  }
  param_stream.close();
  for( const auto& elem : ic_param ){
    std::cout << elem.first << " = " << elem.second << "\n";
  }

}

int run_debug_solver()
{
  bool param_from_file = true;

  int nx;        
  int ny;        
  double delta_x;
  double delta_y;
  double dt;     
  double tmax;   
  int n_tsteps;  
  int save_rate; 
  double c;      
  double m_e;    
  double q_e;    
  double v_thermal;
  int n_particles;
  std::map<std::string, double> ic_param;


  if( param_from_file ) {
  ic_param = Load_parameters( "./parameters.txt" );
  nx =        std::lrint( ic_param["nx"] );       
  ny =        std::lrint( ic_param["ny"] );       
  delta_x =   ic_param["dx"];  
  delta_y =   ic_param["dy"];  
  dt =        ic_param["dt"];       
  tmax =      ic_param["tmax"];     
  n_tsteps =  std::lrint( ic_param["n_tsteps"] ); 
  save_rate = std::lrint( ic_param["save_rate"] );
  c =         ic_param["c"];        
  m_e =       ic_param["m_e"];      
  q_e =       ic_param["q_e"];      
               
  v_thermal = ic_param["v_thermal"];


  n_particles = ic_param["n_particles"];  
  } else {
  // Spatial simulation param
  nx = 128;
  ny = 128;

  delta_x = 3.91e-7;
  delta_y = 3.91e-7;
  double plasma_wavelen = nx*delta_x;

  double c_cgs = 2.99792458 * 1e10;

  // Temporal simulation param
  n_tsteps = 200;
  dt = delta_x / ( 2 * c_cgs);
  tmax = n_tsteps * dt;
  save_rate = 5;
  int weight = 800;

  // Physical constants and parameters
  std::vector<double>::size_type n_particles = 100000;
  printf("particle per cell: %lf \n", n_particles / (double) (nx*ny) );

  c = 1.0;
  v_thermal = 0.047;
  double q_e_cgs_LH = std::sqrt(4*PI) *4.80320425*1e-10;
  double q_e_cgs = 4.80320425*1e-10;
  q_e = 0.30282212088;
  m_e = 511; // with c=1, in KeV ( 511 keV)
  double m_e_cgs = 9.1093819*1e-28;
  double q_by_m_cgs = q_e_cgs / m_e_cgs;
  double q_by_m = q_e / m_e;

  // conversion factors
  double len_KeV_to_cm = 0.19732705*1e-7;
  double sq_len_to_cm2 = len_KeV_to_cm * len_KeV_to_cm;
  double nat_KeV_to_Hz = 6.58*1e18;
  
  
  double density_cgs = weight * n_particles 
    / (nx * ny * delta_x * delta_y);
  
  double plasma_freq_cgs = std::sqrt(  4 * PI * density_cgs * q_e_cgs * q_e_cgs / m_e_cgs );
  
  
  printf("density cgs: %2.2e \n", density_cgs );
  printf("plasma freq cgs: %2.2e \n", plasma_freq_cgs );
  printf("plasma period cgs: %2.2e \n", 2*PI/plasma_freq_cgs );
  printf("plasma freq Hz: %2.2e \n", plasma_freq_cgs/(2*PI) );

  double nat_temp = v_thermal*v_thermal*m_e;
  double k_boltz_Kev = 8.617333262*1e-8;
  double k_boltz_erg = 1.38064 * 1e-16;
  double kelvin_temp = nat_temp / k_boltz_Kev;

  double debye_length_cgs = std::sqrt((v_thermal * c_cgs )*(v_thermal * c_cgs))  
    / std::sqrt( 2 * plasma_freq_cgs*plasma_freq_cgs );
  
  printf("debye length cgs = %2.2e \n", debye_length_cgs );
  

  double wavenum = 2*PI / ( plasma_wavelen );
  printf( "Plasma wavelen: %2.2e \n", plasma_wavelen );
  printf( "Plasma wavenum: %2.2e \n", wavenum );


  double dispersion_freq = plasma_freq_cgs 
  * std::sqrt( 1 + 3 * debye_length_cgs * debye_length_cgs
  * wavenum * wavenum);
  
  printf( "Check this is <<1: %2.2e \n", 
    wavenum*wavenum*debye_length_cgs*debye_length_cgs );
  
  
  
  // Put param into map for easy access and manipulation
  ic_param["nx"] = nx;
  ic_param["ny"] = ny;
  ic_param["delta_x"] = delta_x;
  ic_param["delta_y"] = delta_y;
  ic_param["dt"] = dt;
  ic_param["tmax"] = tmax;
  ic_param["n_tsteps"] = n_tsteps;
  ic_param["save_rate"] = save_rate;
  ic_param["c"] = c_cgs;
  ic_param["m_e"] = weight * m_e_cgs;
  ic_param["q_e"] = weight * q_e_cgs;

  ic_param["q_by_m"] = q_by_m_cgs;
  ic_param["density"] = density_cgs;
  ic_param["v_thermal"] = v_thermal * c_cgs;
  ic_param["kelvin_temp"] = kelvin_temp;
  ic_param["actual_freq"] = dispersion_freq;

  }


  double num_megabytes = (n_tsteps / save_rate * (nx * ny * 9.0 * 8.0 + n_particles * 2.0 * 12.0 * 8.0)) / 1e6;
  printf("Simulation will require %lf megabytes of harddrive space! \n", num_megabytes);

  for( const auto& elem : ic_param ){
    std::cout << elem.first << " = " << elem.second << "\n";
  }
  
  std::ofstream myStream("log.txt");
  for( const auto& elem : ic_param ){
    myStream << elem.first << " = " << elem.second << "\n";
  }

  printf( "Continue? Press any key. ");
  std::cin.get();

  std::string EM_filename("./data/EM");
  std::string RFD_filename("./data/RFD");
  std::string particle_filename("./data/particle");
  std::string current_filename("./data/J");
  std::string charge_filename("./data/rho_q");

  // INITIAL CONDITIONS FOR THE EM FIELD
  EM_field_matrix EM_IC(nx, ny);
  
  if( std::lrint( ic_param["set_wave_ic"]  ) ) {
    Set_EM_field(EM_IC, ic_param );
  } else {
    printf( "Note: no wave IC set. \n" );
  }
  

  Solver mySolver(nx, ny, n_particles, dt, n_tsteps, save_rate, delta_x,
                  delta_y, ic_param);
  std::string filename("./config.csv");

  mySolver.Initialize(EM_IC);
  mySolver.Save_parameters_to_text(filename, 0);
  mySolver.Save_current_state(EM_filename, particle_filename,
                              RFD_filename, current_filename, charge_filename);

  for (int tx = 0; tx < n_tsteps; tx++)
  {
    mySolver.Iterate_RFD();
    
    if( tx % ( n_tsteps / 10 ) == 0 ) { printf("tx = %d \n", tx); }

    if (tx % save_rate == 0) // && tx != 0)
    {
      mySolver.Append_current_state(EM_filename, particle_filename,
                                    RFD_filename, current_filename, charge_filename);
    }
  }

  return 0;
}

void Set_EM_field(EM_field_matrix &EM_IC, std::map<std::string, double> &ic_param )
{

  int nx = ic_param["nx"];
  int ny = ic_param["ny"];
  double delta_x = ic_param["dx"];
  double delta_y = ic_param["dy"];
  double wave1_A = ic_param["wave1_amplitude"];
  double wave1_k = ic_param["wave1_wavevect"];
  double wave2_A = ic_param["wave2_amplitude"];
  double wave2_k = ic_param["wave2_wavevect"];
  double Ex_A = ic_param["Ex_raw"];
  double Ex_k = ic_param["Ex_wavevect"];
  printf("Ex_A %lf \n", Ex_A);
  printf("delta_x %2.2e \n", delta_x);

  std::vector<double>::size_type num_waves = 2;
  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};

  wave_config_init[0][0] = wave1_A;      // amplitude
  wave_config_init[0][1] = wave1_k;      // wavevect = ang_freq, ok for c=1 or t=0
  wave_config_init[0][2] = 0.0; // prop angle
  wave_config_init[0][3] = 0.0;      // phase

  wave_config_init[1][0] = wave2_A;
  wave_config_init[1][1] = wave2_k;
  wave_config_init[1][2] = 0.0;
  wave_config_init[1][3] = 0.0;
  EM_wave_config config(wave_config_init);

  for (int iy = 0; iy < ny; iy++)
  {
    for (int ix = 0; ix < nx; ix++)
    {
      double x = ix * delta_x;
      double y = iy * delta_y;
      // printf( "(%.2lf, %.2lf)", x, y );
      
      EM_IC.E_x[ix + iy * nx] = 0.0; 
      EM_IC.E_y[ix + iy * nx] = Ex_A*std::sin( Ex_k * x );
      EM_IC.E_z[ix + iy * nx] = 0.0; //Ex_A*std::cos( Ex_k * x );

      EM_IC.B_x[ix + iy * nx] = 0;
      EM_IC.B_y[ix + iy * nx] = 0.0; //Ex_A*std::sin( Ex_k * x );
      EM_IC.B_z[ix + iy * nx] = Ex_A*std::sin( Ex_k * x );
      
      EM_IC.E_x[ix + iy * nx] += Get_EM_wave_component(0, config, x, y, 0);
      EM_IC.E_y[ix + iy * nx] += Get_EM_wave_component(1, config, x, y, 0);
      EM_IC.E_z[ix + iy * nx] += Get_EM_wave_component(2, config, x, y, 0);
      //printf( "%lf, ", EM_IC.E_z[ix + iy * nx] );

      EM_IC.B_x[ix + iy * nx] += Get_EM_wave_component(3, config, x, y, 0);
      EM_IC.B_y[ix + iy * nx] += Get_EM_wave_component(4, config, x, y, 0);
      EM_IC.B_z[ix + iy * nx] += Get_EM_wave_component(5, config, x, y, 0);
            
    }
    // printf( "\n" );
  }
}


int write_vector_to_binary(std::string filename, std::vector<double> vect)
{

  std::ofstream filestream(filename, std::ios::out | std::ios::binary);

  filestream.write((char *)&vect[0], vect.size() * sizeof(double));
  filestream.close();

  return 0;
}

// Takes EM matrix object and writes content to file,
// nx columns, ny rows
// adapted from old code, TODO: Refresh memory on how it works
int Write_EM_to_binary(std::string filename_E, std::string filename_B,
                       EM_field_matrix EM_field)
{

  // write all components of EM_field to file
  write_vector_to_binary(filename_E + "_x.dat", EM_field.E_x);
  write_vector_to_binary(filename_E + "_y.dat", EM_field.E_y);
  write_vector_to_binary(filename_E + "_z.dat", EM_field.E_z);

  write_vector_to_binary(filename_B + "_x.dat", EM_field.B_x);
  write_vector_to_binary(filename_B + "_y.dat", EM_field.B_y);
  write_vector_to_binary(filename_B + "_z.dat", EM_field.B_z);

  return 0;
}

double Gaussian(double x, double y) { return 1.0 * std::exp(-x * x - y * y); }

double Point_dist(double x, double y)
{
  if (x * x + y * y < 1.0)
    return 0.0;
  else
    return 1.0 / (x * x + y * y);
}

int run_debug_FDTD()
{
  int nx = 101;
  int ny = nx;
  int save_rate = 10;
  int tmax = 1000;
  std::vector<double>::size_type num_waves = 2;

  double x_max = 3.0;
  double delta_x = (2 * x_max) / (nx - 1);
  double delta_t = 0.005000;
  printf("delta_x, delta_t: %lf, %lf \n", delta_x, delta_t);
  double t = 0;

  double max_ampl = 1.0;
  // double max_ang_freq = 1.0;
  // double max_phi = 2.0 * 3.14;
  // double max_phase = 3.0;

  std::ofstream header_filestream("./data/header.csv");
  std::ofstream point_filestream("./data/time_ev_at_pt.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  // Setup rng
  // unsigned seed =
  // std::chrono::system_clock::now().time_since_epoch().count();
  // std::default_random_engine generator(seed);
  // std::uniform_real_distribution<double> distribution;
  // auto random = std::bind(distribution, generator);

  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};

  /*
  std::ofstream config_filestream("./config.csv");
  config_filestream
      << "Wave, amplitude, ang_freq, propagation angle from +x, phase \n";
  for (size_t ix = 0; ix < num_waves; ix++) {
    wave_config_init[ix][0] = max_ampl;
    wave_config_init[ix][1] = 2.0;
    wave_config_init[ix][2] = 0.78;
    wave_config_init[ix][3] = 0.0;
    config_filestream << ix << ", " << wave_config_init[ix][0] << ", "
                      << wave_config_init[ix][1] << ", "
                      << wave_config_init[ix][2] << ", "
                      << wave_config_init[ix][3] << "\n";
  }
  config_filestream.close();
  */
  wave_config_init[0][0] = max_ampl;
  wave_config_init[0][1] = 2.0;
  wave_config_init[0][2] = -PI / 2.0;
  wave_config_init[0][3] = 0.0;
  wave_config_init[1][0] = max_ampl;
  wave_config_init[1][1] = 2.0;
  wave_config_init[1][2] = 0.0;
  wave_config_init[1][3] = 0.0;
  EM_wave_config config(wave_config_init);
  // FDTD_mode TMz( nx, ny );
  // FDTD_mode TEz( nx, ny );

  // TMz mode is ( Ez, Hz, Hy )
  // TEz mode is ( Hz, Ex, Ey )

  EM_field_matrix all_modes(nx, ny);
  for (int iy = 0; iy < ny; iy++)
  {
    for (int ix = 0; ix < nx; ix++)
    {
      double x = 2.0 * ix / ((double)nx - 1) * x_max - x_max;
      double y = 2.0 * iy / ((double)ny - 1) * x_max - x_max;
      // printf( "(%.2lf, %.2lf)", x, y );

      
      all_modes.E_x[ix + iy * nx] = 0.01;
      all_modes.E_y[ix + iy * nx] = 0.0;
      all_modes.E_z[ix + iy * nx] = 0.0;

      all_modes.B_x[ix + iy * nx] = 0.0;
      all_modes.B_y[ix + iy * nx] = 0.0;
      all_modes.B_z[ix + iy * nx] = 0.0;

      /*
      all_modes.E_x[ix + iy * nx] = Gaussian(x, y);
      all_modes.E_y[ix + iy * nx] = Gaussian(x, y);
      all_modes.E_z[ix + iy * nx] = Gaussian(x, y);

      all_modes.B_x[ix + iy * nx] = Gaussian(x, y);
      all_modes.B_y[ix + iy * nx] = Gaussian(x, y);
      all_modes.B_z[ix + iy * nx] = Gaussian(x, y);
      */

      all_modes.E_x[ix + iy * nx] = Get_EM_wave_component(0, config, x, y, t);
      all_modes.E_y[ix + iy * nx] = Get_EM_wave_component(1, config, x, y, t);
      all_modes.E_z[ix + iy * nx] = Get_EM_wave_component(2, config, x, y, t);

      all_modes.B_x[ix + iy * nx] = Get_EM_wave_component(3, config, x, y, t);
      all_modes.B_y[ix + iy * nx] = Get_EM_wave_component(4, config, x, y, t);
      all_modes.B_z[ix + iy * nx] = Get_EM_wave_component(5, config, x, y, t);
    }
    // printf( "\n" );
  }
  const std::string E_filename("./data/E_data");
  const std::string B_filename("./data/B_data");

  const std::string E_history("./data/time/E_data");
  const std::string B_history("./data/time/B_data");

  // all_modes.E_z[Get_index( nx/2, ny/2, nx, ny )] = std::cos( PI/10.0 * t );

  Write_EM_to_binary(E_filename + "before", B_filename + "before", all_modes);
  for (int tx = 0; tx < tmax; tx++)
  {
    // all_modes.E_z[Get_index( nx/2, ny/2, nx, ny )] = std::cos( PI/10.0 * t );
    Update_mode(all_modes, delta_t, delta_x);
    t = t + delta_t;
    if (tx % save_rate == 0)
    {
      Write_EM_to_binary(E_history + std::to_string(t),
                         B_history + std::to_string(t), all_modes);
      point_filestream << all_modes.E_z[nx / 2 + nx * ny / 2] << ", " << t
                       << "\n";
    }
  }
  Write_EM_to_binary(E_filename + "after", B_filename + "after", all_modes);

  return 0;
}

int run_debug_particle_propagation()
{

  int nx = 64;
  int ny = 64;
  std::vector<double>::size_type n_particles = 100;
  int pp = 1;
  int sign_plus = 1;

  std::vector<double>::size_type num_waves = 4;
  double max_ampl = 1.0;
  double max_ang_freq = 1.0;
  double max_phi = 2.0 * 3.14;
  double max_phase = 3.0;

  // double c = 1.0;

  double x_max = 12.0;
  double y_max = 12.0;
  double z_max = 2.0;
  const double tmax = 20.0;
  const double dt = 0.01;

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("./data/header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  const std::string E_filename("./data/time/E_data");
  const std::string B_filename("./data/time/B_data");
  const std::string xpos_filename("./data/time/pos_x");
  const std::string ypos_filename("./data/time/pos_y");
  const std::string zpos_filename("./data/time/pos_z");
  const std::string filename_RFD("./data/time/RFD");

  std::vector<std::vector<double>> particles{n_particles,
                                             std::vector<double>(3)};
  // Setup rng
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution;
  auto random = std::bind(distribution, generator);

  std::vector<std::vector<double>> wave_config_init{num_waves,
                                                    std::vector<double>(4)};

  std::ofstream config_filestream("./config.csv");
  config_filestream
      << "Wave nr, amplitude, ang_freq, propagation angle from +x, phase \n";
  for (size_t ix = 0; ix < num_waves; ix++)
  {
    wave_config_init[ix][0] = random() * max_ampl;
    wave_config_init[ix][1] = random() * max_ang_freq;
    wave_config_init[ix][2] = random() * max_phi;
    wave_config_init[ix][3] = random() * max_phase;
    config_filestream << ix << ", " << wave_config_init[ix][0] << ", "
                      << wave_config_init[ix][1] << ", "
                      << wave_config_init[ix][2] << ", "
                      << wave_config_init[ix][3] << "\n";
  }
  config_filestream.close();
  EM_wave_config config(wave_config_init);

  for (size_t ix = 0; ix < n_particles; ix++)
  {
    particles[ix][0] = 2.0 * random() * x_max - x_max;
    particles[ix][1] = 2.0 * random() * y_max - y_max;
    particles[ix][2] = 2.0 * random() * z_max - z_max;
  }

  EM_field_matrix field_at_particles(n_particles, pp);
  EM_field_matrix EM_field(nx, ny);

  double t = 0.0;
  // Main simulation loop
  for (int tx = 0; t < tmax; t = t + dt, tx++)
  {

    // For visualizing the fields
    for (int ix = 0; ix < nx; ix++)
    {
      for (int iy = 0; iy < ny; iy++)
      {

        std::vector<double> coords = {
            2.0 * iy / (double)(ny - 1) * x_max - x_max,
            2.0 * ix / (double)(nx - 1) * y_max - y_max, 0};
        // printf( "%lf, %lf, %lf \n", coords[0] - coords[1],
        // coords[1], coords[2] );

        EM_field.E_x[ix * ny + iy] =
            Get_EM_wave_component(0, config, coords[0], coords[1], t);
        EM_field.E_y[ix * ny + iy] =
            Get_EM_wave_component(1, config, coords[0], coords[1], t);
        EM_field.E_z[ix * ny + iy] =
            Get_EM_wave_component(2, config, coords[0], coords[1], t);

        EM_field.B_x[ix * ny + iy] =
            Get_EM_wave_component(3, config, coords[0], coords[1], t);
        EM_field.B_y[ix * ny + iy] =
            Get_EM_wave_component(4, config, coords[0], coords[1], t);
        EM_field.B_z[ix * ny + iy] =
            Get_EM_wave_component(5, config, coords[0], coords[1], t);
      }
    }
    Get_EM_field_from_positions(n_particles, particles, field_at_particles, t,
                                config);
    RFD_matrix RFD_at_particles(field_at_particles, sign_plus);

    Propagate_particles(particles, RFD_at_particles, dt);

    if (tx % 10 == 0)
    {
      // Temporary solution for storing particle positions
      std::vector<double> vect(n_particles);
      for (size_t ix = 0; ix < n_particles; ix++)
      {
        vect[ix] = particles[ix][0];
      }
      write_vector_to_binary(xpos_filename + std::to_string(t), vect);
      for (size_t ix = 0; ix < n_particles; ix++)
      {
        vect[ix] = particles[ix][1];
      }
      write_vector_to_binary(ypos_filename + std::to_string(t), vect);
      for (size_t ix = 0; ix < n_particles; ix++)
      {
        vect[ix] = particles[ix][2];
      }
      write_vector_to_binary(zpos_filename + std::to_string(t), vect);
      Write_EM_to_binary(E_filename + std::to_string(t),
                         B_filename + std::to_string(t), EM_field);
      std::string pdense_filename("./data/time/power_dense");
      std::vector<double> p_dense(nx * ny);
      for (int kx = 0; kx < nx * ny; kx++)
      {
        double E_squared = EM_field.E_x[kx] * EM_field.E_x[kx] +
                           EM_field.E_y[kx] * EM_field.E_y[kx] +
                           EM_field.E_z[kx] * EM_field.E_z[kx];

        double B_squared = EM_field.B_x[kx] * EM_field.B_x[kx] +
                           EM_field.B_y[kx] * EM_field.B_y[kx] +
                           EM_field.B_z[kx] * EM_field.B_z[kx];

        p_dense[kx] = E_squared + B_squared;
      }

      write_vector_to_binary(pdense_filename + std::to_string(t), p_dense);

      write_vector_to_binary(filename_RFD + "_x" + std::to_string(t),
                             RFD_at_particles.RFD_x);
      write_vector_to_binary(filename_RFD + "_y" + std::to_string(t),
                             RFD_at_particles.RFD_y);
      write_vector_to_binary(filename_RFD + "_z" + std::to_string(t),
                             RFD_at_particles.RFD_z);
    }
  }

  return 0;
}

int run_debug_RFD_function()
{
  int nx = 37;
  int ny = 37;
  double EB_max = 3.0;
  int sign = 1;
  EM_field_matrix EM_field(nx, ny);

  // Write nx, ny, ... to "header"
  std::ofstream header_filestream("./data/header.csv");
  header_filestream << nx << ", " << ny;
  header_filestream.close();

  std::string E_filename("./data/E_data");
  std::string B_filename("./data/B_data");

  for (int ix = 0; ix < nx; ix++)
  {
    for (int iy = 0; iy < ny; iy++)
    {
      EM_field.E_x[ix * ny + iy] =
          2.0 * (iy / (double)(ny - 1)) * EB_max - EB_max;
      EM_field.E_y[ix * ny + iy] =
          2.0 * (ix / (double)(nx - 1)) * EB_max - EB_max;
      EM_field.B_x[ix * ny + iy] = 1.0;
    }
  }

  Write_EM_to_binary(E_filename, B_filename, EM_field);

  RFD_matrix RFD(EM_field, sign);

  std::vector<double> magnitudes(nx * ny);

  for (int ix = 0; ix < nx * ny; ix++)
  {
    magnitudes[ix] = std::sqrt(RFD.RFD_x[ix] * RFD.RFD_x[ix] +
                               RFD.RFD_y[ix] * RFD.RFD_y[ix] +
                               RFD.RFD_z[ix] * RFD.RFD_z[ix]);
  }

  std::string filename_RFD = "./data/RFD";
  write_vector_to_binary(filename_RFD + "_x.dat", RFD.RFD_x);
  write_vector_to_binary(filename_RFD + "_y.dat", RFD.RFD_y);
  write_vector_to_binary(filename_RFD + "_z.dat", RFD.RFD_z);
  write_vector_to_binary("./data/magnitudes.dat", magnitudes);

  return 0;
}
