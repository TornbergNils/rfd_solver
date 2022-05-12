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

#include "RFD.hpp"
#include "EM.hpp"
#include "solver.hpp"
#include "Experiment_slab.hpp"
#include "Experiment_wave.hpp"
#include "Experiment_gauss.hpp"
#include "Experiment_Ez_pinch.hpp"
#include "Experiment_Jz_pinch.hpp"

int run_experiment( std::string&, IC_struct& );
int run_specific( std::string& );
int run_all();


int main( int argc, char** argv )
{
  std::string exp_string( argv[1] );
  if( exp_string == "-all" ) {
    run_all();
  } else {
    run_specific( exp_string );
  }
  return 0;
}

int run_specific( std::string &exp_string ) {
  
  if( exp_string == "-Jz" ) {
    Experiment_Jz_pinch IC;
    std::string data_dir( "./data/Jz" );
    run_experiment( data_dir, IC );
  }
  else if( exp_string == "-Ez" ) {
    Experiment_Ez_pinch IC;
    std::string data_dir( "./data/Ez" );
    run_experiment( data_dir, IC );
  }
  else if( exp_string == "-wave" ) {
    Experiment_wave IC;
    std::string data_dir( "./data/wave" );
    run_experiment( data_dir, IC );
  }
  else if( exp_string == "-gauss" ) {
    Experiment_Jz_pinch IC;
    std::string data_dir( "./data/gauss" );
    run_experiment( data_dir, IC );
  }
  else if( exp_string == "-slab" ) {
    Experiment_Jz_pinch IC;
    std::string data_dir( "./data/slab" );
    run_experiment( data_dir, IC );
  }
  else if( exp_string == "-langm" ) {
    Experiment_Jz_pinch IC;
    std::string data_dir( "./data/langm" );
    run_experiment( data_dir, IC );
  }
  
  return 0;
}

int run_all(){
  
  // Plane wave experiment
  std::string data_dir("./data/wave");
  Experiment_wave IC_w;
  run_experiment( data_dir, IC_w );
  
  // Gaussian experiment
  data_dir = "./data/gauss";
  Experiment_gauss IC_g;
  run_experiment( data_dir, IC_g );
  
  // Jz experiment
  data_dir = "./data/Jz";
  Experiment_Jz_pinch IC_jz;
  run_experiment( data_dir, IC_jz );

  return 0;
}

int run_experiment( 
    std::string &data_dir,
    IC_struct &IC
    )
{

  // denote if data is from boris or RFD solver
  if( IC.use_RFD == 1 ) {
    data_dir + "_R";
  } else {
    data_dir + "_B";
  }

  std::string EM_filename( data_dir + "/EM");
  std::string RFD_filename(data_dir + "/RFD");
  std::string particle_filename(data_dir + "/particle");
  std::string current_filename(data_dir + "/J");
  std::string charge_filename(data_dir + "/rho_q");

  Solver mySolver( IC );
  std::string filename( data_dir + "/config.csv");

  mySolver.Initialize( IC );
  mySolver.Save_parameters_to_text(filename, 0);
  mySolver.Save_current_state(EM_filename, particle_filename,
                              RFD_filename, current_filename, charge_filename);
  if( IC.use_RFD == 1 ) {
    printf("Using RFD! \n" );
  } else {
    printf("Using Boris solver! \n" );
  }

  int n_tsteps = IC.n_tsteps;
  int save_rate = IC.save_rate;
  for (int tx = 0; tx < n_tsteps; tx++)
  {
    if( IC.use_RFD == 1 ) {
      mySolver.Iterate_RFD();
    } else {
      mySolver.Iterate_boris();
    }
    
    if( tx % ( n_tsteps / 10 ) == 0 ) { printf("tx = %d \n", tx); }

    if (tx % save_rate == 0) // && tx != 0)
    {
      mySolver.Append_current_state(EM_filename, particle_filename,
                                    RFD_filename, current_filename, charge_filename);
    }
  }

  return 0;
}
