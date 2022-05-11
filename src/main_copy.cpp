//#include <chrono>
//#include <cmath>
//#include <fstream>
//#include <iostream>
//#include <random>
//#include <string>
//#include <vector>
//#include <map>
//#include <sstream>
//
//
//const double PI = 3.14159265358979;
//
//#include "classes.hpp"
//#include "RFD.hpp"
//#include "propagation.hpp"
//#include "FDTD.hpp"
//#include "solver.hpp"
//#include "Experiment_slab.hpp"
//#include "Experiment_wave.hpp"
//#include "Experiment_gauss.hpp"
//#include "Experiment_Ez_pinch.hpp"
//#include "Experiment_Jz_pinch.hpp"
//
//int run_experiment(std::string data_dirname);
//
//// electron_temp += (gamma_e -1 ) * (c_SI * c_SI ) * ( Me_by_Kb ) / (3*n_particles);
//void run_temp_test();
//
//int main( int argc, char** argv )
//{
//  for( int ix = 0; ix <= argc; ix++ ) {
//    printf( "Argument %d = %s \n", ix, argv[ix] );
//  }
//
//  std::string data_dirname;
//  if( argc < 2 ) {
//    printf("No argument 2, using ./data as directory name for data \n");
//    data_dirname = std::string( "./data" );
//  } else{
//    data_dirname = std::string( argv[1] );
//  }
//
//  if( argc < 3 ) {
//    printf("No argument 3, using default sim \n");
//    std::string IC_string( "wave" );
//    run_experiment( data_dirname, IC_string );
//    
//  } else {
//
//  }
//  
//  
//  
//  return 0;
//}
//
//int run_experiment( std::string data_dirname, std::string IC_string )
//{
//
//
//  std::string EM_filename("/EM");
//  std::string RFD_filename("/RFD");
//  std::string particle_filename("/particle");
//  std::string current_filename("/J");
//  std::string charge_filename("/rho_q");
//  
//  // INITIAL CONDITIONS
//
//  Solver mySolver( IC );
//  std::string filename("./config.csv");
//
//  mySolver.Initialize( IC );
//  mySolver.Save_parameters_to_text(filename, 0);
//  mySolver.Save_current_state(EM_filename, particle_filename,
//                              RFD_filename, current_filename, charge_filename);
//  if( IC.use_RFD == 1 ) {
//    printf("Using RFD! \n" );
//  } else {
//    printf("Using Boris solver! \n" );
//  }
//
//  int n_tsteps = IC.n_tsteps;
//  int save_rate = IC.save_rate;
//  for (int tx = 0; tx < n_tsteps; tx++)
//  {
//    if( IC.use_RFD == 1 ) {
//      mySolver.Iterate_RFD();
//    } else {
//      mySolver.Iterate_boris();
//    }
//    
//    if( tx % ( n_tsteps / 10 ) == 0 ) { printf("tx = %d \n", tx); }
//
//    if (tx % save_rate == 0) // && tx != 0)
//    {
//      mySolver.Append_current_state(EM_filename, particle_filename,
//                                    RFD_filename, current_filename, charge_filename);
//    }
//  }
//
//  return 0;
//}