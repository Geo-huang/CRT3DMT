/******************************************************************************
 *    Function       : Forward modelling and inversion for 3D MT problems     *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2019                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2019.06.10                                             *
 *    Last revision  : 2024.04.26                                             *
 ******************************************************************************/
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <assert.h>
#include "timer.h"
#include "struct_defs.h"
#include "fwd_sens_comp.h"
#include "adaptive.h"
#include "mesh3d.h"
#include "lbfgs_inv.h"

/********I. construction of struct types for reading input files *************/
// Please see the contents in file 'struct_defs.h'

/***************II. definition of global variables ***************************/
// save the output log
std::ofstream screen_output("screen_output.log");
std::ofstream para_output("input_parameter_list.log");
// the global variables will be used in the file reading functions
// for skipping useless characters when reading file
std::string skip;
// struct definition for reading startup file
Startup            startup; 
// struct definition for reading inversion parameters
Inv_Para           inv_para;
// !!! for MT3D
// struct definition for storing parameters of 3D forward modelling
// which will be directly used by external function void electric_parameters()
// and used to initialize the member variable of Fwd_Sens_Comp class
Fwd_Para_3D        fwd_para_3d;
// struct definition for storing 3D MT data
Data_3D            data_3d;
// struct definition for storing initial and priori model parameters
Init_Prior_Mod_3D  init_prior_mod_3d;

/*****III. declarations of function for structure variable initialization*****/
void read_startup_file(std::string startup_file, Startup& startup);
void input_inv_parameters(std::string inv_control_file, Inv_Para& inv_para); 

// !!!for MT3D
// for "Fwd_only"
void input_forward_parameters(std::string fwd_para_file, 
                              Fwd_Para_3D& parameter);  
// for inversion
void input_data_parameters(std::string data_para_file, Data_3D& data_3d); 
void input_ini_pri_mod_parameters(std::string ini_pri_mod_file, 
                                  Init_Prior_Mod_3D &init_prior_mod_3d);
void input_fwd_control_parameters(std::string fwd_control_file, 
                                  Fwd_Para_3D& parameter);
void check_MT3D_input_data(Data_3D& data_3d, Init_Prior_Mod_3D &init_prior_mod_3d,
                           int site_marker);


///////////////////////////////////////////////////////////////////////////////
//=========================entrance of CRT3DMT program=======================//
///////////////////////////////////////////////////////////////////////////////
void main(int argc, char** argv)
{  

  // for timekeeping
  Timer time;
  time.start();
  double elapsedtime = 0.0;
  // check the way of execution 
  if (argc <2 ) 
  {
    printf("Usage: %s .start_file_name\n", argv[0]); 
    return;    
  }  

/***************IV. initialize struct variables and objects ******************/
// only choose the needed variables to initialize
// read startup_file_name from screen
  std::string startup_file(argv[1]);
  read_startup_file(startup_file, startup);

  // check input approach
  if(startup.inv_algorithm != "Fwd_only" && startup.inv_algorithm != "L-BFGS")
  {
    std::cout<< "The input approach is incorrect, please input"
              << " Fwd_only or L-BFGS in .startup file!\n";
    screen_output<< "The input approach is incorrect, please input"
                << " Fwd_only or L-BFGS in .startup file!\n";
    para_output.close();
    std::abort();
  }

  // read different file(s) according to different input approach
  if (startup.data_method == "MT1D"){
    std::cout<< "This module is to be constructed." << std::endl;
    exit(1);
  }
  else if(startup.data_method == "MT2D"){
    std::cout<< "This module is to be constructed." << std::endl;
    exit(1);
  }
  else if(startup.data_method == "MT3D"){
    if(startup.inv_algorithm == "Fwd_only"){
      // read 3D MT forward modelling parameter related files
      input_forward_parameters(startup.fwd_para_file, fwd_para_3d);
    }
    else{
      // read inversion parameter file
      input_inv_parameters(startup.inv_control_file, inv_para);
      // read parameters in the forward control file, which is used for controling 
      // 3D MT forward modelling during inversion
      input_fwd_control_parameters(startup.fwd_control_file, fwd_para_3d);
      // read data parameter file
      input_data_parameters(startup.data_file, data_3d);
      // read initial_prior model file
      input_ini_pri_mod_parameters(startup.ini_pri_model_file, init_prior_mod_3d);
      check_MT3D_input_data(data_3d, init_prior_mod_3d, fwd_para_3d.marker);
    }
  }
  else{
    std::cout << "Input error, please check the data method in startup file!" 
              << std::endl;
    screen_output << "Input error, please check the data method in startup file!" 
                  << std::endl;
    para_output.close();
    std::exit(1);
  }
  para_output.close();
  time.stop();
  elapsedtime = time.getElapsedTime();
  //std::cout << "The duration time for reading the input files is: " 
  //          << elapsedtime << " s" << '\n';

/****************************V. choose method to proceed**********************/
  // !!!performs only a forward modelling to produce the forward data 
  if (startup.inv_algorithm == "Fwd_only"){
    if (startup.data_method == "MT1D"){
      std::cout << "This module is under construction, please wait!" << std::endl; 
      exit(1);
    }
    else if (startup.data_method == "MT2D"){
      std::cout << "This module is under construction, please wait!" << std::endl; 
      exit(1);
    }
    else if (startup.data_method == "MT3D"){
      std::cout << "The 3D MT forward data sets of "
                << fwd_para_3d.n_f << " frequencies are producing ... \n";  
      screen_output << "The 3D MT forward data sets of "
                    << fwd_para_3d.n_f << " frequencies are producing ... \n";     
      // object definition 
      Fwd_Sens_Comp fwd_sens_comp(fwd_para_3d, startup.data_method);                                        
      // compute and contaminate 3D MT forward impedance and tipper data
      fwd_sens_comp.Produce_MT3d_Synthetic_Data(startup.rel_err,
                                                startup.abs_err,
                                                startup.rel_tipper_err,
                                                startup.abs_tipper_err,
                                                startup.inv_algorithm);          
    }

    // get the duration time of forward modelling
    time.stop();
    elapsedtime = time.getElapsedTime() - elapsedtime; 
    std::cout << "Producing MT synthetic data took: " 
              << elapsedtime << " s" << " in total." << '\n';\
    screen_output << "Producing MT synthetic data took: " 
                  << elapsedtime << " s" << " in total." << '\n';
  }
  // !!!perform inversion
  else{
    // !!perform traditional inversion
    if (inv_para.adap_method == 0){
      // !do traditional Gauss-Newton inversion in model space
      if (startup.inv_algorithm == "GN"){
        std::cout << "This module is under construction\n";
        exit(1);
      } 
      // !do traditional L-BFGS inversion
      else if (startup.inv_algorithm == "L-BFGS"){
        std::cout << "Perform traditional L-BFGS inversion\n"; 
        screen_output << "Perform traditional L-BFGS inversion\n"; 
        if(startup.data_method == "MT1D"){
          std::cout << "L-BFGS_MT1D module is under construction\n";
          exit(1);
        }
        else if(startup.data_method == "MT2D"){
          std::cout << "L-BFGS_MT2D module is under construction\n";
          exit(1);
        }
        else if(startup.data_method == "MT3D"){
          // define a Meshe3D object to read the inversion mesh
          std::stringstream inv_counter;
          inv_counter << init_prior_mod_3d.counter;
          std::string inv_s = inv_counter.str();
          std::string inv_mesh_name = init_prior_mod_3d.model_name
                                    + std::string(".") + inv_s;
          Mesh3D MT3D_inv_mesh(inv_mesh_name); 
          MT3D_inv_mesh.prepare_mesh_by_reading_files();
          assert(MT3D_inv_mesh.n_elems() > 0);
          // define a LBFGS_Inv object
          LBFGS_Inv LBFGS(startup, data_3d, init_prior_mod_3d, inv_para,
                          fwd_para_3d, MT3D_inv_mesh);
          // test computation of gradient of data misfit term \phi_d
          //LBFGS.test_grad_phi_d();
          // Do 3D MT L-BFGS inversion
          LBFGS.Do_LBFGS_Inv();
          // get the time consumption of the 3D L-BFGS inversion
          time.stop();
          elapsedtime = time.getElapsedTime() - elapsedtime; 
          std::cout << "L-BFGS iterations took " 
                    << elapsedtime << " s" << " in total." << '\n';
          screen_output << "L-BFGS iterations took " 
                        << elapsedtime << " s" << " in total." << '\n';
        }
      }
      else{     
        std::cout<< "The input inversion algorithm is incorrect, please input"
                 << " L-BFGS in .startup file!\n";
       screen_output<< "The input inversion algorithm is incorrect, please input"
                    << " L-BFGS in .startup file!\n";
      }
    }
    // !!perform adaptive inversion
    // optimal algorithm is chosen in the Adaptive class
    else{
      if (startup.data_method == "MT1D"){
        std::cout << "MT1D adaptive module is under construction" 
                  << std::endl;
        exit(1);
      }
      else if (startup.data_method == "MT2D"){
        std::cout << "MT2D adaptive module is under construction" 
                  << std::endl;
        exit(1);
      }
      else if (startup.data_method == "MT3D"){
        std::cout << "Perform adaptive inversion iteration\n";
        screen_output << "Perform adaptive inversion iteration\n";
          // define a Meshe3D object to read the inversion mesh
          std::stringstream inv_counter;
          inv_counter << init_prior_mod_3d.counter;
          std::string inv_s = inv_counter.str();
          std::string inv_mesh_name = init_prior_mod_3d.model_name
                                    + std::string(".") + inv_s;
          // check whether the additional node file are given,
          // which will be used for local refinemet at the observing sites
          // for the forward mesh
          std::string a_node_file_name = init_prior_mod_3d.model_name
                                       + std::string(".") + inv_s + ".a.node";
          std::ifstream check_additional_node_file(a_node_file_name.c_str());
          if( !(check_additional_node_file.good()) ){
            std::cout << "Error, the additional node file for local refinements "
                      << "at the observing sites\nwas not found or has a "
                      << "wrong name, please add it or revise it's name!\n"
                      << "The format of the file name: 'inv_model_name'.'n'.a.node\n";
            std::exit(0);
          }

          // mesh object definition and initialization
          Mesh3D MT3D_inv_mesh(inv_mesh_name); 
          MT3D_inv_mesh.prepare_mesh_by_reading_files();
          assert(MT3D_inv_mesh.n_elems() > 0);
          // define a Adaptive object
        Adaptive adapt(startup, data_3d, init_prior_mod_3d, inv_para,
                       fwd_para_3d, MT3D_inv_mesh);
        adapt.Do_Adapt_Inv();
        time.stop();
        elapsedtime = time.getElapsedTime() - elapsedtime; 
        std::cout << "Adaptive 3D MT inversion took "
                  << elapsedtime << " s" << " in total." << '\n'; 
        screen_output << "Adaptive 3D MT inversion took "
                      << elapsedtime << " s" << " in total." << '\n'; 
      }
    }
  }

screen_output.close();
return; // exit program
}

/*****************************************************************************/
/********VI. Function definition for initializing structure variable**********/
/*****************by reading input file and post-processing*******************/
/*****************************************************************************/

void read_startup_file(std::string startup_file, Startup& startup)
{
  // function for startup file reading
  // change the string variable to a const* char 
  std::ifstream in_stream(startup_file.c_str());
  assert(in_stream.good());
  // read parameters from startup file
  in_stream >> skip >> startup.proj_name;
  para_output << "*************************Parameters from startup file*****"
              << "*********************\n" << "Project_name:\t\t" 
              << startup.proj_name << '\n';
  // data method
  in_stream >> skip >> startup.data_method;
  para_output << "Data_method:\t\t" << startup.data_method << '\n';
  if ( startup.data_method != "MT3D" ){
    std::cout << "Data method input error!\n"
              << "The keyword of data method in startup file must be 'MT3D'!\n";
    para_output.close();
    std::abort();
  }
  // function
  in_stream >> skip >> startup.inv_algorithm;
  para_output << "Function:\t\t" << startup.inv_algorithm << '\n';
  if( startup.inv_algorithm == "Fwd_only" ){
    // for forward modeling
    in_stream >> skip >> startup.fwd_para_file;
    para_output << "Fwd_para_file:\t\t" << startup.fwd_para_file << '\n';
    in_stream >> skip >> startup.rel_err;
    para_output << "Rel_err:\t" << startup.rel_err << '\n';
    if((startup.rel_err < 0.0) || (startup.rel_err > 0.5) ){
      std::cout << "Input error, please set 0 =< rel_err <= 0.5 " 
                << "in .startup file!\n";
      para_output.close();
      std::abort();
    }
    in_stream >> skip >> startup.abs_err;
    para_output << "Abs_err:\t\t" << startup.abs_err << '\n';
    if(startup.abs_err < 0.0){
      std::cout 
      << "Input error, please set abs_err > 0 in .startup file!\n";
      para_output.close();
      std::abort();
    }
    in_stream >> skip >> startup.rel_tipper_err;
    para_output << "Rel_tipper_err:\t" 
                << startup.rel_tipper_err << '\n';
    if( (startup.rel_tipper_err < 0.0) || (startup.rel_tipper_err > 0.5) ){
      std::cout 
      << "Input error, please set 0 =< rel_tipper_err <= 0.5 "
      <<  "(commonly 0 for on-shore data) in .startup file!\n";
      para_output.close();
      std::abort();
    }
    in_stream >> skip >> startup.abs_tipper_err;
    para_output << "Abs_tipper_err:\t\t" << startup.abs_tipper_err << '\n';
    if( (startup.abs_tipper_err < 0.0) || (startup.abs_tipper_err > 0.1) ){
      std::cout 
      << "Input error, please set 0 =< abs_tipper_err <= 0.1 "
      <<  "(commonly 0.01-0.03) in .startup file!\n";
      para_output.close();
      std::abort();
    }
  }
  else if( startup.inv_algorithm == "L-BFGS" ){
    // for inversion
    in_stream >> skip >> startup.data_file;
    para_output << "Data_file:\t\t" << startup.data_file << '\n';
    // for model section
    // initilizing for MT3D problem
    in_stream >> skip >> startup.ini_pri_model_file;
    para_output << "Ini_pri_model_file:\t"
		<< startup.ini_pri_model_file << '\n';
    // for inversion method 
    in_stream >> skip >> startup.inv_control_file;
    para_output << "Inv_control_file:\t" << startup.inv_control_file << '\n';
    in_stream >> skip >> startup.fwd_control_file;
    para_output << "Fwd_control_file:\t" << startup.fwd_control_file << '\n';				
  }

  // reading finished
  in_stream.close(); 
  
  return;    
}

void input_inv_parameters(std::string inv_control_file, Inv_Para& inv_para)
{
  // function for inversion parameter file reading
  para_output << "\n******************Parameters from inversion parameter"
              << " file*********************\n";
  std::ifstream in_inv_para( inv_control_file.c_str() );
  if( !(in_inv_para.good()) ){
    std::cout << "Input error: " << inv_control_file << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  // input the inversion parameters
  // !!!for adaptive inversion
  in_inv_para >> skip >> inv_para.adap_method;
  if( (inv_para.adap_method != 0) && (inv_para.adap_method != 1) ){
    std::cout << "Input error, please set adap_method to 0 or 1 in .inv file!\n";
    screen_output << "Input error, please set adap_method to 0 or 1 in .inv file!\n";
    std::abort();
  }
  in_inv_para >> skip >> inv_para.max_No_inv_unknowns;
  in_inv_para >> skip >> inv_para.max_adjust_times;
  para_output << "Adap_method:\t\t" << inv_para.adap_method << '\n';
  para_output << "Max_inv_unknowns:\t" << inv_para.max_No_inv_unknowns << '\n';
  para_output << "Max_adjust_times:\t" << inv_para.max_adjust_times << '\n';
  if(inv_para.adap_method)
    if(inv_para.max_adjust_times < 0){
      std::cout 
      << "Input error, please set max_adjust_times >= 0 in .inv file!\n";
      para_output.close();
      std::abort();
    }
  in_inv_para >> skip >> inv_para.max_iter_times;
  para_output << "Max_iter_times_p_mesh:\t" << inv_para.max_iter_times << '\n';
  if(inv_para.adap_method !=0){
    if( inv_para.max_iter_times < 2){
      std::cout << "Warining, it's better to set max_inv_iter_times_per_mesh "
                << ">= 2 in .inv file!\n";
      screen_output << "Warining, it's better to set max_inv_iter_times_per_mesh "
                    << ">= 2 in .inv file!\n";
    }
  }
  in_inv_para >> skip >> inv_para.frac_m_r;
  para_output << "Frac_of_grad_m:\t\t" << inv_para.frac_m_r << '\n';
  if((inv_para.frac_m_r < 0.0) || (inv_para.frac_m_r >1.0)){
    std::cout 
    << "Input error, please set 0=< inv_para.frac_m_r <= 1.0 in .inv file!\n";
    para_output.close();
    std::abort();
  }
  in_inv_para >> skip >> inv_para.frac_others_r;
  para_output << "Frac_of_others:\t\t" << inv_para.frac_others_r << '\n';
  if((inv_para.frac_others_r < 0.0) || (inv_para.frac_others_r > 1.0)){
    std::cout 
    << "Input error, please set 0=< inv_para.frac_others_r <= 1.0 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  // !!!for standard inversion
  // for lower and upper bounding constraints
  in_inv_para >> skip >> inv_para.a;
  para_output << "Cond_lower_bound:\t" << inv_para.a << '\n';
  if( !(inv_para.a > 0) )
  {
    std::cout << "Input error, please set inv_para.a > 0 in .inv file!\n";
    para_output.close();
    std::abort();
  }
  in_inv_para >> skip >> inv_para.b;
  para_output << "Cond_upper_bound:\t" << inv_para.b << '\n';
  if( !(inv_para.b > 0) )
  {
    std::cout << "Input error, please set inv_para.b > 0 in .inv file!\n";
    para_output.close();
    std::abort();
  }
  if(inv_para.a > inv_para.b){
    std::cout << "Input error, please set inv_para.a < inv_para.b "
              << "in .inv file!\n";
    para_output.close();
    std::abort();
  }
  in_inv_para >> skip >> inv_para.n;
  para_output << "Transformation_para:\t" << inv_para.n << '\n';
  if(inv_para.n < 0 || inv_para.n > 5){
    std::cout << "Input error, please set 0 < inv_para.n < 5 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.lambda0;
  para_output << "Lambda0:\t\t" << inv_para.lambda0 << '\n';
  if(inv_para.lambda0 < 0){
    std::cout 
    << "Input error, please set inv_para.lambda0 > 0 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.tol_rms_reduc_frac;
  para_output << "Cool_condition:\t\t" << inv_para.tol_rms_reduc_frac << '\n';
  if( (1.0e-4 > inv_para.tol_rms_reduc_frac) || (inv_para.tol_rms_reduc_frac > 0.2) ){
    std::cout 
    << "Input error, please set 1.0e-4 =< inv_para.tol_rms_reduc_frac <= 0.2 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.cool_factor;
  para_output << "Cool_factor:\t\t" << inv_para.cool_factor << '\n';
  if(inv_para.cool_factor < 1){
    std::cout 
    << "Input error, please set inv_para.cool_factor > 1 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.tol_times_cont_cool_lambda;
  para_output << "Tol_contin_cool_lambda:\t" 
              << inv_para.tol_times_cont_cool_lambda << '\n';
  if(inv_para.tol_times_cont_cool_lambda < 2){
    std::cout << "Input error, please set inv_para.tol_times_cont_cool_lambda"
    << " > 2 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.tol_lambda;
  para_output << "Tol_lambda:\t\t" << inv_para.tol_lambda << '\n';
  if( !( (0 < inv_para.tol_lambda) && ( inv_para.tol_lambda < 1 ) ) ){
    std::cout 
    << "Input error, please set inv_para.tol_lambda < 1 in .inv file!\n";
    para_output.close();
    std::abort();
  }

  in_inv_para >> skip >> inv_para.tol_iter_times;
  para_output << "Max_iteration_No:\t" << inv_para.tol_iter_times << '\n';
  if(inv_para.tol_iter_times < 1){
    std::cout 
    << "Input error, please set inv_para.tol_iter_times > 1 in .inv file!\n";
    para_output.close();
    std::abort();
  }
  if(inv_para.adap_method != 0){
    if(inv_para.tol_iter_times < 
      ((inv_para.max_adjust_times + 1) * inv_para.max_iter_times)){
      std::cout << "Input error, please set tol_iter_times > (max_adjust_times+1)"
                << "*max_iter_times in .inv file for adaptive inversion!\n";
      para_output.close();
      std::abort();
    }
  }
  in_inv_para >> skip >> inv_para.tol_rms;
  para_output << "TOl_RMS:\t\t" << inv_para.tol_rms << std::endl;

  return;
}

/*********************************!!!for MT3D*********************************/
void input_data_parameters(std::string data_para_file, Data_3D& data_3d)
{
  // function for 3d data file reading
  para_output << "\n************************Parameters from MT3D data file"
              << "*************************\n";
  std::ifstream input_data(data_para_file.c_str());
  if( !(input_data.good()) )
  {
    std::cout << "Input error: " << data_para_file << " file was not found\n";
    para_output.close();
    std::abort();
  } 
  bool using_i_omega_t_or_not = false;
  std::string line;
  std::string s_code = "Hello", s_comp = "Hello";
  // line number
  unsigned int k = 0;
  // for skipping number when reading the data file
  double skip_num = 0.0;
  // for inputting periods and frequencies
  double period = 0.0;
  double lat = 0.0, lon = 0.0, x = 0.0, y = 0.0, z = 0.0;
  while(getline(input_data, line)){
    k++;
    std::stringstream input(line);
    if( (k==1) || (k==2)){
      // 'find' returns the position of the first character of the first match;
      // if no matches were found, the function returns std::string::npos.
      if(line.find("#") == 0)
      continue;
      else{
        std::cout << "Error, the first two lines of MT3D data file should be"
                  << " the comment lines (starting with #)!";
        para_output.close();
        std::abort();
      }
    }else if(k==3) {
      input >> skip;
      input >> data_3d.data_type;
      para_output << "Data_type:\t\t" << data_3d.data_type << '\n';
      if( (data_3d.data_type != "Full_Impedance") &&
          (data_3d.data_type != "Full_Impedance_plus_Tipper") &&
          (data_3d.data_type != "Off_Diagonal_Impedance") &&
          (data_3d.data_type != "Off_Diagonal_Impedance_plus_Tipper") ) {
        std::cout<<"Error, please use full impedance data"
                 << " (with Keyword: Full_Impedance in data file) or\n"
                 << "full impedance and tipper data"
                 << " (with Keyword: Full_Impedance_plus_Tipper in data file) or\n"
                 << "off-diagonal impedance data"
                 << " (with Keyword: Off_Diagonal_Impedance in data file) or\n"
                 << "off-diagonal impedance and tipper data"
                 << " (with Keyword: Off_Diagonal_Impedance_plus_Tipper in data file)!\n";
        para_output.close();
        std::abort();
      }
      if(data_3d.data_type == "Full_Impedance") 
        data_3d.n_comps = 4; // Zxx, Zxy, Zyx, Zyy
      else if(data_3d.data_type == "Full_Impedance_plus_Tipper")
        data_3d.n_comps = 6; // Zxx, Zxy, Zyx, Zyy, Tx, Ty
      else if(data_3d.data_type == "Off_Diagonal_Impedance")
        data_3d.n_comps = 2; // Zxy, Zyx
      else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
        data_3d.n_comps = 4; // Zxy, Zyx, Tx, Ty
    }else if(k==4) {
      input >> skip;
      input >> data_3d.time_h_f; 
      para_output << "Time_factor:\t\t" << data_3d.time_h_f << '\n';
      std::size_t found = line.find("exp");
      if( found == std::string::npos ) 
      {
        std::cout<<"Error, please input time harmonic factor (with exp letters) "
                 << "adopted in data, e.g., exp(-i_omega_t) in the data file!\n";
        para_output.close();
        std::abort();
      }  
      found = line.find("-");
      if( found == std::string::npos ){
        using_i_omega_t_or_not = true;
        std::cout<< "Note: data with time dependence of +iwt is input!\n";
        screen_output<< "Note: data with time dependence of +iwt is input!\n";
      }
    }else if(k==5) {
      input >> skip;
      input >> data_3d.unit;   
      para_output << "Unit_of_data:\t\t" << data_3d.unit << '\n';
      if( !( (data_3d.unit == "[V/m]/[A/m]") ||
          (data_3d.unit == "Ohm") || (data_3d.unit == "[V/m]/[T]") 
          || (data_3d.unit == "[mV/km]/[nT]") ) ) {
        std::cout<<"Error, the data unit should be Ohm or [V/m]/[A/m] for E/H"
                << " data and [V/m]/[T] or [mV/km]/[nT] for E/B data!\n";
        para_output.close();
        std::abort();
      }   
      if((data_3d.unit == "[V/m]/[T]") || (data_3d.unit == "[mV/km]/[nT]")){
        std::cout<< "Note: E/B data in " << data_3d.unit << " is input!\n";
        screen_output<< "Note: E/B data in " << data_3d.unit << " is input!\n";
      }
    }else if(k==6) {
      input >> skip;
      input >> data_3d.orient;       
      para_output << "Orientation:\t\t" << data_3d.orient << '\n';
    }else if(k==7) {
      input >> skip;
      input >> data_3d.lat_origin;         
      input >> data_3d.lon_origin; 
      para_output << "Lat_origion:\t\t" << data_3d.lat_origin << '\n';
      para_output << "Lon_origin:\t\t"  << data_3d.lon_origin << '\n';
    }else if(k==8) {
      input >> skip;
      input >> data_3d.n_f;         
      input >> data_3d.n_sites;
      para_output << "n_frequencies:\t\t" << data_3d.n_f << '\n';
      para_output << "n_sites:\t\t"  << data_3d.n_sites << '\n';
      if( !((data_3d.n_f > 0) && (data_3d.n_sites > 0)) ) {
        std::cout<<"Error, please set n_f > 0 and n_sites > 0!\n";
        para_output.close();
        std::abort();
      }
      // the number of the data, including real and imaginary parts
      data_3d.N = data_3d.n_f * data_3d.n_sites * data_3d.n_comps * 2;
      // 2 indicates the data sets have two parts, 
      // which are imaginary and real parts
      data_3d.d.resize(data_3d.N);
      data_3d.d.setZero();
      data_3d.std_dev.resize(data_3d.N);
      data_3d.std_dev.setZero();
    }else{
      if(k < 9 + data_3d.N / 2){// avoid reading the last blank line and
        // convienient for both Z, Z & T input data
        // initilize data_3d.f and data_3d.periods with period
        input >> period;
        // set variable, automatically ranking from small to large
        (data_3d.periods).insert(period); 
        // input the observing frequencies of the first site in order of those
        // ranked in the data file
        // here, k = 9,...,8 + data_3d.N / 2
        if(k == 9)
          data_3d.f.push_back(period);
        // It means we only check one site for full impedance related data type,
        // but check two sites for off-diagonal impedace related data type
        else if(k < 9 + data_3d.n_f * 4) {  
          bool true_or_not = 0;
          for(unsigned int i = 0; i < (data_3d.f).size(); i++){
            // judging whether the period has already been pushed into the vector 
            true_or_not = true_or_not || 
                          (std::abs(data_3d.f[i] - period) < EM::TOL);
          }
          if(!true_or_not)
            data_3d.f.push_back(period); 
        }
        // input the code (i.e., name) of the site
        input >> s_code;
        (data_3d.site_code).push_back(s_code);

        // !!! input coordinates
        // Lat, Lon
        input >> lat; 
        (data_3d.lat).push_back(lat);
        input >> lon;
        (data_3d.lon).push_back(lon);
        // x, y in meter
        input >> x; 
        (data_3d.x).push_back(x);
        input >> y; 
        (data_3d.y).push_back(y);
        if(k == 9){
          (data_3d.site_x).push_back(x);
          (data_3d.site_y).push_back(y);
        }
        else{
          bool having_site = false;
          for(unsigned int i = 0; i <  (data_3d.site_x).size(); i++){
            // judging whether the site has already been recorded
            double diff_x = x - data_3d.site_x[i];
            double diff_y = y - data_3d.site_y[i];
            // The two points are approximated to be the the same 
            // in allowable around 1m error range
            if( (std::abs(diff_x) < 1.0) && (std::abs(diff_y) < 1.0) )
              having_site = true;
          }
          if(!having_site){
            (data_3d.site_x).push_back(x);
            (data_3d.site_y).push_back(y);
          }
        }
        // z in meter
        input >> z;  
        (data_3d.z).push_back(z); 
        // name of the data component, i.e., ZXX
        input >> s_comp;   
        (data_3d.data_comp).push_back(s_comp);
        // check for off-diagonal impedance data
        if ( (data_3d.data_type == "Off_Diagonal_Impedance") ||
             (data_3d.data_type == "Off_Diagonal_Impedance_plus_Tipper") ){
               if( (s_comp == "ZXX") || (s_comp == "ZYY") ){
                std::cout << "Error, the none off-diagonal component of impedance"
                          << " is input, please check the data file!\n";
                para_output.close();
                std::abort();
               }
        }

        // !!! input data
        // !! the most important parts 
        input >> data_3d.d[k - 9];
        input >> data_3d.d[k - 9 + data_3d.N / 2];
        input >> data_3d.std_dev[k - 9];
        // std_dev_Real = std_dev_Imag
        data_3d.std_dev[k - 9 + data_3d.N / 2] = data_3d.std_dev[k - 9];
      }
    }
  } // end while

  /////////////////////////////////////////////////////////////////////////////
  /***************************Post postprocessing*****************************/
  /////////////////////////////////////////////////////////////////////////////
  // check frequency information
  if(data_3d.periods.size() != data_3d.n_f){
    std::cout << "Error, please check the info of observing"
              << " frequencies in the data file!\n";
    para_output.close();
    std::abort();
  }

  if(data_3d.f.size() != data_3d.n_f){
    std::cout << "Error, please check the info of observing"
              << " frequencies in the data file!\n";
    para_output.close();
    std::abort();
  }
  else{
    for(unsigned int i = 0; i < data_3d.n_f; i++){
      // changing from period to frequency
      data_3d.f[i] = 1.0 / data_3d.f[i];
    }
  }
  // check whether having zero-valued input data
  for(unsigned int i = 0; i < data_3d.N; i++)
    if( std::abs(data_3d.d[i]) < EM::TOL ){
      std::cout << "Error: found zero-valued input data, which maybe caused"
                << " by missing\n(Tipper) data or space between lines, "
                << "please check the data file!\n"
                << "Program stopped due to input error!\n";
      screen_output << "Error: found zero-valued input data, which maybe caused"
                    << " by missing\n(Tipper) data or space between lines, "
                    << "please check the data file!\n"
                    << "Program stopped due to input error!\n";
      exit(0);
    }

  // variable definitions
  unsigned int n_f = data_3d.n_f;
  unsigned int n_s = data_3d.n_sites;
  unsigned int n_c = data_3d.n_comps;
  unsigned int half_N = data_3d.N / 2;

  // !! change +iwt data to -iwt data if needed, checked!
  // take the opposite of the imaginary parts, i.e., take the conjugate complex
  // number for the input complex data, while the errors don't need to be changed
  // This approach is only correct for MT data and can't use it for RMT data
  if(using_i_omega_t_or_not){
    double max_f = *std::max_element( (data_3d.f).begin(), (data_3d.f).end() );
    if(max_f > 1.0E5){
      std::cout << "Input error, please input RMT data with time harmonic factor"
                << " of -iwt!\n";
      para_output.close();
      std::abort();
    }
    for(unsigned int p = half_N; p < data_3d.N; p++)
      // take the conjugate
      data_3d.d[p] = - data_3d.d[p];
  }

  // !! change impedance data to Ohm ([V/m]/[A/m]) data if needed, checked!
  //    assume: mu_r = 1
  if(data_3d.data_type == "Full_Impedance"){
    if(data_3d.unit == "[V/m]/[T]"){
      // change E/B data in [V/m]/[T] to E/H data in Ohm
      data_3d.d       = EM::MU0 * data_3d.d;
      data_3d.std_dev = EM::MU0 * data_3d.std_dev;
    }
    else if(data_3d.unit == "[mV/km]/[nT]"){
      // change E/B data in [mV/km]/[nT] to E/H data in Ohm
      data_3d.d       = EM::MU0 * 1000.0 * data_3d.d;
      data_3d.std_dev = EM::MU0 * 1000.0 * data_3d.std_dev;
    }
  }
  else if(data_3d.data_type == "Full_Impedance_plus_Tipper"){
    // please note: the tipper data are dimesionless, so we only 
    // need to do the unit conversion for impedance data
    // the number of impedance components per site and per f
    unsigned int n_c_Z    = 4;
    if(data_3d.unit == "[V/m]/[T]"){
      // change E/B data in [V/m]/[T] to E/H data in Ohm
      for(unsigned int i = 0; i <  n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c_Z; k++){
            unsigned int temp_index = i*n_f*n_c_Z + j*n_c_Z + k;
            // real parts
            data_3d.d[temp_index]      
            = EM::MU0 * data_3d.d[temp_index];
            data_3d.std_dev[temp_index]
            = EM::MU0 * data_3d.std_dev[temp_index];
            // imaginary parts
            data_3d.d[temp_index + half_N]      
            = EM::MU0 * data_3d.d[temp_index + half_N];
            data_3d.std_dev[temp_index + half_N]
            = EM::MU0 * data_3d.std_dev[temp_index + half_N];
          }
    }
    else if(data_3d.unit == "[mV/km]/[nT]"){
      // change E/B data in [mV/km]/[nT] to E/H data in Ohm
      for(unsigned int i = 0; i <  n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c_Z; k++){
            unsigned int temp_index = i*n_f*n_c_Z + j*n_c_Z + k;
            // real parts
            data_3d.d[temp_index]      
            = EM::MU0 * 1000.0 * data_3d.d[temp_index];
            data_3d.std_dev[temp_index]
            = EM::MU0 * 1000.0 * data_3d.std_dev[temp_index];
            // imaginary parts
            data_3d.d[temp_index + half_N]      
            = EM::MU0 * 1000.0 * data_3d.d[temp_index + half_N];
            data_3d.std_dev[temp_index + half_N]
            = EM::MU0 * 1000.0 * data_3d.std_dev[temp_index + half_N];
          }
    }
  }
  else if(data_3d.data_type == "Off_Diagonal_Impedance"){
    if(data_3d.unit == "[V/m]/[T]"){
      // change E/B data in [V/m]/[T] to E/H data in Ohm
      data_3d.d       = EM::MU0 * data_3d.d;
      data_3d.std_dev = EM::MU0 * data_3d.std_dev;
    }
    else if(data_3d.unit == "[mV/km]/[nT]"){
      // change E/B data in [mV/km]/[nT] to E/H data in Ohm
      data_3d.d       = EM::MU0 * 1000.0 * data_3d.d;
      data_3d.std_dev = EM::MU0 * 1000.0 * data_3d.std_dev;
    }
  }
  else{// for data type of "Off_Diagonal_Impedance_plus_Tipper"
    // please note: the tipper data are dimesionless, so we only 
    // need to do the unit conversion for diagonal impedance data
    // the number of diagonal impedance components per site and per f
    unsigned int n_c_Z    = 2;
    if(data_3d.unit == "[V/m]/[T]"){
      // change E/B data in [V/m]/[T] to E/H data in Ohm
      for(unsigned int i = 0; i <  n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c_Z; k++){
            unsigned int temp_index = i*n_f*n_c_Z + j*n_c_Z + k;
            // real parts
            data_3d.d[temp_index]      
            = EM::MU0 * data_3d.d[temp_index];
            data_3d.std_dev[temp_index]
            = EM::MU0 * data_3d.std_dev[temp_index];
            // imaginary parts
            data_3d.d[temp_index + half_N]      
            = EM::MU0 * data_3d.d[temp_index + half_N];
            data_3d.std_dev[temp_index + half_N]
            = EM::MU0 * data_3d.std_dev[temp_index + half_N];
          }
    }
    else if(data_3d.unit == "[mV/km]/[nT]"){
      // change E/B data in [mV/km]/[nT] to E/H data in Ohm
      for(unsigned int i = 0; i <  n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c_Z; k++){
            unsigned int temp_index = i*n_f*n_c_Z + j*n_c_Z + k;
            // real parts
            data_3d.d[temp_index]      
            = EM::MU0 * 1000.0 * data_3d.d[temp_index];
            data_3d.std_dev[temp_index]
            = EM::MU0 * 1000.0 * data_3d.std_dev[temp_index];
            // imaginary parts
            data_3d.d[temp_index + half_N]      
            = EM::MU0 * 1000.0 * data_3d.d[temp_index + half_N];
            data_3d.std_dev[temp_index + half_N]
            = EM::MU0 * 1000.0 * data_3d.std_dev[temp_index + half_N];
          }
    }
  }

  // !! rearranging data and standard deviation: 
  // loop frequency --> loop sites --> loop components
  EM::VectorXD d_temp = data_3d.d, std_dev_temp = data_3d.std_dev;
  if(data_3d.data_type == "Full_Impedance"){
    // ! full impedance data
    for(unsigned int i = 0; i <  n_f; i++)
      for(unsigned int j = 0; j < n_s; j++)
        for(unsigned int k = 0; k < n_c; k++)
        {
          data_3d.d[i*n_s*n_c + j*n_c + k]                
          = d_temp[j*n_f*n_c + i*n_c + k];

          data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
          = d_temp[j*n_f*n_c + i*n_c + k + half_N];

          data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
          = std_dev_temp[j*n_f*n_c + i*n_c + k];

          data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
          = std_dev_temp[j*n_f*n_c + i*n_c + k + half_N];
        }
  }
  else if(data_3d.data_type == "Full_Impedance_plus_Tipper"){
    // ! impedance plus tipper data
    // the number of all the real parts of impedance components
    // = the number of all the complex-valued impedance components
    unsigned int N_Z_comp = n_f * n_s * 4;
    // the number of impedance components per site and per frequency
    unsigned int n_c_Z    = 4;
    // the number of tipper components per site and per frequency
    unsigned int n_c_T    = 2;
    for(unsigned int i = 0; i <  n_f; i++)
      for(unsigned int j = 0; j < n_s; j++)
        for(unsigned int k = 0; k < n_c; k++)
        {
          if(k < 4){
            // ! impedance responses
            // real part
            data_3d.d[i*n_s*n_c + j*n_c + k]                
            = d_temp[j*n_f*n_c_Z + i*n_c_Z + k];
            
            // imaginary part
            data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
            = d_temp[j*n_f*n_c_Z + i*n_c_Z + k + half_N];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
            = std_dev_temp[j*n_f*n_c_Z + i*n_c_Z + k];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
            = std_dev_temp[j*n_f*n_c_Z + i*n_c_Z + k + half_N];
          }
          else{
            // tipper responses, (k-4) denotes the index of the 
            // component of the tipper
            data_3d.d[i*n_s*n_c + j*n_c + k]                
            = d_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-4];

            data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
            = d_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-4 + half_N];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
            = std_dev_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-4];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
            = std_dev_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-4 + half_N];
          }
        }
  }
  else if(data_3d.data_type == "Off_Diagonal_Impedance"){
    // ! off-diagonal impedance data
    for(unsigned int i = 0; i <  n_f; i++)
      for(unsigned int j = 0; j < n_s; j++)
        for(unsigned int k = 0; k < n_c; k++)
        {
          data_3d.d[i*n_s*n_c + j*n_c + k]                
          = d_temp[j*n_f*n_c + i*n_c + k];

          data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
          = d_temp[j*n_f*n_c + i*n_c + k + half_N];

          data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
          = std_dev_temp[j*n_f*n_c + i*n_c + k];

          data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
          = std_dev_temp[j*n_f*n_c + i*n_c + k + half_N];
        }
  }
  else{ // for data type of "Off_Diagonal_Impedance_plus_Tipper"
    // ! off-diagonal impedance plus tipper data
    // the number of all the real parts of off-diagonal impedance components
    // = the number of all the complex-valued impedance off-diagonal components
    unsigned int N_Z_comp = n_f * n_s * 2;
    // the number of off-diagonal impedance components per site and per frequency
    unsigned int n_c_Z    = 2;
    // the number of tipper components per site and per frequency
    unsigned int n_c_T    = 2;
    for(unsigned int i = 0; i <  n_f; i++)
      for(unsigned int j = 0; j < n_s; j++)
        for(unsigned int k = 0; k < n_c; k++)
        {
          if(k < 2){
            // ! impedance responses
            // real part
            data_3d.d[i*n_s*n_c + j*n_c + k]                
            = d_temp[j*n_f*n_c_Z + i*n_c_Z + k];
            
            // imaginary part
            data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
            = d_temp[j*n_f*n_c_Z + i*n_c_Z + k + half_N];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
            = std_dev_temp[j*n_f*n_c_Z + i*n_c_Z + k];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
            = std_dev_temp[j*n_f*n_c_Z + i*n_c_Z + k + half_N];
          }
          else{
            // tipper responses, (k-2) denotes the index of the 
            // component of the tipper
            data_3d.d[i*n_s*n_c + j*n_c + k]                
            = d_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-2];

            data_3d.d[i*n_s*n_c + j*n_c + k + half_N]       
            = d_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-2 + half_N];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k]          
            = std_dev_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-2];

            data_3d.std_dev[i*n_s*n_c + j*n_c + k + half_N] 
            = std_dev_temp[N_Z_comp + j*n_f*n_c_T + i*n_c_T + k-2 + half_N];
          }
        }
  }

  /////////////////////////////////////////////////////////////////////////////
  /*******************Output data with +iwt and unit of ohm*******************/
  /////////////////////////////////////////////////////////////////////////////
  //if( (data_3d.unit == "[V/m]/[T]") || (data_3d.unit == "[mV/km]/[nT]")
  //    ||  using_i_omega_t_or_not )
  {
    std::string file_name = "iwt_ohm_" + data_para_file;
    std::ofstream output_data(file_name.c_str());
    // ! output header part
    output_data << "# 3D MT input data with time factor of exp(+i_omega_t)"
                << " and unit of [V/m]/[A/m] (ohm)\n"
                << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m)"
                << " Component Real Imag Error \n"
                << "> MT3D \n"
                << "> " << data_3d.data_type << " \n"
                << "> exp(+i_omega_t) \n"
                << "> [V/m]/[A/m] \n"
                << "> " << data_3d.orient << "\n"
                << "> " << data_3d.lat_origin 
                << "  " << data_3d.lon_origin << "\n"
                << "> " << data_3d.n_f << "  " 
                << data_3d.n_sites << '\n';

  // ! output data part
    if( (data_3d.data_type == "Full_Impedance") ||
        (data_3d.data_type == "Off_Diagonal_Impedance") ){
      for(unsigned int i = 0; i < n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c; k++){
            unsigned int index = j * n_s * n_c + i * n_c + k;
            unsigned int other_index = i * n_f * n_c + j * n_c + k;
            // output the (changed) data
            output_data << setiosflags(std::ios::scientific)
                        << 1.0 / data_3d.f[j] << '\t' 
                        // Code 
                        << data_3d.site_code[other_index] << '\t'
                        // CG_Lat CG_Lon X(m) Y(m) Z(m)
                        << resetiosflags(std::ios::scientific)
                        << std::setprecision(8)
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << data_3d.lat[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << data_3d.lon[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << data_3d.x[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << data_3d.y[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << data_3d.z[other_index] << '\t'
                        // component
                        << data_3d.data_comp[other_index] << '\t'
                        << setiosflags(std::ios::scientific)
                        << std::setprecision(6)
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // real part
                        << data_3d.d(index) << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // imaginary part
                        << -1.0 * data_3d.d(index + half_N) << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // error
                        << data_3d.std_dev(index) << '\n';
            assert( std::abs(data_3d.std_dev(index) - 
                             data_3d.std_dev(index+half_N)) < 1.0e-7);
          }
    }
    else // for data type of "Full_Impedance_plus_Tipper"
    {    // and "Off_Diagonal_Impedance_plus_Tipper"
      // ! store data into respective vectors
      // real and imaginary part of Z and T (changed) data
      std::vector< double > real_d_Z, real_d_T, imag_d_Z, imag_d_T;
      // error and period
      std::vector< double > err_Z, err_T, period_Z, period_T;
      for(unsigned int i = 0; i < n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c; k++)
          {
            unsigned int index = j * n_s * n_c + i * n_c + k;
            unsigned int other_index = i * n_f * n_c + j * n_c + k;
            unsigned int remainder = other_index % n_c;
            assert( std::abs(data_3d.std_dev(index) - 
                             data_3d.std_dev(index+half_N)) < 1.0e-7);
            // n_c - 3 = 3 for data type of "Full_Impedance_plus_Tipper"
            // n_c - 3 = 1 for data type of "Off_Diagonal_Impedance_plus_Tipper"
            if(remainder <= (n_c - 3))
            {
              // (changed) data
              real_d_Z.push_back( data_3d.d(index) );
              imag_d_Z.push_back( data_3d.d(index + half_N) );
              // error
              err_Z.push_back( data_3d.std_dev(index) );
              // period
              period_Z.push_back(1.0 / data_3d.f[j]);
            }
            else
            {
              // (changed) data
              real_d_T.push_back( data_3d.d(index) );
              imag_d_T.push_back( data_3d.d(index + half_N) );
              // error
              err_T.push_back( data_3d.std_dev(index) );
              // period
              period_T.push_back(1.0 / data_3d.f[j]);
            }
          }
      // ! merge Z and T related vectors into one vector
      std::vector< double > real_d, imag_d, err, period;
      for(unsigned int i = 0; i < period_Z.size(); i++){
        real_d.push_back( real_d_Z[i] );
        imag_d.push_back( imag_d_Z[i] );
        err.push_back( err_Z[i] );
        period.push_back( period_Z[i] );
      }
      for(unsigned int i = 0; i < period_T.size(); i++){
        real_d.push_back( real_d_T[i] );
        imag_d.push_back( imag_d_T[i] );
        err.push_back( err_T[i] );
        period.push_back( period_T[i] );
      }      

      // ! output data
      for(unsigned int i = 0; i < half_N; i++){
        // output the (changed) data
        output_data << setiosflags(std::ios::scientific)
                    << period[i] << '\t'
                    << data_3d.site_code[i] << '\t'
                    // CG_Lat CG_Lon X(m) Y(m) Z(m)
                    << resetiosflags(std::ios::scientific)
                    << std::setprecision(8)
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << data_3d.lat[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << data_3d.lon[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << data_3d.x[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << data_3d.y[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << data_3d.z[i] << '\t'
                    // component
                    << data_3d.data_comp[i] << '\t'
                    << setiosflags(std::ios::scientific)
                    << std::setprecision(6)
                    << setiosflags(std::ios::right)
                    << std::setw(15)
                    // real part
                    << real_d[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(15)
                    // imaginary part
                    << -1.0 * imag_d[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(15)
                    // error
                    << err[i] << '\n';
      }
    }
    output_data.close();
  }

  return; //3D MT data file reading is finished
} 

void input_ini_pri_mod_parameters(std::string ini_pri_mod_file, 
                                  Init_Prior_Mod_3D &init_prior_mod_3d)
{
  // function for 3D MT initial and priori model related file reading
  // and some processing
  para_output << "\n*************************Parameters from ini_ref"
              << " file**************************\n";
  bool having_marker_9999999 = 0;
  std::ifstream input(ini_pri_mod_file.c_str());
  if( !(input.good()) )
  {
    std::cout << "Input error: " << ini_pri_mod_file 
              << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  // skip comment line
  input >> skip;
  input >> init_prior_mod_3d.model_name;
  input >> init_prior_mod_3d.counter;
  para_output << "Inv_model_name:\t\t" << init_prior_mod_3d.model_name << '\n';
  para_output << "Counter:\t\t" << init_prior_mod_3d.counter << '\n';
  if(init_prior_mod_3d.counter < 1){
    std::cout<<"Error, please set counter = 1,2,..,\n";
    para_output.close();
    std::abort();
  }
  input >> init_prior_mod_3d.n_regions;
  para_output << "n_regions:\t\t" << init_prior_mod_3d.n_regions << '\n';
  if(init_prior_mod_3d.n_regions < 1){
    std::cout<<"Error, please set n_regions >= 1!\n";
    para_output.close();
    std::abort();
  }

  // initialize map variable region_table_ini by reading corresponding file
  // skip comment line
  input >> skip;
  if(skip == "1"){
    std::cout 
    << "Input error: please input a string for comment before inputting\n"
    << "the parameters of the initial model in .ini_ref file! "
    << "Program stopped!\n";
    exit(0);
  }
  for(int i = 0; i < init_prior_mod_3d.n_regions; i++){
    int marker;
    std::vector<double> temp(3);
    input >> marker;
    for(int j=0; j<3; j++) 
      input >> temp[j];
    init_prior_mod_3d.region_table_ini[marker] = temp;
    if(marker == 9999999){
      having_marker_9999999 = 1;
      if( init_prior_mod_3d.region_table_ini[marker][0] > 1.e-10 )
      {
        std::cout << "Input error: the conductivity of air > 1.e-10, "
                  << "please reset it!\n";
        para_output.close();
        exit(0);
      }
    }
  }
  if(!having_marker_9999999){
    std::cout << "Error, the regional marker of air space should be 9999999,"
              << " please check the initial model file!\n";
    para_output.close();
    exit(0);
  }

  // initialize map variable region_table_ref reading corresponding file
  // skip comment line
  input >> skip;
  if(skip == "1"){
    std::cout 
    << "Input error: please input a string for comment before inputting\n"
    << "the parameters of the reference model in .ini_ref file! "
    << "Program stopped!\n";
    exit(0);
  }
  for(int i = 0; i < init_prior_mod_3d.n_regions; i++){
    std::string using_default;
    int marker;
    // choose one of the two way to initialize region_table_ini map
    if(i == 0){
      input >> using_default;
      if(using_default == "default"){
        init_prior_mod_3d.region_table_ref = init_prior_mod_3d.region_table_ini;
        break;
      }
      else if(using_default == "1"){
        having_marker_9999999 = 0;
        marker = std::stoi(using_default);
        assert(marker == 1);
      }
      else{
        std::cout << "Input error, please check the .ini_.ref file!\n"
                  << " Program stopped!\n";
        exit(0);
      }
    }
    else
      input >> marker;
    std::vector<double> temp(3);
    for(int j=0; j<3; j++) 
      input >> temp[j];
    init_prior_mod_3d.region_table_ref[marker] = temp;
    if(marker == 9999999){
      having_marker_9999999 = 1;
      if( init_prior_mod_3d.region_table_ref[marker][0] > 1.e-10 )
      {
        std::cout << "Input error: the conductivity of air > 1.e-10, "
                  << "please reset it!\n";
        para_output.close();
        exit(0);
      }
    }
  }
  if(!having_marker_9999999){
    std::cout << "Error, the regional marker of air space should be 9999999,"
              << " please check the priori model file!\n";
    para_output.close();
    exit(0);
  }

  // for reading inversion mesh
  std::stringstream inv_counter;
  inv_counter << init_prior_mod_3d.counter;
  std::string inv_s = inv_counter.str();
  std::string inv_mesh_name = init_prior_mod_3d.model_name
                              + std::string(".") + inv_s;
  Mesh3D inv_mesh(inv_mesh_name); 
  inv_mesh.prepare_mesh_by_reading_files();
  assert(inv_mesh.n_elems() > 0);

  assert(!init_prior_mod_3d.region_table_ini.empty());
  assert(!init_prior_mod_3d.region_table_ref.empty());
  typedef std::map<int, std::vector<double> >::iterator IT;
  unsigned int k = 0;
  for (int e = 0; e < inv_mesh.n_elems(); e++){
    // get e-th tetrahedral element 
    Tet* inv_tet= static_cast<Tet*>(&(inv_mesh.get_elem(e)));
    assert(e == inv_tet->get_id());
    const unsigned int inv_tet_marker = inv_tet->get_tet_marker();
    // make sure each tetrahedral element in the subsuface inversion 
    // domain (having fee parameter) has a unique marker, ranked from
    // 1 to M, M represents the number of the inversion unknowns
    assert(inv_tet_marker != 0);

    // make sure each inv_tet_marker has a value corresponding to it
    IT it_ini = init_prior_mod_3d.region_table_ini.find(inv_tet_marker);
    // found the inv_tet_marker in init_prior_mod_3d.region_table_ini map
    assert(it_ini != init_prior_mod_3d.region_table_ini.end());
    // make sure each fwd_tet_marker has a value corresponding to it
    IT it_ref = init_prior_mod_3d.region_table_ref.find(inv_tet_marker);
    // found the inv_tet_marker in init_prior_mod_3d.region_table_ref map
    assert(it_ref != init_prior_mod_3d.region_table_ref.end());

    // dealing with the tetrahedral elements in the inversion domain
    if( (inv_tet_marker != 9999999) && (inv_tet_marker != 6666666) ){
      // initialize set variable sub_node_id, where the node id
      // will be automatically ranked from small to large
      for(unsigned int i = 0; i < 4; i++){
        unsigned int node_id = inv_tet->get_node_id(i);
        (init_prior_mod_3d.sub_node_id).insert(node_id);
      }

      // initialize vector m_ini and m_ini_id
      std::vector<double>& parameter_ini = (*it_ini).second;
      assert(parameter_ini.size()==3);  
      (init_prior_mod_3d.m_ini).push_back(parameter_ini[0]);
      (init_prior_mod_3d.m_ini_id).push_back( e );

      // initialize vector m_ref
      std::vector<double>& parameter_ref = (*it_ref).second;
      assert(parameter_ref.size()==3); 
      init_prior_mod_3d.m_ref.push_back(parameter_ref[0]);

      // guarantee (1 + index of the entry in m_ini vector) equals 
      // to the marker of the corresponding tetrahedral element
      // That relationship will be used to do the parameter passing
      // from inversion mesh to forward mesh when doing inversion
      if( (k + 1) != inv_tet_marker ){
        std::cout << "Input error: duplicated markers are used for tets in inversion domain!\n"
                  << "Note: each tet in the inversion domain need a unique marker and the\n"
                  << "marker need to be ranked from 1 to M, please check the inversion mesh!\n"
                  << "Program stopped due to input error!\n"; 
        para_output.close();
        exit(0);
      }
      k++;
    }
  }

  // ! initialize vector fwd_tet_id_in_mj
  // size definition
  (init_prior_mod_3d.fwd_tet_id_in_mj).resize( init_prior_mod_3d.m_ini.size() );
  // read forward mesh,               
  std::stringstream fwd_counter;
  fwd_counter << fwd_para_3d.starting_counter;
  std::string fwd_s = fwd_counter.str();
  std::string fwd_mesh_name = fwd_para_3d.starting_model[0] 
                              + std::string(".") + fwd_s;
  Mesh3D fwd_mesh(fwd_mesh_name);
  fwd_mesh.prepare_mesh_by_reading_files();
  assert(fwd_mesh.n_elems() > 0);
  // At this version, all frequencies share one forward mesh
  assert(fwd_para_3d.starting_model.size() == 1);
  for(unsigned int i = 0; i < fwd_mesh.n_elems(); i++){
    Tet* fwd_tet= static_cast<Tet*>(&(fwd_mesh.get_elem(i)));
    assert(i == fwd_tet->get_id());
    const unsigned int fwd_tet_marker = fwd_tet->get_tet_marker();
    if( (fwd_tet_marker != 9999999) && (fwd_tet_marker != 6666666) ){
      if( fwd_tet_marker > init_prior_mod_3d.m_ini.size() ){
        std::cout << "Input error, please check forward mesh (the nested forward mesh\n"
                  << "should have the same regional markers as inversion mesh)!\n";
        para_output.close();
        std::abort();
      }
      (init_prior_mod_3d.fwd_tet_id_in_mj[fwd_tet_marker - 1]).push_back( i );
    }
  }
  // guarantee the markers of the forward and those of inversion mesh are the same
  for(unsigned int j = 0; j < (init_prior_mod_3d.fwd_tet_id_in_mj).size(); j++)
    if( (init_prior_mod_3d.fwd_tet_id_in_mj[j]).size() == 0 ){
      std::cout << "Input error, please check forward mesh (the nested forward mesh\n"
                << "should have the same regional markers as inversion mesh)!\n";
      para_output.close();
      std::abort();
    }
  
  return;
}

void input_fwd_control_parameters(std::string fwd_control_file, 
                                  Fwd_Para_3D& parameter)
{
  // !!! read the parameter info offered by forward control file into parameter
  // function for control parameter of 3D forward modeling related file reading
  para_output << "\n*******************Parameters from MT3D_fwd_control_para "
              << "file******************\n";
  std::ifstream in_stream(fwd_control_file.c_str());
  if( !(in_stream.good()) ){
    std::cout << "Input error: " << fwd_control_file << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  
  // default using standard fem
  parameter.algorithm = 0;
  // for convenient, only one model name is neeeded to put in fwd_control_file
  // In the current version, all the frequencies share the same forward mesh,
  // hence only one model name is needed.
  // Later, the size of the 'starting_model' vector will be set to 
  // the number of observing frequencies shown in data, when we want to use Openmp
  // and mesh refinement or different mesh of different frequencies, also its elements 
  // should be assigned with different names (string type)
  (parameter.starting_model).resize(1);
  // skip comment line
  in_stream >> skip;
  in_stream >> parameter.starting_model[0];
  in_stream  >> parameter.starting_counter; 
  para_output << "\n" << "starting_model:\t" << parameter.starting_model[0];
  para_output << "\n" << "starting_counter:\t" << parameter.starting_counter;
  if(parameter.starting_counter < 1) {
    std::cout<<"Error, please set starting_counter = 1,2,..,\n";
    para_output.close();
    std::abort();
  }
  in_stream >> parameter.theta; 
  para_output << "\n" << "theta:\t" << parameter.theta <<"\n";
  if(parameter.theta<0||parameter.theta>90) {
    std::cout<<"Error, please set 0 < theta < 90!\n";
    para_output.close();
    std::abort();
  }

  // read parameters for boundary conditions
  // skip comment line
  in_stream >> skip;
  in_stream >> parameter.n_layer; 
  para_output << "n_layers:\t" << parameter.n_layer <<"\n";
  if(parameter.n_layer<2) {
    std::cout<<"Error, please set n_layer >= 2!\n";
    para_output.close();
    std::abort();
  }
  parameter.EP.resize(parameter.n_layer);
  for(int i=0; i<parameter.n_layer; i++) {
    parameter.EP[i].resize(4);
    for(int j=0; j<4; j++) {
      in_stream >> parameter.EP[i][j];
      para_output << parameter.EP[i][j] <<"\t";
    }
    para_output<<"\n";
  }
  
  // read sites information
    // skip comment line
  in_stream >> skip;
  in_stream >> parameter.site_file;
  para_output << "site_file:\t" << parameter.site_file <<"\n";
  para_output << "\n**************************Parameters from site file"
              << "****************************\n";
  std::ifstream site_stream (parameter.site_file.c_str());
  if( !(site_stream.good()) )
  {
    std::cout << "Input error: " << parameter.site_file 
              << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  site_stream >> parameter.n_sites; 
  site_stream >> parameter.marker; 
  para_output << "n_sites:\t" << parameter.n_sites << "\n";  
  para_output << "site_marker:\t" << parameter.marker << "\n";  
  if(parameter.n_sites<0) {
    std::cout << "n_sites:\t"<<parameter.n_sites<<"\n";  
    std::cout<<"Error, please set n_sites>0 in file:\t"<<parameter.site_file;
    std::cout<<"!\n";
    para_output.close();
    std::abort();
  }
  parameter.u.resize(parameter.n_sites);
  parameter.v.resize(parameter.n_sites);
  parameter.mu.resize(parameter.n_sites);
  for(int i=0; i<parameter.n_sites; i++) {
    double x1 = 0., x2 = 0., x3 = 0.;
    site_stream >> x1 >> x2 >>x3; 
    parameter.u[i] = Point(x1,x2,x3);  // u

    // u can be only equals to (1., 0., 0.) for this version!
    Point u_temp = Point(1., 0., 0.) - Point(x1,x2,x3);
    if(u_temp.size() > 0.){
      std::cout<< "Note: the local coordinate system with u = ("
               << Point(x1,x2,x3) << ") for site " << i +1 << " is used!\n";
      screen_output<< "Note: the local coordinate system with u = ("
                   << Point(x1,x2,x3) << ") for site " << i +1 << " is used!\n";
    }
  }
  for(int i=0; i<parameter.n_sites; i++) {
    double x1 = 0., x2 = 0., x3 = 0.;
    site_stream >> x1 >> x2 >>x3; 
    parameter.v[i] = Point(x1,x2,x3);  // v

    // v can be only equals to (0., 1., 0.) for this version!
    Point v_temp = Point(0., 1., 0.) - Point(x1,x2,x3);
    if(v_temp.size() > 0.){
      std::cout<< "Note: the local coordinate system with v = ("
               << Point(x1,x2,x3) << ") for site " << i+1 << " is used!\n";
      screen_output<< "Note: the local coordinate system with v = ("
                   << Point(x1,x2,x3) << ") for site " << i+1 << " is used!\n";
    }
  }
  for(int i=0; i<parameter.n_sites; i++) {
    double temp = 0.;
    site_stream >> temp;
    if( !(std::abs(temp) > 0.) ){
      std::cout << "Error, negative or zero mu_r emerges!\n";
      para_output.close();
      std::abort();
    }
    parameter.mu[i] = temp;            // mu
  }
  // output u, v and mu_r belonging to the sites into file
  for(int i=0; i<parameter.n_sites; i++) {
    para_output << "u:\t" << parameter.u[i](0) << "\t"
                        << parameter.u[i](1) << "\t"
                        << parameter.u[i](2) << "\t";
    para_output << "v:\t" << parameter.v[i](0) << "\t"
                        << parameter.v[i](1) << "\t"
                        << parameter.v[i](2) << "\t";
    para_output << "mu:\t" << parameter.mu[i] << "\n";
  }
  
  return;
}

void check_MT3D_input_data(Data_3D& data_3d, Init_Prior_Mod_3D &init_prior_mod_3d,
                           int site_marker)
{
  // !!! read inversion mesh
  std::stringstream inv_counter;
  inv_counter << init_prior_mod_3d.counter;
  std::string inv_s = inv_counter.str();
  std::string inv_mesh_name = init_prior_mod_3d.model_name
                              + std::string(".") + inv_s;
  Mesh3D inv_mesh(inv_mesh_name); 
  inv_mesh.prepare_mesh_by_reading_files();
  assert(inv_mesh.n_elems() > 0);

  // !!! store the observing nodes into site vector in 
  //     the order of those ranked in the poly file
  std::vector< Node > sites;
  const unsigned int n_nodes = inv_mesh.n_nodes();
  for(unsigned int i = 0; i < n_nodes; i++){
    Node *a = &(inv_mesh).get_node(i);
    if(a->get_marker() == site_marker)
      sites.push_back(*a); 
  }

  if(sites.size() != (data_3d.site_x).size())
  {
    std::cout << "The number of sites shown in poly file is: " << sites.size() 
              << ",\nthe number of different (x, y) pairs shown in data file is:"
              << (data_3d.site_x).size() << std::endl;
    std::cout << "The two number must be the same, please check the poly and "
              << "data files!\nProgram stopped with input error!\n";
    exit(0);
  }

  // check the order 
  for(unsigned int i = 0; i < sites.size(); i++){
    Point s = sites[i];
    double x = s(0);
    double y = s(1);
    double diff_x = x - data_3d.site_x[i];
    double diff_y = y - data_3d.site_y[i];
    // The two points are approximated to be the the same 
    // in allowable around 5m error range
    if( (std::abs(diff_x) > 5.0) || (std::abs(diff_y) > 5.0) ){
      std::cout << "The horizontal location of " << i + 1 << "-th sites shown"
                << " in poly file\nand data file are not the same!\n";
      std::cout << "The horizontal location and ranked order of the sites"
                << " shown in\npoly file and data file must be the same,"
                << " please check them!\n";
      std::cout << "Program stopped with input error!\n";
      exit(0);
    }
  }

  return;
}

void input_forward_parameters(std::string fwd_para_file, Fwd_Para_3D& parameter)
{
  // function for 3D forward modelling parameter related file reading
  // this function is slightly different from input_fwd_control_parameters
  // function in starting_model and region info reading.
  bool having_marker_9999999 = 0;
  para_output << "\n***********************Parameters from MT3D_fwd_para file"
              << "**********************\n";
  std::ifstream in_stream(fwd_para_file.c_str());
  if( !(in_stream.good()) ){
    std::cout << "Input error: " << fwd_para_file << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  
  // default
  parameter.algorithm = 0;
  in_stream >> parameter.n_f;
  para_output << "n_frequencies:\t"<<parameter.n_f <<"\n";
  if(parameter.n_f<1){
    std::cout<<"Error, please set n_f > 1!\n";
    para_output.close();
    std::abort();
  }
  parameter.f.resize(parameter.n_f);
  para_output << "frequencies:\t";
  for(unsigned int i=0; i<parameter.n_f; i++){
    in_stream >> parameter.f[i]; 
    para_output << parameter.f[i] <<"\t";
    if(parameter.f[i]<0) {
    std::cout<<"Error, please set f > 0!\n";
    para_output.close();
    std::abort();
    }
  }
  parameter.starting_model.resize(parameter.n_f);
  para_output <<"\n" << "starting_model:\t";
  for(unsigned int i=0; i<parameter.n_f; i++){
  in_stream >> parameter.starting_model[i];
  para_output << parameter.starting_model[i] <<"\t";
  }
  in_stream  >> parameter.starting_counter; 
  para_output <<"\n" << "starting_counter:\t"
            << parameter.starting_counter <<"\n";
  if(parameter.starting_counter < 1) {
    std::cout<<"Error, please set starting_counter = 1,2,..,\n";
    para_output.close();
    std::abort();
  }
  in_stream >> parameter.theta; 
  para_output << "theta:\t" << parameter.theta <<"\n";
  if(parameter.theta<0||parameter.theta>90) {
    std::cout<<"Error, please set 0 < theta < 90!\n";
    para_output.close();
    std::abort();
  }

  // read parameters for boundary conditions
  in_stream >> parameter.n_layer; 
  para_output << "n_layers:\t" << parameter.n_layer <<"\n";  
  if(parameter.n_layer<2) {
    std::cout<<"Error, please set n_layer >= 2!\n";
    para_output.close();
    std::abort();
  }
  parameter.EP.resize(parameter.n_layer);
  for(int i=0; i<parameter.n_layer; i++) {
    parameter.EP[i].resize(4);
    for(int j=0; j<4; j++){
      in_stream >> parameter.EP[i][j];
      para_output << parameter.EP[i][j] <<"\t";
    }
    para_output<<"\n";
  }

  in_stream >> parameter.n_regions;
  para_output << "n_regions:\t" << parameter.n_regions <<"\n"; 
  if(parameter.n_regions<1) {
    std::cout<<"Error, please set n_regions >= 1!\n";
    para_output.close();
    std::abort();
  }
  for(int i=0; i<parameter.n_regions; i++)  {
    int marker;
    std::vector<double> temp(3);
    in_stream >> marker;
    para_output << marker <<"\t";
    for(int j=0; j<3; j++) {
      in_stream >> temp[j];
      para_output << temp[j] <<"\t";
    }
    para_output << "\n";
    parameter.region_table[marker] = temp;
    if(marker == 9999999){
      having_marker_9999999 = 1;
      if( parameter.region_table[marker][0] > 1.e-10 )
      {
        std::cout << "Input error: the conductivity of air > 1.e-10, "
                  << "please reset it!\n";
        para_output.close();
        std::abort();
      }
    }
  }
  if(!having_marker_9999999){
    std::cout << "Error, the regional marker of air space should be 9999999,"
              << " please check the forward parameter file!\n";
    para_output.close();
    std::abort();
  }

  // read sites information
  in_stream >> parameter.site_file;
  para_output << "site_file:\t" << parameter.site_file <<"\n";  
  std::ifstream site_stream (parameter.site_file.c_str());
  if( !(site_stream.good()) )
  {
    std::cout << "Input error: " << parameter.site_file 
              << " file was not found!\n";
    para_output.close();
    std::abort();
  }
  site_stream >> parameter.n_sites; 
  site_stream >> parameter.marker;   
  para_output << "n_sites:\t" << parameter.n_sites << "\n";  
  para_output << "site_marker:\t" << parameter.marker << "\n";
  if(parameter.n_sites<0) {
    std::cout << "n_sites:\t"<<parameter.n_sites<<"\n";  
    std::cout<<"Error, please set n_sites>0 in file:\t"<<parameter.site_file;
    std::cout<<"!\n";
    para_output.close();
    std::abort();
  }
  parameter.u.resize(parameter.n_sites);
  parameter.v.resize(parameter.n_sites);
  parameter.mu.resize(parameter.n_sites);
  for(int i=0; i<parameter.n_sites; i++) {
    double x1 = 0., x2 = 0., x3 = 0.;
    site_stream >> x1 >> x2 >>x3; 
    parameter.u[i] = Point(x1,x2,x3);  // u
  }
  for(int i=0; i<parameter.n_sites; i++) {
    double x1 = 0., x2 = 0., x3 = 0.;
    site_stream >> x1 >> x2 >>x3; 
    parameter.v[i] = Point(x1,x2,x3);  // v
  }
  for(int i=0; i<parameter.n_sites; i++) {
    double temp = 0.;
    site_stream >> temp;
    if( !(std::abs(temp) > 0.) ){
      std::cout << "Error, negative or zero mu_r emerges!\n";
      para_output.close();
      std::abort();
    }
    parameter.mu[i] = temp;            // mu
  }
  // output u, v and mu_r belonging to the sites into file
  for(int i=0; i<parameter.n_sites; i++) {
    para_output << "u:\t" << parameter.u[i](0) << "\t"
                        << parameter.u[i](1) << "\t"
                        << parameter.u[i](2) << "\t";
    para_output << "v:\t" << parameter.v[i](0) << "\t"
                        << parameter.v[i](1) << "\t"
                        << parameter.v[i](2) << "\t";
    para_output << "mu:\t" << parameter.mu[i] << "\n";
  }

  return;
}
