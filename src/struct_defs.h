/******************************************************************************
 *    Function       : Definition of struct types for reading input files     *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2019                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2019.07.22                                             *
 *    Last revision  : 2023.09.16                                             *
 ******************************************************************************/
#include "em.h"
#include "point.h"
#include <set>

#ifndef STRUCT_DEFS_H
#define STRUCT_DEFS_H

//=======================for parameters in startup file========================
// struct for .startup file
struct Startup 
{  
  // name of the project
  std::string proj_name;

  // !!!parameters in data section part
  // method of the data
  std::string data_method;
  // name of the data file
  std::string data_file;

  // !!!file_names in model section part
  // only be initialized and used for MT1D problem
  std::string initial_model_file;
  // only be initialized and used for MT1D problem
  std::string prior_model_file;
  // only be initialized and used for MT2D and MT3D problem
  // including the name of the mesh and the respective physical properties
  std::string ini_pri_model_file;

  // !!!parameters in inversion method section part
  // name of the parameter file of inversion
  std::string inv_control_file;
  // file name for 3D MT forward modelling control parameters
  std::string fwd_control_file;
  // function variable for forward modeling or L-BFGS inversion
  std::string inv_algorithm;

  // !!!forward solver related parameters in forward section part
  // name of the parameter file used for forward modelling
  std::string fwd_para_file;
  // !!the following variables are only useful 
  //   for Fwd_only  
  // relative error usual in range of (0, 0.1) for 
  // contaminating the synthetic data except the tipper data
  double rel_err;
  // absolute error (the value depends on the amplitude of the data) 
  // for contaminating the synthetic data except the tipper data
  double abs_err;

  // !the following variables are only used for MT3D
  // relative tipper error (usually set to 0)
  // for contaminating the synthetic tipper data of MT3D
  double rel_tipper_err;
  // absolute tipper error 
  // for contaminating the synthetic tipper data of MT3D
  double abs_tipper_err;
};

//========================for parametres in data file==========================
// for 1D MT problem
struct Data_1D
{

};

// for 2D MT problem
struct Data_2D
{

};

// for 3D MT problem
struct Data_3D
{
  std::string               data_type;
  std::string               time_h_f;
  std::string               unit;
  // geographic orientation of all data components relative to geographic north
  double                    orient;
  // geographic_origin, latitude and lontitude for the (0,0,0) 
  // point in the data file
  double                    lat_origin;
  double                    lon_origin;
  // number of observing frequencies
  unsigned int               n_f;
  // number of observing sites 
  unsigned int               n_sites;
  // The following 7 vectors are used for outputting the predicted
  // and residual data and write the data file (optional)
  std::vector<std::string>   site_code, data_comp;
  std::vector<double>        lat, lon, x, y, z;

  // used for checking data
  std::vector<double>        site_x;
  std::vector<double>        site_y;

  // the following two variables are initialized by assignment
  // the number of components: if the data_type is Full_impedance,
  // it will be set to be 4 (ZXX, ZXY, ZYX, ZYY), while if that is 
  // Full_impedance_Tipper, it will be set to be 6
  // (ZXX, ZXY, ZYX, ZYY, TX, TY)
  unsigned int               n_comps;
  // N = n_f * n_sites * n_comps * 2 (including both real and imaginary parts)
  unsigned int               N;
  // used for checking the input number of frequencies
  std::set<double>           periods;
  // used to initialize fwd_para_3d.f
  // the order of the frequencies in fwd_para_3d.f vector should be the same as
  // that in data file
  std::vector<double>        f;
  EM::VectorXD               d;
  // stanard deviation of the data
  EM::VectorXD               std_dev;            
};


//==========================for inversion parameters===========================
struct Inv_Para
{
  // !!!parameters for adaptive inversion
  // adap_method = 0: performing standard inversion, others: adaptive inversion
  // = 1: with mesh adjustment strategy of gradient \phi_d and gradient of m
  unsigned int adap_method;
  // tol of maximum number of the inversion unknowns. Once the number of the
  // inversion unknowns of the refined mesh reaches to this number, the 
  // refinement will stop and the optimization will be end with this mesh
  unsigned int max_No_inv_unknowns;
  // tol of maximum adjustment times of inversion mesh
  unsigned int max_adjust_times;
  // tol of maximum inversion iteration times per mesh
  unsigned int max_iter_times;
  // tol of fraction of posterior error grad_m for refining inversion mesh 
  // frac_m_r is also used to multiply by the max resistivity in each pair for 1D MT
  double frac_m_r;
  // fraction of other posterior errors for refining inversion mesh  
  double frac_others_r;

  // !!!parameters for standard inversion
  // lower bound of all free parameters (conductivities), 
  // i.e., a in Kim et al., 2011
  double a;
  // upper bound of all free parameters (conductivities)
  // i.e., b in Kim et al., 2011
  double b;
  // positive constant for parameter trasformation, 
  // i.e., n in Kim et al., 2011
  double n;
  // trade-off parameter
  double lambda0;
  // tolerance of the reduction fraction of RMS between two iterations
  // for cooling lambda
  double tol_rms_reduc_frac;
  double cool_factor;
  // !! variables used for terminating the inversion
  // tolerance of the times of continuously cooling lambda
  double tol_times_cont_cool_lambda;
  // tolerance of lambda
  double tol_lambda; 
  // maximum times of inversion iterations
  int tol_iter_times;
  // tolerance of rms
  double tol_rms;

};

//================for parametres in initial or prior model file================
struct Init_Prior_Mod_1D
{
  // used for guaranteeing the model is compatible with the data method
  std::string  data_method;
  // the number of layers;
  unsigned int n_layer;
  // !!!parameters of individual layer, including 
  // conductivity in linear domain
  EM::VectorXD sigma;
  // relative permittivity
  EM::VectorXD epsilon_r;
  // relative permeability
  EM::VectorXD mu_r;
  // depth of the bottom of each layer to ground
  EM::VectorXD depth;
};

struct Init_Prior_Mod_2D
{

};

struct Init_Prior_Mod_3D
{
  // !!!for inversion mesh, which will be changed during adaptive inversion
  // mesh name = model_name.counter.*
  // * are ele, face, neigh, poly, node etc. 
  std::string model_name;  
  // counter for the inversion mesh   
  unsigned int counter;  
  // the number of the regions
  unsigned int n_regions;  

  // !!!initial model related 
  // marker to conductivity (sigma),    
  // relative_mu (mu_r), relative_epsilon (epsilon_r) 
  // used for initializing vector m_ini and assigning 
  // the electric parameter of the tetrahedral elements     
  std::map< int, std::vector< double > > region_table_ini;   
  // for storing the id of the nodes in the inversion domain
  // which will be used for calculating of the inverse of Cm and the
  // model gradient-based inversion mesh refinement indicator 
  std::set< unsigned int > sub_node_id;
  // m_ini only stores the conductivities of the tet elements in inversion domain
  // m_ini = m0 
  std::vector< double > m_ini; 
  // m_ini_id[i] storing the element id of the parameter m_ini[i]
  // which will be used for calculating the inverse of Cm and the
  // model gradient-based inversion mesh refinement indicator
  std::vector< unsigned int > m_ini_id;
  // fwd_tet_id_in_mj[j], j = 0 ... M-1 storing the ids of the forward elements  
  // nested in j-th parameter element, which will be used for calculating of  
  // sensitivity　and the gradient of the misfit term
  std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj;

  // !!!priori model related   
  // marker to conductivity (sigma),    
  // relative_mu (mu_r), relative_epsilon (epsilon_r)   
  // only used for initialization of vector m_ref     
  std::map<int, std::vector< double > > region_table_ref; 
  // m_ref only stores the conductivities of the tet elements in inversion domain
  std::vector< double > m_ref;   
};

//=======================for forward modeling parameters=======================
//   the forward modeling struct variable contains the infomation of the     //
//   geoelectrical model (geometry and physical property) and the source     //
//   for 1D-3D, and the control parameter of forward algorithms for 2D & 3D .//
//===========================================================================//
// for 1D MT problem, used for both producing synthetic data and inversion 
// (but respective ways of assignment are different, for inversion, n_f and f 
//  are initialized by using Data_1D variable and others are initialized by 
//  Init_Prior_Mod_1D variable
struct Fwd_Para_1D 
{
  // used for guaranteeing the fwd parameter is compatible with the data method
  std::string data_method;
  // the number of layers
  unsigned int n_layer;
  // parameters of individual layer, including
  // conductivity, relative permittivity, relative permeability
  // and depth from ground to the bottom of the layer
  EM::MatrixXD EP;
  // the number of the observing frequencies (in Hz)
  unsigned int n_f;
  // observing frequencies
  EM::VectorXD f;
};

struct Fwd_Para_2D
{

};

// for 3D MT problem, used for both synthetic data producing and data inversion 
// (but respective ways of assignment are different, for inversion, n_f and f 
//  are initialized by Data_3D variable, n_regions and region_table are  
//  initialized by Init_Prior_Mod_3D variable and others are intialized by .fwd file)
struct Fwd_Para_3D 
{
  // !!! 1. algorithm. 
  // 0---standard fem; 
  int algorithm;

  // !!! 2. source information
  // the number of frequencies
  unsigned int n_f;  
  // frequencis (Hz)
  std::vector< double > f; 
  // incident angle of plane waves                            
  double theta;       

  // !!! 3. for forward mesh, 
  // when using Openmp with mesh refinement, the size of the 'starting_model' 
  // vector should be set to the number of observing frequencies
  // mesh name = starting_model.starting_counter.*
  // * are ele, face, neigh, poly, node etc. 
  std::vector< std::string > starting_model;  
  // counter for the initial mesh for all frequencies    
  unsigned int starting_counter;    

  // !!! 4. 1D boundary conditions
  // numbers of 1D layers
  unsigned int n_layer;
  // sigma mu_r epsilon_r  d   
  // conductivity(sigma), relative mu (mu_r),
  // (d) depth is from air-earth
  // relative epsilon (epsilon_r)
  // interface(z=0) to bottom 
  // layer in each layer                 
  std::vector< std::vector< double > > EP;  

  // !!! 5. map of regional marker and corresponding \sigma, \mu_r, \epsilon_r
  //        or conductivity vector for element physical property assignment
  // numbers of regions
  unsigned int n_regions;     
  // map to conductivity (sigma),    
  // relative_mu (mu_r), relative_epsilon (epsilon_r)        
  std::map<int, std::vector< double > > region_table;      
                    
  // !!! 6. sites information
  // file for sites
  std::string site_file;     
  // numbers of sites            
  int n_sites;     
  // site marker                         
  int marker;     
  // Eu, Hu, Ev, Hv are calculated                                               
  std::vector< Point > u, v;        
  // relative_mu at sites for calculating impedance        
  std::vector< double > mu;                 
};

// Parameters to contro the line search procedure
struct Line_Search_Para
{
  /*****************Parameters for line searching of step size****************/
  //
  // The maximum number of trials for the line search.
  // This parameter controls the number of function and gradients evaluations
  // per iteration for the line search routine. The default value is \c 20.
  //
  int    max_linesearch;
  //
  // The minimum step length allowed in the line search.
  // The default value is \c 1e-20. Usually this value does not need to be
  // modified. This member is only used for Wolfe condition.
  //
  double min_step;
  //
  // The maximum step length allowed in the line search.
  // The default value is \c 1e+20. Usually this value does not need to be
  // modified.   // modified. This member is only used for Wolfe condition.
  //
  double max_step;
  //
  // A parameter to control the accuracy of the line search routine.
  // The default value is \c 1e-4. This parameter should be greater
  // than zero and smaller than \c 0.5.
  // ftol = c1 in book "Numerical Optimization"
  //
  double ftol;
  //
  // The coefficient for the Wolfe condition.
  // This parameter is valid only when the line-search
  // algorithm is used with the Wolfe condition.
  // The default value is \c 0.9. This parameter should be greater
  // the \ref ftol parameter and smaller than \c 1.0.
  // wolfe = c2 in book "Numerical Optimization"
  //
  double wolfe;
};

#endif // STRUCT_DEFS_H
