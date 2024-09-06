/******************************************************************************
 *    Function       : Definition of Adaptive class                           *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2020                                       *
 *    Email          : chenhuang@cqu.edu.cn                                *
 *    Created time   : 2020.07.07                                             *
 *    Last revision  : 2023.10.12                                             *
 ******************************************************************************/

#ifndef ADAPTIVE_H
#define ADAPTIVE_H
#include <iostream>
#include <fstream>
#include "struct_defs.h"
#include "lbfgs_inv.h"
// used for calling the sleep function to suspend the execution temporarily
#include <unistd.h>
#include <cstdlib>

class Adaptive
{
  public:
  // constructor used for object defining without initialization
  Adaptive();
  // other reloaded constructors, _param is only used for L-BFGS inversion
  Adaptive(Startup _startup, Data_1D _data_1d, Init_Prior_Mod_1D _init_mod_1d,
           Init_Prior_Mod_1D _prior_mod_1d, Inv_Para _inv_para);
  Adaptive(Startup _startup, Data_2D _data_2d, Init_Prior_Mod_2D _init_prior_mod_2d,
           Inv_Para _inv_para, Fwd_Para_2D& _fwd_para_2d);
  Adaptive(Startup _startup, Data_3D _data_3d, Init_Prior_Mod_3D _init_prior_mod_3d,
           Inv_Para _inv_para, Fwd_Para_3D& _fwd_para_3d, Mesh3D& _MT3D_inv_mesh);

  // destructor
  ~Adaptive();

  // template function for sorting error vector
  // used for doing mesh refinement
  template <typename T> 
  std::vector<size_t>  sort_indexes(const std::vector<T>& v);

  // perform adaptive inversion
  void Do_Adapt_Inv();
  // adjust the inversion mesh and update the corresponding members 
  // of this class
  void Adjust_Inv_Mesh();
  // get the number of dofs of forward modelling on the present fwd mesh for MT3D
  unsigned int Get_n_dofs(Mesh3D& mesh3d);


  private:
  //=============================== commonly used ================================
  // !!!for "data inputting"
  // used for determining the method of the data and the method of adaptive
  Startup                    startup;
  // struct variable of inversion parameter 
  Inv_Para                   inv_para;
  // used for performing L-BFGS inversion
  LBFGS_Inv                  LBFGS; 

  // !!!for "data processing"
  // for the gradient of m based inversion mesh refinement strategy
  EM::VectorXD               abs_gra_m;
  // for the gradient of \phi_d based inversion mesh refinement strategy 
  // the abs(gradient) of \phi_d of the last updated model on the present mesh
  EM::VectorXD               abs_gra_phi_d; 
  // for the sensitivity based inversion mesh refinement strategy 
  std::vector<double>        weig_sens;
  // for the model resolution matrix based inversion mesh refinement strategy
  std::vector<double>        diag_resol;
  // data-driven inversion mesh adjustment indicator
  // adjust_indicator == 'abs_gra_phi_d', weig_sens, and diag_resol 
  // corresponding to adap_method == 1, ...
  std::vector<double>        adjust_indicator;
  // the number of the inversion unknowns
  unsigned int               M;
  // counter of the times of adjustment of inversion mesh;
  unsigned int               k;
  // counter of the total number of the inversion iteration
  unsigned int              total_iter_number;
  // used for L-BFGS iteration on the last mesh
  unsigned int              left_iter_number;
  // RMS of the last L-BFGS iteration, used for stopping inversion 
  // mesh adjustment 
  double                     RMS_last;
  // the lambda of the last L-BFGS iteration
  // used for updating the initial lambda for changed inversion mesh
  double                     lambda_last;
  // used for storing the inversion result on the previous mesh
  EM::VectorXD               m_inv; 
  // the number of dofs of forward modelling on the present fwd mesh for MT3D
  unsigned int               n_dofs;

  // used for storing the parameters of originally initial and priori models
  // on previous mesh. Here, epsilon_r_ini = epsilon_r_ref = epsilon_r; 
  // mu_r_ini = mu_r_ref = mu_r, which keep constant during the inversion
  EM::VectorXD m_ini, m_ref, epsilon_r, mu_r;
  // used for storing the parameters of originally initial and priori models
  // on the updated mesh
  std::vector<double> m_ini_updated, m_ref_updated, epsilon_r_updated, mu_r_updated;
  // used for sotring the conductivity inversion result on the updated mesh
  std::vector<double> m_inv_updated;

  //=============================only for MT3D ================================
  Data_3D                    data_3d; 
  Init_Prior_Mod_3D          init_prior_mod_3d;
  Fwd_Para_3D*               fwd_para_3d;      
  std::string                init_fwd_mod_name;    
};

inline Adaptive::Adaptive()
{

}

inline Adaptive::~Adaptive()
{

}

#endif ADAPTIVE_H
