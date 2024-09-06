/******************************************************************************
 *    Function       : Class Inv_Basis contains some basic routines, which    *
 *                   : are used in a variety of inverse schemes               *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2019                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2019.07.23                                             *
 *    Last revision  : 2023.06.30                                             *
 ******************************************************************************/
#ifndef INV_BASIS_H
#define INV_BASIS_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <stdexcept>
#include <Eigen/Core>
#include "struct_defs.h"
#include "fwd_sens_comp.h"
#include "timer.h"
#include "em.h"

class Inv_Basis
{
  public:
  // default constructor used for object defining without initialization
  Inv_Basis();
  // other reloaded constructors
  Inv_Basis(Startup _startup, Data_1D _data_1d, Init_Prior_Mod_1D _init_mod_1d,
            Init_Prior_Mod_1D _prior_mod_1d, Inv_Para _inv_para);
  Inv_Basis(Startup _startup, Data_2D _data_2d, Init_Prior_Mod_2D _init_prior_mod_2d,
            Inv_Para _inv_para, Fwd_Para_2D _fwd_para_2d);
  Inv_Basis(Startup _startup, Data_3D _data_3d, Init_Prior_Mod_3D _init_prior_mod_3d,
            Inv_Para _inv_para, Fwd_Para_3D& _fwd_para_3d, Mesh3D& _MT3D_inv_mesh);
  // destructor
  ~Inv_Basis();

  /////////////////////////////////////////////////////////////////////////////
  //                               PLEASE NOTE:                              //
  // THE PRESENT MODEL mk USED BY THE FOLLOWING FUNCTIONS IS A DATA MEMBER,  //
  // WHICH WILL BE UPDATED WITH THE INVERSION ITERATIONS!                    //
  /////////////////////////////////////////////////////////////////////////////

  // !!! commonly used for L-BFGS algorithm
  // initialize the inverse of model covariance matrix for MT2D (later) and MT3D
  void init_inv_mod_covar_matrix(EM::EigenSparseMatrixXD& Cm_inv);
  // only compute forward data F_mk of present model mk
  void comp_F_mk();
  void compute_misfit_RMS_roughness(const double& phi, const EM::VectorXD& m,
                                    const double& lambda);
  // output the predicted data and its residual compared with the input data
  // formal parameters n_adjust and k are used for file name generating
  void output_data_and_residual(unsigned int n_adjust, unsigned int k);                                
  // output the inversion results (inverted resistivity model)
  // formal parameters n_adjust and k are used for file name generating
  void Output_Inv_Mod(unsigned int n_adjust, unsigned int k);
  // !! for MT1D (not present), MT2D (not present) and MT3D
  // compute the forward responses at model m (arbitrary)
  void comp_F_m(const EM::VectorXD& m, EM::VectorXD& F_m);

  // !! for MT1D, MT2D and MT3D
  // compute the value of the objective function \phi and
  // the gradient of \phi by using sensitivity equation method at model m
  void comp_phi_g_m_sen_eq(const EM::VectorXD& m, const double& lambda,
                           double& phi, EM::VectorXD& grad_phi);
  // do inexact line search algorithm with strong Wolfe condition and the to 
  // update the value of objective function and it's gradient at the updated x
  void line_search_strong_Wolfe(double& fx, EM::VectorXD& x, EM::VectorXD& grad,
                                double& step, const EM::VectorXD& drt, 
                                const EM::VectorXD& xp, const Line_Search_Para& param,
                                unsigned int& iter);

  // !!!Only used for adaptive inversion 
  // compute the gradient of m
  void Comp_Abs_Gra_M(EM::VectorXD& abs_gra_m);
  // compute the abs(gradient) of \phi_d of the last updated model 
  void Comp_Abs_Gra_Phi_D(EM::VectorXD& abs_gra_phi_d);
  // return the last updated model based on the present mesh
  void get_m_last(EM::VectorXD& m_ini);


  protected:
//=============================== commonly used =============================
  // !!!used for 'data' inputting  
  // used for getting the project name, data method  and inversion algorithm
  Startup                                                     startup;
  // struct variable of inversion parameter 
  Inv_Para                                                    inv_para;
  // abn[0], ~[1], ~[2] represents a, b and n 
  // variables, respectively, only used in inv_basis.cpp
  double                                                      abn[3];
  // the tolerance of the max times of inversion iteration
  unsigned int                                                tol_times;
  // !!!used for 'data' processing
  // trade-off parameter, used for each inversion iteration
  double                                                      lambda;
  // the number of the model parameters
  unsigned int                                                M; 
  // the number of the data 
  unsigned int                                                N;
  EM::VectorXD                                                d;
  // data weighted matrix (CRS stored sparse matrix)
  EM::EigenSparseMatrixXD                                     Wd;
  // inverse of model covariance matrix, 
  // which is only used for MT2D and MT3D                            
  EM::EigenSparseMatrixXD                                     Cm_inv;
  EM::VectorXD                                                m_ref;
  // at the beginning, m_k = m_init (initial model)
  EM::VectorXD                                                m_k;
  EM::VectorXD                                                m_kp1;
  // forward data vector of present model mk F_mk = F(mk) 
  // F_mk has form of [(rho_a0,...,rho_aN), (phi0,...,phiN)]
  // at the beginning, F_mk = F(m0)
  EM::VectorXD                                                F_mk;

  // !!!used for "data" outputting
  // used for stopping inversion iteration and / or outputting
  std::vector<double>                                         data_misfit;
  std::vector<double>                                         RMS;
  std::vector<double>                                         roughness;
  // for storing the history of lambda, used for terminating  
  // the L-BFGS iteration on the present inversion mesh
  std::vector<double>                                         lambda_history;

  //=============================only for MT3D ================================
  // !!!used for 'data' inputting                                                      
  Data_3D                                                     data_3d;
  // for initial and prior model
  Init_Prior_Mod_3D                                           init_prior_mod_3d;
  // !!!used for 'data' inputting and processing
  // A pointer pointed to the global variable fwd_para_3d
  // defined in CRT3DMT.cpp
  Fwd_Para_3D*                                                fwd_para_3d;
  // inversion mesh of MT3D
  Mesh3D*                                                     MT3D_inv_mesh;  
  // map variable, mapping the tet id in inversion-domain to 
  // the id set of its surrounding tetrahedra, which will be used for calculating  
  // Cm_inv and the gradient-based inversion mesh refinement indicator
  std::map< unsigned int, std::set< unsigned int > >          tet_id_to_tet_indexs;                                                   
};

inline Inv_Basis::Inv_Basis()
{

}

inline Inv_Basis::~Inv_Basis()
{
  //std::cout << "Hello Inv_Basis" << std::endl;
}

#endif //INV_BASIS_H
