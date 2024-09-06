/******************************************************************************
 *    Function       : Forwad modeling interface to produce synthetic data,   * 
                       compute forward data and calculate sensitivity and     *
                       gradient of data misfit for 1D, 2D and 3D MT problems  *
 *    Author         : Huang Chen and Zhengyong Ren                           *
 *    Copyright      : Huang Chen, 2019 and Zhengyong Ren, 2017               *
 *    Email          : chenhuang@cqu.edu.cn; renzhengyong@csu.edu.cn       *
 *    Created time   : 2019.07.22                                             *
 *    Last revision  : 2023.08.25                                              *
 ******************************************************************************/
#ifndef FWD_COMP_H
#define FWD_COMP_H

#include <iomanip>
#include <cstdlib>
#include "struct_defs.h"
#include "em.h"
// for random number generation
#include <chrono>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
// for parallelization
#include <omp.h>
#include "mkl_service.h"
// for Pardiso solver
#include <Eigen/PardisoSupport>
#include "mesh3d.h"
#include "fem.h"
#include "bc.h"
#include "post_process.h"
#include "node.h"
#include "timer.h"

class Fwd_Sens_Comp
{
  public:
  // overloaded constructors
  // !!! for Fwd_only (produce synthetic data)
  // for 1D MT problems;  
  Fwd_Sens_Comp(Fwd_Para_1D _fwd_para_1d, std::string _data_method);
  // for 2D MT problems
  Fwd_Sens_Comp(Fwd_Para_2D _fwd_para_2d, std::string _data_method);
  // for 3D MT problems
  Fwd_Sens_Comp(Fwd_Para_3D _fwd_para_3d, std::string _data_method);
  // !!! for inversion
  // for 1D MT problems;  
  Fwd_Sens_Comp(Fwd_Para_1D _fwd_para_1d, Data_1D _data_1d);
  // for 2D MT problems
  Fwd_Sens_Comp(Fwd_Para_2D _fwd_para_2d, Data_2D _data_2d);
  // for 3D MT problems
  Fwd_Sens_Comp(Fwd_Para_3D _fwd_para_3d, Data_3D _data_3d);

  // destructor
  ~Fwd_Sens_Comp();

  // !! for FWd_only
  // !!! for 3D MT forward modelling and inversion  
  // !! for FWd_only                         
  // compute and contaminate 3D MT forward impedance data
  void Produce_MT3d_Synthetic_Data(double rel_err, double abs_err, 
                                   double rel_tipper_err, double abs_tipper_err,
                                   std::string inv_approach);
  // !! for inversion
  // compute forward responses F_m (full impedance data) at model m
  void Comp_MT3d_Fwd_Responses(EM::VectorXD& F_m, const EM::VectorXD& m,
                  const std::string& inv_approach, const double abn[3]);
  // compute sensitivity of individual frequency by using sensitivity 
  // equation method
  void compute_sensitivity(FEM& fem, PostProcess& post_p, 
                           EM::MatrixXD& J, unsigned int i,
                 std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj);

  // ! for L-BFGS inversion
  // compute forward responses F[m] and the gradient of the data misfit term 
  // at model m
  void comp_MT3d_fwd_grad_phi_d(EM::VectorXD& F_m, EM::VectorXD& grad_phi_d,
                                const EM::VectorXD& m, const EM::VectorXD& d,
                                EM::EigenSparseMatrixXD& Wd,
                                const std::string& inv_approach, 
                   std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj,
                                const double abn[3]);
  // compute the gradient of the data misfit term of individual frequency  
  // at model m by using sensitivity equation method
  void compute_gradient_phi_d_resp(EM::VectorXC& grad_phi_d_resp,
                                   FEM& fem, PostProcess& post_p,
                                   const int& index_f, const EM::MatrixXC& F_m_resp,
                                   const EM::VectorXD& d, 
                                   EM::EigenSparseMatrixXD& Wd,
                          std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj);
  
  private:
  // !!! for common use
  // forward parameter struct variables
  // which will be used to control individual forward modelling
  Fwd_Para_1D            fwd_para_1d;
  Fwd_Para_2D            fwd_para_2d;  
  Fwd_Para_3D            fwd_para_3d;

  // !!! for Fwd_only (produce synthetic data)
  // method of the data, i.e., MT1D, MT2D or MT3D
  std::string            data_method;   

  // !!! for inversion (including forward modelling)
  // data struct variables, which will be used to offer the info about 
  // the number of sites, frequencies, data components, and so on
  // and will Only be used for sensitivity calculating
  // just read, no writing!
  Data_1D                data_1d;
  Data_2D                data_2d;
  Data_3D                data_3d;
};

inline Fwd_Sens_Comp::~Fwd_Sens_Comp()
{

}
#endif // FWD_COMP_H
