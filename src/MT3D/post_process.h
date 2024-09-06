/*****************************************************************************/
/*                                                                           */
/*  Copyright 2021                                                           */
/*  Huang Chen and Zhengyong Ren                                             */
/*  chenhuang@cqu.edu.cn and renzhengyong@csu.edu.cn                         */
/*                                                                           */
/*****************************************************************************/

#ifndef  _POST_PROCESS_H
#define  _POST_PROCESS_H

#include <vector>
#include <ostream>
#include <fstream>
#include <iomanip>
#include "point.h"
#include "node.h"
#include "em.h"

class FEM;   // E-formula
class PostProcess
{
 public: 
 // constructor for Fwd_only approach
  PostProcess(FEM& fem, int marker,
              std::vector<Point>& u,
              std::vector<Point>& v,
              std::vector<Real> & mu,
              std::ofstream& out_file);

 // overloaded constructor for other algorithms (including LBFGS, ...)
  PostProcess(FEM& fem, int marker,
              std::vector<Real> & mu,
              std::string& data_type,
              bool calculate_L1_L2_or_not);
  // destructor
  ~PostProcess() {}

  // compute Z, App.Resis., Phases, E, H, and output 
  // which will be used for producing synthetic data
  void post_processing(FEM& fem, int marker, 
                       std::vector<Point>& u,
                       std::vector<Point>& v,
                       std::vector<Real> & mu,
                       std::ofstream& OUT);

  // compute Z, App.Resis., Phases, E, H, and no output 
  // which will be used when doing inversion
  // the bool variable is used for judging whether to calculate the global
  // interpolation operator or not (when we call comp_F_m or comp_F_mk function
  // , we don't need to calculate the global interpolation operator L_1,2_Z_XXX)
  void post_processing_without_output(FEM& fem, int marker, 
                                      std::vector<Real> & mu,
                                      std::string& data_type,
                                      bool calculate_L1_L2_or_not);

  // used for producing synthetic data and forward responses when doing inversion
  void Return_Sites(std::vector< Node >& _sites);
  void Return_Impedances(EM::MatrixXC& _Z);
  void Return_Impedances_plus_Tippers(EM::MatrixXC& _ZpT);
  void Return_Diagonal_Impedances(EM::MatrixXC& _diagonal_Z);
  void Return_Diagonal_Impedances_plus_Tippers(EM::MatrixXC& _diagonal_ZpT);
  
  // only used for producing synthetic data
  void Return_Impedances_Tippers(EM::MatrixXC& _Z, EM::MatrixXC& _T);
  // used for sensitivity calculation
  void get_n_sites(unsigned int& n_sites);

  // !!!used for calculation of sensitivity and gradient
  // discrete global interpolation operator for
  // full impedance components of all sites
  EM::EigenSparseMatrixXC L_1_Z, L_2_Z;
  // off-diagonal impedance components of all sites
  EM::EigenSparseMatrixXC L_1_diagonal_Z, L_2_diagonal_Z;
  // full impedance and tipper components of all sites
  EM::EigenSparseMatrixXC L_1_ZpT, L_2_ZpT;
  // off-diagonal impedance and tipper components of all sites
  EM::EigenSparseMatrixXC L_1_diagonal_ZpT, L_2_diagonal_ZpT;

  private:
  // !!!used for compute forward responses
  // for returning sites and impedances
  std::vector< Node >                           sites_;
  // the number of sites
  int                                           n_sites_;
  // for storing full impedance data
  EM::MatrixXC                                  Z_;
  // for storing off-diagonal impedance data
  EM::MatrixXC                                  diagonal_Z_;
  // for storing full tipper data
  EM::MatrixXC                                  T_;
  // observing frequencies
  double                                        f_;
};

#endif // POST_PROCESS
