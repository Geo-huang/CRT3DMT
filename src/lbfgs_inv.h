/******************************************************************************
 *    Function       : class LBFGS_Inv is used to perform inversion by using  *
 *                   : L-BFGS algorithm in model space. The algorithm are     *
 *                   : from Algorithm 9.1 and 9.2 showm in book Nocedal, J.,  *
 *                   : & Wright, S. (2006). Numerical optimization.           *
 *    Author         : Huang Chen,                                            *
 *    Copyright      : Huang Chen, 2020                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2020.11.23                                             *
 *    Last revision  : 2023.09.30                                             *
 ******************************************************************************/
#ifndef LBFGS_INV_H
#define LBFGS_INV_H
#include <iostream>
#include <Eigen/Core>
#include <stdexcept> 
#include <assert.h>
#include "struct_defs.h"
#include "inv_basis.h"

class LBFGS_Inv: public Inv_Basis
{
  public:
  // the default constructor used for defining object without initialization
  LBFGS_Inv();
  // other reloaded constructors
  LBFGS_Inv(Startup _startup, Data_1D _data_1d, Init_Prior_Mod_1D _init_mod_1d,
            Init_Prior_Mod_1D _prior_mod_1d, Inv_Para _inv_para, 
            unsigned int _n_adjust = 0);
  LBFGS_Inv(Startup _startup, Data_2D _data_2d, Init_Prior_Mod_2D _init_prior_mod_2d,
            Inv_Para _inv_para, Fwd_Para_2D& _fwd_para_2d, unsigned int _n_adjust = 0);
  LBFGS_Inv(Startup _startup, Data_3D _data_3d, Init_Prior_Mod_3D _init_prior_mod_3d,
            Inv_Para _inv_para, Fwd_Para_3D& _fwd_para_3d, Mesh3D& _MT3D_inv_mesh,
            unsigned int _n_adjust = 0);

  // default destructor
  ~LBFGS_Inv();

  // checke the validity of line search parameters.
  void check_param() const;
  // define the size of vector or matrix
  void define_size();
  // run L-BFGS inversion
  void Do_LBFGS_Inv();

  // add correction vectors to the BFGS matrix
  void add_correction(const EM::RefConstVec& s, const EM::RefConstVec& y,
                      const EM::RefConstVec& d);
  // calculate a * Hk * v, which then stores in res vector
  void apply_Hv(const EM::VectorXD& v, const double& a, EM::VectorXD& res);
  // correct all the recent y vectors and the coorsponding variables
  // ys and theta by using the updated lambda
  void correct_y_ys_theta(double lambda, double lambda_p);

  // !!! used for adaptive inversion
  // get the RMS of the last L-BFGS iteration,
  // which is possibly less than tolerance of RMS
  void get_RMS_last(double& RMS_last);
  void get_lambda_last(double& lambda_last);
  void get_iter_number(unsigned int& n_iter);
  // update the initial lambda, i.e., change the value of  
  // the inherited variable Inv_Basis::lambda
  void update_initial_lambda(unsigned int k);
  void update_initial_lambda(double lambda0);
  void update_LBFGS_iter_times(unsigned int tol_LBFGS_times);

  // !!! used for testing calculation of the gradient of the \phi_d wsp model m
  void test_grad_phi_d();
  
  private:
  // derived private member variables
  // the times of inversion mesh adjustment
  int                                n_adjust;
  // counter of inversion iteration on the present mesh;
  unsigned int                       k;
  // the following three variables are used to correct the
  // stored y vectors when lambda was changed
  // diff_grad_phi_m = gradp_phi_m_{k} - gradp_phi_m_{k-1}, where
  // gradp_phi_m_{k} denotes the gradient of the model roughness of m_{k}
  EM::VectorXD                       diff_grad_phi_m;    
  // history of the diff_grad_phi_m vectors, similar to m_s and m_y
  EM::MatrixXD                       m_d;                

  // !!! variables usd in Do_LBFGS_Inv() function
  // struct parameters to control the line search procedure
  Line_Search_Para                   m_param;   
  // old x  
  EM::VectorXD                       m_xp;   
  // new gradient   
  EM::VectorXD                       m_grad; 
  // old gradient  
  EM::VectorXD                       m_gradp;  
  // moving direction, p_k
  EM::VectorXD                       m_drt;    

  // !! variables used in add_correction() and apply_Hv() functions
  // The number of corrections to approximate the inverse Hessian matrix.
   // The L-BFGS routine stores the computation results of previous \ref m
  // iterations to approximate the inverse Hessian matrix of the current
  // iteration. This parameter controls the size of the limited memories
  // (corrections). The default value is \c 6. Values m < 3 or m > 20 are
  // not recommended. Large values will result in excessive computing time.
  int                                m_m;      
  // 1./m_theta * I is the initial approximation to the Hessian matrix
  double                             m_theta; 
  // history of the s vectors 
  EM::MatrixXD                       m_s;      
  // history of the y vectors
  EM::MatrixXD                       m_y;   
  // history of the y's values, m_ys = 1.0 / \rho 
  EM::VectorXD                       m_ys;     
  // temporary values used in computing H * v, where H and v represent the
  // approximate inverse of Hessian matrix and a column vector, respectively
  EM::VectorXD                       m_alpha;  
  // number of correction vectors in the history, m_ncorr <= m_param.m
  int                                m_ncorr;  
  // a Pointer to locate the most recent history, 1 <= m_ptr <= m_param.m
  // details: s, y vectors and ys value are stored in cyclic order.
  // for example, if the current s-vector is stored in m_s[, m-1],
  // then in the next iteration m_s[, 0] will be overwritten.
  // m_s[, m_ptr-1] points to the most recent history (i.e., k-1 in book),
  // and m_s[, m_ptr % m] points to the most distant one (i.e., k-m).
  int                                m_ptr;   
};

inline LBFGS_Inv::LBFGS_Inv()
{

}

inline void LBFGS_Inv::check_param() const
{
  // Checking the validity of line search parameters.
  // An `std::invalid_argument` exception will be thrown if some parameter
  // is invalid.
  if(this->m_param.max_linesearch <= 0)
      throw std::invalid_argument("'max_linesearch' must be positive");
  if(this->m_param.min_step < 0)
      throw std::invalid_argument("'min_step' must be positive");
  if(this->m_param.max_step < this->m_param.min_step )
      throw std::invalid_argument("'max_step' must be greater than 'min_step'");
  if(this->m_param.ftol <= 0 || this->m_param.ftol >= 0.5)
      throw std::invalid_argument("'ftol' must satisfy 0 < ftol < 0.5");
  if(this->m_param.wolfe <= this->m_param.ftol || this->m_param.wolfe >= 1)
      throw std::invalid_argument("'wolfe' must satisfy ftol < wolfe < 1");
}

inline LBFGS_Inv::~LBFGS_Inv()
{

}

#endif // LBFGS_INV_H
