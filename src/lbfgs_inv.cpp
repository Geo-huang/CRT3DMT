/******************************************************************************
 *    Function       : class LBFGS_Inv is used to perform inversion by using  *
 *                   : L-BFGS algorithm in model space. The algorithm are     *
 *                   : from Algorithm 9.1 and 9.2 showm in book Nocedal, J.,  *
 *                   : & Wright, S. (1999). Numerical optimization.           *
 *    Author         : Huang Chen,                                            *
 *    Copyright      : Huang Chen, 2020                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2020.11.23                                             *
 *    Last revision  : 2023.09.30                                             *
 ******************************************************************************/
 #include "lbfgs_inv.h"
 // used for generating screen output log
 extern std::ofstream screen_output;

LBFGS_Inv::LBFGS_Inv(Startup _startup, Data_3D _data_3d, 
                     Init_Prior_Mod_3D _init_prior_mod_3d,
                     Inv_Para _inv_para, Fwd_Para_3D& _fwd_para_3d,
                     Mesh3D& _MT3D_inv_mesh, unsigned int _n_adjust):
               Inv_Basis(_startup, _data_3d, _init_prior_mod_3d,
                         _inv_para, _fwd_para_3d, _MT3D_inv_mesh)
{
  this->n_adjust               = _n_adjust;
  this->k                      = 0;

  // !!! initialize sturct variable m_param
  this->m_param.max_linesearch = 20;
  this->m_param.min_step       = 1e-20;
  this->m_param.max_step       = 1e+20;
  this->m_param.ftol           = 1e-4;
  this->m_param.wolfe          = 0.9;
  // checking the validity of line search parameters
  this->check_param();

  // The convergence rate increases with m but so with storage requirements
  // Hence, m will be problem and machine dependent
  // Nocedal and Wright (1999, p. 225) advocated m between 3 and 20
  // Large values will result in excessive computing time
  // default m = 6, which is also used in (Avdeev and Avdeeva, 2009)
  this->m_m                    = 6;
  this->m_theta                = 1.0;
  this->m_ncorr                = 0;
  // this makes sure that m_ptr % m == 0 in the first step
  this->m_ptr                  = this->m_m; 
  this->define_size();
}

void LBFGS_Inv::define_size()
{
  // !!! size definition of some vector and matrix variables
  // n: the number of the inversion unknowns;
  const int M = Inv_Basis::M;
  this->diff_grad_phi_m.resize(M);
  this->m_xp.resize(M);
  this->m_grad.resize(M);
  this->m_gradp.resize(M);
  this->m_drt.resize(M);

  // m: the number of corrections to approximate the inverse Hessian matrix
  const int m = this->m_m;
  this->m_d.resize(M, m);
  this->m_s.resize(M, m);
  this->m_y.resize(M, m);
  this->m_ys.resize(m);
  this->m_alpha.resize(m);
}

void LBFGS_Inv::Do_LBFGS_Inv()
{
  // x is a alternative name of the unknown earth model,
  // which will be iteratively solved
  EM::VectorXD& x = Inv_Basis::m_kp1; 
  // initialized by using initial guess
  // m_kp1 = m_k
  x = Inv_Basis::m_k;
  // the value of the objective function, fx = \phi
  double fx = 0.0;    
  // lambda is an alternative name of Inv_Basis::lambda
  double& lambda = Inv_Basis::lambda;
  Inv_Basis::lambda_history.push_back(lambda);
  // for judging whether lambda is continuously cooled by 
  // n (pre-defined) times or not
  double lambda_continuously_cooled_factor = 0.0;
  
  // for outputting RMS into file
  std::string RMS_file;
  // for outputting roughness into file
  std::string roughness_file;
  // for outputting lambda into file 
  std::string lambda_file;
  // for outputting step length into file
  std::string alfa_file;
  // for outputting search times into file
  std::string search_times_file;
  if(Inv_Basis::inv_para.adap_method == 0){
    RMS_file = Inv_Basis::startup.proj_name + ".RMS";
    roughness_file = Inv_Basis::startup.proj_name + ".roughness";
    lambda_file = Inv_Basis::startup.proj_name + ".lambda";
    alfa_file = Inv_Basis::startup.proj_name + ".alpha";
    search_times_file = Inv_Basis::startup.proj_name + ".search_times";
  }
  else
  {
    std::stringstream number;
    number << (this->n_adjust + 1);
    std::string n;
    n = number.str();
    RMS_file = Inv_Basis::startup.proj_name + std::string(".") + n + ".RMS";
    roughness_file = Inv_Basis::startup.proj_name + std::string(".") + n 
                     + ".roughness";
    lambda_file = Inv_Basis::startup.proj_name + std::string(".") + n 
                  + ".lambda";
    alfa_file = Inv_Basis::startup.proj_name + std::string(".") + n + ".alpha";
    search_times_file = Inv_Basis::startup.proj_name + std::string(".") + n
                        + ".search_times";
  }
  std::ofstream output_RMS(RMS_file.c_str());
  std::ofstream output_roughness(roughness_file.c_str());
  std::ofstream output_lambda( lambda_file.c_str() );
  std::ofstream output_alfa( alfa_file.c_str() );
  std::ofstream output_iter( search_times_file.c_str() );

  // calculate the gradient m_grad at m0
  if(Inv_Basis::startup.data_method == "MT1D"){
    std::cout << "L-BFGS inversion for MT1D is under construction" << std::endl;
    exit(1);
  }
  else if(Inv_Basis::startup.data_method == "MT2D"){
    std::cout << "L-BFGS inversion for MT2D is under construction" << std::endl;
    exit(1);
  }
  else{ // for MT3D
    Inv_Basis::comp_phi_g_m_sen_eq(x, lambda, fx, this->m_grad);
  } 
  Inv_Basis::compute_misfit_RMS_roughness(fx, x, lambda);

  std::cout << "RMS of the initial model is "
            << Inv_Basis::RMS[0]<< std::endl; 
  screen_output << "RMS of the initial model is "
                << Inv_Basis::RMS[0]<< std::endl; 
  if(Inv_Basis::startup.data_method == "MT1D"){
    std::cout << "L-BFGS inversion for MT1D is under construction" << std::endl;
    exit(1);
  }
  else if(Inv_Basis::startup.data_method == "MT2D"){
    std::cout << "L-BFGS inversion for MT2D is under construction" << std::endl;
    exit(1);
  }
  else // for MT3D
  {
    /*
    std::cout << "The initial resistivity model m_" << this->k 
            <<" in linear domain is:" << std::endl;
    for(unsigned int i = 0; i < Inv_Basis::m_k.size(); i++){
      double nomin = 0, denomin = 0., cond = 0.;
      nomin = Inv_Basis::inv_para.a + ( Inv_Basis::inv_para.b 
              * std::exp( Inv_Basis::inv_para.n * Inv_Basis::m_k(i) ) );
      denomin = 1.0 + std::exp( Inv_Basis::inv_para.n * Inv_Basis::m_k(i) );
      // conductivity in linear space
      cond = nomin / denomin;
      std::cout << 1.0 / cond << '\n';
    }
    */
    std::cout << '\n' << "Inversion iteration " << 1 << ":" << std::endl;
    screen_output << '\n' << "Inversion iteration " << 1 << ":" << std::endl;
  }
  // write the RMS of the initial model into the RMS file
  output_RMS << Inv_Basis::RMS[0] << std::endl;
  // write the roughness of the initial model into the roughness file
  output_roughness << Inv_Basis::roughness[0] << std::endl;
  // write the initial lambda corresponding to the initial model
  // into the lambda file, which is meaningless
  output_lambda << lambda << std::endl;
  // write the initial model into file
  Inv_Basis::Output_Inv_Mod(this->n_adjust, 0);
  Inv_Basis::output_data_and_residual(this->n_adjust, 0);

//////////////////////// !!! begin the core of L-BFGS /////////////////////////

  // use the normalized gradient direction as the search direction for the 
  // first iteration to avoid generating too rough model after the first 
  // iteration
  this->m_drt = -this->m_grad / this->m_grad.norm();
  
  // initial step length for the first L-BFGS iteration
  unsigned int iter = 0;
  double step = 1.0;
  double lambda_temp = lambda;

  while( (Inv_Basis::RMS[this->k] > Inv_Basis::inv_para.tol_rms) &&
         (this->k < Inv_Basis::tol_times) &&
         (lambda_temp > Inv_Basis::inv_para.tol_lambda) &&
         (std::abs(lambda_continuously_cooled_factor) < EM::TOL) ){
    
    // save the current model and gradient
    this->m_xp = x;
    this->m_gradp = this->m_grad;

    // search step length by using inexact line search algorithm  
    // and then to update fx, x and gradient (m_grad)
    Inv_Basis::line_search_strong_Wolfe(fx, x, this->m_grad, step, this->m_drt,
                                        this->m_xp, this->m_param, iter);
    std::cout << "Step length is: " << step << std::endl;
    screen_output << "Step length is: " << step << std::endl;
    output_alfa << step << std::endl;
    output_iter << iter << std::endl;

    // !!! for next interation
    Inv_Basis::compute_misfit_RMS_roughness(fx, x, lambda);
    this->k++;
    lambda_temp = lambda;
    // !! cool lambda
    // if the RMS is increased, we cool lambda as well
    double RMS_reduction = 0.0;
    RMS_reduction = Inv_Basis::RMS[this->k-1] - Inv_Basis::RMS[this->k];
    // the tolerance of the reduction fraction of RMS
    double tol_rms_reduc_frac = Inv_Basis::inv_para.tol_rms_reduc_frac;
    if( RMS_reduction < (tol_rms_reduc_frac * Inv_Basis::RMS[this->k-1]) )
      // cool lambda
      lambda /= Inv_Basis::inv_para.cool_factor;
    Inv_Basis::lambda_history.push_back(lambda);
    
    // calculate vector diff_grad_phi_m = gradp_phi_m_{k} - gradp_phi_m_{k-1} 
    if(Inv_Basis::startup.data_method == "MT1D"){
      std::cout << "L-BFGS inversion for MT1D is under construction" << std::endl;
      exit(1);
    }
    else if(Inv_Basis::startup.data_method == "MT2D"){
      std::cout << "L-BFGS inversion for MT2D is under construction" << std::endl;
      exit(1);
    }
    else // for MT3D
      diff_grad_phi_m = 2.0 * Inv_Basis::Cm_inv * (x - m_xp);  

    // !update vector pairs s, y and d
    // s_{k} = x_{k+1} - x_k, y_{k} = g_{k+1} - g_k 
    // d[k] = gradp_phi_m_{k} - gradp_phi_m_{k-1}
    this->add_correction(x - this->m_xp, m_grad - this->m_gradp, diff_grad_phi_m);
    
    // !correct the past y vectors and the corresponding ys and theta
    //  once lambda was changed
    if(std::abs(lambda - lambda_temp) > EM::TOL){ // lambda was changed 
      // do correction
      correct_y_ys_theta(lambda, lambda_temp);
    }
    
    // recursive formula to compute search direction of the present L-BFGS iter
    // m_drt = -H * m_grad
    this->apply_Hv(this->m_grad, -1.0, this->m_drt);
    // reset step = 1.0, initial guess of step for the next line search
    step = 1.0;
    iter = 0;

///////////////////////// !!! end the core of L-BFGS //////////////////////////

    // used for judging whether continuously cool lambda too many times
    unsigned int temp_tol_times = Inv_Basis::inv_para.tol_times_cont_cool_lambda;
    // (this->k+1) is the size of lambda_history vector
    if( (this->k+1) > temp_tol_times ){
      lambda_continuously_cooled_factor = 1.0;
      for(unsigned i = this->k - temp_tol_times; i < this->k; i++){
        lambda_continuously_cooled_factor *= lambda_history[i] 
                                           - lambda_history[i+1];
      }
    }
    // !! write model into file or print cooresponding info on screen and
    // update m_k for the next optimal iteration except the last iteration,
    // the last iteration is used for generating the initial model
    // of the next adjusted mesh for adaptive inversion
    if( (this->k <= Inv_Basis::tol_times) &&
        (lambda_temp > Inv_Basis::inv_para.tol_lambda) &&
        (std::abs(lambda_continuously_cooled_factor) < EM::TOL) ) {
      // update m_k, which will be used as the initial model of the 
      // next adjusted mesh
      Inv_Basis::m_k = Inv_Basis::m_kp1;
      //print the updated model and the corresponding RMS misfit on screen
      std::cout << "Iter = " << this->k << ", lambda = " << lambda_temp
                << ", RMS = " << Inv_Basis::RMS[this->k]
                << ", roughness = " << Inv_Basis::roughness[this->k] << std::endl;
      screen_output << "Iter = " << this->k << ", lambda = " << lambda_temp
                    << ", RMS = " << Inv_Basis::RMS[this->k]
                    << ", roughness = " << Inv_Basis::roughness[this->k] << std::endl;

      if(Inv_Basis::startup.data_method == "MT1D"){
        std::cout << "This module for MT1D is under construction" << std::endl;
        exit(1);
      }  
      else{ // for MT2D and MT3D
        /*
        std::cout << "The current resistivity model m_" << this->k 
                  <<" in linear domain is:" << std::endl;
        for(unsigned int i = 0; i < Inv_Basis::m_k.size(); i++){
          double nomin = 0, denomin = 0., cond = 0.;
          nomin = Inv_Basis::inv_para.a + ( Inv_Basis::inv_para.b 
                  * std::exp( Inv_Basis::inv_para.n * Inv_Basis::m_k(i) ) );
          denomin = 1.0 + std::exp( Inv_Basis::inv_para.n * Inv_Basis::m_k(i) );
          // conductivity in linear space
          cond = nomin / denomin;
          if( std::isnan( cond ) )
            // we set sigma = b (upper bound), if it is equals to NAN
            cond =  Inv_Basis::inv_para.b;
          std::cout << 1.0 / cond << '\n';
        }
        */
        if(this->k != Inv_Basis::tol_times){
          std::cout << '\n' << "Inversion iteration " << this->k + 1 << ":\n";
          screen_output << '\n' << "Inversion iteration " << this->k + 1 << ":\n";
        }
      }
      // write the RMS, model roughness and inversion results into the files
      output_RMS << Inv_Basis::RMS[this->k] << std::endl;
      output_roughness << Inv_Basis::roughness[this->k] << std::endl;
      output_lambda << lambda_temp << std::endl;
      Inv_Basis::Output_Inv_Mod(this->n_adjust, this->k);
      Inv_Basis::output_data_and_residual(this->n_adjust, this->k);
      
      if(this->k == Inv_Basis::tol_times){
        std::cout << "L-BFGS iteration stopped due to tol of max iteration times\n";
        screen_output << "L-BFGS iteration stopped due to tol of max iteration times\n";
      }
      if(!(Inv_Basis::RMS[this->k] > Inv_Basis::inv_para.tol_rms)){
        std::cout << "L-BFGS iteration stopped due to tol of RMS\n";
        screen_output << "L-BFGS iteration stopped due to tol of RMS\n";
      }
      else if(Inv_Basis::RMS[this->k] > Inv_Basis::inv_para.tol_rms){
        ; // do nothing
      }
      else{
        std::cout << "RMS of the last iteration is " << Inv_Basis::RMS[this->k] 
                  << std::endl;
        screen_output << "RMS of the last iteration is " << Inv_Basis::RMS[this->k] 
                      << std::endl;
        std::cout << "Abnormally stopped, possibly '(-)nan' emerged during the " 
                  << "inversion process!" << std::endl;
        screen_output << "Abnormally stopped, possibly '(-)nan' emerged during the " 
                      << "inversion process!" << std::endl;
        std::abort();
      }
    } // end if
    else if( !(lambda_temp > Inv_Basis::inv_para.tol_lambda) ){
      std::cout << "L-BFGS iteration stopped due to tol of lambda\n";
      screen_output << "L-BFGS iteration stopped due to tol of lambda\n";
      this->k--;
    }
    else if(std::abs(lambda_continuously_cooled_factor) > 0.){
      std::cout << "L-BFGS iteration stopped due to tol times of cooling" 
                << " lambda continuously\n";
      screen_output << "L-BFGS iteration stopped due to tol times of cooling" 
                    << " lambda continuously\n";
      this->k--;
    }
    else{
      std::cout << "RMS of the last iteration is " << Inv_Basis::RMS[this->k] 
                << std::endl;
      screen_output << "RMS of the last iteration is " << Inv_Basis::RMS[this->k] 
                    << std::endl;
      std::cout << "Abnormally stopped, possibly '(-)nan' emerged during the " 
                << "inversion process!" << std::endl;
      screen_output << "Abnormally stopped, possibly '(-)nan' emerged during the " 
                    << "inversion process!" << std::endl;
      std::abort();
    }
  } // end while

/*  
  output_RMS << '\n' 
             << "The total number of RMS is: " << this->k << std::endl;
  output_RMS << "The total number of L-BFGS iterations is: " 
             << this->k-1 << std::endl;
*/
  output_RMS.close();
  output_roughness.close();
  output_lambda.close();
  output_alfa.close();
  output_iter.close();

  return;
} 

void LBFGS_Inv::add_correction(const EM::RefConstVec& s, const EM::RefConstVec& y,
                               const EM::RefConstVec& d)
{
  /////////////////////////////////////////////////////////////////////////////
  //              Add correction vectors to the L-BFGS matrix                //
  ///////////////////////////////////////////////////////////////////////////// 
  // m_ptr:
  // a Pointer to locate the most recent history, 1 <= m_ptr <= m_param.m
  // details: s and y vectors are stored in cyclic order. for example, 
  // if the current s-vector is stored in m_s[, m-1], then in the next 
  // iteration m_s[, 0] will be overwritten. m_s[, m_ptr-1] points to the
  // most recent history (i.e., k-1 in book), and m_s[, m_ptr % m] points to 
  // the most distant one (i.e., k-m).
  const int loc = this->m_ptr % this->m_m;

  this->m_d.col(loc).noalias() = d;
  this->m_s.col(loc).noalias() = s;
  this->m_y.col(loc).noalias() = y;

  // ys = y's = 1/rho
  const double ys = this->m_s.col(loc).dot(this->m_y.col(loc));
  this->m_ys[loc] = ys;
  this->m_theta = this->m_y.col(loc).squaredNorm() / ys;

  if(this->m_ncorr < this->m_m)
      this->m_ncorr++;

  this->m_ptr = loc + 1;
  
  return;
}

void LBFGS_Inv::correct_y_ys_theta(double lambda_c, double lambda_p)
{
  EM::VectorXD corrected_y;
  double corrected_ys;
  // !correct the most recent y vector and corresponding ys and theta
  int j = this->m_ptr - 1;
  corrected_y = this->m_y.col(j) + (lambda_c - lambda_p) * this->m_d.col(j);
  corrected_ys = this->m_s.col(j).dot(corrected_y);
  if(corrected_ys > 0)
  {
    // DO CORRECTION!
    this->m_y.col(j).noalias() = corrected_y;
    const double ys = corrected_ys;
    this->m_ys[j] = ys;
    this->m_theta = corrected_y.squaredNorm() / ys;
  }                                                
                                                                                       
  // !correct other y vectors and other ys values
  j = this->m_ptr % this->m_m;
  for(int i = 0; i < this->m_ncorr; i++)
  {
    // loop clockwise from the most recent one to the most distant one
    // i.e., from k-1 to k-m, shown in Algorithm 9.1 of that book
    j = (j + this->m_m - 1) % this->m_m;
    if(i != 0){ // i = 0 has been dealed
      corrected_y = this->m_y.col(j) + (lambda_c - lambda_p) * this->m_d.col(j);
      corrected_ys = this->m_s.col(j).dot(corrected_y);
      if(corrected_ys > 0)
      {
        // DO CORRECTION!
        this->m_y.col(j).noalias() = corrected_y;
        this->m_ys[j] = corrected_ys;
      }
    }
  }
  
  return;
}

void LBFGS_Inv::apply_Hv(const EM::VectorXD& v, const double& a, EM::VectorXD& res)
{
/////////////////////////////////////////////////////////////////////////////////
//Recursive formula to compute a * H * v, where a is a scalar, and v is [n x 1]//
//   H0 = (1/theta) * I is the initial approximation to H of each iteration    //
//   cf. Algorithm 9.1 (L-BFGS two-loop recursion) of J., & Wright, S. (1999). //
//                     Numerical optimization Nocedal                          //
/////////////////////////////////////////////////////////////////////////////////
  res.resize(v.size());
  // Loop 1
  res.noalias() = a * v;
  // this->m_ptr will be updated after calling the add_correction function
  // (j + this->m_m - 1) % this->m_m represent the most recent history 
  int j = this->m_ptr % this->m_m;
  for(int i = 0; i < this->m_ncorr; i++)
  {
      // loop clockwise from the most recent one to the most distant one
      // i.e., from k-1 to k-m, shown in Algorithm 9.1 of that book
      j = (j + this->m_m - 1) % this->m_m;
      this->m_alpha[j] = this->m_s.col(j).dot(res) / this->m_ys[j];
      res.noalias() -= this->m_alpha[j] * this->m_y.col(j);
  }

  // Apply initial H0 by using scaled identity matrix
  // H_k^{0} = sk'yk / yk'yk, k = 0, 1, ... 
  // cf. P226 in second edition of "Numerical Optimization" book
  res /= this->m_theta; 

  // Loop 2
  // here the value of loop variable is slightly different from
  // that in second edition of "Numerical Optimization" book
  for(int i = 0; i < this->m_ncorr; i++)
  {
      const double beta = this->m_y.col(j).dot(res) / this->m_ys[j];
      res.noalias() += (this->m_alpha[j] - beta) * this->m_s.col(j);
      // loop anticlockwise from the most distant one to the most recent one
      // i.e., from k-m to k-1, shown in Algorithm 9.1 of that book
      j = (j + 1) % this->m_m;
  }
  
  return;
}

/**************************used for adaptive inversion************************/
void LBFGS_Inv::get_RMS_last(double& RMS_last)
{
  unsigned int LBFGS_iter_times = Inv_Basis::RMS.size();
  RMS_last = Inv_Basis::RMS[LBFGS_iter_times - 1];

  return;
}

void LBFGS_Inv::get_lambda_last(double& lambda_last)
{
  lambda_last = Inv_Basis::lambda;
  
  return;
}

void LBFGS_Inv::get_iter_number(unsigned int& n_iter)
{
  n_iter = this->k;

  return;
}

void LBFGS_Inv::update_initial_lambda(unsigned int k)
{
  double total_cool_factor = std::pow(Inv_Basis::inv_para.cool_factor, k);
  Inv_Basis::lambda = Inv_Basis::inv_para.lambda0 / total_cool_factor;
}

void LBFGS_Inv::update_initial_lambda(double lambda0)
{
  Inv_Basis::lambda = lambda0;

  return;
}

void LBFGS_Inv::update_LBFGS_iter_times(unsigned int tol_LBFGS_times)
{
  Inv_Basis::tol_times = tol_LBFGS_times;

  return;
}

/****************used for test of gradient calculation for MT3D***************/
void LBFGS_Inv::test_grad_phi_d()
{
  // !!! test of gradient calculation for MT3D
  // directional derivative computed by using sensitivity equation equation
  double dir_deri_sens_eq = 0.0;
  // directional derivative computed by using pertubation method
  double dir_deri_pert = 0.0;
  // for storing the gradient of \phi_d at model mk
  double phi_d_mk       = 0.0;
  // for storing the gradient of \phi_d at model mk + m_pertubation
  double phi_d_mK_pert  = 0.0;

  // unit directional vector
  EM::VectorXD p;
  p.resize(Inv_Basis::M);
  double norm_p = 0.0;
  for(unsigned int i = 0; i < Inv_Basis::M; i++){
    //p(i) = std::log10(i + 1);
    p(i) = 2.0;
  }
  norm_p = p.norm();
  // unitization
  p = p / norm_p;
  
  // !!! compute the project of gradient of \phi_d_mk at director p by using 
  // sensitivity equation method, i.e., compute the directional derivative
  EM::VectorXD grad_phi_d;
  grad_phi_d.resize(Inv_Basis::M);
  // !! compute the gradient of \phi_d_mk
  Inv_Basis::comp_phi_g_m_sen_eq(Inv_Basis::m_k, 0, phi_d_mk, grad_phi_d);
  // !! compute the directional derivative
  dir_deri_sens_eq = grad_phi_d.transpose() * p;

  std::cout << "The directional derivative computed by using sensitivity "
            << "equation approach is:\n" 
            << dir_deri_sens_eq << std::endl;
  screen_output << "The directional derivative computed by using sensitivity "
                << "equation approach is:\n" 
                << dir_deri_sens_eq << std::endl;
  
  // !!! compute directional derivative by using pertubation method
  double h = 1.0e-4;
  EM::VectorXD m_perted = Inv_Basis::m_k + h * p;
  // !! compute the gradient of \phi_d_mk_pertubation
  Inv_Basis::comp_phi_g_m_sen_eq(m_perted, 0, phi_d_mK_pert, grad_phi_d);
  // !! compute the directional derivative
  dir_deri_pert = (phi_d_mK_pert - phi_d_mk) / h;

  std::cout << "The directional derivative computed by using pertubation "
            << "approach is:\n" 
            << dir_deri_pert << std::endl;
  screen_output << "The directional derivative computed by using pertubation "
                << "approach is:\n" 
                << dir_deri_pert << std::endl;
  
  double re_error = std::abs( (dir_deri_sens_eq - dir_deri_pert) / dir_deri_pert );
  std::cout << "Absolute relative error: " << re_error * 100 << "%\n";
  screen_output << "Absolute relative error: " << re_error * 100 << "%\n";

  return;
}
