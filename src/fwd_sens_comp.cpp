/******************************************************************************
 *    Function       : Forwad modeling interface to produce synthetic data,   * 
                       compute forward data and calculate sensitivity and     *
                       gradient of data misfit for 1D, 2D and 3D MT problems  *
 *    Author         : Huang Chen and Zhengyong Ren                           *
 *    Copyright      : Huang Chen, 2019 and Zhengyong Ren, 2017               *
 *    Email          : chenhuang@cqu.edu.cn; renzhengyong@csu.edu.cn       *
 *    Created time   : 2019.07.22                                             *
 *    Last revision  : 2023.8.25                                              *
 ******************************************************************************/
#include "fwd_sens_comp.h"
 // used for generating screen output log
 extern std::ofstream screen_output;

// !!! constructors for Fwd_only (produce synthetic data)
Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_1D _fwd_para_1d, std::string _data_method)
{
  this->fwd_para_1d = _fwd_para_1d;
  this->data_method = _data_method;
}

Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_2D _fwd_para_2d, std::string _data_method)
{
  this->fwd_para_2d   = _fwd_para_2d;
  this->data_method = _data_method;
}

Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_3D _fwd_para_3d, std::string _data_method)
{
  // initialize member varialbles
  this->fwd_para_3d = _fwd_para_3d;
  this->data_method = _data_method;
}

// !!! overloaded constructors for inversion
Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_1D _fwd_para_1d, Data_1D _data_1d)
{
  this->fwd_para_1d = _fwd_para_1d;
  this->data_1d     = _data_1d;
}

Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_2D _fwd_para_2d, Data_2D _data_2d)
{
  this->fwd_para_2d = _fwd_para_2d;
  this->data_2d     = _data_2d;
}

Fwd_Sens_Comp::Fwd_Sens_Comp(Fwd_Para_3D _fwd_para_3d, Data_3D _data_3d)
{
  // initialize member varialbles
  this->fwd_para_3d = _fwd_para_3d;
  this->data_3d     = _data_3d;
}


///////////////////////////////////////////////////////////////////////////////
/*********************for 3D MT forward modelling (fwd_only)******************/
///////////////////////////////////////////////////////////////////////////////
void Fwd_Sens_Comp::Produce_MT3d_Synthetic_Data(double rel_err, double abs_err,
                                                double rel_tipper_err,
                                                double abs_tipper_err,
                                                std::string inv_approach)
{
/*****************compute and contaminate 3D MT forward data******************/
  Fwd_Para_3D& parameter = this->fwd_para_3d;
  // for storing site information;
  std::vector<Node> sites;

  if(parameter.algorithm==1){ // using adaptive fem
    std::cout << "This module is under construction\n";
    exit(1);
  }
  else{ // using standard fem
    // for storing impedances of all frequencies
    std::vector< EM::MatrixXC > Z;
    Z.resize(parameter.n_f);
    // for storing tipper response of all frequencies
    std::vector< EM::MatrixXC > T;
    T.resize(parameter.n_f);
    // for data comtaminating
    // construct a trivial random generator engine from a time-based seed:
    // Each run will get different random values
    // We can also use constant seed to allow for reproducibility, but here
    // we use a dynamic seed
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // step1: define a clock based random engine object
    std::default_random_engine generator(seed);
    // temporary variable for error and standard deviation of each data
    double err = 0.0, std_dev = 0.0;
    double abs_Zxx = 0.0, abs_Zxy = 0.0, abs_Zyx = 0.0, abs_Zyy = 0.0;
    double max_abs_Tzx_Tzy = 0.0;

    // !!!produce MT3d forward impedance data
    // parallelization over frequencies
    // loop frequencies
    #ifdef _OPENMP
      #pragma omp parallel for
    #endif
    for(int i = 0; i < parameter.n_f; i++){
      #ifdef _OPENMP
      ///////////////////////////////////////////////////////////////////// ///
      //!!! added to solve the multi-threads closed problem of MKL (Pardiso) //
      //                under OpenMP parallel region.                        //
      /////////////////////////////////////////////////////////////////////////
        // set the number of threads for the nested mkl parallel region
        //mkl_set_num_threads_local(12); //the same as omp_set_num_threads();
        // disable adjustment of the number of threads
        mkl_set_dynamic(0); 
        // has been deprecated by the latest compiler, has the same effect as 
        // environment variable OMP_NESTED=true,
        //omp_set_nested(1);
        // used to control MKL PARDISO to spawn only the first level of 
        // nested threads and avoid significant overhead 
        // If the number of active levels requested exceeds the number of active
        // levels of parallelism supported by the implementation, the value of the 
        // max-active-levels-var ICV will be set to the number of active levels 
        // supported by the implementation.
        omp_set_max_active_levels(2);
      #endif
      // linear source polarization TE
      BC bc_te(parameter.EP, parameter.f[i], parameter.theta, EM::TE);
      // linear source polarization TM
      BC bc_tm(parameter.EP, parameter.f[i], parameter.theta, EM::TM); 
      // FEM object
      std::stringstream counter;
      counter << parameter.starting_counter;
      std::string s = counter.str();
      std::string model_name = parameter.starting_model[i]+std::string(".")+s;
      Mesh3D mesh(model_name); 
      mesh.prepare_mesh_by_reading_files();
      FEM fem(mesh, parameter.f[i], bc_te, bc_tm, inv_approach, parameter.marker);           
      fem.solve();   
      std::stringstream counter_f;
      counter_f << i + 1;
      std::string s_f = counter_f.str();
      std::string solution_name = "f" + s_f + "_FEM_solution";
      std::ofstream output(solution_name.c_str()); 
      
      PostProcess postp(fem, 
                        parameter.marker, 
                        parameter.u, parameter.v, 
                        parameter.mu, 
                        output);
      // return sites only once for all the frequencies,
      // as all the frequencies share the same sites
      if(i == 0)
        postp.Return_Sites(sites);
      // return impedances for each frequency
      postp.Return_Impedances_Tippers(Z[i], T[i]);
      output.close();
      // delete the solution file of each frequency
      std::string command = "rm " + solution_name;
      std::system( command.c_str() );
    }
   
    // writing file
    std::string data_file_name;
    data_file_name = "MT3D_Synthetic.data";
    std::ofstream out(data_file_name.c_str());
    std::string Z_comp[4];
    Z_comp[0] = "ZXX"; Z_comp[1] = "ZXY";
    Z_comp[2] = "ZYX"; Z_comp[3] = "ZYY";
    std::string T_comp[2];
    T_comp[0] = "TX"; T_comp[1] = "TY";
    std::cout << "\nMT3D synthetic impedance and tipper data were"
              << " wroten into file.\n";
    screen_output << "\nMT3D synthetic impedance and tipper data were"
                  << " wroten into file.\n";
    out << "# 3D MT synthetic data produced by standard FEM \n";
    out << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error \n";
    out << "> Full_Impedance_Tipper \n";
    out << "> exp(-i_omega_t) \n";
    out << "> [V/m]/[A/m] \n";
    out << "> " << 0.0 << "\n";
    out << "> " << 0.0 << "  " << 0.0 << "\n";
    out << "> " << parameter.n_f << "  " << sites.size() << '\n';

    // !!!contaminate impedance data and output into data file
    for(unsigned int i = 0; i < sites.size(); i++)
      for(unsigned int j = 0; j < parameter.n_f; j++){ 
        abs_Zxx = std::abs( Z[j](i,0) );
        abs_Zxy = std::abs( Z[j](i,1) );
        abs_Zyx = std::abs( Z[j](i,2) );
        abs_Zyy = std::abs( Z[j](i,3) );
        if( (abs_Zxy < EM::TOL) && (abs_Zyx < EM::TOL) )
          std_dev = std::max(abs_Zxx, abs_Zyy) * rel_err + abs_err;
        else 
        {
          if( abs_Zxy < EM::TOL )
            abs_Zxy = abs_Zyx;
          if( abs_Zyx < EM::TOL )
            abs_Zyx = abs_Zxy;
          // the stdandard deviation is lager than or equal to zero (i.e., noise-free)
          // this formula is from (Alexander Grayver, 2013) GJI paper
          // abs(Zxy) * abs(Zyx) = abs(Zxy * Zyx)
          std_dev = std::sqrt(abs_Zxy * abs_Zyx) * rel_err + abs_err;
        }
        // step2: define a normal distribution object with 
        // E = 0 and \sigma = std_dev
        std::normal_distribution< double > distribution (0.0, std_dev);

        // !! contaminate 3D forward impedance data with standard deviation of
        // of n% * sqrt(abs(Zxy*Zyx)) for all components icluding both real
        // and imaginary parts
        for(unsigned int k = 0; k < 4; k++){
          //if( (k == 1) || (k == 2) ) // only for off-diagonal components
          { 
            // step3: generate random number by calling member function operator()
            err = distribution( generator );
            Z[j](i,k).real( Z[j](i,k).real() + err );
            err = distribution( generator );
            Z[j](i,k).imag( Z[j](i,k).imag() + err );
            //  write impedance data into file
            out << setiosflags(std::ios::scientific)
                << 1.0 / parameter.f[j] << '\t' 
                // Code CG_Lat CG_Lon
                << sites[i].get_id() << '\t' 
                << 0.0 << '\t' 
                << 0.0 << '\t'
                // X(m) Y(m) Z(m)
                << resetiosflags(std::ios::scientific) 
                << std::setprecision(8)
                << setiosflags(std::ios::right)
                << std::setw(10)
                << sites[i].get_xyz()[0] << '\t' 
                << setiosflags(std::ios::right)
                << std::setw(10)
                << sites[i].get_xyz()[1] << '\t'
                << setiosflags(std::ios::right)
                << std::setw(10)
                << sites[i].get_xyz()[2] << '\t'
                // Component Real Imag
                << Z_comp[k] << '\t' 
                << setiosflags(std::ios::scientific)
                << std::setprecision(6)
                << setiosflags(std::ios::right)
                << std::setw(15)
                << std::real( Z[j](i,k) ) << '\t'
                << setiosflags(std::ios::right)
                << std::setw(15) 
                << std::imag( Z[j](i,k) ) << '\t'
                // Error
                << setiosflags(std::ios::right)
                << std::setw(15) 
                << std_dev << '\n';
          }
        }
      }

    // !!!contaminate tipper data and output into data file
    for(unsigned int i = 0; i < sites.size(); i++)
      for(unsigned int j = 0; j < parameter.n_f; j++){ 
        max_abs_Tzx_Tzy = std::max( std::abs( T[j](i,0) ), 
                                    std::abs( T[j](i,1) ) );
        // the stdandard deviation is lager than or equal to zero (i.e., noise-free)
        std_dev = max_abs_Tzx_Tzy * rel_tipper_err + abs_tipper_err;
        // step2: define a normal distribution objects with 
        // E = 0 and \sigma = std_dev
        std::normal_distribution< double > distribution (0.0, std_dev);

        // contaminate 3D synthetic tipper data with standard deviation of 
        // of n% * max( abs(TX), abs(TY) ) + abs_tipper_err
        for(unsigned int k = 0; k < 2; k++){
          // step3: generate random number by calling member function operator()
          err = distribution( generator );
          T[j](i,k).real( T[j](i,k).real() + err );
          err = distribution( generator );
          T[j](i,k).imag( T[j](i,k).imag() + err );
          // write tipper data into file
          out << setiosflags(std::ios::scientific)
              << 1.0 / parameter.f[j] << '\t' 
              // Code CG_Lat CG_Lon
              << sites[i].get_id() << '\t' 
              << 0.0 << '\t' 
              << 0.0 << '\t'
              // X(m) Y(m) Z(m)
              << resetiosflags(std::ios::scientific) 
              << std::setprecision(8)
              << sites[i].get_xyz()[0] << '\t' 
              << setiosflags(std::ios::right)
              << std::setw(10)
              << sites[i].get_xyz()[1] << '\t'
              << setiosflags(std::ios::right)
              << std::setw(10)
              << sites[i].get_xyz()[2] << '\t'
              // Component Real Imag
              << T_comp[k] << '\t' 
              << setiosflags(std::ios::scientific)
              << std::setprecision(6)
              << setiosflags(std::ios::right)
              << std::setw(15)
              << std::real( T[j](i,k) ) << '\t'
              << setiosflags(std::ios::right)
              << std::setw(15) 
              << std::imag( T[j](i,k) ) << '\t'
              // Error bar or error floor
              << setiosflags(std::ios::right)
              << std::setw(15) 
              << std_dev << '\n';
        }
      }
    out.close();
    std::cout << "File writing finished. \n";
    screen_output << "File writing finished. \n";
  }

  return;
}

void Fwd_Sens_Comp::Comp_MT3d_Fwd_Responses(EM::VectorXD& F_m, const EM::VectorXD& m,
                                const std::string& inv_approach, const double abn[3])
{
  /************compute forward data of initial or inverted models*************/
  /*****************calculation of misfit and data outputting*****************/
  Fwd_Para_3D& parameter = this->fwd_para_3d;
  // for storing site information;
  std::vector<Node> sites;
  assert(parameter.n_f  = this->data_3d.n_f);
  unsigned int n_s      = this->data_3d.n_sites;
  unsigned int n_f      = this->data_3d.n_f;
  unsigned int n_c      = this->data_3d.n_comps;
  unsigned int N        = this->data_3d.N;
  std::string data_type = this->data_3d.data_type;
  // half number of the data
  // which is used for storing the imaginary part of impedance
  unsigned int half_N  = N / 2;

  // checking the size
  assert(F_m.size() == N);
  F_m.setZero();

  if(parameter.algorithm == 1){
    // to be done for maping of forward mesh and inversion mesh
    std::cout << "This module is under construction, please turn to "
              << " standard FEM" << std::endl;
    std::abort();
  }
  else{
    // for storing (off-diagonal) impedances and (tippers) of all frequencies
    std::vector< EM::MatrixXC > Z, ZpT, diagonal_Z, diagonal_ZpT;
    if( data_type == "Full_Impedance" )
      // for storing full impedances of all frequencies
      Z.resize(parameter.n_f);
    else if( data_type == "Full_Impedance_plus_Tipper" )
      ZpT.resize(parameter.n_f);
    else if( data_type == "Off_Diagonal_Impedance" )
      diagonal_Z.resize(parameter.n_f);
    else // for data_type of "Off_Diagonal_Impedance_plus_Tipper"
      diagonal_ZpT.resize(parameter.n_f);
    
    // loop frequencies
    #ifdef _OPENMP
      #pragma omp parallel for
    #endif
    for(int i = 0; i < parameter.n_f; i++){
      #ifdef _OPENMP
      ///////////////////////////////////////////////////////////////////// ///
      //!!! added to solve the multi-threads closed problem of MKL (Pardiso) //
      //                under OpenMP parallel region.                        //
      /////////////////////////////////////////////////////////////////////////
        // set the number of threads for the nested mkl parallel region
        //mkl_set_num_threads_local(12); //the same as omp_set_num_threads();
        // disable adjustment of the number of threads
        mkl_set_dynamic(0); 
        // has been deprecated by the latest compiler, has the same effect as 
        // environment variable OMP_NESTED=true,
        //omp_set_nested(1);
        // used to control MKL PARDISO to spawn only the first level of 
        // nested threads and avoid significant overhead 
        // If the number of active levels requested exceeds the number of active
        // levels of parallelism supported by the implementation, the value of the 
        // max-active-levels-var ICV will be set to the number of active levels 
        // supported by the implementation.
        omp_set_max_active_levels(2);
      #endif
      // linear source polarization TE
      BC bc_te(parameter.EP, parameter.f[i], parameter.theta, EM::TE);
      // linear source polarization TM
      BC bc_tm(parameter.EP, parameter.f[i], parameter.theta, EM::TM);
      // FEM object
      std::stringstream counter;
      counter << parameter.starting_counter;
      std::string s = counter.str();
      std::string model_name = parameter.starting_model[0]+std::string(".")+s;
      Mesh3D mesh(model_name); 
      mesh.prepare_mesh_by_reading_files();
      FEM fem(mesh, parameter.f[i], bc_te, bc_tm, inv_approach,
              parameter.marker, parameter.u, parameter.v, m, abn);
      fem.solve();   
      // false represents we don't need to calculate the global interpolator operators
      PostProcess postp(fem, parameter.marker, parameter.mu, data_type, false);
      // return sites only once for all the frequencies,
      // as all the frequencies share the same sites
      if(i == 0)
        postp.Return_Sites(sites);
      
      if(data_type == "Full_Impedance")
      // return impedances for each frequency
        postp.Return_Impedances(Z[i]);
      else if(data_type == "Full_Impedance_plus_Tipper")
        postp.Return_Impedances_plus_Tippers(ZpT[i]);
      else if(data_type == "Off_Diagonal_Impedance")
        postp.Return_Diagonal_Impedances(diagonal_Z[i]);
      else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
        postp.Return_Diagonal_Impedances_plus_Tippers(diagonal_ZpT[i]);
    }

    // check site number
    unsigned int n_sites = sites.size();
    assert(n_s == n_sites);
    // !!!compute Forward data F[m]
    if(data_type == "Full_Impedance"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // full impedance responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( Z[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( Z[i](j,k) );
          }
    }
    else if(data_type == "Full_Impedance_plus_Tipper"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // impedance plus tipper responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( ZpT[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( ZpT[i](j,k) );
          }
    }
    else if(data_type == "Off_Diagonal_Impedance"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // off-diagonal impedance responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( diagonal_Z[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( diagonal_Z[i](j,k) );
          }
    }
    else{ // for data type of "Off_Diagonal_Impedance_plus_Tipper"
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // off-diagonal impedance plus tipper responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( diagonal_ZpT[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( diagonal_ZpT[i](j,k) );
          }
    }
  } // end if (algorithm == 0)

  return;
}

void Fwd_Sens_Comp::compute_sensitivity(FEM& fem, PostProcess& post_p,
                                        EM::MatrixXD& J, unsigned int i,
               std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj)
{
  // the number of the measuring sites
  unsigned int n_s      = this->data_3d.n_sites;
  // the number of the components
  unsigned int n_c      = this->data_3d.n_comps; 
  // the number of the tetrahedral elements of inversion mesh
  unsigned int n_m      = fwd_tet_id_in_mj.size(); 
  // half number of the data
  unsigned int half_N   = this->data_3d.N / 2;
  std::string data_type = this->data_3d.data_type;

  EM::MatrixXC X_1, X_2, L_1_2_T;
  // calculate X_1 and X_2
  // A * X_1,2 = (L_1,2)^T
  // transpose of the discrete global interpolation operators for source 1
  // changing from sparse matrix to dense matrix will accelerate the following
  // adjoint forward modelling, i.e., (fem.solver).solve( ( L_1_2_T ), 
  // the speed-up ratio had ever reached to 4; but need much more memory, i.e.
  // the memory costs are 1.8G for dense matrix and 1.3G for sparse matrix
  if(data_type == "Full_Impedance")
    L_1_2_T = (post_p.L_1_Z).transpose();
  else if(data_type == "Full_Impedance_plus_Tipper")
    L_1_2_T = (post_p.L_1_ZpT).transpose();
  else if(data_type == "Off_Diagonal_Impedance")
    L_1_2_T = (post_p.L_1_diagonal_Z).transpose();
  else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
    L_1_2_T = (post_p.L_1_diagonal_ZpT).transpose();
  X_1 = (fem.solver).solve( L_1_2_T );

  // transpose of the discrete global interpolation operators for source 2
  if(data_type == "Full_Impedance")
    L_1_2_T = (post_p.L_2_Z).transpose();
  else if(data_type == "Full_Impedance_plus_Tipper")
    L_1_2_T = (post_p.L_2_ZpT).transpose();
  else if(data_type == "Off_Diagonal_Impedance")
    L_1_2_T = (post_p.L_2_diagonal_Z).transpose();
  else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
    L_1_2_T = (post_p.L_2_diagonal_ZpT).transpose();
  X_2 = (fem.solver).solve( L_1_2_T );
  // P_1 = - \partial A / \partial m * fem.E0
  // P_2 = - \partial A / \partial m * fem.E1
  // the j-th column represents the vector of \partial A / \partial m_j * fem.E0,1
  EM::EigenSparseMatrixXC P_1, P_2;
  // calculate P_1 and P_2
  fem.comp_partial_A_partial_m(P_1, P_2, fwd_tet_id_in_mj);
  P_1 = - P_1;
  P_2 = - P_2;
  // directly calculate J_i by using matrix computation
  // real part of J_i
  J.block(i * n_s * n_c, 0, n_s * n_c, n_m) = 
  (X_1.transpose() * P_1 + X_2.transpose() * P_2).real();
  // imaginary part of J_i
  J.block(i * n_s * n_c + half_N, 0, n_s * n_c, n_m) =
  (X_1.transpose() * P_1 + X_2.transpose() * P_2).imag();

  return;
}

void Fwd_Sens_Comp::comp_MT3d_fwd_grad_phi_d(EM::VectorXD& F_m, EM::VectorXD& grad_phi_d,
                                             const EM::VectorXD& m, const EM::VectorXD& d,
                                             EM::EigenSparseMatrixXD& Wd,
                                             const std::string& inv_approach,
                               std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj,
                                             const double abn[3]) 
{
/*************Compute forward data and gradient of the misfit term \phi_d***************/
  Fwd_Para_3D& parameter = this->fwd_para_3d;
  // for storing site information;
  std::vector<Node> sites;
  assert(parameter.n_f  = this->data_3d.n_f);
  unsigned int n_s      = this->data_3d.n_sites;
  unsigned int n_f      = this->data_3d.n_f;
  unsigned int n_c      = this->data_3d.n_comps;
  unsigned int N        = this->data_3d.N;
  std::string data_type = this->data_3d.data_type;
  // half number of the data
  // which is used for storing the imaginary part of impedance
  unsigned int half_N  = N / 2;

  // size double checking
  assert(F_m.size() == N);
  F_m.setZero();
  if(parameter.algorithm == 1){
    // using adaptive fem
    std::cout << "This module is under construction, please turn to "
              << " standard FEM" << std::endl;
    std::abort();
  }
  else{
    // for storing (off-diagonal) impedances and (tippers) of all frequencies
    std::vector< EM::MatrixXC > Z, ZpT, diagonal_Z, diagonal_ZpT;
    if( data_type == "Full_Impedance" )
      // for storing full impedances of all frequencies
      Z.resize(parameter.n_f);
    else if( data_type == "Full_Impedance_plus_Tipper" )
      ZpT.resize(parameter.n_f);
    else if( data_type == "Off_Diagonal_Impedance" )
      diagonal_Z.resize(parameter.n_f);
    else // for data_type of "Off_Diagonal_Impedance_plus_Tipper"
      diagonal_ZpT.resize(parameter.n_f);

    // for sotring the gradient of the misfit term of respective frequency
    // i.e. \nabla \phi_d = \nabla \phi_d_1 + ... + \nabla \phi_d_(n_f)
    std::vector< EM::VectorXC > g_d_respect;
    g_d_respect.resize(parameter.n_f);

    // loop frequencies
    #ifdef _OPENMP
      #pragma omp parallel for
    #endif
    for(int i = 0; i < parameter.n_f; i++){
      #ifdef _OPENMP
      ///////////////////////////////////////////////////////////////////// ///
      //!!! added to solve the multi-threads closed problem of MKL (Pardiso) //
      //                under OpenMP parallel region.                        //
      /////////////////////////////////////////////////////////////////////////
        // set the number of threads for the nested mkl parallel region
        //mkl_set_num_threads_local(12); //the same as omp_set_num_threads();
        // disable adjustment of the number of threads
        mkl_set_dynamic(0); 
        // has been deprecated by the latest compiler, has the same effect as 
        // environment variable OMP_NESTED=true,
        //omp_set_nested(1);
        // used to control MKL PARDISO to spawn only the first level of 
        // nested threads and avoid significant overhead 
        // If the number of active levels requested exceeds the number of active
        // levels of parallelism supported by the implementation, the value of the 
        // max-active-levels-var ICV will be set to the number of active levels 
        // supported by the implementation.
        omp_set_max_active_levels(2);
      #endif
      // linear source polarization TE
      BC bc_te(parameter.EP, parameter.f[i], parameter.theta, EM::TE); 
      // linear source polarization TM
      BC bc_tm(parameter.EP, parameter.f[i], parameter.theta, EM::TM); 
      // FEM object
      std::stringstream counter;
      counter << parameter.starting_counter;
      std::string s = counter.str();
      std::string model_name = parameter.starting_model[0]+std::string(".")+s;
      Mesh3D mesh(model_name); 
      mesh.prepare_mesh_by_reading_files();
      FEM fem(mesh, parameter.f[i], bc_te, bc_tm, inv_approach,
              parameter.marker, parameter.u, parameter.v, m, abn);
      fem.solve();  

      // true represents we need to calculate the global interpolator operators
      PostProcess postp(fem, parameter.marker, parameter.mu, data_type, true);
      // return sites only once for all the frequencies,
      // as all the frequencies share the same sites
      if(i == 0)
        postp.Return_Sites(sites);

      if(data_type == "Full_Impedance")
      // return impedances for each frequency
        postp.Return_Impedances(Z[i]);
      else if(data_type == "Full_Impedance_plus_Tipper")
        postp.Return_Impedances_plus_Tippers(ZpT[i]);
      else if(data_type == "Off_Diagonal_Impedance")
        postp.Return_Diagonal_Impedances(diagonal_Z[i]);
      else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
        postp.Return_Diagonal_Impedances_plus_Tippers(diagonal_ZpT[i]);

      // compute the gradient of the misfit term of respective frequency
      g_d_respect[i].resize( fwd_tet_id_in_mj.size() );
      g_d_respect[i].setZero();

      if(data_type == "Full_Impedance")
        this->compute_gradient_phi_d_resp(g_d_respect[i], fem, postp, i, Z[i], 
                                          d, Wd, fwd_tet_id_in_mj);
      else if(data_type == "Full_Impedance_plus_Tipper")
        this->compute_gradient_phi_d_resp(g_d_respect[i], fem, postp, i, ZpT[i], 
                                          d, Wd, fwd_tet_id_in_mj);    
      else if(data_type == "Off_Diagonal_Impedance") 
        this->compute_gradient_phi_d_resp(g_d_respect[i], fem, postp, i, diagonal_Z[i], 
                                          d, Wd, fwd_tet_id_in_mj); 
      else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
        this->compute_gradient_phi_d_resp(g_d_respect[i], fem, postp, i, diagonal_ZpT[i], 
                                          d, Wd, fwd_tet_id_in_mj); 
    }

    // for checking site number
    unsigned int n_sites = sites.size();
    assert(n_s == n_sites);

    // !!!compute Forward data F[m]
    if(data_type == "Full_Impedance"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            F_m[i*n_s*n_c + j*n_c + k] = std::real( Z[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( Z[i](j,k) );
          }
    }
    else if(data_type == "Full_Impedance_plus_Tipper"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // impedance plus tipper responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( ZpT[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( ZpT[i](j,k) );
          }
    }
    else if(data_type == "Off_Diagonal_Impedance"){
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // off-diagonal impedance responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( diagonal_Z[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( diagonal_Z[i](j,k) );
          }
    }
    else{ // for data type of "Off_Diagonal_Impedance_plus_Tipper"
      for(unsigned int i = 0; i < n_f; i++)
        for(unsigned int j = 0; j < n_s; j++)
          for(unsigned int k = 0; k < n_c; k++){
            // off-diagonal impedance plus tipper responses
            F_m[i*n_s*n_c + j*n_c + k] = std::real( diagonal_ZpT[i](j,k) );
            F_m[i*n_s*n_c + j*n_c + k + half_N] = std::imag( diagonal_ZpT[i](j,k) );
          }
    }
    
    // !!! Compute the gradient of the misfit term \phi_d
    // the number of the tetrahedral elements of inversion mesh
    unsigned int M = fwd_tet_id_in_mj.size();
    assert(grad_phi_d.size() == M);
    grad_phi_d.setZero();
    for(unsigned int i = 0; i < n_f; i++){
      grad_phi_d = grad_phi_d - 2 * (g_d_respect[i]).real();
    }
  }

  return;
}

void Fwd_Sens_Comp::compute_gradient_phi_d_resp(EM::VectorXC& grad_phi_d_resp,
                                                FEM& fem, 
                                                PostProcess& post_p, 
                                                const int& index_f,
                                                const EM::MatrixXC& F_m_resp,
                                                const EM::VectorXD& d,
                                                EM::EigenSparseMatrixXD& Wd,
                    std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj)
{
  // the number of the sites
  unsigned int n_s = this->data_3d.n_sites;
  // for checking site number
  unsigned int n_sites = 0;
  post_p.get_n_sites( n_sites );
  assert(n_s == n_sites);
  // the number of the observing frequencies
  unsigned int n_f = this->data_3d.n_f;
  // the number of the components
  unsigned int n_c = this->data_3d.n_comps; 
  // the number of the half datum
  unsigned int half_N = this->data_3d.N / 2;
  assert( half_N == n_s * n_f * n_c );
  std::string data_type = this->data_3d.data_type;
  // the number of the tetrahedral elements of inversion mesh
  unsigned int M = fwd_tet_id_in_mj.size(); 
  assert(grad_phi_d_resp.size() == M);;
  grad_phi_d_resp.setZero();

  // !!! for calculating g_1 and g_2
  // !! calculate Wd_f
  EM::EigenSparseMatrixXD Wd_f;
  unsigned int i = index_f * n_s * n_c;
  unsigned int p = n_s * n_c;
  // matrix block of size (p, p), starting at (i, i)
  Wd_f = Wd.block(i, i, p, p);
  // !! calculate d_f
  EM::VectorXC d_f; 
  // initialize d_f by using vector block of size p, starting at i 
  // and vector block of size p, starting at i+half_N
  d_f = d.segment(i, p) + EM::II * d.segment(i + half_N, p);
  // !! calculate F_f
  EM::VectorXC F_f;
  F_f.resize(p);
  F_f.setZero();
  for(unsigned int i = 0; i < n_s; i++)
    for(unsigned int j = 0; j < n_c; j++)
      F_f[i * n_c + j] = F_m_resp(i, j);
  // !! calculate V
  // V = (Wd_f)^T * Wd_f * [d_f -F_f]^*
  EM::VectorXC V;
  V.noalias() = (Wd_f.transpose() * Wd_f) * (d_f - F_f).conjugate();
  // !! calculate g_1 and g_2
  EM::VectorXC g_1, g_2;

  if(data_type == "Full_Impedance"){
    g_1.noalias() = (post_p.L_1_Z).transpose() * V;
    g_2.noalias() = (post_p.L_2_Z).transpose() * V;
  }
  else if(data_type == "Full_Impedance_plus_Tipper"){
    g_1.noalias() = (post_p.L_1_ZpT).transpose() * V;
    g_2.noalias() = (post_p.L_2_ZpT).transpose() * V;
  }
  else if(data_type == "Off_Diagonal_Impedance"){
    g_1.noalias() = (post_p.L_1_diagonal_Z).transpose() * V;
    g_2.noalias() = (post_p.L_2_diagonal_Z).transpose() * V;
  }
  else // for data type of "Off_Diagonal_Impedance_plus_Tipper"
  {
    g_1.noalias() = (post_p.L_1_diagonal_ZpT).transpose() * V;
    g_2.noalias() = (post_p.L_2_diagonal_ZpT).transpose() * V;
  }

  // !!! calculate x_1 and x_2
  EM::VectorXC x_1, x_2;
  // A * x_1,2 = g_1,2
  x_1 = (fem.solver).solve( g_1 );
  x_2 = (fem.solver).solve( g_2 );

  // !!! calculate P_1 and P_2
  // P_1 = - \partial A / \partial m * fem.E0
  // P_2 = - \partial A / \partial m * fem.E1
  // the j-th column represents the vector of \partial A / \partial m_j * fem.E0,1
  EM::EigenSparseMatrixXC P_1, P_2;
  fem.comp_partial_A_partial_m(P_1, P_2, fwd_tet_id_in_mj);
  P_1 = - P_1;
  P_2 = - P_2;

  // !!! compute grad_phi_d_resp
  grad_phi_d_resp = P_1.transpose() * x_1 + P_2.transpose() * x_2;

  return;
}
