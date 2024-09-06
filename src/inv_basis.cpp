/******************************************************************************
 *    Function       : Class Inv_Basis contains some basic routines, which    *
 *                   : are used in various inversion schemes                  *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2019                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2019.07.23                                             *
 *    Last revision  : 2023.06.30                                             *
 ******************************************************************************/
#include "inv_basis.h"
// used for generating screen output log
extern std::ofstream screen_output;

Inv_Basis::Inv_Basis(Startup _startup, Data_3D _data_3d,  
                     Init_Prior_Mod_3D _init_prior_mod_3d, Inv_Para _inv_para, 
                     Fwd_Para_3D& _fwd_para_3d, Mesh3D& _MT3D_inv_mesh)
{
  // !!!initialize the variables used only for MT3D
  this->data_3d  = _data_3d;
  this->init_prior_mod_3d = _init_prior_mod_3d;
  // pointer fwd_para_3d points to the global variable fwd_para_3d
  this->fwd_para_3d = &_fwd_para_3d;
  // further initilization of global struct variable fwd_para_3d
  // source info
  (this->fwd_para_3d)->n_f = _data_3d.n_f;
  (this->fwd_para_3d)->f   = _data_3d.f;
  // regional info, the forwar and inverse mesh share the same region,
  // which is only validated for the same and the nested forward, inverse mesh
  (this->fwd_para_3d)->n_regions    = _init_prior_mod_3d.n_regions;
  (this->fwd_para_3d)->region_table = _init_prior_mod_3d.region_table_ini;
  if((this->fwd_para_3d)->n_sites != this->data_3d.n_sites)
  {
    std::cout << "Error, the settings of n_sites in the .data file and .site file"
              << " are different, \nplease check the two files and the"
              << " corresponding poly file\n"; 
    screen_output << "Error, the settings of n_sites in the .data file and .site file"
                 << " are different, \nplease check the two files and the"
                 << " corresponding poly file\n"; 
    exit(1);
  }
  // initialize 3D inversion mesh related pointer variable MT3D_inv_mesh
  this->MT3D_inv_mesh = &_MT3D_inv_mesh;
  assert((*this->MT3D_inv_mesh).n_elems() > 0);

  // !!!initialize the commonly used variables
  this->startup  = _startup;
  this->inv_para = _inv_para;
  // initialize array abn
  this->abn[0] = this->inv_para.a;
  this->abn[1] = this->inv_para.b;
  this->abn[2] = this->inv_para.n;
  // initialize variable tol_times
  if(this->inv_para.adap_method == 0) // traditional inversion
    this->tol_times = this->inv_para.tol_iter_times;
  else // adaptive inversion
    this->tol_times = this->inv_para.max_iter_times;
  // initialize lambda
  this->lambda = this->inv_para.lambda0;

  // !!!initialize the variable used for 'data' processing (involed in computing)
  // note! we only invert the conductivity in inversion domain, so this->M only represents
  // the number of tetrahedral elements in the inversion domain.
  this->M = _init_prior_mod_3d.m_ini.size();
  this->N = _data_3d.N;
  this->d = _data_3d.d;

  // !!!initialize the data weighted matrix Wd
  // define the dimensions of Wd 
  this->Wd.resize(this->N, this->N);
  // estimated number of the non-zero entries of the sparse matrix
  unsigned int Wd_estimated_non_zero_entries = this->N;
  // a list of triplets, which store the info of all the non-zero entries
  std::vector< EM::T > Wd_tripletList;
  // reserve space for vector
  Wd_tripletList.reserve(Wd_estimated_non_zero_entries);
  // for judging the data is with noise or not
  double tot_std_dev = 0.0;
  // !!initialize Wd_tripletList vector
  for(unsigned int i = 0; i< this->N; i++)
    tot_std_dev += this->data_3d.std_dev(i);
  // for setting Wd = I for noise-free data
  if ( tot_std_dev < 1e-20){
    // a reminder for user to check the inversion parameter file,
    // which is not related to the following contents
    if(this->inv_para.tol_rms > 0.5){
      std::cout << "Error, please set the tol of RMS < 0.5 for noise-free data"
                << " in .inv file\n";
      screen_output << "Error, please set the tol of RMS < 0.5 for noise-free data"
                    << " in .inv file\n";
      std::abort();
    }
    // for setting Wd = I for noise-free data
    for(unsigned int i = 0; i < this->N; i++)
      for(unsigned int j = 0; j < this->N; j++)
        if(i == j)
        Wd_tripletList.push_back( EM::T(i, j, 1.0) );
  }
  else {
    // a reminder for user to check the inversion parameter file,
    // which is not related to the following contents
    if(this->inv_para.tol_rms < 0.5){
      std::cout << "Error, please set the tol of RMS >= 0.5 for noisy data"
                << " in .inv file\n";
      screen_output << "Error, please set the tol of RMS >= 0.5 for noisy data"
                    << " in .inv file\n";
      std::abort();
    }
    // for seting Wd = diag{1/sigma1, 1/sigma2, ..., 1/sigmaN] for noisy data
    // where, sigmaj are the associated standard deviations
    double temp = 0.0;
    for(unsigned int i = 0; i < this->N; i++)
      for(unsigned int j = 0; j < this->N; j++) {
        if(i == j){
          // a reminder for user to check the data file,
          // when the standard deviation is too small
          if( this->data_3d.std_dev(i) < 1.0e-10){
            std::cout << "Warning, the standard deviation of the " << i
            << "-th data < 1.0e-10, please check the data in the data file!\n";
            screen_output << "Warning, the standard deviation of the " << i
            << "-th data < 1.0e-10, please check the data in the data file!\n";
            assert(this->data_3d.std_dev(i) > 0.);
         }
         temp = 1.0 / this->data_3d.std_dev(i);
         Wd_tripletList.push_back( EM::T(i, j, temp));
        }
      }
  }
  // !!initilize Wd by using Wd_triplet vector
  Wd.setFromTriplets(Wd_tripletList.begin(), Wd_tripletList.end());
  Wd_tripletList.clear();

  // initialize prior model covariance matrix Cm_inv
  std::cout << "The number of tetrahedral elements in inversion domain is: "
            << this->M << std::endl;
  screen_output << "The number of tetrahedral elements in inversion domain is: "
                << this->M << std::endl;
  // define the dimensions of matrix Cm_inv 
  this->Cm_inv.resize(this->M, this->M);
  this->init_inv_mod_covar_matrix(this->Cm_inv);
  this->m_ref.resize(this->M);
  this->m_k.resize(this->M);
  for(unsigned int i = 0; i < this->M; i++){
    // initialize m_ref
    double nomin = _init_prior_mod_3d.m_ref[i] - this->inv_para.a;
    double denomin = this->inv_para.b - _init_prior_mod_3d.m_ref[i];
    if( !((nomin > 0) && (denomin > 0)) ){
      std::cout << "Error, reference conductivities should be (set) in the range of\n"
                << "(inv_para.a, inv_para.b) shown in .inv file!\n";
      screen_output << "Error, reference conductivities should be (set) in the range of\n"
                    << "(inv_para.a, inv_para.b) shown in .inv file!\n";
      std::abort(); 
    }
    // parametrized in logarithmic space
    this->m_ref(i) = 1.0 / this->inv_para.n * std::log(nomin / denomin);

    // initialize m_k
    nomin = _init_prior_mod_3d.m_ini[i] - this->inv_para.a;
    denomin = this->inv_para.b - _init_prior_mod_3d.m_ini[i];
    if( !((nomin > 0) && (denomin > 0)) ){
      std::cout << "Error, initial conductivities should be (set) in the range of\n"
                << "(inv_para.a, inv_para.b) shown in .inv file!\n";
      screen_output << "Error, initial conductivities should be (set) in the range of\n"
                    << "(inv_para.a, inv_para.b) shown in .inv file!\n";
      std::abort(); 
    }
    // parametrized in logarithmic space, m_k = m_ini (initial model)
    this->m_k(i) = 1.0 / this->inv_para.n * std::log(nomin / denomin);
  }

  // !!!initialize dynamic marices (including vector) to zero matrices
  this->m_kp1.resize(this->M);
  this->m_kp1.setZero();
  this->F_mk.resize(this->N);
  this->F_mk.setZero();
  this->data_misfit.reserve(this->tol_times + 2);
  this->data_misfit.clear();
  this->RMS.reserve(this->tol_times + 2);
  this->RMS.clear();
  this->roughness.reserve( this->tol_times + 2 );
  this->roughness.clear();
  this->lambda_history.reserve(this->tol_times + 2);
  this->lambda_history.clear();
}

void Inv_Basis::init_inv_mod_covar_matrix(EM::EigenSparseMatrixXD &Cm_inv)
{
  /* compute the inverse of model covariance matrix for 3D MT problem */
  std::cout << "CRT3DMT is computing the inverse of model covariance matrix\n"
            << "for 3D MT problem ..." << std::endl;
  screen_output << "CRT3DMT is computing the inverse of model covariance matrix\n"
                << "for 3D MT problem ..." << std::endl;
  assert( (Cm_inv.rows()== this->M) && (Cm_inv.cols()== this->M));
  // for timekeeping
  Timer time;
  time.start();
  double elapsedtime = 0.0;
  unsigned int Cm_inv_estimated_non_zero_entries = 20 * this-> M;
  // !!!preparing for initializing of sparse matrix
  // a list of triplets, which store the info of all the non-zero entries
  std::vector< EM::T > Cm_inv_tripletList;
  // reserve space for Cm_inv_tripletList vector
  Cm_inv_tripletList.reserve(Cm_inv_estimated_non_zero_entries);

  // !!! map the node id to the index vector of its enbracing tetrahedra,
  //     where the index means the local id of tet in all tets of inversion domain
  //     the number of the nodes in the inversion domain
  unsigned int n_sub_nodes = (this->init_prior_mod_3d.sub_node_id).size();
  // map variable definization, used for mapping the inversion-domain node global
  // id to the index vector of its enbracing tetrahedra
  std::map< unsigned int, std::vector < unsigned int > > nod_id_to_tet_indexs;
  // take alternative shorter name
  std::set< unsigned int >& sub_n_id = this->init_prior_mod_3d.sub_node_id;
  std::set< unsigned int >::iterator it;
  for(it = sub_n_id.begin(); it != sub_n_id.end(); it++ ){
    // get node global id
    unsigned int it_node_id = *it;
    for(unsigned int j = 0; j < this->M; j++){
      // get the global id of the j-th tet in inversion domain
      unsigned int mj_id = this->init_prior_mod_3d.m_ini_id[j];
      // define the object of the mj_id-th tetrahedral element 
      Tet* tet_j = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(mj_id)));
      for(unsigned int k = 0; k < 4; k++){
        unsigned int nk_id = tet_j->get_node_id(k);
        // judging whether node it_node is one vertex of tet j
        if(nk_id == it_node_id)
          nod_id_to_tet_indexs[it_node_id].push_back( j );
      }
    }
  }

  // !!! map the inversion-domain tet global id to the index set of its surrounding 
  //     tetrahedra by using nod_id_to_tet_indexs map. The same as before,
  //     the index means the local id of the tet in all the tets in inversion domain
  for(unsigned int i = 0; i < this->M; i++){
    // get the global id of the i-th tet in inversion domain
    unsigned int mi_id = this->init_prior_mod_3d.m_ini_id[i];
    // define the object of the mi_id-th tetrahedral element 
    Tet* tet_i = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(mi_id)));
    for(unsigned int j = 0; j < 4; j++){
      unsigned int nj_id = tet_i->get_node_id(j);
      std::vector< unsigned int > tet_indexs = nod_id_to_tet_indexs[nj_id];
      for(unsigned int k = 0; k < tet_indexs.size(); k++){
        unsigned int tet_index = tet_indexs[k];
        // exclude itself
        if(tet_index != i)
          this->tet_id_to_tet_indexs[mi_id].insert( tet_index );
      }
    }
  }

  // !!! calculate all the componenets of matrix Cm_inv
  for(unsigned int i = 0; i < this->M; i++){
    // get the global id of the i-th tet in inversion domain
    unsigned int mi_id = this->init_prior_mod_3d.m_ini_id[i];
    // define the object of the mi_id-th tetrahedral element  
    Tet* tet_i = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(mi_id)));
    // get the volume of the tetrahedral element i in inversion domain
    double vi = tet_i->get_size();
    // get the centroid(center of gravity) of inversion-domain tetrahedron i
    Point ri = tet_i->get_gpoint();

    // !! calculate the relevant info of the i-th tetrahedral element
    //    (the main purpose is to calculate sum_v_i)
    // the sum of the volumes of the tetrahedron which share at least one
    // vertex with element i (referred to 'those tetrahedron' in the following)
    double sum_v_i = 0.0;
    // for storing respective volume of those tetrahedra in inversion domain
    std::vector < double > vol_i;
    // for storing the distances from gpoint of tetrahedron i to that of
    // its surrounding tetrahedron
    std::vector < double > distance_i;
    std::set< unsigned int > tet_indexs_i = this->tet_id_to_tet_indexs[mi_id];
    std::set< unsigned int >::iterator it_i;
    for(it_i = tet_indexs_i.begin(); it_i != tet_indexs_i.end(); it_i++){
      // get the global id of the *it_i-th tetrahedral element in inversion domain
      unsigned int tet_id = this->init_prior_mod_3d.m_ini_id[*it_i];
      // define the object of the tet_id-th tetrahedral element
      Tet* tet = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(tet_id)));
      // get the volume of tet_id-th tet
      double v_tet = tet->get_size();
      // get the centroid of tet_id-th tet
      Point gpoint = tet->get_gpoint();
      // calcualte the distance from gpoint to ri 
      double distance = (ri - gpoint).size();
      // assignment of sum_v_i, vol_i and distance_i
      sum_v_i += v_tet;
      vol_i.push_back(v_tet);
      distance_i.push_back(distance);
    }

    // !! Cm_inv calculating
    // another index of the surrongding tetrahedron of tet i
    unsigned int j = 0;
    for(it_i = tet_indexs_i.begin(); it_i != tet_indexs_i.end(); it_i++){
      // get the index, i.e., local id, of the tetrahedral element in inversion domain
      unsigned int tet_index = *it_i;
      assert( (this->M - tet_index) >=0 );
      double wj = vol_i[j] / sum_v_i;
      double temp = (vi * wj) / (distance_i[j] * distance_i[j]);
      Cm_inv_tripletList.push_back( EM::T(i, i, temp) );
      Cm_inv_tripletList.push_back( EM::T(tet_index, tet_index, temp) );
      assert(i != tet_index);
      Cm_inv_tripletList.push_back( EM::T(tet_index, i, -1.0*temp) );
      Cm_inv_tripletList.push_back( EM::T(i, tet_index, -1.0*temp) );
      j++;
    }
  }
  // initilize Cm_inv by using Cm_inv_triplet vector
  Cm_inv.setFromTriplets(Cm_inv_tripletList.begin(), Cm_inv_tripletList.end());
  Cm_inv_tripletList.clear();

  time.stop();
  elapsedtime = time.getElapsedTime();
  std::cout << "Computing the inverse of the model covariance matrix took " 
            << elapsedtime << " s" << '\n';
  screen_output << "Computing the inverse of the model covariance matrix took " 
                << elapsedtime << " s" << '\n';

  return;
}

void Inv_Basis::comp_F_mk()
{
  // vector definition for storing total forward data of model mk
  EM::VectorXD F_data; 
  if(this->startup.data_method == "MT1D"){
    std::cout << "The module is under construction" << std::endl;
    exit(1);
  }
  else if (this->startup.data_method == "MT2D"){
    std::cout << "The module is under construction" << std::endl;
    exit(1);
  }
  else{ // MT3D
    std::cout << "CRT3DMT is computing F[m_k] ...\n";
    screen_output << "CRT3DMT is computing F[m_k] ...\n";
    Timer time;
    double elapsedtime = 0.0;
    time.start();
    Fwd_Sens_Comp fwd_sens_comp(*(this->fwd_para_3d), this->data_3d);
    fwd_sens_comp.Comp_MT3d_Fwd_Responses(this->F_mk, this->m_k,  
                                          this->startup.inv_algorithm, this->abn);
    time.stop();
    elapsedtime = time.getElapsedTime();
    std::cout << "It took " << elapsedtime << " s" << std::endl;
    screen_output << "It took " << elapsedtime << " s" << std::endl;
  }

  return;
}

void Inv_Basis::compute_misfit_RMS_roughness(const double& phi, const EM::VectorXD& m,
                                             const double& lambda)
{
  double temp_roughness = 0.0;
  if(this->startup.data_method == "MT1D")
    std::cout << "The MT1D_sen_eq module is under construction" << std::endl;
  else if(this->startup.data_method == "MT2D")
    std::cout << "The MT2D_sen_eq module is under construction" << std::endl;
  else // for MT3D
    temp_roughness = (m - this->m_ref).transpose() * this->Cm_inv * (m - this->m_ref);           
  this->roughness.push_back(temp_roughness);

  double temp_data_misfit = phi - lambda * temp_roughness;
  this->data_misfit.push_back(temp_data_misfit);

  double temp_RMS = std::sqrt(temp_data_misfit / this->N);
  this->RMS.push_back(temp_RMS);

  return;
}

void Inv_Basis::output_data_and_residual(unsigned int n_adjust, unsigned int k)
{
  // the formal parameters n_adjust and k represent the counter of the times of
  // inversion mesh adjustment, and the optimal iteration number at the present
  // mesh, respectively
  double residual = 0.0;
  std::stringstream number_n, number_k;
  std::string n, s, dat_filename, res_filename;
  number_n << (n_adjust + 1);
  n = number_n.str();
  number_k << k;
  s = number_k.str();
  if(this->inv_para.adap_method == 0)
  {
    dat_filename = this->startup.proj_name
                    + std::string(".") + s + std::string(".dat");
    res_filename = this->startup.proj_name
                    + std::string(".") + s + std::string(".res");
  }
  else
  {
    dat_filename = this->startup.proj_name + std::string(".") + n
                    + std::string(".") + s + std::string(".dat");
    res_filename = this->startup.proj_name + std::string(".") + n
                    + std::string(".") + s + std::string(".res");
  }
  std::ofstream output_data(dat_filename.c_str());
  std::ofstream output_res(res_filename.c_str());

  if(this->startup.data_method == "MT1D")
    std::cout << "The module is under construction." << std::endl;
  else if (this->startup.data_method == "MT2D")
    std::cout << "The module is under construction." << std::endl;
  else if (this->startup.data_method == "MT3D"){
    unsigned int n_s    = this->data_3d.n_sites;
    unsigned int n_f    = this->data_3d.n_f;
    unsigned int n_c    = this->data_3d.n_comps;
    unsigned int half_N = this->data_3d.N / 2;
   
    // ! output header part
    output_data << "# 3D MT predicted data produced by FEM \n"
                << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m)"
                << " Component Real Imag Error \n"
                << "> " << this->data_3d.data_type << " \n"
                << "> exp(+i_omega_t) \n"
                << "> [V/m]/[A/m] \n"
                << "> " << this->data_3d.orient << "\n"
                << "> " << this->data_3d.lat_origin 
                << "  " << this->data_3d.lon_origin << "\n"
                << "> " << this->data_3d.n_f << "  " 
                << this->data_3d.n_sites << '\n';
    output_res << "# 3D MT residual data \n"
               << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m)"
               << " Component Real Imag Error \n"
               << "> " << this->data_3d.data_type << " \n"
               << "> exp(+i_omega_t) \n"
               << "> [V/m]/[A/m] \n"
               << "> " << this->data_3d.orient << "\n"
               << "> " << this->data_3d.lat_origin 
               << "  " << this->data_3d.lon_origin << "\n"
               << "> " << this->data_3d.n_f << "  " 
               << this->data_3d.n_sites << '\n';

    // ! output data part
    if( (this->data_3d.data_type == "Full_Impedance") ||
        (this->data_3d.data_type == "Off_Diagonal_Impedance") ){
      for(unsigned int i = 0; i < n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c; k++){
            unsigned int index = j * n_s * n_c + i * n_c + k;
            unsigned int other_index = i * n_f * n_c + j * n_c + k;
            // output the predicted data
            output_data << setiosflags(std::ios::scientific)
                        << 1.0 / this->data_3d.f[j] << '\t' 
                        // Code 
                        << this->data_3d.site_code[other_index] << '\t'
                        // CG_Lat CG_Lon X(m) Y(m) Z(m)
                        << resetiosflags(std::ios::scientific)
                        << std::setprecision(8)
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.lat[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.lon[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.x[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.y[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.z[other_index] << '\t'
                        // component
                        << this->data_3d.data_comp[other_index] << '\t'
                        << setiosflags(std::ios::scientific)
                        << std::setprecision(6)
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // real part
                        << this->F_mk(index) << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // imaginary part
                        << -1.0 * this->F_mk(index + half_N) << '\t'
                        // error
                        << 2.0e+15 << '\n';
            // output residual between input data and the predicted data 
            output_res  << setiosflags(std::ios::scientific)
                        << 1.0 / this->data_3d.f[j] << '\t'
                        // Code 
                        << this->data_3d.site_code[other_index] << '\t'
                        // CG_Lat CG_Lon X(m) Y(m) Z(m)
                        << resetiosflags(std::ios::scientific)
                        << std::setprecision(8)
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.lat[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.lon[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.x[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.y[other_index] << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(10)
                        << this->data_3d.z[other_index] << '\t'
                        // component
                        << this->data_3d.data_comp[other_index] << '\t'
                        << setiosflags(std::ios::scientific)
                        << std::setprecision(6)
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // real part
                        << std::abs( this->F_mk(index) - this->d(index) ) << '\t'
                        << setiosflags(std::ios::right)
                        << std::setw(15)
                        // imaginary part
                        << std::abs( this->d(index + half_N) - this->F_mk(index + half_N) )
                        // error
                        << '\t' << this->data_3d.std_dev(index) << '\n';
            assert( std::abs(this->data_3d.std_dev(index) - 
                             this->data_3d.std_dev(index+half_N)) < 1.0e-7);
          }
    }
    else // for data type of "Full_Impedance_plus_Tipper"
    {    // and "Off_Diagonal_Impedance_plus_Tipper"
      // ! store data into respective vectors
      // real and imaginary part of Z and T predicted data
      std::vector< double > real_d_Z, real_d_T, imag_d_Z, imag_d_T;
      // real and imaginary part of the residual of Z and T predicted data
      std::vector< double > real_r_Z, real_r_T, imag_r_Z, imag_r_T;
      // error and period
      std::vector< double > err_Z, err_T, period_Z, period_T;
      for(unsigned int i = 0; i < n_s; i++)
        for(unsigned int j = 0; j < n_f; j++)
          for(unsigned int k = 0; k < n_c; k++)
          {
            unsigned int index = j * n_s * n_c + i * n_c + k;
            unsigned int other_index = i * n_f * n_c + j * n_c + k;
            unsigned remainder = other_index % n_c;
            assert( std::abs(this->data_3d.std_dev(index) - 
                             this->data_3d.std_dev(index+half_N)) < 1.0e-7);
            // n_c - 3 = 3 for data type of "Full_Impedance_plus_Tipper"
            // n_c - 3 = 1 for data type of "Off_Diagonal_Impedance_plus_Tipper"
            if(remainder <= (n_c -3))
            {
              // predicted data
              real_d_Z.push_back( this->F_mk(index) );
              imag_d_Z.push_back( this->F_mk(index + half_N) );
              // residal data
              real_r_Z.push_back( this->F_mk(index) - this->d(index) );
              imag_r_Z.push_back( this->F_mk(index + half_N) - 
                                  this->d(index + half_N) );
              // error
              err_Z.push_back( this->data_3d.std_dev(index) );
              // period
              period_Z.push_back(1.0 / this->data_3d.f[j]);
            }
            else
            {
              // predicted data
              real_d_T.push_back( this->F_mk(index) );
              imag_d_T.push_back( this->F_mk(index + half_N) );
              // residal data
              real_r_T.push_back( this->F_mk(index) - this->d(index) );
              imag_r_T.push_back( this->F_mk(index + half_N) - 
                                  this->d(index + half_N) );
              // error
              err_T.push_back( this->data_3d.std_dev(index) );
              // period
              period_T.push_back(1.0 / this->data_3d.f[j]);
            }
          }
      // ! merge Z and T related vectors into one vector
      std::vector< double > real_d, imag_d, real_r, imag_r, err, period;
      for(unsigned int i = 0; i < period_Z.size(); i++){
        real_d.push_back( real_d_Z[i] );
        imag_d.push_back( imag_d_Z[i] );
        real_r.push_back( real_r_Z[i] );
        imag_r.push_back( imag_r_Z[i] );
        err.push_back( err_Z[i] );
        period.push_back( period_Z[i] );
      }
      for(unsigned int i = 0; i < period_T.size(); i++){
        real_d.push_back( real_d_T[i] );
        imag_d.push_back( imag_d_T[i] );
        real_r.push_back( real_r_T[i] );
        imag_r.push_back( imag_r_T[i] );
        err.push_back( err_T[i] );
        period.push_back( period_T[i] );
      }      

      // ! output data part
      for(unsigned int i = 0; i < half_N; i++){
        // output the predicted data
        output_data << setiosflags(std::ios::scientific)
                    << period[i] << '\t'
                    << this->data_3d.site_code[i] << '\t'
                    // CG_Lat CG_Lon X(m) Y(m) Z(m)
                    << resetiosflags(std::ios::scientific)
                    << std::setprecision(8)
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.lat[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.lon[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.x[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.y[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.z[i] << '\t'
                    // component
                    << this->data_3d.data_comp[i] << '\t'
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
                    // error
                    << 2.0e+15 << '\n';
        // output residual between input data and the predicted data 
        output_res  << setiosflags(std::ios::scientific)
                    << period[i] << '\t'
                    << this->data_3d.site_code[i] << '\t'
                    // CG_Lat CG_Lon X(m) Y(m) Z(m)
                    << resetiosflags(std::ios::scientific)
                    << std::setprecision(8)
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.lat[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.lon[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.x[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.y[i] << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(10)
                    << this->data_3d.z[i] << '\t'
                    // component
                    << this->data_3d.data_comp[i] << '\t'
                    << setiosflags(std::ios::scientific)
                    << std::setprecision(6)
                    << setiosflags(std::ios::right)
                    << std::setw(15)
                    // real part
                    << std::abs(real_r[i]) << '\t'
                    << setiosflags(std::ios::right)
                    << std::setw(15)
                    // imaginary part
                    << std::abs(imag_r[i]) << '\t'
                    // error
                    << err[i] << '\n';
      }
    }
  }
  else {
    std::cout << "Input error, please check the data method"
              << " given in the startup file!\n";
    screen_output << "Input error, please check the data method"
                  << " given in the startup file!\n";
    std::exit(1);
  }
  output_data.close();
  output_res.close();

  return;
}

void Inv_Basis::Output_Inv_Mod(unsigned int n_adjust, unsigned int k)
{
  // the formal parameters n_adjust and k represent the counter of the times of
  // inversion mesh adjustment, and the optimal iteration number at the present
  // mesh, respectively
  // changing type int to type string of iteration number n
  std::stringstream number_n, number_k;
  std::string n, s, filename;
  number_n << (n_adjust + 1);
  n = number_n.str();
  number_k << k;
  s = number_k.str();

  if (this->startup.data_method == "MT1D")
    std::cout << "The module is under construction" << std::endl;
  else if (this->startup.data_method == "MT2D")
    std::cout << "The module is under construction" << std::endl;
  else if (this->startup.data_method == "MT3D"){
    if(this->inv_para.adap_method == 0)
      filename = this->startup.proj_name
                + std::string(".") + s + std::string(".rho");
    else
      filename = this->startup.proj_name + std::string(".") + n
                + std::string(".") + s + std::string(".rho");
    // vector definition and initilization for storing inverted resistivity model 
    std::vector < double > inv_res_mod;
    unsigned int n_elements = (*this->MT3D_inv_mesh).n_elems();
    inv_res_mod.resize( n_elements );
    typedef std::map<int, std::vector<double> >::iterator IT;
    for(unsigned int e = 0; e < n_elements; e++){
      // get e-th tetrahedral element 
      Tet* tet= static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(e)));
      assert(e == tet->get_id());
      const unsigned int tet_marker = tet->get_tet_marker();
      if( (tet_marker == 9999999) || (tet_marker == 6666666) ){
        // tets in air-space
        std::vector<double>& ini_elec_para = 
        this->init_prior_mod_3d.region_table_ini[tet_marker];
        // resistivity of air
        inv_res_mod[e] = 1.0 / ini_elec_para[0];
      }
      else{ 
        // tets in inversion domain （having free parameters)
        // index of the elemental parameter in m_k vector 
        // equals to （its marker - 1）
        double nomin = 0, denomin = 0., cond = 0.;
        nomin = this->inv_para.a + ( this->inv_para.b 
                * std::exp( this->inv_para.n * this->m_k( tet_marker -1 ) ) );
        denomin = 1.0 + std::exp( this->inv_para.n * this->m_k( tet_marker -1 ) );
        // conductivity in linear space
        cond = nomin / denomin;
          if( std::isnan( cond ) )
            // we set sigma = b (upper bound), if it is equals to NAN
            cond =  this->inv_para.b;
        // resistivity of individual tet elements in inversion domain
        inv_res_mod[e] = 1.0 / cond;
      }
    }
    // output inverted resistivity model into .vtk file
    (*this->MT3D_inv_mesh).write_out_vtk(filename, inv_res_mod);
  }
  else {
    std::cout << "Input error, please check the data method"
              << " given in the startup file!\n";
    screen_output << "Input error, please check the data method"
                  << " given in the startup file!\n";
    std::exit(1);
  }

  return;
}


/***************************Only for L-BFGS inversion*************************/
void Inv_Basis::comp_F_m(const EM::VectorXD& m, EM::VectorXD& F_m)
{
  assert(F_m.size() == this->N);
  // vector definition for storing total forward data of model m
  EM::VectorXD F_data; 
  if(this->startup.data_method == "MT1D"){
    std::cout << "The MT1D module is under construction" << std::endl;
    exit(0);
  }
  else if (this->startup.data_method == "MT2D"){
    std::cout << "The MT2D module is under construction" << std::endl;
    exit(0);
  }
  else{ // MT3D
    std::cout << "CRT3DMT is computing F[m] ...\n";
    screen_output << "CRT3DMT is computing F[m] ...\n";
    Timer time;
    double elapsedtime = 0.0;
    time.start();
    Fwd_Sens_Comp fwd_sens_comp(*(this->fwd_para_3d), this->data_3d);
    fwd_sens_comp.Comp_MT3d_Fwd_Responses(F_m, m, this->startup.inv_algorithm,
                                          this->abn);
    time.stop();
    elapsedtime = time.getElapsedTime();
    std::cout << "It took " << elapsedtime << " s" << std::endl;
    screen_output << "It took " << elapsedtime << " s" << std::endl;
  }

  return;
}

void Inv_Basis::comp_phi_g_m_sen_eq(const EM::VectorXD& m, const double& lambda,
                                    double& phi, EM::VectorXD& grad_phi)
{
  /////////////////////////////////////////////////////////////////////////////
  //  Compute the value of the objective function \phi and the gradient of   //
  //          \phi by using sensitivity equation method  at model m          //
  /////////////////////////////////////////////////////////////////////////////
  // dimension checking of vector grad_phi (formal argument)
  assert(grad_phi.size() == this->M);
  // !!! Define and compute vector F_m and dev which will be used for 
  //     calculating of phi and grad_phi
  // F_m is an alias of this->F_mk, which will be used for outputting
  EM::VectorXD& F_m = this->F_mk;
  F_m.resize(this->N);
  F_m.setZero();
  // define deviation vector at m, dev = d - F_m;
  EM::VectorXD dev;
  dev.resize(this->N);
  dev.setZero();
  EM::VectorXD grad_phi_d;
  grad_phi_d.resize(this->M);
  grad_phi_d.setZero();

  // !!! compute F[m] and \nabla \phi_d(m)
  if(this->startup.data_method == "MT1D"){
    std::cout << "The MT1D_sen_eq module is under construction" << std::endl;
  }
  else if (this->startup.data_method == "MT2D"){
    std::cout << "The MT2D_sen_eq module is under construction" << std::endl;
  }
  else{ // MT3D
    std::cout << "CRT3DMT is computing F[m] and g_d[m] ...\n";
    screen_output << "CRT3DMT is computing F[m] and g_d[m] ...\n";
    Timer time;
    double elapsedtime = 0.0;
    time.start();
    Fwd_Sens_Comp fwd_sens_comp(*(this->fwd_para_3d), this->data_3d);
    fwd_sens_comp.comp_MT3d_fwd_grad_phi_d(F_m, grad_phi_d, m, this->d,
                                        this->Wd, this->startup.inv_algorithm, 
                          this->init_prior_mod_3d.fwd_tet_id_in_mj, this->abn);
    time.stop();
    elapsedtime = time.getElapsedTime();
    std::cout << "It took " << elapsedtime << " s" << std::endl;
    screen_output << "It took " << elapsedtime << " s" << std::endl;
  }

  // !!! compute \phi_d(m)
  // compute deviation vector
  dev = this->d - F_m;
  // compute data misfit \phi_d(m)
  double phi_d = dev.transpose() * this->Wd.transpose() * this->Wd * dev;

  assert( !(lambda < 0.) );
  if(lambda > 0.){ 
    // !!! compute the objective function \phi and its gradient
    // model roughness \phi_m
    double phi_m = 0.0;
    if(this->startup.data_method == "MT1D")
      std::cout << "The MT1D_sen_eq module is under construction" << std::endl;
    else if(this->startup.data_method == "MT2D")
      std::cout << "The MT2D_sen_eq module is under construction" << std::endl;
    else{ // for MT3D
      phi_m = (m - this->m_ref).transpose() * this->Cm_inv * (m - this->m_ref);
      grad_phi = grad_phi_d + 2.0 * lambda * this->Cm_inv * (m - this->m_ref);
    }
    phi = phi_d + lambda * phi_m;
  }
  else{ 
    // lambda = 0, meaning we only need to calculate phi_d and grad_phi_d
    // i.e., phi = phi_d and grad_phi = grad_phi_d
    assert(lambda < EM::TOL);
    phi = phi_d;
    grad_phi = grad_phi_d;
  }

  return;
}

void Inv_Basis::line_search_strong_Wolfe(double& fx, EM::VectorXD& x, 
                     EM::VectorXD& grad, double& step, const EM::VectorXD& drt, 
                         const EM::VectorXD& xp, const Line_Search_Para& param,
                         unsigned int& iter)
{
  // Do inexact line search algorithm with strong Wolfe condition and then to 
  // update the value of objective function and it's gradient at the updated x
  /////////////////////////////////////////////////////////////////////////////
  // Description of the function:                                            //
  // A line search algorithm with the strong Wolfe condition.                //
  // Implementation based on:                                                //
  // "Numerical Optimization" 2nd Edition,                                   //
  // Jorge Nocedal Stephen J. Wright,                                        //
  // Chapter 3. Line Search Methods, page 59.                                //
  //                                                                         //
  // Introduction to the input parameters:                                   //
  // \param fx   In: The objective function value at the current point.      //
  //             Out: The function value at the new point.                   //
  // \param x    Out: The new point moved to.                                //
  // \param grad In: The current gradient vector. Out: The gradient at the   //
  //                 new point.                                              //
  // \param step  In: The initial step length. Out: The calculated step length.
  // \param drt   The current moving direction.                              //
  // \param xp    The current point.                                         //
  // \param param  Parameters for the line search procedure                  //
  /////////////////////////////////////////////////////////////////////////////

  // dimension double checking
  assert( (x.size() == this->M) && (grad.size() == this->M) ); 
  assert( (drt.size() == this->M) && (xp.size() == this->M) );
  std::cout << "CRT3DMT is doing line search for choosing step length ...\n";
  screen_output << "CRT3DMT is doing line search for choosing step length ...\n";
  // Check the initial value of step
  if(step <= 0)
    throw std::invalid_argument("'step' must be positive");

  // To make this implementation more similar to the other line search
  // methods, the symbol names from the literature
  // ("Numerical Optimizations") have been changed.
  //
  // Literature | this code
  // -----------|--------
  // alpha      | step
  // phi        | fx
  // phi'       | dg

  // the rate, by which the 
  const double expansion = 2.0;      

  // Save the function value at the current model x
  const double fx_init = fx;
  // Projection of gradient on the search direction, \phi'(0)
  const double dg_init = grad.dot(drt);
  // Make sure d points to a descent direction
  if(dg_init > 0)
    throw std::logic_error("the moving direction increases the objective function value");
    
    const double test_decr = param.ftol * dg_init,    // Sufficient decrease
                 test_curv = -param.wolfe * dg_init;  // Curvature

    // Ends of the line search range (step_lo > step_hi is allowed)
    double step_hi, step_lo = 0,
             fx_hi,   fx_lo = fx_init,
             dg_hi,   dg_lo = dg_init;

    // STEP 1: Bracketing Phase
    //   Find a range guaranteed to contain a step satisfying strong Wolfe.
    //
    //   See also:
    //     "Numerical Optimization", "Algorithm 3.5 (Line Search Algorithm)".
    iter = 0;
    for(;;)
    {
      x = xp + step * drt;

      if(this->startup.data_method == "MT1D"){
        std::cout << "The MT1D module is under construction" 
                  << std::endl;
        exit(1);
      }
      else if(this->startup.data_method == "MT2D"){
        std::cout << "The MT2D module is under construction" 
                  << std::endl;
        exit(1);
      }
      else // for MT3D
        this->comp_phi_g_m_sen_eq(x, this->lambda, fx, grad);

      if(iter++ >= param.max_linesearch)
        return;

      const double dg = grad.dot(drt);

      if( fx - fx_init > step * test_decr || (0 < step_lo && fx >= fx_lo) )
      {
        step_hi = step;
          fx_hi = fx;
          dg_hi = dg;
        break;
      }

      if( std::abs(dg) <= test_curv )
        return;

      step_hi = step_lo;
        fx_hi = fx_lo;   
        dg_hi = dg_lo;
      step_lo = step;
        fx_lo = fx;
        dg_lo = dg;

      if( dg >= 0 )
        break;

      step *= expansion;
    }

    // STEP 2: Zoom Phase
    //   Given a range (step_lo,step_hi) that is guaranteed to
    //   contain a valid strong Wolfe step value, this method
    //   finds such a value.
    //
    //   See also:
    //   "Numerical Optimization", "Algorithm 3.6 (Zoom)".
    for(;;)
    {
      // use {fx_lo, fx_hi, dg_lo} to make a quadric interpolation of
      // the function said interpolation is used to estimate the minimum
      //
      // polynomial: p (x) = c0*(x - step)² + c1
      // conditions: p (step_hi) = fx_hi
      //             p (step_lo) = fx_lo
      //             p'(step_lo) = dg_lo
      step  = (fx_hi-fx_lo)*step_lo - (step_hi*step_hi - step_lo*step_lo)*dg_lo/2;
      step /= (fx_hi-fx_lo)         - (step_hi         - step_lo        )*dg_lo;

      // if interpolation fails, bisection is used
      if( step <= std::min(step_lo,step_hi) ||
          step >= std::max(step_lo,step_hi) )
          step  = step_lo/2 + step_hi/2;

      x = xp + step * drt;

      if(this->startup.data_method == "MT1D"){
        std::cout << "The MT1D module is under construction" 
                  << std::endl;
        exit(1);
      }
      else if(this->startup.data_method == "MT2D"){
        std::cout << "The MT2D module is under construction" 
                  << std::endl;
        exit(1);
      }
      else // for MT3D
        this->comp_phi_g_m_sen_eq(x, this->lambda, fx, grad);

      if(iter++ >= param.max_linesearch)
        return;

      const double dg = grad.dot(drt);

      if( fx - fx_init > step * test_decr || fx >= fx_lo )
      {
        if( step == step_hi )
          throw std::runtime_error("the line search routine failed, possibly due to insufficient numeric precision");

        step_hi = step;
          fx_hi = fx;
          dg_hi = dg;
      }
      else
      {
        if( std::abs(dg) <= test_curv )
          return;

        if( dg * (step_hi - step_lo) >= 0 )
        {
          step_hi = step_lo;
            fx_hi = fx_lo;
            dg_hi = dg_lo;
        }

        if( step == step_lo )
          throw std::runtime_error("the line search routine failed, possibly due to insufficient numeric precision");

        step_lo = step;
          fx_lo =   fx;
          dg_lo =   dg;
      }
    }
  std::cout << "Line search finished for current inversion iteration. \n";
  screen_output << "Line search finished for current inversion iteration. \n";

  return;
}

/***************************For adaptive inversion****************************/
void Inv_Basis::Comp_Abs_Gra_M(EM::VectorXD& abs_gra_m)
{
  abs_gra_m.resize(this->M);

  if(this->startup.data_method == "MT1D"){
    std::cout << "The function Comp_Abs_Gra_M of MT2D is under construction."
              << std::endl;
    screen_output << "The function Comp_Abs_Gra_M of MT2D is under construction."
              << std::endl;
    exit(1);    
  }
  else if(this->startup.data_method == "MT2D"){
    std::cout << "The function Comp_Abs_Gra_M of MT2D is under construction."
              << std::endl;
    screen_output << "The function Comp_Abs_Gra_M of MT2D is under construction."
              << std::endl;
    exit(1);
  }
  else{ // for MT3D
    std::cout << "\nCRT3DMT is computing model-gradient-based indicator used for\n"
              << "refinement of inversion mesh ..." << std::endl;
    screen_output << "\nCRT3DMT is computing model-gradient-based indicator used for\n"
                  << "refinement of inversion mesh ..." << std::endl;

    // calculate the absolute value of gradient m (conductivity in logarithmic space),
    // m is the last L-BFGS results based on the present inversion mesh
    for(unsigned int i = 0; i < this->M; i++){
      // get the global id of the i-th tet in the inversion domain
      unsigned int mi_id = this->init_prior_mod_3d.m_ini_id[i];
        if(mi_id >= (*this->MT3D_inv_mesh).n_elems()){
          std::cout << "tet_id: " << mi_id << "\n n_elems: " 
                    <<  (*this->MT3D_inv_mesh).n_elems() << std::endl;
          screen_output << "tet_id: " << mi_id << "\n n_elems: " 
                        <<  (*this->MT3D_inv_mesh).n_elems() << std::endl;
          assert(mi_id < (*this->MT3D_inv_mesh).n_elems());
        }
      // define the object of the mi_id-th tetrahedral element
      Tet* tet_i = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(mi_id)));
      // get the volume of the tetrahedral element i in inversion domain
      double vi = tet_i->get_size();
      // get the centroid(center of gravity) of inversion-domain tetrahedron i
      Point ri = tet_i->get_gpoint();
      
      // !! calculate the relevant info of the i-th tetrahedral element 
      //    (the main purpose is to calculate sum_v_i)
      // the sum of the volumes of the tetrahedron which share at least one
      // vertex with element i (referred to 'those tetrahedron' in the following)
      double sum_v_i = 0.0;
      // for storing respective volume of those tetrahedra in inversion domain
      std::vector < double > vol_i;
      // for storing the distances from gpoint of tetrahedron i to that of
      // its surrounding tetrahedron
      std::vector < double > distance_i;
      std::set< unsigned int > tet_indexs_i = this->tet_id_to_tet_indexs[mi_id];
      std::set< unsigned int >::iterator it_i;
      for(it_i = tet_indexs_i.begin(); it_i != tet_indexs_i.end(); it_i++){
        // get the global id of the *it_i-th tetrahedral element in inversion domain
        unsigned int tet_id = this->init_prior_mod_3d.m_ini_id[*it_i];
        // define the object of the tet_id-th tetrahedral element
        Tet* tet = static_cast<Tet*>(&((*this->MT3D_inv_mesh).get_elem(tet_id)));
        // get the volume of tet_id-th tet
        double v_tet = tet->get_size();
        // get the centroid of tet_id-th tet
        Point gpoint = tet->get_gpoint();
        // calcualte the distance from gpoint to ri 
        double distance = (ri - gpoint).size();
        // assignment of sum_v_i, vol_i and distance_i
        sum_v_i += v_tet;
        vol_i.push_back(v_tet);
        distance_i.push_back(distance);
      }

      // !! calculate the absolute value of the gradient m
      unsigned int j = 0;
      double temp = 0.0;
      for(it_i = tet_indexs_i.begin(); it_i != tet_indexs_i.end(); it_i++){
        // get the index, i.e., local id, of the tetrahedral element 
        // in the inversion domain
        unsigned int tet_index = *it_i;
        assert(i != tet_index);
        assert( (this->M - tet_index) >=0 );
        double wj = vol_i[j] / sum_v_i;
        double delta_m_ij = this->m_k( i ) - this->m_k( tet_index );
        temp += wj * std::abs( delta_m_ij / distance_i[j] );
        j++;
      }
      abs_gra_m(i) = vi * temp;
      assert( !(abs_gra_m(i) < 0) );
    }
  }

  return;
}

void Inv_Basis::Comp_Abs_Gra_Phi_D(EM::VectorXD& abs_gra_phi_d)
{
  // size definition
  abs_gra_phi_d.resize(this->M);

  if(this->startup.data_method == "MT1D"){
    std::cout << "The function Comp_Abs_Gra_Phi_D of MT1D is under construction."
              << std::endl;
    exit(1);
  }
  else if(this->startup.data_method == "MT2D"){
    std::cout << "The function Comp_Abs_Gra_Phi_D of MT2D is under construction."
              << std::endl;
    exit(1);
  }
  else{ // for MT3D
    std::cout << "CRT3DMT is computing gradient-phi-d-based indicator used for\n"
              << "refinement of inversion mesh ..." << std::endl;
    screen_output << "CRT3DMT is computing gradient-phi-d-based indicator used for\n"
                  << "refinement of inversion mesh ..." << std::endl;
    // just used to call function comp_phi_g_m_sen_eq
    double useless_variable = 0.0;
    this->comp_phi_g_m_sen_eq(this->m_k, 0., useless_variable, abs_gra_phi_d);
  }

  // obtain the absolute value of gra_phi_d
  for(unsigned int i = 0; i < this->M; i++)
    abs_gra_phi_d(i) = std::abs( abs_gra_phi_d(i) );

  return;
}

void Inv_Basis::get_m_last(EM::VectorXD& m_ini)
{
  m_ini = this->m_k;

  return;
}
