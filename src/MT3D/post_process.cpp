/*****************************************************************************/
/*                                                                           */
/*  Copyright 2021                                                           */
/*  Huang Chen and Zhengyong Ren                                             */
/*  chenhuang@cqu.edu.cn and renzhengyong@csu.edu.cn                         */
/*                                                                           */
/*****************************************************************************/
//  Please note: all the symbol x, y, z (global coordinate) used in this     
//  file should be u, v, n (local coordinate), respectively!    

#include <iostream>
#include <fstream>
#include "post_process.h"
#include "fem.h"
#include "struct_defs.h"
// used for generating screen output log
extern std::ofstream screen_output;

// ---------------------------------------------------------------------------
PostProcess::PostProcess(FEM& fem, int marker, 
                         std::vector<Point>& u,
                         std::vector<Point>& v,
	                       std::vector<Real>&  mu,
                         std::ofstream& out_file)
{ 
  this->f_  = fem._f;
  this->post_processing(fem, marker, u, v, mu, out_file);
}

PostProcess::PostProcess(FEM& fem, int marker,
	                       std::vector<Real>&  mu,
                         std::string& data_type,
                         bool calculate_L1_L2_or_not)
{ 
  this->f_ = fem._f;
  this->post_processing_without_output(fem, marker, mu, data_type,
                                       calculate_L1_L2_or_not);
}


// this->_EH to Z, app and phase, and output
void PostProcess::post_processing(FEM& fem, int marker,
                                  std::vector<Point>& u,
                                  std::vector<Point>& v,
				                          std::vector<Real>& mu,
                                  std::ofstream& OUT)
{
  // compute Z, App.Resis., Phases, E, H, and output 
  // which will be used for producing synthetic data
  std::map<unsigned int, std::vector<ComplexPoint> >& EH= fem._EH;
  const Real f       = fem._f;
  const Real omega   = 2.0*EM::PI*f; 

  // !!! store the observing nodes into site vector in 
  //     the order of those ranked in the poly file
  std::vector< Node > sites;
  const unsigned int n_nodes = (fem._mesh3d).n_nodes();
  for(unsigned int i = 0; i < n_nodes; i++){
    Node *a = &(fem._mesh3d).get_node(i);
    if(a->get_marker() == marker)
      sites.push_back(*a); 
  }

  int n_sites = sites.size();
  // mu.size() = parameter.n_sites;
  if(n_sites != mu.size()){
    std::cout << "The input site number is wrong, please check polygon or site file!" 
              << std::endl; 
    std::abort();
  }
  // Z_ij, E_i,j, H_i,j, i,j = x,y
  std::vector< std::vector<Dcomplex> > Z(n_sites), E(n_sites), H(n_sites);
  // HZ and Ez at the earth side, which will be only used for outputting
  std::vector< std::vector<Dcomplex> > HZZ(n_sites), EZZ(n_sites); 
  // rho_ij, phi_ij, i,j = x,y
  std::vector< std::vector<Real> >     App(n_sites), Phase(n_sites);
  // A denotes Tx, and B denotes Ty
  std::vector<Dcomplex>                A(n_sites), B(n_sites); 

  for(unsigned int i=0; i<n_sites; i++) {
    unsigned int s= sites[i].get_id();
    std::map<unsigned int, std::vector<ComplexPoint> >::iterator it= EH.find(s);
    assert(it!=EH.end());
    // Ex,y and Hx,y
    Dcomplex e[2][2], h[2][2]; 
    // Ez and Hz on earth side
    Dcomplex Hz[2], Ez[2];     
    std::vector<ComplexPoint>& value= (*it).second;
    Point n= (u[i].cross(v[i])).unit();
    // Ex1: e[0][0]; Ex2: e[0][1]
    // Ey1: e[1][0]; Ey2: e[1][1]
    // Hx1: h[0][0]; Hx2: h[0][1]
    // Hy1: h[1][0]; Hy2: h[1][1]
    // value[1,3,5,7] are the fields at the earth side, 
    // while value[0,2,4,6] are the fields at the air side
    e[0][0]= value[1]*u[i]; e[1][0]= value[1]*v[i];
    e[0][1]= value[5]*u[i]; e[1][1]= value[5]*v[i];
    h[0][0]= value[3]*u[i]; h[1][0]= value[3]*v[i];
    h[0][1]= value[7]*u[i]; h[1][1]= value[7]*v[i];
    // Hz1: Hz[0]; Hz2: Hz[1]
    // Ez1: Ez[0]; Ez2: Ez[1]
    Hz[0]  = value[3]*n;    Hz[1]  = value[7]*n;
    Ez[0]  = value[1]*n;    Ez[1]  = value[5]*n;

    /*
    // !!! for checking local interpolation operator

    // PLEASE NOTE: when using the following code, we need 
    // to let gx,y_e,h map be calculated in fem.cpp file, 
    // i.e., where this->_inv_algorithm != "Fwd_only" should be 
    // changed into this->_inv_algorithm == "Fwd_only" except
    // the one used to change ex,y,z to u,v,n

    // assumed!!! u = (1, 0, 0); v = (0, 1, 0);
    std::map<unsigned int, EM::EigenSparseVectorXD>::iterator 
    it_x_e = fem.gx_e.find(s);
    std::map<unsigned int, EM::EigenSparseVectorXD>::iterator 
    it_y_e = fem.gy_e.find(s);
    std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
    it_x_h = fem.gx_h.find(s);
    std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
    it_y_h = fem.gy_h.find(s);
    assert(it_x_e != fem.gx_e.end());
    EM::EigenSparseVectorXD& x_e = (*it_x_e).second;
    EM::EigenSparseVectorXD& y_e = (*it_y_e).second;       
    EM::EigenSparseVectorXC& x_h = (*it_x_h).second;
    EM::EigenSparseVectorXC& y_h = (*it_y_h).second;
    e[0][0] = x_e * fem.E0;
    e[0][1] = x_e * fem.E1;
    e[1][0] = y_e * fem.E0;
    e[1][1] = y_e * fem.E1;
    h[0][0] = x_h * fem.E0;
    h[0][1] = x_h * fem.E1;
    h[1][0] = y_h * fem.E0;
    h[1][1] = y_h * fem.E1;
    std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
    it_z_h = fem.gz_h.find(s);
    EM::EigenSparseVectorXC& z_h = (*it_z_h).second;
    std::cout << "deviation Hz1: " << Hz[0] - z_h * fem.E0 << "\t"
              << "deviation Hz2: " << Hz[1] - z_h * fem.E1 << "\n";
    Hz[0] = z_h * fem.E0;
    Hz[1] = z_h * fem.E1;
    */
  
    // Ty
    B[i]= (Hz[0]*h[0][1]-Hz[1]*h[0][0])/(h[1][0]*h[0][1]-h[1][1]*h[0][0]);
    // Tx = (Hz1 -Hy1 * Ty) / Hx1 checked!
    A[i]= (Hz[0]-B[i]*h[1][0])/h[0][0];
    HZZ[i].resize(2); EZZ[i].resize(2);
    HZZ[i][0]= Hz[0]; // Hz1
    HZZ[i][1]= Hz[1]; // Hz2
    EZZ[i][0]= Ez[0]; // Ez1
    EZZ[i][1]= Ez[1]; // Ez2

    // used for calculating of impedance 
    Dcomplex deth= (h[0][0]*h[1][1]-h[0][1]*h[1][0]);
    if(std::isnan(std::abs(deth))){
      std::cout << "abs(deth) = NAN, program stopped!\n";
      exit(0);
    }    
    assert(std::abs(deth) > 0.);
    Dcomplex inv_deth= 1.0/deth;

    // calculate Zxx, Zyx, Zxy and Zyy
    Z[i].resize(4);
    Z[i][0]= inv_deth*( e[0][0]*h[1][1]-e[0][1]*h[1][0]); // Zxx
    Z[i][2]= inv_deth*(-e[0][0]*h[0][1]+e[0][1]*h[0][0]); // Zxy
    Z[i][1]= inv_deth*( e[1][0]*h[1][1]-e[1][1]*h[1][0]); // Zyx
    Z[i][3]= inv_deth*(-e[1][0]*h[0][1]+e[1][1]*h[0][0]); // Zyy
    App[i].resize(4); Phase[i].resize(4);

    Real mu_v=mu[i]*EM::MU0;
 
    for(int j=0; j<4; j++) {
      App[i][j] = (1./(omega*mu_v))*(std::abs(Z[i][j]))*(std::abs(Z[i][j]));
      Phase[i][j]= std::abs(std::atan(std::imag(Z[i][j])/
                   std::real(Z[i][j]))*180.0/EM::PI);
	    }
    E[i].resize(4); H[i].resize(4);
    E[i][0]= e[0][0]; E[i][1]= e[0][1];
    E[i][2]= e[1][0]; E[i][3]= e[1][1];
    H[i][0]= h[0][0]; H[i][1]= h[0][1];
    H[i][2]= h[1][0]; H[i][3]= h[1][1];
  }

  // output-------------------
  // p,x,y,z,marker,E,H,Z,app,phase,AB,2HZ and 2Ez (5+8*3+4*5=49 terms per row)
  for(int n=0; n<n_sites; n++) {
    OUT<<sites[n].get_id()<<"\t"
       <<sites[n].get_xyz()[0]<<"\t"<<sites[n].get_xyz()[1]<<"\t"
       <<sites[n].get_xyz()[2]<<"\t"<<sites[n].get_marker()<<"\t";
    for(int i=0; i<4; i++) 
      OUT<<std::real(E[n][i])<<"\t"<<std::imag(E[n][i])<<"\t";
    for(int i=0; i<4; i++) 
      OUT<<std::real(H[n][i])<<"\t"<<std::imag(H[n][i])<<"\t";
    for(int i=0; i<4; i++) {
      OUT<<std::real(Z[n][i])<<"\t"<<std::imag(Z[n][i])<<"\t";
      //std::cout <<std::real(Z[n][i])<<"\t"<<std::imag(Z[n][i])<<"\t";
      }
    for(int i=0; i<4; i++)  OUT<<App[n][i]<<"\t";
    for(int i=0; i<4; i++)  OUT<<Phase[n][i]<<"\t";

    // print frequency info on screen
    if(n == 0){
      std::cout << "Frequency = " << this->f_ << " Hz:\n";
      screen_output << "Frequency = " << this->f_ << " Hz:\n";
    }
    //std::cout << std::setprecision(20) << App[n][1] << std::endl;
    // rho_yx, checked by comparing our results with those from Ren's paper
    std::cout <<  App[n][1] << std::endl;
    screen_output <<  App[n][1] << std::endl;

    // A,B of VMT (4)
    OUT<<std::real(A[n])<<"\t"<<std::imag(A[n])<<"\t"
       <<std::real(B[n])<<"\t"<<std::imag(B[n])<<"\t";

    // HZ (4)
    OUT<<std::real(HZZ[n][0])<<"\t"<<std::imag(HZZ[n][0])<<"\t"
       <<std::real(HZZ[n][1])<<"\t"<<std::imag(HZZ[n][1])<<"\t";

    // Ez on the earth side (4)
    OUT<<std::real(EZZ[n][0])<<"\t"<<std::imag(EZZ[n][0])<<"\t"
       <<std::real(EZZ[n][1])<<"\t"<<std::imag(EZZ[n][1])<<"\t";

    OUT<<"\n";
  }

  //initialize member variables for producing synthetic data
  this->sites_ = sites;
  this->n_sites_ = n_sites;
  // Z_i,j; i = 0, ..., (n_sites-1); j = 0,...,3, represents Zxx, Zxy, Zyx, Zyy
  this->Z_.resize(n_sites, 4);
  for(unsigned int i = 0; i < n_sites; i++){
    this->Z_(i,0) = Z[i][0]; // Zxx
    this->Z_(i,1) = Z[i][2]; // Zxy
    this->Z_(i,2) = Z[i][1]; // Zyx
    this->Z_(i,3) = Z[i][3]; // Zyy
  }

  // T_i,j; i = 0, ..., (n_sites-1); j = 0, 1, represents Tx, Ty
  this->T_.resize(n_sites, 2);
  for(unsigned int i = 0; i < n_sites; i++){
    this->T_(i, 0) = A[i]; // Tx
    this->T_(i, 1) = B[i]; // Ty
  }

  return;
} // end postprocessing with output

// this->_EH to Z, app and phase (no output)
void PostProcess::post_processing_without_output(FEM& fem, int marker, 
				                                         std::vector<Real>& mu,
                                                 std::string& data_type,
                                                 bool calculate_L1_L2_or_not)
{
  // !!! compute E, H, Z, T (optional) and no both file and screen outputting,
  //     which will be called when doing inversion
  // !! prepare to do computing
  std::vector<Point>& u = fem.u;
  std::vector<Point>& v = fem.v;
  std::map<unsigned int, std::vector<ComplexPoint> >& EH= fem._EH; 
  const Real f       = fem._f;
  const Real omega   = 2.0*EM::PI*f; 
  // !!! store the observing nodes into site vector in 
  //     the order of those ranked in the poly file
  std::vector< Node > sites;
  const unsigned int n_nodes = (fem._mesh3d).n_nodes();
  for(unsigned int i = 0; i < n_nodes; i++){
    if(sites.size() < mu.size() + 1){
      Node *a = &(fem._mesh3d).get_node(i);
      if(a->get_marker() == marker)
      sites.push_back(*a); 
    }
    else
    break;  
  }
  int n_sites = sites.size();
  assert(n_sites == mu.size());
  std::vector<std::vector<Dcomplex> > E(n_sites), H(n_sites), Z(n_sites);
  // A denotes Tx, B denotes Ty
  std::vector<Dcomplex>  A, B;
  unsigned int n_comps;
  if(data_type == "Full_Impedance")
    // the number of the components of the full impedance data
    n_comps = 4;
  else if(data_type == "Full_Impedance_plus_Tipper"){
    // the total number of the components of Z and T(tipper)
    n_comps = 6;
    A.resize(n_sites);
    B.resize(n_sites);
  }
  else if(data_type == "Off_Diagonal_Impedance")
    n_comps = 2;
  else{ // for data type of "Off_Diagonal_Impedance_plus_Tipper"
    n_comps = 4;
    A.resize(n_sites);
    B.resize(n_sites);
  }

  if(calculate_L1_L2_or_not){
    // the number of the unknowns of the forward modelling
    unsigned int n_dofs = (fem._dofs).get_n_dofs();
    unsigned n_rows = n_sites * n_comps;
    unsigned n_cols = n_dofs;
    // estimated number of non-zero entries in each row,
    // which is the same as that used for gx,y_e,h and gz_h in FEM class
    unsigned int n_estimated = 72;
    // resize the sparse matrices to a n_rows x n_cols matrix and initialize 
    // them to zero.
    if(data_type == "Full_Impedance"){
      this->L_1_Z.resize(n_rows, n_cols);
      this->L_1_Z.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
      this->L_2_Z.resize(n_rows, n_cols);
      this->L_2_Z.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
    }
    else if(data_type == "Full_Impedance_plus_Tipper"){
      this->L_1_ZpT.resize(n_rows, n_cols);
      this->L_1_ZpT.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
      this->L_2_ZpT.resize(n_rows, n_cols);
      this->L_2_ZpT.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );    
    }
    else if(data_type == "Off_Diagonal_Impedance"){
      this->L_1_diagonal_Z.resize(n_rows, n_cols);
      this->L_1_diagonal_Z.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
      this->L_2_diagonal_Z.resize(n_rows, n_cols);
      this->L_2_diagonal_Z.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
    }
    else{ // for data type of "Off_Diagonal_Impedance_plus_Tipper"
      this->L_1_diagonal_ZpT.resize(n_rows, n_cols);
      this->L_1_diagonal_ZpT.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
      this->L_2_diagonal_ZpT.resize(n_rows, n_cols);
      this->L_2_diagonal_ZpT.reserve( Eigen::VectorXi::Constant(n_rows, n_estimated) );
    }
  }

  // !! calculate Z, T (tipper, if needed) and the corresponding global  
  //    interpolation operators (if needed)
  for(unsigned int i=0; i<n_sites; i++) {
    unsigned int s = sites[i].get_id();
    std::map<unsigned int, std::vector<ComplexPoint> >::iterator it=
  	EH.find(s);
    assert(it!=EH.end());
    Dcomplex e[2][2], h[2][2]; // E and H on (u,v)
    Dcomplex Hz[2];            // Hz on earth side
    std::vector<ComplexPoint>& value= (*it).second;
    Point n = (u[i].cross(v[i])).unit();
    // Ex1: e[0][0]; Ex2: e[0][1]
    // Ey1: e[1][0]; Ey2: e[1][1]
    // Hx1: h[0][0]; Hx2: h[0][1]
    // Hy1: h[1][0]; Hy2: h[1][1]
    // x components,            y components
    e[0][0]= value[1]*u[i]; e[1][0]= value[1]*v[i];  // TE
    e[0][1]= value[5]*u[i]; e[1][1]= value[5]*v[i];  // TM
    h[0][0]= value[3]*u[i]; h[1][0]= value[3]*v[i];  // TE
    h[0][1]= value[7]*u[i]; h[1][1]= value[7]*v[i];  // TM  
    // used for both impedance and global interpolation operator calculations
    Dcomplex deth = (h[0][0]*h[1][1] - h[0][1]*h[1][0]);
    if(std::isnan(std::abs(deth))){
      std::cout << "abs(deth) = NAN, program stopped!\n";
      exit(0);
    }
    assert(std::abs(deth) > 0.);
    Dcomplex inv_deth= 1.0 / deth;
    // calculate Zxx, Zyx, Zxy and Zyy
    Z[i].resize(4);
    if( (data_type == "Full_Impedance") ||
        (data_type == "Full_Impedance_plus_Tipper") ){
      Z[i][0]= inv_deth*( e[0][0]*h[1][1]-e[0][1]*h[1][0]); // Zxx
      Z[i][2]= inv_deth*(-e[0][0]*h[0][1]+e[0][1]*h[0][0]); // Zxy
      Z[i][1]= inv_deth*( e[1][0]*h[1][1]-e[1][1]*h[1][0]); // Zyx
      Z[i][3]= inv_deth*(-e[1][0]*h[0][1]+e[1][1]*h[0][0]); // Zyy
    }
    else{
      Z[i][2]= inv_deth*(-e[0][0]*h[0][1]+e[0][1]*h[0][0]); // Zxy
      Z[i][1]= inv_deth*( e[1][0]*h[1][1]-e[1][1]*h[1][0]); // Zyx
    }
    if( (data_type == "Full_Impedance_plus_Tipper") ||
        (data_type == "Off_Diagonal_Impedance_plus_Tipper") ){
      // calculate Tipper data
      //     TE (Hz1)                TM (Hz2)
      // here, EH[s][3,7] = value[3,7]
      Hz[0]  = value[3]*n;     Hz[1]  = value[7]*n;
      B[i]= (Hz[0]*h[0][1]-Hz[1]*h[0][0])/(h[1][0]*h[0][1]-h[1][1]*h[0][0]);// Ty
      A[i]= (Hz[0]-B[i]*h[1][0])/h[0][0];// Tx
    }

    // ! calculate global interpolation operator L_1,2_Z or L_1,2_ZpT
    if(calculate_L1_L2_or_not){
      Dcomplex Ex1 = e[0][0]; 
      Dcomplex Ex2 = e[0][1];
      Dcomplex Ey1 = e[1][0];
      Dcomplex Ey2 = e[1][1];
      Dcomplex Hx1 = h[0][0];
      Dcomplex Hx2 = h[0][1];
      Dcomplex Hy1 = h[1][0];
      Dcomplex Hy2 = h[1][1];
      Dcomplex Hz1 = Hz[0];
      Dcomplex Hz2 = Hz[1];
      EM::EigenSparseVectorXC temp1, temp3;
      Dcomplex temp2 = 0.0; 
      temp1.setZero(); temp3.setZero();
      std::map<unsigned int, EM::EigenSparseVectorXD>::iterator 
      it_x_e = (fem.gx_e).find(s);
      std::map<unsigned int, EM::EigenSparseVectorXD>::iterator 
      it_y_e = (fem.gy_e).find(s);
      std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
      it_x_h = (fem.gx_h).find(s);
      std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
      it_y_h = (fem.gy_h).find(s);
      std::map<unsigned int, EM::EigenSparseVectorXC>::iterator 
      it_z_h = (fem.gz_h).find(s);
      assert(it_x_e != (fem.gx_e).end());
      EM::EigenSparseVectorXD& x_e = (*it_x_e).second;
      EM::EigenSparseVectorXD& y_e = (*it_y_e).second;       
      EM::EigenSparseVectorXC& x_h = (*it_x_h).second;
      EM::EigenSparseVectorXC& y_h = (*it_y_h).second;      
      EM::EigenSparseVectorXC& z_h = (*it_z_h).second;

      if(data_type == "Full_Impedance"){
        // ! calculate L_1_Z and L_2_Z
        // gxx_1
        temp1 = -Hy2 * x_e + Ex2 * y_h;
        temp2 =  Ex1 * Hy2 - Ex2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_Z.row(i * n_comps + 0) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxx_2
        temp1 = -Ex1 * y_h + Hy1 * x_e;
        temp2 =  Ex1 * Hy2 - Ex2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_Z.row(i * n_comps + 0) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gxy_1
        temp1 = -Ex2 * x_h + Hx2 * x_e;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_Z.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxy_2
        temp1 = -Hx1 * x_e + Ex1 * x_h;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_Z.row(i * n_comps + 1) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyx_1
        temp1 = -Hy2 * y_e + Ey2 * y_h;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_Z.row(i * n_comps + 2) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyx_2
        temp1 = -Ey1 * y_h + Hy1 * y_e;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_Z.row(i * n_comps + 2) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyy_1 
        temp1 = -Ey2 * x_h + Hx2 * y_e;
        temp2 =  Ey2 * Hx1 - Ey1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_Z.row(i * n_comps + 3) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyy_2
        temp1 = -Hx1 * y_e + Ey1 * x_h;
        temp2 =  Ey2 * Hx1 - Ey1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_Z.row(i * n_comps + 3) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3); 
      }
      else if(data_type == "Full_Impedance_plus_Tipper"){
        // ! calculate L_1_ZpT and L_2_ZpT
        // interpolation operators for impedance component
        // gxx_1
        temp1 = -Hy2 * x_e + Ex2 * y_h;
        temp2 =  Ex1 * Hy2 - Ex2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 0) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxx_2
        temp1 = -Ex1 * y_h + Hy1 * x_e;
        temp2 =  Ex1 * Hy2 - Ex2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 0) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gxy_1
        temp1 = -Ex2 * x_h + Hx2 * x_e;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxy_2
        temp1 = -Hx1 * x_e + Ex1 * x_h;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 1) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyx_1
        temp1 = -Hy2 * y_e + Ey2 * y_h;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 2) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyx_2
        temp1 = -Ey1 * y_h + Hy1 * y_e;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 2) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyy_1 
        temp1 = -Ey2 * x_h + Hx2 * y_e;
        temp2 =  Ey2 * Hx1 - Ey1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 3) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyy_2
        temp1 = -Hx1 * y_e + Ey1 * x_h;
        temp2 =  Ey2 * Hx1 - Ey1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 3) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3); 

        // interpolation operators for tipper component
        // gx_1 for Tx1
        temp1 = Hy2 * z_h - Hz2 * y_h;
        temp2 = Hy2 * Hz1  - Hy1 * Hz2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 4) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gx_2 for Tx2
        temp1 = Hz1 * y_h - Hy1 * z_h;
        temp2 = Hy2 * Hz1  - Hy1 * Hz2; 
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 4) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gy_1 for Ty1
        temp1 = Hz2 * x_h - Hx2 * z_h;
        temp2 = Hx1 * Hz2  - Hx2 * Hz1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_ZpT.row(i * n_comps + 5) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gy_2 for Ty2
        temp1 = Hx1 * z_h - Hz1 * x_h;
        temp2 = Hx1 * Hz2  - Hx2 * Hz1; 
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_ZpT.row(i * n_comps + 5) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
      }
      else if(data_type == "Off_Diagonal_Impedance"){
        // ! calculate L_1_diagonal_Z and L_2_diagonal_Z
        // gxy_1
        temp1 = -Ex2 * x_h + Hx2 * x_e;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_Z.row(i * n_comps + 0) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxy_2
        temp1 = -Hx1 * x_e + Ex1 * x_h;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_Z.row(i * n_comps + 0) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyx_1
        temp1 = -Hy2 * y_e + Ey2 * y_h;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_Z.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyx_2
        temp1 = -Ey1 * y_h + Hy1 * y_e;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_Z.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
      }
      else { // for data type of "Off_Diagonal_Impedance_plus_Tipper"
        // ! calculate L_1_diagonal_ZpT and L_2_diagonal_ZpT
        // interpolation operators for off-diagonal impedance component
        // gxy_1
        temp1 = -Ex2 * x_h + Hx2 * x_e;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_ZpT.row(i * n_comps + 0) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gxy_2
        temp1 = -Hx1 * x_e + Ex1 * x_h;
        temp2 =  Ex2 * Hx1 - Ex1 * Hx2;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_ZpT.row(i * n_comps + 0) = 
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // gyx_1
        temp1 = -Hy2 * y_e + Ey2 * y_h;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_ZpT.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);
        // gyx_2
        temp1 = -Ey1 * y_h + Hy1 * y_e;
        temp2 =  Ey1 * Hy2 - Ey2 * Hy1;
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_ZpT.row(i * n_comps + 1) =
        -1.0*(inv_deth * temp1 + inv_deth * inv_deth * temp2 * temp3);

        // interpolation operators for tipper component
        // gx_1 for Tx1
        temp1 = Hy2 * z_h - Hz2 * y_h;
        temp2 = Hy2 * Hz1  - Hy1 * Hz2;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_ZpT.row(i * n_comps + 2) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gx_2 for Tx2
        temp1 = Hz1 * y_h - Hy1 * z_h;
        temp2 = Hy2 * Hz1  - Hy1 * Hz2; 
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_ZpT.row(i * n_comps + 2) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gy_1 for Ty1
        temp1 = Hz2 * x_h - Hx2 * z_h;
        temp2 = Hx1 * Hz2  - Hx2 * Hz1;
        temp3 = -Hx2 * y_h + Hy2 * x_h;
        this->L_1_diagonal_ZpT.row(i * n_comps + 3) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
        // gy_2 for Ty2
        temp1 = Hx1 * z_h - Hz1 * x_h;
        temp2 = Hx1 * Hz2  - Hx2 * Hz1; 
        temp3 = -Hy1 * x_h + Hx1 * y_h;
        this->L_2_diagonal_ZpT.row(i * n_comps + 3) = 
        inv_deth * temp1 -  inv_deth * inv_deth * temp2 * temp3;
      }
    }
  }

  //initialize member variables for producing synthetic data
  this->sites_ = sites;
  this->n_sites_ = n_sites;

  // for impedance data
  if( (data_type == "Full_Impedance") ||
      (data_type == "Full_Impedance_plus_Tipper") ){
        // Z_i,j; j = 0,1,2,3, represents Zxx, Zxy, Zyx, Zyy, respectively
        this->Z_.resize(n_sites, 4);
        for(unsigned int i = 0; i < n_sites; i++){
          this->Z_(i,0) = Z[i][0]; // Zxx
          this->Z_(i,1) = Z[i][2]; // Zxy
          this->Z_(i,2) = Z[i][1]; // Zyx
          this->Z_(i,3) = Z[i][3]; // Zyy
        }
      }
  else{
    // diagonal_Z_i,j; j = 0, 1, respresents Zxy, Zyx, respectively
    this->diagonal_Z_.resize(n_sites, 2);
    for(unsigned int i = 0; i < n_sites; i++){
      this->diagonal_Z_(i,0) = Z[i][2]; // Zxy
      this->diagonal_Z_(i,1) = Z[i][1]; // Zyx
    }
  }

  // for tipper data
  if( (data_type == "Full_Impedance_plus_Tipper") ||
      (data_type == "Off_Diagonal_Impedance_plus_Tipper") ){
    // T_i,j; i = 0, ..., (n_sites-1); j = 0, 1, represents Tx, Ty
    this->T_.resize(n_sites, 2);
    for(unsigned int i = 0; i < n_sites; i++){
      this->T_(i, 0) = A[i]; // Tx
      this->T_(i, 1) = B[i]; // Ty
    }
  }

/*
  // !!!claculate apparent resistivities and phases
  //std::vector<std::vector<Real> >     App(n_sites), Phase(n_sites);
  App[i].resize(4); Phase[i].resize(4);
  Real mu_v=mu[i]*EM::MU0;
  // calculate \rho_{uv} and \phi_{uv}, u,v = x,y
  for(int j=0; j<4; j++) {
    App[i][j]= (1./(omega*mu_v))*(std::abs(Z[i][j]))*(std::abs(Z[i][j]));
    Phase[i][j]= 
	std::abs(std::atan(std::imag(Z[i][j])/std::real(Z[i][j]))*180.0/EM::PI);
    }
  // !!!output for vertification
  // print frequency info and apparent resistivitis \rho_a^{yx} on screen
  for(int n=0; n<n_sites; n++) {
    if(n == 0)
      std::cout << "Frequency = " << this->f_ << " Hz:\n";
    //std::cout << std::setprecision(20) << App[n][1] << std::endl;
    // rho_yx, checked by comparing our results with those from Ren's paper
    std::cout << App[n][1] << std::endl;
  }
*/

  return;
}

void PostProcess::Return_Sites(std::vector< Node >& _sites)
{
  _sites  =  this->sites_;
  
  return;
}

void PostProcess::Return_Impedances(EM::MatrixXC& _Z)
{
  _Z  =  this->Z_;

  return;
}

void PostProcess::Return_Impedances_plus_Tippers(EM::MatrixXC& _ZpT)
{
  _ZpT.resize(this->n_sites_, 6);
  for(unsigned int i = 0; i < this->n_sites_; i++)
    for(unsigned int j = 0; j < 6; j++){
      if(j < 4)
        _ZpT(i, j) = this->Z_(i, j);
      else
        _ZpT(i, j) = this->T_(i, j-4);
    }

  return;
}

void PostProcess::Return_Diagonal_Impedances(EM::MatrixXC& _diagonal_Z)
{
  _diagonal_Z = this->diagonal_Z_;

  return;
}

void PostProcess::Return_Diagonal_Impedances_plus_Tippers(EM::MatrixXC& _diagonal_ZpT)
{
  _diagonal_ZpT.resize(this->n_sites_, 4);
  for(unsigned int i = 0; i < this->n_sites_; i++)
    for(unsigned int j = 0; j < 4; j++){
      if(j < 2)
        _diagonal_ZpT(i, j) = this->diagonal_Z_(i, j);
      else
        _diagonal_ZpT(i, j) = this->T_(i, j-2);
    }

  return;
}

// only used for producing synthetic data
void PostProcess::Return_Impedances_Tippers(EM::MatrixXC& _Z, EM::MatrixXC& _T)
{
  _Z = this->Z_;
  _T = this->T_;

  return;
}

void PostProcess::get_n_sites(unsigned int& n_sites)
{
  n_sites = this->n_sites_;

  return;
}
