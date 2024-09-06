/*****************************************************************************/
/*                                                                           */
/*  Copyright 2019                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/

#include "bc.h"


//------------------------Implementation-uninline functions--------------------
void BC::compute_E_H(const Point& p)
{
   // Clear 
   this->E.clear();  this->H.clear();
   // two modes
   if(mode==TE)  TE_oblique_E_H(p);
   else if(mode==TM) TM_oblique_E_H(p);
}

int BC::p_in_which_layer(const Point& p)
{
   const Real z=p(2);
   double h_up[this->N+1], h_down[this->N+1];
   for(int i=0; i<this->N+1; i++) {
      if(i==0) { h_up[i]= -1E50; h_down[i]=this->depth[0]; } 
      else { h_up[i]= this->depth[i-1]; h_down[i]=this->depth[i]; } 
   }
   int located_n = -1; // z must be in [-1E50, 1E50], so we must find its layer number
   for(int i=0; i<this->N+1; i++) {
     if(z>h_up[i]&&z<h_down[i]) {
       located_n=i;
       break;
     }else if(std::abs(z-h_up[i])<1e-10) { // on h_up[i]
       located_n=i;
       break;
     }else if(std::abs(z-h_down[i])<1e-10) {
       located_n=i;
       break;
     }
  }

  assert(located_n!=-1);
  assert(located_n>=0&&located_n<this->N+1);
  return located_n;
}

void BC::TE_oblique_E_H(const Point& p)
{    
    //    x  (Ex)
    //    /
    //   / air  
    //   o----------y         TE-mode
    //   | earth
    //   |
    //   |
    //   z

   // to locate which layer z belongs to
   E.resize(4);    // E=(Ex, 0, 0, 0)
   H.resize(4);    // H=(0, Hy, Hz_-, Hz_+);
   int located_n = p_in_which_layer(p);
   const Real y=p(1); 
   const Real z=p(2);

  if(located_n==N) {

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   EM::Dcomplex t1 = EM::II*ky*y - EM::II*kz[located_n]*z;
   EM::Dcomplex t2 = EM::II*ky*y + EM::II*kz[located_n]*z;
   // re_XX and im_XX represent the real and imaginary part respectively
   double re_t1 = std::real(t1), im_t1 = std::imag(t1);
   double re_t2 = std::real(t2), im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   double exp_re_t1 = std::exp(re_t1), exp_re_t2 = std::exp(re_t2);
   EM::Dcomplex exp_im_t1 = std::exp(EM::II * im_t1), 
                exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(1.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(1.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////

   E[0]= A[located_n]*exp_re_t1*exp_im_t1 + B[located_n]*exp_re_t2*exp_im_t2; //Ex
   E[1]= 0.; // Ey
   E[2]= 0.; // Ez_-
   E[3]= 0.; // Ez_+
   H[0]= 0.; // Hx
   H[1]= (A[located_n]*exp_re_t1*exp_im_t1 - B[located_n]*exp_re_t2*exp_im_t2) 
         * Y[located_n];//Hy
   H[2]= E[0]*(EM::II*ky)/z_hat[located_n]; // Hz_-
   H[3]= H[2]; // Hz_+	

  }else {

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   EM::Dcomplex t1 = EM::II*ky*y - EM::II*kz[located_n]*(z-depth[located_n]);
   EM::Dcomplex t2 = EM::II*ky*y + EM::II*kz[located_n]*(z-depth[located_n]);
   // re_XX and im_XX represent the real and imaginary part respectively
   double re_t1 = std::real(t1), im_t1 = std::imag(t1);
   double re_t2 = std::real(t2), im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   double exp_re_t1 = std::exp(re_t1), exp_re_t2 = std::exp(re_t2);
   EM::Dcomplex exp_im_t1 = std::exp(EM::II * im_t1), 
                exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(2.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(2.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////

   E[0]= A[located_n]*exp_re_t1*exp_im_t1 + B[located_n]*exp_re_t2*exp_im_t2; //Ex
   E[1]= 0.; // Ey
   E[2]= 0.; // Ez_-
   E[3]= 0.; // Ez_+
   H[0]= 0.; // Hx
   H[1]= (A[located_n]*exp_re_t1*exp_im_t1 - B[located_n]*exp_re_t2*exp_im_t2) 
         * Y[located_n]; // Hy          
   H[2]= E[0]*(EM::II*ky)/z_hat[located_n]; // Hz_-
   H[3]= H[2]; // Hz_+	

  }
  // E and H are calculated.    
  return;

}

void BC::TM_oblique_E_H(const Point& p)
{
   //    x  (Hx)
   //	  /
   //	  / air  
   //	  o----------y         TM-mode
   //	  | earth
   //	  |
   //	  |
   //	  z

   // to locate which layer z belongs to
   E.resize(4);    // E=(Ex, 0, 0, 0)
   H.resize(4);    // H=(0, Hy, Hz_-, Hz_+);
   int located_n = p_in_which_layer(p);
   const Real y=p(1); 
   const Real z=p(2);

  if(located_n!=N) {

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   EM::Dcomplex t1 = EM::II*ky*y - EM::II*kz[located_n]*(z-depth[located_n]);
   EM::Dcomplex t2 = EM::II*ky*y + EM::II*kz[located_n]*(z-depth[located_n]);
   // re_XX and im_XX represent the real and imaginary part respectively
   double re_t1 = std::real(t1), im_t1 = std::imag(t1);
   double re_t2 = std::real(t2), im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   double exp_re_t1 = std::exp(re_t1), exp_re_t2 = std::exp(re_t2);
   EM::Dcomplex exp_im_t1 = std::exp(EM::II * im_t1), 
                exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(2.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(2.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////

    H[0]= A[located_n]*exp_re_t1*exp_im_t1 + B[located_n]*exp_re_t2*exp_im_t2;
    H[1]= 0.0; // Hy
    H[2]= 0.0; // Hz_-
    H[3]= 0.0; // Hz_+
    E[0]= 0.0; // Ex
    E[1]= (A[located_n]*exp_re_t1*exp_im_t1 - B[located_n]*exp_re_t2*exp_im_t2) 
           * Z[located_n];
    E[2]= H[0]*(-EM::II*ky)/(y_hat[located_n]); // Ez_-
    E[3]= E[2]; // Ez_+  

  }else {

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   EM::Dcomplex t1 = EM::II*ky*y - EM::II*kz[located_n]*z;
   EM::Dcomplex t2 = EM::II*ky*y + EM::II*kz[located_n]*z;
   // re_XX and im_XX represent the real and imaginary part respectively
   double re_t1 = std::real(t1), im_t1 = std::imag(t1);
   double re_t2 = std::real(t2), im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   double exp_re_t1 = std::exp(re_t1), exp_re_t2 = std::exp(re_t2);
   EM::Dcomplex exp_im_t1 = std::exp(EM::II * im_t1), 
                exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(1.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(1.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////

    H[0]= A[located_n]*exp_re_t1*exp_im_t1 + B[located_n]*exp_re_t2*exp_im_t2;
    H[1]= 0.0; // Hy
    H[2]= 0.0; // Hz_-
    H[3]= 0.0; // Hz_+
    E[0]= 0.0; // Ex
    E[1]= (A[located_n]*exp_re_t1*exp_im_t1 - B[located_n]*exp_re_t2*exp_im_t2) 
           * Z[located_n];
    E[2]= H[0]*(-EM::II*ky)/(y_hat[located_n]); // Ez_-
    E[3]= E[2]; // Ez_+  
  }
  // E and H are calculated.   

  return;
}

std::ostream& operator <<(std::ostream& os, const BC& loe)
{
  // impedance Z
  Dcomplex Z=0.0;
  // app_resistivity,
  Real app_resistivity=0.0;
  // phase,
  Real phase=0.0;
  
  // print info about wavenumber
  //std::cout<<loe.E[0]<<"\t"<<loe.H[1]<<"\n";
  //std::cout<<loe.E[1]<<"\t"<<loe.H[0]<<"\n";

  if(loe.mode==EM::TE)   Z = loe.E[0]/loe.H[1]; //Z=Ex/Hy
  else if(loe.mode==EM::TM)  Z = loe.E[1]/loe.H[0]; //Z=Ey/Hx  
  // Calculate the apparent resistivity of Z 
  app_resistivity = std::abs(Z)*std::abs(Z)/std::real((II*loe.z_hat[1]));
  // Calculate the phase of impandence c
  phase = std::atan(std::imag(Z)/std::real(Z))*180/EM::PI;  

  // print 
  std::cout<<"\nk:\t"<<loe.k[0]<<"\t"<<loe.k[1]
           <<"\nzhat:\t"<<loe.z_hat[0]*-1.0<<"\t"<<loe.z_hat[1]*-1.0
           <<"\nyhat:\t"<<loe.y_hat[0]<<"\t"<<loe.y_hat[1]
           <<"\nr:\t"<<loe.ky
           <<"\nu:\t"<<loe.kz[0]<<"\t"<<loe.kz[1];
  // print wavelength
  std::cout<<"\nwavelength: \n"
           //<<"horerantial:(m) "
           //<<1.0/std::abs(std::max(std::real(loe.ky),std::imag(loe.ky)))
           <<"vertical:(m) " 
           <<1.0/std::abs(std::max(std::real(loe.kz[0]),std::imag(loe.kz[0])))<<"\t"
           <<1.0/std::abs(std::max(std::real(loe.kz[1]),std::imag(loe.kz[1])))<<"\n";

  // print for debug
  if(loe.mode==EM::TE) {
    os<<"\nTE\n"
      <<loe.E[0]<<"\t"<<loe.E[1]<<"\t"<<loe.E[2]<<"\t"<<loe.E[3]<<"\n"
      <<loe.H[0]<<"\t"<<loe.H[1]<<"\t"<<loe.H[2]<<"\t"<<loe.H[3]<<"\n"
      <<Z<<"\t"<<app_resistivity<<"\t"<<phase<<"\n";
   } else if(loe.mode==EM::TM) {
    os<<"\nTM\n"
      <<loe.E[0]<<"\t"<<loe.E[1]<<"\t"<<loe.E[2]<<"\t"<<loe.E[3]<<"\n"
      <<loe.H[0]<<"\t"<<loe.H[1]<<"\t"<<loe.H[2]<<"\t"<<loe.H[3]<<"\n"
      <<Z<<"\t"<<app_resistivity<<"\t"<<phase<<"\n";
  }
   
  return os;
}

void BC::init()
{
  // In this function,we compuate the Coefficients in each layer 
  assert(this->N>=1);
  this->A.resize(this->N+1);
  this->B.resize(this->N+1);
  this->T.resize(this->N+1);
  this->S.setZero();

  // variable definitions for tackling overflow problem
  EM::Dcomplex t1 = 0., t2 = 0.;
  double re_t1 = 0., im_t1 = 0., re_t2 = 0., im_t2 = 0.;
  double exp_re_t1 = 0., exp_re_t2 = 0.;
  EM::Dcomplex exp_im_t1 = 0., exp_im_t2 = 0.;

  // compute Tn, n=1,2,...,N-1
  for(int i=1; i<=N-1; i++) {
    T[i].setZero();

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   t1 =  EM::II*kz[i]*h[i];
   t2 = -EM::II*kz[i]*h[i];
   // re_XX and im_XX represent the real and imaginary part respectively
   re_t1 = std::real(t1); im_t1 = std::imag(t1);
   re_t2 = std::real(t2); im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   exp_re_t1 = std::exp(re_t1); exp_im_t1 = std::exp(EM::II * im_t1);
   exp_re_t2 = std::exp(re_t2); exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(3.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(3.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////

    if(this->mode==EM::TM) {
      T[i](0,0) = exp_re_t1 * exp_im_t1 * ( Z[i]/Z[i-1]+1.0)*0.5;
      T[i](0,1) = exp_re_t2 * exp_im_t2 * (-Z[i]/Z[i-1]+1.0)*0.5;
      T[i](1,0) = exp_re_t1 * exp_im_t1 * (-Z[i]/Z[i-1]+1.0)*0.5;
      T[i](1,1) = exp_re_t2 * exp_im_t2 * ( Z[i]/Z[i-1]+1.0)*0.5;
    }else if (this->mode==EM::TE) {
      T[i](0,0) = exp_re_t1 * exp_im_t1 * ( Y[i]/Y[i-1]+1.0)*0.5;
      T[i](0,1) = exp_re_t2 * exp_im_t2 * (-Y[i]/Y[i-1]+1.0)*0.5;
      T[i](1,0) = exp_re_t1 * exp_im_t1 * (-Y[i]/Y[i-1]+1.0)*0.5;
      T[i](1,1) = exp_re_t2 * exp_im_t2 * ( Y[i]/Y[i-1]+1.0)*0.5;
    }
    //std::cout<<i<<"\n"<<T[i]<<"\n";
    //S*=T[i]; 
  } 
  

  /////////////////////////////////////////////////////////////////////////////////
   // !!! tend to tackle the overflow problem brought by exp function
   t1 = -EM::II*kz[N]*depth[N-1];
   t2 =  EM::II*kz[N]*depth[N-1];
   // re_XX and im_XX represent the real and imaginary part respectively
   re_t1 = std::real(t1); im_t1 = std::imag(t1);
   re_t2 = std::real(t2); im_t2 = std::imag(t2);
   // std::exp(t1,2,3) = exp_re_t1,2,3  * exp_im_t1,2,3
   // exp_re_t1,2,3 = std::exp(re_t1,2,3), exp_im_t1,2,3 = std::exp(imag_t1,2,3*i)
   exp_re_t1 = std::exp(re_t1); exp_im_t1 = std::exp(EM::II * im_t1);
   exp_re_t2 = std::exp(re_t2); exp_im_t2 = std::exp(EM::II * im_t2);
   if( std::isinf(exp_re_t1) ){
     exp_re_t1 = EM::MAX_DB;
     //std::cout << "Note(3.0): +inf was substituded by max value of double-type data!\n";
   }
   if( std::isinf(exp_re_t2) ){
     exp_re_t2 = EM::MAX_DB;
     //std::cout << "Note(3.1): +inf was substituded by max value of double-type data!\n";
   }
  ////////////////////////////////////////////////////////////////////////////////
  // TN
  if(this->mode==EM::TM) {
      T[N](0,0) = exp_re_t1 * exp_im_t1 * ( Z[N]/Z[N-1]+1.0)*0.5;
      T[N](0,1) = exp_re_t2 * exp_im_t2 * (-Z[N]/Z[N-1]+1.0)*0.5;
      T[N](1,0) = exp_re_t1 * exp_im_t1 * (-Z[N]/Z[N-1]+1.0)*0.5;
      T[N](1,1) = exp_re_t2 * exp_im_t2 * ( Z[N]/Z[N-1]+1.0)*0.5;
    }else if (this->mode==EM::TE) {
      T[N](0,0) = exp_re_t1 * exp_im_t1 * ( Y[N]/Y[N-1]+1.0)*0.5;
      T[N](0,1) = exp_re_t2 * exp_im_t2 * (-Y[N]/Y[N-1]+1.0)*0.5;
      T[N](1,0) = exp_re_t1 * exp_im_t1 * (-Y[N]/Y[N-1]+1.0)*0.5;
      T[N](1,1) = exp_re_t2 * exp_im_t2 * ( Y[N]/Y[N-1]+1.0)*0.5;
  }
    //std::cout<<N<<"\n"<<T[N]<<"\n";

  S=T[1]; 
  for(int i=2; i<=N; i++)  S = S*T[i];
  //S*=T[N]; 

    //std::cout<<S<<"\n";

  const Dcomplex S22 = S(1,1);
  const Dcomplex S12 = S(0,1);
  //assert(std::abs(S22)>1E-20); //S22!=0.0
  assert(std::abs(S22)> 0.0);
  const Dcomplex BN =  Dcomplex(1.0,0.0)*this->B0/S22;
  const Dcomplex A0 =  S12*this->B0/S22;

  this->A[0] = A0;
  this->B[0] = this->B0;
  this->A[N] = 0.0;
  this->B[N] = BN;

  for(int j=N-1; j>=(int)1; j--) {
    Vector2C A_B_j_1; 
    A_B_j_1(0)  = this->A[j+1];
    A_B_j_1(1)  = this->B[j+1];
    //std::cout<<T[j+1]<<"\n*********************************"<<A_B_j_1<<"\n";
    Vector2C A_B_j = this->T[j+1]*A_B_j_1; 
    this->A[j]= A_B_j(0);
    this->B[j]= A_B_j(1);
  }

  //for(int i=0; i<=N; i++) {
    //std::cout<<i<<"\t"<<A[i]<<"\t\t"<<B[i]<<"\n";
  //}
  return;
}

void BC::inner_parameters(const std::vector<std::vector<Real> >& ep,
                          const Real frequency)
{
  // exp(-iwt) used
  // zhat=-i*u*omega; yhat=cond-i*omega*epsilon
  // k2= -zhat*yhat
  this->N = ep.size()-1; // ep contains air layer, N does not acount for air layer
  assert(N>=1);      // at least, we are doing half-space model
  k.resize(N+1);
  z_hat.resize(N+1);
  y_hat.resize(N+1);
  depth.resize(N+1); 
  // Compute the wavenumber, impedivity, admittivity 
  for(unsigned int n=0; n<N+1; n++) {
    const Real cond=ep[n][0];
    const Real epsilon=ep[n][1]*EM::EPSILON0;
    const Real mu=ep[n][2]*EM::MU0;
    const Real omega=2.*EM::PI*frequency;
    assert(omega>0.);
    z_hat[n]= -EM::II*omega*mu;
    y_hat[n]= cond-EM::II*omega*epsilon;
    k[n]=std::sqrt(-z_hat[n]*y_hat[n]);
    //std::cout<<k[n]<<"\n";
    assert(std::real(k[n])>-1E-50); //a>0
    assert(std::imag(k[n])>-1E-50); //b>0
    depth[n]=ep[n][3];  // depth to botton of each layer
    // std::cout<<n<<"\t"<<depth[n]<<"\n";
    assert(depth[0]<1e-12);   // h0 is zero at the air-earth interface
  }

  // kz,ky
  kz.resize(k.size());   
  const Dcomplex k_air=k[0];  
  this->ky=k_air*std::sin(theta);   
  for(unsigned int i=0; i<kz.size(); i++) { 
     kz[i]= std::sqrt(k[i]*k[i]-ky*ky);
  }
 
  // Z and Y
  //std::cout<<N+1<<"\n";
  this->Z.clear(); this->Y.clear();
  this->Z.resize(N+1);
  this->Y.resize(N+1);
  for(unsigned int n=0; n<N+1; n++) {
    this->Z[n] = -EM::II*kz[n]/this->y_hat[n]; // Ai**-Bi**
    this->Y[n] =  EM::II*kz[n]/this->z_hat[n]; // Ai**-Bi**
  }
  
  // h
  this->h.resize(N+1);
  h[0] = 1E50;
  h[N] = 1E50;
  for(int i=1; i<=N-1; i++) {
    //std::cout<<i<<"\t"<<depth[i]<<"\t"<<depth[i-1]<<"\n";
    assert(depth[i]>depth[i-1]);
    h[i]= this->depth[i] - this->depth[i-1];
  }

  return ; 
}

