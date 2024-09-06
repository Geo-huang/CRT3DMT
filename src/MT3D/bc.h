/*****************************************************************************/
/*                                                                           */
/*  Copyright 2019                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/

// This class offers electricmagnetic fields (E and H) at a given point P 
// for an n-layers Earth model with air at TE and TM modes
// with an oblique incident plance wave (angle, theta0) at a given frequency (f).
// N is the number of layers below the air space (N>=1) 
// N=1, the half-space model

// usage:  
// exp(-i*w*t)
// +z downward, z=0 is the flat air-earth-interface
// For TE-mode,  E=(Ex,0,0),  H=(0,Hy,Hz);
// For TM-model, E=(0,Ey,Ez), H=(Hx,0,0);
// written by Zhengyong Ren, 10/18/2012

// !!!!!!!!!!!!!
// Please note, we define z_hat = -i*omega*mu, please do not use the sub-routines
// outside this class BC!

#ifndef _BC_H
#define _BC_H


#include <vector>
#include <iostream>
#include <limits>
#include "em.h"
#include "point.h" 

//----------------------------------------------------------------------
class BC
{
 public:
  // Constructor 
  BC(const std::vector<std::vector<Real> >& electric_parameters,
                   const Real frequency,
		   const Real theta_in, // incident angle
		   const Mode working_mode, // TM-mode
		   const Real  B0_EM=1E0); // ampltitute of incident Hx-TM, Ex-TE at z=0
  // Destructor
  ~BC();
  
  // init
  void init();
  // Compuate E and H at a given point.
  void compute_E_H(const Point& p);
  // Compuate the E and H at Point p in TM mode
  void TM_oblique_E_H(const Point& p);
  // Compuate the E and H at Point p in TE mode
  void TE_oblique_E_H(const Point& p);
  // A routine to calculate and print the phase and apparent resistivity on the surface
  friend std::ostream& operator <<(std::ostream& os, const BC& loe);
    // transform electric_parameters to wavenumber, admittivity, impedivity etc.
  void inner_parameters(const std::vector<std::vector<Real> >& ep,
                        const Real frequency);
   // return analytical solutions E and H
   // E=Ex,Ey, Ez_-, Ez_+; H=Hx,Hy,Hz_-,Hz_+; 
   void get_E_H(std::vector<Dcomplex>& e, std::vector<Dcomplex>& h) 
       { e=this->E;  h=this->H; } 
   // get z0 and y0 in air layer
   Dcomplex get_z0() { return z_hat[0]*-1.0; } /*Please note that in other codes, we use*/
                                               /*z_hat =  i*omega*\mu */
                                               /*in this class, we use z_hat = -1*omega*\mu */
                                               /*so, there is a -1.0 in this return */
                                               /*which will be used by other class's subroutines*/
   Dcomplex get_y0() { return y_hat[0]; }
   Dcomplex get_k0() { return k[0]; }
   int p_in_which_layer(const Point& p);  
		   
 public:
  // incident angle 
  const Real                   theta;
  // TE or TM, x is the direction norml to the incident plance (y,z)
  const EM::Mode               mode;
  // ampltitute of H and E at z=0;
  const Real                   B0;  

  unsigned int                 N; /*N-layer earth model*/
  std::vector<Dcomplex>        k; /*k=a+b*i,a>0,b>0, k^2 = -z_hat*y_hat*/
  std::vector<Dcomplex>        kz; /*vertical wavenumber in each layer*/
  Dcomplex                     ky; /*horizontal wavenumber in each layer (which is same in all layers*/
  std::vector<Dcomplex>        z_hat; /*z_hat = -i*omega*\mu*/
  std::vector<Dcomplex>        y_hat; /*y_hat = sigma-i*omega*epsilon*/
  std::vector<double>          h; /*thickness of each layer, h>0*/
  std::vector<double>          depth;/*depth to bottom of each layer*/
  std::vector<Dcomplex>        Z; /*Z_{n} = -i*kz_{n}/y_hat_{n}*/ 
  std::vector<Dcomplex>        Y; /*Y_{n} = i*kz_{n}/z_{hat}_{n}*/
  std::vector<Matrix2C>        T; /*Transfer matrix of each layer*/
  Matrix2C                     S; /*T1*T2*...*TN-1*/
  // Coefficients in each layers
  std::vector<Dcomplex>        A, B; /*a pair of coefficients in each layer*/

  //---------speical for a given p----------
  // electric fields E, E=(Ex, Ey, Ez_-, Ez_+)
  std::vector<Dcomplex>        E;
  // magentic fields H, H=(Hx, Hy, Hz_-, Hz_+)
  std::vector<Dcomplex>        H;

};


//------------------------Implementation-------------------------------
inline
BC:: 
BC(const std::vector<std::vector<Real> >& electric_parameters,
   const Real frequency, const Real theta_in, 
	 const Mode working_mode, const Real  C_EM):
	theta                (theta_in*EM::PI/180.0),
	mode                 (working_mode),
	B0                   (C_EM)
{	
  // Prepare k, z_hat, y_hat
  this->inner_parameters(electric_parameters,frequency);
  // to compuate u,r,D,Q
  this->init();  
}


inline
BC::~BC()
{
  E.clear();
  H.clear();
}


#endif //_BC_H



