/******************************************************************************
 *    Function       : Namespace definition for adinvem software              *
 *    Author         : Zhengyong Ren, Huang Chen                              *
      Copyright      : (C) Zhengyong Ren, 2017; Huang Chen, 2019              *
 *    Email          : renzhengyong@csu.edu.cn; csuchenhuang@csu.edu.cn       *
 *    Created        : in 2017 by Zhengyong Ren                               *
 *    Revised        : on 2019.07.10 by Huang Chen                            *
 ******************************************************************************/

/** compilation under condition **/
#ifndef NAMESPACE_ADINVEM_H
#define NAMESPACE_ADINVEM_H
// remove issues with STL containers
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_USE_MKL_ALL

/** files including **/
#include <complex>
class ComplexPoint;

//----------------------------------------------------------------------------
namespace EM
{  
  // !!!aternative name of some data-types 
  // here, Real can be changable from double to float
  typedef double Real; 
  typedef std::complex<Real> Dcomplex;  

  // !!!declarations and definitions of the
  // global variables 
  static const Real           PI           = 3.1415926535897932384626433832795;
  static const Dcomplex       ZERO         = Dcomplex(0., 0.);
  static const Dcomplex       II           = Dcomplex(0., 1.);
  static const unsigned int   INVALID_UNIT = static_cast <unsigned int> (-1);
  static const Real           TOL          = 1e-14;  
  static const Real           MU0          = 4.0 * PI * 1e-7;
  static const Real           EPSILON0     = 8.854187817 * 1e-12;
  // used in MT3D forward modelling by Zhengyong, hence added here 
  enum Mode {TE, TM, INVALID_MODE};
  // + or - signs 
  double sgn(const double m);
}

inline 
double EM::sgn(const double m) 
{
  if (std::abs(m)<EM::TOL)
    return 0.;
  else 
    return m/std::abs(m);
}

using namespace EM;
#endif // NAMESPACE_ADINVEM_H
