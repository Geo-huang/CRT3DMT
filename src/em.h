/******************************************************************************
 *    Function       : definition of EM Namespace used by CRT3DMT             *
 *    Author         : Zhengyong Ren, Huang Chen                              *
 *    Copyright      : (C) Zhengyong Ren, 2017; Huang Chen, 2019              *
 *    Email          : renzhengyong@csu.edu.cn; chenhuang@cqu.edu.cn          *
 *    Created        : in 2017 by Zhengyong Ren                               *
 *    Revised        : on 2019.07.10 by Huang Chen                            *
 ******************************************************************************/

/** compilation under condition **/
#ifndef NAMESPACE_CRT3DMT_H
#define NAMESPACE_CRT3DMT_H
// remove issues with STL containers
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_USE_MKL_ALL

/** files including **/
#include <complex>
// !!!Eigen class template including
// Linear Algebra Lib.
#include <Eigen/Dense>  
// Sparse Lib
#include <Eigen/Sparse> 
#include <Eigen/StdVector>
class ComplexPoint;

//----------------------------------------------------------------------------
namespace EM
{  
  // !!!aternative name of some data-types 
  // here, Real can be changable from double to float
  typedef double Real; 
  typedef std::complex<Real> Dcomplex;  
  // dynamic matrix by using Eigen::Matrix
  typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>        MatrixXI;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>       MatrixXD;
  typedef Eigen::Matrix<Dcomplex, Eigen::Dynamic, Eigen::Dynamic>   MatrixXC;  
  // dynamic column-vector by using Eigen:Matrix  
  typedef Eigen::Matrix<int, Eigen::Dynamic, 1>                     VectorXI;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>                    VectorXD;
  typedef Eigen::Matrix<Dcomplex, Eigen::Dynamic, 1>                VectorXC;  
  typedef Eigen::Matrix<Dcomplex, 3, 1>                             Grad;
  typedef Eigen::Matrix<Real,     6, 6>                             Matrix6D;
  typedef Eigen::Matrix<Dcomplex, 6, 6>                             Matrix6C;
  typedef Eigen::Matrix<Dcomplex, 2, 2>                             Matrix2C;
  typedef Eigen::Matrix<Dcomplex, 2, 1>                             Vector2C; 
  typedef Eigen::Matrix<Real,     4, 1>                             Vector4D;  
  // used for initializing sparse matrix, which stores three 
  // elements (row number,column number, value(type is variable))
  typedef Eigen::Triplet<double>                                    T;
  typedef Eigen::SparseMatrix<Real, Eigen::RowMajor>                EigenSparseMatrixXD;
  typedef Eigen::SparseMatrix<Dcomplex, Eigen::RowMajor>            EigenSparseMatrixXC; 
  typedef Eigen::SparseVector<Real, Eigen::RowMajor>                EigenSparseVectorXD;
  typedef Eigen::SparseVector<Dcomplex, Eigen::RowMajor>            EigenSparseVectorXC;
  typedef Eigen::VectorXcd                            	      EigenDenseVector;  
  // A column vector expression mapping an existing expression
  // tested!
  typedef Eigen::Ref<const VectorXD>                                RefConstVec;   

  // !!!declarations and definitions of the
  // global variables 
  static const Real           PI           = 3.1415926535897932384626433832795;
  static const Dcomplex       ZERO         = Dcomplex(0., 0.);
  static const Dcomplex       II           = Dcomplex(0., 1.);
  // max value of unsigned-int-type data
  static const unsigned int   INVALID_UNIT = static_cast <unsigned int> (-1);
  // max value of double-type data
  static const unsigned int   MAX_DB       = std::numeric_limits<double>::max();
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
#endif // NAMESPACE_CRT3DMT_H
