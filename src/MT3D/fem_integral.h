/*****************************************************************************/
/*                                                                           */
/*  Copyright 2020                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/
#ifndef  _FEM_Integral_H
#define  _FEM_Integral_H

// C++ includes
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <string>
#include <utility>
#include "em.h"
#include "bc.h"
#include "tet.h"
#include "tri.h"
#include "gauss_tri.h"
#include "gauss_tet.h"

class FEM_Integral
{
 public: 
  FEM_Integral(Tet& tet);
  FEM_Integral();
  ~FEM_Integral();
 void init(); //checked

 public:

  // the value of Vector shape functions (N) at point p
  // address passing as the formal parameter is an array
  void N(Point p, Point n[6]); //checked
  // the value of the curl of vector shape functions (curl_N) at point t;
  // address passing as the formal parameter is an array
  void Curl_N(Point p, Point curl_n[6]); //checked
  // Volume integral of Int_{v}(Curl_Ni_dot_Curl_Nj) and Int_{v}(Ni_dot_Nj)
  void Curl_N_dot_N(Matrix6D& E,  Matrix6D& F);  // checked
  // Surface integral of Int_{Ni_dot_(n_x_A) A is on T1
  void N_dot_n_x_A(unsigned int face, BC& bc, Grad& C); // checked
  // Volume integral of Int_{v}(Ni_dot_Nj),
  // added by Huang Chen for calculating of sensitivity matrix
  void N_dot_N(Matrix6D& F);
   
 private:
  // A determinate of 3x3 small square matrix. //checked
  double det(double a1, double a2, double a3,
	     double b1, double b2, double b3,
	     double c1, double c2, double c3);

  
 private:
  // tet
  Tet*                          tet;
  // a surface integral rule
  GaussTri*                     g_tri;
  // a volume integral rule
  GaussTet*                     g_tet;
  // Volume of tet
  double                        V;
  // four parameters used to defined vector shapes functions on
  // each edge of "tet, see Jin, FEM for EM, p169  
  double                       a[4], b[4], c[4], d[4]; 
  // length of 6 edges
  double                       l[6];

  
  // An array to modify direction of local vector shape functions 
  // (e2n) in a unique way in the xyz global system
  // rule: direction of vector shape function on an edge is
  // from i_1 to i_2 where id(i_1)<id(i_2);
  double                       directions[6];
  //------------------------------------------------------------
  // local indexing rule for our tet
  // see "tet.h" for details
  // a static table connecting edge and its nodes
  int                          e2n[6][2];
  int                          f2n[4][3];
  int                          f2e[4][3];
  //--------------------------------------------------------------
  
};


inline
double FEM_Integral::det(double a1, double a2, double a3,
			 double b1, double b2, double b3,
	                 double c1, double c2, double c3)
{
  //NOTE: checking me
  //http://mathworld.wolfram.com/Determinant.html
  return (a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1);
}

#endif // _FEM_Integral_H


