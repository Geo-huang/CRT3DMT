/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



// Gauss_Legendre quadrature for unit triangle [0,1]x[0,1]
// Inverse algorithm from 1DX1D Gauss rule 
#ifndef _GaussTRI_H
#define _GaussTRI_H

// C++ includes
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
// Local includes
#include "gauss_12D.h"

class GaussTri
{
 public:
  GaussTri (const unsigned int p, const bool using_2D_rule=false);
  ~GaussTri() {}  
  // Compute points and weights.
  void init_2D (const unsigned int p);
  void init_tri(const unsigned int p);
  // Print
  friend std::ostream& operator <<(std::ostream& os, 
  				     const GaussTri& GaussTri);

  // Get the number of quadrature points.
  unsigned int n_points() { return _points.size();}
  // Get a const reference to quadrature points.
  const std::vector<Point>& get_points()  { return _points;}
  // Get a const reference to quadrature weights.
  const std::vector<Real>&  get_weights() { return _weights;} 

  // These three functions are modifed from the DUNAVANT libarary
  // DUNAVANT is a C++ library which defines the weights and abscisass
  // for a sequence of 20 quadrature rules on a triangle, 
  // which are exact for polynomials up to degree 20. 
  // http://people.sc.fsu.edu/~burkardt/cpp_src/dunavant/dunavant.html
  // Author: John Burkardt with the GNU LGPL license. 
  // Get the suborder vector of gauss rule for p
  void get_suborder(const unsigned int p,
		    std::vector<unsigned int>& suborder);
  void get_suborder_L0L1L2_w(const unsigned int p,
			     std::vector<Real>& w,
			     std::vector< std::vector<Real> >& L0L1L2);
  void expand_sub_to_full_order(const unsigned int p);
  
 protected:
  // The quadrature points.
  std::vector<Point>           _points;
  // The quadrature weights.
  std::vector<Real>            _weights;  
};


//--------------------------------------------------------------------
inline
GaussTri::GaussTri(const unsigned int p,const bool using_2D_rule)
{
  // init.
  if(using_2D_rule==true)
    this->init_2D(p);
  else
    this->init_tri(p);
}


inline  
void GaussTri::init_2D(const unsigned int p)
{
  Gauss1D g1d(p);
  const std::vector<Real>&  q=g1d.get_points();
  const std::vector<Real>&  w=g1d.get_weights();
  const unsigned int n=q.size();
  this->_points.resize(n*n);
  this->_weights.resize(n*n);
  for(unsigned int i=0; i<n; i++)
   for(unsigned int j=0; j<n; j++) {
     _weights[i*n+j]= ((double)1.0+q[i])*w[i]*w[j]*(double)0.125*2;
     _points [i*n+j](0)= ((double)1.0+q[i])*((double)1.0+q[j])*(double)0.25;
     _points [i*n+j](1)= ((double)1.0+q[i])*((double)1.0-q[j])*(double)0.25; 
     _points [i*n+j](2)= double(1.0)-_points[i*n+j](0)-_points[i*n+j](1);
   }
  //checking all weights ==1
  double all_weights=0.;
  for(unsigned int i=0; i<_weights.size(); i++) 
    all_weights+= _weights[i];
  assert(std::abs(all_weights-1.0)<1e-10);
  
  // checking all points inside triangle
  for(unsigned int i=0; i<_points.size(); i++) {
    for(unsigned int j=0; j<3; j++) {
      if(_points[i](j)<0.) {
        std::cout<<"WARMING: points outside triangle at"
                 <<" order "<<p<<"\n";
      }
    }
  }
  return ;
}



inline  
std::ostream& operator <<(std::ostream& os, const GaussTri& GaussTri)
{
  os<<"Gauss rule for [0,1]x[0,1] triangle with "
    << GaussTri._points.size()<<" points\n";

  assert(GaussTri._points.size()==GaussTri._weights.size());
  for(unsigned int i=0; i<GaussTri._points.size(); i++)
    {
      os<<std::scientific<<std::setprecision(12);
      os<<GaussTri._weights[i]<<"\t"
	<<GaussTri._points[i](0)<<"\t"
	<<GaussTri._points[i](1)<<"\t"
	<<GaussTri._points[i](2)<<"\n";
    }
  return os;
}

  
#endif //_GaussTriTRI_H


