/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/


#ifndef _GUASS_12D_H
#define _GUASS_12D_H


// C++ includes
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "em.h"
#include "point.h"
// -----------------------------------------------------------------------
// 1D Gauss class definition [-1,+1]
class Gauss1D
{
 public:
  // Constructor. 
  Gauss1D (const unsigned int p); // p-point rule for 2p-1 order.
  ~Gauss1D() {}  
  // Initialize points and weights.
  void init (const unsigned int p);  

  // Get the numbers of quadrature points.
  unsigned int n_points() { return _points.size();}
  // Get the const reference to the quadrature points.
  const std::vector<Real>&  get_points()  { return _points;}
  // Get the const reference to the quadrature weights.
  const std::vector<Real>&  get_weights() { return _weights;} 

  // print for debugging
  friend std::ostream& operator <<(std::ostream& os, const Gauss1D& Gauss1D);
  
 protected:
  
  // The quadrature points.
  std::vector<Real>            _points;
  // The quadrature weights. 
  std::vector<Real>            _weights;  
};

// -----------------------------------------------------------------------
// 2D Gauss class definition [-1,1]x[-1,1]
class Gauss2D
{
 public:
  // Constructor. 
  Gauss2D (const unsigned int p); // p-point rule for 2p-1 order.
  Gauss2D () {}
  // Destructor.
  ~Gauss2D() {}  
  // Initialize points and weights.
  void init (const unsigned int p);  

  // Get the numbers of quadrature points.
  unsigned int n_points() { return _points.size();}
  // Get the const reference to the quadrature points.
  const std::vector<Point>&  get_points()  { return _points;}
  // Get the const reference to the quadrature weights.
  const std::vector<Real>&   get_weights() { return _weights;} 

  // print for debugging
  friend std::ostream& operator <<(std::ostream& os, const Gauss2D& Gauss2D);
  
 protected:  
  // The quadrature points.
  std::vector<Point>           _points;
  // The quadrature weights. 
  std::vector<Real>            _weights;  
};


#endif // _GUASS_12D_H
