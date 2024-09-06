/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



// Gauss_Legendre quadrature for unit tet [0,1]x[0,1]x[0,1]

#ifndef _GaussTet_H
#define _GaussTet_H

// C++ includes
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

class GaussTet
{
 public:
  GaussTet (const unsigned int p);
  ~GaussTet() {}  
  // Compute points and weights.
  void init (const unsigned int p);
  // Print
  friend std::ostream& operator <<(std::ostream& os, const GaussTet& GaussTet);

  // Get the number of quadrature points.
  unsigned int n_points() { return _points.size();}
  // Get a const reference to quadrature points.
  const std::vector<Point>& get_points()  {  return _points;}
  // Get a const reference to quadrature weights.
  const std::vector<Real>&  get_weights() {  return _weights;} 

 protected:
  // The quadrature points.
  std::vector<Point>           _points;
  // The quadrature weights.
  std::vector<Real>            _weights;  
};


//--------------------------------------------------------------------
inline
GaussTet::GaussTet(const unsigned int p)
{
  // init.
  this->init(p);


}


inline  
void GaussTet::init(const unsigned int p)
{
  switch (p) {
    case 1:
      {
	// 1-order
	_points.resize(1);
	_weights.resize(1);	
	_points[0](0) = .25;
	_points[0](1) = .25;
	_points[0](2) = .25;	
	_weights[0]   = 1.0;
        return;
      }
    case 2:
      {
	// 2-order
	_points.resize(4);
	_weights.resize(4);		    
	const Real a = .585410196624969;
	const Real b = .138196601125011;
	    
	_points[0](0) = a;
	_points[0](1) = b;
	_points[0](2) = b;   
	_points[1](0) = b;
	_points[1](1) = a;
	_points[1](2) = b;	    
  	_points[2](0) = b;
	_points[2](1) = b;
	_points[2](2) = a;	    
	_points[3](0) = b;
	_points[3](1) = b;
	_points[3](2) = b;
		    
	_weights[0] = 0.25;
	_weights[1] = _weights[0];
	_weights[2] = _weights[0];
	_weights[3] = _weights[0];

	return;
      }  
    case 3:
      {
       // 3-order
        _points.resize(5);
	_weights.resize(5);		    
		    
	_points[0](0) = .25;
	_points[0](1) = .25;
	_points[0](2) = .25;	  
	_points[1](0) = .5;
	_points[1](1) = .16666666666666666666666666666666666666666667;
	_points[1](2) = .16666666666666666666666666666666666666666667;		    
	_points[2](0) = .16666666666666666666666666666666666666666667;
	_points[2](1) = .5;
	_points[2](2) = .16666666666666666666666666666666666666666667;		    
	_points[3](0) = .16666666666666666666666666666666666666666667;
	_points[3](1) = .16666666666666666666666666666666666666666667;
	_points[3](2) = .5;		    
	_points[4](0) = .16666666666666666666666666666666666666666667;
	_points[4](1) = .16666666666666666666666666666666666666666667;
	_points[4](2) = .16666666666666666666666666666666666666666667;	    
		    
	_weights[0] = -0.8000000;
	_weights[1] = .45;
	_weights[2] = _weights[1];
	_weights[3] = _weights[1];
	_weights[4] = _weights[1];
	    
	return;
    }
    case 4:
    case 5:
      {
	_points.resize(14);
	_weights.resize(14);

	const Real a[3] = {0.31088591926330060980,    // a1 from the paper
		            0.092735250310891226402,   // a2 from the paper
			    0.045503704125649649492};  // a3 from the paper

	      // weights.  a[] and w[] are the only floating-point inputs required
	      // for this rule.
	      const Real w[3] = {0.018781320953002641800,    // w1 from the paper
				 0.012248840519393658257,    // w2 from the paper
				 0.0070910034628469110730};  // w3 from the paper

	      // The first two sets of 4 points are formed in a similar manner
	      for (unsigned int i=0; i<2; ++i)
		{
		  // Where we will insert values into _points and _weights
		  const unsigned int offset=4*i;

		  // Stuff points and weights values into their arrays
		  const Real b = 1. - 3.*a[i];

		  // Here are the permutations.  Order of these is not important,
		  // all have the same weight
		  _points[offset + 0] = Point(a[i], a[i], a[i]);
		  _points[offset + 1] = Point(a[i],    b, a[i]);
		  _points[offset + 2] = Point(   b, a[i], a[i]);
		  _points[offset + 3] = Point(a[i], a[i],    b);
		      			    
		  // These 4 points all have the same weights 
		  for (unsigned int j=0; j<4; ++j)
		    _weights[offset + j] = w[i]*6.0;
		} // end for


	      {
		// The third set contains 6 points and is formed a little differently
		const unsigned int offset = 8;
		const Real b = 0.5*(1. - 2.*a[2]);

		// Here are the permutations.  Order of these is not important,
		// all have the same weight
		_points[offset + 0] = Point(b   ,    b, a[2]);
		_points[offset + 1] = Point(b   , a[2], a[2]);
		_points[offset + 2] = Point(a[2], a[2],    b);
		_points[offset + 3] = Point(a[2],    b, a[2]);
		_points[offset + 4] = Point(   b, a[2],    b);
		_points[offset + 5] = Point(a[2],    b,    b);
		  
		// These 6 points all have the same weights 
		for (unsigned int j=0; j<6; ++j)
		  _weights[offset + j] = w[2]*6.0;
	      }
      return;
      }
    default:
      {
	std::cout << "\n";
	std::cout << "Only 5-order of Tet - Fatal error!\n";
	std::cout << "Illegal RULE = " << p << "\n";
        std::abort();
       }
  }// switch end

  // all weights =1
  double all_weights=0.;
  for(unsigned int i=0; i<_weights.size(); i++) 
    all_weights+= _weights[i];
  assert(std::abs(all_weights-1.0)<1e-10);

  return ;
}



inline  
std::ostream& operator <<(std::ostream& os, const GaussTet& g)
{
  os<<"Gauss rule for [0,1]x[0,1]x[0,1] tet with "
    << g._points.size()<<" points\n";

  assert(g._points.size()==g._weights.size());
  for(unsigned int i=0; i<g._points.size(); i++)
    {
      os<<std::fixed<<std::setprecision(8);
      os<<g._weights[i]<<"\t"
	<<g._points[i](0)<<"\t"
	<<g._points[i](1)<<"\t"
	<<g._points[i](2)<<"\n";
    }
  return os;
}

  
#endif //_GaussTetTRI_H


