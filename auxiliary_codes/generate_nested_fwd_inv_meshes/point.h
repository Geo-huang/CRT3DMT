/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



#ifndef    _POINT_COMPLEX_POINT_H
#define    _POINT_COMPLEX_POINT_H

#include <cmath>
#include <cassert>
#include "em.h"

class ComplexPoint;

class Point
{
public:  
  Point  (const Real x=0., const Real y=0., const Real z=0.);
  Point (const Point& p); 
  virtual ~Point ();
  
  friend class ComplexPoint;
  friend std::ostream& operator <<(std::ostream& os, const Point& point);

  // Overload operators.
  Point& operator=(const Point& p);
  Real operator () (const unsigned int i) const; 
  Real & operator () (const unsigned int i);
  Point operator + (const Point& v) const;
  ComplexPoint operator + (const ComplexPoint& v) const;
  Point operator - (const Point& v) const;
  Point operator * (const Real a) const;
  ComplexPoint operator * (const Dcomplex a) const;
  Real operator *(const Point& v) const;
  Dcomplex operator *(const ComplexPoint& v) const;
  bool operator == (const Point& v) const;  
  bool operator < (const Point& v)  const; 
  void operator = (const Real a);

  // Math functions
  Point cross(const Point& v) const;
  ComplexPoint cross(const ComplexPoint& v) const;
  Point unit() const;
  Real size() const; 
  void zero();
  Real* get_xyz() { return _coords; }


 protected:
  Real                                _coords[3];
};


class ComplexPoint
{
public:  
  ComplexPoint  (const Dcomplex x=EM::ZERO,
                 const Dcomplex y=EM::ZERO,
                 const Dcomplex z=EM::ZERO);
  ComplexPoint (const ComplexPoint& cp); 
  virtual ~ComplexPoint ();
  
  friend class Point;
  friend std::ostream& operator <<(std::ostream& os, const ComplexPoint& point);

  // Overload operators
  ComplexPoint& operator =(const ComplexPoint& cp);
  Dcomplex operator () (const unsigned int i) const; 
  Dcomplex & operator () (const unsigned int i);
  ComplexPoint operator + (const ComplexPoint& v) const;
  ComplexPoint operator - (const ComplexPoint& v) const;
  ComplexPoint operator * (const Real a) const;
  ComplexPoint operator * (const Dcomplex a) const;
  Dcomplex operator *(const Point& v) const;
  Dcomplex operator *(const ComplexPoint& v) const;
  bool operator == (const ComplexPoint& v) const;  
  void operator = (const Dcomplex a);

  // Math functions
  ComplexPoint cross(const ComplexPoint& v) const;
  ComplexPoint cross(const Point& v) const;
  ComplexPoint unit() const;
  Dcomplex size() const; 
  void zero();  
  Dcomplex* get_xyz() { return _coords; }
 protected:
  Dcomplex                            _coords[3];
};


//-----------------------------------------------------------------------
inline
Point::Point (const Real x,
	      const Real y,
	      const Real z)
{
  _coords[0] = x;  _coords[1] = y;  _coords[2] = z;
}

inline
Point::Point (const Point &p)
{
  for (unsigned int i=0; i<3; i++)
    _coords[i] = p._coords[i];
}

inline
Point& Point::operator=(const Point& p)
{
  _coords[0] = p._coords[0];  
  _coords[1] = p._coords[1];  
  _coords[2] = p._coords[2];

  return (*this);
}


inline
Point::~Point ()
{
  //no space is allocated by the new operator.
  //so,there is no delete [].
}


inline
Real Point::operator ()(const unsigned int i) const
{
  assert (i<3);
  return _coords[i];
}


inline
Real & Point::operator () (const unsigned int i)
{
  assert (i<3);  
  return _coords[i];
}


inline
Point Point::operator + (const Point &p) const
{
  return Point(_coords[0] + p._coords[0],
	       _coords[1] + p._coords[1],
	       _coords[2] + p._coords[2]);	       
}



inline
ComplexPoint Point::operator + (const ComplexPoint& p) const
{
  return ComplexPoint(_coords[0] + p._coords[0],
                      _coords[1] + p._coords[1],
                      _coords[2] + p._coords[2]);
}



inline
Point Point::operator - (const Point& p) const
{
  return Point(_coords[0] - p._coords[0],
	       _coords[1] - p._coords[1],
	       _coords[2] - p._coords[2]);
}


inline
Point Point::operator * (const Real factor) const
{
  return Point(_coords[0]*factor, _coords[1]*factor,_coords[2]*factor);
}



inline
ComplexPoint Point::operator * (const Dcomplex factor) const
{
  return ComplexPoint(_coords[0]*factor, _coords[1]*factor,_coords[2]*factor);
}


inline
Real Point::operator * (const Point&p) const
{
  return (_coords[0]*p(0) + _coords[1]*p(1) + _coords[2]*p(2));
}



inline
Dcomplex Point::operator * (const ComplexPoint&p) const
{
  return (_coords[0]*p(0) + _coords[1]*p(1) + _coords[2]*p(2));
}



inline
bool Point::operator == (const Point& rhs) const
{
  return ((fabs(_coords[0] - rhs._coords[0]) +
	   fabs(_coords[1] - rhs._coords[1]) +
	   fabs(_coords[2] - rhs._coords[2])) < EM::TOL );
}


inline
bool Point::operator < (const Point& rhs) const
{
  //First we assume (this)<rhs true
  if (*this == rhs) return false;  
  if ((*this)(0) < rhs(0)) return true;      //  <
  else if ((*this)(0) > rhs(0)) return false;//  >
  else if (std::abs((*this)(0)-rhs(0))<EM::TOL ){//vx=rhsx
    if ((*this)(1) < rhs(1)) return true;
    else if ((*this)(1) > rhs(1)) return false;
    else if (std::abs((*this)(1)-rhs(1))<EM::TOL ){//vy=rhsy
       if ((*this)(2) < rhs(2)) return true;
       else if ((*this)(2) > rhs(2)) return false;
    }
  }
  return false;
}


inline
void Point::operator = (const Real a)
{
  _coords[0] = a;  _coords[1] = a;  _coords[2] = a;
}


inline
Point Point::cross(const Point& p) const
{
  return Point(_coords[1]*p._coords[2] - _coords[2]*p._coords[1],
	       -_coords[0]*p._coords[2] + _coords[2]*p._coords[0],
	       _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}


inline
ComplexPoint Point::cross(const ComplexPoint& p) const
{
  return ComplexPoint(_coords[1]*p._coords[2] - _coords[2]*p._coords[1],
	             -_coords[0]*p._coords[2] + _coords[2]*p._coords[0],
	              _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}


inline
Point Point::unit() const
{
  const Real length = size();
  return Point(_coords[0]/length,
	       _coords[1]/length, 
	       _coords[2]/length);
}

inline
Real Point::size() const
{
  Real value= std::pow(_coords[0],2)+
    std::pow(_coords[1],2)+
    std::pow(_coords[2],2);

  return std::sqrt(value);    
}


inline
void Point::zero()
{
  _coords[0]=0.;
  _coords[1]=0.;
  _coords[2]=0.;

}

inline
std::ostream& operator <<(std::ostream& os, const Point& point)
{
  os<<point(0)<<"\t"<<point(1)<<"\t"<<point(2);
  return os;
}

//-----------------------------------------------------------------------
inline
ComplexPoint::ComplexPoint (const Dcomplex x,
	                    const Dcomplex y,
	                    const Dcomplex z)
{
  _coords[0] = x;  _coords[1] = y;  _coords[2] = z;
}

inline
ComplexPoint::ComplexPoint (const ComplexPoint &p)
{
  for (unsigned int i=0; i<3; i++)
    _coords[i] = p._coords[i];
}

inline
ComplexPoint& ComplexPoint::operator =(const ComplexPoint& cp)
{
   _coords[0] = cp._coords[0]; 
   _coords[1] = cp._coords[1];  
   _coords[2] = cp._coords[2];
 
   return (*this);
}


inline
ComplexPoint::~ComplexPoint ()
{
  //no space is allocated by the new operator.
  //so,there is no delete [].
}


inline
Dcomplex ComplexPoint::operator ()(const unsigned int i) const
{
  assert (i<3);
  return _coords[i];
}


inline
Dcomplex & ComplexPoint::operator () (const unsigned int i)
{
  assert (i<3);  
  return _coords[i];
}


inline
ComplexPoint ComplexPoint::operator + (const ComplexPoint &p) const
{
  return ComplexPoint(_coords[0] + p._coords[0],
	              _coords[1] + p._coords[1],
	              _coords[2] + p._coords[2]);	       
}



inline
ComplexPoint ComplexPoint::operator - (const ComplexPoint& p) const
{
  return ComplexPoint(_coords[0] - p._coords[0],
	              _coords[1] - p._coords[1],
	              _coords[2] - p._coords[2]);
}


inline
ComplexPoint ComplexPoint::operator * (const Real factor) const
{
  return ComplexPoint(_coords[0]*factor, _coords[1]*factor,_coords[2]*factor);
}



inline
ComplexPoint ComplexPoint::operator * (const Dcomplex factor) const
{
  return ComplexPoint(_coords[0]*factor, _coords[1]*factor,_coords[2]*factor);
}



inline
Dcomplex ComplexPoint::operator * (const Point&p) const
{
  return (_coords[0]*p(0) + _coords[1]*p(1) + _coords[2]*p(2));
}


inline
Dcomplex ComplexPoint::operator * (const ComplexPoint&p) const
{
  return (_coords[0]*p(0) + _coords[1]*p(1) + _coords[2]*p(2));
}


inline
bool ComplexPoint::operator == (const ComplexPoint& rhs) const
{
  return ((std::abs(_coords[0] - rhs._coords[0]) +
	   std::abs(_coords[1] - rhs._coords[1]) +
	   std::abs(_coords[2] - rhs._coords[2])) < EM::TOL );
}


inline
void ComplexPoint::operator = (const Dcomplex a)
{
  _coords[0] = a;  _coords[1] = a;  _coords[2] = a;
}


inline
ComplexPoint ComplexPoint::cross(const ComplexPoint& p) const
{
  return ComplexPoint(_coords[1]*p._coords[2] - _coords[2]*p._coords[1],
	              _coords[2]*p._coords[0] - _coords[0]*p._coords[2],
	              _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}


inline
ComplexPoint ComplexPoint::cross(const Point& p) const
{
  return ComplexPoint(_coords[1]*p._coords[2] - _coords[2]*p._coords[1],
	              _coords[2]*p._coords[0] - _coords[0]*p._coords[2],
	              _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
}


inline
ComplexPoint ComplexPoint::unit() const
{
  const Dcomplex length = size();
  return ComplexPoint(_coords[0]/length,
	              _coords[1]/length, 
	              _coords[2]/length);
}

inline
Dcomplex ComplexPoint::size() const
{
  Dcomplex value= _coords[0]*_coords[0]+
                  _coords[1]*_coords[1]+
                  _coords[2]*_coords[2];

  return std::sqrt(value);    
}


inline
void ComplexPoint::zero()
{
  _coords[0]=EM::ZERO;
  _coords[1]=EM::ZERO;
  _coords[2]=EM::ZERO;

}


inline
std::ostream& operator <<(std::ostream& os, const ComplexPoint& point)
{
  os<<point(0)<<"\t"<<point(1)<<"\t"<<point(2);
  return os;
}


#endif // _POINT_COMPLEX_POINT_H
