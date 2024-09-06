/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/


#ifndef _NODE_H 
#define _NODE_H

// C++   includes
#include <iomanip>
#include <string>
// local includes
#include "point.h"

class Node : public Point
{  
 public:
  Node (const Real x=.0,
	const Real y=.0,
	const Real z=.0,
	const unsigned int id=INVALID_UNIT,
	const int marker=0);
  Node (const Node& n);
  Node (const Point p,
	const unsigned int id = INVALID_UNIT);
  ~Node ();

  friend std::ostream& operator <<(std::ostream& os, const Node& node);

  // Overload operators.
  Node& operator=  (const Node& node);
  bool operator == (const Node& node);
  bool operator < (const Node& v) const; 
 
  void set_id(const unsigned int id) { _id=id; }
  unsigned int get_id() { return _id; }
  void set_marker(const int marker) { _marker=marker; }
  int get_marker() { return _marker; }  

 public:
  unsigned int               _id;
  int                        _marker;
};


//-----------------------------------------------------------------------
inline
Node::Node (const Real x, const Real y,const Real z, 
	    const unsigned int id,
	    const int marker) :
  Point    (x,y,z),
     _id    (id),
     _marker(marker)
{ 
}

inline
Node::Node (const Node& n):
  Point (n),
     _id   (n._id),
     _marker(n._marker)
{
}

inline
Node::Node (const Point p, 
	    const unsigned int id ):
  Point(p),
  _id  (id)
{  
}

inline
Node::~Node ()
{
}



inline
Node & Node::operator= (const Node& node)
{
  (*this)(0) = node(0);
  (*this)(1) = node(1);
  (*this)(2) = node(2); 
  _id        = node._id;

  return *this;
}

inline
bool Node::operator== (const Node& node) 
{
  if(_id!=node._id)
    return false;
  else 
    {
      //compare the coordinates.
      return Point(*this)==Point(node);
    }
}


inline
bool Node::operator < (const Node& rhs) const
{
  // Only compare xyz location.
  return (Point)(*this)<(Point)(rhs);
}



inline
std::ostream& operator <<(std::ostream& os, const Node& node)
{
  os<<node._id
    <<"\t"<<node._coords[0]<<"\t"<<node._coords[1]<<"\t"<<node._coords[2]
    <<"\t"<<node._marker;
  return os;
}

#endif // NODE_H
