/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/




/*
	    edge 1
         0---------2
          \       /
           \     /edge 0
    edge 2  \   /
             \ /              
              1

*/

#ifndef  _TRI_H
#define  _TRI_H

// C++ includes
#include <cstdlib>
// Local includes
#include "elem.h"


class Tri: public Element
{
 public:
  Tri (const unsigned int id=INVALID_UNIT);
  ~Tri(){}

  unsigned int n_vertices() { return 3;}
  unsigned int n_edges ()   { return 3;}
  Real get_size();
  // sides of the triangle are edges,
  // each edge is descripted by two nodes
  void get_side(const unsigned int n,std::vector<Node*>& edge); 
  Point get_gpoint();
  const Point get_normal();
  std::auto_ptr<Element> build_side(const unsigned int side);
  Real get_diameter();
};



inline
Tri::Tri(const unsigned int id):
Element(3, id, 3)
{
}



inline
Real Tri::get_diameter() 
{
  std::vector<Real> l(3);
  l[0]= ((*get_node(0))-(*get_node(1))).size();
  l[1]= ((*get_node(1))-(*get_node(2))).size();
  l[2]= ((*get_node(0))-(*get_node(2))).size();

  const Real maxl= *(std::max_element(l.begin(), l.end()));

  return maxl;
}



inline
Real Tri::get_size() 
{
  // The area of a triangle.
  // Algorithm is from, http://mathworld.wolfram.com/Triangle.html
  // edge a,b,C's sizes.
  const Real a=((*get_node(0))+(*get_node(1))*-1.0).size(); 
  const Real b=((*get_node(1))+(*get_node(2))*-1.0).size();
  const Real c=((*get_node(2))+(*get_node(0))*-1.0).size();
  assert(a*b*c!=0);

  const Real p=(a+b+c)*0.5;
  return std::sqrt(p*(p-a)*(p-b)*(p-c));
}




inline
Point Tri::get_gpoint() 
{
  return (*get_node(0)+*get_node(1)+*get_node(2))*(1./3.);
}



inline
void Tri::get_side(const unsigned int n,std::vector<Node*>& edge) 
{
  // clear edge, two nodes of edge of tri3
  edge.clear(); edge.resize(2);
  
  switch (n)
    {
    case 0:
      {
	edge[0]=this->get_node(2);
	edge[1]=this->get_node(1);
	return;
      }
    case 1:
      {
	edge[0]=this->get_node(0);
	edge[1]=this->get_node(2);
	return;
      }  
    case 2:
      {
	edge[0]=this->get_node(1);
	edge[1]=this->get_node(0);
	return;
      }
      
    default:
      {
	std::cout << "\n";
	std::cout << "Only 3 edges of tri3 - Fatal error!\n";
	std::cout << "Illegal RULE = " << n << "\n";
        std::abort();
       }
    }// switch end

  assert(edge.size()==2);
  return;
}



inline
const Point Tri::get_normal() 
{
  //v0-v1-v2
  const Node& v0= *this->get_node(0);
  const Node& v1= *this->get_node(1);
  const Node& v2= *this->get_node(2);
  const Point temp_normal= ((v1-v0).cross(v2-v0)).unit();
  return temp_normal;
}

inline
std::auto_ptr<Element> Tri::build_side (const unsigned int side) 
{
  assert(side<3);
  // Run here,something is wrong.
  std::cout << "Calling build_side function in tri.h file is wrong!" 
            << std::endl;
  std::abort();
  return std::auto_ptr<Element>(NULL); 
}

#endif  // _Tri_h


