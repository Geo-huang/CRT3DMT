/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/
/*
	      0
   TET        
             /|\
          0 / | \
           /  |  \2
          /  1|   \
       1  ..4.|....3
          \   |   /
           \  |  /
          3 \ | /5
             \|/
              
              2
 ------------------------
  edge(i)    i1   i2
   0	       0     1
   1         0     2
   2         0     3
   3         1     2
   4         3     1
   5         2     3
 ------------------------
  (1) edge to nodes table

 ------------------------
  face (i)  i1   i2   i3
  0         1     3    2
  1         0     2    3
  2         0     3    1
  3         0     1    2
 ------------------------
  (2) face to nodes table

 ------------------------
  face (i)  i1   i2   i3
  0         5     3    4
  1         5     2    1
  2         4     0    2
  3         3     1    0
 ------------------------
  (3) face to edges table

*/


#ifndef   _TET_H
#define   _TET_H

// C++ includes
// Local includes
#include "tri.h"
#include <algorithm>    // for max_element

//-------------------------------------------------------------------------
// class Tet declares
class Tet: public Element
{
 public:
  // constructor
  Tet (const unsigned int id=INVALID_UNIT);
  // Destructor
  ~Tet(){}

  unsigned int n_vertices() { return 4;}
  unsigned int n_edges ()   { return 6;}
  Real get_size();
  // sides of the tetrahedron are triangles,
  // each triangle is descripted by Class Tri
  void get_side(const unsigned int face, Element& tri); 
  // get one edge of the tetrahedron, which is descripted by two nodes
  void get_edge(const unsigned int edge, std::vector<Node*>& nodes);
  Point get_gpoint();
  std::auto_ptr<Element> build_side(const unsigned int side); 
  Real get_diameter();
  
  
 public:

};


//-----------------------------------------------------------------------
//Edge class implementations.
inline
Tet::Tet(const unsigned int id):
Element(4, id, 4)
{
}


inline Real Tet::get_size() 
{
  // Calculate the volume of a tetrahedron.
  // for the algorithm in details,please refer to 
  // http://mathworld.wolfram.com/Tetrahedron.html
  // edge vectors a ,b c. down casting.
  const Point & a=(*get_node(1))+((*get_node(0))*-1.0);
  const Point & b=(*get_node(1))+((*get_node(2))*-1.0);
  const Point & c=(*get_node(1))+((*get_node(3))*-1.0);
  return std::abs( a*(c.cross(b)) )/6.;
}


inline Real Tet::get_diameter() 
{
  // get the length of longest edge
  std::vector<Real> l(6);
  l[0]= ((*get_node(0))-(*get_node(1))).size();
  l[1]= ((*get_node(0))-(*get_node(2))).size();
  l[2]= ((*get_node(0))-(*get_node(3))).size();
  l[3]= ((*get_node(1))-(*get_node(3))).size();
  l[4]= ((*get_node(1))-(*get_node(2))).size();
  l[5]= ((*get_node(2))-(*get_node(3))).size();

  const Real maxl= *(std::max_element(l.begin(), l.end()));

  return maxl;
}



inline
Point Tet::get_gpoint() 
{
  // the gpoint is called centroid in English
  return (*get_node(0)+*get_node(1)+*get_node(2)+*get_node(3))*(1./4.);
}


inline
void Tet::get_side(const unsigned int face, Element& tri)
{
/*
 ------------------------
  face (i)  i1   i2   i3
  0         1     3    2
  1         0     2    3
  2         0     3    1
  3         0     1    2
 ------------------------
   (2) face to nodes table
*/
  assert(tri.n_nodes()==3);
  // clear edge, 3 nodes of tri of tet
    switch (face)
    {
    case 0:
      {
	tri.set_node(0,this->get_node(1));
	tri.set_node(1,this->get_node(3));
	tri.set_node(2,this->get_node(2));
	return;
      }
    case 1:
      {
	tri.set_node(0,this->get_node(0));
	tri.set_node(1,this->get_node(2));
	tri.set_node(2,this->get_node(3));
	return;
      }  
    case 2:
      {
	tri.set_node(0,this->get_node(0));
	tri.set_node(1,this->get_node(3));
	tri.set_node(2,this->get_node(1));
	return;
      }
    case 3:
      {
	tri.set_node(0,this->get_node(0));
	tri.set_node(1,this->get_node(1));
	tri.set_node(2,this->get_node(2));
	return;
      }

    default:
      {
	std::cout << "\n";
	std::cout << "Only 3 edges of tri3 - Fatal error!\n";
	std::cout << "Illegal RULE = " << face << "\n";
        std::abort();
       }
    }// switch end

  return;

}


inline
std::auto_ptr<Element> Tet::build_side (const unsigned int side) 
{
  // New a Tri3 object for the side element.
  std::auto_ptr<Element> face(new Tri);  
  // Check the side.
  assert(side<this->n_sides());
  this->get_side(side,*face); 

  return face;  

}



inline
void Tet::get_edge(const unsigned int edge_id, std::vector<Node*>& edge)
{
    // clear edge, two nodes of edge of tri3
  edge.clear(); edge.resize(2);
  assert(edge_id<6);
  
  switch (edge_id)
    {
    case 0:
      {
	edge[0]=this->get_node(0);
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
	edge[0]=this->get_node(0);
	edge[1]=this->get_node(3);
	return;
      }
    case 3:
      {
	edge[0]=this->get_node(1);
	edge[1]=this->get_node(2);
	return;
      }  
    case 4:
      {
	edge[0]=this->get_node(3);
	edge[1]=this->get_node(1);
	return;
      }
    case 5:
      {
	edge[0]=this->get_node(2);
	edge[1]=this->get_node(3);
	return;
      }
      
    default:
      {
	std::cout << "\n";
	std::cout << "Only 6 edges of Tet - Fatal error!\n";
	std::cout << "Illegal RULE = " << edge_id << "\n";
        std::abort();
       }
    }// switch end

  assert(edge.size()==2);
  return;

}

#endif     // _Tet_h
