/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/


#ifndef  _ELEMENT_H
#define  _ELEMENT_H

// C++ includes
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
// Local includes   
#include "node.h"

class Element
{   
 public:  
  Element(const unsigned int n_nodes, const unsigned int id,
	  const unsigned int n_side);
  virtual ~Element();

  friend std::ostream& operator <<(std::ostream&os, Element& elem);   

  // Return the id of the element.A const version.
  unsigned int get_id () const;
  // Return the id of the element.A non const version.
  unsigned int get_id();
  // Assign the id to the element.
  void set_id (const unsigned int id);
  // Return true, if the element has a valid id, 
  bool valid_id () const;

  // Return the ith node const pointer.
  const Node*  get_node(const unsigned int i) const; 
  // Return the ith node non-const pointer.
  Node* get_node(const unsigned int i); 
  // Assign the node n at the ith local position.
  void  set_node(const unsigned int i,  Node* n);
  // Return the id of ith node.
  unsigned int get_node_id(const unsigned int i) const;

  // Get the attribute as a const version.
  const std::vector<Real>& get_attribute() const;
  // Get the i-th attribute as a non-const version.
  std::vector<Real>& get_attribute();
  // Assign the value of attribute.
  void set_attribute(const std::vector<Real>& value);
  // Get tet marker
  unsigned int get_tet_marker() const
    { return static_cast<unsigned int>(_attribute[0]); }
  unsigned int get_tet_marker() 
    { return static_cast<unsigned int>(_attribute[0]); }
  
  // Get the jth neighbor element const pointer.
  const Element*  get_neighbor(const unsigned int j) const;
  // Get the jth neighbor element non-const pointer.
  Element* get_neighbor(const unsigned int j);
  // Assign the element e at the jth neighbor position.
  void set_neighbor(const unsigned int j, Element* e);
  // Return true,if the element e is an element neighbor.
  bool is_neighbor(const Element* e) const;
  // Returns true, if the element on boundary.
  bool on_boundary()const;
  // Returns its local face on boundary
  void get_boundary_side(std::vector<unsigned short int>& side);

  // Get the numbers of nodes as a const version.
  unsigned int n_nodes() { return _nodes.size();}
  // Get the numbers of nodes as a non-const verstion.
  unsigned int n_nodes() const { return _nodes.size();}
  // Get the numbers of neighbors or sides as a const.
  unsigned int n_sides() { return _neighbors.size();}
  // Get the numbers of neighbors or sides as a non-const.
  unsigned int n_sides() const { return _neighbors.size();}

  //----------------------------------------------------
  // Get the numbers of nodes.
  virtual unsigned int n_vertices()=0;
  // Get the numbers of edges.
  virtual unsigned int n_edges()=0;
  // Get the element volume.
  virtual Real  get_size()=0;
  // Build the side_Th side element on the master element.
  virtual std::auto_ptr<Element> build_side(const unsigned int side); 
  // Build the edge of element
  virtual void get_edge(const unsigned int edge, std::vector<Node*>& nodes) 
  {  
   std::cout << "Calling get_edge function in elem.h file is wrong!\n";
   std::abort(); 
  }
  // Build the edge of element
  virtual Point get_gpoint()
  {  
   std::cout << "Calling get_gpoint function in elem.h file is wrong!\n";
   std::abort(); 
  }
  // Build the normla of element
  virtual const Point get_normal()
  {  
   std::cout << "Calling get_normal function in elem.h file is wrong!\n";
   std::abort(); 
  }
  virtual Real get_diameter() 
  {  
   std::cout << "Calling get_diameter function in elem.h file is wrong!\n";
   std::abort(); 
  }


 protected:
  unsigned int              	     _id;
  std::vector<Node*>	             _nodes;
  std::vector<Element*>		     _neighbors;
  std::vector<Real>                  _attribute;
};


// ------------------------------------------------------------------------
inline
Element::Element(const unsigned int n_nodes, const unsigned int id,
                 const unsigned int n_side):
  _id          (id),
  _nodes       (n_nodes,NULL),
  _neighbors   (n_side,NULL)
{
}   

inline
Element::~Element() 
{
}


inline
unsigned int Element::get_id () const
{
  return this->_id;
}

inline
unsigned int Element::get_id () 
{
  return this->_id;
}

inline
void Element::set_id (const unsigned int id)
{
  this->_id=id;
}


inline
bool Element::valid_id () const
{
  return (this->_id!= INVALID_UNIT);
}

inline
const Node* Element::get_node(const unsigned int i) const
{
  assert (i < n_nodes());
  assert (this->_nodes[i] != NULL);
  return this->_nodes[i];
}

inline
Node* Element::get_node(const unsigned int i)
{
  assert (i < n_nodes());
  assert (this->_nodes[i] != NULL);
  return  this->_nodes[i];
}

inline
void Element::set_node(const unsigned int i, Node* n)
{
  assert (i < n_nodes());;
  assert (n!=NULL);
  this->_nodes[i]=n;  
}

inline
unsigned int Element::get_node_id(const unsigned int i) const
{
  assert (i < n_nodes());
  assert (this->_nodes[i] != NULL);
  assert (this->_nodes[i]->get_id() !=INVALID_UNIT);
  return this->_nodes[i]->get_id();  
}


inline
std::ostream& operator <<(std::ostream&os, Element& elem)
{
  os<<elem._id<<"\t";
  for(unsigned int i=0; i<elem.n_nodes();i++)
    os<<(elem.get_node_id(i))<<"\t";
  os<<"\n";

  return os;
}


inline
const Element* Element::get_neighbor (const unsigned int j) const
{
  assert (j < n_sides());
  assert( j< this->_neighbors.size());
  return  this->_neighbors[j];
}


inline
Element*  Element::get_neighbor (const unsigned int j) 
{
  assert (j < n_sides());
  assert( j< this->_neighbors.size());
  return  this->_neighbors[j];
}

inline
void Element::set_neighbor (const unsigned int j, Element* e)
{
  assert (j < n_sides());
  this->_neighbors[j] = e;
}

inline
bool Element::is_neighbor (const  Element* e) const
{
  for (unsigned int n=0; n<this->n_sides(); n++)
    if (this->get_neighbor(n) == e)
      return true;
  return false;
}




inline
bool Element::on_boundary () const
{
  // The element is on the boundary,if it has a NULL neighbor.
  return this->is_neighbor(NULL);
}



inline
void  Element::get_boundary_side(std::vector<unsigned short int>& side)
{
  assert(on_boundary()==true); // must on boundary
  side.clear();
  for (unsigned int n=0; n<this->n_sides(); n++) {
    if(this->get_neighbor(n) == NULL) {
      side.push_back(n);
    }
  }
  std::sort(side.begin(), side.end());
  assert(side.size()>=1);

}




inline
const std::vector<Real>& Element::get_attribute() const
{
  return this->_attribute;
}


inline
std::vector<Real>& Element::get_attribute()
{
  return this->_attribute;
}

inline
void Element::set_attribute(const std::vector<Real>& value)
{
  this->_attribute=value;
}


inline
std::auto_ptr<Element> Element::build_side (const unsigned int side) 
{
  assert(side<3);
  // Run here,something is wrong.
  std::cout << "Calling build_side function in elem.h file is wrong!" 
            << std::endl;
  std::abort();
  return std::auto_ptr<Element>(NULL); 
}

#endif // _ELEMENT_H
