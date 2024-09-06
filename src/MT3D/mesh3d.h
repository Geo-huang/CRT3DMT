/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



/*
  Read Tetgen files-------->(Mesh3D BCInfo)
*/


#ifndef  _MESH3D_H
#define  _MESH3D_H

#include <vector>
#include <string>
#include <cassert>
#include "node.h"
#include "elem.h"
#include "bc_info.h"
#include "tri.h"
#include "tet.h"
#include "surface.h"


//---------------------------------------------------------------------------
class Mesh3D
{  
 public:
  Mesh3D(const std::string name);
  ~Mesh3D();
  void init();
  void clear();
 

  typedef std::vector<Node*>::iterator           node_iterator;
  typedef std::vector<Node*>::const_iterator     const_node_iterator; 
  typedef std::vector<Element*>::iterator        elem_iterator;
  typedef std::vector<Element*>::const_iterator  const_elem_iterator;

  friend              class Element;
  friend              class BCInfo;

  unsigned int mesh_dimension()const { return _dim; }

  unsigned int n_nodes() const  { return _nodes.size(); }
  const Node&  get_node (const unsigned int i) const; 
  Node& get_node(const unsigned int i);
  void add_node (Node* other_node);

  unsigned int n_elems()  const { return _elements.size(); }
  const Element&  get_elem (const unsigned int i) const;
  Element&  get_elem (const unsigned int i);
  void add_elem (Element*  other_elem);

  const std::vector<Element*>& get_elem() const { return _elements;}
  std::vector<Element*>& get_elem() { return _elements; }

  const BCInfo&  get_boundary_info()const;
  BCInfo&  get_boundary_info();

  bool initialized() const { return _initialized; }

  friend std::ostream& operator <<(std::ostream& os, const Mesh3D& mesh);
  void write_out_vtk(const std::string name, 
                     std::vector<Real> elem_data= std::vector<Real>(0));
  void write_out_vtk(const std::string name,  std::vector<unsigned int> elem);
  void write_out_T_vtk(const std::string name);
  void write_out_T1_vtk(const std::string name);


  void prepare_mesh_by_reading_files();  
  void read_tetgen_files();
  void skip_comment_lines(std::istream &in, const char comment_start);
  void node_in (std::istream& node_stream);
  void element_in (std::istream& ele_stream);
  void neighbor_face_in(std::istream& nei_stream,std::istream& face_stream);

  void debugging();  

 public:
  std::string                          _name_input_file;
  const unsigned int	               _dim;  // _dim==3
  std::vector<Node*>                   _nodes;
  std::vector<Element*>	               _elements;
  bool                                 _initialized;
  BCInfo			                         _bc_info;
};    

#endif // _MESH3D_H


