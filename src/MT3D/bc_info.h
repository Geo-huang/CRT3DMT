/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/

// BCInfo contains information about marked points and marked surfaces of mesh

#ifndef _MESH_2D_H
#define _MESH_2D_H

// C++ includes
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <cassert>
#include <vector>
#include "elem.h"

class Mesh3D;

//----------------------------------------------------------------------------
class BCInfo
{
  public:
  BCInfo (const Mesh3D& m);
  ~BCInfo ();
  void init ();
  void clear ();
  friend class Mesh3D;
  friend class DOFs;
  
  
 public:
  void add_node (const unsigned int n,const short int marker);
  void add_node (const Node* n,       const short int marker);
  short int get_node_marker(const Node* node) const;
  void build_node_list (std::vector<const Node*>& node_list,
			std::vector<short int>  & marker_list) const;  

  void add_side (const unsigned int elem,
		 const unsigned short int side,
		 const short int marker);
  void add_side (const Element* elem,
		 const unsigned short int side,
		 const short int marker);
  short int get_side_marker(const Element* elem,
			    const unsigned short int side) const;
  unsigned int n_side_markers() const { return _side_markers.size();}
  void build_side_marker_list(std::set<short int>& side_marker_list);


  bool initialized()  const { return _initialized; }

 public:
  // Input info
  const Mesh3D&                                _mesh3d;
  std::map<const Node*, short int>             _boundary_node_marker;
  std::multimap<const Element*,
    std::pair<unsigned short int, short int> > _boundary_side_marker;
  std::set<short int>                          _side_markers;
  bool                                         _initialized;
};

#endif // MESH_2D_H


