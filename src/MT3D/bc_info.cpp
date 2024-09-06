/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/

#include "bc_info.h"
#include "mesh3d.h"

//---------------------------------------------------------------------------
BCInfo::BCInfo(const Mesh3D& m) :
  _mesh3d (m),
_initialized (false)
{  
  this->init();
}

BCInfo::~BCInfo()
{
  this->clear();
}

void BCInfo::init()
{
  
}

void BCInfo::clear()
{
  // clear marker vectors.
  this->_boundary_node_marker.clear();
  this->_boundary_side_marker.clear();
  this->_side_markers.clear();
}


void BCInfo::add_node(const unsigned int n,
	              const short int marker)
{
  this->add_node (&(this->_mesh3d.get_node(n)), marker);
}

void BCInfo::add_node(const Node* n, 
		      const short int marker)
{
  assert(n!=NULL);
  /*
  // the judgement is not necessary from my opinion (added on 2021/05/06)
  if (marker == static_cast<short int>(INVALID_UNIT))
    // Output the error message when using tetgen 1.5.0, 1.5.1 or 1.6.0
    // as their default marker is 0 and -1 for node, edge and face.
    std::cerr<<"Error, a node with marker: "<< marker
	     <<" ,please check your .POLY file\n";  
  //otherwise.
  */
  this->_boundary_node_marker[n] = marker;
}
  
short int BCInfo::get_node_marker(const Node* node) const
{ 
  assert(node!=NULL);
  // search the node.
  std::map<const Node*, short int>::const_iterator
    n = _boundary_node_marker.find(node);  
  // node not in the data structure
  if (n == _boundary_node_marker.end())
    return static_cast<short int>(INVALID_UNIT);  
  return n->second;
}

void 
BCInfo::build_node_list(std::vector<const Node*>& node_list,
			std::vector<short int>&marker_list) const
{
  // Reserve the size, then use push_back
  node_list.reserve     (_boundary_node_marker.size());
  marker_list.reserve   (_boundary_node_marker.size());
  // loop the map.
  std::map<const Node*, short int>::const_iterator pos;  
  for (pos=_boundary_node_marker.begin(); pos != _boundary_node_marker.end();
       ++pos)
    {
      node_list.push_back (pos->first);
      marker_list.push_back (pos->second);
    }
}

void BCInfo::add_side(const unsigned int elem,
			    const unsigned short int side,
			    const short int marker)
{
  this->add_side ( &(this->_mesh3d.get_elem(elem)), side, marker);
}

void BCInfo::add_side(const Element* elem,
			    const unsigned short int side,
			    const short int marker)
{
  assert (elem != NULL); 
  /*
  // the judgement is not necessary from my opinion (added on 2021/05/06)
  if (marker == static_cast<short int>(INVALID_UNIT))
    // Output the error message when using tetgen 1.5.0, 1.5.1 or 1.6.0
    // as their default marker is 0 and -1 for node, edge and face.
    std::cerr<<"Error, a side with marker: "<<marker
	     <<" ,please check your .POLY file\n";  
  */
  std::pair<unsigned short int, short int> p(side,marker);
  std::pair<const Element*, std::pair<unsigned short int, short int> >
  kv (elem, p);
  //adds ids into sides_elem_id map.
  this->_boundary_side_marker.insert(kv);
  this->_side_markers.insert(marker);
}

short int BCInfo::get_side_marker(const Element* elem,
		                  const unsigned short int side) const
{
  // get_side_marker of tetrahedral elements
  assert (elem != NULL);
  // search the side on the elem.
  std::pair<std::multimap<const Element*,
    std::pair<unsigned short int, short int> >::const_iterator,
    std::multimap<const Element*,
    std::pair<unsigned short int, short int> >::const_iterator > 
    e = _boundary_side_marker.equal_range(elem);  
  // elem not in the data structure
  if (e.first == e.second)
    return static_cast<short int>(INVALID_UNIT);  
  // elem is there, maybe multiple occurrences
  while (e.first != e.second)
    {
      // if this is true we found the requested side
      // of the element and want to return the id
      if (e.first->second.first == side)
	return e.first->second.second;      
      ++e.first;
    }  
  // if we get here,the requested side index is not achieved.
  return static_cast<short int>(INVALID_UNIT);  
}

void BCInfo::build_side_marker_list(std::set<short int>&side_marker_list)
{
  // clear the previous set.
  side_marker_list.clear();
  // simple copy.
  side_marker_list=_side_markers;
}
