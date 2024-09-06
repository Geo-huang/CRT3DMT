/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/

// a set of triangular surface with a same marker
// which is extracted from bc_info class
// tasks:
// 1. build T1 surface which enclose the air and earth domain
// 2. bouid T  surface which is the air-earth interface

#ifndef  _SURFACE
#define  _SURFACE

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "tri.h"
#include "bc_info.h"

class Surface
{
 public:
  Surface() {}
  ~Surface();

  friend class BCInfo;
  friend std::ostream& operator <<(std::ostream& os, Surface& surface);

  void build_surface(BCInfo& bc_info, 
                     short int facet_marker, // 1-----air-earth-interface (T)
                                             // 2-----\partial\Omega (T1)
                     int exclude_tet_marker=-1111); // each facet is shared by 
                                                       // 1 or 2 tets
                                                       // -1111 means no exclusion
                                                       //  9999999 means no tets in air

  void build_n();  // normal vector on surface;

  std::vector<Tri*>& get_surface() { return _T;}
  Point get_n_surface(unsigned int id) { return _n[id]; }
  // Get tri's father "tet" and its local number in tet
  std::pair<const Element*, unsigned short int> get_tet_face(Tri* tri);
  int get_T_marker()  { return _T_marker;  }
  int get_T1_marker() { return _T1_marker; }

 public:     
  std::vector<Tri*>	               _T;           // triangles
  std::vector<Point>               _n;           // normal vector
  std::map<Tri*, 
    std::pair<const Element*,
              unsigned short int> >    _tri_2_tet;   // tri to its tets
  bool                                 _initialized; 
  static const short int               _T_marker =1; // air-earth-interface
  static const short int               _T1_marker=2; // \partial\Omega (outer boundary) 
};    

#endif  // _SURFACE
