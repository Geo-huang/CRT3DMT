/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/


#ifndef _DOFS_H
#define _DOFS_H

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "surface.h"
#include "mesh3d.h"
#include "tet.h"


class DOFs 
{
 public:
  DOFs(Mesh3D& mesh_3d);
  ~DOFs() {}
  friend std::ostream& operator <<(std::ostream& os, const DOFs& dofs);

  unsigned int distribute_dofs(Mesh3D& mesh_3d);

  // return N
  unsigned int get_n_dofs()  { return _N; }
  // global dofs(changing edge local number in tet element to the global number  
  // (global dofs) among the edges in the computational domain)
  void global_dofs_E(Tet* elem, 
                     std::vector<unsigned int>& gdofs);
  // global dofs at side  
  void global_dofs_E(Tet* elem, 
                     unsigned int side,
	             std::vector<unsigned int>& gdofs);
  // consistant checking
  void debugging();


 public:
  /*
   FEM case,  dimenion of _global_dofs = M*6
   M is the number of tets in mesh3d
  */
  std::vector< std::vector<unsigned int> > _global_dofs;
  unsigned int                             _N; 
  const Mesh3D*                            _mesh_3d;
};

#endif  // ifndef _DOFS_H
