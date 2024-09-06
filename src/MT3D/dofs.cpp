/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/


#include "dofs.h"

DOFs::DOFs(Mesh3D& mesh_3d):
 _mesh_3d  (&mesh_3d)
{  
  assert(_mesh_3d->initialized()==true);
  this->distribute_dofs(mesh_3d);
}


void DOFs::global_dofs_E(Tet* elem,   
		         std::vector<unsigned int>& gdofs)
{
  assert(elem!=NULL);
  assert(elem->get_id()<this->_global_dofs.size());  
 
  gdofs.clear();
  gdofs = this->_global_dofs[elem->get_id()];  
  assert(gdofs.size()==6);
  for(unsigned int i=0; i<6; i++) 
    assert(gdofs[i]!=INVALID_UNIT);
  
  return;
}



void DOFs::global_dofs_E(Tet* elem,  unsigned int side,
		         std::vector<unsigned int>& gdofs)
{
/*
 ------------------------
  face (i)  i1   i2   i3
  0         5     3    4
  1         5     2    1
  2         4     0    2
  3         3     1    0
 ------------------------
   (3) face to edges table
 see "tet.h"

*/
  // check 
  assert(elem != NULL);
  gdofs.clear();
  assert(elem->get_id()<this->_global_dofs.size());  
  gdofs.resize(3);
  std::vector<unsigned int>& GDOFS= this->_global_dofs[elem->get_id()];
  // we get it
  switch (side)  {
   case 0:
    {
      gdofs[0]= GDOFS[5];
      gdofs[1]= GDOFS[3];
      gdofs[2]= GDOFS[4];
      return;
    }
   case 1:
    {
      gdofs[0]= GDOFS[5];
      gdofs[1]= GDOFS[2];
      gdofs[2]= GDOFS[1];
      return;
    }
   case 2:
    {
      gdofs[0]= GDOFS[4];
      gdofs[1]= GDOFS[0];
      gdofs[2]= GDOFS[2];
      return;
    }
   case 3:
    {
      gdofs[0]= GDOFS[3];
      gdofs[1]= GDOFS[1];
      gdofs[2]= GDOFS[0];
      return;
    }
   default:
     {
	std::cout << "\n";
	std::cout << "Only 4 sides of Tet - Fatal error!\n";
	std::cout << "Illegal RULE = " << side << "\n";
        std::abort();
     }
  }// switch end

  assert(gdofs.size()==3);
  for(unsigned int i=0; i<3; i++)
    assert(gdofs[i]!=INVALID_UNIT);
  
  return;
}


std::ostream& operator <<(std::ostream& os, const DOFs& dofs)
{
  os<<dofs._N<<"\tDOFs\n";
  // loop and print
  for(unsigned int n=0; n<dofs._global_dofs.size(); n++) 
    {
      os<<n<<"\t:";
      for(unsigned int m=0; m<dofs._global_dofs[n].size(); m++)
	os<<dofs._global_dofs[n][m]<<"\t";
      os<<"\n";
  }
  
  return os;
  
}   


/*
 E in tet for FEM is numbered continously
 global_dofs=[0]....[5] 
                Tet       
*/
unsigned int DOFs::distribute_dofs(Mesh3D& mesh3d)
{
  // CLEAR
  this->_global_dofs.clear();
  // to store edge on where DOFs are numbered already
  std::map<std::vector<unsigned int>, 
           unsigned int> numbered_edges;
  // Reserve size of global_dofs vector
  this->_global_dofs.resize(mesh3d.n_elems());

  // The number of global_dofs
  unsigned int next_dofs=0;
  // First each edge of tet in mesh for FEM part (E)
  for (unsigned int i=0; i<mesh3d.n_elems(); i++) {
    // each tet at least have 6 DOFs, one for one edge
    Element* T= &mesh3d.get_elem(i);
    unsigned int T_id= T->get_id();
    /*
      Init with INVALID_UNIT
    */
    this->_global_dofs[T_id].resize(6);
    for(unsigned int n=0; n<6; n++)
      this->_global_dofs[T_id][n]= INVALID_UNIT;

    std::vector<unsigned int> elem_global_dofs;
    // loop each edge
    for(unsigned int e=0; e<T->n_edges(); e++) {
      std::vector<Node*> edge;
      T->get_edge(e, edge);
      assert(edge.size()==2);
      std::vector<unsigned int> edge_node_id(2);
      edge_node_id[0]= edge[0]->get_id();
      edge_node_id[1]= edge[1]->get_id(); 
      // sort for comparision
      std::sort(edge_node_id.begin(), edge_node_id.end());
      // check whether edge is numbered or not?
      typedef std::map<std::vector<unsigned int>,
			unsigned int >::iterator IT;
      IT it= numbered_edges.find(edge_node_id);
      if(it==numbered_edges.end()) { // not numbered yet
        //using current next_dofs,++
        numbered_edges[edge_node_id] = next_dofs;
        this->_global_dofs[T_id][e] = next_dofs; 
        next_dofs++;
      }else { // already numbered
        this->_global_dofs[T_id][e]= (*it).second;
      }
    }
  }
  this->_N = next_dofs;

  /*
    Checking subroutine
  */
  // debugging();

  // return the numbers of global dofs
  return _N;
}



void DOFs::debugging() 
{
  std::cout<<"debugging of DOFs class began\n";
  // All data of DOFs are prepared
  std::set<std::vector<unsigned int> > edges_tet; // E 
  std::set<unsigned int> id_tet;
  const Mesh3D& mesh3d = *(this->_mesh_3d);

  for (unsigned int i=0; i<mesh3d.n_elems(); i++) {
    Element* T= const_cast<Element*>(&mesh3d.get_elem(i));
    for(unsigned int e=0; e<T->n_edges(); e++) {
      std::vector<Node*> edge;
      T->get_edge(e, edge);
      assert(edge.size()==2);
      std::vector<unsigned int> edge_node_id(2);
      edge_node_id[0]= edge[0]->get_id();
      edge_node_id[1]= edge[1]->get_id(); 
      std::sort(edge_node_id.begin(), edge_node_id.end());
      edges_tet.insert(edge_node_id);
    }
    std::vector<unsigned int> E_dofs;
    global_dofs_E(static_cast<Tet*>(T), E_dofs);
    assert(E_dofs.size()==6);
    for(unsigned int n=0; n<E_dofs.size(); n++)
      id_tet.insert(E_dofs[n]);
  }
  std::cout<<"edges on all tet \t"<<edges_tet.size()<<"\t"<<_N<<"\n";
  assert(edges_tet.size()==this->_N);
  assert(id_tet.size()==this->_N); 
  edges_tet.clear(); id_tet.clear();
    
  std::cout<<"debugging of DOFs class done\n";
  
}



