/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/
#include "surface.h"

Surface::~Surface()
{
  // Delete the elements.
  for(unsigned int i=0;i<(_T.size());++i)
    if(_T[i]!=NULL) {
      delete _T[i]; _T[i]=NULL;
    }
  this->_T.clear();
}

std::ostream& operator <<(std::ostream& os, Surface& surface)
{
  os<<"\n"<<surface._T.size()<<"\n";
  for(unsigned int i=0; i<surface._T.size(); i++) {
    os<<surface._T[i]->get_id()<<"\t";
    for(unsigned int j=0; j<surface._T[i]->n_nodes(); j++)
      os<<(surface._T[i]->get_node(j)->get_id())<<"\t";
    os<<"\n";
  }
    
  os<<"\n";

  return os;
}

void Surface::build_surface(BCInfo& bc_info, 
                            short int facet_marker, 
                            int exclude_tet_marker) 
{
  // write info to _tri_2_tet and _T data strucutres

  this->_tri_2_tet.clear();
  this->_T.clear();
  assert(facet_marker!=0);
  assert(facet_marker == this->_T_marker || facet_marker == this->_T1_marker);

  // facet_marker must be in bc_info._side_markers
  std::set<short int>::iterator it_set= bc_info._side_markers.find(facet_marker);
  assert(it_set!=bc_info._side_markers.end());

  unsigned int n_triangles =0;

  // All boundary tets are in bc_info._boundary_side_marker
  typedef std::multimap<const Element*,
    std::pair<unsigned short int, short int> >::iterator IT;
  IT it_map= bc_info._boundary_side_marker.begin();
  for(; it_map!=bc_info._boundary_side_marker.end(); it_map++) { // a boundary tet
    // get face's marker
    short int marker= (*it_map).second.second;
    if(facet_marker== marker) { 
      //get mother "tet
      Element& elem= const_cast<Element&>(*(*it_map).first);
      if(elem.get_tet_marker()!=exclude_tet_marker) { // no tets with exclude_tet_marker
          // which face
          unsigned short int face= (*it_map).second.first;
          //build this face "tri"
          Tri* tri= static_cast<Tri*>(elem.build_side(face).release()); 
          tri->set_id(n_triangles);
          n_triangles++;
          this->_T.push_back(tri);
          this->_tri_2_tet[tri]=std::make_pair(&elem, face);
      }
    }
  }

  // call this->build_n();
  this->build_n();


  // ready
  this->_initialized= true;


  return;
}




std::pair<const Element*, 
  unsigned short int>  Surface::get_tet_face(Tri* tri)
{
  assert(tri!=NULL);
  typedef  std::map<Tri*, 
    std::pair<const Element*,unsigned short int> >::iterator IT;
  IT it= this->_tri_2_tet.find(tri);
  assert(it!=this->_tri_2_tet.end());
  
  const Element* tet= (*it).second.first;
  unsigned short int face=(*it).second.second;

  assert(tet!=NULL);
  assert(face<4);
  return std::make_pair(tet, face);

}



void Surface::build_n()
{
  unsigned int n_elems= this->_T.size();
  this->_n.clear();
  this->_n.resize(n_elems);

  for(unsigned int e=0; e<n_elems; e++) 
    {
      // get an element
      Tri* tri= this->get_surface()[e];
      const unsigned int ID= tri->get_id();
      assert(ID!= INVALID_UNIT);
      Point          n= tri->get_normal();
      // n
      std::pair<const Element*, unsigned short int> 
	 tet_face= this->get_tet_face(tri);  // where am I?
      Element* tet= const_cast<Element*>(tet_face.first);
      Point gpoint= tet->get_gpoint();
      Point gpoint_tri= tri->get_gpoint();
      if((gpoint-gpoint_tri)*n<0.) { // n is outgoing of surface
                                     // on \partial\Omega, n is outgoing
                                     // on air-earth-interface, n is from earth to air
        n= n*-1.0;
      }
      this->_n[ID]=n;
   }

  return;

}


