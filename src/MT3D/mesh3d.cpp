/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



#include <iostream>
#include <fstream>
#include <map>
#include "mesh3d.h"

// --------------------------------------------------------------------------
Mesh3D::Mesh3D (const std::string  name):
     _dim           (3),
     _nodes         (0),
     _elements      (0),
     _bc_info       (*this),
     _initialized   (false),
     _name_input_file(name)
{
  assert(_dim==3);
  this->init();
}

void Mesh3D::init()
{
  //Nothing to do.
  return;
}

Mesh3D::~Mesh3D()
{
  // Free all allocated space.
  this->clear();     
}


void Mesh3D::clear ()
{
  // Flag we be not initialized.
  this->_initialized=false;
  // Clear Boundary information.
  this->_bc_info.clear();

  // Delete the elements.
  for(unsigned int i=0;i<(_elements.size());++i)
    if(_elements[i]!=NULL) {
      delete _elements[i]; _elements[i]=NULL;
    }
  this->_elements.clear();

  // Delete the nodes.
  for(unsigned int i=0;i<(_nodes.size());++i) {
    if(_nodes[i]!=NULL) {
      delete _nodes[i];    _nodes[i]=NULL;
    }
  }
  this->_nodes.clear();  
}


const Node& Mesh3D::get_node (const unsigned int n) const
{
  assert (n < this->n_nodes());
  assert (this->_nodes[n] != NULL);
  assert (this->_nodes[n]->get_id() != INVALID_UNIT);    
  return (*_nodes[n]);
}


Node& Mesh3D::get_node(const unsigned int n)
{
  assert (n < this->n_nodes());
  assert (this->_nodes[n] != NULL);
  assert (this->_nodes[n]->get_id() != INVALID_UNIT);   
  return (*_nodes[n]);  
}


void  Mesh3D::add_node (Node*  other_node)
{  
  assert(other_node!=NULL);
  // set the id
  other_node->set_id(this->n_nodes());
  // push back.
  this-> _nodes.push_back(other_node); 
}


const Element& Mesh3D::get_elem (const unsigned int i) const
{
  assert (i < this->n_elems());
  assert (_elements[i] != NULL);  
  return *this->_elements[i];
}


Element& Mesh3D::get_elem (const unsigned int i) 
{
  assert (i < this->n_elems());
  assert (_elements[i] != NULL);
  return *this->_elements[i];
}


void Mesh3D::add_elem (Element* other_elem)
{
  assert(other_elem!=NULL);
  // set the id
  other_elem->set_id(this->n_elems());
  // push back.
  this->_elements.push_back(other_elem);
}

const BCInfo& Mesh3D::get_boundary_info()const
{
  return this->_bc_info;
}


BCInfo& Mesh3D::get_boundary_info()
{
  return this->_bc_info;
}



//--------------------------------------------------------------
// Format prints.
std::ostream& operator <<(std::ostream& os, const Mesh3D& mesh)
{
  os << "n_nodes()="    << mesh.n_nodes()
     << " \nn_elems()=" << mesh.n_elems()
     << " \n"; 
  return os;
}

void Mesh3D::debugging()
{
  //std::cout<<"debugging of Mesh3D class began\n";
  //std::cout<<(*this); // No. of nodes and elems
  write_out_vtk("tet_debugging");
  // show the surface mesh of the computational domain
  write_out_T1_vtk("T1_debugging");
  // show the mesh of the air-Earth interface
  write_out_T_vtk("T_debugging");
  //std::cout<<"debugging of Mesh3D class done\n";
  return;
}

void Mesh3D::write_out_vtk(const std::string name,std::vector<Real> elem_data)
{
  // Open the of stream
  std::ofstream vtk_mesh((name+".vtk").c_str()); 
  if(!vtk_mesh.good()) {
    std::cerr<<"Can not open file:\t"<<name+".vtk"
	     <<std::endl;
  } else {
    // Parts 1-2-3, mandatory 
    vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
               // Header info, Really cool data
	    <<"For visulization of ele data\n"
	    <<"ASCII\n"; //ASCII data (not BINARY)
    // Part 4, Geometry/topology, unstructured mesh
    vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n"; // topography and geometry
    // Write out POINTS info (0-->n-1)
    vtk_mesh<<"\nPOINTS\t"<<this->n_nodes()<<"\tdouble\n"; 
    // Loop POINTS to write out coordinates
    for(unsigned int i=0; i<this->n_nodes(); i++) {
      const Node& temp_node= this->get_node(i);
      vtk_mesh<<temp_node(0)<<"\t"  //x-coordinate
	      <<temp_node(1)<<"\t"  //y-coordinate
	      <<temp_node(2)<<"\n"; //z-coordinate
    }
    // Write out CELL info (0-->m-1)
    vtk_mesh<<"\nCELLS\t"<<this->n_elems()<<"\t"<<this->n_elems()*5<<"\n";
    // Loop ELEMS to write out connectivity indices over elem
    for(unsigned int i=0; i<this->n_elems(); i++) {
      const Element& temp_elem= this->get_elem(i);
      vtk_mesh<<(unsigned int)4<<"\t" 
	      <<temp_elem.get_node_id(0)<<"\t" // 0th node id
	      <<temp_elem.get_node_id(1)<<"\t" // 1st node id
	      <<temp_elem.get_node_id(2)<<"\t" // 2nd node id
	      <<temp_elem.get_node_id(3)<<"\n";// 3th node id
    }
    // Write out CELL types (m)
    vtk_mesh<<"\nCELL_TYPES\t"<<this->n_elems()<<"\n";
    for(unsigned int i=0; i<this->n_elems(); i++) {
      vtk_mesh<<(unsigned int)10<<"\n"; // 10--tet
      //For more details please check the figure 2 and 3 in vtk format file
    }

    // Part 5, dataset attributes
    vtk_mesh<<"\nCELL_DATA\t"<<this->n_elems()<<"\n"
            <<"SCALARS tet_marker int 1\n"
	    <<"LOOKUP_TABLE default\n";
    for(unsigned int i=0; i<this->n_elems(); i++) {
       vtk_mesh<<this->get_elem(i).get_tet_marker()<<"\n";
    }
    if(elem_data.size()==this->n_elems()) {
       vtk_mesh<<"\nSCALARS elem_data double 1\n"
               <<"LOOKUP_TABLE default\n";
       for(unsigned int i=0; i<this->n_elems(); i++) 
          vtk_mesh<<elem_data[i]<<"\n";
    }
               
  }// file opened successfully

  vtk_mesh.close();
  // done
}




void Mesh3D::write_out_vtk(const std::string name,
			   std::vector<unsigned int> selected_elem)
{
  // Open the of stream
  std::ofstream vtk_mesh((name+".vtk").c_str()); 
  if(!vtk_mesh.good()) {
    std::cerr<<"Can not open file:\t"<<name+".vtk"
	     <<std::endl;
  } else {
    // Parts 1-2-3, mandatory 
    vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
               // Header info, Really cool data
	    <<"Surface mesh generated for"
	    <<"RMT modeling by Zhengyong Ren,ETHZurich,2009\n"
	    <<"ASCII\n"; //ASCII data (not BINARY)
    // Part 4, Geometry/topology, unstructured mesh
    vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n"; // topography and geometry
    // Write out POINTS info (0-->n-1)
    vtk_mesh<<"\nPOINTS\t"<<this->n_nodes()<<"\tdouble\n"; 
    // Loop POINTS to write out coordinates
    for(unsigned int i=0; i<this->n_nodes(); i++) {
      const Node& temp_node= this->get_node(i);
      vtk_mesh<<temp_node(0)<<"\t"  //x-coordinate
	      <<temp_node(1)<<"\t"  //y-coordinate
	      <<temp_node(2)<<"\n"; //z-coordinate
    }
    // Write out CELL info (0-->m-1)
    vtk_mesh<<"\nCELLS\t"<<selected_elem.size()<<"\t"<<selected_elem.size()*5<<"\n";
    // Loop ELEMS to write out connectivity indices over elem
    for(unsigned int i=0; i<selected_elem.size(); i++) {
      const Element& temp_elem= this->get_elem(selected_elem[i]);
      vtk_mesh<<(unsigned int)4<<"\t" 
	      <<temp_elem.get_node_id(0)<<"\t" // 0th node id
	      <<temp_elem.get_node_id(1)<<"\t" // 1st node id
	      <<temp_elem.get_node_id(2)<<"\t" // 2nd node id
	      <<temp_elem.get_node_id(3)<<"\n";// 3th node id
    }
    // Write out CELL types (m)
    vtk_mesh<<"\nCELL_TYPES\t"<<selected_elem.size()<<"\n";
    for(unsigned int i=0; i<selected_elem.size(); i++) {
      vtk_mesh<<(unsigned int)10<<"\n"; // 10--tet
      //For more details please check the figure 2 and 3 in vtk format file.
    }
    vtk_mesh<<"\n";
               
  }// file opened successfully

  vtk_mesh.close();
  // done
}



void Mesh3D::write_out_T1_vtk(const std::string name)
{
  // used for showing the surface mesh of the computational domain
  // build \partial \Omega
  BCInfo& bc_info = this->get_boundary_info();
  const short int T1_marker= 2; 
  Surface T1;
  T1.build_surface(bc_info, T1_marker);

  // Open the ofstream
  std::ofstream vtk_mesh((name+".vtk").c_str());
  if(!vtk_mesh.good()) {
    std::cerr<<"Can not open file:\t"<<name+".vtk"
	     <<std::endl;
  } else {
    // Write out the fixed header information.
    // Parts 1-2-3, mandatory 
    vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
               // Header info, Really cool data
	    <<"Surface mesh generated for"
	    <<"RMT modeling by Zhengyong Ren,ETHZurich,2009\n"
	    <<"ASCII\n"; //ASCII data (not BINARY)

    // Part 4, Geometry/topology, unstructured mesh
    vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n"; // topography and geometry
    // Write out POINTS info (0-->n-1)
    vtk_mesh<<"\nPOINTS\t"<<this->n_nodes()<<"\tdouble\n"; 
    // Loop POINTS to write out coordinates
    for(unsigned int i=0; i<this->n_nodes(); i++) {
      const Node& temp_node= this->get_node(i);
      vtk_mesh<<temp_node(0)<<"\t"  //x-coordinate
	      <<temp_node(1)<<"\t"  //y-coordinate
	      <<temp_node(2)<<"\n"; //z-coordinate
    }
    // Write out CELL info (0-->m-1)
    vtk_mesh<<"\nCELLS\t"<<T1.get_surface().size()<<"\t"
            <<T1.get_surface().size()*4<<"\n";
    // Loop ELEMS to write out connectivity indices over elem
    for(unsigned int i=0; i<T1.get_surface().size(); i++) {
      const Element& tri= *(T1.get_surface()[i]);
      vtk_mesh<<(unsigned int)3<<"\t" // linear triangle has only 3 vertexes
	      <<tri.get_node_id(0)<<"\t" // 0th node id
	      <<tri.get_node_id(1)<<"\t" // 1st node id
	      <<tri.get_node_id(2)<<"\n";// 2nd node id
    }
    // Write out CELL types (m)
    vtk_mesh<<"\nCELL_TYPES\t"<<T1.get_surface().size()<<"\n";
    for(unsigned int i=0; i<T1.get_surface().size(); i++) {
      vtk_mesh<<(unsigned int)5<<"\n"; // 5--triangle
      //For more details please check the figure 2 and 3 in vtk format file.
    }
  } // if,else

  vtk_mesh.close();
  // done...
}



void Mesh3D::write_out_T_vtk(const std::string name)
{
  // used for showing the mesh of the air-Earth interface
  // build air-earth-interface
  BCInfo& bc_info = this->get_boundary_info();
  const short int T_marker= 1; 
  const int air_tet_marker = 9999999;
  Surface T;
  T.build_surface(bc_info, T_marker, air_tet_marker);

  // SurfMesh* T= (_bc_info._T);
  // assert((T)!=NULL);
  
  // Open the of stream
  std::ofstream vtk_mesh((name+".vtk").c_str());
  if(!vtk_mesh.good()) {
    std::cerr<<"Can not open file:\t"<<name+".vtk"
	     <<std::endl;
  } else {
    // Write out the fixed header information.
    // Parts 1-2-3, mandatory 
    vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
               // Header info, Really cool data
	    <<"Surface mesh generated for"
	    <<"RMT modeling by Zhengyong Ren,ETHZurich,2009\n"
	    <<"ASCII\n"; //ASCII data (not BINARY)

    // Part 4, Geometry/topology, unstructured mesh
    vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n"; // topography and geometry
    // Write out POINTS info (0-->n-1)
    vtk_mesh<<"\nPOINTS\t"<<this->n_nodes()<<"\tdouble\n"; 
    // Loop POINTS to write out coordinates
    for(unsigned int i=0; i<this->n_nodes(); i++) {
      const Node& temp_node= this->get_node(i);
      vtk_mesh<<temp_node(0)<<"\t"  //x-coordinate
	      <<temp_node(1)<<"\t"  //y-coordinate
	      <<temp_node(2)<<"\n"; //z-coordinate
    }
    // Write out CELL info (0-->m-1)
    vtk_mesh<<"\nCELLS\t"<<T.get_surface().size()<<"\t"
            <<T.get_surface().size()*4<<"\n";
    // Loop ELEMS to write out connectivity indices over elem
    for(unsigned int i=0; i<T.get_surface().size(); i++) {
      const Element& tri= *(T.get_surface()[i]);
      vtk_mesh<<(unsigned int)3<<"\t" // linear triangle has only 3 vertexes
	      <<tri.get_node_id(0)<<"\t" // 0th node id
	      <<tri.get_node_id(1)<<"\t" // 1st node id
	      <<tri.get_node_id(2)<<"\n";// 2nd node id
    }
    // Write out CELL types (m)
    vtk_mesh<<"\nCELL_TYPES\t"<<T.get_surface().size()<<"\n";
    for(unsigned int i=0; i<T.get_surface().size(); i++) {
      vtk_mesh<<(unsigned int)5<<"\n"; // 5--triangle
      //For more details please check the figure 2 and 3 in vtk format file.
    }
  } // if,else

  vtk_mesh.close();
  // done...
}




void Mesh3D::skip_comment_lines (std::istream &in, const char comment_start)
{    
  char c, line[256];
  
  while (in.get(c), c==comment_start) 
    in.getline (line, 255);
  
  // put back first character of
  // first non-comment line
  in.putback (c);
}


//----------------------------------------------------------------------
void Mesh3D::node_in (std::istream& node_stream)
{
  BCInfo&          bc_info= this->_bc_info;
  // Check input buffer
  assert (node_stream.good());
  
  unsigned int dimension=0, nAttributes=0, BoundaryMarkers=0;
  unsigned int n_nodes=0;
  node_stream >> n_nodes    // Read the number of nodes from the stream
	      >> dimension        // Read the dimension from the stream
	      >> nAttributes      // Read the number of attributes from stream
	      >> BoundaryMarkers; // Read if or not boundary markers are 
	      			  // included in *.node (0 or 1)
  assert(dimension==3); 
  assert(n_nodes>0);	      
  assert(nAttributes==0);	
  assert(BoundaryMarkers==1);		  
  this->_nodes.clear();
  this->_nodes.resize(n_nodes);
  // Read the nodal coordinates from the node_stream (*.node file).
  unsigned int node_lab=0;
  Point xyz;
  int marker;
  for (unsigned int i=0; i<n_nodes; i++)
    {
      // Check input buffer
      assert (node_stream.good());      
      node_stream >> node_lab  // node number
		  >> xyz(0)    // x-coordinate value
		  >> xyz(1)    // y-coordinate value
		  >> xyz(2)    // z-coordinate value    
		  >> marker;   // node marker  
      // new a node, using tetgen's numbering rule
      Node* newnode = new Node(xyz, node_lab);
      this->_nodes[i]= newnode;
      this->_nodes[i]->set_marker(marker);
      if(marker!=INVALID_UNIT) bc_info.add_node(this->_nodes[i],  marker);
    }
  
  return;
}

//----------------------------------------------------------------------
void Mesh3D::element_in (std::istream& ele_stream)
{
  bool having_marker_9999999 = 0;
  // Check input buffer
  assert (ele_stream.good());

  // Read the elements from the ele_stream (*.ele file). 
  unsigned int element_lab=0, n_nodes_tet=0, nAttri=0;
  unsigned int n_elems=0;  unsigned int tet_marker=0;

  ele_stream >> n_elems     // Read the number of tetrahedrons from the stream.
	           >> n_nodes_tet // Linear or second Tet(defaults to 4).
	           >> nAttri;     // Read the number of attributes from stream.
  assert(n_elems>0);
  assert(n_nodes_tet==4);
  assert(nAttri==1);	     
  this->_elements.clear();
  this->_elements.resize(n_elems);
  
  for (unsigned int i=0; i<n_elems; i++)
  {
    assert (ele_stream.good());      
    // TetGen only supports Tet4 and Tet10 elements.
    this->_elements[i]= new Tet;
    // Read the element label
    ele_stream >> element_lab;
    this->_elements[i]->set_id(element_lab);
    // Read node labels
    for (unsigned int j=0; j<n_nodes_tet; j++)
      {
        unsigned long int node_label;
        ele_stream >> node_label;
        // Assign node to element
        this->_elements[i]->set_node(j,this->_nodes[node_label]);
      }
    // Read attributes from the stream.
    ele_stream >> tet_marker;
    std::vector<Real> attri(nAttri); 
    attri[0]= tet_marker;
    this->_elements[i]->set_attribute(attri);// Only 1 atrri.
    if(tet_marker == 9999999)
      having_marker_9999999 = 1;
  }
  if(!having_marker_9999999){
    std::cout << "Error, the regional marker of air space should be 9999999,\n"
              << "please check the poly file!\n";
    exit(0);
  }
    
  return;
}


//----------------------------------------------------------------------
void Mesh3D::neighbor_face_in(std::istream& nei_stream,
			      std::istream& face_stream)
{
  bool having_marker_1 = 0, having_marker_2 = 0;
  // Check input buffer
  assert (nei_stream.good());
  assert (face_stream.good());
  assert (this->_elements.size()>0);
  const unsigned int n_elems= this->_elements.size(); 

  // Read the elements from the nei_stream (*.neigh file). 
  unsigned int num_elems=0;  unsigned int n_neighs=0;
  unsigned int elem_lab=0;
  nei_stream >> num_elems     // Read the number of tetrahedrons 
	     >> n_neighs;     // Read the number of neighbors(4)
  assert(n_elems==num_elems);
  assert(n_neighs==4);      
  int minius_one_counter=0;
  for (unsigned int i=0; i<num_elems; i++)
    {
      assert (nei_stream.good());      
      Element* tet= this->_elements[i];
      bool have_minus=false;
      int have_minus_neighbr[4];
      nei_stream >> elem_lab;
      assert(elem_lab==i);
      // Read its neighbors' id
      for (unsigned int j=0; j<n_neighs; j++)
	{
	  Element* neighbor= static_cast<Element*>(NULL);
	  int neighbor_label;
	  nei_stream >> neighbor_label;
          have_minus_neighbr[j]= neighbor_label;
	  // Find& Assign neighbor to element
          if(neighbor_label==-1) {
	     minius_one_counter++;
	     have_minus=true;
          }
	  if(neighbor_label!=-1) {
	  	neighbor= this->_elements[neighbor_label];
	  }
	  tet->set_neighbor(j,neighbor);
	}
/*
       if(have_minus) {
         std::cout<<i<<"\t"<<have_minus_neighbr[0]<<"\t"
                  <<have_minus_neighbr[1]<<"\t"
                  <<have_minus_neighbr[2]<<"\t"
                  <<have_minus_neighbr[3]<<"\n";
       }
*/
    }
  // std::cout<<"minius_one_counter  "<<minius_one_counter<<"\n";
  /*
    Give data to BCINFO
  */
  
  // Node: only trifaces with non-zero markers are considered.
  // we first build a map of boundary element, side and triangle nodes list.
  // therefore,the interior boundaries (with marker=0) are not considered.
  // a searching map contains faces of each elements
  std::multimap<std::vector<unsigned int>,
    std::pair<Element*,unsigned int> > element_side_dictionary;
  const unsigned int n_side=4;
  std::vector<Element*>&  elems = this->_elements;
  // loop all elements to full fill this map.
  for(unsigned int elem=0; elem< n_elems; ++elem)  {
      Element* tet= elems[elem];
      for(unsigned int side=0;side<n_side;++side) {
	  Element* tri=tet->build_side(side).release();
	  std::vector<unsigned int> side_node_ids(3);
	  assert(tri->n_nodes()==3);
	  for(unsigned int i=0;i<tri->n_nodes();++i)
	    side_node_ids[i]=tri->get_node(i)->get_id();
	  // now inset the node_id,and element,side into map with sort.
	  std::sort(side_node_ids.begin(),side_node_ids.end());
	  element_side_dictionary.insert(std::make_pair(side_node_ids,
	     std::pair<Element*,unsigned int>(tet,side)));
	  //delete the side_face.
	  delete tri; tri=NULL;
      }
  }
  
  // Read the elements from the face_stream (*.face file). 
  unsigned int num_faces=0;  unsigned int n_markers=0;
  unsigned int face_lab=0; unsigned int face_marker=0;
  face_stream >> num_faces     // Read the number of faces 
	      >> n_markers;     // Read the number of marker(1)
  assert(num_faces>0);
  assert(n_markers==1);    
  const unsigned int M=3; // 3 faces
  std::vector<unsigned int> triface_nodes(3);  
  typedef std::multimap<std::vector<unsigned int>,
	std::pair<Element*,unsigned int> >::iterator MIT;
  BCInfo&          bc_info= this->_bc_info;
  for (unsigned int i=0; i<num_faces; i++)
  {
    assert (face_stream.good());   
    face_stream >> face_lab
          >> triface_nodes[0]
          >> triface_nodes[1]
          >> triface_nodes[2]
          >> face_marker;
    std::sort(triface_nodes.begin(),triface_nodes.end());
    if(face_marker!=0) {
      std::pair<MIT, MIT> tri_iter=
      element_side_dictionary.equal_range(triface_nodes);
      assert(tri_iter.first!=tri_iter.second);
      for(MIT mit= tri_iter.first; mit!= tri_iter.second; mit++) {
        assert((*mit).second.first!=NULL);
        // set the id of side on element.
        bc_info.add_side((*mit).second.first,
                      (*mit).second.second,  face_marker);
      }  
    }
    if(face_marker == 1)
      having_marker_1 = 1;
    else if(face_marker == 2)
      having_marker_2 = 1;
  }
  
  bc_info._initialized = true;
  if(!having_marker_1){
    std::cout << "Error, the facet marker of the air-Earth interface should be 1,\n"
              << "please check the poly file!\n";
    exit(0);
  }
  if(!having_marker_2){
    std::cout << "Error, the facet marker of the computational boundary should be 2,\n"
              << "please check the poly file!\n";
    exit(0);
  }

  return;
}



void Mesh3D::read_tetgen_files()
{
  std::string name_node, name_ele, name_nei, name_face;
  std::string dummy= this->_name_input_file;
  // Check name for *.node or *.ele *.neigh *.face extension.
  name_node= dummy+ std::string(".node");
  name_ele = dummy+ std::string(".ele");
  name_nei=  dummy+ std::string(".neigh");
  name_face= dummy+ std::string(".face");
 
  // Set the streams from which to read in
  std::ifstream node_stream (name_node.c_str());
  std::ifstream ele_stream  (name_ele.c_str());
  std::ifstream nei_stream  (name_nei.c_str());
  std::ifstream face_stream (name_face.c_str());
    
  assert(node_stream.good());
  assert(ele_stream.good());
  assert(nei_stream.good());
  assert(face_stream.good());

  // Skip the comment lines at the beginning
  this->skip_comment_lines (node_stream, '#');
  this->skip_comment_lines (ele_stream, '#');
  this->skip_comment_lines (nei_stream, '#');
  this->skip_comment_lines (face_stream, '#');

  // Read from the streams
  this->node_in (node_stream);
  this->element_in (ele_stream);
  this->neighbor_face_in(nei_stream, face_stream);

  //std::cout << "Mesh3D read in nodes and elements, neighbors and faces " <<std::endl; 

  return;
  
}



void Mesh3D::prepare_mesh_by_reading_files()
{
  //Step 1: clear the mesh member data.
  this->clear();
  //Step 2: read mesh generated by tetgen
  this->read_tetgen_files();

  // Tell we have did the initializing.
  this->_initialized= true;
  /*
    Checking subroutine
  */
  //debugging();
  // Finish.
  return;  

}







