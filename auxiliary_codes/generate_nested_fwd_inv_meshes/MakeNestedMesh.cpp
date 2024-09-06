#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iomanip>
#include <unistd.h>
#include <cmath>
#include "point.h"

int main(int argc, const char **argv)
{
  // Make sure the program is used correctly
  if(argc < 2){
  printf("Usage: %s model_name_without_extension\n", argv[0]);
  return 1;
  }

  // !!!read the regional electric parameter from the 'elec.para' file
  // !define a map from marker to conductivity (sigma), relative_mu (mu_r)  
  //  and relative_epsilon (epsilon_r)    
  std::map<int, std::vector< double > > region_table; 
  typedef std::map<int, std::vector<double> >::iterator IT;
  // the total number of regions shown in the input poly file
  unsigned int n_region;
  bool having_marker_9999999 = 0;
  std::string elec_para_file = "elec.para";
  std::ifstream read_elec_para(elec_para_file.c_str());
  if(!(read_elec_para.good()))
  {
    std::cout << "Input error, didn't find the electrical parameter file"
              << " named 'elec.para', please offer it!\n";
    std::cout << "Program stopped\n";
    exit(1);
  }
  read_elec_para >> n_region;
  // initialize region_table map
  for(unsigned int i = 0; i < n_region; i++)
  {
    int marker;
    std::vector<double> temp(3);
    read_elec_para >> marker;
    for(unsigned int j = 0; j < 3; j++)
      read_elec_para >> temp[j];
    region_table[marker] = temp;

    if(marker == 9999999)
      having_marker_9999999 = 1;
  }
  if(!having_marker_9999999){
    std::cout << "Error, the regional marker of air space should be 9999999,"
              << " please check the poly and elec.para files!\n";
    exit(0);
  }
  else
  {
    if(region_table[9999999][0] > 1.e-10 )
      {
        std::cout << "Input error: the conductivity of air should be larger"
                  << " than 1.e-10, please reset it!\n";
        exit(0);
      }
  }
  read_elec_para.close();

  // !!!definition of golobal variables 
  // name of the generated para file containing the electrical preperty of each 
  // inversion element and other regions
  std::string regional_para_file = "elec_regional_property.para";
  std::ofstream output_regional_para(regional_para_file.c_str());
  std::map <int, int> map_old_marker_2_new_marker;
  typedef std::map<int, int>::iterator IT_int;
  bool having_marker_6666666 = 0;
  std::string model_name(argv[1]);
  std::string line;
  // line counter
  unsigned int k = 0;
  int t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, marker = 0;
  unsigned int t1_temp = 0, tm = 0;
  double tx = 0.0, ty = 0.0, tz = 0.0;
  // used for storing the coordinates of nodes
  std::vector< double > x, y, z;
  x.clear(); y.clear(); z.clear();
  // used for storing the original marker of each tet of inversion mesh
  std::vector< int > original_marker;
  // used for storing the resistivity of each tet of inversion mesh
  std::vector< double > resistivity;
  // !! generate the vtk file contains the info of marker, resistivity
  //    of each inversion tet
  std::ofstream vtk_mesh("resistivity.vtk");
  // Parts 1-2-3 of the vtk file, mandatory 
  vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
              // Header info, Really cool data
    <<"For visulization of ele data\n"
    <<"ASCII\n"; //ASCII data (not BINARY)
  // Part 4 of the vtk file, Geometry/topology including info of node, ele, 
  // type of the ele
  vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n";

  // !!! generating the inversion mesh, electric property file and the vtk file 
  //     about the marker, resistivity of the model
  std::string inv_poly_name = model_name + ".poly";
  std::string command, command1, command2;
  std::cout<< "Please input tetgen command, i.e., ./tetgen -pAnzkKiq1.6a\n";
  std::cin >> command1 >> command2;
  command = command1 + " " + command2 + " " + inv_poly_name;
  // call tetgen to discretize the inversion poly file
  std::cout << "Call tetgen to discretize the inversion poly file ..." << "\n";
  int status = std::system (command.c_str());
  if(status!=0)
    std::cout << "Usage error, please make sure only the model name instead of "
              << "the whole name of the poly file is input!"
              << " Or check the correctness of the poly file and the existence of tetgen!";
  assert(status==0); 

  std::cout << "Generate files of the inversion mesh ...\n";
  std::string old_node_name = model_name + ".1.node";
  std::string old_face_name = model_name + ".1.face";
  std::string old_neigh_name = model_name + ".1.neigh";
  std::string old_vtk_name = model_name + ".1.vtk";
  std::string new_node_name = model_name + "_inv.1.node";
  std::string new_face_name = model_name + "_inv.1.face";
  std::string new_neigh_name = model_name + "_inv.1.neigh";
  std::string new_vtk_name = model_name + "_inv.1.vtk";
  // change the names of .node, .face and .neigh file
  std::rename(old_node_name.c_str(), new_node_name.c_str());
  std::rename(old_face_name.c_str(), new_face_name.c_str());
  std::rename(old_neigh_name.c_str(), new_neigh_name.c_str());
  std::rename(old_vtk_name.c_str(), new_vtk_name.c_str());
  // changing the elemental marker in .ele file
  // name of the generated ele file
  std::string old_ele_file = model_name + ".1.ele";
  std::ifstream input_old_ele(old_ele_file.c_str());
  assert(input_old_ele.good());
  std::string new_ele_name = model_name + "_inv.1.ele";
  std::ofstream output_new_ele(new_ele_name.c_str());
  std::vector<double> temp(3);
  // t1 represents the number of the tetrahedral elements
  input_old_ele >> t1 >> t2 >> t3;
  t1_temp = t1;
  output_new_ele << t1 << '\t' << t2 << '\t' << t3 << '\n';
  while( getline(input_old_ele, line) ){
    // get line number
    k++;
    if(k <= t1_temp){
      input_old_ele >> t1 >> t2 >> t3 >> t4 >> t5 >> marker;
      original_marker.push_back( marker );
      IT it = region_table.find(marker);
      assert(it != region_table.end());
      temp = region_table[marker];
      resistivity.push_back(1. / temp[0]); 
      if( (marker != 9999999) && (marker != 6666666) ){
        tm++;
        map_old_marker_2_new_marker[tm] = marker;
        marker = tm;
      }
      else{
        if(marker == 6666666)
          having_marker_6666666 = 1;
      }
      output_new_ele << t1 << '\t' << t2 << '\t' << t3 << '\t' << t4
                     << '\t' << t5 << '\t' << marker << '\n';
      assert( k == (t1 + 1) );
    }
    // write the comment at the last line of the .ele file
    // to the new .ele file
    if ( line.find("#") == 0 )
      output_new_ele << line;
  }
  std::cout << "The total number of tets in inversion domain is: " 
            << tm << "\n";

  // ! output electrical info of each region into file
  // output the total number of regions
  if(having_marker_6666666)
    output_regional_para << tm + 2 << '\n';
  else
    output_regional_para << tm + 1 << '\n';    
  for(IT_int it = map_old_marker_2_new_marker.begin(); 
            it != map_old_marker_2_new_marker.end(); it++)
  {
    unsigned old_marker = (*it).second;
    output_regional_para << (*it).first << '\t';
    temp = region_table[old_marker];
    for(unsigned int i = 0; i < 3; i++)
      output_regional_para << temp[i] << '\t';
    output_regional_para << '\n';
  }
  if(having_marker_6666666){
    temp = region_table[6666666];
    output_regional_para << 6666666 << '\t';
    for(unsigned int i = 0; i < 3; i++)
      output_regional_para << temp[i] << '\t';
    output_regional_para << '\n';
  }
  temp = region_table[9999999];
  output_regional_para << 9999999 << '\t' << temp[0] << '\t'  
                       << temp[1] << '\t' << temp[2];

  input_old_ele.close();
  output_new_ele.close();
  output_regional_para.close();
  // delete the old .ele file
  command = "rm " + old_ele_file;
  std::system( command.c_str() );

  // !!! generate the volume constraint file for mesh refinement
  //     by using -ra command to generate the nested forward mesh
  std::cout << "Generate the volume constraint file\n";
  // !! for inputting 
  // name of the generated node file
  std::string node_file = model_name + "_inv.1.node";
  std::ifstream input_node(node_file.c_str());
  assert(input_node.good());
  // name of the generated ele file
  std::string ele_file = model_name + "_inv.1.ele";
  std::ifstream input_ele(ele_file.c_str());
  assert(input_ele.good());

  // !! for outputting 
  std::string vol_file = model_name + "_inv.1.vol";
  std::ofstream output_vol(vol_file.c_str());
  // ! read .node file
  // t1 represents the number of the nodes
  input_node >> t1;
  // output point info of the vtk file
  vtk_mesh<<"\nPOINTS\t"<< t1 <<"\tdouble\n";
  t1_temp = t1;
  k = 0;
  while( getline(input_node, line) ){
    // get line number
    k++;
    if(k <= t1_temp){
      input_node >> t1 >> tx >> ty >> tz;
      //storing the node coordinates for calculating the coordinates of the centroids
      //of tetrahedral elements, which will be used in the region part
      x.push_back(tx);
      y.push_back(ty);
      z.push_back(tz);
      // output point info of the vtk file
      vtk_mesh << tx << '\t' << ty << '\t' << tz << '\n';
      assert(k == (t1 + 1));
    }
    // skip the comment at the last line of the .node file
    if ( line.find("#") == 0 )
      continue;
  }
  assert( x.size() == t1_temp );

  // ! write the volume constraint file by using .ele file
  // t1 represents the number of the tetrahedral elements, has been revised
  input_ele >> t1 >> t2 >> t3;
  // output ele info of the vtk
  vtk_mesh << "\nCELLS\t" << t1 << '\t' << 5*t1 << '\n'; 
  output_vol << t1 << '\n';
  t1_temp = t1;
  k = 0;
  assert(t1_temp < 6666666);
  std::cout << "Please input the average frequency, conductivity"
            << " (formal: f sigma):\n";
  double f_aver = 0., sigma_aver = 0.;
  std::cin >> f_aver >> sigma_aver;
  if( (f_aver < 1.0e-14) || (sigma_aver < 1.0e-14)){
    std::cout << "Input error, please type a space between f and sigma!\n";
    std::abort();
  }
  while( getline(input_ele, line) ){
    k++;
    if(k <= t1_temp){
      input_ele >> t1 >> t2 >> t3 >> t4 >> t5 >> marker;
      // output ele info of the vtk
      vtk_mesh << (unsigned int)4 << '\t' << t2 << '\t' << t3 << '\t'
                                          << t4 << '\t' << t5 << '\n';
      // calculate the volume of t1-th element
      Point p0(x[t2], y[t2], z[t2]); Point p1(x[t3], y[t3], z[t3]);
      Point p2(x[t4], y[t4], z[t4]); Point p3(x[t5], y[t5], z[t5]);
      Point a = p1 + (p0 * -1.0);
      Point b = p1 + (p2 * -1.0);
      Point c = p1 + (p3 * -1.0);
      double volume = std::abs( a*(c.cross(b)) )/6.;
      
      // !!! calculate the volume constraints of each region
      double vol_constraint = 0.;
      double sigma_cond = 0., epsilon = 0, mu = 0., omega = 0.; 
      double temp1 = 0., temp2 = 0., lamda = 0., half_lamda = 0.;
      if(marker == 9999999)
        sigma_cond = 1.0e-16;
      else
        sigma_cond = sigma_aver;
      epsilon =  8.854187817 * 1e-12;
      mu = 4.0 * EM::PI * 1e-7;
      assert(std::abs(mu) > 0.);
      omega = 2.0 * EM::PI * f_aver;
      temp1 = std::sqrt(1.0 + std::pow( sigma_cond / (omega*epsilon), 2.0));
      temp2 = std::sqrt(0.5 * mu * epsilon * (temp1 + 1.0));
      lamda = 2.0 * EM::PI / (omega * temp2);
      // From (Ren, 2013), the mesh is generated by enforcing the spatial  
      // constraints of side lengths of 0.5 wavelength in each tetrahedron.
      half_lamda = lamda / 2.0;
      vol_constraint = std::sqrt(2.0) / 12.0 * std::pow(half_lamda, 3.0);

      if(marker == 9999999)
        vol_constraint = 1.22E19;
      else if(volume < vol_constraint){
        vol_constraint = -1.0e+100;
        // for global refinement
        //vol_constraint = volume / divisor;
      }

      output_vol << t1 << "\t" << vol_constraint << "\n";       

      assert( k == (t1 + 1) );
    }
    // skip the comment at the last line of the .ele file
    if (line.find("#") == 0)
      continue;
  }
  // output ele info of the vtk, the element type
  vtk_mesh<< "\nCELL_TYPES\t" << t1_temp <<"\n";
  for(unsigned int i = 0; i < t1_temp; i++){
    // 10--tet
    vtk_mesh << (unsigned int)10 << '\n'; 
    // for more details please check the figure 2 and 3 in vtk format file
  }
  // Part 5 of the vtk file, dataset attributes
  vtk_mesh << "\nCELL_DATA\t"<< t1_temp <<"\n" << "SCALARS tet_marker int 1\n"
           << "LOOKUP_TABLE default\n";
  for(unsigned int i = 0; i < t1_temp; i++) {
    vtk_mesh<< original_marker[i] <<"\n";
  }
  vtk_mesh << "\nSCALARS elem_data double 1\n" << "LOOKUP_TABLE default\n";
  for(unsigned int i = 0; i < t1_temp; i++){
    vtk_mesh << resistivity[i] << '\n';
  }
  input_node.close();
  input_ele.close();
  output_vol.close();
  vtk_mesh.close();

  // !! re-generate the forward mesh by using global refinement
  std::cout << "Call tetgen to re-generate forward mesh:" << "\n";
  std::cout << "Please input tetgen command, i.e.,\n"
            << "./tetgen -AnzkKraiq1.4 or\n"
            << "./tetgen -AnzkKriq1.2 (without lambda constraints, recommended!)\n";
  std::cin >> command1 >> command2;
  std::string mesh_name =  model_name + "_inv.1";
  command = command1 + " " + command2 + " " + mesh_name;
    status = std::system(command.c_str());
  if(status!=0)
    std::cout << " Please check the existence of input tetgen!";
    // exit normally
    assert(status==0);
  // rename newly generated .node, .face, .neigh, .ele and .vtk files
  // old names
  std::string tot_refined_mesh_name =  model_name + "_inv.2";
  old_node_name = tot_refined_mesh_name + ".node";
  old_face_name = tot_refined_mesh_name + ".face";
  old_neigh_name = tot_refined_mesh_name + ".neigh";
  old_vtk_name = tot_refined_mesh_name + ".vtk"; 
  std::string old_ele_name = tot_refined_mesh_name + ".ele"; 
  // new names
  std::string new_fwd_mesh_name = model_name + "_fwd";
  new_node_name = new_fwd_mesh_name + ".1.node";
  new_face_name = new_fwd_mesh_name + ".1.face";
  new_ele_name = new_fwd_mesh_name + ".1.ele";
  new_neigh_name = new_fwd_mesh_name + ".1.neigh";
  new_vtk_name = new_fwd_mesh_name + ".1.vtk";  
  // change the file names
  std::rename(old_node_name.c_str(), new_node_name.c_str());
  std::rename(old_face_name.c_str(), new_face_name.c_str());
  std::rename(old_ele_name.c_str(), new_ele_name.c_str());
  std::rename(old_neigh_name.c_str(), new_neigh_name.c_str());
  std::rename(old_vtk_name.c_str(), new_vtk_name.c_str());

  return 0;
}




