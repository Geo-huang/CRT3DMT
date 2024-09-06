/******************************************************************************
 *    Function       : Adaptively adjust the inversion mesh                   *   
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2020                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2020.07.07                                             *
 *    Last revision  : 2023.10.12                                             *
 ******************************************************************************/
#include "adaptive.h"
 // used for generating screen output log
 extern std::ofstream screen_output;

Adaptive::Adaptive(Startup _startup, Data_3D _data_3d, 
                   Init_Prior_Mod_3D _init_prior_mod_3d, Inv_Para _inv_para,
                   Fwd_Para_3D& _fwd_para_3d, Mesh3D& _MT3D_inv_mesh)
{
  // to be checked 
  this->startup           = _startup;
  this->inv_para          = _inv_para;
  this->data_3d           = _data_3d;
  this->init_prior_mod_3d = _init_prior_mod_3d;
  // _fwd_para_3d is the reference of the global variable fwd_para_3d
  // pointer fwd_para_3d points to the global variable fwd_para_3d
  this->fwd_para_3d       = &_fwd_para_3d;
  this->init_fwd_mod_name = _fwd_para_3d.starting_model[0];

  if(this->startup.inv_algorithm == "GN"){
    std::cout << "This module is under construction\n";
    exit(1);
  }
  else if(this->startup.inv_algorithm == "L-BFGS")
    // call defaulted implicitly-declared copy assignment operator
    this->LBFGS = LBFGS_Inv(this->startup, this->data_3d, this->init_prior_mod_3d,
                    this->inv_para, *(this->fwd_para_3d), _MT3D_inv_mesh);
  else{
    std::cout<< "The input inversion algorithm is incorrect, please input"
              << " L-BFGS in .startup file!\n";
    screen_output<< "The input inversion algorithm is incorrect, please input"
                 << " L-BFGS in .startup file!\n";
    exit(1);
  }         
}

template <typename T> 
std::vector<size_t>  Adaptive::sort_indexes(const std::vector<T>& v)
{
  // initialize original index locations
  std::vector<size_t>  idx(v.size());
  for (size_t i=0; i!=idx.size(); ++i) 
    idx[i] = i;

  // sort the element index according to the relative 
  // value of the element
  // lambda expressions are a C++11 feature
  sort(idx.begin(), idx.end(),
      [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

void Adaptive::Do_Adapt_Inv()
{
  this->k = 0;
  this->total_iter_number = 0;
  std::string unkowns_filename = this->startup.proj_name + ".unknowns";
  std::ofstream out_unknowns(unkowns_filename.c_str());
  std::string fwd_unkowns_filename = this->startup.proj_name + ".fwd_unknowns";
  std::ofstream out_fwd_unknowns(fwd_unkowns_filename.c_str());

  if(this->startup.inv_algorithm == "GN")
  {
    std::cout << "This module is under construction\n";
    exit(1);
  } 
  else if(this->startup.inv_algorithm == "L-BFGS")
  {
    // performing L-BFGS inversion on the initial inversion mesh
    this->LBFGS.Do_LBFGS_Inv();
    this->LBFGS.get_RMS_last( this->RMS_last );
    this->LBFGS.get_lambda_last( this->lambda_last );
    this->LBFGS.get_iter_number( this->total_iter_number );

    if(this->startup.data_method == "MT1D"){
      std::cout<< "This module is to be constructed." << std::endl;
      exit(1);
    }
    else if(this->startup.data_method == "MT3D"){
      // get the n_dofs of forward modelling
      std::stringstream counter;
      counter << (*this->fwd_para_3d).starting_counter;
      std::string s = counter.str();
      std::string model_name = (*this->fwd_para_3d).starting_model[0]+std::string(".")+s;
      Mesh3D mesh(model_name); 
      mesh.prepare_mesh_by_reading_files();
      std::cout << "The number of the dofs of the initial forward "
                << "mesh is: " << this->Get_n_dofs(mesh) << std::endl;
      screen_output << "The number of the dofs of the initial forward "
                    << "mesh is: " << this->Get_n_dofs(mesh) << std::endl;
      out_fwd_unknowns << this->Get_n_dofs(mesh) << "\n";
                       
      std::cout << "The number of the free parameters of the initial inversion "
                << "mesh is: " << (this->init_prior_mod_3d.m_ini).size() 
                << std::endl; 
      screen_output << "The number of the free parameters of the initial inversion "
                    << "mesh is: " << (this->init_prior_mod_3d.m_ini).size() 
                    << std::endl; 
      out_unknowns << (this->init_prior_mod_3d.m_ini).size() << "\n";
      
      std::cout<< "\n\nAdjust inversion mesh:";
      screen_output<< "\n\nAdjust inversion mesh:";
    }

    // !!! for refining the initial inversion mesh
    // calculate the gradient of m in logarithmic space for MT3D
    this->LBFGS.Inv_Basis::Comp_Abs_Gra_M(this->abs_gra_m);
    if(this->inv_para.adap_method == 1){
      // calclulate the gradient of \phi_d
      this->LBFGS.Inv_Basis::Comp_Abs_Gra_Phi_D(this->abs_gra_phi_d);
    }
    else if(this->inv_para.adap_method == 2){
      std::cout << "This module is under construction\n";
      exit(1);
    }
    else if(this->inv_para.adap_method == 3){
      std::cout << "This module is under construction\n";
      exit(1);
    }
    else{
      std::cout << '\n' 
                << "The input adaptive method is incorrect!"
                << " Please check the cooresponding part in the .inv file!" 
                << '\n' << '\n';
      screen_output << '\n' 
                    << "The input adaptive method is incorrect!"
                    << " Please check the cooresponding part in the .inv file!" 
                    << '\n' << '\n';
      exit(1);
    }
    // get the last updated model on the old mesh, which
    // will be used for assignment of the intial model 
    // and/or the prior model on the updated mesh      
    this->LBFGS.Inv_Basis::get_m_last(this->m_inv);
    // adjust the inversion mesh and update the L-BFGS related data members
    this->k++;
    this->Adjust_Inv_Mesh();

    while(true)
    {
      // Assign the object of class LBFGS_Inv again
      if(this->startup.data_method == "MT1D"){
        std::cout<< "This module is to be constructed." << std::endl;
        exit(1);
      }
      else if(this->startup.data_method == "MT3D"){
        // define a Mesh3D object
        std::stringstream inv_counter;
        inv_counter << this->init_prior_mod_3d.counter;
        std::string inv_s = inv_counter.str();
        std::string inv_mesh_name = this->init_prior_mod_3d.model_name
                                    + std::string(".") + inv_s;
        Mesh3D MT3D_inv_mesh(inv_mesh_name);
        MT3D_inv_mesh.prepare_mesh_by_reading_files();
        assert(MT3D_inv_mesh.n_elems() > 0);
        // update object this->LBFGS
        this->LBFGS = LBFGS_Inv(this->startup, this->data_3d, this->init_prior_mod_3d,
                                this->inv_para, *(this->fwd_para_3d),
                                MT3D_inv_mesh, this->k);
        // update the initial lambda for the new L-BFGS iterations
        //this->LBFGS.update_initial_lambda( this->k );
        this->LBFGS.update_initial_lambda(this->lambda_last);
        this->M = (this->init_prior_mod_3d.m_ini).size();
        if( (this->k == this->inv_para.max_adjust_times) || 
             (this->M >= this->inv_para.max_No_inv_unknowns) ){
          this->left_iter_number = this->inv_para.tol_iter_times
                                 - this-> total_iter_number;
          this->LBFGS.update_LBFGS_iter_times(this->left_iter_number);
        }

        std::cout << "\nNew L-BFGS iteration is in progress:\n"; 
        screen_output << "\nNew L-BFGS iteration is in progress:\n"; 
        // NOTE:  FOR THE FOLLOWING OPERATION, WE NEED TO REUSE LOCAL 
        // VARIABLE MT3D_inv_mesh!
        this->LBFGS.Do_LBFGS_Inv();
        
        // get the n_dofs of forward modelling
        std::stringstream counter;
        counter << (*this->fwd_para_3d).starting_counter;
        std::string s = counter.str();
        std::string model_name = (*this->fwd_para_3d).starting_model[0]+std::string(".")+s;
        Mesh3D mesh(model_name); 
        mesh.prepare_mesh_by_reading_files();
        std::cout << "The number of the dofs of the present forward "
                  << "mesh is: " << this->Get_n_dofs(mesh) << std::endl;
        screen_output << "The number of the dofs of the present forward "
                      << "mesh is: " << this->Get_n_dofs(mesh) << std::endl;
        out_fwd_unknowns << this->Get_n_dofs(mesh) << "\n";
                         
        std::cout << "The number of the free parameters of the present inversion "
                  << "mesh is: " << (this->init_prior_mod_3d.m_ini).size() 
                  << std::endl;
        screen_output << "The number of the free parameters of the present inversion "
                      << "mesh is: " << (this->init_prior_mod_3d.m_ini).size() 
                      << std::endl;
        out_unknowns << (this->init_prior_mod_3d.m_ini).size() << "\n";

        // calculating the gradient of m
        this->LBFGS.get_RMS_last( this->RMS_last );
        this->LBFGS.get_lambda_last( this->lambda_last );
        unsigned int n_iter_temp = 0;
        this->LBFGS.get_iter_number( n_iter_temp );
        this->total_iter_number += n_iter_temp;
        this->k++;

        if(this->k == (this->inv_para.max_adjust_times + 1)){
          std::cout << "Inversion mesh adjustment stopped due to"
                    << " tol of max adjustment times\n"; 
          std::cout << "The times of the inversion mesh adjustment is: " 
                    << this->k - 1 << std::endl;
          screen_output << "Inversion mesh adjustment stopped due to"
                        << " tol of max adjustment times\n"; 
          screen_output << "The times of the inversion mesh adjustment is: " 
                        << this->k - 1 << std::endl;
          break;
        }
        else if( this->M >= this->inv_para.max_No_inv_unknowns ){
          std::cout << "Inversion mesh adjustment stopped due to tol of inv unkowns\n";
          std::cout << "The times of the inversion mesh adjustment is: " 
                    << this->k - 1 << std::endl;
          screen_output << "Inversion mesh adjustment stopped due to tol of inv unkowns\n";
          screen_output << "The times of the inversion mesh adjustment is: " 
                        << this->k - 1 << std::endl;
          break;
        }
        else if( !(this->RMS_last > this->inv_para.tol_rms) ){
          std::cout << "Inversion mesh adjustment stopped due to tol of RMS\n";
          std::cout << "The times of the inversion mesh adjustment is: " 
                    << this->k - 1 << std::endl;
          screen_output << "Inversion mesh adjustment stopped due to tol of RMS\n";
          screen_output << "The times of the inversion mesh adjustment is: " 
                        << this->k - 1 << std::endl;
        break;
        }
        else if( !(this->lambda_last > this->inv_para.tol_lambda) ){
          std::cout << "Inversion mesh adjustment stopped due to tol of lambda\n";
          std::cout << "The times of the inversion mesh adjustment is: " 
                    << this->k - 1 << std::endl;
          screen_output << "Inversion mesh adjustment stopped due to tol of lambda\n";
          screen_output << "The times of the inversion mesh adjustment is: " 
                        << this->k - 1 << std::endl;
          break;
        }
        // NOTE:  FOR THE FOLLOWING OPERATION, WE NEED TO REUSE LOCAL 
        // VARIABLE MT3D_inv_mesh!
        // calculate the gradient of m in logarithmic space for MT3D
        std::cout<< "\n\nAdjust inversion mesh:";
        screen_output<< "\n\nAdjust inversion mesh:";
        this->LBFGS.Inv_Basis::Comp_Abs_Gra_M(this->abs_gra_m);
      }

      // !!! refine subsequent inversion mesh
      if(this->inv_para.adap_method == 1){
        // calclulate the gradient of \phi_d
        this->LBFGS.Inv_Basis::Comp_Abs_Gra_Phi_D(this->abs_gra_phi_d);
      }
      else if(this->inv_para.adap_method == 2){
        std::cout << "This module is under construction\n";
        exit(1);
      }
      else if(this->inv_para.adap_method == 3){
        std::cout << "This module is under construction\n";
        exit(1);
      }
      else{
        std::cout<< '\n' 
                  << "The input adaptive method is incorrect!"
                  << " Please check the cooresponding part in the .inv file!" 
                  << '\n' << '\n';
        screen_output<< '\n' 
                     << "The input adaptive method is incorrect!"
                     << " Please check the cooresponding part in the .inv file!" 
                     << '\n' << '\n';
        exit(1);
      }
      // get the last updated model on the old mesh, which
      // will be used for assignment of the intial model 
      // and/or the prior model on the updated mesh      
      this->LBFGS.Inv_Basis::get_m_last(this->m_inv);
      // adjust the inversion mesh and update the L-BFGS related data members
      this->Adjust_Inv_Mesh();
    }// end while
  }
  else
  {
    std::cout<< "The input inversion algorithm is incorrect, please input"
              << " L-BFGS in .startup file!\n";
    screen_output<< "The input inversion algorithm is incorrect, please input"
                 << " L-BFGS in .startup file!\n";
    exit(1);
  }
  
  return;
}

void Adaptive::Adjust_Inv_Mesh()
{
  if(this->startup.data_method == "MT1D"){
    std::cout << "The function Adjust_Inv_Mesh of MT1D is under construction."
              << std::endl;
    exit(1);
  }
  else if(this->startup.data_method == "MT2D"){
    std::cout << "The function Adjust_Inv_Mesh of MT2D is under construction."
              << std::endl;
    exit(1);
  }
  else{ // for MT3D
    // !!! refine 3D inversion mesh
    this->M = (this->init_prior_mod_3d.m_ini).size();
    this->adjust_indicator.resize(this->M);
    if(this->inv_para.adap_method == 1){
    // this->adjust_indicator == abs_gra_phi_d, but in different data structure
      for (unsigned int i = 0; i < this->M; i++)
        this->adjust_indicator[i] = this->abs_gra_phi_d(i);
    }
    else if(this->inv_para.adap_method == 2){
      std::cout << "This module is under construction\n";
      exit(1);
    }
    else {
      std::cout << "This module is under construction\n";
      exit(1);
    }
    // ensure input refinement-control parameter is right
    assert( (this->inv_para.frac_m_r) && (this->inv_para.frac_m_r < 1.0) );

    /*
    // !!! do refinement by using the fixed fraction of the largest error
    // maximum component of this->adjust_indicator vector
    double max_adjust_indicator = 0.0;
    max_adjust_indicator = *std::max_element( this->adjust_indicator.begin(),
                                              this->adjust_indicator.end() );
    // maximum component of this->abs_gra_m vector
    double max_abs_gra_m = 0.0;
    // store the absolute value of the gradient m_k by using another data type
    std::vector<double> abs_grad_m( this->abs_gra_m.size() );
    for(unsigned int i = 0; i < abs_grad_m.size(); i++)
      abs_grad_m[i] = this->abs_gra_m(i);
    max_abs_gra_m = *std::max_element( abs_grad_m.begin(), abs_grad_m.end() );

    // read the present inversion mesh
    std::stringstream counter;
    counter << this->init_prior_mod_3d.counter;
    std::string s = counter.str();
    std::string inv_mesh_name = this->init_prior_mod_3d.model_name
                                + std::string(".") + s;
    Mesh3D Present_inv_mesh(inv_mesh_name);
    Present_inv_mesh.prepare_mesh_by_reading_files();
    unsigned int tot_tet_number = Present_inv_mesh.n_elems();
    assert( tot_tet_number > this->M );
    // used for generating the volume constrained file .vol file
    std::vector< double > ele_constrained_vol(tot_tet_number, -1.0E100);
    // stores the id of the tetrahedral element to be refined
    std::vector< unsigned int > id_records_of_refined_ele;
    for(unsigned int e = 0; e < tot_tet_number; e++){
      Tet* inv_tet = static_cast< Tet* >( &(Present_inv_mesh.get_elem(e)) );
      assert(e == inv_tet->get_id()); 
      unsigned int inv_tet_marker = inv_tet->get_tet_marker();

      // only refine the tetrahedral elements in the inversion domain,
      // where the parameters are free
      if( (inv_tet_marker != 9999999) && (inv_tet_marker != 6666666) ){
        //local id of the tet element in the inversion domain
        unsigned int tet_index = inv_tet_marker - 1;
        if( (this->abs_gra_m(tet_index) > this->inv_para.frac_m_r * max_abs_gra_m) 
             && (this->adjust_indicator[tet_index] > this->inv_para.frac_others_r 
                 * max_adjust_indicator) ){ 
          double tet_volume = inv_tet->get_size();
          ele_constrained_vol[e] = tet_volume / 2.0;
          id_records_of_refined_ele.push_back(e);
        }
      }
    }
    std::cout << id_records_of_refined_ele.size()
              << " marked inv elements for next refinement\n";
    // assignment of vector ele_constrained_vol was finished!
    */

    // !!! refinement by using the fixed fraction of the number of the inversion tets
    // store the absolute value of the gradient m_k by using another data type
    std::vector<double> abs_grad_m( this->abs_gra_m.size() );
    for(unsigned int i = 0; i < abs_grad_m.size(); i++)
      abs_grad_m[i] = this->abs_gra_m(i);
    // the planed number of elements needed to be refined by using indicator of 
    // the gradient of model
    unsigned int num_ele_grad_m = this->M * this->inv_para.frac_m_r;
    // the planed number of elements needed to be refined by using other indicator
    unsigned int num_ele_other  = this->M * this->inv_para.frac_others_r;
    std::vector<size_t> index_grad_m = this->sort_indexes(abs_grad_m);
    std::vector<size_t> index_other = this->sort_indexes(this->adjust_indicator);
    // read the present inversion mesh
    std::stringstream counter;
    counter << this->init_prior_mod_3d.counter;
    std::string s = counter.str();
    std::string inv_mesh_name = this->init_prior_mod_3d.model_name
                                + std::string(".") + s;
    Mesh3D Present_inv_mesh(inv_mesh_name);
    Present_inv_mesh.prepare_mesh_by_reading_files();
    unsigned int tot_tet_number = Present_inv_mesh.n_elems();
    assert( tot_tet_number > this->M );
    // used for generating the volume constrained file .vol file
    std::vector< double > ele_constrained_vol(tot_tet_number, -1.0E100);
    // !! refine by using estimator gradient of m
    for(unsigned int i = this->M - num_ele_grad_m; i < this->M; i++){
      unsigned int index = index_grad_m[i];
      // get the global id of the index-th tet in the inversion domain
      unsigned int global_id = this->init_prior_mod_3d.m_ini_id[index];
      // define the object of the global_id-th tetrahedral element
      Tet* inv_tet = static_cast< Tet* >( &(Present_inv_mesh.get_elem(global_id)) );
      // get the volume of the tetrahedral element i in inversion domain
      double tet_volume = inv_tet->get_size();
      ele_constrained_vol[global_id] = tet_volume / 2.0;
    }
    // !! undo the unnecessary refinements by using other indicator, i.e. 
    // this->adjust_indicator
    for(unsigned int j = 0; j < this->M - num_ele_other; j++){
      unsigned int index = index_other[j];
      // get the global id of the index-th tet in the inversion domain
      unsigned int global_id = this->init_prior_mod_3d.m_ini_id[index];
      ele_constrained_vol[global_id] = -1.0E100;
    }
    unsigned int counter_refined_ele = 0;
    for(unsigned int k = 0; k < ele_constrained_vol.size(); k++){
      double vol = ele_constrained_vol[k];
      if(std::abs(vol) < 1.0E100)
        counter_refined_ele++;
    }
    std::cout << counter_refined_ele
              << " marked inv elements for next refinement\n";
    screen_output << counter_refined_ele
                  << " marked inv elements for next refinement\n";
    // assignment of vector ele_constrained_vol was finished!

    // write new volume constraints into .vol file
    std::string vol_name = this->init_prior_mod_3d.model_name
                           + std::string(".") + s + std::string(".vol");
    std::ofstream output_vol(vol_name.c_str());
    assert(output_vol.good());
    // <# of tetrahedra>
    output_vol << tot_tet_number << "\n";
    for(unsigned int e = 0; e < tot_tet_number; e++){
      // <tetrahedron #> <maximum volume>
      output_vol << e << '\t' << ele_constrained_vol[e] << '\n';
    }
    output_vol.close();

    std::cout << "Call tetgen to re-generate inversion mesh:" << "\n";
    screen_output << "Call tetgen to re-generate inversion mesh:" << "\n";
    // -q1.4 is necessary to avoid generating the bad grid?
    // which one is better -q1.4 or -q1.6?
    std::string command("./tetgen -Anzarq1.6KQ ");
    command = command + inv_mesh_name;
    int status = std::system(command.c_str());
    if(status!=0)
    {
      std::string command1, command2;
      std::cout << "\nPlease put tetgen program into"
                << " the current work directory first!\n";
      std::cout << "Then input the tetgen command,"
                << " i.e., ./tetgen -Anzarq1.6KQ\n";
      std::cin >> command1 >> command2;
      command = command1 + " " + command2 + " " + inv_mesh_name;
      int status = std::system(command.c_str());
      // exit normally
      assert(status==0);
    }

    // preparation for initialization of variable this->init_prior_mod_3d
    (this->init_prior_mod_3d.sub_node_id).clear();
    (this->init_prior_mod_3d.m_ini).clear();
    (this->init_prior_mod_3d.m_ini_id).clear();
    (this->init_prior_mod_3d.fwd_tet_id_in_mj).clear();

    // variable definitions for mesh operation
    std::string line;
    // line counter
    unsigned int l = 0;
    int t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, marker = 0;
    unsigned int t1_temp = 0, tm = 0;
    double tx = 0.0, ty = 0.0, tz = 0.0;
    // used for storing the coordinates of nodes
    std::vector< double > x, y, z;
    x.clear(); y.clear(); z.clear();

    // !!! rewrite the info of elemental marker in .ele file
    std::stringstream counter_p_1;
    counter_p_1 << (this->init_prior_mod_3d.counter + 1);
    std::string s_p_1 = counter_p_1.str();
    std::string ele_file_name = this->init_prior_mod_3d.model_name
                              + std::string(".") + s_p_1 + std::string(".ele");
    std::string re_ele_file_name = this->init_prior_mod_3d.model_name
                        + std::string(".") + s_p_1 + std::string("_temp.ele");
    std::rename(ele_file_name.c_str(), re_ele_file_name.c_str());
    std::ifstream input_old_ele(re_ele_file_name.c_str());
    assert(input_old_ele.good());
    std::ofstream output_new_ele(ele_file_name.c_str());
    input_old_ele >> t1 >> t2 >> t3;
    t1_temp = t1;
    output_new_ele << t1 << '\t' << t2 << '\t' << t3 << '\n';
    // used for checking tet_marker is a integer or not
    double tet_marker = 0.0;
    while( getline(input_old_ele, line) ){
      // get line number
      l++;
      if(l <= t1_temp){
        input_old_ele >> t1 >> t2 >> t3 >> t4 >> t5 >> tet_marker;
        // ensure tet_marker is a integer
        assert(std::abs(int(tet_marker) - tet_marker) < EM::TOL);
        marker = tet_marker;
        if( (marker != 9999999) && (marker != 6666666) ){
          // ensure no extra marker was produced by mesh refinement 
          assert( (marker > 0) && (marker < this->init_prior_mod_3d.n_regions) );
          // update member variable init_prior_mod_3d
          double nomin = 0, denomin = 0., linear_sigma = 0.;
          nomin = this->inv_para.a + ( this->inv_para.b 
                  * std::exp( this->inv_para.n * this->m_inv( marker - 1 ) ) );
          denomin = 1.0 + std::exp( this->inv_para.n * this->m_inv( marker - 1 ) );
          // conductivity in lienar space
          linear_sigma = nomin / denomin;
          if( std::isnan( linear_sigma ) )
            // we set sigma = b (upper bound), if it is equals to NAN
            linear_sigma =  this->inv_para.b;

          // update m_ini
          (this->init_prior_mod_3d.m_ini).push_back( linear_sigma );
          // update m_ini_id
          (this->init_prior_mod_3d.m_ini_id).push_back( t1 );
          // update region_table_ini
          std::vector< double > para(3);
          para = this->init_prior_mod_3d.region_table_ini[marker];
          tm++;
          marker = tm;
          assert(this->init_prior_mod_3d.m_ini.size() == marker);
          // mu_r and epsilon_r are kept constant
          para[0] = linear_sigma;
          this->init_prior_mod_3d.region_table_ini[marker] = para;
          // update sub_node_id
          (this->init_prior_mod_3d.sub_node_id).insert(t2);
          (this->init_prior_mod_3d.sub_node_id).insert(t3);
          (this->init_prior_mod_3d.sub_node_id).insert(t4);
          (this->init_prior_mod_3d.sub_node_id).insert(t5);
        }
        output_new_ele << t1 << '\t' << t2 << '\t' << t3 << '\t' << t4
                       << '\t' << t5 << '\t' << marker << '\n';
        assert( l == (t1 + 1) );
      }
      // write the comment at the last line of the .ele file
      // to the new .ele file
      if ( line.find("#") == 0 )
       output_new_ele << line;
    }
    input_old_ele.close();
    output_new_ele.close();
    // delete the old .ele file
    command = "rm " + re_ele_file_name;
    std::system( command.c_str() );

    // !!! re-generate forward mesh
    // !! generating the volume constraint file for mesh refinement
    //    by using -ra command to generate the nested forward mesh
    std::cout << "Generate the volume constraint file for local mesh refinement\n"
              << "by using -ra command to generate the nested fwd mesh\n";
    screen_output << "Generate the volume constraint file for local mesh refinement\n"
                  << "by using -ra command to generate the nested fwd mesh\n";
    // !! for inputting 
    // name of the generated node file
    // name of the refined inversion node file
    std::string node_file = this->init_prior_mod_3d.model_name 
                            + std::string(".") + s_p_1 + std::string(".node"); 
    std::ifstream input_node(node_file.c_str());
    assert(input_node.good());
    // name of the refined inversion ele file
    std::string ele_file = this->init_prior_mod_3d.model_name 
                           + std::string(".") + s_p_1 + std::string(".ele");
    std::ifstream input_ele(ele_file.c_str());
    assert(input_ele.good());

    // ! for outputting
    std::string vol_file = this->init_prior_mod_3d.model_name
                           + std::string(".") + s_p_1 + std::string(".vol");
    std::ofstream output_volume(vol_file.c_str());
    
    // ! read .node file
    // t1 represents the number of the nodes
    input_node >> t1;
    t1_temp = t1;
    l = 0;
    while( getline(input_node, line) ){
      // get line number
      l++;
      if(l <= t1_temp){
        input_node >> t1 >> tx >> ty >> tz;
        //storing the node coordinates for calculating the coordinates of the centroids
        //of tetrahedral elements, which will be used in the region part
        x.push_back(tx);
        y.push_back(ty);
        z.push_back(tz);
        assert(l == (t1 + 1));
      }
      // skip the comment at the last line of the .node file
      if ( line.find("#") == 0 )
        continue;
    }
    assert( x.size() == t1_temp );

    // ! write the volume constraint file by using .ele file
    // t1 represents the number of the tetrahedral elements
    input_ele >> t1;
    output_volume << t1 << '\n';
    t1_temp = t1; 
    if(t1_temp >= 6666666){
      std::cout << "The refined inversion mesh is too dense!" << std::endl;
      screen_output << "The refined inversion mesh is too dense!" << std::endl;
      std::abort();
    }
    l = 0;
    bool have_marker_6666666_or_not = false;
    while( getline(input_ele, line) ){
      l++;
      if(l <= t1_temp){
        input_ele >> t1 >> t2 >> t3 >> t4 >> t5 >> marker;
        if(marker == 6666666)
          have_marker_6666666_or_not = true;
          
        // calculate the volume of t1-th element
        Point p0(x[t2], y[t2], z[t2]); Point p1(x[t3], y[t3], z[t3]);
        Point p2(x[t4], y[t4], z[t4]); Point p3(x[t5], y[t5], z[t5]);
        Point a = p1 + (p0 * -1.0);
        Point b = p1 + (p2 * -1.0);
        Point c = p1 + (p3 * -1.0);
        double volume = std::abs( a*(c.cross(b)) )/6.;
        
        // !!! calculate the volume constraints of each region
        // f_highest represents the highest frequency in log10 domain
        double vol_constraint = 0., f_sum = 0., f_aver = 0.;
        double sigma_cond = 0., epsilon = 0, mu = 0., omega = 0.; 
        double temp1 = 0., temp2 = 0., lamda = 0., half_lamda = 0.;
        std::vector<double> parameter;

        // ! calculate the average frequency in log10 domain
        for(unsigned int i = 0; i < this->data_3d.n_f; i++){
          f_sum += std::log10(this->data_3d.f[i]);
        }
        f_aver = std::pow(10., f_sum / this->data_3d.n_f);

        /*
        f_aver = *std::max_element( (this->data_3d.f).begin(), 
                                    (this->data_3d.f).end() );
        */
        parameter = this->init_prior_mod_3d.region_table_ini[marker];
        sigma_cond = parameter[0];
        epsilon = parameter[1] * EM::EPSILON0;
        if( std::abs(epsilon) < EM::TOL)
          epsilon =  EM::EPSILON0;
        mu = parameter[2] * EM::MU0;
        assert(std::abs(mu) > 0.);
        omega = 2.0 * EM::PI * f_aver;
        temp1 = std::sqrt(1.0 + std::pow( sigma_cond / (omega*epsilon), 2.0));
        temp2 = std::sqrt(0.5 * mu * epsilon * (temp1 + 1.0));
        lamda = 2.0 * EM::PI / (omega * temp2);
        // From (Ren, 2013), the mesh is generated by enforcing the spatial  
        // constraints of side lengths of 0.5 wavelength in each tetrahedron.
        half_lamda = lamda / 2.0;
        vol_constraint = std::sqrt(2.0) / 12.0 * std::pow(half_lamda, 3.0);

        if(marker == 9999999){
          vol_constraint = -1.0e+100;
          //vol_constraint = 1.22e19;
        }
        else if(volume < vol_constraint){
          vol_constraint = -1.0e+100;
          // for global refinement
          //vol_constraint = volume / 2.0;
        }

        output_volume << t1 << "\t" << vol_constraint << "\n";

        assert( l == (t1 + 1) );
      }
      // skip the comment at the last line of the .ele file
      if (line.find("#") == 0)
        continue;
    }
    input_node.close();
    input_ele.close();
    output_volume.close();
    std::cout << "Done, the new volume constraint file has been generated!\n\n";
    screen_output << "Done, the new volume constraint file has been generated!\n\n";

    // !! re-generate the forward mesh by using local refinement
    std::string refined_mesh_name = this->init_prior_mod_3d.model_name
                                    + std::string(".") + s_p_1;
    // ! create a new additional-node file and copy the content in the old one 
    //   into it for local refinement at the observing sites when generating
    //   the forward mesh
    // create a new node file
    std::string new_node_file_name = refined_mesh_name + ".a.node";
    std::ofstream create_a_new_node_file(new_node_file_name.c_str());
    assert(create_a_new_node_file.good());
    std::string old_node_file_name = inv_mesh_name + ".a.node";
    // copy the content in the old node file into the new node file
    std::string cp_command = "cp -f " + old_node_file_name + " "
                                      + new_node_file_name;
    std::system(cp_command.c_str());
                                    
    std::cout << "Call tetgen to re-generate forward mesh:" << "\n";
    screen_output << "Call tetgen to re-generate forward mesh:" << "\n";
    // ! -i switch is used to insert additional nodes for local refinement 
    //   at the observing sites
    // do local refinement with only local nodes, -q1.4 is used to improve mesh quality
    command = std::string("./tetgen -Anzrq1.4Ki  ");
    // do local refinement with only wavelength-based volume constraints
    //command = std::string("./tetgen -Anzrq1.4Ka  ");
    // do local refinement with local nodes and wavelength-based volume constraints
    //command = std::string("./tetgen -Anzrq1.4Kai  ");

    command = command + refined_mesh_name;
    status = std::system(command.c_str());
    // exit normally
    assert(status==0);

    // ! rename newly generated .node, .face, .neigh, .ele and .vtk files
    // the old names
    std::stringstream counter_p_2;
    counter_p_2 << (this->init_prior_mod_3d.counter + 2);
    std::string s_p_2 = counter_p_2.str();
    std::string tot_refined_mesh_name = this->init_prior_mod_3d.model_name
                                       + std::string(".") + s_p_2;
    std::string old_node_name = tot_refined_mesh_name + ".node";
    std::string old_face_name = tot_refined_mesh_name + ".face";
    std::string old_ele_name = tot_refined_mesh_name + ".ele"; 
    std::string old_neigh_name = tot_refined_mesh_name + ".neigh";
    std::string old_vtk_name = tot_refined_mesh_name + ".vtk";  
    // the new names
    std::stringstream counter_fwd;
    counter_fwd << (this->k + 1);
    std::string s_fwd = counter_fwd.str();
    std::string new_fwd_mesh_name = this->init_fwd_mod_name
                                    + std::string(".") + s_fwd;
    std::string new_node_name = new_fwd_mesh_name + ".node";
    std::string new_face_name = new_fwd_mesh_name + ".face";
    std::string new_ele_name = new_fwd_mesh_name + ".ele";
    std::string new_neigh_name = new_fwd_mesh_name + ".neigh";
    std::string new_vtk_name = new_fwd_mesh_name + ".vtk";
    // change the file name
    std::rename(old_node_name.c_str(), new_node_name.c_str());
    std::rename(old_face_name.c_str(), new_face_name.c_str());
    std::rename(old_ele_name.c_str(), new_ele_name.c_str());
    std::rename(old_neigh_name.c_str(), new_neigh_name.c_str());
    std::rename(old_vtk_name.c_str(), new_vtk_name.c_str());

    // !!! update struct variables fwd_para_3d and init_prior_mod_3d
    // !! update fwd_para_3d
    // At this version, all frequencies share one forward mesh
    assert( (this->fwd_para_3d->starting_model).size() == 1);
    (*this->fwd_para_3d).starting_model[0] = this->init_fwd_mod_name;
    (*this->fwd_para_3d).starting_counter = this->k + 1;
    // !! update init_prior_mod_3d
    (this->init_prior_mod_3d.counter)++;
    if(have_marker_6666666_or_not){
      // ensure the inversion mesh are refined
      assert( (tm + 2) >= this->init_prior_mod_3d.n_regions );
      this->init_prior_mod_3d.n_regions = tm + 2;
    }
    else{
      // ensure the inversion mesh are refined
      assert( (tm + 1) >= this->init_prior_mod_3d.n_regions );
      this->init_prior_mod_3d.n_regions = tm + 1;
    }
    // ! update this->init_prior_mod_3d.fwd_tet_id_in_mj
    assert( (this->init_prior_mod_3d.m_ini).size() == tm );
    (this->init_prior_mod_3d.fwd_tet_id_in_mj).resize( tm );
    // read forward mesh
    std::stringstream fwd_counter;
    fwd_counter << (*this->fwd_para_3d).starting_counter;
    std::string fwd_s = fwd_counter.str();
    std::string fwd_mesh_name = (*this->fwd_para_3d).starting_model[0]
                                + std::string(".") + fwd_s;
    Mesh3D fwd_mesh(fwd_mesh_name);
    fwd_mesh.prepare_mesh_by_reading_files();
    assert(fwd_mesh.n_elems() > 0);
    for(unsigned int i = 0; i < fwd_mesh.n_elems(); i++){
      Tet* fwd_tet= static_cast<Tet*>(&(fwd_mesh.get_elem(i)));
      assert(i == fwd_tet->get_id());
      const unsigned int fwd_tet_marker = fwd_tet->get_tet_marker();
      if( (fwd_tet_marker != 9999999) && (fwd_tet_marker != 6666666) ){
        if( fwd_tet_marker > tm ){
          std::cout << "fwd_tet_marker: " << fwd_tet_marker << "tm: " << tm << std::endl;
          std::cout << "Input error, please check forward mesh (the nested forward mesh\n"
                    << "should have the same regional markers as inversion mesh)!\n";
          screen_output << "fwd_tet_marker: " << fwd_tet_marker << "tm: " << tm << std::endl;
          screen_output << "Input error, please check forward mesh (the nested forward mesh\n"
                        << "should have the same regional markers as inversion mesh)!\n";
          std::abort();
        }
        (this->init_prior_mod_3d.fwd_tet_id_in_mj[fwd_tet_marker - 1]).push_back( i );
      }
    }
    // guarantee the markers of the forward and those of inversion mesh are the same
    for(unsigned int j = 0; j < tm; j++)
      if( (this->init_prior_mod_3d.fwd_tet_id_in_mj[j]).size() == 0 ){
        std::cout << "j: " << j << std::endl;
        std::cout << "Input error, please check forward mesh (the nested forward mesh\n"
                  << "should have the same regional markers as inversion mesh)!\n";
        screen_output << "j: " << j << std::endl;
        screen_output << "Input error, please check forward mesh (the nested forward mesh\n"
                      << "should have the same regional markers as inversion mesh)!\n";
        std::abort();
      }
    // ! update this->init_prior_mod_3d.m_ref
    // there is no need to update region_table_ref variable
    /*choosing the inverted model as the reference model is better than 
    choosing the original homogeneous half-space reference model,
    which is because the homogeneous model is more far way from 
    the true model than the inverted model*/
    this->init_prior_mod_3d.m_ref = this->init_prior_mod_3d.m_ini;
  }

  return;
}

unsigned int Adaptive::Get_n_dofs(Mesh3D& mesh3d)
{
  this->n_dofs = 0;
  std::set<std::vector<unsigned int> > edges_tet; // E 
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
  }
  
  this->n_dofs = edges_tet.size();
  // return the numbers of global dofs
  return this->n_dofs;
}
