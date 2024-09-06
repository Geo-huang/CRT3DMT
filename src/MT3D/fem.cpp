/*****************************************************************************/
/*                                                                           */
/*  Copyright 2020                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/
#include <numeric> 
#include "fem.h"
#include "fem_integral.h" // subroutins for all integrals 
#include "gauss_12D.h"
#include "gauss_tet.h"
#include "gauss_tri.h"
#include "timer.h"
 // used for generating screen output log
 extern std::ofstream screen_output;

// user_defined routines to return conductivity and 
// dielectric permittivity and magnetic permeability
// and wavenumber at a given point
extern void electric_parameters(Point p, Real& cond, Real& epsilon, 
				Real& mu, Dcomplex& k,const Real f,
				const unsigned int tet_marker);

extern void electric_parameters(Point p, Real& cond, Real& epsilon, 
		                     Real& mu, Dcomplex& k, const Real f, 
                         const unsigned int tet_marker,
                         EM::VectorXD m_k, double abn[3]);

// ---------------------------------------------------------------------------
// used for producing synthetic data 
FEM::FEM(Mesh3D&  mesh3d, 
	       const Real frequency, 
       	 BC&  bc_te, // boundary condtions
       	 BC&  bc_tm, // boundary condtions
         std::string inv_algorithm,
         const int sites_marker):
  _mesh3d          (mesh3d),  
  _dofs            (mesh3d),
  _f               (frequency),
  _inv_algorithm   (inv_algorithm),
  _sites_marker    (sites_marker)
{ 
   this->_bc[0]= &bc_te;
   this->_bc[1]= &bc_tm;
   this->_X.resize(2);
}

// used for computing the forward responses and sensitivities  
// when doing inversion
FEM::FEM(Mesh3D&  mesh3d, 
	       const Real frequency, 
       	 BC&  bc_te, // boundary condtions
       	 BC&  bc_tm, // boundary condtions
         std::string inv_algorithm,
         const int sites_marker,
         std::vector<Point>& _u,
         std::vector<Point>& _v,
         EM::VectorXD m_k,
         const double abn[3]):
  _mesh3d          (mesh3d),  
  _dofs            (mesh3d),
  _f               (frequency),
  _inv_algorithm   (inv_algorithm),
  _sites_marker    (sites_marker),
  u                (_u),
  v                (_v),
  _m_k             (m_k)
{ 
  // initialize array _abn;
  this->_abn[0] = abn[0];
  this->_abn[1] = abn[1];
  this->_abn[2] = abn[2];

  this->_bc[0]= &bc_te;
  this->_bc[1]= &bc_tm;
  this->_X.resize(2);

  // ! initialize site_ID_to_local_index map
  //   which will be used to calculate the local interpolator
  unsigned int local_index = 0;
  const unsigned int n_nodes = (this->_mesh3d).n_nodes();
  for(unsigned int i = 0; i < n_nodes; i++){
    Node *a = &(this->_mesh3d).get_node(i);
    if(a->get_marker() == this->_sites_marker){
      // to be checked!
      site_ID_to_local_index[i] = local_index;
      local_index++;
    }
  }
}

void FEM::solve()
{
  Timer t;  
  t.start();
  this->compute_A();  
  this->compute_B();  
  t.stop();  
  //std::cout<<"Compute A, B0, B1 USING: \t" << t.getElapsedTimeInSec()<<" (s)\n";
  t.start(); 
  // Performs a symbolic decomposition on the sparcity of matrix
  this->solver.analyzePattern(A); 
  // Performs a numeric decomposition of matrix
  this->solver.factorize(A);
  t.stop(); 
  //std::cout<<"Time of factorizing A is:\t" << t.getElapsedTime()<<"\n";
  t.start(); 
  // the solution x of $ A x = B0 $ using the current decomposition of A
  E0 = this->solver.solve(B0);
  E1 = this->solver.solve(B1);
  //std::cout << E0 << std::endl;
  t.stop(); 
  //std::cout<<"Time of solving AE0=B0, AE1=B1 is:\t" <<t.getElapsedTime()<<"\n";
  // calculate E and H at on the air-earth surface
  this->compute_electromagnetic_fields(); 
  //std::cout<<"EM computed on the air-earth surface\n";

  return;
}

void FEM::compute_A()
{
  typedef Eigen::Triplet<Dcomplex> T;
  std::vector<T> coefficients;    // list of non-zeros coefficients
  unsigned int n_E = this->_dofs.get_n_dofs();
  coefficients.reserve(n_E * 100);// space reservation
  this->A.resize(n_E, n_E);       // space reserved
  this->A.setZero();              // A=0
  // loop all tets
  assert(_mesh3d.get_elem().size()>0);
  //#pragma omp parallel for
  for (int e=0; e<_mesh3d.n_elems(); e++) {
    Tet* tet= static_cast<Tet*>(&(_mesh3d.get_elem(e)));
    std::vector<unsigned int> N_dofs;
    this->_dofs.global_dofs_E(tet, N_dofs); // for <6x6> submatrix for E
    assert(N_dofs.size()==6);
    Matrix6D E;   // E = integral of curl_Ni*curl_Nj over tet, i,j = 1,...,6
    Matrix6D F;   // F = integral of Ni*Nj ovet tet, i,j = 1,...,6
    E.setZero();  F.setZero(); // VERY IMPORTANT E=0 F=0
    FEM_Integral fem_integral(*tet);
    fem_integral.Curl_N_dot_N(E,F);
    Point gpoint = tet->get_gpoint(); // center point of tet
    Real cond_tet=0., epsilon_tet=0., mu_tet=0.; //relative mu and epsilon
    Dcomplex k_tet=0.;  // k in tet
    const unsigned int tet_marker= tet->get_tet_marker();
    // !!!get parameter in each tet
    if(this->_inv_algorithm == "Fwd_only")
      electric_parameters(gpoint, cond_tet, epsilon_tet, mu_tet, k_tet, 
			                    this->_f, tet_marker);
    else{
      electric_parameters(gpoint, cond_tet, epsilon_tet, mu_tet, k_tet, 
			                    this->_f, tet_marker, this->_m_k, this->_abn);
    }
    assert(cond_tet>0.);
    assert(mu_tet>0);
    Real mu_r_tet= mu_tet/EM::MU0; 
    assert(mu_r_tet>1e-3); //relative mu
    assert(std::abs(k_tet)>0.);
    Dcomplex k1= k_tet*k_tet;
    Matrix6C EF;  
    EF.setZero(); // VERY IMPORTANT EF=0
    for(unsigned int i=0; i<6; i++)
      for(unsigned int j=0; j<6; j++) {
       //#pragma omp critical 
        {
	        EF(i,j) = (E(i,j)- F(i,j)*k1)*(1./(mu_r_tet));	
          coefficients.push_back(T(N_dofs[i], N_dofs[j], EF(i,j)));
        }
      }
  } // o(m) loop
  
  this->A.setFromTriplets(coefficients.begin(), coefficients.end());
  coefficients.clear();

  return;
}

void FEM::compute_B()
{
  unsigned int n_E = this->_dofs.get_n_dofs();
  B0.resize(n_E); B0.setZero();
  B1.resize(n_E); B1.setZero();
  BCInfo& bc_info = _mesh3d.get_boundary_info();
  // 2 is the default marker of the boundary of the computational domain
  const short int T1_marker= 2;  // \partial \Omega
  // loop all tets
  assert(_mesh3d.get_elem().size()>0);
  // std::cout<<_mesh3d.n_elems()<<"\n";
  for (unsigned int e=0; e<_mesh3d.n_elems(); e++) {
    Tet* tet= static_cast<Tet*>(&(_mesh3d.get_elem(e)));       
    if(tet->on_boundary()) {
      std::vector<unsigned short int> side;
      tet->get_boundary_side(side);      
      for(unsigned int i=0; i<side.size(); i++) {
	      if(bc_info.get_side_marker(tet, side[i])==T1_marker) {
	        FEM_Integral fem_integral(*tet);
	        std::vector<unsigned int> E_dofs;
	        this->_dofs.global_dofs_E(tet, side[i], E_dofs);    
	        assert(E_dofs.size()==3); 
	        Grad N_dot_n_x_A[2];    
          N_dot_n_x_A[0].setZero();
          N_dot_n_x_A[1].setZero();
          // A = H_0*z0   
          fem_integral.N_dot_n_x_A(side[i], *this->_bc[0], N_dot_n_x_A[0]);
          fem_integral.N_dot_n_x_A(side[i], *this->_bc[1], N_dot_n_x_A[1]);
          // A = -H_0*z0  
          N_dot_n_x_A[0]*=-1.0;
          N_dot_n_x_A[1]*=-1.0;
          for(int m=0; m<3; m++) {
	          this->B0.coeffRef(E_dofs[m]) += N_dot_n_x_A[0](m);
	          this->B1.coeffRef(E_dofs[m]) += N_dot_n_x_A[1](m);
          }
	      }
      }
    }
  } // o(m) loop
  
  return;
}

void FEM::compute_electromagnetic_fields()
{
/*****************************************************************************/
// This is for edge-based finite element methods.
// Calculating E/H at all nodes on air-earth-surface
// E0, E1, and Mesh3D are ready
// Using E(p) = sum_{i}{Ne}E_{i}N_{i}(p) to calculate E at point p
// point p is on air-earth-interface
// H = \nabla \times E * (1/z0) which is constant in each tet
// E and H are computed by simply averaging their values from all tets
// these tets should share the point p
// exp(-iwt)
/*****************************************************************************/
  // get \hat y and \hat z in air layer
  const Dcomplex y0= this->_bc[0]->get_y0();  // sigma -i omega epsilon
  const Dcomplex z0= this->_bc[0]->get_z0();  // +i omega mu
  this->_X.clear();
  this->_X.resize(2);
  for(unsigned int i=0; i<(E0).size(); i++) { // _X[0]=E0 _X[1]=E1
    this->_X[0].push_back((E0)(i));
    this->_X[1].push_back((E1)(i));
  }
  // build air-earth-interface (T with marker == 1 in our code)
  BCInfo& bc_info = _mesh3d.get_boundary_info();
  // 1 represents the facet marker of the air-Earth interface
  const short int T_marker= 1; 
  const int air_tet_marker = 9999999;
  Surface T;
  T.build_surface(bc_info, T_marker, air_tet_marker);
  //          E00,E01,H00,H01        &         E10,E11,H10,H11,
  // ------TE source polarization----  -----TM source polarization----
  // map from the node id of site on the air-Earth interface to the number of 
  // the tets which share a common vertex at that node 
  std::map<unsigned int, unsigned int> Count;
  // map from the node id of site on the air-Earth interface to the vector of
  // the id of tets which have a common vertex at site
  std::map<unsigned int, std::vector<unsigned int> > node_tet;
  // loop tri of "T", which are triangles on the air-Earth surface
  for(unsigned int i=0; i < T.get_surface().size(); i++) {
    Tri* tri_i = T.get_surface()[i]; // get a tri on T
    // normal vector of the triangle surface tri_i
    Point n    = T.get_n_surface(tri_i->get_id());
    std::vector<unsigned int> DOFs_E;
    std::pair<const Element*, unsigned short int> 
    // !!! get the triangular face on the air-Earth interface
    tet_face_i = T.get_tet_face(tri_i);  // where am I?
    // !!!get the tetrahedral element tet_i which has a face of tet_face_i
    Element* tet_i= const_cast<Element*>(tet_face_i.first);
    // make sure no air-tets are involved!
    assert(tet_i->get_tet_marker()!=9999999); // 9999999--air layer
    unsigned short int face_i= tet_face_i.second;
    // changing edge local number in tet element to the global number  
    // (global dofs) among the edges in the computational domain)
    this->_dofs.global_dofs_E(static_cast<Tet*>(tet_i), DOFs_E);
    assert(DOFs_E.size() == 6);
    FEM_Integral fem_integral(*(static_cast<Tet*>(tet_i)));
    const unsigned int tet_marker= tet_i->get_tet_marker();
    unsigned int tet_id= tet_i->get_id();  
    // !!!loop nodes of tri on "tet"
    for(unsigned int k=0; k<3; k++) {
      Node& v= *(tri_i->get_node(k));
      unsigned int v_marker = v.get_marker();
      if(v_marker == this->_sites_marker)
      {
        unsigned int ID = v.get_id();
        // "real" point
        Point n_E[6]; Point curl_n[6];
        fem_integral.N(v, n_E); // Shape function n_E at "v"
        fem_integral.Curl_N(v, curl_n);
        node_tet[v.get_id()].push_back(tet_i->get_id());
        Real cond_v=0., epsilon_v=0., mu_v=0.; Dcomplex k_v=0.;  
        if(this->_inv_algorithm == "Fwd_only")
          electric_parameters(v, cond_v, epsilon_v, mu_v, k_v, this->_f,
                              tet_marker);
        else
          electric_parameters(v, cond_v, epsilon_v, mu_v, k_v, this->_f,
                              tet_marker, this->_m_k, this->_abn);
        const double omega = 2*EM::PI*(this->_f);
        const Dcomplex y1 = Dcomplex(cond_v,omega*epsilon_v*-1.0);
        const Dcomplex z1 = Dcomplex(0., omega*mu_v);
        // D and Q are used for claculating e and h fields at air side
        const Dcomplex D  = y1/y0;
        const Dcomplex Q  = z1/z0;
        assert(std::abs(Q)>0.9); 
        // e, h fields at earth side
        ComplexPoint E1[2], H1[2];

        unsigned int local_index = site_ID_to_local_index[ID];
        // define and initialize three orthogonal bases of the local coordinate
        // system at (local_index)-th site, which will be  used for computing the
        // local interpolator 
        // unit vector basis at x direction
        Point e_x = Point(1.0, 0.0, 0.0);
        // unit vector basis at y direction
        Point e_y = Point(0.0, 1.0, 0.0);
        // unit vector basis at z direction
        Point e_z = Point(0.0, 0.0, 1.0);
        if(this->_inv_algorithm != "Fwd_only"){
          // rotate the x,y,z coordiante system to u,v,n coordiante system 
          e_x = this->u[local_index];
          e_y = this->v[local_index];
          e_z = (e_x.cross(e_y)).unit();
        }

        //E(p) = sum_{i}{Ne}E_{i}N_{i}(p)
        for(unsigned int ns=0; ns<6; ns++) {
          E1[0] = E1[0] + n_E[ns]*_X[0][DOFs_E[ns]]; 
          E1[1] = E1[1] + n_E[ns]*_X[1][DOFs_E[ns]];
          H1[0] = H1[0] + (curl_n[ns]*_X[0][DOFs_E[ns]])*(1./z1);
	        H1[1] = H1[1] + (curl_n[ns]*_X[1][DOFs_E[ns]])*(1./z1);
        }
        // e, h fields at air side, these fields are not used 
        // computing via continuity condition
        ComplexPoint E0[2], H0[2];
        E0[0]= (n.cross(E1[0])).cross(n)  + n*((n*E1[0])*D);
        E0[1]= (n.cross(E1[1])).cross(n)  + n*((n*E1[1])*D);
        H0[0] = (n.cross(H1[0])).cross(n) + n*((n*H1[0])*Q);
        H0[1] = (n.cross(H1[1])).cross(n) + n*((n*H1[1])*Q);     
        // iterator of the map from node id to ComplexPoint
        typedef std::map<unsigned int, 
	      std::vector<ComplexPoint> >::iterator  IT;
        IT it = _EH.find(ID);
        if(it!= _EH.end()) { // for existed node id!
          std::vector<ComplexPoint>& value= (*it).second;
          assert(value.size()==8);
          // please note: the E0,1 and H0,1 is a local variable, not this->E0, E1
          // the different source polarizations are differentiated by index, 0 and 1
          // which represent source polarizations TE and TM, respectively 
          value[0]= value[0]+ E0[0];  // air    E  TE
          value[1]= value[1]+ E1[0];  // earth  E  TE
          value[2]= value[2]+ H0[0];  // air    H  TE
          value[3]= value[3]+ H1[0];  // earth  H  TE
          value[4]= value[4]+ E0[1];  // air    E  TM
          value[5]= value[5]+ E1[1];  // earth  E  TM
          value[6]= value[6]+ H0[1];  // air    H  TM
          value[7]= value[7]+ H1[1];  // earth  H  TM
          std::map<unsigned int, unsigned int>::iterator cit=
	        Count.find(ID);
          assert(cit!=Count.end());
          (*cit).second = (*cit).second + 1;

          // !!!added for calculation of local interpolator 
          if(this->_inv_algorithm != "Fwd_only"){
            // used for assigning of the local interpolation operator  
            // iterator of the map from node id to local interpolator 
            // (type:EigenSparseVectorXD)
            typedef std::map<unsigned int, EM::EigenSparseVectorXD>::iterator IT_e;
            IT_e it_x_e = gx_e.find(ID);
            IT_e it_y_e = gy_e.find(ID);
            // iterator of the map from node id to local interpolator 
            // (type:EigenSparseVectorXC)
            typedef std::map<unsigned int, EM::EigenSparseVectorXC>::iterator IT_h;
            IT_h it_x_h = gx_h.find(ID);
            IT_h it_y_h = gy_h.find(ID);
            IT_h it_z_h = gz_h.find(ID);
            assert(it_x_e != gx_e.end());
            EM::EigenSparseVectorXD& x_e = (*it_x_e).second;
            EM::EigenSparseVectorXD& y_e = (*it_y_e).second;       
            EM::EigenSparseVectorXC& x_h = (*it_x_h).second;
            EM::EigenSparseVectorXC& y_h = (*it_y_h).second;
            EM::EigenSparseVectorXC& z_h = (*it_z_h).second;
            for(unsigned int ns=0; ns<6; ns++){
              x_e.coeffRef( DOFs_E[ns] ) += n_E[ns]     * e_x;
              y_e.coeffRef( DOFs_E[ns] ) += n_E[ns]     * e_y;
              x_h.coeffRef( DOFs_E[ns] ) += (curl_n[ns] * e_x) * (1./z1);
              y_h.coeffRef( DOFs_E[ns] ) += (curl_n[ns] * e_y) * (1./z1);
              z_h.coeffRef( DOFs_E[ns] ) += (curl_n[ns] * e_z) * (1./z1);
            }
          }

        }else { // for new node id
          std::vector<ComplexPoint> value(8);
          // Ei[j], Hi[i], i=0 denotes the value at air side,
          // i=1 ... at earth side, j=0 denotes TE, j=1 denotes TM
          value[0]= E0[0]; 
          value[1]= E1[0]; 
          value[2]= H0[0]; 
          value[3]= H1[0];  
          value[4]= E0[1];
          value[5]= E1[1];
          value[6]= H0[1];
          value[7]= H1[1];
          _EH[ID] = value;
          std::map<unsigned int, unsigned int>::iterator cit=
  	      Count.find(ID);
          assert(cit==Count.end());
          Count[ID]=1;

          // !!!added for local interpolator calculation
          if(this->_inv_algorithm != "Fwd_only"){
            EM::EigenSparseVectorXD x_e;
            EM::EigenSparseVectorXD y_e;
            EM::EigenSparseVectorXC x_h;
            EM::EigenSparseVectorXC y_h;
            EM::EigenSparseVectorXC z_h;
            // size defining and initializing of the vector before using
            // is mandatory!
            // for right tetrahedron: the solid angle at a vertex subtended
            // by  a face is arccos(23/27), approx. 0.55129 steradians
            unsigned int estimation_of_non_zero_entries = 72;
            x_e.resize( this->_dofs.get_n_dofs() );
            x_e.setZero(); 
            x_e.reserve(estimation_of_non_zero_entries);
            y_e.resize( this->_dofs.get_n_dofs() );
            y_e.setZero();
            y_e.reserve(estimation_of_non_zero_entries);
            x_h.resize( this->_dofs.get_n_dofs() );
            x_h.setZero();
            x_h.reserve(estimation_of_non_zero_entries);
            y_h.resize( this->_dofs.get_n_dofs() );
            y_h.setZero();
            y_h.reserve(estimation_of_non_zero_entries);
            z_h.resize( this->_dofs.get_n_dofs() );
            z_h.setZero();
            z_h.reserve(estimation_of_non_zero_entries);
            for(unsigned int ns=0; ns<6; ns++) {
              x_e.insert( DOFs_E[ns] ) = n_E[ns] * e_x;
              y_e.insert( DOFs_E[ns] ) = n_E[ns] * e_y;
              x_h.insert( DOFs_E[ns] ) = (curl_n[ns] * e_x) * (1./z1);
              y_h.insert( DOFs_E[ns] ) = (curl_n[ns] * e_y) * (1./z1);
              z_h.insert( DOFs_E[ns] ) = (curl_n[ns] * e_z) * (1./z1);
            }
            gx_e[ID] = x_e;
            gy_e[ID] = y_e;
            gx_h[ID] = x_h;
            gy_h[ID] = y_h; 
            gz_h[ID] = z_h;
          } // if method != Fwd_only

        } // judging whether node is new 
      } // if v_marker == marker
    } // loop all 3 verteces
  } // loop all tri on "T"
  // electromagnetic fields averaging at sides on air-Earth interface...
  typedef std::map<unsigned int, 
    std::vector<ComplexPoint> >::iterator IT;
  for(IT it=_EH.begin(); it!=_EH.end(); it++) {
    // get the number of the tets which share a common vertex at node (*it).first 
    std::map<unsigned int, unsigned int>::iterator cit=
      Count.find((*it).first);
    assert(cit!=Count.end());
    const unsigned int n= (*cit).second;
    std::vector<ComplexPoint>& value= (*it).second;
    assert(value.size()==8);
    for(unsigned int i=0; i<value.size(); i++)
      value[i]= value[i]*(1./(n*1.0));
  }

  // !!!added for local interpolator calculation
  if(this->_inv_algorithm != "Fwd_only"){
    // local interpolation operator averaging at sides on air-Earth interface
    typedef std::map<unsigned int, EM::EigenSparseVectorXD>::iterator IT_e;
    // averaging gx_e
    for(IT_e it_x_e = gx_e.begin(); it_x_e != gx_e.end(); it_x_e++){
      // get the number of the tets which share a common vertex at node (*it).first 
      std::map<unsigned int, unsigned int>::iterator cit=
      Count.find((*it_x_e).first);
      assert(cit!=Count.end());
      const unsigned int n = (*cit).second;
      EM::EigenSparseVectorXD& x_e = (*it_x_e).second;
      x_e *= 1.0 / n;
    }
    // averaging gy_e
    for(IT_e it_y_e = gy_e.begin(); it_y_e != gy_e.end(); it_y_e++){
      // get the number of the tets which share a common vertex at node (*it).first 
      std::map<unsigned int, unsigned int>::iterator cit=
      Count.find((*it_y_e).first);
      assert(cit!=Count.end());
      const unsigned int n = (*cit).second;
      EM::EigenSparseVectorXD& y_e = (*it_y_e).second;
      y_e *= 1.0 / n;
    }
    typedef std::map<unsigned int, EM::EigenSparseVectorXC>::iterator IT_h;
    // averaging gx_h
    for(IT_h it_x_h = gx_h.begin(); it_x_h != gx_h.end(); it_x_h++){
      // get the number of the tets which share a common vertex at node (*it).first 
      std::map<unsigned int, unsigned int>::iterator cit=
      Count.find((*it_x_h).first);
      assert(cit!=Count.end());
      const unsigned int n = (*cit).second;
      EM::EigenSparseVectorXC& x_h = (*it_x_h).second;
      x_h *= 1.0 / n;
    }
    // averaging gy_h
    for(IT_h it_y_h = gy_h.begin(); it_y_h != gy_h.end(); it_y_h++){
      // get the number of the tets which share a common vertex at node (*it).first 
      std::map<unsigned int, unsigned int>::iterator cit=
        Count.find((*it_y_h).first);
      assert(cit!=Count.end());
      const unsigned int n = (*cit).second;
      EM::EigenSparseVectorXC& y_h = (*it_y_h).second;
      y_h *= 1.0 / n;
    }
    // averaging gz_h
    for(IT_h it_z_h = gz_h.begin(); it_z_h != gz_h.end(); it_z_h++){
      // get the number of the tets which share a common vertex at node (*it).first 
      std::map<unsigned int, unsigned int>::iterator cit=
        Count.find((*it_z_h).first);
      assert(cit!=Count.end());
      const unsigned int n = (*cit).second;
      EM::EigenSparseVectorXC& z_h = (*it_z_h).second;
      z_h *= 1.0 / n;
    }
  } // if method != Fwd_only

  return;
}

void FEM::comp_partial_A_partial_m(EigenSparseMatrixXC& dA_t_E0,
                                   EigenSparseMatrixXC& dA_t_E1,
      std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj)
{
  /////////////////////////////////////////////////////////////////////////////
  // This function is used for calculating of senstivity matrix and gradient //
  //     of the misfit term of the inversion objective function, where,      //
  //                 dA_t_E0 = \partial A / \partial m * E0,                 //
  //                 dA_t_E1 = \partial A / \partial m * E1                  //
  /////////////////////////////////////////////////////////////////////////////
  // guaranting the E0 and E1 has been solved by FEM when calling this function 
  assert(std::abs(this->E0.norm()) > 0.0);
  assert(std::abs(this->E1.norm()) > 0.0);
  // the number of the dofs corresponding to the forward mesh
  unsigned int n_E = this->_dofs.get_n_dofs();
  // the number of the tetrahedral elements in the inversion domain, i.e.,
  // the dimension of the model space
  unsigned int n_m = this->_m_k.size();
  // dimension definitions
  dA_t_E0.resize(n_E, n_m); dA_t_E1.resize(n_E, n_m);
  dA_t_E0.setZero();        dA_t_E1.setZero();

  typedef Eigen::Triplet< Dcomplex > T;
  // lists of non-zero coefficients
  std::vector< T > coefficients_0, coefficients_1; 
  // space reservation for coefficients vector
  coefficients_0.reserve(36 * n_m); coefficients_1.reserve(36 * n_m);

  // calculate all the non-zero coefficients in the j-th column of 
  // dA_t_E0 and dA_t_E1 matrices
  for(unsigned int j = 0; j < n_m; j++)
  {
    // get id vector of the forward tetrahedral elements nested in 
    // the j-th parameter element (belonging to inversion mesh)
    std::vector<unsigned int> fwd_tet_ids = fwd_tet_id_in_mj[j];

    // ! for calculationg of partial \sigma_j / partial m_j
    // linear conductivity of the j-th tet element in the inversion domain
    double cond = 0.;
    double nomin = this->_abn[0] + ( this->_abn[1] 
                   * std::exp( this->_abn[2] * this->_m_k(j) ) );
    double denomin = 1.0 + std::exp( this->_abn[2] * this->_m_k(j) );
    cond = nomin / denomin;
    if( std::isnan( cond ) )
      // we set sigma = b (upper bound), if it is equals to NAN
      cond =  this->_abn[1];
    // partial \sigma_j / partial m_j
    double p_sj_p_mj = 0.;
    nomin = this->_abn[2] * (this->_abn[1] - cond) * (cond - this->_abn[0]);
    if( nomin < 0 ){
      std::cout << "WARNING: CONDUCTIVITY VALE OUT OF BOUNDS EMERGENED!\n";
      screen_output << "WARNING: CONDUCTIVITY VALE OUT OF BOUNDS EMERGENED!\n";
    }
    denomin = this->_abn[1] - this->_abn[0];
    p_sj_p_mj = nomin / denomin;

    // used for storing all the maps of local edge number to global edge number
    // of all the tetrahedra in j-th parameter element
    std::vector < std::vector<unsigned int> > all_N_dofs;
    // used for storing all the F matrices beloning to all the tets
    // in j-th parameter element; F = integral of Ni*Nj over element tet
    std::vector< Matrix6D > all_F;
    for(unsigned int p = 0; p < fwd_tet_ids.size(); p++){
      Tet* tet= static_cast<Tet*>( &(_mesh3d.get_elem( fwd_tet_ids[p] )) );
      // !! initialize all_N_dofs vector
      std::vector<unsigned int> N_dofs;
      // maping of the local edge number to the global edge number
      this->_dofs.global_dofs_E(tet, N_dofs);
      assert(N_dofs.size()==6);
      all_N_dofs.push_back(N_dofs);

      // !! initialize all_F vector
      // F = integral of Ni*Nj over element tet
      Matrix6D F; 
      F.setZero(); 
      FEM_Integral fem_integral(*tet);
      // compute F 
      fem_integral.N_dot_N(F);
      all_F.push_back(F);
    }

    double omega = 2.0 * EM::PI * this->_f;
    Dcomplex zhat0 = EM::II * omega * EM::MU0;
    for(unsigned int e = 0; e < fwd_tet_ids.size(); e++){
      for(unsigned int i = 0; i < 6; i++){
        Dcomplex temp0, temp1 = 0.0; 
        for(unsigned int k = 0; k < 6; k++){
          temp0 -= all_F[e](i,k) * (zhat0 * p_sj_p_mj) 
                   * this->E0( all_N_dofs[e][k] );
          temp1 -= all_F[e](i,k) * (zhat0 * p_sj_p_mj) 
                   * this->E1( all_N_dofs[e][k] );
        }
        coefficients_0.push_back( T(all_N_dofs[e][i], j, temp0) );
        coefficients_1.push_back( T(all_N_dofs[e][i], j, temp1) );
      }
    }
  }

  dA_t_E0.setFromTriplets(coefficients_0.begin(), coefficients_0.end());
  dA_t_E1.setFromTriplets(coefficients_1.begin(), coefficients_1.end());
  coefficients_0.clear(); 
  coefficients_1.clear();

  return;
}

