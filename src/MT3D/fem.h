/*****************************************************************************/
/*                                                                           */
/*  Copyright 2020                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/

/* 
  a standard FEM approach to compute EM response based on E-field formulation
  The final linear system looks like
  AE0=B0  TE mode source
  AE1=B1  TM mode source
*/

#ifndef  _FEM_H
#define  _FEM_H

// Local includes
#include "mesh3d.h"
#include "bc.h"
#include "dofs.h"
#include <Eigen/PardisoSupport>
// ----------------------------------------------------------------------------
class FEM
{
  public: 
  // construction function, which will only be used for inv_algorithm = Fwd_only
  FEM( Mesh3D&  mesh3d,               // input mesh
       const Real frequency,          // Hz
       BC&  bc_te,                    // boundary condtions
       BC&  bc_tm,
       std::string inv_algorithm,     // including Fwd_only, L-BFGS
       const int sites_marker);                              

  // reloaded construction function, which will only be used for 
  // inv_algorithm = others (i.e., != Fwd_only)
  FEM( Mesh3D&  mesh3d,               // input mesh
       const Real frequency,          // Hz
       BC&  bc_te,                    // boundary condtions
       BC&  bc_tm,
       std::string inv_algorithm,     // including Fwd_only, L-BFGS
       const int sites_marker,        // marker of the sites 
       std::vector<Point>& _u,        // u-direction orthogonal bases at sites
       std::vector<Point>& _v,        // v-direction orthogonal bases at sites
       EM::VectorXD m_k,             // k-th inverted model, the present model  
       const double abn[3]);          // parameters controlling transformation
                                                   
  ~FEM() {}
  
  friend class PostProcess;
 
  // step 1: FEM on earth and air domains
  void compute_A();
  void compute_B(); 
  void solve();  // checked
  // compute and ouput E, H at each node 
  void compute_electromagnetic_fields(); // checked

  // added for senstivity and gradient of the misfit term calculation
  // dA_t_E0 = \partial A / \partial m * E0,
  // dA_t_E1 = \partial A / \partial m * E1,
  void comp_partial_A_partial_m(EigenSparseMatrixXC& dA_t_E0,
                                EigenSparseMatrixXC& dA_t_E1,
    std::vector< std::vector<unsigned int> > fwd_tet_id_in_mj);

  public:
  Mesh3D&                                _mesh3d;
  DOFs                                   _dofs;
  const Real                             _f;
  const int                              _sites_marker;
  // the following three variables are added for passing the parameter 
  // on the inversion mesh to the forward mesh when doing inversion
  std::string                            _inv_algorithm;
  EM::VectorXD                           _m_k;
  // abn[0], ~[1], ~[2] represents the variables a(lower bound), 
  // b(upper bound) and n used for parameter transformation, respectively. 
  double                                 _abn[3];

  BC*                                    _bc[2];
  EigenSparseMatrixXC                     A;
  EigenDenseVector                        B0, B1;
  EigenDenseVector                        E0, E1;
  // added for sensitivity calculation
  // call parsidso from eigen library to solve Ax=B
  Eigen::PardisoLU<EigenSparseMatrixXC>   solver;
  std::vector< std::vector<Dcomplex> >    _X;  // _X[0]=E0 _X[1]=E1
  //--------------------------------------------
  // store the E and H at each node on the ground
  // E00,E10,H00,H10 & E01,E11,H01,H11,
  // ------TE----        -----TM----
  std::map<unsigned int, 
	std::vector<ComplexPoint>  >           _EH; // E, H at nodes on ground
  //--------------------------------------------

  /////////////////////////////////////////////////////////////////////////////
  //             !!! used for calculating the local interpolators,
  //     which will be used to calculate the sensitivities and gradients
  //                         when doing inversion
  // map from the global id of the site (node) to its local id (the order
  // shown in the poly file)
  std::map<unsigned int, unsigned int>            site_ID_to_local_index;
  // u[i], v[i], n[i] = (u[i].cross(v[i])).unit() represent three orthonormal 
  // bases of the local coordinate system at the i-th site
  std::vector<Point>                              u;
  std::vector<Point>                              v;
  // map from the node id of the site on the air-Earth interface to its discrete 
  // local interpolation operator, which interpolate the electric
  // and magnetic fields from the forward modelling grid to the receiver
  // which will be initialized when inv_algorithm != Fwd_only
  // Note: we use gx,y,z_e,h to denote gu,v,n_e,h in this code!!!
  std::map<unsigned int, EM::EigenSparseVectorXD>  gx_e, gy_e;
  std::map<unsigned int, EM::EigenSparseVectorXC>  gx_h, gy_h, gz_h;
  /////////////////////////////////////////////////////////////////////////////
};

#endif // _FEM_H
