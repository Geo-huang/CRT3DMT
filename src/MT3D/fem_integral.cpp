/*****************************************************************************/
/*                                                                           */
/*  Copyright 2020                                                           */
/*  Zhengyong Ren and Huang Chen                                             */
/*  renzhengyong@csu.edu.cn, chenhuang@cqu.edu.cn                            */
/*                                                                           */
/*****************************************************************************/



// exp(-iwt) 
#include <cassert>
#include "fem_integral.h"
#include "bc_info.h"
#include "tet.h"
#include "gauss_12D.h"
#include "gauss_tri.h"
#include "gauss_tet.h"

FEM_Integral::FEM_Integral(Tet& tet_):
  tet  (&tet_),
  g_tri (new GaussTri(5)), // 2-order enough
  g_tet (new GaussTet(3))  // 2-order enough
{

  int e2n_[6][2]={{0,1},{0,2},{0,3},{1,2},{3,1},{2,3}};
  int f2n_[4][3]={{1,3,2},{0,2,3},{0,3,1},{0,1,2}};
  int f2e_[4][3]={{5,3,4},{5,2,1},{4,0,2},{3,1,0}};
  for(int i=0; i<6; i++)
    for(int j=0; j<2; j++) {
      this->e2n[i][j]=e2n_[i][j];
    }
  for(int i=0; i<4; i++)
    for(int j=0; j<3; j++) {
      this->f2n[i][j]=f2n_[i][j];
      this->f2e[i][j]=f2e_[i][j];
    }

  this->init();

}


FEM_Integral::FEM_Integral():
  tet  (NULL),
  g_tri (NULL), 
  g_tet (NULL)  
{

}


FEM_Integral::~FEM_Integral()
{
  if(g_tri!=NULL) delete g_tri;
  if(g_tet!=NULL) delete g_tet;
}



void FEM_Integral::init()
{
  // volume 
  Point& v0= *tet->get_node(0);
  Point& v1= *tet->get_node(1);
  Point& v2= *tet->get_node(2);
  Point& v3= *tet->get_node(3);
  this->V= tet->get_size();


  // length of edges
  Point v[4]={v0, v1, v2, v3};
  for(int i=0; i<6; i++)   
     this->l[i]= (v[e2n[i][0]]-v[e2n[i][1]]).size();

  // a,b,c,d
  a[0]= det(v1(0),v2(0),v3(0),v1(1),v2(1),v3(1),v1(2),v2(2),v3(2));
  a[1]= det(v0(0),v2(0),v3(0),v0(1),v2(1),v3(1),v0(2),v2(2),v3(2))*-1.0;
  a[2]= det(v0(0),v1(0),v3(0),v0(1),v1(1),v3(1),v0(2),v1(2),v3(2));
  a[3]= det(v0(0),v1(0),v2(0),v0(1),v1(1),v2(1),v0(2),v1(2),v2(2))*-1.0;
  b[0]= det(1., 1., 1., v1(1),v2(1),v3(1),v1(2),v2(2),v3(2))*-1.0;
  b[1]= det(1., 1., 1., v0(1),v2(1),v3(1),v0(2),v2(2),v3(2));
  b[2]= det(1., 1., 1., v0(1),v1(1),v3(1),v0(2),v1(2),v3(2))*-1.0;
  b[3]= det(1., 1., 1., v0(1),v1(1),v2(1),v0(2),v1(2),v2(2));
  c[0]= det(v1(0),v2(0),v3(0), 1., 1., 1.,v1(2),v2(2),v3(2))*-1.0;
  c[1]= det(v0(0),v2(0),v3(0), 1., 1., 1.,v0(2),v2(2),v3(2));
  c[2]= det(v0(0),v1(0),v3(0), 1., 1., 1.,v0(2),v1(2),v3(2))*-1.0;
  c[3]= det(v0(0),v1(0),v2(0), 1., 1., 1.,v0(2),v1(2),v2(2));
  d[0]= det(v1(0),v2(0),v3(0),v1(1),v2(1),v3(1), 1., 1., 1.)*-1.0;
  d[1]= det(v0(0),v2(0),v3(0),v0(1),v2(1),v3(1), 1., 1., 1.);
  d[2]= det(v0(0),v1(0),v3(0),v0(1),v1(1),v3(1), 1., 1., 1.)*-1.0;
  d[3]= det(v0(0),v1(0),v2(0),v0(1),v1(1),v2(1), 1., 1., 1.);

  // modify directions of vector shape functions 
  for(int e=0; e<6; e++) {
    Node& i_1=  *tet->get_node(e2n[e][0]);
    Node& i_2=  *tet->get_node(e2n[e][1]);
    assert(i_1.get_id()!=i_2.get_id());
    // we assume that id(i_1)<id(i_2)
    if(i_1.get_id()>i_2.get_id()) {
      directions[e]= -1.0;
    }
    else {
      directions[e]=  1.0;
    }
  }


  return;

}

void FEM_Integral::N(Point p, Point n[6])
{
  // Here, we use shape function Ni=Li=Ni(x,y,z), 
  // see Jin,2002, P168, eqn.5.11, FEM for electromagentics
  // geometrical mapping functions and its xyz gradients
  assert(this->V>0.);
  const double factor= 1.0/(6.0*this->V);
  double L[4];
  Point grad_L[4];
  for(int i=0; i<4; i++) L[i]= (a[i]+b[i]*p(0)+c[i]*p(1)+d[i]*p(2))*factor;
  for(int i=0; i<4; i++) {
    grad_L[i](0)= b[i]*factor;
    grad_L[i](1)= c[i]*factor;
    grad_L[i](2)= d[i]*factor;
  }
/*
  // checking------------
  for(int i=0; i<4; i++) {
    Point& v= *tet->get_node(i);

    double t=(a[i]+b[i]*v(0)+c[i]*v(1)+d[i]*v(2))*factor;
    std::cout<<t<<"\t";
  }
*/  
  // vector shape function on each edge
  for(int i=0; i<6; i++) {
    n[i]= (grad_L[e2n[i][1]]*L[e2n[i][0]]-
	   grad_L[e2n[i][0]]*L[e2n[i][1]])*l[i];
  }
  
  // checking directions
  for(int i=0; i<6; i++) {
    n[i] = n[i]* directions[i];
  }

  return;

}



void FEM_Integral::Curl_N(Point p, Point curl_n[6])
{
  // p should be within "tet"
  assert(V>0.);
  for(int i=0; i<6; i++) {
    double factor= 2.0*l[i]/std::pow(6.0*V,2);
    curl_n[i](0)= (c[e2n[i][0]]*d[e2n[i][1]]-d[e2n[i][0]]*c[e2n[i][1]])*factor;
    curl_n[i](1)= (d[e2n[i][0]]*b[e2n[i][1]]-b[e2n[i][0]]*d[e2n[i][1]])*factor;
    curl_n[i](2)= (b[e2n[i][0]]*c[e2n[i][1]]-c[e2n[i][0]]*b[e2n[i][1]])*factor;
  }

  // checking directions
  for(int i=0; i<6; i++) {
    curl_n[i] = curl_n[i]*directions[i];
  }

  return;
}


void FEM_Integral::Curl_N_dot_N(Matrix6D& E,  Matrix6D& F)
{
  E.setZero(); 
  F.setZero();

  // fullfill E(6x6)
  for(int i=0; i<6; i++)
   for(int j=0; j<6; j++) {
     E(i,j)= 4.0*l[i]*l[j]*V*(1.0/std::pow(6.*V,4))*(
	(c[e2n[i][0]]*d[e2n[i][1]]-d[e2n[i][0]]*c[e2n[i][1]])*
        (c[e2n[j][0]]*d[e2n[j][1]]-d[e2n[j][0]]*c[e2n[j][1]])
      + (d[e2n[i][0]]*b[e2n[i][1]]-b[e2n[i][0]]*d[e2n[i][1]])*
        (d[e2n[j][0]]*b[e2n[j][1]]-b[e2n[j][0]]*d[e2n[j][1]])
      +	(b[e2n[i][0]]*c[e2n[i][1]]-c[e2n[i][0]]*b[e2n[i][1]])*
        (b[e2n[j][0]]*c[e2n[j][1]]-c[e2n[j][0]]*b[e2n[j][1]]) );
   }
  long double f[4][4];
  for(int i=0; i<4; i++)
   for(int j=0; j<4; j++) 
     f[i][j]= b[i]*b[j]+c[i]*c[j]+d[i]*d[j];
  // fullfill F(6x6)
  // Jin's book is wrong for F(1,2),F(0,0),F(1,1),F(2,2),F(3,3),F(4,4),F(5,5) (Ren)!
  // to be checked for the present edition, those problems may have been corrected!
  F(0,0)= (l[0]*l[0]/(360.*V))*(f[1][1]-f[0][1] +f[0][0]);
  F(0,1)= (l[0]*l[1]/(720.*V))*(2.*f[1][2]-f[1][0]-f[0][2] +f[0][0]);
  F(0,2)= (l[0]*l[2]/(720.*V))*(2.*f[1][3]-f[1][0]-f[0][3] +f[0][0]);
  F(0,3)= (l[0]*l[3]/(720.*V))*(f[1][2]-f[1][1]-2.*f[0][2] +f[0][1]);
  F(0,4)= (l[0]*l[4]/(720.*V))*(f[1][1]-f[1][3]-f[0][1] +2.*f[0][3]);
  F(0,5)= (l[0]*l[5]/(720.*V))*(f[1][3]-f[1][2]-f[0][3] +f[0][2]);
  F(1,1)= (l[1]*l[1]/(360.*V))*(f[2][2]-f[0][2] +f[0][0]);
  F(1,2)= (l[1]*l[2]/(720.*V))*(2.*f[2][3]-f[2][0]-f[0][3] +f[0][0]);
  F(1,3)= (l[1]*l[3]/(720.*V))*(f[2][2]-f[1][2]-f[0][2] +2.*f[0][1]);
  F(1,4)= (l[1]*l[4]/(720.*V))*(f[1][2]-f[2][3]-f[0][1] +f[0][3]);
  F(1,5)= (l[1]*l[5]/(720.*V))*(f[0][2]-f[2][2]-2.*f[0][3] +f[2][3]);
  F(2,2)= (l[2]*l[2]/(360.*V))*(f[3][3]-f[0][3] +f[0][0]);
  F(2,3)= (l[2]*l[3]/(720.*V))*(f[2][3]-f[1][3]-f[0][2]+f[0][1]);
  F(2,4)= (l[2]*l[4]/(720.*V))*(f[1][3]-f[3][3]-2.*f[0][1]+f[0][3]);
  F(2,5)= (l[2]*l[5]/(720.*V))*(f[3][3]-f[2][3]-f[0][3]+2.*f[0][2]);
  F(3,3)= (l[3]*l[3]/(360.*V))*(f[2][2]-f[1][2] +f[1][1]);
  F(3,4)= (l[3]*l[4]/(720.*V))*(f[1][2]-2.*f[2][3] -f[1][1] +f[1][3]);
  F(3,5)= (l[3]*l[5]/(720.*V))*(f[2][3]-f[2][2] -2.*f[1][3] +f[1][2]);
  F(4,4)= (l[4]*l[4]/(360.*V))*(f[1][1]-f[1][3] +f[3][3]);
  F(4,5)= (l[4]*l[5]/(720.*V))*(f[1][3]-2.*f[1][2] -f[3][3]+f[2][3]);
  F(5,5)= (l[5]*l[5]/(360.*V))*(f[3][3]-f[2][3] +f[2][2]);
  for(int i=0; i<6; i++)
   for(int j=0; j<i; j++)
     F(i,j)=F(j,i);

  // direction correctness
  for(int i=0; i<6; i++) 
    for(int j=0; j<6; j++) {
      E(i,j)*= directions[i]*directions[j];
      F(i,j)*= directions[i]*directions[j];
    }
    
  return;
}


// Surface integral of Int_{Ni_dot_(n_x_A) A is on T1
// A = H_0*z0
void FEM_Integral::N_dot_n_x_A(unsigned int face, BC& bc, Grad& C)
{
  // bc is well tested!
  C.setZero();
  assert(tet->get_neighbor(face)==NULL);
  // get normal of "face"
  // make sure normal is correct (pointing out)
  std::auto_ptr<Element> tri= tet->build_side(face);
  Point normal= tri->get_normal();
  Point gpoint= tet->get_gpoint();
  Point gpoint_tri= tri->get_gpoint();
  if((gpoint_tri-gpoint)*normal<0.) {
    normal= normal*-1.0;
  }
  double size= tri->get_size();
  assert(size>0.);
  // nodes on "face"
  Point& v0= *tet->get_node(f2n[face][0]);
  Point& v1= *tet->get_node(f2n[face][1]);
  Point& v2= *tet->get_node(f2n[face][2]);

  GaussTri& g= *this->g_tri;
  const std::vector<Point>& q=g.get_points();
  const std::vector<Real>&  w=g.get_weights();
  for(unsigned int i=0; i<q.size(); i++) {
     Point r=v0*q[i](0)+ v1*q[i](1)+v2*q[i](2);
     Point N_tet[6]; // N's values on each edge of tet at "r"
     this->N(r, N_tet);
     Point n[3]; // N's values on face 
     for(int j=0; j<3; j++) n[j]= N_tet[f2e[face][j]];
     std::vector<Dcomplex> e, h; // E and H at gauss points
     bc.compute_E_H(r);
     bc.get_E_H(e,h);  assert(e.size()==h.size()&& e.size()==4);
     Dcomplex z0= bc.get_z0(); // z0=i*omega*u0
     ComplexPoint A0=ComplexPoint(h[0], h[1], h[2])*z0;  
     Real WJ=w[i]*size;  // Jacobi*Weights over "face"
   
     for(int j=0; j<3; j++) {
       // "+" will be changed to "-" in the function FEM::compute_B()
       C(j) += n[j]*(normal.cross(A0))*WJ;
     }
  }

  return;
}


void FEM_Integral::N_dot_N(Matrix6D& F)
{
  F.setZero();

  long double f[4][4];
  for(int i=0; i<4; i++)
   for(int j=0; j<4; j++) 
     f[i][j]= b[i]*b[j]+c[i]*c[j]+d[i]*d[j];
  // fullfill F(6x6)
  // Jin's book is wrong for F(1,2),F(0,0),F(1,1),F(2,2),F(3,3),F(4,4),F(5,5) (Ren)!
  // to be checked for the present edition, those problems may have been corrected!
  F(0,0)= (l[0]*l[0]/(360.*V))*(f[1][1]-f[0][1] +f[0][0]);
  F(0,1)= (l[0]*l[1]/(720.*V))*(2.*f[1][2]-f[1][0]-f[0][2] +f[0][0]);
  F(0,2)= (l[0]*l[2]/(720.*V))*(2.*f[1][3]-f[1][0]-f[0][3] +f[0][0]);
  F(0,3)= (l[0]*l[3]/(720.*V))*(f[1][2]-f[1][1]-2.*f[0][2] +f[0][1]);
  F(0,4)= (l[0]*l[4]/(720.*V))*(f[1][1]-f[1][3]-f[0][1] +2.*f[0][3]);
  F(0,5)= (l[0]*l[5]/(720.*V))*(f[1][3]-f[1][2]-f[0][3] +f[0][2]);
  F(1,1)= (l[1]*l[1]/(360.*V))*(f[2][2]-f[0][2] +f[0][0]);
  F(1,2)= (l[1]*l[2]/(720.*V))*(2.*f[2][3]-f[2][0]-f[0][3] +f[0][0]);
  F(1,3)= (l[1]*l[3]/(720.*V))*(f[2][2]-f[1][2]-f[0][2] +2.*f[0][1]);
  F(1,4)= (l[1]*l[4]/(720.*V))*(f[1][2]-f[2][3]-f[0][1] +f[0][3]);
  F(1,5)= (l[1]*l[5]/(720.*V))*(f[0][2]-f[2][2]-2.*f[0][3] +f[2][3]);
  F(2,2)= (l[2]*l[2]/(360.*V))*(f[3][3]-f[0][3] +f[0][0]);
  F(2,3)= (l[2]*l[3]/(720.*V))*(f[2][3]-f[1][3]-f[0][2]+f[0][1]);
  F(2,4)= (l[2]*l[4]/(720.*V))*(f[1][3]-f[3][3]-2.*f[0][1]+f[0][3]);
  F(2,5)= (l[2]*l[5]/(720.*V))*(f[3][3]-f[2][3]-f[0][3]+2.*f[0][2]);
  F(3,3)= (l[3]*l[3]/(360.*V))*(f[2][2]-f[1][2] +f[1][1]);
  F(3,4)= (l[3]*l[4]/(720.*V))*(f[1][2]-2.*f[2][3] -f[1][1] +f[1][3]);
  F(3,5)= (l[3]*l[5]/(720.*V))*(f[2][3]-f[2][2] -2.*f[1][3] +f[1][2]);
  F(4,4)= (l[4]*l[4]/(360.*V))*(f[1][1]-f[1][3] +f[3][3]);
  F(4,5)= (l[4]*l[5]/(720.*V))*(f[1][3]-2.*f[1][2] -f[3][3]+f[2][3]);
  F(5,5)= (l[5]*l[5]/(360.*V))*(f[3][3]-f[2][3] +f[2][2]);
  
  for(int i=0; i<6; i++)
   for(int j=0; j<i; j++)
     F(i,j)=F(j,i);

  // direction correctness
  for(int i=0; i<6; i++) 
    for(int j=0; j<6; j++)
      F(i,j)*= directions[i]*directions[j];
    
  return;
}
