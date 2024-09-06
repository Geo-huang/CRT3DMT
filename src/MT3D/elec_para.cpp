/******************************************************************************
 *    Function       : external function                                      *
 *    Author         : Huang Chen                                             *
 *    Copyright      : Huang Chen, 2020                                       *
 *    Email          : chenhuang@cqu.edu.cn                                   *
 *    Created time   : 2020.07.27                                             *
 *    Last revision  :                                                        *
 ******************************************************************************/
#include <iostream>
#include <vector>
#include "struct_defs.h"
#include "point.h"

extern Fwd_Para_3D fwd_para_3d;  
// external function, which will be frequently used by classes FEM in file fem.cpp
// the formal parameter Point p of the following two functions is useless 
void electric_parameters(Point p, Real& cond, Real& epsilon, 
		         Real& mu, Dcomplex& k, const Real f, 
                         const unsigned int tet_marker)
{
  /********assigning electric parameter when producing synthetic data*********/
  /****************************by using tet marker****************************/
  assert(!fwd_para_3d.region_table.empty());
  typedef std::map<int, std::vector<double> >::iterator IT;
  IT it= fwd_para_3d.region_table.find(tet_marker);
  // found the tet_marker in fwd_para_3d.region_table map variable
  assert(it!=fwd_para_3d.region_table.end());
  // assignment of electric property parameter vector
  std::vector<double>& parameter = (*it).second;
  assert(parameter.size()==3);  
  cond    = parameter[0];    
  epsilon = parameter[1]*EM::EPSILON0;
  mu      = parameter[2]*EM::MU0;
  double omega = 2.0*EM::PI*f;
  Dcomplex zhat = EM::II*omega*mu;
  Dcomplex yhat = cond - EM::II*omega*epsilon;
  k = std::sqrt(zhat*yhat);
  return;
}

void electric_parameters(Point p, Real& cond, Real& epsilon, 
		                     Real& mu, Dcomplex& k, const Real f, 
                         const unsigned int tet_marker, 
                         EM::VectorXD m_k, double abn[3])
{
  /************electric parameter assignement when doing inversion************/
  /*************************by using tet marker and tet ID *******************/
  // please note: the conductivity of the tet in the inversion domain is  
  // initialized by using tet marker, vector m_k and abn array
  // abn[0], abn[1], abn[2] represents the variables a(lower bound),  
  // b(upper bound) and n used for parameter transformation, respectively. 
  std::vector<double>& parameter = fwd_para_3d.region_table[tet_marker];
  assert(parameter.size()==3);  
  epsilon = parameter[1]*EM::EPSILON0;
  mu      = parameter[2]*EM::MU0;
  double omega = 2.0*EM::PI*f;
  Dcomplex zhat = EM::II*omega*mu;
  
  if( (tet_marker == 9999999) || (tet_marker == 6666666) ){ 
    // !!! for tet out of the inversion domain (having free parameters)
    cond = parameter[0];    
  }
  else{ 
    // !!! for tet in the inversion domain
    // lower bound a < upper bound b
    assert(abn[0] < abn[1]);
    // transformed from logarithmic space to linear domain
    double nomin = 0, denomin = 0.;
    nomin = abn[0] + ( abn[1] * std::exp( abn[2] * m_k(tet_marker -1) ) );
    denomin = 1.0 + std::exp( abn[2] * m_k(tet_marker -1) );
    cond = nomin / denomin;
    if( std::isnan( cond ) )
      // we set sigma = b (upper bound), if it is equals to NAN
      cond =  abn[1];
  }

  Dcomplex yhat = cond - EM::II*omega*epsilon;
  k = std::sqrt(zhat*yhat);

  return;
}

