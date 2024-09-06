/*****************************************************************************/
/*                                                                           */
/*  Copyright 2017                                                           */
/*  Zhengyong Ren                                                            */
/*  renzhengyong@csu.edu.cn                                                  */
/*                                                                           */
/*****************************************************************************/



#include "gauss_12D.h"
#include "gauss_legendre.h" //1D library for any order p

//----------------------------------------------------------------------
Gauss1D::Gauss1D(const unsigned int p)
{
  // init.
  this->init(p);
  // assert
  assert(this->_points.size()>0);
  assert(this->_weights.size()>0);
  assert(this->_points.size()==this->_weights.size());
}


 
void Gauss1D::init(const unsigned int p)
{
  // inites points and weights
  this->_points.clear(); this->_points.resize(p);
  this->_weights.clear();this->_weights.resize(p);

  // p-points accuracte enough for 2p-1 order
  switch(p)
    {
    case 1:
      {
	_points[0]  =0.;   _weights[0] =2.;
	
	return;
      }
    case 2:
      {
	_points[0]  =-0.57735026919; 	_weights[0] =1.0;
	_points[1]  = 0.57735026919; 	_weights[1] =1.0;
	
	return;
      }
    case 3:
      {
	_points[0]  =-0.774596669241;   _weights[0] =0.555555555556;
	_points[1]  = 0.;               _weights[1] =0.888888888889;
	_points[2]  = 0.774596669241; 	_weights[2] =0.555555555556;

	return;
      }
    case 4:
      {
	_points[0]  =-0.861136311594;   _weights[0] =0.347854845137;	
	_points[1]  =-0.339981043585;   _weights[1] =0.652145154863;
	_points[2]  = 0.339981043585;   _weights[2] =0.652145154863;
	_points[3]  = 0.861136311594;   _weights[3] =0.347854845137;	
	
	return;
      }
    case 5:
      {
	_points[0]  =-0.906179845939;   _weights[0] =0.236926885056;
	_points[1]  =-0.538469310106;   _weights[1] =0.478628670499;
	_points[2]  = 0.;               _weights[2] =0.568888888889;
	_points[3]  = 0.538469310106;   _weights[3] =0.478628670499;
	_points[4]  = 0.906179845939;   _weights[4] =0.236926885056;

	return;	
      }
    case 6:
      {
	_points[0]  =-0.932469514203;   _weights[0] =0.171324492379;
	_points[1]  =-0.661209386466;   _weights[1] =0.360761573048;
	_points[2]  =-0.238619186083;   _weights[2] =0.467913934573;
	_points[3]  = 0.238619186083;   _weights[3] =0.467913934573;
	_points[4]  = 0.661209386466;   _weights[4] =0.360761573048;
	_points[5]  = 0.932469514203;   _weights[5] =0.171324492379;
	
	return;
      }
    case 7:
      {
	_points[0]  =-0.949107912343;   _weights[0] =0.129484966169;
	_points[1]  =-0.741531185599;   _weights[1] =0.279705391489;
	_points[2]  =-0.405845151377;   _weights[2] =0.381830050505;
	_points[3]  = 0.;               _weights[3] =0.417959183673;
	_points[4]  = 0.405845151377;   _weights[4] =0.381830050505;
	_points[5]  = 0.741531185599;   _weights[5] =0.279705391489;
	_points[6]  = 0.949107912343;   _weights[6] =0.129484966169;
	
	return;
      }
    case 8:
      {
	_points[0]  =-0.960289856498;   _weights[0] =0.10122853629;
	_points[1]  =-0.796666477414;   _weights[1] =0.222381034453;
	_points[2]  =-0.525532409916;   _weights[2] =0.313706645878;
	_points[3]  =-0.183434642496;   _weights[3] =0.362683783378;
	_points[4]  = 0.183434642496;   _weights[4] =0.362683783378;
	_points[5]  = 0.525532409916;   _weights[5] =0.313706645878;
	_points[6]  = 0.796666477414;   _weights[6] =0.222381034453;
	_points[7]  = 0.960289856498;   _weights[7] =0.10122853629;
	
	return;
      }
    case 9:
      {
	_points[0]  =-0.968160239508;   _weights[0] =0.0812743883616;
	_points[1]  =-0.836031107327;   _weights[1] =0.180648160695;
	_points[2]  =-0.613371432701;   _weights[2] =0.260610696403;
	_points[3]  =-0.324253423404;   _weights[3] =0.31234707704;
	_points[4]  = 0.;               _weights[4] =0.330239355001;
	_points[5]  = 0.324253423404;   _weights[5] =0.31234707704;
	_points[6]  = 0.613371432701;   _weights[6] =0.260610696403;
	_points[7]  = 0.836031107327;   _weights[7] =0.180648160695;
	_points[8]  = 0.968160239508;   _weights[8] =0.0812743883616;

	return;
      }
    case 10:
      {
	_points[0]  =-0.973906528517;   _weights[0] =0.0666713443087;
	_points[1]  =-0.865063366689;   _weights[1] =0.149451349151;
	_points[2]  =-0.679409568299;   _weights[2] =0.219086362516;
	_points[3]  =-0.433395394129;   _weights[3] =0.26926671931;
	_points[4]  =-0.148874338982;   _weights[4] =0.295524224715;
	_points[5]  = 0.148874338982;   _weights[5] =0.295524224715;
	_points[6]  = 0.433395394129;   _weights[6] =0.26926671931;
	_points[7]  = 0.679409568299;   _weights[7] =0.219086362516;
	_points[8]  = 0.865063366689;   _weights[8] =0.149451349151;
	_points[9]  = 0.973906528517;   _weights[9] =0.0666713443087;

	return;
      }
    case 11:
      {
	_points[0]  =-0.978228658146;   _weights[0] =0.0556685671162;
	_points[1]  =-0.887062599768;   _weights[1] =0.125580369465;
	_points[2]  =-0.730152005574;   _weights[2] =0.186290210928;
	_points[3]  =-0.519096129207;   _weights[3] =0.233193764592;
	_points[4]  =-0.269543155952;   _weights[4] =0.26280454451;
	_points[5]  = 0.;               _weights[5] =0.272925086778;
	_points[6]  = 0.269543155952;   _weights[6] =0.26280454451;
	_points[7]  = 0.519096129207;   _weights[7] =0.233193764592;
	_points[8]  = 0.730152005574;   _weights[8] =0.186290210928;
	_points[9]  = 0.887062599768;   _weights[9] =0.125580369465;
	_points[10] = 0.978228658146;   _weights[10]=0.0556685671162;

	return;
      }
    case 12:
      {
	_points[0]  =-0.981560634247;   _weights[0] =0.0471753363866;
	_points[1]  =-0.90411725637;    _weights[1] =0.106939325995;
	_points[2]  =-0.769902674194;   _weights[2] =0.160078328543;
	_points[3]  =-0.587317954287;   _weights[3] =0.203167426723;
	_points[4]  =-0.367831498998;   _weights[4] =0.233492536538;
	_points[5]  =-0.125233408511;   _weights[5] =0.249147045813;
	_points[6]  = 0.125233408511;   _weights[6] =0.249147045813;
	_points[7]  = 0.367831498998;   _weights[7] =0.233492536538;
	_points[8]  = 0.587317954287;   _weights[8] =0.203167426723;
	_points[9]  = 0.769902674194;   _weights[9] =0.160078328543;
	_points[10] = 0.90411725637;    _weights[10]=0.106939325995;
	_points[11] = 0.981560634247;   _weights[11]=0.0471753363866;

	return;
      }
    case 13:
      {
	_points[0]  =-0.984183054719;   _weights[0] =0.0404840047656;
	_points[1]  =-0.917598399223;   _weights[1] =0.0921214998379;
	_points[2]  =-0.801578090733;   _weights[2] =0.13887351022;
	_points[3]  =-0.64234933944;    _weights[3] =0.178145980762;
	_points[4]  =-0.448492751036;   _weights[4] =0.207816047537;
	_points[5]  =-0.230458315955;   _weights[5] =0.226283180263;
	_points[6]  = 0.0;              _weights[6] =0.232551553231;
	_points[7]  = 0.984183054719;   _weights[7] =0.0404840047656;
	_points[8]  = 0.917598399223;   _weights[8] =0.0921214998379;
	_points[9]  = 0.801578090733;   _weights[9] =0.13887351022;
	_points[10]  = 0.64234933944;   _weights[10] =0.178145980762;
	_points[11]  = 0.448492751036;  _weights[11] =0.207816047537;
	_points[12]  = 0.230458315955;   _weights[12] =0.226283180263;	
	
	return;
      }
    case 14:
      {
	_points[0]  =-0.986283808697;   _weights[0] =0.0351194603328;
	_points[1]  =-0.928434883664;   _weights[1] =0.0801580871593;
	_points[2]  =-0.82720131507;    _weights[2] =0.121518570688;
	_points[3]  =-0.687292904812;   _weights[3] =0.157203167158;
	_points[4]  =-0.515248636358;   _weights[4] =0.185538397478;
	_points[5]  =-0.319112368928;   _weights[5] =0.205198463721;
	_points[6]  =-0.108054948707;   _weights[6] =0.215263853463;
	_points[7]  = 0.986283808697;   _weights[7] =0.0351194603328;
	_points[8]  = 0.928434883664;   _weights[8] =0.0801580871593;
	_points[9]  = 0.82720131507;    _weights[9] =0.121518570688;
	_points[10]  = 0.687292904812;   _weights[10] =0.157203167158;
	_points[11]  = 0.515248636358;   _weights[11] =0.185538397478;
	_points[12]  = 0.319112368928;   _weights[12] =0.205198463721;
	_points[13]  = 0.108054948707;   _weights[13] =0.215263853463;

	return;
      }
    case 15:
      {
	_points[0]  =-0.98799251802;    _weights[0] =0.030753241995;
	_points[1]  =-0.937273392401;   _weights[1] =0.0703660474883;
	_points[2]  =-0.84820658341;    _weights[2] =0.107159220467;
	_points[3]  =-0.72441773136 ;   _weights[3] =0.139570677926;
	_points[4]  =-0.570972172609;   _weights[4] =0.166269205817;
	_points[5]  =-0.394151347078;   _weights[5] =0.186161000016;
	_points[6]  =-0.201194093997;   _weights[6] =0.198431485327;
	_points[7]  = 0.0;              _weights[7] =0.202578241926;
	_points[8]  = 0.98799251802;    _weights[8] =0.030753241995;
	_points[9]  = 0.937273392401;   _weights[9] =0.0703660474883;
	_points[10]  = 0.84820658341;    _weights[10] =0.107159220467;
	_points[11]  = 0.72441773136 ;   _weights[11] =0.139570677926;
	_points[12]  = 0.570972172609;   _weights[12] =0.166269205817;
	_points[13]  = 0.394151347078;   _weights[13] =0.186161000016;
	_points[14]  = 0.201194093997;   _weights[14] =0.198431485327;

	return;
      }
    case 16:
      {
	_points[0]  =-0.989400934992;   _weights[0] =0.027152459411;
	_points[1]  =-0.944575023073;   _weights[1] =0.0622535239372;
	_points[2]  =-0.865631202388;   _weights[2] =0.0951585116838;
	_points[3]  =-0.755404408355;   _weights[3] =0.124628971256;
	_points[4]  =-0.617876244403;   _weights[4] =0.149595988817;
	_points[5]  =-0.458016777657;   _weights[5] =0.169156519395;
	_points[6]  =-0.281603550779;   _weights[6] =0.182603415045;
	_points[7]  =-0.0950125098376;  _weights[7] =0.189450610455;
	_points[8]  = 0.989400934992;   _weights[8] =0.027152459411;
	_points[9]  = 0.944575023073;   _weights[9] =0.0622535239372;
	_points[10]  = 0.865631202388;   _weights[10] =0.0951585116838;
	_points[11]  = 0.755404408355;   _weights[11] =0.124628971256;
	_points[12]  = 0.617876244403;   _weights[12] =0.149595988817;
	_points[13]  = 0.458016777657;   _weights[13] =0.169156519395;
	_points[14]  = 0.281603550779;   _weights[14] =0.182603415045;
	_points[15]  = 0.0950125098376;  _weights[15] =0.189450610455;
	
	return;
      }
    case 17:
      {
	_points[0]  =-0.990575475314;   _weights[0] =0.0241483028649;
	_points[1]  =-0.950675521769;   _weights[1] =0.055459529376;
	_points[2]  =-0.880239153727;   _weights[2] =0.0850361483129;
	_points[3]  =-0.781514003897;   _weights[3] =0.111883847193;
	_points[4]  =-0.657671159217;   _weights[4] =0.13513636847;
	_points[5]  =-0.512690537086;   _weights[5] =0.154045761078;
	_points[6]  =-0.351231763454;   _weights[6] =0.168004102158;
	_points[7]  =-0.178484181496;   _weights[7] =0.176562705368;
	_points[8]  = 0.0;              _weights[8] =0.179446470358;
	_points[9]  =  0.990575475314;   _weights[9] =0.0241483028649;
	_points[10]  = 0.950675521769;   _weights[10] =0.055459529376;
	_points[11]  = 0.880239153727;   _weights[11] =0.0850361483129;
	_points[12]  = 0.781514003897;   _weights[12] =0.111883847193;
	_points[13]  = 0.657671159217;   _weights[13] =0.13513636847;
	_points[14]  = 0.512690537086;   _weights[14] =0.154045761078;
	_points[15]  = 0.351231763454;   _weights[15] =0.168004102158;
	_points[16]  = 0.178484181496;   _weights[16] =0.176562705368;
	
	return;

      }
    case 18:
      {
	_points[0] = -0.9915651684209310;    
	_points[1] =-0.9558239495713977;    
	_points[2] =-0.8926024664975558;    
	_points[3] =-0.8037049589725230;    
	_points[4] =-0.6916870430603533;    
	_points[5] =-0.5597708310739477;    
	_points[6] =-0.4117511614628424;    
	_points[7] =-0.2518862256915056;    
	_points[8] = -0.8477501304173471E-01;
	_points[9] =  0.8477501304173529E-01;
	_points[10] = 0.2518862256915056;    
	_points[11] =  0.4117511614628425;    
	_points[12] =  0.5597708310739474;    
	_points[13] =  0.6916870430603530;    
	_points[14] =  0.8037049589725234;    
	_points[15] =  0.8926024664975557;    
	_points[16] = 0.9558239495713978;    
	_points[17] =  0.9915651684209308;

	_weights[0] =0.2161601352648338E-01;
	_weights[1] =0.4971454889496975E-01;
	_weights[2] =0.7642573025488915E-01;
	_weights[3] =0.1009420441062867;    
	_weights[4] =0.1225552067114792;    
	_weights[5] =0.1406429146706502;    
	_weights[6] =0.1546846751262653;    
	_weights[7] =0.1642764837458330;    
	_weights[8] =0.1691423829631437;    
	_weights[9] =0.1691423829631434 ;   
	_weights[10] =0.1642764837458328;    
	_weights[11] =0.1546846751262653;    
	_weights[12] =0.1406429146706505 ;   
	_weights[13] =0.1225552067114788 ;   
	_weights[14] =0.1009420441062872 ;  
	_weights[15] =0.7642573025488915E-01;
	_weights[16] =0.4971454889496916E-01;
	_weights[17] =0.2161601352648266E-01;
	
	return;
      }
    case 19:
      {
	_points[0] =  -0.9924068438435842;    
	_points[1] =  -0.9602081521348304;    
	_points[2] =  -0.9031559036148182;   
	_points[3] =  -0.8227146565371429;   
	_points[4] =  -0.7209661773352291;   
	_points[5] =  -0.6005453046616805; 
	_points[6] =  -0.4645707413759609;  
	_points[7] =  -0.3165640999636295;    
	_points[8] =  -0.1603586456402254;   
	_points[9] =   0.1667949964511336E-15;
	_points[10] =   0.1603586456402253;    
	_points[11] =   0.3165640999636297;   
	_points[12] =   0.4645707413759609;   
	_points[13] =   0.6005453046616809;   
	_points[14] =   0.7209661773352291;   
	_points[15] =   0.8227146565371429;   
	_points[16] =   0.9031559036148178;   
	_points[17] =   0.9602081521348298;  
	_points[18] =   0.9924068438435842; 
	
	_weights[0] =  0.1946178822972707E-01;
	_weights[1] = 0.4481422676569900E-01;
	_weights[2] = 0.6904454273764107E-01;
	_weights[3] = 0.9149002162244987E-01;
	_weights[4] = 0.1115666455473343 ;   
	_weights[5] =  0.1287539625393365 ;   
	_weights[6] =  0.1426067021736069  ;  
	_weights[7] =  0.1527660420658594   ; 
	_weights[8] = 0.1589688433939545;    
	_weights[9] = 0.1610544498487837;    
	_weights[10] =  0.1589688433939540;    
	_weights[11] =  0.1527660420658593;    
	_weights[12] =  0.1426067021736076;    
	_weights[13] = 0.1287539625393368;    
	_weights[14] =  0.1115666455473341;    
	_weights[15] =  0.9149002162244975E-01;
	_weights[16] = 0.6904454273764120E-01;
	_weights[17] = 0.4481422676569921E-01;
	_weights[18] = 0.1946178822972656E-01;
	
	return;
      }
    case 20:
      {
	_points[0] =  -0.9931285991850950;    
	_points[1] =  -0.9639719272779136;    
	_points[2] = -0.9122344282513257;    
	_points[3] = -0.8391169718222189;    
	_points[4] = -0.7463319064601508;    
	_points[5] = -0.6360536807265149;    
	_points[6] = -0.5108670019508270;    
	_points[7] = -0.3737060887154196;    
	_points[8] = -0.2277858511416451;    
	_points[9] = -0.7652652113349742E-01;
	_points[10] = 0.7652652113349727E-01;
	_points[11] = 0.2277858511416456;   
	_points[12] = 0.3737060887154196;    
	_points[13] = 0.5108670019508271;    
	_points[14] = 0.6360536807265151;    
	_points[15] = 0.7463319064601509;    
	_points[16] = 0.8391169718222189;    
	_points[17] = 0.9122344282513254;    
	_points[18] = 0.9639719272779137;    
	_points[19] = 0.9931285991850947;    

	_weights[0] =   0.1761400713915209E-01;
	_weights[1] =   0.4060142980038740E-01;
	_weights[2] =  0.6267204833410854E-01;
	_weights[3] =  0.8327674157670478E-01;
	_weights[4] =  0.1019301198172402;    
	_weights[5] =  0.1181945319615188;    
	_weights[6] =  0.1316886384491762;    
	_weights[7] =  0.1420961093183826;    
	_weights[8] =  0.1491729864726034;    
	_weights[9] =   0.1527533871307258;    
	_weights[10] =   0.1527533871307260;    
	_weights[11] =  0.1491729864726043;    
	_weights[12] =   0.1420961093183828;    
	_weights[13] =  0.1316886384491771;    
	_weights[14] =  0.1181945319615193 ;   
	_weights[15] = 0.1019301198172400 ;   
	_weights[16] =  0.8327674157670395E-01;
	_weights[17] =  0.6267204833410935E-01;
	_weights[18] =  0.4060142980038609E-01;
	_weights[19] =   0.1761400713915237E-01;
	
	return;
      }
      
    default:
      {
        // for case p>20, call third library p[21, infinite]
        assert(p>20);
        unsigned int n=p;
        double eps=1e-16; unsigned int m= (n+1)>>1;
        double x[m], w[m];
        gauss_legendre_tbl(n, x, w,eps); // call third-library
        this->_points.clear();  this->_points.resize(n);
        this->_weights.clear(); this->_weights.resize(n);
        if((n%2)==0) { //even n
           assert(m==n/2);
           for(unsigned int i=0; i<m; i++) {
             _points[i+n/2]= x[i];  
             _weights[i+n/2]= w[i]; 
           }
           for(unsigned int i=0; i<n-m; i++) {
             _points[n-m-i-1]= x[i]*-1.0;  
             _weights[n-m-i-1]= w[i]; 
           }
        }else { // old
           assert(m==n/2+1);
           for(unsigned int i=0; i<m; i++) {
             _points[i+n/2]= x[i];  
             _weights[i+n/2]= w[i]; 
           }
           for(unsigned int i=0; i<n-m; i++) {
             _points[n-m-i-1]= x[i+1]*-1.0;  
             _weights[n-m-i-1]= w[i+1]; 
           }
        }
             
        //std::cout<<(*this);
      }
    } // switch end

  return;  
}




 
std::ostream& operator <<(std::ostream& os, const Gauss1D& gauss)
{
  os<<"1D gauss rule for [-1,1] with"
    <<gauss._points.size()<<"points\n";

  assert(gauss._points.size()==gauss._weights.size());
  for(unsigned int i=0; i<gauss._points.size(); i++)
    {
      os<<std::setprecision(8);
      os<<gauss._weights[i]<<"\t"
	<<gauss._points[i]<<"\n";
    }
  os<<std::endl;
  return os;
}



//----------------------------------------------------------------------
Gauss2D::Gauss2D(const unsigned int p)
{
  // init.
  this->init(p);
  // assert
  assert(this->_points.size()>0);
  assert(this->_weights.size()>0);
  assert(this->_points.size()==this->_weights.size());
}


void Gauss2D::init(const unsigned int p)
{
  Gauss1D g1d(p);
  const std::vector<Real>&  q=g1d.get_points();
  const std::vector<Real>&  w=g1d.get_weights();
  unsigned int N=q.size();
  _points.resize(N*N);
  _weights.resize(N*N);
  for(unsigned int i=0; i<N; i++)
   for(unsigned int j=0; j<N; j++) {
     _weights[i*N+j]= w[i]*w[j];
     _points [i*N+j](0)= q[i]; //epsilon
     _points [i*N+j](1)= q[j];
     _points [i*N+j](2)= 0.;
   }

  return;
}

 
std::ostream& operator <<(std::ostream& os, const Gauss2D& gauss)
{
  os<<"2D gauss rule [-1,1]x[-1,1] with "
    <<gauss._points.size()<<" points\n";

  assert(gauss._points.size()==gauss._weights.size());
  for(unsigned int i=0; i<gauss._points.size(); i++)
    {
      os<<std::setprecision(8);
      os<<gauss._weights[i]<<"\t"
	<<gauss._points[i]<<"\n";
    }
  os<<std::endl;
  return os;
}

