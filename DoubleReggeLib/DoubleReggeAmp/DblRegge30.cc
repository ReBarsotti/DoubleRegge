#include <cassert>
#include <complex>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "DblRegge30.h"


DblRegge30::DblRegge30( const vector< string >& args ) :
  UserAmplitude< DblRegge30 >( args )
{       
  assert( args.size() == 18 );
	
  m_b_par = AmpParameter( args[0] );
  m_c0 = AmpParameter( args[1]);
  m_c1 = AmpParameter( args[2]);
  m_c2 = AmpParameter( args[3]);
  m_n0 = AmpParameter( args[4]);
  m_n1 = AmpParameter( args[5]);
  m_n2 = AmpParameter( args[6]);
  m_d10 = AmpParameter(args[7]);
  m_d11 = AmpParameter( args[8]); 
  m_d12 = AmpParameter( args[9]);
  m_d20 = AmpParameter( args[10] );
  m_d21 = AmpParameter( args[11] );
  m_d22 = AmpParameter( args[12] );	
  m_aPrime = AmpParameter( args[13] );
  m_a0 = AmpParameter( args[14] );
  m_s0 = AmpParameter( args[15] );
  m_fastParticle = atoi( args[16].c_str() );
  m_charge = atoi( args[17].c_str() );

  registerParameter( m_b_par );
  registerParameter( m_c0 );
  registerParameter( m_c1 );
  registerParameter( m_c2 );
  registerParameter( m_n0 );
  registerParameter( m_n1 );
  registerParameter( m_n2 );
  registerParameter( m_d10 );
  registerParameter( m_d11 );
  registerParameter( m_d12 );
  registerParameter( m_d20);
  registerParameter( m_d21 );
  registerParameter( m_d22 );
  registerParameter( m_aPrime );
  registerParameter( m_a0 );
  registerParameter( m_s0 );
}

complex< GDouble >
DblRegge30::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  // use a short variable name
  GDouble* u = userVars;
  complex<double> i(0,1);

  double alpha1 = m_aPrime*u[k_t1] + m_a0;
  double alpha2 = m_aPrime*u[k_t]  + m_a0;

  alpha2 = ( alpha2 < 0 ? 0 : alpha2 );

  double v1 = m_d10 + m_d11*u[k_s1] + m_d12*u[k_s2];
  double v2 = m_d20 + m_d21*u[k_s1] + m_d22*u[k_s2];

  complex<double> xi1   = ( complex<double>( -1, 0 ) + exp( -i*M_PI*alpha1 ) ) / 2.;
  complex<double> xi2   = ( complex<double>( -1, 0 ) + exp( -i*M_PI*alpha2 ) ) / 2.;
  complex<double> xi12  = ( complex<double>(  1, 0 ) + exp( -i*M_PI*(alpha1-alpha2) ) ) / 2.;
  complex<double> xi21  = ( complex<double>(  1, 0 ) + exp( -i*M_PI*(alpha2-alpha1) ) ) / 2.;

  complex<double> ss2 = pow( u[k_s]/m_s0, alpha1 ) * pow( u[k_s2]/m_s0, alpha2-alpha1 );
  complex<double> ss1 = pow( u[k_s]/m_s0, alpha2 ) * pow( u[k_s1]/m_s0, alpha1-alpha2 );
  
  // this is the t dependence of the cross section -- we'll multiply the
  // amplitude by the sqrt of this since it the amplitude is squared in
  // in the intensity -- note the arguments of the exponentials are always
  // negative and the expreission os contructed to be strictly positive
  double tDep = sqrt( m_n0*m_n0*exp( m_c0*m_c0*u[k_t] ) +
		      m_n1*m_n1*exp( m_c1*m_c1*u[k_t] )*fabs(u[k_t]) +  // fabs to keep this term positive
		      m_n2*m_n2*exp( m_c2*m_c2*u[k_t] )*u[k_t]*u[k_t] );

  return ( ss2*xi1*xi21*v1 + ss1*xi2*xi12*v2 )*tDep;
}

void DblRegge30::calcUserVars( GDouble** pKin, GDouble* userVars ) const{
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
  TLorentzVector eta    ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
  TLorentzVector pi     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );

  userVars[k_s]  = ( proton + eta + pi ).M2();
  userVars[k_s1] = ( eta + pi ).M2();
  userVars[k_t]  = ( eta + pi - beam ).M2();

  // m_fastParticle is set to 2 for eta and 3 for pi
  userVars[k_s2] = ( m_fastParticle == 2 ? ( proton + pi ).M2() : ( proton + eta ).M2() );
  userVars[k_t1] = ( m_fastParticle == 2 ? ( eta - beam ).M2()  : ( pi - beam ).M2()    );
}
