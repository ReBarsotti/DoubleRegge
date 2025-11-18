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
#include "IUAmpTools/report.h"

#include "DblRegge30.h"

const char* DblRegge30::kModule = "DblRegge30";

DblRegge30::DblRegge30( const vector< string >& args ) :
	UserAmplitude< DblRegge30 >( args )
{       
	assert( args.size() == 14 );
	
	m_c0 = AmpParameter(args[0]);
	m_c1 = AmpParameter(args[1]);
	m_c2 = AmpParameter(args[2]);
	m_c3 = AmpParameter(args[3]);
	m_d0 = AmpParameter(args[4]);
	m_d1 = AmpParameter(args[5]);
	m_d2 = AmpParameter(args[6]);
	m_d3 = AmpParameter(args[7]);
	m_a0 = AmpParameter( args[8] );
	m_a0_2 = AmpParameter( args[9] );
	m_aPrime = AmpParameter( args[10] );
	m_aPrime_2 = AmpParameter( args[11] );
	m_s0 = AmpParameter( args[12] );
	m_fastParticle = atoi( args[13].c_str() );


	registerParameter(m_c0);
	registerParameter(m_c1);
	registerParameter(m_c2);
	registerParameter(m_c3);
	registerParameter(m_d0);
	registerParameter(m_d1);
	registerParameter(m_d2);
	registerParameter(m_d3);
	registerParameter( m_aPrime );
	registerParameter( m_aPrime_2 );
	registerParameter( m_a0 );
	registerParameter( m_a0_2 );
	registerParameter( m_s0 );
}

complex< GDouble >
DblRegge30::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

	// use a short variable name
	GDouble* u = userVars;
	complex<double> i(0,1);


	double alpha1 = m_aPrime*m_aPrime*u[k_t1] + m_a0;
	double alpha2 = m_aPrime_2*m_aPrime_2*u[k_t]  + m_a0_2;

	double eta = m_s0*u[k_s] / (u[k_s2]*u[k_s1]);


	complex< double> v1 = (m_c0 + m_c1*(1./eta) + m_c2*(1./eta)*(1./eta) + m_c3*(1./eta)*(1./eta)*(1./eta));
	complex< double> v2 = (m_d0 + m_d1*(1./eta) + m_d2*(1./eta)*(1./eta) + m_d3*(1./eta)*(1./eta)*(1./eta));

	complex<double> xi1   = ( complex<double>( -1, 0 ) + exp( -i*M_PI*alpha1 ) ) / 2.;
	complex<double> xi2   = ( complex<double>( -1, 0 ) + exp( -i*M_PI*alpha2 ) ) / 2.;
	complex<double> xi12  = ( complex<double>(  1, 0 ) + exp( -i*M_PI*(alpha1-alpha2) ) ) / 2.;
	complex<double> xi21  = ( complex<double>(  1, 0 ) + exp( -i*M_PI*(alpha2-alpha1) ) ) / 2.;


	complex<double> ss2 = pow( u[k_s]/m_s0, alpha1 ) * pow( u[k_s2]/m_s0, alpha2-alpha1 );
	complex<double> ss1 = pow( u[k_s]/m_s0, alpha2 ) * pow( u[k_s1]/m_s0, alpha1-alpha2 );

	double barrier = sqrt(-u[k_t1]); //*(1+m_c0*cos(u[k_phiGJ]));

	// this is the t dependence of the cross section -- we'll multiply the
	// amplitude by the sqrt of this since it the amplitude is squared in
	// in the intensity -- note the arguments of the exponentials are always
	// negative and the expreission os contructed to be strictly positive

	return barrier*((ss2*xi1*xi21*v1 + ss1*xi2*xi12*v2 ));

}

void
DblRegge30::updatePar( const AmpParameter& par ) {

	if( par == m_a0 || par == m_aPrime ){

		report( DEBUG, kModule ) << "in updatePar(): a0, aPrime = "
			<< m_a0 << ", " << m_aPrime << endl;
	}
}

void
DblRegge30::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector target (          0,          0,          0,      0.938 );
	TLorentzVector proton ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector eta    ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
	TLorentzVector pi     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );

	userVars[k_s]  = ( proton + eta + pi ).M2();
	userVars[k_s1] = ( eta + pi ).M2();
	userVars[k_t]  = ( eta + pi - beam ).M2();
	userVars[k_pi] = pi.M();

	// m_fastParticle is set to 2 for eta and 3 for pi
	userVars[k_s2] = ( m_fastParticle == 2 ? ( proton + pi ).M2() : ( proton + eta ).M2() );
	userVars[k_t1] = ( m_fastParticle == 2 ? ( eta - beam ).M2()  : ( pi - beam ).M2()    );

	// get the scattering angle
	TLorentzRotation cmBoost = -( target + beam ).BoostVector();
	TLorentzVector protonCM = cmBoost * proton;
	userVars[k_thetaCM] = acos( fabs( protonCM.CosTheta() ) );

	TLorentzVector beamCM = cmBoost * beam;
	TLorentzVector etaCM = cmBoost * eta;


	userVars[k_tmin] = etaCM.M2() - 2*etaCM.E()*beamCM.E() + 2*etaCM.Vect().Mag()*beamCM.Vect().Mag();

	// now boost to the GJ frame -- everything about is frame invaraint anyway...
	TLorentzVector etaPi = eta + pi;
	eta.Boost( -1 * etaPi.BoostVector() );
	beam.Boost( -1 * etaPi.BoostVector() );
	proton.Boost( -1 * etaPi.BoostVector() );

	TVector3 z = beam.Vect().Unit();
	TVector3 y = (beam.Vect().Unit().Cross(proton.Vect().Unit())).Unit();
	TVector3 x = y.Cross(z);

	TVector3 angles( (eta.Vect()).Dot(x),
			(eta.Vect()).Dot(y),
			(eta.Vect()).Dot(z) );

	userVars[k_phiGJ] = angles.Phi();
	userVars[k_thetaGJ] = angles.Theta();

	double alpha1 =  m_aPrime*m_aPrime*userVars[k_t1] + m_a0;
	double alpha2 =  m_aPrime_2*m_aPrime_2*userVars[k_t]  + m_a0_2;

}
