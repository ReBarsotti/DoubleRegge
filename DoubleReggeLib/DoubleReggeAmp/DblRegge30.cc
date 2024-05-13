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
	assert( args.size() == 22 );

	m_f0 = AmpParameter(args[0]);
	m_f1 = AmpParameter(args[1]);
	m_f2 = AmpParameter(args[2]);
	m_f3 = AmpParameter(args[3]);
	m_f4 = AmpParameter(args[4]);
	m_f5 = AmpParameter(args[5]);
	m_c0 = AmpParameter( args[6]);
	m_c1 = AmpParameter( args[7]);
	m_c2 = AmpParameter( args[8]);
	m_c3 = AmpParameter( args[9]);
	m_g0 = AmpParameter( args[10]);
	m_g1 = AmpParameter( args[11]);
	m_g2 = AmpParameter( args[12]);
	m_g3 = AmpParameter( args[13]);
	m_d1 = AmpParameter(args[14]);
	m_d2 = AmpParameter( args[15] );
	m_a0 = AmpParameter( args[16] );
	m_a0_2 = AmpParameter( args[17] );	
	m_aPrime = AmpParameter( args[18] );
	m_aPrime_2 = AmpParameter( args[19] );
	m_s0 = AmpParameter( args[20] );
	m_fastParticle = atoi( args[21].c_str() );


	registerParameter( m_f0 );
	registerParameter( m_f1 );
	registerParameter( m_f2 );
	registerParameter( m_f3 );
	registerParameter( m_f4 );
	registerParameter( m_f5 );
	registerParameter( m_c0 );
	registerParameter( m_c1 );
	registerParameter( m_c2 );
	registerParameter( m_c3 );

	registerParameter( m_g0 );
	registerParameter( m_g1 );
	registerParameter( m_g2 );
	registerParameter( m_g3 );
	registerParameter( m_d1 );
	registerParameter( m_d2);
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



	double alpha1 = m_aPrime*u[k_t1] + m_a0;
	double alpha2 = m_aPrime_2*u[k_t2] + m_a0_2;

	complex<double> xfunc = (m_c0 + m_c1*u[k_t1] + m_c2*u[k_t1]*u[k_t1])*exp(m_c3*m_c3*u[k_t1]);
	complex<double> yfunc = (m_g0 + m_g1*u[k_t2] + m_g2*u[k_t2]*u[k_t2])*exp(m_g3*m_g3*u[k_t2]);


	complex<double> ss2 = pow( u[k_s]/m_s0, alpha1 ) * pow( u[k_s2]/m_s0, alpha2-alpha1 );
	complex<double> ss1 = pow( u[k_s]/m_s0, alpha2 ) * pow( u[k_s1]/m_s0, alpha1-alpha2 );

	// careful here!!  if you decide to float a0 or a0Prime in the fit then
	// these precalculated values cannot be used
	//complex<double> gam1( u[k_gam1_re], u[k_gam1_im] );
	//complex<double> gam2( u[k_gam2_re], u[k_gam2_im] );
	//  these two lines would be appropriate if a0 and/or a0Prime and hence
	//  the values of alpha are floating in the fit...
	//	complex< double > gam1 = cgamma( -alpha1 );
	//	complex< double > gam2 = cgamma( -alpha2 );
	double spinFactor =  sin( u[k_thetaGJ] );
	double MASSPI2 = 0.01948;
	double phi_approx = (u[k_t1]+ u[k_t2] - MASSPI2 + (u[k_s2]*u[k_s1]/u[k_s])) / (2*sqrt(-u[k_t1])*sqrt(-u[k_t2])); 

	//	double angles = sqrt(u[k_s]*u[k_s]*(-u[k_t1])*sin(u[k_phiGJ])*sin(u[k_phiGJ]));
	double explicitAngle = 1 + m_f0*cos( u[k_phiGJ] ) + m_f1*cos( u[k_phiGJ])*cos( u[k_phiGJ]) + m_f2*cos(u[k_phiGJ])*cos(u[k_phiGJ])*cos( u[k_phiGJ]);
	//double explicitAngle = 1+ m_f0*phi_approx + m_f1*phi_approx*phi_approx + m_f2 * phi_approx*phi_approx*phi_approx;

	double explicitAngle2 =  1 +m_f3*cos( u[k_phiGJ] ) + m_f4*cos( u[k_phiGJ])*cos( u[k_phiGJ]) + m_f5*cos(u[k_phiGJ])*cos(u[k_phiGJ])*cos( u[k_phiGJ]);
	double v1 = m_d1;
	double v2 = m_d2;
	//double spinFactor =  u[k_s]*u[k_s1]*sqrt(abs((u[k_t1] - u[k_tmin] )));  
	//	double spinFactor= 1 + m_d22*cos(u[k_phiGJ]);
	return ( ss2*xfunc*v1*explicitAngle +yfunc*ss1*v2*explicitAngle2 )*spinFactor;
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
	userVars[k_t2]  = ( eta + pi - beam ).M2();

	//	userVars[k_tmin] = eta.M2() - 2*eta.E()*beam.E() + 2*eta.Vect().Mag()*beam.Vect().Mag();

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





	// the gamma functions are expensive -- as long as aPrime and a0
	// aren't varying, then we can cache them
	//	double alpha1 =  m_aPrime*userVars[k_t1] + m_a0;
	//	double alpha2 =  m_aPrime_2*userVars[k_t]  + m_a0_2;

	//	complex< double > gam1 = cgamma( -alpha1 );
	//	complex< double > gam2 = cgamma( -alpha2 );

	//	userVars[k_gam1_re] = real( gam1 );
	//	userVars[k_gam1_im] = imag( gam1 );

	//	userVars[k_gam2_re] = real( gam2 );
	//	userVars[k_gam2_im] = imag( gam2 );
}


complex<double>
DblRegge30::cgamma( complex<double> z ) const {

	complex<double> ui (0,1);
	complex<double> g, infini= 1e308+ 0.0*ui; // z0,z1
	double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
	double na=0.0,t,x1 = 1,y1=0.0,sr,si;
	int j,k;


	static double a[] = {
		8.333333333333333e-02,
		-2.777777777777778e-03,
		7.936507936507937e-04,
		-5.952380952380952e-04,
		8.417508417508418e-04,
		-1.917526917526918e-03,
		6.410256410256410e-03,
		-2.955065359477124e-02,
		1.796443723688307e-01,
		-1.39243221690590};

	x = real(z);
	y = imag(z);

	if (x > 171) return infini;
	if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
		return infini;
	else if (x < 0.0) {
		x1 = x;
		y1 = y;
		x = -x;
		y = -y;
	}
	x0 = x;
	if (x <= 7.0) {
		na = (int)(7.0-x);
		x0 = x+na;
	}
	q1 = sqrt(x0*x0+y*y);
	th = atan(y/x0);
	gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
	gi = th*(x0-0.5)+y*log(q1)-y;
	for (k=0;k<10;k++){
		t = pow(q1,-1.0-2.0*k);
		gr += (a[k]*t*cos((2.0*k+1.0)*th));
		gi -= (a[k]*t*sin((2.0*k+1.0)*th));
	}
	if (x <= 7.0) {
		gr1 = 0.0;
		gi1 = 0.0;
		for (j=0;j<na;j++) {
			gr1 += (0.5*log((x+j)*(x+j)+y*y));
			gi1 += atan(y/(x+j));
		}
		gr -= gr1;
		gi -= gi1;
	}


	if (x1 <= 0.0) {
		q1 = sqrt(x*x+y*y);
		th1 = atan(y/x);
		sr = -sin(M_PI*x)*cosh(M_PI*y);
		si = -cos(M_PI*x)*sinh(M_PI*y);
		q2 = sqrt(sr*sr+si*si);
		th2 = atan(si/sr);
		if (sr < 0.0) th2 += M_PI;
		gr = log(M_PI/(q1*q2))-gr;
		gi = -th1-th2-gi;
		x = x1;
		y = y1;
	}

	g0 = exp(gr);
	gr = g0*cos(gi);
	gi = g0*sin(gi);

	g = gr + ui*gi;

	return g;
}

void DblRegge30::updatePar( const AmpParameter& par) {

}
