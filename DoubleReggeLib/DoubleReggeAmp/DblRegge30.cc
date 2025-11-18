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


//	complex< double > gam1 = cgamma( -alpha1 );
//	complex< double > gam2 = cgamma( -alpha2 );

	double eta = m_s0*u[k_s] / (u[k_s2]*u[k_s1]);



//	complex< double > v1 = (cgamma(alpha1-alpha2)/gam2)*CHGM(-alpha1 , 1-alpha1+alpha2,-1./eta);
//	complex< double > v2 =( cgamma(alpha2-alpha1)/gam1)*CHGM(-alpha2 , 1-alpha2+alpha1, -1./eta);

	complex< double> v1 = (m_c0 + m_c1*(1./eta) + m_c2*(1./eta)*(1./eta) + m_c3*(1./eta)*(1./eta)*(1./eta));
	complex< double> v2 = (m_d0 + m_d1*(1./eta) + m_d2*(1./eta)*(1./eta) + m_d3*(1./eta)*(1./eta)*(1./eta));

//	complex< double> v1 = 1 + m_c0*cos(u[k_phiGJ]);
//	complex< double> v2 = 1 + m_d0*cos(u[k_phiGJ]);

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


	// careful here!!  if you decide to float a0 or a0Prime in the fit then
	// these precalculated values cannot be used
	//complex<double> gam1( u[k_gam1_re], u[k_gam1_im] );
	//	complex<double> gam2( u[k_gam2_re], u[k_gam2_im] );

	//  these two lines would be appropriate if a0 and/or a0Prime and hence
	//  the values of alpha are floating in the fit...


	//GAMMAS REMOVED
	return barrier*((ss2*xi1*xi21*v1 + ss1*xi2*xi12*v2 ));
//	return gam1*gam2*barrier*((ss2*xi1*xi21*v1 + ss1*xi2*xi12*v2 ));
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
	double alpha1 =  m_aPrime*m_aPrime*userVars[k_t1] + m_a0;
	double alpha2 =  m_aPrime_2*m_aPrime_2*userVars[k_t]  + m_a0_2;

//	complex< double > gam1 = cgamma( -alpha1 );
//	complex< double > gam2 = cgamma( -alpha2 );

//	userVars[k_gam1_re] = real( gam1 );
//	userVars[k_gam1_im] = imag( gam1 );

//	userVars[k_gam2_re] = real( gam2 );
//	userVars[k_gam2_re] = imag( gam2 );
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
double DblRegge30::CHGM(double A, double B, double X) const{
	double A0=A, A1=A, X0=X, HG = 0.0;
	double TBA, TB, TA, Y0=0.0, Y1=0.0, RG, LA = (int) A, NL, R, M, INF = pow(10,300);
	double sum1, sum2, R1, R2, HG1, HG2;
	if (B == 0.0 || B == -abs( (int) B)){
		HG = INF;
	} else if(A == 0.0 || X == 0.0) {
		HG = 1.0;
	} else if(A == -1.0){
		HG = 1.0 - X/B;
	} else if(A == B){
		HG = exp(X);
	} else if (A-B == 1.0){
		HG = (1.0+X/B)*exp(X);
	} else if (A == 1.0 && B == 2.0){
		HG = (exp(X)-1.0)/X;
	} else if(A == (int)A && A < 0.0){
		M = (int) -A;
		R = 1.0;
		HG = 1.0;
		for (int k = 1; k<= M ; k++) {
			R = R*(A+k-1.0)/k/(B+k-1.0)*X;
			HG+=R;
		}
	}
	if(HG != 0){return HG;}

	if(X<0.0){
		A = B-A;
		A0 = A;
		X = fabs(X);
	}
	if(A<2.0) {NL = 0;}
	else{
		NL = 1;
		LA = (int) A;
		A  = A-LA-1.0;
	}
	for (int n = 0; n<= NL; n++) {
		if(A0 >= 2.0 ) { A+=1.0; }
		if(X <= 30.0 + fabs(B) || A < 0.0){
			HG = 1.0;
			RG = 1.0;
			for (int j = 1; j<= 500; j++) {
				RG = RG*(A+j-1)/(j*(B+j-1))*X;
				HG += RG;
				if(fabs(RG/HG) < pow(10.,-15.)) {
					if(n==0) {Y0 = HG;}
					if(n==1) {Y1 = HG;}
				}
				continue;
			}
		} else {
			TA = tgamma(A);
			TB = tgamma(B);
			TBA = tgamma(B-A);
			sum1 = 1.0;
			sum2 = 1.0;
			R1 = 1.0;
			R2 = 1.0;
			for (int i = 1; i<=8; i++) {
				R1 = - R1*(A+i-1)*(A-B+i)/(X*i);
				R2 = - R2*(B-A+i-1)*(A-i)/(X*i);
				sum1+=R1;
				sum2+=R2;
			}
			HG1 = TB/TBA*pow(X,-A)*cos(M_PI*A)*sum1;
			HG2 = TB/TA*exp(X)*pow(X,A-B)*sum2;
			HG = HG1+HG2;
		}
		if(n==0) {Y0 = HG;}
		if(n==1) {Y1 = HG;}
	}
	if(A0 >= 2.0){
		for (int i=1; i<=LA-1; i++) {
			HG = ((2.*A-B+X)*Y1+(B-A)*Y0)/A;
			Y0 = Y1;
			Y1 = HG;
			A += 1.;
		}
	}
	if(X0<0.0) {HG = HG*exp(X0);}
	A = A1;
	X = X0;

	return HG;
}
