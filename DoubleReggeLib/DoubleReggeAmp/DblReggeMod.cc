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
#include "DblReggeMod.h"


double si[2];
double alp[2];


DblReggeMod::DblReggeMod( const vector< string >& args ) :
	UserAmplitude< DblReggeMod >( args )
{       
	assert( args.size() == 18 );
	b_par = AmpParameter( args[0] );
	c0 = AmpParameter( args[1]);
	c1 = AmpParameter( args[2]);
	c2 = AmpParameter( args[3]);
	n0 = AmpParameter( args[4]);
	n1 = AmpParameter( args[5]);
	n2 = AmpParameter( args[6]);
	d10 = AmpParameter(args[7]);
	d11 = AmpParameter( args[8]); 
	d12 = AmpParameter( args[9]);
	d20 = AmpParameter( args[10] );
	d21 = AmpParameter( args[11] );
	d22 = AmpParameter( args[12] );	
	aPrime = AmpParameter(args[13] );
	a0 = AmpParameter(args[14]);
	S0 = AmpParameter( args[15] );
	fastParticle = atoi(args[16].c_str());
	charge = atoi( args[17].c_str() );

	registerParameter( b_par );
	registerParameter( c0 );
	registerParameter( c1 );
	registerParameter( c2 );
	registerParameter( n0 );
	registerParameter( n1 );
	registerParameter( n2 );
	registerParameter(d10);
	registerParameter( d11 );
	registerParameter( d12 );
	registerParameter(d20);
	registerParameter( d21 );
	registerParameter( d22 );
	registerParameter( aPrime);
	registerParameter( a0 );
	registerParameter( S0 );

}


complex< GDouble >
DblReggeMod::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

	GDouble s12 = userVars[u_s12];
	GDouble s23 = userVars[u_s23];
	GDouble t1 = userVars[u_t1];
	GDouble s = userVars[u_s];
	GDouble u3 = userVars[u_u3];

	double param[17] = {b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0,S0};
	double inv[5] = {s, s12, s23, t1 ,u3};

	double mass2[4] = { userVars[u_beamM2], userVars[u_p1M2], userVars[u_p2M2], userVars[u_recoilM2]};

	int hel[3] = {1,-1,-1};
	std::complex<double> amp = ampEtaPi0(param, hel, inv, mass2);

	///if(abs(amp) > 35)
	//{
	//cout << "amp: " << amp << endl;
	//cout << "s12: " << s12 << endl;
	//cout << "s23: " << s23 << endl;
	//cout << "t1: " << t1 << endl;
	//cout << "u3: " << u3 << endl;
	//}
	return amp;
}

void DblReggeMod::calcUserVars( GDouble** pKin, GDouble* userVars ) const{
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );
	TLorentzVector resonance = p1 + p2;

	userVars[u_s12] = resonance.M2();
	userVars[u_s23] =  (p2 + recoil).M2();
	userVars[u_t1] = (beam - p1).M2();
	userVars[u_t2] = (beam - p2).M2();
	userVars[u_s] = (recoil + p1 + p2).M2();
	userVars[u_u3] = userVars[u_t1] + userVars[u_t2] + userVars[u_s12] - (beam.M2() + p1.M2() + p2.M2());

	userVars[u_beamM2] = beam.M2();
	userVars[u_p1M2] = p1.M2();
	userVars[u_p2M2] = p2.M2();
	userVars[u_recoilM2] = recoil.M2();

}

std::complex<double> DblReggeMod::ampEtaPi0(double par[11], int hel[3], double inv[5], double mass2[4]) const{

	std::complex<double> zero (0,0);

	//      if(abs(hel[0]*hel[1]*hel[2]) != 1 || hel[0]==0){return zero;}

	double s, s12,s23,t1,u3;
	double m12, m22, m32, ma2;
	s   = inv[0];   s12 = inv[1];   s23 = inv[2];   t1  = inv[3];   u3  = inv[4];
	ma2 = mass2[0]; m12 = mass2[1]; m22 = mass2[2]; m32 = mass2[3];
	double t2  = -t1+u3-s12+ma2+m12+m22;
	double s13 = s-s12-s23+m12+m22+m32;

		// scalar part
     // slope of Regge trajectories alpha'
	double alp0eta = aPrime*t1 + a0;
	double alp0pi0 = aPrime*t2 + a0;
	double alp1    = aPrime*u3 + a0;
	int tau[2];

	if(charge ==0){
		tau[0] = -1;
		tau[1] = -1;    // only vector exchange
	}
	else if(charge == 1){
		tau[0] = 1; 
		tau[1]= -1; //for charged channel, a2 exchange?
	}	

	si[0] = s12; si[1] = s23;
	alp[0] = alp0eta; alp[1] = alp1;

	std::complex<double> ADR1, ADR2;
	double Bot1, Bot2;

	if(fastParticle == 2){
	
		ADR1 = DoubleRegge(tau, s, si, alp); // fast eta
		Bot1 = exp(b_par*b_par*t1 )* sqrt(n0*exp(c0*c0*u3) + n1*exp(c1*c1*u3)*u3 + n2*exp(c2*c2*u3)*u3*u3) ;
	}
	else{
		ADR1 = (0,0);
		Bot1 = 0;
	}
	si[1] = s13; alp[0] = alp0pi0;
	
	if(fastParticle == 3){
		ADR2 = DoubleRegge(tau, s, si, alp); // fast pi0
		Bot2 = exp(b_par*b_par*t2) * sqrt(n0*exp(c0*c0*u3) + n1*exp(c1*c1*u3)*u3 + n2*exp(c2*c2*u3)*u3*u3);
	}
	else{
		Bot2 = 0;
		ADR2 = (0,0);
	}
	// helicity part
	double fac1 =  sqrt(-t1/mass2[2]);
	double fac2 = sqrt(-t2/mass2[2]);    // use the pion mass in both fac1 and fac2
	double fac3 = pow(-u3/4./mass2[0],abs((hel[1]-hel[2])/4.)); // hel[1,2] are twice the nucleon helicities!
	double parity = pow(-1,(hel[1]-hel[2])/2.);
	if(hel[1] == -1){fac3 = fac3*parity;}

	return fac3*(Bot1*fac1*ADR1 + Bot2*fac2*ADR2 );
}

double DblReggeMod::V12(double d0, double d1, double d2,double si[2]) const{
	double res = d0 + si[0]*d1 + si[1]*d2;
	return res;
}

std::complex<double> DblReggeMod::DoubleRegge( int tau[2], double s, double si[2], double alp[2]) const{
	std::complex<double> ui (0,1);
	// signature factors:

	std::complex<double> x0  = 1/2.*((double)tau[0] + exp(-ui*M_PI*alp[0]));
	std::complex<double> x1  = 1/2.*((double)tau[1] + exp(-ui*M_PI*alp[1]));
	std::complex<double> x01 = 1/2.*((double)tau[0]*tau[1] + exp(-ui*M_PI*(alp[0]-alp[1])));
	std::complex<double> x10 = 1/2.*((double)tau[1]*tau[0] + exp(-ui*M_PI*(alp[1]-alp[0])));
	// double Regge vertices:

	double eta = S0*s/(si[0]*si[1]);
	std::complex<double> V0 = V12(d10,d11,d12, si);
	std::complex<double> V1 = V12(d20,d21,d22, si);
	std::complex<double> up1 = pow(s/S0,alp[1])*pow(si[0]/S0,alp[0]-alp[1]);
	std::complex<double> up2 = pow(s/S0,alp[0])*pow(si[1]/S0,alp[1]-alp[0]);

	// combine pieces:


	std::complex<double> t1 =up1*x1*x01*V1;
	std::complex<double> t0 = up2*x0*x10*V0;
	//return (t0+t1)*cgamma(-alp[0],0)*cgamma(-alp[1],0);;
	return (t0+t1);
}


void
DblReggeMod::updatePar( const AmpParameter& par ){

	// could do expensive calculations here on parameter updates
}

#ifdef GPU_ACCELERATION
void DblReggeMod::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const{

	GPUDblReggeMod_exec( dimGrid, dimBlock, GPU_AMP_ARGS, b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0 ,S0, fastParticle, charge);

}
#endif

