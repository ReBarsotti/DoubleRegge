#ifndef CUDA_DBLREGGEHELPER
#define CUDA_DBLREGGEHELPER

#include "GPUManager/GPUCustomTypes.h"
#include <stdio.h>


static __device__ GDouble V12(GDouble d0, GDouble d1, GDouble d2,GDouble si[2]) {

	GDouble res =d0 + si[1]*d1 + si[0]*d2;

	return res;
}



static __device__ WCUComplex GPU_DoubleRegge(int tau[2], GDouble s, GDouble si[2], GDouble alp[2], GDouble S0, GDouble d10, GDouble d11,GDouble d12, GDouble d20, GDouble d21, GDouble d22){

	WCUComplex ui=  {0,1};
	// signature factors:

	GDouble real,realE,imagE, imag;
	real = (-1.*ui*M_PI*alp[0]).Re();
	imag = (-1.*ui*M_PI*alp[0]).Im();

	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step = {realE, imagE};
	WCUComplex x0  = 1/2.*((double)tau[0] +step);

	real = (-1.*ui*M_PI*alp[1]).Re();
	imag = (-1.*ui*M_PI*alp[1]).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step1 = {realE, imagE};

	WCUComplex x1  = 1/2.*((double)tau[1] + step1);

	real = (-1.*ui*M_PI*(alp[0] - alp[1])).Re();
	imag = (-1.*ui*M_PI*(alp[0]-alp[1])).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step2 = {realE, imagE};

	WCUComplex x01 = 1/2.*((double)tau[0]*tau[1] + step2);

	real = (-1.*ui*M_PI*(alp[1] - alp[0])).Re();
	imag = (-1.*ui*M_PI*(alp[1]-alp[0])).Im();
	realE = exp(real)*cos(imag);
	imagE = exp(real)*sin(imag);
	WCUComplex step3 = {realE, imagE};
	WCUComplex x10 = 1/2.*((double)tau[1]*tau[0] + step3);
	// double Regge vertices:

	double eta = S0*s/(si[0]*si[1]);
	GDouble V0 = V12(d10,d11,d12, si);
	GDouble V1 = V12(d20,d21,d22, si);

	GDouble up1 = pow(s/S0,alp[1])*pow(si[0]/S0,alp[0]-alp[1]);
	GDouble up2 = pow(s/S0,alp[0])*pow(si[1]/S0,alp[1]-alp[0]);

	// combine pieces:


	WCUComplex  t1 =up1*x1*x01*V1;
	WCUComplex t0 = up2*x0*x10*V0;

	return (t0+t1);


}


static __device__ WCUComplex GPU_ampEtaPi0(GDouble par[11], int hel[3], GDouble inv[5], GDouble mass2[4], GDouble b_par,GDouble c0, GDouble c1,GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble d10, GDouble d11, GDouble d12, GDouble d20,GDouble d21, GDouble d22, GDouble aPrime,GDouble a0, GDouble S0, int fastParticle, int charge ){

	WCUComplex zero =  {0,0};

	if(abs(hel[0]*hel[1]*hel[2]) != 1 || hel[0]==0){return zero;}

	GDouble s, s12,s23,t1,u3;
	GDouble m12, m22, m32, ma2;
	s   = inv[0];   s12 = inv[1];   s23 = inv[2];   t1  = inv[3];   u3  = inv[4];
	ma2 = mass2[0]; m12 = mass2[1]; m22 = mass2[2]; m32 = mass2[3];
	GDouble t2  = -t1+u3-s12+ma2+m12+m22;
	GDouble s13 = s-s12-s23+m12+m22+m32;

	// scalar part
	GDouble app = aPrime;     // slope of Regge trajectories alpha'
	GDouble alp0eta = app*t1 + a0;
	GDouble alp0pi0 = app*t2 + a0;
	GDouble alp1    = app*u3 + a0;

	//	int tau[2] = {-1, -1};    // only vector exchange
	int tau[2];

	if(charge ==0){
		tau[0] = -1;
		tau[1] = -1;    // only vector exchange
	}
	else if(charge == 1){
		tau[0] = 1;
		tau[1]= -1; //for charged channel, a2 exchange?
	}



	GDouble si[2]; 
	GDouble alp[2];
	si[0] = s12; si[1] = s23;
	alp[0] = alp0eta; alp[1] = alp1;

	WCUComplex ADR1,ADR2;
	GDouble Bot1, Bot2;

	if(fastParticle==2){
		ADR1 = GPU_DoubleRegge(tau, s, si, alp,S0, d10,d11,d12,d20,d21,d22); // fast eta
		Bot1 = exp(b_par*b_par*t1 )*sqrt(n0*exp(c0*c0*u3) + n1*exp(c1*c1*u3)*u3 + n2*exp(c2*c2*u3)*u3*u3);
	}
	else{
		Bot1 = 0;
		ADR1 = {0,0};
	}

	si[1] = s13; alp[0] = alp0pi0;

	if(fastParticle==3){
		ADR2 = GPU_DoubleRegge(tau, s, si, alp, S0, d10, d11, d12, d20,d21,d22); // fast pi0
		Bot2 = exp(b_par*b_par*t2) * sqrt(n0*exp(c0*u3) + n1*exp(c1*u3)*u3 + n2*exp(c2*u3)*u3*u3);
	}
	else{
		Bot2 = 0;
		ADR2 = {0,0};
	}

	// helicity part
	GDouble fac1 =  sqrt(-t1/mass2[2]);
	GDouble fac2 = sqrt(-t2/mass2[2]);    // use the pion mass in both fac1 and fac2
	GDouble fac3 = pow(-u3/4./mass2[0],abs((hel[1]-hel[2])/4.)); // hel[1,2] are twice the nucleon helicities!
	GDouble parity = pow(-1.0,(hel[1]-hel[2])/2.);
	if(hel[1] == -1){fac3 = fac3*parity;}

	WCUComplex finalFactor = fac3*(Bot1*fac1*ADR1 + Bot2*fac2*ADR2 );
	return finalFactor;
}


static __device__ WCUComplex GPU_calcAmplitude(GDouble s, GDouble s12, GDouble s23, GDouble t1, GDouble u3, GDouble b_par,GDouble c0,GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble d10, GDouble d11,GDouble d12,GDouble d20,GDouble d21, GDouble d22,GDouble aPrime,GDouble a0, GDouble S0, GDouble beamM2, GDouble p1M2, GDouble p2M2, GDouble recoilM2, int fastParticle, int charge){

	GDouble param[17] = {b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0, S0 };
	GDouble inv[5] = {s, s12, s23, t1 ,u3};
	GDouble mass2[4] = {beamM2, p1M2, p2M2, recoilM2};

	int hel[3] = {1,-1,-1};

	WCUComplex amp = GPU_ampEtaPi0(param, hel, inv, mass2, b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0, S0, fastParticle, charge);


	return amp;
}


#endif
