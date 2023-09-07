
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
#include "DblReggeHelper.cuh"
//#include "AMPTOOLS_AMPS/DblReggeHelper.cuh"

__global__ void
DblReggeMod_kernel(GPU_AMP_PROTO, GDouble b_par, GDouble c0, GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble d10, GDouble d11, GDouble d12,GDouble d20, GDouble d21,GDouble d22, GDouble aPrime, GDouble a0, GDouble S0, int fastParticle, int charge){

	int iEvent = GPU_THIS_EVENT;

	// here we need to be careful to index the user-defined
	// data with the proper integer corresponding to the
	// enumeration in the C++ header file

	//user vars as defined in enum in header:


	GDouble s12 = GPU_UVARS(0);
	GDouble s23 = GPU_UVARS(1);
	GDouble t1 = GPU_UVARS(2);
//	GDouble t2 = GPU_UVARS(3);
	GDouble s = GPU_UVARS(4);
	GDouble u3 = GPU_UVARS(5);
	GDouble beamM2 = GPU_UVARS(6);
	GDouble p1M2 = GPU_UVARS(7);
	GDouble p2M2 = GPU_UVARS(8);
	GDouble recoilM2 = GPU_UVARS(9);
//	GDouble up1 = GPU_UVARS(10);
//	GDouble up2 = GPU_UVARS(11);



	WCUComplex amp =  GPU_calcAmplitude(s, s12, s23, t1, u3, b_par, c0,c1,c2,n0,n1,n2, d10, d11,d12,d20,d21,d22,aPrime,a0, S0,beamM2, p1M2, p2M2, recoilM2, fastParticle, charge);

//if((amp.Re() + amp.Im()) > 35)
//{
//printf( "amp: %f \n", amp);
//printf( "u3: %f \n", u3);
//printf( "t1: %f \n", t1);
//printf( "s23: %f \n", s23);
//printf( "s12: %f \n", s12);
//}	

	pcDevAmp[iEvent] = amp;
}



void GPUDblReggeMod_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		GDouble b_par,GDouble c0, GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble d10, GDouble d11,GDouble d12,GDouble d20, GDouble d21, GDouble d22, GDouble aPrime, GDouble a0, GDouble S0, int fastParticle, int charge  )
{

	DblReggeMod_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS, b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0,S0,fastParticle, charge );
}

