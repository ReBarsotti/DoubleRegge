
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"
#include "DblReggeHelper.cuh"
//#include "AMPTOOLS_AMPS/DblReggeHelper.cuh"

__global__ void
DblReggeMod_kernel(GPU_AMP_PROTO, GDouble S0, GDouble b_par, GDouble c0, GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble a0, GDouble a1, GDouble a2, GDouble a3, GDouble a4, int fastParticle, int charge){

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
	GDouble phi = GPU_UVARS(12);
//	GDouble up2 = GPU_UVARS(11);



	WCUComplex amp =  GPU_calcAmplitude(phi,s, s12, s23, t1, u3, S0, b_par,c0,c1,c2,n0,n1,n2,a0,a1,a2,a3,a4,beamM2, p1M2, p2M2, recoilM2, fastParticle, charge);

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
		GDouble S0, GDouble b_par, GDouble c0, GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble a0, GDouble a1, GDouble a2, GDouble a3, GDouble a4,int fastParticle, int charge  )
{

	DblReggeMod_kernel<<< dimGrid, dimBlock >>>( GPU_AMP_ARGS, S0, b_par, c0,c1,c2,n0,n1,n2, a0,a1,a2,a3,a4,fastParticle, charge );
}

