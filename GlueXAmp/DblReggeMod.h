#if !defined(DBLREGGEMOD)
#define DBLREGGEMOD

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPUDblReggeMod_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, GDouble b_par,GDouble c0, GDouble c1, GDouble c2, GDouble n0, GDouble n1, GDouble n2, GDouble d10, GDouble d11,GDouble d12, GDouble d20, GDouble d21, GDouble d22, GDouble aPrime, GDouble a0, GDouble S0, int fastParticle, int charge );
#endif //GPU_ACCELERATION


using std::complex;
using namespace std;

class Kinematics;

class DblReggeMod : public UserAmplitude< DblReggeMod >
{

public:

        DblReggeMod() : UserAmplitude< DblReggeMod >() { };
        DblReggeMod( const vector< string >& args );

        string name() const { return "DblReggeMod"; }

        enum UserVars {u_s12=0,u_s23=1,u_t1=2,u_t2=3,u_s=4,u_u3=5,u_beamM2=6, u_p1M2=7, u_p2M2=8, u_recoilM2=9,u_up1=10, u_up2=11, kNumUserVars };
        
        unsigned int numUserVars() const {return kNumUserVars; }
	

	complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
        void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
        void updatePar( const AmpParameter& par );
        double V12(double d0, double d1, double d2,double si[2]) const;
        std::complex<double> DoubleRegge(int tau[2], double s, double si[2], double alp[2]) const;
        std::complex<double> ampEtaPi0(double par[13], int hel[3],  double inv[5], double mass2[4]) const;

        bool needsUserVarsOnly() const { return true; }
	bool areUserVarsStatic() const { return true; }
       // void init();

#ifdef GPU_ACCELERATION
	void launchGPUKernel ( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
	bool isGPUEnabled() const { return true; }
#endif // GPU_ACCELERATION

private:
  	int j;
        int fastParticle; //2 for eta, 3 for pion
	int charge; // 0 for neutral, 1 for charged
        AmpParameter b_par,c0,c1,c2,n0,n1,n2,d10,d11,d12,d20,d21,d22,aPrime,a0;
        AmpParameter S0;
};

#endif
                    
