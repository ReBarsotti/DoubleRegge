#if !defined(DBLREGGE30)
#define DBLREGGE30

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class DblRegge30 : public UserAmplitude< DblRegge30 >
{

 public:

 DblRegge30() : UserAmplitude< DblRegge30 >() {}
  DblRegge30( const vector< string >& args );

  ~DblRegge30(){}

  string name() const { return "DblRegge30"; }

  enum UserVars { k_s = 0, k_s1, k_s2, k_t1, k_t2, k_tmin, k_phiGJ, k_thetaGJ, k_thetaCM,
		  k_gam1_re, k_gam1_im, k_gam2_re, k_gam2_im, kNumUserVars };     
  unsigned int numUserVars() const {return kNumUserVars; }
   
  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;

  bool needsUserVarsOnly() const { return true; }
  bool areUserVarsStatic() const { return false; }

  void updatePar( const AmpParameter& par );
  
#ifdef GPU_ACCELERATION
  void launchGPUKernel ( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  bool isGPUEnabled() const { return false; }
#endif // GPU_ACCELERATION

 private:

  complex<double> cgamma( complex<double> z ) const;
  int m_fastParticle; //2 for eta, 3 for pion
  AmpParameter m_c0, m_c1, m_c2, m_c3, m_d1, m_d2, m_aPrime, m_a0,m_aPrime_2,m_a0_2,m_f0,m_f1,m_f2, m_g0,m_g1,m_g2,m_g3,m_f3,m_f4,m_f5;
  AmpParameter m_s0;

  static const char* kModule;
};

#endif
