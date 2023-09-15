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

  enum UserVars { k_s = 0, k_s1, k_s2, k_t1, k_t, k_phiGJ, k_thetaGJ, k_thetaCM,
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
  
  int m_j;
  int m_fastParticle; //2 for eta, 3 for pion
  int m_charge; // 0 for neutral, 1 for charged
  AmpParameter m_b_par, m_c0, m_c1, m_c2, m_n0, m_n1, m_n2, m_d10, m_d11, m_d12, m_d20, m_d21, m_d22, m_aPrime, m_a0;
  AmpParameter m_s0;

  static const char* kModule;
};

#endif
