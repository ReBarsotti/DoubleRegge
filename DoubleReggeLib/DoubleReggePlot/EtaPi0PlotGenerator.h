#if !(defined ETAPI0PLOTGENERATOR)
#define ETAPI0PLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

class Kinematics;

class EtaPi0PlotGenerator : public PlotGenerator
{

public:

  EtaPi0PlotGenerator( const FitResults& results, Option opt = kDefault );
	

  enum {
    khm12 = 0, khm13, khm23, kdltz, cosT, phiAng, PhiT, cosT_m23, Omega, cosT_phi, cosT_Phi, cosT_lab, phiAng_lab, cosT_m23_lab, phi_m23_lab, t_eta, t_pi,t2, s12, s13, s, s23,eta,cosT_cm, pion,recoil, beamE,phi_s23,phi_s12,phi_s13,
    kNumHists
  };

private:

  void projectEvent( Kinematics* kin );
  void projectEvent( Kinematics* kin, const string& reactionName );
};

#endif
