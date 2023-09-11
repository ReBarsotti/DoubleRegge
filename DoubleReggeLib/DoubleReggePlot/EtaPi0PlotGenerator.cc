#include "EtaPi0PlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

EtaPi0PlotGenerator::EtaPi0PlotGenerator( const FitResults& results, Option opt ) :
	PlotGenerator( results, opt )
{
	bookHistogram( khm12, new Histogram1D( 60, 1.0, 3.0, "hm12", "Mass( 1 2 )" ) );
	bookHistogram( khm13, new Histogram1D( 60, 1.0, 3.0, "hm13", "Mass( 1 3 )" ) );
	bookHistogram( khm23, new Histogram1D( 60, 1.0, 3.0, "hm23", "Mass( 2 3 )" ) );
	bookHistogram( kdltz, new Histogram2D( 80, 0.0, 25.0, 80, 0.0, 9.0, "dltz", "Dalitz Plot" ) );

	bookHistogram( eta, new Histogram1D(60, 0.0 , 1.0,"eta", "Mass of #eta" ));
	bookHistogram( pion, new Histogram1D(60, 0.0 , 1.0,"pion", "Mass of #pi" ));
	bookHistogram( recoil, new Histogram1D(60, 0.0 , 1.0,"recoil", "Mass of recoil" ));


	bookHistogram( s12, new Histogram1D(100, 0., 14., "s12", "S(1 2) (Vincent's notation)") );
	bookHistogram( s13, new Histogram1D(100, 0., 14., "s13", "S(1 3) (Vincent's notation)") );
	bookHistogram( s23, new Histogram1D(100, 0., 14., "s23", "S(2 3) (Vincent's notation)") );
	bookHistogram( t1, new Histogram1D(60, -10., 1., "t1", "t1 (Vincent's notation)"));
	bookHistogram( t2, new Histogram1D(60, -10., 1., "t2", "t2 (Vincent's notation)"));
	bookHistogram( u3, new Histogram1D(60, -6., 1., "u3", "u3"));
	bookHistogram(beamE, new Histogram1D(60,8.0, 9.0, "beamE", "Beam Energy"));

	bookHistogram( cosT, new Histogram1D( 60, -1.1, 1.1, "cosT", "CosTheta") );
	bookHistogram( phiAng, new Histogram1D(40, -3.2, 3.2, "phiAng", "#phi") );
	bookHistogram( cosT_lab, new Histogram1D( 60, -1.1, 1.1, "cosT_lab", "CosThetaLab") );
	bookHistogram( cosT_cm, new Histogram1D( 60, -1.1, 1.1, "cosT_cm", "CosThetaCM") );
	bookHistogram( phiAng_lab, new Histogram1D(40, -3.2, 3.2, "phiAng_lab", "#phi_{lab}") );
	bookHistogram( cosT_m23_lab, new Histogram2D(200, 0.6, 2.6, 100, -1.0, 1.0, "cosTLab_m23", "cos(#theta_{lab}) vs. Mass(#eta #pi^{0})" ) );
	bookHistogram( phi_m23_lab, new Histogram2D(200, 0.6, 2.6, 100, -3.2, 3.2, "PhiLab_m23", "#phi_{lab} vs. Mass(#eta #pi^{0})" ) );
	bookHistogram( PhiT, new Histogram1D(80, -3.2, 3.2, "PhiT", "#Phi") );
	bookHistogram( Omega, new Histogram1D(100, -3.2, 3.2, "Omega", "#Omega") );
	bookHistogram( cosT_m23, new Histogram2D(100, 1.8, 3.2, 100, -1.0, 1.0, "cosT_m23", "cos(#theta) vs. Mass(#eta #pi^{0})" ) );
	bookHistogram( cosT_phi, new Histogram2D(40, -3.2, 3.2, 100, -1.0, 1.0, "cosT_phi", "cos(#theta) vs. #phi" ) );
	bookHistogram( cosT_Phi, new Histogram2D(40, -3.2, 3.2, 100, -1.0, 1.0, "cosT_Phi", "cos(#theta) vs. #Phi" ) );
}

void EtaPi0PlotGenerator::projectEvent( Kinematics* kin ){

	// this function will make this class backwards-compatible with older versions
	// (v0.10.x and prior) of AmpTools, but will not be able to properly obtain
	// the polariation plane in the lab when multiple orientations are used

	projectEvent( kin, "" );
}

void EtaPi0PlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){

	// obtain the polarzation angle for this event by getting the list of amplitudes
	// associated with this reaction -- we know all are Zlm amplitudes
	// take the sixth argument of the first factor of the first amplitude in the first sum
	//  double polAngle = stod(cfgInfo()->amplitudeList( reactionName, "", "" ).at(0)->factors().at(0).at(5));
	double polAngle = 45.;
	double s1,s2,teta,tpi;
	TLorentzVector P0 = kin->particle(0); //beam
	TLorentzVector P1 = kin->particle(1); //proton
	TLorentzVector P2 = kin->particle(2); //eta
	TLorentzVector P3 = kin->particle(3); //pi0
	//	TLorentzVector recoil = kin->particle(4); //recoil proton


	fillHistogram( khm12, (P1+P2).M() );
	fillHistogram( khm13, (P1+P3).M() );
	fillHistogram( khm23, (P2+P3).M() );
	fillHistogram( kdltz, (P1+P2).M2(), (P2+P3).M2() );
	fillHistogram( beamE, P0.E());

	TLorentzVector proton = {0,0,0,0};
	proton.SetE(.93827231);

	fillHistogram(u3, (proton - P1).M2() );


	fillHistogram( s12, (P2 + P3).M2() );
	fillHistogram( s23, (P3 + P1).M2() );
	fillHistogram( s13, (P2 + P1).M2() );
	fillHistogram( t1, (P0 - P2).M2() );
	fillHistogram( t2, (P0 - P3).M2() );


//MODEL TESTING CODE, COMMENT OUT FOR REGULAR USE

s1 = (P2+P3).M();
s2 = (P3+P1).M2();
teta = (P0-P2).M();
tpi = (P0-P3).M();

if((s1<=2.6 && s1>=2.5) && (teta>= 0.36 && teta<= 0.49) && (tpi>= 16.3 && tpi<=17.3)){
        fillHistogram(eta, (1.1*2*proton.M()*P0.E())/(s1*s2));
        fillHistogram( s23, (P3 + P1).M2() );

} 

//END OF MODEL TESTING CODE	


//	fillHistogram(eta, P2.M());
	fillHistogram(pion, P3.M());
	fillHistogram(recoil, P1.M());

	TLorentzVector cm = proton +  P0;
	//cout << cm.Px() << ","<<cm.Py() <<"," <<cm.Pz() <<"," <<cm.E() <<endl;
	TLorentzRotation cmBoost( -cm.BoostVector() );
	TLorentzVector recoil_cm = cmBoost*P1;
	TLorentzVector target_cm = cmBoost*proton;
	GDouble cm_theta = recoil_cm.Angle(target_cm.Vect());
	GDouble costheta_cm = TMath::Cos(cm_theta);
	fillHistogram(cosT_cm, costheta_cm);


	TLorentzVector resonance = P2 + P3;
	TLorentzRotation resRestBoost( -resonance.BoostVector() );

	TLorentzVector beam_res   = resRestBoost * P0;
	TLorentzVector recoil_res = resRestBoost * P1;
	TLorentzVector p2_res = resRestBoost * P2;

	// Angles of etapi system in LAB frame:
	GDouble locCosT_lab = resonance.CosTheta();
	GDouble locPhi_lab = resonance.Phi();


	fillHistogram( cosT_lab, locCosT_lab );
	fillHistogram( phiAng_lab, locPhi_lab );
	fillHistogram( cosT_m23_lab, resonance.M(), locCosT_lab );
	fillHistogram( phi_m23_lab, resonance.M(), locPhi_lab );

	// Helicity Frame:
	// TVector3 z = -1. * recoil_res.Vect().Unit();
	// TVector3 y = (P0.Vect().Unit().Cross(-P1.Vect().Unit())).Unit();
	TVector3 z = beam_res.Vect().Unit();
	TVector3 y = recoil_res.Vect().Cross(-z).Unit();
	TVector3 x = y.Cross(z);

	TVector3 angles( (p2_res.Vect()).Dot(x),
			(p2_res.Vect()).Dot(y),
			(p2_res.Vect()).Dot(z) );

	Double_t cosTheta = angles.CosTheta();
	Double_t phi = angles.Phi();

	TVector3 eps(cos(polAngle*TMath::DegToRad()), sin(polAngle*TMath::DegToRad()), 0.0); // beam polarization vector
	GDouble Phi = atan2(y.Dot(eps), P0.Vect().Unit().Dot(eps.Cross(y)));



	TVector3 eta = P2.Vect();
	GDouble omega = atan2(y.Dot(eta), P0.Vect().Unit().Dot(eta.Cross(y)));


	fillHistogram( PhiT, Phi); 
	fillHistogram( cosT, cosTheta);
	fillHistogram( phiAng, phi);
	fillHistogram( Omega, omega);
	fillHistogram( cosT_m23, (P2+P3).M(), cosTheta);
	fillHistogram( cosT_phi, phi, cosTheta);
	fillHistogram( cosT_Phi, Phi, cosTheta);
}
