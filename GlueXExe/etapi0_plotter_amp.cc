#include <iostream>
#include <string>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"


#include "GlueXPlot/EtaPi0PlotGenerator.h"
#include "GlueXDataIO/FSRootDataReader.h"
#include "GlueXDataIO/ROOTDataReader.h"
//#include "AMPTOOLS_AMPS/TwoPiAngles.h"
//#include "GlueXAmp/Zlm.h"
//#include "AMPTOOLS_AMPS/dblRegge.h"
//#include "AMPTOOLS_AMPS/TwoPiAngles_amp.h"
//#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "GlueXAmp/DblReggeMod.h"
//#include "GlueXAmp/DblRegge_FastEta.h"
//#include "GlueXAmp/DblRegge_FastPi.h"
//#include "GlueXAmp/Uniform.h"
typedef EtaPi0PlotGenerator PlotGen;
PlotGen::Option option = PlotGen::kDefault;
//PlotGen::Option option = PlotGen::kNoGenMC;

void atiSetup(){
  
//  AmpToolsInterface::registerAmplitude( dblRegge() );
 // AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerAmplitude( DblReggeMod() );  // Modification to go into dblRegge.cc ? Will check...
 // AmpToolsInterface::registerAmplitude( Uniform() );  
  //AmpToolsInterface::registerAmplitude( DblRegge_FastEta() ); 
  //AmpToolsInterface::registerAmplitude( DblRegge_FastPi() ); 
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );
}


//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter *** " << endl << endl;

  if (argc <= 1){
    
    cout << "Usage:" << endl << endl;
    cout << "\tampPlotter <fit results name>" << endl << endl;
    return 0;
  }

    // ************************
    // parse the command line parameters
    // ************************

  string resultsName(argv[1]);
  FitResults results( resultsName );
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }

    // ************************
    // set up the plot generator
    // ************************

  atiSetup();  
  PlotGen plotGen( results, option );

    // ************************
    // start the GUI
    // ************************

  cout << ">> Plot generator ready, starting GUI..." << endl;

  int dummy_argc = 0;
  char* dummy_argv[] = {};  
  TApplication app( "app", &dummy_argc, dummy_argv );
  
  gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetFillStyle(1001);
  gStyle->SetPalette(1);
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameFillStyle(1001);
  
  PlotFactory factory( plotGen );	
  PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
	
  app.Run();
    
  return 0;

}

