#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TFile.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"

#include "DoubleReggePlot/EtaPi0PlotGenerator.h"
#include "DoubleReggeDataIO/FSRootDataReader.h"
#include "DoubleReggeAmp/DblReggeMod.h"
#include "DoubleReggeAmp/DblRegge30.h"

typedef EtaPi0PlotGenerator PlotGen;

#include "IUAmpTools/report.h"
static const char* kModule = "plotResults";

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){


  // ************************
  // usage
  // ************************

  if (argc < 3){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\tplotResults <fit results name> <output file name>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Plotting Results from the Fit *** " << endl << endl;


  // ************************
  // parse the command line parameters
  // ************************

  string resultsname(argv[1]);
  string outname(argv[2]);

  report( INFO, kModule ) << "Fit results file name    = " << resultsname << endl;
  report( INFO, kModule ) << "Output file name    = " << outname << endl << endl;


  // ************************
  // load the results and display the configuration info
  // ************************

  FitResults results( resultsname );  
  string reactionName = results.reactionList()[0];

  // ************************
  // set up an output Root file
  // ************************

  TFile* plotfile = new TFile( outname.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  // ************************
  // set up a PlotGenerator and make plots
  // ************************

  AmpToolsInterface::registerAmplitude(DblReggeMod());
  AmpToolsInterface::registerAmplitude(DblRegge30());
  AmpToolsInterface::registerDataReader(FSRootDataReader());
 
  PlotGen plotGenerator( results );
  plotGenerator.enableReaction( reactionName );
  vector<string> amps = plotGenerator.uniqueAmplitudes();

  for (unsigned int iamp = 0; iamp < amps.size(); iamp++){
    
    plotGenerator.enableAmp(iamp);
  }
  
  // loop over background, data, accMC, and genMC

  for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){

    // loop over different variables

    for (unsigned int ivar  = 0; ivar  < PlotGen::kNumHists; ivar++){

      string histname =  plotGenerator.getHistogram( ivar )->name();
      histname += "_";

      if (iplot == PlotGenerator::kData) histname += "dat";
      if (iplot == PlotGenerator::kBkgnd) histname += "bkg";
      if (iplot == PlotGenerator::kAccMC) histname += "acc";
      if (iplot == PlotGenerator::kGenMC) histname += "gen";
   
      Histogram* hist = plotGenerator.projection(ivar, reactionName, iplot);
      TH1* thist = hist->toRoot();
      thist->SetName(histname.c_str());
      plotfile->cd();
      thist->Write();
    }
  }

  plotfile->Close();

  // ************************
  // print results to the screen
  // ************************
  
  report( NOTICE, kModule ) << "TOTAL EVENTS = " << results.intensity().first << " +- "
			    << results.intensity().second << endl;
  vector<string> fullamps = plotGenerator.fullAmplitudes();
  for (unsigned int i = 0; i < fullamps.size(); i++){
    vector<string> useamp;  useamp.push_back(fullamps[i]);
    report( NOTICE, kModule ) << "FIT FRACTION " << i+1 << " = "
			      << results.intensity(useamp).first /
      results.intensity().first <<  " +- "
			      << results.intensity(useamp).second /
      results.intensity().first <<  endl;
  }

  return 0;

}

