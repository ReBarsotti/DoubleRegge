#include <iostream>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"

#include "DoubleReggeDataIO/FSRootDataReader.h"
#include "DoubleReggeAmp/DblReggeMod.h"
#include "DoubleReggeAmp/DblRegge30.h"

using namespace std;

#include "IUAmpTools/report.h"
static const char* kModule = "printAmplitudes";

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  if (argc <= 1){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\tprintAmplitudes <config file name>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Printing Amplitudes *** " << endl << endl;

    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  report( INFO, kModule ) << "Config file name = " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();


    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterface::registerAmplitude(DblReggeMod());
  AmpToolsInterface::registerAmplitude(DblRegge30());
  AmpToolsInterface::registerDataReader(FSRootDataReader());

  AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kPlotGeneration);

  DataReader* dataReader = ATI.genMCReader(cfgInfo->reactionList()[0]->reactionName());
  for (int i = 0; i < 5; i++){
    Kinematics* kin = dataReader->getEvent();
    ATI.printEventDetails(cfgInfo->reactionList()[0]->reactionName(),kin);
    delete kin;
  }

}
