#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "TSystem.h"

#include "DoubleReggeDataIO/ROOTDataReader.h"
#include "DoubleReggeDataIO/ROOTDataReaderBootstrap.h"
#include "DoubleReggeDataIO/ROOTDataReaderWithTCut.h"
//#include "DoubleReggeDataIO/ROOTDataReaderTEM.h"
#include "DoubleReggeDataIO/FSRootDataReader.h"
#include "DoubleReggeAmp/DblReggeMod.h"
#include "DoubleReggeAmp/DblRegge30.h"
//#include "DoubleReggeAmp/DblRegge_FastEta.h"
//#include "DoubleReggeAmp/DblRegge_FastPi.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/report.h"

using std::complex;
using namespace std;

string kModule = "fit";

double runSingleFit(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile) {
  AmpToolsInterface ati( cfgInfo );


  double lik = ati.likelihood();
  report( NOTICE, kModule ) << "LIKELIHOOD BEFORE MINIMIZATION:  " << lik << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);

  if( useMinos ){

    fitManager->minosMinimization();
  }
  else{

    fitManager->migradMinimization();
  }

  if( hesse )
     fitManager->hesseEvaluation();

  bool fitFailed =
    ( fitManager->status() != 0 || fitManager->eMatrixStatus() != 3 );

  if( fitFailed ){
    report( ERROR, kModule ) << "ERROR: fit failed use results with caution..." << endl;
    return 1e6;
  }

  lik = ati.likelihood();
  report( NOTICE, kModule ) << "LIKELIHOOD AFTER MINIMIZATION:  " << lik << endl;

  ati.finalizeFit();

  if( seedfile.size() != 0 && !fitFailed ){
    ati.fitResults()->writeSeed( seedfile );
  }

  return ati.likelihood();
}

void runRndFits(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, int numRnd, double maxFraction) {
  AmpToolsInterface ati( cfgInfo );
  string fitName = cfgInfo->fitName();

  double lik = ati.likelihood();
  report( NOTICE, kModule ) << "LIKELIHOOD BEFORE MINIMIZATION:  " << lik << endl;

  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);

  vector< vector<string> > parRangeKeywords = cfgInfo->userKeywordArguments("parRange");

  // keep track of best fit (mininum log-likelihood)
  double minLL = 0;
  int minFitTag = -1;

  for(int i=0; i<numRnd; i++) {
    report( NOTICE, kModule ) << endl << "###############################" << endl;
    report( NOTICE, kModule ) << "FIT " << i << " OF " << numRnd << endl;
    report( NOTICE, kModule ) << endl << "###############################" << endl;

    // randomize parameters
    ati.randomizeProductionPars(maxFraction);
    for(size_t ipar=0; ipar<parRangeKeywords.size(); ipar++) {
      ati.randomizeParameter(parRangeKeywords[ipar][0], atof(parRangeKeywords[ipar][1].c_str()), atof(parRangeKeywords[ipar][2].c_str()));
    }

    if(useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if(hesse)
       fitManager->hesseEvaluation();

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

    if( fitFailed )
      report( ERROR, kModule ) << "ERROR: fit failed use results with caution..." << endl;

    lik = ati.likelihood();
    report( NOTICE, kModule ) << "LIKELIHOOD AFTER MINIMIZATION:  " << lik << endl;

    ati.finalizeFit(to_string(i));

    if( seedfile.size() != 0 && !fitFailed ){
      string seedfile_rand = seedfile + Form("_%d.txt", i);
      ati.fitResults()->writeSeed( seedfile_rand );
    }

    // update best fit
    if( !fitFailed && ati.likelihood() < minLL ) {
      minLL = ati.likelihood();
      minFitTag = i;
    }
  }

  // print best fit results
  if(minFitTag < 0) report( ERROR, kModule ) << "ALL FITS FAILED!" << endl;
  else {
    report( NOTICE, kModule ) << "MINIMUM LIKELIHOOD FROM " << minFitTag << " of " << numRnd << " RANDOM PRODUCTION PARS = " << minLL << endl;
    gSystem->Exec(Form("cp %s_%d.fit %s.fit", fitName.data(), minFitTag, fitName.data()));
    if( seedfile.size() != 0 )
      gSystem->Exec(Form("cp %s_%d.txt %s.txt", seedfile.data(), minFitTag, seedfile.data()));
  }
}

void runParScan(ConfigurationInfo* cfgInfo, bool useMinos, bool hesse, int maxIter, string seedfile, string parScan) {
  double minVal=0, maxVal=0, stepSize=0;
  int steps=0;

  vector< vector<string> > parScanKeywords = cfgInfo->userKeywordArguments("parScan");

  if(parScanKeywords.size()==0) {
    report( ERROR, kModule ) << "No parScan keyword found in configuration file. Set up at least one parameter for scanning! Aborting." << endl;
    return;
  } else {
    for(size_t ipar=0; ipar<parScanKeywords.size(); ipar++) {
      if(parScanKeywords[ipar][0]==parScan) {
	minVal = atof(parScanKeywords[ipar][1].c_str());
	maxVal = atof(parScanKeywords[ipar][2].c_str());
	stepSize = atof(parScanKeywords[ipar][3].c_str());
	steps = trunc((maxVal-minVal)/stepSize)+1;
	break;
      } else
	report( NOTICE, kModule ) << "Skipping configuration to scan " << parScanKeywords[ipar][0] << "since scanning of " << parScan << " was requested..." << endl;
    }
  }

  AmpToolsInterface ati( cfgInfo );

  string fitName = cfgInfo->fitName();
  double lik = ati.likelihood();
  report( NOTICE, kModule ) << "LIKELIHOOD BEFORE MINIMIZATION:  " << lik << endl;

  ParameterManager* parMgr = ati.parameterManager();
  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  fitManager->setMaxIterations(maxIter);


  for(int i=0; i<steps; i++) {
    report( NOTICE, kModule ) << endl << "###############################" << endl;
    report( NOTICE, kModule ) << "FIT " << i << " OF " << steps << endl;
    report( NOTICE, kModule ) << endl << "###############################" << endl;

    // reinitialize production parameters from seed file
    ati.reinitializePars();

    // set parameter to be scanned
    vector<ParameterInfo*> parInfoVec = cfgInfo->parameterList();

    auto parItr = parInfoVec.begin();
    for( ; parItr != parInfoVec.end(); ++parItr ) {
      if( (**parItr).parName() == parScan ) break;
    }

    if( parItr == parInfoVec.end() ){
      report( ERROR, kModule ) << "ERROR:  request to scan nonexistent parameter:  " << parScan << endl;
      return;
    }

    // set and fix parameter for scan
    double value = minVal + i*stepSize;
    parMgr->setAmpParameter( parScan, value );

    cfgInfo->setFitName(fitName + "_scan");

    if(useMinos)
      fitManager->minosMinimization();
    else
      fitManager->migradMinimization();

    if(hesse)
       fitManager->hesseEvaluation();

    bool fitFailed = (fitManager->status() != 0 || fitManager->eMatrixStatus() != 3);

    if( fitFailed )
      report( ERROR, kModule ) << "ERROR: fit failed use results with caution..." << endl;

    lik = ati.likelihood();
    report( NOTICE, kModule ) << "LIKELIHOOD AFTER MINIMIZATION:  " << lik << endl;

    ati.finalizeFit(to_string(i));

    if( seedfile.size() != 0 && !fitFailed ){
      string seedfile_scan = seedfile + Form("_scan_%d.txt", i);
      ati.fitResults()->writeSeed( seedfile_scan );
    }
  }
}

int main( int argc, char* argv[] ){

   // set default parameters

   bool useMinos = false;
   bool hesse = false;

   string configfile;
   string seedfile;
   string scanPar;
   int numRnd = 0;
   int maxIter = 10000;

   // parse command line

   for (int i = 1; i < argc; i++){

      string arg(argv[i]);

      if (arg == "-c"){  
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  configfile = argv[++i]; }
      if (arg == "-s"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  seedfile = argv[++i]; }
      if (arg == "-r"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  numRnd = atoi(argv[++i]); }
      if (arg == "-m"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  maxIter = atoi(argv[++i]); }
      if (arg == "-n") useMinos = true;
      if (arg == "-H") hesse = true;
      if (arg == "-p"){
         if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
         else  scanPar = argv[++i]; }
      if (arg == "-h"){
         report( ERROR, kModule ) << endl << " Usage for: " << argv[0] << endl << endl;
         report( ERROR, kModule ) << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
         report( ERROR, kModule ) << "   -H \t\t\t\t\t evaluate HESSE matrix after minimization" << endl;
         report( ERROR, kModule ) << "   -c <file>\t\t\t\t config file" << endl;
         report( ERROR, kModule ) << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
         report( ERROR, kModule ) << "   -r <int>\t\t\t Perform <int> fits each seeded with random parameters" << endl;
         report( ERROR, kModule ) << "   -p <parameter> \t\t\t\t Perform a scan of given parameter. Stepsize, min, max are to be set in cfg file" << endl;
         report( ERROR, kModule ) << "   -m <int>\t\t\t Maximum number of fit iterations" << endl; 
         exit(1);}
   }

   if (configfile.size() == 0){
      report( ERROR, kModule ) << "No config file specified" << endl;
            exit(1);
   }

   ConfigFileParser parser(configfile);
   ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
   cfgInfo->display();

 /*  AmpToolsInterface::registerAmplitude( BreitWigner() );
   AmpToolsInterface::registerAmplitude( BreitWigner3body() );
   AmpToolsInterface::registerAmplitude( TwoPSAngles() );
   AmpToolsInterface::registerAmplitude( TwoPSHelicity() );
   AmpToolsInterface::registerAmplitude( TwoPiAngles() );
   AmpToolsInterface::registerAmplitude( TwoPiAngles_amp() );
   AmpToolsInterface::registerAmplitude( TwoPiAngles_primakoff() );
   AmpToolsInterface::registerAmplitude( TwoPiWt_primakoff() );
   AmpToolsInterface::registerAmplitude( TwoPiWt_sigma() );
   AmpToolsInterface::registerAmplitude( TwoPitdist() );
   AmpToolsInterface::registerAmplitude( TwoPiNC_tdist() );
   AmpToolsInterface::registerAmplitude( ThreePiAngles() );
   AmpToolsInterface::registerAmplitude( ThreePiAnglesSchilling() );
   AmpToolsInterface::registerAmplitude( VecRadiative_SDME() );
   AmpToolsInterface::registerAmplitude( TwoLeptonAngles() );
   AmpToolsInterface::registerAmplitude( TwoLeptonAnglesGJ() );
   AmpToolsInterface::registerAmplitude( Zlm() );
   AmpToolsInterface::registerAmplitude( b1piAngAmp() );
   AmpToolsInterface::registerAmplitude( polCoef() );
   AmpToolsInterface::registerAmplitude( Uniform() );
*/  
 
   AmpToolsInterface::registerAmplitude( DblReggeMod() );
   AmpToolsInterface::registerAmplitude( DblRegge30() );
//AmpToolsInterface::registerAmplitude( DblRegge_FastEta() );
//AmpToolsInterface::registerAmplitude( DblRegge_FastPi() );
  /* AmpToolsInterface::registerAmplitude( omegapi_amplitude() );
   AmpToolsInterface::registerAmplitude( Vec_ps_refl() );
   AmpToolsInterface::registerAmplitude( PhaseOffset() );
   AmpToolsInterface::registerAmplitude( Piecewise() );
*/
   AmpToolsInterface::registerDataReader( ROOTDataReader() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
   AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
 //  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
   AmpToolsInterface::registerDataReader( FSRootDataReader() );

   if(numRnd==0){
      if(scanPar=="")
         runSingleFit(cfgInfo, useMinos, hesse, maxIter, seedfile);
      else
         runParScan(cfgInfo, useMinos, hesse, maxIter, seedfile, scanPar);
   } else {
      runRndFits(cfgInfo, useMinos, hesse, maxIter, seedfile, numRnd, 0.5);
   }

  return 0;
}


