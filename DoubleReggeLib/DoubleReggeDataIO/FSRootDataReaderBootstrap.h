#if !defined(FSROOTDATAREADERBOOTSTRAP)
#define FSROOTDATAREADERBOOTSTRAP

#include <set>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"
#include "TRandom2.h"
using namespace std;

class FSRootDataReaderBootstrap : public UserDataReader< FSRootDataReaderBootstrap >{

	public:

		FSRootDataReaderBootstrap() : UserDataReader< FSRootDataReaderBootstrap >(), 
		m_inFile(NULL),
		m_randGenerator( NULL ) { }

		FSRootDataReaderBootstrap( const vector< string >& args );

		string name() const { return "FSRootDataReaderBootstrap"; }

		virtual Kinematics* getEvent();

		virtual void resetSource();

		virtual unsigned int numEvents() const;

		int eventCounter() const { return m_eventCounter; }


	private:

		TFile* m_inFile;
		TTree* m_inTree;
		TTree* m_inFriendTree;
		int m_nPart;     
		int m_eventCounter;
		unsigned int m_numParticles;


		TRandom2* m_randGenerator;

		double m_EnPB;
		double m_PxPB;
		double m_PyPB;
		double m_PzPB;
		double m_EnP[50];
		double m_PxP[50];
		double m_PyP[50];
		double m_PzP[50];

		double m_weight;

		multiset< unsigned int > m_entryOrder;
		mutable multiset< unsigned int >::const_iterator m_nextEntry;

};

#endif
