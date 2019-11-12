//! VPlotOptimalCut plot results from optical cut search

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMath.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VPlotOptimalCut
{
    private:
    
        TFile* fFile;
        TTree* fData;
        vector< string > fListOfVariables;
        
        double getCutValue( string iVariable, int iEntryNumber );
        
    public:
    
        VPlotOptimalCut( string iFile );
        ~VPlotOptimalCut() {}
        int  findOptimalCut( int iSourceStrength = 4, double size = -99 );
        void listVariables();
        void plotAll( int iSourceStrength = 4, bool bPrint = false );
        void plotHistograms( string iVariable, int i_opt = -1, bool bPrint = false );
        void plotOptimalCuts( string iVariable, int iSourceStrength = 4, bool iMax = true, double size = -999, double iMaxObs = 1.e20 );
        
        void plotScanParameter( string iVar1, string iVar2 );
};
