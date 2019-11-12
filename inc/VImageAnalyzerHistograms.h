//! VImageAnalyzerHistograms  histogramming of run parameter, etc.

#ifndef VImageAnalyzerHistograms_H
#define VImageAnalyzerHistograms_H

#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VImageAnalyzerHistograms
{
    private:
        unsigned int fTelescopeID;
        
    public:
        VImageAnalyzerHistograms( unsigned int iTel );
        ~VImageAnalyzerHistograms();
        void init();                              //!< book histograms
        void fillL2DiagnosticTree( int rN, int eN, int iMJD, double it, vector< double >& iF, vector< double >& iS );
        void terminate( TFile* );                 //!< write results to same file as VAnalyzer class
        
        TList* hisList;
        
        TTree* fdiagno;
        int runNumber;
        int eventNumber;
        int MJD;
        float time;
        float fFADCstopTZero[4];
        float fFADCstopSum[4];
        
};
#endif                                            // VImageAnalyzerHistograms_H
