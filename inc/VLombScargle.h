//! VLombScargle spectral analysis of unevenly sampled data (after Lomb/Scargle)

#ifndef VLombScargle_H
#define VLombScargle_H

#include "VFluxAndLightCurveUtilities.h"
#include "VFluxDataPoint.h"
#include "VLightCurveAnalyzer.h"
#include "VHistogramUtilities.h"
#include "VPlotUtilities.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TText.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VLombScargle : public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        bool fDebug;
        
        // central data object
        vector< VFluxDataPoint > fFluxDataVector;
        
        vector< double > fVFrequency;
        vector< double > fVPeriodigram;
        
        TGraph*      fPeriodigramGraph;
        TH1D*        fPeriodigramHisto;
        TCanvas*     fPeriodigramCanvas;
        
        unsigned int fNFrequencies;
        double       fFrequency_min;
        double       fFrequency_max;
        
        TRandom3*    fRandom;
        
        vector< double > fProbabilityLevels;
        vector< int >    fProbabilityLevelDigits;
        
        void   reset();
        
    public:
    
        VLombScargle();
        VLombScargle( vector< VFluxDataPoint > iDataVector );
        ~VLombScargle() {}
        
        void    fillPeriodigram( bool iShuffle = false );
        TGraph* getPeriodigramGraph();
        TH1D*   getPeriodigramHistogram( string iName );
        void    plotFrequencyLine( double iFrequencyLine_plot = -99., int iColor = 4 );
        void    plotPeriodigram( string iXTitle = "", string iYTitle = "", bool bLogX = true );
        void    plotProbabilityLevels( bool iPlotinColor = false );
        void    plotProbabilityLevelsFromToyMC( unsigned int iMCCycles = 500, unsigned int iSeed = 0, bool iPlotinColor = false );
        void    setDataVector( vector< VFluxDataPoint > iDataVector )
        {
            fFluxDataVector = iDataVector;
        }
        void    setFrequencyRange( unsigned int iNFrequencies = 1000, double iFrequency_min = 1. / 1000., double iFrequency_max = 1. / 10. );
        void    setProbabilityLevels( vector< double > iProbabilityLevels );
        void    setProbabilityLevels( vector< double > iProbabilityLevels, vector< int > iProbabilityLevelDigits );
};


#endif
