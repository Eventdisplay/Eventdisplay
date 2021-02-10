//!  VEffectiveAreaCalculatorMCHistograms filling, reading, writing of MC histograms for effective area calculation

#ifndef VEffectiveAreaCalculatorMCHistograms_H
#define VEffectiveAreaCalculatorMCHistograms_H

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

#include "VSpectralWeight.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VEffectiveAreaCalculatorMCHistograms : public TNamed
{
    private:
    
        bool     fDebug;
        
        // cuts
        bool     fMCCuts;
        double   fArrayxyoff_MC_min;
        double   fArrayxyoff_MC_max;
        
        // spectral weight calculator
        VSpectralWeight*           fSpectralWeight;                 //! backwards compatibility
        vector< VSpectralWeight* > fVSpectralWeight;                //!
        double   fMCEnergyRange_TeV_min;
        double   fMCEnergyRange_TeV_max;
        double   fMCSpectralIndex;
        
        // range in log(E) for the EffArea hisotrams
        int      fEnergyAxisBins_log10;
        double   fEnergyAxisMin_log10;
        double   fEnergyAxisMax_log10;
        
        int      checkParameters( const VEffectiveAreaCalculatorMCHistograms* );
        
    public:
    
        // azimuth and spectral index bins
        vector< double > fVMinAz;
        vector< double > fVMaxAz;
        vector< double > fVSpectralIndex;
        
        // MC histograms
        vector< vector< TH1D* > > hVEmc;                     // [spectral index][az]
        vector< vector< TProfile* > > hVEmcSWeight;          // [spectral index][az]
        vector< TH1D* > hVEmcUnWeighted;                     // [spectral index]
        
        VEffectiveAreaCalculatorMCHistograms();
        ~VEffectiveAreaCalculatorMCHistograms() {}
        
        bool      add( const VEffectiveAreaCalculatorMCHistograms* );
        bool      fill( double i_ze, TTree* i_MCData, bool iBAzimuthBins );
        double    getEnergyAxisMin_log10()
        {
            return fEnergyAxisMin_log10;
        }
        double    getEnergyAxisMax_log10()
        {
            return fEnergyAxisMax_log10;
        }
        TH1D*     getHistogram_Emc( unsigned int iAz, unsigned int iIndex );
        TProfile* getHistogram_EmcWeight( unsigned int iAz, unsigned int iIndex );
        TH1D*     getHistogram_EmcUnweighted( unsigned int iAz );
        void      initializeHistograms();
        void      initializeHistograms( vector< double > iAzMin, vector< double > iAzMax,
                                        vector< double > iSpectralIndex,
                                        int nbins = 60, double iMin = -2., double iMax = 4. );
        bool      isMCCuts()
        {
            return fMCCuts;
        }
        void      listEntries()
        {
            print();
        }
        void      print();
        bool      readFromEffectiveAreaTree( string iFile );
        bool      readFromEffectiveAreaFile( string iFile );
        bool      matchDataVectors( vector< double > iAzMin, vector< double > iAzMax, vector< double > iSpectralIndex );
        void      setCuts( double iArrayxyoff_MC_min = -99., double iArrayxyoff_MC_max = -99. );
        void      setDebug( bool iB = false )
        {
            fDebug = iB;
        }
        void      setDefaultValues( bool b90DegIntervalls = true );
        bool      setMonteCarloEnergyRange( double iMin, double iMax, double iMCIndex = 2. );
        
        ClassDef( VEffectiveAreaCalculatorMCHistograms, 13 );
};

#endif
