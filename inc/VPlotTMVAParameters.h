//! VPlotTMVAParameters plot results from TMVA optimization

#ifndef VPlotTMVAParameters_H
#define VPlotTMVAParameters_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"

#include "VTMVAEvaluator.h"

using namespace std;

class VPlotTMVAParameters
{
    private:
    
        vector< string > fSubArrays;
        string  fDataDirectory;
        
        vector< TH1D* > hSignalEfficiency;
        vector< TH1D* > hBackgroundEfficiency;
        vector< TH1D* > hMVA;
        
        bool initializeHistograms( unsigned int iEnergyWeightFileIndex_min, unsigned int iEnergyWeightFileIndex_max, unsigned int iZenithWeightFileIndex_min, unsigned int iZenithWeightFileIndex_max );
        
    public:
        VPlotTMVAParameters();
        ~VPlotTMVAParameters() {}
        
        void initializeWeightFiles( string iDirectory, string iTMVADirectory, unsigned int iEnergyWeightFileIndex_min, unsigned int iEnergyWeightFileIndex_max, unsigned int iZenithWeightFileIndex_min, unsigned int iZenithWeightFileIndex_max, double iParticleNumberFile_Conversion_Rate_to_seconds = 60. );
        void plot( bool iPrint = false );
        void setDirectories( string iDataDirectory )
        {
            fDataDirectory = iDataDirectory;
        }
        bool setSubArrays( string iSubarrayFile );
};



#endif
