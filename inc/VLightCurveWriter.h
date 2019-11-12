//! VLightCurveWriter write / print light curves in different formats
#ifndef VLightCurveWriter_H
#define VLightCurveWriter_H

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include "VFluxAndLightCurveUtilities.h"
#include "VFluxDataPoint.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

class VLightCurveWriter
{
    private:
    
        bool fDebug;
        
        // central data object
        vector< VFluxDataPoint > fFluxDataVector;
        
    public:
    
        VLightCurveWriter();
        VLightCurveWriter( vector< VFluxDataPoint > iDataVector );
        ~VLightCurveWriter() {}
        
        void setDataVector( vector< VFluxDataPoint > iDataVector );
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        bool writeDataVectorToRootFile( string i_root_file );
        
        bool writeASCIIt_SimpleLongTable( string ASCIIFile, double iMultiplier = 1. );
        void writeDCFFormat();
        void writeWikiFormat();
        void writeTexFormat_TexTableRow( double iSigmaMinFluxLimits = -99., double iFluxMultiplicator = 1., bool iPrintPhaseValues = false );
        bool writeTexFormat_SimpleLongTable( string iTexFile, double iMultiplier = 1. );
        
};

#endif
