//! VDispTable fill and read tables for angular reconstruction with disp method

#ifndef VDispTable_H
#define VDispTable_H

#include "Cshowerpars.h"
#include "Ctpars.h"

#include "VDispTableReader.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VDispTable
{
    private:
        bool fDebug;
        
        unsigned int fNTel;
        vector< float > fAz_min;
        vector< float > fAz_max;
        // 2D histograms
        vector< TProfile2D* > h2D_AzDispTable;
        vector< TH2F* >       h2D_AzDispTableN;
        vector< TProfile2D* > h2D_AzDispPhiTable;
        vector< TProfile2D* > h2D_AzDispMissTable;
        // 3D histograms
        vector< TProfile3D* > h3D_AzDispTable;
        vector< TH3F* >       h3D_AzDispTableN;
        vector< TProfile3D* > h3D_AzDispPhiTable;
        vector< TProfile3D* > h3D_AzDispMissTable;
        
        TFile* fTableFile;
        VDispTableReader* fData;
        
        // scaling parameter
        double fWidthScaleParameter;
        double fLengthScaleParameter;
        
        // quality cuts
        int    fQC_Ntubes_min;
        double fQC_Size_min;
        double fQC_Length_min;
        double fQC_Loss_max;
        
        bool prepareTraining( string );
        bool isGoodEvent( Ctpars* );                              // check quality cuts
        
    public:
    
        VDispTable();
        VDispTable( unsigned int iNtel, string iOutFile );
        ~VDispTable() {}
        void addAzBin( float iMin, float iMax );
        bool fillTable( string iMCfile, float i_ze, float i_woff, int iNentries = -1 );
        void setNoise( unsigned int k, float iNoise )
        {
            if( fData )
            {
                fData->setNoise( k, iNoise );
            }
        }
        void setDataVectors( vector< string > i_ze, vector< string > i_woff, vector< string > i_noise );
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        void setQualityCuts( int iNtubes_min = 4, double iSize_min = 0., double iLength_min = 0., double iLoss_max = 0.5 );
        void setWidthLengthScalingParameters( double iWS = 0.08, double iLS = 0.21 )
        {
            fWidthScaleParameter = iWS;
            fLengthScaleParameter = iLS;
        }
        bool terminate();
        
};
#endif
