//! VDispTableReader reading disp tables for angular reconstruction

#ifndef VDispTableReader_H
#define VDispTableReader_H

#include "TCanvas.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TTree.h"

#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

class VDispTableReader : public TNamed
{
    private:
    
        TList* hisList;
        TTree* fData;
        
        void         getIndexBoundary( unsigned int* ib, unsigned int* il, vector< float >& iV, float x );
        unsigned int getTreeEntryID( int i_found_ze, int i_found_az, int i_found_woff, int i_found_noise );
        bool         isHistoBinningSet( string iVarName );
        
    public:
    
        map<unsigned int, int> fTreeEntry_Map;
        vector< float > f_ze;
        vector< float > f_az_min;
        vector< float > f_az_max;
        vector< float > f_woff;
        vector< float > f_noise;
        
        float ze;
        unsigned int az_bin;
        float az_min;
        float az_max;
        float woff;
        float pedvar;
        // 2D histograms
        TH2F* h2D_DispTable;
        TH2F* h2D_DispPhiTable;
        TH2F* h2D_DispTableN;
        TH2F* h2D_DispMissTable;
        // 3D histograms
        TH3F* h3D_DispTable;
        TH3F* h3D_DispPhiTable;
        TH3F* h3D_DispTableN;
        TH3F* h3D_DispMissTable;
        // histogram binning
        vector< string > fHisto_ListOfVariables;
        map< string, int >    fHisto_binning;
        map< string, double > fHisto_min;
        map< string, double > fHisto_max;
        
        
        // disp table scaling
        double fWidthScaleParameter;
        double fLengthScaleParameter;
        
        VDispTableReader();
        ~VDispTableReader() {}
        
        bool fill( float i_ze, unsigned int i_az, float i_az_min, float i_az_max, float i_woff, float i_meanPedvars, TH2* iH2D, TH2* iH2DN, TH2* h2DPhi, TH2* h2DMiss, TH3* iH3D, TH3* iH3DN, TH3* h3DPhi, TH3* h3DMiss );
        unsigned int getAzBin( float az );
        float getLowerZe( float );
        float getUpperZe( float );
        TTree* getTree()
        {
            return fData;
        }
        int  getTreeEntryFinder( unsigned int );
        int  getTreeEntryFinder( float iZe, float iAz, float iWoff, float iPedvar, int iZe_Inter );
        bool initialize( bool iRead = true );
        bool plot( int iEntry = 0 );
        void print( bool bDetailed = false );
        void reset();
        double scaleParameter( double iPara, double iScale );
        double scaleLengthParameter( double iLength )
        {
            return scaleParameter( iLength, fLengthScaleParameter );
        }
        double scaleWidthParameter( double iWidth )
        {
            return scaleParameter( iWidth, fWidthScaleParameter );
        }
        void setDataVectors( vector< string > i_ze, vector< float > i_az_min, vector< float > i_az_max, vector< string > i_woff, vector< string > i_noise );
        bool setHistoBinning( string iVariableName, int iBin, double iMin, double iMax );
        void setNoise( unsigned int k, float iNoise )
        {
            if( k < f_noise.size() )
            {
                f_noise[k] = iNoise;
            }
        }
        void setWidthLengthScalingParameters( double iWS = 0.08, double iLS = 0.21 )
        {
            fWidthScaleParameter = iWS;
            fLengthScaleParameter = iLS;
        }
        void terminate();
        
        ClassDef( VDispTableReader, 6 );
};
#endif
