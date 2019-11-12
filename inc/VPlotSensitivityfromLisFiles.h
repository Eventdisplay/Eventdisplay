//! VPlotSensitivityfromLisFiles plot sensitivity parameters from lis files

#ifndef VPlotSensitivityfromLisFiles_H
#define VPlotSensitivityfromLisFiles_H

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMath.h"

#include "VPlotUtilities.h"

using namespace std;

class VLisFileData
{
    public:
    
        unsigned int fID;
        string  fName;
        string  fFileName;
        
        // data from list file
        map< string, vector< double > > fVar;
        
        VLisFileData();
        ~VLisFileData() {};
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

class VPlotSensitivityfromLisFiles : public VPlotUtilities
{
    private:
    
        bool fDebug;
        
        // Variables
        vector< string > fVarName;
        
        // min/max for plotting
        map< string, double > fVarMin;
        map< string, double > fVarMax;
        
        // data
        vector< VLisFileData* > fData;
        
    public:
    
        VPlotSensitivityfromLisFiles();
        ~VPlotSensitivityfromLisFiles() {}
        
        bool addLisFile( string iFile, string iCut = "" );
        bool applycuts( double amp = 1., double NTel = 1., double NPix = 1. );
        bool checkVarName( string V );
        void compare_cuts();
        unsigned int getID_Index( unsigned int iID );
        void listDataSets();
        void listVariableNames();
        bool printDataSet( unsigned int iID );
        bool removeDataSet( unsigned int iDataSetID );
        void setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        
        TCanvas* plot( string iVName, unsigned int iDataSetID = 0, TCanvas* c = 0, Style_t iLineStyle = 1, Style_t iMarkerStyle = 20 );
        TCanvas* plot_AllIDs( string iVName, Style_t iLineStyle = 1, Style_t iMarkerStyle = 20, TCanvas* c = 0 );
};

#endif
