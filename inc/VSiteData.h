//! VSiteData  data class for sites

#ifndef VSiteData_H
#define VSiteData_H

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VSiteData
{
    private:
    
        bool             fDebug;
        
    public:
    
        string           fSiteName;
        string           fHemisphere;
        float            fSite_asl;
        float            fSite_B_N;
        float            fSite_B_S;
        float            fSite_B_dB;
        float            fSiteRequirementID;
        string           fReferenceSiteName;
        string           fSiteFileType;
        vector< string > fSiteFileName;
        vector< bool   > fSiteFile_exists;
        vector< float >  fSiteFile_Emin;
        vector< float >  fSiteFile_Emax;
        vector< float >  fObservationTime_s;
        vector< int >    fTelCutMSTs;
        vector< int >    fTelCutSSTs;
        vector< float >  fCameraOffset_deg;
        vector< string > fArray;
        vector< TGraphAsymmErrors* > fGraphSensitivity;
        TGraphAsymmErrors* fGraphSensitivityInterPolated;
        
        vector< int >    fPlottingColor;
        vector< int >    fPlottingLineStyle;
        vector< int >    fPlottingFillStyle;
        vector< int >    fPlottingMarkerStyle;
        vector< string > fLegend;
        
        
        VSiteData();
        ~VSiteData();
        bool  addDataSet( string iDataList, unsigned int iSiteCounter, string iDirectionString = "" );
        bool  checkIntegrity();
        TGraphAsymmErrors* getCombinedSensitivityGraph( bool iInterpolate = false, string iDirectionString = "", bool iIntegratedSensitivity = false );
        void  print();
        void  reset();
        void  setDebug( bool iB = false )
        {
            fDebug = iB;
        }
};


#endif
