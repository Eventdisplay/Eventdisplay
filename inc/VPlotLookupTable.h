//! VPlotLookupTable plot and compare lookup tables

#ifndef VPlotLookupTable_H
#define VPlotLookupTable_H

#include <iostream>
#include <string>
#include <vector>

#include "VPlotUtilities.h"
#include "VInterpolate2DHistos.h"
#include "VHistogramUtilities.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"

using namespace std;

class VPlotLookupTableData
{
    public:
    
        string fLookupTable;
        string fLookupTableFileName;
        TFile* fLookupTableFile;
        int    fZe;
        int    fAz;
        int    fTelID;
        int    fNoise;
        int    fWobbleOffset;
        
        TH2F*  hmedian;
        TH2F*  hmean;
        TH2F*  hmpv;
        TH2F*  hsigma;
        TH2F*  hnevents;
        
        VPlotLookupTableData();
        ~VPlotLookupTableData() {}
};

class VPlotLookupTable : public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        vector< string > fListOfTableNames;
        
        vector< VPlotLookupTableData* > fLookupTableData;
        
        double  fLogSizeAxis_min;
        double  fLogSizeAxis_max;
        double  fLogEnergyAxis_min;
        double  fLogEnergyAxis_max;
        double  fDistanceAxis_min;
        double  fDistanceAxis_max;
        
        bool    checkTableName( string iTableName );
        TH2F*   divide2DHistograms( TH2F* h1, TH2F* h2, char* hname );
        void    plot2DHistogram( TH2F* h, unsigned int iSetID, string htitle, int iCanvasX,
                                 double i_min = -999., double i_max = -999., bool iLogZ = false );
                                 
                                 
    public:
    
        VPlotLookupTable();
        ~VPlotLookupTable() {}
        
        bool  addLookupTable( string iLookupTableFile, string iTable = "mscw", int ze = 20, int az = 0,
                              int telID = 1, int noise = 455, int woff = 500 );
                              
        void  printLookupTables();
        void  plotLookupTables( unsigned int iSetID = 0, double i_ymin = 1.e-3, bool iMedianOnly = false );
        void  plotRelativeTables( unsigned int iSetID1, unsigned int iSetID2, double iMin = 0.95, double iMax = 1.05 );
        void  plotLookupTableSlice( double iLogSize = 5., double iR = 0. );
        void  setPlottingDistanceAxis( double imin = 0., double imax = 800. )
        {
            fDistanceAxis_min = imin;
            fDistanceAxis_max = imax;
        }
        void  setPlottingLogEnergyAxis( double imin = -2., double imax = 2.5 )
        {
            fLogEnergyAxis_min = imin;
            fLogEnergyAxis_max = imax;
        }
        void  setPlottingLogSizeAxis( double imin = 1.0, double imax = 7. )
        {
            fLogSizeAxis_min = imin;
            fLogSizeAxis_max = imax;
        }
        bool smoothLookupTables( unsigned int iSetID = 0, string iMethod = "interpolate", int iMinEvents = 0 );
        
};

#endif
