//! VExposure calculate VERITAS exposure

#ifndef VEXPOSURE_H
#define VEXPOSURE_H

#include "VAstronometry.h"
#include "VDBTools.h"
#include "VStarCatalogue.h"
#include "VGlobalRunParameter.h"
#include "VUtilities.h"
#include "VDB_Connection.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <bitset>

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TEventList.h>
#include <TF1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TMath.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>
#include <TString.h>

using namespace std;

class VExposure : public TObject, public VGlobalRunParameter
{
    private:
    
        bool   fDebug;
        
        bool bPlotElevationPlots;
        
        bool fMakeRunList;
        int fSelectLaser;
        int fDataStartTime; // Start Date
        vector< unsigned int > fLaserRunID; // Laser Runs
        bool bPrintVerbose;
        bool bPrintTimeMask;
        
        // Type of observing mode
        string fObservingMode;
        
        // Minimium Duration
        double fMinDuration;
        
        // Telscope Min Elevation
        double fTelMinElevation;
        
        // Target Name
        string fTargetSourceName;
        
        // date range
        string fStartDate_SQL;
        string fStopDate_SQL;
        
        // acceptance curves
        TF1* fAcceptance;
        double fAcceptance_MaxDistance;
        double fMaximumIntegrationRadius;
        
        bool fPlotExtendedSources;
        bool fPlotSourceNames;
        
        bool fDoCheckSums;
        vector<int> fRunsNoChecksum;
        vector<int> fRunsGoodChecksum;
        vector<int> fRunsBadChecksum;
        
        // list of sources from DB
        VDB_ObservingSources* fVDB_ObservingSources;
        
        vector< int > fRun;
        vector< unsigned int> fRunConfigMask;
        vector< string > fRunStatus;
        vector< double > fRunStartMJD;
        vector< double > fRunStopMJD;
        vector< double > fRunDuration;
        vector< string > fRunSourceID;
        vector< string > fRunObsMode;
        vector< double > fRunRA;
        vector< double > fRunDec;
        vector< double > fRunoffsetRA;
        vector< double > fRunoffsetDec;
        vector< double > fRunGalLong1958;
        vector< double > fRunGalLat1958;
        vector< double > fRunTelElevation;
        vector< double > fRunTelAzimuth;
        vector< double > fWobbleN;
        vector< double > fWobbleE;
        vector< int >    fRunDate;
        vector< vector < unsigned int > > fRunLaserList;
        vector< unsigned int > fDateLaserList;
        
        vector< int > fRunDownload;
        vector< int > fRunDownloadDate;
        vector< unsigned int > fLaserDownload;
        vector< int > fLaserDownloadDate;
        vector< int > fRunDownloadList;
        vector< string > fDataCat;
        vector< string > fStatus;
        vector< string > fStatReason;
        vector< string > fTelCutMask;
        vector< string > fUsable;
        vector< string > fTimeCutMask;
        vector< string > fLightLevel;
        vector< string > fVPMcon;
        vector< string > fAuthor;
        vector< string > fComment;
        
        TH2D* fMapGal2D;
        TH2D* fRadAccMapGal2D;
        TH2D* fMapGal2D_aitoff;
        TH2D* fRadAccMapGal2D_aitoff;
        
        TH1D* fTimeDifferencesBetweenRuns;
        
        vector< TCanvas* > fPlottingCanvas;
        int  fCanvasSize_x;
        int  fCanvasSize_y;
        
        bool   fPlotVTSObjects;
        vector< string > fCatalogue;
        vector< int >    fCatalogueMarkerColor;
        vector< int >    fCatalogueMarkerStyle;
        vector< double > fCatalogueTextAngle;
        
        vector< string > fTexTable;
        double fTexTable_EFlux_min;
        double fTexTable_EFlux_max;
        
        void   aitoff2xy( Double_t l, Double_t b, Double_t& Al, Double_t& Ab );
        bool   doDQM( unsigned int iIndex, double iMinDuration = 600. );
        void   drawAitoffCoordinateSystem();
        double getAcceptance( double r );
        void   getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip );
        bool   getDBSourceCoordinates( TSQLServer* f_db, string iSource, double& iEVNTargetDec, double& iEVNTargetRA );
        void   analyseCatalogue( string iCatalogue = "../../eventdisplay/astro/tevcat.dat",
                                 double ibmin = -90., double ibmax = 90., double ilmin = -180., double ilmax = 180.,
                                 TH2* h = 0, bool bAitoff = false, int iMarkerStyle = 5, int iMarkerColor = 1, double iTextAngle = 45. );
        void   plotObject( double l, double b, string iname, double iextension,
                           double ibmin = -90., double ibmax = 90., double ilmin = -180., double ilmax = 180.,
                           TH2* h = 0, bool bAitoff = false, int iMarkerStyle = 5, int iMarkerColor = 1, double iTextAngle = 45. );
        TCanvas* plot2DGalactic( string iName, string iTitle, int ix, int iy, int iwx, int iwy, TH2* h,
                                 double ibmin = -90., double ibmax = 90., double ilmin = -180., double ilmax = 180.,
                                 bool bAitoff = false, int iPalette = 1, TString opt = "A colz" );
        void plotVTSObjects( bool bAitoff = false, double ibmin = -90., double ibmax = 90., double ilmin = -180., double ilmax = 180.,
                             int iMarkerStyle = 5, int iMarkerColor = 1, double iTextAngle = 45., TH2* h = 0 );
        void set_plot_style();
        void resetDataVectors();
    public:
    
        VExposure( int nBinsL = 5000, int nBinB = 2000 );
        
        bool readFromDB();
        bool readFromDBList();
        bool setPlannedObservation( vector<double> ra, vector<double> dec, vector<double> t );
        void fillElevationPlot( int iYear = 2014, int iMonth = -1, int ze_max_deg = 40 );
        void fillExposureMap();
        TH2D* getGalacticMap()
        {
            return fMapGal2D;
        }
        TH2D* getGalacticMapAcceptanceCorrected()
        {
            return fRadAccMapGal2D;
        }
        TH2D* getGalacticMapAitoff()
        {
            return fMapGal2D_aitoff;
        }
        TH2D* getGalacticMapAcceptanceCorrectedAitoff()
        {
            return fRadAccMapGal2D_aitoff;
        }
        void  plotMarker( double l, double b, double r = 0., string iText = "", int iMarkerStyle = 5, int iMarkerColor = 1, int iMarkerSize = 1, double iTextAngle = 45. );
        TCanvas* plot( double ilmin = -180., double ilmax = 180., double ibmin = -90., double ibmax = 90., unsigned int iReturnCanvas = 0, int iPalette = 55, TString opt = "A colz" );
        void     plot_HESSSkySurvey( TCanvas* c );
        void plotTimeDifferencesbetweenRuns();
        void printListOfRuns( string iCatalogue, double iR = 2.5, double iMinDuration = 600., string iTeVCatalogue = "",
                              double r_min = 0.1, string iEventListFile = "" );
        void printListOfRuns( double il, double ib, double iR = 2.5, double iMinDuration = 600., string iDQMfileList = "",
                              string ofile = "", unsigned int iVerbose = 0 );
        void printListOfRuns();
        void printShortRunList();
        void outputAnasumRunlist( string fAnasumFile );
        void printTexTable();
        bool readAcceptanceCurveFromFile( string iAcc, double iAcceptance_MaxDistance = 1.e9 );
        bool readRootFile( string iname = "2006_2008.root", double iMinMJD = -99., double iMaxMJD = -99. );
        void setPlotExtendedSourcesFromCataloge( bool iB = false )
        {
            fPlotExtendedSources = iB;
        }
        void setPlotSourceNames( bool iB = false )
        {
            fPlotSourceNames = iB;
        }
        void setDoCheckSums( bool iB = true )
        {
            fDoCheckSums = iB;
        }
        void setMaximumIntegrationRadius( double iR = 1.5 )
        {
            fMaximumIntegrationRadius = iR;    // FOV (VERITAS is 3.5 deg, this is the (optimistic) region of constant radial acceptance)
        }
        void setTimeRange( string iStart = "2011-01-01", string iStopp = "2011-02-01" );
        void setSourceName( string sourceName = "Crab" );
        void setTelMinElevation( double iElevation = 0. );
        void setSelectLaser( int iSelectLaser );
        void setMinDuration( double iDuration = 0. );
        void setPrintTimeMask( int iPrintTimeMask );
        void setPrintVerbose( int iPrintVerbose );
        bool writeRootFile( string );
        void setMakeRunList( bool iSet );
        void getLaserList();
        unsigned int  getLaserDate( unsigned int iRunNumber );
        void readRunListFromFile( string runlist );
        void readLaserRunDateListFromFile( string runlist );
        void readLaserRunListFromFile( string runlist );
        void downloadRunList();
        void checkRunList();
        void readRunCommentsFromDB();
        void setRunNumber( unsigned int number );
        void setLaserNumber( unsigned int number );
        void setObservingMode( bool bObs );
        void setCanvasSize( double x = 600, double y = 400 );
        void setPlotVTSObjects( bool iVTS = true )
        {
            fPlotVTSObjects = iVTS;
        }
        
        vector< unsigned int > getLaserRun( unsigned int iRunNumber, unsigned int iNTel );
        //	TSQLServer* connectToSQLServer( string iServer );
        
        void addCatalogue( string, int iMarker = 5, int iColor = 50, double iAngle = 45. );
        void listCatalogues();
        bool removeCataloge( unsigned int iB );
        
        TString getArchiveMD5sum( int date, int run, bool force_download = false );
        TString calcMD5sum( int date, int run );
        TString readMD5sumFromFile( TString filename, int run, bool warn = true );
        int checkMD5sum( int date, int run, bool force_download = false ) ;
        void printChecksumSummary();
        
        ClassDef( VExposure, 10 );
};
#endif
