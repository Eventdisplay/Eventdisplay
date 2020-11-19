//! VEvndispData  central event data class

#ifndef VDATA_H
#define VDATA_H

#include "VImageAnalyzerData.h"
#include "VEvndispReconstructionParameter.h"
#include "VCalibrationData.h"
#include "VDeadChannelFinder.h"
#include "VDetectorGeometry.h"
#include "VDetectorTree.h"
#include "VDSTReader.h"
#include "VGrIsuReader.h"
#include "VMCParameters.h"
#include "VMultipleGrIsuReader.h"
#ifndef NOVBF
#include "VBaseRawDataReader.h"
#endif
#include "VPEReader.h"
#include "VDB_PixelDataReader.h"
#include "VEvndispRunParameter.h"
#include "VStarCatalogue.h"
#include "VShowerParameters.h"
#include "VFrogsParameters.h"
//#include "VFrogsImageData.h"
#include "VPointing.h"
#include "VArrayPointing.h"
#include "VTraceHandler.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"

#include <bitset>
#include <iostream>
#include <string>
#include <vector>


using namespace std;

class VEvndispData
{
    private:
    
    protected:
        static bool fDebug;
        static int  fDebugLevel;
        static int  fNDebugMessages;
        
        // global run parameters
        static int fRunNumber;                    //!< run number
        static VEvndispRunParameter* fRunPar;          //!< all command line/configuration parameters
        
        // telescope data
        static unsigned int fNTel;                //!< total number of telescopes
        static unsigned int fTelID;               //!< telescope number of current telescope
        static vector< unsigned int > fTeltoAna;  //!< analyze only this subset of telescopes (this is dynamic and can change from event to event)
        // telescope pointing (one per telescope)
        static VArrayPointing* fArrayPointing;
        static vector< VPointing* > fPointing;
        // no pointing information
        static bool fNoTelescopePointing;
        // cameras
        static VDetectorGeometry* fDetectorGeo;
        static VDetectorTree* fDetectorTree;
        
        // reader
        VVirtualDataReader*  fReader;
        static VGrIsuReader* fGrIsuReader;
        static VMultipleGrIsuReader* fMultipleGrIsuReader;
#ifndef NOVBF
        static VBaseRawDataReader* fRawDataReader;
#endif
        static VDSTReader* fDSTReader;
        static VPEReader*  fPEReader;
        
        // DB pixel data
        static VDB_PixelDataReader* fDB_PixelDataReader;
        
        // event data
        static unsigned int fEventNumber;         //!< current event number (array event)
        //!< event number of telescope event
        static vector< unsigned int > fTelescopeEventNumber;
        static unsigned int fEventType;           //!< current event type
        static int fArrayEventMJD;                //!< MJD of current event
        static int fArrayPreviousEventMJD;        //!< MJD of previous event
        static double fArrayEventTime;            //!< time of current event
        static vector< int > fEventMJD;           //!< MJD of current event (per telescope)
        static vector< double > fEventTime;       //!< time of current event (per telescope)
        
        static vector< vector< int > > fTriggeredTel;
        static vector< int > fTriggeredTelN;
        
        // event status from data reader
        static unsigned long int fExpectedEventStatus;
        static unsigned int fNumberofGoodEvents;
        static unsigned int fNumberofIncompleteEvents;
        // event status from analysis
        //!< 0: good event
        static unsigned int fAnalysisArrayEventStatus;
        //!< 0: good event
        static vector< unsigned int > fAnalysisTelescopeEventStatus;
        
        // global trace handler
        static VTraceHandler* fTraceHandler;
        
        // calibrator and calibration data
        static vector< bool > fCalibrated;        //!< this telescope is calibrated
        //! data class for calibration data [ntel]
        static vector< VCalibrationData* > fCalData;
        // dead channel finder
        static vector< VDeadChannelFinder* > fDeadChannelDefinition_HG;
        static vector< VDeadChannelFinder* > fDeadChannelDefinition_LG;
        
        //  analysis cuts
        //!< cuts for array analysis
        static VEvndispReconstructionParameter* fEvndispReconstructionParameter;
        
        // analysis results
        static TFile* fOutputfile;                //!< root output file for image parameter trees, histograms, etc.
        static vector< TDirectory* > fAnaDir;     //! directories in root file
        static vector< VImageAnalyzerData* > fAnaData; //!< data class with analysis results for each telescope
        //!< data class with analysis results from all telescopes
        static VShowerParameters* fShowerParameters;
        static VFrogsParameters* fFrogsParameters;
        //	static vector< VFrogImageData* > fFrogsData;    //!< frogs Template tube information
        static VMCParameters* fMCParameters;      //!< data class with MC parameters
        
        // timing results
        static vector< TGraphErrors* > fXGraph;   //!< Long axis timing graph
        static vector< TGraphErrors* > fXGraphDP2;   //!< Long axis timing graph (double pass, pass 2)
        
        // default pedestals for plotraw option
        static valarray<double> fPlotRawPedestals;
        
        // temp
        static valarray<double> fTempValArray;
        
        // set detector geometry
        unsigned int        checkSummationWindow( unsigned int iTelID, unsigned int iSumWindow );
        bool                setDetectorGeometry( unsigned int iNTel, vector< string > icamera , string idir );
        
        // names of dead channels
        static vector< string > fDeadChannelText;
        
        // star catalogue
        static VStarCatalogue* fStarCatalogue;
        
        // dummy vector
        static vector< float > fDummyVector_float;
        
    public:
        VEvndispData();
        ~VEvndispData() {}
        void                dumpTreeData();       //!< print all tree data to stdout
        void                endOfRunInfo();       //!< print some statistics at end of run
        bool                get_reconstruction_parameters( string ifile, bool iMakeNeighbourList );
        
        // getters apply always to current telescope (fTelID) if telID is not a function argument
        VImageAnalyzerData*      getAnaData( unsigned int iTel )
        {
            if( iTel < fAnaData.size() )
            {
                return fAnaData[iTel];
            }
            else
            {
                return 0;
            }
        }
        VImageAnalyzerData*      getAnaData()
        {
            return fAnaData[fTelID];
        }
        vector< TDirectory* > getAnaDirectories()
        {
            return fAnaDir;
        }
        VImageAnalyzerHistograms*     getAnaHistos()
        {
            return fAnaData[fTelID]->fAnaHistos;
        }
        VImageAnalyzerHistograms*     getAnaHistos( unsigned int itelID )
        {
            return fAnaData[itelID]->fAnaHistos;
        }
        VEvndispReconstructionParameter* getEvndispReconstructionParameter()
        {
            return fEvndispReconstructionParameter;
        }
        VEvndispReconstructionParameterData* getEvndispReconstructionParameter( unsigned int iMethod )
        {
            if( fEvndispReconstructionParameter )
            {
                return fEvndispReconstructionParameter->getReconstructionParameterData( iMethod );
            }
            return 0;
        }
        unsigned int        getAnalysisArrayEventStatus()
        {
            return fAnalysisArrayEventStatus;
        }
        vector< unsigned int >& getAnalysisTelescopeEventStatus()
        {
            return fAnalysisTelescopeEventStatus;
        }
        double getAverageElevation();
        vector<bool>&       getBorder()
        {
            return fAnaData[fTelID]->fBorder;
        }
        VCalibrationData*   getCalData()
        {
            return fCalData[fTelID];
        }
        VCalibrationData*   getCalData( unsigned int iTel )
        {
            if( iTel < fCalData.size() )
            {
                return fCalData[iTel];
            }
            else
            {
                return 0;
            }
        }
        VCalibrationData*   getCalibrationData()
        {
            return fCalData[fTelID];
        }
        VCalibrationData*   getCalibrationData( unsigned int iTel )
        {
            if( iTel < fCalData.size() )
            {
                return fCalData[iTel];
            }
            else
            {
                return 0;
            }
        }
        vector<int>&        getChannelStatus()
        {
            return fCalData[fTelID]->fChannelStatus;
        }
        bool                getCalibrated()
        {
            return fCalibrated[fTelID];
        }
        vector<double>&     getBorderCorrelationCoefficient()
        {
            return fAnaData[fTelID]->fCorrelationCoefficient;
        }
        VEvndispData*           getData()
        {
            return this;
        }
        unsigned int        getDead( unsigned int iChannel, bool iLowGain );
        vector<unsigned int>&   getDead( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fAnaData[fTelID]->fDead;
            }
            else
            {
                return fAnaData[fTelID]->fLowGainDead;
            }
        }
        vector<unsigned int>&   getMasked()
        {
            return fAnaData[fTelID]->fMasked;
        }
        vector<bool>&       getDeadRecovered( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fAnaData[fTelID]->fDeadRecovered;
            }
            else
            {
                return fAnaData[fTelID]->fLowGainDeadRecovered;
            }
        }
        VDeadChannelFinder* getDeadChannelFinder( bool iLowGain = false );
        unsigned int        getNDead( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fAnaData[fTelID]->fNDead;
            }
            else
            {
                return fAnaData[fTelID]->fLowGainNDead;
            }
        }
        vector<unsigned int>&       getDeadUI( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fAnaData[fTelID]->fDeadUI;
            }
            else
            {
                return fAnaData[fTelID]->fLowGainDeadUI;
            }
        }
        vector<string>&     getDeadChannelText()
        {
            return fDeadChannelText;
        }
        bool                getDebugFlag()
        {
            return fDebug;
        }
        int                 getDebugLevel();
        VDetectorGeometry*  getDetectorGeo() const
        {
            return fDetectorGeo;
        }
        VDetectorGeometry*  getDetectorGeometry() const
        {
            return fDetectorGeo;
        }
        TTree*              getDetectorTree();
        bool                getDoublePassErrorWeighting2005()
        {
            if( fTelID < fRunPar->fDoublePassErrorWeighting2005.size() )
            {
                return fRunPar->fDoublePassErrorWeighting2005[fTelID];
            }
            return false;
        }
        int                 getEventMJD()
        {
            return fArrayEventMJD;
        }
        vector< int >&      getEventMJDVector()
        {
            return fEventMJD;
        }
        unsigned int        getEventNumber()
        {
            return fEventNumber;
        }
        string              getEventDisplayVersion()
        {
            return getRunParameter()->getEVNDISP_VERSION();
        }
        double              getEventTime()
        {
            return fArrayEventTime;
        }
        vector< double >&   getEventTimeVector()
        {
            return fEventTime;
        }
        unsigned int        getEventType()
        {
            return fEventType;
        }
        unsigned long int   getExpectedEventStatus()
        {
            return fExpectedEventStatus;
        }
        valarray<double>&   getFADCStopOffsets()
        {
            return fCalData[fTelID]->fFADCStopOffsets;
        }
        vector< double >&   getFADCstopSums()
        {
            return fAnaData[fTelID]->fFADCstopSum;
        }
        vector< unsigned int >&  getFADCstopTrig()
        {
            return fAnaData[fTelID]->getFADCstopTrigChannelID();
        }
        vector< double >&   getFADCstopTZero()
        {
            return fAnaData[fTelID]->fFADCstopTZero;
        }
        double getTelescopeAverageFADCtoPhe( bool iLowGain = false )
        {
            return fCalData[fTelID]->getTelescopeAverageFADCtoPhe( iLowGain );
        }
        valarray<double>& getFADCtoPhe( bool iLowGain = false )
        {
            if( iLowGain )
            {
                return fCalData[fTelID]->fLowGainFADCtoPhe;
            }
            return fCalData[fTelID]->fFADCtoPhe;
        }
        bool                getFillMeanTraces()
        {
            return fAnaData[fTelID]->fFillMeanTraces;
        }
        bool                getFillPulseSum()
        {
            return fAnaData[fTelID]->fFillPulseSum;
        }
        TH1F*               getGainDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getGainDist( iLowGain );
        }
        TH1F*               getGainVarsDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getGainVarsDist( iLowGain );
        }
        valarray<double>&   getGains( bool iLowGain = false )
        {
            if( iLowGain )
            {
                return fCalData[fTelID]->fLowGainGains;
            }
            return fCalData[fTelID]->fGains;
        }
        valarray< bool >&   getGains_DefaultValue( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fGains_DefaultSetting;
            }
            else
            {
                return fCalData[fTelID]->fLowGainGains_DefaultSetting;
            }
        }
        valarray<double>&   getGainvars( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fGainvars;
            }
            else
            {
                return fCalData[fTelID]->fLowGainGainvars;
            }
        }
        double              getHIGHQE_gainfactor( unsigned int iChannel )
        {
            return fAnaData[fTelID]->getHIGHQE_gainfactor( iChannel );
        }
        vector<bool>&       getHiLo()
        {
            return fAnaData[fTelID]->fHiLo;
        }
        VImageAnalyzerHistograms*     getHistograms()
        {
            return fAnaData[fTelID]->fAnaHistos;
        }
        vector<bool>&       getImage()
        {
            return fAnaData[fTelID]->fImage;
        }
        VImageCleaningRunParameter* getImageCleaningParameter( bool iDoublePassParameters = false, unsigned int iTelID = 99999 );
        vector<bool>&       getImageBorderNeighbour()
        {
            return fAnaData[fTelID]->fImageBorderNeighbour;
        }
        vector<bool>&       getBorderBorderNeighbour()
        {
            return fAnaData[fTelID]->fBorderBorderNeighbour;
        }
        VImageParameter*    getImageParameters()
        {
            return fAnaData[fTelID]->fImageParameter;
        }
        VImageParameter*    getImageParameters( int );
        VImageParameter*    getImageParametersLogL()
        {
            return fAnaData[fTelID]->fImageParameterLogL;
        }
        vector<int>&        getImageUser()
        {
            return fAnaData[fTelID]->fImageUser;
        }
        vector<bool>&       getLLEst()
        {
            return fAnaData[fTelID]->fLLEst;
        }
        TList*              getIntegratedChargeHistograms()
        {
            return fAnaData[fTelID]->getIntegratedChargeHistograms();
        }
        TGraphErrors*  getIPRGraph()
        {
            return fCalData[fTelID]->getIPRGraph( getSumWindow(), false );
        }
        TGraphErrors*  getIPRGraph( unsigned int iSumWindow, bool iMakeNewGraph = false )
        {
            return fCalData[fTelID]->getIPRGraph( iSumWindow, iMakeNewGraph );
        }
        float               getL1Rate( unsigned int iChannel )
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getL1Rate( getTelID(), iChannel, getEventMJD(), getEventTime() );
            }
            else
            {
                return 0.;
            }
        }
        float               getHV( unsigned int iChannel )
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getHV( getTelID(), iChannel, getEventMJD(), getEventTime() );
            }
            else
            {
                return 0.;
            }
        }
        float               getCurrent( unsigned int iChannel )
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getCurrent( getTelID(), iChannel, getEventMJD(), getEventTime() );
            }
            else
            {
                return 0.;
            }
        }
        vector< float >     getL1Rates()
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getL1Rates( getTelID(), getEventMJD(), getEventTime() );
            }
            else
            {
                return fDummyVector_float;
            }
        }
        vector< float >     getHV()
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getHV( getTelID(), getEventMJD(), getEventTime() );
            }
            else
            {
                return fDummyVector_float;
            }
        }
        vector< float >     getCurrents()
        {
            if( fDB_PixelDataReader )
            {
                return fDB_PixelDataReader->getCurrents( getTelID(), getEventMJD(), getEventTime() );
            }
            else
            {
                return fDummyVector_float;
            }
        }
        valarray<double>&   getLowGainMultiplier_Camera()
        {
            return fCalData[fTelID]->getLowGainMultiplier_Camera() ;
        }
        vector< pair<unsigned int, int> >&        getLowGainDefaultSumWindows()
        {
            return fCalData[fTelID]->fLowGainDefaultSumWindows ;
        }
        double              getLowGainMultiplier_Trace( )
        {
            return fCalData[fTelID]->getLowGainMultiplier_Trace() ;
        }
        double              getLowGainMultiplier_Sum( unsigned int iMethod, int iWindow, int jWindow )
        {
            return fCalData[fTelID]->getLowGainMultiplier_Sum( iMethod, iWindow, jWindow );
        }
        double              getLowGainSumCorrection( unsigned int iMethod, int iSumWindow, int jSumWindow , bool HiLo = true )
        {
            return fCalData[fTelID]->getLowGainSumCorrection( iMethod, iSumWindow, jSumWindow , HiLo ) ;
        }
        bool                getLowGainPedestals()
        {
            return fCalData[fTelID]->fBoolLowGainPedestals;
        }
        bool                getLowGainGains()
        {
            return fCalData[fTelID]->fBoolLowGainGains;
        }
        bool                getLowGainTOff()
        {
            return fCalData[fTelID]->fBoolLowGainTOff;
        }
        double              getSumWindowMaxTimeDifferenceLGtoHG()
        {
            if( fTelID < fRunPar->fSumWindowMaxTimeDifferenceLGtoHG.size() )
            {
                return fRunPar->fSumWindowMaxTimeDifferenceLGtoHG[fTelID];
            }
            return -999.;
        }
        VMCParameters*      getMCParameters()
        {
            return fMCParameters;
        }
        unsigned int        getMC_FADC_TraceStart()
        {
            return fRunPar->fMC_FADCTraceStart;
        }
        
        double              getmeanPedvars( bool iLowGain = false, unsigned int iSumWindow = 0 )
        {
            return getCalData()->getmeanPedvars( iLowGain, iSumWindow );
        }
        double              getmeanRMSPedvars( bool iLowGain, unsigned int iSumWindow )
        {
            return getCalData()->getmeanRMSPedvars( iLowGain, iSumWindow );
        }
        void                getmeanPedvars( double& imean, double& irms, bool iLowGain = false, unsigned int iSumWindow = 0,  double iTime = -99. )
        {
            if( iTime < -90. )
            {
                getCalData()->getmeanPedvars( imean, irms, iLowGain, iSumWindow, getEventTime() );
            }
            else
            {
                getCalData()->getmeanPedvars( imean, irms, iLowGain, iSumWindow, getEventTime() );
            }
        }
        
        vector< double >&   getmeanPedvarsAllSumWindow( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fVmeanPedvars;
            }
            else
            {
                return fCalData[fTelID]->fVmeanLowGainPedvars;
            }
        }
        vector< double >&   getmeanRMSPedvarsAllSumWindow( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fVmeanRMSPedvars;
            }
            else
            {
                return fCalData[fTelID]->fVmeanRMSLowGainPedvars;
            }
        }
        TList*              getMeanPulseHistograms()
        {
            return fAnaData[fTelID]->getMeanPulseHistograms();
        }
        unsigned int        getNChannels()
        {
            return fDetectorGeo->getNChannels( fTelID );
        }
        vector<bool>&       getBrightNonImage()
        {
            return fAnaData[fTelID]->fBrightNonImage;
        }
        bool                getNoPointing()
        {
            return fNoTelescopePointing;
        }
        unsigned int        getNSamples( unsigned int iTelID = 9999 )
        {
            if( iTelID == 9999 )
            {
                return fDetectorGeo->getNSamples( fTelID );
            }
            else
            {
                return fDetectorGeo->getNSamples( iTelID );
            }
        }
        unsigned int        getNSamplesAnalysis( unsigned int iTelID = 9999 );
        unsigned int        getNTel() const
        {
            return fNTel;
        }
        TFile*              getOutputFile()
        {
            return fOutputfile;
        }
        valarray<double>&   getPE()
        {
            if( getReader() )
            {
                return getReader()->getPE();
            }
            return fTempValArray;
        }
        bool                getPedsFromPLine()
        {
            return fCalData[fTelID]->fPedFromPLine;
        }
        double              getPed_min( bool iLowGain = false )
        {
            return fCalData[fTelID]->getPed_min( iLowGain );
        }
        double              getPed_max( bool iLowGain = false )
        {
            return fCalData[fTelID]->getPed_max( iLowGain );
        }
        TH1F*               getPedDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getPedDist( iLowGain );
        }
        TH1F*               getPedvarsDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getPedvarsDist( iLowGain );
        }
        TH1F*               getPedLowGainDist()
        {
            return fCalData[fTelID]->getPedDist( true );
        }
        TH1F*               getPedvarsLowGainDist()
        {
            return fCalData[fTelID]->getPedvarsDist( true );
        }
        
        ///////////////// pedestals /////////////////////////////////
        // getters for pedestals
        valarray<double>&   getPeds( bool iLowGain = false, double iTime = -99. );
        valarray<double>&   getPedsLowGain( double iTime = -99. )
        {
            return getPeds( true, iTime );
        }
        
        // getters for pedestal variation
        valarray<double>&   getPedvars( bool iLowGain = false, unsigned int iSW = 0, double iTime = -99. );
        valarray<double>&   getPedvars( unsigned int iSW, bool iLowGain = false )
        {
            return getPedvars( iLowGain, iSW );
        }
        valarray<double>&   getPedvarsLowGain()
        {
            return getPedvars( true );
        }
        
        vector< valarray<double> >& getPedvarsAllSumWindows( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fVPedvars;
            }
            else
            {
                return fCalData[fTelID]->fVLowGainPedvars;
            }
        }
        // getter for pedestal rms
        valarray<double>&   getPedrms( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fPedrms;
            }
            else
            {
                return fCalData[fTelID]->fLowGainPedsrms;
            }
        }
        valarray<double>&   getPedsLowGainrms()
        {
            return fCalData[fTelID]->fLowGainPedsrms;
        }
        VDB_PixelDataReader* getDBPixelDataReader()
        {
            return fDB_PixelDataReader;
        }
        // padding stuff (probably out of date)
        /////////////// end pedestals //////////////////////////////
        
        valarray<double>&   getRawTZeros()
        {
            return fAnaData[fTelID]->getTZeros( false );
        }
        VVirtualDataReader* getReader()
        {
            return fReader;
        }
        int                 getRunNumber()
        {
            return fRunNumber;
        }
        VEvndispRunParameter*    getRunParameter()
        {
            return fRunPar;
        }
        VShowerParameters*  getShowerParameters()
        {
            return fShowerParameters;
        }
        VFrogsParameters*    getFrogsParameters()
        {
            return fFrogsParameters;
        }
        int                 getSumFirst()
        {
            return fRunPar->fsumfirst[fTelID];
        }
        unsigned int getSearchWindowLast()
        {
            return fRunPar->fSearchWindowLast[fTelID];
        }
        valarray<double>&   getSums()
        {
            return fAnaData[fTelID]->fSums;
        }
        valarray<double>&   getSums2()
        {
            return fAnaData[fTelID]->fSums2;
        }
        valarray<double>&   getTemplateMu()
        {
            return fAnaData[fTelID]->fTemplateMu;
        }
        double              getTemplateMuMin()
        {
            return fAnaData[fTelID]->fTemplateMu.min();
        }
        double              getTemplateMuMax()
        {
            return fAnaData[fTelID]->fTemplateMu.max();
        }
        unsigned int        getLargestSumWindow();
        unsigned int        getLargestSumWindow( unsigned int iTelID );
        VStarCatalogue*     getStarCatalogue()
        {
            return fStarCatalogue;
        }
        unsigned int        getSumWindow()
        {
            return checkSummationWindow( fTelID, fRunPar->fsumwindow_1[fTelID] );
        }
        unsigned int        getSumWindow_2()
        {
            return checkSummationWindow( fTelID, fRunPar->fsumwindow_2[fTelID] );
        }
        unsigned int        getSumWindow_Pass1()
        {
            return checkSummationWindow( fTelID, fRunPar->fsumwindow_pass1[fTelID] );
        }
        unsigned int        getSumWindow( unsigned int iTelID )
        {
            if( iTelID < fRunPar->fsumwindow_1.size() )
            {
                return checkSummationWindow( iTelID, fRunPar->fsumwindow_1[iTelID] );
            }
            else
            {
                return 0;
            }
        }
        unsigned int        getSumWindow_2( unsigned int iTelID )
        {
            if( iTelID < fRunPar->fsumwindow_2.size() )
            {
                return checkSummationWindow( iTelID, fRunPar->fsumwindow_2[iTelID] );
            }
            else
            {
                return 0;
            }
        }
        unsigned int        getSumWindow_Pass1( unsigned int iTelID )
        {
            if( iTelID < fRunPar->fsumwindow_pass1.size() )
            {
                return checkSummationWindow( iTelID, fRunPar->fsumwindow_pass1[iTelID] );
            }
            else
            {
                return 0;
            }
        }
        valarray< unsigned int >& getCurrentSumWindow()
        {
            return fAnaData[fTelID]->fCurrentSummationWindow;
        }
        valarray< unsigned int >& getCurrentSumWindow_2()
        {
            return fAnaData[fTelID]->fCurrentSummationWindow_2;
        }
        float                    getSumWindowShift()
        {
            return fRunPar->fTraceWindowShift[fTelID];
        }
        float         getSumwWindowStart_T_maxT0Diff()
        {
            return fRunPar->fsumfirst_maxT0startDiff[fTelID];
        }
        bool getSumWindow_searchmaxreverse()
        {
            return fRunPar->fSumWindow_searchmaxreverse[fTelID];
        }
        unsigned int  getSumWindowStart_T_method()
        {
            return fRunPar->fsumfirst_startingMethod[fTelID];
        }
        double              getSumWindowMaxTimedifferenceToDoublePassPosition()
        {
            return fRunPar->fSumWindowMaxTimedifferenceToDoublePassPosition[fTelID];
        }
        valarray<unsigned int>& getTCorrectedSumFirst()
        {
            return fAnaData[fTelID]->fTCorrectedSumFirst;
        }
        valarray<unsigned int>& getTCorrectedSumLast()
        {
            return fAnaData[fTelID]->fTCorrectedSumLast;
        }
        valarray<unsigned int>& getTCorrectedSum2First()
        {
            return fAnaData[fTelID]->fTCorrectedSum2First;
        }
        valarray<unsigned int>& getTCorrectedSum2Last()
        {
            return fAnaData[fTelID]->fTCorrectedSum2Last;
        }
        unsigned int        getTelescopeEventNumber( unsigned int iTelID )
        {
            if( iTelID < fTelescopeEventNumber.size() )
            {
                return fTelescopeEventNumber[iTelID];
            }
            else
            {
                return 0;
            }
        }
        vector< unsigned int >& getTelescopeEventNumber()
        {
            return fTelescopeEventNumber;
        }
        bool                getTelescopeStatus( unsigned int iTelID );
        unsigned int        getTelID()
        {
            return fTelID;
        }
        unsigned int        getTeltoAnaID()
        {
            return getTeltoAnaID( fTelID );
        }
        unsigned int        getTeltoAnaID( unsigned int iTelID );
        ULong64_t           getTelType( unsigned int iTelID );
        unsigned int        getTelType_Counter( ULong64_t iTelType );
        vector< unsigned int>& getTeltoAna()
        {
            return fTeltoAna;
        }
        double              getTimeSinceRunStart()
        {
            return fAnaData[fTelID]->fTimeSinceRunStart;
        }
        double              getMeanAverageTZero( bool iLowGain = false )
        {
            return fCalData[fTelID]->getAverageTZero( iLowGain );
        }
        TH1F*               getAverageTZeroDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getAverageTzerosetDist( iLowGain );
        }
        TH1F*               getToffsetDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getToffsetDist( iLowGain );
        }
        TH1F*               getToffsetVarsDist( bool iLowGain = false )
        {
            return fCalData[fTelID]->getToffsetVarsDist( iLowGain );
        }
        valarray<double>&   getAverageTZeros( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fAverageTzero;
            }
            else
            {
                return fCalData[fTelID]->fLowGainAverageTzero ;
            }
        }
        valarray<double>&   getAverageTZerosvars( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fAverageTzerovars;
            }
            else
            {
                return fCalData[fTelID]->fLowGainAverageTzerovars;
            }
        }
        valarray<double>&   getTOffsets( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fTOffsets;
            }
            else
            {
                return fCalData[fTelID]->fLowGainTOffsets;
            }
        }
        valarray<double>&   getTOffsetvars( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return fCalData[fTelID]->fTOffsetvars;
            }
            else
            {
                return fCalData[fTelID]->fLowGainTOffsetvars;
            }
        }
        valarray<double>&   getTraceAverageTime( bool iCorrected = true )
        {
            if( iCorrected )
            {
                return fAnaData[fTelID]->fPulseTimingAverageTimeCorrected;
            }
            return fAnaData[fTelID]->fPulseTimingAverageTime;
        }
        unsigned int        getTraceIntegrationMethod()
        {
            return fRunPar->fTraceIntegrationMethod[fTelID];
        }
        unsigned int        getDigitalFilterMethod()
        {
            return fRunPar->fDF_DigitalFilter[fTelID];
        }
        unsigned int        getDigitalFilterUpSample()
        {
            return fRunPar->fDF_UpSample[fTelID];
        }
        float               getDigitalFilterPoleZero()
        {
            return fRunPar->fDF_PoleZero[fTelID];
        }
        bool                hasFADCData()
        {
            return ( bool )fRunPar->fTraceIntegrationMethod[fTelID];
        }
        unsigned int        getTraceIntegrationMethod_pass1()
        {
            return fRunPar->fTraceIntegrationMethod_pass1[fTelID];
        }
        valarray<double>&   getTraceMax()
        {
            return fAnaData[fTelID]->fTraceMax;
        }
        valarray<unsigned int>& getTraceN255()
        {
            return fAnaData[fTelID]->fTraceN255;
        }
        valarray<double>&   getTraceRawMax()
        {
            return fAnaData[fTelID]->fRawTraceMax;
        }
        valarray<double>&   getTraceWidth()
        {
            return fAnaData[fTelID]->getTraceWidth( true );
        }
        VTraceHandler*      getTraceHandler()
        {
            return fTraceHandler;
        }
        vector<bool>&       getTrigger()          // MS
        {
            return fAnaData[fTelID]->fTrigger;
        }
        vector< vector< int > >& getTriggeredTel()
        {
            return fTriggeredTel;
        }
        vector< int >&      getTriggeredTelN()
        {
            return fTriggeredTelN;
        }
        vector< valarray< double > >& getPulseTiming( bool iCorrected = true );
        valarray<double>&   getPulseTime( bool iCorrected = true );
		valarray<double>&   getTTrigger()
		{
			return fAnaData[fTelID]->getTTrigger();
		}
        valarray<double>&   getTZeros( bool iCorrected = true )
        {
            return fAnaData[fTelID]->getTZeros( iCorrected );
        }
        TGraphErrors*       getXGraph( bool iIsSecondPass )
        {
            if( iIsSecondPass )
            {
                return fXGraphDP2[fTelID];
            }
            return fXGraph[fTelID];
        }
        VArrayPointing*     getArrayPointing()
        {
            return fArrayPointing;
        }
        vector< VPointing* > getPointing()
        {
            return fPointing;
        }
        vector<unsigned short int>&       getZeroSuppressed()
        {
            return fAnaData[fTelID]->fZeroSuppressed;
        }
        void                incrementNumberofIncompleteEvents()
        {
            fNumberofIncompleteEvents++;
        }
        void                incrementNumberofGoodEvents()
        {
            fNumberofGoodEvents++;
        }
        //!< set pointer to data reader
        bool                initializeDataReader();
        bool                initializeDeadChannelFinder();
        bool                initializeStarCatalogue( int iMJD, double iTime );
        bool                isDoublePass()
        {
            if( fTelID < fRunPar->fDoublePass.size() )
            {
                return fRunPar->fDoublePass[fTelID];
            }
            return false;
        }
        bool                isEqualSummationWindows();
        bool                isMC()
        {
            if( getReader() )
            {
                return getReader()->isMC();
            }
            
            return false;
        }
        bool                isDST_MC()
        {
            if( isMC() && ( fRunPar->fsourcetype == 4 || fRunPar->fsourcetype == 7 ) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        bool                isTeltoAna( unsigned int iTel );
        void                printDeadChannels( bool iLowGain = false, bool iGrepAble = false );
        void                resetAnaData();
        void                setAnalysisArrayEventStatus( unsigned int i )
        {
            fAnalysisArrayEventStatus = i;
        }
        void                setAnaData()
        {
			fAnaData[fTelID]->initialize( fDetectorGeo->getNChannels( fTelID ), getReader()->getMaxChannels(),
                                                      getDebugFlag(), getRunParameter()->fMCNdeadSeed, getNSamples(), 
                                                      getRunParameter()->fpulsetiminglevels.size(), getRunParameter()->fpulsetiming_tzero_index, 
                                                      getRunParameter()->fpulsetiming_width_index, getRunParameter()->fpulsetiming_triggertime_index  );
        }
        void                setBorder( bool iBo )
        {
            fAnaData[fTelID]->fBorder.assign( fDetectorGeo->getNChannels( fTelID ), iBo );
        }
        void                setBorder( unsigned int iChannel, bool iBo )
        {
            fAnaData[fTelID]->fBorder[iChannel] = iBo;
        }
        void                setBorderCorrelationCoefficient( double iC )
        {
            fAnaData[fTelID]->fCorrelationCoefficient.assign( fDetectorGeo->getNChannels( fTelID ), iC );
        }
        void                setBorderCorrelationCoefficient( unsigned int iChannel, double iC )
        {
            fAnaData[fTelID]->fCorrelationCoefficient[iChannel] = iC;
        }
        void                setBorderThresh( double ithresh )
        {
            if( getImageCleaningParameter() )
            {
                getImageCleaningParameter()->fborderthresh = ithresh;
            }
        }
        void                setBrightNonImageThresh( double ithresh )
        {
            if( getImageCleaningParameter() )
            {
                getImageCleaningParameter()->fbrightnonimagetresh = ithresh;
            }
        }
        void                setCalData()
        {
            fCalData[fTelID]->initialize( fDetectorGeo->getNChannels( fTelID ), getDebugFlag() );
        }
        void                setCalibrated()
        {
            if( fTelID < fCalibrated.size() )
            {
                fCalibrated[fTelID] = true;
            }
        }
        void                setCalibrated( bool iCal )
        {
            if( fTelID < fCalibrated.size() )
            {
                fCalibrated[fTelID] = iCal;
            }
        }
        void                setCurrentSummationWindow( unsigned int iw, bool iSecondWindow )
        {
            if( !iSecondWindow )
            {
                fAnaData[fTelID]->fCurrentSummationWindow = iw;
            }
            else
            {
                fAnaData[fTelID]->fCurrentSummationWindow_2 = iw;
            }
        }
        void                setCurrentSummationWindow( unsigned int imin, unsigned int imax, bool iSecondWindow );
        void                setCurrentSummationWindow( unsigned int iChannel, unsigned int imin, unsigned int imax, bool iSecondWindow );
        void                setDead( unsigned int iDead, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fAnaData[fTelID]->fDead.assign( fDetectorGeo->getNChannels( fTelID ), iDead );
            }
            else
            {
                fAnaData[fTelID]->fLowGainDead.assign( fDetectorGeo->getNChannels( fTelID ), iDead );
            }
        }
        void                setDead( unsigned int iChannel, unsigned int iDead, bool iLowGain = false, bool iFullSet = false, bool iReset = false );
        void                setNDead( unsigned int iN, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fAnaData[fTelID]->fNDead = iN;
            }
            else
            {
                fAnaData[fTelID]->fLowGainNDead = iN;
            }
        }
        void                setDeadChannelText();
        void                setFADCStopOffsets( double iOffset )
        {
            fCalData[fTelID]->fFADCStopOffsets = iOffset;
        }
        void                setFADCStopOffsets( unsigned int iChannel, double iOffset )
        {
            fCalData[fTelID]->fFADCStopOffsets[iChannel] = iOffset;
        }
        void                setFADCtoPhe( double iFADCtoPhe = 1., bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fFADCtoPhe = iFADCtoPhe;
            }
            else
            {
                fCalData[fTelID]->fLowGainFADCtoPhe = iFADCtoPhe;
            }
        }
        void                setFADCtoPhe( unsigned int iChannel, double iFADCtoPhe = 1., bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fFADCtoPhe[iChannel] = iFADCtoPhe;
            }
            else
            {
                fCalData[fTelID]->fLowGainFADCtoPhe[iChannel] = iFADCtoPhe;
            }
        }
        void                setGains( double iGain, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fGains = iGain;
            }
            else
            {
                fCalData[fTelID]->fLowGainGains = iGain;
            }
        }
        void                setGains( unsigned int iChannel, double iGain, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fGains[iChannel] = iGain;
            }
            else
            {
                fCalData[fTelID]->fLowGainGains[iChannel] = iGain;
            }
        }
        void                setGains_DefaultValue( bool iV, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fGains_DefaultSetting = iV;
            }
            else
            {
                fCalData[fTelID]->fLowGainGains_DefaultSetting = iV;
            }
        }
        void                setGainvars( unsigned int iChannel, double iGainvar, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fGainvars[iChannel] = iGainvar;
            }
            else
            {
                fCalData[fTelID]->fLowGainGainvars[iChannel] = iGainvar;
            }
        }
        void                setGainvars( double iGainvar, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fGainvars = iGainvar;
            }
            else
            {
                fCalData[fTelID]->fLowGainGainvars = iGainvar;
            }
        }
        void                setHiLo( bool iHL )
        {
            fAnaData[fTelID]->fHiLo.assign( fDetectorGeo->getNChannels( fTelID ), iHL );
        }
        void                setHiLo( unsigned int iChannel, bool iHL )
        {
            fAnaData[fTelID]->fHiLo[iChannel] = iHL;
        }
        void                setHistoFilling( bool ifill )
        {
            fRunPar->ffillhistos = ifill;
        }
        void                setImage( bool iIm )
        {
            fAnaData[fTelID]->fImage.assign( fDetectorGeo->getNChannels( fTelID ), iIm );
        }
        void                setImageBorderNeighbour( bool iIm )
        {
            fAnaData[fTelID]->fImageBorderNeighbour.assign( fDetectorGeo->getNChannels( fTelID ), iIm );
        }
        void                setBorderBorderNeighbour( bool iIm )
        {
            fAnaData[fTelID]->fBorderBorderNeighbour.assign( fDetectorGeo->getNChannels( fTelID ), iIm );
        }
        void                setImage( unsigned int iChannel, bool iIm )
        {
            fAnaData[fTelID]->fImage[iChannel] = iIm;
        }
        void                setImageThresh( double ithresh )
        {
            if( getImageCleaningParameter() )
            {
                getImageCleaningParameter()->fimagethresh = ithresh;
            }
        }
        void                setImageUser( int iIu )
        {
            fAnaData[fTelID]->fImageUser.assign( fDetectorGeo->getNChannels( fTelID ), iIu );
        }
        void                setImageUser( unsigned int iChannel, int iIu )
        {
            fAnaData[fTelID]->fImageUser[iChannel] = iIu;
        }
        void                setLowGainPedestals()
        {
            fCalData[fTelID]->fBoolLowGainPedestals = true;
        }
        void                setLowGainGains()
        {
            fCalData[fTelID]->fBoolLowGainGains = true;
        }
        void                setLowGainTOff()
        {
            fCalData[fTelID]->fBoolLowGainTOff = true;
        }
        void                setLLEst( vector<bool> iEst )
        {
            fAnaData[fTelID]->fLLEst = iEst;
        }
        void	setLowGainMultiplier_Trace( double lmult )
        {
            fCalData[fTelID]->setLowGainMultiplier_Trace( lmult );
        }
        void	setLowGainMultiplier_Trace( unsigned int iTelID, double lmult )
        {
            if( iTelID < fCalData.size() )
            {
                fCalData[iTelID]->setLowGainMultiplier_Trace( lmult );
            }
        }
        void	setLowGainPedestalFile( string file )
        {
            fCalData[fTelID]->setLowGainPedestalFile( file );
        }
        void	setLowGainMultiplier_Sum( int iMethod, int iSumWindow, int jSumWindow, double lmult )
        {
            fCalData[fTelID]->setLowGainMultiplier_Sum( iMethod, iSumWindow, jSumWindow, lmult );
        }
        
        void                setNChannels( unsigned int iChan )
        {
            fDetectorGeo->setNChannels( fTelID, iChan );
        }
        void                setBrightNonImage( bool iIm )
        {
            fAnaData[fTelID]->fBrightNonImage.assign( fDetectorGeo->getNChannels( fTelID ), iIm );
        }
        void                setBrightNonImage( unsigned int iChannel, bool iIm )
        {
            fAnaData[fTelID]->fBrightNonImage[iChannel] = iIm;
        }
        void                setNoPointing( bool iP )
        {
            fNoTelescopePointing = iP;
        }
        void                setNSamples( unsigned int iSamp )
        {
            fDetectorGeo->setNSamples( fTelID, iSamp, fRunPar->fUseVBFSampleLength );
        }
        void                setNSamples( unsigned int iTelID, unsigned int iSamp )
        {
            fDetectorGeo->setNSamples( iTelID, iSamp, fRunPar->fUseVBFSampleLength );
        }
        void                setDebugLevel( int i );
        void                setZeroSuppressed( unsigned short int iZ )
        {
            fAnaData[fTelID]->fZeroSuppressed.assign( fDetectorGeo->getNChannels( fTelID ), iZ );
        }
        void                setZeroSuppressed( unsigned int iChannel, unsigned short int iZ )
        {
            if( iChannel < fAnaData[fTelID]->fZeroSuppressed.size() )
            {
                fAnaData[fTelID]->fZeroSuppressed[iChannel] = iZ;
            }
        }
        
        void             setClusterNpix( int iID, int clusterNpix )
        {
            fAnaData[fTelID]->fClusterNpix[iID] = clusterNpix;
        }
        vector<int>&     getClusterNpix()
        {
            return fAnaData[fTelID]->fClusterNpix;
        }
        void             setClusterID( unsigned int iChannel, int iID )
        {
            fAnaData[fTelID]->fClusterID[iChannel] = iID;
        }
        vector<unsigned int>&     getClusterID()
        {
            return fAnaData[fTelID]->fClusterID;
        }
        void             setMainClusterID( int iID )
        {
            fAnaData[fTelID]->fMainClusterID = iID;
        }
        int              getMainClusterID()
        {
            return fAnaData[fTelID]->fMainClusterID;
        }
        
        void             setClusterSize( int iID, double clustersize )
        {
            fAnaData[fTelID]->fClusterSize[iID] = clustersize;
        }
        vector<double>&  getClusterSize()
        {
            return fAnaData[fTelID]->fClusterSize;
        }
        void             setClusterTime( int iID, double clustertime )
        {
            fAnaData[fTelID]->fClusterTime[iID] = clustertime;
        }
        vector<double>&  getClusterTime()
        {
            return fAnaData[fTelID]->fClusterTime;
        }
        
        void             setClusterCenx( int iID, double clustercenx )
        {
            fAnaData[fTelID]->fClusterCenx[iID] = clustercenx;
        }
        vector<double>&  getClusterCenx()
        {
            return fAnaData[fTelID]->fClusterCenx;
        }
        void             setClusterCeny( int iID, double clusterceny )
        {
            fAnaData[fTelID]->fClusterCeny[iID] = clusterceny;
        }
        vector<double>&  getClusterCeny()
        {
            return fAnaData[fTelID]->fClusterCeny;
        }
        
        void             setNcluster_cleaned( int i_Ncluster )
        {
            fAnaData[fTelID]->fncluster_cleaned = i_Ncluster;
        };
        int              getNcluster_cleaned()
        {
            return fAnaData[fTelID]->fncluster_cleaned;
        };
        void             setNcluster_uncleaned( int i_Ncluster )
        {
            fAnaData[fTelID]->fncluster_uncleaned = i_Ncluster;
        };
        int              getNcluster_uncleaned()
        {
            return fAnaData[fTelID]->fncluster_uncleaned;
        };
        void  setIPRGraph( unsigned int iSumWindow, TGraphErrors* g )
        {
            return fCalData[fTelID]->setIPRGraph( iSumWindow, g );
        }
        /////////////// pedestals /////////////////////
        void                setPeds( unsigned int iChannel, double iPed, bool iLowGain = false )
        {
            fCalData[fTelID]->setPeds( iChannel, iPed, iLowGain );
        }
        void                setPedsFromPLine()
        {
            fCalData[fTelID]->fPedFromPLine = true;
        }
        /////////////// end pedestals /////////////////////
        void                setRootDir( unsigned int iTel, TDirectory* iDir )
        {
            fAnaDir[iTel] = iDir;
        }
        void                setRunNumber( int iRunN )
        {
            fRunNumber = iRunN;
        }
        void                setSumFirst( int iSum )
        {
            fRunPar->fsumfirst[fTelID] = iSum;
        }
        void                setSums( double iSum )
        {
            fAnaData[fTelID]->fSums = iSum;
        }
        void                setSums( unsigned int iChannel, double iSum )
        {
            fAnaData[fTelID]->fSums[iChannel] = iSum;
        }
        bool                setSums( valarray< double > iVSum );
        void                setSums2( double iSum )
        {
            fAnaData[fTelID]->fSums2 = iSum;
        }
        void                setSums2( unsigned int iChannel, double iSum )
        {
            fAnaData[fTelID]->fSums2[iChannel] = iSum;
        }
        void                setSums2( valarray< double > iVSum )
        {
            fAnaData[fTelID]->fSums2 = iVSum;
        }
        void                setTemplateMu( valarray< double > iVTemplateMu )
        {
            fAnaData[fTelID]->fTemplateMu = iVTemplateMu;
        }
        void                setTCorrectedSumFirst( unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSumFirst = iT;
        }
        void                setTCorrectedSumFirst( unsigned int iChannel, unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSumFirst[iChannel] = iT;
        }
        void                setTCorrectedSumLast( unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSumLast = iT;
        }
        void                setTCorrectedSumLast( unsigned int iChannel, unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSumLast[iChannel] = iT;
        }
        void                setTCorrectedSum2First( unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSum2First = iT;
        }
        void                setTCorrectedSum2First( unsigned int iChannel, unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSum2First[iChannel] = iT;
        }
        void                setTCorrectedSum2Last( unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSum2Last = iT;
        }
        void                setTCorrectedSum2Last( unsigned int iChannel, unsigned int iT )
        {
            fAnaData[fTelID]->fTCorrectedSum2Last[iChannel] = iT;
        }
        void                setTelID( unsigned int iTel );
        void                setTeltoAna( vector< unsigned int > iT );
        void                setTOffsets( double iToff, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fTOffsets = iToff;
            }
            else
            {
                fCalData[fTelID]->fLowGainTOffsets = iToff;
            }
        }
        void                setTOffsets( unsigned int iChannel, double iToff, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fTOffsets[iChannel] = iToff;
            }
            else
            {
                fCalData[fTelID]->fLowGainTOffsets[iChannel] = iToff;
            }
        }
        void                setTOffsetvars( double iToffv, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fTOffsetvars = iToffv;
            }
            else
            {
                fCalData[fTelID]->fLowGainTOffsetvars = iToffv;
            }
        }
        void                setTOffsetvars( unsigned int iChannel, double iToffv, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fTOffsetvars[iChannel] = iToffv;
            }
            else
            {
                fCalData[fTelID]->fLowGainTOffsetvars[iChannel] = iToffv;
            }
        }
        void                setAverageTZero( double iTZero, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fAverageTzero = iTZero;
            }
            else
            {
                fCalData[fTelID]->fLowGainAverageTzero = iTZero;
            }
        }
        bool                setAverageTZero( unsigned int iChannel, double iTZero, bool iLowGain = false );
        void                setAverageTZerovars( double iTZerovars, bool iLowGain = false )
        {
            if( !iLowGain )
            {
                fCalData[fTelID]->fAverageTzerovars = iTZerovars;
            }
            else
            {
                fCalData[fTelID]->fLowGainAverageTzerovars = iTZerovars;
            }
        }
        bool                setAverageTZerovars( unsigned int iChannel, double iTZero, bool iLowGain = false );
        void                setMeanAverageTZero( double iTZero, bool iLowGain = false )
        {
            fCalData[fTelID]->setAverageTZero( iTZero, iLowGain );
        }
        void                setTrace( unsigned int iChannel, vector< double > fT, bool iHiLo, double iPeds )
        {
            fAnaData[fTelID]->setTrace( iChannel, fT, iHiLo, iPeds );
        }
        void                setTraceAverageTime( double iT )
        {
            fAnaData[fTelID]->fPulseTimingAverageTime = iT;
        }
        void                setTraceAverageTime( unsigned int iChannel, double iT )
        {
            fAnaData[fTelID]->fPulseTimingAverageTime[iChannel] = iT;
        }
        void                setTraceMax( unsigned int iChannel, double iS )
        {
            fAnaData[fTelID]->fTraceMax[iChannel] = iS;
        }
        void                setTraceMax( double iV )
        {
            fAnaData[fTelID]->fTraceMax = iV;
        }
        void                setTraceMax( valarray< double > iV )
        {
            fAnaData[fTelID]->fTraceMax = iV;
        }
        void                setTraceN255( unsigned int iS )
        {
            fAnaData[fTelID]->fTraceN255 = iS;
        }
        void                setTraceN255( unsigned int iChannel, unsigned int iS )
        {
            fAnaData[fTelID]->fTraceN255[iChannel] = iS;
        }
        void                setTraceRawMax( unsigned int iChannel, double iS )
        {
            fAnaData[fTelID]->fRawTraceMax[iChannel] = iS;
        }
        void                setTraceRawMax( double iV )
        {
            fAnaData[fTelID]->fRawTraceMax = iV;
        }
        void                setTraceRawMax( valarray< double > iV )
        {
            fAnaData[fTelID]->fRawTraceMax = iV;
        }
        void                setTraceWidth( double iV )
        {
            fAnaData[fTelID]->fTraceWidth = iV;
        }
        void                setTraceWidth( unsigned int iChannel, double iS )
        {
            fAnaData[fTelID]->fTraceWidth[iChannel] = iS;
        }
        void                setTraceWidth( valarray< double > iV )
        {
            fAnaData[fTelID]->fTraceWidth = iV;
        }
        void                setTrigger( bool iIm )
        {
            fAnaData[fTelID]->fTrigger.assign( fDetectorGeo->getNChannels( fTelID ), iIm );
        }
        void                setTrigger( unsigned int iChannel, bool iIm )
        {
            fAnaData[fTelID]->fTrigger[iChannel] = iIm;
        }
        void                setPulseTiming( vector< valarray< double > > iPulseTiming, bool iCorrected );
        void                setPulseTiming( float iTZero, bool iCorrected );
        void                setPulseTiming( unsigned int iChannel, vector< float > iTZero, bool iCorrected );
        void                setPulseTiming( unsigned int iChannel, bool iCorrection, double iValue );
        void                setPulseTimingCorrection( unsigned int iChannel, double iCorrection );
        void                setXGraph( TGraphErrors* igraph, bool iSecondPass )
        {
            if( iSecondPass )
            {
                fXGraphDP2[fTelID] = igraph;
            }
            else
            {
                fXGraph[fTelID] = igraph;
            }
        }
        bool                usePedestalsInTimeSlices( bool iLowGain = false )
        {
            if( !iLowGain )
            {
                return getRunParameter()->fUsePedestalsInTimeSlices;
            }
            else
            {
                return getRunParameter()->fLowGainUsePedestalsInTimeSlices;
            }
        }
        
        void                testDataReader();     //!< check if reader is available, set pointers
};
#endif
