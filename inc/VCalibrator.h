//! VCalibrator calibration class (calculation of pedestals/gains/toffsets), reading and writing of all calibration data

#ifndef VCALIBRATOR_H
#define VCALIBRATOR_H

#include "VImageBaseAnalyzer.h"
#include "VPedestalCalculator.h"
#include "VDB_CalibrationInfo.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TTree.h"

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VCalibrator : public VImageBaseAnalyzer
{
    private:
        int fCalibrationfileVersion;

        map< ULong64_t, int > fNumberPedestalEvents;        //!< number of events used in pedestal analysis
        vector<int> fNumberGainEvents;            //!< number of events used in gain and toffset analysis
        vector<int> fNumberTZeroEvents;

        TFile* fPedSingleOutFile;
        map< ULong64_t, TFile* > fPedOutFile;
        map< ULong64_t, vector< vector<TH1F* > > > hpedPerTelescopeType;  //<! one histogram per teltype/channel/sumwindow
        map< ULong64_t, vector< vector<TH1F* > > > hped_vec;     //<! one histogram per telescope/channel/sumwindow
        map< ULong64_t, TClonesArray* > fPedestalsHistoClonesArray;
        TFile* opfgain;
        TFile* opftoff;
        vector<TH1F* > hgain;
        vector<TProfile* > hpulse;
        vector<TProfile* > htcpulse;
        vector<TH1F* > htoff;
        vector<TProfile* > htoff_vs_sum;
        int fPedPerTelescopeTypeMinCnt;                         // statistical limit for IPR calculation

        //Extra calib output.
        TTree* tExtra_ChargeTree;
        vector<double>* fExtra_sum;
        vector<double>* fExtra_tzero;
        vector<short>* fExtra_HiLo;
        vector<short>* fExtra_sumfirst;
        vector<short>* fExtra_sumwindow;
        vector<short>* fExtra_dead;
        vector<short>* fExtra_use;
        vector<double>* fExtra_ped;
        vector<double>* fExtra_pedVar;
        double fExtra_QMon;
        double fExtra_TZeroMon;
        int fExtra_nHiLo;
        int fExtra_nMon;
        int fExtra_nPix;
        int fExtra_eventNumber;

        // average Tzero calculation
        vector< TFile* > fTZeroOutFile;
        // one histogram per telescope and channel
        vector< vector< TH1F* > > htzero;
        vector< vector< TH1F* > > htaverage;

        vector< string > fPedFileNameC;
        vector< string > fGainFileNameC;
        vector< string > fToffFileNameC;
        vector< string > fPixFileNameC;
        vector< string > fTZeroFileNameC;
        vector< bool > fBlockTel;
        vector< string > fNewLowGainPedFileNameC;
        vector< string > fLowGainPedFileNameC;
        vector< string > fLowGainGainFileNameC;
        vector< string > fLowGainToffFileNameC;
        vector< string > fLowGainMultiplierNameC;
        vector< string > fLowGainTZeroFileNameC;

        TTree* fillCalibrationSummaryTree( unsigned int itel, string iName, vector<TH1F* > h );
        bool   fillPedestalTree( unsigned int tel, VPedestalCalculator* iP );
        bool   initializePedestalHistograms( ULong64_t iTelType, bool iLowGain,
                                             vector< double > minSumPerSumWindow,
                                             vector< double > maxSumPerSumWindow );
        void getCalibrationRunNumbers();
        int  getCalibrationRunNumbers_fromCalibFile();
        unsigned int getNumberOfEventsUsedInCalibration( vector< int > iE, int iTelID );
        unsigned int getNumberOfEventsUsedInCalibration( map< ULong64_t, int > iE, int iTelID );
        TFile* getPedestalRootFile( ULong64_t iTel );
        int  readLowGainCalibrationValues_fromCalibFile( string iVariable = "LOWGAINPED", unsigned int iTel = 9999 );
        string getCalibrationFileName( int iTel, int irun, string iSuffix, string name = "" );
        void readCalibrationData();
        bool readCalibrationDatafromDSTFiles( string iSourceFile, bool iPedOnly = false );
        void readfromVOFFLINE_DB( int gain_or_toff, string& iFile, vector< unsigned int >& VchannelList, vector< double >& Vmean, vector< double >& Vrms );
        void readGains( bool iLowGain = false );
        bool readIPRGraph_from_DSTFile( string iDSTFile, unsigned int iSummationWindow, ULong64_t iTelType );
        bool calculateIPRGraphs();
        bool calculateIPRGraphs( string iPedFileName, unsigned int iSummationWindow, ULong64_t iTelType, unsigned int i_tel );
        bool readLowGainMultiplier( );
        bool readPeds( string iFile, bool, unsigned int );
        bool readPeds_from_grisufile( bool, unsigned int );
        bool readPeds_from_rootfile( string iFile, bool, unsigned int );
        bool readPeds_from_textfile( string iFile, bool iLowGain, unsigned int i_SumWindow );
        bool readPeds_from_combinedfile( string iFile, bool iLowGain, unsigned int i_SumWindow );
        void readPixelstatus();
        void readTOffsets( bool iLowGain = false );
        bool readAverageTZeros( bool iLowGain = false );
        void setCalibrationFileNames();

        void writeGains( bool iLowGain = false );
        void writePeds( bool iLowGain, VPedestalCalculator* iP = 0, bool iWriteAsciiFile = true );
        void writeTOffsets( bool iLowGain = false );
        void writeAverageTZeros( bool iLowGain = false );
        bool writeIPRgraphs( string iFile = "" );


    public:
        VCalibrator();
        ~VCalibrator() {}

        void calculateAverageTZero( bool iLowGain = false );
        void calculatePedestals( bool iLowGain = false );
        void calculateGainsAndTOffsets( bool iLowGain = false );
        unsigned int getNumberOfEventsUsedInCalibration( int iTelID, int iType );
        void initialize();
        void terminate( VPedestalCalculator* );
};
#endif
