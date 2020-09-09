//! VEvndispRunParameter  input parameter storage class
//
//  IMPORTANT   IMPORTANT    IMPORTANT   IMPORTANT   IMPORTANT    IMPORTANT
//
//   ClassDef VERSION NUMBER HAS TO BE INCREASED BY ONE AFTER ANY CHANGES HERE
//
//  IMPORTANT   IMPORTANT    IMPORTANT   IMPORTANT   IMPORTANT    IMPORTANT

#ifndef VEvndispRunParameter_H
#define VEvndispRunParameter_H

#include <TDatime.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#include "VGlobalRunParameter.h"
#include "VImageCleaningRunParameter.h"

using namespace std;

class VEvndispRunParameter : public TNamed, public VGlobalRunParameter
{
    public:
    
        // local host parameters
        string fEventDisplayUser;
        string fEventDisplayHost;
        string fEventDisplayDate;
        
        // computer system parameters
        string fEventDisplayBuildCompiler;
        string fEventDisplayBuildCompilerVersion;
        string fEventDisplayBuildArch;
        string fEventDisplayBuildNode;
        string fEventDisplayBuildDir;
        SysInfo_t* fEventDisplaySystemInfo;
        
        // root parameters
        string fEventDisplayBuildROOTVersion;
        int    fEventDisplayBuildROOTVersionInt;
        
        // JOB ids
        int    fSGE_TASK_ID;
        
        // DB parameters
        bool   fuseDB;
        string fDBRunType;                        // run type according to DB
        string fDBRunStartTimeSQL;                // run start (SQL type)
        string fDBRunStoppTimeSQL;                // run stopp (SQL type)
        double fDBDataStartTimeMJD;               // run start (mjd)
        double fDBDataStoppTimeMJD;               // run stopp (mjd)
        double fDBDataStartTimeSecOfDay;          // run start (sec of day)
        double fDBDataStoppTimeSecOfDay;          // run stopp (sec of day)
        
        // run parameters
        int    frunnumber;                        // runnumber
        string fRunTitle;                         // run title (usually empty, useful for debugging)
        int    frunmode;                          // run mode
        //    (0=analysis, 1=pedestal calculation, 2=gain/toffset calculation,
        //    5 = gain/toffset calculation for low gain channels, 6=pedestals for low gain channels)
        int    fsourcetype;                       // source type (0=rawdata,1=GrIsu,2=MC in vbf format,3=Rawdata in vbf,4=DST,5=multiple GrIsu file)
        bool   fRunIsZeroSuppressed;              // run is zero suppressed
        string fsourcefile;                       // name of data file
        int    fnevents;                          // total number of events to be analyzed
        int    fFirstEvent;                       // skip up till this event
        int    fTimeCutsMin_min;                  // start to analyse run at this min
        int    fTimeCutsMin_max;                  // stop to analyse this run at this min
        
        bool fprintdeadpixelinfo ;       // DEADCHAN if true, will print list of dead pixels
        // at end of run to evndisp.log
        
        bool fSaveDeadPixelRegistry ; // will save an extra tree containing info about the dead pixels
        
        bool fForceLLImageFit ;          // FORCELL if true, will use log-likelihood image fitting
        // on all images, irregardless of if they are near the
        // edge of the camera or not.  Set in
        // EVNDISP.reconstruction.runparameter
        bool fSquaredImageCalculation;   // use unusual definition of image calculation (weighted by square)
        
        string fInstrumentEpoch;                  // Instrumental epoch (e.g. for VTS V5 or V6)
        unsigned int fAtmosphereID;               // corsika ID of atmosphere
        vector< float > fEpochGain;               // absolute gain (dc/pe) during epoch (e.g. for VTS V5 or V6)
        string fEpochFile;                        // file with list of epochs and atmospheres
        
        float  fRunDuration;                      // duration of runs in [s]
        
        // output parameters
        bool   ffillhistos;                       // fill some analysis histograms
        string foutputfileName;                   // file name for analysis output (root file), -1 = no output
        
        // debugging
        bool  fDebug;                             // print debug output
        bool  fPrintSmallArray;                   // some printout for small arrays only
        unsigned int fPrintAnalysisProgress;      // print a line each time this number of events have been processed
        unsigned int fPrintGrisuHeader;           // Print full grisu header, includeing detector config file used in
        //   grisu simulations (sims only, from VBF header), OR name of
        //   config file for detector simulation if available
        
        // array/telescope geometry parameters
        unsigned int fNTelescopes;                // number of telescopes
        vector<string> fcamera;                   // name of camera configuration files
        vector< unsigned int > fTelToAnalyze;     // analyze only this telescope (Telescope 1 = 0!! )
        bool   fUseVBFSampleLength;               // use number of samples from VBF file (ignore .cfg file)
        
        // pointing parameters
        string fTargetName;                       // target name
        double fTargetDec;                        // target declination [deg] (J2000)
        double fTargetRA;                         // target RA [deg] (J2000)
        double fTargetDecOffset;                  // target offset in dec [deg]
        double fTargetRAOffset;                   // target offset in RA [deg]
        double fWobbleNorth;                      // wobble offset NORTH (deg)
        double fWobbleEast;                       // wobble offset EAST (deg)
        double felevation;                        // elevation in [deg] (preli)
        double fazimuth;                          // azimuth in [deg] (preli)
        bool fDBTracking;                         // use tracking from DB
        string fDBTrackingCorrections;            // apply tracking corrections from this period (SQL time string), empty for use of measured values in db
        string fPMTextFileDirectory;              // pointing monitor text file directory
        bool fDBVPM;                              // use calibrated VPM tracking from database
        bool fDBUncalibratedVPM;                  // use uncalibrated VPM tracking from database
        vector<double> fPointingErrorX;           // pointing error, in camera coordinates [deg]
        vector<double> fPointingErrorY;           // pointing error, in camera coordinates [deg]
        double fCheckPointing;                    // stop run if pointing difference between calculated pointing direction and vbf is larger than this value
        bool fDBCameraRotationMeasurements;       // read camera rotations from DB
        
        // calibration (pedestals, gains, etc.)
        bool   fcalibrationrun;                   // true if this is a pedestal/gain/toff calculation run
        string fcalibrationfile;                  // file with file names for calibration
        bool freadCalibfromDB;                    // if true, calibration information are read in VOFFLINE DB
        int freadCalibfromDB_versionquery;        // require a given version of calibration
        bool freadCalibfromDB_save_file;          // calibration information read from the DB are stored in
        // VGlobalRunParameter::getDirectory_EVNDISPCalibrationData() +/Tel_?
        bool fNoCalibNoPb;                        // if true, when no information for gain and toff can be found, the analysis is done filling thenm with 1 and 0 respectively (in VCalibrator)
        bool fNextDayGainHack;            //if true, and > 100 channels in one telescope have gain=0, all gains in that tel will be set to 1; gains won't be tested in the dead channel finder.
        bool fWriteExtraCalibTree;        // write additional tree into .gain.root file with channel charges/monitor charge/nHiLo for each event
        string fLowGainCalibrationFile;           // file with file name for low-gain calibration
        int fNCalibrationEvents;                  // events to be used for calibration
        float faverageTZeroFiducialRadius;        // fiducial radius for average tzero calculation (DST), in fraction of FOV
        unsigned int fCombineChannelsForPedestalCalculation; // combine all channels per telescope type for the pedestal calculation
        vector< int > fGainFileNumber;
        vector< int > fTOffFileNumber;
        vector< int > fPedFileNumber;
        vector< int > fTZeroFileNumber;
        vector< int > fPedLowGainFileNumber;
        string fPedLowGainFile;
        vector< int > fGainLowGainFileNumber;
        vector< int > fTOffLowGainFileNumber;
        vector< int > fLowGainMultiplierFileNumber;
        vector< int > fTZeroLowGainFileNumber;
        vector< int > fPixFileNumber;
        string fDeadChannelFile;
        float fCameraCoordinateTransformX;        // multiply all X coordinates in camera by this value
        float fCameraCoordinateTransformY;        // multiply all Y coordinates in camera by this value
        vector< double > fGainCorrection;         // all gains are divided by this number
        bool   fIgnoreDSTGains;
        double fLaserSumMin;                      // minimal total charge sum for a event to be a laser event
        bool   fUsePedEvents;                     // use only true pedestal events (evttype=2) for pedestal analysis (default: off )
        bool   fLowGainPeds;                      // calculate pedestals from low gain channels only
        bool   fUsePedestalsInTimeSlices;         // use pedestal in time slices for image calculation (high gain)
        bool   fLowGainUsePedestalsInTimeSlices;  // use pedestals in time slices for image calculation (low gain)
        bool   fPedestalsInTimeSlices;            // calculating time pedendent pedestals
        double fPedestalsLengthOfTimeSlice;       // length of a time slice for pedestal calculations (tracking test)
        bool   fPedestalSingleRootFile;           // write pedestal trees and histograms into a single root file
        int    fCalibrationSumWindow;             // sumwindow for all calibration calculation
        int    fCalibrationSumFirst;              // starting point all calibration calculation
        int    fCalibrationSumWindowAverageTime;  // sumwindow for average arrival time calculation
        float  fCalibrationIntSumMin;             // minimum integrated charge in a channel and event to be taken into account in gain or tzero calibration runs
        string fsetSpecialChannels;               // set channels with L2 channels to correct for FADC crate time offsets (file name of channel settings)
        bool   fL2TimeCorrect;                    // use L2 pulses to correct FADC times (default: on )
        unsigned fCalibrationDataType;            // for DSTs: kind of calibration data available: 1: full (peds, pedvars, etc). 0: (no calibration data)
        
        //////////////////////////////////////////////////
        // FADC integration
        string  fFADCChargeUnit;                  // FADC charge unit (DC or PE)
        vector< unsigned int > fTraceIntegrationMethod;   // trace integration method
        vector< unsigned int > fTraceIntegrationMethod_pass1;   // trace integration method for pass 1 (doublepass only)
        vector< unsigned int > fDF_DigitalFilter; // use a digital filter in the trace analysis
        vector< unsigned int > fDF_UpSample;      // digital filter: up sample parameter
        vector< float >        fDF_PoleZero;      // digital filter: pole-zero cancellation parameter
        vector<int> fsumfirst;                    // parameter for window summation start (window 1)
        vector<unsigned int>   fSearchWindowLast; // last sample searched for trace integration
        vector<int> fsumwindow_1;                 // parameter for window summation (window 1)
        vector<int> fsumwindow_2;                 // parameter for window summation (window 2)
        vector<int> fsumwindow_pass1;             // parameter for window summation (double pass - pass 1)
        bool   fFixWindowStart_sumwindow2;        // fixed window start for summation window 2
        vector<bool> fDoublePass;                       // double pass image cleaning (default: on)
        vector<bool> fDoublePassErrorWeighting2005;     // use error weighting from 2004 or today
        vector< float >    fTraceWindowShift;       // shift the summation window by value (in doublepass: low gain channels only, default: 0 )
        vector< unsigned int >   fsumfirst_startingMethod;   // start the summation window at T0, TAverage, ... (+shift; not for doublepass)
        vector< int >  fsumfirst_maxT0startDiff;              // maximum difference between T0 and trace integration window start
        vector< double > fSumWindowMaxTimeDifferenceLGtoHG;   // maximum difference between lg and hg window in doublepass method
        vector< double > fSumWindowMaxTimedifferenceToDoublePassPosition; // maximum difference between doublepass calculated window start and t0 (in samples, default: 10 )
        
        float  fFADCtoPhe[VDST_MAXTELTYPES];                  //! default conversion factor c[phes/fadc]: [phes]=c*[fadc] for certain integ. window (4slices)
        
        // FADC timing parameters
        vector< float > fpulsetiminglevels;       // levels at which timing of FADC pulses is calculated
        unsigned int fpulsetiming_tzero_index;
        unsigned int fpulsetiming_width_index;
        unsigned int fpulsetiming_max_index;
        unsigned int fpulsetiming_triggertime_index;
        
        // image cleaning
        vector< VImageCleaningRunParameter* >  fImageCleaningParameters;
        vector< VImageCleaningRunParameter* >  fImageCleaningParameters_DB_Pass1;      // image cleaning parameters for double pass, pass1
        
        // image analysis
        int    fImageLL;                          // loglikelihood image parameterisation 0=off/1=on/2=verbose mode (default: 0=off )
        // neighbour scaling
        vector< float > fNeighbourDistanceFactor;
        vector< bool > fSquarePixels;
        vector< unsigned int > fNeighbourMultiplicity;
        // image analysis
        vector< float > fLogLikelihoodLoss_min;   // do loglikelihood image parameterisation if loss is larger than this value
        vector< float > fLogLikelihoodLoss_max;   // do loglikelihood image parameterisation if loss is smaller than this value
        vector< int > fLogLikelihood_Ntubes_min; //  do loglikelihood image parameterisation if ntubes is larger than this value
        float  fImageAnalysisFUIFactor;           // FUI factor for image analysis
        bool   frecoverImagePixelNearDeadPixel;
        bool   fFillImageBorderNeighbours;
        bool   fSmoothDead;                       // smooth dead pixels (default: off )
        // NN cleaning
        bool   ifWriteGraphsToFile;               // flag to write run-time NN image cleaning settings to root file (read from NN image cleaning input card)
        bool   ifReadIPRfromDatabase;             // flag to read IPR from external IPR database
        bool   ifCreateIPRdatabase;               // flag to create external IPR database
        bool   ifReadIPRfromDSTFile;              // flat to read IPR graph from DST file
        
        string freconstructionparameterfile;      // reconstruction parameter file
        // MC parameters
        string fsimu_pedestalfile;                // use external pedestal file for MC
        int    fsimu_noiselevel;                  // noise level used for this external file
        double fsimu_pedestalfile_DefaultPed;     // default pedestal for Grisu MC (must be the same is in grisu .cfg file)
        double fsimu_lowgain_pedestal_DefaultPed; // default low-gain pedestal level for MC
        bool   fsimu_HILO_from_simFile;           // read hilo multiplier from MC header in vbf file
        int    fIsMC;                             // source data is Monte Carlo (==1 write all MC events, ==2 write only triggered events)
        bool   fIgnoreCFGversions;                // ignore configuration file versions
        double finjectGaussianNoise;              // add Gaussian noise to trace (standard deviation in units of DC)
        UInt_t finjectGaussianNoiseSeed;          // seed for random Gaussian noise to trace adding
        
        int fgrisuseed;
        int ftelescopeNOffset;                    // offset in telescope numbering (default: first tele = 0 ), for GrIsu Reader
        int fsampleoffset;                        // offset in samples
        unsigned int fMC_FADCTraceStart;          // start of FADC trace (in case the simulated trace is longer than needed)
        string ftracefile;                        // background trace library file for MC
        bool fMCnoDead;                           // ignore all dead channel labels in camera configuration files
        int fMCNdead;                             // number of pixels set randomly dead
        int fMCNdeadSeed;                         // seed for setting pixels randomly dead
        int fMCNdeadboard;                        // number of boards set randomly dead (10 dead pixels in a row)
        double fMCScale;                          // scale factor for MC data
        
        // tree filling
        unsigned int fShortTree;                  // 0: full tree; 1: short tree
        unsigned int fwriteMCtree;                // 0: do not write MC tree
        bool fWriteTriggerOnly;                   // true: write triggered events for simulation only
        bool fFillMCHistos;                       // true: fill MC histograms with thrown events
        
        // display parameters
        bool   fdisplaymode;                      // display mode or command line mode
        bool   floopmode;                         // infinite event loop
        unsigned int fw;                          // main window width
        unsigned int fh;                          // main window height
        bool fPlotRaw;                            // plot raw values only, used fpeds = 0
        bool fPlotPaper;                          // clean up plots for papers and talks (no dead channels, no small text, ...)
        unsigned int fPlotAllInOneMethod;         // from which method are the angular reconstruction results to taken to plot in all in one window
        
        // star catalogue
        string fStarCatalogueName;
        float  fMinStarBrightness_B;
        float  fMinStarPixelDistance_deg;        // closest distance of a pixel to a bright star
        int    fMinStarNTubes;                   // closest distance analysis applies only for image smaller than this number
        
        // muon parameters
        bool fmuonmode;                           // search for complete muon rings, Martin
        // Hough transform muon parameters
        bool fhoughmuonmode;                      // Use hough transform muon analysis
        
        // Frogs parameters
        bool ffrogsmode;                          // for Frogs template Analysis, GH
        string ffrogsmscwfile;                      // frogs file for getting table energy
        int ffrogsRecID;                              // RecID or Cut_ID Frogs Uses - combine with table
        string ffrogstemplatelist;                // text file containing the elevation and
        //   elevation-dependant template file names
        string ffrogsparameterfile;               // parameter file for frogs settings
        
        // write pulse histograms to gain files
        int  fwriteLaserPulseN;                    // number of pulse histogram written to gain file
        bool fwriteAverageLaserPulse;              // write average laser pulse to file
        
        // dst parameters
        string fdstfile;                          // dst output file name (root file)
        int fdstminntubes;                        // write only events with more than fdstminntubes ntubes into dst file
        bool fdstwriteallpixel;                   // write all information of all pixel into dst output files
        bool fdstcalibration;                     // write only the branches required for calibration to DST file
        
        // trigsim parameters (note: different telescope IDs!)
        map< unsigned int, int >   fTrigSim_referenceTrigger;         //!
        map< unsigned int, float > fTrigSim_threshold;                //!
        TString  fTrigSimInputcard;                      // input card for trigsim and next-neighbour image cleaning
        TString  fNNGraphsFile;
        TString  fIPRdatabase;                    // file to read IPRs from external database
        TString  fIPRdatabaseFile;                // file to write the IPR database
        
        // Parallaxwidth
        int fPWmethod;                            // how to make the trigger-map to calculate the trigger-level image parameters
        int fPWcleanNeighbors;                    // number of neighbors required for a center pixel to survive the cleaning procedure
        float fPWcleanThreshold;                  // cleaning threshold to use to determine hit pixels from the summed FADC charge (dc)
        int fPWlimit;                             // limits the number of pixels transmitted per sector, if =0, then the function is ignored and no cut is applied on the generation of the trigger map
        
        // functions
        void print();
        void print( int iEV );
        void printCTA_DST();
        
        VEvndispRunParameter( bool bSetGlobalParameter = true );
        ~VEvndispRunParameter();
        
        bool         isMC()
        {
            return fIsMC;
        }
        void         setDefaultParameters( string iObservatory );
        void         setPulseZeroIndex();
        void         setSystemParameters();
        bool         useDB()
        {
            return fuseDB;
        }
        
        ClassDef( VEvndispRunParameter, 186 ); //(increase this number)
};
#endif
