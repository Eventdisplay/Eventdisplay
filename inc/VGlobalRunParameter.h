//! VGlobalRunParameter global definitions for eventdisplay package

#ifndef VGlobalRunParameter_H
#define VGlobalRunParameter_H

#include <TChain.h>
#include <TSystem.h>
#include <TTree.h>
#include <TROOT.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>

//////////////////////////////////////////////////////////////////
// HARDWIRED MAXIMUM NUMBER OF TELESCOPES AND CHANNELS, etc.
// (changing anything here means that you have to rerun all stages
//  of the analysis again)
#define VDST_MAXTELESCOPES  125    // maximum number of telescopes
#define VDST_MAXTELTYPES      7    // maximum number of telescope types
#define VDST_MAXNNGROUPTYPES  6    // maximum number of NN-group types searched in NN-image cleaning
#ifndef CTA_PROD3_SC
#define VDST_MAXCHANNELS  12000    // maximum number of channels per telescopes
#else
#define VDST_MAXCHANNELS  12000    // maximum number of channels per telescopes
#endif
#define VDST_MAXSUMWINDOW   130    // maximum number of summation windows (=maximum number of samples per FADC trace)
#define VDST_PEDTIMESLICES 5000    // maximum number of time slices for pedestal calculation
#define VDST_MAXRECMETHODS  100    // maximum number of arrayreconstruction method
#define VDST_MAXTIMINGLEVELS 10    // maximum number of timing levels
//////////////////////////////////////////////////////////////////

using namespace std;

class VGlobalRunParameter
{
    private:
    
        static bool      bReadRunParameter;
        static bool      bDebug;                                         // print debug output
        
        static string       fObservatory;
        static double       fObservatory_Latitude_deg;
        static double       fObservatory_Longitude_deg;
        static double       fObservatory_Height_m;
        
        
        // OUTPUT TREE VERSION
        //
        // changes from 5 to 6: LTrig  now ULong64_t
        //                      ImgSel now ULong64_t
        // changes from 6 to 7: introduced list of selected telescopes
        //
        // changes from 7 to 8: add MC primary to showerpars
        static unsigned int fEVNDISP_TREE_VERSION;
        
        static string    fDBServer;                                         // database location (VTS)
        static string    fRawDataServer;                                    // location of raw data (VTS)
        
        // DIRECTORIES
        static string fEVNDISPAnaDataDirectory;          // directory where all data (detectorgeometry, ...) is expected and written to (output file)
        static string fVBFRawDataDirectory;              // directory with VERITAS vbf data (vbf files)
        static string fEVNDISPCalibrationDataDirectory;  // directory where calibration data is expected and written to
        static string fEVNDISPOutputDirectory;           // output- and result files are written into this directory
        
    public:
    
        static string       fEVNDISP_VERSION;                             // EVNDISPLAY VERSION
        
        VGlobalRunParameter( bool bSetGlobalParameter = true );
        virtual ~VGlobalRunParameter();
        
        string       getDBServer() const
        {
            return fDBServer;
        }
        static string getDirectory_EVNDISPAnaData()
        {
            return fEVNDISPAnaDataDirectory;
        }
        string       getDirectory_EVNDISPCalibrationData()
        {
            return fEVNDISPAnaDataDirectory + "/Calibration/";
        }
        string       getDirectory_EVNDISPCalibrationData_perRun()
        {
            return fEVNDISPCalibrationDataDirectory + "/Calibration/";
        }
        string       getDirectory_EVNDISPDetectorGeometry()
        {
            return fEVNDISPAnaDataDirectory + "/DetectorGeometry/";
        }
        string       getDirectory_EVNDISPParameterFiles()
        {
            return fEVNDISPAnaDataDirectory + "/ParameterFiles/";
        }
        string       getDirectory_VBFRawData()
        {
            return fVBFRawDataDirectory;
        }
        string       getDirectory_EVNDISPOutput()
        {
            return fEVNDISPOutputDirectory;
        }
        static string getEVNDISP_VERSION()
        {
            return fEVNDISP_VERSION;
        }
        static float  getEVNDISP_VERSION_FL()
        {
            return atof( fEVNDISP_VERSION.substr( 2, 5 ).c_str() );
        }
        static unsigned int getEVNDISP_VERSION_UI()
        {
            return ( unsigned int )( 100.*getEVNDISP_VERSION_FL() );
        }
        static unsigned int getEVNDISP_TREE_VERSION()
        {
            return fEVNDISP_TREE_VERSION;
        }
        static unsigned int getEVNDISP_TREE_VERSION( TTree* );
        static bool         getEVNDISP_TREE_isShort( TTree* );
        string       getObservatory()
        {
            return fObservatory;
        }
        double       getObservatory_Height_m()
        {
            return fObservatory_Height_m;
        }
        static double       getObservatory_Latitude_deg()
        {
            return fObservatory_Latitude_deg;
        }
        static double       getObservatory_Longitude_deg()
        {
            return fObservatory_Longitude_deg;
        }
        string       getRawDataServer() const
        {
            return fRawDataServer;
        }
        void         printGlobalRunParameter();
        bool         readRunparameterFile( string iFile );
        bool         setDirectories();
        void         setDirectory_EVNDISPOutput( string iDir )
        {
            fEVNDISPOutputDirectory = iDir;
        }
        bool         setDirectory_EVNDISPCalibrationData( string iDir );
        bool         update( TChain* ic );
        
        ClassDef( VGlobalRunParameter, 10 );
};

#endif
