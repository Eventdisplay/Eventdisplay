/*! \class VReadRunParameter
    \brief read command line parameters and do some tests for consistency

     evndisp -help gives an overview of the most important parameters
     all parameters are listed in README/README.EVNDISP.commandline

*/

#include <VReadRunParameter.h>

VReadRunParameter::VReadRunParameter()
{
    fDebug = false;
    if( fDebug )
    {
        cout << "VReadRunParameter::VReadRunParameter()" << endl;
    }
    fRunPara = new VEvndispRunParameter();
    fRunPara->printGlobalRunParameter();
    fusercamera = false;
    f_boolCommandline = false;
    f_boolConfigfile = false;
    
    fTelToAna = 0;
    fTelToAnaString = "";
    
    fPrintOutputFile = false;
    
    fTargetName_overwriteDB = "";
    fRA_overwriteDB = -9999.;
    fDec_overwriteDB = -9999.;
    fWobbleNorth_overwriteDB = -9999.;
    fWobbleEast_overwriteDB  = -9999.;
}

VReadRunParameter::~VReadRunParameter()
{

    if( fRunPara )
    {
        delete fRunPara;
    }
    
}


/*!
   \attention
      reading from the command line assumes same parameters for all telescopes (preliminary)
*/
bool VReadRunParameter::readCommandline( int argc, char* argv[] )
{
    if( fDebug )
    {
        cout << "VReadRunParameter::readCommandline()" << endl;
    }
    f_boolCommandline = true;
    int i = 1;
    // no command line parameters
    if( argc == 1 )
    {
        printHelp();
        return false;
    }
    // read all command line parameters
    while( i++ < argc )
    {
        string iTemp  = argv[i - 1];
        string iTemp1 = argv[i - 1];              // this is to get the camera name correctly
        string iTemp2 = "";
        if( i < argc )
        {
            iTemp2 = argv[i];
        }
        iTemp = VUtilities::lowerCase( iTemp );
        // print help text
        if( iTemp.find( "help" ) < iTemp.size() )
        {
            printHelp();
            return false;
        }
        if( iTemp.find( "debug" ) < iTemp.size() )
        {
            fRunPara->fDebug = true;
            fDebug = true;
        }
        // nice plots for talks/papers
        else if( iTemp.find( "plotpape" ) != string::npos )
        {
            fRunPara->fPlotPaper = true;
        }
        // camera
        else if( iTemp.find( "camera" ) < iTemp.size() && !( iTemp.find( "cameradirectory" ) < iTemp.size() )
                 && !( iTemp.find( "cameracoordinatetransformx" ) < iTemp.size() )  && !( iTemp.find( "cameracoordinatetransformy" ) < iTemp.size() ) )
        {
            // reading iTemp1, to get upper/lower cases right
            fRunPara->fcamera[0] = iTemp1.substr( iTemp1.rfind( "=" ) + 1, iTemp1.size() );
            fusercamera = true;
        }
        else if( iTemp.find( "vbfnsamples" ) < iTemp.size() )
        {
            fRunPara->fUseVBFSampleLength = true;
        }
        // number of telescopes (>0)
        else if( iTemp.find( "ntel" ) < iTemp.size() )
        {
            fRunPara->fNTelescopes = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fRunPara->fNTelescopes <= 0 )
            {
                cout << "can only work with more than one telescope" << endl;
                return false;
            }
        }
        // display mode
        else if( iTemp.find( "disp" ) < iTemp.size() )
        {
            fRunPara->fdisplaymode = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // muon mode
        else if( iTemp.find( "muon" ) < iTemp.size() )
        {
            fRunPara->fmuonmode = 1;
        }
        else if( iTemp.find( "hough" ) < iTemp.size() )
        {
            fRunPara->fhoughmuonmode = 1;
        }
        /////////////////////////////////////////////////////////
        // Frogs mode
        else if( iTemp1 == "-frogs" )
        {
            fRunPara->ffrogsmode = 1;
            if( iTemp2.size() > 0 )
            {
                fRunPara->ffrogsmscwfile = iTemp2;
                i++;
                // setting default parameters for frogs
                // (might be overwritten by the following command lines)
                if( fRunPara->ffrogsparameterfile.size() == 0 )
                {
                    fRunPara->ffrogsparameterfile = "FROGS.runparameter";
                }
                if( fRunPara->ffrogsRecID < 0 )
                {
                    fRunPara->ffrogsRecID = 0;
                }
            }
            else
            {
                fRunPara->ffrogsmscwfile = "";
            }
        }
        else if( iTemp1 == "-frogsid" )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->ffrogsRecID = atoi( iTemp2.c_str() );
                i++;
            }
            else
            {
                fRunPara->ffrogsRecID = 0;
            }
        }
        else if( iTemp1 == "-templatelistforfrogs" )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->ffrogstemplatelist = iTemp2 ;
                i++;
            }
            else
            {
                fRunPara->ffrogstemplatelist = "" ;
            }
        }
        else if( iTemp1 == "-frogsparameterfile" )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->ffrogsparameterfile = iTemp2 ;
                i++;
            }
            else
            {
                fRunPara->ffrogsparameterfile = "" ;
            }
        }
        /////////////////////////////////////////////////////////
        // Model3D
        else if( iTemp.find( "model3d" ) < iTemp.size() )
        {
            fRunPara->fUseModel3D = true;
            fRunPara->fLnLTableFile = "table_LnL.root";
        }
        else if( iTemp.find( "plot3d" ) < iTemp.size() )
        {
            fRunPara->fUseModel3D = true;
            fRunPara->fUseDisplayModel3D = true;
            fRunPara->fLnLTableFile = "table_LnL.root";
            fRunPara->fdisplaymode = 1;
        }
        else if( iTemp.find( "lnlfile" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fLnLTableFile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fLnLTableFile = "table_LnL.root";
            }
        }
        else if( iTemp.find( "createlnltablefile" ) < iTemp.size() )
        {
            fRunPara->fUseModel3D = true;
            fRunPara->fCreateLnLTable = true;
        }
        // source file
        else if( iTemp.find( "sourcefi" ) < iTemp.size() )
        {
            checkSecondArgument( iTemp1, iTemp2, true );
            if( iTemp2.size() > 0 )
            {
                fRunPara->fsourcefile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fsourcefile = "";
            }
        }
        // pedestal file for grisu simulations
        else if( iTemp.find( "pedestalfile" ) < iTemp.size() )
        {
            checkSecondArgument( iTemp1, iTemp2, true );
            if( iTemp2.size() > 0 )
            {
                fRunPara->fsimu_pedestalfile = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "combine_pedestal_channels" ) < iTemp.size() )
        {
            if( iTemp.find( "=" ) != string::npos )
            {
                fRunPara->fCombineChannelsForPedestalCalculation = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            }
            else
            {
                fRunPara->fCombineChannelsForPedestalCalculation = 1;
            }
        }
        // for pedestal calculation: write results into a single root file
        else if( iTemp.find( "singlepedestalrootfile" ) < iTemp.size() )
        {
            fRunPara->fPedestalSingleRootFile = ( bool )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // user name
        else if( iTemp.find( "user" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fEventDisplayUser = iTemp2;
                i++;
            }
        }
        // file with dead channel definition
        else if( iTemp.find( "deadchannelfile" ) < iTemp.size() )
        {
            checkSecondArgument( iTemp1, iTemp2, true );
            if( iTemp2.size() > 0 )
            {
                fRunPara->fDeadChannelFile = iTemp2;
                i++;
            }
        }
        // hilo from MC run her
        else if( iTemp.find( "simu_hilo_from_simfile" ) < iTemp.size() )
        {
            fRunPara->fsimu_HILO_from_simFile = true;
        }
        // noise level of external grisu pedestal file
        else if( iTemp.find( "pedestalnoiselevel" ) < iTemp.size() )
        {
            fRunPara->fsimu_noiselevel = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "pedestaldefaultpedestal" ) < iTemp.size() )
        {
            fRunPara->fsimu_pedestalfile_DefaultPed = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "lowgainpedestallevel" ) < iTemp.size() )
        {
            fRunPara->fsimu_lowgain_pedestal_DefaultPed = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // calibration file
        else if( iTemp.find( "calibrationfi" ) < iTemp.size() && !( iTemp.find( "lowgain" ) < iTemp.size() ) )
        {
            checkSecondArgument( iTemp1, iTemp2, true );
            if( iTemp2.size() > 0 )
            {
                fRunPara->fcalibrationfile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fcalibrationfile = "";
            }
            // no reading of DB in case of external calibration file
            fRunPara->freadCalibfromDB = false;
            
        }
        else if( iTemp.find( "readcalibdb" ) < iTemp.size() )
        {
            fRunPara->freadCalibfromDB = true;
            
            if( iTemp2.size() > 0 && VUtilities::isInteger( iTemp2 ) )
            {
                fRunPara->freadCalibfromDB_versionquery = atoi( iTemp2.c_str() );
                i++;
            }
            else
            {
                fRunPara->freadCalibfromDB_versionquery = -111;
            }
            
        }
        else if( iTemp.find( "readandsavecalibdb" ) < iTemp.size() )
        {
            fRunPara->freadCalibfromDB_save_file = true;
            
            fRunPara->freadCalibfromDB = true;
            
            if( iTemp2.size() > 0 )
            {
                fRunPara->freadCalibfromDB_versionquery = atoi( iTemp2.c_str() );
                i++;
            }
            else
            {
                fRunPara->freadCalibfromDB_versionquery = -111;
            }
            
        }
        else if( iTemp.find( "nocalibnoproblem" ) < iTemp.size() )
        {
        
            if( !fRunPara->freadCalibfromDB )
            {
                fRunPara->fNoCalibNoPb = true;    //Security needed for the step in VImageBaseAnalyzer::findDeadChans, where the Gains and TOffsets are not checked if fNoCalibNoPb = true. This should not happen when reading information from the VOFFLine DB (LG)
            }
            
        }
        else if( iTemp.find( "nextdaygainhack" ) < iTemp.size() )
        {
            fRunPara->fNextDayGainHack = true;
        }
        else if( iTemp.find( "ignoredstgains" ) < iTemp.size() )
        {
            fRunPara->fIgnoreDSTGains = true;
        }
        else if( iTemp.find( "lowgaincalibrationfile" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fLowGainCalibrationFile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fLowGainCalibrationFile = "";
            }
        }
        else if( iTemp.find( "arraycuts" ) < iTemp.size() || iTemp.find( "recopara" ) < iTemp.size() || iTemp.find( "reconstructionparameter" ) < iTemp.size() )
        {
            checkSecondArgument( iTemp1, iTemp2, true );
            if( iTemp2.size() > 0 )
            {
                fRunPara->freconstructionparameterfile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->freconstructionparameterfile = "";
            }
        }
        // dst output file
        else if( iTemp.find( "dstfile" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fdstfile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fdstfile = "dstfile.root";
            }
        }
        else if( iTemp.find( "dstcalibration" ) < iTemp.size() )
        {
            fRunPara->fdstcalibration = ( bool )( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        // star catalogue
        else if( iTemp.find( "starcatalogue" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fStarCatalogueName = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "starbrightness" ) < iTemp.size() )
        {
            fRunPara->fMinStarBrightness_B = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // write all pixel to dst file
        else if( iTemp.find( "dstallpixel" ) < iTemp.size() )
        {
            fRunPara->fdstwriteallpixel  = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // minimal number of tubes for dst event
        else if( iTemp.find( "dstntubes" ) < iTemp.size() )
        {
            fRunPara->fdstminntubes = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // calculate pedestal variations on a short time scale
        else if( iTemp.find( "pedestalsintimeslices" ) < iTemp.size() && !( iTemp.find( "pedestalsintimeslicessumwindow" ) < iTemp.size() )
                 && !( iTemp.find( "pedestalsintimeslicessumfirst" ) < iTemp.size() )
                 && !( iTemp.find( "usepedestalsintimeslices" ) < iTemp.size() )
                 && !( iTemp.find( "usepedestalsintimesliceslowgain" ) < iTemp.size() )
                 && !( iTemp.find( "nopedestalsintimeslices" ) < iTemp.size() ) )
        {
            fRunPara->fPedestalsInTimeSlices = true;
            fRunPara->fUsePedestalsInTimeSlices = true;
            fRunPara->fLowGainUsePedestalsInTimeSlices = true;
        }
        // inject Gaussian noise
        else if( iTemp.find( "injectgaussiannoise" ) < iTemp.size() && !( iTemp.find( "injectgaussiannoiseseed" ) < iTemp.size() ) )
        {
            fRunPara->finjectGaussianNoise = double( atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        else if( iTemp.find( "injectgaussiannoiseseed" ) < iTemp.size() )
        {
            fRunPara->finjectGaussianNoiseSeed = UInt_t( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        else if( iTemp.find( "nopedestalsintimeslices" ) < iTemp.size() )
        {
            fRunPara->fPedestalsInTimeSlices = false;
        }
        else if( iTemp.find( "usepedestalsintimeslices" ) < iTemp.size() && !( iTemp.find( "usepedestalsintimesliceslowgain" ) < iTemp.size() ) )
        {
            fRunPara->fUsePedestalsInTimeSlices = bool( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        else if( iTemp.find( "usepedestalsintimesliceslowgain" ) < iTemp.size() )
        {
            fRunPara->fLowGainUsePedestalsInTimeSlices = bool( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        // parameter for pedestal variations on a short time scale for tracking tests
        else if( iTemp.find( "pedestalslengthoftimeslice" ) < iTemp.size() )
        {
            fRunPara->fPedestalsLengthOfTimeSlice = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            fRunPara->fUsePedestalsInTimeSlices = true;
            fRunPara->fLowGainUsePedestalsInTimeSlices = true;
        }
        else if( iTemp.find( "sumwindowaveragetime" ) < iTemp.size() )
        {
            fRunPara->fCalibrationSumWindowAverageTime = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "calibrationsumwindow" ) < iTemp.size() && iTemp != "pedestalsintimeslices" )
        {
            fRunPara->fCalibrationSumWindow = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "calibrationsumfirst" ) < iTemp.size() && iTemp != "pedestalsintimeslices" )
        {
            fRunPara->fCalibrationSumFirst = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "calibrationsummin" ) < iTemp.size() && iTemp != "pedestalsintimeslices" )
        {
            fRunPara->fCalibrationIntSumMin = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // fast plotting without pedestals
        else if( iTemp.find( "plotraw" ) < iTemp.size() )
        {
            fRunPara->fPlotRaw = true;
            fRunPara->fuseDB = false;
            fRunPara->fDBTracking = false;
            fRunPara->fTargetName = "laser";
            fRunPara->fdisplaymode = 1;
        }
        else if( iTemp.find( "plotmethod" ) < iTemp.size() )
        {
            fRunPara->fPlotAllInOneMethod = ( unsigned int )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "writeallmc" ) < iTemp.size() )
        {
            fRunPara->fWriteTriggerOnly = false;
        }
        else if( iTemp.find( "writenomctree" ) < iTemp.size() )
        {
            fRunPara->fwriteMCtree = 0;
        }
        else if( iTemp.find( "writemctree" ) < iTemp.size() )
        {
            fRunPara->fwriteMCtree = 1;
        }
        else if( iTemp.find( "fillmchisto" ) < iTemp.size() )
        {
            fRunPara->fFillMCHistos = true;
        }
        // use db for run infos
        else if( iTemp.find( "usedbinfo" ) < iTemp.size() && !( iTemp.find( "donotusedbinfo" ) < iTemp.size() ) )
        {
            fRunPara->fuseDB = true;
            isCompiledWithDB();
        }
        else if( iTemp.find( "donotusedbinfo" ) < iTemp.size() )
        {
            fRunPara->fuseDB = false;
        }
        else if( iTemp.find( "-savedeadpixelregistry" ) < iTemp.size() )
        {
            fRunPara->fSaveDeadPixelRegistry = true ;
            cout << "fSaveDeadPixelRegistry = true" << endl;
        }
        // ignore configuration file versions
        else if( iTemp.find( "ignorecfgversions" ) < iTemp.size() )
        {
            fRunPara->fIgnoreCFGversions = true;
        }
        // print analysis progress
        else if( iTemp.find( "printanalysisprogress" ) < iTemp.size() || iTemp.find( "pap" ) < iTemp.size() )
        {
            fRunPara->fPrintAnalysisProgress = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // set run number from command line (default: get it from source file name)
        else if( iTemp.rfind( "runnum" ) < iTemp.size() || iTemp == "run" )
        {
            fRunPara->frunnumber = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "runtitle" ) < iTemp.size() )
        {
            fRunPara->fRunTitle = iTemp.substr( iTemp.find( "=" ) + 1, iTemp.size() );
        }
        else if( iTemp.rfind( "epochfile" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fEpochFile = iTemp2;
                i++;
            }
        }
        else if( iTemp.rfind( "epoch" ) < iTemp.size() && !( iTemp.find( "epochfile" ) < iTemp.size() ) )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fInstrumentEpoch = iTemp2;
                i++;
            }
        }
        else if( iTemp.rfind( "atmosphereid" ) < iTemp.size() )
        {
            fRunPara->fAtmosphereID = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.rfind( "gaincorrection" ) < iTemp.size() )
        {
            fRunPara->fGainCorrection[0] = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // which telescope should be analyzed
        else if( iTemp.rfind( "teltoana" ) < iTemp.size() )
        {
            fTelToAnaString = iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() );
            // comma separated list
            if( fTelToAnaString.find( "," ) < fTelToAnaString.size() )
            {
                fTelToAna = 0;
            }
            else
            {
                fTelToAna = ( unsigned int )( atoi( fTelToAnaString.c_str() ) );
            }
        }
        // elevation
        else if( iTemp.rfind( "elevat" ) < iTemp.size() )
        {
            fRunPara->felevation = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // azimuth
        else if( iTemp.rfind( "azimuth" ) < iTemp.size() )
        {
            fRunPara->fazimuth = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // target declination
        else if( iTemp.rfind( "declination" ) < iTemp.size() )
        {
            fRunPara->fTargetDec = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            fDec_overwriteDB = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // target ra [J2000]
        else if( iTemp.rfind( "rightascension" ) < iTemp.size() )
        {
            fRunPara->fTargetRA = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            fRA_overwriteDB = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // target declination offset [J2000]
        else if( iTemp.rfind( "decoffset" ) < iTemp.size() )
        {
            fRunPara->fTargetDecOffset = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // target ra offset
        else if( iTemp.rfind( "raoffset" ) < iTemp.size() )
        {
            fRunPara->fTargetRAOffset = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // target name
        else if( iTemp.rfind( "target" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fTargetName = iTemp2;
                fTargetName_overwriteDB = iTemp2;
                i++;
            }
            else
            {
                fRunPara->fTargetName = "";
            }
        }
        // these two command line settings might be overwritten by values read from the data base
        // wobble offset NORTH
        else if( iTemp.rfind( "wobblenorth" ) < iTemp.size() && !( iTemp.rfind( "overwritedb_wobblenorth" ) < iTemp.size() ) )
        {
            fRunPara->fWobbleNorth = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // wobble offset EAST
        else if( iTemp.rfind( "wobbleeast" ) < iTemp.size() && !( iTemp.rfind( "overwritedb_wobbleeast" ) < iTemp.size() ) )
        {
            fRunPara->fWobbleEast = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // these two command line settings are never overwritten by values read from the data base
        else if( iTemp.rfind( "overwritedb_wobblenorth" ) < iTemp.size() )
        {
            fWobbleNorth_overwriteDB = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // wobble offset EAST
        else if( iTemp.rfind( "overwritedb_wobbleeast" ) < iTemp.size() )
        {
            fWobbleEast_overwriteDB = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // pointing check
        else if( iTemp.rfind( "checkpointing" ) < iTemp.size() )
        {
            fRunPara->fCheckPointing  = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // pointing error x
        else if( iTemp.rfind( "-pointingerrorx" ) < iTemp.size() )
        {
            unsigned int iTelID = 0;
            double ix = 0.;
            iTelID = ( unsigned int ) atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.rfind( ":" ) ).c_str() ) - 1;
            ix = atof( iTemp.substr( iTemp.rfind( ":" ) + 1, iTemp.size() ).c_str() );
            f_pointingErrorX[iTelID] = ix;
            fRunPara->fDBTracking = false;
        }
        // pointing error y
        else if( iTemp.rfind( "-pointingerrory" ) < iTemp.size() )
        {
            unsigned int iTelID = 0;
            double ix = 0.;
            iTelID = ( unsigned int ) atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.rfind( ":" ) ).c_str() ) - 1;
            ix = atof( iTemp.substr( iTemp.rfind( ":" ) + 1, iTemp.size() ).c_str() );
            f_pointingErrorY[iTelID] = ix;
            fRunPara->fDBTracking = false;
        }
        else if( iTemp.rfind( "-usedbtracking" ) < iTemp.size() )
        {
            fRunPara->fDBTracking = true;
            isCompiledWithDB();
        }
        else if( iTemp.rfind( "-usenodbtracking" ) < iTemp.size() )
        {
            fRunPara->fDBTracking = false;
        }
        else if( iTemp.rfind( "-usedbvpm" ) < iTemp.size() )
        {
            fRunPara->fDBVPM = true;
            fRunPara->fDBUncalibratedVPM = false;
            fRunPara->fDBTracking = true;
        }
        else if( iTemp.rfind( "-usenodbvpm" ) < iTemp.size() )
        {
            fRunPara->fDBVPM = false;
            fRunPara->fDBUncalibratedVPM = false;
            fRunPara->fDBTracking = true;
        }
        else if( iTemp.rfind( "-useuncalibratedvpm" ) < iTemp.size() )
        {
            fRunPara->fDBVPM = false;
            fRunPara->fDBUncalibratedVPM = true;
            fRunPara->fDBTracking = true;
        }
        else if( iTemp.rfind( "-usedbrotations" ) < iTemp.size() )
        {
            fRunPara->fDBCameraRotationMeasurements = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.rfind( "-usetcorrectionfrom" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fDBTrackingCorrections = iTemp2;
                fRunPara->fDBTracking = true;
                i++;
            }
            else
            {
                fRunPara->fDBTrackingCorrections = "";
                cout << "no date for T-Point correction given" << endl;
                return false;
            }
        }
        else if( iTemp.rfind( "-pointingmonitortxt" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fPMTextFileDirectory = iTemp2;
                fRunPara->fDBTracking = true;
                i++;
            }
            else
            {
                fRunPara->fPMTextFileDirectory = "";
                cout << "no pointing monitor text file directory give" << endl;
                return false;
            }
        }
        // min laser charge
        else if( iTemp.rfind( "lasermin" ) < iTemp.size() )
        {
            fRunPara->fLaserSumMin = ( double )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // run mode
        else if( iTemp.rfind( "mod" ) < iTemp.size() )
        {
            fRunPara->frunmode = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fRunPara->frunmode < 0 || fRunPara->frunmode > 8 )
            {
                cout << "Error invalid run mode (0=analysis,1=pedestal,2=gains/toffsets,3=trace library,4=dst output,";
                cout << "5=gains/toffsets low gain,6=low gain pedestals,7=tzero,8=low gain tzero)" << endl;
                return false;
            }
        }
        // number of events
        else if( iTemp.rfind( "neve" ) < iTemp.size() && !( iTemp.find( "calibrationnevents" ) < iTemp.size() ) )
        {
            fRunPara->fnevents = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fRunPara->fnevents <= 0 )
            {
                fRunPara->fnevents = -10000;
            }
        }
        else if( iTemp.find( "calibrationnevents" ) < iTemp.size() )
        {
            fRunPara->fNCalibrationEvents  = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // first event number (skip to this point)
        else if( iTemp.rfind( "firstevent" ) < iTemp.size() )
        {
            fRunPara->fFirstEvent = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fRunPara->fFirstEvent < 0 )
            {
                fRunPara->fFirstEvent = -10000;
            }
        }
        // start analyzing at this minute
        else if( iTemp.find( "timecutmin" ) < iTemp.size() )
        {
            fRunPara->fTimeCutsMin_min = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // stopp analyzing at this minute
        else if( iTemp.find( "timecutmax" ) < iTemp.size() )
        {
            fRunPara->fTimeCutsMin_max = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        
        // check if the user wants to print the list of dead pixels for this run
        else if( iTemp.rfind( "printdeadpixelinfo" ) < iTemp.size() ) // DEADCHAN
        {
            fRunPara->fprintdeadpixelinfo = true ;
        }
        
        
        // source type
        else if( iTemp.find( "type" ) < iTemp.size() )
        {
            fRunPara->fsourcetype = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fRunPara->fsourcetype < 0 || fRunPara->fsourcetype > 7 )
            {
                cout << "unknown sourcetype " << fRunPara->fsourcetype << endl;
                return false;
            }
            // MC grisu file
            if( fRunPara->fsourcetype == 1 || fRunPara->fsourcetype == 2
                    || fRunPara->fsourcetype == 5 || fRunPara->fsourcetype == 6
                    || fRunPara->fsourcetype == 7 )
            {
                fRunPara->fIsMC = 1;
            }
        }
        else if( iTemp.find( "tracewindowshift" ) < iTemp.size() )
        {
            fRunPara->fTraceWindowShift[0] = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "mc_fadctracestart" ) < iTemp.size() )
        {
            fRunPara->fMC_FADCTraceStart = ( unsigned int )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "averagetzerofiducialradius" ) < iTemp.size() )
        {
            fRunPara->faverageTZeroFiducialRadius = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "cameracoordinatetransformx" ) < iTemp.size() )
        {
            fRunPara->fCameraCoordinateTransformX = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "cameracoordinatetransformy" ) < iTemp.size() )
        {
            fRunPara->fCameraCoordinateTransformY = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "loop" ) < iTemp.size() )
        {
            fRunPara->floopmode = true;
        }
        else if( iTemp.find( "printoutputfile" ) < iTemp.size() )
        {
            fPrintOutputFile = true;
        }
        else if( iTemp.find( "outputdirectory" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->setDirectory_EVNDISPOutput( iTemp2 );
                i++;
            }
        }
        else if( iTemp.find( "calibrationdirectory" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                if( !fRunPara->setDirectory_EVNDISPCalibrationData( iTemp2 ) )
                {
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
                i++;
            }
        }
        // analysis output file
        else if( iTemp.find( "output" ) < iTemp.size() && iTemp != "outputdirectory" && iTemp != "printoutputfile" )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->foutputfileName = iTemp2;
                i++;
            }
            else
            {
                fRunPara->foutputfileName = "";
            }
        }
        else if( iTemp.find( "highres" ) < iTemp.size() )
        {
            fRunPara->fh = ( unsigned int )( 833 );
            fRunPara->fw = ( unsigned int )( 1250 );
        }
        else if( iTemp.find( "hdhres" ) < iTemp.size() )
        {
            fRunPara->fh = ( unsigned int )( 833 * 1.5 );
            fRunPara->fw = ( unsigned int )( 1250 * 1.5 );
        }
        else if( iTemp.find( "telenoff" ) < iTemp.size() )
        {
            fRunPara->ftelescopeNOffset = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "sampleoff" ) < iTemp.size() )
        {
            fRunPara->fsampleoffset = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "logl" ) < iTemp.size() && !( iTemp.find( "loglminloss" ) < iTemp.size() ) )
        {
            fRunPara->fImageLL = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "loglminloss" ) < iTemp.size() )
        {
            fRunPara->fLogLikelihoodLoss_min[0] = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "fuifactor" ) < iTemp.size() )
        {
            fRunPara->fImageAnalysisFUIFactor  = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "fixwindow2" ) < iTemp.size() )
        {
            fRunPara->fFixWindowStart_sumwindow2 = ( bool )( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        else if( iTemp.find( "recoverimagesneardeadpixel" ) < iTemp.size() )
        {
            fRunPara->frecoverImagePixelNearDeadPixel = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            i++;
        }
        else if( iTemp.find( "smoothdead" ) < iTemp.size() )
        {
            fRunPara->fSmoothDead = true;
        }
        else if( iTemp.find( "usepeds" ) < iTemp.size() && !( iTemp.find( "donotusepeds" ) < iTemp.size() ) )
        {
            fRunPara->fUsePedEvents = true;
        }
        else if( iTemp.find( "donotusepeds" ) < iTemp.size() )
        {
            fRunPara->fUsePedEvents = false;
        }
        else if( iTemp.find( "l2timecorrect" ) < iTemp.size() )
        {
            fRunPara->fL2TimeCorrect = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "l2setspecialchannels" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->fsetSpecialChannels = iTemp2;
                if( fRunPara->fsetSpecialChannels == "nofile" )
                {
                    fRunPara->fsetSpecialChannels = "";
                }
                i++;
            }
        }
        else if( iTemp.find( "tracelib" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunPara->ftracefile = iTemp2;
                i++;
            }
            else
            {
                fRunPara->ftracefile = "";
            }
        }
        else if( iTemp.find( "ignoredead" ) < iTemp.size() )
        {
            fRunPara->fMCnoDead = true;
        }
        else if( iTemp.find( "ndeadchannel" ) < iTemp.size() )
        {
            fRunPara->fMCNdead = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "seedndead" ) < iTemp.size() )
        {
            fRunPara->fMCNdeadSeed = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "ndeadboard" ) < iTemp.size() )
        {
            fRunPara->fMCNdeadboard = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "writelaserpulsen" ) < iTemp.size() )
        {
            fRunPara->fwriteLaserPulseN = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "writeaveragelaserpulse" ) < iTemp.size() )
        {
            fRunPara->fwriteAverageLaserPulse = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "fadcscale" ) < iTemp.size() )
        {
            fRunPara->fMCScale = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "pedestalseed" ) < iTemp.size() )
        {
            fRunPara->fgrisuseed = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "ismc" ) < iTemp.size() )
        {
            fRunPara->fIsMC = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "shorttree" ) < iTemp.size() && !( iTemp.find( "noshorttree" ) < iTemp.size() ) )
        {
            fRunPara->fShortTree = 1;
        }
        else if( iTemp.find( "noshorttree" ) < iTemp.size() )
        {
            fRunPara->fShortTree = 0;
        }
        else if( iTemp.find( "pwmethod" ) < iTemp.size() )
        {
            fRunPara->fPWmethod = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "pwcleanneighbors" ) < iTemp.size() )
        {
            fRunPara->fPWcleanNeighbors = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "pwcleanthreshold" ) < iTemp.size() )
        {
            fRunPara->fPWcleanThreshold = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "pwlimit" ) < iTemp.size() )
        {
            fRunPara->fPWlimit = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.rfind( "printgrisuheader" ) < iTemp.size() )
        {
            fRunPara->fPrintGrisuHeader = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.rfind( "writeextracalibtree" ) < iTemp.size() )
        {
            fRunPara->fWriteExtraCalibTree = true;
        }
        else if( i > 1 )
        {
            cout << "unknown command line parameter: " << iTemp << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // test and adjust some parameters
    test_and_adjustParams();
    
    return true;
}


/*!

     test command line parameters given, compare them is requreseted with the DB results
     fill values for all telescopes

*/
void VReadRunParameter::test_and_adjustParams()
{
    if( fDebug )
    {
        cout << "VReadRunParameter::test_and_adjustParams()" << endl;
    }
    
    // test if this is a DST file
    getRunParametersFromDST();
    // standard VTS analysis
    if( (fRunPara->frunmode == 0 || fRunPara->frunmode == 4 )
      && fRunPara->getObservatory().find( "VERITAS" ) != string::npos
            && fRunPara->fIsMC == 0 )
    {
        if( !fRunPara->fNoCalibNoPb && fRunPara->fcalibrationfile == "" )
        {
            fRunPara->freadCalibfromDB = true;
        }
    }
    // set pulse timing for DST case
    if( fRunPara->frunmode == 4 && fRunPara->getObservatory().find( "VERITAS" ) == string::npos )
    {
        fRunPara->fpulsetiminglevels.clear();
        fRunPara->fpulsetiminglevels.push_back( 0.2 );
        fRunPara->fpulsetiminglevels.push_back( 0.5 );
        fRunPara->fpulsetiminglevels.push_back( 0.8 );
        fRunPara->fpulsetiminglevels.push_back( 1.0 );
        unsigned int i_fps = fRunPara->fpulsetiminglevels.size();
        if( i_fps > 1 )
        {
            for( unsigned int i = 0; i < i_fps - 1; i++ )
            {
                fRunPara->fpulsetiminglevels.push_back( fRunPara->fpulsetiminglevels[i_fps - i - 2] );
            }
        }
        // calculate tzero index and width index
        fRunPara->setPulseZeroIndex();
    }
    
    fRunPara->setSystemParameters();
    
    /////////////////////////////////////////////////////////////////
    // MC adjustments
    if( fRunPara->fIsMC > 0 )
    {
        fRunPara->fDBTracking = false;
        fRunPara->fuseDB = false;
        fRunPara->fL2TimeCorrect = false;
        fRunPara->fDBCameraRotationMeasurements = false;
        // this is for example required for Grisudet sims, when low gains are read from
        // the cfg file
        if( fRunPara->fLowGainCalibrationFile == "NOFILE" || ( fRunPara->fsourcetype == 4 || fRunPara->fsourcetype == 7 ) )
        {
            fRunPara->fLowGainCalibrationFile = "";
        }
        if( !fRunPara->fUseModel3D )  //gains needed for Model3D
        {
            fRunPara->fEpochFile = "";
        }
    }
    
    // CTA/AGIS adjustments
    if( fRunPara->getObservatory().find( "cta" ) != string::npos || fRunPara->getObservatory().find( "CTA" ) != string::npos
            || fRunPara->getObservatory().find( "agis" ) != string::npos || fRunPara->getObservatory().find( "AGIS" ) != string::npos )
    {
        // no special channels allowed (e.g. L2 timing channels)
        fRunPara->fsetSpecialChannels = "";
        // no dead channels allowed
        fRunPara->fDeadChannelFile = "";
        // no epochs (e.g. summer and winter) for CTA
        fRunPara->fEpochFile = "";
    }
    
    ////////////////////////////////////////////////////
    // get runnumbers from file name (if source file is given)
    if( fRunPara->frunnumber <= 0 )
    {
        if( fRunPara->fsourcefile.find( "pure" ) < fRunPara->fsourcefile.size() )
        {
            fRunPara->frunnumber = atoi( fRunPara->fsourcefile.substr( fRunPara->fsourcefile.rfind( "/" ) + 1, fRunPara->fsourcefile.rfind( "." ) ).c_str() );
        }
        // MC (grisu)
        else if( fRunPara->fsourcefile.find( "grisu" ) < fRunPara->fsourcefile.size() )
        {
            fRunPara->frunnumber = atoi( fRunPara->fsourcefile.substr( fRunPara->fsourcefile.rfind( "/" ) + 4, 6 ).c_str() );
        }
        // vbf
        else if( fRunPara->fsourcefile.find( "vbf" ) < fRunPara->fsourcefile.size() || fRunPara->fsourcefile.find( "cvbf" ) < fRunPara->fsourcefile.size() )
        {
            fRunPara->frunnumber = atoi( fRunPara->fsourcefile.substr( fRunPara->fsourcefile.rfind( "/" ) + 1, fRunPara->fsourcefile.rfind( "." ) ).c_str() );
            if( fRunPara->fsourcetype == 0 )
            {
                fRunPara->fsourcetype = 3;
            }
        }
        // dst
        else if( fRunPara->fsourcefile.find( "dst" ) < fRunPara->fsourcefile.size() )
        {
            fRunPara->frunnumber = atoi( fRunPara->fsourcefile.substr( fRunPara->fsourcefile.rfind( "/" ) + 4, 6 ).c_str() );
        }
        else
        {
            fRunPara->frunnumber = atoi( fRunPara->fsourcefile.substr( fRunPara->fsourcefile.rfind( "/" ) + 1, fRunPara->fsourcefile.size() ).c_str() );
        }
    }
    /////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    // settings for pedestal and gain calculations
    if( fRunPara->frunmode == 1 || fRunPara->frunmode == 2 || fRunPara->frunmode == 5 || fRunPara->frunmode == 6 || fRunPara->frunmode == 7 )
    {
        if( fRunPara->frunmode != 7 )
        {
            fRunPara->fcalibrationfile = "";
        }
        fRunPara->fDBTracking = false;
        //        fRunPara->fuseDB = false;
        fRunPara->fcalibrationrun = true;
    }
    if( fRunPara->frunmode != 1 && fRunPara->frunmode != 6 )
    {
        fRunPara->fPedestalsInTimeSlices = false;
    }
    // low gain pedestals: ignore sample length in cfg file
    if( fRunPara->frunmode == 6 )
    {
        fRunPara->fUseVBFSampleLength = true;
        fRunPara->fPedestalsInTimeSlices = false;
        fRunPara->fUsePedestalsInTimeSlices = false;
        fRunPara->fLowGainUsePedestalsInTimeSlices = false;
    }
    
    /////////////////////////////////////////////////////////////////
    // set vector sizes for calibration numbers
    for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
    {
        if( i >= fRunPara->fGainFileNumber.size() )
        {
            fRunPara->fGainFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fTOffFileNumber.size() )
        {
            fRunPara->fTOffFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fPedFileNumber.size() )
        {
            fRunPara->fPedFileNumber.push_back( fRunPara->frunnumber );   // default: take peds from same run
            fRunPara->fTZeroFileNumber.push_back( fRunPara->frunnumber );   // default: take tzeros from same run
        }
        else
        {
            fRunPara->fPedFileNumber[i] = fRunPara->frunnumber;
            fRunPara->fTZeroFileNumber[i] = fRunPara->frunnumber;
        }
        if( !fRunPara->fPlotRaw )
        {
            // writing low gain pedestals
            if( fRunPara->frunmode == 6 )
            {
                if( i >= fRunPara->fPedLowGainFileNumber.size() )
                {
                    fRunPara->fPedLowGainFileNumber.push_back( fRunPara->frunnumber );
                }
                else
                {
                    fRunPara->fPedLowGainFileNumber[i] = fRunPara->frunnumber;
                }
            }
            // reading low gain pedestals
            else
            {
                if( i >= fRunPara->fPedLowGainFileNumber.size() )
                {
                    fRunPara->fPedLowGainFileNumber.push_back( 0 );
                }
                else
                {
                    fRunPara->fPedLowGainFileNumber[i] = 0;
                }
            }
        }
        if( i >= fRunPara->fLowGainMultiplierFileNumber.size() )
        {
            fRunPara->fLowGainMultiplierFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fGainLowGainFileNumber.size() )
        {
            fRunPara->fGainLowGainFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fTOffLowGainFileNumber.size() )
        {
            fRunPara->fTOffLowGainFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fPixFileNumber.size() )
        {
            fRunPara->fPixFileNumber.push_back( 0 );
        }
        if( i >= fRunPara->fTZeroLowGainFileNumber.size() )
        {
            fRunPara->fTZeroLowGainFileNumber.push_back( 0 );
        }
    }
    // gain and toff calculation: set output run numbers
    if( fRunPara->frunmode == 2 || fRunPara->frunmode == 5 )
    {
        for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
        {
            fRunPara->fGainFileNumber[i] = fRunPara->frunnumber;
            fRunPara->fGainLowGainFileNumber[i] = fRunPara->frunnumber;
            fRunPara->fTOffFileNumber[i] = fRunPara->frunnumber;
            fRunPara->fTOffLowGainFileNumber[i] = fRunPara->frunnumber;
        }
    }
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    // read run info from database
    // (some values can be overwritten by the command line)
    if( fRunPara->fuseDB )
    {
        VDBRunInfo i_DBinfo( fRunPara->frunnumber, fRunPara->getDBServer(), fRunPara->fNTelescopes );
        if( i_DBinfo.isGood() )
        {
            i_DBinfo.print();
            
            if( fTargetName_overwriteDB.size() > 1 )
            {
                fRunPara->fTargetName = fTargetName_overwriteDB;
                cout << "Overwrite DB value for target name, set to " << fRunPara->fTargetName << endl;
            }
            else
            {
                fRunPara->fTargetName = i_DBinfo.getTargetName();
            }
            // DB coordinates are in J2000
            if( fDec_overwriteDB < -9998. )
            {
                fRunPara->fTargetDec = i_DBinfo.getTargetDec();
            }
            else
            {
                fRunPara->fTargetDec = fDec_overwriteDB;
                cout << "Overwrite DB value for target dec, set to " << fRunPara->fTargetDec << " deg" << endl;
            }
            if( fRA_overwriteDB < -9998. )
            {
                fRunPara->fTargetRA = i_DBinfo.getTargetRA();
            }
            else
            {
                fRunPara->fTargetRA = fRA_overwriteDB;
                cout << "Overwrite DB value for target ra, set to " << fRunPara->fTargetRA << " deg" << endl;
            }
            if( fWobbleNorth_overwriteDB < -9998. )
            {
                fRunPara->fWobbleNorth = i_DBinfo.getWobbleNorth();
            }
            else
            {
                fRunPara->fWobbleNorth = fWobbleNorth_overwriteDB;
                cout << "VReadRunParameter::test_and_adjustParams() info: overwriting DB wobble north (";
                cout << i_DBinfo.getWobbleNorth() << " deg)";
                cout << "with command line value: " << fWobbleNorth_overwriteDB << " deg" << endl;
            }
            if( fWobbleEast_overwriteDB < -9998. )
            {
                fRunPara->fWobbleEast = i_DBinfo.getWobbleEast();
            }
            else
            {
                fRunPara->fWobbleEast = fWobbleEast_overwriteDB;
                cout << "VReadRunParameter::test_and_adjustParams() info: overwriting DB wobble east (";
                cout << i_DBinfo.getWobbleEast() << " deg)";
                cout << "with command line value: " << fWobbleEast_overwriteDB << " deg" << endl;
            }
            fRunPara->fDBRunType = i_DBinfo.getRunType();
            fRunPara->fDBRunStartTimeSQL = i_DBinfo.getDataStartTimeSQL();
            fRunPara->fDBRunStoppTimeSQL = i_DBinfo.getDataStoppTimeSQL();
            fRunPara->fDBDataStartTimeMJD = i_DBinfo.getDataStartTimeMJD();
            fRunPara->fDBDataStoppTimeMJD = i_DBinfo.getDataStoppTimeMJD();
            fRunPara->fDBDataStartTimeSecOfDay = i_DBinfo.getDataStartTime();
            fRunPara->fDBDataStoppTimeSecOfDay = i_DBinfo.getDataStoppTime();
            fRunPara->fRunDuration = ( float )i_DBinfo.getDuration();
            if( fTelToAna == 0 )
            {
                fTelToAna = i_DBinfo.getTelToAna();
            }
            // get source file from run number (if run number is given and no sourcefile)
            if( fRunPara->fsourcefile.size() < 1 )
            {
                char iname[5000];
                // set raw data directories
                if( fRunPara->getDirectory_VBFRawData().size() > 0 )
                {
                    sprintf( iname, "%s/d%d/%d.cvbf", fRunPara->getDirectory_VBFRawData().c_str(), i_DBinfo.getRunDate(), fRunPara->frunnumber );
                    // (note: returns false if path name can be accessed)
                    if( gSystem->AccessPathName( iname ) )
                    {
                        sprintf( iname, "%s/data/d%d/%d.cvbf", fRunPara->getDirectory_VBFRawData().c_str(), i_DBinfo.getRunDate(), fRunPara->frunnumber );
                    }
                }
                else
                {
                    sprintf( iname, "data/d%d/%d.cvbf", i_DBinfo.getRunDate(), fRunPara->frunnumber );
                    cout << iname << endl;
                }
                fRunPara->fsourcefile = iname;
            }
            
            // get laser runs
            if( fRunPara->frunmode != 2 && fRunPara->frunmode != 5 )
            {
                vector< unsigned int > iL = i_DBinfo.getLaserRun();
                if( iL.size() != fRunPara->fNTelescopes )
                {
                    cout << "VReadRunParameter::test_and_adjustParams() error: list of laser file has wrong length ";
                    cout << iL.size() << "\t" << fRunPara->fNTelescopes << endl;
                    exit( EXIT_FAILURE );
                }
                else
                {
                    for( unsigned int i = 0; i < iL.size(); i++ )
                    {
                        fRunPara->fGainFileNumber[i] = ( int )iL[i];
                        fRunPara->fTOffFileNumber[i] = ( int )iL[i];
                    }
                }
            }
            
        }
        else
        {
            cout << endl;
            cout << "FATAL ERROR: cannot connect to VERITAS database" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    if( !readEpochsAndAtmospheres() )
    {
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    /////////////////////////////////////////////////////////////////
    // check if sourcefile is given
    if( fRunPara->fsourcefile.size() < 1 )
    {
        cout << "error: no sourcefile (test failed in VReadRunParameter)" << endl;
        cout << "possible reasons:" << endl;
        cout << "\t - source file does not exist" << endl;
        cout << "\t - eventdisplay was compiled without database support and couldn't find corresponding data directory" << endl;
        cout << "exit...." << endl;
        exit( EXIT_FAILURE );
    }
    
    if( fRunPara->fsourcefile.find( "vbf" ) < fRunPara->fsourcefile.size() && fRunPara->fsourcetype == 0 )
    {
        fRunPara->fsourcetype = 3;
    }
    
    // bz2 and gz reading does not work for vbf files; dont allow
    if( fRunPara->fsourcefile.find( "vbf" ) < fRunPara->fsourcefile.size() )
    {
        if( fRunPara->fsourcefile.find( ".gz" ) < fRunPara->fsourcefile.size() || fRunPara->fsourcefile.find( ".bz2" ) < fRunPara->fsourcefile.size() )
        {
            cout << "error: cannot read gzipped or bzipped files" << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // can't apply T-point corrections to pointing monitor data
    if( fRunPara->fDBTrackingCorrections.size() > 0 && fRunPara->fPMTextFileDirectory.size() > 0 )
    {
        cout << "error: can't apply T-point corrections to pointing monitor data" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    // switch of display in calibration and dst mode
    if( fRunPara->fdisplaymode == true && fRunPara->frunmode != 0 )
    {
        fRunPara->fdisplaymode = 0;
    }
    if( fRunPara->fdisplaymode )
    {
        fRunPara->fWriteTriggerOnly = false;
    }
    
    // fRunPara->fImageLL: values between 0-2
    if( fRunPara->fImageLL < 0 || fRunPara->fImageLL > 2 )
    {
        cout << "warning: logl parameter out of range, setting it to 1" << endl;
        fRunPara->fImageLL = 1;
    }
    // check fadc trace scaling
    if( fRunPara->fMCScale <= 0. )
    {
        cout << "errors: fadcscale <= 0" << endl;
        exit( EXIT_FAILURE );
    }
    
    
    // set target name for calibration runs
    if( fRunPara->frunmode == 1 && fRunPara->fTargetName == "NONAME" )
    {
        fRunPara->fTargetName = "wy";
    }
    if( fRunPara->frunmode == 6 && fRunPara->fTargetName == "NONAME" )
    {
        fRunPara->fTargetName = "wy";
    }
    if( fRunPara->frunmode == 2 && fRunPara->fTargetName == "NONAME" )
    {
        fRunPara->fTargetName = "laser";
    }
    if( fRunPara->frunmode == 5 && fRunPara->fTargetName == "NONAME" )
    {
        fRunPara->fTargetName = "laser";
    }
    
    // remove .cam if present
    for( unsigned int i = 0; i < fRunPara->fcamera.size(); i++ )
    {
        if( fRunPara->fcamera[i].find( ".cam" ) < fRunPara->fcamera[i].size() )
        {
            fRunPara->fcamera[i] = fRunPara->fcamera[i].substr( 0, fRunPara->fcamera[i].find( ".cam" ) - 1 );
            cout << "adjusted camera name: " << fRunPara->fcamera[i] << endl;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // set default camera configurations
    if( !fusercamera )
    {
        ////////////////////////////////////////////
        // configuration from Sept 2009 to July 2012: V5
        // (after the T1 move)
        if( fRunPara->frunnumber > 46642 && fRunPara->frunnumber < 63409 )
        {
            fRunPara->fcamera[0] = "EVN_V5_Oct2012_newArrayConfig_20121027_v420.txt";
        }
        ////////////////////////////////////////////
        // configuration until August 2009: V4
        // (before the T1 move)
        else if( fRunPara->frunnumber <= 46642 )
        {
            fRunPara->fcamera[0] = "EVN_V4_Oct2012_oldArrayConfig_20130428_v420.txt";
        }
        ////////////////////////////////////////////
        // configuration from Sep 2012: V6
        // (after the PMT and camera upgrade
        else
        {
            fRunPara->fcamera[0] = "EVN_V6_Upgrade_20121127_v420.txt";
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    
    // set camera file name to dstfile for dst reading
    if( fRunPara->fsourcetype == 4 || fRunPara->fsourcetype == 7 )
    {
        for( unsigned int t = 0; t < fRunPara->fcamera.size(); t++ )
        {
            fRunPara->fcamera[t] = "dstfile";
        }
    }
    
    setDirectories();
    
    // for command line case, all parameters are the same for all telescopes
    if( f_boolCommandline )
    {
        // start at 1, parameters for first telescope are filled by default
        for( unsigned int i = 1; i < fRunPara->fNTelescopes; i++ )
        {
            fRunPara->fcamera.push_back( fRunPara->fcamera[0] );
            fRunPara->fImageCleaningParameters.push_back( new VImageCleaningRunParameter( fRunPara->fImageCleaningParameters[0]->getName() ) );
            fRunPara->fImageCleaningParameters.back()->setTelID( i );
            // use fixed image/border thresholds for the case of no pedvars in the DSTs
            if( fRunPara->fCalibrationDataType == 0 )
            {
                fRunPara->fImageCleaningParameters.back()->fUseFixedThresholds = true;
            }
            fRunPara->fImageCleaningParameters_DB_Pass1.push_back( new VImageCleaningRunParameter( fRunPara->fImageCleaningParameters_DB_Pass1[0]->getName() ) );
            
            fRunPara->fDoublePass.push_back( fRunPara->fDoublePass[0] );
            fRunPara->fDoublePassErrorWeighting2005.push_back( fRunPara->fDoublePassErrorWeighting2005[0] );
            fRunPara->fGainCorrection.push_back( fRunPara->fGainCorrection[0] );
            fRunPara->fsumwindow_1.push_back( fRunPara->fsumwindow_1[0] );
            fRunPara->fsumwindow_2.push_back( fRunPara->fsumwindow_2[0] );
            fRunPara->fsumwindow_pass1.push_back( fRunPara->fsumwindow_pass1[0] );
            fRunPara->fsumfirst.push_back( fRunPara->fsumfirst[0] );
            fRunPara->fSearchWindowLast.push_back( fRunPara->fSearchWindowLast[0] );
            fRunPara->fTraceWindowShift.push_back( fRunPara->fTraceWindowShift[0] );
            fRunPara->fsumfirst_startingMethod.push_back( fRunPara->fsumfirst_startingMethod[0] );
            fRunPara->fsumfirst_maxT0startDiff.push_back( fRunPara->fsumfirst_maxT0startDiff[0] );
            fRunPara->fTraceIntegrationMethod.push_back( fRunPara->fTraceIntegrationMethod[0] );
            fRunPara->fTraceIntegrationMethod_pass1.push_back( fRunPara->fTraceIntegrationMethod_pass1[0] );
            fRunPara->fDF_DigitalFilter.push_back( fRunPara->fDF_DigitalFilter[0] );
            fRunPara->fDF_UpSample.push_back( fRunPara->fDF_UpSample[0] );
            fRunPara->fDF_PoleZero.push_back( fRunPara->fDF_PoleZero[0] );
            fRunPara->fLogLikelihoodLoss_min.push_back( fRunPara->fLogLikelihoodLoss_min[0] );
            fRunPara->fLogLikelihoodLoss_max.push_back( fRunPara->fLogLikelihoodLoss_max[0] );
            fRunPara->fLogLikelihood_Ntubes_min.push_back( fRunPara->fLogLikelihood_Ntubes_min[0] );
            fRunPara->fSumWindowMaxTimedifferenceToDoublePassPosition.push_back( fRunPara->fSumWindowMaxTimedifferenceToDoublePassPosition[0] );
            fRunPara->fSumWindowMaxTimeDifferenceLGtoHG.push_back( fRunPara->fSumWindowMaxTimeDifferenceLGtoHG[0] );
            fRunPara->fNeighbourDistanceFactor.push_back( fRunPara->fNeighbourDistanceFactor[0] );
            fRunPara->fSquarePixels.push_back( fRunPara->fSquarePixels[0] );
            fRunPara->fNeighbourMultiplicity.push_back( fRunPara->fNeighbourMultiplicity[0] );
        }
    }
    
    // set pointing errors
    fRunPara->fPointingErrorX.clear();
    fRunPara->fPointingErrorY.clear();
    for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
    {
        if( f_pointingErrorX.find( i ) != f_pointingErrorX.end() )
        {
            fRunPara->fPointingErrorX.push_back( f_pointingErrorX[i] );
        }
        else
        {
            fRunPara->fPointingErrorX.push_back( 0. );
        }
        if( f_pointingErrorY.find( i ) != f_pointingErrorY.end() )
        {
            fRunPara->fPointingErrorY.push_back( f_pointingErrorY[i] );
        }
        else
        {
            fRunPara->fPointingErrorY.push_back( 0. );
        }
    }
    
    // set telescope to analyse
    //   1,2,3 = Telescope 1,2,3
    //   12    = Telescope 1,2
    //   23    = Telescope 2,3
    //
    fRunPara->fTelToAnalyze.clear();
    if( fTelToAna != 0 && fRunPara->fNTelescopes < 10 && fTelToAna > 10 )
    {
        unsigned int iM = 10;
        for( int i = 0; i < 10; i++ )
        {
            unsigned int tt = ( fTelToAna % ( iM ) ) / ( iM / 10 );
            if( tt == 0 )
            {
                break;
            }
            fRunPara->fTelToAnalyze.push_back( tt - 1 );
            iM *= 10;
        }
        sort( fRunPara->fTelToAnalyze.begin(), fRunPara->fTelToAnalyze.end() );
    }
    else if( fTelToAna != 0 && fTelToAna < 10 )
    {
        fRunPara->fTelToAnalyze.push_back( fTelToAna - 1 );
    }
    else if( fTelToAna == 0 && fTelToAnaString.size() > 0 )
    {
        int pos_prev = -1;
        size_t pos_curr = 0;
        vector< int > comma_positions;
        comma_positions.push_back( -1 );
        while( ( pos_curr = fTelToAnaString.find( ",", pos_prev + 1 ) ) < fTelToAnaString.size() )
        {
            comma_positions.push_back( pos_curr );
            pos_prev = pos_curr;
        }
        for( size_t i = 0; i < comma_positions.size(); i++ )
        {
            string tmpstr = fTelToAnaString.substr( comma_positions[i] + 1, ( i == comma_positions.size() - 1 ) ? fTelToAnaString.size() : comma_positions[i + 1] - comma_positions[i] - 1 );
            fRunPara->fTelToAnalyze.push_back( atoi( tmpstr.c_str() ) - 1 );
        }
        sort( fRunPara->fTelToAnalyze.begin(), fRunPara->fTelToAnalyze.end() );
    }
    // nothing set: analyse all telescope
    else
    {
        for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
        {
            fRunPara->fTelToAnalyze.push_back( i );
        }
    }
    
    // check if it is possible to analyze the requested telescope
    if( fRunPara->fTelToAnalyze.size() > fRunPara->fNTelescopes )
    {
        cout << "error: can't analyze this telescope ";
        cout << fRunPara->fTelToAnalyze.size() << "\t";
        cout << "(total number of telescopes: " << fRunPara->fNTelescopes << ")" << endl;
    }
    
    for( unsigned int i = 0; i < fRunPara->fTelToAnalyze.size(); i++ )
    {
        if( fRunPara->fTelToAnalyze[i] > fRunPara->fNTelescopes )
        {
            cout << "telescope number to analyze out of range " << fRunPara->fTelToAnalyze[i] << "\t" << fRunPara->fNTelescopes << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // gain/toff only per telescope
    if( ( fRunPara->fTelToAnalyze.size() > 1 || fRunPara->fTelToAnalyze.size() == 0 ) && ( fRunPara->frunmode == 2 || fRunPara->frunmode == 5 ) )
    {
        cout << fRunPara->fTelToAnalyze.size() << "\t" << fRunPara->frunmode << endl;
        cout << "gain/toff calibration only possible for one telescope, use teltoana command line parameter" << endl;
        exit( EXIT_FAILURE );
    }
    
    // set trace window shift to zero for fixed window start
    for( unsigned int t = 0; t < fRunPara->fsumfirst_startingMethod.size(); t++ )
    {
        if( fRunPara->fsumfirst_startingMethod[t] == 0 )
        {
            fRunPara->fTraceWindowShift[t] = 0;
        }
    }
    
    // throws an error if its a simulation file and it was asked to calculate the PW parameters from the CFD hits, since CFDs don't exist in the simulation record
    if( fRunPara->fIsMC > 0 && ( fRunPara->fPWmethod == 0 || fRunPara->fPWmethod == 1 ) )
    {
        cout << " GrISU simulations don't have true CFD hits! PWmethod=" << fRunPara->fPWmethod << " can only be used with data with CFD hits." << endl;
        cout << " Switching to FADC signals cleaned with PWcleanThreshold!!  " << endl << endl;
    }
    
    if( fRunPara->fPWmethod > 3 )
    {
        cout << "Incorrect method for calculating binary-image trigger map: " << fRunPara->fPWmethod << endl;
        exit( EXIT_FAILURE );
    }
    if( fRunPara->fPWcleanNeighbors < 0 )
    {
        cout << "Incorrect number of nearest neighbors required for cleaning: " << fRunPara->fPWcleanNeighbors << endl;
        exit( EXIT_FAILURE );
    }
    if( fRunPara->fPWlimit < 0 )
    {
        cout << "Can't send less than zero pixels to moment-generating function: " << fRunPara->fPWlimit << endl;
        exit( EXIT_FAILURE );
    }
    
    // (TMPTMPTMP ) don't start frogs analysis unless it's a 4-Telescope run
    if( fRunPara->ffrogsmode && fRunPara->fTelToAnalyze.size() != 4 )
    {
        cout << "Error: You requested a frogs analysis for " << fRunPara->fTelToAnalyze.size() << " telescopes; will only do frogs analysis for exactly 4 telescopes!" << endl;
        exit( EXIT_FAILURE );
    }
    
    
    if( fPrintOutputFile )
    {
        cout << fRunPara->foutputfileName << endl;
        exit( EXIT_SUCCESS );
    }
    
    printStartMessage();
}

void VReadRunParameter::printShortHelp()
{
    cout << endl;
    cout << "Type 'evndisp help' for full list of command line parameters";
    cout << endl;
    cout << "VERITAS analysis ";
    cout << " (please check scripts in directory ./scripts/VTS/)" << endl;
    cout << endl;
    cout << "  How to display events: " << endl;
    cout << endl;
    cout << "  How to analyse a run" << endl;
    cout << endl;
    cout << endl;
    cout << "CTA analysis: please check scripts in directory ./scripts/CTA" << endl;
    cout << endl;
}


void VReadRunParameter::printHelp()
{
    if( fDebug )
    {
        cout << "VReadRunParameter::printHelp()" << endl;
    }
    
    printStartMessage();
    cout << endl;
    if( gSystem->Getenv( "EVNDISPSYS" ) )
    {
        int syst_ret = system( "cat $EVNDISPSYS/README/README.EVNDISP" );
        if( syst_ret == -1 )
        {
            cout << "VReadRunParameter::printHelp() error: could not find helper file in README directory" << endl;
        }
    }
    else
    {
        cout << "VReadRunParameter::printHelp(): no help available (environmental variable EVNDISPSYS not set)" << endl;
    }
}


void VReadRunParameter::setDirectories()
{
    char i_text[600];
    
    // outputfilename
    //suppress output file for peds			gain/toffset		low gain/toffset		lpeds			DSTs
    if( fRunPara->frunmode == 1 || fRunPara->frunmode == 2 || fRunPara->frunmode == 5 || fRunPara->frunmode == 6 || fRunPara->frunmode == 4 )
    {
        fRunPara->foutputfileName = "-1";
    }
    if( fRunPara->foutputfileName != "-1" && fRunPara->foutputfileName.size() == 0 )
    {
        sprintf( i_text, "%s/%d.root", fRunPara->getDirectory_EVNDISPOutput().c_str(), fRunPara->frunnumber );
        fRunPara->foutputfileName = i_text;
        // check if output directory exists, otherwise create it
        string i_worDir;
        i_worDir = gSystem->WorkingDirectory();
        if( !gSystem->cd( fRunPara->getDirectory_EVNDISPOutput().c_str() ) )
        {
            if( gSystem->mkdir( fRunPara->getDirectory_EVNDISPOutput().c_str() ) != 0 )
            {
                // check again if the directory exists - might have been created
                // by another job in the meanwhile
                if( !gSystem->cd( fRunPara->getDirectory_EVNDISPOutput().c_str() ) )
                {
                    cout << "VReadRunParameter::test_and_adjustParams() warning: unable to create output directory: " << endl;
                    cout <<  fRunPara->getDirectory_EVNDISPOutput() << endl;
                    exit( EXIT_FAILURE );
                }
            }
            cout << "VReadRunParameter::test_and_adjustParams() info: found / created output directory" << endl;
        }
        gSystem->cd( i_worDir.c_str() );
    }
    
    // check if calibration directories exist, otherwise crate
    if( fRunPara->frunmode == 1 || fRunPara->frunmode == 2 || fRunPara->frunmode == 5 || fRunPara->frunmode == 6 || fRunPara->frunmode == 7 )
    {
        for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
        {
            sprintf( i_text, "%s/Tel_%d/", fRunPara->getDirectory_EVNDISPCalibrationData_perRun().c_str(), i + 1 );
            if( gSystem->AccessPathName( i_text ) )
            {
                cout << "\t creating calibration directory for Telescope " << i + 1 << " : " << i_text << endl;
                if( gSystem->mkdir( i_text, kTRUE ) != 0 )
                {
                    cout << "VReadRunParameter::test_and_adjustParams() error: unable to create calibration directory for Telescope ";
                    cout << i + 1 << ": " << endl;
                    cout << i_text << endl;
                    exit( EXIT_FAILURE );
                }
            }
        }
    }
}


void VReadRunParameter::isCompiledWithDB()
{
#ifndef RUNWITHDB
    cout << endl;
    cout << "#########################################################################################" << endl;
    cout << "ERROR: eventdisplay is not compiled with mysql support" << endl;
    cout << "       unable to use database" << endl;
    cout << "       do not use command line parameters for database usage" << endl;
    cout << endl;
    cout << "       uncomment line 'DBFLAG=-DRUNWITHDB' in Makefile to compile eventdisplay with mysql support" << endl;
    cout << "#########################################################################################" << endl;
    exit( EXIT_FAILURE );
#endif
}


/*!

   open DST file and read run parameters

*/
bool VReadRunParameter::getRunParametersFromDST()
{
    // check if suffix is .root (required for dst files)
    if( fRunPara->fsourcefile.size() < 6 || fRunPara->fsourcefile.substr( fRunPara->fsourcefile.size() - 5, 5 ) != ".root" )
    {
        return false;
    }
    if( fDebug )
    {
        cout << "VReadRunParameter::getRunParametersFromDST() " << endl;
    }
    
    cout << "reading run parameters from dst file: " << fRunPara->fsourcefile << endl;
    // open dst file
    TFile iF( fRunPara->fsourcefile.c_str() );
    if( iF.IsZombie() )
    {
        cout << "VReadRunParameter::getRunParametersFromDST error opening DST file: " << fRunPara->fsourcefile << endl;
        exit( EXIT_FAILURE );
    }
    // check if this is a MC or data DST
    if( ( TTree* )iF.Get( "mc" ) )
    {
        fRunPara->fsourcetype = 7;
        fRunPara->fIsMC = 1;
    }
    else
    {
        fRunPara->fsourcetype = 4;
        fRunPara->fIsMC = 0;
    }
    fRunPara->fuseDB = false;
    fRunPara->fDBTracking = false;
    VEvndispRunParameter* iV = ( VEvndispRunParameter* )iF.Get( "runparameterDST" );
    if( !iV )
    {
        iV = ( VEvndispRunParameter* )iF.Get( "runparameterV2" );
    }
    if( iV )
    {
        fRunPara->fTargetName = iV->fTargetName;
        if( fRunPara->fTargetName == "NOSOURCE" )
        {
            fRunPara->fTargetName = "laser";
        }
        fRunPara->frunnumber  = iV->frunnumber;
        fRunPara->fNTelescopes = iV->fNTelescopes;
        fRunPara->fIsMC = iV->fIsMC;
        if( fRunPara->fIsMC )
        {
            fRunPara->fsourcetype = 7;
        }
        else
        {
            fRunPara->fsourcetype = 4;
        }
        fRunPara->getObservatory() = iV->getObservatory();
        fRunPara->fCalibrationDataType = iV->fCalibrationDataType;
        fRunPara->fFADCChargeUnit = iV->fFADCChargeUnit;
        fRunPara->fpulsetiminglevels = iV->fpulsetiminglevels;
        fRunPara->setPulseZeroIndex();
    }
    // no run parameters found
    else
    {
        // get dst tree and read run number and number of telescopes
        TTree* t = ( TTree* )iF.Get( "dst" );
        if( !t )
        {
            cout << "VReadRunParameter::getRunParametersFromDST error: DST tree not found" << endl;
            exit( EXIT_FAILURE );
        }
        unsigned int iR = 0;
        unsigned int iNTel = 0;
        t->SetBranchAddress( "runNumber", &iR );
        t->SetBranchAddress( "ntel", &iNTel );
        
        if( t->GetEntries() > 0 )
        {
            t->GetEntry( 0 );
            fRunPara->frunnumber = iR;
            fRunPara->fNTelescopes = iNTel;
        }
        else
        {
            return false;
        }
    }
    
    return true;
}

void VReadRunParameter::printStartMessage()
{
    cout << endl;
    cout << "\t---------------------------------------------------" << endl;
    cout << "\t|                                                 | " << endl;
    cout << "\t|        IACT event analysis and display          |" << endl;
    cout << "\t|                                                 | " << endl;
    cout << "\t|     \t\t Version " << fRunPara->getEVNDISP_VERSION() << "                   |" << endl;
    cout << "\t|                                                 | " << endl;
    cout << "\t---------------------------------------------------" << endl;
}

bool VReadRunParameter::checkSecondArgument( string iPara1, string iPara2, bool iExitIfFails )
{
    if( ( iPara2.size() > 0 && iPara2.substr( 0, 1 ) == "-" ) || iPara2.size() == 0 )
    {
        if( iExitIfFails )
        {
            cout << "command line parameter <" << iPara1 << "> requires second argument" << endl;
            exit( EXIT_FAILURE );
        }
        else
        {
            return false;
        }
    }
    
    return true;
}

/*
 *  read instrument epoch and atmospheric ID from a parameter file
 *
 *  (this might be extended in future to read the values from the db)
 *
 *  an epoch might be summer or winter
 *
 */
bool VReadRunParameter::readEpochsAndAtmospheres()
{
    if( !fRunPara )
    {
        return false;
    }
    
    if( fRunPara->fEpochGain.size() != fRunPara->fNTelescopes )
    {
        fRunPara->fEpochGain.resize( fRunPara->fNTelescopes, 5.5 );
    }
    
    // no epoch file given - return
    if( fRunPara->fEpochFile.size() == 0 )
    {
        return true;
    }
    
    ifstream is;
    is.open( fRunPara->fEpochFile.c_str(), ifstream::in );
    if( !is )
    {
        string iTemp = fRunPara->getDirectory_EVNDISPParameterFiles() + fRunPara->fEpochFile;
        is.open( iTemp.c_str(), ifstream::in );
        if( !is )
        {
            cout << "error opening epoch parameter file " << fRunPara->fEpochFile << endl;
            cout << iTemp << endl;
            exit( EXIT_FAILURE );
        }
    }
    string is_line;
    string temp;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            
            int run_min = 0;
            int run_max = 0;
            
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            if( is_stream.eof() )
            {
                continue;
            }
            is_stream >> temp;
            // EPOCH (e.g. V4, V5 and V6)
            // * EPOCH V6 <run min> <run max>
            if( temp == "EPOCH" && fRunPara->fInstrumentEpoch == "noepoch" )
            {
                string iTemp = "";
                if( !is_stream.eof() )
                {
                    is_stream >> iTemp;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> run_min;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> run_max;
                }
                
                if( fRunPara->frunnumber >= run_min && fRunPara->frunnumber <= run_max )
                {
                    fRunPara->fInstrumentEpoch = iTemp;
                }
            }
            // atmosphere (e.g. summer or winter)
            //    expect input string in sql format without hours: 2014-06-16)
            //  * ATMOSPHERE 21 2014-06-01 2014-06-16
            else if( temp == "ATMOSPHERE" && fRunPara->fAtmosphereID == 0 )
            {
                string iTemp = "";
                string date_min = "";
                string date_max = "";
                double imjd_min = 0.;
                double isec_min = 0.;
                double imjd_max = 0.;
                double isec_max = 0.;
                if( !is_stream.eof() )
                {
                    is_stream >> iTemp;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> date_min;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> date_max;
                }
                date_min += " 12:00:00";
                date_max += " 12:00:00";
                if( VSkyCoordinatesUtilities::getMJD_from_SQLstring( date_min, imjd_min, isec_min ) != 0 )
                {
                    cout << "VReadRunParameter::readEpochsAndAtmospheres() error: invalid date string: " << date_min << endl;
                    return false;
                }
                if( VSkyCoordinatesUtilities::getMJD_from_SQLstring( date_max, imjd_max, isec_max ) != 0 )
                {
                    cout << "VReadRunParameter::readEpochsAndAtmospheres() error: invalid date string: " << date_max << endl;
                    return false;
                }
                if( fRunPara->fDBDataStartTimeMJD >= imjd_min && fRunPara->fDBDataStartTimeMJD <= imjd_max )
                {
                    fRunPara->fAtmosphereID = atoi( iTemp.c_str() );
                }
            }
            // absolute gains (dc/pe) for each VERITAS epoch
            //  * GAIN V4 1 5.11
            else if( temp == "GAIN" )
            {
                string iEpoch = "";
                string iTel = "";
                string iGain = "";
                if( !is_stream.eof() )
                {
                    is_stream >> iEpoch;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> iTel;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> iGain;
                }
                unsigned int iTelNum = 0;
                
                if( iEpoch == fRunPara->fInstrumentEpoch )
                {
                    iTelNum = atoi( iTel.c_str() );
                    for( unsigned int i = 0; i < fRunPara->fNTelescopes; i++ )
                    {
                        if( iTelNum == i + 1 )
                        {
                            fRunPara->fEpochGain[i] = atof( iGain.c_str() );
                        }
                    }
                }
            }
        }
    }
    is.close();
    
    // print warnings / errors if no correct atmosphere is read from epoch file
    if( fRunPara->fAtmosphereID <= 0 )
    {
        cout << "Error file reading epoch for summer/winter atmosphere from " << fRunPara->fEpochFile << endl;
        cout << "\t no corresponding epoch fine. Exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    return true;
}
