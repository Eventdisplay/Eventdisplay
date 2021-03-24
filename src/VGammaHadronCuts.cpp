/*! \class VGammaHadronCuts
  \brief class containing parameter cut definitions


  List of cut selectors:

  gamma/hadron cut selector values consist of two digits: ID1+ID2*10

  ID2:

     0: apply gamma/hadron cuts on parameters in given data tree
     1: apply gamma/hadron cuts on probabilities given by a friend to the data tree (e.g. random forest analysis)
     2: same as 2
     3: apply cuts on probabilities given by a friend to the data tree already at the level of
        the event quality level
     4: TMVA gamma/hadron separation

  ID1:

     0: apply cuts on MSCW/MSCL (mean reduced scaled width/length)
     1: apply cuts on mean width/length (no lookup tables necessary)
     2: no cut applied (always passed)
     3: apply cuts on MWR/MLR (mean scaled width/length)

  Example:

    cut selector = 0 : apply MSCW/MSCL cuts (default)
    cut selector = 22 : apply event probability cuts
    cut selector = 10 : apply cuts from a tree AND apply MSCW/MSCL cuts

  direction cut selector value consist of one digit:

     0: energy independent direction cut (theta2 cut; set fArrayTheta2_max and fArrayTheta2_min)
        (in cut file: *theta2cut <theta2cut_min> <theta2cut_max>
     1: energy dependent TF1 function
        (in cut file: *theta2file ...)
     2: graph from instrument response function calculation
        (in cut file: *theta2file ... IRF)
     3: TMVA: use gamma/hadron part evaluation (no direction cuts applied in VGammaHadronCuts)
     4: TMVA: get direction cut from TMVA evaluator for each event
     5: TMVA: get direction cut from theta2 graph obtained during initialization of TMVA evaluator


*/

#include "VGammaHadronCuts.h"

VGammaHadronCuts::VGammaHadronCuts( string iCutID )
{
    fCutID = iCutID;
    setDebug( false );
    
    resetCutValues();
    
    fStats = 0;
    
    fReconstructionType = NOT_SET;
    fGammaHadronCutSelector = 0;
    fDirectionCutSelector = 0;
    bMCCuts = false;
    
    fData = 0;
    setNTel( 0 );
    fNLTrigs = 0;
    fDataDirectory = "";
    fInstrumentEpoch = "NOT_SET";
    
    // mean width/length/distance
    fMeanImageWidth = 0.;
    fMeanImageLength = 0.;
    fMeanImageDistance = 0.;
    
    // probabilities cuts
    fProbabilityCut_File = 0;
    fProbabilityCut_Tree = 0;
    fProbabilityCut_QualityFlag = 0;
    fProbabilityCut_NSelectors = VANACUTS_PROBSELECTIONCUTS_MAX;
    fProbabilityCut_ProbID = 0;
    for( unsigned int i = 0; i < fProbabilityCut_NSelectors; i++ )
    {
        fProbabilityCut_SelectionCut[i] = -1.;
    }
    
    // TMVA evaluator
    fTMVAEvaluator = 0;
    fTMVA_MVAMethod = "";
    fTMVA_MVAMethodCounter = 0;
    fTMVAWeightFileIndex_Emin = 0;
    fTMVAWeightFileIndex_Emax = 0;
    fTMVAWeightFileIndex_Zmin = 0;
    fTMVAWeightFileIndex_Zmax = 0;
    fTMVAEnergyStepSize = 0.2;
    fTMVAWeightFile = "";
    fTMVASignalEfficiency.clear();
    fTMVA_MVACut.clear();
    fTMVA_MVACutGraphFileName = "";
    fTMVA_MVACutGraphSmoothing = 0.2;
    fTMVA_MVACutGraphSmoothingMax = 1.e5;
    fTMVA_MVACutGraphConstantCutEnergy_TeV = 2.e10;
    // Note: for TMVA is this not the probability threshold but the MVA cut value
    fTMVAProbabilityThreshold = -99.;
    fTMVAOptimizeSignalEfficiencyParticleNumberFile = "";
    fTMVAParticleNumberFile_Conversion_Rate_to_seconds = 60.;
    fTMVAOptimizeSignalEfficiencySignificance_Min = 5.;
    fTMVAOptimizeSignalEfficiencySignalEvents_Min = 10.;
    fTMVAOptimizeSignalEfficiencyObservationTime_h = 50.;
    fTMVAFixedSignalEfficiencyMax = 1.;
    fTMVAFixedThetaCutMin = 0.;
    fTMVA_EvaluationResult = -99.;
    fTMVAEvaluatorResults = 0;
    fTMVAMinSourceStrength = 0.001;
    
    // energy dependent theta2 cut
    fFileNameAngRes = "";
    fFileAngRes = 0;
    fF1AngResName = "";
    fF1AngRes = 0;
    fAngRes_ScalingFactor = 1.;
    fAngRes_AbsoluteMinimum = 0.;
    fAngRes_AbsoluteMaximum = 1.e10;
    fAngRes_FixedAboveEnergy_TeV = 1.e30;
    fAngResContainmentProbability = 0;

    fCutCharacteristicsMCAZ = -999.;
    fCutCharacteristicsMCAZ_tolerance = 60.;
}

void VGammaHadronCuts::initialize()
{
    // statistics
    fStats = new VGammaHadronCutsStatistics();
    fStats->initialize( fCutID );
}

VGammaHadronCuts::~VGammaHadronCuts()
{
    if( fStats )
    {
        delete fStats;
    }
    for( unsigned int i = 0; i < fNTelTypeCut.size(); i++ )
    {
        if( fNTelTypeCut[i] )
        {
            delete fNTelTypeCut[i];
        }
    }
    if( fTMVAEvaluator )
    {
        delete fTMVAEvaluator;
    }
    if( fProbabilityCut_File && !fProbabilityCut_File->IsZombie() )
    {
        fProbabilityCut_File->Close();
    }
    map< string, TGraph* >::iterator it;
    for( it = fEnergyDependentCut.begin(); it != fEnergyDependentCut.end(); ++it )
    {
        if( it->second )
        {
            delete it->second;
        }
    }
}

void VGammaHadronCuts::resetCutStatistics()
{
    if( fStats )
    {
        fStats->reset();
    }
}

void VGammaHadronCuts::resetCutValues()
{
    // array cuts
    fCut_Theta2_min = -100.;
    fCut_Theta2_max = 100.;
    fCut_Chi2_min = -1.e6;
    fCut_Chi2_max = 1.e10;
    fCut_Dir_max = 1.e10;
    fCut_MeanImageDistance_min = -100;
    fCut_MeanImageDistance_max = 100;
    fCut_MeanImageLength_min = -100;
    fCut_MeanImageLength_max = 100;
    fCut_MeanImageWidth_min = -100;
    fCut_MeanImageWidth_max = 100;
    fCut_Size_min = -1000;
    fCut_Size_max = 1e10;
    fCut_MSCW_min   = -1.;
    fCut_MSCW_max   = 1.;
    fCut_MSCL_min   = -1.;
    fCut_MSCL_max   = 1.;
    fCut_MSW_min    = -1.;
    fCut_MSW_max    = 10.;
    fCut_MSL_min    = -1.;
    fCut_MSL_max    = 10.;
    fCut_CameraFiducialSize_min = -1.;
    fCut_CameraFiducialSize_max = 100.;
    fCut_CameraFiducialSize_MC_min = -100.;
    fCut_CameraFiducialSize_MC_max = 100.;
    fCut_dE_min = -99.;
    fCut_dE_max = 1.e12;
    fCut_EmultiplicityRatio_min = 0.;
    fCut_EChi2_min = 0.;
    fCut_EChi2_max = 99999.;
    fCut_Erec_min = 0.;
    fCut_Erec_max = 1.e10;
    fCut_Emmission_min = 0.;
    fCut_Emmission_max = 1.e12;
    fCut_NImages_min = 0;
    fCut_NImages_max = 100000;
    fCut_DispNImages_min = 0;
    fCut_DispNImages_max = 100000;
    
    fCut_AverageCoreDistanceToTelescopes_nimages = -99;
    fCut_AverageCoreDistanceToTelescopes_min = -99.;
    fCut_AverageCoreDistanceToTelescopes_max = 99999.;
    fCut_MinimumCoreDistanceToTelescopes_max = 1.e10;
    
    fCut_SizeSecondMax_min = -1000;
    fCut_SizeSecondMax_max = 1.e10;
    
    fProbabilityCut = 0.5;
}


/*

    read cuts from a text file

    iPrint == 0 : surpress printing of any statements to screen

*/
bool VGammaHadronCuts::readCuts( string i_cutfilename, int iPrint )
{
    // reset trigger vector
    fNLTrigs = 0;
    fCut_ImgSelect.clear();
    
    // open text file
    ifstream is;
    i_cutfilename = VUtilities::testFileLocation( i_cutfilename, "GammaHadronCutFiles", true );
    if( iPrint == 1 )
    {
        cout << "\t reading analysis cuts from " << i_cutfilename << endl;
    }
    else if( iPrint == 2 )
    {
        cout << "reading analysis cuts from " << i_cutfilename << endl;
    }
    
    TString iPath( i_cutfilename.c_str() );
    gSystem->ExpandPathName( iPath );
    is.open( iPath.Data(), ifstream::in );
    if( !is )
    {
        cout << "VGammaHadronCuts::readCuts: cut input file not found, " << i_cutfilename << endl;
        return false;
    }
    string is_line;
    string temp;
    string iCutVariable;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            is_stream >> iCutVariable;
            //////////////////////////////////////
            // choose gamma/hadron cut selectors and direction cut selector
            if( iCutVariable == "cutselection" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fGammaHadronCutSelector;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fDirectionCutSelector;
                }
            }
            /////////////////////////////////////
            // stereo ('array') telescope cuts
            else if( iCutVariable == "arraywidth" )
            {
                is_stream >> temp;
                fCut_MeanImageWidth_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MeanImageWidth_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "arraylength" )
            {
                is_stream >> temp;
                fCut_MeanImageLength_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MeanImageLength_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "arrayntel" || iCutVariable == "nnimages" )
            {
                is_stream >> temp;
                fCut_NImages_min = ( atoi( temp.c_str() ) );
                is_stream >> temp;
                fCut_NImages_max = ( atoi( temp.c_str() ) );
            }
            else if( iCutVariable == "dispnimages" )
            {
                is_stream >> temp;
                fCut_DispNImages_min = ( atoi( temp.c_str() ) );
                is_stream >> temp;
                fCut_DispNImages_max = ( atoi( temp.c_str() ) );
            }
            else if( iCutVariable == "arraysize" )
            {
                is_stream >> temp;
                fCut_Size_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_Size_max = ( atof( temp.c_str() ) );
                // write warning out about size cuts and number of telescopes
                if( fNTel > 2 && fCut_Size_min > 0. )
                {
                    cout << "--------- VGammaHadronCuts warning: ignoring size cut for data with more than 2 telescopes ------" << endl;
                }
            }
            else if( iCutVariable == "arraychi2" )
            {
                is_stream >> temp;
                fCut_Chi2_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_Chi2_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "dirdiff" )
            {
                is_stream >> temp;
                fCut_Dir_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "arraydist" )
            {
                is_stream >> temp;
                fCut_MeanImageDistance_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MeanImageDistance_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "mscw" || iCutVariable == "arraymscw" )
            {
                is_stream >> temp;
                fCut_MSCW_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MSCW_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "mscl" || iCutVariable == "arraymscl" )
            {
                is_stream >> temp;
                fCut_MSCL_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MSCL_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "msw" || iCutVariable == "arraymsw" )
            {
                is_stream >> temp;
                fCut_MSW_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MSW_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "msl" || iCutVariable == "arraymsl" )
            {
                is_stream >> temp;
                fCut_MSL_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_MSL_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "mc_xy_off" || iCutVariable == "arrayxyoff_mc" || iCutVariable == "mc_xyoff" || iCutVariable == "cameraedge_mc" )
            {
                is_stream >> temp;
                fCut_CameraFiducialSize_MC_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_CameraFiducialSize_MC_max = ( atof( temp.c_str() ) );
                bMCCuts = true;
            }
            
            else if( iCutVariable == "xy_off" || iCutVariable == "arrayxyoff" || iCutVariable == "xyoff" || iCutVariable == "cameraedge" )
            {
                is_stream >> temp;
                fCut_CameraFiducialSize_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_CameraFiducialSize_max = ( atof( temp.c_str() ) );
            }
            else if( iCutVariable == "telcoredistance" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCut_MinimumCoreDistanceToTelescopes_max;
                }
            }
            else if( iCutVariable == "arraycore" )
            {
                is_stream >> temp;
                fCut_AverageCoreDistanceToTelescopes_min = atof( temp.c_str() );
                is_stream >> temp;
                fCut_AverageCoreDistanceToTelescopes_max = atof( temp.c_str() );
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fCut_AverageCoreDistanceToTelescopes_nimages = atoi( temp.c_str() );
                }
            }
            // these allow to deselect certain telescope combinations (changed from ltrig to ImgSel)
            else if( iCutVariable == "arrayltrig" || iCutVariable == "allowedImageCombinations" )
            {
                if( fNTel == 0 || fNTel > 10 )
                {
                    cout << "VGammaHadronCuts::readCuts warning: cut identifier " << temp << " ignored for ";
                    cout << "current telescope configuration ( " << fNTel << "telescopes)" << endl;
                    continue;
                }
                // calculate how many possible telescope combinations exist (16 for 4 telescopes)
                if( fNLTrigs == 0 )
                {
                    int num_ltrigs = 0;
                    for( unsigned int i = 0; i < fNTel; i++ )
                    {
                        num_ltrigs += int( pow( 2., double( i ) ) );
                    }
                    num_ltrigs += 1;
                    fNLTrigs = num_ltrigs;
                    for( unsigned int i = 0; i < fNLTrigs; i++ )
                    {
                        fCut_ImgSelect.push_back( 1 );
                    }
                }
                // telescope combination
                is_stream >> temp;
                int index = ( atoi( temp.c_str() ) );
                // on = 0 or off = 1
                is_stream >> temp;
                int i_select = atoi( temp.c_str() );
                // check epoch
                bool i_useTheseCuts = true;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    if( temp != fInstrumentEpoch )
                    {
                        i_useTheseCuts = false;
                    }
                }
                // check telescope combinations
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    if( fTelToAnalyze.size() > 0 && temp != getTelToAnalyzeString() )
                    {
                        i_useTheseCuts = false;
                    }
                }
                // set telescope combinations on or off
                if( i_useTheseCuts )
                {
                    if( index < 0 )
                    {
                        for( unsigned int i = 0; i < fCut_ImgSelect.size(); i++ )
                        {
                            fCut_ImgSelect[i] = i_select;
                        }
                    }
                    else if( index < ( int )fCut_ImgSelect.size() )
                    {
                        fCut_ImgSelect[index] = i_select;
                    }
                }
            }
            // probability cut variables (e.g. random forest)
            else if( iCutVariable == "RFthresh" || iCutVariable == "Probthresh" )
            {
                is_stream >> temp;
                fProbabilityCut = ( atof( temp.c_str() ) );
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fProbabilityCut_ProbID = atoi( temp.c_str() );
                }
                else
                {
                    fProbabilityCut_ProbID = 0;
                }
            }
            
            // to define the lower bounds in probablity cut ranges  (e.g. random forest)
            else if( iCutVariable == "RFCutLowerVals" )
            {
                while( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fProbabilityCutRangeLower.push_back( atof( temp.c_str() ) );
                }
            }
            
            // to define the upper bounds in probablity cut ranges  (e.g. random forest)
            else if( iCutVariable == "RFCutUpperVals" )
            {
                while( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fProbabilityCutRangeUpper.push_back( atof( temp.c_str() ) );
                }
            }
            
            // to define the upper bounds in probablity cut ranges  (e.g. random forest)
            else if( iCutVariable == "RFProbID" )
            {
                is_stream >> temp;
                fProbabilityCut_ProbID = atoi( temp.c_str() );
            }
            
            
            // energy reconstruction cuts
            else if( iCutVariable == "arrayechi2" )
            {
                is_stream >> temp;
                fCut_EChi2_min = atof( temp.c_str() );
                is_stream >> temp;
                fCut_EChi2_max = atof( temp.c_str() );
            }
            else if( iCutVariable == "arraydE" )
            {
                is_stream >> temp;
                fCut_dE_min = atof( temp.c_str() );
                is_stream >> temp;
                fCut_dE_max = atof( temp.c_str() );
            }
            else if( iCutVariable == "arrayemultiplicityratio" )
            {
                is_stream >> temp;
                fCut_EmultiplicityRatio_min = atof( temp.c_str() );
                // no max value (as it cannot be larger than 1)
            }
            else if( iCutVariable == "arrayerec" )
            {
                is_stream >> temp;
                fCut_Erec_min = atof( temp.c_str() );
                is_stream >> temp;
                fCut_Erec_max = atof( temp.c_str() );
            }
            else if( iCutVariable == "arrayemission" || iCutVariable == "emissionheight" )
            {
                is_stream >> temp;
                fCut_Emmission_min = atof( temp.c_str() );
                is_stream >> temp;
                fCut_Emmission_max = atof( temp.c_str() );
            }
            else if( iCutVariable == "sizesecondmax" )
            {
                float isize_min = -1000.;
                float isize_max = 1.e10;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    isize_min = atof( temp.c_str() );
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    isize_max = atof( temp.c_str() );
                }
                // check instrument epoch
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    if( temp == fInstrumentEpoch )
                    {
                        fCut_SizeSecondMax_min = isize_min;
                        fCut_SizeSecondMax_max = isize_max;
                    }
                }
                else
                {
                    fCut_SizeSecondMax_min = isize_min;
                    fCut_SizeSecondMax_max = isize_max;
                }
            }
            // telescope type dependent cut on number of images (implementation 1)
            // syntax:  teltype_nnimages <min images> <tel type counter1> ..<tel type counter 2>
            else if( iCutVariable == "teltypegroup" )
            {
                fNTelTypeCut.push_back( new VNTelTypeCut() );
                is_stream >> temp;
                fNTelTypeCut.back()->fNTelType_min = atoi( temp.c_str() );
                while( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fNTelTypeCut.back()->fTelType_counter.push_back( atoi( temp.c_str() ) );
                }
            }
            // telescope type dependent cut on number of images (implementation 2)
            // syntax:  teltype_nnimages <min images > <tel type1 > <tel type2> ...
            else if( iCutVariable == "nimages_teltype" )
            {
                fNTelTypeCut.push_back( new VNTelTypeCut() );
                is_stream >> temp;
                fNTelTypeCut.back()->fNTelType_min = atoi( temp.c_str() );
                while( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fNTelTypeCut.back()->fTelTypeID.push_back( atol( temp.c_str() ) );
                }
                fNTelTypeCut.back()->purgeTelTypeIDs();
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // TMVA values
            else if( iCutVariable == "TMVAPARAMETER" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    if( temp == fInstrumentEpoch )
                    {
                        while( !(is_stream>>std::ws).eof() )
                        {
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVA_MVAMethod;
                            }
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVA_MVAMethodCounter;
                            }
                            // files should have endings _fTMVAWeightFileIndex_min to _fTMVAWeightFileIndex_max
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVAWeightFileIndex_Emin;
                            }
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVAWeightFileIndex_Emax;
                            }
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVAWeightFileIndex_Zmin;
                            }
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> fTMVAWeightFileIndex_Zmax;
                            }
                            if( !(is_stream>>std::ws).eof() )
                            {
                                if( !( is_stream >> fTMVAEnergyStepSize ) )
                                {
                                    cout << "VGammaHadronCuts::readCuts: missing TMVAPARAMETER energy step size  " << endl;
                                    break;
                                }
                            }
                            string iWeightFileDirectory;
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> iWeightFileDirectory;
                            }
                            string iWeightFileName;
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> iWeightFileName;
                            }
                            fTMVAWeightFile = gSystem->ExpandPathName( iWeightFileDirectory.c_str() );
                            // check if path name is complete
                            // note: this method returns FALSE if one **can** access the file
                            if( gSystem->AccessPathName( fTMVAWeightFile.c_str() ) )
                            {
                                fTMVAWeightFile = VGlobalRunParameter::getDirectory_EVNDISPAnaData() + "/" + fTMVAWeightFile;
                                if( gSystem->AccessPathName( fTMVAWeightFile.c_str() ) )
                                {
                                    cout << "VGammaHadronCuts::readCuts error,";
                                    cout << " weight file directory not found: ";
                                    cout << fTMVAWeightFile << endl;
                                    cout << "exiting..." << endl;
                                    exit( EXIT_FAILURE );
                                }
                            }
                            fTMVAWeightFile += "/" + iWeightFileName;
                            break;
                        }
                    }
                    else if( iPrint != 0 )
                    {
                        cout << "VGammaHadronCuts::readCuts: ignoring ";
                        cout << temp << " TMVAPARAMETER because instrument epoch is: " << fInstrumentEpoch << endl;
                        //cout << "VGammaHadronCuts::readCuts: skipping TMVAPARAMETER due to epoch mismatch:";
                        //cout << " required: " << fInstrumentEpoch << ", is: " << temp << endl;
                    }
                }
            }
            else if( iCutVariable == "TMVACUTS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAOptimizeSignalEfficiencySignificance_Min;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAOptimizeSignalEfficiencySignalEvents_Min;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    // observing time is given as "50h", or "5m", "5s"
                    if( temp.find( "h" ) != temp.npos )
                    {
                        fTMVAOptimizeSignalEfficiencyObservationTime_h = atof( temp.substr( 0, temp.find( "h" ) ).c_str() );
                    }
                    else if( temp.find( "m" ) != temp.npos )
                    {
                        fTMVAOptimizeSignalEfficiencyObservationTime_h = atof( temp.substr( 0, temp.find( "m" ) ).c_str() );
                        fTMVAOptimizeSignalEfficiencyObservationTime_h /= 60.;
                    }
                    else if( temp.find( "s" ) != temp.npos )
                    {
                        fTMVAOptimizeSignalEfficiencyObservationTime_h = atof( temp.substr( 0, temp.find( "s" ) ).c_str() );
                        fTMVAOptimizeSignalEfficiencyObservationTime_h /= 3600.;
                    }
                    else
                    {
                        fTMVAOptimizeSignalEfficiencyObservationTime_h = atof( temp.c_str() );
                    }
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAOptimizeSignalEfficiencyParticleNumberFile;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAFixedSignalEfficiencyMax;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAMinSourceStrength;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAFixedThetaCutMin;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVAParticleNumberFile_Conversion_Rate_to_seconds;
                }
            }
            else if( iCutVariable == "TMVASignalEfficiency" )
            {
                while( !(is_stream>>std::ws).eof() )
                {
                    unsigned int iKey = 0;
                    double iS = 0.;
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iKey;
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iS;
                    }
                    fTMVASignalEfficiency[iKey] = iS;
                }
            }
            else if( iCutVariable == "TMVA_MVACut" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    if( temp == fInstrumentEpoch )
                    {
                        while( !(is_stream>>std::ws).eof() )
                        {
                            unsigned int iKey = 0;
                            double iS = 0.;
                            if( is_stream >> iKey && is_stream >> iS )
                            {
                                fTMVA_MVACut[iKey] = iS;
                            }
                        }
                    }
                }
            }
            else if( iCutVariable == "TMVA_MVACutGraph" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVA_MVACutGraphFileName;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVA_MVACutGraphSmoothing;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVA_MVACutGraphSmoothingMax;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTMVA_MVACutGraphConstantCutEnergy_TeV;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // direction cut values
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            // fixed theta2 cut
            else if( iCutVariable == "arraytheta2" || iCutVariable == "theta2cut" )
            {
                is_stream >> temp;
                fCut_Theta2_min = ( atof( temp.c_str() ) );
                is_stream >> temp;
                fCut_Theta2_max = ( atof( temp.c_str() ) );
            }
            // read in values for energy dependent theta2 cut
            // * theta2file <root file> <function name>
            // (note that fF1AngResName == "IRF" means that the graph from the IRF file is extrapolated)
            else if( iCutVariable == "theta2file" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    string iFileNameAngRes;
                    is_stream >> iFileNameAngRes;
                    if( fFileNameAngRes.find( "$" ) != string::npos )
                    {
                        fFileNameAngRes = gSystem->ExpandPathName( iFileNameAngRes.c_str() );
                    }
                    else
                    {
                        fFileNameAngRes = iFileNameAngRes;
                    }
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fF1AngResName;
                }
                else
                {
                    fF1AngResName = "fitAngRes";
                }
            }
            
            // use angular resolution calculated for example from same data with makeEffectiveArea
            else if( iCutVariable == "angres" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fAngResContainmentProbability;    // should be an integer; probability x 100
                }
            }
            // theta2 scaling
            // * theta2scaling <scale factor> <minimum theta> <maximum theta>
            else if( iCutVariable == "theta2scaling" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fAngRes_ScalingFactor = atof( temp.c_str() );
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fAngRes_AbsoluteMinimum = atof( temp.c_str() );
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fAngRes_AbsoluteMaximum = atof( temp.c_str() );
                    fCut_Theta2_max = fAngRes_AbsoluteMaximum * fAngRes_AbsoluteMaximum;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fAngRes_FixedAboveEnergy_TeV = atof( temp.c_str() );
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }
    // check cut selection
    if( fReconstructionType == NOT_SET )
    {
        if( fGammaHadronCutSelector / 10 == 7 )
        {
            fReconstructionType = DEEPLEARNER;
        }
        else
        {
            fReconstructionType = GEO;
        }
    }
    if( fData )
    {
        fData->setReconstructionType( fReconstructionType );
    }
    
    // initialize theta2 cut (energy dependent)
    if( fFileNameAngRes.size() > 0 )
    {
        if( !initAngularResolutionFile() )
        {
            cout << "VGammaHadronCuts::readCuts error: error initializing angular resolution file" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    return true;
}

void VGammaHadronCuts::printDirectionCuts()
{
    cout << "Direction cuts (cut selector " << fDirectionCutSelector << ")" << endl;
    
    // theta2 cut (no energy dependence)
    if( fDirectionCutSelector == 0 )
    {
        cout << "Theta2 cut: ";
        if( fCut_Theta2_min > 0. )
        {
            cout << fCut_Theta2_min << " < Theta^2 [deg^2]";
            if( fCut_Theta2_max > 0. )
            {
                cout << "< " << fCut_Theta2_max;
            }
        }
        else if( fCut_Theta2_max > 0. )
        {
            cout << "Theta^2 [deg^2] < " << fCut_Theta2_max;
        }
        else
        {
            cout << "no Theta^2 cut set";
        }
        cout << endl;
    }
    // theta cut from TF1 function (energy dependent)
    else if( fFileAngRes && fF1AngRes )
    {
        cout << "Direction cut from angular resolution function in file: ";
        cout << fFileAngRes->GetName();
        cout << "(" << fF1AngRes->GetName() << ")";
    }
    // theta cut from IRF graph (energy depedendent)
    else if( getTheta2Cut_IRF_Max() || fDirectionCutSelector == 2 )
    {
        cout << "Direction cut from IRF graph " << endl;
        if( getTheta2Cut_IRF_Max() )
        {
            cout << "Number of points: " << getTheta2Cut_IRF_Max()->GetN() << endl;
            getTheta2Cut_IRF_Max()->Print();
        }
        else if( getAngularResolutionContainmentRadius() > 0. )
        {
            cout << " (calculated from same file, containment probability " << ( ( double )getAngularResolutionContainmentRadius() ) / 100. << ")" << endl;
        }
        else
        {
            cout << endl;
            cout << "VGammaHadronCuts::printDirectionCuts WARNING: no function found" << endl;
        }
    }
    // theta cut using TMVA
    else if( fDirectionCutSelector == 4 )
    {
        cout << "Direction cut from TMVA angular containment optimizer" << endl;
    }
    else if( fDirectionCutSelector == 3 || fDirectionCutSelector == 5 )
    {
        cout << "Direction cut from TMVA " << endl;
        if( fDirectionCutSelector == 3 )
        {
            cout << "using TMVA evaluator (part of gamma/hadron cuts)" << endl;
        }
        if( fDirectionCutSelector == 5 )
        {
            if( fDirectionCutSelector == 5 && getTheta2Cut_TMVA_max() )
            {
                getTheta2Cut_TMVA_max()->Print();
            }
        }
    }
    cout << "Direction cut scale factor " << fAngRes_ScalingFactor;
    cout << ", minimum : " << fAngRes_AbsoluteMinimum << " [deg] ";
    cout << ", maximum : " << fAngRes_AbsoluteMaximum << " [deg]";
    cout << " (constant above " << fAngRes_FixedAboveEnergy_TeV << " TeV)";
    cout << endl;
    
}

void VGammaHadronCuts::printCutSummary()
{
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "VGammaHadronCuts::printCutSummary()";
    cout << " (ntel=" << fTelToAnalyze.size() << ": T" << getTelToAnalyzeString() << ") ";
    cout << endl;
    cout << "Gamma/hadron cut selector: " << fGammaHadronCutSelector;
    cout << " (" << fReconstructionType << ")" << endl;
    if( fInstrumentEpoch != "NOT_SET" )
    {
        cout << "Instrument epoch selected: " << fInstrumentEpoch << endl;
    }
    else
    {
        cout << "Instrument epoch not set" << endl;
    }
    
    // direction cuts
    printDirectionCuts();
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // gamma/hadron cuts
    
    // mean reduced scaled cuts
    if( fGammaHadronCutSelector % 10 < 1 )
    {
        cout << "Shape cuts: ";
        cout << fCut_MSCW_min << " < MSCW < " << fCut_MSCW_max;
        cout << ", " << fCut_MSCL_min << " < MSCL < " << fCut_MSCL_max << ", ";
    }
    // mean cuts
    else if( fGammaHadronCutSelector % 10 == 1 )
    {
        cout << "Shape cuts: ";
        cout << fCut_MeanImageWidth_min  << " < mean width < " << fCut_MeanImageWidth_max;
        cout << ", " << fCut_MeanImageLength_min << " < mean length < " << fCut_MeanImageLength_max << ", ";
    }
    // mean scaled cuts
    else if( fGammaHadronCutSelector % 10 == 3 )
    {
        cout << "Shape cuts: ";
        cout << fCut_MSW_min << " < MWR < " << fCut_MSW_max;
        cout << ", " << fCut_MSL_min << " < MLR < " << fCut_MSL_max << ", ";
    }
    cout << endl;
    cout << "Average core distance < " << fCut_AverageCoreDistanceToTelescopes_max << " m";
    cout << " (max distance to telescopes (mintel) " << fCut_MinimumCoreDistanceToTelescopes_max << " m)";
    if( fCut_AverageCoreDistanceToTelescopes_nimages > 0 )
    {
        cout << " (NImages = " << fCut_AverageCoreDistanceToTelescopes_nimages << ")";
    }
    cout << endl;
    // probability cuts
    if( fGammaHadronCutSelector / 10 >= 1 && fGammaHadronCutSelector / 10 <= 3 )
    {
        cout << endl;
        cout << "event probability threshold: " << fProbabilityCut << " (element " << fProbabilityCut_ProbID << ")";
        if( fGammaHadronCutSelector / 10 == 3 )
        {
            cout << " (applied as quality cut)";
        }
        if( fProbabilityCut_File )
        {
            cout << ", read from " << fProbabilityCut_File->GetName();
        }
    }
    // TMVA cuts
    if( useTMVACuts() )
    {
        cout << endl;
        cout << "TMVA gamma/hadron separation with MVA method " << fTMVA_MVAMethod;
        cout << " (counter " << fTMVA_MVAMethodCounter << ")";
        cout << endl;
        cout << "weight files: " << fTMVAWeightFile;
        cout << " (" << fTMVAWeightFileIndex_Emin << "," << fTMVAWeightFileIndex_Emax << ")";
        cout << " (" << fTMVAWeightFileIndex_Zmin << "," << fTMVAWeightFileIndex_Zmax << ")" << endl;
        if( fTMVAOptimizeSignalEfficiencyParticleNumberFile.size() > 0. )
        {
            cout << "using optimal signal efficiency with a requirement of ";
            cout << fTMVAOptimizeSignalEfficiencySignificance_Min << " sigma and at least ";
            cout << fTMVAOptimizeSignalEfficiencySignalEvents_Min << " signal events in ";
            cout << fTMVAOptimizeSignalEfficiencyObservationTime_h << " h observing time" << endl;
            cout << "reading particle counts from " << fTMVAOptimizeSignalEfficiencyParticleNumberFile;
            cout << " (unit of rate graphs: " << fTMVAParticleNumberFile_Conversion_Rate_to_seconds << ")" << endl;
            cout << "   (max signal efficiency: " << fTMVAFixedSignalEfficiencyMax << ")" << endl;
            cout << "   (min source strength: " << fTMVAMinSourceStrength << ")" << endl;
        }
        else
        {
            if( fDebug )
            {
                printSignalEfficiency();
                printTMVA_MVACut();
            }
        }
    }
    // deep learner cuts
    if( useDeepLearnerCuts() )
    {
        cout << "Deep learning cuts (hard wired)" << endl;
    }
    // other cut parameters
    if( fNTel == 2 )
    {
        cout << ", size > " << fCut_Size_min << endl;
    }
    cout << "Fiducial area (camera) < " << fCut_CameraFiducialSize_max << " deg, ";
    cout << " stereo reconstruction quality: " << fCut_Chi2_min << " <= sChi2 <= " << fCut_Chi2_max;
    cout << " (direction method difference < " << fCut_Dir_max << ")" << endl;
    cout << "Energy reconstruction: ";
    cout << fCut_EChi2_min << " <= EChi2 <= " << fCut_EChi2_max;
    cout << ", " << fCut_dE_min << " < dE < " << fCut_dE_max;
    cout << ", " << fCut_Erec_min << " < Erec < " << fCut_Erec_max;
    cout << ", multiplicity ratio for E-reconstruction > " << fCut_EmultiplicityRatio_min;
    cout << endl;
    cout << "SizeSecondMax: " << fCut_SizeSecondMax_min << " < SizeSecondMax < " << fCut_SizeSecondMax_max << endl;
    cout << "Emission height: " << fCut_Emmission_min << " < Emission height [km] < " << fCut_Emmission_max;
    cout << endl;
    cout << "NImage cut: " << fCut_NImages_min << " <= Ntel <= " << fCut_NImages_max;
    cout << endl;
    if( fCut_DispNImages_min > 0 || fCut_DispNImages_max < 1000 )
    {
         cout << "Disp-NImage cut: " << fCut_DispNImages_min << " <= Ntel <= " << fCut_DispNImages_max;
         cout << endl;
    }
    if( bMCCuts )
    {
        cout << "MC cuts: " << fCut_CameraFiducialSize_MC_min << " < fiducial area (camera) < " << fCut_CameraFiducialSize_MC_max << " deg ";
        cout << endl;
    }
    if( fCut_ImgSelect.size() > 0 )
    {
        cout << "Allowed telescope combinations: ";
        for( unsigned int i = 0; i < fCut_ImgSelect.size(); i++ )
        {
            cout << i << ": " << fCut_ImgSelect[i];
            if( i < fCut_ImgSelect.size() - 1 )
            {
                cout << ", ";
            }
        }
        cout << endl;
    }
    if( fNTelTypeCut.size() > 0 )
    {
        for( unsigned int j = 0; j < fNTelTypeCut.size(); j++ )
        {
            fNTelTypeCut[j]->print();
        }
    }
    printEnergyDependentCuts();
    cout << "-----------------------------------------------------------------------------------------" << endl;
}


/////////////////////////////////////////////////////////////////////////////////
// require certain quality in stereo reconstruction
bool VGammaHadronCuts::applyStereoQualityCuts_Chi2( bool bCount )
{
    if( fData->getChi2() < fCut_Chi2_min || fData->getChi2() > fCut_Chi2_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eArrayChi2 );
        }
        return false;
    }
    // cut on differences between e.g. DISP methods and line intersection
    // method
    if( fData->getDirectionReconstructionDifference() > fCut_Dir_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eArrayChi2 );
        }
        return false;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////
// apply number of images cut
bool VGammaHadronCuts::applyStereoQualityCuts_NImages( bool bCount )
{
    if( fData->getNImages() < fCut_NImages_min || fData->getNImages() > fCut_NImages_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eNImages );
        }
        return false;
    }
    // disp nimages cut (applies only for disp stereo reconstruction)
    if( fData->DispNImages < fCut_DispNImages_min || fData->DispNImages > fCut_DispNImages_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eNImages );
        }
        return false;
    }
    return true;
}

/*
 * check quality of variables from lookup tables
 *
 * iEnergyReconstructionMethod == 99: apply no cuts
 *
 */
bool VGammaHadronCuts::applyStereoQualityCuts_TableVariables( bool bCount, unsigned int iEnergyReconstructionMethod )
{
    if( iEnergyReconstructionMethod != 99 )
    {
        // quality cut for MSCW/L reconstruction cuts
        if( fGammaHadronCutSelector % 10 < 1 && ( fData->MSCW < -50. || fData->MSCL < -50. ) )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eMSC_Quality );
            }
            return false;
        }
        // quality cut for MWR/MLR reconstruction cuts
        if( fGammaHadronCutSelector % 10 == 3 && ( fData->MWR < -50. || fData->MLR < -50. ) )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eMSC_Quality );
            }
            return false;
        }
        
        /////////////////////////////////////////////////////////////////////////////////
        // check energy reconstruction
        double iErec  = fData->getEnergy_TeV();
        
        if( iErec < 0. || iErec < fCut_Erec_min )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eErec );
            }
            return false;
        }
        if( iErec > fCut_Erec_max )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eErec );
            }
            return false;
        }
    }
    
    return true;
    
}

/*
 * apply quality cut on core position / distance to telescopes
 *
 */
bool VGammaHadronCuts::applyStereoQualityCuts_TelescopeCoreDistance( bool bCount )
{
    // closest distance between telescope and core
    double iR_min = 1.e10;
    // average distance to telescopes with images
    double iR     = 0.;
    double iNTR   = 0.;
    
    // loop over all images
    for( int i = 0; i < fData->getNImages(); i++ )
    {
        if( fData->getR()[i] > 0. )
        {
            iR += fData->getR()[i];
            iNTR++;
            if( fData->getR()[i] < iR_min )
            {
                iR_min = fData->getR()[i];
            }
        }
    }
    if( iNTR > 0. )
    {
        iR /= iNTR;
    }
    
    // apply cut on distance
    if( iR < fCut_AverageCoreDistanceToTelescopes_min || iR > fCut_AverageCoreDistanceToTelescopes_max )
    {
        if( fCut_AverageCoreDistanceToTelescopes_nimages < 0
                || fCut_AverageCoreDistanceToTelescopes_nimages == fData->getNImages() )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eCorePos );
            }
            return false;
        }
    }
    if( iR_min > fCut_MinimumCoreDistanceToTelescopes_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eCorePos );
        }
        return false;
    }
    return true;
}

bool VGammaHadronCuts::applyStereoQualityCuts_ImageSelection( bool bCount )
{
    /////////////////////////////////////////////////////////////////////////////////
    // apply image selection cuts (check which telescopes were used in the reconstruction)
    if( fCut_ImgSelect.size() > 0 )
    {
        if( fData->getImgSel() < fCut_ImgSelect.size() )
        {
            if( !fCut_ImgSelect[fData->getImgSel() ] )
            {
                if( bCount && fStats )
                {
                    fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                    fStats->updateCutCounter( VGammaHadronCutsStatistics::eLTrig );
                }
                return false;
            }
        }
        else
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eLTrig );
            }
            return false;
        }
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////
// apply cuts on second max
bool VGammaHadronCuts::applyStereoQualityCuts_SecondMaxCut( bool bCount )
{
    if( fData->SizeSecondMax < fCut_SizeSecondMax_min || fData->SizeSecondMax > fCut_SizeSecondMax_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eSizeSecondMax );
        }
        return false;
    }
    return true;
}

/*
  ensure event quality, reasonable output from the table variables, etc
*/
bool VGammaHadronCuts::applyStereoQualityCuts( unsigned int iEnergyReconstructionMethod, bool bCount, int iEntry, bool iIsOn )
{

    if( !applyStereoQualityCuts_Chi2( bCount ) )
    {
        return false;
    }
    if( !applyStereoQualityCuts_NImages( bCount ) )
    {
        return false;
    }
    if( !applyStereoQualityCuts_TableVariables( bCount, iEnergyReconstructionMethod ) )
    {
        return false;
    } 
    if( !applyStereoQualityCuts_TelescopeCoreDistance( bCount ) )
    {
        return false;
    }
    if( !applyStereoQualityCuts_ImageSelection( bCount ) )
    {
        return false;
    }
    if( !applyStereoQualityCuts_SecondMaxCut( bCount ) )
    {
        return false;
    } 
    
    /////////////////////////////////////////////////////////////////////////
    // apply cut selector from probability tree
    if( fGammaHadronCutSelector / 10 == 3 )
    {
        if( !applyProbabilityCut( iEntry , iIsOn ) )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eStereoQuality );
            }
            return false;
        }
    }
    
    return true;
}


/*!
  cuts apply in energies in linear scale!

  iEnergyReconstructionMethod == 100: return always true
*/
bool VGammaHadronCuts::applyEnergyReconstructionQualityCuts( unsigned int iEnergyReconstructionMethod, bool bCount )
{
    if( iEnergyReconstructionMethod == 100 )
    {
        return true;
    }
    
    double iErec = fData->getEnergy_TeV( );
    double iErecChi2 = fData->getEnergyChi2( );
    double idE = fData->getEnergyDelta( );
    
    if( iErecChi2 < fCut_EChi2_min || iErecChi2 > fCut_EChi2_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eEnergyRec );
        }
        return false;
    }
    if( iErec < fCut_Erec_min || iErec > fCut_Erec_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eEnergyRec );
        }
        return false;
    }
    if( idE < fCut_dE_min || idE > fCut_dE_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eEnergyRec );
        }
        return false;
    }
    
    ////////////////////////////////////////////////////////
    // test fraction of images used in energy reconstruction
    // (good events should use all images)
    unsigned int i_NErec = fData->getNImages_used_in_EnergyReconstruction();
    
    if( fData->NImages > 0 && ( double )i_NErec / ( double )fData->NImages < fCut_EmultiplicityRatio_min )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eEnergyRec );
        }
        return false;
    }
    
    
    return true;
}

/*

  apply gamma/hadron separation cut

*/

bool VGammaHadronCuts::isGamma( int i, bool bCount, bool fIsOn )
{

    /////////////////////////////////////////////////////////////////////////////
    // apply box cuts  (e.g. MSCW/MSCL or MWR/MLR)
    if( fGammaHadronCutSelector % 10 <= 3 )
    {
        if( fDebug )
        {
            cout << "VGammaHadronCuts::isGamma: applyStereoShapeCuts" << endl;
        }
        if( !applyStereoShapeCuts() )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eIsGamma );
            }
            return false;
        }
        // all cut selectors >= 10 are different
        if( fGammaHadronCutSelector < 10 )
        {
            return true;
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    // apply probability threshold cut (e.g. random forest cuts)
    if( fGammaHadronCutSelector / 10 == 1 || fGammaHadronCutSelector / 10 == 2 )
    {
        if( fDebug )
        {
            cout << "VGammaHadronCuts::isGamma: applyProbabilityCut" << endl;
        }
        if( !applyProbabilityCut( i, fIsOn ) )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eIsGamma );
            }
            return false;
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    // apply cut using TMVA reader
    else if( useTMVACuts() )
    {
        if( fDebug )
        {
            cout << "VGammaHadronCuts::isGamma: applyTMVACut" << endl;
        }
        if( !applyTMVACut( i ) )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eIsGamma );
            }
            return false;
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    // apply deep learner cuts 
    // (using data_DL tree for cuts)
    else if( useDeepLearnerCuts() )
    {
        if( !applyDeepLearnerCut() )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eIsGamma );
            }
            return false;
        }
    }
    
    return true;
}

/*

   use TMVA reader and apply cuts

*/
bool VGammaHadronCuts::applyTMVACut( int i )
{
    if( fDebug )
    {
        cout << "VGammaHadronCuts::applyTMVACut event " << i;
        cout << ", signal efficiency " << fTMVASignalEfficiency.size();
        cout << ", probability threshold/MVA cut " << fTMVA_MVACut.size();
        cout << " (" << fTMVAEvaluator << ")";
        cout << endl;
    }
    fTMVA_EvaluationResult = -99.;
    if( fTMVAEvaluator )
    {
        bool i_TMVA_Evaluation = fTMVAEvaluator->evaluate();
        fTMVA_EvaluationResult = fTMVAEvaluator->getTMVA_EvaluationResult();
        return i_TMVA_Evaluation;
    }
    
    return false;
}

/*
   apply deep learner cut
*/
bool VGammaHadronCuts::applyDeepLearnerCut()
{
    if( !fData->isDeepLearner() )
    {
        cout << "VGammaHadronCuts::applyDeepLearnerCut error: input data (mscw file) does not contain DL parameters" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }

    // cut on boolean gammaness
    // if( !fData->dl_isGamma )
    
    // cut on gammmaness
    // (hardwired)
    if( fData->dl_gammaness < 0.2 )
    {
        return false;
    }
    
    return true;
}

/*

  use cut selector from a friend tree to the given data tree

*/
bool VGammaHadronCuts::applyProbabilityCut( int i, bool fIsOn )
{
    // tree exists and has the requested entry
    if( fProbabilityCut_Tree && fProbabilityCut_Tree->GetEntry( i ) )
    {
        // check cut quality
        if( fProbabilityCut_QualityFlag > 0 )
        {
            if( fProbabilityCutRangeLower.size() != fProbabilityCutRangeUpper.size() )
            {
                cout << "Error in definitions of RF probablity ranges" << endl
                     << "RFCutLowerVals and  RFCutLowerVals have different numbers of entries in cut file" << endl;
                exit( -1 );
            }
            else
            {
                if( fProbabilityCut_ProbID < fProbabilityCut_NSelectors && fProbabilityCut_NSelectors < VANACUTS_PROBSELECTIONCUTS_MAX )
                {
                    for( unsigned int dex = 0; dex < fProbabilityCutRangeLower.size(); dex++ )
                    {
                        if( fIsOn && fProbabilityCut_SelectionCut[fProbabilityCut_ProbID] >= fProbabilityCutRangeLower[dex]
                                && fProbabilityCut_SelectionCut[fProbabilityCut_ProbID] <  fProbabilityCutRangeUpper[dex] )
                        {
                            return true;
                        }
                        if( !fIsOn && fProbabilityCut_SelectionCut[fProbabilityCut_ProbID] >= -fProbabilityCutRangeLower[dex]
                                && fProbabilityCut_SelectionCut[fProbabilityCut_ProbID] <  -fProbabilityCutRangeUpper[dex] )
                        {
                            return true;
                        }
                    }
                }
            }
        }
    }
    
    return false;
}

/*!

    apply mean reduced scaled cuts

*/
bool VGammaHadronCuts::applyMeanReducedScaledStereoShapeCuts()
{
    if( fData->MSCW > fCut_MSCW_max )
    {
        return false;
    }
    if( fData->MSCW < fCut_MSCW_min )
    {
        return false;
    }
    if( fData->MSCL > fCut_MSCL_max )
    {
        return false;
    }
    if( fData->MSCL < fCut_MSCL_min )
    {
        return false;
    }
    
    return true;
}

/*!

    apply mean scaled cuts

*/
bool VGammaHadronCuts::applyMeanScaledStereoShapeCuts()
{
    if( fData->MWR > fCut_MSW_max )
    {
        return false;
    }
    if( fData->MWR < fCut_MSW_min )
    {
        return false;
    }
    if( fData->MLR > fCut_MSL_max )
    {
        return false;
    }
    if( fData->MLR < fCut_MSL_min )
    {
        return false;
    }
    
    return true;
}

/*

   apply mean shape cuts

   (not scaled width, this calculation is MC independent)

*/
bool VGammaHadronCuts::applyMeanStereoShapeCuts()
{
    fMeanImageWidth = 0.;
    fMeanImageLength = 0.;
    fMeanImageDistance = 0.;
    int intel = 0;
    // loop over all telescopes and calculate mean values
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        if( fData->size[i] <= fCut_Size_min || fData->size[i] >= fCut_Size_max )
        {
            continue;
        }
        
        fMeanImageWidth += fData->width[i];
        fMeanImageLength += fData->length[i];
        fMeanImageDistance += fData->dist[i];
        intel++;
    }
    if( intel > 0 )
    {
        fMeanImageWidth      /= ( double )intel;
        fMeanImageLength     /= ( double )intel;
        fMeanImageDistance   /= ( double )intel;
    }
    else
    {
        return false;
    }
    // apply cuts
    if( fMeanImageDistance <= fCut_MeanImageDistance_min || fMeanImageDistance >= fCut_MeanImageDistance_max )
    {
        return false;
    }
    if( fMeanImageLength <= fCut_MeanImageLength_min || fMeanImageLength >= fCut_MeanImageLength_max )
    {
        return false;
    }
    if( fMeanImageWidth <= fCut_MeanImageWidth_min || fMeanImageWidth >= fCut_MeanImageWidth_max )
    {
        return false;
    }
    
    return true;
}

/*

   stereo shape cuts

*/
bool VGammaHadronCuts::applyStereoShapeCuts()
{
    /////////////////////////////////////////////////////////////////////////////
    // cut selector 2 is always true
    if( fGammaHadronCutSelector % 10 == 2 )
    {
        return true;
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // apply cuts
    /////////////////////////////////////////////////////////////////////////////
    // MSCW/MSCL cuts
    if( fGammaHadronCutSelector % 10 < 1 )
    {
        if( !applyMeanReducedScaledStereoShapeCuts() )
        {
            return false;
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    // mean width/length cuts
    // (not scaled width, this calculation is MC independent)
    else if( fGammaHadronCutSelector % 10 == 1 )
    {
        if( !applyMeanStereoShapeCuts() )
        {
            return false;
        }
    }
    /////////////////////////////////////////////////////////////////////////////
    // MWR/MLR cuts
    else if( fGammaHadronCutSelector % 10 == 3 )
    {
        if( !applyMeanScaledStereoShapeCuts() )
        {
            return false;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // emission height cuts
    if( fData->EmissionHeight > 0. && fData->EmissionHeight > fCut_Emmission_max )
    {
        return false;
    }
    if( fData->EmissionHeight > 0. && fData->EmissionHeight < fCut_Emmission_min )
    {
        return false;
    }
    
    return true;
}


/*
  initialize gamma/hadron cuts

  (see description of cuts at the top of this file)
*/
void VGammaHadronCuts::initializeCuts( int irun, string iFile )
{
    // probability cuts
    if( fGammaHadronCutSelector / 10 >= 1 && fGammaHadronCutSelector / 10 <= 3 )
    {
        if( irun > 0 )
        {
            initProbabilityCuts( irun );
        }
        else if( iFile.size() > 0 )
        {
            initProbabilityCuts( iFile );
        }
        else
        {
            cout << "VGammaHadronCuts::initializeCuts: failed setting probability cuts for " << irun << " " << iFile << endl;
            cout << "exiting..." << endl;
            exit( -1 );
        }
    }
    // TMVA cuts
    else if( useTMVACuts() )
    {
        if( !initTMVAEvaluator( fTMVAWeightFile, fTMVAWeightFileIndex_Emin, fTMVAWeightFileIndex_Emax,
                                fTMVAWeightFileIndex_Zmin, fTMVAWeightFileIndex_Zmax, fTMVAEnergyStepSize ) )
        {
            cout << "VGammaHadronCuts::initializeCuts: failed setting TMVA reader for " << fTMVAWeightFile;
            cout << "(" << fTMVAWeightFileIndex_Emin << "," << fTMVAWeightFileIndex_Emax << ")" << endl;
            cout << "(" << fTMVAWeightFileIndex_Zmin << "," << fTMVAWeightFileIndex_Zmax << ")" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
}

/*
 * initialize TMVA evaluator
 *
 * take into account number of energy and zenith bins
*/
bool VGammaHadronCuts::initTMVAEvaluator( string iTMVAFile, 
        unsigned int iTMVAWeightFileIndex_Emin, unsigned int iTMVAWeightFileIndex_Emax,
        unsigned int iTMVAWeightFileIndex_Zmin, unsigned int iTMVAWeightFileIndex_Zmax,
        double iTMVAEnergy_StepSize )
{
    TDirectory* cDir = gDirectory;
    
    fTMVAEvaluator = new VTMVAEvaluator();
    
    // turn off theta2 optimization except for TMVA direction cut selector flags
    if( fDirectionCutSelector == 3 || fDirectionCutSelector == 4 || fDirectionCutSelector == 5 )
    {
        fTMVAEvaluator->setOptimizeAngularContainment( true );
    }
    else
    {
        fTMVAEvaluator->setOptimizeAngularContainment( false );
    }
    
    fTMVAEvaluator->setDebug( fDebug );
    // smoothing of MVA values
    fTMVAEvaluator->setSmoothAndInterpolateMVAValues( true );
    // set parameters for optimal MVA cut value search
    // (always assume an alpha value of 0.2)
    if( fTMVAOptimizeSignalEfficiencyParticleNumberFile.size() > 0. )
    {
        fTMVAEvaluator->setSensitivityOptimizationParameters( fTMVAOptimizeSignalEfficiencySignificance_Min,
                fTMVAOptimizeSignalEfficiencySignalEvents_Min,
                fTMVAOptimizeSignalEfficiencyObservationTime_h,
                1. / 5. );
        fTMVAEvaluator->setSensitivityOptimizationFixedSignalEfficiency( fTMVAFixedSignalEfficiencyMax );
        fTMVAEvaluator->setParticleNumberFile( fTMVAOptimizeSignalEfficiencyParticleNumberFile,
                                               fTMVAParticleNumberFile_Conversion_Rate_to_seconds );
        fTMVAEvaluator->setSensitivityOptimizationSourceStrength( fTMVAMinSourceStrength );
    }
    // set a constant signal efficiency
    else if( fTMVASignalEfficiency.size() > 0 )
    {
        fTMVAEvaluator->setSignalEfficiency( fTMVASignalEfficiency );
    }
    // set a fixed probability threshold or (for TMVA) a fixed MVA cut value
    else if( fTMVA_MVACut.size() > 0 )
    {
        fTMVAEvaluator->setTMVACutValue( fTMVA_MVACut );
    }
    else if( fTMVA_MVACutGraphFileName.size() == 0 )
    {
        cout << "VGammaHadronCuts::initTMVAEvaluator error: unclear TMVA cut settings" << endl;
        cout << "\t fTMVASignalEfficiency: " << fTMVASignalEfficiency.size() << endl;
        cout << "\t fTMVAProbabilityThreshold: " << fTMVA_MVACut.size() << endl;
        cout << "exiting... " << endl;
        exit( EXIT_FAILURE );
    }
    fTMVAEvaluator->setTMVAMethod( fTMVA_MVAMethod, fTMVA_MVAMethodCounter );
    fTMVAEvaluator->setTMVAAngularContainmentThetaFixedMinRadius( fTMVAFixedThetaCutMin );
    // read MVA weight files; set MVA cut values (e.g. find optimal values)
    if( !fTMVAEvaluator->initializeWeightFiles( iTMVAFile,
                                                iTMVAWeightFileIndex_Emin, iTMVAWeightFileIndex_Emax,
                                                iTMVAWeightFileIndex_Zmin, iTMVAWeightFileIndex_Zmax,
                                                iTMVAEnergy_StepSize, fInstrumentEpoch,
                                                "UseInterpolatedCounts",
                                                fCutID ) )
    {
        cout << "VGammaHadronCuts::initTMVAEvaluator: error while initializing TMVA weight files" << endl;
        cout << "exiting... " << endl;
        exit( EXIT_FAILURE );
    }
    // set MVA cut files from a list of graphs in a root file
    if( fTMVA_MVACutGraphFileName.size() > 0 )
    {
        fTMVAEvaluator->setTMVACutValueFromGraph( fTMVA_MVACutGraphFileName,
                                                  fTMVA_MVACutGraphSmoothing, 
                                                  fTMVA_MVACutGraphSmoothingMax,
                                                  fTMVA_MVACutGraphConstantCutEnergy_TeV, false );
        fMVACutGraphs = fTMVAEvaluator->getTMVACutValueGraphs();
    }

    
    fTMVAEvaluatorResults = fTMVAEvaluator->getTMVAEvaluatorResults();
    fTMVAEvaluator->printSignalEfficiency();
    fTMVAEvaluator->printAngularContainmentRadius();
    
    if( fDirectionCutSelector == 3 )
    {
        fTMVAEvaluator->setIgnoreTheta2Cut( false );
    }
    else
    {
        fTMVAEvaluator->setIgnoreTheta2Cut( true );
    }
    fEnergyDependentCut[ "TMVABoxCut_Theta2_max" ] = fTMVAEvaluator->getOptimalTheta2Cut_Graph();
    if( getTheta2Cut_TMVA_max() )
    {
        cout << "VGammaHadronCuts::initTMVAEvaluator: found theta2_max graph from TMVA" << endl;
        getTheta2Cut_TMVA_max()->SetName( "TMVAtheta2" );
    }
    
    if( cDir )
    {
        cDir->cd();
    }
    
    return !fTMVAEvaluator->IsZombie();
}

bool VGammaHadronCuts::setDataTree( CData* idata )
{
    fData = idata;
    
    if( !fData )
    {
        return false;
    }
    
    if( fTMVAEvaluator )
    {
        fTMVAEvaluator->initializeDataStrutures( fData );
    }
    return true;
}


bool VGammaHadronCuts::initProbabilityCuts( int irun )
{
    ostringstream iFile;
    iFile << fDataDirectory << "/" << irun << ".mscw.rf.root";
    
    return initProbabilityCuts( iFile.str() );
}


bool VGammaHadronCuts::initProbabilityCuts( string iFile )
{
    TDirectory* cDir = gDirectory;
    
    fProbabilityCut_File = new TFile( iFile.c_str() );
    if( fProbabilityCut_File->IsZombie() )
    {
        cout << "Error while opening file with probability cuts: " << iFile << endl;
        exit( 0 );
    }
    cout << "\t opening file with probability cuts: " << fProbabilityCut_File->GetName() << endl;
    
    fProbabilityCut_Tree = ( TTree* )gDirectory->Get( "rf" );
    if( !fProbabilityCut_Tree )
    {
        cout << "Error: could not find tree with probability cuts" << endl;
        exit( 0 );
    }
    fProbabilityCut_Tree->SetBranchAddress( "cut", &fProbabilityCut_QualityFlag );
    if( fProbabilityCut_Tree->GetBranchStatus( "Ng" ) )
    {
        fProbabilityCut_Tree->SetBranchAddress( "Ng", &fProbabilityCut_NSelectors );
    }
    else
    {
        fProbabilityCut_NSelectors = 2;
    }
    fProbabilityCut_Tree->SetBranchAddress( "g", fProbabilityCut_SelectionCut );
    
    if( cDir )
    {
        cDir->cd();
    }
    
    return true;
}


/*

   check that event is inside a certain region in the camera

*/
bool VGammaHadronCuts::applyInsideFiducialAreaCut( bool bCount )
{

    return applyInsideFiducialAreaCut( fData->getXoff(), fData->getYoff(), bCount );
}

bool VGammaHadronCuts::applyInsideFiducialAreaCut( float Xoff, float Yoff, bool bCount )
{
    double xy = Xoff * Xoff + Yoff * Yoff;
    
    if( xy > fCut_CameraFiducialSize_max * fCut_CameraFiducialSize_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eXYoff );
        }
        return false;
    }
    
    if( fCut_CameraFiducialSize_min >= 0. )
    {
        if( xy < fCut_CameraFiducialSize_min * fCut_CameraFiducialSize_min )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eXYoff );
            }
            return false;
        }
    }
    
    return true;
}

/*

    cut on Monte Carlo source position
    (camera coordinates)

*/
bool VGammaHadronCuts::applyMCXYoffCut( double xoff, double yoff, bool bCount )
{
    if( !fData->isMC() )
    {
        return true;
    }
    
    double xy = xoff * xoff + yoff * yoff;
    
    if( xy > fCut_CameraFiducialSize_MC_max * fCut_CameraFiducialSize_MC_max )
    {
        if( bCount && fStats )
        {
            fStats->updateCutCounter( VGammaHadronCutsStatistics::eMC_XYoff );
        }
        return false;
    }
    
    if( fCut_CameraFiducialSize_MC_min >= 0. )
    {
        if( xy < fCut_CameraFiducialSize_MC_min * fCut_CameraFiducialSize_MC_min )
        {
            if( bCount && fStats )
            {
                fStats->updateCutCounter( VGammaHadronCutsStatistics::eMC_XYoff );
            }
            return false;
        }
    }
    
    return true;
}

/*

   check telescope type (e.g. remove all LSTs)

   returns true if current event fails the fulfil the conditions

*/
bool VGammaHadronCuts::applyTelTypeTest( bool bCount )
{
    bool icut = false;
    
    if( fNTelTypeCut.size() == 0 )
    {
        return true;
    }
    
    for( unsigned int i = 0; i < fNTelTypeCut.size(); i++ )
    {
        // test: true means this event has enough telescopes
        icut = ( icut || fNTelTypeCut[i]->test( fData ) );
    }
    
    if( bCount && fStats && !icut )
    {
        fStats->updateCutCounter( VGammaHadronCutsStatistics::eTelType );
    }
    
    return icut;
}


/*
   apply cut on event direction (theta2 cut)

  * function is called for effective areas calculations only (MC)
  * the cut might be energy dependent (if a function is defined)

  x0, y0:   calculate theta2 relative to these points (-99999. if relative to MCx/yoff)

*/
bool VGammaHadronCuts::applyDirectionCuts( bool bCount, double x0, double y0 )
{
    double theta2 = 0.;
    
    // define reference direction
    // (default is simulated direction (MCxoff/MCyoff)
    if( x0 < -99990. )
    {
        if( fData->isMC() )
        {
            x0 = fData->MCxoff;
        }
        else
        {
            x0 = 0.;
        }
    }
    if( y0 < -99990. )
    {
        if( fData->isMC() )
        {
            y0 = fData->MCyoff;
        }
        else
        {
            y0 = 0.;
        }
    }
    
    // calculate theta2
    theta2 = ( fData->getXoff() - x0 ) * ( fData->getXoff() - x0 ) + ( fData->getYoff() - y0 ) * ( fData->getYoff() - y0 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // fetch theta2 cut (max) (might be energy dependent)
    double i_theta2_cut_max = getTheta2Cut_max( fData->getEnergy_TeV( ) );
    
    // direction cut (ring around center of camera)
    if( theta2 < i_theta2_cut_max && theta2 > fCut_Theta2_min )
    {
        return true;
    }
    
    if( bCount && fStats )
    {
        fStats->updateCutCounter( VGammaHadronCutsStatistics::eDirection );
    }
    
    return false;
}

double VGammaHadronCuts::getEnergyDependentCut( double energy_TeV, TGraph* iG, bool bUseEvalue, bool bMaxCut )
{
    double i_cut_value = 99999;
    if( !bMaxCut )
    {
        i_cut_value = -1;
    }
    if( !iG )
    {
        return i_cut_value;
    }
    // use TGraph::Eval
    if( bUseEvalue )
    {
        i_cut_value  = iG->Eval( energy_TeV );
    }
    // return value at this point without interpolation
    else
    {
        for( int i = 0; i < iG->GetN() - 1; i++ )
        {
            if( energy_TeV >= iG->GetX()[i] && energy_TeV < iG->GetX()[i + 1 ] )
            {
                i_cut_value = iG->GetY()[i + 1 ];
            }
        }
    }
    
    // for e outside of graph range, return edge values
    if( iG->GetN() > 0 && iG->GetX() && iG->GetY() )
    {
        if( energy_TeV < iG->GetX()[0] )
        {
            i_cut_value = iG->GetY()[0];
        }
        if( energy_TeV > iG->GetX()[iG->GetN() - 1] )
        {
            i_cut_value = iG->GetY()[iG->GetN() - 1];
        }
    }
    return i_cut_value;
}

/*
   fetch theta2 cut (might be energy dependent)

   e      :   [TeV] energy (linear)
*/
double VGammaHadronCuts::getTheta2Cut_max( double e )
{
    double theta_cut_max = -1.;
    // theta2 at fAngRes_FixedAboveEnergy_TeV
    double theta_cut_fixed = -1.;
    
    //////////////////////////////////////////////
    // energy independent theta2 cut
    //////////////////////////////////////////////
    if( fDirectionCutSelector == 0 )
    {
        theta_cut_max = TMath::Sqrt( fCut_Theta2_max );   // will be squared later
        theta_cut_fixed = theta_cut_max;
    }
    //////////////////////////////////////////////
    // energy dependent theta2 cut
    //////////////////////////////////////////////
    if( e > 0. )
    {
        // all graphs/functions are vs log10(energy [TeV])
        e = log10( e );
        
        /////////////////////////////////////////////
        // use a function to get the angular resolution
        if( fDirectionCutSelector == 1 && fF1AngRes )
        {
            // energy outside of functions range:, return edge values
            if( e < fF1AngRes->GetXmin() )
            {
                e = fF1AngRes->GetXmin();
            }
            else if( e > fF1AngRes->GetXmax() )
            {
                e = fF1AngRes->GetXmax();
            }
            
            // get angular resolution and apply scaling factor
            theta_cut_max  = fF1AngRes->Eval( e );
            if( fAngRes_FixedAboveEnergy_TeV > 0. )
            {
                theta_cut_fixed = fF1AngRes->Eval( log10( fAngRes_FixedAboveEnergy_TeV ) );
            }
        }
        /////////////////////////////////////////////
        // use IRF graph for angular resolution
        else if( ( fDirectionCutSelector == 1 || fDirectionCutSelector == 2 ) && getTheta2Cut_IRF_Max() )
        {
            // get theta2 cut
            theta_cut_max  = getEnergyDependentCut( e, getTheta2Cut_IRF_Max(), true );
            if( fAngRes_FixedAboveEnergy_TeV > 0. )
            {
                theta_cut_fixed = getEnergyDependentCut( log10( fAngRes_FixedAboveEnergy_TeV ), 
                                                         getTheta2Cut_IRF_Max(), 
                                                         true );
            }
        }
        /////////////////////////////////////////////
        // use TMVA determined cut
        else if( fDirectionCutSelector == 3 )
        {
            theta_cut_max = -1.;
            theta_cut_fixed = -1.;
        }
        // optimal theta2 cut
        else if( fDirectionCutSelector == 4 && fTMVAEvaluator )
        {
            theta_cut_max = fTMVAEvaluator->getOptimalTheta2Cut( e );
            if( theta_cut_max > 0. )
            {
                theta_cut_max = TMath::Sqrt( theta_cut_max );
            }
            else
            {
                theta_cut_max = 0.;
            }
            if( fAngRes_FixedAboveEnergy_TeV > 0. )
            {
                theta_cut_fixed = fTMVAEvaluator->getOptimalTheta2Cut( log10( fAngRes_FixedAboveEnergy_TeV ) );
            }
        }
        /////////////////////////////////////////////
        // use a graph with theta2 cuts
        else if( fDirectionCutSelector == 5 && getTheta2Cut_TMVA_max() )
        {
            theta_cut_max = getEnergyDependentCut( e, getTheta2Cut_TMVA_max(), true, true );
            if( fAngRes_FixedAboveEnergy_TeV > 0. )
            {
                theta_cut_fixed = getEnergyDependentCut( log10(fAngRes_FixedAboveEnergy_TeV),
                                                         getTheta2Cut_TMVA_max(), true, true );
            }
            if( theta_cut_max < 0. )
            {
                theta_cut_max = 0.;
            }
            // set minimum theta2 cut
            if( theta_cut_max < fTMVAFixedThetaCutMin )
            {
                theta_cut_max = fTMVAFixedThetaCutMin;
            }
        }
    }
    
    // apply scale factors
    theta_cut_max *= fAngRes_ScalingFactor;

    // check if a fixed theta should be applied above
    // a certain energy
    if( fAngRes_FixedAboveEnergy_TeV > 0. && e > log10(fAngRes_FixedAboveEnergy_TeV) && theta_cut_fixed > 0. )
    {
        theta_cut_max = theta_cut_fixed;
    }

    // check if theta2 is below/above absolute min/max
    if( theta_cut_max < fAngRes_AbsoluteMinimum )
    {
        return fAngRes_AbsoluteMinimum * fAngRes_AbsoluteMinimum;
    }
    if( theta_cut_max > fAngRes_AbsoluteMaximum )
    {
        return fAngRes_AbsoluteMaximum * fAngRes_AbsoluteMaximum;
    }
    return theta_cut_max * theta_cut_max;
}

/*

   read angular resolution from root file

*/

bool VGammaHadronCuts::initAngularResolutionFile()
{

    // open angular resolution file
    fFileAngRes = new TFile( fFileNameAngRes.c_str() );
    if( fFileAngRes->IsZombie() )
    {
        // try it at the default location
        string iTempFileL = VGlobalRunParameter::getDirectory_EVNDISPAnaData();
        iTempFileL += "/GammaHadronCutFiles/" + fFileNameAngRes;
        fFileAngRes = new TFile( iTempFileL.c_str() );
        if( fFileAngRes->IsZombie() )
        {
            cout << "VGammaHadronCuts::initAngularResolutionFile: error open angular resolution file: " << fFileNameAngRes << endl;
            cout << " (tried also " << iTempFileL << ")" << endl;
            return false;
        }
    }
    
    /////////////////////////
    // get energy dependent theta values (probably angular resolution)
    
    // read angular resolution as function (TF1) from a fit
    if( fF1AngResName != "IRF" && fF1AngResName != "GRAPH" )
    {
        // get energy dependent theta values (probably angular resolution)
        if( fFileAngRes->Get( fF1AngResName.c_str() ) )
        {
            fF1AngRes = ( TF1* )fFileAngRes->Get( fF1AngResName.c_str() )->Clone();
            cout << "VGammaHadronCuts::initAngularResolutionFile: read angular resolution fit function from file " << fF1AngResName << " : " << endl;
            fF1AngRes->Print();
        }
        else
        {
            cout << "VGammaHadronCuts::initAngularResolutionFile: error finding angular resolution function with name ";
            cout << fF1AngResName << " in " << fFileAngRes->GetName() << endl;
            return false;
        }
    }
    // read angular resolution from instrument response function tree
    else if( fF1AngResName == "IRF" )
    {
        char iTreeName[200];
        if( getAngularResolutionContainmentRadius() - 68 != 0 )
        {
            sprintf( iTreeName, "t_angular_resolution_%03dp", getAngularResolutionContainmentRadius() );

        }
        else
        {
           sprintf( iTreeName, "t_angular_resolution" );
        }
        
        cout << "VGammaHadronCuts::initAngularResolutionFile: reading angular resolution graph from file (" << iTreeName << "):" << endl;
        cout << fFileAngRes->GetName() << endl;
        TTree* t = ( TTree* )fFileAngRes->Get( iTreeName );
        if( !t )
        {
            cout << "VGammaHadronCuts::initAngularResolutionFile: error finding tree with instrument response function for angular resolution " << endl;
            return false;
        }
        // TODO: this should go in the future
        // (currently there should be only one entry in the angular resolution tree, which is the case in the CTA analysis)
        if( t->GetEntries() != 1 )
        {
            cout << "VGammaHadronCuts::initAngularResolutionFile error: invalid number of entries in instrument response function tree ";
            cout << " (should be 1, is " << t->GetEntries() << " )" << endl;
            cout << "   (as we do not know which tree entry to select...)" << endl;
            return false;
        }
        // read the tree data and clone the angular resolution graph
        VInstrumentResponseFunctionData* c = 0;
        TBranch* br = t->GetBranch( "IRF" );
        br->SetAddress( &c );
        t->GetEntry( 0 );
        if( c && c->fResolutionGraph.size() > VInstrumentResponseFunctionData::E_DIFF + 1
                && c->fResolutionGraph[VInstrumentResponseFunctionData::E_DIFF] )
        {
            TGraph* iG = ( TGraph* )c->fResolutionGraph[VInstrumentResponseFunctionData::E_DIFF]->Clone();
            if( iG )
            {
                iG->SetName( "IRFAngRes" );
                fEnergyDependentCut[ "IRFAngRes" ] = iG;
            }
            else
            {
                cout << "VGammaHadronCuts::initAngularResolutionFile error: reading angular resolution graph from ";
                cout << fF1AngResName << endl;
                return false;
            }
        }
    }
    // read a single graph with angular resolution from file
    else if( fF1AngResName == "GRAPH" )
    {
        TGraphErrors* g = ( TGraphErrors* )fFileAngRes->Get( "gAngularResolution" );
        if( !g )
        {
            cout << "VGammaHadronCuts::initAngularResolutionFile error: reading angular resolution graph from " << fF1AngResName << endl;
            return false;
        }
        setIRFGraph( g );
    }
    if( fFileAngRes )
    {
        fFileAngRes->Close();
    }
    
    return true;
}

bool VGammaHadronCuts::setIRFGraph( TGraphErrors* g )
{
    if( !g )
    {
        cout << "VGammaHadronCuts::setIRFGraph warning: IRF pointer is zero" << endl;
        return false;
    }
    
    
    TGraphErrors* iG = new TGraphErrors( 1 );
    double x = 0.;
    double y = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        iG->SetPoint( i, x, y );
        iG->SetPointError( i, g->GetErrorX( i ), g->GetErrorY( i ) );
    }
    iG->SetName( "IRFAngRes" );
    fEnergyDependentCut[ "IRFAngRes" ] = iG;
    
    // print results
    cout << "replaced IRF graph for direction cut" << endl;
    printDirectionCuts();
    
    return true;
}


/*!

     update statistics

*/
void VGammaHadronCuts::newEvent( bool iFillTree )
{
    // fill previous event
    if( iFillTree && fStats->getCounterValue( VGammaHadronCutsStatistics::eTot > 0 ) )
    {
        fStats->fill();
    }
    
    fStats->updateCutCounter( VGammaHadronCutsStatistics::eTot );
}


double VGammaHadronCuts::getProbabilityCutAlpha( bool fIsOn )
{
    if( fProbabilityCutRangeLower.size() != fProbabilityCutRangeUpper.size() )
    {
        cout << "Error in definitions of RF probablity ranges" << endl;
        cout << "RFCutLowerVals and RFCutLowerVals have different numbers of entries in cut file" << endl;
        exit( -1 );
    }
    //////////////////////////////////////////////////////////////
    // return 1 if probability cuts are not set
    if( fProbabilityCutRangeLower.size() == 0 && fProbabilityCutRangeUpper.size() == 0 )
    {
        return 1.;
    }
    //////////////////////////////////////////////////////////////
    
    double on_size = 0;
    double off_size = 0;
    for( unsigned int i = 0; i < fProbabilityCutRangeLower.size(); i++ )
    {
        if( fProbabilityCutRangeLower[i] >= 0 && fProbabilityCutRangeUpper[i] >= 0 )
        {
            on_size = on_size + ( fProbabilityCutRangeUpper[i] - fProbabilityCutRangeLower[i] );
        }
        else if( fProbabilityCutRangeLower[i] <= 0 && fProbabilityCutRangeUpper[i] <= 0 )
        {
            off_size = off_size - ( fProbabilityCutRangeUpper[i] - fProbabilityCutRangeLower[i] );
        }
        else
        {
            cout << "Error in definitions of RF probablity ranges" << endl;
            cout << "One pair of RFCutLowerVals and RFCutLowerVals values have oposite sign" << endl;
        }
    }
    
    if( fIsOn )
    {
        return on_size;
    }
    else
    {
        return off_size;
    }
    
}

/*
 * check if this cut should be used
 * in IRF or effective area calculation
 */
bool VGammaHadronCuts::useThisCut( CData *c )
{
    if( !c ) return false;

    if( fCutCharacteristicsMCAZ < -998. ) return true;

    if( TMath::Abs( c->MCaz - fCutCharacteristicsMCAZ )
       < fCutCharacteristicsMCAZ_tolerance )
    {
        return true;
    }
    if( TMath::Abs( c->MCaz - fCutCharacteristicsMCAZ - 360. )
       < fCutCharacteristicsMCAZ_tolerance )
    {
        return true;
    }

    return false;
}

void VGammaHadronCuts::terminate( bool iShort,
                                  string iObjectName )
{
    SetName( iObjectName.c_str() );
    
    if( fStats->getDataTree() && !iShort )
    {
        fStats->terminate();
        fStats->getDataTree()->Write();
    }
    if( fTMVAEvaluatorResults )
    {
        fTMVAEvaluatorResults->Write();
    }
    else
    {
        cout << "No TMVAEvaluator Results." << endl;
    }
    
    Write();
}

void VGammaHadronCuts::printSignalEfficiency()
{
    if( fTMVASignalEfficiency.size() == 0 )
    {
        cout << "no signal efficiency set, already optimized?" << endl;
        return;
    }
    
    map< unsigned int, double >::iterator iIter;
    for( iIter = fTMVASignalEfficiency.begin(); iIter != fTMVASignalEfficiency.end(); iIter++ )
    {
        cout << "signal efficiency for energy bin " << iIter->first << ": ";
        cout << iIter->second << endl;
    }
}

void VGammaHadronCuts::printTMVA_MVACut()
{
    if( fTMVA_MVACut.size() == 0 )
    {
        cout << "no MVA cut set" << endl;
        return;
    }
    
    map< unsigned int, double >::iterator iIter;
    for( iIter = fTMVA_MVACut.begin(); iIter != fTMVA_MVACut.end(); ++iIter )
    {
        cout << "MVA cut for energy/zenith bin " << iIter->first << ": ";
        cout << iIter->second << endl;
    }
}

void VGammaHadronCuts::printEnergyDependentCuts()
{
    map< string, TGraph* >::iterator it;
    for( it = fEnergyDependentCut.begin(); it != fEnergyDependentCut.end(); ++it )
    {
        if( it->second )
        {
            double x = 0.;
            double y = 0.;
            cout << "Energy dependent cut variable " << it->first << ": " << endl;
            for( int i = 0; i < it->second->GetN(); i++ )
            {
                it->second->GetPoint( i, x, y );
                cout <<  "[" << x << ",>" << y << "] ";
            }
            cout << endl;
        }
    }
}

string VGammaHadronCuts::getTelToAnalyzeString()
{
    stringstream iTemp;
    sort( fTelToAnalyze.begin(), fTelToAnalyze.end() );
    for( unsigned int i = 0; i < fTelToAnalyze.size(); i++ )
    {
        iTemp << fTelToAnalyze[i] + 1;
    }
    
    return iTemp.str();
}


TGraph* VGammaHadronCuts::getEnergyDependentCut( string iCutName )
{
    if( fEnergyDependentCut.find( iCutName ) != fEnergyDependentCut.end() )
    {
        return fEnergyDependentCut[iCutName];
    }
    
    return 0;
}

bool VGammaHadronCuts::getEnergyDependentCutFromFile( string iFileName, string iVariable )
{
    string iTemp = "g" + iVariable;
    if( gSystem->AccessPathName( iFileName.c_str() ) )
    {
        string iTempFileL = VGlobalRunParameter::getDirectory_EVNDISPAnaData();
        iFileName = iTempFileL + "/GammaHadronCutFiles/" + iFileName;
    }
    TFile* i_f = new TFile( iFileName.c_str(), "READ" );
    if( i_f->IsZombie() )
    {
        cout << "VGammaHadronCuts::getEnergyDependentCutFromFile error opening file with energy dependent cut: " << iFileName << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    TGraph* g = ( TGraph* )i_f->Get( iTemp.c_str() );
    if( !g )
    {
        cout << "VGammaHadronCuts::getEnergyDependentCutFromFile error while reading energy dependent cut: " << iTemp << endl;
        cout << "\t tried to read from " << iFileName << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else
    {
        TGraph* iG = ( TGraph* )g->Clone();
        iG->SetName( iVariable.c_str() );
        fEnergyDependentCut[iVariable] = iG;
    }
    i_f->Close();
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cut on number of images per telescope type depend
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VNTelTypeCut::VNTelTypeCut()
{
    fNTelType_min = 0;
}


/*
 *   count number of images per telescope type(s)
 *   test if this is above the given threshold cut
 *
 */
bool VNTelTypeCut::test( CData* c )
{
    if( !c )
    {
        return false;
    }
    
    unsigned int ntel_type = 0;
    // identification by teltype ID (preferred)
    if( fTelTypeID.size() > 0 )
    {
        for( unsigned int i = 0; i < fTelTypeID.size(); i++ )
        {
            for( int t = 0; t < c->NTtype; t++ )
            {
                if( fTelTypeID[i] == c->TtypeID[t] )
                {
                    ntel_type += c->NImages_Ttype[t];
                }
            }
        }
    }
    // identification by teltype counter (dangerous, as sequence might change from
    // one array type to another)
    else
    {
        for( unsigned int i = 0; i < fTelType_counter.size(); i++ )
        {
            if( ( int )fTelType_counter[i] < c->NTtype )
            {
                ntel_type += c->NImages_Ttype[fTelType_counter[i]];
            }
        }
    }
    // OBS! >=
    if( ntel_type >= fNTelType_min && fNTelType_min > 0 )
    {
        return true;
    }
    
    return false;
}

/*
 * make sure that the list of telescope IDs is unique 
 * remove duplicates
 */
void VNTelTypeCut::purgeTelTypeIDs()
{
       sort( fTelTypeID.begin(), fTelTypeID.end() );
       fTelTypeID.erase( unique( fTelTypeID.begin(), fTelTypeID.end() ), fTelTypeID.end() );
}

void VNTelTypeCut::print()
{
    cout << "telescope type cut: mintel >= " << fNTelType_min << " for type(s) ";
    if( fTelTypeID.size() > 0 )
    {
        for( unsigned int i = 0; i < fTelTypeID.size(); i++ )
        {
            cout << fTelTypeID[i] << " ";
        }
    }
    else
    {
        for( unsigned int i = 0; i < fTelType_counter.size(); i++ )
        {
            cout << fTelType_counter[i] << " ";
        }
    }
    cout << endl;
}

