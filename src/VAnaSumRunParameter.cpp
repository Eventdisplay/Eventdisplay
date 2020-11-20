/*! \class VAnaSumRunParameter
 *
 *
*/

#include "VAnaSumRunParameter.h"


/*

        definition of data class for anasum run parameters

*/
VAnaSumRunParameterDataClass::VAnaSumRunParameterDataClass()
{
    fEventDisplayVersion = "";
    
    fRunOn = 0;
    fRunOnFileName = "";
    fRunOff = 0;
    fRunOffFileName = "";
    
    fMJDOn = 0.;
    fMJDOff = 0.;
    
    fTarget = "";
    fTargetRAJ2000 = 0.;
    fTargetDecJ2000 = -90.;
    fTargetRA = 0.;
    fTargetDec = 0.;
    
    fPairOffset = 0.;
    fOff_Target         = "";
    fOff_TargetRAJ2000  = 0.;
    fOff_TargetDecJ2000 = -90.;
    fOff_WobbleNorth    = 0.;
    fOff_WobbleWest     = 0.;
    
    fWobbleNorth = 0.;                            // [deg]
    fWobbleWest = 0.;                             // [deg]
    fWobbleNorthMod = 0.;                         // [deg]
    fWobbleWestMod = 0.;                          // [deg]
    
    fSkyMapCentreNorth = 0.;
    fSkyMapCentreWest = 0.;
    fSkyMapCentreRAJ2000 = 0.;
    fSkyMapCentreDecJ2000 = 0.;
    
    fTargetShiftNorth = 0.;
    fTargetShiftWest = 0.;
    fTargetShiftRAJ2000 = 0.;
    fTargetShiftDecJ2000 = 0.;
    
    fNTel = 4;
    fTelToAna = "";
    fMaxTelID = fNTel;
    
    fBackgroundModel = 0;
    fSourceRadius = 0.;                           // actually radius^2
    fmaxradius = 0.;                              // maximum accepted distance from camera center [deg]
    
    fCutFile = "";
    
    fAcceptanceFile = "";
    fRunWiseRadialAcceptance = false;
    fListOfExclusionRegions = 0;
    
    fEffectiveAreaFile = "";                      // file with effective areas, use NOFILE if not avaible
    
    // ON/OFF MODEL
    fOO_alpha = 0.;
    
    // RING BACKGROUND MODEL
    fRM_RingRadius = 0.;                          // ring radius [deg]
    fRM_RingWidth = 0.;                           // ring width [deg]
    
    // REFLECTED REGION MODEL
    fRE_distanceSourceOff = 0.2;                  // minimal distance of off source regions in number of background regions from the source region
    fRE_nMinoffsource = 3;                        // minmum number of off source regions (default 3)
    fRE_nMaxoffsource = 7;                        // maximum number of off source regions (default 7)
    
    // TEMPLATE MODEL
    fTE_mscw_min = 0.;
    fTE_mscw_max = 0.;
    fTE_mscl_min = 0.;
    fTE_mscl_max = 0.;
}

VAnaSumRunParameterDataClass::~VAnaSumRunParameterDataClass()
{
/*        if( fListOfExclusionRegions )
        {
             delete fListOfExclusionRegions;
        } */
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


VAnaSumRunParameter::VAnaSumRunParameter()
{
    // default version number (important for reading of run parameters from root file)
    fVersion = 6;
    
    // bin size for sky maps [deg]
    //   must be the same for all runs (add up sky maps)
    fSkyMapBinSize   = 0.01;
    fSkyMapBinSizeUC = 0.05;
    fSkyMapSizeXmin = -2.;
    fSkyMapSizeXmax =  2.;
    fSkyMapSizeYmin = -2.;
    fSkyMapSizeYmax =  2.;
    
    // sky maps are centred around this point
    fSkyMapCentreNorth = 0.;
    fSkyMapCentreWest = 0.;
    fSkyMapCentreRAJ2000 = 0.;
    fSkyMapCentreDecJ2000 = 0.;
    
    // position relative to which 1D histograms are filled
    fTargetShiftNorth = 0.;                       // [deg]
    fTargetShiftWest = 0.;                        // [deg];
    fTargetShiftRAJ2000 = 0.;
    fTargetShiftDecJ2000 = 0.;
    
    // parameter for energy spectra (in log E)
    fEnergyReconstructionSpectralIndex = 2.5;
    fReconstructionType = GEO;
    fEnergyReconstructionMethod = 1;
    fEffectiveAreaVsEnergyMC = 1;             // default: use effective areas vs reconstructed energy (accurate method)
    fEnergySpectrumBinSize = 0.05;
    fEnergyEffectiveAreaSmoothingIterations = -1;
    fEnergyEffectiveAreaSmoothingThreshold = -1.;
    fDeadTimeCalculationMethod = 0;
    
    // background model
    fTMPL_fBackgroundModel = 0;
    fTMPL_RM_RingRadius = 0.;
    fTMPL_RM_RingWidth = 0.;
    fTMPL_RE_distanceSourceOff = 0.;
    fTMPL_RE_nMinoffsource = 0;
    fTMPL_RE_nMaxoffsource = 0;
    fTMPL_RE_RemoveOffRegionsRandomly = false;
    
    // cut, effective areas and acceptance files
    fTMPL_SourceRadius = 0.;
    fTMPL_maxradius = 0.;
    fTMPL_CutFile = "";
    fTMPL_AcceptanceFile = "";
    fTMPL_EffectiveAreaFile = "";
    
    // length of time intervals in seconds for rate plots and short term histograms
    fTimeIntervall = 4. * 60.;
    
    fWriteEventTree = 2;
    
    // Binned Likelihood
    fLikelihoodAnalysis = false; 

    // if 0, use default 1D radial acceptance
    // if >0, use alternate 2D-dependent acceptance
    f2DAcceptanceMode = 0 ; // USE2DACCEPTANCE
    // use run-wise radial acceptance curve
    fRunWiseRadialAcceptance = false;
    
    // for deadtime fraction storage
    fScalarDeadTimeFrac = 0.0 ;
    
    // set monte carlo zenith angles
    setMCZenith();
    
    // exclusion regions
    // (this is just a placeholder, full list of
    //  exclusion region is kept in
    // VAnaSumRunParameterDataClass )
    fExclusionRegions = new VExclusionRegions();
}

VAnaSumRunParameter::~VAnaSumRunParameter()
{
    if( fExclusionRegions )
    {
        delete fExclusionRegions;
    }
}


int VAnaSumRunParameter::returnWithError( string iL, string iM, string iC )
{
    cout << endl;
    cout << iL << endl;
    cout << iM << endl;
    if( iC.size() > 0 )
    {
        cout << "correct writing: " << endl;
        cout << iC << endl;
    }
    return 0;
}

/*

    read run parameters from an ascii file

*/
int VAnaSumRunParameter::readRunParameter( string i_filename, bool fIgnoreZeroExclusionRegion )
{
    ifstream is;
    is.open( i_filename.c_str(), ifstream::in );
    if( !is )
    {
        string itemp = getDirectory_EVNDISPParameterFiles();
        itemp += "/" + i_filename;
        is.open( itemp.c_str(), ifstream::in );
        if( !is )
        {
            cout << "no file found to read run parameters: " << itemp << endl;
            exit( EXIT_FAILURE );
        }
        i_filename = itemp;
    }
    cout << "Reading anasum parameters from " << i_filename << " :" << endl;
    cout << endl;
    string is_line;
    string temp;
    string temp2;
    
    // loop over all lines in the run parameter file
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
            // print runparameter to stdout
            cout << is_line << endl;
            if( (is_stream>>std::ws).eof() )
            {
                return returnWithError( "VAnaSumRunParameter::readRunParameter: not enough parameters", is_line );
            }
            is_stream >> temp;
            if( (is_stream>>std::ws).eof() )
            {
                return returnWithError( "VAnaSumRunParameter::readRunParameter: not enough parameters", is_line );
            }
            is_stream >> temp2;
            if( temp == "TIMEMASKFILE" )
            {
                fTimeMaskFile = temp2;
                // check if timemask file needs an additional path
                ifstream is_test;
                is_test.open( fTimeMaskFile.c_str(), ifstream::in );
                if( !is_test )
                {
                    string iDIR_temp = i_filename.substr( 0, i_filename.rfind( "/" ) );
                    iDIR_temp += "/" + fTimeMaskFile;
                    is_test.open( iDIR_temp.c_str(), ifstream::in );
                    if( !is_test )
                    {
                        cout << "Error opening time mask file: " << fTimeMaskFile << endl;
                        cout << "exiting..." << endl;
                        exit( EXIT_FAILURE );
                    }
                    else
                    {
                        fTimeMaskFile = iDIR_temp;
                    }
                }
                is_test.close();
            }
            else if( temp == "GAMMAHADRONCUT" )
            {
                fTMPL_CutFile = temp2;
            }
            else if( temp == "RADIALACCEPTANCEFILE" )
            {
                fTMPL_AcceptanceFile = temp2;
            }
            else if( temp == "EFFECTIVEAREAFILE" )
            {
                fTMPL_EffectiveAreaFile = temp2;
            }
            else if( temp == "REFLECTEDREGION" )
            {
                fTMPL_fBackgroundModel = eREFLECTEDREGION;
                fTMPL_RE_distanceSourceOff = atof( temp2.c_str() );
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp2;
                    fTMPL_RE_nMinoffsource = atoi( temp2.c_str() );
                }
                else
                {
                    returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* REFLECTEDREGION dist noff_min noff_max" );
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp2;
                    fTMPL_RE_nMaxoffsource = atoi( temp2.c_str() );
                }
                else
                {
                    returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* REFLECTEDREGION dist noff_min noff_max" );
                }
            }
            else if( temp == "REFLECTEDREGION_OFFREMOVAL" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp2;
                    fTMPL_RE_RemoveOffRegionsRandomly = bool( atoi( temp2.c_str() ) );
                }
            }
            else if( temp == "RINGBACKGROUND" )
            {
                fTMPL_fBackgroundModel = eRINGMODEL;
                fTMPL_RM_RingRadius = atof( temp2.c_str() );
                // important: filling here temporary the
                // area ratio of off-to-on regions
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp2;
                    fTMPL_RM_RingWidth = atof( temp2.c_str() );
                }
                else
                {
                    returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* RINGBACKGROUND ring_radius ring_area" );
                }
            }
            else if( temp == "SKYMAPBINSIZE" )
            {
                fSkyMapBinSize = atof( temp2.c_str() );
            }
            else if( temp == "SKYMAPBINSIZEUC" )
            {
                fSkyMapBinSizeUC = atof( temp2.c_str() );
            }
            else if( temp == "SKYMAPSIZEX" )
            {
                fSkyMapSizeXmin = atof( temp2.c_str() );
                fSkyMapSizeXmin = -1. *TMath::Abs( fSkyMapSizeXmin );
                fSkyMapSizeXmax = TMath::Abs( fSkyMapSizeXmin );
            }
            else if( temp == "SKYMAPSIZEY" )
            {
                fSkyMapSizeYmin = atof( temp2.c_str() );
                fSkyMapSizeYmin = -1. *TMath::Abs( fSkyMapSizeYmin );
                fSkyMapSizeYmax = TMath::Abs( fSkyMapSizeYmin );
            }
            else if( temp == "BRIGHTSTARCATALOGUE" )
            {
                if( fExclusionRegions )
                {
                    fExclusionRegions->setStarCatalogue( temp2 );
                }
            }
            else if( temp == "BRIGHTSTARSETTINGS" )
            {
                double iMinBrightness = atof( temp2.c_str() );
                double iExclusionRadiusDeg = -1.;
                string iStarBand = "";
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iExclusionRadiusDeg;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iStarBand;
                }
                if( fExclusionRegions )
                {
                    fExclusionRegions->addBrightStarSettings( iMinBrightness, iExclusionRadiusDeg, iStarBand );
                }
            }
            else if( temp == "SKYMAPCENTRE_XY" )
            {
                if( checkNumberOfArguments( is_line ) != 4 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* SKYMAPCENTRE_XY x y" );
                }
                
                fSkyMapCentreWest = -1.*atof( temp2.c_str() );
                is_stream >> temp2;
                fSkyMapCentreNorth = -1.*atof( temp2.c_str() );
            }
            else if( temp == "SKYMAPCENTRE_RADECJ2000_DEG" )
            {
                if( checkNumberOfArguments( is_line ) != 4 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line,
                                            "* SKYMAPCENTRE_RADECJ2000_DEG (RA(deg) DEC(deg)" );
                }
                fSkyMapCentreRAJ2000 = atof( temp2.c_str() );
                is_stream >> temp2;
                fSkyMapCentreDecJ2000 = atof( temp2.c_str() );
            }
            else if( temp == "SKYMAPCENTRE_RADECJ2000_HOUR" )
            {
                if( checkNumberOfArguments( is_line ) != 8 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line,
                                            "* SKYMAPCENTRE_RADECJ2000_HOUR RA(Hour Min Sec)  DEC(Deg Min Sec)" );
                }
                
                // ra
                double h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double s = ( double )atof( temp2.c_str() );
                fSkyMapCentreRAJ2000 = VSkyCoordinatesUtilities::getRightAscension_inDegrees_fromHour( h, m, s );
                // dec
                is_stream >> temp2;
                h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                s = ( double )atof( temp2.c_str() );
                fSkyMapCentreDecJ2000 = VSkyCoordinatesUtilities::getDeclination_inDegrees_fromHour( h, m, s );
            }
            else if( temp == "TARGETXYSHIFT" )
            {
                if( checkNumberOfArguments( is_line ) != 4 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* TARGETXYSHIFT x y" );
                }
                
                fTargetShiftWest = -1.*atof( temp2.c_str() );
                is_stream >> temp2;
                fTargetShiftNorth = -1.*atof( temp2.c_str() );
            }
            else if( temp == "TARGETPOSITION_RADECJ2000_DEG" )
            {
                if( checkNumberOfArguments( is_line ) != 4 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* TARGETPOSITION_RADECJ2000_DEG (RA(deg) DEC(deg)" );
                }
                
                fTargetShiftRAJ2000 = atof( temp2.c_str() );
                is_stream >> temp2;
                fTargetShiftDecJ2000 = atof( temp2.c_str() );
            }
            else if( temp == "TARGETPOSITION_RADECJ2000_HOUR" )
            {
                if( checkNumberOfArguments( is_line ) != 8 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* TARGETPOSITION_RADECJ2000_HOUR RA(Hour Min Sec)  DEC(Deg Min Sec)" );
                }
                // ra
                double h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double s = ( double )atof( temp2.c_str() );
                fTargetShiftRAJ2000 = VSkyCoordinatesUtilities::getRightAscension_inDegrees_fromHour( h, m, s );
                // dec
                is_stream >> temp2;
                h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                s = ( double )atof( temp2.c_str() );
                fTargetShiftDecJ2000 = VSkyCoordinatesUtilities::getDeclination_inDegrees_fromHour( h, m, s );
            }
            
            else if( temp == "REGIONTOEXCLUDE" || temp == "REGIONTOEXCLUDE_RADECJ2000_DEG" )
            {
                if( checkNumberOfArguments( is_line ) < 5 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* REGIONTOEXCLUDE (West(deg)  North(deg)  Radius(deg)) (or * REGIONTOEXCLUDE_RADECJ2000_DEG (RA(deg) DEC(deg) Radius(deg)) OR * REGIONTOEXCLUDE (West(deg)  North(deg)  Radius1(deg) Radius2(deg) RotAngle(deg)) (or * REGIONTOEXCLUDE_RADECJ2000_DEG (RA(deg) DEC(deg) Radius1(deg) Radius2(deg) RotAngle(deg)). Check if you want extended or point source!" );
                }
                double iExcludeFromBackground_West = -1.;
                double iExcludeFromBackground_North = -1.;
                double iExcludeFromBackground_RAJ2000 = -99.;
                double iExcludeFromBackground_DecJ2000 = -99.;
                double iExcludeFromBackground_Radius1 = -1.;
                double iExcludeFromBackground_Radius2 = -1.;
                double iExcludeFromBackground_RotAngle = -1.;
                string iExcludeFromBackground_Name = "";
                
                if( temp == "REGIONTOEXCLUDE" )
                {
                    iExcludeFromBackground_West = -1.* ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_North = -1.* ( double )atof( temp2.c_str() );
                }
                else if( temp == "REGIONTOEXCLUDE_RADECJ2000_DEG" )
                {
                    iExcludeFromBackground_RAJ2000 = ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_DecJ2000 = ( double )atof( temp2.c_str() );
                }
                // for circular region
                if( checkNumberOfArguments( is_line ) == 5 || checkNumberOfArguments( is_line ) == 6 )
                {
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius1 = ( double )atof( temp2.c_str() );
                    iExcludeFromBackground_Radius2 = iExcludeFromBackground_Radius1;
                    iExcludeFromBackground_RotAngle = 0.;
                }
                // for ellipsoidal region
                else if( checkNumberOfArguments( is_line ) == 7 || checkNumberOfArguments( is_line ) == 8 )
                {
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius1  = ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius2 = ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_RotAngle = ( double ) atof( temp2.c_str() );
                }
                else
                {
                    return returnWithError( "VAnaSumRunparameter: wrong number of parameters: ", is_line, "Check if you want point or extended source in AnasumRunParameter file!" );
                }
                // read name in
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iExcludeFromBackground_Name;
                }
                bool i_this_is_a_zero_exclusionregion = false;
                if( fIgnoreZeroExclusionRegion )
                {
                    if( TMath::Abs( iExcludeFromBackground_North ) < 1.e-2
                            && TMath::Abs( iExcludeFromBackground_West ) < 1.e-2 )
                    {
                        i_this_is_a_zero_exclusionregion = true;
                    }
                }
                
                if( fExclusionRegions && !i_this_is_a_zero_exclusionregion )
                {
                    fExclusionRegions->addExclusionRegion( iExcludeFromBackground_North, iExcludeFromBackground_West,
                                                           iExcludeFromBackground_RAJ2000, iExcludeFromBackground_DecJ2000,
                                                           iExcludeFromBackground_Radius1, iExcludeFromBackground_Radius2,
                                                           iExcludeFromBackground_RotAngle, iExcludeFromBackground_Name );
                }
            }
            else if( temp == "REGIONTOEXCLUDE_RADECJ2000_HOUR" )
            {
                if( checkNumberOfArguments( is_line ) < 9 )
                {
                    return returnWithError( "VAnaSumRunparameter: not enough parameters: ", is_line, "* REGIONTOEXCLUDE_RADECJ2000_HOUR (RA(Hour Min Sec)  DEC(Deg Min Sec)  Radius(deg)) OR * REGIONTOEXCLUDE_RADECJ2000_HOUR (RA(Hour Min Sec)  DEC(Deg Min Sec)  Radius1(deg) Radius2(deg) RotAngle(deg)). Check if you want extended or point source!" );
                }
                double iExcludeFromBackground_RAJ2000 = -1.;
                double iExcludeFromBackground_DecJ2000 = -1.;
                double iExcludeFromBackground_Radius1 = -1.;
                double iExcludeFromBackground_Radius2 = -1.;
                double iExcludeFromBackground_RotAngle = -1.;
                string iExcludeFromBackground_Name = "";
                
                // ra
                double h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                double s = ( double )atof( temp2.c_str() );
                iExcludeFromBackground_RAJ2000 = VSkyCoordinatesUtilities::getRightAscension_inDegrees_fromHour( h, m, s );
                
                // dec
                is_stream >> temp2;
                h = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                m = ( double )atof( temp2.c_str() );
                is_stream >> temp2;
                s = ( double )atof( temp2.c_str() );
                iExcludeFromBackground_DecJ2000 = VSkyCoordinatesUtilities::getDeclination_inDegrees_fromHour( h, m, s );
                
                // for circular exclusion region
                if( checkNumberOfArguments( is_line ) == 9 || checkNumberOfArguments( is_line ) == 10 )
                {
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius1 = ( double )atof( temp2.c_str() );
                    iExcludeFromBackground_Radius2 = iExcludeFromBackground_Radius1;
                    iExcludeFromBackground_RotAngle = 0.;
                }
                // for ellipsoidal exclusion region
                else if( checkNumberOfArguments( is_line ) == 11 || checkNumberOfArguments( is_line ) == 12 )
                {
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius1 = ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_Radius2 = ( double )atof( temp2.c_str() );
                    is_stream >> temp2;
                    iExcludeFromBackground_RotAngle = ( double )atof( temp2.c_str() );
                }
                else
                {
                    return returnWithError( "VAnaSumRunparameter: wrong number of parameters: ", is_line, "Check if you want point or extended source in AnasumRunParameter file!" );
                }
                // read name in
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iExcludeFromBackground_Name;
                }
                if( fExclusionRegions )
                {
                    fExclusionRegions->addExclusionRegion( -1., -1.,
                                                           iExcludeFromBackground_RAJ2000, iExcludeFromBackground_DecJ2000,
                                                           iExcludeFromBackground_Radius1, iExcludeFromBackground_Radius2,
                                                           iExcludeFromBackground_RotAngle, iExcludeFromBackground_Name );
                }
            }
            else if( temp == "ENERGYBINSIZE" )
            {
                fEnergySpectrumBinSize = atof( temp2.c_str() );
            }
            else if( temp == "ENERGYEFFECTIVEAREAS" )
            {
                if( temp2 == "MC" )
                {
                    fEffectiveAreaVsEnergyMC = 0;
                }
                else if( temp2 == "REC" )
                {
                    fEffectiveAreaVsEnergyMC = 1;
                }
                else
                {
                    cout << "Unknown parameter for ENERGYEFFECTIVEAREAS in parameter file " << i_filename << ": " << temp2 << endl;
                    cout << "use MC or REC (default)" << endl;
                    return 0;
                }
            }
            else if( temp == "ENERGYRECONSTRUCTIONMETHOD" )
            {
                fEnergyReconstructionMethod = ( unsigned int )atoi( temp2.c_str() );
                // print a clear warning if method 0 is selected
                if( fEnergyReconstructionMethod == 0 )
                {
                    cout << endl;
                    cout << "WARNING: energy reconstruction 0 is no longer valid. For any standard analysise, please use method 1 by:" << endl;
                    cout << "  open your anasum run parameter file and replace " << endl;
                    cout << "* ENERGYRECONSTRUCTIONMETHOD 0" << endl;
                    cout << "   by " << endl;
                    cout << "* ENERGYRECONSTRUCTIONMETHOD 1" << endl;
                    cout << "(if you really want to use method 0, you will have to look into the code to find the detour)" << endl;
                    return 0;
                }
                // horrible detour to make sure that users don't use the wrong method
                else if( fEnergyReconstructionMethod == 100 )
                {
                    if( fReconstructionType == NOT_SET )
                    {
                        cout << "Warning: using energy reconstruction method 0" << endl;
                        fEnergyReconstructionMethod = 0;
                        fReconstructionType = ENERGY_ER;
                    }
                    else
                    {
                        cout << "Warning: Analysis type has already been set, this line has no effect." << endl;
                    }
                }
                else if( fEnergyReconstructionMethod > 1 )
                {
                    cout << "Unknown parameter for ENERGYRECONSTRUCTIONMETHOD in parameter file " << i_filename << ": " << temp2 << endl;
                    cout << "allowed values are 0 and 1" << endl;
                    return 0;
                }
            }
            else if( temp == "ReconstructionType" )
            {
                if( temp2 == "GEO" )
                {
                    fReconstructionType = GEO;
                }
                else if( temp2 == "ENERGY_ER" )
                {
                    cout << "Warning: using energy reconstruction method 0" << endl;
                    fReconstructionType = ENERGY_ER;
                    fEnergyReconstructionMethod = 0;
                }
                else if( temp2 == "DEEPLEARNER" )
                {
                    fReconstructionType = DEEPLEARNER;
                }
                else
                {
                    fReconstructionType = NOT_SET;
                    cout << "VAnasumRunParameter::readRunParameter warning: Unknown analysis type ";
                    cout << temp2 << ", will use GammaHadron cut file to decide." << endl;
                }
            }
            else if( temp == "DEADTIMECALCULATIONMETHOD" )
            {
                fDeadTimeCalculationMethod = atoi( temp2.c_str() );
                if( fDeadTimeCalculationMethod != 0 && fDeadTimeCalculationMethod != 1 )
                {
                    cout << "Unknown dead time calculation method (0=scalar method, 1=time difference)" << endl;
                    return 0;
                }
            }
            else if( temp == "RATEINTERVALLLENGTH" )
            {
                fTimeIntervall = atof( temp2.c_str() ) * 60.;
            }
            // expect spectral index positive
            else if( temp == "ENERGYSPECTRALINDEX" )
            {
                fEnergyReconstructionSpectralIndex = fabs( atof( temp2.c_str() ) );
            }
            // effective area smoothing
            else if( temp == "ENERGYEFFAREASMOOTHITER" )
            {
                fEnergyEffectiveAreaSmoothingIterations = atoi( temp2.c_str() );
            }
            else if( temp == "ENERGYEFFAREASMOOTHTHRESH" )
            {
                fEnergyEffectiveAreaSmoothingThreshold = atof( temp2.c_str() );
            }
            
            ////////////////////////////////////////////
            // Option USE2DACCEPTANCE within ANASUM.runparameter
            // * USE2DACCEPTANCE 0
            //     use normal radial acceptance
            // * USE2DACCEPTANCE 1
            //     use simple 2d acceptance model
            else if( temp == "USE2DACCEPTANCE" )
            {
                f2DAcceptanceMode = ( unsigned int )atoi( temp2.c_str() ) ;
            }
            else if( temp == "USERUNWISERADIALACCEPTANCE" )
            {
                fRunWiseRadialAcceptance = true;
            }
            ///////////////////////////////////////////////////////////
            // WRITEEVENTTREE
            // write tree tree with on/off into anasum output file
            else if( temp == "WRITEEVENTTREE" )
            {
                fWriteEventTree = ( unsigned int )atoi( temp2.c_str() );
            }
            /// enable likelihood analysis ///
            else if (temp == "ENABLEBINNEDLIKELIHOOD")
            {
                 unsigned int tmpLikelihood = ( unsigned int )atoi( temp2.c_str() ) ;
                 if( tmpLikelihood == 1)
                 {
                      fLikelihoodAnalysis = true;
                 }
            }
            else
            {
                cout << "Warning: unknown line in parameter file " << i_filename << ": " << endl;
                cout << is_line << endl;
            }
        }
    }
    if( fTMPL_CutFile.size() > 0 )
    {
        fTMPL_SourceRadius = readSourceRadius( fTMPL_CutFile );
        if( fTMPL_SourceRadius <= 0. )
        {
            cout << "error in reading run parameters: ";
            cout << "invalid source radius " << fTMPL_SourceRadius << endl;
            exit( EXIT_FAILURE );
        }
        // note: fTMPL_RM_RingWidth is overwritten to carry to correct variable
        fTMPL_RM_RingWidth = getRingWidth( TMath::Pi() * fTMPL_SourceRadius, fTMPL_RM_RingRadius, fTMPL_RM_RingWidth );
        fTMPL_maxradius = readMaximumDistance( fTMPL_CutFile );
        if( fTMPL_maxradius <  0. )
        {
            cout << "error in reading run parameters: ";
            cout << "invalid maximum distance " << fTMPL_maxradius << endl;
            exit( EXIT_FAILURE );
        }
    }
    else
    {
        fTMPL_SourceRadius = 0.1;
        fTMPL_maxradius = 2.0;
    }
    // prelimary: require same extension in x and y
    if( fabs( fSkyMapSizeXmax - fSkyMapSizeYmax ) > 1.e-3 )
    {
        return returnWithError( "VAnaSumRunParameter::readRunParameter: x and y extension of the sky map should be the same (preliminary)", "" );
    }
    is.close();
    cout << "========================================================" << endl;
    cout << "        end reading run parameters                      " << endl;
    cout << "========================================================" << endl;
    cout << endl;
    
    return 1;
}

/*
 * this is used by the radial acceptance code
 */
int VAnaSumRunParameter::loadSimpleFileList( string i_listfilename )
{
    int i_nline = 0;
    ifstream is;
    is.open( i_listfilename.c_str(), ifstream::in );
    if( !is )
    {
        cout << " VAnaSumRunParameter:::loadSimpleFileList error: file with list of runs not found : " << i_listfilename << endl;
        cout << "exiting..." << endl;
        exit( -1 );
    }
    string is_line;
    string temp;
    VAnaSumRunParameterDataClass i_sT;
    reset( i_sT );
    
    cout << "Reading simple run list from: " << i_listfilename << endl;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            // read run list
            i_sT.fRunOn = atoi( temp.c_str() );
            i_sT.fRunOff = atoi( temp.c_str() );
            // fill the runlist vector
            i_sT.f2DAcceptanceMode = f2DAcceptanceMode ; // USE2DACCEPTANCE
            i_sT.fRunWiseRadialAcceptance = fRunWiseRadialAcceptance;
            fRunList.push_back( i_sT );
            // fill the runlist map
            fMapRunList[i_sT.fRunOn] = fRunList.back();
            ++i_nline;
        }
    }
    return i_nline;
}

/*

   read long run list

*/
int VAnaSumRunParameter::loadLongFileList( string i_listfilename, bool bShortList, bool bTotalAnalysisOnly )
{
    int i_nline = 0;
    ifstream is;
    is.open( i_listfilename.c_str(), ifstream::in );
    if( !is )
    {
        cout << " VAnaSumRunParameter::loadLongFileList error: file with list of runs not found : " << i_listfilename << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    string temp;
    VAnaSumRunParameterDataClass i_sT;
    reset( i_sT );
    
    cout << "Reading long run list from (S" << bShortList << ", TA" << bTotalAnalysisOnly << "): " << i_listfilename << endl;
    
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
            
            int narg = checkNumberOfArguments( is_line );
            // check version number
            if( narg == 3 )
            {
                is_stream >> temp;
                
                if( temp == "VERSION" || temp == "Version" || temp == "version" )
                {
                    is_stream >> temp;
                    fVersion = atoi( temp.c_str() );
                    continue;
                }
            }
            checkNumberOfArguments( -1, narg, i_listfilename, is_line, fVersion, bShortList );
            is_stream >> temp;
            // read run list
            i_sT.fRunOn = atoi( temp.c_str() );
            is_stream >> temp;
            i_sT.fRunOff = atoi( temp.c_str() );
            
            // short list, only read run numbers and target name
            if( bShortList )
            {
                // fill the runlist vector
                i_sT.f2DAcceptanceMode = f2DAcceptanceMode ; // USE2DACCEPTANCE
                i_sT.fRunWiseRadialAcceptance = fRunWiseRadialAcceptance;
                fRunList.push_back( i_sT );
                // fill the runlist map
                fMapRunList[i_sT.fRunOn] = fRunList.back();
                ++i_nline;
                continue;
            }
            
            // offset in min between on and off run (positive if off run after on run)
            // (now read from VEvndispRunParameter)
            if( fVersion < 6 )
            {
                is_stream >> temp;
                i_sT.fPairOffset = atof( temp.c_str() );
            }
            // cut selector (now in cut file, therefore ignored)
            is_stream >> temp;
            ////////////
            // cut file
            // (in >=7, cuts are read from the effective area file
            if( fVersion < 7 )
            {
                is_stream >> temp;
                i_sT.fCutFile = temp;

                // source radius (actually (source radius)^2 )
                // (read theta2 cut from cut file)
                if( !bTotalAnalysisOnly )
                {
                    i_sT.fSourceRadius = readSourceRadius( i_sT.fCutFile );
                }
                else
                {
                    i_sT.fSourceRadius = 0.1;
                }
                if( i_sT.fSourceRadius <= 0. )
                {
                    cout << "VAnaSumRunParameter::loadLongFileList: error in run list: " << endl;
                    cout << is_line << endl;
                    cout << "invalid source radius " << i_sT.fSourceRadius << endl;
                    exit( EXIT_FAILURE );
                }
            }
            // background model
            is_stream >> temp;
            if( temp == "RE" )
            {
                i_sT.fBackgroundModel = eREFLECTEDREGION;
            }
            else if( temp == "RB" )
            {
                i_sT.fBackgroundModel = eRINGMODEL;
            }
            else if( temp == "OO" || temp == "ONOFF" )
            {
                i_sT.fBackgroundModel = eONOFF;
            }
            else if( temp == "TML" )
            {
                i_sT.fBackgroundModel = eTEMPLATE;
            }
            else
            {
                i_sT.fBackgroundModel = atoi( temp.c_str() );
            }
            checkNumberOfArguments( i_sT.fBackgroundModel, narg, i_listfilename, is_line, fVersion, bShortList );
            // maximum distance for events from camera center
            // (read maximum distance from cut file)
            if( fVersion < 2 )
            {
                is_stream >> temp;
                i_sT.fmaxradius = atof( temp.c_str() );
            }
            else if( fVersion < 7 )
            {
                if( !bTotalAnalysisOnly )
                {
                    i_sT.fmaxradius =  readMaximumDistance( i_sT.fCutFile );
                }
                else
                {
                    i_sT.fmaxradius = 2.0;
                }
                if( i_sT.fmaxradius < 0. )
                {
                    cout << "VAnaSumRunParameter::loadLongFileList: error in run list: " << endl;
                    cout << is_line << endl;
                    cout << "invalid maximum distance " << i_sT.fmaxradius << endl;
                    exit( EXIT_FAILURE );
                }
            }
            // file for effective areas
            is_stream >> temp;
            i_sT.fEffectiveAreaFile = temp;
            // cuts are in the effective area files
            if( fVersion >= 7 )
            {
                if( i_sT.fEffectiveAreaFile.find( "IGNOREEFFECTIVEAREA" ) != string::npos )
                {
                    cout << "VAnaSumRunParameter::loadLongFileList warning: ";
                    cout << "ignore effective areas - cannot read cuts from effective area files" << endl;
                }
                else
                {
                    // check if IRF runparameters are consistent with ANASUM.runparameter file
                    checkAnasumParameter( i_sT.fEffectiveAreaFile );
                    
                    i_sT.fCutFile = temp;
                    // source radius (actually (source radius)^2 )
                    // (read theta2 cut from cut file)
                    if( !bTotalAnalysisOnly )
                    {
                        readCutParameter( i_sT.fCutFile, i_sT.fSourceRadius, i_sT.fmaxradius );
                    }
                    else
                    {
                        i_sT.fSourceRadius = 0.1;
                        i_sT.fmaxradius = 2.0;
                    }
                    if( i_sT.fSourceRadius <= 0. )
                    {
                        cout << "VAnaSumRunParameter::loadLongFileList: error in run list: " << endl;
                        cout << is_line << endl;
                        cout << "invalid source radius " << i_sT.fSourceRadius << endl;
                        exit( EXIT_FAILURE );
                    }
                }
            }
            // background model dependend parameters
            //
            //	  if( i_sT.fBackgroundModel == eONOFF )
            //	  {
            // nothing here
            //	  }
            if( i_sT.fBackgroundModel == eRINGMODEL )
            {
                is_stream >> temp;
                i_sT.fRM_RingRadius = atof( temp.c_str() );
                is_stream >> temp;
                i_sT.fRM_RingWidth  = getRingWidth( TMath::Pi() * i_sT.fSourceRadius, i_sT.fRM_RingRadius, atof( temp.c_str() ) );
                is_stream >> temp;
                i_sT.fAcceptanceFile = temp;
            }
            
            else if( i_sT.fBackgroundModel == eREFLECTEDREGION )
            {
                is_stream >> temp;
                i_sT.fRE_distanceSourceOff = atof( temp.c_str() );
                is_stream >> temp;
                i_sT.fRE_nMinoffsource = atoi( temp.c_str() );
                is_stream >> temp;
                i_sT.fRE_nMaxoffsource = atoi( temp.c_str() );
                if( i_sT.fRE_nMaxoffsource > 10 )
                {
                    cout << "VAnaSumRunParameter::loadLongFileList():";
                    cout << "warning, a large number of reflection regions might introduce a gradient into sky maps (consider values <10)" << endl;
                }
                is_stream >> temp;
                i_sT.fAcceptanceFile = temp;
            }
            
            /////////////////
            
            else if( i_sT.fBackgroundModel == eTEMPLATE )
            {
                is_stream >> temp;
                i_sT.fTE_mscw_min = atof( temp.c_str() );
                is_stream >> temp;
                i_sT.fTE_mscw_max = atof( temp.c_str() );
                is_stream >> temp;
                i_sT.fTE_mscl_min = atof( temp.c_str() );
                is_stream >> temp;
                i_sT.fTE_mscl_max = atof( temp.c_str() );
                cout << "DO NOT USE " << endl;
                exit( 0 );
            }
            else if( i_sT.fBackgroundModel == eONOFF )
            {
                // off runs are weighted the same as on runs
                i_sT.fOO_alpha = 1.;
            }
            // fill the runlist vector
            i_sT.f2DAcceptanceMode = f2DAcceptanceMode ; // USE2DACCEPTANCE
            i_sT.fRunWiseRadialAcceptance = fRunWiseRadialAcceptance;
            fRunList.push_back( i_sT );
            // fill the runlist map
            fMapRunList[i_sT.fRunOn] = fRunList.back();
            ++i_nline;
        }
    }
    cout << "\t finished reading " << i_nline << " lines from run list " << endl;
    
    return i_nline;
}


void VAnaSumRunParameter::printStereoParameter( int ion )
{
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        if( fRunList[i].fRunOn == ion )
        {
            printStereoParameter( i );
        }
    }
}


void VAnaSumRunParameter::printStereoParameter( unsigned int i )
{
    if( i < fRunList.size() )
    {
        int ioff = fRunList[i].fRunOff;
        
        cout << "Stereo analysis for run: " << fRunList[i].fRunOn << "\t" << ioff;
        cout << "  (run " << i + 1 << " out of " << fRunList.size() << ")";
        if( fRunList[i].fEventDisplayVersion.size() > 0 )
        {
            cout << ", eventdisplay version " << fRunList[i].fEventDisplayVersion;
        }
        cout << endl;
        cout << "\t Object: " << fRunList[i].fTarget;
        cout << " (background: " << fRunList[i].fOff_Target << ")" << endl;
        cout << "\t Wobble: (N" << fRunList[i].fWobbleNorth << ", W" << fRunList[i].fWobbleWest << ")";
        cout << ", sky maps centred at (ra,dec) (" << fSkyMapCentreRAJ2000 << ", " << fSkyMapCentreDecJ2000 << ")";
        cout << ", target shift: (N" << fRunList[i].fTargetShiftNorth << ", W" << fRunList[i].fTargetShiftWest << ")";
        cout << " (RA/DEC)_J2000 [" << fRunList[i].fTargetShiftRAJ2000 << ", " << fRunList[i].fTargetShiftDecJ2000 << "]" <<  endl;
        cout << "\t user defined exclusion regions:" << endl;
        if( fExclusionRegions )
        {
            fExclusionRegions->printExclusionRegions();
        }
        cout << "\t number of telescopes: " << fRunList[i].fNTel << endl;
        cout << "\t time interval for rate plots: " << fTimeIntervall << " s (" << fTimeIntervall / 60. << " min)" << endl;
        cout << "\t effective areas from " << fRunList[i].fEffectiveAreaFile << endl;
        cout << "\t sky plot binning [deg] " << fSkyMapBinSize << "\t" << fSkyMapBinSizeUC << endl;
        cout << "\t sky plot size [deg]: " << fSkyMapSizeXmin << " < X < " << fSkyMapSizeXmax;
        cout << ", " << fSkyMapSizeYmin << " < Y < " << fSkyMapSizeYmax << endl;
        cout << "\t energy spectra parameters (binsize, log10): " << fEnergySpectrumBinSize;
        if( fEffectiveAreaVsEnergyMC == 0 )
        {
            cout << " (use effective area A_MC)";
        }
        else if( fEffectiveAreaVsEnergyMC == 1 )
        {
            cout << " (use effective area A_REC)";
        }
        cout << ", Method " << fEnergyReconstructionMethod << endl;
        cout << "\t dead time calculation method: ";
        if( fDeadTimeCalculationMethod == 0 )
        {
            cout << "scalar method" << endl;
        }
        else
        {
            cout << "tdiff method" << endl;
        }
        
        cout << "\t background model: ";
        if( fRunList[i].fBackgroundModel == eONOFF )
        {
            cout << "ON/OFF BACKGROUND MODEL" << endl;
            cout << "\t theta2 cut: " << fRunList[i].fSourceRadius << " deg2" << endl;
            cout << "\t maximum distance to camera center: " << fRunList[i].fmaxradius << endl;
        }
        else if( fRunList[i].fBackgroundModel == eRINGMODEL )
        {
            cout << "RING BACKROUND MODEL" << endl;
            cout << "\t theta2 cut: " << fRunList[i].fSourceRadius << " deg2" << endl;
            cout << "\t ring radius: " << fRunList[i].fRM_RingRadius << " deg" << endl;
            cout << "\t ring width: " << fRunList[i].fRM_RingWidth << " deg" << endl;
            cout << "\t area ratio source region to ring: " << 2 * fRunList[i].fRM_RingRadius* fRunList[i].fRM_RingWidth / fRunList[i].fSourceRadius << endl;
            cout << "\t acceptance file: " << fRunList[i].fAcceptanceFile << endl;
            cout << "\t maximum distance to camera center: " << fRunList[i].fmaxradius << " deg" << endl;
        }
        else if( fRunList[i].fBackgroundModel == eREFLECTEDREGION )
        {
            cout << "REFLECTED REGIONS BACKGROUND MODEL";
            if( fTMPL_RE_RemoveOffRegionsRandomly )
            {
                cout << " (excess regions are removed randomly)";
            }
            cout << endl;
            cout << "\t theta2 cut: " << fRunList[i].fSourceRadius << " deg2" << endl;
            cout << "\t distance to source region: " << fRunList[i].fRE_distanceSourceOff << endl;
            cout << "\t minimum number of off source regions: " << fRunList[i].fRE_nMinoffsource << endl;
            cout << "\t maximum number of off source regions: " << fRunList[i].fRE_nMaxoffsource << endl;
            cout << "\t acceptance file (for theta2 plots): " << fRunList[i].fAcceptanceFile << endl;
            cout << "\t maximum distance to camera center: " << fRunList[i].fmaxradius << " deg" << endl;
        }
        if( fRunWiseRadialAcceptance )
        {
            cout << "\t Using run-wise radial acceptance curves" << endl;
        }
    }
    cout << endl;
}


int VAnaSumRunParameter::checkNumberOfArguments( string is )
{
    // get rid of trailing spaces
    while( is.rfind( " " ) > is.size() - 2 )
    {
        is = is.substr( 0, is.size() - 1 );
    }
	// Need to remove newline character from the string
	// since it is counted as an additional parameter 
	is.erase(std::remove(is.begin(), is.end(), '\n'), is.end());
	is.erase(std::remove(is.begin(), is.end(), '\r'), is.end());
    istringstream is_stream( is );
    string itemp;
    int z = 0;
    while( !(is_stream>>std::ws).eof() )
    {
        is_stream >> itemp;
        z++;
    }
    return z;
}


void VAnaSumRunParameter::checkNumberOfArguments( int im, int narg, string i_listfilename, string is_line, int iversion, bool bShortList )
{
    if( bShortList && narg > 3 )
    {
        return;
    }
    
    int n_tot = 0;
    if( im == -1 )
    {
        n_tot = 10;
    }
    else if( im == 0 )
    {
        n_tot = 12;
    }
    else if( im == 1 )
    {
        n_tot = 16;
    }
    else if( im == 2 )
    {
        n_tot = 16;
    }
    else if( im == 3 )
    {
        n_tot = 15;
    }
    else
    {
        cout << "VAnaSumRunParameter::checkNumberOfArguments error: unknown background model" << endl;
        cout << "exiting..." << endl;
        exit( -1 );
    }
    
    // wobble offsets removed with version >=3
    if( iversion > 2 )
    {
        n_tot -= 2;
    }
    if( iversion > 3 )
    {
        n_tot -= 2;
    }
    if( iversion > 4 && ( im == 1 || im == 3 ) )
    {
        n_tot -= 1;
    }
    if( iversion > 5 )
    {
        n_tot -= 1;    // no more RA offset for off runs
    }
    if( iversion > 6 )
    {
        n_tot -= 1;    // no more cut file
    }
    
    if( ( im == -1 && narg < n_tot ) || ( im >= 0 && narg != n_tot ) )
    {
        cout << "error: not enough/too many parameter in " << i_listfilename << ": " << endl;
        cout << is_line << endl;
        cout << "expected " << n_tot << " parameter, found " << narg << " parameters" << endl;
        cout << "(" << im << ", " << narg << ", run list version " << iversion;
        if( bShortList )
        {
            cout << ", shortlist";
        }
        cout << ")" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
}


void VAnaSumRunParameter::reset( VAnaSumRunParameterDataClass it )
{
    it.fRunOn = 0;
    it.fRunOff = 0;
    it.fTarget = "target";
    it.fPairOffset = 0.;
    it.fWobbleNorth = 0.;
    it.fWobbleWest = 0.;
    it.fWobbleNorthMod = 0.;
    it.fWobbleWestMod = 0.;
    it.fOff_WobbleNorth = 0.;
    it.fOff_WobbleWest = 0.;
    it.fNTel = 4;
    it.fMaxTelID = it.fNTel;
    it.fTelToAnalyze.clear();
    it.fBackgroundModel = 0;
    it.fSourceRadius = 0.;
    it.fmaxradius = 0.;
    
    it.fAcceptanceFile = "";
    
    it.fOO_alpha = 0.;
    
    it.fRM_RingRadius = 0.;
    it.fRM_RingWidth = 0.;
    it.fAcceptanceFile = "";
    
    it.fRE_distanceSourceOff = 0.;
    it.fRE_nMinoffsource = 0;
    it.fRE_nMaxoffsource = 0;
    
    it.fAcceptanceFile = "";
    
    it.fTE_mscw_min = 0.;
    it.fTE_mscw_max = 0.;
    it.fTE_mscl_min = 0.;
    it.fTE_mscl_max = 0.;
}

/* 
 * calculate ring width for the ring background model
 *
 * a_on = Area of on region
 * rr   = ring radius
 * rat = area ratio of off-to-on regions
 *
 * inner radius of ring: rr - ring width / 2
 * outer radius of ring: rr + ring width / 2
 *
 */
double VAnaSumRunParameter::getRingWidth( double a_on, double rr, double rat )
{
    if( rr == 0. )
    {
        return 0.;
    }
    
    return rat / 4. / TMath::Pi() / rr * a_on * 2.;
}

double VAnaSumRunParameter::readSourceRadius( string ifile )
{
    VGammaHadronCuts iC;
    iC.setNTel( 1 );  // irrelevant - but suppresses some warnings
    if( !iC.readCuts( ifile, 0 ) )
    {
        return -1;
    };
    
    if( iC.getTheta2Cut_max() < 0. && iC.getDirectionCutSelector() == 2 )
    {
        return iC.getAngularResolutionAbsoluteMaximum();
    }
    return iC.getTheta2Cut_max();
}

/*
 * read some parameters from gamma/hadron cuts (root file)
 *
 */
bool VAnaSumRunParameter::readCutParameter( string ifile, double& iSourceRadius, double& iMaximumDistance )
{
    iSourceRadius = -1.;
    iMaximumDistance = -1.;
    
    string iEffFile = VUtilities::testFileLocation( ifile, "EffectiveAreas", true );
    
    TFile* iF  = new TFile( iEffFile.c_str() );
    if( iF->IsZombie() )
    {
        cout << "VAnaSumRunParameter::readCutParameter error opening file to read direction cuts: " << endl;
        cout << "\t" << ifile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    VGammaHadronCuts* iC = ( VGammaHadronCuts* )iF->Get( "GammaHadronCuts" );
    if( !iC )
    {
        cout << "VAnaSumRunParameter::readCutParameter error reading direction cut from file: " << endl;
        cout << "\t" << ifile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( iC->getTheta2Cut_max() < 0. && iC->getDirectionCutSelector() == 2 )
    {
        iSourceRadius = iC->getAngularResolutionAbsoluteMaximum();
    }
    else
    {
        iSourceRadius = iC->getTheta2Cut_max();
    }
    iMaximumDistance = iC->fCut_CameraFiducialSize_max;
    
    if( iF )
    {
        iF->Close();
    }
    
    return true;
}

/*
 *  get gamma/hadron parameters from a file
 *
 *  Note: file stays open
 */
VGammaHadronCuts* VAnaSumRunParameter::getGammaHadronCuts( string ifile )
{
    string iEffFile = VUtilities::testFileLocation( ifile, "EffectiveAreas", true );
    
    TFile* iF  = new TFile( iEffFile.c_str() );
    if( iF->IsZombie() )
    {
        cout << "VAnaSumRunParameter::getGammaHadronCuts: error opening file to read direction cuts: " << endl;
        cout << "\t" << ifile << endl;
        return 0;
    }
    VGammaHadronCuts* iC = ( VGammaHadronCuts* )iF->Get( "GammaHadronCuts" );
    if( !iC )
    {
        cout << "VAnaSumRunParameter::getGammaHadronCuts error reading direction cut from file: " << endl;
        cout << "\t" << ifile << endl;
        return 0;
    }
    return iC;
}

/*
 *  Check if requested anasum runparameter make sense, i.e.,
 *  - energy reconstruction method should be the same as in effective area file
 *  - spectral index should be in the range of simulated index range
 *
 */
bool VAnaSumRunParameter::checkAnasumParameter( string ifile )
{
    string iEffFile = VUtilities::testFileLocation( ifile, "EffectiveAreas", true );
    TFile* iF  = new TFile( iEffFile.c_str() );
    if( iF->IsZombie() )
    {
        cout << "VAnaSumRunParameter::checkAnasumParameter error opening file to read IRF parameters: " << endl;
        cout << "\t" << ifile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    VInstrumentResponseFunctionRunParameter* iIRF = ( VInstrumentResponseFunctionRunParameter* )iF->Get( "makeEffectiveArea_runparameter" );
    if( !iIRF )
    {
        cout << "VAnaSumRunParameter::checkAnasumParameter error reading IRF parameter from file: " << endl;
        cout << "\t" << ifile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else
    {
        // check energy reconstruction method
        if( iIRF->fEnergyReconstructionMethod != fEnergyReconstructionMethod )
        {
            cout << "VAnaSumRunParameter::checkAnasumParameter error in energy reconstruction method specified in runparameter file. " << endl;
            cout << "\t Effective area file (" << ifile << ") uses energy reconstruction method " << iIRF->fEnergyReconstructionMethod << endl;
            cout << "\t but energy reconstruction method " << fEnergyReconstructionMethod << " is requested in the anasum runparameter file. " << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        // check spectral index range
        double iIndexMin = iIRF->fSpectralIndexMin;
        double iIndexMax = iIRF->fSpectralIndexMin + iIRF->fNSpectralIndex * iIRF->fSpectralIndexStep;
        if( fEnergyReconstructionSpectralIndex < iIndexMin || fEnergyReconstructionSpectralIndex > iIndexMax )
        {
            cout << "VAnaSumRunParameter::checkAnasumParameter warning: spectral index out of range. " << endl;
            cout << "\t Requested spectral index (" << fEnergyReconstructionSpectralIndex << ") is outsided the range simulated ["
                 << iIndexMin << "-" << iIndexMax << "]." << endl << endl;
        }
    }
    
    if( iF )
    {
        iF->Close();
    }
    
    return true;
}



double VAnaSumRunParameter::readMaximumDistance( string ifile )
{
    VGammaHadronCuts iC;
    iC.setNTel( 1 );  // irrelevant - but suppressed printing of warnings to screen
    if( !iC.readCuts( ifile, 0 ) )
    {
        return -1;
    }
    
    return iC.fCut_CameraFiducialSize_max;
}


unsigned int VAnaSumRunParameter::getMaxNumberofTelescopes()
{
    unsigned int iMax = 0;
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        if( fRunList[i].fNTel > iMax )
        {
            iMax = fRunList[i].fNTel;
        }
    }
    return iMax;
}


bool VAnaSumRunParameter::setTargetShifts( unsigned int i, double west, double north, double ra, double dec )
{
    if( i < fRunList.size() )
    {
        if( fMapRunList.find( fRunList[i].fRunOn ) != fMapRunList.end() )
        {
            fMapRunList[fRunList[i].fRunOn].fTargetShiftWest = west;
            fMapRunList[fRunList[i].fRunOn].fTargetShiftNorth = north;
            fMapRunList[fRunList[i].fRunOn].fTargetShiftRAJ2000 = ra;
            fMapRunList[fRunList[i].fRunOn].fTargetShiftDecJ2000 = dec;
        }
        return true;
    }
    return false;
}

bool VAnaSumRunParameter::setSkyMapCentreJ2000( unsigned int i, double ra, double dec )
{
    if( i < fRunList.size() )
    {
        fRunList[i].fSkyMapCentreRAJ2000  = ra;
        fRunList[i].fSkyMapCentreDecJ2000 = dec;
        if( fMapRunList.find( fRunList[i].fRunOn ) != fMapRunList.end() )
        {
            fMapRunList[fRunList[i].fRunOn].fSkyMapCentreRAJ2000  = ra;
            fMapRunList[fRunList[i].fRunOn].fSkyMapCentreDecJ2000 = dec;
        }
        return true;
    }
    return false;
}


bool VAnaSumRunParameter::setTargetRADecJ2000( unsigned int i, double ra, double dec, string iTargetName )
{
    if( i < fRunList.size() )
    {
        fRunList[i].fTargetRAJ2000  = ra;
        fRunList[i].fTargetDecJ2000 = dec;
        fRunList[i].fTarget = iTargetName;
        if( fMapRunList.find( fRunList[i].fRunOn ) != fMapRunList.end() )
        {
            fMapRunList[fRunList[i].fRunOn].fTargetRAJ2000  = ra;
            fMapRunList[fRunList[i].fRunOn].fTargetDecJ2000 = dec;
            fMapRunList[fRunList[i].fRunOn].fTarget = iTargetName;
        }
        // set centre of stereo maps (if this parameter is not set in the file runparameter.dat)
        if( TMath::Abs( fSkyMapCentreNorth ) < 1.e-8 && TMath::Abs( fSkyMapCentreWest ) < 1.e-8
                && TMath::Abs( fSkyMapCentreRAJ2000 ) < 1.e-8 && TMath::Abs( fSkyMapCentreDecJ2000 ) < 1.e-8 )
        {
            fRunList[i].fSkyMapCentreNorth    = 0.;
            fRunList[i].fSkyMapCentreWest     = 0.;
            fRunList[i].fSkyMapCentreRAJ2000  = ra;
            fRunList[i].fSkyMapCentreDecJ2000 = dec;
            if( fMapRunList.find( fRunList[i].fRunOn ) != fMapRunList.end() )
            {
                fMapRunList[fRunList[i].fRunOn].fSkyMapCentreNorth    = 0.;
                fMapRunList[fRunList[i].fRunOn].fSkyMapCentreWest     = 0.;
                fMapRunList[fRunList[i].fRunOn].fSkyMapCentreRAJ2000  = ra;
                fMapRunList[fRunList[i].fRunOn].fSkyMapCentreDecJ2000 = dec;
            }
        }
        return true;
    }
    return false;
}


bool VAnaSumRunParameter::setTargetRADec_currentEpoch( unsigned int i, double ra, double dec )
{
    if( i < fRunList.size() )
    {
        fRunList[i].fTargetRA = ra;
        fRunList[i].fTargetDec = dec;
        if( fMapRunList.find( fRunList[i].fRunOn ) != fMapRunList.end() )
        {
            fMapRunList[fRunList[i].fRunOn].fTargetRA = ra;
            fMapRunList[fRunList[i].fRunOn].fTargetDec = dec;
        }
        return true;
    }
    return false;
}


void VAnaSumRunParameter::getEventdisplayRunParameter( string fDatadir )
{
    cout << "\t reading run parameter from mscw file (may take a second...)" << endl;
    char i_temp[200];
    int i_run;
    string i_treename = "data";
    string fPrefix = "";
    string fSuffix = ".mscw.root";
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        i_run = fRunList[i].fRunOn;
        
        sprintf( i_temp, "%s/%s%d%s", fDatadir.c_str(), fPrefix.c_str(), i_run, fSuffix.c_str() );
        TFile* i_f = new TFile( i_temp );
        if( i_f->IsZombie() )
        {
            cout << "VAnaSumRunParameter::getEventdisplayRunParameter fatal error: file not found, " << i_temp << endl;
            exit( -1 );
        }
        VEvndispRunParameter* iParV2 = ( VEvndispRunParameter* )i_f->Get( "runparameterV2" );
        if( iParV2 )
        {
            fRunList[i].fEventDisplayVersion = iParV2->getEVNDISP_VERSION();
            fRunList[i].fTarget              = iParV2->fTargetName;
            fRunList[i].fTargetRAJ2000       = iParV2->fTargetRA;
            fRunList[i].fTargetDecJ2000      = iParV2->fTargetDec;
            fRunList[i].fWobbleNorth         = iParV2->fWobbleNorth;
            fRunList[i].fWobbleWest          = -1.*iParV2->fWobbleEast;
            fRunList[i].fWobbleNorthMod      = iParV2->fWobbleNorth;
            fRunList[i].fWobbleWestMod       = -1.*iParV2->fWobbleEast;
            fRunList[i].fNTel                = ( int )iParV2->fTelToAnalyze.size();
            fRunList[i].fTelToAnalyze        = iParV2->fTelToAnalyze;
        }
        i_f->Close();
        // get maximum telescope ID
        fRunList[i].fMaxTelID = 0;
        for( unsigned int t = 0; t < fRunList[i].fTelToAnalyze.size(); t++ )
        {
            if( fRunList[i].fTelToAnalyze[t] > fRunList[i].fMaxTelID )
            {
                fRunList[i].fMaxTelID = fRunList[i].fTelToAnalyze[t];
            }
        }
        // go from T1 = 0 to T1 = 1
        fRunList[i].fMaxTelID += 1;
    }
    //////////////////////////////////////////////////////////////////////////////
    // background (off) runs
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        i_run = fRunList[i].fRunOff;
        
        sprintf( i_temp, "%s/%s%d%s", fDatadir.c_str(), fPrefix.c_str(), i_run, fSuffix.c_str() );
        TFile* i_f = new TFile( i_temp );
        if( i_f->IsZombie() )
        {
            cout << "VAnaSumRunParameter::getEventdisplayRunParameter fatal error: off file not found, " << i_temp << endl;
            exit( EXIT_FAILURE );
        }
        VEvndispRunParameter* iParV2 = ( VEvndispRunParameter* )i_f->Get( "runparameterV2" );
        if( iParV2 )
        {
            fRunList[i].fPairOffset         = iParV2->fTargetRAOffset * 24. * 60. / 360.;
            fRunList[i].fOff_Target         = iParV2->fTargetName;
            fRunList[i].fOff_TargetRAJ2000  = iParV2->fTargetRA;
            fRunList[i].fOff_TargetDecJ2000 = iParV2->fTargetDec;
            fRunList[i].fOff_WobbleNorth    = iParV2->fWobbleNorth;
            fRunList[i].fOff_WobbleWest     = -1.*iParV2->fWobbleEast;
            if( TMath::Abs( fRunList[i].fPairOffset ) > 1.e-3 )
            {
                cout << "\t on/off pair offset in RA: " << fRunList[i].fPairOffset << endl;
                cout << "\t Warning: pair offset for on/off analysis not taken into account";
                cout << "- this part of the analysis needs development work" << endl;
            }
        }
        i_f->Close();
    }
}


void VAnaSumRunParameter::getWobbleOffsets( string fDatadir )
{
    cout << "\t read wobble offsets from root files (may take a second...)" << endl;
    char i_temp[200];
    int i_run;
    string i_treename = "data";
    string fPrefix = "";
    string fSuffix = ".mscw.root";
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        i_run = fRunList[i].fRunOn;
        
        sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), i_run, fSuffix.c_str() );
        TFile* i_f = new TFile( i_temp );
        if( i_f->IsZombie() )
        {
            cout << "VAnaSumRunParameter::getWobbleOffset fatal error: file not found, " << i_temp << endl;
            exit( -1 );
        }
        
        TTree* i_tree = ( TTree* )i_f->Get( i_treename.c_str() );
        if( !i_tree )
        {
            cout << "VAnaSumRunParameter::getWobbleOffset tree not found " << i_treename << endl;
        }
        i_tree->GetEntry( 1 );
        
        if( i_tree->GetBranch( "WobbleN" ) && i_tree->GetBranch( "WobbleE" ) )
        {
            fRunList[i].fWobbleNorth = i_tree->GetBranch( "WobbleN" )->GetLeaf( "WobbleN" )->GetValue();
            fRunList[i].fWobbleWest = -1.*i_tree->GetBranch( "WobbleE" )->GetLeaf( "WobbleE" )->GetValue();
            fRunList[i].fWobbleNorthMod = i_tree->GetBranch( "WobbleN" )->GetLeaf( "WobbleN" )->GetValue();
            fRunList[i].fWobbleWestMod = -1.*i_tree->GetBranch( "WobbleE" )->GetLeaf( "WobbleE" )->GetValue();
        }
        else
        {
            cout << endl;
            cout << "VAnaSumRunParameter::getWobbleOffset error: cannot determine wobble offset for run " << fRunList[i].fRunOn << " from mscw_energy output file" << endl;
            cout << "\t old file version?" << endl;
            exit( 0 );
        }
        i_f->Close();
    }
}


/*!
    observe that these are hardwired values according to the VERITAS simulation sets
*/
void VAnaSumRunParameter::setMCZenith()
{
    fMCZe.clear();
    
    fMCZe.push_back( 0.0 );
    fMCZe.push_back( 20.0 );
    fMCZe.push_back( 30.0 );
    fMCZe.push_back( 35.0 );
    fMCZe.push_back( 40.0 );
    fMCZe.push_back( 45.0 );
    fMCZe.push_back( 50.0 );
    fMCZe.push_back( 55.0 );
    fMCZe.push_back( 60.0 );
    fMCZe.push_back( 65.0 );
}

/*
    write exclusion regions into root file (into the current directory)
*/
bool VAnaSumRunParameter::writeListOfExcludedSkyRegions( int ionRun )
{
    if( ionRun < 0 && fExclusionRegions )
    {
        return fExclusionRegions->writeExclusionRegionTree();
    }
    else
    {
        for( unsigned int i = 0; i < fRunList.size(); i++ )
        {
            if( fRunList[i].fRunOn == ionRun
                    && fRunList[i].fListOfExclusionRegions )
            {
                fRunList[i].fListOfExclusionRegions->writeExclusionRegionTree( ionRun );
            }
        }
    }
    
    return false;
}

/*

    read list of exclusion regions from a anasum file

*/
bool VAnaSumRunParameter::getListOfExcludedSkyRegions( TFile* f, int iRunOn )
{
    if( fExclusionRegions )
    {
        return fExclusionRegions->readExclusionRegionTree( f, iRunOn );
    }
    
    return false;
}

/*
 * return largest star exclusion radius
 *
 */
double VAnaSumRunParameter::getLargestStarExlusionRadius()
{
    if( fExclusionRegions )
    {
        return fExclusionRegions->getLargestStarExlusionRadius();
    }
    return 0.;
}

vector< VListOfExclusionRegions* > VAnaSumRunParameter::getExclusionRegions( unsigned int iRunCounter )
{
    if( iRunCounter < fRunList.size() && fRunList[iRunCounter].fListOfExclusionRegions )
    {
        return fRunList[iRunCounter].fListOfExclusionRegions->getListOfExclusionRegions();
    }
    vector< VListOfExclusionRegions* > a;
    return a;
}

/*
    initialize exclusion regions

*/
bool VAnaSumRunParameter::initializeExclusionRegions( unsigned int iRunCounter, VStarCatalogue* iCatalogue,
        double iSkyMapCentre_ra_deg, double iSkyMapCentre_dec_deg,
        double iCameraCentre_ra_deg, double iCameraCentre_dec_deg )
{
    if( iRunCounter >= fRunList.size() )
    {
        return false;
    }
    if( !fRunList[iRunCounter].fListOfExclusionRegions )
    {
        fRunList[iRunCounter].fListOfExclusionRegions = new VExclusionRegions();
    }
    // copy stuff from general list
    if( fExclusionRegions )
    {
        fRunList[iRunCounter].fListOfExclusionRegions->setStarCatalogue( fExclusionRegions->fStarCatalogue );
        fRunList[iRunCounter].fListOfExclusionRegions->addBrightStarSettings( fExclusionRegions->getBrightStarSettings() );
        for( unsigned int i = 0; i < fExclusionRegions->fExclusionRegions.size(); i++ )
        {
            if( fExclusionRegions->fExclusionRegions[i] )
            {
                fRunList[iRunCounter].fListOfExclusionRegions->addExclusionRegion( fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_North,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_West,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_RAJ2000,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_DecJ2000,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_Radius1,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_Radius2,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_RotAngle,
                        fExclusionRegions->fExclusionRegions[i]->fExcludeFromBackground_StarName );
            }
        }
        // for consistency: make sure that that general list if correct (not used)
        fExclusionRegions->initializeExclusionRegions( iCatalogue,
                iSkyMapCentre_ra_deg, iSkyMapCentre_dec_deg,
                iCameraCentre_ra_deg, iCameraCentre_dec_deg );
    }
    // initialize exclusion regions
    fRunList[iRunCounter].fListOfExclusionRegions->initializeExclusionRegions( iCatalogue,
            iSkyMapCentre_ra_deg, iSkyMapCentre_dec_deg,
            iCameraCentre_ra_deg, iCameraCentre_dec_deg );
            
            
    // print list to screen
    fRunList[iRunCounter].fListOfExclusionRegions->printExclusionRegions();
    
    return true;
}
