/*! \class VEvndispReconstructionParameter
    \brief reading and storage class for eventdisplay reconstruction parameters

    note: due to historical reasons, the term 'method' is not used in the right way
    in this routine and in VArrayAnalyzer.cpp

    one value per telescope

    there should be one instance per reconstruction method of this class

*/

#include "VEvndispReconstructionParameter.h"

ClassImp( VEvndispReconstructionParameter )
ClassImp( VEvndispReconstructionParameterData )
ClassImp( VEvndispReconstructionCut )

VEvndispReconstructionParameter::VEvndispReconstructionParameter()
{
    fDebug = false;
    fNTel_type = 0;
    fRunPara = 0;
}


VEvndispReconstructionParameter::VEvndispReconstructionParameter( vector< ULong64_t > i_telType, VEvndispRunParameter* iRunPara )
{
    reset();
    
    fRunPara = iRunPara;
    
    fTel_type_perTelescope = i_telType;
    
    // get set with telescope types
    for( unsigned int i = 0; i < i_telType.size(); i++ )
    {
        fTel_type.insert( fTel_type_perTelescope[i] );
    }
    fNTel_type = fTel_type.size();
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        fTel_typeUniqueVector.push_back( *fTel_type_iter );
    }
    
}

void VEvndispReconstructionParameter::reset()
{
    fDebug = false;
    fNTel_type = 0;
    fRunPara = 0;
}


/*

     apply array analysis cuts for this set of image parameters

*/
bool VEvndispReconstructionParameter::applyArrayAnalysisCuts( unsigned int iMeth, unsigned int iTel, unsigned int iTelType,
        VImageParameter* iImageParameter, unsigned short int iLocalTriggerType,
        VStarCatalogue* iStarCatalogue )
{
    // sanity checks
    if( iMeth >= fReconstructionParameterData.size() || !fReconstructionParameterData[iMeth] )
    {
        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCuts error: invalid method number " << iMeth << "\t" << fReconstructionParameterData.size() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( iTelType >= fReconstructionParameterData[iMeth]->fNTel_type )
    {
        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCuts error: invalid telescope type " << iTelType << "\t" << fNTel_type << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( !iImageParameter )
    {
        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCuts error: no image parameters given" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    if( fDebug )
    {
        cout << "APPLY ARRAY ANALYSIS CUTS FOR METHOD " << iMeth << " AND TELESCOPE " << iTel + 1 << ", TYPE " << iTelType << endl;
    }
    
    // return value
    bool iArrayCut = true;
    
    // eventstatus
    if( iImageParameter->eventStatus > 0 )
    {
        iArrayCut = false;
        if( fDebug )
        {
            cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut: event status > 0: " << iImageParameter->eventStatus << endl;
        }
    }
    
    ////////////////////////////////////////////
    // L2 trigger type (mainly for CTA prod2)
    if( !fReconstructionParameterData[iMeth]->testL2TriggerType( iTel, iTelType, iLocalTriggerType ) )
    {
        iArrayCut = false;
    }
    ////////////////////////////////////////////
    // general reconstruction quality cuts
    
    // user set: remove image
    iArrayCut = iArrayCut && fReconstructionParameterData[iMeth]->testUserImage( iTelType );
    
    // image size
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "SIZE", iImageParameter->size, iImageParameter->ntubes );
    // number of ntubes (<=!!!)
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "NTUBES", ( int )iImageParameter->ntubes, iImageParameter->ntubes );
    // number of saturated channels
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "NLOWGAIN", ( int )iImageParameter->nlowgain, iImageParameter->ntubes );
    // image width
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "WIDTH", iImageParameter->width, iImageParameter->ntubes );
    // image distance to camera centre
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "DISTANCE", iImageParameter->dist, iImageParameter->ntubes );
    // loss cut
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "LOSS1", iImageParameter->loss, iImageParameter->ntubes );
    // loss cut
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "LOSS2", iImageParameter->loss, iImageParameter->ntubes );
    // width over length cut
    if( iImageParameter->length )
    {
        iArrayCut = iArrayCut
                    && fReconstructionParameterData[iMeth]->test( iTelType, "WIDTHLENGTH", iImageParameter->width / iImageParameter->length, iImageParameter->ntubes );
    }
    // fui cut
    iArrayCut = iArrayCut
                && fReconstructionParameterData[iMeth]->test( iTelType, "FUI", iImageParameter->fui, iImageParameter->ntubes );
    // MC only: cut on MC energy (use with care!)
    if( fRunPara->isMC() )
    {
        iArrayCut = iArrayCut
                    && fReconstructionParameterData[iMeth]->test( iTelType, "MCENERGY_LINTEV", iImageParameter->MCenergy, iImageParameter->ntubes );
    }
    
    ////////////////////////////////////////////
    // cut on successfull LL reconstruction on the edge of the FOV
    if( fRunPara )
    {
        if( iTel < fRunPara->fLogLikelihoodLoss_min.size() && iTel < fRunPara->fLogLikelihoodLoss_max.size() )
        {
            if( iImageParameter->loss > fRunPara->fLogLikelihoodLoss_min[iTel] && iImageParameter->Fitstat < 1
                    && iImageParameter->loss > fRunPara->fLogLikelihoodLoss_max[iTel] )
            {
                iArrayCut = false;
                if( fDebug )
                {
                    cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut: fit stat cut: ";
                    cout << iImageParameter->Fitstat << endl;
                }
            }
            // check number of events at the edge of the FOV
            if( iTel < fRunPara->fLogLikelihoodLoss_min.size()
                    && iTel < fRunPara->fLogLikelihoodLoss_max.size()
                    && iTel < fRunPara->fLogLikelihood_Ntubes_min.size() )
            {
                if( iImageParameter->loss   >  fRunPara->fLogLikelihoodLoss_min[iTel]
                        && iImageParameter->loss  <  fRunPara->fLogLikelihoodLoss_max[iTel]
                        && iImageParameter->ntubes <= fRunPara->fLogLikelihood_Ntubes_min[iTel] )
                {
                    iArrayCut = false;
                    if( fDebug )
                    {
                        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut: LL ntubes cut: ";
                        cout << iImageParameter->ntubes;
                        cout << " (" <<  fRunPara->fLogLikelihood_Ntubes_min[iTel] << ")" << endl;
                    }
                }
            }
        }
    }
    
    
    ////////////////////////////////////////////
    // remove image which is too close to a bright star
    // (use list of image and border pixels)
    // __this cut is disabled__
    if( iStarCatalogue && fRunPara && iImageParameter->ntubes < fRunPara->fMinStarNTubes )
    {
        for( unsigned int i = 0; i < iImageParameter->fImageBorderPixelPosition_x.size(); i++ )
        {
            if( i < iImageParameter->fImageBorderPixelPosition_y.size() )
            {
                if( iStarCatalogue->getDistanceToClosestStar( iImageParameter->fImageBorderPixelPosition_x[i],
                        iImageParameter->fImageBorderPixelPosition_y[i] ) < fRunPara->fMinStarPixelDistance_deg )
                {
                    iArrayCut = false;
                    if( fDebug )
                    {
                        cout << "Telescope " << iTel + 1 << endl;
                        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut: bright star cut: ";
                        cout << iStarCatalogue->getDistanceToClosestStar( iImageParameter->cen_x, iImageParameter->cen_y );
                        cout << " (" << fRunPara->fMinStarPixelDistance_deg << " deg )" << endl;
                    }
                    if( !iArrayCut )
                    {
                        break;
                    }
                }
            }
        }
    }
    
    if( fDebug )
    {
        cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut: cut: " << iArrayCut << endl;
    }
    
    return iArrayCut;
}

/*

    add a new cut method

*/
void VEvndispReconstructionParameter::addNewMethod( unsigned int iRecordID, unsigned int iMethodID )
{
    // new reconstruction parameter data element and fill with default values
    fReconstructionParameterData.push_back( new VEvndispReconstructionParameterData( fNTel_type, fTel_type ) );
    fReconstructionParameterData.back()->setDebug( fDebug );
    fReconstructionParameterData.back()->setMethodID( iMethodID );
    
}


void VEvndispReconstructionParameter::print_arrayAnalysisCuts()
{
    //////////////////////////////////////////////////////////////////////
    // print FADC integration parameters
    
    for( unsigned int i = 0; i < fRunPara->fTelToAnalyze.size(); i++ )
    {
        if( fRunPara->fTelToAnalyze[i] < fRunPara->fDoublePass.size() && fRunPara->fDoublePass[fRunPara->fTelToAnalyze[i]] )
        {
            cout << "double pass cleaning for Telescope " << fRunPara->fTelToAnalyze[i] + 1;
            cout << " (uncertainties: ";
            if( fRunPara->fTelToAnalyze[i] < fRunPara->fDoublePassErrorWeighting2005.size() && fRunPara->fDoublePassErrorWeighting2005[fRunPara->fTelToAnalyze[i]] )
            {
                cout << "2005) ";
            }
            else
            {
                cout << "2013) ";
            }
            cout << " (low gain window shift: " << fRunPara->fTraceWindowShift[fRunPara->fTelToAnalyze[i]] << ")" << endl;
        }
    }
    
    //////////////////////////////////////////////////////////////////////
    // print array analysis cuts
    cout << endl;
    cout << "------------------------------" << endl;
    cout << "----- Array Analysis Cuts ----" << endl;
    if( fReconstructionParameterData.size() > 1 )
    {
        cout << "------(" << fReconstructionParameterData.size() << " methods)------" << endl;
    }
    else
    {
        cout << "------(" << fReconstructionParameterData.size() << " method)------" << endl;
    }
    cout << endl;
    for( unsigned m = 0; m < fReconstructionParameterData.size(); m++ )
    {
        cout << "\t set number: " << m << ", array reconstruction method: " << fReconstructionParameterData[m]->fMethodID;
        if( fReconstructionParameterData[m]->fUseEventdisplayPointing )
        {
            cout << ", use eventdisplay pointing (no pointing correction from DB or pointing monitor)";
        }
        cout << endl;
        cout << "\t\t minimum angle between image axes [deg]: " << fReconstructionParameterData[m]->fAxesAngles_min << endl;
        // loop over all telescope types
        for( unsigned int t = 0; t < fReconstructionParameterData[m]->fTelescopeType.size(); t++ )
        {
            // is this imaged used?
            cout << "\t\t TelType " << fReconstructionParameterData[m]->fTelescopeType[t];
            if( t < fReconstructionParameterData[m]->fLocalUseImage.size()
                    &&      fReconstructionParameterData[m]->fLocalUseImage[t] )
            {
                cout << " image(s) used in reconstruction";
            }
            else
            {
                cout << " image(s) not used in reconstruction";
            }
            cout << endl;
            // L2 trigger type
            if( t < fReconstructionParameterData[m]->fL2TriggerType.size()
                    && fReconstructionParameterData[m]->fL2TriggerType[t] < 9999 )
            {
                cout << "\t\t TelType " << fReconstructionParameterData[m]->fTelescopeType[t];
                cout << ": L2 trigger type ";
                cout << fReconstructionParameterData[m]->fL2TriggerType[t] << endl;
            }
            // loop over all reconstruction / quality cuts and print them
            map< string, VEvndispReconstructionCut* >::iterator i_iterator_intMap;
            for( i_iterator_intMap = fReconstructionParameterData[m]->fTelescopeTypeCut[t].begin();
                    i_iterator_intMap != fReconstructionParameterData[m]->fTelescopeTypeCut[t].end();
                    ++i_iterator_intMap )
            {
                i_iterator_intMap->second->print( fReconstructionParameterData[m]->fTelescopeType[t],
                                                  i_iterator_intMap->first );
            }
        }
        if( fReconstructionParameterData[m]->fDISP_TMVAFileNameVector.size() > 0 )
        {
            for( unsigned int ze = 0; ze < fReconstructionParameterData[m]->fDISP_TMVAFileNameVector.size(); ze++ )
            {
                if( fReconstructionParameterData[m]->fDISP_TMVAFileNameVector[ze].size() > 0
                        && fReconstructionParameterData[m]->fDISP_TMVAFileNameVector[ze].find( "USE_BDT_METHOD" ) != string::npos )
                {
                    cout << "\t\t TMVA (BDT) file used from method ";
                    cout << fReconstructionParameterData[m]->fDISP_TMVAFileNameVector[ze] << endl;
                }
                else
                {
                    cout << "\t\t TMVA (BDT) file for zenith angle ";
                    if( ze < fReconstructionParameterData[m]->fDISP_TMVAZenithBin.size() )
                    {
                        cout << fReconstructionParameterData[m]->fDISP_TMVAZenithBin[ze] << " deg: ";
                    }
                    cout << fReconstructionParameterData[m]->fDISP_TMVAFileNameVector[ze] << endl;
                }
                if( ze < fReconstructionParameterData[m]->fDISP_MAXTelescopes.size()
                        && ze < fReconstructionParameterData[m]->fDISP_MAXMethodID.size() )
                {
                    cout << "\t\t (use DISP BDT method for up to ";
                    cout << fReconstructionParameterData[m]->fDISP_MAXTelescopes[ze];
                    cout << " images; fallback method is reconstruction method ";
                    cout << fReconstructionParameterData[m]->fDISP_MAXMethodID[ze];
                    cout << ")" << endl;
                }
            }
        }
    }
    
    if( fRunPara && fRunPara->fStarCatalogueName.size() > 0 )
    {
        cout << endl;
        cout << "reading star catalogue from: " << fRunPara->fStarCatalogueName << endl;
        cout << "\t minimum brightness (B): " << fRunPara->fMinStarBrightness_B;
        if( fRunPara->fMinStarPixelDistance_deg < 1.e20 )
        {
            cout << " (max pixel distance " << fRunPara->fMinStarPixelDistance_deg << " deg)";
        }
        if( fRunPara->fMinStarNTubes < 100000 )
        {
            cout << ", (max number of pixels: " << fRunPara->fMinStarNTubes << ")";
        }
        cout << endl;
    }
    cout << "------------------------------" << endl;
    cout << endl;
}

/*

    read all key words for FADC analysis from runparameter file

*/
bool VEvndispReconstructionParameter::readKeyWord_FADCANALYSIS( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
    {
        if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
        {
            // trace integration
            if( i < fRunPara->fTraceIntegrationMethod.size() )
            {
                fRunPara->fTraceIntegrationMethod[i] = ( unsigned int )atoi( iTemp[1].c_str() );
            }
            // digital filter parameters
            if( i < fRunPara->fDF_DigitalFilter.size() && iTemp.size() > 2 )
            {
                fRunPara->fDF_DigitalFilter[i] = ( unsigned int )atoi( iTemp[2].c_str() );
            }
            if( i < fRunPara->fDF_UpSample.size() && iTemp.size() > 3 )
            {
                fRunPara->fDF_UpSample[i] = ( unsigned int )atoi( iTemp[3].c_str() );
            }
            if( i < fRunPara->fDF_PoleZero.size() && iTemp.size() > 4 )
            {
                fRunPara->fDF_PoleZero[i] = atof( iTemp[4].c_str() );
            }
        }
    }
    return true;
}

/*

    read all key words for double pass trace integration from runparameter file

*/
bool VEvndispReconstructionParameter::readKeyWord_FADCDOUBLEPASS( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 6 )
    {
        return false;
    }
    // on/off switch for double pass method
    if( iTemp[1].size() > 0 )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fDoublePass.size() )
                {
                    fRunPara->fDoublePass[i] = atoi( iTemp[1].c_str() );
                }
            }
        }
    }
    // length of summation window for pass 1
    if( iTemp[2].size() > 0 )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fsumwindow_pass1.size() )
                {
                    fRunPara->fsumwindow_pass1[i] = atoi( iTemp[2].c_str() );
                }
            }
        }
    }
    // trace integration method for pass 1
    if( iTemp[3].size() > 0 )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fTraceIntegrationMethod_pass1.size() )
                {
                    fRunPara->fTraceIntegrationMethod_pass1[i] = ( unsigned int )atoi( iTemp[3].c_str() );
                }
            }
        }
    }
    // set double pass error option
    if( iTemp[4].size() > 0 )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fDoublePassErrorWeighting2005.size() )
                {
                    fRunPara->fDoublePassErrorWeighting2005[i] = ( unsigned int )atoi( iTemp[4].c_str() );
                }
            }
        }
    }
    // set maximal allowed time difference (in samples) between lg and hg T0
    if( iTemp[5].size() > 0 )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fSumWindowMaxTimeDifferenceLGtoHG.size() )
                {
                    fRunPara->fSumWindowMaxTimeDifferenceLGtoHG[i] = atof( iTemp[5].c_str() );
                }
            }
        }
    }
    
    return true;
}

/*

      read key words for summation window setting from runparameter file

*/
bool VEvndispReconstructionParameter::readKeyWord_FADCSUMMATIONWINDOW( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 4 )
    {
        return false;
    }
    for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
    {
        if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
        {
            if( i < fRunPara->fsumwindow_1.size() && iTemp[1].size() > 0 )
            {
                fRunPara->fsumwindow_1[i] = atoi( iTemp[1].c_str() );
            }
            if( iTemp[2].size() > 0 )
            {
                if( i < fRunPara->fsumwindow_2.size() )
                {
                    fRunPara->fsumwindow_2[i] = atoi( iTemp[2].c_str() );
                }
            }
            else
            {
                // window 1 and 2 are the same unless stated differently
                if( i < fRunPara->fsumwindow_2.size() && iTemp[3].size() > 0 )
                {
                    fRunPara->fsumwindow_2[i] = atoi( iTemp[3].c_str() );
                }
            }
        }
    }
    
    return true;
}

/*

      read key words for summation window start from runparameter file

*/
bool VEvndispReconstructionParameter::readKeyWord_FADCSUMMATIONSTART( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 4 )
    {
        return false;
    }
    for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
    {
        if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
        {
            if( i < fRunPara->fsumfirst.size() && iTemp[1].size() > 0 )
            {
                fRunPara->fsumfirst[i] = atoi( iTemp[1].c_str() );
            }
            if( iTemp[2].size() > 0 )
            {
                if( i < fRunPara->fTraceWindowShift.size() )
                {
                    fRunPara->fTraceWindowShift[i] = atof( iTemp[2].c_str() );
                }
            }
            if( iTemp[3].size() > 0 )
            {
                if( i < fRunPara->fsumfirst_startingMethod.size() )
                {
                    // use T0 for timing of pulse integration
                    if( iTemp[3] == "T0" || iTemp[3] == "TZERO" || atoi( iTemp[3].c_str() ) == 1 )
                    {
                        fRunPara->fsumfirst_startingMethod[i] = 1;
                    }
                    // use average pulse arrival time for timing of pulse integration
                    else if( iTemp[3] == "TAVERAGE" || atoi( iTemp[3].c_str() ) == 2 )
                    {
                        fRunPara->fsumfirst_startingMethod[i] = 2;
                    }
                    // pixel trigger time
                    else if( iTemp[3] == "TTRIGGER" || atoi( iTemp[3].c_str() ) == 3 )
                    {
                           fRunPara->fsumfirst_startingMethod[i] = 3;
                    }
                    // fixed window start
                    else if( iTemp[3] == "FIXED" || atoi( iTemp[3].c_str() ) == 0 )
                    {
                        fRunPara->fsumfirst_startingMethod[i] = 0;
                    }
                    // default
                    else
                    {
                        cout << "VEvndispReconstructionParameter::read_arrayAnalysisCuts error:";
                        cout << " unknown timing method used for calculation of window start";
                        cout << " (valid parameters are TZERO/TAVERAGE/FIXED/TTRIGGER): ";
                        cout << iTemp[3] << endl;
                        fRunPara->fsumfirst_startingMethod[i] = 1;
                        cout << "...exiting" << endl;
                        exit( EXIT_FAILURE );
                    }
                }
            }
            if( iTemp[4].size() > 0 )
            {
                if( i < fRunPara->fsumfirst_maxT0startDiff.size() )
                {
                    fRunPara->fsumfirst_maxT0startDiff[i] = atof( iTemp[4].c_str() );
                }
            }
            if( iTemp[5].size() > 0 )
            {
                if( i < fRunPara->fSearchWindowLast.size() )
                {
                    fRunPara->fSearchWindowLast[i] = ( unsigned int )( atoi( iTemp[5].c_str() ) );
                }
            }
        }
    }
    return true;
}

/*

      read key words for cleaning from runparameter file

*/
bool VEvndispReconstructionParameter::readKeyWord_CLEANING( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 4 )
    {
        return false;
    }
    // image cleaning for double pass, pass1 (expext _DP in key word)
    if( iTemp[0].find( "_DP" ) != string::npos )
    {
        return fillImageCleaningParameter( iTemp, t_temp, fRunPara->fImageCleaningParameters_DB_Pass1 );
    }
    // image cleaning for main pass
    else
    {
        return fillImageCleaningParameter( iTemp, t_temp, fRunPara->fImageCleaningParameters );
    }
    return true;
}

bool VEvndispReconstructionParameter::readKeyWord_BRIGHTSTARS( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 5 )
    {
        return false;
    }
    fRunPara->fStarCatalogueName = iTemp[1];
    if( iTemp[2].size() > 0 )
    {
        fRunPara->fMinStarBrightness_B = atof( iTemp[2].c_str() );
    }
    if( iTemp[3].size() > 0 )
    {
        fRunPara->fMinStarPixelDistance_deg = atof( iTemp[3].c_str() );
    }
    if( iTemp[4].size() > 0 )
    {
        fRunPara->fMinStarNTubes = atoi( iTemp[4].c_str() );
    }
    return true;
}

bool VEvndispReconstructionParameter::readKeyWord_NEIGHBOURS( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
    {
        if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
        {
            if( i < fRunPara->fNeighbourDistanceFactor.size() && iTemp[1].size() > 0 )
            {
                fRunPara->fNeighbourDistanceFactor[i] = atof( iTemp[1].c_str() );
            }
            if( iTemp.size() >= 3 )
            {
                if( i < fRunPara->fNeighbourMultiplicity.size() && iTemp[2].size() > 0 )
                {
                    fRunPara->fNeighbourMultiplicity[i] = atoi( iTemp[2].c_str() );
                }
            }
        }
        if( iTemp.size() == 2 )
        {
            if( i < fRunPara->fNeighbourMultiplicity.size() && iTemp[2].size() > 0 )
            {
                fRunPara->fNeighbourMultiplicity[i] = atoi( iTemp[2].c_str() );
            }
        }
    }
    return true;
}

bool VEvndispReconstructionParameter::readKeyWord_SQUARE( vector< string > iTemp, int t_temp )
{
        if( !fRunPara || iTemp.size() < 2 )
        {
        return false;
        }  
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
            {
                if( i < fRunPara->fSquarePixels.size() && iTemp[1].size() > 0 )
                {
                    if (iTemp[1] == "TRUE")
                    {
                        fRunPara->fSquarePixels[i] = true;
                    }
                    else
                    {
                        fRunPara->fSquarePixels[i] = false;
                    }
                }
            }
        }
        return true;
}

bool VEvndispReconstructionParameter::readKeyWord_LLEDGEFIT( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 3 )
    {
        return false;
    }
    for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
    {
        if( t_temp < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_temp )
        {
            if( i < fRunPara->fLogLikelihoodLoss_min.size() && iTemp[1].size() > 0 )
            {
                fRunPara->fLogLikelihoodLoss_min[i] = atof( iTemp[1].c_str() );
            }
            if( iTemp[2].size() > 0 )
            {
                if( i < fRunPara->fLogLikelihood_Ntubes_min.size() )
                {
                    fRunPara->fLogLikelihood_Ntubes_min[i] = atoi( iTemp[2].c_str() );
                }
            }
            if( iTemp.size() > 3 && iTemp[3].size() > 0
                    && i < fRunPara->fLogLikelihoodLoss_max.size() )
            {
                fRunPara->fLogLikelihoodLoss_max[i] = atof( iTemp[3].c_str() );
            }
            
        }
    }
    return true;
}

/* FORCELL
   Set in EVNDISP.reconstruction.runparameter , for example:
   * -1 FORCELL 1
   if set to 1, the imagefitting will use log-likelihood image fitting for all
   images, regardless of whether the image is on the edge or not.
   This will ignore the options specified by 'LLEDGEFIT'
   If option is set to 0 or is not present, will behave normally
  (i.e. LL image fit on edges, hillas ellipse in middle of camera)
*/
bool VEvndispReconstructionParameter::readKeyWord_FORCELL( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    if( iTemp[1].size() > 0 && atoi( iTemp[1].c_str() ) == 1 )
    {
        fRunPara->fForceLLImageFit = true ;
    }
    else
    {
        fRunPara->fForceLLImageFit = false ;
    }
    cout << endl;
    cout << "FORCELL set to " << fRunPara->fForceLLImageFit << endl;
    cout << endl;
    return true;
}

/* CreateIPRdatabase

   Set in EVNDISP.reconstruction.runparameter, for example:
   * -1 CreateIPRdatabase TRUE
   If set to TRUE,
   create IPR database from DST file containing ped charge histos (hpedPerTelescopeType) and write it to the IPRdatabaseFile
*/
bool VEvndispReconstructionParameter::readKeyWord_CreateIPRdatabase( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    if( iTemp[1].size() > 0 && iTemp[1] == "TRUE" )
    {
        fRunPara->ifCreateIPRdatabase = true ;
    }
    else
    {
        fRunPara->ifCreateIPRdatabase = false ;
    }
    
    cout << std::boolalpha << "\n------------------------------" << endl;
    cout << "NN cleaning general settings:" << endl;
    cout << "\nCreateIPRdatabase set to " << fRunPara->ifCreateIPRdatabase << endl;
    
    return true;
}
/*  IPRdatabaseFile

    Set IPRdatabaseFile in EVNDISP.reconstruction.runparameter
    Usually set to EVNDISP_ANALYSIS_DIRECTORY so it will be written to this directory
    with output file <runnumber>.IPR.root)
*/
bool VEvndispReconstructionParameter::readKeyWord_IPRdatabaseFile( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    
    if( iTemp[1].size() > 0 )
    {
        fRunPara->fIPRdatabaseFile = iTemp[1];
        if( fRunPara->fIPRdatabaseFile == "EVNDISP_ANALYSIS_DIRECTORY" )
        {
            stringstream i_ssIPR;
            i_ssIPR << fRunPara->getDirectory_EVNDISPCalibrationData_perRun();
            i_ssIPR << fRunPara->frunnumber;
            i_ssIPR << ".IPR.root";
            fRunPara->fIPRdatabaseFile = i_ssIPR.str();
        }
    }
    else
    {
        cout << "VEvndispReconstructionParameter:: error:readKeyWord_IPRdatabaseFile";
        cout << " Empty value given to IPRdatabaseFile";
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    // (not entirely correct, as it depends on the sequences of
    // occurence of the parameters in the parameter file)
    if( fRunPara->ifCreateIPRdatabase )
    {
        cout << "IPRdatabaseFile set to " << fRunPara->fIPRdatabaseFile << endl;
    }
    
    return true;
}
/* ReadIPRfromDST

   Set ReadIPRfromDST in EVNDISP.reconstruction.runparameter, for example:
   * -1 ReadIPRfromDST TRUE
   If set to TRUE,
   read IPRs from DST file
*/
bool VEvndispReconstructionParameter::readKeyWord_ReadIPRfromDST( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    if( iTemp[1].size() > 0 && iTemp[1] == "TRUE" )
    {
        fRunPara->ifReadIPRfromDSTFile = true ;
    }
    else
    {
        fRunPara->ifReadIPRfromDSTFile = false ;
    }
    
    cout << "ReadIPRfromDST set to " << fRunPara->ifReadIPRfromDSTFile << endl;
    
    return true;
}
/* ReadIPRfromDatabase

   Set ReadIPRfromDatabase in EVNDISP.reconstruction.runparameter, for example:
   * -1 ReadIPRfromDatabase TRUE
   If set to TRUE,
   read IPRs from IPRdatabase (overrides CreateIPRdatabase)
*/
bool VEvndispReconstructionParameter::readKeyWord_ReadIPRfromDatabase( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    if( iTemp[1].size() > 0 && iTemp[1] == "TRUE" )
    {
        fRunPara->ifReadIPRfromDatabase = true ;
    }
    else
    {
        fRunPara->ifReadIPRfromDatabase = false ;
    }
    
    cout << "ReadIPRfromDatabase set to " << fRunPara->ifReadIPRfromDatabase << endl;
    
    return true;
}
/* IPRdatabase

   Set IPRdatabase in EVNDISP.reconstruction.runparameter
   Usually set to EVNDISP_ANALYSIS_DIRECTORY so it will be written to this directory
   with output file <runnumber>.IPR.root)
*/
bool VEvndispReconstructionParameter::readKeyWord_IPRdatabase( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    
    if( iTemp[1].size() > 0 )
    {
        fRunPara->fIPRdatabase = iTemp[1];
        if( fRunPara->fIPRdatabase == "EVNDISP_ANALYSIS_DIRECTORY" )
        {
            stringstream i_ssIPR;
            i_ssIPR << fRunPara->getDirectory_EVNDISPCalibrationData_perRun();
            i_ssIPR << fRunPara->frunnumber;
            i_ssIPR << ".IPR.root";
            fRunPara->fIPRdatabase = i_ssIPR.str();
        }
    }
    else
    {
        cout << "VEvndispReconstructionParameter:: error:readKeyWord_IPRdatabase";
        cout << " Empty value given to IPRdatabase";
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    cout << "IPRdatabase set to " << fRunPara->fIPRdatabase << endl;
    
    return true;
}
/* WriteGraphsToFile

   Set WriteGraphsToFile in EVNDISP.reconstruction.runparameter, for example:
   * -1 WriteGraphsToFile TRUE
   If set to TRUE, write NN Image cleaning graphs (prob. curves, IPR graphs, etc..) to a GraphsFile file
*/
bool VEvndispReconstructionParameter::readKeyWord_WriteGraphsToFile( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    if( iTemp[1].size() > 0 && iTemp[1] == "TRUE" )
    {
        fRunPara->ifWriteGraphsToFile = true ;
    }
    else
    {
        fRunPara->ifWriteGraphsToFile = false ;
    }
    
    cout << "WriteGraphsToFile set to " << fRunPara->ifWriteGraphsToFile << endl;
    
    return true;
}
/*  GraphsFile

    Set GraphsFile in EVNDISP.reconstruction.runparameter
    Usually set to EVNDISP_ANALYSIS_DIRECTORY so it will be written to this directory
    with output file <runnumber>.IPR.root)
*/
bool VEvndispReconstructionParameter::readKeyWord_GraphsFile( vector< string > iTemp, int t_temp )
{
    if( !fRunPara || iTemp.size() < 2 )
    {
        return false;
    }
    
    if( iTemp[1].size() > 0 )
    {
        fRunPara->fNNGraphsFile = iTemp[1];
        if( fRunPara->fNNGraphsFile == "EVNDISP_ANALYSIS_DIRECTORY" )
        {
            stringstream i_ssIPR;
            i_ssIPR << fRunPara->getDirectory_EVNDISPCalibrationData_perRun();
            i_ssIPR << fRunPara->frunnumber;
            i_ssIPR << ".IPRcontours.root";
            fRunPara->fNNGraphsFile = i_ssIPR.str();
        }
    }
    else
    {
        cout << "VEvndispReconstructionParameter:: error:readKeyWord_GraphsFile";
        cout << " Empty value given to GraphsFile";
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    cout << "GraphsFile set to " << fRunPara->fNNGraphsFile << endl;
    cout << std::noboolalpha << "------------------------------\n" << endl;
    
    return true;
}

/*

    read reconstruction (quality) cuts

*/
bool VEvndispReconstructionParameter::readKeyWord_RECMETHOD( vector< string > iTemp, int t_temp, ULong64_t t_type )
{
    if( !fRunPara || iTemp.size() < 4 )
    {
        return false;
    }
    if( iTemp[0] == "USEEVNPOINTING" )
    {
        fReconstructionParameterData.back()->fUseEventdisplayPointing = true;
        return true;
    }
    else if( iTemp[0] == "MNIMAGE" && iTemp[1].size() > 0 )
    {
        fReconstructionParameterData.back()->fNImages_min = atoi( iTemp[1].c_str() );
        return true;
    }
    else if( iTemp[0] == "MINANGLE" && iTemp[1].size() > 0 )
    {
        fReconstructionParameterData.back()->fAxesAngles_min = atof( iTemp[1].c_str() );
        return true;
    }
    else if( iTemp[0] == "TMVABDTFILE" && iTemp[1].size() > 0 )
    {
        fReconstructionParameterData.back()->fDISP_TMVAZenithBin.push_back( atof( iTemp[1].c_str() ) );
        if( iTemp[2].size() > 0 )
        {
            fReconstructionParameterData.back()->fDISP_TMVAFileNameVector.push_back( iTemp[2] );
        }
        // use disp up to this number of images
        if( iTemp[3].size() > 0 )
        {
            fReconstructionParameterData.back()->fDISP_MAXTelescopes.push_back( atoi( iTemp[3].c_str() ) );
        }
        // (default value is 4)
        else
        {
            fReconstructionParameterData.back()->fDISP_MAXTelescopes.push_back( 4 );
        }
        // above the maximum number of images: use this reconstruction method
        if( iTemp.size() > 3 && iTemp[4].size() > 0 )
        {
            fReconstructionParameterData.back()->fDISP_MAXMethodID.push_back( atoi( iTemp[4].c_str() ) );
        }
        // (default value is 4)
        else
        {
            fReconstructionParameterData.back()->fDISP_MAXMethodID.push_back( 4 );
        }
        
    }
    else if( iTemp[0] == "USEIMAGE" )
    {
        if( t_temp >= 0 )                 // use defaults for telescope number < 0
        {
            if( fReconstructionParameterData.back()->isValidTelescopeType( t_type ) )
            {
                unsigned int t = fReconstructionParameterData.back()->getTelescopeType_counter( t_type );
                fReconstructionParameterData.back()->fLocalUseImage[t] = ( bool )atoi( iTemp[1].c_str() );
            }
        }
    }
    // all min/max cuts
    // note: strict matching of key words
    // always: KEYWORD MIN MAX
    else if( iTemp[1].size() > 0 && iTemp[2].size() > 0 )
    {
        // is this a valid key word for a double value?
        if( fReconstructionParameterData.back()->isValidKeyword( iTemp[0] ) )
        {
            int i_ntubes_min = -99999;
            int i_ntubes_max = -99999;
            if( iTemp[3].size() > 0 )
            {
                i_ntubes_min = atoi( iTemp[3].c_str() );
            }
            if( iTemp[4].size() > 0 )
            {
                i_ntubes_max = atoi( iTemp[4].c_str() );
            }
            // set this cut for all telescope types
            if( t_temp < 0 )
            {
                for( unsigned int t = 0;
                        t < fReconstructionParameterData.back()->fTelescopeTypeCut.size();
                        t++ )
                {
                    // double
                    if( iTemp[1].find( "." ) != string::npos || iTemp[2].find( "." ) != string::npos )
                    {
                        fReconstructionParameterData.back()->fTelescopeTypeCut[t][iTemp[0]]
                        ->setCutValues( atof( iTemp[1].c_str() ),
                                        atof( iTemp[2].c_str() ),
                                        i_ntubes_min, i_ntubes_max );
                    }
                    // integer cut
                    else
                    {
                        fReconstructionParameterData.back()->fTelescopeTypeCut[t][iTemp[0]]
                        ->setCutValues( atoi( iTemp[1].c_str() ),
                                        atoi( iTemp[2].c_str() ),
                                        i_ntubes_min, i_ntubes_max );
                    }
                }
            }
            //////////////////////////////////////////////////////////
            // set values for a particular telescope type
            else if( fReconstructionParameterData.back()->isValidTelescopeType( t_type ) )
            {
                unsigned int t = fReconstructionParameterData.back()->getTelescopeType_counter( t_type );
                // double
                if( iTemp[1].find( "." ) != string::npos || iTemp[2].find( "." ) != string::npos )
                {
                    fReconstructionParameterData.back()->fTelescopeTypeCut[t][iTemp[0]]
                    ->setCutValues( atof( iTemp[1].c_str() ),
                                    atof( iTemp[2].c_str() ),
                                    i_ntubes_min, i_ntubes_max );
                }
                // integer cut
                else
                {
                    fReconstructionParameterData.back()->fTelescopeTypeCut[t][iTemp[0]]
                    ->setCutValues( atoi( iTemp[1].c_str() ),
                                    atoi( iTemp[2].c_str() ),
                                    i_ntubes_min, i_ntubes_max );
                }
            }
            else
            {
                cout << "\t VEvndispReconstructionParameter::read_arrayAnalysisCuts warning:";
                cout << " unknown telescope type " << t_type << endl;
            }
        }
    }  // if( iTemp[1].size() > 0 && iTemp[2].size() > 0 )
    return true;
}

/*!

   A RECMETHOD LINE STARTS ALWAYS A NEW RECORD

   line without '*' in the beginning are ignored

*/
unsigned int VEvndispReconstructionParameter::read_arrayAnalysisCuts( string ifile )
{
    if( ifile.size() == 0 )
    {
        return 0;
    }
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VEvndispReconstructionParameter::read_arrayAnalysisCuts error while opening array analysis cut file: " << ifile << endl;
        return 0;
    }
    cout << endl;
    cout << "reading reconstruction parameters (e.g. image quality cuts) from: " << endl;
    cout << ifile << endl;
    string iLine;
    vector< string > iTemp( 6, "" );
    ULong64_t t_type = 0;
    int t_temp = 0;
    vector< int > v_temp;
    int m_temp = -1;
    
    ////////////////////////////////////////////////////////////
    // loop over all files in the reconstruction parameter file
    //
    // (lines without '*' in the first column are ignored)
    while( getline( is, iLine ) )
    {
        // line without '*' in the beginning are ignored
        if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*" )
        {
            iTemp.assign( 6, "" );
            istringstream is_stream( iLine );
            is_stream >> iTemp[0];
            // proceed
            is_stream >> iTemp[0];
            ////////////////////////////////////////////////////////////////////////
            // telescope types (horrible...)
            v_temp.clear();
            if( atoi( iTemp[0].c_str() ) >= 0 )
            {
                t_type = ULong64_t( atoi( iTemp[0].c_str() ) );
                t_temp = getTelescopeType_counter( t_type );
                v_temp = getTelescopeType_counterVector( t_type );
            }
            else if( atoi( iTemp[0].c_str() ) < -10 &&  atoi( iTemp[0].c_str() ) > -1000 )
            {
                // get telescope type counter (0,1,2...) for this type
                t_temp = getTelescopeType_counter_from_MirrorArea( ULong64_t( -1 * atoi( iTemp[0].c_str() ) ) );
                v_temp = getTelescopeType_counter_from_MirrorAreaVector( ULong64_t( -1 * atoi( iTemp[0].c_str() ) ) );
                if( t_temp >= 0 && t_temp < ( int )fTel_typeUniqueVector.size() )
                {
                    t_type = fTel_typeUniqueVector[t_temp];
                }
                
            }
            else if( atoi( iTemp[0].c_str() ) < -1000 )
            {
                t_temp = getTelescopeType_counter_from_MirrorArea_and_PixelSize( ULong64_t( -1 * atoi( iTemp[0].c_str() ) ) );
                v_temp = getTelescopeType_counter_from_MirrorArea_and_PixelSizeVector( ULong64_t( -1 * atoi( iTemp[0].c_str() ) ) );
                if( t_temp >= 0 && t_temp < ( int )fTel_typeUniqueVector.size() )
                {
                    t_type = fTel_typeUniqueVector[t_temp];
                }
            }
            else
            {
                t_temp = -1;
            }
            // unknown telescope types
            if( t_temp == -2 )
            {
                continue;
            }
            
            ///////////////////////////////////
            // read variable identifieres
            for( unsigned int i = 0; i < iTemp.size(); i++ )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> iTemp[i];
                }
                else
                {
                    break;
                }
                if( i == 0 )
                {
                    iTemp[0] = VUtilities::upperCase( iTemp[0] );
                }
            }
            //////////////////////////////////////////////////////////////////////////////////////////////
            // fadc trace analysis
            if( iTemp[0] == "FADCANALYSIS" && fRunPara )
            {
                readKeyWord_FADCANALYSIS( iTemp, t_temp );
                continue;
            }
            // double pass options
            else if( iTemp[0] == "FADCDOUBLEPASS" && fRunPara )
            {
                readKeyWord_FADCDOUBLEPASS( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "FADCSUMMATIONWINDOW" && fRunPara )
            {
                readKeyWord_FADCSUMMATIONWINDOW( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "FADCSUMMATIONSTART" && fRunPara )
            {
                readKeyWord_FADCSUMMATIONSTART( iTemp, t_temp );
                continue;
            }
            // image cleaning
            else if( iTemp[0].find( "CLEANING" ) != string::npos || iTemp[0].find( "TIMETWOLEVELPARAMETERS" ) != string::npos )
            {
                readKeyWord_CLEANING( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "BRIGHTSTARS" && fRunPara )
            {
                readKeyWord_BRIGHTSTARS( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "LLEDGEFIT" && fRunPara )
            {
                readKeyWord_LLEDGEFIT( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "NEIGHBOURDISTANCEFACTOR" && fRunPara )
            {
                readKeyWord_NEIGHBOURS( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "SQUAREPIXELS" && fRunPara )
            {
                readKeyWord_SQUARE( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "FORCELL" && fRunPara )
            {
                readKeyWord_FORCELL( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "CREATEIPRDATABASE" && fRunPara )
            {
                readKeyWord_CreateIPRdatabase( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "IPRDATABASEFILE" && fRunPara )
            {
                readKeyWord_IPRdatabaseFile( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "READIPRFROMDST" && fRunPara )
            {
                readKeyWord_ReadIPRfromDST( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "READIPRFROMDATABASE" && fRunPara )
            {
                readKeyWord_ReadIPRfromDatabase( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "IPRDATABASE" && fRunPara )
            {
                readKeyWord_IPRdatabase( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "WRITEGRAPHSTOFILE" && fRunPara )
            {
                readKeyWord_WriteGraphsToFile( iTemp, t_temp );
                continue;
            }
            else if( iTemp[0] == "GRAPHSFILE" && fRunPara )
            {
                readKeyWord_GraphsFile( iTemp, t_temp );
                continue;
            }
            // Model3D: reconstruction ID for starting values
            else if( iTemp[0] == "MODEL3DSTARTID" && fRunPara && iTemp[1].size() > 0 )
            {
                fRunPara->fIDstartDirectionModel3D = atoi( iTemp[1].c_str() );
                continue;
            }
            else if( iTemp[0] == "CreateIPRdatabase" && fRunPara )
            {
                readKeyWord_FORCELL( iTemp, t_temp );
                continue;
            }
            
            /////////////////////////////////////////////////
            // check for exit statement
            if( iTemp[0] == "EXIT" )
            {
                return fReconstructionParameterData.size();
            }
            // check for non-MC exit statement
            if( iTemp[0] == "MCONLY" )
            {
                if( !fRunPara->isMC() )
                {
                    return fReconstructionParameterData.size();
                }
                else
                {
                    continue;
                }
            }
            /////////////////////////////////////////////////
            /////////////////////////////////////////////////
            /////////////////////////////////////////////////
            // A NEW RECORD OF QUALITY CUTS USED IN THE
            // STEREO RECONSTRUCTION STARTS ALWAYS WITH 'RECMETHOD'
            /////////////////////////////////////////////////
            /////////////////////////////////////////////////
            if( iTemp[0] == "RECMETHOD" )
            {
                m_temp = fReconstructionParameterData.size();
                // reset all parameters for this method number
                addNewMethod( m_temp, ( unsigned int )atoi( iTemp[1].c_str() ) );
                continue;
            }
            // check if a first record was created
            if( m_temp < 0 || fReconstructionParameterData.size() == 0 )
            {
                cout << "VEvndispReconstructionParameter::read_arrayAnalysisCuts error: no valid set of cuts found " << endl;
                cout << "invalid line is: " << endl;
                cout << "\t " << iLine << endl;
                cout << "start set of cuts with the line: " << endl;
                cout << "* -1 RECMETHOD 0" << endl;
                exit( EXIT_FAILURE );
            }
            if( fReconstructionParameterData.size() > 0 )
            {
                readKeyWord_RECMETHOD( iTemp, t_temp, t_type );
                continue;
            }
            
            else
            {
                cout << "\t VEvndispReconstructionParameter::read_arrayAnalysisCuts warning: unknown identifier: " << iTemp[0] << endl;
            }
        }                                         // if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*" )
    }                                             // while( getline( is, iLine ) )
    return fReconstructionParameterData.size();
}


int VEvndispReconstructionParameter::getTelescopeType_counter( ULong64_t t )
{
    unsigned int z = 0;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        if( *fTel_type_iter == t )
        {
            return z;
        }
        z++;
    }
    return -2;
}

/*
   use mirror area flag to identify array analysis cut
*/
int VEvndispReconstructionParameter::getTelescopeType_counter_from_MirrorArea( ULong64_t t )
{
    unsigned int z = 0;
    ULong64_t v = 0;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        v = *fTel_type_iter;
        v /= 100000;
        
        if( v > 2000 )
        {
            v -= 2000;
        }
        if( v > 1000 )
        {
            v -= 1000;
        }
        if( v == t )
        {
            return z;
        }
        z++;
    }
    return -2;
}

/*
   use mirror area and pixel size to identify array analysis cut
*/
int VEvndispReconstructionParameter::getTelescopeType_counter_from_MirrorArea_and_PixelSize( ULong64_t t )
{
    unsigned int z = 0;
    ULong64_t v = 0;
    ULong64_t v2 = 0;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        v = *fTel_type_iter;
        
        ULong64_t v1 = v / 100000;
        if( v1 > 1000 )
        {
            v1 -= 1000;
        }
        v2 = v1 * 100 + ( v % 100 );
        
        if( v2 == t )
        {
            return z;
        }
        z++;
    }
    return -2;
}


vector< int > VEvndispReconstructionParameter::getTelescopeType_counterVector( ULong64_t t )
{
    unsigned int z = 0;
    vector< int > v;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        if( *fTel_type_iter == t )
        {
            v.push_back( ( int )z );
        }
        z++;
    }
    return v;
}

/*
   use mirror area flag to identify array analysis cut
*/
vector< int > VEvndispReconstructionParameter::getTelescopeType_counter_from_MirrorAreaVector( ULong64_t t )
{
    unsigned int z = 0;
    ULong64_t v = 0;
    vector< int > x;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); fTel_type_iter++ )
    {
        v = *fTel_type_iter;
        v /= 100000;
        
        if( v > 2000 )
        {
            v -= 2000;
        }
        if( v > 1000 )
        {
            v -= 1000;
        }
        if( v == t )
        {
            x.push_back( z );
        }
        z++;
    }
    return x;
}

/*
   use mirror area and pixel size to identify array analysis cut
*/
vector< int > VEvndispReconstructionParameter::getTelescopeType_counter_from_MirrorArea_and_PixelSizeVector( ULong64_t t )
{
    unsigned int z = 0;
    vector< int > x;
    ULong64_t v = 0;
    ULong64_t v2 = 0;
    set< ULong64_t >::iterator fTel_type_iter;
    for( fTel_type_iter = fTel_type.begin(); fTel_type_iter != fTel_type.end(); ++fTel_type_iter )
    {
        v = *fTel_type_iter;
        
        ULong64_t v1 = v / 100000;
        if( v1 > 1000 )
        {
            v1 -= 1000;
        }
        v2 = v1 * 100 + ( v % 100 );
        
        if( v2 == t )
        {
            x.push_back( z );
        }
        z++;
    }
    return x;
}


bool VEvndispReconstructionParameter::fillImageCleaningParameter( vector< string > iTemp, int t_telescopeCounter,
        vector< VImageCleaningRunParameter* >& iImageCleaningParameters )
{
    if( iTemp[0].find( "IMAGECLEANINGMETHOD" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() )
                {
                    if( !iImageCleaningParameters[i]->setImageCleaningMethod( iTemp[1] ) )
                    {
                        cout << "VEvndispReconstructionParameter: unknown image cleaning method: " << iTemp[1] << endl;
                    }
                    if( iTemp[2].size() > 0 )
                    {
                        if( iTemp[2] == "FIXED" )
                        {
                            iImageCleaningParameters[i]->fUseFixedThresholds = true;
                        }
                        else
                        {
                            iImageCleaningParameters[i]->fUseFixedThresholds = false;
                        }
                    }
                }
            }
        }
        return true;
    }
    else if( iTemp[0].find( "IMAGECLEANINGTHRESHOLDS" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->fimagethresh = atof( iTemp[1].c_str() );
                }
                if( iTemp[2].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->fborderthresh = atof( iTemp[2].c_str() );
                    }
                }
            }
        }
        return true;
    }
    //////////////////////////////////////////////////////
    // NN cleaning - fake image cleaning probability
    else if( iTemp[0].find( "IMAGECLEANING_FAKEPROBABILITY" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->fNNOpt_FakeImageProb = atof( iTemp[1].c_str() );
                }
            }
        }
        return true;
    }
    //////////////////////////////////////////////////////
    // NN cleaning - active multiplicities
    else if( iTemp[0].find( "IMAGECLEANING_ACTIVEMULTIPLICITIES" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() )
                {
                    for( unsigned int f = 0; f < iTemp.size(); f++ )
                    {
                        if( iTemp[f].size() > 0 )
                        {
                            for( unsigned int m = 0; m < iImageCleaningParameters[i]->fNNOpt_Multiplicities.size(); m++ )
                            {
                                if( iTemp[f] == iImageCleaningParameters[i]->fNNOpt_Multiplicities[m] )
                                {
                                    iImageCleaningParameters[i]->fNNOpt_ActiveNN[m] = true;
                                }
                            }
                        }
                    }
                    // (not used; always true)
                    iImageCleaningParameters[i]->fNNOpt_ActiveNN[4] = true;
                }
            }
        }
    }
    //
    else if( iTemp[0].find( "TIMECLEANINGPARAMETERS" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->ftimecutpixel = atof( iTemp[1].c_str() );
                }
                if( iTemp[2].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->ftimecutcluster = atof( iTemp[2].c_str() );
                    }
                }
                if( iTemp[3].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->fminpixelcluster = atoi( iTemp[3].c_str() );
                    }
                }
                if( iTemp[4].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->floops = atoi( iTemp[4].c_str() );
                    }
                }
            }
        }
        return true;
    }
    else if( iTemp[0].find( "CLUSTERCLEANINGPARAMETERS" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->fnmaxcluster = atof( iTemp[1].c_str() );
                }
                if( iTemp[2].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->fminsizecluster = atof( iTemp[2].c_str() );
                    }
                }
                if( iTemp[3].size() > 0 )
                {
                    if( i < iImageCleaningParameters.size() )
                    {
                        iImageCleaningParameters[i]->fminpixelcluster = atoi( iTemp[3].c_str() );
                    }
                }
            }
        }
        return true;
    }
    else if( iTemp[0].find( "CORRELATIONCLEANINGPARAMETER" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->fCorrelationCleanBoardThresh = atof( iTemp[1].c_str() );
                }
                if( iTemp[2].size() > 0 && i < iImageCleaningParameters.size() )
                {
                    iImageCleaningParameters[i]->fCorrelationCleanCorrelThresh = atof( iTemp[2].c_str() );
                }
                if( iTemp[3].size() > 0 && i < iImageCleaningParameters.size() )
                {
                    iImageCleaningParameters[i]->fCorrelationCleanNpixThresh = atoi( iTemp[3].c_str() );
                }
            }
        }
        return true;
    }
    else if( iTemp[0].find( "TIMETWOLEVELPARAMETERS" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    iImageCleaningParameters[i]->ftimediff = atof( iTemp[1].c_str() );
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Maximum number of rings to be searched for boundary pixels
    else if( iTemp[0].find( "IMAGECLEANING_NRINGS" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    unsigned int nRings = atoi( iTemp[1].c_str() );
                    if( 0 < nRings && nRings < 10 ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_nRings = nRings;
                    }
                    else
                    {
                        cout << "\t Error: Value of the maximum number of rings to be searched ";
                        cout << "for boundary pixels is out of scope, ";
                        cout << "defaulting to nRings = 3" << endl;
                        iImageCleaningParameters[i]->fNNOpt_nRings = 3;
                    }
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Maximal coincidence window allowed between neighbouring pixels for any NN group (in nano seconds)
    else if( iTemp[0].find( "IMAGECLEANING_COINCWINLIMIT" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    float CoincWinLimit = atof( iTemp[1].c_str() );
                    if( 0 < CoincWinLimit && CoincWinLimit < 100 ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_CoincWinLimit = CoincWinLimit;
                    }
                    else
                    {
                        cout << "\t Error: The maximum coincidence window cannot be " << CoincWinLimit;
                        cout << ", defaulting to CoincWinLimit = 8" << endl;
                        iImageCleaningParameters[i]->fNNOpt_CoincWinLimit = 8;
                    }
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Ignore timing information in NN image cleaning
    else if( iTemp[0].find( "IMAGECLEANING_NNOPTNOTIMEING" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    if( iTemp[1] == "TRUE" ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_ifNNoptNoTimeing = true;
                    }
                    else
                    {
                        iImageCleaningParameters[i]->fNNOpt_ifNNoptNoTimeing = false;
                    }
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Set the sample time slice and number of ADC bins to be read explicitly
    else if( iTemp[0].find( "IMAGECLEANING_SETEXPLICITSAMPLETIMESLICE" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    if( iTemp[1] == "TRUE" ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_ifExplicitSampleTimeSlice = true;
                    }
                    else
                    {
                        iImageCleaningParameters[i]->fNNOpt_ifExplicitSampleTimeSlice = false;
                    }
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Size of time slice in ns (usually 1 or 2 ns, only used if setExplicitSampleTimeSlice is TRUE)
    else if( iTemp[0].find( "IMAGECLEANING_SAMPLETIMESLICE" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    float sampleTimeSlice = atof( iTemp[1].c_str() );
                    if( 0 < sampleTimeSlice && sampleTimeSlice < 20 ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_sampleTimeSlice = sampleTimeSlice;
                    }
                    else
                    {
                        cout << "\t Error: The time slice interval cannot be " << sampleTimeSlice;
                        cout << ", defaulting to sampleTimeSlice = 1 ns" << endl;
                        iImageCleaningParameters[i]->fNNOpt_sampleTimeSlice = 1;
                    }
                }
            }
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // NN cleaning - Number of ADC bins summed up, each bin the size of sampleTimeSlice (only used if setExplicitSampleTimeSlice is TRUE)
    else if( iTemp[0].find( "IMAGECLEANING_NBINSADC" ) != string::npos )
    {
        for( unsigned int i = 0; i < fTel_type_perTelescope.size(); i++ )
        {
            if( t_telescopeCounter < 0 || getTelescopeType_counter( fTel_type_perTelescope[i] ) == t_telescopeCounter )
            {
                if( i < iImageCleaningParameters.size() && iTemp[1].size() > 0 )
                {
                    unsigned int nBinsADC = atoi( iTemp[1].c_str() );
                    if( 0 < nBinsADC && nBinsADC < 150 ) // Sanity check
                    {
                        iImageCleaningParameters[i]->fNNOpt_nBinsADC = nBinsADC;
                    }
                    else
                    {
                        cout << "\t Error: Value of the number of ADC bins to read is out of scope, ";
                        cout << "defaulting to nBinsADC = 50" << endl;
                        iImageCleaningParameters[i]->fNNOpt_nBinsADC = 50;
                    }
                }
            }
        }
        return true;
    }
    
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// data class for reconstruction parameters
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

VEvndispReconstructionParameterData::VEvndispReconstructionParameterData( unsigned int iNTel_type, set< ULong64_t > fTel_type )
{
    fDebug          = false;
    fNTel_type      = iNTel_type;
    
    for( set<ULong64_t>::iterator it = fTel_type.begin(); it != fTel_type.end(); ++it )
    {
        fTelescopeType.push_back( *it );
    }
    
    fMethodID       = 0;
    
    fNImages_min    = 0.;
    fAxesAngles_min = 0.;
    
    // cut parameters: set default parameters
    for( unsigned int i = 0; i < fNTel_type; i++ )
    {
        map< string, VEvndispReconstructionCut* > iTempM_C;
        iTempM_C.insert( make_pair( "NTUBES", new VEvndispReconstructionCut() ) );
        iTempM_C[ "NTUBES" ]->setCutValues( 2, 100000 );
        iTempM_C.insert( make_pair( "NLOWGAIN", new VEvndispReconstructionCut() ) );
        iTempM_C[ "NLOWGAIN" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "DISTANCE", new VEvndispReconstructionCut() ) );
        iTempM_C[ "DISTANCE" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "SIZE", new VEvndispReconstructionCut() ) );
        iTempM_C[ "SIZE" ]->setCutValues( 0., 1.e10 );
        iTempM_C.insert( make_pair( "WIDTH", new VEvndispReconstructionCut() ) );
        iTempM_C[ "WIDTH" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "LOSS1", new VEvndispReconstructionCut() ) );
        iTempM_C[ "LOSS1" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "LOSS2", new VEvndispReconstructionCut() ) );
        iTempM_C[ "LOSS2" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "WIDTHLENGTH", new VEvndispReconstructionCut() ) );
        iTempM_C[ "WIDTHLENGTH" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "FUI", new VEvndispReconstructionCut() ) );
        iTempM_C[ "FUI" ]->setCutValues( false );
        iTempM_C.insert( make_pair( "MCENERGY_LINTEV", new VEvndispReconstructionCut() ) );
        iTempM_C[ "MCENERGY_LINTEV" ]->setCutValues( false );
        
        fTelescopeTypeCut.push_back( iTempM_C );
        
        fLocalUseImage.push_back( true );
        fL2TriggerType.push_back( 9999 );
    }
    
    fUseEventdisplayPointing = false;
    
    // MLP parameters for array reconstruction
    fDISP_MLPFileName = "";
}

/*
 *
 *  note: cut on teltype!
 */
bool VEvndispReconstructionParameterData::testUserImage( unsigned int iTelType )
{
    if( iTelType < fLocalUseImage.size() && fLocalUseImage[iTelType] )
    {
        return true;
    }
    if( fDebug )
    {
        cout << "VEvndispReconstructionParameterData::testUserImage: image removed by user selection" << endl;
    }
    
    return false;
}

/*

       test a cut value for a given cut (iVarName)
       (double version)

*/
bool VEvndispReconstructionParameterData::test( unsigned int iTelType, string iVarName, double iVarD, int iNtubes )
{
    // first requirement: teltype must exist
    if( iTelType < fTelescopeTypeCut.size() )
    {
        // second requirement: cut must exist
        if( fTelescopeTypeCut[iTelType].find( iVarName ) != fTelescopeTypeCut[iTelType].end() )
        {
            // third requirement: value must pass cuts
            if( fTelescopeTypeCut[iTelType][iVarName]->test( iVarD, iNtubes ) )
            {
                if( fDebug )
                {
                    cout << "VEvndispReconstructionParameterData::test " << iVarName << " cut ";
                    cout << " (teltype " << iTelType << "): ";
                    cout << iVarD;
                    cout << " (" << fTelescopeTypeCut[iTelType][iVarName]->fCut_double_min << ",";
                    cout << fTelescopeTypeCut[iTelType][iVarName]->fCut_double_max;
                    cout << " ntubes [" << fTelescopeTypeCut[iTelType][iVarName]->fCut_ntubes_min;
                    cout << ", " << fTelescopeTypeCut[iTelType][iVarName]->fCut_ntubes_max << "] )";
                    cout << endl;
                }
                return true;
            }
        }
        else
        {
            cout << "VEvndispReconstructionParameterData::test error: unknown variable name " << iVarName << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    return false;
}

/*

       test a cut value for a given cut (iVarName)
       (int version)

*/
bool VEvndispReconstructionParameterData::test( unsigned int iTelType, string iVarName, int iVarI, int iNtubes )
{
    // first requirement: teltype must exist
    if( iTelType < fTelescopeTypeCut.size() )
    {
        // second requirement: cut must exist
        if( fTelescopeTypeCut[iTelType].find( iVarName ) != fTelescopeTypeCut[iTelType].end() )
        {
            // third requirement: value must pass cuts
            if( fTelescopeTypeCut[iTelType][iVarName]->test( iVarI, iNtubes ) )
            {
                if( fDebug )
                {
                    cout << "VEvndispReconstructionParameterData::test " << iVarName << " cut ";
                    cout << " (teltype " << iTelType << "): ";
                    cout << iVarI;
                    cout << " (" << fTelescopeTypeCut[iTelType][iVarName]->fCut_int_min << ",";
                    cout << fTelescopeTypeCut[iTelType][iVarName]->fCut_int_max;
                    cout << "ntubes [" << fTelescopeTypeCut[iTelType][iVarName]->fCut_ntubes_min;
                    cout << endl;
                }
                return true;
            }
        }
        else
        {
            cout << "VEvndispReconstructionParameterData::test error: unknown variable name " << iVarName << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    return false;
}

bool VEvndispReconstructionParameterData::testL2TriggerType( unsigned int iTel, unsigned int iTelType, unsigned short int iLocalTriggerType )
{
    // 9999: any trigger allowed
    if( iTelType < fL2TriggerType.size() && fL2TriggerType[iTelType] != 9999 )
    {
        bitset< 8 > i_L2TrigType( iLocalTriggerType );
        if( !i_L2TrigType.test( fL2TriggerType[iTelType] ) )
        {
            return false;
        }
        if( fDebug )
        {
            cout << "VEvndispReconstructionParameter::applyArrayAnalysisCut Tel " << iTel + 1 << ", type " << iTelType;
            cout << ": L2 trigger type " << iLocalTriggerType;
            cout << " test: " << i_L2TrigType.test( fL2TriggerType[iTelType] );
            cout << " (require " << fL2TriggerType[iTelType] << ")" << endl;
        }
    }
    return true;
}

void VEvndispReconstructionParameterData::setMethodID( unsigned int iMethodID )
{
    // hardwired: allowed array reconstruction numbers
    if( iMethodID != 0 && iMethodID != 3 && iMethodID != 4
            && iMethodID != 5 && iMethodID != 6 && iMethodID != 7 && iMethodID != 8 && iMethodID != 9 )
    {
        cout << "VEvndispReconstructionParameterData::setMethodID: invalid array reconstruction method: " << iMethodID << endl;
        cout << "(allowed is 0,3,4,5,6,7,8,9)" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    fMethodID = iMethodID;
}

unsigned int VEvndispReconstructionParameterData::getTelescopeType_counter( unsigned int iTelType )
{
    for( unsigned int i = 0; i < fTelescopeType.size(); i++ )
    {
        if( iTelType == fTelescopeType[i] )
        {
            return i;
        }
    }
    // hmm, is this right?
    return 0;
}

bool VEvndispReconstructionParameterData::isValidTelescopeType( unsigned int iTelType )
{
    for( unsigned int i = 0; i < fTelescopeType.size(); i++ )
    {
        if( iTelType == fTelescopeType[i] )
        {
            return true;
        }
    }
    return false;
}

bool VEvndispReconstructionParameterData::isValidKeyword( string iKeyword )
{
    if( fTelescopeTypeCut.size() > 0 )
    {
        map< string, VEvndispReconstructionCut* >::iterator i_iterator_intMap;
        for( i_iterator_intMap = fTelescopeTypeCut[0].begin();
                i_iterator_intMap != fTelescopeTypeCut[0].end();
                ++i_iterator_intMap )
        {
            if( iKeyword == i_iterator_intMap->first )
            {
                return true;
            }
        }
    }
    return false;
}

VEvndispReconstructionParameterData* VEvndispReconstructionParameter::getReconstructionParameterData( unsigned int iMethod )
{
    if( iMethod < fReconstructionParameterData.size() )
    {
        return fReconstructionParameterData[iMethod];
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

VEvndispReconstructionCut::VEvndispReconstructionCut()
{
    fCutSet = false;
    
    fCut_ntubes_min = -99999;
    fCut_ntubes_max =  99999;
    
    fCut_double_min = -9999.;
    fCut_double_max = -9999.;
    fCut_int_min    = -9999;
    fCut_int_max    = -9999;
}

void  VEvndispReconstructionCut::print( unsigned int telType, string iCutName )
{
    if( !fCutSet )
    {
        return;
    }
    
    if( fCut_double_min > -9998. && fCut_double_max > -9998. )
    {
        cout << "\t\t TelType " << telType;
        cout << " cut: " << iCutName << " [";
        cout << fCut_double_min << ", " << fCut_double_max << "]";
        if( fCut_ntubes_min > 0 )
        {
            cout << ", ntubes > " << fCut_ntubes_min;
        }
        if( fCut_ntubes_max > 0 )
        {
            cout << ", ntubes <= " << fCut_ntubes_max;
        }
        cout << endl;
    }
    else if( fCut_int_min > -9998. && fCut_int_max > -9998. )
    {
        cout << "\t\t TelType " << telType;
        cout << " cut: " << iCutName << " [";
        cout << fCut_int_min << ", " << fCut_int_max << "]";
        if( fCut_ntubes_min > 0 )
        {
            cout << ", ntubes > " << fCut_ntubes_min;
        }
        if( fCut_ntubes_max > 0 )
        {
            cout << ", ntubes <= " << fCut_ntubes_max;
        }
        cout << endl;
    }
}

void  VEvndispReconstructionCut::setCutValues( bool iCutSet )
{
    fCutSet = false;
}

/*
    set the cut values (with or without ntubes cut)
    (int version)
*/
void  VEvndispReconstructionCut::setCutValues( int iMin, int iMax, int ntubes_min, int ntubes_max )
{
    fCutSet = true;
    fCut_int_min = iMin;
    fCut_int_max = iMax;
    fCut_ntubes_min = ntubes_min;
    fCut_ntubes_max = ntubes_max;
}

/*
    set the cut values (with or without ntubes cut)
    (double version)
*/
void  VEvndispReconstructionCut::setCutValues( double iMin, double iMax, int ntubes_min, int ntubes_max )
{
    fCutSet = true;
    fCut_double_min = iMin;
    fCut_double_max = iMax;
    fCut_ntubes_min = ntubes_min;
    fCut_ntubes_max = ntubes_max;
}

/*
    apply a cut on the number of tubes per image

    return true if cut was passed

*/
bool  VEvndispReconstructionCut::testNtubes( int ntubes )
{
    if( ntubes > 0 && fCut_ntubes_min > 0 )
    {
        if( ntubes < fCut_ntubes_min )
        {
            return false;
        }
    }
    if( ntubes > 0 && fCut_ntubes_max > 0 )
    {
        if( ntubes >= fCut_ntubes_max )
        {
            return false;
        }
    }
    return true;
}

/*
      test this particular cut (int version)

      take also number of ntubes into account

      return true if cut was passed (or not applicable)

*/
bool  VEvndispReconstructionCut::test( int iValue, int ntubes )
{
    if( !fCutSet )
    {
        return true;
    }
    
    if( !testNtubes( ntubes ) )
    {
        return true;
    }
    if( iValue > fCut_int_min && iValue < fCut_int_max )
    {
        return true;
    }
    return false;
}

/*
      test this particular cut (double version)

      take also number of ntubes into account

      return true if cut was passed  (or not applicable)

*/
bool  VEvndispReconstructionCut::test( double iValue, int ntubes )
{
    if( !fCutSet )
    {
        return true;
    }
    
    if( !testNtubes( ntubes ) )
    {
        return true;
    }
    if( iValue > fCut_double_min && iValue < fCut_double_max )
    {
        return true;
    }
    return false;
}
