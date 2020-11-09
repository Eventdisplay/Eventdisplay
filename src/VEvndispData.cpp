/*! \class VEvndispData
    \brief central data class

   \attention

     fReader is not static! Call of VEvndispData::initializeDataReader() necessary  in constructors of all inherit classes


*/

#include "VEvndispData.h"

VEvndispData::VEvndispData()
{
    fReader = 0;
}


/*!
 * this function should be called only once in the initialization
 */
void VEvndispData::setTeltoAna( vector< unsigned int > iT )
{
    fTeltoAna = iT;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // fExpectedEventStatus is used in VDisplay -> don't allow to display more than 8*sizeof(size_t) telescope (64?)
    //
    bitset<8 * sizeof( size_t )> ib;
    if( getRunParameter()->fdisplaymode && getNTel() >= ib.size() )
    {
        cout << "Warning: maximum telescopes ID allowed in display mode: " << ib.size() << endl;
    }
    for( unsigned int i = 0; i < iT.size(); i++ )
    {
        if( iT[i] >= ib.size() )
        {
            if( getRunParameter()->fdisplaymode )
            {
                if( getNTel() >= ib.size() )
                {
                    fExpectedEventStatus = 0;
                }
                else
                {
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
            }
            fExpectedEventStatus = 0;
        }
        else
        {
            ib.set( iT[i], 1 );
        }
    }
    if( fExpectedEventStatus != 0 )
    {
        fExpectedEventStatus = ( unsigned long int ) ib.to_ulong();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // init event counting statistics
    vector< int > itemp;
    itemp.assign( fNTel + 1, 0 );
    fTriggeredTel.assign( fNTel, itemp );
    fTriggeredTelN.assign( getTeltoAna().size() + 1, 0 );
    
    // analysis event status
    fAnalysisArrayEventStatus = 0;
    fAnalysisTelescopeEventStatus.assign( fNTel, 0 );
}


void VEvndispData::setTelID( unsigned int iTel )
{
    if( iTel < fNTel )
    {
        fTelID = iTel;
        if( fReader != 0 )
        {
            fReader->setTelescopeID( iTel );
        }
        else
        {
            cout << "VEvndispData::setTelID: fatal error: fReader == 0" << endl;
            exit( EXIT_FAILURE );
        }
        if( fDetectorGeo != 0 )
        {
            fDetectorGeo->setTelID( iTel );
        }
        else
        {
            cout << "VEvndispData::setTelID warning: fDetectorGeo == 0" << endl;
        }
    }
    else
    {
        cout << "VEvndispData::setTelescope: error: invalid telescope number " << iTel << "\t" << fTelID << endl;
    }
}


bool VEvndispData::initializeDataReader()
{
    if( getDebugFlag() )
    {
        cout << "VEvndispData::initializeDataReader()" << endl;
    }
    
    if( fGrIsuReader != 0 )
    {
        fReader = ( VVirtualDataReader* )fGrIsuReader;
    }
#ifndef NOVBF
    else if( fRawDataReader != 0 )
    {
        fReader = ( VVirtualDataReader* )fRawDataReader;
    }
#endif
    else if( fDSTReader != 0 )
    {
        fReader = ( VVirtualDataReader* )fDSTReader;
    }
    else if( fMultipleGrIsuReader != 0 )
    {
        fReader = ( VVirtualDataReader* )fMultipleGrIsuReader;
    }
    else if( fPEReader != 0 )
    {
        fReader = ( VVirtualDataReader* )fPEReader;
    }
    
    if( !fReader )
    {
        return false;
    }
    
    // set trace integration method
    fReader->setTeltoAna( fTeltoAna );
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fReader->setTelescopeID( i );
        if( i < getRunParameter()->fTraceIntegrationMethod.size()
                && getRunParameter()->fTraceIntegrationMethod[i] == 0 )
        {
            fReader->setPerformFADCAnalysis( i, false );
        }
        else
        {
            fReader->setPerformFADCAnalysis( i, true );
        }
        
    }
    
    return true;
}


void VEvndispData::testDataReader()
{
    if( fReader == 0 )
    {
        if( fRunPar->fsourcetype == 1 && fGrIsuReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fGrIsuReader;
        }
#ifndef NOVBF
        else if( fRunPar->fsourcetype == 0 && fRawDataReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fRawDataReader;
        }
#endif
        else if( fRunPar->fsourcetype == 4  && fDSTReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fDSTReader;
        }
        else if( fRunPar->fsourcetype == 5  && fMultipleGrIsuReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fMultipleGrIsuReader;
        }
        else if( fRunPar->fsourcetype == 6 && fPEReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fPEReader;
        }
        else if( fRunPar->fsourcetype == 7  && fDSTReader != 0 )
        {
            fReader = ( VVirtualDataReader* )fDSTReader;
        }
        else
        {
            cout << "VEvndispData::testDataReader() error: no reader found" << endl;
            exit( EXIT_FAILURE );
        }
    }
}


void VEvndispData::resetAnaData()
{
    if( fTelID < fAnaData.size() )
    {
        fAnaData[fTelID]->fSums = 0.;
        fAnaData[fTelID]->fTCorrectedSumFirst = fRunPar->fsumfirst[fTelID];
        fAnaData[fTelID]->fTCorrectedSumLast = fRunPar->fsumfirst[fTelID] + fRunPar->fsumwindow_1[fTelID];
        fAnaData[fTelID]->fCurrentSummationWindow = fRunPar->fsumwindow_1[fTelID];
        fAnaData[fTelID]->fCurrentSummationWindow_2 = fRunPar->fsumwindow_2[fTelID];
        
        fAnaData[fTelID]->fTemplateMu = 0;
    }
}

/*

   define detector geometry (telescope positions, pixel position, pixel numbering, etc)

*/
bool VEvndispData::setDetectorGeometry( unsigned int iNTel, vector< string > iCamera, string iDir )
{
    if( fDebug )
    {
        cout << "VEvndispData::setDetectorGeometry " << iNTel << "\t" << iDir << endl;
        for( unsigned int i = 0; i < iCamera.size(); i++ )
        {
            cout << "\t T" << i + 1 << "\t" << iCamera[i] << endl;
        }
    }
    bool iMakeNeighbourList = false;
    //////////////////////////////////////////////////////////////////////////////////////////
    // read detector geometry from a configuration file and/or DB
    // (all cases but DSTs)
    if( getRunParameter()->fsourcetype != 7 && getRunParameter()->fsourcetype != 4 )
    {
        fDetectorGeo = new VDetectorGeometry( iNTel, iCamera, iDir, fDebug,
                                              getRunParameter()->fCameraCoordinateTransformX, getRunParameter()->fCameraCoordinateTransformY,
                                              getRunParameter()->fsourcetype );
        // get camera rotations from the DB
        if( getRunParameter()->fDBCameraRotationMeasurements )
        {
            fDetectorGeo->readDetectorGeometryFromDB( getRunParameter()->fDBRunStartTimeSQL, getRunParameter()->fDBCameraRotationMeasurements );
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    // for DST files: read detector geometry from DST file
    // (telconfig tree)
    else
    {
        fDetectorGeo = new VDetectorGeometry( iNTel, fDebug );
        fDetectorGeo->setSourceType( getRunParameter()->fsourcetype );
        TFile iDetectorFile( getRunParameter()->fsourcefile.c_str() );
        if( iDetectorFile.IsZombie() )
        {
            cout << "VEvndispData::setDetectorGeometry error opening detector file: " << getRunParameter()->fsourcefile << endl;
            exit( EXIT_FAILURE );
        }
        TTree* iTree = ( TTree* )iDetectorFile.Get( "telconfig" );
        if( !iTree )
        {
            cout << "VEvndispData::setDetectorGeometry error: cannot find detector tree (telconfig) in " << getRunParameter()->fsourcefile << endl;
            exit( EXIT_FAILURE );
        }
        VDetectorTree iDetectorTree;
        iDetectorTree.readDetectorTree( fDetectorGeo, iTree, ( getRunParameter()->getObservatory() != "VERITAS" ) );
        // neighbour list will be set up later, after reading the reconstruction runparameter files
        iMakeNeighbourList = true;
        if( fDebug )
        {
            cout << "VEvndispData::setDetectorGeometry reading detector geometry from DST file" << endl;
        }
    }
    // print most important parameters of the detector
    if( fDetectorGeo )
    {
        fDetectorGeo->print( false );
    }
    return iMakeNeighbourList;
}

TTree* VEvndispData::getDetectorTree()
{
    if( fDetectorTree )
    {
        return fDetectorTree->getTree();
    }
    
    fDetectorTree = new VDetectorTree();
    fDetectorTree->fillDetectorTree( fDetectorGeo );
    
    return fDetectorTree->getTree();
}




void VEvndispData::setDeadChannelText()
{
    if( fDebug )
    {
        cout << "VEvndispData::setDeadChannelText()" << endl;
    }
    
    fDeadChannelText.clear();
    
    fDeadChannelText.push_back( "working" );
    fDeadChannelText.push_back( "dead: out of pedestal range" );
    fDeadChannelText.push_back( "dead: small absolute pedvars" );
    fDeadChannelText.push_back( "dead: small relative pedvars" );
    fDeadChannelText.push_back( "dead: large relative pedvars" );
    fDeadChannelText.push_back( "dead: outside gain range" );
    fDeadChannelText.push_back( "dead: small/large gain vars" );
    fDeadChannelText.push_back( "dead: large gain deviation" );
    fDeadChannelText.push_back( "dead: large time offset" );
    fDeadChannelText.push_back( "disabled: FADC stop signal" );
    fDeadChannelText.push_back( "disabled: masked" );
    fDeadChannelText.push_back( "disabled: user set" );
    fDeadChannelText.push_back( "disabled: MC set" );
    fDeadChannelText.push_back( "dead: L1 rates out of range" );
    fDeadChannelText.push_back( "dead: HV values out of range" );
    fDeadChannelText.push_back( "dead: currents out of range" );   // not used (yet)
}

bool VEvndispData::get_reconstruction_parameters( string ifile, bool iMakeNeighbourList )
{
    fEvndispReconstructionParameter = new VEvndispReconstructionParameter( getDetectorGeometry()->getTelType(), getRunParameter() );
    fEvndispReconstructionParameter->SetName( "EvndispReconstructionParameter" );
    
    unsigned int iNMethods = 0;
    if( ifile.size() > 0 )
    {
        // read array analysis cuts
        if( ifile.find( "/" ) != string::npos )
        {
            iNMethods = fEvndispReconstructionParameter->read_arrayAnalysisCuts( ifile );
        }
        else
        {
            iNMethods = fEvndispReconstructionParameter->read_arrayAnalysisCuts( getRunParameter()->getDirectory_EVNDISPParameterFiles() + "/" + ifile );
        }
        if( !iNMethods )
        {
            return false;
        }
        // make list of neighbours in detector geometry element (scaling factors are only now known )
        if( iMakeNeighbourList )
        {
            cout << "\t making list of neighbours for image cleaning";
            getDetectorGeometry()->makeNeighbourList( getRunParameter()->fNeighbourDistanceFactor,
                    getRunParameter()->fSquarePixels );
            cout << "...done" << endl;
        }
        // print cuts to screen (only done for the relevant run modes)
        if( fRunPar->frunmode != 1 && fRunPar->frunmode != 2 && fRunPar->frunmode != 5 && fRunPar->frunmode != 6 )
        {
            fEvndispReconstructionParameter->print_arrayAnalysisCuts();
        }
        // notify array analysis of the number of reconstruction methods
        if( fShowerParameters )
        {
            fShowerParameters->setNArrayReconstructionMethods( iNMethods );
        }
    }
    else
    {
        return false;
    }
    
    return true;
}


/*!
    dump all tree data for current event to stdout
*/
void VEvndispData::dumpTreeData()
{
    //////////////////////////////////////////////////////////////////////////////
    //
    // THIS PRODUCES A SEGMENTATION FAULT
    //
    //////////////////////////////////////////////////////////////////////////////
    /*    if( fDetectorTree )
        {
            cout << endl << "telescopes" << endl;
    cout << fDetectorTree->getTree()->GetName() << endl;
    fDetectorTree->getTree()->Print();
    
            fDetectorTree->getTree()->Scan();
            cout << endl;
        } */
    //////////////////////////////////////////////////////////////////////////////
    if( fShowerParameters )
    {
        fShowerParameters->getTree()->Show( -1 );
    }
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        setTelID( i );
        cout << "===============================================================" << endl;
        cout << "Telescope " << i + 1 << ":" << endl;
        if( getImageParameters() )
        {
            if( getImageParameters()->getTree() )
            {
                getImageParameters()->getTree()->Show( -1 );
            }
            else
            {
                cout << "\t no image parameter tree found..." << endl;
            }
        }
        else
        {
            cout << "\t no image parameters found..." << endl;
        }
        if( getRunParameter()->fImageLL )
        {
            cout << "Log-likelihood tree: " << endl;
            if( getImageParameters( getRunParameter()->fImageLL )->getTree() )
            {
                getImageParameters( getRunParameter()->fImageLL )->getTree()->Show( -1 );
            }
        }
    }
}


void VEvndispData::setDead( unsigned int iChannel, unsigned int iDead, bool iLowGain, bool iFullSet, bool iReset )
{
    if( iFullSet )
    {
        if( iChannel < getDead( iLowGain ).size() )
        {
            getDead( iLowGain )[iChannel] = iDead;
        }
        bitset<8 * sizeof( uint32_t )> idead = getDead( iLowGain )[iChannel];
        for( unsigned int i = 0; i < idead.size(); i++ )
        {
            if( idead.test( i ) )
            {
                getDeadUI( iLowGain )[iChannel] = i;
                break;
            }
        }
    }
    else
    {
        if( iDead == 0 )
        {
            return;
        }
        
        if( iChannel < getDead( iLowGain ).size() )
        {
            bitset<8 * sizeof( uint32_t )> ibit_dead = getDead( iLowGain )[iChannel];
            if( !iReset )
            {
                ibit_dead.set( iDead, 1 );
            }
            else
            {
                ibit_dead.set( iDead, 0 );
            }
            getDead( iLowGain )[iChannel] = ibit_dead.to_ulong();
        }
        if( iChannel < getDeadUI( iLowGain ).size() )
        {
            getDeadUI( iLowGain )[iChannel] = iDead;
        }
    }
}


void VEvndispData::endOfRunInfo()
{
    cout << endl;
    cout << "Event statistics:" << endl;
    cout << "-----------------" << endl;
    // print this only for a small number of telescopes
    if( fTriggeredTel.size() < 10 )
    {
        cout << "\t Multiplicity:  | ";
        for( unsigned int i = 0; i < fTriggeredTel[0].size(); i++ )
        {
            cout << i << "\t";
        }
        cout << endl;
        cout << "\t ----------------------------------------------" << endl;
        for( unsigned int i = 0; i < fTriggeredTel.size(); i++ )
        {
            cout << "\t Telescope   " << i + 1 << ": | ";
            for( unsigned int j = 0; j < fTriggeredTel[i].size(); j++ )
            {
                cout << fTriggeredTel[i][j] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    for( unsigned int i = 0; i < fTriggeredTelN.size(); i++ )
    {
        cout << "\t number of " << i << "-fold events: " << fTriggeredTelN[i] << endl;
    }
    cout << "-----------------------------------------------" << endl;
}

//if iGrepAble == true, prints DEADCHAN at beginning of line. Prints list of dead pixels, some properties, and the reason why they are dead.
void VEvndispData::printDeadChannels( bool iLowGain, bool iGrepAble )
{
    TString mode = iLowGain ? " low" : "high" ;
    
    unsigned int idead = 0;
    cout << endl;
    
    cout << "Dead channel list (" << mode << " gain) for Telescope " << getTelID() + 1 << endl;
    cout << "==================================" << endl;
    
    for( unsigned int i = 0; i < getDead( iLowGain ).size(); i++ )
    {
        if( getDead( iLowGain )[i] > 0 )
        {
            bitset<16 * sizeof( uint32_t )> i_dead = getDead( iLowGain )[i];
            TString linestart, extrainfo ;
            if( iGrepAble )
            {
                linestart.Form( "DEADCHAN Tel %d, Channel %3d, %s gain, run %d: ", getTelID() + 1, i, mode.Data(), getRunNumber() ) ;
            }
            else
            {
                linestart.Form( "\t%3d: ", i );
            }
            
            if( i_dead.test( 9 ) ) //L2 channel
            {
                cout << linestart  << fDeadChannelText[9] << endl;
            }
            else if( i_dead.test( 11 ) ) //user set
            {
                cout << linestart  << fDeadChannelText[11] << endl;
            }
            //set dead by eventdisplay
            else
            {
                extrainfo.Form( "(pedestal %5.1f, pedvar %5.1f, rel. gain %3.2f, gainvar %3.2f, L1 rate %.2e Hz, HV %4d V, I %4.1f muA)", getPeds( iLowGain )[i] , getPedvars( ( bool )iLowGain, getSumWindow() )[ i ], getGains( iLowGain )[i] , getGainvars( iLowGain )[i] , getL1Rate( i ), ( int )getHV( i ), getCurrent( i ) );
                
                cout << linestart << extrainfo ;
                for( unsigned j = 0; j < i_dead.size(); j++ )
                {
                    if( i_dead.test( j ) && j < fDeadChannelText.size() )
                    {
                        cout << " - " << fDeadChannelText[j] ;
                    }
                }
                cout << endl;
                idead++;
                
            }
            
        }
    }
    cout << "Total number of dead channels (" << mode << " gain) for Telescope " << getTelID() + 1 << ": " << idead << endl;
    cout << "==================================" << endl;
}



VImageParameter* VEvndispData::getImageParameters( int iselect )
{
    if( iselect > 0 )
    {
        return fAnaData[fTelID]->fImageParameterLogL;
    }
    
    return fAnaData[fTelID]->fImageParameter;
}


/*
 * read parameters for dead channel search
 *
 */
bool VEvndispData::initializeDeadChannelFinder()
{
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        bool itelStats = getTelescopeStatus( getTeltoAna()[i] );
        
        fDeadChannelDefinition_HG.push_back( new VDeadChannelFinder( getRunParameter()->frunmode, getTeltoAna()[i], false, fReader->isMC() ) );
        fDeadChannelDefinition_LG.push_back( new VDeadChannelFinder( getRunParameter()->frunmode, getTeltoAna()[i], true, fReader->isMC() ) );
        if( itelStats && getRunParameter()->fsourcetype != 6
                && getRunParameter()->fsourcetype != 4
                && getRunParameter()->fsourcetype != 7 )
        {
            if( getRunParameter()->fDeadChannelFile.size() > 0 )
            {
                fDeadChannelDefinition_HG.back()->readDeadChannelFile( getRunParameter()->getDirectory_EVNDISPParameterFiles() + "/" + getRunParameter()->fDeadChannelFile );
            }
            fDeadChannelDefinition_HG.back()->printDeadChannelDefinition();
            if( getRunParameter()->fDeadChannelFile.size() > 0 )
            {
                fDeadChannelDefinition_LG.back()->readDeadChannelFile( getRunParameter()->getDirectory_EVNDISPParameterFiles() + "/" + getRunParameter()->fDeadChannelFile );
            }
            fDeadChannelDefinition_LG.back()->printDeadChannelDefinition();
        }
        else if( getRunParameter()->fsourcetype != 7 && getRunParameter()->fsourcetype != 4 )
        {
            cout << "Ignoring dead channel settings for this run mode (";
            cout << getRunParameter()->fsourcetype << ") and ";
            cout << " telescope " << getTeltoAna()[i] + 1 << endl;
        }
    }
    return true;
}


bool VEvndispData::getTelescopeStatus( unsigned int iTelID )
{
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getTeltoAna()[i] == iTelID )
        {
            return true;
        }
    }
    return false;
}


valarray<double>& VEvndispData::getPeds( bool iLowGain, double iTime )
{
    if( !getRunParameter()->fPlotRaw )
    {
        if( iTime < 90. )
        {
            return fCalData[fTelID]->getPeds( iLowGain, getEventTime() );
        }
        else
        {
            return fCalData[fTelID]->getPeds( iLowGain, iTime );
        }
    }
    fPlotRawPedestals.resize( getNChannels(), getDetectorGeometry()->getDefaultPedestal() );
    
    return fPlotRawPedestals;
}

valarray<double>& VEvndispData::getPedvars( bool iLowGain, unsigned int iSW, double iTime )
{
    if( iSW == 0 )
    {
        iSW = getSumWindow();
    }
    if( iTime < -90. )
    {
        return fCalData[fTelID]->getPedvars( iLowGain, iSW, getEventTime() );
    }
    
    return fCalData[fTelID]->getPedvars( iLowGain, iSW, iTime );
}


void VEvndispData::setDebugLevel( int i )
{
    fDebugLevel = i;
    if( i > 0 )
    {
        fNDebugMessages++;
    }
}


int VEvndispData::getDebugLevel()
{
    int nDebugMessagesMax = 100;
    
    // allow a maximum of 100 debug messages
    if( fNDebugMessages > nDebugMessagesMax )
    {
        return 2;
    }
    else if( fNDebugMessages == nDebugMessagesMax )
    {
        cout << "VEvndispData:: - more than " << nDebugMessagesMax << " error messages or warnings" << endl;
        cout << " ---------- stopped printout of error messages ---------------------------" << endl;
        fNDebugMessages++;
        return 2;
    }
    
    return fDebugLevel;
}


void VEvndispData::setCurrentSummationWindow( unsigned int imin, unsigned int imax, bool iSecondWindow )
{
    unsigned int iT = 0;
    if( imax > imin )
    {
        iT = imax - imin;
    }
    
    // this should never happen
    if( iT > getNSamples() )
    {
        cout << "VEvndispData::setCurrentSummationWindow (a) error: summation window too large ";
        if( iSecondWindow )
        {
            cout << " (2nd window): ";
        }
        else
        {
            cout << ": ";
        }
        cout << imin << "\t" << imax << "\t" << iT << endl;
        iT = 0;
    }
    if( iSecondWindow )
    {
        fAnaData[fTelID]->fCurrentSummationWindow = iT;
    }
    else
    {
        fAnaData[fTelID]->fCurrentSummationWindow_2 = iT;
    }
    
}


void VEvndispData::setCurrentSummationWindow( unsigned int iChannel, unsigned int imin, unsigned int imax, bool iSecondWindow )
{
    // first summation window
    if( !iSecondWindow )
    {
        if( imax > imin )
        {
            fAnaData[fTelID]->fCurrentSummationWindow[iChannel] = imax - imin;
        }
        else
        {
            fAnaData[fTelID]->fCurrentSummationWindow[iChannel] = 0;
        }
        
        // this should never happen
        if( fAnaData[fTelID]->fCurrentSummationWindow[iChannel] > getNSamples() )
        {
            cout << "VEvndispData::setCurrentSummationWindow (b) error: summation window too large: ";
            cout << fTelID << "\t" << iChannel << "\t" << imin << "\t" << imax << "\t" << fAnaData[fTelID]->fCurrentSummationWindow[iChannel] << endl;
            fAnaData[fTelID]->fCurrentSummationWindow[iChannel] = 0;
        }
    }
    // second summation window
    else
    {
        if( imax > imin )
        {
            fAnaData[fTelID]->fCurrentSummationWindow_2[iChannel] = imax - imin;
        }
        else
        {
            fAnaData[fTelID]->fCurrentSummationWindow_2[iChannel] = 0;
        }
        
        // this should never happen
        if( fAnaData[fTelID]->fCurrentSummationWindow_2[iChannel] > getNSamples() )
        {
            cout << "VEvndispData::setCurrentSummationWindow (2nd window) (b) error: summation window too large: ";
            cout << fTelID << "\t" << iChannel << "\t" << imin << "\t" << imax << "\t" << fAnaData[fTelID]->fCurrentSummationWindow_2[iChannel] << endl;
            fAnaData[fTelID]->fCurrentSummationWindow_2[iChannel] = 0;
        }
    }
}


unsigned int VEvndispData::getTeltoAnaID( unsigned int iTelID )
{
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getTeltoAna()[i] == iTelID )
        {
            return i;
        }
    }
    return 0;
}


VDeadChannelFinder* VEvndispData::getDeadChannelFinder( bool iLowGain )
{
    if( !iLowGain )
    {
        return fDeadChannelDefinition_HG[getTeltoAnaID()];
    }
    
    return fDeadChannelDefinition_LG[getTeltoAnaID()];
}

bool VEvndispData::isTeltoAna( unsigned int iTel )
{
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getTeltoAna()[i] == iTel )
        {
            return true;
        }
    }
    
    return false;
}

/*
 * return typical (arrival) time of FADC pulse
 *
 * a) T0
 * b) average pulse time
 */
valarray<double>& VEvndispData::getPulseTime( bool iCorrected )
{
     // in run parameter file: FADCSUMMATIONSTART set to TZERO
    if( getSumWindowStart_T_method() == 1 )
    {
        return fAnaData[fTelID]->getTZeros( iCorrected );
    }
     // in run parameter file: FADCSUMMATIONSTART set to TTRIGGER
     else if( getSumWindowStart_T_method() == 3 )
     {
         return fAnaData[fTelID]->getTTrigger();
     }
    
    // default: return average pulse time
     // in run parameter file: FADCSUMMATIONSTART set to TAVERAGE
    return fAnaData[fTelID]->getTraceAverageTime( iCorrected );
}


/*
 *  return timing values characterizing a pulse
 *  time of max, tzero, etc. (0.5, 1., 0.5 values)
 *
 */
vector< valarray<double> >& VEvndispData::getPulseTiming( bool iCorrected )
{
    if( iCorrected )
    {
        return fAnaData[fTelID]->fPulseTimingCorrected;
    }
    
    return fAnaData[fTelID]->fPulseTimingUncorrected;
}

void VEvndispData::setPulseTiming( vector< valarray< double > > iPulseTiming, bool iCorrected )
{
    if( iCorrected )
    {
        fAnaData[fTelID]->fPulseTimingCorrected   = iPulseTiming;
    }
    else
    {
        fAnaData[fTelID]->fPulseTimingUncorrected = iPulseTiming;
    }
}

/*!
    reset pulse timing
*/
void VEvndispData::setPulseTiming( float iP, bool iCorrected )
{
    for( unsigned int i = 0; i < getPulseTiming( iCorrected ).size(); i++ )
    {
        getPulseTiming( iCorrected )[i] = iP;
    }
}

void VEvndispData::setPulseTiming( unsigned int iChannel, vector< float > iP, bool iCorrected )
{
    // check pulse timing vector sizes
    if( iP.size() != getPulseTiming( iCorrected ).size() )
    {
        cout << "VEvndispData::setPulseTiming( unsigned int iChannel, vector< float > iP ) vector size problem: ";
        cout << iP.size() << "\t" << getPulseTiming( iCorrected ).size() << endl;
        return;
    }
    // fill all values
    for( unsigned int i = 0; i < getPulseTiming( iCorrected ).size(); i++ )
    {
        if( iChannel < getPulseTiming( iCorrected )[i].size() )
        {
            getPulseTiming( iCorrected )[i][iChannel] = iP[i];
        }
    }
}

void VEvndispData::setPulseTimingCorrection( unsigned int iChannel, double iCorrection )
{
    // pulse timing parameters <= fpulsetiming_tzero_index are corrected
    // (all times later are pulse widths)
    if( getRunParameter()->fpulsetiming_tzero_index < getPulseTiming( true ).size() )
    {
        for( unsigned int i = 0; i <= getRunParameter()->fpulsetiming_tzero_index; i++ )
        {
            if( iChannel < getPulseTiming( true )[i].size() )
            {
                getPulseTiming( true )[i][iChannel] += iCorrection;
            }
        }
    }
    // average trace time corrections
    if( iChannel < getTraceAverageTime( false ).size() && iChannel < getTraceAverageTime( false ).size() )
    {
        getTraceAverageTime( true )[iChannel] = getTraceAverageTime( false )[iChannel] + iCorrection;
    }
}

/*
 * set all pulse timing values to the same vale
 */
void VEvndispData::setPulseTiming( unsigned int iChannel, bool iCorrection, double iValue )
{
    if( getRunParameter()->fpulsetiming_tzero_index < getPulseTiming( iCorrection ).size() )
    {
        for( unsigned int i = 0; i <= getRunParameter()->fpulsetiming_tzero_index; i++ )
        {
            if( iChannel < getPulseTiming( iCorrection )[i].size() )
            {
                getPulseTiming( iCorrection )[i][iChannel] = iValue;
            }
        }
    }
}

unsigned int VEvndispData::getLargestSumWindow()
{
    unsigned int iSW = 0;
    
    if( getSumWindow() > iSW )
    {
        iSW = getSumWindow();
    }
    if( getSumWindow_2() > iSW )
    {
        iSW = getSumWindow_2();
    }
    if( getSumWindow_Pass1() > iSW )
    {
        iSW = getSumWindow_Pass1();
    }
    
    return iSW;
}

unsigned int VEvndispData::getLargestSumWindow( unsigned int iTelID )
{
    unsigned int iSW = 0;
    
    if( getSumWindow() > iSW )
    {
        iSW = getSumWindow( iTelID );
    }
    if( getSumWindow_2() > iSW )
    {
        iSW = getSumWindow_2( iTelID );
    }
    if( getSumWindow_Pass1() > iSW )
    {
        iSW = getSumWindow_Pass1( iTelID );
    }
    
    return iSW;
}

/*
   make sure that summation window is not larger than readout window

*/
unsigned int VEvndispData::checkSummationWindow( unsigned int iTelID, unsigned int iSumWindow )
{
    if( iSumWindow <= getNSamples( iTelID ) || getNSamples( iTelID ) == 0 )
    {
        return iSumWindow;
    }
    
    return getNSamples( iTelID );
}

bool VEvndispData::setSums( valarray< double > iVSum )
{
    if( iVSum.size() != getNChannels() )
    {
        cout << "VEvndispData::setSums() error: setting wrong vector size for integrated charges: " << getNChannels() << "\t" << iVSum.size() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fAnaData[fTelID]->fSums = iVSum;
    return true;
}

bool VEvndispData::setAverageTZero( unsigned int iChannel, double iTZero, bool iLowGain )
{
    if( !iLowGain )
    {
        if( iChannel < fCalData[fTelID]->fAverageTzero.size() )
        {
            fCalData[fTelID]->fAverageTzero[iChannel] = iTZero;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if( iChannel < fCalData[fTelID]->fLowGainAverageTzero.size() )
        {
            fCalData[fTelID]->fLowGainAverageTzero[iChannel] = iTZero;
        }
        else
        {
            return false;
        }
    }
    return true;
}

bool VEvndispData::setAverageTZerovars( unsigned int iChannel, double iTZero, bool iLowGain )
{
    if( !iLowGain )
    {
        if( iChannel < fCalData[fTelID]->fAverageTzerovars.size() )
        {
            fCalData[fTelID]->fAverageTzerovars[iChannel] = iTZero;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if( iChannel < fCalData[fTelID]->fLowGainAverageTzerovars.size() )
        {
            fCalData[fTelID]->fLowGainAverageTzerovars[iChannel] = iTZero;
        }
        else
        {
            return false;
        }
    }
    return true;
}

ULong64_t VEvndispData::getTelType( unsigned int iTelID )
{
    if( getDetectorGeometry() )
    {
        if( iTelID < getDetectorGeometry()->getTelType().size() )
        {
            return getDetectorGeometry()->getTelType()[iTelID];
        }
    }
    return 99999;
}

unsigned int VEvndispData::getTelType_Counter( ULong64_t iTelType )
{
    if( getDetectorGeometry() )
    {
        return getDetectorGeometry()->getTelType_Counter( iTelType );
    }
    return 99999;
}


/*

    called for first events

    assume that MJD is not changing significantly
*/
bool VEvndispData::initializeStarCatalogue( int iMJD, double iTime )
{
    if( getRunParameter()->fStarCatalogueName.size() > 0 )
    {
        double i_MJD = ( double )iMJD + iTime / 86400.;
        fStarCatalogue = new VStarCatalogue();
        if( !fStarCatalogue->init( i_MJD, getRunParameter()->fStarCatalogueName ) )
        {
            return false;
        }
        double iMaxFOV = 0.;
        for( unsigned int i = 0; i < getDetectorGeometry()->getFieldofView().size(); i++ )
        {
            if( getDetectorGeometry()->getFieldofView()[i] > iMaxFOV )
            {
                iMaxFOV = getDetectorGeometry()->getFieldofView()[i];
            }
        }
        iMaxFOV *= 1.5;
        fStarCatalogue->setFOV( getArrayPointing()->getTelRA()*TMath::RadToDeg(), getArrayPointing()->getTelDec()*TMath::RadToDeg(),
                                iMaxFOV / 2., iMaxFOV / 2., false, getRunParameter()->fMinStarBrightness_B, "B" );
    }
    return true;
}

/*
 *  check if a channel is dead
 *
 *  for low gain channel: return != 0 if low gain is ok, but high gain is dead
 */
unsigned int VEvndispData::getDead( unsigned int iChannel, bool iLowGain = false )
{
    // high gain channels:
    if( !iLowGain )
    {
        return getDead( iLowGain )[iChannel];
    }
    
    // for low gain channels only:
    
    // low gain is dead
    if( getDead( iLowGain )[iChannel] )
    {
        return getDead( iLowGain )[iChannel];
    }
    // low gain is ok, but high gain is dead
    return getDead( false )[iChannel];
}

double VEvndispData::getAverageElevation()
{
    double iAverageElevation = 0.;
    // for MC: return average elevation from CORSIKA
    if( isMC() )
    {
        if( getReader() && getReader()->getMonteCarloHeader() )
        {
            iAverageElevation = 0.5 * ( getReader()->getMonteCarloHeader()->alt_range[0]
                                        + getReader()->getMonteCarloHeader()->alt_range[1] );
            iAverageElevation *= TMath::RadToDeg();
        }
    }
    // for data: loop over run and calculate average elevation
    else
    {
        double iN = 0.;
        double ze  = 0.;
        double az = 0.;
        float i_start = getRunParameter()->fDBDataStartTimeSecOfDay;
        if( getRunParameter()->fTimeCutsMin_min > 0 )
        {
            i_start += ( float )( getRunParameter()->fTimeCutsMin_min ) * 60.;
        }
        float i_stopp = getRunParameter()->fDBDataStoppTimeSecOfDay;
        if( getRunParameter()->fTimeCutsMin_max > 0 && getRunParameter()->fDBDataStartTimeSecOfDay + getRunParameter()->fTimeCutsMin_max * 60 < getRunParameter()->fDBDataStoppTimeSecOfDay )
        {
            i_stopp = getRunParameter()->fDBDataStartTimeSecOfDay + ( float )( getRunParameter()->fTimeCutsMin_max ) * 60.;
        }
        
        for( float i = i_start; i < i_stopp; i++ )
        {
            VSkyCoordinatesUtilities::getHorizontalCoordinates( getRunParameter()->fDBDataStartTimeMJD, i,
                    getRunParameter()->fTargetDec, getRunParameter()->fTargetRA,
                    az, ze );
            iAverageElevation += 90. - ze;
            iN++;
        }
        if( iN > 0. )
        {
            iAverageElevation /= iN;
        }
        else
        {
            iAverageElevation = 0.;
        }
    }
    
    return iAverageElevation;
}

VImageCleaningRunParameter* VEvndispData::getImageCleaningParameter( bool iDoublePassParameters, unsigned int iTelID )
{
    if( iTelID == 99999 )
    {
        iTelID = fTelID;
    }
    if( iDoublePassParameters )
    {
        if( iTelID < getRunParameter()->fImageCleaningParameters_DB_Pass1.size() )
        {
            return getRunParameter()->fImageCleaningParameters_DB_Pass1[iTelID];
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if( iTelID < getRunParameter()->fImageCleaningParameters.size() )
        {
            return getRunParameter()->fImageCleaningParameters[iTelID];
        }
        else
        {
            return 0;
        }
    }
    return 0;
}

/*
 * check if summation window 1 and 2 are of equal length
 *
 */
bool VEvndispData::isEqualSummationWindows()
{
    if( getSumWindow() == getSumWindow_2() )
    {
        return true;
    }
    return false;
}

////////////////////////////////
// initialize static variables

bool VEvndispData::fDebug = false;
int VEvndispData::fDebugLevel = 0;
int VEvndispData::fNDebugMessages = 0;

// run data
int VEvndispData::fRunNumber = 0;
VEvndispRunParameter* VEvndispData::fRunPar = 0;

// telescope data
unsigned int VEvndispData::fNTel = 1;
unsigned int VEvndispData::fTelID = 0;
vector< unsigned int > VEvndispData::fTeltoAna;
VDetectorGeometry* VEvndispData::fDetectorGeo = 0;
VDetectorTree* VEvndispData::fDetectorTree = 0;

// pointing
VArrayPointing* VEvndispData::fArrayPointing = 0;
vector<VPointing*> VEvndispData::fPointing;
bool VEvndispData::fNoTelescopePointing = false;

// reader
VGrIsuReader* VEvndispData::fGrIsuReader = 0;
VMultipleGrIsuReader* VEvndispData::fMultipleGrIsuReader = 0;
#ifndef NOVBF
VBaseRawDataReader* VEvndispData::fRawDataReader = 0;
#endif
VDSTReader* VEvndispData::fDSTReader = 0;
VPEReader* VEvndispData::fPEReader = 0;

// event data
unsigned int VEvndispData::fEventNumber = 0;
vector< unsigned int > VEvndispData::fTelescopeEventNumber;
unsigned int VEvndispData::fEventType = 0;
unsigned int VEvndispData::fNumberofGoodEvents = 0;
unsigned int VEvndispData::fNumberofIncompleteEvents = 0;
unsigned long int VEvndispData::fExpectedEventStatus = 99;
unsigned int VEvndispData::fAnalysisArrayEventStatus = 0;
vector< unsigned int > VEvndispData::fAnalysisTelescopeEventStatus;

vector< vector< int > > VEvndispData::fTriggeredTel;
vector< int > VEvndispData::fTriggeredTelN;
int VEvndispData::fArrayEventMJD = 0;
int VEvndispData::fArrayPreviousEventMJD = 0;
double VEvndispData::fArrayEventTime = 0.;
vector< int > VEvndispData::fEventMJD;
vector< double > VEvndispData::fEventTime;

// trace handler
VTraceHandler* VEvndispData::fTraceHandler = 0;

//calibration data
vector< bool > VEvndispData::fCalibrated;
vector< VCalibrationData* > VEvndispData::fCalData;

// DB pixel data
VDB_PixelDataReader* VEvndispData::fDB_PixelDataReader = 0;

// dead channel finders
vector< VDeadChannelFinder* > VEvndispData::fDeadChannelDefinition_HG;
vector< VDeadChannelFinder* > VEvndispData::fDeadChannelDefinition_LG;

// analyzer data
TFile* VEvndispData::fOutputfile = 0;
vector< TDirectory* > VEvndispData::fAnaDir;
vector< VImageAnalyzerData* > VEvndispData::fAnaData;
VShowerParameters* VEvndispData::fShowerParameters = 0;
VMCParameters* VEvndispData::fMCParameters = 0;
VEvndispReconstructionParameter* VEvndispData::fEvndispReconstructionParameter = 0;
VFrogsParameters* VEvndispData::fFrogsParameters = 0;
//vector< VFrogImageData* > VEvndispData::fFrogData;

// timing graphs
vector< TGraphErrors* > VEvndispData::fXGraph;
vector< TGraphErrors* > VEvndispData::fXGraphDP2;

// dead channel text
vector< string > VEvndispData::fDeadChannelText;

//  default pedestals for plotraw option
valarray<double> VEvndispData::fPlotRawPedestals;
valarray<double> VEvndispData::fTempValArray;

// star catalogue
VStarCatalogue* VEvndispData::fStarCatalogue = 0;

// dummy vectors
vector< float > VEvndispData::fDummyVector_float;
