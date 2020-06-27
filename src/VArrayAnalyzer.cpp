/*! \class VArrayAnalyzer

     VArrayAnalyzer class for shower direction and core reconstruction

     different analysis methods are provided

     this class is also taking care of filling the showerpars tree correctly

*/

#include "VArrayAnalyzer.h"

VArrayAnalyzer::VArrayAnalyzer()
{
    fDebug = getDebugFlag();
    if( fDebug )
    {
        cout << "VArrayAnalyzer::VArrayAnalyzer()" << endl;
    }
    
    fInitialized = false;
    
    // set up data storage class (all analysis results are stored here)
    fShowerParameters = new VShowerParameters( getNTel(), getRunParameter()->fShortTree,
            getEvndispReconstructionParameter()->getNReconstructionCuts() );
    // set up MC data storage class
    fMCParameters = new VMCParameters( fDebug );
    // test if number of telescopes exceeds value in fShowerParameters
    if( getNTel() > fShowerParameters->getMaxNTelescopes() )
    {
        cout << "VArrayAnalyzer::VArrayAnalyzer(): error, to many telescopes ";
        cout << getNTel() << "\t" << fShowerParameters->getMaxNTelescopes() << endl;
        cout << "\t adjust VDST_MAXTELESCOPES in VGlobalRunParameter" << endl;
        exit( EXIT_FAILURE );
    }
}


VArrayAnalyzer::~VArrayAnalyzer()
{

}


/*!

     stereo analysis

     called for every event (for standard analysis)
*/
void VArrayAnalyzer::doAnalysis()
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::doAnalysis()" << endl;
    }
    
    // only at first call in the analysis run: initialize data class, set trees
    if( !fInitialized )
    {
        initAnalysis();
        fInitialized = true;
    }
    
    // fill simulation data
    if( isMC() )
    {
        fillSimulationEvent( fReader->hasArrayTrigger() );
    }
    
    // return if triggered only events are written to output
    if( getRunParameter()->fWriteTriggerOnly && !fReader->hasArrayTrigger() )
    {
        return;
    }
    
    // reset data vectors, etc.
    initEvent();
    
    if( isMC() )
    {
        // fill shower parameter tree
        if( getShowerParameters() )
        {
            getShowerParameters()->setMC();
            getShowerParameters()->MCprimary = fReader->getMC_primary();
            getShowerParameters()->MCenergy  = fReader->getMC_energy();
            getShowerParameters()->MCxcore   = fReader->getMC_X();
            getShowerParameters()->MCycore   = fReader->getMC_Y();
            getShowerParameters()->MCxcos    = fReader->getMC_Xcos();
            getShowerParameters()->MCycos    = fReader->getMC_Ycos();
            getShowerParameters()->MCTel_Xoff = fReader->getMC_Xoffset();
            getShowerParameters()->MCTel_Yoff = fReader->getMC_Yoffset();
            getShowerParameters()->MCze = fReader->getMC_Ze();
            getShowerParameters()->MCaz = fReader->getMC_Az();
            getShowerParameters()->MCCorsikaRunID = fReader->getMCCorsikaRunID();
            getShowerParameters()->MCCorsikaShowerID = fReader->getMCCorsikaShowerID();
            getShowerParameters()->MCFirstInteractionHeight = fReader->getMCFirstInteractionHeight();
            getShowerParameters()->MCFirstInteractionDepth = fReader->getMCFirstInteractionDepth();
            
        }
    }
    
    // get reference pointing
    if( getArrayPointing() )
    {
        getShowerParameters()->fArrayPointing_Elevation = getArrayPointing()->getTelElevation();
        getShowerParameters()->fArrayPointing_Azimuth   = getArrayPointing()->getTelAzimuth();
        if( isMC() )
        {
            getShowerParameters()->fArrayPointing_deRotationAngle_deg = 0.;
        }
        else
        {
            getShowerParameters()->fArrayPointing_deRotationAngle_deg = 
                         getArrayPointing()->getDerotationAngle( getShowerParameters()->MJD, getShowerParameters()->time )
                         *  TMath::RadToDeg();

        }
        getShowerParameters()->fWobbleNorth             = getArrayPointing()->getWobbleNorth();
        getShowerParameters()->fWobbleEast              = getArrayPointing()->getWobbleEast();
    }
    else
    {
        getShowerParameters()->fArrayPointing_Elevation = 0.;
        getShowerParameters()->fArrayPointing_Azimuth   = 0.;
        getShowerParameters()->fArrayPointing_deRotationAngle_deg = 0.;
        getShowerParameters()->fWobbleNorth = 0.;
        getShowerParameters()->fWobbleEast = 0.;
    }
    
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        // get pointing direction
        if( i < getPointing().size() && getPointing()[i] )
        {
            getShowerParameters()->fTelElevation[i] = getPointing()[i]->getTelElevation();
            getShowerParameters()->fTelAzimuth[i]   = getPointing()[i]->getTelAzimuth();
            // dec/ra of current epoch
            getShowerParameters()->fTelDec[i]       = getPointing()[i]->getTelDec();
            getShowerParameters()->fTelRA[i]        = getPointing()[i]->getTelRA();
        }
        else
        {
            getShowerParameters()->fTelElevation[i] = 0.;
            getShowerParameters()->fTelAzimuth[i]   = 0.;
            getShowerParameters()->fTelDec[i]       = 0.;
            getShowerParameters()->fTelRA[i]        = 0.;
        }
        // compare calculated pointing with vbf pointing
        checkPointing();
    }
    
    ////////////////////////////////////////////////////////
    // do array analysis only if there is an array trigger
    if( getReader()->hasArrayTrigger() )
    {
        // transform telescope positions into shower coordinates (direction is telescope orientation)
        for( unsigned int i = 0; i < getNTel(); i++ )
        {
            // (use only telescopes with valid pointings, otherwise use array pointing)
            if( i < getPointing().size() && getPointing()[i] && getPointing()[i]->getTelElevation() > 0. )
            {
                transformTelescopePosition( i, 90. - getPointing()[i]->getTelElevation(),
                                            getPointing()[i]->getTelAzimuth(),
                                            fReader->isMC() && ( i == getTeltoAna()[0] ) );
            }
            else
            {
                transformTelescopePosition( i, 90. - getShowerParameters()->fArrayPointing_Elevation,
                                            getShowerParameters()->fArrayPointing_Azimuth,
                                            fReader->isMC() && ( i == getTeltoAna()[0] ) );
            }
        }
        //////////////////////////////////////////////////////////////////////////////////////////////
        // calculate shower direction and shower core position
        calcShowerDirection_and_Core();
        //////////////////////////////////////////////////////////////////////////////////////////////
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // calculate RA and dec of shower direction
    // (these are current epoch values)
    for( unsigned int i = 0; i < getShowerParameters()->fNMethods; i++ )
    {
        // test for successfull reconstruction (zenith angle > 0.)
        if( getShowerParameters()->fShowerZe[i] < 0. )
        {
            continue;
        }
        if( getArrayPointing() )
        {
            getShowerParameters()->fDec[i] = getArrayPointing()->getTelDec() + getShowerParameters()->fShower_YoffsetDeRot[i];
            if( cos( getShowerParameters()->fDec[i] / TMath::RadToDeg() ) != 0. )
            {
                getShowerParameters()->fRA[i] = getArrayPointing()->getTelRA()
                                                + getShowerParameters()->fShower_XoffsetDeRot[i] / cos( getShowerParameters()->fDec[i] / TMath::RadToDeg() );
            }
            else
            {
                getShowerParameters()->fRA[i] = getArrayPointing()->getTelRA();
            }
        }
        else
        {
            getShowerParameters()->fDec[i] = 0.;
            getShowerParameters()->fRA[i] = 0.;
        }
    }
    getShowerParameters()->eventStatus = getAnalysisArrayEventStatus();
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // fill shower parameter tree with results
    getShowerParameters()->getTree()->Fill();
    
}


/*!
     this routine is called for each event before the event analysis
*/
void VArrayAnalyzer::initEvent()
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::initEvent()" << endl;
    }
    
    setAnalysisArrayEventStatus( 0 );
    
    // reset shower parameters
    getShowerParameters()->reset( getNTel() );
    
    // get basic infos for this event
    getShowerParameters()->runNumber = getRunNumber();
    getShowerParameters()->eventNumber = getEventNumber();
    getShowerParameters()->fNTelescopes = getNTel();
    
    getShowerParameters()->time = getEventTime();
    getShowerParameters()->MJD  = getEventMJD();
    
    calcTelescopePointing() ;
    getArrayPointing()->fillPointingTree( isMC() );
    
    // set telescope pointing for data (again) and fill pointing tree
    if( !fReader->isMC() )
    {
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            if( getTeltoAna()[i] < getPointing().size() && getPointing()[getTeltoAna()[i]] )
            {
                getPointing()[getTeltoAna()[i]]->setTelPointing( getShowerParameters()->MJD, getShowerParameters()->time,
                        getRunParameter()->fDBTracking, true );
            }
        }
    }
    /////////////////////////////////////////////////////////////////////
    
    // get number of telescopes with local trigger
    // (works only for relatively small number of telescopes)
    bitset<8 * sizeof( unsigned long )> i_ntrigger;
    if( fNTel < i_ntrigger.size() )
    {
        if( fReader->getLocalTrigger().size() >= fNTel )
        {
            for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
            {
                if( fReader->getLocalTrigger()[getTeltoAna()[i]] )
                {
                    if( getTeltoAna()[i] < i_ntrigger.size() )
                    {
                        i_ntrigger.set( getTeltoAna()[i], 1 );
                    }
                }
            }
        }
        getShowerParameters()->fLTrig = i_ntrigger.to_ulong();
    }
    else
    {
        getShowerParameters()->fLTrig = 0;
    }
    // list of telescopes with local trigger
    // (note: ntrig ignores different trigger types)
    if( fReader->getLocalTrigger().size() >= fNTel )
    {
        unsigned int ztrig = 0;
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            if( fReader->getLocalTrigger()[getTeltoAna()[i]] )
            {
                getShowerParameters()->fTrig_list[ztrig] = getTeltoAna()[i];
                getShowerParameters()->fTrig_type[ztrig] = fReader->getLocalTriggerType( i );
                ztrig++;
            }
        }
        getShowerParameters()->fNTrig = ztrig;
    }
}


// do all the calculations just before we write data to the
// pointingData TTree, so we can do it per-event or per-time
void VArrayAnalyzer::calcTelescopePointing()
{

    //////////////////////////////////////////////////////////////////////////////////////////////
    // get the telescope pointing
    //////////////////////////////////////////////////////////////////////////////////////////////
    // MC readers
    // source file is Monte Carlo
    if( fReader->isMC() )
    {
        getArrayPointing()->setMC();
        // calculate mean el and az from all telescopes (should probably be the same for all)
        double i_el = 0.;
        double i_az = 0.;
        double i_n = 0.;
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            // require non-zero elevation
            if( getTeltoAna()[i] < fReader->getTelElevation().size()
                    && getTeltoAna()[i] < fReader->getTelAzimuth().size()
                    && fReader->getTelElevation()[getTeltoAna()[i]] > 0. )
            {
                i_el += fReader->getTelElevation()[getTeltoAna()[i]];
                i_az += fReader->getTelAzimuth()[getTeltoAna()[i]];
                i_n++;
            }
        }
        if( i_n > 1.e-2 )
        {
            getArrayPointing()->setTelElevation( i_el / i_n );
            getArrayPointing()->setTelAzimuth( i_az / i_n );
        }
        else
        {
            getArrayPointing()->setTelElevation( 0. );
            getArrayPointing()->setTelAzimuth( 0. );
        }
    }
    // set pointing direction from command line
    // (very unusual use case)
    else if( getRunParameter()->felevation > 0. && getRunParameter()->fazimuth > 0. )
    {
        getArrayPointing()->setTelElevation( getRunParameter()->felevation );
        getArrayPointing()->setTelAzimuth( getRunParameter()->fazimuth );
    }
    // set pointing direction with target
    // target is set via command line, getPointing() is initiated in run/VEventLoop::initEventLoop()
    else
    {
        if( getArrayPointing()->isSet() && getArrayPointing()->getTargetName() != "laser" )
        {
            getArrayPointing()->updatePointing( getImageParameters()->MJD, getImageParameters()->time );
            if( !getArrayPointing()->isPrecessed() )
            {
                getArrayPointing()->precessTarget( getImageParameters()->MJD, -1 );
                // set wobble offsets
                getArrayPointing()->setWobbleOffset( getRunParameter()->fWobbleNorth, getRunParameter()->fWobbleEast, -1, getImageParameters()->MJD );
                getArrayPointing()->updatePointing( getImageParameters()->MJD, getImageParameters()->time );
            }
        }
    }
}

// fill reduced pointing data (1 row per second,
// instead of 1 per line like VArrayPointing::fillPointingTree() )
void VArrayAnalyzer::generateReducedPointingTreeData()
{
    double MJDStart = 0 ; // start MJD of run
    double MJDStopp = 0 ; // end   MJD of run
    double TimeStart = 0.0 ; // start Time of run
    double TimeStopp = 0.0 ; // end   Time of run
    int i_stat_1 = VSkyCoordinatesUtilities::getMJD_from_SQLstring( getRunParameter()->fDBRunStartTimeSQL, MJDStart, TimeStart );
    int i_stat_2 = VSkyCoordinatesUtilities::getMJD_from_SQLstring( getRunParameter()->fDBRunStoppTimeSQL, MJDStopp, TimeStopp );
    if( i_stat_1 != 0 || i_stat_2 != 0 )
    {
        return;
    }
    
    int    iMJD  = TMath::Nint( MJDStart );
    double iTime = TimeStart ;
    while( iMJD <= TMath::Nint( MJDStopp )  &&  iTime < TimeStopp )
    {
        // loop over run duration, 1 loop per second
        iTime += 1.0 ;
        
        // cycle to next MJD if needed
        if( iTime >= 86400.0 )
        {
            iTime = iTime - 86400.0 ;
            iMJD += 1 ;
        }
        
        // set time, init telescope pointing
        setAnalysisArrayEventStatus( 0 );
        getShowerParameters()->reset( getNTel() );
        getShowerParameters()->runNumber    = getRunNumber();
        getShowerParameters()->eventNumber  = 103 ;
        getShowerParameters()->fNTelescopes = 4 ;
        getShowerParameters()->time = iTime ;
        getShowerParameters()->MJD  = iMJD  ;
        //calcTelescopePointing() ;
        updatePointingToArbitraryTime( iMJD, iTime ) ;
        
        // fill pointing to tree in VArrayPointing
        getArrayPointing()->fillPntReduced() ;
    }
}

/*!
    this routine is called once at the beginning of the analysis
*/
void VArrayAnalyzer::initAnalysis()
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::initAnalysis()" << endl;
    }
    
    fMeanPointingMismatch.assign( getNTel(), 0. );
    fNMeanPointingMismatch.assign( getNTel(), 0. );
    
    // initialize star catalogue
    initializeStarCatalogue( getEventMJD(), getEventTime() );
    
    // loop over all methods and read in necessary MLPs, TMVAs, tables, etc....
    for( unsigned int i = 0; i < getEvndispReconstructionParameter()->getNReconstructionCuts(); i++ )
    {
        // set reconstruction method
        if( getEvndispReconstructionParameter( i ) )
        {
            if( getEvndispReconstructionParameter( i )->fMethodID == 5 )
            {
                fDispAnalyzer.push_back( new VDispAnalyzer() );
                fDispAnalyzer.back()->setTelescopeTypeList( getDetectorGeometry()->getTelType_list() );
                if( !fDispAnalyzer.back()->initialize( getEvndispReconstructionParameter( i )->fDISP_MLPFileName, "MLP" ) )
                {
                    exit( EXIT_FAILURE );
                }
            }
            else if( getEvndispReconstructionParameter( i )->fMethodID == 7 )
            {
                initializeDispAnalyzer( i );
            }
            else
            {
                fDispAnalyzer.push_back( 0 );
            }
        }
    }
}


/*!
   open root outputfile (to be called once before starting the array analysis)
*/
void VArrayAnalyzer::initOutput()
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::initOuput()" << endl;
    }
    // check if root outputfile exist
    if( fOutputfile != 0 )
    {
        return;
    }
    // otherwise create it
    if( fRunPar->foutputfileName != "-1" )
    {
        // tree versioning numbers used in mscw_energy
        stringstream i_sst;
        i_sst << "VERSION " << getRunParameter()->getEVNDISP_TREE_VERSION();
        if( getRunParameter()->fShortTree )
        {
            i_sst << " (short tree)";
        }
        fOutputfile = new TFile( fRunPar->foutputfileName.c_str(), "RECREATE", i_sst.str().c_str() );
    }
}


/*!
   create a new tree in the top directory of the analysis output file
   (to be called once before starting the array analysis)
*/
void VArrayAnalyzer::initTree()
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::initTrees() " <<  endl;
    }
    if( !fOutputfile )
    {
        cout <<  "VArrayAnalyzer::initTrees() warning: No output file, will not initialize trees." <<  endl;
        return;
    }
    fOutputfile->cd();
    // now book the tree (for MC with additional MC block)
    // tree versioning numbers used in mscw_energy
    stringstream i_sst;
    i_sst << "Shower Parameters (VERSION " << getRunParameter()->getEVNDISP_TREE_VERSION() << ")";
    if( getRunParameter()->fShortTree )
    {
        i_sst << "(short tree)";
    }
    fShowerParameters->initTree( "showerpars", i_sst.str(), fReader->isMC() );
    if( isMC() && fMCParameters )
    {
        fMCParameters->initTree();
    }
}


/*!
    write result tree into ROOT output file
    (called once after finishing the analysis)
*/
void VArrayAnalyzer::terminate( bool iWriteDebug )
{
    if( fDebug )
    {
        cout << "void VArrayAnalyzer::terminate()" << endl;
    }
    
    generateReducedPointingTreeData() ;
    
    for( unsigned int i = 0; i < fDispAnalyzer.size(); i++ )
    {
        if( fDispAnalyzer[i] )
        {
            fDispAnalyzer[i]->terminate();
        }
    }
    
    if( fOutputfile != 0 && fRunPar->foutputfileName != "-1" )
    {
        // write pointing mismatch
        if( !getNoPointing() && getArrayPointing() && getArrayPointing()->getTargetName() != "laser" && !fReader->isMC() )
        {
            cout << endl;
            cout << "Pointing check results (using VBF data):" << endl;
            for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
            {
                if( getMeanPointingMismatch( getTeltoAna()[i] ) > 0.1 )
                {
                    cout << endl;
                    cout << "WARNING: LARGE MISMATCH BETWEEN EVENTDISPLAY AND VBF POINTING DATA FOR TELESCOPE " << getTeltoAna()[i] + 1;
                    cout << " (mean mismatch is " << getMeanPointingMismatch( getTeltoAna()[i] ) << " deg)" << endl;
                    cout << "\t IS THE WOBBLE OFFSET RIGHT?" << endl;
                }
                else
                {
                    cout << "\t mean pointing mismatch between pointing data (analysis and VBF) for Telescope " << getTeltoAna()[i] + 1 << ": ";
                    cout << getMeanPointingMismatch( getTeltoAna()[i] ) << " deg" << endl;
                }
            }
        }
        fOutputfile->cd();
        
        if( !getShowerParameters()->getTree() )
        {
            return;
        }
        cout << "---------------------------------------------------------------------------------------------------------" << endl;
        cout << "writing analysis tree (";
        cout << getShowerParameters()->getTree()->GetEntries();
        cout << " entries) to : " << fOutputfile->GetName() << endl;
        int i_nbytes = getShowerParameters()->getTree()->Write();
        if( iWriteDebug )
        {
            cout << "WRITEDEBUG: showerpars tree (nbytes " << i_nbytes << "): ";
            if( fOutputfile )
            {
                cout << fOutputfile->Get( "showerpars" );
            }
            cout << endl;
        }
        ///////////////////////////////////////////////////////////////
        // MC tree and histograms
        if( isMC() )
        {
            // get MC tree
            TTree* i_tMC = 0;
            if( getRunParameter()->fsourcetype == 7 )
            {
                if( getReader()->getMCTree() )
                {
                    if( getRunParameter()->fwriteMCtree )
                    {
                        if( getRunParameter()->fwriteMCtree )
                        {
                            i_tMC = ( TTree* )getReader()->getMCTree()->CloneTree( -1, "fast" );
                        }
                        if( i_tMC )
                        {
                            i_tMC->SetName( "MCpars" );
                        }
                    }
                    else
                    {
                        i_tMC = ( TTree* )getReader()->getMCTree();
                    }
                }
            }
            else
            {
                if( getMCParameters()->getTree() )
                {
                    i_tMC = getMCParameters()->getTree();
                }
            }
            // MC tree is not written by default to output file (takes a while)
            if( getRunParameter()->fwriteMCtree && i_tMC )
            {
                cout << "writing MC tree " << endl;
                i_tMC->Write();
                cout << "(number of MC events in tree " << i_tMC->GetName() << ": " << i_tMC->GetEntries() << ")" << endl;
                cout << "...done" << endl;
            }
            if( getRunParameter()->fFillMCHistos && i_tMC )
            {
                ////////////////////////////////////////////////////////////
                // filling of MC histograms -> needed for effective area calculations
                cout << "filling MC histograms" << endl;
                VEffectiveAreaCalculatorMCHistograms iMC_histos;
                double i_ze = 0.;
                // set energy range for spectral weighting
                if( !getReader()->getMonteCarloHeader() )
                {
                    cout << "VArrayAnalyzer::terminate warning: no simulation run header in vbf file found" << endl;
                    cout << "Please check the consistency of your MC vbf file" << endl;
                }
                else
                {
                    iMC_histos.setMonteCarloEnergyRange( getReader()->getMonteCarloHeader()->E_range[0],
                                                         getReader()->getMonteCarloHeader()->E_range[1],
                                                         TMath::Abs( getReader()->getMonteCarloHeader()->spectral_index ) );
                    i_ze = getReader()->getMonteCarloHeader()->getMeanZenithAngle_Deg();
                }
                iMC_histos.setDefaultValues();
                iMC_histos.initializeHistograms();
                // backwards compatibility: no ze in MC run header for old files
                if( TMath::Abs( i_ze - 90. ) < 1.e-3 )
                {
                    cout << "\t MC run header seems to be of old type without correct zenith/azimuth range" << endl;
                    i_ze = getShowerParameters()->MCze;
                }
                iMC_histos.fill( i_ze, i_tMC, true );
                iMC_histos.print();
                int i_nbytes = iMC_histos.Write();
                if( iWriteDebug )
                {
                    cout << "WRITEDEBUG: MC histograms (nbytes " << i_nbytes << "):";
                    if( fOutputfile )
                    {
                        cout << fOutputfile->Get( "MChistos" );
                    }
                    cout << endl;
                }
                // END: filling of MC histograms
                ////////////////////////////////////////////////////////////
            }
        }
        if( getEvndispReconstructionParameter() )
        {
            int i_nbytes = getEvndispReconstructionParameter()->Write();
            if( iWriteDebug )
            {
                cout << "WRITEDEBUG: evndisp reconstruction parameters (nbytes " << i_nbytes << "): ";
                if( fOutputfile )
                {
                    cout << fOutputfile->Get( "EvndispReconstructionParameter" );
                }
                cout << endl;
            }
        }
        if( getArrayPointing() )
        {
            getArrayPointing()->terminate( iWriteDebug, isMC() );
        }
        if( fOutputfile )
        {
            fOutputfile->Flush();
        }
        cout << "---------------------------------------------------------------------------------------------------------" << endl;
    }
}


/*!
          use GrIsu routines for the actual transformation

      \param iTel telescope number
      \param i_ze zenith angle [deg]
      \param i_az azimuth angle [deg]
      \param i_MC data is MC data (true)
*/
void VArrayAnalyzer::transformTelescopePosition( int iTel, float i_ze, float i_az, bool i_MC )
{
    // transform telescope positions from ground into shower coordinates
    float i_xrot, i_yrot, i_zrot;
    float i_xcos = 0.;
    float i_ycos = 0.;
    // calculate direction cosinii
    i_xcos = sin( i_ze / TMath::RadToDeg() ) * sin( ( i_az - 180. ) / TMath::RadToDeg() );
    i_ycos = sin( i_ze / TMath::RadToDeg() ) * cos( ( i_az - 180. ) / TMath::RadToDeg() );
    
    setTelID( iTel );
    // call to GrIsu routine
    tel_impact( i_xcos, i_ycos, getDetectorGeo()->getTelXpos()[iTel], getDetectorGeo()->getTelYpos()[iTel], getDetectorGeo()->getTelZpos()[iTel], &i_xrot, &i_yrot, &i_zrot, false );
    
    getImageParameters()->Tel_x_SC = i_xrot;
    getImageParameters()->Tel_y_SC = i_yrot;
    getImageParameters()->Tel_z_SC = i_zrot;
    
    if( getRunParameter()->fImageLL )
    {
        getImageParametersLogL()->Tel_x_SC = i_xrot;
        getImageParametersLogL()->Tel_y_SC = i_yrot;
        getImageParametersLogL()->Tel_z_SC = i_zrot;
    }
    
    getShowerParameters()->fTel_x_SC[iTel] = i_xrot;
    getShowerParameters()->fTel_y_SC[iTel] = i_yrot;
    getShowerParameters()->fTel_z_SC[iTel] = i_zrot;
    
    if( i_MC )
    {
        // transformation of Monte Carlo shower core into shower coordinates
        // call to GrIsu routine
        tel_impact( i_xcos, i_ycos, getShowerParameters()->MCxcore, getShowerParameters()->MCycore, getShowerParameters()->MCzcore, &i_xrot, &i_yrot, &i_zrot, false );
        getShowerParameters()->MCxcore_SC = i_xrot;
        getShowerParameters()->MCycore_SC = i_yrot;
        getShowerParameters()->MCzcore_SC = i_zrot;
    }
}


/*!
     select images used in shower reconstruction

     cuts possible on ntubes, dist, alpha (set with in files with arraycuts,
     see example array_analysis_cuts.txt)

     results stored in boolean vector with length of number of telescopes

     parameter: iMeth shower reconstruction method number
*/
void VArrayAnalyzer::selectShowerImages( unsigned int iMeth )
{
    bool fSelectDebug = false;
    
    getShowerParameters()->fTelIDImageSelected[iMeth].clear();
    getShowerParameters()->fTelIDImageSelected_bitcode[iMeth] = 0;
    getShowerParameters()->fShowerNumImages[iMeth] = 0;
    
    // make sure that reconstruction parameters exist
    if( !getEvndispReconstructionParameter() )
    {
        return;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    // loop over all telescopes and check which image is suitable for reconstruction
    for( unsigned int t = 0; t < getNTel(); t++ )
    {
        setTelID( t );
        
        // reset list with selected images
        getShowerParameters()->fTelIDImageSelected_list[iMeth][t] = 0;
        getShowerParameters()->fTelIDImageSelected[iMeth].push_back( true );
        
        if( fSelectDebug )
        {
            cout << "VArrayAnalyzer::selectShowerImages " << getTelID() + 1 << "\t" << getEventNumber() << "\t" << iMeth << endl;
        }
        if( fSelectDebug )
        {
            cout << "VArrayAnalyzer::selectShowerImages eventstatus " << getImageParameters( getRunParameter()->fImageLL )->eventStatus << endl;
        }
        
        // get telescope type for this telescope
        int iTelType = getEvndispReconstructionParameter()->getTelescopeType_counter( getDetectorGeometry()->getTelType()[t] );
        if( iTelType < 0 )
        {
            cout << "VArrayAnalyzer::selectShowerImages error: invalid telescope counter: " << t << "\t" << iTelType << endl;
            exit( EXIT_FAILURE );
        }
        
        // apply array analysis cuts
        updatePointingToStarCatalogue( t );
        getShowerParameters()->fTelIDImageSelected[iMeth].back() =
            getEvndispReconstructionParameter()->applyArrayAnalysisCuts( iMeth, t, iTelType,
                    getImageParameters( getRunParameter()->fImageLL ),
                    getReader()->getLocalTriggerType( t ),
                    getStarCatalogue() );
        // apply cut on image distance to stars
        
        ///////////////////////////
        
        // check if fit was successfull
        if( getRunParameter()->fImageLL && getImageParametersLogL()->Fitstat < 3 )
        {
            getShowerParameters()->fTelIDImageSelected[iMeth].back() = false;
        }
        
        // list of selected images
        if( getShowerParameters()->fTelIDImageSelected[iMeth].back() )
        {
            getShowerParameters()->fTelIDImageSelected_list[iMeth][getTelID()] = 1;
            getShowerParameters()->fShowerNumImages[iMeth]++;
        }
        
        bitset<8 * sizeof( unsigned long )> i_nimage;
        if( fNTel < i_nimage.size() )
        {
            for( unsigned int i = 0; i < getNTel(); i++ )
            {
                if( i < i_nimage.size() )
                {
                    if( getShowerParameters()->fTelIDImageSelected[iMeth][i] )
                    {
                        i_nimage.set( i, 1 );
                    }
                    else
                    {
                        i_nimage.set( i, 0 );
                    }
                }
            }
            getShowerParameters()->fTelIDImageSelected_bitcode[iMeth] = i_nimage.to_ulong();
        }
        else
        {
            getShowerParameters()->fTelIDImageSelected_bitcode[iMeth] = 0;
        }
    }
}


void VArrayAnalyzer::calcShowerDirection_and_Core()
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::calcShowerDirection_and_Core()" << endl;
    }
    
    // loop over all methods
    for( unsigned int i = 0; i < getEvndispReconstructionParameter()->getNReconstructionCuts(); i++ )
    {
        // set reconstruction method
        if( getEvndispReconstructionParameter( i ) )
        {
            getShowerParameters()->fMethodID[i] = getEvndispReconstructionParameter( i )->fMethodID;
            // select shower images to be used to determinate shower coordinates
            selectShowerImages( i );
            
            unsigned int iReconstructionMethodID = getEvndispReconstructionParameter( i )->fMethodID;
            // DISP method: use only up to a certain number of images
            // (used for one zenith angle only)
            if( getEvndispReconstructionParameter( i )->fDISP_MAXTelescopes.size() == 1
                    && getEvndispReconstructionParameter( i )->fDISP_MAXMethodID.size() == 1
                    && getShowerParameters()->fShowerNumImages[i] > getEvndispReconstructionParameter( i )->fDISP_MAXTelescopes[0] )
            {
                iReconstructionMethodID = getEvndispReconstructionParameter( i )->fDISP_MAXMethodID[0];
            }
            
            // call reconstruction method
            if( iReconstructionMethodID == 0 )
            {
                rcs_method_0( i );
            }
            else if( iReconstructionMethodID == 3 )
            {
                rcs_method_3( i );
            }
            else if( iReconstructionMethodID == 4 )
            {
                rcs_method_4( i );
            }
            else if( iReconstructionMethodID == 5 )
            {
                rcs_method_5( i, 5 );
            }
            else if( iReconstructionMethodID == 6 )
            {
                rcs_method_5( i, 6 );
            }
            else if( iReconstructionMethodID == 7 )
            {
                rcs_method_5( i, 7 );
            }
            else if( iReconstructionMethodID == 8 )
            {
                rcs_method_8( i );
            }
            else
            {
                cout << "VArrayAnalyzer::calcShowerDirection_and_Core(): unknown array reconstruction method ";
                cout << iReconstructionMethodID;
                cout << endl;
                continue;
            }
        }
        else
        {
            cout << "VArrayAnalyzer::calcShowerDirection_and_Core(): error, analysis vectors with different size" << endl;
            cout << "(nmethods: " << getEvndispReconstructionParameter()->getNReconstructionCuts() << ")" << endl;
            cout << "...fatal error" << endl;
            exit( EXIT_FAILURE );
        }
    }
}


/*!
      reconstruction of shower direction and core

      shower direction by intersection of image axes
      shower core by intersection of lines connecting reconstruced shower
      direction and image centroids

      basic code is a copy from Charlie Duke's routines in GrIsu

      adjustments:
      * technical stuff to run in the eventdisplay frame
* weighting
*/
int VArrayAnalyzer::rcs_method_0( unsigned int iMethod )
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::rcs_method_0 " << iMethod << endl;
    }
    
    /* Requires phi, xc, yc, and sum image parameters
    
       use image parameters, phi, xc, and yc to determine image
       axes. Use rcs_perpendicular_fit to find source point.
    
       find impact point by rotation through delta to position
       source point at camera center, then find impact point
       using ground lines connecting the reconstructed shower
    direction with the image centroids. Impact
       point found using rcs_perpendicular_fit.
    */
    int num_images = 0;                           /* number of images included in the fit */
    
    float xs = 0.;
    float ys = 0.;
    float stds = 0.;
    
    float ximp = 0.;                              /* impact point, from perpen. minimization */
    float yimp = 0.;
    float stdp = 0.;
    
    if( !getEvndispReconstructionParameter( iMethod ) )
    {
        return 0;
    }
    
    num_images = getShowerParameters()->fShowerNumImages[iMethod];
    
    // are there enough images the run an array analysis
    if( num_images >= ( int )getEvndispReconstructionParameter( iMethod )->fNImages_min )
    {
        prepareforDirectionReconstruction( iMethod, 0 );
    }
    else
    {
        fillShowerDirection( iMethod, 0., 0., -1. );
        return 0;
    }
    
    // don't do anything if angle between image axis is too small (for 2 images only)
    if( num_images == 2 && !testMinAngleBetweenImages( iMethod, fabs( atan( m[0] ) - atan( m[1] ) ) ) )
    {
        return 0;
    }
    else
    {
        getShowerParameters()->fiangdiff[iMethod] = 0.;
    }
    
    ///////////////////////////////
    // direction reconstruction
    ///////////////////////////////
    // calculate offset in camera coordinates
    
    rcs_perpendicular_fit( x, y, w, m, num_images, &xs, &ys, &stds );
    
    // assume all focal lengths of all telescopes are the same
    xs = atan( xs / ( 1000. * getDetectorGeo()->getFocalLength()[0] ) ) * TMath::RadToDeg();
    ys = atan( ys / ( 1000. * getDetectorGeo()->getFocalLength()[0] ) ) * TMath::RadToDeg();
    stds = atan( stds / ( 1000. * getDetectorGeo()->getFocalLength()[0] ) ) * TMath::RadToDeg();
    
    // calculate and fill directions
    if( !fillShowerDirection( iMethod, xs, ys, stds ) )
    {
        return 0;
    }
    
    // end of direction reconstruction
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    
    /* - - - - - now find the impact location in the plane
        perpendicular to the array. Since we know the source location,
        we can rotate the array to point in this direction thus moving
        the source point to the center of the camera.  Source location
        will be at (0, 0).
    */
    prepareforCoreReconstruction( iMethod, xs, ys );
    
    /* Now call perpendicular_distance for the fit, returning ximp and yimp */
    
    rcs_perpendicular_fit( x, y, w, m, num_images, &ximp, &yimp, &stdp );
    
    if( isnormal( ximp ) && isnormal( yimp ) )
    {
        fillShowerCore( iMethod, ximp, yimp );
        getShowerParameters()->fShower_stdP[iMethod] = stdp;
        /* indicates fit ok */
        getShowerParameters()->fShower_Chi2[iMethod] = 0.0;
    }
    // this happens sometimes for two-telescope events when image axis are exactly parallel
    else
    {
        getShowerParameters()->fShower_stdP[iMethod] = 0.;
        getShowerParameters()->fShower_Chi2[iMethod] = -2.;
    }
    
    return 0;
}


/***************** end of rcs_method_0 *********************************/

/*!
    \par iMeth reconstruction method number
    \par ximp shower core position in shower coordinates
    \par yimp shower core position in shower coordinates
*/
bool VArrayAnalyzer::fillShowerCore( unsigned int iMeth, float ximp, float yimp )
{
    // check validity
    if( !isnormal( ximp ) || !isnormal( yimp ) )
    {
        ximp = -99999.;
        yimp = -99999.;
        return false;
    }
    // reconstructed shower core in shower coordinates
    getShowerParameters()->fShowerXcore_SC[iMeth] = ximp;
    getShowerParameters()->fShowerYcore_SC[iMeth] = yimp;
    
    // reconstructed shower core in ground coordinates
    float i_xcos = 0.;
    float i_ycos = 0.;
    float zimp = 0.;
    float igx, igy, igz;
    // calculate z in shower coordinates (for z=0 in ground coordinates)
    if( getShowerParameters()->fShowerZe[iMeth] != 0. )
    {
        zimp = yimp / tan( ( 90. - getShowerParameters()->fShowerZe[iMeth] ) / TMath::RadToDeg() );
    }
    else
    {
        zimp = 0.;
    }
    // calculate direction cosinii
    // taking telescope plane as reference plane. This is not exactly correct for method 0
    // taking pointing direction of first telescope in teltoana vector
    if( getArrayPointing() )
    {
        i_xcos = sin( ( 90. - getArrayPointing()->getTelElevation() ) / TMath::RadToDeg() )
                 * sin( ( getArrayPointing()->getTelAzimuth() - 180. ) / TMath::RadToDeg() );
        if( fabs( i_xcos ) < 1.e-7 )
        {
            i_xcos = 0.;
        }
        i_ycos = sin( ( 90. - getArrayPointing()->getTelElevation() ) / TMath::RadToDeg() )
                 * cos( ( getArrayPointing()->getTelAzimuth() - 180. ) / TMath::RadToDeg() );
        if( fabs( i_ycos ) < 1.e-7 )
        {
            i_ycos = 0.;
        }
    }
    else
    {
        i_xcos = 0.;
        i_ycos = 0.;
    }
    tel_impact( i_xcos, i_ycos, ximp, yimp, zimp, &igx, &igy, &igz, true );
    if( isinf( igx ) || isinf( igy ) || TMath::IsNaN( igx ) || TMath::IsNaN( igy ) )
    {
        igx = -99999;
        igy = -99999;
        getShowerParameters()->fShower_Chi2[iMeth] = -2.;
    }
    getShowerParameters()->fShowerXcore[iMeth] = igx;
    getShowerParameters()->fShowerYcore[iMeth] = igy;
    
    // z coordinate should be zero in ground coordinate system, if not -> problem
    if( fabs( igz ) > 1.e3 )
    {
        return false;
    }
    
    return true;
}


void VArrayAnalyzer::checkPointing()
{
    // there is no pointing available
    if( getNoPointing() )
    {
        return;
    }
    
    // calculate difference between calculated pointing direction and vbf pointing direction
    if( getReader()->getArrayTrigger() )
    {
        for( unsigned int j = 0; j < getTeltoAna().size(); j++ )
        {
            unsigned int i = getTeltoAna()[j];
            // get vbf telescope index
            unsigned int ivbf = 9999;
            for( int t = 0; t < ( int )getReader()->getArrayTrigger()->getNumSubarrayTelescopes(); t++ )
            {
                if( getTeltoAna()[j] == getReader()->getArrayTrigger()->getSubarrayTelescopeId( t ) )
                {
                    ivbf = t;
                    break;
                }
            }
            if( ivbf == 9999 )
            {
                continue;
            }
            
            getShowerParameters()->fTelElevationVBF[i] = getReader()->getArrayTrigger()->getAltitude( ivbf );
            getShowerParameters()->fTelAzimuthVBF[i]   = getReader()->getArrayTrigger()->getAzimuth( ivbf );
            if( i < getPointing().size() && getPointing()[i] )
            {
                float iPointingDiff = VAstronometry::vlaDsep( getPointing()[i]->getTelAzimuth() * TMath::DegToRad(),
                                                             ( 90. - getPointing()[i]->getTelElevation() ) * TMath::DegToRad(),
                                                             getReader()->getArrayTrigger()->getAzimuth( ivbf )* TMath::DegToRad() ,
                                                             ( 90. - getReader()->getArrayTrigger()->getAltitude( ivbf ) ) * TMath::DegToRad() )
                                                             * TMath::RadToDeg();
                getShowerParameters()->fTelPointingMismatch[i] = iPointingDiff;
                getShowerParameters()->fTelPointingErrorX[i] = getPointing()[i]->getPointingErrorX();
                getShowerParameters()->fTelPointingErrorY[i] = getPointing()[i]->getPointingErrorY();
                
                fMeanPointingMismatch[i] += iPointingDiff;
                fNMeanPointingMismatch[i]++;
                
                // check pointing difference, abort if too large
                if( getRunParameter()->fCheckPointing < 900. && iPointingDiff > getRunParameter()->fCheckPointing )
                {
                    cout << "VArrayAnalyzer::checkPointing() large mismatch between calculated telescope pointing direction and VBF pointing: ";
                    cout << iPointingDiff << " deg" << endl;
                    cout << "\t calculated telescope pointing direction: ";
                    cout << getPointing()[i]->getTelAzimuth() << "\t" << getPointing()[i]->getTelElevation() << endl;
                    cout << "\t VBF telescope pointing direction:        ";
                    cout << getReader()->getArrayTrigger()->getAzimuth( ivbf ) << "\t" << getReader()->getArrayTrigger()->getAltitude( ivbf ) << endl;
                    cout << "ABORT due to large pointing error (>" << getRunParameter()->fCheckPointing << " deg, event " << getEventNumber() << ")" << endl;
                    exit( EXIT_FAILURE );
                }
            }
            else
            {
                getShowerParameters()->fTelPointingMismatch[i] = 0.;
                getShowerParameters()->fTelPointingErrorX[i] = 0.;
                getShowerParameters()->fTelPointingErrorY[i] = 0.;
            }
        }
    }
}


double VArrayAnalyzer::getMeanPointingMismatch( unsigned int iTel )
{
    if( getNoPointing() )
    {
        return -1.;
    }
    if( iTel > fMeanPointingMismatch.size() )
    {
        return -2.;
    }
    if( TMath::IsNaN( fMeanPointingMismatch[iTel] ) == 0 && TMath::IsNaN( fNMeanPointingMismatch[iTel] ) == 0 )
    {
        if( fNMeanPointingMismatch[iTel] > 0. )
        {
            if( TMath::Abs( fMeanPointingMismatch[iTel] / fNMeanPointingMismatch[iTel] ) < 1.e-5 )
            {
                return 0.;
            }
            else
            {
                return fMeanPointingMismatch[iTel] / fNMeanPointingMismatch[iTel];
            }
        }
    }
    
    return -3.;
}


/*
    add a small pointing error to centroids and recalculate image angle phi
*/
float VArrayAnalyzer::recalculateImagePhi( double iDeltaX, double iDeltaY )
{
    float i_phi = 0.;
    // LL: fits without good error matrix
    // LL: f_d, s_s and f_sdevxy depend on the state of the error matrix
    if( getImageParameters( getRunParameter()->fImageLL )->Fitstat == 1 )
    {
        i_phi = getImageParameters( getRunParameter()->fImageLL )->phi;
    }
    else
    {
        double xmean = getImageParameters( getRunParameter()->fImageLL )->cen_x + iDeltaX;
        double ymean = getImageParameters( getRunParameter()->fImageLL )->cen_y + iDeltaY;
        double d = getImageParameters( getRunParameter()->fImageLL )->f_d;
        double z = getImageParameters( getRunParameter()->fImageLL )->f_s;
        double sdevxy = getImageParameters( getRunParameter()->fImageLL )->f_sdevxy;
        const double ac = ( d + z ) * ymean + 2.0 * sdevxy * xmean;
        const double bc = 2.0 * sdevxy * ymean - ( d - z ) * xmean;
        const double cc = sqrt( ac * ac + bc * bc );
        i_phi = atan2( ac / cc, bc / cc );
    }
    
    return i_phi;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
/*!
      reconstruction of shower direction and core

      Hofmann et al 1999, Method 1 (HEGRA method)

      shower direction by intersection of image axes
      shower core by intersection of lines connecting reconstruced shower
      direction and image centroids

      todo: core reconstruction

*/
int VArrayAnalyzer::rcs_method_3( unsigned int iMethod )
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::rcs_method_3 " << iMethod << endl;
    }
    
    int num_images = 0;                           /* number of images included in the fit */
    
    float xs = 0.;
    float ys = 0.;
    float stds = 0.;
    
    float ximp = 0.;                              /* impact point, from perpen. minimization */
    float yimp = 0.;
    float stdp = 0.;
    
    num_images = getShowerParameters()->fShowerNumImages[iMethod];
    
    if( !getEvndispReconstructionParameter( iMethod ) )
    {
        return 0;
    }
    
    // are there enough images the run an array analysis
    if( num_images >= ( int )getEvndispReconstructionParameter( iMethod )->fNImages_min )
    {
        prepareforDirectionReconstruction( iMethod, 3 );
    }
    else
    {
        fillShowerDirection( iMethod, 0., 0., -1. );
        return 0;
    }
    // don't do anything if angle between image axis is too small (for 2 images only)
    if( num_images == 2 && !testMinAngleBetweenImages( iMethod, fabs( atan( m[0] ) - atan( m[1] ) ) ) )
    {
        return 0;
    }
    else
    {
        getShowerParameters()->fiangdiff[iMethod] = 0.;
    }
    
    ///////////////////////////////
    // direction reconstruction
    ///////////////////////////////
    /* Now call perpendicular_distance for the fit, returning xs and ys */
    
    // Hofmann et al 1999, Method 1 (HEGRA method)
    
    vector< float > xx( 2, 0. );
    vector< float > yy( 2, 0. );
    vector< float > ww( 2, 0. );
    vector< float > mm( 2, 0. );
    float itotweight = 0.;
    float iweight = 1.;
    float ixs = 0.;
    float iys = 0.;
    float iangdiff = 0.;
    for( unsigned int ii = 0; ii < m.size(); ii++ )
    {
        for( unsigned int jj = 1; jj < m.size(); jj++ )
        {
            if( ii == jj || ii > jj )
            {
                continue;
            }
            xx[0] = x[ii];
            yy[0] = y[ii];
            ww[0] = 1.;
            mm[0] = m[ii];
            xx[1] = x[jj];
            yy[1] = y[jj];
            ww[1] = 1.;
            mm[1] = m[jj];
            
            rcs_perpendicular_fit( xx, yy, ww, mm, 2, &xs, &ys, &stds );
            
            iangdiff = sin( fabs( atan( m[jj] ) - atan( m[ii] ) ) );
            
            // discard all pairs with almost parallel lines
            float i_diff =  fabs( atan( m[0] ) - atan( m[1] ) );
            if( i_diff < getEvndispReconstructionParameter( iMethod )->fAxesAngles_min / TMath::RadToDeg() ||
                    fabs( 180. / TMath::RadToDeg() - i_diff ) < getEvndispReconstructionParameter( iMethod )->fAxesAngles_min / TMath::RadToDeg() )
            {
                continue;
            }
            
            iweight  = 1. / ( 1. / w[ii] + 1. / w[jj] );
            iweight *= 1. / ( l[ii] + l[jj] );
            iweight *= iangdiff;
            
            ixs += atan( xs / ( 1000. * getDetectorGeo()->getFocalLength()[0] ) ) * iweight * TMath::RadToDeg();
            iys += atan( ys / ( 1000. * getDetectorGeo()->getFocalLength()[0] ) ) * iweight * TMath::RadToDeg();
            itotweight += iweight;
        }
    }
    if( itotweight > 0. )
    {
        xs = ixs / itotweight;
        ys = iys / itotweight;
    }
    else
    {
        xs = -99999.;
        ys = -99999.;
    }
    
    if( !fillShowerDirection( iMethod, xs, ys, stds ) )
    {
        return 0;
    }
    
    // core reconstruction
    
    prepareforCoreReconstruction( iMethod, xs, ys );
    
    // Now call perpendicular_distance for the fit, returning ximp and yimp
    
    rcs_perpendicular_fit( x, y, w, m, num_images, &ximp, &yimp, &stdp );
    
    // convert xs,ys into [deg]
    
    if( isnormal( ximp ) && isnormal( yimp ) )
    {
        fillShowerCore( iMethod, ximp, yimp );
        getShowerParameters()->fShower_stdP[iMethod] = stdp;
        /* indicates fit ok */
        getShowerParameters()->fShower_Chi2[iMethod] = 0.0;
    }
    else
    {
        getShowerParameters()->fShower_Chi2[iMethod] = -2.;
    }
    
    return 0;
}


/***************** end of rcs_method_3 *********************************/

/*!
      reconstruction of shower direction and core

      Hofmann et al 1999, Method 1 (HEGRA method)

      shower direction by intersection of image axes
      shower core by intersection of lines connecting reconstruced shower
      direction and image centroids

      (difference to rcs_method_3: direction reconstruction independent of CD routines)

todo: core construction
*/
int VArrayAnalyzer::rcs_method_4( unsigned int iMethod )
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::rcs_method_4 " << iMethod << endl;
    }
    int num_images = 0;                           /* number of images included in the fit */
    
    float xs = 0.;
    float ys = 0.;
    
    float ximp = 0.;                              /* impact point, from perpen. minimization */
    float yimp = 0.;
    float stdp = 0.;
    
    num_images = getShowerParameters()->fShowerNumImages[iMethod];
    
    if( !getEvndispReconstructionParameter( iMethod ) )
    {
        return 0;
    }
    
    // are there enough images the run an array analysis
    if( num_images >= ( int )getEvndispReconstructionParameter( iMethod )->fNImages_min )
    {
        prepareforDirectionReconstruction( iMethod, 4 );
    }
    else
    {
        fillShowerDirection( iMethod, 0., 0., -1. );
        return 0;
    }
    
    // don't do anything if angle between image axis is too small (for 2 images only)
    if( num_images == 2 && !testMinAngleBetweenImages( iMethod, fabs( atan( m[0] ) - atan( m[1] ) ) ) )
    {
        return 0;
    }
    else
    {
        getShowerParameters()->fiangdiff[iMethod] = 0.;
    }
    
    ///////////////////////////////
    // direction reconstruction
    ////////////////////////////////////////////////
    // Hofmann et al 1999, Method 1 (HEGRA method)
    // (modified weights)
    
    float itotweight = 0.;
    float iweight = 1.;
    float ixs = 0.;
    float iys = 0.;
    float iangdiff = 0.;
    float b1 = 0.;
    float b2 = 0.;
    vector< float > v_xs;
    vector< float > v_ys;
    
    double i_weight_max = 0.;
    
    for( unsigned int ii = 0; ii < m.size(); ii++ )
    {
        for( unsigned int jj = 1; jj < m.size(); jj++ )
        {
            if( ii == jj || ii > jj )
            {
                continue;
            }
            
            // check minimum angle between image lines; ignore if too small
            iangdiff = fabs( atan( m[jj] ) - atan( m[ii] ) );
            if( iangdiff < getEvndispReconstructionParameter( iMethod )->fAxesAngles_min / TMath::RadToDeg() ||
                    fabs( 180. / TMath::RadToDeg() - iangdiff ) < getEvndispReconstructionParameter( iMethod )->fAxesAngles_min / TMath::RadToDeg() )
            {
                continue;
            }
            // weight is sin of angle between image lines
            iangdiff = fabs( sin( fabs( atan( m[jj] ) - atan( m[ii] ) ) ) );
            
            b1 = y[ii] - m[ii] * x[ii];
            b2 = y[jj] - m[jj] * x[jj];
            
            // line intersection
            if( m[ii] != m[jj] )
            {
                xs = ( b2 - b1 )  / ( m[ii] - m[jj] );
            }
            else
            {
                xs = 0.;
            }
            ys = m[ii] * xs + b1;
            
            
            iweight  = 1. / ( 1. / w[ii] + 1. / w[jj] ); // weight 1: size of images
            iweight *= ( 1. - l[ii] ) * ( 1. - l[jj] ); // weight 2: elongation of images (width/length)
            iweight *= iangdiff;                      // weight 3: angular differences between the two image axis
            iweight *= iweight;                       // use squared value
            
            if( iweight > i_weight_max )
            {
                i_weight_max = iweight;
            }
            
            ixs += xs * iweight;
            iys += ys * iweight;
            itotweight += iweight;
            
            v_xs.push_back( xs );
            v_ys.push_back( ys );
        }
    }
    // check validity of weight
    if( itotweight > 0. )
    {
        xs = ixs / itotweight;
        ys = iys / itotweight;
        // calculate dispdiff
        // (this is not exactly dispdiff, but
        //  an equivalent measure comparable to dispdiff)
        float fShower_DispDiff = 0.;
        for( unsigned n = 0; n < v_xs.size(); n++ )
        {
            for( unsigned m = n; m < v_xs.size(); m++ )
            {
                fShower_DispDiff += ( v_xs[n] - v_xs[m] ) * ( v_xs[n] - v_xs[m] );
                fShower_DispDiff += ( v_ys[n] - v_ys[m] ) * ( v_ys[n] - v_ys[m] );
            }
        }
        getShowerParameters()->fDispDiff[iMethod] = fShower_DispDiff;
    }
    else
    {
        xs = -99999.;
        ys = -99999.;
        getShowerParameters()->fDispDiff[iMethod] = -99999.;
    }
    
    // fill correct shower direction
    // do not continue with core reconstruction in case there is no valid
    // direction reconstruction
    if( !fillShowerDirection( iMethod, xs, ys, 0. ) )
    {
        return 0;
    }
    
    ////////////////////////////////////////////////
    
    // core reconstruction
    
    prepareforCoreReconstruction( iMethod, xs, ys );
    
    // Now call perpendicular_distance for the fit, returning ximp and yimp
    
    rcs_perpendicular_fit( x, y, w, m, num_images, &ximp, &yimp, &stdp );
    
    // convert xs,ys into [deg]
    
    fillShowerCore( iMethod, ximp, yimp );
    getShowerParameters()->fShower_stdP[iMethod] = stdp;
    /* indicates fit ok */
    getShowerParameters()->fShower_Chi2[iMethod] = 0.0;
    
    return 0;
}


/***************** end of rcs_method_4 *********************************/

/*
 * fill variables related to shower direction in camera, ground and sky coordinates
 *
 * check validity of reconstruction
 *
 */
bool VArrayAnalyzer::fillShowerDirection( unsigned int iMethod, float xs, float ys, float stds )
{
    // required successfull reconstruction
    if( TMath::IsNaN( xs ) || TMath::IsNaN( ys ) || xs < -9998. || ys < -9998. )
    {
        getShowerParameters()->fShower_Xoffset[iMethod] = -99999.;
        getShowerParameters()->fShower_Yoffset[iMethod] = -99999.;
        getShowerParameters()->fShowerZe[iMethod] = -99999.;
        getShowerParameters()->fShowerAz[iMethod] = -99999.;
        getShowerParameters()->fShower_stdS[iMethod] = -99999.;
        getShowerParameters()->fShower_Chi2[iMethod] = -99999.;
        getShowerParameters()->fShower_XoffsetDeRot[iMethod] = -99999.;
        getShowerParameters()->fShower_YoffsetDeRot[iMethod] = -99999.;
        
        return false;
    }
    
    getShowerParameters()->fShower_Xoffset[iMethod] = xs;
    // convention in coordinate system
    getShowerParameters()->fShower_Yoffset[iMethod] = -1. * ys;
    
    double ze = 0.;
    double az = 0.;
    if( getArrayPointing() )
    {
        getArrayPointing()->getRotatedShowerDirection( -1.*getShowerParameters()->fShower_Yoffset[iMethod],
                -1.*getShowerParameters()->fShower_Xoffset[iMethod], ze, az );
    }
    if( TMath::IsNaN( ze ) )
    {
        ze = -99999.;
    }
    if( TMath::IsNaN( az ) )
    {
        az = -99999.;
    }
    getShowerParameters()->fShowerZe[iMethod] = ze;
    getShowerParameters()->fShowerAz[iMethod] = VSkyCoordinatesUtilities::adjustAzimuthToRange( az );
    getShowerParameters()->fShower_stdS[iMethod] = stds;
    
    // calculate derotated shower directions in camera coordinates
    // (current epoch values)
    if( !fReader->isMC() )
    {
        double xrot = 0.;
        double yrot = 0.;
        if( getArrayPointing() )
        {
            double iUTC = VSkyCoordinatesUtilities::getUTC( getShowerParameters()->MJD, getShowerParameters()->time );
            getArrayPointing()->derotateCoords( iUTC, getShowerParameters()->fShower_Xoffset[iMethod],
                                                -1.*getShowerParameters()->fShower_Yoffset[iMethod], xrot, yrot );
        }
        
        getShowerParameters()->fShower_XoffsetDeRot[iMethod] = xrot;
        getShowerParameters()->fShower_YoffsetDeRot[iMethod] = -1.*yrot;
    }
    else
    {
        getShowerParameters()->fShower_XoffsetDeRot[iMethod] = getShowerParameters()->fShower_Xoffset[iMethod];
        getShowerParameters()->fShower_YoffsetDeRot[iMethod] = getShowerParameters()->fShower_Yoffset[iMethod];
    }
    
    return true;
}


void VArrayAnalyzer::prepareforDirectionReconstruction( unsigned int iMethodIndex, unsigned iReconstructionMethod )
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::prepareforDirectionReconstruction; preparing method " << iMethodIndex << endl;
    }
    
    double iPointingErrorX = 0.;
    double iPointingErrorY = 0.;
    float i_cen_x = 0.;
    float i_cen_y = 0.;
    float i_phi = 0.;
    
    // reset data vectors
    telID.clear();
    x.clear();
    y.clear();
    l.clear();
    m.clear();
    phi.clear();
    sinphi.clear();
    cosphi.clear();
    w.clear();
    length.clear();
    width.clear();
    asym.clear();
    loss.clear();
    dist.clear();
    pedvar.clear();
    teltype.clear();
    tgrad.clear();
    xcore.clear();
    ycore.clear();
    Xoff.clear();
    Yoff.clear();
    cen_x.clear();
    cen_y.clear();
    ltrig.clear();
    ze.clear();
    az.clear();
    
    ///////////////////////////////////////////////
    // fill the x, y, w, and m arrays for the fit
    for( unsigned int tel = 0; tel < getNTel(); tel++ )
    {
        setTelID( tel );
        if( getShowerParameters()->fTelIDImageSelected[iMethodIndex][tel]
                && getEvndispReconstructionParameter( iMethodIndex ) )
        {
            telID.push_back( tel );
            // get pointing difference between expected pointing towards source and measured pointing
            // (by command line, tracking program or pointing monitors)
            if( !getEvndispReconstructionParameter( iMethodIndex )->fUseEventdisplayPointing && tel < getPointing().size() && getPointing()[tel] )
            {
                iPointingErrorX = getPointing()[tel]->getPointingErrorX();
                iPointingErrorY = getPointing()[tel]->getPointingErrorY();
            }
            // do not use pointing corrections
            else
            {
                iPointingErrorX = 0.;
                iPointingErrorY = 0.;
            }
            // get image centroids corrected for pointing errors
            i_cen_x = getImageParameters( getRunParameter()->fImageLL )->cen_x + iPointingErrorX;
            i_cen_y = getImageParameters( getRunParameter()->fImageLL )->cen_y + iPointingErrorY;
            // centroid locations in getImageParameters( getRunParameter()->fImageLL ) are in [deg]
            // (in contrary to centroids in GrIsu ([mm]))
            if( iReconstructionMethod == 4 || iReconstructionMethod == 5 )
            {
                x.push_back( i_cen_x );
                y.push_back( i_cen_y );
            }
            else
            {
                x.push_back( tan( i_cen_x * TMath::DegToRad() )*getDetectorGeo()->getFocalLength()[tel] * 1000. );
                y.push_back( tan( i_cen_y * TMath::DegToRad() )*getDetectorGeo()->getFocalLength()[tel] * 1000. );
            }
            // weight is size
            w.push_back( getImageParameters( getRunParameter()->fImageLL )->size );
            // calculate new 'phi' with pointing errors taken into account
            i_phi = recalculateImagePhi( iPointingErrorX, iPointingErrorY );
            if( cos( i_phi ) != 0. )
            {
                m.push_back( sin( i_phi ) / cos( i_phi ) );
            }
            else
            {
                m.push_back( 1.e9 );
            }
            
            phi.push_back( i_phi );
            sinphi.push_back( sin( i_phi ) );
            cosphi.push_back( cos( i_phi ) );
            if( getImageParameters( getRunParameter()->fImageLL )->length > 0. )
            {
                l.push_back( getImageParameters( getRunParameter()->fImageLL )->width /
                             getImageParameters( getRunParameter()->fImageLL )->length );
            }
            else
            {
                l.push_back( 1. );
            }
            length.push_back( getImageParameters( getRunParameter()->fImageLL )->length );
            width.push_back( getImageParameters( getRunParameter()->fImageLL )->width );
            loss.push_back( getImageParameters( getRunParameter()->fImageLL )->loss );
            asym.push_back( getImageParameters( getRunParameter()->fImageLL )->asymmetry );
            dist.push_back( getImageParameters( getRunParameter()->fImageLL )->dist );
            pedvar.push_back( getImageParameters( getRunParameter()->fImageLL )->fmeanPedvar_Image );
            tgrad.push_back( getImageParameters( getRunParameter()->fImageLL )->tgrad_x );
            cen_x.push_back( getImageParameters( getRunParameter()->fImageLL )->cen_x );
            cen_y.push_back( getImageParameters( getRunParameter()->fImageLL )->cen_y );
            xcore.push_back( getShowerParameters()->fShowerXcore[0] );
            ycore.push_back( getShowerParameters()->fShowerYcore[0] );
            Xoff.push_back( getShowerParameters()->fShower_Xoffset[0] );
            Yoff.push_back( getShowerParameters()->fShower_Xoffset[0] );
            ltrig.push_back( getShowerParameters()->fLTrig );
            ze.push_back( 90. - getShowerParameters()->fTelElevation[tel] );
            az.push_back( getShowerParameters()->fTelAzimuth[tel] );
            teltype.push_back( getDetectorGeometry()->getTelType()[tel] );
        }
    }
}


void VArrayAnalyzer::prepareforCoreReconstruction( unsigned int iMethodIndex, float xs, float ys )
{
    vector< float > xtel;                         /* arrays for telescope locations */
    vector< float > ytel;
    vector< float > ztel;
    xtelnew.clear();
    ytelnew.clear();
    ztelnew.clear();
    for( unsigned int tel = 0; tel < getNTel(); tel++ )
    {
        setTelID( tel );
        xtel.push_back( getImageParameters( getRunParameter()->fImageLL )->Tel_x_SC );
        xtelnew.push_back( 0. );
        ytel.push_back( getImageParameters( getRunParameter()->fImageLL )->Tel_y_SC );
        ytelnew.push_back( 0. );
        ztel.push_back( getImageParameters( getRunParameter()->fImageLL )->Tel_z_SC );
        ztelnew.push_back( 0. );
    }
    rcs_rotate_delta( xtel, ytel, ztel, xtelnew, ytelnew, ztelnew, xs / TMath::RadToDeg(), ys / TMath::RadToDeg(), getNTel() );
    
    /* fill the x, y, w, and m arrays  */
    x.clear();
    y.clear();
    m.clear();
    w.clear();
    float i_cen_x = 0.;
    float i_cen_y = 0.;
    float i_weight = 0.;
    for( unsigned int tel = 0; tel < getNTel(); tel++ )
    {
        setTelID( tel );
        if( getShowerParameters()->fTelIDImageSelected[iMethodIndex][tel] )
        {
            x.push_back( xtelnew[tel] );          /* telescope locations */
            y.push_back( ytelnew[tel] );
            i_weight  = getImageParameters( getRunParameter()->fImageLL )->size;
            i_weight *= ( 1. - getImageParameters( getRunParameter()->fImageLL )->width /
                          getImageParameters( getRunParameter()->fImageLL )->length );
            w.push_back( i_weight * i_weight );
            i_cen_x = getImageParameters( getRunParameter()->fImageLL )->cen_x - xs;
            i_cen_y = getImageParameters( getRunParameter()->fImageLL )->cen_y - ys;
            m.push_back( -1.*i_cen_y / i_cen_x );
        }
    }
}


/*!
      reconstruction of shower direction and core

      modified disp method

      method 5: MLP disp analysis
      method 6: table based disp analysis
      method 7: TMVA BDT based anlaysis

      todo: core construction
*/
int VArrayAnalyzer::rcs_method_5( unsigned int iMethod, unsigned int iDisp )
{
    if( fDebug )
    {
        cout << "VArrayAnalyzer::rcs_method_5/6/7 " << iMethod << ", DISP " << iDisp << endl;
    }
    
    // check that disp analyzer exists
    if( iMethod >= fDispAnalyzer.size() || !fDispAnalyzer[iMethod] || fDispAnalyzer[iMethod]->isZombie() )
    {
        fillShowerDirection( iMethod, -9999., -9999., -2. );
        return 0;
    }
    
    // number of images included in the fit
    int num_images = 0;
    
    // impact point, from perpen. minimization
    float xs = 0.;
    float ys = 0.;
    
    // dispdiff gamma/hadron parameter
    float dispdiff = 0.;
    
    float ximp = 0.;
    float yimp = 0.;
    float stdp = 0.;
    
    num_images = getShowerParameters()->fShowerNumImages[iMethod];
    
    if( !getEvndispReconstructionParameter( iMethod ) )
    {
        return 0;
    }
    
    ////////////////////////////////////////////////////////////////////
    // run a geometrical reconstruction
    // (mean value of pairs of intersecting lines)
    rcs_method_4( iMethod );
    
    float xoff_4 = getShowerParameters()->fShower_Xoffset[iMethod];
    float yoff_4 = getShowerParameters()->fShower_Yoffset[iMethod];
    
    // are there enough images the run an array analysis
    if( num_images >= ( int )getEvndispReconstructionParameter( iMethod )->fNImages_min )
    {
        prepareforDirectionReconstruction( iMethod, 5 );
    }
    else
    {
        fillShowerDirection( iMethod, -9999., -9999., -1. );
        return 0;
    }
    
    // core direction
    
    
    /////////////////////////////////////////////////////
    // direction reconstruction with MLP/TMVA or disp tables
    ////////////////////////////////////////////////////
    prepareforDirectionReconstruction( iMethod, 5 );
    
    vector< float > v_disp;
    vector< float > v_weight;
    
    // loop over all valid images
    for( unsigned int ii = 0; ii < m.size(); ii++ )
    {
        // calculate displacement from image values
        float disp = fDispAnalyzer[iMethod]->evaluate( width[ii], length[ii], asym[ii], dist[ii], w[ii], pedvar[ii], tgrad[ii],
                loss[ii], cen_x[ii], cen_y[ii], xoff_4, yoff_4,
                teltype[ii], ze[ii], az[ii], true );
        v_disp.push_back( disp );
        
        // weigth for averaging
        if( fDispAnalyzer[iMethod]->getDispE() > 0. )
        {
            v_weight.push_back( 1. / fDispAnalyzer[iMethod]->getDispE() );
        }
        else
        {
            setTelID( telID[ii] );
            v_weight.push_back( getImageParameters( getRunParameter()->fImageLL )->ntubes *
                                ( 1. - getImageParameters( getRunParameter()->fImageLL )->width
                                  / getImageParameters( getRunParameter()->fImageLL )->length ) );
        }
    }
    
    // calculate average direction
    // (this is the direction used later in the analysis)
    fDispAnalyzer[iMethod]->calculateMeanDirection( xs, ys, x, y, cosphi, sinphi, v_disp, v_weight, dispdiff );
    
    // dispdiff is a measure how close the directions reconstructed from the individual
    // images are
    getShowerParameters()->fDispDiff[iMethod] = dispdiff;
    
    /*
     * calculate shower directions for each image
     */
    for( unsigned int ii = 0; ii < m.size(); ii++ )
    {
        getShowerParameters()->addDISPPoint( telID[ii], iMethod,
                                             fDispAnalyzer[iMethod]->getXcoordinate_disp( ii, x[ii], cosphi[ii] ),
                                             fDispAnalyzer[iMethod]->getYcoordinate_disp( ii, y[ii], sinphi[ii] ), 1. );
    }
    if( !fillShowerDirection( iMethod, xs, ys, 0. ) )
    {
        return 0;
    }
    
    ////////////////////////////////////////////////
    
    // core reconstruction (still to be replaced by DISP method)
    prepareforCoreReconstruction( iMethod, xs, ys );
    
    // Now call perpendicular_distance for the fit, returning ximp and yimp
    rcs_perpendicular_fit( x, y, w, m, num_images, &ximp, &yimp, &stdp );
    
    // convert xs,ys into [deg]
    fillShowerCore( iMethod, ximp, yimp );
    
    getShowerParameters()->fShower_stdP[iMethod] = stdp;
    /* indicates fit ok */
    if( xs > -9998. && ys > -9998. )
    {
        getShowerParameters()->fShower_Chi2[iMethod] = 0.0;
    }
    
    return 0;
}


// get Monte Carlo shower parameters and fill them into showerpars tree
bool VArrayAnalyzer::fillSimulationEvent( bool iHasArrayTrigger )
{
#ifndef NOVBF
    // ignore pedestal events
    // (assume that vbf files only have pedestal events)
    if( fReader->getATEventType() == VEventType::PED_TRIGGER )
    {
        return false;
    }
#endif
    
    if( fReader->isMC() )
    {
        if( getMCParameters() )
        {
            // fill MC parameter tree
            getMCParameters()->runNumber = getRunNumber();
            getMCParameters()->eventNumber = getEventNumber();
            getMCParameters()->MCprimary = fReader->getMC_primary();
            getMCParameters()->MCenergy  = fReader->getMC_energy();
            getMCParameters()->MCxcore   = fReader->getMC_X();
            getMCParameters()->MCycore   = fReader->getMC_Y();
            getMCParameters()->MCxcos    = fReader->getMC_Xcos();
            getMCParameters()->MCycos    = fReader->getMC_Ycos();
            getMCParameters()->MCTel_Xoff = fReader->getMC_Xoffset();
            getMCParameters()->MCTel_Yoff = fReader->getMC_Yoffset();
            getMCParameters()->MCze = fReader->getMC_Ze();
            getMCParameters()->MCaz = fReader->getMC_Az();
            getMCParameters()->MCCorsikaRunID = fReader->getMCCorsikaRunID();
            getMCParameters()->MCCorsikaShowerID = fReader->getMCCorsikaShowerID();
            getMCParameters()->MCFirstInteractionHeight = fReader->getMCFirstInteractionHeight();
            getMCParameters()->MCFirstInteractionDepth = fReader->getMCFirstInteractionDepth();
            getMCParameters()->ArrayTrigger = ( short int )iHasArrayTrigger;
            getMCParameters()->fill();
        }
    }
    return true;
}



/*!
      reconstruction of shower direction and core

      combination of method 4 and 6

      weighted by f(cos(ze)) (after M.Beilicke 2010)

*/
int VArrayAnalyzer::rcs_method_8( unsigned int iMethod )
{
    double cosze = 0.;
    if( getArrayPointing() )
    {
        cosze = cos( ( 90. - getArrayPointing()->getTelElevation() ) / TMath::RadToDeg() );
    }
    // weight calculation according to Beilicke 2010
    double weight = 0.;
    if( cosze < 0.4 )
    {
        weight = 1.;
    }
    else
    {
        weight = TMath::Exp( -12.5 * ( cosze - 0.4 ) * ( cosze - 0.4 ) );
    }
    
    // first call intersection of lines
    rcs_method_4( iMethod );
    
    float xoff_4 = getShowerParameters()->fShower_Xoffset[iMethod];
    float yoff_4 = getShowerParameters()->fShower_Yoffset[iMethod];
    float stds_4 = getShowerParameters()->fShower_stdS[iMethod];
    
    //  disp method
    rcs_method_5( iMethod, 6 );
    
    float xoff_6 = getShowerParameters()->fShower_Xoffset[iMethod];
    float yoff_6 = getShowerParameters()->fShower_Yoffset[iMethod];
    
    // calculate weighted mean between the two methods
    float xoff_7 = xoff_4 * ( 1. - weight ) + xoff_6 * weight;
    float yoff_7 = -1.* ( yoff_4 * ( 1. - weight ) + yoff_6 * weight );
    
    if( !fillShowerDirection( iMethod, xoff_7, yoff_7, stds_4 ) )
    {
        return 0;
    }
    
    // (core reconstruction: same in method 4 and 7)
    
    return 0;
}


/*

   pass telescope pointing to star catalogue for calcuation of
   distance to bright stars

*/
bool VArrayAnalyzer::updatePointingToStarCatalogue( unsigned int iTelescope )
{
    if( !getStarCatalogue() )
    {
        return false;
    }
    
    // get pointing of telescope
    float iTel_dec = 0.;
    float iTel_ra  = 0.;
    if( iTelescope < getPointing().size() && getPointing()[iTelescope] )
    {
        iTel_dec = getPointing()[iTelescope]->getTelDec() * TMath::RadToDeg();
        iTel_ra  = getPointing()[iTelescope]->getTelRA() * TMath::RadToDeg();
    }
    else if( getArrayPointing() )
    {
        iTel_dec = getArrayPointing()->getTelDec() * TMath::RadToDeg();
        iTel_ra  = getArrayPointing()->getTelRA() * TMath::RadToDeg();
    }
    else
    {
        return false;
    }
    double iScale = 1.;
    if( iTelescope < getDetectorGeometry()->getCameraScaleFactor().size() )
    {
        iScale = getDetectorGeometry()->getCameraScaleFactor()[iTelescope];
    }
    
    getStarCatalogue()->setTelescopePointing( iTelescope,
            getArrayPointing()->getDerotationAngle( getEventMJD(), getEventTime() ) * TMath::RadToDeg(),
            iTel_ra, iTel_dec, iScale );
            
    return true;
}

// for getting the interpolated pointing at a specific MJD and Time
// for filling the pointingDataReduced TTree in the evndisp.<>.root file
void VArrayAnalyzer::updatePointingToArbitraryTime( int iMJD, double iTime )
{

    // same conditions as above
    if( !fReader->isMC() &&
            !( getRunParameter()->felevation > 0.0 ) &&
            !( getRunParameter()->fazimuth > 0.0 ) &&
            getArrayPointing()->isSet() &&
            ( getArrayPointing()->getTargetName() != "laser" ) )
    {
    
        getArrayPointing()->updatePointing( iMJD, iTime );
        if( !getArrayPointing()->isPrecessed() )
        {
            getArrayPointing()->precessTarget( iMJD, -1 );
            // set wobble offsets
            getArrayPointing()->setWobbleOffset( getRunParameter()->fWobbleNorth, getRunParameter()->fWobbleEast, -1, iMJD );
            getArrayPointing()->updatePointing( iMJD, iTime );
        }
        
    }
    
}

/*
 * return TMVA file for angular reconstruction
 *
 * return file with zenith angle closest to the simulated ones
 *
 */
string VArrayAnalyzer::getTMVAFileNameForAngularReconstruction( unsigned int iStereoCutCounter, string iBDTFileName )
{
    string iName = "";
    
    // use cos(ze) to determine the closest ze bin
    double iAverageZenith = cos( ( 90. - getAverageElevation() ) * TMath::DegToRad() );
    
    // loop over all zenith angles and files given in the reconstructionparameter file
    if( getEvndispReconstructionParameter( iStereoCutCounter ) )
    {
        unsigned int iBinSelected = 0;
        double iDiff = 1.e20;
        for( unsigned int i = 0; i < getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAZenithBin.size(); i++ )
        {
            double iZe = cos( getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAZenithBin[i] * TMath::DegToRad() );
            if( TMath::Abs( iZe - iAverageZenith ) < iDiff )
            {
                iDiff = TMath::Abs( iZe - iAverageZenith );
                iBinSelected = i;
            }
        }
        // this is the closest bin
        if( iBinSelected < getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAFileNameVector.size() )
        {
            iName = getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAFileNameVector[iBinSelected];
        }
    }
    //////////////
    // check if full file name is given, otherwise adjust it
    
    // full path given
    if( iName.find( "/" ) != string::npos )
    {
        iName += "/" + iBDTFileName;
    }
    // no full path given - expect file to be at the default location
    // $OBS_EVNDISP_AUX_DIR/DISP_BDTs/VX/ze...
    else
    {
    
        // With a data file, fInstrumentEpoch is "V<epoch>", but with a
        // simulation file, fInstrumentEpoch is "<epoch>", so we have to
        // only get the epoch number as a string
        string epostr = "V" ;
        if( getRunParameter()->fInstrumentEpoch.find( "V" ) != string::npos )
        {
            // if the fInstrumentEpoch has the format "V<epoch>", just use it as the epostr
            epostr = getRunParameter()->fInstrumentEpoch ;
        }
        else
        {
            // if the fInstrumentEpoch has the format "<epoch>", add it to the existing "V"
            epostr.append( getRunParameter()->fInstrumentEpoch ) ;
        }
        
        cout << "scanned fInstrumentEpoch : converted '" << getRunParameter()->fInstrumentEpoch << "' to '" << epostr << "' ..." << endl;
        
        //string iFullFileName = getRunParameter()->getDirectory_EVNDISPAnaData() + "/DISP_BDTs/V" + getRunParameter()->fInstrumentEpoch + "/";
        string iFullFileName = getRunParameter()->getDirectory_EVNDISPAnaData() + "/DISP_BDTs/" + epostr + "/";
        
        iName = iFullFileName + iName + "/" + iBDTFileName;
    }
    cout << "initializing TMVA disp analyzer for average zenith angle in current run: " << 90. - getAverageElevation() << endl;
    if( fDebug )
    {
        cout << "TMVA BDTs available for: " << endl;
        if( getEvndispReconstructionParameter( iStereoCutCounter ) )
        {
            for( unsigned int i = 0; i < getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAZenithBin.size(); i++ )
            {
                cout << "\t" << getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAZenithBin[i] << endl;
            }
        }
    }
    
    return iName;
}

void VArrayAnalyzer::initializeDispAnalyzer( unsigned int iStereoCutCounter )
{
    // first check if we can 'reuse' a MVA from another, previously defined method
    //   this is indicated by a single TMVABDTFILE BDT file with the file name USE_BDT_METHOD_<methodID>
    if( getEvndispReconstructionParameter( iStereoCutCounter )
            && getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAFileNameVector.size() == 1
            && getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAFileNameVector[0].find( "USE_BDT_METHOD_" ) != string::npos )
    {
        string iTemp = getEvndispReconstructionParameter( iStereoCutCounter )->fDISP_TMVAFileNameVector[0];
        unsigned int iMethodID = ( unsigned int )atoi( iTemp.substr( iTemp.rfind( "_" ) + 1, iTemp.size() ).c_str() );
        cout << "initializing TMVA disp analyzer for array reconstruction method " << iStereoCutCounter;
        cout << " using disp analyser from method " << iMethodID << endl;
        if( iMethodID < fDispAnalyzer.size() && fDispAnalyzer[iMethodID] )
        {
            fDispAnalyzer.push_back( fDispAnalyzer[iMethodID] );
        }
        else
        {
            cout << "VArrayAnalyzer::initAnalysis() error initializing MVA-BDT (method " << iStereoCutCounter << ")" << endl;
            cout << "\t could not find disp analyzer from previous method " << iMethodID << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    // initialize disp analyzer
    else
    {
        fDispAnalyzer.push_back( new VDispAnalyzer() );
        fDispAnalyzer.back()->setTelescopeTypeList( getDetectorGeometry()->getTelType_list() );
        // initialize disp with BDTs from closest zenith angle
        if( !fDispAnalyzer.back()->initialize( getTMVAFileNameForAngularReconstruction( iStereoCutCounter ), "TMVABDT" ) )
        {
            cout << "VArrayAnalyzer::initAnalysis() error initializing MVA-BDT (method " << iStereoCutCounter << ")" << endl;
            cout << "\t file " << getTMVAFileNameForAngularReconstruction( iStereoCutCounter ) << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
}

/*
   test angle between two image lines

   return false if angle is smaller then the angle given in the reconstruction parameters

*/
bool VArrayAnalyzer::testMinAngleBetweenImages( unsigned int iMethod, float iangdiff )
{
    if( !getEvndispReconstructionParameter( iMethod ) || iMethod >= VDST_MAXRECMETHODS )
    {
        getShowerParameters()->fShower_Chi2[iMethod] = -9999.;
        return false;
    }
    getShowerParameters()->fiangdiff[iMethod] = iangdiff * TMath::RadToDeg();
    if( iangdiff < getEvndispReconstructionParameter( iMethod ) ->fAxesAngles_min / TMath::RadToDeg() ||
            fabs( 180. / TMath::RadToDeg() - iangdiff ) < getEvndispReconstructionParameter( iMethod )->fAxesAngles_min / TMath::RadToDeg() )
    {
        getShowerParameters()->fShower_Chi2[iMethod] = -1.*fabs( atan( m[0] ) - atan( m[1] ) ) * TMath::RadToDeg();
        return false;
    }
    return true;
}
