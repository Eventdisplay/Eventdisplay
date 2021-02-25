/*! \class VDLWriter
 *  \brief writes DL2 cut event lists
 *
 */

#include "VDL2Writer.h"

/*!
 *
 */
VDL2Writer::VDL2Writer( VInstrumentResponseFunctionRunParameter* iRunPara, 
                        VGammaHadronCuts* icuts )
{
    fRunPara = iRunPara;
    if( !fRunPara )
    {
        cout << "VDL2Writer: no run parameters given" << endl;
        cout << "...exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // cuts
    fCuts = icuts;
    setIgnoreEnergyReconstructionCuts( fRunPara->fIgnoreEnergyReconstructionQuality );
    setIsotropicArrivalDirections( fRunPara->fIsotropicArrivalDirections );
    setTelescopeTypeCuts( fRunPara->fTelescopeTypeCuts );
    
    ////////////////////////////////////
    // tree with cut results for each event
    fCut_Class = 0;
    fCut_MVA = 0.;
    fEventTreeCuts = new TTree( "fEventTreeCuts", "event cuts" );
    fEventTreeCuts->Branch( "CutClass", &fCut_Class, "Class/I" );
    fEventTreeCuts->Branch( "MVA", &fCut_MVA, "MVA/F" );
}

/////////////////////////////////////////////////////////////////////////////////////

VDL2Writer::~VDL2Writer()
{
}

/*
 *
 *  event loop 
 *
 */
bool VDL2Writer::fill( CData* d, 
                       unsigned int iMethod )
{
    // lots of debug output
    bool bDebugCuts = false;
    
    // do not require successfull energy reconstruction
    if( fIgnoreEnergyReconstruction )
    {
        iMethod = 100;
    }
    
    // reset unique event counter
    //	fUniqueEventCounter.clear();
    int iSuccessfullEventStatistics = 0;
    
    //////////////////////////////////////////////////////////////////
    // print some run information
    cout << endl;
    cout << "DL2 Eventlist filling" << endl;
    if( fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
    {
        cout << "\t ignore first " << fRunPara->fIgnoreFractionOfEvents * 100. << " % of events" << endl;
    }
    cout << endl;
    
    // make sure that all data pointers exist
    if( !d )
    {
        cout << "VDL2Writer::fill error: no data tree" << endl;
        return false;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // reset cut statistics
    fCuts->resetCutStatistics();
    
    ///////////////////////////////////////////////////////
    // get full data set and loop over all entries
    ///////////////////////////////////////////////////////
    Long64_t d_nentries = d->fChain->GetEntries();
    Long64_t i_start = 0;
    if( fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
    {
        i_start = ( Long64_t )( fRunPara->fIgnoreFractionOfEvents * d_nentries );
    }
    cout << "\t total number of data events: " << d_nentries;
    cout << " (start at event " << i_start << ")" << endl;
    
    // loop over all events
    for( Long64_t i = i_start; i < d_nentries; i++ )
    {
        d->GetEntry( i );

        // update cut statistics
        fCuts->newEvent();
        
        if( bDebugCuts )
        {
            cout << "============================== " << endl;
            cout << "EVENT entry number " << i << endl;
        }
        
        // apply MC cuts
        if( bDebugCuts )
        {
            cout << "#0 CUT MC " << fCuts->applyMCXYoffCut( d->MCxoff, d->MCyoff, false ) << endl;
        }
        
        if( !fCuts->applyMCXYoffCut( d->MCxoff, d->MCyoff, true ) )
        {
            fillEventDataTree( VGammaHadronCutsStatistics::eMC_XYoff, -1. );
            continue;
        }
        
        ////////////////////////////////
        // apply general quality and gamma/hadron separation cuts
        
        // apply reconstruction cuts
        if( bDebugCuts )
        {
            cout << "#1 CUT applyInsideFiducialAreaCut ";
            cout << fCuts->applyInsideFiducialAreaCut();
            cout << "\t" << fCuts->applyStereoQualityCuts( iMethod, false, i, true ) << endl;
        }
        
        // apply fiducial area cuts
        if( !fCuts->applyInsideFiducialAreaCut( true ) )
        {
            fillEventDataTree( 2, -1. );
            continue;
        }
        
        // apply reconstruction quality cuts
        if( !fCuts->applyStereoQualityCuts( iMethod, true, i , true ) )
        {
            fillEventDataTree( 3, -1. );
            continue;
        }
        
        // apply telescope type cut (e.g. for CTA simulations)
        if( fTelescopeTypeCutsSet )
        {
            if( bDebugCuts )
            {
                cout << "#2 Cut NTELType " << fCuts->applyTelTypeTest( false ) << endl;
            }
            if( !fCuts->applyTelTypeTest( true ) )
            {
                fillEventDataTree( 4, -1. );
                continue;
            }
        }
        
        //////////////////////////////////////
        // apply direction cut
        //
        // bDirectionCut = false: if direction is inside
        // theta_min and theta_max
        //
        // point source cut; use MC shower direction as reference direction
        bool bDirectionCut = false;
        if( !fIsotropicArrivalDirections )
        {
            if( !fCuts->applyDirectionCuts( true ) )
            {
                bDirectionCut = true;
            }
        }
        // background cut; use (0,0) as reference direction
        // (command line option -d)
        else
        {
            if( !fCuts->applyDirectionCuts( true, 0., 0. ) )
            {
                bDirectionCut = true;
            }
        }
        
        //////////////////////////////////////
        // apply energy reconstruction quality cut
        if( !fIgnoreEnergyReconstruction )
        {
            if( bDebugCuts )
            {
                cout << "#4 EnergyReconstructionQualityCuts";
                cout << fCuts->applyEnergyReconstructionQualityCuts( iMethod ) << endl;
            }
            if( !fCuts->applyEnergyReconstructionQualityCuts( iMethod, true ) )
            {
                fillEventDataTree( 6, -1. );
                continue;
            }
        }
        //////////////////////////////////////
        // apply gamma hadron cuts
        if( bDebugCuts )
        {
            cout << "#3 CUT ISGAMMA " << fCuts->isGamma( i ) << endl;
        }
        if( !fCuts->isGamma( i, true ) )
        {
            fillEventDataTree( 7, fCuts->getTMVA_EvaluationResult() );
            continue;
        }
        if( !bDirectionCut )
        {
            fillEventDataTree( 5, fCuts->getTMVA_EvaluationResult() );
        }
        // remaining events
        else
        {
            fillEventDataTree( 0, fCuts->getTMVA_EvaluationResult() );
        }
        
        // unique event counter
        // (make sure that map doesn't get too big)
        if( !bDirectionCut && iSuccessfullEventStatistics >= 0 )
        {
            iSuccessfullEventStatistics++;
        }
    }  // end of loop
    /////////////////////////////////////////////////////////////////////////////
    fCuts->printCutStatistics();
    if( iSuccessfullEventStatistics < 0 )
    {
        iSuccessfullEventStatistics *= -1;
    }
    cout << "\t total number of events after cuts: " << iSuccessfullEventStatistics << endl;
    
    return true;
}


/* 
 * fill tree with cut numbers and MVA values per event
 *
 * 1    hhEcutTrigger_R (color: 1, marker: 20)
 * 2    hhEcutFiducialArea_R (color: 2, marker: 20)
 * 3    hhEcutStereoQuality_R (color: 3, marker: 20)
 * 4    hhEcutTelType_R (color: 4, marker: 20)
 * 5    hhEcutDirection_R (color: 5, marker: 20)
 * 6    hhEcutEnergyReconstruction_R (color: 6, marker: 20)
 * 7    hhEcutGammaHadron_R (color: 7, marker: 20)
 *
 *  1. Events passing gamma/hadron separation cut and direction cut
 *  fEventTreeCuts->Draw("MVA", "Class==5" );
 *
 *  2. Events passing gamma/hadron separation cut and not direction cut
 *  fEventTreeCuts->Draw("MVA", "Class==0" );
 *
 *  3. Events before gamma/hadron separation cut and before direction cut
 *   fEventTreeCuts->Draw("MVA", "Class==0||Class==7||Class==5", "");
 *
 */
void VDL2Writer::fillEventDataTree( int iCutClass, float iMVA )
{
      if( fEventTreeCuts )
      {
          fCut_Class = iCutClass;
          fCut_MVA = iMVA;
          fEventTreeCuts->Fill();
      }
}

