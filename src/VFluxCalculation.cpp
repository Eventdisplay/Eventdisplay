/*! \class VFluxCalculation
    \brief calculate integral flux for a given time bin / run

    Questions:

    - get rid of run-wise fluxes?
    - change the way the effective areas are filled in anasum, right now there are filled if events pass gamma/hadron cuts
    - check how time masks are taken into account
    	- calculation of mean time of bin
	- exposure calculation
    - use VLightCurveBinner in calculateCombinedFluxes

    Limitations:

    - hLinerecCounts histograms define min/max energy: currently in bins of 10 GeV. Change to 5 GeV to allow more energy thresholds?

*/

#include "VFluxCalculation.h"

/*
   emtpy constructor
*/

VFluxCalculation::VFluxCalculation()
{
    reset();
    
    // no files are loaded
    bZombie = true;
}

/*

   calculate fluxes using a single anasum result file

*/
VFluxCalculation::VFluxCalculation( string iDataFile,
                                    int iRunMin, int iRunMax,
                                    double iMJDMin, double iMJDMax,
                                    double iFluxMultiplier,
                                    bool iDebug )
{
    reset();
    fDebug = iDebug;
    
    if( iDataFile.find( ".root" ) < iDataFile.size() )
    {
        bZombie = openAnasumDataFile( iDataFile );
        loadRunListFromAnasumDataFile( iRunMin, iRunMax, iMJDMin, iMJDMax );
    }
    else
    {
        bZombie = loadFluxDataVectorFromAsciiFile( iDataFile, iFluxMultiplier, iMJDMin, iMJDMax );
    }
    
}


/*
   calculate fluxes using multiple anasum result files
*/
VFluxCalculation::VFluxCalculation( vector< string > iAnasumFile_vector,
                                    int iRunMin, int iRunMax,
                                    double iMJDMin, double iMJDMax,
                                    bool iDebug )
{
    reset();
    fDebug = iDebug;
    
    bZombie = openAnasumDataFile( iAnasumFile_vector );
    loadRunListFromAnasumDataFile( iRunMin, iRunMax, iMJDMin, iMJDMax );
}

VFluxCalculation::~VFluxCalculation()
{
    closeFiles();
}

/*

   close all anasum result files

*/
void VFluxCalculation::closeFiles()
{
    for( unsigned int i = 0; i < fFile.size(); i++ )
    {
        if( fFile[i] && !fFile[i]->IsZombie() )
        {
            fFile[i]->Close();
        }
    }
}


void VFluxCalculation::reset()
{
    fFile.clear();
    fData = 0;
    
    clearTimeBinVectors();
    clearPhaseBinVectors();
    
    // long (debug) printout
    fDebug = false;
    
    // default maximum usable energy [TeV]
    setMaxSaveMC_Energy_TeV();
    // assumed spectral parameters
    setSpectralParameters();
    
    // significance parameters: decide when to calculate fluxes and when upper limits
    setSignificanceParameters( 3., 5., 0.99, 5, 17, true );
    
    fCalculateExpectedLimitsN = 0;
    
    fUseRunWiseBins = true;
    fUseIntraRunBins = false;
}

void VFluxCalculation::clearTimeBinVectors()
{
    fRequestedTimeBins_MJD_start.clear();
    fRequestedTimeBins_MJD_stopp.clear();
}

void VFluxCalculation::clearPhaseBinVectors()
{
    fRequestedPhaseBins_start.clear();
    fRequestedPhaseBins_stopp.clear();
}

/*
   clear all data vectors
*/
void VFluxCalculation::resetRunList()
{
    fFluxData_perRun.clear();
    fFluxData_perTimeBin.clear();
    fFluxDataVector.clear();
    fFluxDataCombined.reset();
}


/*
    get a run list from tree tRunSummary in anasum output file

    get run duration, zenith angle, excess events from this tree (per run)

    apply cuts on run number and/or MJD range
*/
unsigned int VFluxCalculation::loadRunListFromAnasumDataFile( int iRunMin, int iRunMax, double iMJDMin, double iMJDMax )
{
    if( fDebug )
    {
        cout << "VFluxCalculation::loadRunListFromAnasumDataFile" << endl;
    }
    if( bZombie )
    {
        return 0;
    }
    
    resetRunList();
    
    /////////////////////////////////////////////////////////////////
    // loop over all files and chain the run summary trees together
    if( !fData )
    {
        char hname[2000];
        TChain* c = new TChain( "tRunSummary" );
        for( unsigned int i = 0; i < fFile.size(); i++ )
        {
            sprintf( hname, "%s/total_1/stereo/tRunSummary", fFile[i]->GetName() );
            c->Add( hname );
            if( fDebug )
            {
                cout << "\t chaining " << hname << endl;
            }
        }
        fData = new CRunSummary( c );
    }
    
    int nentries = fData->fChain->GetEntries();
    if( fDebug )
    {
        cout << "total number of runs: " << nentries - 1;
        cout << " (run number interval [" << iRunMin << ", " << iRunMax << "])" << endl;
    }
    
    ///////////////////////////////////////////////////////
    // loop over all run summary trees and read data
    for( int i = 0; i < nentries; i++ )
    {
        fData->GetEntry( i );
        
        // skip combined data element (is calculated by combining all individual data elements)
        if( fData->runOn < 0 )
        {
            continue;
        }
        // check run min/max requirements
        if( fData->runOn > 0 && fData->runOn < iRunMin )
        {
            continue;
        }
        if( fData->runOn > 0 && fData->runOn > iRunMax )
        {
            continue;
        }
        // check MJD min/max requirements
        if( iMJDMin > 0. && fData->MJDOn < iMJDMin && fData->runOn != -1 )
        {
            continue;
        }
        if( iMJDMax > 0. && fData->MJDOn > iMJDMax && fData->runOn != -1 )
        {
            continue;
        }
        
        // skip runs with zero elevation (no events)
        if( fData->elevationOn < 5. && fData->elevationOff < 5 )
        {
            if( fDebug )
            {
                cout << "\t skip run " << fData->runOn << ": no entries" << endl;
            }
            continue;
        }
        
        /////////////////////////////////////////////////////
        // new flux data element
        VFluxDataPoint i_FluxData;
        
        i_FluxData.fRunNumber    = fData->runOn;
        i_FluxData.fMJD_RunStart = fData->MJDOn_runStart;
        i_FluxData.fMJD_RunStop  = fData->MJDOn_runStart;
        i_FluxData.fMJD          = fData->MJDOn;
        i_FluxData.fMJD_Start    = fData->MJDOn_runStart;
        i_FluxData.fMJD_Stop     = fData->MJDOn_runStopp;
        i_FluxData.fTimeBinDuration_sec = fData->RunDurationOn;
        i_FluxData.fSecondsIntoRun = 0.5 * fData->RunDurationOn;
        
        i_FluxData.fExposure   = fData->tOn;
        i_FluxData.fDeadTime   = fData->DeadTimeFracOn;
        i_FluxData.fExposure_deadTimeCorrected = i_FluxData.fExposure * ( 1. - i_FluxData.fDeadTime );
        
        // mean zenith angle for this run
        if( fData->elevationOn > 1. )
        {
            i_FluxData.fZe = 90. - fData->elevationOn;
        }
        else
        {
            i_FluxData.fZe = 90. - fData->elevationOff;
        }
        i_FluxData.fWobbleOffset = sqrt( fData->WobbleNorth * fData->WobbleNorth + fData->WobbleWest * fData->WobbleWest );
        i_FluxData.fPedvars      = fData->pedvarsOn;
        i_FluxData.fNon          = fData->NOn;
        i_FluxData.fNoff         = fData->NOff;
        i_FluxData.fAlpha        = fData->OffNorm;
        
        fFluxData_perRun.push_back( i_FluxData );
    }
    
    if( fFluxData_perRun.size() == 0 )
    {
        cout << "run list empty - no flux calculation / analysis possible" << endl;
        return 0;
    }
    
    if( fDebug )
    {
        cout << "total number of runs in runlist: ";
        cout << fFluxData_perRun.size() << endl;
    }
    
    return fFluxData_perRun.size();
}

void VFluxCalculation::printRunList( ostream& output )
{
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        fFluxData_perRun[i].printRunData( output );
    }
}


bool VFluxCalculation::openAnasumDataFile( string ifile )
{
    vector< string > iTemp;
    iTemp.push_back( ifile );
    return openAnasumDataFile( iTemp );
}

/*
   open a list of root data files

   return true: failure and at least one file is a Zombie
*/
bool VFluxCalculation::openAnasumDataFile( vector< string > ifile )
{
    fFile.clear();
    for( unsigned int i = 0; i < ifile.size(); i++ )
    {
        fFile.push_back( new TFile( ifile[i].c_str() ) );
        if( fFile.back()->IsZombie() )
        {
            return true;
        }
    }
    
    return false;
}

/*

   calculate integral flux for each run, each time bin and `for the merged data vector

*/
bool VFluxCalculation::calculateIntegralFlux( double iMinEnergy_TeV )
{
    fMinEnergy_TeV = iMinEnergy_TeV;
    
    // get number of events above threshold energy (and define time dependent data vectors)
    if( !getNumberOfEventsinEnergyInterval() )
    {
        cout << "VFluxCalculation::calculateIntegralFlux error calculating events above energy " << fMinEnergy_TeV << " TeV" << endl;
        return false;
    }
    
    // set spectral and orbital parameters for each time bin
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        fFluxData_perRun[i].setSpectralParameters( fMinEnergy_TeV, fE0_TeV, fSpectralIndex, fMaxEnergy_TeV );
        fFluxData_perRun[i].calculateOrbitalPhaseData( fOrbitalPhaseData );
    }
    for( unsigned int t = 0; t < fFluxData_perTimeBin.size(); t++ )
    {
        fFluxData_perTimeBin[t].setSpectralParameters( fMinEnergy_TeV, fE0_TeV, fSpectralIndex, fMaxEnergy_TeV );
        fFluxData_perTimeBin[t].calculateOrbitalPhaseData( fOrbitalPhaseData );
    }
    
    // calculate significances and upper limits
    if( !calculateSignificancesAndUpperLimits() )
    {
        cout << "VFluxCalculation::calculateIntegralFlux error calculating significances and upper limits" << endl;
        return false;
    }
    
    // calculate integrated spectral weighted effective areas
    if( !getIntegralEffectiveArea() )
    {
        cout << "VFluxCalculation::calculateIntegralFlux error calculating spectral weighted effective areas" << endl;
        return false;
    }
    
    // calculate fluxes and upper flux limits
    if( !calculateFluxes() )
    {
        cout << "VFluxCalculation::calculateIntegralFlux error calculating fluxes" << endl;
        return false;
    }
    
    // combine / merge flux points to final results
    if( !calculateCombinedFluxes() )
    {
        cout << "VFluxCalculation::calculateIntegralFlux error calculating combined fluxes (combined in time bins)" << endl;
        return false;
    }
    
    return true;
}

/*

   merge/combine fluxes to calculate flux point in the requested time bins

   also calculate total flux here

   treatment of time bins should probably be moved to an extra class (VLightCurveBinner)

*/
bool VFluxCalculation::calculateCombinedFluxes()
{
    if( fDebug )
    {
        cout << "calculating combined fluxes" << endl;
        cout << "------------------" << endl;
    }
    
    
    //////////////////////////////////////////////////
    // time bins / run wise analysis results
    fFluxDataVector.clear();
    
    //////////////////////////////////////////////////////////////////////////
    // combined results for all runs
    // (use run-wise results; should be equivalent to summing up time bins)
    fFluxDataCombined.reset();
    fFluxDataCombined.setCombinedDataElement();
    fFluxDataCombined.setSignificanceParameters( fThresholdSignificance, fMinEvents, fUpperLimit, fUpperLimitMethod, fLiMaEqu, fBoundedLimits );
    fFluxDataCombined.setSpectralParameters( fMinEnergy_TeV, fE0_TeV, fSpectralIndex, fMaxEnergy_TeV );
    
    //////////////////////////////////////////////////
    // sort internal data vectors
    sort( fFluxData_perRun.begin(), fFluxData_perRun.end() );
    sort( fFluxData_perTimeBin.begin(), fFluxData_perTimeBin.end() );
    
    //////////////////////////////////////////////////
    // first check if a vector of periods is given
    if( fRequestedTimeBins_MJD_start.size() > 0 )
    {
        // check that number of start and stop times are consistent
        if( fRequestedTimeBins_MJD_start.size() != fRequestedTimeBins_MJD_stopp.size() )
        {
            cout << "VFluxCalculation::calculateCombinedFluxes() error: MJD start and stop vectors do not match" << endl;
            cout << "\t start: " << fRequestedTimeBins_MJD_start.size() << ", stop: " << fRequestedTimeBins_MJD_stopp.size() << endl;
            return false;
        }
        // loop over all periods
        for( unsigned int i = 0; i < fRequestedTimeBins_MJD_start.size(); i++ )
        {
            bool i_new_period = true;
            for( unsigned int j = 0; j < fFluxData_perTimeBin.size(); j++ )
            {
                // check if MJD is in given period
                double iMiddleOfTimeBin = ( fFluxData_perTimeBin[j].fMJD_Start + fFluxData_perTimeBin[j].fMJD_Stop ) * 0.5;
                if( iMiddleOfTimeBin >= fRequestedTimeBins_MJD_start[i]
                        && iMiddleOfTimeBin < fRequestedTimeBins_MJD_stopp[i] )
                {
                    // period-wise data vector
                    if( i_new_period )
                    {
                        VFluxDataPoint iTempDataElement;
                        fFluxDataVector.push_back( iTempDataElement );
                        fFluxDataVector.back().setSignificanceParameters( fThresholdSignificance, fMinEvents, fUpperLimit,
                                fUpperLimitMethod, fLiMaEqu, fBoundedLimits );
                        fFluxDataVector.back().setSpectralParameters( fMinEnergy_TeV, fE0_TeV, fSpectralIndex, fMaxEnergy_TeV );
                        i_new_period = false;
                    }
                    if( fFluxDataVector.size() > 0 )
                    {
                        fFluxDataVector.back() = fFluxDataVector.back() + fFluxData_perTimeBin[j];
                        // combined data vector
                        fFluxDataCombined = fFluxDataCombined + fFluxData_perTimeBin[j];
                    }
                }
            }
        }
        // sort analysis results
        sort( fFluxDataVector.begin(), fFluxDataVector.end() );
        
        
    }
    
    ///////////////////////////////////
    // Runwise binning
    else if( fUseRunWiseBins )
    {
        if( fDebug )
        {
            cout << "Doing run-wise analysis." << endl;
        }
        for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
        {
            if( fDebug )
            {
                cout << "Run " << i << endl;
            }
            
            // run-wise data vector
            fFluxDataVector.push_back( fFluxData_perRun[i] );
            // combined data vector
            fFluxDataCombined = fFluxDataCombined + fFluxData_perRun[i];
        }
        // sort analysis results
        sort( fFluxDataVector.begin(), fFluxDataVector.end() );
        
        
    }
    ///////////////////////////////////
    // Intra-run binning
    else if( fUseIntraRunBins )
    {
        if( fDebug )
        {
            cout << "Keeping anasum time bins." << endl;
        }
        for( unsigned int i = 0; i < fFluxData_perTimeBin.size(); i++ )
        {
            // small time-bin wise data vector
            fFluxDataVector.push_back( fFluxData_perTimeBin[i] );
            // combined data vector
            fFluxDataCombined = fFluxDataCombined + fFluxData_perTimeBin[i];
        }
        // sort analysis results
        sort( fFluxDataVector.begin(), fFluxDataVector.end() );
        
        
    }
    
    //phase binning
    else if( fRequestedPhaseBins_start.size() > 0 )
    {
        if( fRequestedPhaseBins_start.size() != fRequestedPhaseBins_stopp.size() )
        {
            cout << "VFluxCalculation::calculateCombinedFluxes() error: phase start and stop vectors do not match" << endl;
            cout << "\t start: " << fRequestedPhaseBins_start.size() << ", stop: " << fRequestedPhaseBins_stopp.size() << endl;
            return false;
        }
        // loop over all periods
        for( unsigned int i = 0; i < fRequestedPhaseBins_start.size(); i++ )
        {
            bool i_new_period = true;
            for( unsigned int j = 0; j < fFluxData_perTimeBin.size(); j++ )
            {
                // check if phase is in given period
                // note: phase bins may wrap around 0, eg. the user might specify a bin from 0.9 to 0.1.
                if(
                    ( fRequestedPhaseBins_start[i] < fRequestedPhaseBins_stopp[i] &&  fFluxData_perTimeBin[j].fOrbitalPhase >= fRequestedPhaseBins_start[i] && fFluxData_perTimeBin[j].fOrbitalPhase < fRequestedPhaseBins_stopp[i] )
                    || ( fRequestedPhaseBins_start[i] > fRequestedPhaseBins_stopp[i] && ( fFluxData_perTimeBin[j].fOrbitalPhase >= fRequestedPhaseBins_start[i] || fFluxData_perTimeBin[j].fOrbitalPhase < fRequestedPhaseBins_stopp[i] ) ) )
                {
                    // period-wise data vector
                    if( i_new_period )
                    {
                        VFluxDataPoint iTempDataElement;
                        fFluxDataVector.push_back( iTempDataElement );
                        fFluxDataVector.back().setSignificanceParameters( fThresholdSignificance, fMinEvents, fUpperLimit,
                                fUpperLimitMethod, fLiMaEqu, fBoundedLimits );
                        fFluxDataVector.back().setSpectralParameters( fMinEnergy_TeV, fE0_TeV, fSpectralIndex, fMaxEnergy_TeV );
                        i_new_period = false;
                    }
                    if( fFluxDataVector.size() > 0 )
                    {
                        fFluxDataVector.back() = fFluxDataVector.back() + fFluxData_perTimeBin[j];
                        // combined data vector
                        fFluxDataCombined = fFluxDataCombined + fFluxData_perTimeBin[j];
                    }
                }
            }
        }
        
    }
    
    else
    {
        cout << "VFluxCalculation::calculateCombinedFluxes() error: no method of combining fluxes was specified." << endl;
        return false;
    }
    
    
    return true;
}

void VFluxCalculation::calculateExpectedLimitCombinedOnly( int n )
{
    //if(!fFluxDataCombined) return;
    fFluxDataCombined.setCalculateExpectedLimitsN( n );
    fFluxDataCombined.calculateExpectedLimit();
    fFluxDataCombined.calculateFlux();
}



/*!
    read effective areas from anasum output file for each run

    calculate spectral weighted effective area average

    < A (E) > = INT_{E}^{INF} A(E) (E/E_0)^{-fSpectralIndex} dE

*/
bool VFluxCalculation::getIntegralEffectiveArea()
{
    if( fDebug )
    {
        cout << "calculating spectral weighted integral effective areas" << endl;
        cout << "------------------------------------------------------" << endl;
        cout << endl;
    }
    
    /////////////////////////
    // loop over all files
    for( unsigned int f = 0; f < fFile.size(); f++ )
    {
        if( bZombie || !fFile[f] )
        {
            return false;
        }
        
        fFile[f]->cd();
        if( fDebug )
        {
            cout << "now at file " << fFile[f]->GetName() << " (runlist size: " << fFluxData_perRun.size() << ")" << endl;
        }
        char hname[800];
        
        ///////////////////////////////////
        // loop over all runs in this file
        for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
        {
            if( fDebug )
            {
                cout << "Calculating spectral weighted integral effective area for run " << fFluxData_perRun[i].fRunNumber;
                cout << " at ze [deg]: " << fFluxData_perRun[i].fZe << endl;
            }
            
            sprintf( hname, "run_%d/stereo/EffectiveAreas", fFluxData_perRun[i].fRunNumber );
            if( !fFile[f]->Get( hname ) )
            {
                continue;
            }
            if( !fFile[f]->cd( hname ) )
            {
                cout << "directory " << hname << " not found" << endl;
                cout << "continue..." << endl;
                continue;
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////
            // graphs with mean effective areas (mean for given run; filled in anasum)
            /////////////////////////////////////////////////////////////////////////////////////////////////
            TGraphAsymmErrors* g = ( TGraphAsymmErrors* )gDirectory->Get( "gMeanEffectiveArea_on" );
            /////////////////////////////////////////////////////////////////////////////////////////////////
            // found graph with effective areas
            if( !g )
            {
                // try go get off graph
                g = ( TGraphAsymmErrors* )gDirectory->Get( "gMeanEffectiveArea_off" );
                if( !g )
                {
                    cout << "VFluxCalculation::getIntegralEffectiveArea(): error: effective area graph not found" << endl;
                    cout << "continue..." << endl;
                    fFluxData_perRun[i].fEffArea_cm2 = 0.;
                    continue;
                }
            }
            
            // cp graph in vectors
            vector< double > iV_energy;
            vector< double > iV_effArea;
            double i_x = 0.;
            double i_y = 0.;
            for( int p = 0; p < g->GetN(); p++ )
            {
                g->GetPoint( p, i_x, i_y );
                iV_energy.push_back( i_x );
                iV_effArea.push_back( i_y );
            }
            fFluxData_perRun[i].calculateIntegralEffectiveArea( iV_energy, iV_effArea );
            
            /////////////////////////////////////////////////////////////////////
            // TIME dependent effective areas
            /////////////////////////////////////////////////////////////////////
            // read effective area graphs for time binned intra run light curves
            TGraph2DErrors* g_time = ( TGraph2DErrors* )gDirectory->Get( "gTimeBinnedMeanEffectiveArea" );
            // find 2D graph with time dependent effective areas
            if( !g_time )
            {
                // if not try go get off graph
                g_time = ( TGraph2DErrors* )gDirectory->Get( "gTimeBinnedMeanEffectiveArea_off" );
                if( !g_time )
                {
                    cout << "VFluxCalculation::getIntegralEffectiveArea(): error: 2D effective area graph not found" << endl;
                    cout << "continue..." << endl;
                    continue;
                }
            }
            // calculate spectral weighted integral effective area
            int     iTB_bins        = g_time->GetN();
            double* iTB_energy_axis = g_time->GetX();
            double* iTB_time_axis   = g_time->GetZ();
            double* iTB_effArea     = g_time->GetY();
            
            // loop over all time bins
            for( unsigned t = 0; t < fFluxData_perTimeBin.size(); t++ )
            {
                // select time bins for current run and ignore all other runs
                if( fFluxData_perTimeBin[t].fRunNumber == fFluxData_perRun[i].fRunNumber )
                {
                    // reset temporary vectors for effective area calculation
                    iV_energy.clear();
                    iV_effArea.clear();
                    // assume effective areas per time bin are sorted in energy
                    for( int b = 0; b < iTB_bins; b++ )
                    {
                        if( fFluxData_perTimeBin[t].isTimeInsideRun( iTB_time_axis[b] ) )
                        {
                            iV_energy.push_back( iTB_energy_axis[b] );
                            iV_effArea.push_back( iTB_effArea[b] );
                        }
                    }
                    // the last bin might be zero due to time cuts - take second last bin
                    // QQQQ
                    // (this happens as effective area value are filled in VStereoAnalysis only
                    //  when there was an event after gamma / hadron separation cuts)
                    // QQQQ need to discuss a better solution
                    if( iV_energy.size() == 0 && t > 0 )
                    {
                        fFluxData_perTimeBin[t].fEffArea_cm2 = fFluxData_perTimeBin[t - 1].fEffArea_cm2;
                    }
                    // everything went well in the reading of the effective areas
                    else
                    {
                        fFluxData_perTimeBin[t].calculateIntegralEffectiveArea( iV_energy, iV_effArea );
                    }
                }
            }
        }
    }
    
    return true;
}

/*

   print results for each data point

*/
void VFluxCalculation::printResultSummary( ostream& output )
{
    for( unsigned int t = 0; t < fFluxDataVector.size(); t++ )
    {
        fFluxDataVector[t].printResultSummary( output );
    }
    output << "All data combined:" << endl;
    fFluxDataCombined.printResultSummary( output );
    output << endl;
}

/*

   print results for each data point

*/
void VFluxCalculation::printResults( ostream& output )
{
    for( unsigned int t = 0; t < fFluxDataVector.size(); t++ )
    {
        fFluxDataVector[t].printResults( output );
    }
    output << endl << "All data combined:" << endl;
    fFluxDataCombined.printResults( output );
    output << endl;
}

/*

   calculate integral flux for each run and all runs together

   use spectral weighted effective areas from getIntegralEffectiveArea()

*/
bool VFluxCalculation::calculateFluxes()
{
    if( fDebug )
    {
        cout << "calculating fluxes" << endl;
        cout << "------------------" << endl;
    }
    
    // loop over all runs in run list
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        if( fDebug )
        {
            cout << "Run " << i << endl;
        }
        fFluxData_perRun[i].calculateFlux();
    }
    
    // Time binned analysis
    for( unsigned int t = 0; t < fFluxData_perTimeBin.size(); t++ )
    {
        if( fDebug )
        {
            cout << "Bin " << t << endl;
        }
        fFluxData_perTimeBin[t].calculateFlux();
    }
    return true;
}

/*

    set spectral parameters

    these values are used to calculate a spectral weighted effective area

*/
void VFluxCalculation::setSpectralParameters( double iMinEnergy_TeV, double E0, double alpha , double iMaxEnergy_TeV )
{
    fMinEnergy_TeV = iMinEnergy_TeV;
    if( iMaxEnergy_TeV <  fMaxSave_MCEnergy_TeV && iMaxEnergy_TeV > 0. )
    {
        fMaxEnergy_TeV = iMaxEnergy_TeV;
    }
    else
    {
        fMaxEnergy_TeV = fMaxSave_MCEnergy_TeV;
    }
    fE0_TeV = E0;
    fSpectralIndex = alpha;
}

/*

   setting significance parameters

   (used for the decision of when to use fluxes and when upper flux limits)

*/
void VFluxCalculation::setSignificanceParameters( double iThresholdSignificance, double iMinEvents, double iUpperLimit,
        int iUpperlimitMethod, int iLiMaEqu, bool iBoundedLimits )
{
    fThresholdSignificance = iThresholdSignificance;
    fMinEvents = iMinEvents;
    fUpperLimit = iUpperLimit;
    fUpperLimitMethod = iUpperlimitMethod;
    fLiMaEqu = iLiMaEqu;
    fBoundedLimits = iBoundedLimits;
}

/*

     get number of on and off events above a given energy threshold [TeV]

*/
bool VFluxCalculation::getNumberOfEventsinEnergyInterval()
{
    if( fDebug )
    {
        cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval: " << fMinEnergy_TeV << " < E [TeV] < " << fMaxEnergy_TeV << endl;
    }
    
    /////////////////////////////
    // run (file) wise analysis
    for( unsigned int f = 0; f < fFile.size(); f++ )
    {
        if( !fFile[f] || fFile[f]->IsZombie() || !fFile[f]->cd() )
        {
            cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval: error: no input file found" << endl;
            return false;
        }
        if( !getNumberOfEventsinEnergyInterval( fFile[f] ) )
        {
            return false;
        }
        
        fFile[f]->cd();
    }      // end loop over all files
    return true;
}

/*

     get number of on and off events above a given energy threshold [TeV]

     note: also initialize the time-binned data vector

*/
bool VFluxCalculation::getNumberOfEventsinEnergyInterval( TFile* iFile )
{
    if( !iFile )
    {
        return false;
    }
    
    char hname[200];
    fFluxData_perTimeBin.clear();
    
    ///////////////////////////////////////////////////
    // loop over all runs in run list
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        if( fDebug )
        {
            cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval in run " << fFluxData_perRun[i].fRunNumber << endl;
        }
        sprintf( hname, "run_%d/stereo/energyHistograms", fFluxData_perRun[i].fRunNumber );
        
        // check if this run is in the current file
        if( !iFile->Get( hname ) )
        {
            continue;
        }
        if( !iFile->cd( hname ) )
        {
            cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval: error finding directory " << hname;
            cout << " in file " << iFile->GetName() << endl;
            return false;
        }
        
        //////////////////////////////////////////////////////
        // run wise analysis
        //////////////////////////////////////////////////////
        TH1D* hon  = ( TH1D* )gDirectory->Get( "hLinerecCounts_on" );
        TH1D* hoff = ( TH1D* )gDirectory->Get( "hLinerecCounts_off" );
        if( !hon || !hoff )
        {
            cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval error finding counting histograms (energy): ";
            cout << hon << "\t" << hoff << "\t" << fFluxData_perRun[i].fRunNumber << endl;
            return false;
        }
        // get number of on and off events above the given energy threshold
        fFluxData_perRun[i].fNon = 0.;
        fFluxData_perRun[i].fNoff = 0.;
        for( int b = hon->GetNbinsX(); b > 0; b-- )
        {
            if( hon->GetXaxis()->GetBinLowEdge( b ) < fMinEnergy_TeV )
            {
                // adjust energy range to binned values (upper edge of bin)
                fMinEnergy_TeV = hon->GetXaxis()->GetBinUpEdge( b );
                break;
            }
            if( fMaxEnergy_TeV < fMaxSave_MCEnergy_TeV && hon->GetXaxis()->GetBinUpEdge( b ) > fMaxEnergy_TeV )
            {
                // (don't think that adjusting the energy to the binning is important for the upper cut value)
                continue;
            }
            fFluxData_perRun[i].fNon  += hon->GetBinContent( b );
            fFluxData_perRun[i].fNoff += hoff->GetBinContent( b );
        }
        
        // clean up
        gDirectory->RecursiveRemove( hon );
        gDirectory->RecursiveRemove( hoff );
        
        // end run-wise analysis
        //////////////////////////
        
        ////////////////////////////////////////////////////////////////////////////////
        // time binned analysis
        // (time bins shorter than the length of a run)
        //
        // also: declaration of vector data element for time binned analysis
        //
        ////////////////////////////////////////////////////////////////////////////////
        // get number of on events above energy threshold in Time BIN
        // get time binned histograms from anasum file
        // (expect that all histograms have some y-axis binning)
        TH2D* hon2DtimeBinned          = ( TH2D* )gDirectory->Get( "hLinerecCounts2DtimeBinned_on" );
        TH2D* hoff2DtimeBinned         = ( TH2D* )gDirectory->Get( "hLinerecCounts2DtimeBinned_off" );
        TH1D* honRealDuration1DtimeBinned  = ( TH1D* )gDirectory->Get( "hRealDuration1DtimeBinned_on" );
        TH1D* honDuration1DtimeBinned  = ( TH1D* )gDirectory->Get( "hDuration1DtimeBinned_on" );
        if( !hon2DtimeBinned || !hoff2DtimeBinned || !honDuration1DtimeBinned || !honRealDuration1DtimeBinned )
        {
            cout << "VFluxCalculation::getNumberOfEventsinEnergyInterval: error finding time binned 2D counting histogram (energy): ";
            cout << hon2DtimeBinned << "\t" << hoff2DtimeBinned << "\t" << honDuration1DtimeBinned << "\t" << honRealDuration1DtimeBinned;
            cout << "\t" << fFluxData_perRun[i].fRunNumber << endl;
            return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // set up time binnned vector fFluxData_perTimeBin[runnumber][time bin]
        // one flux data vector per bin in time binned histogram
        
        // loop over all time bins and
        //   1. set up data vector for time binned analysis
        //   2. get number of on/off events above energy threshold
        for( int t = 1; t <= hon2DtimeBinned->GetNbinsY(); t++ )
        {
            // initialize a new data element with the values from the current run
            fFluxData_perTimeBin.push_back( fFluxData_perRun[i] );
            fFluxData_perTimeBin.back().resetFluxValues();
            fFluxData_perTimeBin.back().fNon = 0.;
            fFluxData_perTimeBin.back().fNoff = 0.;
            // time bins
            fFluxData_perTimeBin.back().fTimeBinDuration_sec = hon2DtimeBinned->GetYaxis()->GetBinWidth( t );
            fFluxData_perTimeBin.back().fMJD_Start = fFluxData_perRun[i].fMJD_Start + ( t - 1 ) * fFluxData_perTimeBin.back().fTimeBinDuration_sec / 86400.;
            fFluxData_perTimeBin.back().fMJD_Stop  = fFluxData_perTimeBin.back().fMJD_Start + fFluxData_perTimeBin.back().fTimeBinDuration_sec / 86400.;
            // QQQQ time mask should be taken into account??? -> bins within the time mask should be removed
            fFluxData_perTimeBin.back().fMJD = 0.5 * ( fFluxData_perTimeBin.back().fMJD_Start + fFluxData_perTimeBin.back().fMJD_Stop );
            fFluxData_perTimeBin.back().fSecondsIntoRun = ( fFluxData_perTimeBin.back().fMJD - fFluxData_perTimeBin.back().fMJD_RunStart ) * 86400.;
            // exposure = length of time bin, time mask taken into account
            fFluxData_perTimeBin.back().fExposure     = honDuration1DtimeBinned->GetBinContent( t );
            // dead time corrected exposure
            fFluxData_perTimeBin.back().fExposure_deadTimeCorrected = honRealDuration1DtimeBinned->GetBinContent( t );
            // dead time
            fFluxData_perTimeBin.back().fDeadTime = 1.0 - fFluxData_perTimeBin.back().fExposure_deadTimeCorrected / fFluxData_perTimeBin.back().fExposure;
            // check if bin is completely within an open time mask (0 exposure)
            if( fFluxData_perTimeBin.back().fExposure == 0 )
            {
                fFluxData_perTimeBin.back().fTimeMask_open = true;
            }
            
            // 2. get number of events above energy threshold
            for( int b = hon2DtimeBinned->GetNbinsX(); b > 0; b-- )
            {
                if( hon2DtimeBinned->GetXaxis()->GetBinLowEdge( b ) < fMinEnergy_TeV )
                {
                    break;
                }
                if( fMaxEnergy_TeV < fMaxSave_MCEnergy_TeV && hon2DtimeBinned->GetXaxis()->GetBinUpEdge( b ) > fMaxEnergy_TeV )
                {
                    continue;
                }
                fFluxData_perTimeBin.back().fNon  += hon2DtimeBinned->GetBinContent( hon2DtimeBinned->GetBin( b, t ) );
                fFluxData_perTimeBin.back().fNoff += hoff2DtimeBinned->GetBinContent( hoff2DtimeBinned->GetBin( b, t ) );
            }
        } // loop over all time bins
        // end of time-binned analysis
    } // END: loop over all runs
    
    return true;
}

/*

   calculate significance and upper limits (in event numbers) for each run and time bin

*/
bool VFluxCalculation::calculateSignificancesAndUpperLimits()
{
    if( fDebug )
    {
        cout << "calculating significances and upper flux limits" << endl;
        cout << "-----------------------------------------------" << endl;
        cout << "using Li & Ma equation " << fLiMaEqu << endl;
        cout << "using upper flux calculation after ";
        if( fUpperLimitMethod == 0 )
        {
            cout << "Helene";
        }
        else if( fUpperLimitMethod == 3 )
        {
            cout << "Feldman & Cousins";
        }
        else if( fUpperLimitMethod == 4 || fUpperLimitMethod == 5 )
        {
            cout << "Rolke et al";
            if( fBoundedLimits )
            {
                cout << " (bounded limits)";
            }
            else
            {
                cout << " (unbounded limits)";
            }
        }
        cout << ", " << fUpperLimit * 100. << "\% upper limit" << endl;
        cout << "threshold significance for upper limit calculation is " << fThresholdSignificance << " sigma or more than " << fMinEvents << " events" << endl;
        cout << endl;
    }
    
    ////////////////////////////////////////////////////////
    // loop over all runs
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        // calculate significances and upper limits
        fFluxData_perRun[i].setSignificanceParameters( fThresholdSignificance, fMinEvents, fUpperLimit, fUpperLimitMethod, fLiMaEqu, fBoundedLimits );
        fFluxData_perRun[i].setCalculateExpectedLimitsN( fCalculateExpectedLimitsN );
        fFluxData_perRun[i].calculateSignificancesAndUpperLimits();
        if( fDebug )
        {
            cout << "VFluxCalculation::calculateSignificancesAndUpperLimits: " << fFluxData_perRun[i].fRunNumber << ": ";
            cout << fFluxData_perRun[i].fSignificance << " (" << fFluxData_perRun[i].isSignificantDataPoint() << ")" << endl;
        }
    }
    // loop over all times bins
    for( unsigned int t = 0; t < fFluxData_perTimeBin.size(); t++ )
    {
        fFluxData_perTimeBin[t].setSignificanceParameters( fThresholdSignificance, fMinEvents, fUpperLimit, fUpperLimitMethod, fLiMaEqu, fBoundedLimits );
        fFluxData_perTimeBin[t].setCalculateExpectedLimitsN( fCalculateExpectedLimitsN );
        fFluxData_perTimeBin[t].calculateSignificancesAndUpperLimits();
    }
    
    return true;
}


/*

    read time bin vector from ascii file

    Format:

    <MJD_min_1> <MJD_max_1>
    <MJD_min_2> <MJD_max_2>
    <MJD_min_3> <MJD_max_3>
    ...

    line with a '#' in the first column are treated as a comment

*/

bool VFluxCalculation::setTimeBinVector( string iTimeBinFile, bool iPrint )
{
    // open ascii files with time bins
    ifstream is( iTimeBinFile.c_str() );
    if( !is )
    {
        cout << "VFluxCalculation::setTimeBinVector(): error reading " << iTimeBinFile  << endl;
        return false;
    }
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    cout << "VFluxCalculation::setTimeBinVector(): reading time bin vector from " << iTimeBinFile << endl;
    
    string iTemp_string = "";
    
    string is_line;
    
    double iTemp_d = 0.;
    
    //////////////////////
    // loop over all lines
    while( getline( is, is_line ) )
    {
        // empty line
        if( is_line.size() == 0 )
        {
            continue;
        }
        // this is a comment line
        if( is_line.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> iTemp_d;
        fRequestedTimeBins_MJD_start.push_back( iTemp_d );
        is_stream >> iTemp_d;
        fRequestedTimeBins_MJD_stopp.push_back( iTemp_d );
    }
    is.close();
    
    // sort time vectors
    sort( fRequestedTimeBins_MJD_start.begin(), fRequestedTimeBins_MJD_start.end() );
    sort( fRequestedTimeBins_MJD_stopp.begin(), fRequestedTimeBins_MJD_stopp.end() );
    
    if( iPrint )
    {
        for( unsigned int i = 0; i < fRequestedTimeBins_MJD_start.size(); i++ )
        {
            cout << "\t Bin " << i;
            cout << ", MJD: [" << fRequestedTimeBins_MJD_start[i];
            cout << ", " << fRequestedTimeBins_MJD_stopp[i] << "]" << endl;
        }
    }
    
    return true;
}

/*

   set a vector of periods to be used for flux calculation

*/
void VFluxCalculation::setTimeBinVector( vector< double > iRequestedTimeBins_MJD_start, vector< double > iRequestedTimeBins_MJD_stopp )
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    fRequestedTimeBins_MJD_start = iRequestedTimeBins_MJD_start;
    sort( fRequestedTimeBins_MJD_start.begin(), fRequestedTimeBins_MJD_start.end() );
    fRequestedTimeBins_MJD_stopp = iRequestedTimeBins_MJD_stopp;
    sort( fRequestedTimeBins_MJD_stopp.begin(), fRequestedTimeBins_MJD_stopp.end() );
}


void VFluxCalculation::setPhaseBinVector( vector< double > iRequestedPhaseBins_start, vector< double > iRequestedPhaseBins_stopp )
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    fRequestedPhaseBins_start = iRequestedPhaseBins_start;
    fRequestedPhaseBins_stopp = iRequestedPhaseBins_stopp;
}

/*
set orbital phase data. Note: You must call this before calling "calculateIntegralFlux" otherwise the phases won't be right.
*/
void VFluxCalculation::setPhaseFoldingValues( double iZeroPhase_MJD, double iOrbit_Days, double iOrbitError_low_Days, double iOrbitError_up_Days )
{
    fOrbitalPhaseData.fZeroPhase_MJD = iZeroPhase_MJD;
    fOrbitalPhaseData.fOrbit_days = iOrbit_Days;
    fOrbitalPhaseData.fOrbit_days_error_low  = iOrbitError_low_Days;
    fOrbitalPhaseData.fOrbit_days_error_high = iOrbitError_up_Days;
    
    // (re-)calculate orbital phases
    // (not clear if at this point the data vector is already filled)
    // set spectral and orbital parameters for each time bin
    //     for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    //     {
    //	fFluxDataVector[i].calculateOrbitalPhaseData( fOrbitalPhaseData );
    //     }
    for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
    {
        fFluxData_perRun[i].calculateOrbitalPhaseData( fOrbitalPhaseData );
    }
    for( unsigned int i = 0; i < fFluxData_perTimeBin.size(); i++ )
    {
        fFluxData_perTimeBin[i].calculateOrbitalPhaseData( fOrbitalPhaseData );
    }
}

/*

    read light curve data from ASCII file for the given period in MJD

*/
unsigned int VFluxCalculation::loadFluxDataVectorFromAsciiFile( string iAsciiFile,
        double iFluxMultiplier,
        double iMJDMin, double iMJDMax )
{
    if( fDebug )
    {
        cout << "VFluxCalculation::loadFluxDataVectorFromAsciiFile (file: ";
        cout << iAsciiFile << ")" << endl;
    }
    resetRunList();
    
    // open ascii file
    ifstream is( iAsciiFile.c_str() );
    if( !is )
    {
        cout << "VFluxCalculation::loadFluxDataVectorFromAsciiFile(): error reading " << iAsciiFile << endl;
        bZombie = true;
        return false;
    }
    cout << "VFluxCalculation::loadFluxDataVectorFromAsciiFile(): reading " << iAsciiFile << endl;
    
    string iTemp_string = "";
    string iFluxUnitString;
    
    string is_line;
    
    // column definitions
    vector< string > iColumn;
    
    //////////////////////
    // loop over all lines
    //   (take differences in file format into account)
    while( getline( is, is_line ) )
    {
        if( is_line.size() == 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        // first column can be:
        // 1. comment (a '!' or a '#' in the first column)
        // 2. ENERGYRANGE_TEV (energy range is given in [TeV]
        // 3. ENERGYRANGE_KEV (energy range is given in [keV]
        // 4. COLUMN .. ... (column definitions)
        // 5. second since fXRTMissionTimeStart or MJD
        is_stream >> iTemp_string;
        if( iTemp_string.size() == 0
                || iTemp_string.substr( 0, 1 ) == "!"
                || iTemp_string.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        else if( iTemp_string == "ENERGYRANGE_TEV" || iTemp_string == "ENERGYRANGE_KEV" )
        {
            if( !is_stream.eof() )
            {
                is_stream >> fMinEnergy_TeV;
            }
            if( !is_stream.eof() )
            {
                is_stream >> fMaxEnergy_TeV;
            }
            if( iTemp_string == "ENERGYRANGE_KEV" )
            {
                fMinEnergy_TeV *= 1.e-6;
                fMaxEnergy_TeV *= 1.e-6;
            }
            cout << "\t energy range read from ASCII file [" << fMinEnergy_TeV << ", " << fMaxEnergy_TeV << "]" << endl;
            continue;
        }
        else if( iTemp_string == "FLUXUNITS_cm2_s" )
        {
            is_stream >> iTemp_string;
            iFluxMultiplier *= atof( iTemp_string.c_str() );
            continue;
        }
        else if( iTemp_string == "FLUXUNITSTRING" )
        {
            is_stream >> iFluxUnitString;
            iFluxUnitString = is_line.substr( is_stream.tellg(), is_line.size() );
            continue;
        }
        else if( iTemp_string == "COLUMN" )
        {
            do
            {
                is_stream >> iTemp_string;
                // remove all spaces from string
                iTemp_string.erase( std::remove( iTemp_string.begin(), iTemp_string.end(), ' ' ), iTemp_string.end() );
                
                iColumn.push_back( iTemp_string );
                
            }
            while( !is_stream.eof() );
            
            cout << "\t found the following data columns: ";
            for( unsigned int i = 0; i < iColumn.size(); i++ )
            {
                if( iColumn[i] != "IGNORE" )
                {
                    cout << iColumn[i];
                    if( i != iColumn.size() - 1 )
                    {
                        cout << ", ";
                    }
                }
            }
            cout << endl;
            continue;
        }
        
        // new flux data point
        VFluxDataPoint i_FluxData;
        i_FluxData.setSpectralParameters( fMinEnergy_TeV, 1, -2.5, fMaxEnergy_TeV );
        i_FluxData.calculateOrbitalPhaseData( fOrbitalPhaseData );
        i_FluxData.fFluxUnitString = iFluxUnitString;
        
        // reset string stream
        is_stream.clear();
        is_stream.seekg( 0, ios::beg );
        
        // read in all columns
        for( unsigned int i = 0; i < iColumn.size(); i++ )
        {
            if( !is_stream.eof() )
            {
                is_stream >> iTemp_string;
            }
            else
            {
                continue;
            }
            // MJD
            if( iColumn[i] == "MJD" )
            {
                i_FluxData.fMJD = atof( iTemp_string.c_str() );
                i_FluxData.fMJD_Start = i_FluxData.fMJD;
                i_FluxData.fMJD_Stop  = i_FluxData.fMJD;
            }
            // width of MJD bin (must be a column after MJD)
            else if( iColumn[i] == "MJD_WIDTH" )
            {
                i_FluxData.fMJD_Start = i_FluxData.fMJD - atof( iTemp_string.c_str() );
                i_FluxData.fMJD_Stop  = i_FluxData.fMJD + atof( iTemp_string.c_str() );
            }
            else if( iColumn[i] == "MJDSTART" )
            {
                i_FluxData.fMJD_Start = atof( iTemp_string.c_str() );
                i_FluxData.fMJD = 0.5 * ( i_FluxData.fMJD_Start + i_FluxData.fMJD_Stop );
            }
            // end of MJD interval (expect this to be in a column after MJDSTART)
            else if( iColumn[i] == "MJDSTOP" )
            {
                i_FluxData.fMJD_Stop = atof( iTemp_string.c_str() );
                i_FluxData.fMJD = 0.5 * ( i_FluxData.fMJD_Start + i_FluxData.fMJD_Stop );
            }
            // XRT times (seconds since mission start)
            else if( iColumn[i] == "XRTSECONDS" )
            {
                i_FluxData.fMJD = 54857.09977457897 + atof( iTemp_string.c_str() ) / ( 24.0 * 60.0 * 60.0 );
            }
            // expect this to be after XRTSECONDS
            else if( iColumn[i] == "XRTSECONDS_WIDTH" )
            {
                i_FluxData.fMJD_Start = i_FluxData.fMJD - atof( iTemp_string.c_str() ) / ( 24.0 * 60.0 * 60.0 );
                i_FluxData.fMJD_Stop  = i_FluxData.fMJD + atof( iTemp_string.c_str() ) / ( 24.0 * 60.0 * 60.0 );
            }
            // expect this to be after XRTSECONDS
            else if( iColumn[i] == "XRTSECONDS_WIDTH_UP" )
            {
                i_FluxData.fMJD_Stop = i_FluxData.fMJD + atof( iTemp_string.c_str() ) / ( 24.0 * 60.0 * 60.0 );
            }
            // expect this to be after XRTSECONDS
            else if( iColumn[i] == "XRTSECONDS_WIDTH_LO" )
            {
                i_FluxData.fMJD_Start = i_FluxData.fMJD + atof( iTemp_string.c_str() ) / ( 24.0 * 60.0 * 60.0 );
            }
            // exposure
            else if( iColumn[i] == "OBSTIME_MIN" )
            {
                i_FluxData.fExposure = atof( iTemp_string.c_str() ) * 60.;
            }
            // elevation
            else if( iColumn[i] == "ELEVATION_DEG" )
            {
                i_FluxData.fZe = 90. - atof( iTemp_string.c_str() );
            }
            // Non, Noff, alpha
            else if( iColumn[i] == "NON" )
            {
                i_FluxData.fNon = atof( iTemp_string.c_str() );
            }
            else if( iColumn[i] == "NOFF" )
            {
                i_FluxData.fNoff = atof( iTemp_string.c_str() );
            }
            else if( iColumn[i] == "ALPHA" )
            {
                i_FluxData.fAlpha = atof( iTemp_string.c_str() );
            }
            // flux
            else if( iColumn[i] == "FLUX" )
            {
                i_FluxData.fFlux = atof( iTemp_string.c_str() ) * iFluxMultiplier;
                i_FluxData.fFlux_Rolke = atof( iTemp_string.c_str() ) * iFluxMultiplier;
                
                if( i_FluxData.fFluxUL <= 0. )
                {
                    i_FluxData.fSignificantDataPoint = true;
                }
            }
            // flux error  (must be a column after FLUX)
            // (negative value indicates a upper flux limit)
            else if( iColumn[i] == "FLUXERROR" )
            {
                i_FluxData.fFluxE = atof( iTemp_string.c_str() ) * iFluxMultiplier;
                i_FluxData.fFluxCI_lo_1sigma = i_FluxData.fFluxE;
                i_FluxData.fFluxCI_up_1sigma = i_FluxData.fFluxE;
            }
            else if( iColumn[i] == "FLUXUL" )
            {
                i_FluxData.fFluxUL = atof( iTemp_string.c_str() );
                if( i_FluxData.fFluxUL > 0. )
                {
                    i_FluxData.fFluxUL *= iFluxMultiplier;
                    i_FluxData.fSignificantDataPoint = false;
                }
            }
            else if( iColumn[i] == "SIGNIFICANCE" )
            {
                i_FluxData.fSignificance = atof( iTemp_string.c_str() );
            }
            else if( iColumn[i] == "DATAFLAG" )
            {
                // 0 = good time bin
                if( atoi( iTemp_string.c_str() ) != 0 )
                {
                    i_FluxData.fTimeMask_open = true;
                }
            }
        }
        i_FluxData.calculateOrbitalPhaseData();
        // a good data element has a MJD > 1.
        if( i_FluxData.fMJD > 1. )
        {
            // check MJD range:
            if( iMJDMin > 0. && i_FluxData.fMJD_Start < iMJDMin )
            {
                continue;
            }
            if( iMJDMax > 0. && i_FluxData.fMJD_Stop > iMJDMax )
            {
                continue;
            }
            
            fFluxDataVector.push_back( i_FluxData );
        }
    }
    is.close();
    
    // sort vector
    sort( fFluxDataVector.begin(), fFluxDataVector.end() );
    
    cout << "VLightCurve::loadFluxDataVectorFromAsciiFile() total number of light curve data: " << fFluxDataVector.size() << endl;
    
    return true;
}

const VFluxDataPoint* VFluxCalculation::getFluxDataPerRun( int iRunNumber )
{
    if( iRunNumber < 0 )
    {
        return &fFluxDataCombined;
    }
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        if( fFluxDataVector[i].fRunNumber == iRunNumber )
        {
            return &fFluxDataVector[i];
        }
    }
    return 0;
}

/*

     print width of time bin

     assume that all time bins are of equal length (this is tested)

*/
void VFluxCalculation::printTimeBinWidth()
{
    double iT_Bin = 0;
    
    if( fFluxData_perTimeBin.size() > 0 )
    {
        iT_Bin = fFluxData_perTimeBin[0].fTimeBinDuration_sec;
        for( unsigned int i = 1; i < fFluxData_perTimeBin.size(); i++ )
        {
            // don't expect differences in time bin widths larger than 10 s
            if( TMath::Abs( fFluxData_perTimeBin[i].fTimeBinDuration_sec - iT_Bin ) > 10. )
            {
                cout << "warning: not all time bins of equal length" << endl;
                cout << "\t" << i << "\t" << iT_Bin << "\t" << fFluxData_perTimeBin[i].fTimeBinDuration_sec << endl;
                return;
            }
        }
    }
    else
    {
        cout << "data vector of zero length: fill (calculateIntegralFlux) it first" << endl;
        return;
    }
    cout << "Length of shortest time bins: " << iT_Bin << " s" << endl;
}



void VFluxCalculation::setRunwiseLightCurve()
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseIntraRunBins = false;
    fUseRunWiseBins = true;
}

/*
Set up time bins for daily, weekly etc lightcurves. User can specify a start time, otherwise the MJD of the first run will be used.
*/
void VFluxCalculation::setTimeBinsInDays( double nNights, double iMJD_start )
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    
    if( iMJD_start == 0 )
    {
        if( fFluxData_perTimeBin.size() > 0 )
        {
            iMJD_start = ( int )fFluxData_perTimeBin[0].fMJD_Start; //round down to day before
        }
        else if( fFluxData_perRun.size() > 0 )
        {
            iMJD_start = ( int )fFluxData_perRun[0].fMJD_Start; //round down to day before
        }
        else
        {
            cout << "VFluxCalculation::setUpForNightlyLightCurve Error: No flux data available (per time bin or per run)" << endl;
            return;
        }
    }
    
    if( fDebug )
    {
        cout << "Starting bins at MJD " << iMJD_start << ", width "  << nNights * 24.*60.*60. << " s." << endl;
    }
    
    for( double t = iMJD_start; t < fFluxData_perRun.back().fMJD_Stop; t += nNights )
    {
        fRequestedTimeBins_MJD_start.push_back( t );
        fRequestedTimeBins_MJD_stopp.push_back( t + nNights );
    }
}

/*
Set up time bins in seconds, for intra-run light curves. Note: Time bins have to be a multiple of the anasum time bins. Thus, there is no overlap; the binning will 'restart' at the beginning of the run.
seconds=0 will use the anasum time bins (no combination).
*/

void VFluxCalculation::setTimeBinsInSeconds( double seconds )
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    
    if( seconds > 0 )
    {
        for( unsigned int i = 0; i < fFluxData_perRun.size(); i++ )
        {
            if( fDebug )
            {
                cout << "Starting bins at MJD " << fFluxData_perRun[i].fMJD_Start << ", width "  << fRequestedTimeBinWidth_s << " s." << endl;
            }
            for( double t = fFluxData_perRun[i].fMJD_Start; t < fFluxData_perRun[i].fMJD_Stop; t += fRequestedTimeBinWidth_s / ( 24.*60.*60. ) )
            {
                fRequestedTimeBins_MJD_start.push_back( t );
                fRequestedTimeBins_MJD_stopp.push_back( t + fRequestedTimeBinWidth_s / ( 24.*60.*60. ) );
            }
            // last time bin of run should not extend past end of run to avoid overlap
            fRequestedTimeBins_MJD_stopp.back() = fFluxData_perRun[i].fMJD_Stop;
        }
    }
    else
    {
        fUseIntraRunBins = true;
    }
}

/*
Set bins for phase-folded light curve. User can specify a starting phase. For example, nBins = 5 and iStart = 0.1 will give bin edges 0.1, 0.3, 0.5, 0.7, 0.9. The last bin will wrap around 0.
*/

void VFluxCalculation::setPhaseBins( int nBins, double iStart )
{
    clearTimeBinVectors();
    clearPhaseBinVectors();
    fUseRunWiseBins = false;
    fUseIntraRunBins = false;
    
    iStart -= ( int )iStart; //force iStart between 0 and 1
    for( int i = 0; i < nBins; i++ )
    {
        fRequestedPhaseBins_start.push_back( iStart );
        iStart += 1.0 / nBins;
        iStart -= ( int )iStart; //force iStart between 0 and 1
        fRequestedPhaseBins_stopp.push_back( iStart );
    }
}
