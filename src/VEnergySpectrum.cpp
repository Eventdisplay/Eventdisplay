/*! class VEnergySpectrum

    analyse and plot spectral energy distributions


*/

#include "VEnergySpectrum.h"

VEnergySpectrum::VEnergySpectrum()
{
    fDebug = 0;
    
    bAsciiDataFile = false;
    
    fDataSetName = "E";
    fTotalRun = -1;
    
    initializeRunVariables();
}

VEnergySpectrum::VEnergySpectrum( string iFile, string iname, int irun, bool iSourceTypeIsAscii )
{
    fDebug = 0;
    bAsciiDataFile = iSourceTypeIsAscii;
    
    fDataSetName = iname;
    fTotalRun = irun;
    
    initializeRunVariables();
    
    openDataFile( iFile, fTotalRun, bAsciiDataFile );
}

bool VEnergySpectrum::openDataFile( string iFile, int irun, bool iSourceTypeIsAscii )
{
    fTotalRun = irun;
    bAsciiDataFile = iSourceTypeIsAscii;
    
    // open anasum file
    if( !bAsciiDataFile )
    {
        if( !openFile( iFile, fTotalRun ) )
        {
            bZombie = true;
            return false;
        }
    }
    else
    {
        if( !openAsciiFile( iFile ) )
        {
            bZombie = true;
            return false;
        }
    }
    return true;
}


void VEnergySpectrum::initializeRunVariables()
{

    nRebinner = 0;
    bUseRebinner = false;
    
    bCombineRuns = false;
    
    bZombie = false;
    
    hErecCountsOn = 0;
    hErecCountsOff = 0;
    hErecTotalTime = 0;
    hErecTotalTimeDeadTimeCorrected = 0;
    hEffArea = 0;
    
    fAnalysisEnergyBinning = -1.;
    
    fEnergyThresholdFixedValue = 0;
    fEnergyThresholdFileName = "";
    
    // mean normalisation factor (alpha)
    // (calculated from total number ofe events)
    fTotalNormalisationFactor = 0.;
    fTotalObservationTime = 0.;
    fTotalObservationTimeDeadTimeCorrected = 0.;
    
    fAnalysisMinEnergy = 0.;
    fAnalysisMaxEnergy = 1.e10;
    
    // default values
    setAddHistogramParameters();
    setEnergyBinning();
    setSignificanceParameters();
    setEnergyThresholdDefinition();
    setErrorCalculationMethod();
    setOffsetdistance();
    
    // set default fitting parameters
    setSpectralFitFunction();
    setSpectralFitFluxNormalisationEnergy();
    setSpectralFitRangeLin();
    setSpectralFitPlottingStyle();
    
    // set default plotting parameters
    setPlottingSpectralWeightForBinCenter();
    setEnergyInBinDefinition();
    fPlottingCanvas = 0;
    setPlottingYaxis();
    setPlottingMultiplierIndex();
    setPlottingEnergyRangeLog();
    setPlottingLogEnergyAxis();
    setPlottingStyle();
    setPlottingUpperLimits();
    
    gEnergySpectrum = 0;
    fEnergySpectrumFit = 0;
    
}


void VEnergySpectrum::setEnergyBinning( double iB )
{
    fAnalysisEnergyBinning = iB;
}

/*

   energies in TeV

*/
void VEnergySpectrum::setEnergyRangeLinear( double xmin, double xmax )
{
    fAnalysisMinEnergy = xmin;
    fAnalysisMaxEnergy = xmax;
}

/*

   energies in log10 TeV

*/
void VEnergySpectrum::setEnergyRangeLog( double xmin, double xmax )
{
    fAnalysisMinEnergy = TMath::Power( 10., xmin );
    fAnalysisMaxEnergy = TMath::Power( 10., xmax );
}

/*

   combine energy spectra from all runs

*/
bool VEnergySpectrum::combineRuns()
{
    vector< int > i_temp;
    return combineRuns( i_temp );
}

/*

   read differential flux points from an ascii file

*/

bool VEnergySpectrum::openAsciiFile( string iFile )
{
    // assume that upper and lower errors are given
    fErrorCalculationMethod = "UPDOWN";
    
    // open ascii file
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VEnergySpectrum::openAsciiFile:: error reading data from " << iFile << endl;
        return false;
    }
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        // comment
        if( is_line.substr( 0, 1 ) == "*" )
        {
            continue;
        }
        
        // new data vector
        VDifferentialFluxData i_flux;
        
        istringstream is_stream( is_line );
        
        is_stream >> i_flux.Energy;
        i_flux.EnergyWeightedMean = i_flux.Energy;
        is_stream >> i_flux.DifferentialFlux;
        is_stream >> i_flux.DifferentialFluxError_low;
        if( i_flux.DifferentialFluxError_low > 0. )
        {
            i_flux.DifferentialFluxError_low = i_flux.DifferentialFlux - i_flux.DifferentialFluxError_low;
        }
        is_stream >> i_flux.DifferentialFluxError_up;
        if( i_flux.DifferentialFluxError_up > 0. )
        {
            i_flux.DifferentialFluxError_up = i_flux.DifferentialFluxError_up - i_flux.DifferentialFlux;
        }
        
        if( i_flux.DifferentialFluxError_low > 0. && i_flux.DifferentialFluxError_low > 0. )
        {
            i_flux.DifferentialFluxError = sqrt( i_flux.DifferentialFluxError_low * i_flux.DifferentialFluxError_low
                                                 + i_flux.DifferentialFluxError_up * i_flux.DifferentialFluxError_up );
        }
        else
        {
            i_flux.DifferentialFluxError = -99.;
        }
        
        // add this bin to data vector
        fDifferentialFlux.push_back( i_flux );
        
    }
    is.close();
    
    return true;
}

/*

   combine energy spectra from runs in the given run list

*/
bool VEnergySpectrum::combineRuns( vector< int > runlist )
{
    if( isZombie() )
    {
        return false;
    }
    
    if( fDebug == 1 )
    {
        cout << "VEnergySpectrum::combineRuns " << runlist.size() << endl;
    }
    
    if( bUseRebinner )
    {
        cout << "VEnergySpectrum::combineRuns: using rebinner" << endl;
    }
    else
    {
        cout << "VEnergySpectrum::combineRuns: not using rebinner" << endl;
    }
    // clear old data vector
    fDifferentialFlux_perRun.clear();
    
    // check if runs are available and read run values from anasum file
    readRunList( runlist, fTotalRun );
    
    // total observation time
    fTotalObservationTime = 0.;
    fTotalObservationTimeDeadTimeCorrected = 0.;
    
    // total normalisation factor
    // (alpha; averaged over all runs and energies)
    fTotalNormalisationFactor = 0.;
    
    // energy threshold calculator (note different definitions of energy thresholds)
    VEnergyThreshold iEnergyThresholdCalculator( fEnergyThresholdFixedValue, fEnergyThresholdFileName );
    
    ///////////////////////////////////////////////////////////////////////
    // loop over all runs in run list
    //
    // calculate fluxes for each run, after that average over all runs
    //
    int z = 0;
    double i_alpha_U = 0.;
    double i_alpha_L = 0.;
    string hname;
    cout << "combining " << fRunList.size() << " runs" << endl;
    if( fRunList.size() == 0 )
    {
        return false;
    }
    double i_nonCounts_counter  = 0.;
    double i_noffCounts_counter = 0.;
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        // differential counting histogram 'on'
        hname = "herecCounts_on";
        if( fOffsetDistance > -998. )
        {
            hname = "herecCounts2D_vs_distance_on";
        }
        TH1D* i_hErecCountsOn = ( TH1D* )getHistogram( hname, fRunList[i].runnumber, "energyHistograms", fOffsetDistance );
        // differential counting histogram 'off'
        hname = "herecCounts_off";
        if( fOffsetDistance > -998. )
        {
            hname = "herecCounts2D_vs_distance_off";
        }
        TH1D* i_hErecCountsOff = ( TH1D* )getHistogram( hname, fRunList[i].runnumber, "energyHistograms", fOffsetDistance );
        // make sure the counting histograms exist
        if( !i_hErecCountsOn || !i_hErecCountsOff )
        {
            cout << "histograms not found for run " << fRunList[i].runnumber << endl;
            continue;
        }
        /////////////////////////////////////////
        // histogram with systematic errors (same for on and off); used for energy threshold calculation
        hname = "gMeanEnergySystematicError";
        TGraphErrors* i_hEsys = ( TGraphErrors* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", fOffsetDistance );
        if( fAnalysisEnergyThresholdDefinition == 1 && !i_hEsys )
        {
            cout << "WARNING: histogram with systematic error in energy reconstruction not found";
            cout << " (run " << fRunList[i].runnumber << ")" << endl;
        }
        /////////////////////////////////////////
        // effective areas
        TGraphErrors* i_gEff = new TGraphErrors( 0 );
        hname = "gMeanEffectiveArea";
        i_gEff = ( TGraphErrors* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", fOffsetDistance );
        
        if( !i_gEff )
        {
            // second choice (if above fails): try to get off effective areas (might even have better statistics)
            hname = "gMeanEffectiveArea_off";
            i_gEff = ( TGraphErrors* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", fOffsetDistance );
            if( !i_gEff )
            {
                cout << "WARNING: no mean effective area graph found, ignoring run ";
                cout << " (run " << fRunList[i].runnumber << ")" << endl;
                continue;
            }
        }
        
        // total number of signal and background event numbers
        i_nonCounts_counter  += i_hErecCountsOn->GetEntries();
        i_noffCounts_counter += i_hErecCountsOff->GetEntries();
        
        /////////////////////////////////////////////////////
        // get energy threshold (may depend on ze, az, ...)
        // (for energy threshold definitions see VEnergyThreshold )
        /////////////////////////////////////////////////////
        // no energy threshold applied
        if( fAnalysisEnergyThresholdDefinition == 0 )
        {
            fRunList[i].energyThreshold = 0.;
        }
        // energy threshold: systematic smaller then given value
        else if( fAnalysisEnergyThresholdDefinition == 1 )
        {
            fRunList[i].energyThreshold = iEnergyThresholdCalculator.getEnergy_maxSystematic( i_hEsys,
                                          fAnalysisMaxEnergySystematic );
        }
        // energy threshold: effective area > given fraction of maximum effective area (typical value is 10%)
        else if( fAnalysisEnergyThresholdDefinition == 2 )
        {
            fRunList[i].energyThreshold = iEnergyThresholdCalculator.getEnergy_MaxEffectiveAreaFraction( i_gEff,
                                          fAnalysisMaxEffectiveAreaFraction );
        }
        // energy threshold given by fixed value (e.g. 0.5 TeV)
        else if( fAnalysisEnergyThresholdDefinition == 3 )
        {
            fRunList[i].energyThreshold = iEnergyThresholdCalculator.getEnergy_fixedValue();
        }
        // make sure that energy threshold is > 0
        if( fRunList[i].energyThreshold < 1.e-10 )
        {
            fRunList[i].energyThreshold = 1.e-10;
        }
        // print run info
        if( fDebug == 1 )
        {
            fRunList[i].print();
        }
        
        /////////////////////////////////////////////////////
        // rebin energy spectra according to user input
        /////////////////////////////////////////////////////
        if( bUseRebinner && nRebinner == 0 && newBinningGroupings.size() > 0 )
        {
            setOriginalBinner( i_hErecCountsOn );
        }
        
        rebinEnergySpectrum( i_hErecCountsOn, fAnalysisEnergyBinning );
        rebinEnergySpectrum( i_hErecCountsOff, fAnalysisEnergyBinning );
        // clone first histogram and make summary histograms
        // (this ensures that binning is correct)
        if( z == 0 )
        {
            // differential counting histogram 'on'
            hname = "herecCounts_on";
            hErecCountsOn = ( TH1D* )i_hErecCountsOn->Clone( hname.c_str() );
            rebinEnergySpectrum( hErecCountsOn, fAnalysisEnergyBinning );
            hErecCountsOn->Reset();
            hErecCountsOn->SetEntries( 0 );
            // differential counting histogram 'off'
            hname = "herecCounts_off";
            hErecCountsOff = ( TH1D* )i_hErecCountsOff->Clone( hname.c_str() );
            rebinEnergySpectrum( hErecCountsOff, fAnalysisEnergyBinning );
            hErecCountsOff->Reset();
            hErecCountsOff->SetEntries( 0 );
            // histogram with total observation time
            hname = "herecTotalTime";
            hErecTotalTime = ( TH1D* )i_hErecCountsOn->Clone( hname.c_str() );
            rebinEnergySpectrum( hErecTotalTime, fAnalysisEnergyBinning );
            hErecTotalTime->Reset();
            hErecTotalTime->SetEntries( 0 );
            hErecTotalTime->SetYTitle( "observing time [s]" );
            // histogram with total observation time
            hname = "hErecTotalTimeDeadTimeCorrected";
            hErecTotalTimeDeadTimeCorrected = ( TH1D* )i_hErecCountsOn->Clone( hname.c_str() );
            rebinEnergySpectrum( hErecTotalTimeDeadTimeCorrected, fAnalysisEnergyBinning );
            hErecTotalTimeDeadTimeCorrected->Reset();
            hErecTotalTimeDeadTimeCorrected->SetEntries( 0 );
            hErecTotalTimeDeadTimeCorrected->SetYTitle( "observing time [s]" );
            // histogram with mean effective area
            hname = "hEffArea";
            hEffArea = ( TH1D* )i_hErecCountsOn->Clone( hname.c_str() );
            rebinEnergySpectrum( hEffArea, fAnalysisEnergyBinning );
            if( hEffArea )
            {
                hEffArea->Reset();
            }
        }
        
        /////////////////////////////////////////////////
        // combine histograms
        /////////////////////////////////////////////////
        
        // calculate total observation time (take the energy threshold into account)
        addValueToHistogram( hErecTotalTime, fRunList[i].tOn, log10( fRunList[i].energyThreshold ) );
        fTotalObservationTime += fRunList[i].tOn;
        addValueToHistogram( hErecTotalTimeDeadTimeCorrected, fRunList[i].tOn * fRunList[i].deadTimeFraction, log10( fRunList[i].energyThreshold ) );
        fTotalObservationTimeDeadTimeCorrected += fRunList[i].tOn * fRunList[i].deadTimeFraction;
        // add current histogram to combined histogram (take energy threshold into account)
        // (counting histograms are not dead time corrected)
        addHistogram( hErecCountsOn, i_hErecCountsOn, log10( fRunList[i].energyThreshold ) );
        addHistogram( hErecCountsOff, i_hErecCountsOff, log10( fRunList[i].energyThreshold ) );
        
        // calculation of mean normalisation factor
        if( fRunList[i].alpha > 0. )
        {
            i_alpha_U += fRunList[i].alpha / ( 1. + fRunList[i].alpha ) * ( fRunList[i].NOn + fRunList[i].NOff );
            i_alpha_L += 1. / ( 1. + fRunList[i].alpha ) * ( fRunList[i].NOn + fRunList[i].NOff );
        }
        
        // calculation of mean effective area (weighted by observation time)
        addValueToHistogram( hEffArea, i_gEff, fRunList[i].tOn, log10( fRunList[i].energyThreshold ) );
        
        // calculate differential flux for this run
        if( fRunList[i].alpha > 0. && fRunList[i].tOn > 0. )
        {
            fDifferentialFlux_perRun.push_back( calculateDifferentialFluxesPerRun( i_hErecCountsOn,
                                                i_hErecCountsOff,
                                                fRunList[i].alpha,
                                                fRunList[i].tOn * fRunList[i].deadTimeFraction,
                                                fRunList[i].tOn,
                                                i_gEff,
                                                log10( fRunList[i].energyThreshold ) ) );
        }
        
        z++;
    } // end of loop over all runs in runlist
    ////////////////////////////////////////////////
    
    // error: histograms for no runs have been found
    if( z == 0 )
    {
         cout << "VEnergySpectrum::combineRuns(): ";
         cout << "ERROR , no histograms for any run found." << endl;
         return false;
    }
    
    if( fDebug == 3 )
    {
        cout << "COUNTING HISTOGRAMS: ";
        cout << i_nonCounts_counter << "\t" << i_noffCounts_counter << endl;
    }
    if( fDebug == 3 )
    {
        cout << "COUNTING HISTOGRAMS (total) after rebinning: ";
        cout << hErecCountsOn->GetEntries() << "\t" << hErecCountsOff->GetEntries() << endl;
    }
    // normalize effective area
    hEffArea->Divide( hErecTotalTime );
    
    // calculate weighted mean normalisation factor
    // (averaged over all runs and all energy bins)
    if( i_alpha_L > 0. )
    {
        fTotalNormalisationFactor = i_alpha_U / i_alpha_L;
    }
    
    //////////////////////////////
    // combine flux vector
    calculateDifferentialFluxes();
    
    cout << "total time [s]: " << fTotalObservationTime;
    cout << " (dead time corrected [s]: " << fTotalObservationTimeDeadTimeCorrected << ", ";
    cout << "total norm factor: " << fTotalNormalisationFactor << ")" << endl;
    
    bCombineRuns = true;
    
    return true;
}

/*

    set significance parameters

    (decide where to plot an upper limit or a differential flux point)

    change only if you know what you do

*/
void VEnergySpectrum::setSignificanceParameters( double iSig, double iMinEvents, double iUL, bool iForFit, int iLiAndMa, int iULAlgo )
{
    fAnalysisSignificance = iSig;
    fAnalysisMinEvents = iMinEvents;
    fAnalysisUpperLimits = iUL;
    fAnalysisLiAndMaEquation = iLiAndMa;
    fAnalysisUpperLimitAlgorithm = iULAlgo;
    fAnalysisForFit = iForFit;
    
    if( fAnalysisForFit && ( fAnalysisSignificance > -5 || fAnalysisMinEvents > -5 ) )
    {
        cout << "ERROR: poor choice of min. significance or min. events for fitting spectrum. Set negative values to fit all points.  " << endl;
    }
}


/*!
    add observation time from current run

    this calculation takes the energy threshold into account (in log10 TeV)
*/
void VEnergySpectrum::addValueToHistogram( TH1* h, double iTObs, double iEThreshold_log10TeV )
{
    if( fDebug == 1 )
    {
        cout << "VEnergySpectrum::addValueToHistogram " << h << " " << iEThreshold_log10TeV << " " << iTObs << endl;
    }
    if( !h )
    {
        return;
    }
    
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        if( !fAnalysisHistogramAddingUseLowEdge && h->GetBinLowEdge( i ) >= iEThreshold_log10TeV )
        {
            h->SetBinContent( i, h->GetBinContent( i ) + iTObs );
        }
        if( fAnalysisHistogramAddingUseLowEdge && ( h->GetBinLowEdge( i ) + h->GetBinWidth( i ) ) >= iEThreshold_log10TeV )
        {
            h->SetBinContent( i, h->GetBinContent( i ) + iTObs );
        }
    }
}

/*

   add values from a graph to a histogram (evaluated at each bin)

   IMPORTANT: THIS IS DANGEROUS TO USE, AS VALUES BEYOND THE GRAPH'S RANGE ARE SIMPLY INTERPOLATED

    this calculation takes the energy threshold into account (in log10 TeV)
*/
void VEnergySpectrum::addValueToHistogram( TH1* h, TGraph* g, double iTObs, double iEThreshold_log10TeV )
{
    if( fDebug == 1 )
    {
        cout << "VEnergySpectrum::addValueToHistogram " << h << " " << iEThreshold_log10TeV << " " << endl;
    }
    if( !h )
    {
        return;
    }
    if( !g )
    {
        return;
    }
    
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        double x = VMathsandFunctions::getMeanEnergyInBin( fEnergyInBinDefinition, h->GetXaxis()->GetBinLowEdge( i ),
                   h->GetXaxis()->GetBinUpEdge( i ), fPlottingSpectralWeightForBinCenter );
        if( x < -1.e30 )
        {
            continue;
        }
        
        if( g->Eval( x ) > 0. )
        {
            if( !fAnalysisHistogramAddingUseLowEdge && h->GetBinLowEdge( i ) >= iEThreshold_log10TeV )
            {
                h->SetBinContent( i, h->GetBinContent( i ) + g->Eval( x ) * iTObs );
            }
            if( fAnalysisHistogramAddingUseLowEdge && ( h->GetBinLowEdge( i ) + h->GetBinWidth( i ) ) >= iEThreshold_log10TeV )
            {
                h->SetBinContent( i, h->GetBinContent( i ) + g->Eval( x ) * iTObs );
            }
        }
    }
}


/*

    this function is destructive to h2

    this calculation takes the energy threshold into account (in log10 TeV)
*/
void VEnergySpectrum::addHistogram( TH1* h1, TH1* h2, double iEThreshold_log10TeV )
{
    if( fDebug == 1 )
    {
        cout << "VEnergySpectrum::addHistogram " << h1 << " " << h2 << " " << iEThreshold_log10TeV << " " << endl;
    }
    
    if( !h1 || !h2 )
    {
        return;
    }
    
    if( h1->GetNbinsX() != h2->GetNbinsX() )
    {
        cout << "Error in addHistogram( TH1*h1, TH1* h2, iEThreshold_log10TeV ): histograms have different bin definitions" << endl;
        cout << "\t" << h1->GetNbinsX() << " " << h2->GetNbinsX() << endl;
        return;
    }
    
    for( int i = 0; i <= h2->GetNbinsX(); i++ )
    {
        // (default)
        if( !fAnalysisHistogramAddingUseLowEdge && h2->GetBinLowEdge( i ) < iEThreshold_log10TeV )
        {
            h2->SetBinContent( i, 0. );
            h2->SetBinError( i, 0. );
        }
        else if( fAnalysisHistogramAddingUseLowEdge && ( h2->GetBinLowEdge( i ) + h2->GetBinWidth( i ) ) < iEThreshold_log10TeV )
        {
            h2->SetBinContent( i, 0. );
            h2->SetBinError( i, 0. );
        }
    }
    h1->Add( h2 );
}


/*!
    this function is for a logarithmic energy axis

    iFluxes = true: return dN and not dN/dE
*/
void VEnergySpectrum::rebinEnergySpectrum( TH1D* h, double iER )
{
    if( !h || iER < 0. )
    {
        return;
    }
    
    // histogram name
    string itemp = h->GetName();
    
    // get current binning of energy axis
    double iBW = h->GetXaxis()->GetBinWidth( 1 );
    if( iBW - iER > 1.e-5 && !bUseRebinner )
    {
        cout << "VEnergySpectrum::rebinEnergySpectrum: error: cannot rebin to smaller than existing bins" << endl;
        cout << "current binning: " << iBW;
        cout << ", requested binning: " << iER;
        cout << " (" << TMath::Abs( iBW - iER ) << " )" << endl;
        return;
    }
    int ngroup = int( iER / iBW + 0.01 );
    if( fabs( ( double )ngroup * iBW - iER ) > 1.e-5 && !bUseRebinner )
    {
        cout << "VEnergySpectrum::rebinEnergySpectrum error: rebinning only possible in multiples of " << iBW << endl;
        cout << "\t" << ngroup << "\t" << iER << "\t" << iBW << endl;
        return;
    }
    if( fDebug == 2 )
    {
        cout << "REBIN HISTOGRAMS: " << itemp;
        cout << ": " << iER << endl;
        cout << "\t bin width: " << iBW << ", ngroup " << ngroup << endl;
    }
    if( ngroup == 1 && !bUseRebinner )
    {
        return;
    }
    
    // counting and timing histograms are simply rebinned
    if( itemp.find( "Counts" ) < itemp.size() || itemp.find( "TotalTime" ) < itemp.size() )
    {
        if( bUseRebinner )
        {
            TH1D* htemp = ( TH1D* )setVariableBinning( h );
            *h = *htemp;
        }
        else
        {
            h->Rebin( ngroup );
        }
        return;
    }
    return;
}

/*!

    calculate differential fluxes (points in energy spectrum)

*/
vector< VDifferentialFluxData > VEnergySpectrum::calculateDifferentialFluxesPerRun( TH1D* i_hErecCountsOn, TH1D* i_hErecCountsOff, double iAlpha,
        double TimeDeadTimeCorrected, double Time, TGraphErrors* i_gEffectiveArea,
        double iEnergyThreshold_log10TeV )
{
    // data vector (which is returned in this function)
    vector< VDifferentialFluxData > iDifferentialFluxPerRun;
    
    // do not do this for spectra given with ascii files
    if( bAsciiDataFile )
    {
        return iDifferentialFluxPerRun;
    }
    // check that all histograms are available
    if( !i_hErecCountsOn || !i_hErecCountsOff || !i_gEffectiveArea || isZombie() )
    {
        return iDifferentialFluxPerRun;
    }
    if( fDebug == 1 )
    {
        cout << "VEnergySpectrum::calculateDifferentialFluxesPerRun() " << i_hErecCountsOn->GetNbinsX() << endl;
    }
    
    double x = 0.;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // loop over all bins in energy spectrum energy spectrum
    if( fDebug >= 2 )
    {
        cout << "VEnergySpectrum::calculateDifferentialFluxesPerRun(): starting loop over herec bins: " << i_hErecCountsOn->GetNbinsX() << endl;
    }
    for( int i = 1; i <= i_hErecCountsOn->GetNbinsX(); i++ )
    {
        /////////////////////////////////////////////////
        // remove points outside wanted energy range
        // (lower limit in energy)
        if( fAnalysisMinEnergy > 0. && i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ) <= log10( fAnalysisMinEnergy ) )
        {
            if( fDebug == 3 )
            {
                cout << "\t\t failed low energy cuts: " << i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i );
                cout << "\t" << log10( fAnalysisMinEnergy ) << endl;
            }
            continue;
        }
        // (upper limit in energy)
        if( fAnalysisMaxEnergy > 0. && i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i ) > log10( fAnalysisMaxEnergy ) )
        {
            if( fDebug == 3 )
            {
                cout << "\t\t failed high energy cuts: " << i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i );
                cout << "\t" << log10( fAnalysisMaxEnergy ) << endl;
            }
            continue;
        }
        
        /////////////////////////////////////////////////
        // new differential flux vector vector
        VDifferentialFluxData i_flux;
        
        // energy (mean point of a log10 bin)
        x = i_hErecCountsOn->GetXaxis()->GetBinCenter( i );
        
        i_flux.Energy = TMath::Power( 10., x );
        // lower and upper edge of energy bin
        i_flux.Energy_lowEdge = TMath::Power( 10., i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ) );
        i_flux.Energy_upEdge  = TMath::Power( 10., i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i ) );
        // mean energy (spectral weighted energy)
        double x = VMathsandFunctions::getMeanEnergyInBin( fEnergyInBinDefinition, i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ),
                   i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i ),
                   fPlottingSpectralWeightForBinCenter );
        if( fEnergyInBinDefinition > 90 )
        {
            cout << "VEnergySpectrum::calculateDifferentialFluxesPerRun() invalid fEnergyInBinDefinition: " << fEnergyInBinDefinition << endl;
            cout << "  allowed values: 0, 1,2 " << endl;
            return iDifferentialFluxPerRun;
        }
        if( x < -1.e20 )
        {
            cout << "VEnergySpectrum::calculateDifferentialFluxesPerRun() invalid energy : ";
            cout << i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ) << "\t" << i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i ) << endl;
            cout << "error" << endl;
            return iDifferentialFluxPerRun;
        }
        
        i_flux.EnergyWeightedMean = TMath::Power( 10., x );
        
        // energy interval dE (on linear axis)
        i_flux.dE = TMath::Power( 10., i_hErecCountsOn->GetXaxis()->GetBinUpEdge( i ) )
                    - TMath::Power( 10., i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ) );
                    
        // get on and off numbers for this bin
        i_flux.NOn        = i_hErecCountsOn->GetBinContent( i );
        i_flux.NOff       = i_hErecCountsOff->GetBinContent( i );
        i_flux.NOff_alpha = iAlpha;
        
        // energy threshold given for this run
        if( i_hErecCountsOn->GetXaxis()->GetBinLowEdge( i ) > iEnergyThreshold_log10TeV )
        {
            // observation time
            i_flux.ObsTime = TimeDeadTimeCorrected;
            i_flux.ExposureTime = Time;
            
            // calculate flux normalization and flux for non-empty bins only
            double iNormalization = 1. / i_flux.dE;
            iNormalization /= i_flux.ObsTime;
            if( i_gEffectiveArea )
            {
                double iEff = i_gEffectiveArea->Eval( log10( i_flux.Energy ), 0, "S" );
                if( iEff > 0. )
                {
                    iNormalization /= ( iEff * 1.e4 );
                    i_flux.EffectiveArea = iEff * 1.e4;
                }
                else
                {
                    iNormalization = 0.;
                    i_flux.EffectiveArea = 0.;
                }
            }
            else
            {
                continue;
            }
            calculateDifferentialFluxes_ErrorsAndSignificances( i_flux, false, iNormalization );
            
        }
        i_flux.fillEvent( fRunList_MJD_min, fRunList_MJD_max );
        // add this bin to data vector
        iDifferentialFluxPerRun.push_back( i_flux );
        
        // debug output
        if( fDebug >= 2 )
        {
            cout << setprecision( 8 ) << "E " << i_flux.Energy_lowEdge << " - " << i_flux.Energy_upEdge;
            cout << ", NOn " << i_flux.NOn;
            cout << ", NOff " << i_flux.NOff;
            cout << ", NDif " << i_flux.NOn - i_flux.NOff* iAlpha;
            cout << ", Sign " << i_flux.Significance;
            cout << ", Norm " << fTotalNormalisationFactor;
            cout << ", TOn " << i_flux.ObsTime;
            cout << ", Flux " << scientific << i_flux.DifferentialFlux;
            cout << ", (F2 " << ( i_flux.NOn - i_flux.NOff * iAlpha ) / i_flux.dE / i_flux.ObsTime
                 / hEffArea->GetBinContent( hEffArea->FindBin( log10( i_flux.Energy ) ) ) * 1.e-4 << ")";
            cout << fixed << endl;
            cout << endl;
        }
    }  // end of loop over all energy bins
    ///////////////////////////////////////////
    
    return iDifferentialFluxPerRun;
}


/*

      calculate average energy spectrum over all runs

*/
void VEnergySpectrum::calculateDifferentialFluxes()
{
    fDifferentialFlux.clear();
    
    if( fDifferentialFlux_perRun.size() == 0 )
    {
        return;
    }
    
    // define and reset flux vector
    fDifferentialFlux = fDifferentialFlux_perRun[0];
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        fDifferentialFlux[i].resetCountsAndFluxes();
    }
    
    // check first that all run-wise flux vectors have the same size
    // (this should always be the case)
    for( unsigned int i = 1; i < fDifferentialFlux_perRun.size(); i++ )
    {
        if( fDifferentialFlux_perRun[i].size() != fDifferentialFlux_perRun[0].size() )
        {
            cout << "======================================================================================================================" << endl;
            cout << "VEnergySpectrum::calculateDifferentialFluxes error, run-wise flux vectors differ in size (element " << i << ")" << endl;
            cout << "\t expect " << fDifferentialFlux_perRun[0].size() << ", found " << fDifferentialFlux_perRun[i].size() << endl;
            cout << "======================================================================================================================" << endl;
            return;
        }
    }
    
    //////////////////////////////////////////////////////////////////
    // loop over all energy bins (all runs with same energy binnning)
    for( unsigned int j = 0; j < fDifferentialFlux_perRun[0].size(); j++ )
    {
        double i_alpha_U = 0.;
        double i_alpha_L = 0.;
        double i_Sum_weight = 0.;
        // loop over all run-wise flux vectors and add them up
        for( unsigned int i = 0; i < fDifferentialFlux_perRun.size(); i++ )
        {
            if( fDifferentialFlux_perRun[i][j].NOff_alpha > -1. )
            {
                fDifferentialFlux[j].ObsTime      += fDifferentialFlux_perRun[i][j].ObsTime;
                fDifferentialFlux[j].ExposureTime += fDifferentialFlux_perRun[i][j].ExposureTime;
                fDifferentialFlux[j].NOn          += fDifferentialFlux_perRun[i][j].NOn;
                fDifferentialFlux[j].NOff         += fDifferentialFlux_perRun[i][j].NOff;
                // modified Li & Ma
                i_alpha_U += fDifferentialFlux_perRun[i][j].NOff_alpha / ( 1. + fDifferentialFlux_perRun[i][j].NOff_alpha )
                             * ( fDifferentialFlux_perRun[i][j].NOn + fDifferentialFlux_perRun[i][j].NOff );
                i_alpha_L += 1. / ( 1 + fDifferentialFlux_perRun[i][j].NOff_alpha )
                             * ( fDifferentialFlux_perRun[i][j].NOn + fDifferentialFlux_perRun[i][j].NOff );
                             
                // calculate differential fluxes
                fDifferentialFlux[j].DifferentialFlux += fDifferentialFlux_perRun[i][j].DifferentialFlux
                        * fDifferentialFlux_perRun[i][j].dE
                        * fDifferentialFlux_perRun[i][j].DifferentialFluxWeighting;
                i_Sum_weight += fDifferentialFlux_perRun[i][j].DifferentialFluxWeighting;
            }
        } // (end) loop over all runs
        
        // alpha calculation
        if( i_alpha_L > 0. )
        {
            fDifferentialFlux[j].NOff_alpha = i_alpha_U / i_alpha_L;
        }
        
        // average differential flux in this energy bin
        if( fDifferentialFlux[j].dE > 0. && i_Sum_weight > 0. )
        {
            fDifferentialFlux[j].DifferentialFlux /= fDifferentialFlux[j].dE;
            fDifferentialFlux[j].DifferentialFlux /= i_Sum_weight;
        }
        
        calculateDifferentialFluxes_ErrorsAndSignificances( fDifferentialFlux[j], true, 1. );
        
        fDifferentialFlux[j].fillEvent( fRunList_MJD_min, fRunList_MJD_max );
    } // loop over all energy bins
    
    ///////////////////////////////////////////////////////////
    // remove all empty bins in the differential flux vector
    std::vector<VDifferentialFluxData>::iterator i = fDifferentialFlux.begin();
    while( i != fDifferentialFlux.end() )
    {
        if( i->NOn == 0 && i->NOff == 0 )
        {
            i = fDifferentialFlux.erase( i );
        }
        else
        {
            ++i;
        }
    }
}

/*

     calculate significances, errors, differential fluxes, upper limits for a given flux vector

     iCombinedAnalysis = true: combined analysis of several runs
     iNormalization    = normalization calculated from observing time and effective areas

*/
void VEnergySpectrum::calculateDifferentialFluxes_ErrorsAndSignificances( VDifferentialFluxData& iFlux,
        bool iCombinedAnalysis,
        double iNormalization )
{
    /////////////////////
    // calculate significance
    // (using Li&Ma, note: event numbers might be too low here for Li & Ma)
    iFlux.Significance = VStatistics::calcSignificance( iFlux.NOn,
                         iFlux.NOff,
                         iFlux.NOff_alpha,
                         fAnalysisLiAndMaEquation );
                         
    /////////////////////
    // Poisson error
    if( iFlux.NOn > 0. )
    {
        iFlux.NOn_error  = sqrt( iFlux.NOn );
    }
    if( iFlux.NOff > 0. )
    {
        iFlux.NOff_error  = sqrt( iFlux.NOff );
    }
    
    ////////////////////////////////////////////////////////
    // calculate differential flux using Rolke et al method
    ////////////////////////////////////////////////////////
    TRolke i_Rolke;
    i_Rolke.SetCLSigmas( fAnalysisSignificance );
    i_Rolke.SetPoissonBkgKnownEff( ( int )iFlux.NOn, ( int )iFlux.NOff, 1. / iFlux.NOff_alpha, 1. );
    
    // calculate flux and Poissonian flux error from event numbers
    double i_ndiff = iFlux.NOn - iFlux.NOff * iFlux.NOff_alpha;
    if( fErrorCalculationMethod == "Rolke" && i_ndiff > 0. )
    {
        // this is an approximation, but very similar result to above value
        i_Rolke.SetCLSigmas( 1.e-4 );
        i_ndiff = 0.5 * ( i_Rolke.GetLowerLimit() + i_Rolke.GetUpperLimit() );
        
        // scale flux to new Rolke mean
        double i_N = iFlux.DifferentialFlux / ( iFlux.NOn - iFlux.NOff * iFlux.NOff_alpha );
        iFlux.DifferentialFlux = i_ndiff * i_N;
        if( iFlux.NOn - iFlux.NOff * iFlux.NOff_alpha != 0. )
        {
            iFlux.DifferentialFlux *= i_ndiff / ( iFlux.NOn - iFlux.NOff * iFlux.NOff_alpha );
        }
    }
    
    
    // normalize by time/area/etc (only for run-wise analysis)
    if( !iCombinedAnalysis )
    {
        iFlux.DifferentialFlux = i_ndiff * iNormalization;
    }
    
    ///////////////////////////////////////////////////////////////
    // calculate flux errors (check first that this is a significant bin)
    ///////////////////////////////////////////////////////////////
    i_Rolke.SetCLSigmas( fAnalysisSignificance );
    if( !iCombinedAnalysis   // never calculate upper limit for individual runs
            || ( ( ( fErrorCalculationMethod == "Rolke" && i_Rolke.GetLowerLimit() > 0 )
                   || ( fErrorCalculationMethod == "Poisson" && iFlux.Significance > fAnalysisSignificance ) )
                 && iFlux.NOn > fAnalysisMinEvents ) )
    {
        // calculate asymmetric flux errors using TRolke
        i_Rolke.SetCLSigmas( 1. );
        if( i_ndiff != 0. )
        {
            iFlux.DifferentialFluxError_low = TMath::Abs( ( i_ndiff - i_Rolke.GetLowerLimit() ) * iFlux.DifferentialFlux / i_ndiff );
            iFlux.DifferentialFluxError_up  = TMath::Abs( ( i_Rolke.GetUpperLimit() - i_ndiff ) * iFlux.DifferentialFlux  / i_ndiff );
            // recalculate poissonian flux error
            iFlux.DifferentialFluxError =  sqrt( iFlux.NOn + iFlux.NOff_alpha * iFlux.NOff_alpha * iFlux.NOff )
                                           * iFlux.DifferentialFlux / i_ndiff;
        }
        // calculate weights for combination of run-by-run fluxes
        iFlux.DifferentialFluxWeighting = iFlux.EffectiveArea * iFlux.ExposureTime;
    }
    ///////////////////////////////////////
    // calculate upper flux limit
    ///////////////////////////////////////
    else
    {
        // note that the upper limit calculation method should match the error calculation method
        double i_UL = VStatistics::calcUpperLimit( iFlux.NOn, iFlux.NOff, iFlux.NOff_alpha, fAnalysisUpperLimits, fAnalysisUpperLimitAlgorithm );
        // scale  upper flux to right value
        if( i_ndiff != 0. )
        {
            iFlux.DifferentialFlux = iFlux.DifferentialFlux  * i_UL / i_ndiff;
        }
        else
        {
            iFlux.DifferentialFlux = -1.e99;
        }
        // flux error is negativ for upper flux value
        iFlux.DifferentialFluxError = -1.;
    }
    
    // debug output
    if( fDebug >= 2 )
    {
        cout << setprecision( 8 ) << "E " << iFlux.Energy_lowEdge << " - " << iFlux.Energy_upEdge;
        cout << ", NOn " << iFlux.NOn;
        cout << ", NOff " << iFlux.NOff;
        cout << ", NDif " << iFlux.NOn - iFlux.NOff* iFlux.NOff_alpha;
        cout << ", Sign " << iFlux.Significance;
        cout << ", Norm " << fTotalNormalisationFactor;
        cout << ", Flux " << scientific << iFlux.DifferentialFlux;
        cout << fixed << endl;
        cout << endl;
    }  // end of loop over all energy bins
    ///////////////////////////////////////////
    
    if( fDebug >= 2 )
    {
        cout << "total events " << fixed << setprecision( 8 ) << iFlux.NOn << " " << iFlux.NOff;
        cout << " " << iFlux.NOff_alpha << " " << iFlux.NOn - iFlux.NOff_alpha* iFlux.NOff << endl;
    }
}



/*
   print differential fluxes

   bSED = true: print vF_v fluxes to be used in VSpectralEnergyDistribution
*/
void VEnergySpectrum::printDifferentialFluxes( bool bSED )
{
    cout << "    \n#----------------------------------------------------------------------------------------------------------------------------#" << endl;
    cout << "#         <E>     E_min    E_max     dE        dN/dE & Err or upper limit      NOn    Err   NOff  Err   alpha    signif.    T "  << endl;
    cout << "#        [TeV]    [TeV]    [TeV]    [TeV]          [1/cm^2/s/TeV]                                                [sigma]   [s]" << endl;
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        cout << setw( 3 ) <<  i << "   ";
        fDifferentialFlux[i].printClean( bSED, fErrorCalculationMethod != "Poisson" );
    }
    cout << "(Error calculation method " << fErrorCalculationMethod << ")" << endl;
    cout << "#------------------------------------------------------------------------------------------------------------------------------#" << endl;
}

/*
   print differential fluxes (for each run)

   bSED = true: print vF_v fluxes to be used in VSpectralEnergyDistribution
*/
void VEnergySpectrum::printDifferentialFluxes_perRun( unsigned int iEnergyBin, bool bSED )
{
    cout << "\n#----------------------------------------------------------------------------------------------------------------------------#" << endl;
    cout << "#      <E>     E_min    E_max     dE        dN/dE & Err or upper limit      NOn  Err   NOff  Err   alpha    signif.    T "  << endl;
    cout << "#     [TeV]    [TeV]    [TeV]    [TeV]          [1/cm^2/s/TeV]                                              [sigma]   [s]" << endl;
    for( unsigned int i = 0; i < fDifferentialFlux_perRun.size(); i++ )
    {
        if( iEnergyBin < fDifferentialFlux_perRun[i].size() )
        {
            fDifferentialFlux_perRun[i][iEnergyBin].printClean( bSED, fErrorCalculationMethod != "Poisson" );
        }
    }
}

/*

    plot the differential energy spectrum taking significances, upper limits, etc into account

*/
TCanvas* VEnergySpectrum::plot( TCanvas* c )
{
    if( isZombie() )
    {
        return 0;
    }
    
    // combine histograms of energy spectra from all the runs in the current run lists (herec)
    if( !bCombineRuns )
    {
        combineRuns();
    }
    
    TH1D* hNull = 0;
    // setup canvas for plotting
    if( c == 0 )
    {
        char hname[600];
        char htitle[600];
        if( fAnalysisForFit )
        {
            sprintf( hname, "c_%s", fDataSetName.c_str() );
        }
        else
        {
            sprintf( hname, "c_%s_UL", fDataSetName.c_str() );
        }
        sprintf( htitle, "energy spectrum (%s)", fDataSetName.c_str() );
        c = new TCanvas( hname, htitle, 10, 10, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hnull_%s", fDataSetName.c_str() );
        hNull = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNull->SetMinimum( fPlottingYaxisMin );
        hNull->SetMaximum( fPlottingYaxisMax );
        hNull->SetStats( 0 );
        hNull->SetXTitle( "log_{10} energy [TeV]" );
        // y-axis: taking multiplication by E^fPlottingMultiplierIndex into account
        hNull->SetYTitle( "dN/dE [cm^{-2}s^{-1}TeV^{-1}]" );
        if( fPlottingMultiplierIndex > 1. )
        {
            sprintf( hname, "E^{%.2f} dN/dE [cm^{-2}s^{-1}TeV^{%.2f}]", fPlottingMultiplierIndex, fPlottingMultiplierIndex - 1. );
            hNull->SetYTitle( hname );
        }
        hNull->GetYaxis()->SetTitleOffset( 1.6 );
        
        // plot an empty histogram with the right axes
        plot_nullHistogram( c, hNull, fPlottingLogEnergyAxis, true, hNull->GetYaxis()->GetTitleOffset(), fPlottingMinEnergy, fPlottingMaxEnergy );
        c->SetLogy( 1 );
    }
    c->cd();
    
    // plot the spectral energy points
    plot_energySpectrum();
    
    fPlottingCanvas = c;
    
    return c;
}

/*

   fill and plot energy spectrum graph

   plot upper limits

   return pointer to graph with differential flux points

*/
TGraphAsymmErrors* VEnergySpectrum::plot_energySpectrum()
{
    // expect that vector with differential fluxes is filled
    if( fDifferentialFlux.size() == 0 )
    {
        return 0;
    }
    
    // loop over vector with differential fluxes and plot fluxes or upper flux limits
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        // upper flux limits
        // (plot as arrows)
        if( fDifferentialFlux[i].DifferentialFluxError < 0. && fPlottingUpperLimits )
        {
            if( fDifferentialFlux[i].DifferentialFlux > fPlottingYaxisMin && fDifferentialFlux[i].DifferentialFlux < fPlottingYaxisMax )
            {
                double i_energy = log10( fDifferentialFlux[i].Energy );
                double i_flux  = fDifferentialFlux[i].DifferentialFlux * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex );
                TArrow* fUL = new TArrow( i_energy, i_flux, i_energy, i_flux * 0.25, 0.03, "|-|>" );
                fUL->SetLineColor( fPlottingColor );
                fUL->SetFillColor( fPlottingColor );
                fUL->SetLineWidth( 1 );
                fUL->Draw();
            }
        }
    }
    // graph with differential flux points
    gEnergySpectrum = getEnergySpectrumGraph();
    if( gEnergySpectrum )
    {
        gEnergySpectrum->Draw( "p" );
        return gEnergySpectrum;
    }
    
    // not successful (!)
    return 0;
}

/*

   fill graph with differential flux points (ignoring points which will be upper limits)

*/
TGraphAsymmErrors* VEnergySpectrum::getEnergySpectrumGraph()
{
    unsigned int z = 0;        // counter
    
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        // differential flux error <0: this is an upper flux limit; ignore
        if( fDifferentialFlux[i].DifferentialFluxError < 0. )
        {
            continue;
        }
        
        // first point: create graph with spectral energies
        // (isn't there a delete missing?)
        if( z == 0 )
        {
            gEnergySpectrum = new TGraphAsymmErrors( 1 );
            gEnergySpectrum->SetMarkerColor( fPlottingColor );
            gEnergySpectrum->SetLineColor( fPlottingColor );
            gEnergySpectrum->SetMarkerSize( fPlottingMarkerSize );
            gEnergySpectrum->SetMarkerStyle( fPlottingMarkerStyle );
        }
        // spectral flux
        gEnergySpectrum->SetPoint( z, log10( fDifferentialFlux[i].EnergyWeightedMean ),
                                   fDifferentialFlux[i].DifferentialFlux * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex ) );
        // error on differential flux
        if( fErrorCalculationMethod == "Poisson" )
        {
            gEnergySpectrum->SetPointEYhigh( z, fDifferentialFlux[i].DifferentialFluxError
                                             * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex ) );
            gEnergySpectrum->SetPointEYlow( z, fDifferentialFlux[i].DifferentialFluxError
                                            * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex ) );
        }
        else if( fErrorCalculationMethod == "Rolke" || fErrorCalculationMethod == "UPDOWN" )
        {
            gEnergySpectrum->SetPointEYhigh( z, fDifferentialFlux[i].DifferentialFluxError_up
                                             * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex ) );
            gEnergySpectrum->SetPointEYlow( z, fDifferentialFlux[i].DifferentialFluxError_low
                                            * TMath::Power( fDifferentialFlux[i].Energy, fPlottingMultiplierIndex ) );
        }
        z++;
    }
    // return graph
    if( z > 0 )
    {
        return gEnergySpectrum;
    }
    
    // not successful (!)
    return 0;
}

/*

   fit to the energy spectrum

   fit function is set with setSpectralFitFunction()

   fSpectralFitFunction == 0 :  power law
   fSpectralFitFunction == 1 :  power law with exponential cutoff
   fSpectralFitFunction == 2 :  broken power law
   fSpectralFitFunction == 3 :  curved power law

*/
TF1* VEnergySpectrum::fitEnergySpectrum( string iname, bool bDraw )
{
    // new fitter
    // (delete missing?)
    fSpectralFitter = new VSpectralFitter( iname );
    // set fit function (see IDs in function comments)
    fSpectralFitter->setSpectralFitFunction( fSpectralFitFunction );
    // all energies in TeV
    fSpectralFitter->setSpectralFitFluxNormalisationEnergy( fSpectralFitFluxNormalisationEnergy );
    fSpectralFitter->setSpectralFitRangeLin( fSpectralFitEnergy_min, fSpectralFitEnergy_max );
    fSpectralFitter->setPlottingStyle( fPlottingEnergySpectrumFitLineColor, fPlottingEnergySpectrumFitLineStyle, fPlottingEnergySpectrumFitLineWidth );
    
    // get energy spectrum
    gEnergySpectrum = getEnergySpectrumGraph();
    if( gEnergySpectrum )
    {
        // perform the fit
        TF1* f = fSpectralFitter->fit( gEnergySpectrum );
        if( f )
        {
            // draw everything
            if( bDraw && fAnalysisForFit )
            {
                f->Draw( "same" );
            }
            // print results
            fSpectralFitter->print();
            // return pointer to fit function
            return f;
        }
    }
    
    // not successful (!)
    return 0;
}

double VEnergySpectrum::calculateIntegralFluxFromFitFunction( double iMinEnergy_TeV, double iMaxEnergy_TeV )
{
    if( fSpectralFitter )
    {
        printf( "Flux F(E > %.2f TeV) [cm^-2 s^-1]: ", iMinEnergy_TeV );
        printf( "%.2e +- ", fSpectralFitter->getIntegralFlux( iMinEnergy_TeV, iMaxEnergy_TeV ) );
        printf( "%.2e\n",   fSpectralFitter->getIntegralFluxError( iMinEnergy_TeV, iMaxEnergy_TeV ) );
    }
    
    return 0.;
}

double VEnergySpectrum::getIntegralFluxFromFitFunction( double iMinEnergy_TeV, double iMaxEnergy_TeV )
{
    if( fSpectralFitter )
    {
        return fSpectralFitter->getIntegralFlux( iMinEnergy_TeV, iMaxEnergy_TeV );
    }
    
    return 0.;
}

double VEnergySpectrum::getIntegralFluxErrorFromFitFunction( double iMinEnergy_TeV, double iMaxEnergy_TeV )
{
    if( fSpectralFitter )
    {
        return fSpectralFitter->getIntegralFluxError( iMinEnergy_TeV, iMaxEnergy_TeV );
    }
    
    return 0.;
}

double VEnergySpectrum::getIntegralFluxFromHistogram( )
{
    if( gEnergySpectrum )
    {
        double flux = 0;
        for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
        {
            flux += fDifferentialFlux[i].DifferentialFlux * fDifferentialFlux[i].dE ;
        }
        return flux;
    }
    
    return 0.;
}

double VEnergySpectrum::getIntegralFluxErrorFromHistogram( )
{
    if( gEnergySpectrum )
    {
        double fluxE = 0;
        for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
        {
            fluxE += TMath::Power( fDifferentialFlux[i].DifferentialFluxError * fDifferentialFlux[i].dE , 2 ) ;
        }
        return sqrt( fluxE );
    }
    
    return 0.;
}




/*

    set plotting style for fit function

*/
void VEnergySpectrum::setSpectralFitPlottingStyle( int iColor, int iStyle, float iWidth )
{
    fPlottingEnergySpectrumFitLineColor = iColor;
    fPlottingEnergySpectrumFitLineStyle = iStyle;
    fPlottingEnergySpectrumFitLineWidth = iWidth;
}


/*

    write fit results (e.g. index and flux normalization) into the canvas

*/
void VEnergySpectrum::plotFitValues( bool iForce )
{
    TCanvas* c = 0;
    if( iForce )
    {
        c = ( TCanvas* )gROOT->FindObject( "c_E_UL" );
    }
    else
    {
        c = fPlottingCanvas;
    }
    
    if( !c )
    {
        cout << "VEnergySpectrum::plotFitValues() error: no canvas to plot things" << endl;
        return;
    }
    
    c->cd();
    
    if( !fSpectralFitter || !fSpectralFitter->getSpectralFitFunction() )
    {
        cout << "VEnergySpectrum::plotFitValues() error: no fit function" << endl;
        return;
    }
    TF1* fEnergy = fSpectralFitter->getSpectralFitFunction();
    
    char hname[500];
    float iTextSize = 0.030;
    if( fSpectralFitFunction > 1 )
    {
        iTextSize = 0.025;
    }
    TLatex* tL1 = new TLatex();
    sprintf( hname, "dN/dE =" );
    tL1->SetNDC( 1 );
    tL1->SetTextSize( iTextSize );
    tL1->DrawLatex( 0.18, 0.23, hname );
    TLatex* tL2 = new TLatex();
    // get exponent (assume negativ exponent)
    int i_expV = ( int )( TMath::Log10( fEnergy->GetParameter( 0 ) ) - 0.5 );
    // get mantissa
    double i_manV = fEnergy->GetParameter( 0 ) * TMath::Power( 10, -1.*( double )i_expV );
    double i_manE = fEnergy->GetParError( 0 ) * TMath::Power( 10, -1.*( double )i_expV );
    // get spectral index
    double i_indexV = fEnergy->GetParameter( 1 ) - fPlottingMultiplierIndex;
    double i_indexE = fEnergy->GetParError( 1 );
    // cutoff energy
    double i_ecutoffV = 0.;
    double i_ecutoffE = 0.;
    // curvature index
    double i_curvatureV = 0.;
    double i_curvatureE = 0.;
    // break energy and second slope for a BPL fit
    double i_breakenergyV = 0.;
    double i_breakenergyE = 0.;
    double i_index2V = 0.;
    double i_index2E = 0.;
    
    if( fSpectralFitFunction == 0 )
    {
        // 1) simple power law
        sprintf( hname, "(%.2f#pm%.2f)#times 10^{%d} (E/%.2f TeV)^{%.2f#pm%.2f} [cm^{-2}s^{-1}TeV^{-1}]",
                 i_manV, i_manE, i_expV, fSpectralFitter->getSpectralFitNormalisationEnergy(),
                 i_indexV, i_indexE );
    }
    else if( fSpectralFitFunction == 1 )
    {
        // get energy cutoff
        i_ecutoffV = fEnergy->GetParameter( 2 );
        i_ecutoffE = fEnergy->GetParError( 2 );
        
        // 2) power law with exponential cutoff
        sprintf( hname, "(%.2f#pm%.2f)#times 10^{%d} (-E/%.2f TeV)^{%.2f#pm%.2f}e^{E/(%.2f#pm%.2f)} [cm^{-2}s^{-1}TeV^{-1}]",
                 i_manV, i_manE, i_expV, fSpectralFitter->getSpectralFitNormalisationEnergy(),
                 i_indexV, i_indexE, i_ecutoffV, i_ecutoffE );
    }
    else if( fSpectralFitFunction == 2 )
    {
        // get break energy and second index
        i_breakenergyV = fEnergy->GetParameter( 3 );
        i_breakenergyE = fEnergy->GetParError( 3 );
        i_index2V = fEnergy->GetParameter( 2 );
        i_index2E = fEnergy->GetParError( 2 );
        
        // 3) broken power law
        sprintf( hname, "#splitline{(%.2f#pm%.2f)#times 10^{%d} (E/(%.2f#pm%0.02f) TeV)^{#Gamma}[cm^{-2}s^{-1}TeV^{-1}];}{   E<%0.2f, #Gamma = %0.2f#pm%0.2f else #Gamma = %0.2f#pm%0.2f}", i_manV, i_manE, i_expV, i_breakenergyV, i_breakenergyE, i_breakenergyV, i_indexV, i_indexE, i_index2V, i_index2E );
    }
    else if( fSpectralFitFunction == 3 )
    {
        // get curvature index
        i_curvatureV = fEnergy->GetParameter( 2 );
        i_curvatureE = fEnergy->GetParError( 2 );
        
        // 4) curved power law
        sprintf( hname, "(%.2f#pm%.2f)#times 10^{%d} (E/%.2f TeV)^{%.2f#pm%.2f + (%.2f#pm%.2f) log_{10}(E/%.2f TeV)} [cm^{-2}s^{-1}TeV^{-1}]",
                 i_manV, i_manE, i_expV, fSpectralFitter->getSpectralFitNormalisationEnergy(),
                 i_indexV, i_indexE, i_curvatureV, i_curvatureE, fSpectralFitter->getSpectralFitNormalisationEnergy() );
    }
    
    
    tL2->SetNDC( 1 );
    tL2->SetTextSize( iTextSize );
    tL2->DrawLatex( 0.18, 0.19, hname );
    TLatex* tL3 = new TLatex();
    double irc2 = 0.;
    if( fEnergy->GetNDF() > 0. )
    {
        irc2 = fEnergy->GetChisquare() / fEnergy->GetNDF();
    }
    sprintf( hname, "#chi^{2}/dof=%.2f/%d (%.1f)", fEnergy->GetChisquare(), ( int )fEnergy->GetNDF(), irc2 );
    tL3->SetNDC( 1 );
    tL3->SetTextSize( iTextSize );
    tL3->DrawLatex( 0.18, 0.15, hname );
}

/*

    plot residuals between fit function and differential energy spectrum

*/
TCanvas* VEnergySpectrum::plotResiduals( TCanvas* c, TF1* f )
{
    if( !fSpectralFitter )
    {
        return 0;
    }
    if( !gEnergySpectrum )
    {
        return 0;
    }
    
    gEnergySpectrumFitResiduals = new TGraphAsymmErrors( gEnergySpectrum->GetN() );
    gEnergySpectrumFitResiduals->SetMarkerColor( fPlottingColor );
    gEnergySpectrumFitResiduals->SetLineColor( fPlottingColor );
    gEnergySpectrumFitResiduals->SetMarkerSize( fPlottingMarkerSize );
    gEnergySpectrumFitResiduals->SetMarkerStyle( fPlottingMarkerStyle );
    
    if( fSpectralFitter && !f )
    {
        f = fSpectralFitter->getSpectralFitFunction();
        if( !f )
        {
            return 0;
        }
    }
    
    double x, y;
    double ye_l = 0.;
    double ye_u = 0.;
    for( int i = 0; i < gEnergySpectrumFitResiduals->GetN(); i++ )
    {
        gEnergySpectrum->GetPoint( i, x, y );
        if( fErrorCalculationMethod == "Poisson" )
        {
            ye_l = gEnergySpectrum->GetErrorY( i );
            ye_u = gEnergySpectrum->GetErrorY( i );
        }
        else
        {
            ye_l = gEnergySpectrum->GetErrorYlow( i );
            ye_u = gEnergySpectrum->GetErrorYhigh( i );
        }
        if( f->Eval( x ) > 0. )
        {
            gEnergySpectrumFitResiduals->SetPoint( i, x, ( y - f->Eval( x ) ) / f->Eval( x ) );
            gEnergySpectrumFitResiduals->SetPointEYhigh( i, ye_u / f->Eval( x ) );
            gEnergySpectrumFitResiduals->SetPointEYlow( i, ye_l / f->Eval( x ) );
        }
    }
    if( !c )
    {
        char hname[600];
        char htitle[600];
        sprintf( hname, "cRes_%s", fDataSetName.c_str() );
        sprintf( htitle, "energy spectrum, residuals (%s)", fDataSetName.c_str() );
        c = new TCanvas( hname, htitle, 700, 10, 500, 500 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hNullGF_%s", fDataSetName.c_str() );
        TH1D* hNullGF = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNullGF->SetStats( 0 );
        hNullGF->SetXTitle( "log_{10} energy [TeV]" );
        hNullGF->SetYTitle( "(f_{rec} - f_{fit}) / f_{fit}" );
        hNullGF->GetYaxis()->SetTitleOffset( 1.3 );
        hNullGF->SetMaximum( 1. );
        hNullGF->SetMinimum( -1. );
        
        plot_nullHistogram( c, hNullGF, fPlottingLogEnergyAxis, false, hNullGF->GetYaxis()->GetTitleOffset(), fPlottingMinEnergy, fPlottingMaxEnergy );
    }
    c->cd();
    
    gEnergySpectrumFitResiduals->Draw( "p" );
    
    return c;
}

TCanvas* VEnergySpectrum::plotMeanEffectiveArea( TCanvas* c, double i_effMin, double i_effMax )
{
    if( !hEffArea )
    {
        return 0;
    }
    if( !c )
    {
        char hname[600];
        char htitle[600];
        sprintf( hname, "ceff_%s", fDataSetName.c_str() );
        sprintf( htitle, "effective area vs energy (%s)", fDataSetName.c_str() );
        c = new TCanvas( hname, htitle, 1200, 10, 500, 500 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLogy( 0 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hNullEF_%s", fDataSetName.c_str() );
        TH1D* hNullGT = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNullGT->SetStats( 0 );
        hNullGT->SetXTitle( "log_{10} energy_{rec} [TeV]" );
        hNullGT->SetYTitle( "effective area [m^{2}" );
        hNullGT->GetYaxis()->SetTitleOffset( 1.3 );
        hNullGT->SetMinimum( i_effMin );
        if( i_effMax < 0. )
        {
            hNullGT->SetMaximum( hEffArea->GetMaximum() * 1.2 );
        }
        else
        {
            hNullGT->SetMaximum( i_effMax );
        }
        hNullGT->SetMinimum( 0.5 );
        hNullGT->Draw();
    }
    c->cd();
    
    hEffArea->SetLineColor( fPlottingColor );
    hEffArea->Draw( "same" );
    
    return c;
}

/*

    plot histogram with events per energy bin

*/
TCanvas*  VEnergySpectrum::plotCountingHistograms( TCanvas* c )
{
    if( !hErecCountsOn || !hErecCountsOff )
    {
        return 0;
    }
    
    if( !c )
    {
        char hname[600];
        char htitle[600];
        sprintf( hname, "cGT_%s", fDataSetName.c_str() );
        sprintf( htitle, "counts vs energy (%s)", fDataSetName.c_str() );
        c = new TCanvas( hname, htitle, 800, 10, 500, 500 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLogy( 1 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hNullGT_%s", fDataSetName.c_str() );
        TH1D* hNullGT = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNullGT->SetStats( 0 );
        hNullGT->SetXTitle( "log_{10} energy [TeV]" );
        hNullGT->SetYTitle( "counts per bin" );
        hNullGT->GetYaxis()->SetTitleOffset( 1.3 );
        hNullGT->SetMaximum( hErecCountsOn->GetMaximum() * 1.2 );
        hNullGT->SetMinimum( 0.5 );
        hNullGT->Draw();
    }
    c->cd();
    
    hErecCountsOn->Draw( "same" );
    TH1D* hErecCountsOff_Clone = ( TH1D* )hErecCountsOff->Clone();
    hErecCountsOff_Clone->Scale( fTotalNormalisationFactor );
    hErecCountsOff_Clone->Draw( "same" );
    
    // get counting values used for energy reconstruction
    TGraph* gErecOn = new TGraph( 1 );
    gErecOn->SetMarkerStyle( 24 );
    gErecOn->SetMarkerColor( 1 );
    gErecOn->SetLineColor( 1 );
    TGraph* gErecOff = new TGraph( 1 );
    gErecOff->SetMarkerStyle( 25 );
    gErecOff->SetMarkerColor( 2 );
    gErecOff->SetLineColor( 2 );
    cout << "Counting histogram on:  color 1" << endl;
    cout << "Counting histogram off: color 2" << endl;
    
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        gErecOn->SetPoint( i, log10( fDifferentialFlux[i].EnergyWeightedMean ), fDifferentialFlux[i].NOn );
        gErecOff->SetPoint( i, log10( fDifferentialFlux[i].EnergyWeightedMean ), fDifferentialFlux[i].NOff * fTotalNormalisationFactor );
    }
    gErecOn->Draw( "p" );
    gErecOff->Draw( "p" );
    
    return c;
}

/*

    plot life time (dead time corrected) vs energy

*/
TCanvas* VEnergySpectrum::plotLifeTimevsEnergy( TCanvas* c )
{
    if( !getTotalTimeHistogram() )
    {
        return 0;
    }
    
    // total time (not deadtime corrected)
    TH1D* h = ( TH1D* )getTotalTimeHistogram( false )->Clone();
    h->Scale( 1. / 60. );                         //! s -> min
    h->SetFillColor( 9 );
    h->SetFillStyle( 3002 );
    
    if( !c )
    {
        char hname[600];
        char htitle[600];
        sprintf( hname, "cLT_%s", fDataSetName.c_str() );
        sprintf( htitle, "life time vs energy (%s)", fDataSetName.c_str() );
        c = new TCanvas( hname, htitle, 800, 10, 500, 500 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hNullLT_%s", fDataSetName.c_str() );
        TH1D* hNullLT = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNullLT->SetStats( 0 );
        hNullLT->SetXTitle( "log_{10} energy [TeV]" );
        hNullLT->SetYTitle( "life time [min]" );
        hNullLT->GetYaxis()->SetTitleOffset( 1.7 );
        hNullLT->SetMaximum( h->GetMaximum() * 1.2 );
        hNullLT->SetMinimum( 0. );
        
        plot_nullHistogram( c, hNullLT, fPlottingLogEnergyAxis, false, hNullLT->GetYaxis()->GetTitleOffset(), fPlottingMinEnergy, fPlottingMaxEnergy );
    }
    c->cd();
    
    h->Draw( "same" );
    
    // total time (deadtime corrected)
    if( getTotalTimeHistogram( true ) )
    {
        TH1D* h = ( TH1D* )getTotalTimeHistogram( true )->Clone();
        h->Scale( 1. / 60. );                         //! s -> min
        h->SetLineColor( 2 );
        h->Draw( "same" );
    }
    
    return c;
}

/*

   write number of events and significances at each point of the energy spectrum

*/
void VEnergySpectrum::plotEventNumbers( Double_t ts )
{
    if( !fPlottingCanvas )
    {
        cout << "Error: no canvas to plot things" << endl;
        return;
    }
    fPlottingCanvas->cd();
    
    char hnum[500];
    for( unsigned int i = 0; i < fDifferentialFlux.size(); i++ )
    {
        sprintf( hnum, "%.1f events (%.1f #sigma)", ( fDifferentialFlux[i].NOn - fTotalNormalisationFactor * fDifferentialFlux[i].NOff ),
                 fDifferentialFlux[i].Significance );
                 
        double y  = ( fDifferentialFlux[i].DifferentialFlux + fDifferentialFlux[i].DifferentialFluxError ) * 1.75;
        // upper flux limit
        if( fDifferentialFlux[i].DifferentialFluxError < 0. )
        {
            y = fDifferentialFlux[i].DifferentialFlux * 1.75;
        }
        
        TLatex* iT = new TLatex( log10( fDifferentialFlux[i].Energy ), y, hnum );
        iT->SetTextSize( ts );
        iT->SetTextAngle( 90. );
        iT->Draw();
    }
}

/*

    rebinning: set up template histogram

*/
void VEnergySpectrum::setOriginalBinner( TH1* h )
{
    nOriginalBinner = ( TH1D* )h->Clone( "nOriginalBinner" );
    nOriginalBinner->SetTitle( "nOriginalBinner" );
    
    vector < Double_t > z;
    
    for( unsigned int i = 0; i < newBinningGroupings.size(); i++ )
    {
        z.push_back( nOriginalBinner->GetBinLowEdge( newBinningGroupings[i] ) );
    }
    
    unsigned int Nbins = z.size();
    Double_t binsx[Nbins];
    
    for( unsigned int i = 0; i < Nbins; i++ )
    {
        binsx[i] = z[i];
    }
    
    nRebinner = new TH1D( "nRebinner", "nRebinner", Nbins - 1, binsx );
    for( Int_t i = 1; i <= nRebinner->GetNbinsX(); i++ )
    {
        cout << "nRebinner: new lower bin edges [TeV]: " << TMath::Power( 10, nRebinner->GetBinLowEdge( i ) ) << endl;
    }
    cout << "nRebinner: number of bins: " << Nbins << endl;
}

/*

   rebin histograms (variable binning)

*/
TH1* VEnergySpectrum::setVariableBinning( TH1* a )
{
    Double_t xbins[nRebinner->GetNbinsX() + 1];
    
    for( Int_t i = 1; i <= nRebinner->GetNbinsX(); i++ )
    {
        xbins[i - 1] = nRebinner->GetBinLowEdge( i );
    }
    
    xbins[nRebinner->GetNbinsX()] = nRebinner->GetBinLowEdge( nRebinner->GetNbinsX() ) + nRebinner->GetBinWidth( nRebinner->GetNbinsX() );
    
    
    return ( TH1* )a->Rebin( nRebinner->GetNbinsX(), a->GetName(), xbins );
}

void VEnergySpectrum::setRebinBool( bool i_rebin )
{
    bUseRebinner = i_rebin;
}

bool VEnergySpectrum::setEnergyInBinDefinition( unsigned int iEF )
{
    if( iEF < 3 )
    {
        fEnergyInBinDefinition = iEF;
        return true;
    }
    
    cout << "VEnergySpectrum::setEnergyInBinDefinition error: energy bin definition out of range: " << iEF << endl;
    cout << "possible values: 0, 1, 2" << endl;
    
    fEnergyInBinDefinition = 99;
    
    return false;
}

void VEnergySpectrum::printEnergyBins()
{
    if( hErecCountsOn )
    {
        cout << "Printing energy binning: " << endl;
        for( int i = 1; i <= hErecCountsOn->GetNbinsX(); i++ )
        {
            cout << "Energy bin " << i << ":";
            printf( " bin edges in log10(E/TeV): [%.2f,%.2f],", hErecCountsOn->GetBinLowEdge( i ),
                    hErecCountsOn->GetBinLowEdge( i ) + hErecCountsOn->GetBinWidth( i ) );
            printf( " bin edges in E/TeV: [%.2f,%.2f]", TMath::Power( 10., hErecCountsOn->GetBinLowEdge( i ) ),
                    TMath::Power( 10., hErecCountsOn->GetBinLowEdge( i ) + hErecCountsOn->GetBinWidth( i ) ) );
            cout << endl;
        }
    }
}

/*

    set error and signficance calculation method

    Rolke method is more suitable to bins with low number of events, Li & Ma is good for bins with >10 bins

*/
bool VEnergySpectrum::setErrorCalculationMethod( string iMethod )
{
    if( iMethod == "Rolke" || iMethod == "Poisson" )
    {
        fErrorCalculationMethod = iMethod;
        return true;
    }
    cout << "VEnergySpectrum::setErrorCalculationMethod error: unknown error calculation method (should be Rolke/Poisson)" << endl;
    
    return false;
}

/*
    check rebinning parameters and return numbers if bins to combine
*/
int VEnergySpectrum::getRebinningGrouping( TH1* h, double iNewBinWidth )
{
    // get current binning of energy axis
    // (assume all bins are the same: is this ok?)
    double iBW = h->GetXaxis()->GetBinWidth( 1 );
    if( iBW - iNewBinWidth > 1.e-5 && !bUseRebinner )
    {
        cout << "VEnergySpectrum::checkRebinningParameters: error: cannot rebin to smaller than existing bins" << endl;
        cout << "current binning: " << iBW;
        cout << ", requested binning: " << iNewBinWidth;
        cout << " (" << TMath::Abs( iBW - iNewBinWidth ) << " )" << endl;
        return 0;
    }
    int ngroup = int( iNewBinWidth / iBW + 0.01 );
    if( fabs( ( double )ngroup * iBW - iNewBinWidth ) > 1.e-5 && !bUseRebinner )
    {
        cout << "VEnergySpectrum::checkRebinningParameters error: rebinning only possible in multiples of " << iBW << endl;
        cout << "\t" << ngroup << "\t" << iNewBinWidth << "\t" << iBW << endl;
        return 0;
    }
    
    return ngroup;
}

/*
 * plot energy threshold vs run number and vs zenith angle
 */
void VEnergySpectrum::plotEnergyThresholds()
{
      if( fRunList.size() == 0 )
      {
           return;
      }

      TGraph *hGEnergyTresholdRunNumber = new TGraph( 1 );
      hGEnergyTresholdRunNumber->SetMarkerStyle( 21 );
      hGEnergyTresholdRunNumber->SetTitle( "" );

      TGraph *hGEnergyTresholdZe = new TGraph( 1 );
      hGEnergyTresholdZe->SetMarkerStyle( 21 );
      hGEnergyTresholdZe->SetTitle( "" );


      unsigned int z = 0;
      for( unsigned int i = 0; i < fRunList.size(); i++ )
      {
              if( fRunList[i].elevationOn > 0. )
              {
                  hGEnergyTresholdRunNumber->SetPoint( z, fRunList[i].runnumber, fRunList[i].energyThreshold * 1.e3 );
                  hGEnergyTresholdZe->SetPoint( z, 90.-fRunList[i].elevationOn, fRunList[i].energyThreshold * 1.e3 );
                  z++;
              }
              else if( fRunList[i].elevationOff > 0. )
              {
                  hGEnergyTresholdRunNumber->SetPoint( z, fRunList[i].runnumber, fRunList[i].energyThreshold * 1.e3 );
                  hGEnergyTresholdZe->SetPoint( z, 90.-fRunList[i].elevationOff, fRunList[i].energyThreshold * 1.e3 );
                  z++;
              }
              else
              {
                  cout << "ignoring run with 0 elevation: " << fRunList[i].runnumber << endl;
              }
      }
      char hname[600];
      char htitle[600];
      sprintf( hname, "cET_%s", fDataSetName.c_str() );
      sprintf( htitle, "energy thresholds vs run number (%s)", fDataSetName.c_str() );
      TCanvas *cET = new TCanvas( hname, htitle, 700, 10, 500, 500 );
      cET->SetGridx( 0 );
      cET->SetGridy( 0 );
      gPad->SetLeftMargin( 0.13 );

      hGEnergyTresholdRunNumber->Draw( "ap" );
      hGEnergyTresholdRunNumber->GetHistogram()->GetYaxis()->SetTitleOffset( 1.4 );
      hGEnergyTresholdRunNumber->GetHistogram()->SetXTitle( "run number #" );
      hGEnergyTresholdRunNumber->GetHistogram()->SetYTitle( "energy threshold (TeV)" );
    
      sprintf( hname, "cEZE_%s", fDataSetName.c_str() );
      sprintf( htitle, "energy thresholds vs zenith angle (%s)", fDataSetName.c_str() );
      TCanvas *cEZE = new TCanvas( hname, htitle, 900, 10, 500, 500 );
      cEZE->SetGridx( 0 );
      cEZE->SetGridy( 0 );
      gPad->SetLeftMargin( 0.13 );

      hGEnergyTresholdZe->Draw( "ap" );
      hGEnergyTresholdZe->GetHistogram()->GetYaxis()->SetTitleOffset( 1.4 );
      hGEnergyTresholdZe->GetHistogram()->SetXTitle( "zenith angle (deg)" );
      hGEnergyTresholdZe->GetHistogram()->SetYTitle( "energy threshold (TeV)" );

}

/*
 * print energy thresholds per run
*/ 
void VEnergySpectrum::printEnergyThreshold()
{
    cout << "Energy thresholds per run" << endl;
    cout << "-------------------------" << endl;
    cout << "Method: ";
    if( fAnalysisEnergyThresholdDefinition == 2 )
    {
        cout << "fraction of Eff_max (";
        cout << fAnalysisEnergyThresholdDefinition << "): ";
        cout << fAnalysisMaxEffectiveAreaFraction << endl;
    }
    else if( fAnalysisEnergyThresholdDefinition == 1 )
    {
        cout << "Esys_max (";
        cout << fAnalysisEnergyThresholdDefinition << "): ";
        cout << fAnalysisMaxEnergySystematic << endl;
    }
    else if( fAnalysisEnergyThresholdDefinition == 3 )
    {
        cout << "Fixed value (";
        cout << fAnalysisEnergyThresholdDefinition << "): ";
        cout << fEnergyThresholdFixedValue << " TeV" << endl;
    }
    else
    {
        cout << "no energy threshold method defined)" << endl;
    }
    
    printEnergyThresholds();
}

/*

   short cut to plot the Crab Nebula spectra
   (sets axis correctly, add spectra from other instruments, etc)

*/
TCanvas* VEnergySpectrum::plotCrabNebulaSpectrum( double iPlottingMultiplierIndex, double i_FitStart_TevLin, double i_FitStop_TeVLin,
        double i_EnergyBinningLog10, int i_SpectralFitFunction )
{
    // binning and statistics
    setEnergyBinning( i_EnergyBinningLog10 );
    setEnergyRangeLinear( 0.10, 500. );
    setSignificanceParameters( 2., 3., 0.95, false );
    setSignificanceParameters( -2., -3., 0.95, false );
    setPlottingMultiplierIndex( iPlottingMultiplierIndex );
    
    // plotting
    setPlottingEnergyRangeLinear( 0.08, 50. );
    setPlottingYaxis( 1e-16, 1.e-6 );
    if( iPlottingMultiplierIndex > 2. )
    {
        setPlottingYaxis( 9.e-12, 2.e-10 );
    }
    
    // plot & fit
    TCanvas* c = plot();
    
    setSpectralFitFunction( i_SpectralFitFunction );
    setSpectralFitFluxNormalisationEnergy( 1. );
    setSpectralFitRangeLin( i_FitStart_TevLin, i_FitStop_TeVLin );
    fitEnergySpectrum();
    if( iPlottingMultiplierIndex < 2. )
    {
        plotFitValues();
        plotEventNumbers();
    }
    plotResiduals();
    plotLifeTimevsEnergy();
    plotCountingHistograms();
    plotMeanEffectiveArea();
    
    VEnergySpectrumfromLiterature l( "$VERITAS_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat" );
    l.setPlottingMultiplierIndex( iPlottingMultiplierIndex );
    
    l.setPlottingStyle( 2, 2, 2, 25 );
    l.listValues( 1 );
    l.plot( 1, c );
    l.setPlottingStyle( 3, 2, 2, 25 );
    l.listValues( 2 );
    l.plot( 2, c );
    l.setPlottingStyle( 4, 2, 2, 25 );
    l.listValues( 3 );
    l.plot( 3, c );
    l.setPlottingStyle( 5, 2, 2, 25 );
    l.listValues( 4 );
    l.plot( 4, c );
    l.setPlottingStyle( 6, 2, 2, 25 );
    l.listValues( 5 );
    l.plot( 5, c );
    l.setPlottingStyle( 7, 2, 2, 25 );
    l.listValues( 6 );
    l.plot( 6, c );
    l.setPlottingStyle( 8, 2, 2, 25 );
    l.listValues( 7 );
    l.plot( 7, c );
    if( l.isValidID( 10 ) )
    {
        l.setPlottingStyle( 47, 2, 2, 25 );
        l.listValues( 10 );
        l.plot( 10, c );
    }
    
    plot( c );
    
    return c;
}

/*

Write all basic flux information to a simple text file,
one energy bin per line.  Useful for exporting flux
information to other programs.

*/
int VEnergySpectrum::writeFluxInfoToTextFile( char* output_text_file )
{

    FILE* fp = fopen( output_text_file, "w" ) ;
    for( unsigned int j = 0 ; j < fDifferentialFlux.size() ; j++ )
    {
        fprintf( fp, " %12.6f", fDifferentialFlux[j].Energy_lowEdge ) ;            // TeV
        fprintf( fp, " %12.6f", fDifferentialFlux[j].Energy_upEdge ) ;             // TeV
        fprintf( fp, " %15.6e", fDifferentialFlux[j].DifferentialFlux ) ;          // 1/cm2/s/TeV
        fprintf( fp, " %15.6e", fDifferentialFlux[j].DifferentialFluxError_up ) ;  // 1/cm2/s/TeV
        fprintf( fp, " %15.6e", fDifferentialFlux[j].DifferentialFluxError_low ) ; // 1/cm2/s/TeV
        fprintf( fp, " %12.6f", fDifferentialFlux[j].NOn ) ;                       // # events
        fprintf( fp, " %12.6f", fDifferentialFlux[j].NOn_error ) ;                 // # events
        fprintf( fp, " %12.6f", fDifferentialFlux[j].NOff ) ;                      // # events
        fprintf( fp, " %12.6f", fDifferentialFlux[j].NOff_error ) ;                // # events
        fprintf( fp, " %12.6f", fDifferentialFlux[j].NOff_alpha ) ;                // unitless
        fprintf( fp, " %12.6f", fDifferentialFlux[j].Significance ) ;
        fprintf( fp, "\n" ) ;
    }
    return fclose( fp );
    
}

/*
 *
 *
*/
int VEnergySpectrum::writeFitInfoToTextFile( char* output_text_file )
{

    FILE* fp = fopen( output_text_file, "w" ) ;
    
    // these lines must be space-delimited!
    fprintf( fp, "FitEquation dN/dE=I*(E/1TeV)^{-Gamma}\n" ) ;
    fprintf( fp, "I %12.6e "                    , getSpectralFitFunction()->GetParameter( 0 ) ) ;
    fprintf( fp, "+- %12.6e 1/cm2/s/TeV\n"      , getSpectralFitFunction()->GetParError( 0 ) ) ;
    fprintf( fp, "Gamma %.6f "                  , getSpectralFitFunction()->GetParameter( 1 ) ) ;
    fprintf( fp, "+- %.6f  \n"                  , getSpectralFitFunction()->GetParError( 1 ) ) ;
    fprintf( fp, "ChiSquared %.6f\n"            , getSpectralFitFunction()->GetChisquare() ) ;
    fprintf( fp, "NumberOfDegreesOfFreedom %d\n", getSpectralFitFunction()->GetNDF() ) ;
    
    return fclose( fp ) ;
}
