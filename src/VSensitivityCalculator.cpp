/*! \class VSensitivityCalculator
    \brief calculate sensitivity from data rates, Monte Carlo, or energy spectra

################################################################################################
## plot and list observation time necessary for a given source strength
################################################################################################

- input: gamma-ray, background rates from Crab Nebula data, background normalisation parameter

VSensitivityCalculator a;
// set gamma-ray, background rates from Crab Nebula data, background normalisation parameter
a.addDataSet( 4.7, 0.21, 0.1, "A" );
// set range of source strengths to be covered (in Crab Units)
a.setSourceStrengthRange_CU(1.e-3, 40. );
// set source strengths (in Crab Units) to be plotted with guided lines
vector< double > v;
v.push_back( 1. ); v.push_back( 0.1 ); v.push_back( 0.03 ); v.push_back( 0.01 );
a.setSourceStrengthVector_CU( v );

// plot observation time necessary for a given source strength
a.plotObservationTimevsFlux( 0, 0, 2 );
// list observation time necessary for a given source strength
a.list_sensitivity( 0 );

################################################################################################
## plot integral and differential sensitivity vs energy using a reconstructed Crab spectrum
################################################################################################

- input: anasum results file from Crab analysis

VSensitivityCalculator a;
// read Crab spectrum and plot integral sensitivity in Crab Unit (CU)
// (replace "CU" by "PFLUX" for particle flux [1/m2/s] or energy flux "ENERGY' [erg/cm2/s])
// (for units other than CU, you need to specify an expected Crab spectrum:
a.setEnergySpectrumfromLiterature( "$VERITAS_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat", 1 );
TCanvas *c = a.plotIntegralSensitivityvsEnergyFromCrabSpectrum(0, "myCrabFile.root", 1, "CU" )
// plot a sensitivity curve from literature on top of it
a.plotSensitivityvsEnergyFromTextTFile( c, "SensitivityFromTextFile_CU", 2, 2, 2, "CU" );
// plot differential sensitivity (4 bins per decade)
a.plotDifferentialSensitivityvsEnergyFromCrabSpectrum( 0, "myCrabFile.root", 1, "ENERGY", 0.25 );

################################################################################################
## plot integral and  differential sensitivity vs energy using MC
################################################################################################

-input: MC gamma and proton effective area

VSensitivityCalculator a;
// set MC parameter for gamma-ray simulations
a.setMonteCarloParameters(1, "../data/TeV_data/EnergySpectrum_literatureValues.dat", 1, "/Users/Shared/data/AGIS/Set-Ia/effectiveAreas/CT_SignalNoise_StdCuts_NT2/effectiveArea-CFG1-0100m-ID1.root" );
// set MC parameter for proton simulations
a.setMonteCarloParameters(14, "../data/TeV_data/EnergySpectrum_literatureValues_CR.dat", 0, "/Users/Shared/data/AGIS/Set-Ib/effectiveAreas/effectiveArea-CFG1-0100m-ID1.root", 20., 0, 0.5, 150, 2.5, 2. );
// plot integral sensitivity
a.plotIntegralSensitivityvsEnergyFromCrabSpectrum( 0, "MC" );

################################################################################################
## plot differential sensitivity vs energy using CTA WP Phys files
################################################################################################

VSensitivityCalculator a;
a.setMonteCarloParametersCTA_MC( "data/DESY/DESY.d20140105.Erec1.V2.ID0_180degNIM2.prod2-Aar-lowE-NS.S.2a.180000s.root", 0., "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat", 5 );
a.plotDifferentialSensitivityvsEnergyFromCrabSpectrum( 0, "CTA-MC", 1, "ENERGY" );

*/

#include "VSensitivityCalculator.h"

ClassImp( VSensitivityCalculator )

VSensitivityCalculator::VSensitivityCalculator()
{
    fConstant_Flux_To_Ergs = TMath::Qe() / 1.e-7;
    setDebug( false );
    
    fEnergySpectrumfromLiterature = 0;
    fEnergySpectrumfromLiterature_ID = 0;
    
    fTMVAEvaluatorResults = 0;
    
    hnull = 0;
    
    reset();
}


void VSensitivityCalculator::reset()
{
    fCurrentDataSet = 0;
    
    fPlotDebugName = "";
    
    // Li & Ma formula for significance calculation (use 17, 5 or 9)
    fLiAndMaEqu = 17;
    
    fSetEvents_minCutOnly = false;
    
    // source strength values (in Crab units)
    // units: fraction of source strenght
    fSourceStrength_min = 0.0003;
    fSourceStrength_max = 3.5;
    fSourceStrength_step = 0.005;
    setSourceStrengthVector_CU();
    // values of Crab fluxes to be plotted as lines into the sensitivity vs energy graph (in Crab Units)
    fPlottingCrabFlux_CU.push_back( 1.e2 );
    fPlottingCrabFlux_CU.push_back( 1.e1 );
    fPlottingCrabFlux_CU.push_back( 1.e0 );
    fPlottingCrabFlux_CU.push_back( 1.e-1 );
    fPlottingCrabFlux_CU.push_back( 1.e-2 );
    fPlottingCrabFlux_CU.push_back( 1.e-3 );
    fPlottingCrabFlux_CU.push_back( 1.e-4 );
    
    // sensitivity graph
    gSensitivityvsEnergy = 0;
    gBGRate = 0;
    gProtonRate = 0;
    gElectronRate = 0;
    
    setPlotCanvasSize();
    setPlotCanvasName();
    setEnergyRange_Log();
    setFluxRange_PFLUX();
    setFluxRange_ENERG();
    setFluxRange_CU();
    setObservationTimeRange();
    setSignificanceParameter();
    setBackgroundMissingParticleFraction();
    setRequireCutsToBeOptimized();
    setMaximumEnergyBias();
    setPlotCrabLines();
    
    fGraphObsvsTime.clear();
    fData.clear();
    
    fMCCTA_File = "";
    fMCCTA_cameraoffset_deg = 0.;
    fEnergy_dE_log10 = 0.25;
}

/*
   iD = data set ID
   iObservationTime = total observation time [h]

   energy           = energy on linear scale [TeV] (for debug output only)

   return sensitivity as fraction of data set used
*/
double VSensitivityCalculator::getSensitivity( unsigned int iD, double energy, unsigned int iFillStatistics )
{
    if( !checkDataSet( iD, "getSensitivity" ) )
    {
        return 0.;
    }
    
    return getSensitivity( fData[iD].fSignal, fData[iD].fBackground, fData[iD].fAlpha, energy, iFillStatistics );
}


/*
   iD = data set ID
   iObservationTime = total observation time [h]

   energy           = energy on linear scale [TeV] (for debug output only)

   return sensitivity as fraction of data set used
*/
double VSensitivityCalculator::getSensitivity( double iSignal, double iBackground, double iAlpha, double energy, unsigned int iFillStatistics )
{
    double f = 0.;
    double s = 0.;
    double t = fObservationTime_h * 60.;            // h -> min
    
    // hardwired implementation of dead time
    /*        double iDeadTimeFraction = 1.0;
            iSignal     *= iDeadTimeFraction;
            iBackground *= iDeadTimeFraction; */
    
    // fSignal = gamma-ray + background rates in source region
    // fBackground = background counts in off source regions
    // (observe that this is different from the standard definition used in addDataSet() (which are rates, etc))
    
    // calculate signal (non - alpha * noff) rate
    double n_diff = iSignal - iBackground * iAlpha;
    
    if( fDebug && energy > 0. )
    {
        cout << "getSensitivity energy: " << energy << "\t obstime: " << fObservationTime_h;
        cout << "\t sig > " << fSignificance_min << "\t events >= " << fEvents_min;
        cout << "\t min ratio of signal to background events > " << fMinBackgroundRateRatio_min << endl;
        cout << "\t signal+background rate: " << iSignal;
        cout << "\t background events: " << iBackground;
        cout << "\t background rate: " <<  iBackground* iAlpha << "\t alpha: " << iAlpha << endl;
        cout << "\t signal rate: " << n_diff;
        cout << "\t nsourcestrengths " << fSourceStrength.size() << endl;
    }
    
    bool bSuccess = false;
    // minimum number of background events
    if( iFillStatistics == 4 )
    {
        if( iBackground * iAlpha > 0. )
        {
            return -1.;
        }
        // TMP
        else
        {
            return 0.01;
        }
    }
    
    unsigned i_nMax = fSourceStrength.size() - 1;
    // Loop optimization
    // calculate start for iFillStatistics == 0
    if( iFillStatistics == 0 )
    {
        for( unsigned int n = i_nMax; n > 0; n-- )
        {
            if( fSourceStrength[n] * t * n_diff >= fEvents_min )
            {
                i_nMax = n;
                break;
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////
    // loop over the source strength vector and check different sensitivity criteria
    for( unsigned int n = i_nMax; n > 0; n-- )
    {
        f = fSourceStrength[n];
        
        // default significance calculation
        s = VStatistics::calcSignificance( t * ( f * n_diff + iBackground * iAlpha ), t * iBackground, iAlpha, fLiAndMaEqu );
        // significance calculation for Crab flares (don't use!)
        //        s = VStatistics::calcSignificance( t * ( f * n_diff + iBackground * iAlpha + n_diff),
        //	                                   t * ( iBackground + n_diff / iAlpha ), iAlpha, fLiAndMaEqu );
        
        //////////////////////////////////////////////////////////////////////////
        // check if this set of observations passes the significance criteria
        //////////////////////////////////////////////////////////////////////////
        // require a certain significance
        bool bPassed_MinimumSignificance = ( s >= fSignificance_min );
        // require a minimum number of signal events
        bool bPassed_MinimumSignalEvents = ( t * f * n_diff >= fEvents_min );
        // require background events
        // (removes most sensitivity values at large energies, but otherwise transition zone
        //  between signal and background limited zone not well defined)
        // NOTE: this cut depends on your MC statistics, not on the sensitivity of your observatory
        bool bPasses_MinimumNumberofBackGroundEvents = ( iBackground * iAlpha > 0. );
        // require the signal to be larger than a certain fraction of background
        bool bPasses_MinimumSystematicCut = false;
        if( iBackground * iAlpha > 0. )
        {
            bPasses_MinimumSystematicCut = ( f * n_diff / ( iBackground * iAlpha ) >= fMinBackgroundRateRatio_min );
        }
        
        // PRELI: allow calculation of sensitivity in event limited region
        if( bPasses_MinimumNumberofBackGroundEvents )
        {
            fSetEvents_minCutOnly = false;
        }
        if( fSetEvents_minCutOnly )
        {
            bPasses_MinimumNumberofBackGroundEvents = true;
            bPasses_MinimumSystematicCut = true;
            bPassed_MinimumSignificance = true;
        }
        
        // sensitivity limitation histograms
        if( iFillStatistics != 0 )
        {
            if( ( iFillStatistics == 1 && bPassed_MinimumSignificance )
                    || ( iFillStatistics == 2 && bPassed_MinimumSignalEvents )
                    || ( iFillStatistics == 3 && bPasses_MinimumSystematicCut ) )
            {
                f = fSourceStrength[n];
                bSuccess = true;
                break;
            }
            continue;
        }
        
        // standard sensitivity calculation
        if( bPassed_MinimumSignificance && bPassed_MinimumSignalEvents
                && bPasses_MinimumSystematicCut
                && bPasses_MinimumNumberofBackGroundEvents )
        {
            f = fSourceStrength[n];
            bSuccess = true;
            break;
        }
        if( fDebug && bSuccess && energy > 0. )
        {
            cout << fixed;
            cout << "\t SignificanceCalculation: ";
            cout << "\t n: " << n - 1 << "\t f " << scientific << f;
            cout << fixed << "\t significance: " << s;
            cout << "\t min events: " << t* f* iSignal;
            cout << "\t systematics: " << ( f * n_diff / ( iBackground * iAlpha ) );
            cout << endl;
            cout << "\t\t ndiff: " << t * ( f * n_diff );
            cout << "\t non: " << t * ( f * n_diff + iBackground * iAlpha );
            cout << "\t noff: " << t* iBackground;
            cout << "\t alpha: " << iAlpha;
            cout << "\t t: " << t;
            cout << endl;
        }
    }
    
    if( !bSuccess )
    {
        return -1.;
    }
    
    // return flux value in CU that passed significance criteria
    return f;
}


bool VSensitivityCalculator::checkDataSet( unsigned int iD, string iName )
{
    if( iD >= fData.size() )
    {
        cout << "ERROR (" << iName << "): data set ID out of range: should be < " << fData.size() << endl;
        return false;
    }
    return true;
}


unsigned int VSensitivityCalculator::listDataSets()
{
    cout << "Available data sets " << fData.size() << ": " << endl;
    cout << "==================== " << endl;
    cout << "Signal rate [1/min] \t Background rate [1/min] \t Alpha \t Name " << endl;
    cout << "--------------------\t-------------------------\t-------\t------" << endl;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        cout << fData[i].fSignal << "\t\t\t";
        cout << fData[i].fBackground << "\t\t\t";
        cout << fData[i].fAlpha << "\t\t";
        cout << fData[i].fName << endl;
    }
    cout << endl;
    
    return fData.size();
}


bool VSensitivityCalculator::setCurrentDataSet( unsigned int iD )
{
    if( !checkDataSet( iD, "setCurrentDataSet" ) )
    {
        return false;
    }
    
    fCurrentDataSet = iD;
    
    return true;
}


bool VSensitivityCalculator::removeDataSet( unsigned int iD )
{
    if( !checkDataSet( iD, "removeDataSet" ) )
    {
        return false;
    }
    
    fData.erase( fData.begin() + iD );
    
    return true;
}


/*!
     iGammaRayRate = signal rate [1/min]
     iBackGroundRate = background rate [1/min]
     iAlpha = normalization parameter (e.g. 1./5.)
     iName  = name of data set
*/
unsigned int VSensitivityCalculator::addDataSet( double iGammaRayRate, double iBackGroundRate, double iAlpha, string iName )
{
    sSensitivityData t;
    
    t.fSignal = iGammaRayRate;
    t.fBackground = iBackGroundRate;
    t.fAlpha = iAlpha;
    t.fName = iName;
    
    fData.push_back( t );
    
    return fData.size() - 1;
}

/*!
     iOnRate = signal rate [1/min]
     iOffRate = background rate [1/min]
     iAlpha = normalization parameter (e.g. 1./5.)
     iName  = name of data set
*/
unsigned int VSensitivityCalculator::addDataSetOnOff( double iOnRate, double iOffRate, double iAlpha, string iName )
{
    sSensitivityData t;
    
    t.fSignal = iOnRate - iOffRate * iAlpha;
    t.fBackground = iOffRate;
    t.fAlpha = iAlpha;
    t.fName = iName;
    
    fData.push_back( t );
    
    return fData.size() - 1;
}


/*!

    fSourceStrength will be sorted in reverse order

*/
void VSensitivityCalculator::setSourceStrengthRange_CU( double iMin, double iMax, double iStep, bool iLog )
{
    fSourceStrength.clear();
    
    fSourceStrength_min = iMin;
    fSourceStrength_max = iMax;
    fSourceStrength_step = iStep;
    
    unsigned int i_Steps = ( unsigned int )( ( fSourceStrength_max - fSourceStrength_min ) / fSourceStrength_step ) + 1;
    for( unsigned int n = 0; n < i_Steps; n++ )
    {
        fSourceStrength.push_back( fSourceStrength_max - n * fSourceStrength_step );
    }
    
    if( iLog )
    {
        for( unsigned int i = 0; i < fSourceStrength.size(); i++ )
        {
            fSourceStrength[i] = TMath::Power( 10., fSourceStrength[i] );
        }
    }
}

/*

    set the source strength vector and sort it in reverse order

*/
void VSensitivityCalculator::setSourceStrengthVector_CU( vector< double > iF )
{
    fSourceStrength = iF;
    
    // sort in reverse order
    sort( fSourceStrength.begin(), fSourceStrength.end() );
    reverse( fSourceStrength.begin(), fSourceStrength.end() );
}

TGraph* VSensitivityCalculator::getCrabSpectrum( bool bInt, string bUnit, bool bReset, bool iUseMAGICCurvedCrabSpectrum )
{
    vector< double > iTemp;
    iTemp.push_back( 1. );
    
    vector< TGraph* > i_Temp_Graph = getCrabSpectrum( iTemp, bInt, bUnit, bReset, iUseMAGICCurvedCrabSpectrum );
    if( i_Temp_Graph.size() == 1 )
    {
        return i_Temp_Graph[0];
    }
    
    return 0;
}

/*!
     fill graphs with Crab-like spectra and different source strengths

     i_fCrabFlux:  vector of Crab fluxes in CU units [1. = 1 Crab ]
     bInt:         integrated or differential flux
     bUnit:        PFLUX, CU, ENERGY
     bReset:       recalculate the graphs (otherwise: graph is not filled if the number of graph requested is the same as CrabFlux_SourceStrength.size() )

*/
vector< TGraph* > VSensitivityCalculator::getCrabSpectrum( vector< double > i_fCrabFlux, bool bInt, string bUnit, bool bReset, bool iUseMAGICCurvedCrabSpectrum )
{
    if( fDebug )
    {
        cout << "VSensitivityCalculator::getCrabSpectrum " << i_fCrabFlux.size();
        cout << "(INT " << bInt << ", UNIT " << bUnit << ")";
        cout << endl;
    }
    
    // check if Crab spectrum is already defined
    if( i_fCrabFlux.size() == fCrabFlux_SourceStrength.size() && !bReset )
    {
        return fCrabFlux_SourceStrength;
    }
    
    // temporary function describing Crab spectrum
    TF1* i_fFunCrabFlux = 0;
    // temporary graph
    TGraph* i_GraphCrabFlux = new TGraph( 100 );
    
    //////////////////////////////////////
    // use MAGIC curved Crab spectrum
    // - preliminary -
    // (integrated spectra do not work yet)
    if( iUseMAGICCurvedCrabSpectrum && bUnit == "ENERGY" && !bInt )
    {
        //Crab Nebula (MAGIC paper) with [0]=1.
        // code from Tarek H.
        TString specstrDef = "[0]*(6.0e-10)*pow((x/(0.3*1E3)),(-2.31+(-0.26*TMath::Log10((x)/(0.3*1E3)))))";
        TF1* origSpec = new TF1( "origSpec", specstrDef, 1E1, 1E6 );
        origSpec->SetNpx( 10000 );
        origSpec->FixParameter( 0, 1. );
        
        TF1* origSpecE2 = new TF1( "origSpecE2", "origSpec*(x*1E-3)*(x*1E-3)*1.602176487", 1E1, 1E6 );
        origSpecE2->SetNpx( 10000 );
        double iE = 0.;
        double iY = 0.;
        // loop over different Crab fluxes
        for( int p = 0; p < i_GraphCrabFlux->GetN(); p++ )
        {
            // equal interval in logE
            iE = fEnergy_min_Log + p * ( fEnergy_max_Log - fEnergy_min_Log ) / i_GraphCrabFlux->GetN();
            // differential flux
            iY = origSpecE2->Eval( TMath::Power( 10., iE ) * 1.e3 );
            i_GraphCrabFlux->SetPoint( p, iE, iY );
        }
        
    }
    ////////////////////////////////////////////////////////
    // Crab spectra as simple power law for all energies
    // (better: use spectra from external text file)
    else
    {
        if( fEnergySpectrumfromLiterature )
        {
            if( bUnit == "PFLUX" || bUnit == "ENERGY" )
            {
                i_fFunCrabFlux = fEnergySpectrumfromLiterature->getEnergySpectrum( fEnergySpectrumfromLiterature_ID,
                                 false, TMath::Power( 10., fEnergy_min_Log ), 10000. );
            }
            else if( bUnit == "CU" )
            {
                i_fFunCrabFlux = new TF1( "i_fFunCrabFlux" , "1.", TMath::Power( 10., fEnergy_min_Log ), 10000. );
            }
        }
        if( !i_fFunCrabFlux )
        {
            cout << "Error reading Crab Nebula spectrum " << fEnergySpectrumfromLiterature << "\t" << fEnergySpectrumfromLiterature_ID << endl;
            vector< TGraph* > xx;
            return xx;
        }
        /////////////////////////////////////////////////////
        // calculate integral spectrum
        if( bInt && bUnit != "CU" )
        {
            const int np = 1000;
            double* x = new double[np];
            double* y = new double[np];
            double iE = 0.;
            double iY = 0.;
            
            i_fFunCrabFlux->CalcGaussLegendreSamplingPoints( np, x, y, 1.e-15 );
            // loop over different Crab flux
            for( int p = 0; p < i_GraphCrabFlux->GetN(); p++ )
            {
                // equal interval in logE
                iE = fEnergy_min_Log + p * ( fEnergy_max_Log - fEnergy_min_Log ) / i_GraphCrabFlux->GetN();
                // integrate above this energy
                iY = i_fFunCrabFlux->IntegralFast( np, x, y, TMath::Power( 10., iE ), 10000. );
                if( bUnit == "ENERGY" )
                {
                    iY *= 1.e12 * fConstant_Flux_To_Ergs * TMath::Power( 10., iE );
                }
                i_GraphCrabFlux->SetPoint( p, iE, iY );
            }
            delete[] x;
            delete[] y;
        }
        /////////////////////////////////////////////////////
        // calculate differential spectrum
        else
        {
            double iE = 0.;
            double iY = 0.;
            // loop over different Crab fluxes
            for( int p = 0; p < i_GraphCrabFlux->GetN(); p++ )
            {
                // equal interval in logE
                iE = fEnergy_min_Log + p * ( fEnergy_max_Log - fEnergy_min_Log ) / i_GraphCrabFlux->GetN();
                // differential flux
                iY = i_fFunCrabFlux->Eval( TMath::Power( 10., iE ) );
                if( bUnit == "ENERGY" )
                {
                    iY *= 1.e12 * fConstant_Flux_To_Ergs * TMath::Power( 10., iE ) * TMath::Power( 10., iE );
                }
                i_GraphCrabFlux->SetPoint( p, iE, iY );
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////
    // now fill a vector with graphs for the different flux values needed
    double xx = 0.;
    double yy = 0.;
    fCrabFlux_SourceStrength.clear();
    for( unsigned int i = 0; i < i_fCrabFlux.size(); i++ )
    {
        fCrabFlux_SourceStrength.push_back( new TGraph( i_GraphCrabFlux->GetN() ) );
        for( int p = 0; p < i_GraphCrabFlux->GetN(); p++ )
        {
            i_GraphCrabFlux->GetPoint( p, xx, yy );
            fCrabFlux_SourceStrength.back()->SetPoint( p, xx, yy * i_fCrabFlux[i] );
        }
        fCrabFlux_SourceStrength.back()->SetLineStyle( 9 );
    }
    if( fDebug )
    {
        cout << "VSensitivityCalculator::getCrabSpectrum " << fCrabFlux_SourceStrength.size() << endl;
    }
    
    return fCrabFlux_SourceStrength;
}


/*

    calculate and plot integral sensitivity

*/
TCanvas* VSensitivityCalculator::plotIntegralSensitivityvsEnergyFromCrabSpectrum( TCanvas* cSensitivity, string iAnasumCrabFile,
        string bUnit,
        double iEnergyMin_TeV_lin, double iEnergyMax_TeV_lin )
{
    fRequireCutsToBeOptimized = false;
    if( fPlotDebugName.size() > 0 )
    {
        prepareDebugPlots();
    }
    if( !calculateSensitivityvsEnergyFromCrabSpectrum( iAnasumCrabFile, bUnit, -1., iEnergyMin_TeV_lin, iEnergyMax_TeV_lin ) )
    {
        return 0;
    }
    return plotSensitivityvsEnergyFromCrabSpectrum( cSensitivity, bUnit, -1. );
}


/*

    calculate and plot differential sensitivity

*/
TCanvas* VSensitivityCalculator::plotDifferentialSensitivityvsEnergyFromCrabSpectrum( TCanvas* cSensitivity, string iAnasumCrabFile,
        string bUnit,
        double dE_Log10,
        double iEnergyMin_TeV_lin, double iEnergyMax_TeV_lin )
{
    fRequireCutsToBeOptimized = false;
    if( fPlotDebugName.size() > 0 )
    {
        prepareDebugPlots();
    }
    if( !calculateSensitivityvsEnergyFromCrabSpectrum( iAnasumCrabFile, bUnit, dE_Log10, iEnergyMin_TeV_lin, iEnergyMax_TeV_lin ) )
    {
        return 0;
    }
    
    return plotSensitivityvsEnergyFromCrabSpectrum( cSensitivity, bUnit, dE_Log10 );
}

bool VSensitivityCalculator::calculateParticleNumberGraphs_MC( double dE_Log10 )
{
    double alpha = 1.;
    vector< VDifferentialFluxData > i_DifferentialFlux = getDifferentialFluxVectorfromMC( dE_Log10, alpha );
    fillParticleNumbersGraphs( i_DifferentialFlux, alpha, dE_Log10 );
    
    return true;
}


/*

    calculate and plot sensitivities

    integral sensitivity:     dE_Log10 < 0
    differential sensitivity: dE_Log10 = bin size in log10 E


    bNewCanvas:   plot histograms and Crab lines

*/
TCanvas* VSensitivityCalculator::plotSensitivityvsEnergyFromCrabSpectrum( 
        TCanvas* cSensitivity, string bUnit, double dE_Log10, bool bNewCanvas )
{

    // get canvas
    if( cSensitivity == 0 )
    {
        cSensitivity = plotCanvas_SensitivityvsEnergy( bUnit, ( dE_Log10 < 0. ) );
    }
    else if( bNewCanvas )
    {
        cSensitivity->cd();
        cSensitivity = plotCanvas_SensitivityvsEnergy( bUnit, ( dE_Log10 < 0. ), bNewCanvas );
    }
    if( !cSensitivity )
    {
        return 0;
    }
    cSensitivity->cd();
    if( !gSensitivityvsEnergy )
    {
        return cSensitivity;
    }
    
    // plot everything
    setGraphPlottingStyle( gSensitivityvsEnergy, fPlottingColor, fPlottingLineWidth,
                           fPlottingMarkerStyle, fPlottingMarkerSize,
                           fPlottingFillStyle, fPlottingLineStyle );
                           
    gSensitivityvsEnergy->Draw( "p" );
    
    cSensitivity->Update();
    
    return cSensitivity;
    
}


bool VSensitivityCalculator::calculateSensitivityvsEnergyFromCrabSpectrum( string iAnasumCrabFile, string bUnit, double dE_Log10,
        double iEnergyMin_TeV_lin, double iEnergyMax_TeV_lin )
{
    if( !checkUnits( bUnit ) )
    {
        cout << "VSensitivityCalculator::calculateSensitivityvsEnergyFromCrabSpectrum: error while checking units: " << bUnit << endl;
        return 0;
    }
    cout << endl;
    if( dE_Log10 > 0. )
    {
        cout << "calculating or reading differential sensitivity";
    }
    else
    {
        cout << "calculating or reading integral sensitivity";
    }
    cout << " (" << iAnasumCrabFile << ")";
    cout << endl;
    cout << "==============================================" << endl;
    fEnergy_dE_log10 = dE_Log10;
    
    // get vector with integral and differential Crab-like spectra for different flux levels
    // (use Whipple spectrum)
    TGraph* i_fFunCrabFlux = getCrabSpectrum( ( dE_Log10 < 0. ), bUnit, true, false );
    if( !i_fFunCrabFlux )
    {
        cout << "VSensitivityCalculator::calculateSensitivityvsEnergyFromCrabSpectrum: error: no reference spectrum found " << endl;
        return 0;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // get differential flux vector (used for differential and integral flux calculations)
    // differential flux is calculated for 1 Crab Unit
    
    // note: alpha is modified by getDifferentialFlux...
    double alpha = 1.;
    vector< VDifferentialFluxData > fDifferentialFlux;
    // differential flux vector from effectivea areas (gamma/electron/proton...)
    if( iAnasumCrabFile == "MC" )
    {
        fDifferentialFlux = getDifferentialFluxVectorfromMC( dE_Log10, alpha );
    }
    // differential flux vector from effectivea areas / background from from MC phys file
    else if( iAnasumCrabFile == "CTA-MC" )
    {
        fDifferentialFlux = getDifferentialFluxVectorfromCTA_MC( alpha );
    }
    // read sensitivity directly from CTA WP Phys root file and return
    else if( iAnasumCrabFile == "CTA-PHYS" )
    {
        gSensitivityvsEnergy = getSensitivityGraphFromWPPhysFile( bUnit, iEnergyMin_TeV_lin, iEnergyMax_TeV_lin, dE_Log10 );
        if( !gSensitivityvsEnergy )
        {
            return false;
        }
        return true;
    }
    // differential flux vector from anasum file (e.g. a measurement towards the Crab Nebula)
    else
    {
        fDifferentialFlux = getDifferentFluxVectorfromData( iAnasumCrabFile, dE_Log10, alpha );
    }
    
    if( fDifferentialFlux.size() == 0 )
    {
        cout << "Error: no entries in differential flux vector" << endl;
        return 0;
    }
    cout << "evaluating differential flux vector (size " << fDifferentialFlux.size() << ")" << endl;
    
    // this is the range of fluxes which will be searched
    setSourceStrengthRange_CU( -6., 10., 4. / 1000., true );
    
    double non = 0;
    double non_error = 0.;
    double noff = 0;
    double noff_error = 0.;
    double s = 0.;
    double s_error_U = 0.;
    double s_error_L = 0.;
    
    // create sensitivity vs energy graph
    gSensitivityvsEnergy = new TGraphAsymmErrors( 1 );
    gSensitivityvsEnergy->SetTitle( "" );
    
    // often no background simulations (off events) available at high energies
    // don't really need it, since fEvents_min is the important one
    fSetEvents_minCutOnly = true;
    
    /////////////////////////////////////////////////////////////////
    // loop over all on and off counts and calculate sensitivity
    int z = 0;
    for( int i = fDifferentialFlux.size() - 1; i >= 0; i-- )
    {
        cout << fixed << "FLUX RESULTS " << z << ":\t" << fDifferentialFlux[i].EnergyWeightedMean << " [TeV]\t";
        cout << fDifferentialFlux[i].Energy << " [TeV] ";
        if( dE_Log10 > 0. )
        {
            cout << " dE_log10: " << fDifferentialFlux[i].dE << "\t";
        }
        cout << "[" << fDifferentialFlux[i].Energy_lowEdge << ", " << fDifferentialFlux[i].Energy_upEdge << "]\t";
        cout << "[" << log10( fDifferentialFlux[i].Energy_lowEdge ) << ", " << log10( fDifferentialFlux[i].Energy_upEdge ) << "]\t";
        cout << endl;
        
        // integral spectrum: add (integrate over) number of events
        if( dE_Log10 < 0. && fDifferentialFlux[i].Energy < iEnergyMax_TeV_lin )
        {
            non        += fDifferentialFlux[i].NOn;
            noff       += fDifferentialFlux[i].NOff;
            non_error  += fDifferentialFlux[i].NOn_error * fDifferentialFlux[i].NOn_error;
            noff_error += fDifferentialFlux[i].NOff_error * fDifferentialFlux[i].NOff_error;
        }
        // differential spectrum: number of events per differential bin
        else
        {
            non        = fDifferentialFlux[i].NOn;
            noff       = fDifferentialFlux[i].NOff;
            non_error  = fDifferentialFlux[i].NOn_error * fDifferentialFlux[i].NOn_error;
            noff_error = fDifferentialFlux[i].NOff_error * fDifferentialFlux[i].NOff_error;
        }
        non_error  = sqrt( non_error );
        noff_error = sqrt( noff_error );
        
        // skip all zero entries
        if( TMath::Abs( non ) < 1.e-2 && TMath::Abs( noff ) < 1.e-2 )
        {
            continue;
        }
        // set errors to zero for no off events
        if( TMath::Abs( noff ) < 1e-2 )
        {
            noff = 0.;
            noff_error = 0.;
        }
        
        // observe that the signal rate is defined differently for list_sensitivity() etc...
        // get fraction of Crab flux for XX sigma etc...
        // perform a toy MC to calculate errors (i_N_iter times)
        //        const unsigned int i_N_iter = 1000;
        // NOTE: this is not the correct approach:
        // Non is randomized using 100% Crab flux, not f*CU
        const unsigned int i_N_iter = 100;
        double i_s_v[i_N_iter];
        double i_s_x = 0.;
        double i_s_xx = 0.;
        int i_s_z = 0;
        for( unsigned int q = 0; q < i_N_iter; q++ )
        {
            double iN_on  = gRandom->Gaus( non, non_error );
            double iN_off = gRandom->Gaus( noff, noff_error );
            // in case of no background events: set background rate to very very small number
            // (otherwise no sensitivity values is calculated)
            if( iN_off <= 0. )
            {
                iN_off = 1.e-15;
            }
            double i_s = getSensitivity( iN_on  / fDifferentialFlux[i].ExposureTime * 60.,
                                         iN_off / fDifferentialFlux[i].ExposureTime * 60., alpha , -100., 0 );
            if( i_s > 0 )
            {
                i_s_v[i_s_z] = i_s;
                i_s_x   += i_s;
                i_s_xx  += i_s * i_s;
                i_s_z++;
            }
        }
        // take median of distribution (distribution is asymmetric)
        if( i_s_z > 1 )
        {
            double i_a[] = { 0.16, 0.5, 0.84 };
            double i_b[] = { 0.0,  0.0, 0.0  };
            TMath::Quantiles( i_s_z, 3, i_s_v, i_b, i_a, kFALSE );
            cout << "\t Quantiles " << i_b[1] << "\t" << i_b[1] - i_b[0] << "\t" << i_b[2] - i_b[1] << "\t" << i_s_z << endl;
            cout << "\t Mean      " << i_s_x / i_s_z;
            cout << ", RMS: " << sqrt( 1. / ( i_s_z - 1. ) * ( i_s_xx - i_s_x * i_s_x / i_s_z ) );
            if( i_s_x / i_s_z > 0. )
            {
                cout << ", RMS/Mean " << sqrt( 1. / ( i_s_z - 1. ) * ( i_s_xx - i_s_x * i_s_x / i_s_z ) ) / ( i_s_x / i_s_z );
            }
            cout << endl;
            s = i_b[1];
            s_error_L = i_b[0];
            s_error_U = i_b[2];
        }
        else
        {
            s = getSensitivity( non / fDifferentialFlux[i].ExposureTime * 60., noff / fDifferentialFlux[i].ExposureTime * 60., alpha, -100., 0 );
            s_error_L = 0.;
            s_error_U = 0.;
        }
        
        // Preliminary: catch cases were lower sensitivity cannnot be calculated
        if( s_error_L < 0. )
        {
            s_error_L = s;
            s_error_U = s;
        }
        // weighted energy of this bin
        double energy = TMath::Log10( fDifferentialFlux[i].EnergyWeightedMean );
        // (int) of energy for maps
        int energy_intX3 = ( int )( energy * 1.e3 );
        // get sensitivity limititations
        fSignificanceLimited[energy_intX3]        = getSensitivity( non / fDifferentialFlux[i].ExposureTime * 60.,
                noff / fDifferentialFlux[i].ExposureTime * 60., alpha, energy, 1 );
        fMinEventsLimited[energy_intX3]           = getSensitivity( non / fDifferentialFlux[i].ExposureTime * 60.,
                noff / fDifferentialFlux[i].ExposureTime * 60., alpha, energy, 2 );
        fMinBackgroundEventsLimited[energy_intX3] = getSensitivity( non / fDifferentialFlux[i].ExposureTime * 60.,
                noff / fDifferentialFlux[i].ExposureTime * 60., alpha, energy, 3 );
        fMinNoBackground[energy_intX3]            = getSensitivity( non / fDifferentialFlux[i].ExposureTime * 60.,
                noff / fDifferentialFlux[i].ExposureTime * 60., alpha, energy, 4 );
                
        // fill sensitivity graphs
        double f1 = i_fFunCrabFlux->Eval( log10( fDifferentialFlux[i].Energy_lowEdge ) );
        double f2 = i_fFunCrabFlux->Eval( log10( fDifferentialFlux[i].Energy_upEdge ) );
        
        if( i_fFunCrabFlux != 0 && s > 0.
                && fDifferentialFlux[i].Energy > iEnergyMin_TeV_lin && fDifferentialFlux[i].Energy < iEnergyMax_TeV_lin
                && checkCutOptimization( fDifferentialFlux[i].Energy ) )
        {
            // integral sensitivity
            if( dE_Log10 < 0. )
            {
                if( bUnit == "CU" )
                {
                    gSensitivityvsEnergy->SetPoint( z, energy, s );
                    gSensitivityvsEnergy->SetPointEYhigh( z, TMath::Abs( s - s_error_U ) );
                    gSensitivityvsEnergy->SetPointEYlow( z, TMath::Abs( s - s_error_L ) );
                }
                else
                {
                    gSensitivityvsEnergy->SetPoint( z, energy, i_fFunCrabFlux->Eval( energy ) * s );
                    gSensitivityvsEnergy->SetPointEYhigh( z, i_fFunCrabFlux->Eval( energy ) * TMath::Abs( s - s_error_U ) );
                    gSensitivityvsEnergy->SetPointEYlow( z, i_fFunCrabFlux->Eval( energy ) * TMath::Abs( s - s_error_L ) );
                }
            }
            // differential sensitivity
            else
            {
                if( bUnit == "PFLUX" )
                {
                    gSensitivityvsEnergy->SetPoint( z, energy, ( f1 - f2 ) * s / fDifferentialFlux[i].dE );
                }
                else if( bUnit == "CU" )
                {
                    gSensitivityvsEnergy->SetPoint( z, energy, s );
                    gSensitivityvsEnergy->SetPointEYhigh( z, TMath::Abs( s - s_error_U ) );
                    gSensitivityvsEnergy->SetPointEYlow( z, TMath::Abs( s - s_error_L ) );
                }
                else if( bUnit == "ENERGY" )
                {
                    gSensitivityvsEnergy->SetPoint( z, energy, s * i_fFunCrabFlux->Eval( energy ) );
                    gSensitivityvsEnergy->SetPointEYhigh( z, TMath::Abs( s - s_error_U ) * i_fFunCrabFlux->Eval( energy ) );
                    gSensitivityvsEnergy->SetPointEYlow( z, TMath::Abs( s - s_error_L ) * i_fFunCrabFlux->Eval( energy ) );
                }
                else
                {
                    cout << "plotSensitivityvsEnergyFromCrabSpectrum error: unknown unit: " << bUnit << endl;
                    return 0;
                }
            }
            z++;
        }
        else
        {
            cout << "ENERGY DISCARDED: " << fDifferentialFlux[i].Energy << " [TeV], [" << iEnergyMin_TeV_lin;
            cout << ", " << iEnergyMax_TeV_lin << "]";
            cout << ", Sens " << s << " + " << s_error_U << " - " << s_error_L << " (" << i_s_z << ")";
            cout << ", FunCrab: " << i_fFunCrabFlux << endl;
            cout << "\t significance: " << fSignificanceLimited[energy_intX3] << endl;
            cout << "\t events: " << fMinEventsLimited[energy_intX3] << endl;
            cout << "\t systematics: " << fMinBackgroundEventsLimited[energy_intX3] << endl;
            cout << "\t no background: " << fMinNoBackground[energy_intX3] << endl;
            cout << "\t cut optimization: " << checkCutOptimization( fDifferentialFlux[i].Energy, true )  << endl;
        }
        // sensitivity limitiations
        // (TMP differential sensitivity only)
        if( dE_Log10 > 0. )
        {
            if( bUnit == "ENERGY" )
            {
                fSignificanceLimited[energy_intX3]        *= i_fFunCrabFlux->Eval( energy );
                fMinEventsLimited[energy_intX3]           *= i_fFunCrabFlux->Eval( energy );
                fMinBackgroundEventsLimited[energy_intX3] *= i_fFunCrabFlux->Eval( energy );
                fMinNoBackground[energy_intX3]            *= i_fFunCrabFlux->Eval( energy );
            }
            else if( bUnit == "PFLUX" )
            {
                fSignificanceLimited[energy_intX3]        *= ( f1 - f2 ) / fDifferentialFlux[i].dE;
                fMinEventsLimited[energy_intX3]           *= ( f1 - f2 ) / fDifferentialFlux[i].dE;
                fMinBackgroundEventsLimited[energy_intX3] *= ( f1 - f2 ) / fDifferentialFlux[i].dE;
                fMinNoBackground[energy_intX3]            *= ( f1 - f2 ) / fDifferentialFlux[i].dE;
            }
        }
        // print out
        cout << "\t CU RESULTS: NON: " << non << " (+-" << non_error << ")";
        cout << "\t NOFF: " << noff << " (+-" << noff_error << ")";
        cout << "\t NDIFF: " << non - alpha* noff  << ", alpha=" << alpha << "\t";
        cout << fDifferentialFlux[i].ExposureTime << " [s]\t";
        cout << endl;
        float i_ndiff = s * ( non - alpha * noff );
        cout << "\t SENS RESULTS: NON: " << alpha* noff + i_ndiff;
        cout << "\t NOFF: " << noff;
        cout << "\t NDIFF: " << i_ndiff;
        cout << endl;
        
        cout << "\t FLUX: " <<  setprecision( 4 ) << s << " (" << s_error_L << "," << s_error_U << ") [CU],  ";
        if( bUnit == "PFLUX" )
        {
            cout <<  i_fFunCrabFlux->Eval( energy ) * s << " [cm^-2s^-1] (" << energy << ")";
        }
        else if( bUnit == "ENERGY" )
        {
            cout << scientific <<  i_fFunCrabFlux->Eval( energy ) * s;
            cout << " (" << TMath::Abs( s - s_error_U ) * i_fFunCrabFlux->Eval( energy );
            cout << ", " << TMath::Abs( s - s_error_L ) * i_fFunCrabFlux->Eval( energy ) << ")";
            cout << " [erg cm^-2s^-1]";
        }
        cout << " (" << fSetEvents_minCutOnly << ")";
        cout << endl;
        cout << "\t Limits:";
        cout << "\t significance: " << fSignificanceLimited[energy_intX3];
        cout << "\t events: " << fMinEventsLimited[energy_intX3];
        cout << "\t systematics: " << fMinBackgroundEventsLimited[energy_intX3];
        cout << "\t no background: " << fMinNoBackground[energy_intX3] << endl;
        cout << "\t cut optimization: " << checkCutOptimization( fDifferentialFlux[i].Energy, true )  << endl;
    }
    if( fDebug )
    {
        gSensitivityvsEnergy->Print();
    }
    
    fillParticleNumbersGraphs( fDifferentialFlux, alpha );
    if( fPlotDebugName.size() > 0 )
    {
        plotDebugPlotsParticleNumbers();
    }
    
    return true;
}

bool VSensitivityCalculator::printSensitivity()
{
    if( gSensitivityvsEnergy )
    {
        cout << "Contents of sensitivity graph: " << endl;
        gSensitivityvsEnergy->Print();
        return true;
    }
    
    return false;
}


/*!
    plot an empty canvas for the sensitivy vs energy graphs

    add lines corresponding to XXX% of the Crab

*/
TCanvas* VSensitivityCalculator::plotCanvas_SensitivityvsEnergy( string bUnit, bool bIntegralSensitivity, bool bNewCanvas )
{
    char htitle[400];
    TCanvas* iCanvas = 0;
    if( !bNewCanvas )
    {
        iCanvas = new TCanvas( fPlot_CanvasName.c_str(), fPlot_CanvasTitle.c_str(), 10, 10, fPlot_CanvasSize_x, fPlot_CanvasSize_y );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        iCanvas->SetLeftMargin( 0.15 );
    }
    else
    {
        iCanvas = ( TCanvas* )gPad;
    }
    if( !iCanvas )
    {
        return 0;
    }
    
    hnull = new TH1D( "hnullSens", "", 10, fEnergy_min_Log, fEnergy_max_Log );
    hnull->SetStats( 0 );
    // integral sensitivity
    if( bIntegralSensitivity )
    {
        hnull->SetXTitle( "log_{10} energy E_{t} [TeV]" );
        if( bUnit == "PFLUX" )
        {
            hnull->SetYTitle( "Integral Flux Sensitivity (E>E_{t}) [cm^{-2} s^{-1}]" );
            hnull->SetMaximum( fPlot_flux_PFLUX_max );
            hnull->SetMinimum( fPlot_flux_PFLUX_min );
        }
        else if( bUnit == "ENERGY" )
        {
            hnull->SetYTitle( "E_{t} x Integral Flux Sensitivity (E>E_{t}) [erg cm^{-2} s^{-1}]" );
            hnull->SetMaximum( fPlot_flux_ENERG_max );
            hnull->SetMinimum( fPlot_flux_ENERG_min );
        }
        else if( bUnit == "CU" )
        {
            hnull->SetYTitle( "Integral Flux Sensitivity (E>E_{t}) [C.U.]" );
            hnull->SetMaximum( fPlot_flux_CU_max );
            hnull->SetMinimum( fPlot_flux_CU_min );
        }
    }
    // differential sensitivity
    else
    {
    
        hnull->SetXTitle( "log_{10} energy [TeV]" );
        if( bUnit == "PFLUX" )
        {
            hnull->SetYTitle( "Flux Sensitivity [cm^{-2} s^{-1} TeV^{-1}]" );
            hnull->SetMaximum( fPlot_flux_PFLUX_max );
            hnull->SetMinimum( fPlot_flux_PFLUX_min );
        }
        else if( bUnit == "ENERGY" )
        {
            hnull->SetYTitle( "E^{2} x Flux Sensitivity [erg cm^{-2} s^{-1}]" );
            hnull->SetMaximum( fPlot_flux_ENERG_max );
            hnull->SetMinimum( fPlot_flux_ENERG_min );
        }
        else if( bUnit == "CU" )
        {
            hnull->SetYTitle( "Differential Flux Sensitivity  [C.U.]" );
            hnull->SetMaximum( fPlot_flux_CU_max );
            hnull->SetMinimum( fPlot_flux_CU_min );
        }
    }
    
    plot_nullHistogram( iCanvas, hnull, true, true, 1.7, TMath::Power( 10., fEnergy_min_Log ), TMath::Power( 10., fEnergy_max_Log ) );
    iCanvas->SetLogy( 1 );

    // get vector with integral or differential Crab-like spectra for different flux levels
    vector< TGraph* > i_fFunCrabFlux = getCrabSpectrum( fPlottingCrabFlux_CU, bIntegralSensitivity, bUnit, false, true );
    
    for( unsigned int i = 0; i < fPlottingCrabFlux_CU.size(); i++ )
    {
        if( i < i_fFunCrabFlux.size() && i_fFunCrabFlux[i] )
        {
            if( !bPlotCrabLines )
            {
                continue;
            }
            i_fFunCrabFlux[i]->SetLineColor( 16 );
            i_fFunCrabFlux[i]->Draw( "l" );
            
            if( bUnit != "CU" )
            {
                if( fPlottingCrabFlux_CU[i] > 1. )
                {
                    sprintf( htitle, "%dx Crab", ( int )( fPlottingCrabFlux_CU[i] ) );
                }
                else if( fPlottingCrabFlux_CU[i] * 100. >= 1. )
                {
                    sprintf( htitle, "%d%% Crab", ( int )( fPlottingCrabFlux_CU[i] * 100 ) );
                }
                else if( fPlottingCrabFlux_CU[i] * 100. > 0.09 )
                {
                    sprintf( htitle, "%.1f%% Crab", fPlottingCrabFlux_CU[i] * 100 );
                }
                else if( fPlottingCrabFlux_CU[i] * 100. > 0.009 )
                {
                    sprintf( htitle, "%.2f%% Crab", fPlottingCrabFlux_CU[i] * 100 );
                }
                double e = 1. - 0.3 * ( double )i;
                if( hnull && i_fFunCrabFlux[i]->Eval( e ) * 1.1 < hnull->GetMaximum() )
                {
                    TText* iT = new TText( e, i_fFunCrabFlux[i]->Eval( e ) * 1.1, htitle );
                    iT->SetTextColor( 16 );
                    iT->SetTextSize( 0.3 * iT->GetTextSize() );
                    iT->Draw();
                }
            }
        }
    }
    iCanvas->Update();
    
    return iCanvas;
}


/*!
    units used in data files should be encoded into the file name

    CU  = Crab units
    PFLUX = particle flux [1/cm2/s]
    ENERG = energy flux [erg/cm2/s]

    energy in TeV

*/
TCanvas* VSensitivityCalculator::plotSensitivityvsEnergyFromTextTFile( TCanvas* c, string iTextFile, int iColor, double iLineWidth, int iLineStyle,
        string bUnit, bool bIntegralSensitivity, double iMultFactor )
{
    // check if canvas exists, otherwise create a new one
    if( c == 0 )
    {
        c = plotCanvas_SensitivityvsEnergy( bUnit, bIntegralSensitivity );
        if( !c )
        {
            return 0;
        }
    }
    c->cd();
    
    // graph with sensitivity values: read them directly from the text file
    TGraph* g = new TGraph( gSystem->ExpandPathName( iTextFile.c_str() ) );
    cout << "reading " << iTextFile << " with " << g->GetN() << " data points." << endl;
    if( g->IsZombie() || g->GetN() < 1 )
    {
        cout << "VSensitivityCalculator::plotSensitivityvsEnergyFromTextTFile warning: ";
        cout << " no data points found" << endl;
        return c;
    }
    g->SetLineColor( iColor );
    g->SetLineWidth( ( Width_t )iLineWidth );
    g->SetLineStyle( iLineStyle );
    
    // check what units the data is in
    string i_fileUnit = "PFLUX";
    if( iTextFile.find( "ENERGY" ) < iTextFile.size() )
    {
        i_fileUnit = "ENERGY";
    }
    else if( iTextFile.find( "CU" ) < iTextFile.size() )
    {
        i_fileUnit = "CU";
    }
    
    // convert energy to log10( energy )
    double i_Energy = 0.;
    double i_fluxSensitivity = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, i_Energy, i_fluxSensitivity );
        // convert flux
        if( bUnit == "ENERGY" && i_fileUnit != "ENERGY" )
        {
            i_fluxSensitivity *= i_Energy * 1.e12 * fConstant_Flux_To_Ergs;
        }
        i_Energy = log10( i_Energy );
        g->SetPoint( i, i_Energy, i_fluxSensitivity * iMultFactor );
    }
    // draw sensitivity graph
    g->Draw( "cp" );
    
    c->Update();
    
    return c;
}

/*

   list the available units for flux sensitivity calculations

*/
void VSensitivityCalculator::listUnits()
{
    cout << "Available units for sensitivity plot: " << endl;
    cout << "\t CU    = Crab units" << endl;
    cout << "\t PFLUX = particle flux [1/cm2/s]" << endl;
    cout << "\t ENERG = energy flux [erg/cm2/s]" << endl;
}


void VSensitivityCalculator::setEnergyRange_Log( double iN, double iX )
{
    fEnergy_min_Log = iN;
    fEnergy_max_Log = iX;
}


void VSensitivityCalculator::setEnergyRange_Lin( double iN, double iX )
{
    fEnergy_min_Log = log10( iN );
    fEnergy_max_Log = log10( iX );
}


void VSensitivityCalculator::setSignificanceParameter( double iSignificance, double iMinEvents,
        double iObservationTime, double iBackgroundRateRatio_min,
        double iAlpha )
{
    fSignificance_min = iSignificance;
    fEvents_min = iMinEvents;
    fObservationTime_h = iObservationTime;
    fMinBackgroundRateRatio_min = iBackgroundRateRatio_min;
    fAlpha = iAlpha;
}


void VSensitivityCalculator::setObservationTimeRange( double iObs_min, double iObs_max, int iObs_steps )
{
    fObservationTime_min = iObs_min;
    fObservationTime_max = iObs_max;
    fObservationTime_steps = iObs_steps;
}


void VSensitivityCalculator::plotObservationTimevsFluxFromTextFile( TCanvas* c, string iTextFile, int iLineColor, double iLineWidth, int iLineStyle )
{
    if( !c )
    {
        return;
    }
    
    TGraph* g = new TGraph( iTextFile.c_str() );
    g->SetLineColor( iLineColor );
    g->SetLineWidth( ( Width_t )iLineWidth );
    g->SetLineStyle( iLineStyle );
    
    cout << "reading " << iTextFile << " with " << g->GetN() << " data points" << endl;
    
    g->Draw( "c" );
}


TCanvas* VSensitivityCalculator::plotObservationTimevsFlux( unsigned int iD, TCanvas* c, int iLineColor, double iLineWidth, bool bGuidingLines )
{
    if( !checkDataSet( iD, "plotObservationTimevsFlux" ) )
    {
        return 0;
    }
    
    calculateObservationTimevsFlux( iD );
    
    bool bNewCanvas = ( c == 0 );
    if( c == 0 )
    {
        c = new TCanvas( "cF", "flux vs time", 510, 10, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLogx( 1 );
        c->SetLogy( 1 );
        c->SetLeftMargin( 0.12 );
        c->SetRightMargin( 0.07 );
        c->Draw();
    }
    setGraphPlottingStyle( fGraphObsvsTime[iD], iLineColor, iLineWidth );
    
    fGraphObsvsTime[iD]->SetTitle( "" );
    if( bNewCanvas )
    {
        fGraphObsvsTime[iD]->Draw( "ac" );
    }
    else
    {
        fGraphObsvsTime[iD]->Draw( "c" );
    }
    fGraphObsvsTime[iD]->GetHistogram()->SetYTitle( "observation time [h]" );
    fGraphObsvsTime[iD]->GetHistogram()->SetXTitle( "flux [Crab Units]" );
    fGraphObsvsTime[iD]->GetHistogram()->GetYaxis()->SetTitleOffset( 1.3 );
    if( fDebug )
    {
        fGraphObsvsTime[iD]->Print();
    }
    
    if( bNewCanvas && bGuidingLines )
    {
        for( unsigned int i = 0; i < fSourceStrength.size(); i++ )
        {
            plot_guidingLines( fSourceStrength[i], fGraphObsvsTime[iD], ( fSourceStrength[i] < 0.05 ) );
        }
    }
    
    return c;
}


void VSensitivityCalculator::list_sensitivity( unsigned int iD, ostream& terminal )
{
    if( !checkDataSet( iD, "plotObservationTimevsFlux" ) )
    {
        return ;
    }
    
    calculateObservationTimevsFlux( iD );
    
    terminal << " Flux               time           time" << endl;
    terminal << " [Crab Units]       [min]           [h]" << endl;
    terminal << " ========================================" << endl;
    
    for( unsigned int i = 0; i < fSourceStrength.size(); i++ )
    {
        terminal << " " << fSourceStrength[i] << "\t\t" << setw( 8 ) << setprecision( 3 );
        terminal << fGraphObsvsTime[iD]->Eval( fSourceStrength[i] ) * 60. << "\t" << setw( 8 ) << setprecision( 3 );
        terminal << fGraphObsvsTime[iD]->Eval( fSourceStrength[i] );
        terminal << endl;
    }
    
    terminal << "(requiring a significance of at least " << fSignificance_min << " sigma or " << fEvents_min;
    terminal << " events, and using Li & Ma formula " << fLiAndMaEqu << ")" << endl;
}

/*
 * calculate observation time vs flux
 *
 */
double VSensitivityCalculator::calculateObservationTimevsFlux( unsigned int iD )
{
    if( !checkDataSet( iD, "plotObservationTimevsFlux" ) )
    {
        return 0.;
    }
    
    fGraphObsvsTime[iD] = new TGraph( 100 );
    
    double s = 0.;
    double t = 0.;
    double x = 0.;
    double iG = fData[iD].fSignal;
    double iB = fData[iD].fBackground;
    double alpha = fData[iD].fAlpha;
    if( alpha > 0. )
    {
        iB /= alpha;
    }
    else
    {
        return 0.;
    }
    
    if( fDebug )
    {
        cout << "Rates for ID " << iD << ", " << iG << ", " << iB << ", " << alpha << endl;
    }
    
    int z = 0;
    int z_max = fGraphObsvsTime[iD]->GetN();
    for( int i = 0; i < z_max; i++ )
    {
        // take logarithmic steps in flux [log10 CU]
        x =   TMath::Log10( fSourceStrength_min ) + ( TMath::Log10( fSourceStrength_max )
                - TMath::Log10( fSourceStrength_min ) ) / ( double )fGraphObsvsTime[iD]->GetN() * ( double )i;
        // linear flux [CU]
        x =   TMath::Power( 10., x );
        
        // loop over possible observation lengths
        bool bSuccess = false;
        for( int j = 0; j < fObservationTime_steps; j++ )
        {
            // log10 hours
            t = TMath::Log10( fObservationTime_min ) + ( TMath::Log10( fObservationTime_max ) -
                    TMath::Log10( fObservationTime_min ) ) / ( double )fObservationTime_steps * ( double )j;
            // log10 hours to min
            t = TMath::Power( 10., t ) * 60.;
            
            s = VStatistics::calcSignificance( iG * t * x + iB * t * alpha, iB * t, alpha, fLiAndMaEqu );
            
            if( s > fSignificance_min && t * x * iG >= fEvents_min )
            {
                bSuccess = true;
                break;
            }
        }
        if( bSuccess )
        {
            fGraphObsvsTime[iD]->SetPoint( z, x, t / 60. );
            z++;
        }
    }
    
    return fGraphObsvsTime[iD]->Eval( fSignificance_min );
}


void VSensitivityCalculator::plot_guidingLines( double x, TGraph* g, bool iHours )
{
    if( !g )
    {
        return;
    }
    
    double i_y = g->Eval( x );
    TLine* iL_x = new TLine( x, 0., x, i_y );
    iL_x->SetLineStyle( 2 );
    iL_x->SetLineColor( 2 );
    iL_x->Draw();
    TLine* iL_y = new TLine( 0., i_y, x, i_y );
    iL_y->SetLineStyle( 4 );
    iL_y->SetLineColor( 2 );
    iL_y->Draw();
    
    char hname[600];
    if( iHours )
    {
        if( x >= 1.e-2 )
        {
            sprintf( hname, "%d%% Crab in %d h", ( int )( x * 100. ), ( int )( i_y + 0.5 ) );
        }
        else
        {
            sprintf( hname, "%.2f%% Crab in %d h", x * 100., ( int )( i_y + 0.5 ) );
        }
    }
    else
    {
        if( x > 1.e-2 )
        {
            sprintf( hname, "%d%% Crab in %d min", ( int )( x * 100. ), ( int )( i_y * 60. + 0.5 ) );
        }
        else
        {
            sprintf( hname, "%.2f%% Crab in %d h", x * 100., ( int )( i_y + 0.5 ) );
        }
    }
    TText* iT = new TText( x * 1.3, i_y, hname );
    iT->SetTextSize( iT->GetTextSize() * 0.7 );
    iT->SetTextAngle( 0. );
    iT->SetTextColor( 2 );
    iT->Draw();
}


void VSensitivityCalculator::setSourceStrengthVector_CU()
{
    fSourceStrength.clear();
    fSourceStrength.push_back( 1. );   // this is 100% Crab
    fSourceStrength.push_back( 0.30 );
    fSourceStrength.push_back( 0.10 );
    fSourceStrength.push_back( 0.05 );
    fSourceStrength.push_back( 0.03 );
    fSourceStrength.push_back( 0.01 );
}


bool VSensitivityCalculator::checkUnits( string bUnit )
{
    if( bUnit != "CU" && bUnit != "PFLUX" && bUnit != "ENERGY" )
    {
        return false;
    }
    
    return true;
}


vector< VDifferentialFluxData > VSensitivityCalculator::getDifferentialFluxVectorfromMC_ErrorMessage( string i_message )
{
    cout << "VSensitivityCalculator::getDifferentialFluxVectorfromMC error: " << i_message << endl;
    
    vector< VDifferentialFluxData > a;
    
    return a;
}

/*

    read integral/differential sensitivities and backround rates (total, proton, electron) from WP Phys File

    allows to cut for maximum energy bias (default: off)

    returns: 
    - sensitivity graph

    fills:
    - proton rate graph gProtonRate
    - electron rate graph gElectronRate


*/
TGraphAsymmErrors* VSensitivityCalculator::getSensitivityGraphFromWPPhysFile( string bUnit, double iEnergyMin_TeV_lin,
        double iEnergyMax_TeV_lin, double dE_log10 )
{
    TGraphAsymmErrors* g = 0;
    TH1F* h = 0;
    
    TFile f( fMCCTA_File.c_str() );
    if( f.IsZombie() )
    {
        cout << "VSensitivityCalculate::getSensitivityGraphFromWPPhysFile(): error: CTA-MC file not found" << endl;
        cout << fMCCTA_File << endl;
        return 0;
    }
    cout << "reading CTA-MC file (unit=" << bUnit << "): " << fMCCTA_File << endl;
    // sensitivities
    if( bUnit == "ENERGY" )
    {
        if( dE_log10 > 0. )
        {
            h = get_CTA_IRF_Histograms( "DiffSensE2Erg", fMCCTA_cameraoffset_deg );
            if( !h )
            {
                h = get_CTA_IRF_Histograms( "DiffSens", fMCCTA_cameraoffset_deg );
            }
        }
        else
        {
            h = get_CTA_IRF_Histograms( "IntSensE2Erg", fMCCTA_cameraoffset_deg );
            if( !h )
            {
                h = get_CTA_IRF_Histograms( "IntSens", fMCCTA_cameraoffset_deg );
            }
            cout << "getting integral sensitivity histograms" << endl;
        }
    }
    else
    {
        if( dE_log10 > 0. )
        {
            h = get_CTA_IRF_Histograms( "DiffSensCU", fMCCTA_cameraoffset_deg );
        }
        else
        {
            h = get_CTA_IRF_Histograms( "IntSensCU", fMCCTA_cameraoffset_deg );
            cout << "getting integral sensitivity histograms" << endl;
        }
    }
    if( h )
    {
        g = new TGraphAsymmErrors( 1 );
        get_Graph_from_Histogram( h, g, false, false, 1.e3, log10( iEnergyMin_TeV_lin ), log10( iEnergyMax_TeV_lin ) );
        setGraphPlottingStyle( g, 1, 1, 20, 2 );
        smoothenSensitivityGraph( g, get_CTA_IRF_Histograms( "BGRate", fMCCTA_cameraoffset_deg ) );
    }
    
    // background rates
    h = 0;
    h = get_CTA_IRF_Histograms( "BGRate", fMCCTA_cameraoffset_deg );
    if( h )
    {
        gBGRate     = new TGraphAsymmErrors( 1 );
        get_Graph_from_Histogram( h, gBGRate, true, false, -1., log10( iEnergyMin_TeV_lin ), log10( iEnergyMax_TeV_lin ) );
        setGraphPlottingStyle( gBGRate, 1, 2, 20, 2 );
    }
    // proton rates
    h = 0;
    h = get_CTA_IRF_Histograms( "ProtRate", fMCCTA_cameraoffset_deg );
    if( !h )
    {
        h = get_CTA_IRF_Histograms( "hProtRate", fMCCTA_cameraoffset_deg );
    }
    if( h )
    {
        gProtonRate = new TGraphAsymmErrors( 1 );
        get_Graph_from_Histogram( h, gProtonRate, true, 0., log10( iEnergyMin_TeV_lin ), log10( iEnergyMax_TeV_lin ) );
        setGraphPlottingStyle( gProtonRate, 1, 2, 21, 2 );
    }
    // electron rates
    h = 0;
    h = get_CTA_IRF_Histograms( "ElecRate", fMCCTA_cameraoffset_deg );
    if( !h )
    {
        h = get_CTA_IRF_Histograms( "hElecRate", fMCCTA_cameraoffset_deg );
    }
    if( h )
    {
        gElectronRate = new TGraphAsymmErrors( 1 );
        get_Graph_from_Histogram( h, gElectronRate, true, 0., log10( iEnergyMin_TeV_lin ), log10( iEnergyMax_TeV_lin ) );
        setGraphPlottingStyle( gElectronRate, 1, 2, 22, 2 );
    }

    // remove points with large energy bias
    if( fMaxEnergyBias < 1.e4 && fMaxEnergyBias > 0. )
    {
        h = get_CTA_IRF_Histograms( "Ebias", fMCCTA_cameraoffset_deg );
        if( !h )
        {
            cout << "Warning: energy bias histogram not found!";
            cout << " No cut on energy bias applied";
            cout << " (requested was " << fMaxEnergyBias << ")" << endl;
        }
        else
        {
            applyEnergyBiasCut( h, g, fMaxEnergyBias );
        }
    }
    
    f.Close();
    
    return g;
}

/*
 * remove all points of a graph for which a set energy bias critera is violated
 *
 * e.g. ignore sensitivities for which the bias is larger than 10%
 *
 */
void VSensitivityCalculator::applyEnergyBiasCut( TH1F* h, TGraphAsymmErrors* g, double iMaxEBias )
{
    if( !h || !g )
    {
        return;
    }

    int nRemoved = 0;

    double x = 0.;
    double y = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        if( h->FindBin( x ) > 0. 
                && h->GetBinContent( h->FindBin( x ) ) > iMaxEBias )
        {
            g->RemovePoint( i );
            nRemoved++;
        }
    }
    if( nRemoved > 0 )
    {
        cout << "removed " << nRemoved << " point(s) due to maximum energy bias constraint";
        cout << " (<" << iMaxEBias << ")" << endl;
    }

}

/*

   calculate different or integral flux from CTA-MC files

   (WP Phys style)

*/
vector< VDifferentialFluxData > VSensitivityCalculator::getDifferentialFluxVectorfromCTA_MC( double& iNorm )
{
    vector< VDifferentialFluxData > a;
    
    TFile f( fMCCTA_File.c_str() );
    if( f.IsZombie() )
    {
        cout << "VSensitivityCalculate::getDifferentialFluxVectorfromCTA_MC: error: CTA-MC file not found" << endl;
        exit( -1 );
    }
    cout << "reading CTA-MC file: " << fMCCTA_File << endl;
    TH1F* hEff = ( TH1F* )f.Get( "EffectiveArea" );
    if( !hEff )
    {
        cout << "VSensitivityCalculate::getDifferentialFluxVectorfromCTA_MC: error: no effective area histogram found " << endl;
        return a;
    }
    vector< double > iEnergy;
    vector< double > iEff;
    for( int i = 1; i <= hEff->GetNbinsX(); i++ )
    {
        iEnergy.push_back( hEff->GetBinCenter( i ) );
        iEff.push_back( hEff->GetBinContent( i ) );
    }
    
    TH1F* hBck = ( TH1F* )f.Get( "BGRate" );
    if( !hBck )
    {
        cout << "VSensitivityCalculate::getDifferentialFluxVectorfromCTA_MC: error: no background histogram found " << endl;
        return a;
    }
    VMonteCarloRateCalculator iMCR;
    // setting up differential flux vector
    for( int i = 1; i <= hEff->GetNbinsX(); i++ )
    {
        if( hEff->GetBinContent( i ) > 0. && hBck->GetBinContent( i ) > 0. )
        {
            VDifferentialFluxData i_flux;
            i_flux.Energy_lowEdge = hEff->GetBinLowEdge( i );
            i_flux.Energy_upEdge = hEff->GetBinLowEdge( i ) + hEff->GetBinWidth( i );
            i_flux.Energy_lowEdge_bin = i;
            i_flux.Energy_upEdge_bin = i;
            // energies are on linear scale
            i_flux.Energy_lowEdge = TMath::Power( 10., i_flux.Energy_lowEdge );
            i_flux.Energy_upEdge  = TMath::Power( 10., i_flux.Energy_upEdge );
            i_flux.dE = i_flux.Energy_upEdge - i_flux.Energy_lowEdge;
            i_flux.Energy = ( i_flux.Energy_lowEdge + i_flux.Energy_upEdge ) / 2.;
            i_flux.EnergyWeightedMean = i_flux.Energy;
            i_flux.ExposureTime = fObservationTime_h * 60. * 60.;
            
            // number of background events
            i_flux.NOff = hBck->GetBinContent( i ) * i_flux.ExposureTime;
            i_flux.NOff_error = sqrt( i_flux.NOff );
            // number of signal events
            i_flux.NOn = iMCR.getMonteCarloRate( iEnergy, iEff, fEnergySpectrumfromLiterature, fEnergySpectrumfromLiterature_ID,
                                                 i_flux.Energy_lowEdge_bin, i_flux.Energy_upEdge_bin );
            i_flux.NOn *= i_flux.ExposureTime / 60.;
            i_flux.NOn += i_flux.NOff;
            i_flux.NOn_error  = sqrt( i_flux.NOn );
            // multiply off events by alpha (number of background regions)
            i_flux.NOff /= fAlpha;
            i_flux.NOff_error /= fAlpha;
            
            a.push_back( i_flux );
        }
    }
    iNorm = fAlpha;
    
    // this is not necessary when count rates are read from CTA WP Phys file
    fRequireCutsToBeOptimized = false;
    
    f.Close();
    
    return a;
}

/*

    calculate differential or integral flux from gamma and proton effective areas

*/
vector< VDifferentialFluxData > VSensitivityCalculator::getDifferentialFluxVectorfromMC( double dE_Log10, double& iNorm )
{
    vector< VDifferentialFluxData > a;
    iNorm = 0.;
    // iterator with MC data
    char hname[800];
    map< unsigned int, VSensitivityCalculatorDataResponseFunctions* >::iterator i_MCData_iterator;
    
    ///////////////////////////////////////////////////////////////////
    // check if data is complete (need at least gamma-ray data)
    if( fMC_Data.find( 1 ) == fMC_Data.end() )
    {
        return getDifferentialFluxVectorfromMC_ErrorMessage( "no gamma-ray MC data given" );
    }
    if( fMC_Data.size() < 2 )   // (need at least one background data set)
    {
        return getDifferentialFluxVectorfromMC_ErrorMessage( "MC Data vector not large enough (should be >=2)" );
    }
    
    ///////////////////////////////////////////////////////////////////
    // differential flux vector (for gammas and background particles)
    vector< VDifferentialFluxData > v_flux;
    
    ///////////////////////////////////////////////////////////////////
    // get Crab spectrum from literature ( index [1] is gamma-ray)
    VEnergySpectrumfromLiterature i_Crab( fMC_Data[1]->fSpectralParameterFile );
    cout << "\t reading Crab Nebula spectrum with ID " <<  fMC_Data[1]->fSpectralParameterID << endl;
    if( i_Crab.isZombie() )
    {
        return a;
    }
    std::cout << "VSensitivityCalculator::getDifferentialFluxVectorfromMC Crab " << std::endl;
    i_Crab.listValues( fMC_Data[1]->fSpectralParameterID );
    
    ///////////////////////////////////////////////////////////////////
    // get effective areas for gamma-rays and background particles
    
    if( fDebug )
    {
        cout << "get MC effective areas" << endl;
    }
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); i_MCData_iterator++ )
    {
        if( !getMonteCarlo_EffectiveArea( ( *i_MCData_iterator ).second ) )
        {
            return a;
        }
    }
    
    ///////////////////////////////////////////////////////////////////
    // set up energy binning for differential flux vector
    
    ///////////////////////////////////////////////////////////////////
    // binning in gamma and proton effective areas must be the same (number of bins, range, and bin width)
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); i_MCData_iterator++ )
    {
        if( ( *i_MCData_iterator ).first == 1 )
        {
            continue;
        }
        
        if( fMC_Data[1]->effArea_Ebins != ( *i_MCData_iterator ).second->effArea_Ebins )
        {
            sprintf( hname, "diffent number of bins in gamma and background effective area (%d, %d)",
                     fMC_Data[1]->effArea_Ebins, ( *i_MCData_iterator ).second->effArea_Ebins );
            return getDifferentialFluxVectorfromMC_ErrorMessage( hname );
        }
        if( TMath::Abs( fMC_Data[1]->effArea_Emin - ( *i_MCData_iterator ).second->effArea_Emin ) > 0.05
                || TMath::Abs( fMC_Data[1]->effArea_Emax - ( *i_MCData_iterator ).second->effArea_Emax ) > 0.05 )
        {
            return getDifferentialFluxVectorfromMC_ErrorMessage( "different energy axis definition in gamma and background effective areas" );
        }
        if( ( *i_MCData_iterator ).second->energy.size() == 0 )
        {
            sprintf( hname, "warning: effective area (%s) vector with length 0",
                     ( *i_MCData_iterator ).second->fName.c_str() );
            return getDifferentialFluxVectorfromMC_ErrorMessage( hname );
        }
    }
    
    // integral sensitivity: take existing energy binning (from gammas)
    double iBinSize = ( fMC_Data[1]->effArea_Emax - fMC_Data[1]->effArea_Emin ) / fMC_Data[1]->effArea_Ebins;
    // offset between gamma and background effective areas
    // (offset in effective area vector, determine first filled value)
    map< unsigned int, int > iEnergyScaleOffset;
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); i_MCData_iterator++ )
    {
        if( ( *i_MCData_iterator ).first == 1 )
        {
            continue;
        }
        if( fMC_Data[1]->energy.size() > 0 && ( *i_MCData_iterator ).second->energy.size() > 0 )
        {
            iEnergyScaleOffset[( *i_MCData_iterator ).first] = ( int )( ( fMC_Data[1]->energy[0] - ( *i_MCData_iterator ).second->energy[0] ) / iBinSize );
            if( fDebug )
            {
                cout << "Energy scale offset (particle " << ( *i_MCData_iterator ).first << "): " << endl;
                cout << "   first bin gamma: " << fMC_Data[1]->energy[0] << endl;
                cout << "   first bin particle " << ( *i_MCData_iterator ).first << ": " << ( *i_MCData_iterator ).second->energy[0] << endl;
                cout << "   bin size: " << iBinSize << endl;
                cout << "   offset (bin bins): " << ( int )( ( fMC_Data[1]->energy[0] - ( *i_MCData_iterator ).second->energy[0] ) / iBinSize ) << endl;
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////
    // GAMMA RAYS (signal)
    ///////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////
    // INTEGRAL SENSITIVITY
    ////////////////////////////////////////////////////////////////
    if( dE_Log10 < 0. )
    {
        for( unsigned int i = 0; i < fMC_Data[1]->energy.size(); i++ )
        {
            VDifferentialFluxData i_flux;
            // energy bins in sensitivity curve
            i_flux.Energy_lowEdge = fMC_Data[1]->energy[i] - iBinSize / 2.;
            i_flux.Energy_upEdge  = fMC_Data[1]->energy[i] + iBinSize / 2.;
            i_flux.Energy_lowEdge_bin = i;
            i_flux.Energy_upEdge_bin = i;
            // energies are on linear scale
            i_flux.Energy_lowEdge = TMath::Power( 10., i_flux.Energy_lowEdge );
            i_flux.Energy_upEdge  = TMath::Power( 10., i_flux.Energy_upEdge );
            
            i_flux.dE = i_flux.Energy_upEdge - i_flux.Energy_lowEdge;
            // simplified (maybe should somehow be weighted by the shape of the energy spectrum)
            i_flux.Energy = ( i_flux.Energy_lowEdge + i_flux.Energy_upEdge ) / 2.;
            i_flux.EnergyWeightedMean = i_flux.Energy;
            // convert from [h] to [s]
            i_flux.ExposureTime = fObservationTime_h * 60. * 60.;
            
            v_flux.push_back( i_flux );
            if( fDebug )
            {
                cout << "ENERGY: " << v_flux.size() << "\t" << i_flux.Energy_lowEdge << " - " << i_flux.Energy_upEdge;
                cout << "\t" << i_flux.Energy << "\t" << i_flux.EnergyWeightedMean;
                cout << "\t" << i_flux.dE << "\t" << i_flux.ExposureTime << "\t" << fObservationTime_h << "  Off: ";
                map< unsigned int, int >::iterator iEnergyScaleOffset_iter;
                for( iEnergyScaleOffset_iter = iEnergyScaleOffset.begin();
                        iEnergyScaleOffset_iter != iEnergyScaleOffset.end(); ++iEnergyScaleOffset_iter )
                {
                    cout << ( *iEnergyScaleOffset_iter ).second << "\t";
                }
                cout << endl;
            }
        }
    }
    ////////////////////////////////////////////////////////////////
    // DIFFERENTIAL SENSITIVITY
    ////////////////////////////////////////////////////////////////
    else
    {
        // dE_Log10 must be a multiple of effective area energy bin size
        if( dE_Log10 < iBinSize )
        {
            return getDifferentialFluxVectorfromMC_ErrorMessage( "error: dE_log10 smaller than bin size of effective areas" );
        }
        if( TMath::Abs( dE_Log10 / iBinSize - ( int )( dE_Log10 / iBinSize ) ) > 1.e-2 )
        {
            return getDifferentialFluxVectorfromMC_ErrorMessage( "error: dE_Log10 must be a multiple of effective area energy bin size" );
        }
        
        // minimum and maximum energy in differential sensitivity curve defined by extrema of true (MC!) energy
        // (despite there might be reconstructed energies with lower/higher values)
        double iBinEnergyMin = fMC_Data[1]->effArea_Emin;
        if( fDebug )
        {
            cout << "VSensitivityCalculator::getDifferentialFluxVectorfromMC: minimum energy: " << iBinEnergyMin << endl;
        }
        // use first bin low edge of gamma-ray effective area
        if( fMC_Data[1]->energy_lowEdge.size() > 0 )
        {
            iBinEnergyMin = fMC_Data[1]->energy_lowEdge[0];
        }
        // now fill vector with energies and energy bins
        while( iBinEnergyMin < fMC_Data[1]->effArea_Emax )
        {
            VDifferentialFluxData i_flux;
            
            bool i_eff_found = false;
            // lower energy bin
            for( unsigned int i = 0; i < fMC_Data[1]->energy.size(); i++ )
            {
                if( TMath::Abs( fMC_Data[1]->energy[i] - iBinSize / 2. - iBinEnergyMin ) < 1.e-2 )
                {
                    // temporary assignment: Energy_lowEdge should be on lin scale, here it is on log scale
                    i_flux.Energy_lowEdge = fMC_Data[1]->energy[i] - iBinSize / 2.;
                    i_flux.Energy_lowEdge_bin = i;
                    i_eff_found = true;
                    break;
                }
            }
            if( i_eff_found )
            {
                // upper energy bin
                if( i_flux.Energy_lowEdge_bin + int( dE_Log10 / iBinSize ) < fMC_Data[1]->energy.size() )
                {
                    i_flux.Energy_upEdge     = i_flux.Energy_lowEdge + dE_Log10;
                    i_flux.Energy_upEdge_bin = i_flux.Energy_lowEdge_bin + int( dE_Log10 / iBinSize ) - 1;
                    // energies are on linear scale
                    i_flux.Energy_lowEdge = TMath::Power( 10., i_flux.Energy_lowEdge );
                    i_flux.Energy_upEdge  = TMath::Power( 10., i_flux.Energy_upEdge );
                    
                    i_flux.dE = i_flux.Energy_upEdge - i_flux.Energy_lowEdge;
                    // simplified
                    i_flux.Energy = ( i_flux.Energy_lowEdge + i_flux.Energy_upEdge ) / 2.;
                    i_flux.EnergyWeightedMean = i_flux.Energy;
                    // convert from [h] to [s]
                    i_flux.ExposureTime = fObservationTime_h * 60. * 60.;
                    
                    v_flux.push_back( i_flux );
                    if( fDebug )
                    {
                        cout << "ENERGY (diff): " << v_flux.size() << "\t";
                        cout << i_flux.Energy_lowEdge_bin << "-" << i_flux.Energy_upEdge_bin << "  ";
                        cout << i_flux.Energy_lowEdge << " - " << i_flux.Energy_upEdge << " TeV, [";
                        cout << log10( i_flux.Energy_lowEdge ) << ", " << log10( i_flux.Energy_upEdge ) << "] ";
                        cout << "\t" << i_flux.Energy << "\t" << i_flux.EnergyWeightedMean;
                        cout << "\t dE: " << i_flux.dE << "\tObsTime: " << i_flux.ExposureTime << " [s]";
                        cout << ", Offset: ";
                        map< unsigned int, int >::iterator iEnergyScaleOffset_iter;
                        for( iEnergyScaleOffset_iter = iEnergyScaleOffset.begin(); iEnergyScaleOffset_iter != iEnergyScaleOffset.end();
                                ++iEnergyScaleOffset_iter )
                        {
                            cout << ( *iEnergyScaleOffset_iter ).second << "\t";
                        }
                        cout << endl;
                    }
                }
            }
            iBinEnergyMin += dE_Log10;
        }
    }
    
    ///////////////////////////////////////////////////////////////////
    // get gamma rate [1/min] for a certain ze, az, noise, wobble offset
    ///////////////////////////////////////////////////////////////////
    // for integral sensitivity set energy bins according to effective areas
    if( fDebug )
    {
        cout << "Calculate rates for gamma rays" << endl;
    }
    VMonteCarloRateCalculator iMCR;
    for( unsigned int i = 0; i < v_flux.size(); i++ )
    {
        v_flux[i].NOn       = getMonteCarloRateFromWeightedRateHistogram( v_flux[i].Energy_lowEdge, v_flux[i].Energy_upEdge,
                              false, fMC_Data[1]->hWeightedRate );
        v_flux[i].NOn_error = getMonteCarloRateFromWeightedRateHistogram( v_flux[i].Energy_lowEdge, v_flux[i].Energy_upEdge,
                              true, fMC_Data[1]->hWeightedRate );
    }
    
    
    ///////////////////////////////////////////////////////////////////
    // COSMIC RAYS (background)
    ///////////////////////////////////////////////////////////////////
    
    map< unsigned int, vector< double > > v_flux_NOff;
    map< unsigned int, vector< double > > v_flux_NOff_error;
    
    // loop over all background files
    // (i.e. background particle types, usually electrons=2 and protons=14)
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); i_MCData_iterator++ )
    {
        if( ( *i_MCData_iterator ).first == 1 )
        {
            continue;    // ignore gamma rays
        }
        //std::cout<<"i_MCData_iterator "<<i_MCData_iterator<<std::endl;
        ///////////////////////////////////////////////////////////////////
        // get CR spectrum from literature
        VEnergySpectrumfromLiterature i_CR( ( *i_MCData_iterator ).second->fSpectralParameterFile );
        cout << "\t reading CR spectrum with ID" << ( *i_MCData_iterator ).second->fSpectralParameterID << endl;
        if( i_CR.isZombie() )
        {
            return a;
        }
        std::cout << "VSensitivityCalculator::getDifferentialFluxVectorfromMC CR " << std::endl;
        i_CR.listValues( ( *i_MCData_iterator ).second->fSpectralParameterID );
        
        // get CR rate for a certain ze, az, noise, wobble offset
        for( unsigned int i = 0; i < v_flux.size(); i++ )
        {
            v_flux_NOff[( *i_MCData_iterator ).first].push_back( 0. );
            v_flux_NOff_error[( *i_MCData_iterator ).first].push_back( 0. );
        }
        
        VMonteCarloRateCalculator iMCR;
        // loop over all energy bins
        for( unsigned int i = 0; i < v_flux.size(); i++ )
        {
            v_flux_NOff[( *i_MCData_iterator ).first][i]       =  getMonteCarloRateFromWeightedRateHistogram(
                        v_flux[i].Energy_lowEdge,
                        v_flux[i].Energy_upEdge, false,
                        ( *i_MCData_iterator ).second->hWeightedRate );
            v_flux_NOff_error[( *i_MCData_iterator ).first][i] =  getMonteCarloRateFromWeightedRateHistogram(
                        v_flux[i].Energy_lowEdge,
                        v_flux[i].Energy_upEdge, true,
                        ( *i_MCData_iterator ).second->hWeightedRate );
                        
            //////////////////////////////////////////////////////////////////////////////////////////
            // take care of space angle and theta2 cut normalisation
            //////////////////////////////////////////////////////////////////////////////////////////
            
            // CR spectrum is given as dN/dt/dOmega (Omega = space angle)
            // multiply CR rate by space angle used in CORSIKA (scattering angle)
            v_flux_NOff[( *i_MCData_iterator ).first][i]       *= ( *i_MCData_iterator ).second->SolidAngle_MCScatterAngle;
            v_flux_NOff_error[( *i_MCData_iterator ).first][i] *= ( *i_MCData_iterator ).second->SolidAngle_MCScatterAngle;
            // scale direction cut CR solid angle to gamma-ray solid angle
            // (CR theta2 might be larger than gamma-ray, simply to gain statistics under the assumption of a flat angular acceptance)
            double iSolidAngle_Gamma =  fMC_Data[1]->getSolidAngle_DirectionCut( v_flux[i].Energy );
            double iSolidAngle_Bck   = ( *i_MCData_iterator ).second->getSolidAngle_DirectionCut( v_flux[i].Energy );
            if( iSolidAngle_Bck > 0. )
            {
                v_flux_NOff[( *i_MCData_iterator ).first][i]       *= iSolidAngle_Gamma / iSolidAngle_Bck;
                v_flux_NOff_error[( *i_MCData_iterator ).first][i] *= iSolidAngle_Gamma / iSolidAngle_Bck;
                if( fDebug )
                {
                    cout << "SOLID ANGLE GAMMA AND BCK: " << i;
                    cout << "\t" << iSolidAngle_Gamma << "\t" << iSolidAngle_Bck << "\t" << iSolidAngle_Gamma / iSolidAngle_Bck;
                    cout << "\t" << v_flux[i].Energy;
                    cout << endl;
                }
            }
            // missing helium (etc) simulations
            v_flux_NOff[( *i_MCData_iterator ).first][i]       *= ( 1. + fMC_BackgroundMissingParticleFraction );
            v_flux_NOff_error[( *i_MCData_iterator ).first][i] *= ( 1. + fMC_BackgroundMissingParticleFraction );
        }
    }
    
    //////////////////////////////////
    // now sum up all background
    //////////////////////////////////
    
    // add background rate to signal rate (NOn is gamma-ray + background rate in signal region)
    map< unsigned int, vector< double > >::iterator v_flux_NOff_iter;
    for( unsigned int i = 0; i < v_flux.size(); i++ )
    {
        v_flux[i].NOff = 0.;
        for( v_flux_NOff_iter = v_flux_NOff.begin(); v_flux_NOff_iter != v_flux_NOff.end(); ++v_flux_NOff_iter )
        {
            if( i < ( *v_flux_NOff_iter ).second.size() )
            {
                // require at least some proton events
                if( i < v_flux_NOff[14].size() && v_flux_NOff[14][i] > 0. )
                {
                    v_flux[i].NOn  += ( *v_flux_NOff_iter ).second[i];
                    v_flux[i].NOff += ( *v_flux_NOff_iter ).second[i];
                }
            }
        }
    }
    // error calculation
    for( unsigned int i = 0; i < v_flux.size(); i++ )
    {
        v_flux[i].NOn_error = v_flux[i].NOn_error * v_flux[i].NOn_error;
        v_flux[i].NOff_error = 0.;
        for( v_flux_NOff_iter = v_flux_NOff_error.begin(); v_flux_NOff_iter != v_flux_NOff_error.end(); ++v_flux_NOff_iter )
        {
            if( i < ( *v_flux_NOff_iter ).second.size() )
            {
                v_flux[i].NOn_error  += ( *v_flux_NOff_iter ).second[i] * ( *v_flux_NOff_iter ).second[i];
            }
            if( i < ( *v_flux_NOff_iter ).second.size() )
            {
                v_flux[i].NOff_error += ( *v_flux_NOff_iter ).second[i] * ( *v_flux_NOff_iter ).second[i];
            }
        }
        v_flux[i].NOn_error  =  sqrt( v_flux[i].NOn_error );
        v_flux[i].NOff_error =  sqrt( v_flux[i].NOff_error );
    }
    fillBackgroundParticleNumbers( v_flux, v_flux_NOff, v_flux_NOff_error );
    if( fPlotDebugName.size() > 0 )
    {
        plotDebugPlotsBackgroundParticleNumbers( v_flux, v_flux_NOff, v_flux_NOff_error );
    }
    
    // calculate mean alpha
    
    // scale number of off events by alpha
    double alpha = 0.;
    double zz = 0;
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); i_MCData_iterator++ )
    {
        if( ( *i_MCData_iterator ).first == 1 )
        {
            continue;
        }
        
        alpha += ( *i_MCData_iterator ).second->alpha;
        zz++;
    }
    if( zz > 0. )
    {
        alpha /= zz;
    }
    else
    {
        alpha  = 1.;
    }
    
    for( unsigned int i = 0; i < v_flux.size(); i++ )
    {
        if( alpha > 0. )
        {
            v_flux[i].NOff /= alpha;
            v_flux[i].NOff_error /= alpha;
        }
    }
    
    if( fPlotDebugName.size() > 0 )
    {
        plotEffectiveArea();
    }
    
    ////////////////////////////////////////////////////
    // calculate number of on/off events from rates
    ////////////////////////////////////////////////////
    double iTotG = 0.;
    double iTotB = 0.;
    if( v_flux.size() > 0 )
    {
        for( unsigned int i = v_flux.size(); i > 0; i-- )
        {
            unsigned int z = i - 1;
            v_flux[z].NOn  *= fObservationTime_h * 60.;               // [min]
            v_flux[z].NOff *= fObservationTime_h * 60.;
            v_flux[z].NOn_error  *= fObservationTime_h * 60.;
            v_flux[z].NOff_error *= fObservationTime_h * 60.;
            iTotG += ( v_flux[z].NOn - v_flux[z].NOff * alpha ) / ( fObservationTime_h * 60. );
            iTotB += v_flux[z].NOff * alpha / ( fObservationTime_h * 60. );
            if( fDebug )
            {
                cout << "NUMBER OF MC ON/OFF EVENTS: " << z << "\t";
                cout << v_flux[z].Energy_lowEdge << " - " << v_flux[z].Energy_upEdge << " TeV,\t";
                cout << "Non-Noff: " << v_flux[z].NOn - v_flux[z].NOff* alpha  << "\t";
                cout << "NON: " << v_flux[z].NOn << " (" << v_flux[z].NOn_error << "), ";
                cout << "NOFF: " << v_flux[z].NOff << "(" << v_flux[z].NOff_error << "), "  << v_flux[z].NOff* alpha << "\t" << alpha << "\t";
                cout << "RateG: " << iTotG << " [1/min], " << iTotB << " [1/min]";
                cout << endl;
            }
        }
    }
    
    ////////////////////////////////////////////////////
    iNorm = alpha;  // (norm is a return value)
    return v_flux;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////


/*
 * read an anasum file and get a differential flux vector
 *
*/
vector< VDifferentialFluxData > VSensitivityCalculator::getDifferentFluxVectorfromData( string iAnasumCrabFile, double dE_Log10, double& iNorm )
{
    iNorm = 0;
    
    // read energy spectrum from file
    VEnergySpectrum e( iAnasumCrabFile );
    if( e.isZombie() )
    {
        vector< VDifferentialFluxData > a;
        return a;
    }
    
    if( dE_Log10 < 0. )
    {
        e.setEnergyBinning( 0.10 );
    }
    else
    {
        e.setEnergyBinning( dE_Log10 );
    }
    e.setEnergyRangeLinear( 0., 1.e10 );
    e.setSignificanceParameters( -50., -1. );
    // energy threshold: 10% in energy bias
    e.setEnergyThresholdDefinition( 1, 0.1, 0.2 );
    e.combineRuns();
    
    iNorm = e.getTotalNormalisationFactor();
    cout << endl << "Normalisation factor: " << iNorm << endl;
    
    cout << endl << "Differential fluxes and significances per differential bin:" << endl;
    e.printDifferentialFluxes();
    cout << endl;
    
    return e.getDifferentialFlux();
}


/*
   set all Monte Carlo parameters for gammas or protons

   iParticleID                 corsika particle ID: 1  = gamma, 14 = proton
   iSpectralParameterFile      read source spectrum from this file (VEnergySpectrumfromLiterature file)
   iSpectralParameterID        spectrum ID to read from iSpectralParameterFile
   iGammaEffectiveAreaFile     file with effective areas
   ze, az, woff, noise, index  parameters for effective area search (to be read from iGammaEffectiveAreaFile)
   theta2_MCSatterAngle        background scattering angle [deg2]
*/
void VSensitivityCalculator::setMonteCarloParameters( unsigned int iParticleID,
        string iSpectralParameterFile, unsigned int iSpectralParameterID,
        string iGammaEffectiveAreaFile,
        double ze, int az, double woff, int noise, double index,
        double iEnergy_min_log, double iEnergy_max_log, string bUnit )
{
    if( fMC_Data.find( iParticleID ) != fMC_Data.end() )
    {
        cout << "VSensitivityCalculator::setMCParameters: particle with ID " << iParticleID << " already in set of MC parameters" << endl;
        cout << "\t ignoring input" << endl;
        return;
    }
    if( iGammaEffectiveAreaFile.size() == 0 )
    {
        cout << "VSensitivityCalculator::setMCParameters: effective area not set for particle with ID " << iParticleID << endl;
        cout << "\t ignoring input" << endl;
        return;
    }
    
    // create a new MC data object
    VSensitivityCalculatorDataResponseFunctions* f = new VSensitivityCalculatorDataResponseFunctions();;
    fMC_Data[iParticleID] = f;
    
    // fill data
    f->fParticleID = iParticleID;
    f->fSpectralParameterFile = gSystem->ExpandPathName( iSpectralParameterFile.c_str() );
    f->fSpectralParameterID = iSpectralParameterID;
    f->fEffectiveAreaFile = gSystem->ExpandPathName( iGammaEffectiveAreaFile.c_str() );
    f->ze = ze;
    f->az = az;
    f->woff = woff;
    f->noise = noise;
    f->index = index;
    f->alpha = fAlpha;
    f->energy_min_log = iEnergy_min_log;
    f->energy_max_log = iEnergy_max_log;
    
    if( iParticleID == 1 )
    {
        f->fName = "gamma";
    }
    else if( iParticleID == 2 )
    {
        f->fName = "electron";
    }
    else if( iParticleID == 14 )
    {
        f->fName = "proton";
    }
    else if( iParticleID == 402 )
    {
        f->fName = "helium";
    }
    else
    {
        cout << "VSensitivityCalculator::setMCParameters: unknown particle ID: " << iParticleID  << endl;
        return;
    }
    f->initializeHistograms( bUnit );
}


/*!

   return rate (in [1/min]) from weighted rate histograms calculated in effective area code

*/
double VSensitivityCalculator::getMonteCarloRateFromWeightedRateHistogram( double iE_low, double iE_up,
        bool iRateError,
        TH1D* iWeightedRateHistogram )
{
    if( !iWeightedRateHistogram )
    {
        return 0.;
    }
    iE_low = log10( iE_low ) + 0.01;
    iE_up  = log10( iE_up ) - 0.01;
    
    // find lower bin (assume identical binning in weightedratehistogram and energy vector
    int iHis_minBin = iWeightedRateHistogram->GetXaxis()->FindBin( iE_low );
    int iHis_maxBin = iWeightedRateHistogram->GetXaxis()->FindBin( iE_up );
    
    double iR = 0.;
    for( int i = iHis_minBin; i <= iHis_maxBin; i++ )
    {
        if( !iRateError )
        {
            iR += iWeightedRateHistogram->GetBinContent( i );
        }
        else
        {
            iR += iWeightedRateHistogram->GetBinError( i ) * iWeightedRateHistogram->GetBinError( i );
        }
    }
    if( iRateError && iR > 0. )
    {
        return sqrt( iR );
    }
    
    return iR;
}


/*!
    fill effective area vector into VSensitivityCalculatorDataResponseFunctions
*/
bool VSensitivityCalculator::getMonteCarlo_EffectiveArea( VSensitivityCalculatorDataResponseFunctions* iMCPara )
{
    //////////////////////////////////////////////////////////////////////////////////////
    // read effective areas from file
    //////////////////////////////////////////////////////////////////////////////////////
    TFile fEff( iMCPara->fEffectiveAreaFile.c_str() );
    if( fEff.IsZombie() )
    {
        cout << "VSensitivityCalculator::getMonteCarlo_EffectiveArea: error, cannot find effective area file ";
        cout << iMCPara->fEffectiveAreaFile.c_str() << endl;
        exit( EXIT_FAILURE );
    }
    TTree* t = ( TTree* )fEff.Get( "fEffArea" );
    if( !t )
    {
        cout << "VSensitivityCalculator::getMonteCarlo_EffectiveArea: error, cannot find effective area tree in ";
        cout << iMCPara->fEffectiveAreaFile.c_str() << endl;
        return false;
    }
    CEffArea* c = new CEffArea( t );
    if( !c )
    {
        return false;
    }
    
    cout << endl;
    cout << "=================================================================================" << endl;
    cout << "reading effective areas for " << iMCPara->fName << " from " << iMCPara->fEffectiveAreaFile.c_str() << endl;
    cout << "\t total number of effective areas in this data file: " << c->fChain->GetEntries();
    cout << " (reading effective areas vs true energy)";
    cout << endl;
    
    iMCPara->energy.clear();
    iMCPara->effArea.clear();
    
    // loop over all effective areas in effective area tree and read the appropriate one
    bool bFound = false;
    for( unsigned int i = 0; i < c->fChain->GetEntries(); i++ )
    {
        c->GetEntry( i );
        
        // do selection of effective areas only if there are several effective areas
        if( c->fChain->GetEntries() > 1 )
        {
            if( TMath::Abs( c->index - iMCPara->index ) > 1.e-2 )
            {
                continue;
            }
            
            if( c->noise != iMCPara->noise )
            {
                continue;
            }
            
            if( c->az != iMCPara->az )
            {
                continue;
            }
            
            // ignore everything else for non-gammas (wobble offsets and zenith angles)
            if( iMCPara->fParticleID == 1 )
            {
                if( TMath::Abs( c->ze - iMCPara->ze ) > 2. )
                {
                    continue;
                }
                if( TMath::Abs( c->Woff - iMCPara->woff ) > 0.05 )
                {
                    continue;
                }
            }
        }
        
        cout << "\t found effective area (entry " << i << ") with " << c->nbins << " bins" << endl;
        bFound = true;
        
        /////////////////////////////////////////////////////////
        // get weighted rate histograms
        //
        //   used for energy binning
        
        // copy weighted rate histogram
        if( c->hWeightedRate )
        {
            char hname[1000];
            sprintf( hname, "%s_%s_%d_%d_%d", c->hWeightedRate->GetName(), iMCPara->fName.c_str(),
                     ( int )( iMCPara->woff * 1000 ), iMCPara->az, ( int )( iMCPara->index * 100 ) );
            iMCPara->hWeightedRate->SetName( hname );
            iMCPara->hWeightedRate->SetTitle( c->hWeightedRate->GetTitle() );
            iMCPara->hWeightedRate->SetBins( c->hWeightedRate->GetNbinsX(), c->hWeightedRate->GetXaxis()->GetXmin(), c->hWeightedRate->GetXaxis()->GetXmax() );
            for( int i = 0; i <= c->hWeightedRate->GetNbinsX(); i++ )
            {
                iMCPara->hWeightedRate->SetBinContent( i, c->hWeightedRate->GetBinContent( i ) );
                iMCPara->hWeightedRate->SetBinError( i, c->hWeightedRate->GetBinError( i ) );
            }
            
            int i_z_weightedFilledBins = 0;
            // set energy binning (note: this is fixed by hWeightedRate histogram)
            for( int i = 1; i <= iMCPara->hWeightedRate->GetNbinsX(); i++ )
            {
                if( iMCPara->hWeightedRate->GetBinContent( i ) )
                {
                    iMCPara->energy.push_back( iMCPara->hWeightedRate->GetBinCenter( i ) );
                    iMCPara->energy_lowEdge.push_back( iMCPara->hWeightedRate->GetBinLowEdge( i ) );
                    iMCPara->energy_upEdge.push_back( iMCPara->hWeightedRate->GetBinLowEdge( i ) + iMCPara->hWeightedRate->GetBinWidth( i ) );
                    i_z_weightedFilledBins++;
                }
            }
            if( fDebug )
            {
                cout << "\t number of filled bins in weighted rate histogram: " << i_z_weightedFilledBins << endl;
            }
            
            /////////////////////////////////////////////////////////
            // fill effective areas vs true energy (log10!)
            
            for( unsigned int i = 0; i < iMCPara->energy.size(); i++ )
            {
                int z = 0;
                iMCPara->effArea.push_back( 0. );
                iMCPara->effArea_error.push_back( 0. );
                for( int n = 0; n < c->nbins; n++ )
                {
                    if( c->e0[n] > iMCPara->energy_lowEdge[i] && c->e0[n] < iMCPara->energy_upEdge[i] )
                    {
                        iMCPara->effArea.back() += c->eff[n];
                        iMCPara->effArea_error.back() += (c->eff_error[n]*c->eff_error[n] );
                        z++;
                    }
                }
                if( z > 0. )
                {
                    iMCPara->effArea.back() /= z;
                    iMCPara->effArea_error.back() = sqrt( iMCPara->effArea_error.back() ) / z;
                }
            }
            // get global energy binning
            iMCPara->effArea_Ebins = iMCPara->hWeightedRate->GetNbinsX();
            iMCPara->effArea_Emin  = iMCPara->hWeightedRate->GetXaxis()->GetXmin();
            iMCPara->effArea_Emax  = iMCPara->hWeightedRate->GetXaxis()->GetXmax();
            cout << "\t global binning for effective areas: ";
            cout << " nbin: " << iMCPara->effArea_Ebins;
            cout << " emin: " << iMCPara->effArea_Emin;
            cout << " emax: " << iMCPara->effArea_Emax;
            cout << endl;
            cout << "\t effective area vector size: " << iMCPara->energy.size();
            if( iMCPara->energy.size() > 0 )
            {
                cout << " [" << iMCPara->energy[0];
            }
            if( iMCPara->energy.size() > 1 )
            {
                cout << "," << iMCPara->energy[iMCPara->energy.size() - 1] << "]";
            }
            cout << endl;
        }
        else
        {
            cout << "\t no global binning for effective areas found" << endl;
            iMCPara->effArea_Emin = 0.;
        }
        break;
    }
    // remove single filled bins in effective areas
    bool bGoodBin = false;
    if( iMCPara->energy.size() > 1 )
    {
        for( unsigned int i = iMCPara->energy.size() - 1; i > 0; i-- )
        {
            if( TMath::Abs( iMCPara->energy_lowEdge[i] - iMCPara->energy_upEdge[i - 1] ) > 1.e-3 && !bGoodBin )
            {
                iMCPara->energy.erase( iMCPara->energy.begin() + i, iMCPara->energy.end() );
                iMCPara->energy_upEdge.erase( iMCPara->energy_upEdge.begin() + i, iMCPara->energy_upEdge.end() );
                iMCPara->energy_lowEdge.erase( iMCPara->energy_lowEdge.begin() + i, iMCPara->energy_lowEdge.end() );
                iMCPara->effArea.erase( iMCPara->effArea.begin() + i, iMCPara->effArea.end() );
                iMCPara->effArea_error.erase( iMCPara->effArea_error.begin() + i, iMCPara->effArea_error.end() );
            }
            else
            {
                bGoodBin = true;
            }
        }
    }
    if( !bFound )
    {
        cout << "\t no effective area found!" << endl;
    }
    // (end of treatment of effective areas)
    
    ///////////////////////////////////////////////////////////
    // copy response matrix for energy reconstruction
    if( c->hResponseMatrix )
    {
        char hname[1000];
        sprintf( hname, "%s_%s_%d_%d_%d", c->hResponseMatrix->GetName(), iMCPara->fName.c_str(),
                 ( int )( iMCPara->woff * 1000 ), iMCPara->az, ( int )( iMCPara->index * 100 ) );
        iMCPara->hResponseMatrix->SetName( hname );
        iMCPara->hResponseMatrix->SetTitle( c->hResponseMatrix->GetTitle() );
        iMCPara->hResponseMatrix->SetBins( c->hResponseMatrix->GetNbinsX(),
                                           c->hResponseMatrix->GetXaxis()->GetXmin(), c->hResponseMatrix->GetXaxis()->GetXmax(),
                                           c->hResponseMatrix->GetNbinsY(),
                                           c->hResponseMatrix->GetYaxis()->GetXmin(), c->hResponseMatrix->GetYaxis()->GetXmax() );
        TH2D* i_hResponse = c->hResponseMatrix;
        // for cosmic rays: interpolate response matrix
        if( iMCPara->fParticleID != 1 && iMCPara->fParticleID != 2 )
        {
            cout << "\t interpolating response matrix for particle type " << iMCPara->fParticleID << endl;
            i_hResponse = ( TH2D* )VHistogramUtilities::interpolateResponseMatrix( c->hResponseMatrix );
        }
        if( i_hResponse )
        {
            for( int i = 0; i <= i_hResponse->GetNbinsX(); i++ )
            {
                for( int j = 0; j <= i_hResponse->GetNbinsY(); j++ )
                {
                    iMCPara->hResponseMatrix->SetBinContent( i, j, i_hResponse->GetBinContent( i, j ) );
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // read MC parameters in case run header is available
    // (e.g. diffuse scatter angle used in simtel simulations)
    //////////////////////////////////////////////////////////////////////////////////////
    VMonteCarloRunHeader* iMCHeader = ( VMonteCarloRunHeader* )fEff.Get( "MC_runheader" );
    if( iMCHeader )
    {
        cout << "reading Monte Carlo header" << endl;
        // diffuse scatter angle
        if( iMCHeader->diffuse )
        {
            iMCPara->theta2_MCScatterAngle = iMCHeader->viewcone[1] * iMCHeader->viewcone[1];
            cout << "\t setting diffuse scattering angle to " << sqrt( iMCPara->theta2_MCScatterAngle ) << " [deg]" << endl;
        }
    }
    else
    {
        cout << "===============================================================================" << endl;
        cout << "WARNING: no Monte Carlo header found; setting diffuse scattering angle to 4 deg!" << endl;
        cout << "===============================================================================" << endl;
        iMCPara->theta2_MCScatterAngle = 4.*4.;
    }
    iMCPara->SolidAngle_MCScatterAngle = 2. * TMath::Pi() * ( 1. - cos( sqrt( iMCPara->theta2_MCScatterAngle ) * TMath::Pi() / 180. ) );
    
    //////////////////////////////////////////////////////////////////////////////////////
    // calculate solid angle from analysis cuts (or from input parameter theta2)
    // (independent of energy binning of differential sensitivity curves)
    //////////////////////////////////////////////////////////////////////////////////////
    VGammaHadronCuts* iCuts = ( VGammaHadronCuts* )fEff.Get( "GammaHadronCuts" );
    if( iCuts )
    {
        // theta2 might be energy dependent
        cout << "calculating solid angle from analysis cuts (might be energy dependent)" << endl;
        cout << "\t TMVA CUT " << iCuts->getDirectionCutSelector() << endl;
        iMCPara->theta2_min = iCuts->getTheta2Cut_min();       // theta2 min assumed to be energy independent
        if( iMCPara->theta2_min < 0. )
        {
            iMCPara->theta2_min = 0.;
        }
        
        // solid angle vs energy is stored in a TGraph
        // (note that lower cut on theta2 is not a function of energy)
        iMCPara->gSolidAngle_DirectionCut_vs_EnergylgTeV = new TGraph( 1000 );
        iMCPara->gTheta2Cuts_vsEnergylgTeV = new TGraph( 1000 );
        
        double e = 0.;
        double iSolidAngle = 0.;
        double itheta2 = 0.;
        for( int i = 0; i < iMCPara->gSolidAngle_DirectionCut_vs_EnergylgTeV->GetN(); i++ )
        {
            e = fEnergy_min_Log + i * ( fEnergy_max_Log - fEnergy_min_Log ) / iMCPara->gSolidAngle_DirectionCut_vs_EnergylgTeV->GetN();
            
            // VGammaHadronCuts::getTheta2Cut_max work with lin E
            itheta2 = iCuts->getTheta2Cut_max( TMath::Power( 10., e ) );
            if( itheta2 > 0. )
            {
                iSolidAngle  = 2. * TMath::Pi() * ( 1. - cos( sqrt( itheta2 ) * TMath::Pi() / 180. ) );
                if( iMCPara->theta2_min > 0. )
                {
                    iSolidAngle -= 2. * TMath::Pi() * ( 1. - cos( sqrt( iMCPara->theta2_min ) * TMath::Pi() / 180. ) );
                }
            }
            else
            {
                cout << "Error: no theta2 cut given " << endl;
                return false;
            }
            iMCPara->gSolidAngle_DirectionCut_vs_EnergylgTeV->SetPoint( i, e, iSolidAngle );
            iMCPara->gTheta2Cuts_vsEnergylgTeV->SetPoint( i, e, itheta2 );
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // read success of TMVA cut optimization from gamma/hadron cuts
    fTMVAEvaluatorResults = ( VTMVAEvaluatorResults* )fEff.Get( "TMVAEvaluatorResults" );
    if( fTMVAEvaluatorResults )
    {
        cout << "reading TMVAEvaluatorResults from effective area file (";
        cout << fTMVAEvaluatorResults->fTMVAData.size() << " bins)" << endl;
    }
    
    return true;
}

void VSensitivityCalculator::plotSignificanceParameters( TCanvas* cSensitivity )
{
    if( !cSensitivity )
    {
        cout << "warning: no sensitivity canvas" << endl;
        return;
    }
    
    char hname[400];
    sprintf( hname, "%.0f hours, %.0f #sigma, >%.0f events", fObservationTime_h, fSignificance_min, fEvents_min );
    TLatex* iL = new TLatex( 0.17, 0.15, hname );
    iL->SetTextSize( 0.03 );
    iL->SetNDC();
    iL->Draw();
}

/*

    prepare three canvases for debugging plots
*/
void VSensitivityCalculator::prepareDebugPlots()
{
    vector< string > iCanvasTitle;
    iCanvasTitle.push_back( "debug plot: number of on/off (total) " );
    iCanvasTitle.push_back( "debug plot: effective areas" );
    iCanvasTitle.push_back( "debug plot: number of on/off (per particle type) " );
    char hname[800];
    for( unsigned int i = 0; i < iCanvasTitle.size(); i++ )
    {
        sprintf( hname, "%s_%u", fPlotDebugName.c_str(), i );
        cPlotDebug.push_back( new TCanvas( hname, iCanvasTitle[i].c_str(), 20 + i * 350, 650, 900, 600 ) );
        
        cPlotDebug.back()->SetGridx( 0 );
        cPlotDebug.back()->SetGridy( 0 );
        cPlotDebug.back()->SetLogy( 1 );
        cPlotDebug.back()->Draw();
    }
}


void VSensitivityCalculator::plotEffectiveArea()
{
    if( cPlotDebug.size() != 3 )
    {
        return;
    }
    
    if( !cPlotDebug[1] || !cPlotDebug[1]->cd() )
    {
        return;
    }
    
    vector< TGraphErrors* > g;
    TLegend* iL = new TLegend( 0.6, 0.6, 0.85, 0.85 );
    
    int z = 0;
    map< unsigned int, VSensitivityCalculatorDataResponseFunctions* >::iterator i_MCData_iterator;
    for( i_MCData_iterator = fMC_Data.begin(); i_MCData_iterator != fMC_Data.end(); ++i_MCData_iterator )
    {
        g.push_back( new TGraphErrors( 1 ) );
        g.back()->SetMinimum( 0.1 );
        g.back()->SetMaximum( 5.e6 );
        setGraphPlottingStyle( g.back(), z + 1, 1, 20 + z, 2 );
        
        for( unsigned int i = 0; i < ( *i_MCData_iterator ).second->energy.size(); i++ )
        {
            g.back()->SetPoint( i, ( *i_MCData_iterator ).second->energy[i], ( *i_MCData_iterator ).second->effArea[i] );
            g.back()->SetPointError( i, 0., ( *i_MCData_iterator ).second->effArea_error[i] );
        }
        if( g.back() )
        {
            if( z == 0 )
            {
                g.back()->Draw( "ap" );
                g.back()->GetHistogram()->SetXTitle( "log_{10} energy_{true} [TeV]" );
                g.back()->GetHistogram()->SetYTitle( "effective area [m^{2}]" );
            }
            else
            {
                g.back()->Draw( "p" );
            }
            iL->AddEntry( g.back(), ( *i_MCData_iterator ).second->fName.c_str(), "pl" );
        }
        z++;
    }
    if( iL )
    {
        iL->Draw();
    }
    
    cPlotDebug[1]->Update();
    
}

void VSensitivityCalculator::fillBackgroundParticleNumbers( vector< VDifferentialFluxData > iDifferentialFlux,
        map< unsigned int, vector< double > > i_flux_NOff,
        map< unsigned int, vector< double > > i_flux_NOffError )
{
    gProtonRate   = new TGraphAsymmErrors( 1 );
    gElectronRate = new TGraphAsymmErrors( 1 );
    
    if( i_flux_NOff.find( 14 ) != i_flux_NOff.end() && i_flux_NOffError.find( 14 ) != i_flux_NOffError.end() )
    {
        for( unsigned int i = 0; i < iDifferentialFlux.size(); i++ )
        {
            if( i < i_flux_NOff[14].size() )
            {
                gProtonRate->SetPoint( i, log10( iDifferentialFlux[i].Energy ), i_flux_NOff[14][i] );
            }
            if( i < i_flux_NOffError[14].size() )
            {
                gProtonRate->SetPointError( i, 0., 0., i_flux_NOffError[14][i], i_flux_NOffError[14][i] );
            }
        }
    }
    if( i_flux_NOff.find( 2 ) != i_flux_NOff.end() && i_flux_NOffError.find( 2 ) != i_flux_NOffError.end() )
    {
        for( unsigned int i = 0; i < iDifferentialFlux.size(); i++ )
        {
            if( i < i_flux_NOff[2].size() )
            {
                gElectronRate->SetPoint( i, log10( iDifferentialFlux[i].Energy ), i_flux_NOff[2][i] );
            }
            if( i < i_flux_NOffError[2].size() )
            {
                gElectronRate->SetPointError( i, 0., 0., i_flux_NOffError[2][i], i_flux_NOffError[2][i] );
            }
        }
    }
}

/*!

    plot particle numbers

*/
void VSensitivityCalculator::plotDebugPlotsBackgroundParticleNumbers( vector< VDifferentialFluxData > iDifferentialFlux,
        map< unsigned int, vector< double > > i_flux_NOff,
        map< unsigned int, vector< double > > i_flux_NOffError )
{
    if( cPlotDebug.size() != 3 )
    {
        return;
    }
    
    if( !cPlotDebug[2] || !cPlotDebug[2]->cd() )
    {
        return;
    }
    
    TGraphAsymmErrors* gNon = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* gNoff = new TGraphAsymmErrors( 1 );
    gNon->SetMinimum( 0.0001 );
    gNon->SetMaximum( 1.e6 );
    setGraphPlottingStyle( gNon, 1, 1, 20, 2 );
    setGraphPlottingStyle( gNoff, 1, 1, 24, 2 );
    
    unsigned int z = 0;
    for( unsigned int i = 0; i < iDifferentialFlux.size(); i++ )
    {
        if( iDifferentialFlux[i].Energy > 0. )
        {
            gNon->SetPoint( z, log10( iDifferentialFlux[i].Energy ), iDifferentialFlux[i].NOn * fObservationTime_h * 60. );
            gNoff->SetPoint( z, log10( iDifferentialFlux[i].Energy ), iDifferentialFlux[i].NOff * fObservationTime_h * 60. );
            z++;
        }
    }
    if( fDebug )
    {
        cout << "particle numbers: signal (on) region): " << endl;
        gNon->Print();
        cout << "particle numbers: signal (off) region): " << endl;
        gNoff->Print();
    }
    
    gNon->Draw( "ap" );
    gNon->GetHistogram()->SetXTitle( "log_{10} energy [TeV]" );
    gNon->GetHistogram()->SetYTitle( "number of on/off events" );
    gNoff->Draw( "p" );
    
    vector< TGraphErrors* > g;
    
    map< unsigned int, vector< double > >::iterator i_flux_NOff_iter;
    map< unsigned int, vector< double > >::iterator i_flux_NOffError_iter;
    
    z = 0;
    for( i_flux_NOff_iter = i_flux_NOff.begin(); i_flux_NOff_iter != i_flux_NOff.end(); i_flux_NOff_iter++ )
    {
        g.push_back( new TGraphErrors( 1 ) );
        g.back()->SetMinimum( 1.e-4 );
        g.back()->SetMaximum( 1.e6 );
        setGraphPlottingStyle( g.back(), z + 2, 1, 21 + z, 2 );
        z++;
    }
    
    for( unsigned int i = 0; i < iDifferentialFlux.size(); i++ )
    {
        z = 0;
        for( i_flux_NOff_iter = i_flux_NOff.begin(); i_flux_NOff_iter != i_flux_NOff.end(); ++i_flux_NOff_iter )
        {
            if( i < ( *i_flux_NOff_iter ).second.size() )
            {
                g[z]->SetPoint( i, log10( iDifferentialFlux[i].Energy ), ( *i_flux_NOff_iter ).second[i] * fObservationTime_h * 60. );
            }
            z++;
        }
        z = 0;
        for( i_flux_NOffError_iter = i_flux_NOffError.begin(); i_flux_NOffError_iter != i_flux_NOffError.end(); ++i_flux_NOffError_iter )
        {
            if( i < ( *i_flux_NOffError_iter ).second.size() )
            {
                g[z]->SetPointError( i, 0., ( *i_flux_NOffError_iter ).second[i] );
            }
            z++;
        }
    }
    
    for( unsigned int i = 0; i < g.size(); i++ )
    {
        g[i]->Draw( "p" );
    }
    
    cPlotDebug[2]->Update();
}

/*
 *
 *   cp events rates from differential flux vector into a graph and write them to the output file given
 *
 */
void VSensitivityCalculator::fillParticleNumbersGraphs( vector< VDifferentialFluxData > iDifferentialFlux, double alpha, double dE_log10 )
{
    gSignalRate = new TGraphAsymmErrors( 1 );
    gBGRate     = new TGraphAsymmErrors( 1 );
    setGraphPlottingStyle( gSignalRate, 1, 1, 20, 2 );
    setGraphPlottingStyle( gBGRate, 2, 1, 21, 2 );
    
    int z = 0;
    for( unsigned int i = 0; i < iDifferentialFlux.size(); i++ )
    {
        if( iDifferentialFlux[i].Energy > 0. && fObservationTime_h > 0. && dE_log10 > 0. )
        {
            gSignalRate->SetPoint( z, log10( iDifferentialFlux[i].Energy ),
                                   ( iDifferentialFlux[i].NOn - iDifferentialFlux[i].NOff * alpha ) / ( fObservationTime_h * 60. ) / dE_log10 );
            gSignalRate->SetPointEXhigh( z, log10( iDifferentialFlux[i].Energy_upEdge ) - log10( iDifferentialFlux[i].Energy ) );
            gSignalRate->SetPointEXlow( z, log10( iDifferentialFlux[i].Energy ) - log10( iDifferentialFlux[i].Energy_lowEdge ) );
            gBGRate->SetPoint( z, log10( iDifferentialFlux[i].Energy ),
                               iDifferentialFlux[i].NOff * alpha / ( fObservationTime_h * 60. ) / dE_log10 );
            gBGRate->SetPointEXhigh( z, log10( iDifferentialFlux[i].Energy_upEdge ) - log10( iDifferentialFlux[i].Energy ) );
            gBGRate->SetPointEXlow( z, log10( iDifferentialFlux[i].Energy ) - log10( iDifferentialFlux[i].Energy_lowEdge ) );
            z++;
        }
    }
    // write graphs with on and off events to disk
    if( fDebugParticleNumberFile.size() > 0 )
    {
        cout << "writing graphs (npoints=" << z << ") with on/off events to file: " << fDebugParticleNumberFile << endl;
        cout << "  (alpha=" << alpha << ", dE_log10=" << dE_log10;
        cout << ", time conversion units " << fObservationTime_h * 60. << ")" << endl;
        
        TFile iF( fDebugParticleNumberFile.c_str(), "RECREATE" );
        if( iF.IsZombie() )
        {
            cout << "VSensitivityCalculator::fillParticleNumbersGraphs: error opening particle number file (write): " << endl;
            cout << fDebugParticleNumberFile << endl;
            return;
        }
        gSignalRate->SetName( "gSignalRate" );
        gBGRate->SetName( "gBGRate" );
        if( TMath::Abs( dE_log10 - 1. ) > 1.e-4 )
        {
            gSignalRate->SetTitle( "differential signal rate (1/min)" );
            gBGRate->SetTitle( "differential background rate (1/min)" );
        }
        else
        {
            gSignalRate->SetTitle( "signal rate (1/min)" );
            gBGRate->SetTitle( "background rate (1/min)" );
        }
        gSignalRate->Write();
        gBGRate->Write();
        
        iF.Close();
    }
}


void VSensitivityCalculator::plotDebugPlotsParticleNumbers()
{
    if( !gSignalRate || !gBGRate )
    {
        return;
    }
    
    gSignalRate->SetMinimum( 0.0001 );
    gSignalRate->SetMaximum( 1.e2 );
    setGraphPlottingStyle( gSignalRate, 1, 1, 20, 2 );
    setGraphPlottingStyle( gBGRate, 2, 1, 21, 2 );
    
    if( cPlotDebug.size() != 3 )
    {
        return;
    }
    
    if( !cPlotDebug[0] || !cPlotDebug[0]->cd() )
    {
        return;
    }
    gSignalRate->Draw( "ap" );
    gSignalRate->GetHistogram()->SetXTitle( "log_{10} energy [TeV]" );
    gSignalRate->GetHistogram()->SetYTitle( "on/off rate [1/min]" );
    gBGRate->Draw( "p" );
    
    cPlotDebug[0]->Update();
}

/*!
    read values for Crab energy spectra from disk
*/
bool VSensitivityCalculator::setEnergySpectrumfromLiterature( string iFile, unsigned int iID )
{
    fEnergySpectrumfromLiterature = new VEnergySpectrumfromLiterature( iFile, false );
    if( fEnergySpectrumfromLiterature->isZombie() )
    {
        return false;
    }
    
    setEnergySpectrumfromLiterature_ID( iID );
    
    return true;
}

void VSensitivityCalculator::plotSensitivityLimitations( TCanvas* c, double iYValue )
{
    if( !c )
    {
        cout << "VSensitivityCalculator::plotSensitivityLimitations: error, canvas not found" << endl;
        return;
    }
    c->cd();
    
    // get maximum in histogram
    if( iYValue < -100. )
    {
        if( !hnull )
        {
            hnull = ( TH1D* )c->GetListOfPrimitives()->FindObject( "hnullSens" );
        }
        if( hnull )
        {
            iYValue = 0.5 * hnull->GetMaximum();
        }
    }
    
    cout << "plotSensitivityLimitations: ";
    cout << fMinEventsLimited.size() << "\t";
    cout << fSignificanceLimited.size() << "\t" << fMinBackgroundEventsLimited.size();
    cout << "\t" << fMinNoBackground.size() << endl;
    
    // minimum number of events
    if( fMinEventsLimited.size() > 0 )
    {
        cout << "Minimum number of events limited: red " << endl;
        map< int, double >::const_iterator itx;
        for( itx = fMinEventsLimited.begin(); itx != fMinEventsLimited.end(); itx++ )
        {
            double energy = ( double( ( *itx ).first ) ) / 1.e3;
            if( ( *itx ).second > 0. )
            {
                TGraphErrors* g = new TGraphErrors( 1 );
                g->SetPoint( 0, energy, ( *itx ).second );
                g->SetPointError( 0, 0.5 * fEnergy_dE_log10, 0. );
                g->SetLineColor( 2 );
                g->SetLineWidth( 3 );
                g->Draw( "LZ" );
            }
        }
    }
    // background event ratio limited
    if( fMinBackgroundEventsLimited.size() > 0 )
    {
        cout << "Background rate ratio limited: green" << endl;
        map< int, double >::const_iterator itx;
        for( itx = fMinBackgroundEventsLimited.begin(); itx != fMinBackgroundEventsLimited.end(); ++itx )
        {
            TGraphErrors* g = new TGraphErrors( 1 );
            double energy = ( double( ( *itx ).first ) ) / 1.e3;
            if( ( *itx ).second > 0. )
            {
                g->SetPoint( 0, energy, ( *itx ).second );
                g->SetPointError( 0, 0.5 * fEnergy_dE_log10, 0. );
                g->SetLineColor( 3 );
                g->SetLineWidth( 3 );
                g->Draw( "ZL" );
            }
        }
    }
    // significance
    if( fSignificanceLimited.size() > 0 )
    {
        cout << "Significance limited: blue" << endl;
        map< int, double >::const_iterator itx;
        for( itx = fSignificanceLimited.begin(); itx != fSignificanceLimited.end(); ++itx )
        {
            TGraphErrors* g = new TGraphErrors( 1 );
            double energy = ( double( ( *itx ).first ) ) / 1.e3;
            if( ( *itx ).second > 0. )
            {
                g->SetPoint( 0, energy, ( *itx ).second );
                g->SetPointError( 0, 0.5 * fEnergy_dE_log10, 0. );
                g->SetLineColor( 4 );
                g->SetLineWidth( 3 );
                g->Draw( "LZ" );
            }
        }
    }
    // no background events
    if( fMinNoBackground.size() > 0 )
    {
        cout << "No background events limited: purple" << endl;
        map< int, double >::const_iterator itx;
        for( itx = fMinNoBackground.begin(); itx != fMinNoBackground.end(); ++itx )
        {
            TGraphErrors* g = new TGraphErrors( 1 );
            double energy = ( double( ( *itx ).first ) ) / 1.e3;
            if( ( *itx ).second > 0. )
            {
                g->SetPoint( 0, energy, ( *itx ).second );
                g->SetPointError( 0, 0.5 * fEnergy_dE_log10, 0. );
                g->SetLineColor( 6 );
                g->SetLineWidth( 3 );
                g->Draw( "LZ" );
            }
        }
    }
    
}

bool VSensitivityCalculator::setMonteCarloParametersCTA_MC( string iCTA_MCFile, double iMCCTA_cameraoffset_deg,
        string iSpectralParameterFile, unsigned int iSpectralParameterID )
{
    if( !setEnergySpectrumfromLiterature( iSpectralParameterFile, iSpectralParameterID ) )
    {
        return false;
    }
    
    fMCCTA_cameraoffset_deg = iMCCTA_cameraoffset_deg;
    fMCCTA_File = iCTA_MCFile;
    
    return true;
}

/*
 * fill a histogram with background rate per square degree from a background rate graph
 *
 */
bool VSensitivityCalculator::fillBackroundvsSquareDegree( TGraphAsymmErrors* i_R, TH1F* iH_sqDeg )
{
    if( !i_R || !iH_sqDeg )
    {
        return false;
    }
    
    if( fMC_Data.find( 1 ) != fMC_Data.end() && fMC_Data[1]->gTheta2Cuts_vsEnergylgTeV )
    {
        double x = 0.;
        double y = 0.;
        for( int i = 0; i < i_R->GetN(); i++ )
        {
            i_R->GetPoint( i, x, y );
            if( y > 0. )
            {
                double iSolidAngle = fMC_Data[1]->gSolidAngle_DirectionCut_vs_EnergylgTeV->Eval( x );
                if( iSolidAngle > 0. )
                {
                    y /= iSolidAngle * TMath::RadToDeg() * TMath::RadToDeg();
                    iH_sqDeg->SetBinContent( iH_sqDeg->FindBin( x ), y / 60. );
                    iH_sqDeg->SetBinError( iH_sqDeg->FindBin( x ), 0.5 * ( i_R->GetErrorYlow( i ) + i_R->GetErrorYhigh( i ) )
                                           / ( iSolidAngle * TMath::DegToRad() * TMath::DegToRad() ) / 60. );
                }
            }
        }
    }
    else
    {
        return false;
    }
    
    return true;
}

/*

   fill sensitivity histograms for CTA WP Phys IRF files

*/
bool VSensitivityCalculator::fillSensitivityHistograms( TH1F* iSensitivity, TH1F* iBGRate, TH1F* iBGRateSqDeg,
        TH1F* iProtonRate, TH1F* iProtonRateSqDeg,
        TH1F* iElectronRate, TH1F* iElectronRateSqDeg,
        bool iHighEnergyFilling,
        TH1F* iTheta2, TH1F* iTheta )
{
    if( iSensitivity && gSensitivityvsEnergy )
    {
        fillSensitivityHistogramfromGraph( gSensitivityvsEnergy, iSensitivity, 1. );
    }
    // make sure the high-energy filling only is done when sensitivities are filled
    else
    {
        iHighEnergyFilling = false;
    }
    if( iBGRate && gBGRate )
    {
        // rate histogram are filled in 1/s, graphs are in 1/min
        fillSensitivityHistogramfromGraph( gBGRate, iBGRate, 1. / 60. );
    }
    fillBackroundvsSquareDegree( gBGRate, iBGRateSqDeg );
    if( iProtonRate && gProtonRate )
    {
        fillSensitivityHistogramfromGraph( gProtonRate, iProtonRate, 1. / 60. );
        fillBackroundvsSquareDegree( gProtonRate, iProtonRateSqDeg );
    }
    if( iElectronRate && gElectronRate )
    {
        fillSensitivityHistogramfromGraph( gElectronRate, iElectronRate, 1. / 60. );
        fillBackroundvsSquareDegree( gElectronRate, iElectronRateSqDeg );
    }
    if( iTheta2
            && fMC_Data.find( 1 ) != fMC_Data.end() && fMC_Data[1]->gTheta2Cuts_vsEnergylgTeV )
    {
        fillSensitivityHistogramfromGraph( fMC_Data[1]->gTheta2Cuts_vsEnergylgTeV, iTheta2, 1., false );
    }
    if( iTheta
            && fMC_Data.find( 1 ) != fMC_Data.end() && fMC_Data[1]->gTheta2Cuts_vsEnergylgTeV )
    {
        fillSensitivityHistogramfromGraph( fMC_Data[1]->gTheta2Cuts_vsEnergylgTeV, iTheta, 1., false, true );
    }
    
    // go again over the high energy part of the sensitivity curve and check conditions for sensitivity:
    // is signal number critera fullfilled early enough so that it can be used?
    if( iHighEnergyFilling )
    {
        double x = 0.;
        double y = 0.;
        int i_x = 0;
        bool iFillBins = false;
        map< int, double >::const_iterator itx;
        for( itx = fMinEventsLimited.begin(); itx != fMinEventsLimited.end(); ++itx )
        {
            x = ( double( ( *itx ).first ) ) / 1.e3;
            y = ( ( *itx ).second );
            i_x = iSensitivity->FindBin( x );
            if( x > 0. && iSensitivity->GetBinContent( i_x ) > 0. )
            {
                if( iSensitivity->GetBinError( i_x ) == 0.
                        || iSensitivity->GetBinError( i_x ) / iSensitivity->GetBinContent( i_x ) < 0.25 )
                {
                    if( TMath::Abs( iSensitivity->GetBinContent( i_x ) - iSensitivity->GetBinError( i_x ) - y ) / y < 0.10 )
                    {
                        iFillBins = true;
                    }
                }
            }
            // fill (don't allow negative sensitivity)
            if( iFillBins && y > 0. )
            {
                iSensitivity->SetBinContent( i_x, y );
                cout << "Filling HE (signal counts): " << x << "\t" << y << "\t" << i_x << endl;
            }
        }
        
    }
    return true;
}

TCanvas* VSensitivityCalculator::plotSignalBackgroundRates( TCanvas* c, bool bPlotParticleBackground, double iRateMinimum, double iRateMaximum )
{
    if( !c )
    {
        c = new TCanvas( "cSignalBackgroundRates", "background rates", 0, 0, 400, 400 );
        c->SetGridx( 0 );
        c->SetGridx( 0 );
        c->SetLeftMargin( 0.15 );
        c->SetRightMargin( 0.07 );
        c->Draw();
    }
    else
    {
        c->cd();
    }
    
    if( !c->GetListOfPrimitives()->FindObject( "hnullSignalParticleRates" ) )
    {
        TH1D* h = new TH1D( "hnullSignalParticleRates", "", 10, fEnergy_min_Log, fEnergy_max_Log );
        h->SetStats( 0 );
        h->SetMaximum( iRateMaximum );
        h->SetMinimum( iRateMinimum );
        h->SetXTitle( "log_{10} energy [TeV]" );
        h->SetYTitle( "background rate [1/s]" );
        plot_nullHistogram( c, h, true, true, 1.7, TMath::Power( 10., fEnergy_min_Log ), TMath::Power( 10., fEnergy_max_Log ) );
        c->SetLogy( 1 );
    }
    
    if( gBGRate )
    {
        setGraphPlottingStyle( gBGRate );
        gBGRate->SetLineWidth( 1. );
        gBGRate->Draw( "pl" );
        
        if( gProtonRate && bPlotParticleBackground )
        {
            setGraphPlottingStyle( gProtonRate );
            gProtonRate->SetLineStyle( 2 );
            gProtonRate->Draw( "pl" );
        }
        if( gElectronRate && bPlotParticleBackground )
        {
            setGraphPlottingStyle( gElectronRate );
            gElectronRate->SetLineStyle( 3 );
            gElectronRate->Draw( "pl" );
        }
    }
    
    return c;
}

/*

   fill sensitivity histogram in CTA WP Phys IRF file style from the sensitivity graphs

*/
bool VSensitivityCalculator::fillSensitivityHistogramfromGraph( TGraph* g, TH1F* h, double iScale, bool bFillErrors, bool bSQRT )
{
    if( !h || !g )
    {
        return false;
    }
    
    double x = 0.;
    double y = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        if( y > 0. && iScale != 0. )
        {
            if( !bSQRT )
            {
                h->SetBinContent( h->FindBin( x ), y * iScale );
            }
            else
            {
                if( y * iScale > 0. )
                {
                    h->SetBinContent( h->FindBin( x ), sqrt( y * iScale ) );
                }
            }
            if( bFillErrors )
            {
                if( 0.5 * ( g->GetErrorYlow( i ) + g->GetErrorYhigh( i ) )*iScale < y * iScale )
                {
                    h->SetBinError( h->FindBin( x ), 0.5 * ( g->GetErrorYlow( i ) + g->GetErrorYhigh( i ) )*iScale );
                }
                else
                {
                    h->SetBinError( h->FindBin( x ), TMath::Min( g->GetErrorYlow( i ), g->GetErrorYhigh( i ) ) * iScale );
                }
            }
            else
            {
                h->SetBinError( h->FindBin( x ), 0. );
            }
        }
    }
    
    return true;
}

bool VSensitivityCalculator::fillSensitivityHistogramfromMap( map< int, double > m, TH1F* h )
{
    if( !h || m.size() == 0 )
    {
        return false;
    }
    
    double x = 0.;
    double y = 0.;
    map< int, double >::const_iterator itx;
    for( itx = m.begin(); itx != m.end(); ++itx )
    {
        x = ( double( ( *itx ).first ) ) / 1.e3;
        y = ( ( *itx ).second );
        if( y > 0. )
        {
            h->SetBinContent( h->FindBin( x ), y );
        }
    }
    
    return true;
}

/*

   fill histograms with sensitivity limits
   (for CTA-style root files

   0 = significance
   1 = signal number
   2 = background ratio
   3 = no off events

*/
bool VSensitivityCalculator::fillSensitivityLimitsHistograms( vector<TH1F*>& h )
{

    if( h.size() >= 4 )
    {
        fillSensitivityHistogramfromMap( fSignificanceLimited, h[0] );
        fillSensitivityHistogramfromMap( fMinEventsLimited, h[1] );
        fillSensitivityHistogramfromMap( fMinBackgroundEventsLimited, h[2] );
        fillSensitivityHistogramfromMap( fMinNoBackground, h[3] );
    }
    return false;
}

bool VSensitivityCalculator::checkCutOptimization( double iEnergy, bool iPrint )
{
    if( !fRequireCutsToBeOptimized )
    {
        if( iPrint )
        {
            cout << "VSensitivityCalculator::checkCutOptimization: " << fRequireCutsToBeOptimized << endl;
        }
        return true;
    }
    unsigned int no_data_counter = 0;
    
    if( fTMVAEvaluatorResults )
    {
        if( iEnergy > 0. )
        {
            iEnergy = log10( iEnergy );
        }
        else
        {
            if( iPrint )
            {
                cout << "VSensitivityCalculator::checkCutOptimization: negative energy " << endl;
            }
            return false;
        }
        for( unsigned int i = 0; i < fTMVAEvaluatorResults->fTMVAData.size(); i++ )
        {
            if( !fTMVAEvaluatorResults->fTMVAData[i] )
            {
                no_data_counter++;
                continue;
            }
            
            if( iEnergy > fTMVAEvaluatorResults->fTMVAData[i]->fEnergyCut_Log10TeV_min
                    && iEnergy < fTMVAEvaluatorResults->fTMVAData[i]->fEnergyCut_Log10TeV_max )
            {
                if( fTMVAEvaluatorResults->fTMVAData[i]->fTMVAOptimumCutValueFound )
                {
                    return true;
                }
                else
                {
                    cout << "VSensitivityCalculator::checkCutOptimization: ignore energy bin ";
                    cout << "[" << fTMVAEvaluatorResults->fTMVAData[i]->fEnergyCut_Log10TeV_min << ", ";
                    cout << fTMVAEvaluatorResults->fTMVAData[i]->fEnergyCut_Log10TeV_max << "] (";
                    cout << iEnergy << "): cuts not optimized " << endl;
                    return false;
                }
            }
        }
    }
    if( iPrint )
    {
        cout << "VSensitivityCalculator::checkCutOptimization: no fTMVAEvaluatorResults results" <<  endl;
        cout << "\t " << fTMVAEvaluatorResults;
        if( fTMVAEvaluatorResults )
        {
            cout << "\t" << fTMVAEvaluatorResults->fTMVAData.size();
        }
        cout << "\t no data counter: " << no_data_counter;
        cout << endl;
    }
    return false;
}

TGraph* VSensitivityCalculator::getObservationTimevsFluxGraph( unsigned int ID )
{
    if( fGraphObsvsTime.find( ID ) == fGraphObsvsTime.end() )
    {
        return 0;
    }
    
    return fGraphObsvsTime[ID];
}

/*
 * quick fix for single energy bins without sensitivity estimations due
 * to missing background: fill average of neighbouring bins
 *
 * note: fixed to work between 1 and 10 TeV
 *
 */
void VSensitivityCalculator::smoothenSensitivityGraph( TGraphAsymmErrors* g, TH1F* h )
{
    if( g && h )
    {
        double y = 0.;
        double yN = 0.;
        double yL = 0.;
        double yH = 0.;
        double x = 0.;

        for( int i = 1; i < g->GetN()-1; i++ )
        {
            g->GetPoint( i-1, x, yL );
            g->GetPoint( i+1, x, yH );
            g->GetPoint( i, x, y );

            if( x > 0. && x < 1. )
            {
                if( h->GetBinContent( h->FindBin( x ) ) < 1.e-20 )
                {
                    yN = 0.5* ( yL + yH );
                    g->SetPoint( i, x, yN );
                    cout << "Interpolated sensitivity for bin: ";
                    cout << x << "\t" << yN;
                    cout << " (was " << y << ")" << endl;
                }
            }
        }
    }
}
    


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// data class for response functions for a give primary
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VSensitivityCalculatorDataResponseFunctions::VSensitivityCalculatorDataResponseFunctions()
{
    fName = "";
    fParticleID = 0;
    fSpectralParameterFile = "";
    fSpectralParameterID = 0;
    fEffectiveAreaFile = "";
    ze = 0.;
    az = 0;
    woff = 0.;
    noise = 0;
    index = 0.;
    theta2_min = 0.;
    theta2_MCScatterAngle = 0.;
    gSolidAngle_DirectionCut_vs_EnergylgTeV = 0;
    gTheta2Cuts_vsEnergylgTeV = 0;
    SolidAngle_MCScatterAngle = 0.;
    alpha = 0.;
    effArea_Ebins = 0;
    effArea_Emin = 0.;
    effArea_Emax = 0.;
    energy_min_log = 0.;
    energy_max_log = 0.;
    hResponseMatrix = 0;
    hWeightedRate = 0;
}

VSensitivityCalculatorDataResponseFunctions::~VSensitivityCalculatorDataResponseFunctions()
{
    //   if( hResponseMatrix ) delete hResponseMatrix;
}

double VSensitivityCalculatorDataResponseFunctions::getSolidAngle_DirectionCut( double e )
{
    if( gSolidAngle_DirectionCut_vs_EnergylgTeV )
    {
        if( e > 0. )
        {
            return gSolidAngle_DirectionCut_vs_EnergylgTeV->Eval( log10( e ) );
        }
    }
    
    return -1.;
}

void VSensitivityCalculatorDataResponseFunctions::initializeHistograms( string iU )
{
    char hname[1000];
    // preliminary initialization
    sprintf( hname, "%s_%d_%d_%d_%d_%s", fName.c_str(), fParticleID, ( int )( woff * 1000 ), az, ( int )( index * 100 ), iU.c_str() );
    hResponseMatrix = new TH2D( hname, "B", 5, 0., 1., 5, 0., 1. );
    sprintf( hname, "%s1D_%d_%d_%d_%d_%s", fName.c_str(), fParticleID, ( int )( woff * 1000 ), az, ( int )( index * 100 ), iU.c_str() );
    hWeightedRate = new TH1D( hname, "B", 5, 0., 1. );
}

