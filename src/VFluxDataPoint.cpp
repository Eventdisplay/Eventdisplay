/*! \class VFluxDataPoint
    \brief data class for integral fluxes

    Assumption for small time binned analysis:

    - all runs must have same binning in time
    - intra-run modified alpha calculation is used

*/

#include "VFluxDataPoint.h"


VFluxDataPoint::VFluxDataPoint( string iName )
{
    setDebug( false );
    
    fName = iName;
    
    setSignificanceParameters();
    
    reset();
}

void VFluxDataPoint::reset()
{
    fDataFileName = "";
    fFluxUnitString = "not_set";
    fIsCombinedFluxData = false;
    fRunNumber = 0;
    fMJD_RunStart = 0.;
    fMJD_RunStop = 0.;
    fSecondsIntoRun = 0.;
    fMJD = 0.;
    fMJD_Start = 0.;
    fMJD_Stop = 0.;
    fOrbitalPhase_Start = -99;
    fOrbitalPhase_Stop = -99;
    fOrbitalPhase = -99;
    fExposure = 0.;
    fTimeBinDuration_sec = 0.;
    fExposure_deadTimeCorrected = 0.;
    fDeadTime = 0.;
    fTimeMask_open = false;
    fZe = 0.;
    fWobbleOffset = 0.;
    fPedvars = 0.;
    fNdiff = 0.;
    fNdiffE = 0.;
    fNdiff_Rolke = 0.;
    fRate = 0.;
    fRateE = 0.;
    fRate_Rolke = 0.;
    fRate_lo_1sigma = 0.;
    fRate_up_1sigma = 0.;
    fNon = 0.;
    fNoff = 0.;
    fAlpha = 0.;
    fAlpha_U = 0.;
    fAlpha_L = 0.;
    fTotFluxSum = 0.;
    fTotFluxWeight = 0.;
    fTotWobbleOffset = 0.;
    fTotZeSum    = 0.;
    fTotPedvars  = 0.;
    fTotEffArea  = 0.;
    
    fMinEnergy_TeV = -99.;
    fMaxEnergy_TeV = -99.;
    fE0_TeV        = -99.;
    fSpectralIndex = -99.;
    
    resetFluxValues();
    
    fDebugCounter = 0;
    fCalculateExpectedLimitsN = 0;
}

void VFluxDataPoint::resetFluxValues()
{
    fSignificantDataPoint = false;
    fSignificance = -99.;
    fUL = -99.;
    fCI_lo_1sigma = -99.;
    fCI_up_1sigma = -99.;
    fCI_lo_3sigma = -99.;
    fCI_up_3sigma = -99.;
    
    fEffArea_cm2 = 0.;
    fWeightedEffArea_cm2 = 0.;
    fFluxConstant = 0.;
    fFluxConstantE = 0.;
    fFlux = -99.;
    fFluxE = -99.;
    fFlux_Rolke = -99.;
    fFluxCI_1sigma = -99.;
    fFluxCI_lo_1sigma = -99.;
    fFluxCI_up_1sigma = -99.;
    fFluxCI_lo_3sigma = -99.;
    fFluxCI_up_3sigma = -99.;
    fFluxUL = -99.;
    
    fUL_Expected = -99;
    fFluxUL_Expected = -99;
    fFluxConstantUL_Expected = -99;
    
}

bool VFluxDataPoint::hasOrbitalPhases()
{
    return ( fOrbitalPhase_Start >= 0. && fOrbitalPhase_Stop <= 1. );
}

void VFluxDataPoint::printRunData( ostream& output )
{
    output << "RUN ";
    if( isCombinedDataElement() )
    {
        output << fRunNumber;
    }
    else
    {
        output << " (combined data)";
    }
    output << ",  Exposure: " << fExposure;
    output << " (dead time correct: " << fExposure_deadTimeCorrected << ")";
    output << ", Ze [deg]: " << fZe;;
    output << endl;
}

/*

    simple printing of flux calculation results

*/
void VFluxDataPoint::printResultSummary( ostream& output )
{
    if( isCombinedDataElement() )
    {
        output << Form( "MJD [%.4f, %.4f]:", fMJD_Start, fMJD_Stop );
        if( fOrbitalPhase_Start >= 0. && fOrbitalPhase_Stop <= 1. )
        {
            output << Form( "  (orbital phase: [%.4f, %.4f])", fOrbitalPhase_Start, fOrbitalPhase_Stop );
        }
    }
    else
    {
        output << Form( "RUN %d at MJD %.4f [%.4f, %.4f])", fRunNumber, fMJD, fMJD_Start, fMJD_Stop );
        if( fOrbitalPhase >= 0. && fOrbitalPhase <= 1. )
        {
            output << Form( "  (orbital phase: %.4f)", fOrbitalPhase );
        }
    }
    output << Form( ": %.1f sigma, ", fSignificance );
    if( isSignificantDataPoint() )
    {
        //	    output << Form( "Flux F(E > %.3f TeV) [cm^-2 s^-1] (P): %.2e +-%.2e", fMinEnergy_TeV, fFlux, fFluxE );
        output << Form( "Flux F(E > %.3f TeV) [cm^-2 s^-1] (R): %.2e [+%.2e, -%.2e]", fMinEnergy_TeV, fFlux_Rolke, fFluxCI_up_1sigma, fFluxCI_lo_1sigma );
    }
    else if( !fTimeMask_open )
    {
        output << Form( "Upper Flux Limit F(E > %.3f TeV, %d%% CL) [cm^-2 s^-1]: %.2e", fMinEnergy_TeV, ( int )( fUpperLimit * 100. ), fFluxUL );
    }
    else
    {
        output << " (Time bin fully excluded due to time cuts.) " ;
    }
    output << endl;
}

/*

    detailed printing of flux calculation results

*/
void VFluxDataPoint::printResults( ostream& output )
{
    if( isCombinedDataElement() )
    {
        cout << Form( "MJD [%.4f, %.4f]:", fMJD_Start, fMJD_Stop );
        if( fOrbitalPhase_Start >= 0. && fOrbitalPhase_Stop <= 1. )
        {
            cout << Form( "  (orbital phase: [%.4f, %.4f])", fOrbitalPhase_Start, fOrbitalPhase_Stop );
        }
        cout << endl;
    }
    else
    {
        cout << Form( "RUN %d at MJD %.4f [%.4f, %.4f])", fRunNumber, fMJD, fMJD_Start, fMJD_Stop );
        if( fOrbitalPhase >= 0. && fOrbitalPhase <= 1. )
        {
            cout << Form( "  (orbital phase: %.4f)", fOrbitalPhase );
        }
        cout << endl;
    }
    cout << Form( "\t ze [deg]: %.1f, ", fZe );
    cout << Form( "wobble offset [deg]: %.1f, ", fWobbleOffset );
    cout << Form( "pedvars: %.2f, ", fPedvars );
    cout << Form( "t_obs [min]: %.1f", fExposure / 60. );
    cout << Form( " (%.1f [s])", fExposure );
    cout << Form( ", dead time corrected [min]: %.1f", fExposure_deadTimeCorrected / 60. );
    cout << Form( " (%.1f [s])", fExposure_deadTimeCorrected );
    if( fTimeMask_open )
    {
        cout << Form( " (bin excluded from analysis (time mask))" );
    }
    cout << Form( ", dead time fraction: %.4f", fDeadTime );
    if( fWeightedEffArea_cm2 > 1. )
    {
        cout << Form( ", spectral weighted effective area [cm^2]: %.2e", fWeightedEffArea_cm2 );
    }
    output << endl;
    if( fNon > 0. ||  fNoff > 0. || fAlpha > 0. )
    {
        cout << Form( "\t Poisson: NOn: %d, NOff: %d, Norm: %.2f, NDiff: %.1f +-%.1f", ( int )fNon, ( int )fNoff, fAlpha , fNdiff, fNdiffE );
        cout << Form( ", Rate: %.3f +-%.3f gammas/min\n", fRate, fRateE );
        cout << Form( "\t Rolke et al: NOn: %d, NOff: %d, Norm: %.2f, NDiff: %.1f [+%.1f, -%.1f]", ( int )fNon, ( int )fNoff, fAlpha , fNdiff_Rolke, fCI_up_1sigma, fCI_lo_1sigma );
        cout << Form( ", Rate: %.3f [+%.3f, -%.3f] gammas/min\n", fRate_Rolke, fRate_up_1sigma, fRate_lo_1sigma );
    }
    if( fSignificance > -90. )
    {
        cout << Form( "\t Significance: %.1f\n", fSignificance );
    }
    if( isSignificantDataPoint() )
    {
        cout << Form( "\t Flux F(E > %.3f TeV) [cm^-2 s^-1] (P): %.2e +-%.2e\n", fMinEnergy_TeV, fFlux, fFluxE );
        cout << Form( "\t Flux F(E > %.3f TeV) [cm^-2 s^-1] (R): %.2e [+%.2e, -%.2e]\n", fMinEnergy_TeV, fFlux_Rolke, fFluxCI_up_1sigma, fFluxCI_lo_1sigma );
        // calculate flux in CU (only for >10 GeV points)
        // (note: we assume a power law throughout the energy range)
        if( fMinEnergy_TeV > 0.1 && fFlux > 0. )
        {
            cout << Form( "\t Flux F(E > %.3f TeV) [C.U. (Whipple)]: %.3f [+%.3f, -%.3f] \n ", fMinEnergy_TeV,
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy_TeV, fFlux ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy_TeV, fFluxCI_up_1sigma ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy_TeV, fFluxCI_lo_1sigma ) );
            cout << Form( "\t Flux F(E > %e Hz) [erg/s/cm^2]: %.2e [+%.2e, -%.2e]\n",
                          VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( fMinEnergy_TeV ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( fMinEnergy_TeV, fFlux_Rolke, true ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( fMinEnergy_TeV, fFluxCI_up_1sigma, true ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( fMinEnergy_TeV, fFluxCI_lo_1sigma, true ) );
        }
        cout << Form( "\t Differential flux normalisation at %.2f TeV [cm^-2 s^-1 TeV^-1]: %.2e\n", fE0_TeV, fFluxConstant );
    }
    else if( !fTimeMask_open )
    {
        cout << Form( "\t Upper limit on events: %.1f\n", fUL );
        cout << Form( "\t Upper Flux Limit F(E > %.3f TeV, %d%% CL) [cm^-2 s^-1]: %.2e\n", fMinEnergy_TeV, ( int )( fUpperLimit * 100. ), fFluxUL );
        if( fMinEnergy_TeV > 0. )
        {
            cout << Form( "\t Upper Flux Limit F(E > %.3f TeV, %d%% CL) [C.U. (Whipple)]: %.3f\n", fMinEnergy_TeV, ( int )( fUpperLimit * 100. ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy_TeV, fFluxUL ) );
        }
        // convert integral flux limit to a differential flux
        double iDF = fFluxUL * ( fAlpha + 1 ) / fMinEnergy_TeV;
        cout << Form( "\t Upper Flux Limit F(E > %e Hz, %d%% CL) [erg/s/cm^2]: %e\n",
                      VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( fMinEnergy_TeV ), ( int )( fUpperLimit * 100. ),
                      VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( fMinEnergy_TeV, iDF, true ) );
        cout << Form( "\t Upper Limit on flux normalisation at %.2f TeV [cm^-2 s^-1 TeV^-1]: %.2e\n", fE0_TeV, fFluxConstant );
    }
    if( fUL_Expected > -99 )
    {
        cout << Form( "\t Expected upper limit on events:  %.2f [+%.2f, -%.2f]\n", fUL_Expected, fUL_Expected_up_1sigma - fUL_Expected, fUL_Expected - fUL_Expected_lo_1sigma );
        cout << Form( "\t Expected upper Flux Limit F(E > %.3f TeV, %d%% CL) [cm^-2 s^-1]: %.2e [+%.2e, -%.2e]\n", fMinEnergy_TeV, ( int )( fUpperLimit * 100. ), fFluxUL_Expected, fFluxUL_Expected_up_1sigma, fFluxUL_Expected_lo_1sigma );
        if( fMinEnergy_TeV > 0. )
        {
            cout << Form( "\t Expected upper Flux Limit F(E > %.3f TeV, %d%% CL) [C.U. (Whipple)]: %.3f\n", fMinEnergy_TeV, ( int )( fUpperLimit * 100. ),
                          VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy_TeV, fFluxUL_Expected ) );
        }
        // convert integral flux limit to a differential flux
        double iDF = fFluxUL_Expected * ( fAlpha + 1 ) / fMinEnergy_TeV;
        cout << Form( "\t Expected upper Flux Limit F(E > %e Hz, %d%% CL) [erg/s/cm^2]: %e\n",
                      VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( fMinEnergy_TeV ), ( int )( fUpperLimit * 100. ),
                      VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( fMinEnergy_TeV, iDF, true ) );
        cout << Form( "\t Expected upper Limit on Flux Constant at %.2f TeV, %d%% CL) [cm^-2 s^-1 TeV^-1]: %.2e\n", fE0_TeV, ( int )( fUpperLimit * 100. ), fFluxConstantUL_Expected );
    }
}

/*

    calculate significances and upper limits using the on and off counts

*/
void VFluxDataPoint::calculateSignificancesAndUpperLimits()
{
    resetFluxValues();
    
    if( fAlpha <= 0. )
    {
        return;
    }
    
    // calculate N_diff
    fNdiff = fNon - fNoff * fAlpha;
    // assuming poisson errors
    fNdiffE = sqrt( fNon + fNoff * fAlpha * fAlpha );
    
    // calculate significance
    fSignificance = VStatistics::calcSignificance( fNon, fNoff, fAlpha, fLiMaEqu );
    
    TRolke i_Rolke;
    // calculate N_diff using Rolke et al method
    i_Rolke.SetCLSigmas( 1.e-4 );
    i_Rolke.SetPoissonBkgKnownEff( ( int )fNon, ( int )fNoff, 1. / fAlpha, 1. );
    fNdiff_Rolke = 0.5 * ( i_Rolke.GetLowerLimit() + i_Rolke.GetUpperLimit() );
    
    // calculate confidence intervals for fluxes
    i_Rolke.SetCLSigmas( 1. );
    i_Rolke.SetPoissonBkgKnownEff( ( int )fNon, ( int )fNoff, 1. / fAlpha, 1. );
    fCI_lo_1sigma = TMath::Abs( fNdiff_Rolke - i_Rolke.GetLowerLimit() );
    fCI_up_1sigma = TMath::Abs( i_Rolke.GetUpperLimit() - fNdiff_Rolke );
    i_Rolke.SetCLSigmas( 3. );
    i_Rolke.SetPoissonBkgKnownEff( ( int )fNon, ( int )fNoff, 1. / fAlpha, 1. );
    fCI_lo_3sigma = TMath::Abs( fNdiff_Rolke - i_Rolke.GetLowerLimit() );
    fCI_up_3sigma = TMath::Abs( i_Rolke.GetUpperLimit() - fNdiff_Rolke );
    // calculate upper flux
    fUL = VStatistics::calcUpperLimit( fNon, fNoff, fAlpha, fUpperLimit, fUpperLimitMethod, fBoundedLimits );
    
    // check if this is a significant point
    if( fSignificance < fThresholdSignificance || fNon < fMinEvents )
    {
        fSignificantDataPoint = false;
    }
    else
    {
        fSignificantDataPoint = true;
    }
    
    
    // rates
    if( fExposure > 0. )
    {
        fRate = fNdiff / fExposure * 60.;
        fRateE = fNdiffE / fExposure * 60.;
        fRate_Rolke = fNdiff_Rolke / fExposure * 60.;
        fRate_lo_1sigma = fCI_lo_1sigma / fExposure * 60.;
        fRate_up_1sigma = fCI_up_1sigma / fExposure * 60.;
    }
    else
    {
        fRate = 0.;
        fRateE = 0.;
        fRate_Rolke = 0.;
        fRate_lo_1sigma = 0.;
        fRate_up_1sigma = 0.;
    }
    calculateExpectedLimit();
}


/*

  calculate spectral weighted integral effective area

*/
bool VFluxDataPoint::calculateIntegralEffectiveArea( vector< double > energy_axis, vector< double > effArea )
{
    if( energy_axis.size() != effArea.size() )
    {
        cout << "VFluxDataPoint::calculateIntegralEffectiveArea: mismatch between energy and effective area vector" << endl;
        cout << "\t" << energy_axis.size() << "\t" << effArea.size() << endl;
        return false;
    }
    if( energy_axis.size() < 2 )
    {
        cout << "VFluxDataPoint::calculateIntegralEffectiveArea: energy/effective area vector too small" << endl;
        cout << "\t" << energy_axis.size() << "\t" << effArea.size() << endl;
        return false;
    }
    
    fEffArea_cm2 = 0.;
    
    double ieff_int = 0.;
    double ieff_mean = 0.;
    double e0 = 0.;
    double e1 = 0.;
    double e2 = 0.;
    
    // loop over all bins in effective area histogram
    for( unsigned int b = 1; b < energy_axis.size() - 1; b++ )
    {
        // get energies and effective areas
        e0 = energy_axis[b - 1];
        e2 = energy_axis[b + 1];
        e1 = energy_axis[b];
        ieff_mean = effArea[b];
        
        // energy bin (linear scale)
        double iE_low = TMath::Power( 10., ( e0 + e1 ) / 2. );
        double iE_up  = TMath::Power( 10., ( e1 + e2 ) / 2. );
        
        // check energy range
        if( iE_up < fMinEnergy_TeV ||  iE_low > fMaxEnergy_TeV )
        {
            continue;
        }
        if( iE_low < fMinEnergy_TeV )
        {
            iE_low = fMinEnergy_TeV;
        }
        if( iE_up  > fMaxEnergy_TeV )
        {
            iE_up = fMaxEnergy_TeV;
        }
        
        // width of energy bin (on linear scale)
        double dE = iE_up - iE_low;
        
        // integrate spectral weighted effective area
        // using trapezoidal rule
        ieff_int = ieff_mean * dE * 0.5 * ( TMath::Power( iE_up / fE0_TeV, fSpectralIndex ) + TMath::Power( iE_low / fE0_TeV, fSpectralIndex ) );
        // (correction term at E_low in order to maximize second derivative)
        ieff_int -=  1. / 12. * dE * dE * dE * ieff_mean * fSpectralIndex * ( fSpectralIndex - 1. ) / fE0_TeV / fE0_TeV * TMath::Power( iE_low / fE0_TeV, fSpectralIndex - 2. );
        
        fEffArea_cm2 += ieff_int;
    }
    // m^2 to cm^2
    fEffArea_cm2 *= 1.e4;
    
    return true;
}

/*

    calculate flux using the event counts, effective area and dead-time corrected observing time

*/
void VFluxDataPoint::calculateFlux()
{
    double iEffNorm = fEffArea_cm2 * fExposure_deadTimeCorrected;
    
    // calculate integral flux constant
    if( fDebug )
    {
        cout << "Calculating Integral fluxes within the restricted range [" << fMinEnergy_TeV << " TeV, "  << fMaxEnergy_TeV << " TeV]  " << endl;
    }
    double iN_D = -1. / ( fSpectralIndex + 1. );
    iN_D *= ( TMath::Power( fMinEnergy_TeV, fSpectralIndex + 1. ) - TMath::Power( fMaxEnergy_TeV, fSpectralIndex + 1. ) ) / TMath::Power( fE0_TeV, fSpectralIndex );
    
    fWeightedEffArea_cm2 = fEffArea_cm2 / iN_D;
    
    // calculate fluxes
    if( iEffNorm > 0. )
    {
        fFluxConstant     = fNdiff / iEffNorm;
        fFluxConstantE    = fNdiffE / iEffNorm;
        fFlux             = fFluxConstant * iN_D;
        fFluxE            = fFluxConstantE * iN_D;
        fFluxUL           = fUL / iEffNorm * iN_D;
        fFlux_Rolke       = fNdiff_Rolke / iEffNorm * iN_D;
        fFluxCI_lo_1sigma = fCI_lo_1sigma / iEffNorm * iN_D;
        fFluxCI_up_1sigma = fCI_up_1sigma / iEffNorm * iN_D;
        fFluxCI_1sigma    = 0.5 * ( fFluxCI_lo_1sigma + fFluxCI_up_1sigma );
        fFluxCI_lo_3sigma = fCI_lo_3sigma / iEffNorm * iN_D;
        fFluxCI_up_3sigma = fCI_up_3sigma / iEffNorm * iN_D;
        
        if( !isSignificantDataPoint() )
        {
            fFluxConstant = fUL / iEffNorm;
            fFluxConstantE = 0.;
        }
        
        if( fUL_Expected > 0 )
        {
            fFluxUL_Expected = fUL_Expected / iEffNorm * iN_D;
            fFluxUL_Expected_lo_1sigma =  fFluxUL_Expected - fUL_Expected_lo_1sigma / iEffNorm * iN_D ;
            fFluxUL_Expected_up_1sigma = -fFluxUL_Expected + fUL_Expected_up_1sigma / iEffNorm * iN_D ;
            fFluxUL_Expected_lo_2sigma =  fFluxUL_Expected - fUL_Expected_lo_2sigma / iEffNorm * iN_D ;
            fFluxUL_Expected_up_2sigma = -fFluxUL_Expected + fUL_Expected_up_2sigma / iEffNorm * iN_D ;
            
            fFluxConstantUL_Expected = fUL_Expected / iEffNorm;
            fFluxConstantUL_Expected_lo_1sigma =  fFluxConstantUL_Expected - fUL_Expected_lo_1sigma / iEffNorm ;
            fFluxConstantUL_Expected_up_1sigma = -fFluxConstantUL_Expected + fUL_Expected_up_1sigma / iEffNorm ;
            fFluxConstantUL_Expected_lo_2sigma =  fFluxConstantUL_Expected - fUL_Expected_lo_2sigma / iEffNorm ;
            fFluxConstantUL_Expected_up_2sigma = -fFluxConstantUL_Expected + fUL_Expected_up_2sigma / iEffNorm ;
        }
    }
}

/*

    check if a given time stamp is inside this time bin

    (all values in seconds)

*/
bool VFluxDataPoint::isTimeInsideRun( double iT_seconds )
{
    if( iT_seconds > 0. )
    {
        if( iT_seconds >= ( fSecondsIntoRun - 0.5 * fTimeBinDuration_sec ) && iT_seconds < ( fSecondsIntoRun + 0.5 * fTimeBinDuration_sec ) )
        {
            return true;
        }
    }
    return false;
}

/*

   setting significance parameters

   (used for the decision of when to use fluxes and when upper flux limits)

*/
void VFluxDataPoint::setSignificanceParameters( double iThresholdSignificance, double iMinEvents,
        double iUpperLimit, int iUpperlimitMethod,
        int iLiMaEqu, bool iBoundedLimits )
{
    fThresholdSignificance = iThresholdSignificance;
    fMinEvents = iMinEvents;
    fUpperLimit = iUpperLimit;
    fUpperLimitMethod = iUpperlimitMethod;
    fLiMaEqu = iLiMaEqu;
    fBoundedLimits = iBoundedLimits;
}

/*

    set spectral parameters

    these values are used to calculate a spectral weighted effective area

*/
void VFluxDataPoint::setSpectralParameters( double iMinEnergy_TeV, double E0, double alpha, double iMaxEnergy_TeV )
{
    fMinEnergy_TeV = iMinEnergy_TeV;
    fMaxEnergy_TeV = iMaxEnergy_TeV;
    fE0_TeV = E0;
    fSpectralIndex = alpha;
}

void VFluxDataPoint::calculateOrbitalPhaseData( VOrbitalPhaseData iOrbitalPhaseData )
{
    fOrbitalPhaseData = iOrbitalPhaseData;
    calculateOrbitalPhaseData();
}

void VFluxDataPoint::calculateOrbitalPhaseData()
{
    fOrbitalPhase_Start = fOrbitalPhaseData.getOrbitalPhase( fMJD_Start );
    fOrbitalPhase_Stop  = fOrbitalPhaseData.getOrbitalPhase( fMJD_Stop );
    fOrbitalPhase       = fOrbitalPhaseData.getOrbitalPhase( fMJD );
}

void VFluxDataPoint::calculateExpectedLimit( )
{
    if( fCalculateExpectedLimitsN < 1 )
    {
        return;
    }
    TH1D* temp = new TH1D( "temp" , "Counts" , 500, -50, 50 );
    TRandom3* r = new TRandom3();
    for( int i = 0; i < fCalculateExpectedLimitsN; i++ )
    {
        int Noff = r->Poisson( fNoff );
        int Non = r->Poisson( fNoff * fAlpha );
        double n = VStatistics::calcUpperLimit( Non, Noff, fAlpha, fUpperLimit, fUpperLimitMethod, fBoundedLimits );
        temp->Fill( n );
    }
    
    Double_t xq[5];  // position where to compute the quantiles in [0,1]
    Double_t yq[5];  // array to contain the quantiles
    
    xq[0] = 0.025;	// ~ -2 sigma
    xq[1] = 0.16;	// ~ -1 sigma
    xq[2] = 0.5;	// median
    xq[3] = 0.84;	// ~ +1 sigma
    xq[4] = 0.975;	// ~ +2 sigma
    
    temp->GetQuantiles( 5, yq, xq );
    
    fUL_Expected = yq[2];
    fUL_Expected_lo_1sigma = yq[1];
    fUL_Expected_up_1sigma = yq[3];
    fUL_Expected_lo_2sigma = yq[0];
    fUL_Expected_up_2sigma = yq[4];
    
    delete temp;
    
}


/*

     operator used for sorting of flux points (use time (MJD))

*/
bool VFluxDataPoint::operator<( const VFluxDataPoint& s1 ) const
{
    return fMJD < s1.fMJD;
}

/*

    add two flux points

*/
VFluxDataPoint& VFluxDataPoint::operator+=( VFluxDataPoint& s )
{

    fDebug = ( fDebug || s.fDebug );
    
    fIsCombinedFluxData = true;
    
    fCalculateExpectedLimitsN = ( fCalculateExpectedLimitsN > s.fCalculateExpectedLimitsN ? fCalculateExpectedLimitsN : s.fCalculateExpectedLimitsN );
    
    //////////////////////
    // run-wise data
    if( fRunNumber < 1 )
    {
        fRunNumber          = s.fRunNumber;
        fMJD_RunStart       = s.fMJD_RunStart;
        fMJD_RunStop        = s.fMJD_RunStop;
        fName               = s.fName;
        fDataFileName       = s.fDataFileName;
        // assume that this is always the same
        fOrbitalPhaseData   = s.fOrbitalPhaseData;
    }
    else if( fRunNumber != s.fRunNumber )
    {
        fRunNumber    = -1;
        fMJD_RunStart = 0;
        fMJD_RunStop  = 0;
        fName         = "(added data)";
        fDataFileName = "(many files)";
    }
    fRunNumber_list.insert( s.fRunNumber );
    
    //////////////////////
    // time bin data
    if( fMJD_Start > s.fMJD_Start || fMJD_Start < 10. )
    {
        fMJD_Start = s.fMJD_Start;
    }
    if( fMJD_Stop < s.fMJD_Stop )
    {
        fMJD_Stop  = s.fMJD_Stop;
    }
    // QQQQ this is not the best solution
    // need to take time mask into account
    fMJD                = 0.5 * ( fMJD_Start + fMJD_Stop );
    
    //orbital phases are messy. If we combine two runs from the same phase, but different orbits, we cannot use getOrbitalPhase with the average MJD.
    //so, we combine the available phases. Note that you have to call
    
    //also, we don't have orbital phases set if one of the flux points is 'empty'.
    if( fOrbitalPhase_Start == -99 && s.fOrbitalPhase_Start != -99 )
    {
        fOrbitalPhase_Start = s.fOrbitalPhase_Start;
        fOrbitalPhase_Stop  = s.fOrbitalPhase_Stop;
        if( fOrbitalPhase_Stop   <   fOrbitalPhase_Start )
        {
            fOrbitalPhase_Stop += 1.0;
        }
        fOrbitalPhase = 0.5 * ( fOrbitalPhase_Start + fOrbitalPhase_Stop );
        fOrbitalPhase = fOrbitalPhase - ( int )fOrbitalPhase;
    }
    else if( s.fOrbitalPhase_Start != -99 )
    {
    
        //if one of the flux points is wrapped around 0, unwrap it first.
        if( fOrbitalPhase_Stop   <   fOrbitalPhase_Start )
        {
            fOrbitalPhase_Stop += 1.0;
        }
        if( s.fOrbitalPhase_Stop < s.fOrbitalPhase_Start )
        {
            s.fOrbitalPhase_Stop += 1.0;
        }
        
        fOrbitalPhase   = 0.5 * ( fOrbitalPhase_Start +   fOrbitalPhase_Stop );
        s.fOrbitalPhase = 0.5 * ( s.fOrbitalPhase_Start + s.fOrbitalPhase_Stop );
        
        
        //if the phases of the two flux points differ by more than 0.5, we probably want to wrap one of them around 1.
        if( fabs( fOrbitalPhase - s.fOrbitalPhase ) > 0.5 )
        {
        
            //unwrap the lower flux point.
            if( fOrbitalPhase < s.fOrbitalPhase )
            {
                fOrbitalPhase_Start += 1.0;
                fOrbitalPhase_Stop += 1.0;
                fOrbitalPhase += 1.0;
            }
            else
            {
                s.fOrbitalPhase_Start += 1.0;
                s.fOrbitalPhase_Stop += 1.0;
                s.fOrbitalPhase += 1.0;
            }
            
        }
        
        //now we have made sure that the flux points are arrange in 'order'.
        
        //new start value is minimum of old start values etc.
        if( s.fOrbitalPhase_Start < fOrbitalPhase_Start )
        {
            fOrbitalPhase_Start = s.fOrbitalPhase_Start;
        }
        if( s.fOrbitalPhase_Stop > fOrbitalPhase_Stop )
        {
            fOrbitalPhase_Stop = s.fOrbitalPhase_Stop;
        }
        fOrbitalPhase = 0.5 * ( fOrbitalPhase_Start + fOrbitalPhase_Stop );
        
        //wrap phase values if necessary.
        fOrbitalPhase_Start  = fOrbitalPhase_Start - ( int )fOrbitalPhase_Start ;
        s.fOrbitalPhase_Start = s.fOrbitalPhase_Start - ( int )s.fOrbitalPhase_Start ;
        
        fOrbitalPhase_Stop  = fOrbitalPhase_Stop - ( int )fOrbitalPhase_Stop ;
        s.fOrbitalPhase_Stop = s.fOrbitalPhase_Stop - ( int )s.fOrbitalPhase_Stop ;
        
        fOrbitalPhase = fOrbitalPhase - ( int )fOrbitalPhase;
        s.fOrbitalPhase = s.fOrbitalPhase - ( int )s.fOrbitalPhase;
    }
    //else if s.fOrbitalPhase_Start == -99 just keep the old values.
    
    
    //////////////////////
    // timing
    fExposure                   += s.fExposure;
    fExposure_deadTimeCorrected += s.fExposure_deadTimeCorrected;
    fTimeBinDuration_sec        += s.fTimeBinDuration_sec;
    fTimeMask_open               = fTimeMask_open && s.fTimeMask_open;
    fDeadTime                    = 1.0 - fExposure_deadTimeCorrected / fExposure;
    //////////////////////
    // zenith etc.
    fTotZeSum        += s.fZe * s.fExposure_deadTimeCorrected;
    fTotWobbleOffset += s.fWobbleOffset * s.fExposure_deadTimeCorrected;
    fTotPedvars      += s.fPedvars * s.fExposure_deadTimeCorrected;
    fTotEffArea      += s.fEffArea_cm2 * s.fExposure_deadTimeCorrected;
    
    //////////////////////
    // event counting
    fNon  += s.fNon;
    fNoff += s.fNoff;
    
    fAlpha_U += s.fAlpha / ( 1. + s.fAlpha ) * ( s.fNon + s.fNoff );
    fAlpha_L += 1. / ( 1. + s.fAlpha ) * ( s.fNon + s.fNoff );
    if( fAlpha_L > 0. )
    {
        fAlpha = fAlpha_U / fAlpha_L;
    }
    else
    {
        fAlpha = 1.;
    }
    calculateSignificancesAndUpperLimits();
    
    // Note that calculateSignificancesAndUpperLimits() resets the effective area.
    // We fill it with an average value, weighted by the exposure.
    
    if( fExposure_deadTimeCorrected > 0. )
    {
        fZe           = fTotZeSum / fExposure_deadTimeCorrected;
        fWobbleOffset = fTotWobbleOffset / fExposure_deadTimeCorrected;
        fPedvars      = fTotPedvars / fExposure_deadTimeCorrected;
        fEffArea_cm2  = fTotEffArea / fExposure_deadTimeCorrected;
    }
    
    //////////////////////
    // flux calculation
    
    calculateFlux();
    
    fDebugCounter++;
    
    return *this;
}

VFluxDataPoint VFluxDataPoint::operator+( VFluxDataPoint& s )
{
    VFluxDataPoint v = *this;
    v += s;
    return v;
}
