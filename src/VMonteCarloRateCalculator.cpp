/*! \class VMonteCarloRateCalculator
    \brief calculate rates from Monte Carlo

*/

#include "VMonteCarloRateCalculator.h"

VMonteCarloRateCalculator::VMonteCarloRateCalculator()
{
    default_settings();
    setPlottingStyle();
}


VMonteCarloRateCalculator::VMonteCarloRateCalculator( string ifile )
{
    fFile = new TFile( ifile.c_str() );
    if( fFile->IsZombie() )
    {
        cout << "VMonteCarloRateCalculator: error opening file " << ifile << endl;
        return;
    }
    fze = 0.;
    faz = 0;
    fwoff = 0.;
    fnoise = 0;
    fnrates = 0;
    fMCTree = ( TTree* )fFile->Get( "fMCRate" );
    if( !fMCTree )
    {
        cout << "VMonteCarloRateCalculator: no rate tree found" << endl;
        return;
    }
    fMCTree->SetBranchAddress( "ze", &fze );
    fMCTree->SetBranchAddress( "az", &faz );
    fMCTree->SetBranchAddress( "Woff", &fwoff );
    fMCTree->SetBranchAddress( "noise", &fnoise );
    fMCTree->SetBranchAddress( "nrates", &fnrates );
    fMCTree->SetBranchAddress( "MCrate", fMCrate );
    
}


/*
    energy in log10

    expect power law
*/
double VMonteCarloRateCalculator::getMonteCarloRate( int nbins, double* e, double* eff,
        double i_gamma, double i_phi,
        double iE0, double iEMin, double iEMax, double bDebug )
{
    if( nbins <= 2 || e == 0 || eff == 0 )
    {
        return 0.;
    }
    
    fenergy.clear();
    feffectiveArea.clear();
    for( int i = 0; i < nbins; i++ )
    {
        fenergy.push_back( e[i] );
        feffectiveArea.push_back( eff[i] );
    }
    return getMonteCarloRate( fenergy, feffectiveArea, i_gamma, i_phi, iE0, iEMin, iEMax, bDebug );
}


/*
   energy vector in log10
*/
double VMonteCarloRateCalculator::getMonteCarloRate( vector< double >& e, vector< double >& eff,
        double i_gamma, double i_phi,
        double iE0, double iEMin, double iEMax, bool bDebug )
{
    iEMin = log10( iEMin );
    iEMax = log10( iEMax );
    
    i_phi *= TMath::Power( iE0, -1.*i_gamma );
    
    // dN
    double y1 = 0.;
    // rate
    double iTot = 0.;
    double x1 = 0.;
    double x2 = 0.;
    double i_bflux = 0.;
    
    for( unsigned int i = 0; i < e.size(); i++ )
    {
        // get lower and upper energy bin
        if( i == 0 )
        {
            x1 = e[i];
            x2 = e[i] + ( e[i + 1] - e[i] ) / 2.;
        }
        else if( i == e.size() - 1 )
        {
            x1 = e[i] - ( e[i] - e[i - 1] ) / 2.;
            x2 = e[i];
        }
        else
        {
            x1 = e[i] - ( e[i] - e[i - 1] ) / 2.;
            x2 = e[i] + ( e[i + 1] - e[i] ) / 2.;
        }
        // effective area (m2 -> cm2) (fluxes are given in cm^-2)
        y1 = eff[i] * 1.e4;
        
        // calculate number of events per energy bin
        // (energy (x) is in log E [TeV]
        
        // (this is the integral flux above energy of upper bin edge)
        i_bflux  = i_phi / ( i_gamma + 1 ) * pow( pow( 10., x2 ), i_gamma + 1 );
        // (subtract integral flux above energy of lower bin edge)
        i_bflux  -= i_phi / ( i_gamma + 1 ) * pow( pow( 10., x1 ), i_gamma + 1 );
        
        // multiply flux by effective areas
        // (this is the number of gammas per energy bin (dN not dN/dE)
        y1  *= i_bflux;
        
        // this is approximate, iEMin and iEMax can be much smaller than interval x1, x2
        if( x1 > iEMin && x2 < iEMax )
        {
            iTot += y1;
        }
    }
    if( bDebug )
    {
        cout << "VMonteCarloRateCalculator::getMonteCarloRate " << e.size() << "\t" << eff.size() << endl;
        cout << "\t Gamma " << i_gamma << " phi: " << i_phi << " E0: " << iE0;
        cout << " EMin: " << iEMin << " EMax: " << iEMax;
        cout << "Rate [1/s]: " << iTot << " Rate [1/min]: " << iTot * 60. << endl;
    }
    if( e.size() != eff.size() || e.size() <= 2 )
    {
        return 0.;
    }
    
    // convert rate [1/s] to [1/min]
    return iTot * 60.;
}


TCanvas* VMonteCarloRateCalculator::plot_MonteCarloRate_vs_wobbleOffsets( TCanvas* c, double ze, int az, int noise, string iPlottingOption )
{
    if( !c )
    {
        c = new TCanvas( "cRatevsWobbleDirection", "cRatevsWobbleDirection", 10, 10, 400, 400 );
        c->Draw();
    }
    else
    {
        c->cd();
    }
    
    TGraphAsymmErrors* g = getMonteCarloRate_vs_wobbleOffsets( ze, az, noise );
    if( g )
    {
        g->Print();
    }
    if( g )
    {
        g->Draw( iPlottingOption.c_str() );
    }
    
    return c;
}


TGraphAsymmErrors* VMonteCarloRateCalculator::getMonteCarloRate_vs_wobbleOffsets( double ze, int az, int noise )
{
    if( !fFile || !fMCTree )
    {
        return 0;
    }
    
    TGraphAsymmErrors* g = new TGraphAsymmErrors( 1 );
    setGraphPlottingStyle( g );
    
    double i_mean = 0.;
    double i_min = 0.;
    double i_max = 0.;
    
    int z = 0;
    for( int i = 0; i < fMCTree->GetEntries(); i++ )
    {
        fMCTree->GetEntry( i );
        
        if( faz != az )
        {
            continue;
        }
        if( fnoise != noise )
        {
            continue;
        }
        if( TMath::Abs( fze - ze ) > 0.7 )
        {
            continue;
        }
        
        getMinMaxRates( fnrates, fMCrate, i_mean, i_min, i_max );
        
        g->SetPoint( z, fwoff, i_mean );
        g->SetPointEYhigh( z, i_max - i_mean );
        g->SetPointEYlow( z, i_mean - i_min );
        
        z++;
    }
    
    return g;
}


void VMonteCarloRateCalculator::getMinMaxRates( unsigned int n, double* r, double& i_mean, double& i_min, double& i_max )
{
    i_mean = 0.;
    double z = 0;
    i_min = 1.e10;
    i_max = -1.e10;
    
    for( unsigned int i = 0; i < n; i++ )
    {
        if( TMath::Abs( r[i] < 1.e-2 ) )
        {
            continue;
        }
        if( r[i] < i_min )
        {
            i_min = r[i];
        }
        if( r[i] > i_max )
        {
            i_max = r[i];
        }
        i_mean += r[i];
        z++;
    }
    if( z > 0. )
    {
        i_mean /= z;
    }
}

/*
    energy in log10

    expect power law
*/
double VMonteCarloRateCalculator::getMonteCarloRate( int nbins, double* e, double* eff,
        VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
        double iEMin, double iEMax, bool bDebug )
{
    if( nbins <= 2 || e == 0 || eff == 0 )
    {
        return 0.;
    }
    
    vector< double > energy;
    vector< double > effectiveArea;
    for( int i = 0; i < nbins; i++ )
    {
        energy.push_back( e[i] );
        effectiveArea.push_back( eff[i] );
    }
    return getMonteCarloRate( energy, effectiveArea, e_lit, e_lit_ID, 0, energy.size(), iEMin, iEMax, 0, bDebug );
}


/*

   keep this function for backward compatibility

   energy in log10
*/
double VMonteCarloRateCalculator::getMonteCarloRate( vector< double > e, vector< double > eff,
        VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
        unsigned int iEMinBin, unsigned int iEMaxBin, TH2* iResponseMatrix, bool bDebug )
{
    return getMonteCarloRate( e, eff, e_lit, e_lit_ID, iEMinBin, iEMaxBin, 0., 0., iResponseMatrix, bDebug );
}

/*
   energy vector e in log10 [TeV]

   iEMin and iEMax on linear axis [TeV]

   (this function is used in VSensitivityCalculator)
*/
double VMonteCarloRateCalculator::getMonteCarloRate( vector< double > e, vector< double > eff,
        VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
        unsigned int iEMinBin, unsigned int iEMaxBin,
        double iEMin, double iEMax, TH2* iResponseMatrix, bool bDebug )
{
    vector< double > e_gamma;
    
    return getMonteCarloRate( e, eff, e_lit, e_lit_ID, iEMinBin, iEMaxBin, iEMin, iEMax, iResponseMatrix, e_gamma, bDebug );
}

/*

*/
double VMonteCarloRateCalculator::getMonteCarloRate( vector< double > e, vector< double > eff,
        VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
        unsigned int iEMinBin, unsigned int iEMaxBin,
        double iEMin, double iEMax, TH2* iResponseMatrix,
        vector< double > e_gamma, bool bDebug )
{
    if( e.size() == 0 )
    {
        return -99;
    }
    if( e_gamma.size() == 0 )
    {
        e_gamma = e;
    }
    // check bin region
    if( iEMinBin > e_gamma.size() )
    {
        return 0.;
    }
    // adjust maximum bin
    if( iEMaxBin > e_gamma.size() - 1 )
    {
        iEMaxBin = e_gamma.size() - 1;
    }
    
    if( !iResponseMatrix && e_gamma.size() != eff.size() )
    {
        cout <<  "VMonteCarloRateCalculator::getMonteCarloRate error: energy and effective area vector have different length: ";
        cout << e_gamma.size() << "\t" << eff.size() << endl;
        return -99.;
    }
    
    // check function ID of cosmic ray / gamma spectrum
    if( !e_lit || !e_lit->isValidID( e_lit_ID ) )
    {
        return -99.;
    }
    
    // dN
    double y1 = 0.;
    // rate
    double iTot = 0.;
    double x1 = 0.;
    double x2 = 0.;
    double i_bflux = 0.;
    int i_Rec_Bin = 0;
    int j_Rec_Bin = 0;
    
    if( bDebug )
    {
        cout << "VMonteCarloRateCalculator::getMonteCarloRate " << e_lit_ID << ": " << e.size() << "\t" << eff.size();
        cout << " (gamma: " << e_gamma.size() << ")" << endl;
        cout << "\t";
        cout << " EMin: " << e_gamma[iEMinBin] << " (" << iEMinBin << ")";
        cout << " EMax: " << e_gamma[iEMaxBin] << " (" << iEMaxBin << ")";
        cout << endl;
    }
    
    // loop over energy interval in units of effective area vector
    for( unsigned int i = iEMinBin; i <= iEMaxBin; i++ )
    {
        // get energy bin min/max
        if( i == 0 )
        {
            x1 = e_gamma[i];
            x2 = e_gamma[i] + ( e_gamma[i + 1] - e_gamma[i] ) / 2.;
        }
        else if( i == e.size() - 1 )
        {
            x1 = e_gamma[i] - ( e_gamma[i] - e_gamma[i - 1] ) / 2.;
            x2 = e_gamma[i];
        }
        else
        {
            x1 = e_gamma[i] - ( e_gamma[i] - e_gamma[i - 1] ) / 2.;
            x2 = e_gamma[i] + ( e_gamma[i + 1] - e_gamma[i] ) / 2.;
        }
        if( bDebug )
        {
            cout << "\t\t integrate [" << x1 << ", " << x2 << "]" << endl;
        }
        // =================================================================
        // apply response matrix
        // (important: expect same binning in response matrix and effective area vector!)
        if( iResponseMatrix )
        {
            i_Rec_Bin = iResponseMatrix->GetXaxis()->FindBin( e_gamma[i] );
            if( bDebug )
            {
                cout << "\t\t REC BIN " << i_Rec_Bin << "\t" << e_gamma[i] << endl;
            }
            for( unsigned int j = 0; j < e.size(); j++ )
            {
                j_Rec_Bin = iResponseMatrix->GetYaxis()->FindBin( e[j] );
                if( iResponseMatrix->GetBinContent( i_Rec_Bin, j_Rec_Bin ) > 1.e-7 )
                {
                    // effective area (m2 -> cm2) * response matrix element
                    y1 = eff[j] * 1.e4 * iResponseMatrix->GetBinContent( i_Rec_Bin, j_Rec_Bin );
                    // get integral flux for this bin
                    i_bflux = e_lit->getIntegralFlux( pow( 10., iResponseMatrix->GetYaxis()->GetBinLowEdge( j_Rec_Bin ) ),
                                                      pow( 10., iResponseMatrix->GetYaxis()->GetBinUpEdge( j_Rec_Bin ) ), e_lit_ID );
                    // multiply flux by effective areas
                    // (this is the number of gammas per energy bin (dN not dN/dE))
                    y1 *= i_bflux;
                    // total flux
                    iTot += y1;
                }
            }
        }
        // =================================================================
        else
        {
            // effective area (m2 -> cm2)
            y1 = eff[i] * 1.e4;
            
            // calculate number of events per energy bin
            // (energy (x) is in log E [TeV]
            
            // get integral flux for this bin
            i_bflux = e_lit->getIntegralFlux( pow( 10., x1 ), pow( 10., x2 ), e_lit_ID );
            
            // multiply flux by effective areas
            // (this is the number of gammas per energy bin (dN not dN/dE)
            y1  *= i_bflux;
            
            // this is approximate, iEMin and iEMax can be much smaller than interval x1, x2
            if( iEMin > 0. && iEMax > 0. )
            {
                // total flux
                if( x1 - log10( iEMin ) > -1.e-3 && log10( iEMax ) - x2 > -1.e3 )
                {
                    iTot += y1;
                }
            }
            else
            {
                iTot += y1;
            }
        }
    }
    if( bDebug )
    {
        cout << "\tRate [1/s]: " << iTot << " Rate [1/min]: " << iTot * 60. << endl;
    }
    
    // convert rate [1/s] to [1/min]
    return iTot * 60.;
}
