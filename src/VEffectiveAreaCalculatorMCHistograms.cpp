/*! \class VEffectiveAreaCalculatorMCHistograms
    \brief filling, reading, writing of MC histograms for effective area calculation


*/

#include "VEffectiveAreaCalculatorMCHistograms.h"

VEffectiveAreaCalculatorMCHistograms::VEffectiveAreaCalculatorMCHistograms()
{
    // default name
    SetName( "MChistos" );
    
    fSpectralWeight = 0;
    
    fMCEnergyRange_TeV_min = 0.05;
    fMCEnergyRange_TeV_max = 50.0;
    fMCSpectralIndex       = 2.0;
    
    fMCCuts = false;
    fArrayxyoff_MC_min = -1.e5;
    fArrayxyoff_MC_max =  1.e5;
    
    fEnergyAxisBins_log10 = 60;
    fEnergyAxisMin_log10  = -2.;
    fEnergyAxisMax_log10  =  4.;
    
    fDebug = false;
}

void VEffectiveAreaCalculatorMCHistograms::setDefaultValues( bool b90DegIntervalls )
{
    // define azimuth bins
    fVMinAz.clear();
    fVMaxAz.clear();
    if( b90DegIntervalls )
    {
        fVMinAz.push_back( 135.0 );
        fVMaxAz.push_back( -135 );
        fVMinAz.push_back( -180.0 );
        fVMaxAz.push_back( -90. );
        for( unsigned int i = 0; i < 6; i++ )
        {
            fVMinAz.push_back( fVMinAz.back() + 45. );
            fVMaxAz.push_back( fVMaxAz.back() + 45. );
        }
    }
    else
    {
        fVMinAz.push_back( 135.0 );
        fVMaxAz.push_back( -165.0 );
        fVMinAz.push_back( 157.5 );
        fVMaxAz.push_back( -142.5 );
        fVMinAz.push_back( -180. );
        fVMaxAz.push_back( -120. );
        for( int i = 0; i < 13; i++ )
        {
            fVMinAz.push_back( fVMinAz.back() + 22.5 );
            fVMaxAz.push_back( fVMaxAz.back() + 22.5 );
        }
    }
    // (no az cut)
    fVMinAz.push_back( -1.e3 );
    fVMaxAz.push_back( +1.e3 );
    
    /////////////////////////////////////////////////////////////////
    // define  spectral index bins
    fVSpectralIndex.clear();
    for( unsigned int i = 0; i < 40; i++ )
    {
        fVSpectralIndex.push_back( 1.5 + ( double )i * 0.1 );
    }
    
    /////////////////////////////////////////////////////////////////
    fEnergyAxisBins_log10 = 60;
    fEnergyAxisMin_log10  = -2.;
    fEnergyAxisMax_log10  =  4.;
}


void VEffectiveAreaCalculatorMCHistograms::print()
{
    cout << "VEffectiveAreaCalculatorMCHistograms::print(): found ";
    cout << fVSpectralIndex.size() << " spectral index bin(s), ";
    cout << fVMinAz.size() << " azimuth bin(s)" << endl;
    if( fMCCuts )
    {
        cout << "\tMC cuts: " << fArrayxyoff_MC_min << " < MCxy < " << fArrayxyoff_MC_max << endl;
    }
    cout << "\tMC Energy range: [" << fMCEnergyRange_TeV_min << "," << fMCEnergyRange_TeV_max << "] TeV, index " << fMCSpectralIndex << endl;
    
    for( unsigned int i = 0; i < fVSpectralIndex.size(); i++ )
    {
        cout << "\tSpectral index (bin " << i << "): " << fVSpectralIndex[i] << endl;
        
        for( unsigned int j = 0; j < fVMinAz.size(); j++ )
        {
            if( getHistogram_Emc( j, i ) && getHistogram_Emc( j, i )->GetEntries() > 0 )
            {
                cout << "\tAzimuth (bin " << j << "): [" << fVMinAz[j] << ", " << fVMaxAz[j] << "]";
                cout << "\tEntries (MC): " << getHistogram_Emc( j, i )->GetEntries();
                if( getHistogram_EmcWeight( j, i ) && getHistogram_EmcWeight( j, i )->GetEntries() > 0 )
                {
                    cout << "\tEntries (MCweights): ";
                    cout << getHistogram_EmcWeight( j, i )->GetEntries();
                }
                if( getHistogram_EmcUnweighted( j ) && getHistogram_EmcUnweighted( j )->GetEntries() > 0 )
                {
                    cout << "\t Entries (Unweighted): ";
                    cout << getHistogram_EmcUnweighted( j )->GetEntries();
                }
                cout << endl;
            }
        }
    }
}

TH1D* VEffectiveAreaCalculatorMCHistograms::getHistogram_EmcUnweighted( unsigned int iAz )
{
    if( iAz < hVEmcUnWeighted.size() )
    {
        return hVEmcUnWeighted[iAz];
    }
    
    return 0;
}


TH1D* VEffectiveAreaCalculatorMCHistograms::getHistogram_Emc( unsigned int iAz, unsigned int iIndex )
{
    if( iIndex < hVEmc.size() )
    {
        if( iAz < hVEmc[iIndex].size() )
        {
            return hVEmc[iIndex][iAz];
        }
    }
    
    return 0;
}

TProfile* VEffectiveAreaCalculatorMCHistograms::getHistogram_EmcWeight( unsigned int iAz, unsigned int iIndex )
{
    if( iIndex < hVEmcSWeight.size() )
    {
        if( iAz < hVEmcSWeight[iIndex].size() )
        {
            return hVEmcSWeight[iIndex][iAz];
        }
    }
    
    return 0;
}

bool VEffectiveAreaCalculatorMCHistograms::readFromEffectiveAreaTree( string iFile )
{
    iFile = "nofile";
    return false;
}

bool VEffectiveAreaCalculatorMCHistograms::readFromEffectiveAreaFile( string iFile )
{
    iFile = "nofile";
    return false;
}

bool VEffectiveAreaCalculatorMCHistograms::fill( double i_ze, TTree* i_MCData, bool bAzimuthBins )
{
    cout << endl;
    cout << "filling MC histograms for effective area calculation for ze " << i_ze << " [deg]" << endl;
    cout << "=========================================================================================" << endl;
    if( fVSpectralWeight.size() > 0 && fVSpectralWeight[0] )
    {
        fVSpectralWeight[0]->print();
    }
    cout << "=========================================================================================" << endl;
    cout << endl;
    
    if( !i_MCData )
    {
        cout << "VEffectiveAreaCalculatorMCHistograms::fill error: no MC data chain" << endl;
        return false;
    }
    
    float i_fMCE0 = 0.;
    float i_fMCAz = 0.;
    float i_fMCxoff = 0.;
    float i_fMCyoff = 0.;
    i_MCData->SetBranchAddress( "MCe0", &i_fMCE0 );
    if( bAzimuthBins )
    {
        i_MCData->SetBranchAddress( "MCaz", &i_fMCAz );
    }
    if( fMCCuts )
    {
        i_MCData->SetBranchAddress( "MCxoff", &i_fMCxoff );
        i_MCData->SetBranchAddress( "MCyoff", &i_fMCyoff );
    }
    
    // spectral weight
    double i_weight = 1.;
    // MC energy (log10)
    double eMC = 0.;
    
    // array lengths
    unsigned int i_vMinAzSize        = fVMinAz.size();
    unsigned int i_vSpectralIndexSize = fVSpectralWeight.size();
    cout << "\t array lengths az: " << i_vMinAzSize << ", spectral index: " << i_vSpectralIndexSize << endl;
    
    // entries in MC tree (must be long, chain could contain lots of events)
    Long64_t nentries = i_MCData->GetEntries();
    cout << "total number of MC events: " << nentries << endl;
    if( fMCCuts )
    {
        cout << "(apply MC cuts)" << endl;
    }
    //////////////////////////////////////////////////////////
    // now loop over all MC entries
    for( Long64_t i = 0; i < nentries; i++ )
    {
        // read data (subset)
        i_MCData->GetEntry( i );
        
        // apply MC cuts
        if( fMCCuts )
        {
            if( i_fMCxoff * i_fMCxoff + i_fMCyoff * i_fMCyoff > fArrayxyoff_MC_max )
            {
                continue;
            }
            if( i_fMCxoff * i_fMCxoff + i_fMCyoff * i_fMCyoff < fArrayxyoff_MC_min )
            {
                continue;
            }
        }
        
        // log of MC energy
        eMC = log10( i_fMCE0 );
        
        // fill the MC histogram for all az bins
        for( unsigned int i_az = 0; i_az < i_vMinAzSize; i_az++ )
        {
            // check which azimuth bin we are
            if( bAzimuthBins && i_ze > 3. )
            {
                // confine MC az to -180., 180.
                if( i_fMCAz > 180. )
                {
                    i_fMCAz -= 360.;
                }
                // expect bin like [135,-135]
                if( fVMinAz[i_az] > fVMaxAz[i_az] )
                {
                    if( i_fMCAz < fVMinAz[i_az] && i_fMCAz > fVMaxAz[i_az] )
                    {
                        continue;
                    }
                }
                // expect bin like [-135,-45.]
                else
                {
                    if( i_fMCAz < fVMinAz[i_az] || i_fMCAz > fVMaxAz[i_az] )
                    {
                        continue;
                    }
                }
            }
            // loop over all spectral index
            for( unsigned int s = 0; s < i_vSpectralIndexSize; s++ )
            {
                // weight by spectral index
                i_weight = fVSpectralWeight[s]->getSpectralWeight( i_fMCE0 );
                
                // fill MC histograms
                if( hVEmc[s][i_az] )
                {
                    hVEmc[s][i_az]->Fill( eMC, i_weight );
                }
                if( hVEmcSWeight[s][i_az] )
                {
                    hVEmcSWeight[s][i_az]->Fill( eMC, i_weight );
                }
            }
            // fill unweighted histogram
            hVEmcUnWeighted[i_az]->Fill( eMC );
        }
    } // end of loop over all MC entries
    
    return true;
}

void VEffectiveAreaCalculatorMCHistograms::initializeHistograms()
{
    initializeHistograms( fVMinAz, fVMaxAz, fVSpectralIndex, fEnergyAxisBins_log10, fEnergyAxisMin_log10, fEnergyAxisMax_log10 );
}


void VEffectiveAreaCalculatorMCHistograms::initializeHistograms( vector< double > iAzMin, vector< double > iAzMax,
        vector< double > iSpectralIndex,
        int nbins, double xmin, double xmax )
{
    fVMinAz = iAzMin;
    fVMaxAz = iAzMax;
    fVSpectralIndex = iSpectralIndex;
    
    char hname[400];
    
    vector< TH1D* > iT_TH1D;
    vector< TProfile* > iT_TProfile;
    
    for( unsigned int i = 0; i < fVSpectralIndex.size(); i++ )
    {
        iT_TProfile.clear();
        
        // histograms for effective area calculation
        iT_TH1D.clear();
        for( unsigned int j = 0; j < fVMinAz.size(); j++ )
        {
            sprintf( hname, "hVVEmc_%u_%u", i, j );
            iT_TH1D.push_back( new TH1D( hname, "", nbins, xmin, xmax ) );
            iT_TH1D.back()->SetXTitle( "energy_{MC} [TeV]" );
            iT_TH1D.back()->SetYTitle( "entries (weighted)" );
            iT_TH1D.back()->Sumw2();
        }
        hVEmc.push_back( iT_TH1D );
        
        for( unsigned int j = 0; j < fVMinAz.size(); j++ )
        {
            sprintf( hname, "hVVEmcSWeight_%u_%u", i, j );
            iT_TProfile.push_back( new TProfile( hname, "", nbins, xmin, xmax, 0., 1.e12 ) );
            iT_TProfile.back()->SetXTitle( "energy_{MC} [TeV]" );
            iT_TProfile.back()->SetYTitle( "spectral weight" );
        }
        hVEmcSWeight.push_back( iT_TProfile );
        
    }
    // unweighted histogram (for debugging purposes)
    for( unsigned int j = 0; j < fVMinAz.size(); j++ )
    {
        sprintf( hname, "hVVEmc_%u", j );
        hVEmcUnWeighted.push_back( new TH1D( hname, "", nbins, xmin, xmax ) );
        hVEmcUnWeighted.back()->SetXTitle( "energy_{MC} [TeV]" );
        hVEmcUnWeighted.back()->SetYTitle( "entries (not weighted)" );
        hVEmcUnWeighted.back()->Sumw2();
    }
    
    
    // set spectral weight vector
    for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
    {
        fVSpectralWeight.push_back( new VSpectralWeight() );
        // weight by spectral index
        fVSpectralWeight.back()->setSpectralIndex( fVSpectralIndex[s] );
        fVSpectralWeight.back()->setMCParameter( fMCSpectralIndex, fMCEnergyRange_TeV_min, fMCEnergyRange_TeV_max );
    }
    // backwards compatibility
    if( fVSpectralWeight.size() > 0 && fVSpectralWeight[0] )
    {
        fSpectralWeight = fVSpectralWeight[0];
    }
}

bool VEffectiveAreaCalculatorMCHistograms::setMonteCarloEnergyRange( double iMin, double iMax, double iMCIndex )
{
    fMCEnergyRange_TeV_min = iMin;
    fMCEnergyRange_TeV_max = iMax;
    fMCSpectralIndex       = iMCIndex;
    
    for( unsigned int i = 0; i < fVSpectralWeight.size(); i++ )
    {
        if( fVSpectralWeight[i] )
        {
            fVSpectralWeight[i]->setMCParameter( fMCSpectralIndex, fMCEnergyRange_TeV_min, fMCEnergyRange_TeV_max );
        }
    }
    // backwards compatibility
    if( fVSpectralWeight.size() > 0 && fVSpectralWeight[0] )
    {
        fSpectralWeight = fVSpectralWeight[0];
    }
    
    return true;
}

void VEffectiveAreaCalculatorMCHistograms::setCuts( double iArrayxyoff_MC_min, double iArrayxyoff_MC_max )
{
    fArrayxyoff_MC_min = iArrayxyoff_MC_min;
    fArrayxyoff_MC_max = iArrayxyoff_MC_max;
    if( fArrayxyoff_MC_min < 0. || fArrayxyoff_MC_max < 0. )
    {
        fMCCuts = false;
    }
    else
    {
        fMCCuts = true;
    }
}

/*

   add histograms

*/
bool VEffectiveAreaCalculatorMCHistograms::add( const VEffectiveAreaCalculatorMCHistograms* iMChis )
{
    if( !iMChis )
    {
        return false;
    }
    
    // make sure that both instants are similar
    int i_check = checkParameters( iMChis );
    if( i_check != 0 )
    {
        cout << "VEffectiveAreaCalculatorMCHistograms::add() error: parameters differ (" << i_check << ")" << endl;
        return false;
    }
    
    // add histograms
    for( unsigned int i = 0; i < hVEmc.size(); i++ )
    {
        for( unsigned int j = 0; j < hVEmc[i].size(); j++ )
        {
            if( hVEmc[i][j] && iMChis->hVEmc[i][j] )
            {
                hVEmc[i][j]->Add( iMChis->hVEmc[i][j] );
            }
        }
    }
    
    for( unsigned int i = 0; i < hVEmcSWeight.size(); i++ )
    {
        for( unsigned int j = 0; j < hVEmcSWeight[i].size(); j++ )
        {
            if( hVEmcSWeight[i][j] && iMChis->hVEmcSWeight[i][j] )
            {
                hVEmcSWeight[i][j]->Add( iMChis->hVEmcSWeight[i][j] );
            }
        }
    }
    
    for( unsigned int j = 0; j < hVEmcUnWeighted.size(); j++ )
    {
        if( hVEmcUnWeighted[j] && iMChis->hVEmcUnWeighted[j] )
        {
            hVEmcUnWeighted[j]->Add( iMChis->hVEmcUnWeighted[j] );
        }
    }
    
    return true;
}

/*

     note the interesting error coding...

*/
int VEffectiveAreaCalculatorMCHistograms::checkParameters( const VEffectiveAreaCalculatorMCHistograms* iMChis )
{
    if( fDebug )
    {
        cout << "VEffectiveAreaCalculatorMCHistograms::checkParameters" << endl;
    }
    if( !iMChis )
    {
        return 1;
    }
    
    if( iMChis->fMCCuts != fMCCuts )
    {
        return 2;
    }
    if( TMath::Abs( iMChis->fArrayxyoff_MC_min - fArrayxyoff_MC_min ) > 1.e-3 )
    {
        return 3;
    }
    if( TMath::Abs( iMChis->fArrayxyoff_MC_max - fArrayxyoff_MC_max ) > 1.e-3 )
    {
        return 4;
    }
    if( TMath::Abs( iMChis->fMCEnergyRange_TeV_min - fMCEnergyRange_TeV_min ) > 1.e-3 )
    {
        return 5;
    }
    if( TMath::Abs( iMChis->fMCEnergyRange_TeV_max - fMCEnergyRange_TeV_max ) > 1.e-3 )
    {
        return 6;
    }
    if( TMath::Abs( iMChis->fMCSpectralIndex - fMCSpectralIndex ) > 1.e-3 )
    {
        return 7;
    }
    if( iMChis->fVMinAz.size() != fVMinAz.size() )
    {
        return 8;
    }
    for( unsigned int i = 0; i < iMChis->fVMinAz.size(); i++ ) if( TMath::Abs( iMChis->fVMinAz[i] - fVMinAz[i] ) > 1.e-3 )
        {
            return 9;
        }
    if( iMChis->fVMaxAz.size() != fVMaxAz.size() )
    {
        return 10;
    }
    for( unsigned int i = 0; i < iMChis->fVMaxAz.size(); i++ ) if( TMath::Abs( iMChis->fVMaxAz[i] - fVMaxAz[i] ) > 1.e-3 )
        {
            return 11;
        }
    if( iMChis->fVSpectralIndex.size() != fVSpectralIndex.size() )
    {
        return 12;
    }
    for( unsigned int i = 0; i < iMChis->fVSpectralIndex.size(); i++ )
    {
        if( TMath::Abs( iMChis->fVSpectralIndex[i] - fVSpectralIndex[i] ) > 1.e-3 )
        {
            return 13;
        }
    }
    
    if( iMChis->hVEmc.size() != hVEmc.size() )
    {
        return 14;
    }
    for( unsigned int i = 0; i < hVEmc.size(); i++ )
    {
        if( hVEmc[i].size() != iMChis->hVEmc[i].size() )
        {
            return 15;
        }
        for( unsigned int j = 0; j < hVEmc[i].size(); j++ )
        {
            if( hVEmc[i][j] && iMChis->hVEmc[i][j] )
            {
                if( hVEmc[i][j]->GetNbinsX() != iMChis->hVEmc[i][j]->GetNbinsX() )
                {
                    return 16;
                }
                if( hVEmc[i][j]->GetXaxis()->GetXmin() != iMChis->hVEmc[i][j]->GetXaxis()->GetXmin() )
                {
                    return 17;
                }
                if( hVEmc[i][j]->GetXaxis()->GetXmax() != iMChis->hVEmc[i][j]->GetXaxis()->GetXmax() )
                {
                    return 18;
                }
            }
        }
    }
    
    if( iMChis->hVEmcSWeight.size() != hVEmcSWeight.size() )
    {
        return 19;
    }
    for( unsigned int i = 0; i < hVEmcSWeight.size(); i++ )
    {
        if( hVEmcSWeight[i].size() != iMChis->hVEmcSWeight[i].size() )
        {
            return 20;
        }
        for( unsigned int j = 0; j < hVEmcSWeight[i].size(); j++ )
        {
            if( hVEmcSWeight[i][j] && iMChis->hVEmcSWeight[i][j] )
            {
                if( hVEmcSWeight[i][j]->GetNbinsX() != iMChis->hVEmcSWeight[i][j]->GetNbinsX() )
                {
                    return 21;
                }
                if( hVEmcSWeight[i][j]->GetXaxis()->GetXmin() != iMChis->hVEmcSWeight[i][j]->GetXaxis()->GetXmin() )
                {
                    return 22;
                }
                if( hVEmcSWeight[i][j]->GetXaxis()->GetXmax() != iMChis->hVEmcSWeight[i][j]->GetXaxis()->GetXmax() )
                {
                    return 23;
                }
            }
        }
    }
    
    return 0;
}

bool VEffectiveAreaCalculatorMCHistograms::matchDataVectors( vector< double > iAzMin, vector< double > iAzMax, vector< double > iSpectralIndex )
{
    vector< double > iVSpectralIndex_new;
    vector< vector< TH1D* > > ihVEmc_new;
    vector< vector< TProfile* > > iVEmcSWeight_new;
    vector< TH1D* > iH;
    vector< TProfile* > iP;
    
    // match spectral index
    if( fDebug )
    {
        cout << "VEffectiveAreaCalculatorMCHistograms::matchDataVectors: ";
        cout << "matching spectral index;";
        cout << " found: " << fVSpectralIndex.size();
        cout << " requested: " << iSpectralIndex.size() << endl;
    }
    for( unsigned int i = 0; i < fVSpectralIndex.size(); i++ )
    {
        for( unsigned int j = 0; j < iSpectralIndex.size(); j++ )
        {
            if( TMath::Abs( iSpectralIndex[j] - fVSpectralIndex[i] ) < 1.e-3 )
            {
                iVSpectralIndex_new.push_back( fVSpectralIndex[i] );
                ihVEmc_new.push_back( hVEmc[i] );
                iVEmcSWeight_new.push_back( hVEmcSWeight[i] );
            }
        }
    }
    fVSpectralIndex = iVSpectralIndex_new;
    hVEmc = ihVEmc_new;
    hVEmcSWeight = iVEmcSWeight_new;
    
    if( iAzMin.size() == 0 )
    {
        iAzMin = fVMinAz;
    }
    if( iAzMax.size() == 0 )
    {
        iAzMax = fVMaxAz;
    }
    
    if( iAzMin.size() != iAzMax.size() )
    {
        cout << "VEffectiveAreaCalculatorMCHistograms::matchDataVectors error: mismatch in az vector size: ";
        cout << iAzMin.size() << "\t" << iAzMax.size() << endl;
        return false;
    }
    // match azimuth vector
    vector< double > iVMinAz_new;
    vector< double > iVMaxAz_new;
    vector< unsigned int > iVAz_match;
    for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
    {
        iVMinAz_new.clear();
        iVMaxAz_new.clear();
        vector< double > iVMinAz = fVMinAz;
        vector< double > iVMaxAz = fVMaxAz;
        vector< TH1D* > ihVEmc_new;
        vector< TProfile* > iVEmcSWeight_new;
        
        for( unsigned int i = 0; i < iAzMin.size(); i++ )
        {
            for( unsigned int j = 0; j < iVMinAz.size(); j++ )
            {
                if( TMath::Abs( iAzMin[i] - iVMinAz[j] ) < 1.e-3 && TMath::Abs( iAzMax[i] - iVMaxAz[j] ) < 1.e-3 )
                {
                    iVMinAz_new.push_back( iAzMin[i] );
                    iVMaxAz_new.push_back( iAzMax[i] );
                    ihVEmc_new.push_back( hVEmc[s][j] );
                    iVEmcSWeight_new.push_back( hVEmcSWeight[s][j] );
                    iVAz_match.push_back( j );
                }
            }
        }
        hVEmc[s] = ihVEmc_new;
        hVEmcSWeight[s] = iVEmcSWeight_new;
    }
    fVMinAz = iVMinAz_new;
    fVMaxAz = iVMaxAz_new;
    
    // unweighted histogram
    vector< TH1D* > ihVEmcUnWeighted_new;
    for( unsigned int i = 0; i < iVAz_match.size(); i++ )
    {
        if( iVAz_match[i] < hVEmcUnWeighted.size() )
        {
            ihVEmcUnWeighted_new.push_back( hVEmcUnWeighted[iVAz_match[i]] );
        }
    }
    hVEmcUnWeighted = ihVEmcUnWeighted_new;
    
    return true;
}
