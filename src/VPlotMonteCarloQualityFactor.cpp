/*! \file VPlotMonteCarloQualityFactor
    \brief fill and plot quality factors for MC


*/

#include "VPlotMonteCarloQualityFactor.h"


VPlotMonteCarloQualityFactor::VPlotMonteCarloQualityFactor()
{
    fDebug = true;
    
    fSignalChain = 0;
    fBackgroundChain = 0;
    
    setMSCCuts();
    setEnergyRange();
    
    initializeHistograms();
}

void VPlotMonteCarloQualityFactor::setMSCCuts( double iMSCW_min, double iMSCW_max, double iMSCL_min, double iMSCL_max )
{
    fMSCW_min = iMSCW_min;
    fMSCW_max = iMSCW_max;
    fMSCL_min = iMSCL_min;
    fMSCL_max = iMSCL_max;
}

void VPlotMonteCarloQualityFactor::initializeHistograms()
{
    // number of bins
    int iBinning = 100;
    
    // parameters for different variables
    
    // MSCW
    fData["MSCW"] = new VPlotMonteCarloQualityFactorData();
    fData["MSCW"]->fVar_min = -2.;
    fData["MSCW"]->fVar_max =  5.;
    // MSCL
    fData["MSCL"] = new VPlotMonteCarloQualityFactorData();
    fData["MSCL"]->fVar_min = -2.;
    fData["MSCL"]->fVar_max =  5.;
    // dE
    fData["dE"] = new VPlotMonteCarloQualityFactorData();
    fData["dE"]->fVar_min =  0.;
    fData["dE"]->fVar_max =  0.5;
    // EChi2  (actually log10(EChi2)
    fData["EChi2"] = new VPlotMonteCarloQualityFactorData();
    fData["EChi2"]->fVar_min =  -3.;
    fData["EChi2"]->fVar_max =  5.;
    // EmissionHeight
    fData["EmissionHeight"] = new VPlotMonteCarloQualityFactorData();
    fData["EmissionHeight"]->fVar_min =  0.;
    fData["EmissionHeight"]->fVar_max =  50.;
    // EmissionHeightChi2
    fData["EmissionHeightChi2"] = new VPlotMonteCarloQualityFactorData();
    fData["EmissionHeightChi2"]->fVar_min =  -200.;
    fData["EmissionHeightChi2"]->fVar_max =   200.;
    
    // create histograms
    char hname[800];
    //   char htitle[800];
    map< string, VPlotMonteCarloQualityFactorData* >::iterator iData;
    for( iData = fData.begin(); iData != fData.end(); iData++ )
    {
        sprintf( hname, "hSi_%s", ( *iData ).first.c_str() );
        //      sprintf( htitle, "%s (signal)", (*iData).first.c_str() );
        ( *iData ).second->hSignal = new TH1D( hname, "", iBinning, ( *iData ).second->fVar_min, ( *iData ).second->fVar_max );
        ( *iData ).second->hSignal->SetXTitle( ( *iData ).first.c_str() );
        ( *iData ).second->hSignal->SetYTitle( "normalized number of events" );
        setHistogramPlottingStyle( ( *iData ).second->hSignal, 1, 2, 1., 1, 1, 3003 );
        
        sprintf( hname, "hBc_%s", ( *iData ).first.c_str() );
        //      sprintf( htitle, "%s (background)", (*iData).first.c_str() );
        ( *iData ).second->hBackground = new TH1D( hname, "", iBinning, ( *iData ).second->fVar_min, ( *iData ).second->fVar_max );
        ( *iData ).second->hBackground->SetXTitle( ( *iData ).first.c_str() );
        ( *iData ).second->hBackground->SetYTitle( "normalized number of events" );
        setHistogramPlottingStyle( ( *iData ).second->hBackground, 2, 2, 1., 1, 1, 3003 );
        
        sprintf( hname, "hQl_%s", ( *iData ).first.c_str() );
        //      sprintf( htitle, "%s (q-factor, lower bound)", (*iData).first.c_str() );
        ( *iData ).second->hQFactors_LowerCut = new TH1D( hname, "", iBinning, ( *iData ).second->fVar_min, ( *iData ).second->fVar_max );
        ( *iData ).second->hQFactors_LowerCut->SetXTitle( ( *iData ).first.c_str() );
        ( *iData ).second->hQFactors_LowerCut->SetYTitle( "q-factor" );
        setHistogramPlottingStyle( ( *iData ).second->hQFactors_LowerCut, 1, 2, 1., 1, 1, 3003 );
        
        sprintf( hname, "hQu_%s", ( *iData ).first.c_str() );
        //      sprintf( htitle, "%s (q-factor, upper bound)", (*iData).first.c_str() );
        ( *iData ).second->hQFactors_UpperCut = new TH1D( hname, "", iBinning, ( *iData ).second->fVar_min, ( *iData ).second->fVar_max );
        ( *iData ).second->hQFactors_UpperCut->SetXTitle( ( *iData ).first.c_str() );
        ( *iData ).second->hQFactors_UpperCut->SetYTitle( "q-factor" );
        setHistogramPlottingStyle( ( *iData ).second->hQFactors_UpperCut, 2, 2, 1., 1, 1, 3003 );
        
        ( *iData ).second->gQFactor_LowerCutE = new TGraphErrors( 1 );
        setGraphPlottingStyle( ( *iData ).second->gQFactor_LowerCutE, 1, 2, 22, 1. );
        ( *iData ).second->gQFactor_UpperCutE = new TGraphErrors( 1 );
        setGraphPlottingStyle( ( *iData ).second->gQFactor_UpperCutE, 2, 2, 23, 1. );
        ( *iData ).second->gQFactorMax_LowerCutE = new TGraphErrors( 1 );
        setGraphPlottingStyle( ( *iData ).second->gQFactorMax_LowerCutE, 1, 2, 22, 1. );
        ( *iData ).second->gQFactorMax_UpperCutE = new TGraphErrors( 1 );
        setGraphPlottingStyle( ( *iData ).second->gQFactorMax_UpperCutE, 2, 2, 23, 1. );
    }
    
}

void VPlotMonteCarloQualityFactor::resetHistograms()
{
    map< string, VPlotMonteCarloQualityFactorData* >::iterator iData;
    for( iData = fData.begin(); iData != fData.end(); iData++ )
    {
        if( ( *iData ).second->hSignal )
        {
            ( *iData ).second->hSignal->Reset();
        }
        if( ( *iData ).second->hBackground )
        {
            ( *iData ).second->hBackground->Reset();
        }
        if( ( *iData ).second->hQFactors_UpperCut )
        {
            ( *iData ).second->hQFactors_UpperCut->Reset();
        }
        if( ( *iData ).second->hQFactors_LowerCut )
        {
            ( *iData ).second->hQFactors_LowerCut->Reset();
        }
    }
    
}

bool VPlotMonteCarloQualityFactor::setSignalDataChain( string iChain )
{
    return setDataChain( iChain, true );
}

bool VPlotMonteCarloQualityFactor::setBackgroundDataChain( string iChain )
{
    return setDataChain( iChain, false );
}

bool VPlotMonteCarloQualityFactor::setDataChain( string iChain, bool bSignal )
{
    TChain* iC = new TChain( "data" );
    int nFil = iC->Add( iChain.c_str() );
    if( nFil == 0 )
    {
        cout << "VPlotMonteCarloQualityFactor::setDataChain: no Files in ";
        if( bSignal )
        {
            cout << "signal chain";
        }
        else
        {
            cout << "background chain";
        }
        cout << endl;
        return false;
    }
    if( bSignal )
    {
        fSignalChain     = new CData( iC );
    }
    else
    {
        fBackgroundChain = new CData( iC );
    }
    
    cout << "added " << iChain << endl;
    
    return true;
}

void VPlotMonteCarloQualityFactor::fillEnergyDependence( int iMaxNevents, double iEmin, double iEmax, double iEbin )
{
    if( !fSignalChain || !fBackgroundChain )
    {
        cout << "VPlotMonteCarloQualityFactor::fillEnergyDependence: missing data chain(s): " << fSignalChain << "\t" << fBackgroundChain << endl;
        return;
    }
    if( iEbin <= 0. )
    {
        return;
    }
    
    
    map< string, VPlotMonteCarloQualityFactorData* >::iterator iData;
    
    // loop over all energy bins
    int  iNbin = int( ( iEmax - iEmin ) / iEbin );
    cout << "Energy binning: " << iNbin << "\t" << iEmin << "\t" << iEmax << "\t" << iEbin << endl;
    for( int i = 0; i < iNbin; i++ )
    {
        resetHistograms();
        setEnergyRange( iEmin + i * iEbin, iEmin + ( i + 1 )*iEbin );
        cout << endl;
        cout << "==========================================================================" << endl;
        cout << "Energybin: " << i << "\t" << iEmin + i* iEbin << "\t" << iEmin + ( i + 1 )*iEbin << endl;
        fill( iMaxNevents );
        
        // get maximum values from q-factor histograms
        for( iData = fData.begin(); iData != fData.end(); iData++ )
        {
            ( *iData ).second->gQFactor_LowerCutE->SetPoint( i, iEmin + i * iEbin + 0.5 * iEbin, ( *iData ).second->hQFactors_LowerCut->GetBinCenter( ( *iData ).second->hQFactors_LowerCut->GetMaximumBin() ) );
            ( *iData ).second->gQFactor_UpperCutE->SetPoint( i, iEmin + i * iEbin + 0.5 * iEbin, ( *iData ).second->hQFactors_UpperCut->GetBinCenter( ( *iData ).second->hQFactors_UpperCut->GetMaximumBin() ) );
            
            ( *iData ).second->gQFactorMax_LowerCutE->SetPoint( i, iEmin + i * iEbin + 0.5 * iEbin, ( *iData ).second->hQFactors_LowerCut->GetMaximum() );
            ( *iData ).second->gQFactorMax_UpperCutE->SetPoint( i, iEmin + i * iEbin + 0.5 * iEbin, ( *iData ).second->hQFactors_UpperCut->GetMaximum() );
        }
    }
    setEnergyRange();
}


bool VPlotMonteCarloQualityFactor::fill( int iMaxNevents )
{
    if( !fSignalChain || !fBackgroundChain )
    {
        cout << "VPlotMonteCarloQualityFactor::fill: missing data chain(s): " << fSignalChain << "\t" << fBackgroundChain << endl;
        return false;
    }
    
    // fill signal histograms
    fill( iMaxNevents, fSignalChain, true );
    
    // fill background histograms
    fill( iMaxNevents, fBackgroundChain, false );
    
    // calculate q-factors
    calculateQfactors();
    
    return true;
}

void VPlotMonteCarloQualityFactor::calculateQfactors()
{
    cout << "calculating q-factors" << endl;
    
    map< string, VPlotMonteCarloQualityFactorData* >::iterator iData;
    for( iData = fData.begin(); iData != fData.end(); iData++ )
    {
        TH1D* hS = ( *iData ).second->hSignal;
        TH1D* hB = ( *iData ).second->hBackground;
        TH1D* hQL = ( *iData ).second->hQFactors_LowerCut;
        TH1D* hQU = ( *iData ).second->hQFactors_UpperCut;
        double iSTot = 0.;
        for( int i = 1; i <= hS->GetNbinsX(); i++ )
        {
            iSTot += hS->GetBinContent( i );
        }
        double iBTot = 0.;
        for( int i = 1; i <= hB->GetNbinsX(); i++ )
        {
            iBTot += hB->GetBinContent( i );
        }
        
        // calculate quality factor (upper)
        double iS = 0.;
        double iB = 0.;
        for( int i = 1; i <= hS->GetNbinsX(); i++ )
        {
            iS += hS->GetBinContent( i );
            iB += hB->GetBinContent( i );
            
            if( iB > 0. && iBTot > 0. && iSTot > 0. )
            {
                hQU->SetBinContent( i, ( iS / iSTot ) / sqrt( iB / iBTot ) );
            }
            else
            {
                hQU->SetBinContent( i, 0. );
            }
        }
        // calculate quality factor (lower)
        iS = 0.;
        iB = 0.;
        for( int i = hS->GetNbinsX(); i > 0; i-- )
        {
            iS += hS->GetBinContent( i );
            iB += hB->GetBinContent( i );
            
            if( iB > 0. && iBTot > 0. && iSTot > 0. )
            {
                hQL->SetBinContent( i, ( iS / iSTot ) / sqrt( iB / iBTot ) );
            }
            else
            {
                hQL->SetBinContent( i, 0. );
            }
        }
        // normalise distributions to number of entries
        // normalise distributions to number of entries
        if( iSTot > 0. )
        {
            hS->Scale( 1. / iSTot );
        }
        if( iBTot > 0. )
        {
            hB->Scale( 1. / iBTot );
        }
        cout << ( *iData ).first << "\t" << iSTot << "\t" << iBTot << endl;
    }
    
}


void VPlotMonteCarloQualityFactor::fill( int iMaxNevents, CData* c, bool bSignal )
{
    if( !c )
    {
        return;
    }
    
    vector< TH1D* > iHis;
    if( bSignal )
    {
        iHis.push_back( fData["MSCW"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["MSCW"]->hBackground );
    }
    if( bSignal )
    {
        iHis.push_back( fData["MSCL"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["MSCL"]->hBackground );
    }
    if( bSignal )
    {
        iHis.push_back( fData["dE"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["dE"]->hBackground );
    }
    if( bSignal )
    {
        iHis.push_back( fData["EChi2"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["EChi2"]->hBackground );
    }
    if( bSignal )
    {
        iHis.push_back( fData["EmissionHeight"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["EmissionHeight"]->hBackground );
    }
    if( bSignal )
    {
        iHis.push_back( fData["EmissionHeightChi2"]->hSignal );
    }
    else
    {
        iHis.push_back( fData["EmissionHeightChi2"]->hBackground );
    }
    
    int iNevents = c->fChain->GetEntries();
    if( iMaxNevents > 0 && iMaxNevents < iNevents )
    {
        iNevents = iMaxNevents;
    }
    cout << "Filling ";
    if( bSignal )
    {
        cout << " (signal): ";
    }
    else
    {
        cout << " (background): ";
    }
    cout << iNevents << endl;
    for( int i = 0; i < iNevents; i++ )
    {
        c->GetEntry( i );
        
        if( c->getEnergy_TeV() <= 0. )
        {
            continue;
        }
        if( log10( c->getEnergy_TeV() ) < fEnergy_min )
        {
            continue;
        }
        if( log10( c->getEnergy_TeV() ) > fEnergy_max )
        {
            continue;
        }
        
        // quality cuts
        if( c->getChi2() < 0 )
        {
            continue;
        }
        if( c->MSCW < -5. )
        {
            continue;
        }
        if( c->MSCW > 20. )
        {
            continue;
        }
        if( c->MSCL < -5. )
        {
            continue;
        }
        if( c->MSCL > 20. )
        {
            continue;
        }
        
        // fill histograms
        iHis[0]->Fill( c->MSCW );
        iHis[1]->Fill( c->MSCL );
        
        // apply basic shape cut
        if( c->MSCW < fMSCW_min )
        {
            continue;
        }
        if( c->MSCW > fMSCW_max )
        {
            continue;
        }
        if( c->MSCL < fMSCL_min )
        {
            continue;
        }
        if( c->MSCL > fMSCL_max )
        {
            continue;
        }
        
        iHis[2]->Fill( c->getEnergyDelta() );
        if( c->getEnergyChi2() > 0. )
        {
            iHis[3]->Fill( log10( c->getEnergyChi2() ) );
        }
        iHis[4]->Fill( c->EmissionHeight );
        iHis[5]->Fill( c->EmissionHeightChi2 );
    }
    cout << "Total number of entries in MSCW histogram ";
    if( bSignal )
    {
        cout << " (signal): ";
    }
    else
    {
        cout << " (background): ";
    }
    cout << iHis[0]->GetEntries() << endl;
}

void VPlotMonteCarloQualityFactor::plot( bool iPrint )
{
    char hname[800];
    char htitle[800];
    int z = 0;
    
    map< string, VPlotMonteCarloQualityFactorData* >::iterator iData;
    for( iData = fData.begin(); iData != fData.end(); iData++ )
    {
        // parameter histograms
        sprintf( hname, "cQ_%s", ( *iData ).first.c_str() );
        sprintf( htitle, "%s", ( *iData ).first.c_str() );
        TCanvas* cHis = new TCanvas( hname, htitle, 10 + z * 30, 100, 1300, 800 );
        cHis->Divide( 2, 2 );
        
        TPad* iP = ( TPad* )cHis->cd( 1 );
        iP->SetGridx( 0 );
        iP->SetGridy( 0 );
        
        ( *iData ).second->hSignal->Draw();
        ( *iData ).second->hBackground->Draw( "same" );
        
        // upper q factor
        iP = ( TPad* )cHis->cd( 2 );
        iP->SetGridx( 0 );
        iP->SetGridy( 0 );
        
        if( ( *iData ).second->hQFactors_UpperCut->GetMaximum() > ( *iData ).second->hQFactors_LowerCut->GetMaximum() )
        {
            ( *iData ).second->hQFactors_UpperCut->SetMaximum( ( *iData ).second->hQFactors_UpperCut->GetMaximum() * 1.5 );
        }
        else
        {
            ( *iData ).second->hQFactors_UpperCut->SetMaximum( ( *iData ).second->hQFactors_LowerCut->GetMaximum() * 1.5 );
        }
        
        ( *iData ).second->hQFactors_UpperCut->Draw();
        ( *iData ).second->hQFactors_LowerCut->Draw( "same" );
        
        // draw energy dependence of q-factor variable
        if( ( *iData ).second->gQFactor_LowerCutE->GetN() > 2 )
        {
            iP = ( TPad* )cHis->cd( 3 );
            iP->SetGridx( 0 );
            iP->SetGridy( 0 );
            
            ( *iData ).second->gQFactor_LowerCutE->Draw( "alp" );
            ( *iData ).second->gQFactor_LowerCutE->SetMaximum( ( *iData ).second->fVar_max );
            ( *iData ).second->gQFactor_LowerCutE->SetMinimum( ( *iData ).second->fVar_min );
            ( *iData ).second->gQFactor_LowerCutE->GetHistogram()->SetXTitle( "log_{10} energy [TeV]" );
            sprintf( hname, "value at maximum q-factor (%s)", ( *iData ).first.c_str() );
            ( *iData ).second->gQFactor_LowerCutE->GetHistogram()->SetYTitle( hname );
            ( *iData ).second->gQFactor_UpperCutE->Draw( "lp" );
            // plot mean value for all energies
            double iMaxLower = ( *iData ).second->hQFactors_LowerCut->GetBinCenter( ( *iData ).second->hQFactors_LowerCut->GetMaximumBin() );
            double iMaxUpper = ( *iData ).second->hQFactors_UpperCut->GetBinCenter( ( *iData ).second->hQFactors_UpperCut->GetMaximumBin() );
            
            TLine* iLUpper = new TLine( ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmin(), iMaxUpper, ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmax(), iMaxUpper );
            iLUpper->SetLineStyle( 2 );
            iLUpper->SetLineWidth( 2 );
            iLUpper->SetLineColor( ( *iData ).second->gQFactorMax_UpperCutE->GetLineColor() );
            iLUpper->Draw();
            
            TLine* iLLower = new TLine( ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmin(), iMaxLower, ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmax(), iMaxLower );
            iLLower->SetLineStyle( 2 );
            iLLower->SetLineWidth( 2 );
            iLLower->SetLineColor( ( *iData ).second->gQFactorMax_LowerCutE->GetLineColor() );
            iLLower->Draw();
            
        }
        // draw energy dependence of q-factor maximum
        if( ( *iData ).second->gQFactorMax_LowerCutE->GetN() > 2 )
        {
            iP = ( TPad* )cHis->cd( 4 );
            iP->SetGridx( 0 );
            iP->SetGridy( 0 );
            
            ( *iData ).second->gQFactorMax_LowerCutE->Draw( "alp" );
            ( *iData ).second->gQFactorMax_LowerCutE->SetMaximum( 5. );
            ( *iData ).second->gQFactorMax_LowerCutE->SetMinimum( 0. );
            ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->SetXTitle( "log_{10} energy [TeV]" );
            sprintf( hname, "maximum q-factor (%s)", ( *iData ).first.c_str() );
            ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->SetYTitle( hname );
            ( *iData ).second->gQFactorMax_UpperCutE->Draw( "lp" );
            
            // plot mean value for all energies
            double iMaxLower = ( *iData ).second->hQFactors_LowerCut->GetBinContent( ( *iData ).second->hQFactors_LowerCut->GetMaximumBin() );
            double iMaxUpper = ( *iData ).second->hQFactors_UpperCut->GetBinContent( ( *iData ).second->hQFactors_UpperCut->GetMaximumBin() );
            
            TLine* iLUpper = new TLine( ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmin(), iMaxUpper, ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmax(), iMaxUpper );
            iLUpper->SetLineStyle( 2 );
            iLUpper->SetLineWidth( 2 );
            iLUpper->SetLineColor( ( *iData ).second->gQFactorMax_UpperCutE->GetLineColor() );
            iLUpper->Draw();
            
            TLine* iLLower = new TLine( ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmin(), iMaxLower, ( *iData ).second->gQFactorMax_LowerCutE->GetHistogram()->GetXaxis()->GetXmax(), iMaxLower );
            iLLower->SetLineStyle( 2 );
            iLLower->SetLineWidth( 2 );
            iLLower->SetLineColor( ( *iData ).second->gQFactorMax_LowerCutE->GetLineColor() );
            iLLower->Draw();
        }
        if( iPrint )
        {
            sprintf( hname, "QFactor-%s.eps", ( *iData ).first.c_str() );
            cHis->Print( hname );
        }
        
        z++;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

VPlotMonteCarloQualityFactorData::VPlotMonteCarloQualityFactorData()
{
    fVar_max = 1.;
    fVar_min = 0.;
    hSignal = 0;
    hBackground = 0;
    hQFactors_LowerCut = 0;
    hQFactors_UpperCut = 0;
    gQFactor_LowerCutE = 0;
    gQFactorMax_UpperCutE = 0;
    gQFactorMax_LowerCutE = 0;
    gQFactor_UpperCutE = 0;
}
