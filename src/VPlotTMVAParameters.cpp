/*  \class VPlotTMVAParameters

    plot signal and background efficiency, and MVA cut variable for different
    subarrays from TMVA output files

    Example:

    .L lib/libVAnaSum.so
    VPlotTMVAParameters a;
    a.setSubArrays( "scripts/CTA/subArray.prod1-red.list");
    a.setDirectories( "$CTA_USER_DATA_DIR/analysis/AnalysisData/cta-ultra3/EffectiveArea-ID2-d20130415/QualityCuts001CU/" );
    a.initializeWeightFiles( "$CTA_USER_DATA_DIR/analysis/AnalysisData/cta-ultra3/", "TMVA/BDT-ID2-d20130415-0.0/BDT_", 0, 8 );

*/

#include "VPlotTMVAParameters.h"

VPlotTMVAParameters::VPlotTMVAParameters()
{
    fDataDirectory = "";
}

void VPlotTMVAParameters::plot( bool iPrint )
{
    char hname[2000];
    char htitle[2000];
    
    for( unsigned int i = 0; i < hSignalEfficiency.size(); i++ )
    {
    
        // signal and background efficiency
        sprintf( hname, "cTMVA_S_BC_%d", i );
        sprintf( htitle, "signal/background efficiency distribution (energy/zenith bin %d)", i );
        TCanvas* c = new TCanvas( hname, htitle, 100 + i * 20, 100 + i * 20, 400, 400 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        
        if( hSignalEfficiency[i] )
        {
            hSignalEfficiency[i]->Draw();
            cout << "Signal efficiency in energy/zenith bin " << i << ": ";
            cout << hSignalEfficiency[i]->GetMean() << " +- " << hSignalEfficiency[i]->GetRMS() << endl;
        }
        if( i < hBackgroundEfficiency.size() && hBackgroundEfficiency[i] && hBackgroundEfficiency[i]->GetEntries() > 0 )
        {
            hBackgroundEfficiency[i]->Draw( "same" );
            cout << "\t Background efficiency in energy/zenith bin " << i << ": ";
            cout << hBackgroundEfficiency[i]->GetMean() << " +- " << hBackgroundEfficiency[i]->GetRMS() << endl;
        }
        if( iPrint )
        {
            sprintf( hname, "efficiency-%d.eps", i );
            c->Print( hname );
        }
        
        
        // MVA cut variable
        sprintf( hname, "cTMVA_MVA_%d", i );
        sprintf( htitle, "MVA cut variable(energy/zenith bin %d)", i );
        TCanvas* d = new TCanvas( hname, htitle, 600 + i * 20, 100 + i * 20, 400, 400 );
        d->SetGridx( 0 );
        d->SetGridy( 0 );
        
        if( i < hMVA.size() )
        {
            hMVA[i]->Draw();
            cout << "\t MVA cut in energy/zenith bin " << i << ": ";
            cout << hMVA[i]->GetMean() << " +- " << hMVA[i]->GetRMS() << endl;
        }
        if( iPrint )
        {
            sprintf( hname, "mva-%d.eps", i );
            d->Print( hname );
        }
    }
}

bool VPlotTMVAParameters::initializeHistograms( unsigned int iEnergyWeightFileIndex_min, unsigned int iEnergyWeightFileIndex_max, unsigned int iZenithWeightFileIndex_min, unsigned int iZenithWeightFileIndex_max )
{
    char hname[2000];
    
    hSignalEfficiency.clear();
    hBackgroundEfficiency.clear();
    hMVA.clear();
    
    for( unsigned int i = iEnergyWeightFileIndex_min; i <= iEnergyWeightFileIndex_max; i++ )
    {
        for( unsigned int j = iZenithWeightFileIndex_min; j <= iZenithWeightFileIndex_max; j++ )
        {
            sprintf( hname, "hSignalEfficiency_%d_%d", i, j );
            hSignalEfficiency.push_back( new TH1D( hname, "", 100, 0., 1. ) );
            hSignalEfficiency.back()->SetXTitle( "efficiency" );
            hSignalEfficiency.back()->SetLineWidth( 2 );
            
            sprintf( hname, "hBackgroundEfficiency_%d_%d", i, j );
            hBackgroundEfficiency.push_back( new TH1D( hname, "", 100, 0., 1. ) );
            hBackgroundEfficiency.back()->SetXTitle( "efficiency" );
            hBackgroundEfficiency.back()->SetLineWidth( 2 );
            hBackgroundEfficiency.back()->SetLineColor( 2 );
            
            sprintf( hname, "hMVA_%d_%d", i, j );
            hMVA.push_back( new TH1D( hname, "", 100, -1., 1. ) );
            hMVA.back()->SetXTitle( "MVA variable" );
            hMVA.back()->SetLineWidth( 2 );
            hMVA.back()->SetLineColor( 4 );
        }
    }
    
    return true;
}

void VPlotTMVAParameters::initializeWeightFiles( string iDirectory, string iTMVADirectory, unsigned int iEnergyWeightFileIndex_min, unsigned int iEnergyWeightFileIndex_max, unsigned int iZenithWeightFileIndex_min, unsigned int iZenithWeightFileIndex_max, double iParticleNumberFile_Conversion_Rate_to_seconds )
{
    if( !initializeHistograms( iEnergyWeightFileIndex_min, iEnergyWeightFileIndex_max, iZenithWeightFileIndex_min, iZenithWeightFileIndex_max ) )
    {
        cout << "VPlotTMVAParameters::initializeWeightFiles error initializing histograms" << endl;
        return;
    }
    
    // loop over all sub arrays, get efficiency from optimization and fill them into histograms
    char hname[2000];
    for( unsigned int i = 0; i < fSubArrays.size(); i++ )
    {
        VTMVAEvaluator a;
        sprintf( hname, "%s/ParticleNumbers.%s.00.root", fDataDirectory.c_str(), fSubArrays[i].c_str() );
        a.setParticleNumberFile( hname, iParticleNumberFile_Conversion_Rate_to_seconds );
        sprintf( hname, "%s/%s/%s", iDirectory.c_str(), fSubArrays[i].c_str(), iTMVADirectory.c_str() );
        a.initializeWeightFiles( hname, iEnergyWeightFileIndex_min, iEnergyWeightFileIndex_max, iZenithWeightFileIndex_min, iZenithWeightFileIndex_max );
        
        for( unsigned int j = 0; j < a.getOptimumCutValueFound().size(); j++ )
        {
            if( a.getOptimumCutValueFound()[j] )
            {
                if( j < a.getSignalEfficiency().size() && j < hSignalEfficiency.size() && hSignalEfficiency[j] )
                {
                    hSignalEfficiency[j]->Fill( a.getSignalEfficiency()[j] );
                }
                if( j < a.getBackgroundEfficiency().size() && j < hBackgroundEfficiency.size() && hBackgroundEfficiency[j] )
                {
                    hBackgroundEfficiency[j]->Fill( a.getBackgroundEfficiency()[j] );
                }
                if( j < a.getTMVACutValue().size() && j < hMVA.size() && hMVA[j] )
                {
                    hMVA[j]->Fill( a.getTMVACutValue()[j] );
                }
            }
        }
    }
}

bool VPlotTMVAParameters::setSubArrays( string iFileTxt )
{
    fSubArrays.clear();
    
    ifstream is;
    is.open( iFileTxt.c_str(), ifstream::in );
    if( !is )
    {
        return false;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        fSubArrays.push_back( is_line );
    }
    
    if( fSubArrays.size() == 0 )
    {
        cout << "VPlotTMVAParameters::setSubArrays: no arrays found" << endl;
        cout << "\t " << iFileTxt << endl;
        return false;
    }
    
    cout << "found " << fSubArrays.size() << " sub arrays" << endl;
    
    return true;
}

