void plotQualityCuts( string iMSCWFile, string iEnergy = "MCe0" )
{
    TFile* f = new TFile( iMSCWFile.c_str() );
    if( f->IsZombie() )
    {
        return;
    }
    TTree* data = ( TTree* )f->Get( "data" );
    
    vector< string > iCutName;
    vector< string > iCut;
    iCutName.push_back( "noCut" );
    iCut.push_back( "" );
    iCutName.push_back( "XYcut" );
    iCut.push_back( "(Xoff*Xoff+Yoff*Yoff)<2.*2." );
    iCutName.push_back( "Chi2" );
    iCut.push_back( "Chi2>=0." );
    iCutName.push_back( "NImages" );
    iCut.push_back( "NImages>1" );
    iCutName.push_back( "MSC" );
    iCut.push_back( "MSCW>-50.&&MSCL>-50." );
    iCutName.push_back( "EChi2" );
    iCut.push_back( "EChi2S>=0.&&EChi2S<=99999" );
    iCutName.push_back( "dES" );
    iCut.push_back( "dES>-99&dES<1e+12" );
    iCutName.push_back( "ErecS" );
    iCut.push_back( "ErecS>0.&&ErecS<1.e10" );
    iCutName.push_back( "Size" );
    iCut.push_back( "SizeSecondMax>200." );
    
    vector< TH1D* > iCutHistogram;
    
    // histograms in true energy with different cuts
    char hname[200];
    for( unsigned int i = 0; i < iCutName.size(); i++ )
    {
        sprintf( hname, "h_%s", iCutName[i].c_str() );
        iCutHistogram.push_back( new TH1D( hname, "", 60, -2., 4. ) );
        iCutHistogram.back()->SetStats( 0 );
        sprintf( hname, "log_{10} %s (TeV)", iEnergy.c_str() );
        iCutHistogram.back()->GetXaxis()->SetTitle( hname );
        iCutHistogram.back()->SetLineColor( i + 1 );
        iCutHistogram.back()->SetLineWidth( 2 );
        iCutHistogram.back()->SetMinimum( 0.5 );
        
        sprintf( hname, "log10(%s)", iEnergy.c_str() );
        data->Project( iCutHistogram.back()->GetName(), hname, iCut[i].c_str() );
        cout << iCutName[i] << "\t color " << iCutHistogram.back()->GetLineColor() << endl;
    }
    
    // plot everything
    //
    TCanvas* c = new TCanvas( "cE", "energy spectra", 10, 10, 800, 600 );
    c->Draw();
    for( unsigned int i = 0; i < iCutHistogram.size(); i++ )
    {
        if( i == 0 )
        {
            iCutHistogram[i]->DrawCopy();
        }
        else
        {
            iCutHistogram[i]->DrawCopy( "same" );
        }
    }
    TCanvas* cR = new TCanvas( "cERatioe", "energy spectra (ratio)", 100, 10, 800, 600 );
    cR->Draw();
    for( unsigned int i = 1; i < iCutHistogram.size(); i++ )
    {
        iCutHistogram[i]->Divide( iCutHistogram[i], iCutHistogram[0] );
        if( i == 1 )
        {
            iCutHistogram[i]->Draw();
        }
        else
        {
            iCutHistogram[i]->Draw( "same" );
        }
    }
}

