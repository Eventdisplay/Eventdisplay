/*  plot migration matrix for a certain offset
 *
 *  */
TH2F* plotMigrationMatrix( string iFile, double iOffset )
{
    // open file and get migration matrix (3D)
    TFile* f = new TFile( iFile.c_str() );
    if( f->IsZombie() )
    {
        return 0;
    }
    TH3F* hEdisp = ( TH3F* )f->Get( "EestOverEtrueNoTheta2cut_offaxis" );
    if( !hEdisp )
    {
        return 0;
    }
    
    char hname[200];
    sprintf( hname, "Edisp_%d", ( int )( iOffset * 100 ) );
    TH2F* hP = new TH2F( hname, "", hEdisp->GetXaxis()->GetNbins(), hEdisp->GetXaxis()->GetXmin(), hEdisp->GetXaxis()->GetXmax(),
                         hEdisp->GetYaxis()->GetNbins(), hEdisp->GetYaxis()->GetXmin(), hEdisp->GetYaxis()->GetXmax() );
    hP->SetXTitle( hEdisp->GetXaxis()->GetTitle() );
    hP->SetYTitle( hEdisp->GetYaxis()->GetTitle() );
    
    int nZ = hEdisp->GetZaxis()->FindBin( iOffset );
    cout << "nZ " << nZ << endl;
    
    for( int i = 1; i <= hEdisp->GetXaxis()->GetNbins(); i++ )
    {
        for( int j = 1; j <= hEdisp->GetYaxis()->GetNbins(); j++ )
        {
            hP->SetBinContent( i, j, hEdisp->GetBinContent( i, j, nZ ) );
        }
    }
    
    hP->Draw( "colz" );
    
    
    return hP;
}


