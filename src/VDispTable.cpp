/*! \file VDispTable.cpp
    \brief fill and read tables for angular reconstruction with disp method


*/

#include "VDispTable.h"

VDispTable::VDispTable()
{
    fDebug = false;
    
    fNTel = 0;
    fTableFile = 0;
    fData = 0;
    
    setQualityCuts();
    setWidthLengthScalingParameters();
}


VDispTable::VDispTable( unsigned int iNTel, string ioutFile )
{
    fDebug = false;
    
    fNTel = iNTel;
    
    fTableFile = 0;
    fData = 0;
    
    setQualityCuts();
    setWidthLengthScalingParameters();
    
    prepareTraining( ioutFile );
}


bool VDispTable::prepareTraining( string iOutFile )
{
    // create output file
    fTableFile = new TFile( iOutFile.c_str(), "RECREATE" );
    cout << "creating output table file: " << fTableFile->GetName() << endl;
    if( fTableFile->IsZombie() )
    {
        cout << "VDispTable::prepareTraining error creating output file: " << iOutFile << endl;
        cout << "...exiting" << endl;
        exit( 0 );
    }
    
    // create data class
    fData = new VDispTableReader();
    fData->SetName( "dispTable" );
    fData->initialize( false );
    
    return true;
}

/*!
    write data class to file
*/
bool VDispTable::terminate()
{
    if( fTableFile && fData )
    {
        cout << endl;
        cout << "writing data to " << fTableFile->GetName() << endl;
        fTableFile->cd();
        fData->Write();
        fData->print( true );
    }
    if( fTableFile )
    {
        fTableFile->Close();
    }
    
    return true;
}


void VDispTable::addAzBin( float iMin, float iMax )
{
    fAz_min.push_back( iMin );
    fAz_max.push_back( iMax );
    
    char hname[500];
    char htitle[500];
    
    // 2D disp tables
    if( fData->h2D_DispTable )
    {
        sprintf( hname, "h2D_DispTable_%d", ( int )fAz_min.size() );
        sprintf( htitle, "az bin (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h2D_AzDispTable.push_back( new TProfile2D( hname, htitle, fData->h2D_DispTable->GetNbinsX(), fData->h2D_DispTable->GetXaxis()->GetXmin(), fData->h2D_DispTable->GetXaxis()->GetXmax(), fData->h2D_DispTable->GetNbinsY(), fData->h2D_DispTable->GetYaxis()->GetXmin(), fData->h2D_DispTable->GetYaxis()->GetXmax(), 0., 100. ) );
        h2D_AzDispTable.back()->SetXTitle( "width/length" );
        h2D_AzDispTable.back()->SetYTitle( "log_{10} size" );
        h2D_AzDispTable.back()->SetZTitle( "disp [deg]" );
    }
    if( fData->h2D_DispTableN )
    {
        sprintf( hname, "h2D_DispTableN_%d", ( int )fAz_min.size() );
        h2D_AzDispTableN.push_back( ( TH2F* )fData->h2D_DispTableN->Clone( hname ) );
        h2D_AzDispTableN.back()->SetTitle( "" );
        h2D_AzDispTableN.back()->SetXTitle( "width/length" );
        h2D_AzDispTableN.back()->SetYTitle( "log_{10} size" );
        h2D_AzDispTableN.back()->SetZTitle( "number of events/bin" );
    }
    if( fData->h2D_DispPhiTable )
    {
        sprintf( hname, "h2D_DispPhiTable_%d", ( int )fAz_min.size() );
        sprintf( htitle, "az bin (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h2D_AzDispPhiTable.push_back( new TProfile2D( hname, htitle, fData->h2D_DispPhiTable->GetNbinsX(), fData->h2D_DispPhiTable->GetXaxis()->GetXmin(), fData->h2D_DispPhiTable->GetXaxis()->GetXmax(), fData->h2D_DispPhiTable->GetNbinsY(), fData->h2D_DispPhiTable->GetYaxis()->GetXmin(), fData->h2D_DispPhiTable->GetYaxis()->GetXmax(), 0., 100. ) );
        h2D_AzDispPhiTable.back()->SetXTitle( "width/length" );
        h2D_AzDispPhiTable.back()->SetYTitle( "log_{10} size" );
        h2D_AzDispPhiTable.back()->SetZTitle( "phi [deg]" );
    }
    if( fData->h2D_DispMissTable )
    {
        sprintf( hname, "h2D_DispMiss_%d", ( int )fAz_min.size() );
        sprintf( htitle, "error in reconstruction (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h2D_AzDispMissTable.push_back( new TProfile2D( hname, htitle, fData->h2D_DispMissTable->GetNbinsX(), fData->h2D_DispMissTable->GetXaxis()->GetXmin(), fData->h2D_DispMissTable->GetXaxis()->GetXmax(), fData->h2D_DispMissTable->GetNbinsY(), fData->h2D_DispMissTable->GetYaxis()->GetXmin(), fData->h2D_DispMissTable->GetYaxis()->GetXmax(), 0., 100. ) );
        h2D_AzDispMissTable.back()->SetTitle( htitle );
        h2D_AzDispMissTable.back()->SetXTitle( "width/length" );
        h2D_AzDispMissTable.back()->SetYTitle( "log_{10} size" );
        h2D_AzDispMissTable.back()->SetZTitle( "miss" );
    }
    
    // 3D disp tables
    if( fData->h3D_DispTable )
    {
        sprintf( hname, "h3D_DispTable_%d", ( int )fAz_min.size() );
        sprintf( htitle, "az bin (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h3D_AzDispTable.push_back( new TProfile3D( hname, htitle, fData->h3D_DispTable->GetNbinsX(), fData->h3D_DispTable->GetXaxis()->GetXmin(), fData->h3D_DispTable->GetXaxis()->GetXmax(), fData->h3D_DispTable->GetNbinsY(), fData->h3D_DispTable->GetYaxis()->GetXmin(), fData->h3D_DispTable->GetYaxis()->GetXmax(), fData->h3D_DispTable->GetNbinsZ(), fData->h3D_DispTable->GetZaxis()->GetXmin(), fData->h3D_DispTable->GetZaxis()->GetXmax() ) );
        h3D_AzDispTable.back()->SetXTitle( "f(width)" );
        h3D_AzDispTable.back()->SetXTitle( "f(length)" );
        h3D_AzDispTable.back()->SetYTitle( "log_{10} size" );
        h3D_AzDispTable.back()->SetZTitle( "disp [deg]" );
    }
    
    if( fData->h3D_DispTableN )
    {
        sprintf( hname, "h3D_DispTableN_%d", ( int )fAz_min.size() );
        h3D_AzDispTableN.push_back( ( TH3F* )fData->h3D_DispTableN->Clone( hname ) );
        h3D_AzDispTableN.back()->SetTitle( "" );
        h3D_AzDispTableN.back()->SetXTitle( "f(width)" );
        h3D_AzDispTableN.back()->SetXTitle( "f(length)" );
        h3D_AzDispTableN.back()->SetYTitle( "log_{10} size" );
        h3D_AzDispTableN.back()->SetZTitle( "number of events/bin" );
    }
    
    if( fData->h3D_DispPhiTable )
    {
        sprintf( hname, "h3D_DispPhiTable_%d", ( int )fAz_min.size() );
        sprintf( htitle, "az bin (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h3D_AzDispPhiTable.push_back( new TProfile3D( hname, htitle, fData->h3D_DispPhiTable->GetNbinsX(), fData->h3D_DispPhiTable->GetXaxis()->GetXmin(), fData->h3D_DispPhiTable->GetXaxis()->GetXmax(), fData->h3D_DispPhiTable->GetNbinsY(), fData->h3D_DispPhiTable->GetYaxis()->GetXmin(), fData->h3D_DispPhiTable->GetYaxis()->GetXmax(), fData->h3D_DispPhiTable->GetNbinsZ(), fData->h3D_DispPhiTable->GetZaxis()->GetXmin(), fData->h3D_DispPhiTable->GetZaxis()->GetXmax() ) );
        h3D_AzDispPhiTable.back()->SetXTitle( "f(width)" );
        h3D_AzDispPhiTable.back()->SetXTitle( "f(length)" );
        h3D_AzDispPhiTable.back()->SetYTitle( "log_{10} size" );
        h3D_AzDispPhiTable.back()->SetZTitle( "phi [deg]" );
    }
    
    if( fData->h3D_DispMissTable )
    {
        sprintf( hname, "h3D_DispMiss_%d", ( int )fAz_min.size() );
        sprintf( htitle, "error in reconstruction (%f < az < %f; %d)", fAz_min.back(), fAz_max.back(), ( int )fAz_min.size() );
        h3D_AzDispMissTable.push_back( new TProfile3D( hname, htitle, fData->h3D_DispMissTable->GetNbinsX(), fData->h3D_DispMissTable->GetXaxis()->GetXmin(), fData->h3D_DispMissTable->GetXaxis()->GetXmax(), fData->h3D_DispMissTable->GetNbinsY(), fData->h3D_DispMissTable->GetYaxis()->GetXmin(), fData->h3D_DispMissTable->GetYaxis()->GetXmax(), fData->h3D_DispMissTable->GetNbinsZ(), fData->h3D_DispMissTable->GetZaxis()->GetXmin(), fData->h3D_DispMissTable->GetZaxis()->GetXmax() ) );
        //       h3D_AzDispMissTable.push_back( new TProfile3D( hname, htitle, iScaleWidth_bin, 0., 1., iScaleLength_bin, 0., 1., iSize_bin, iSize_min, iSize_max ) );
        h3D_AzDispMissTable.back()->SetXTitle( "f(width)" );
        h3D_AzDispMissTable.back()->SetXTitle( "f(length)" );
        h3D_AzDispMissTable.back()->SetYTitle( "log_{10} size" );
        h3D_AzDispMissTable.back()->SetZTitle( "miss" );
    }
    
}


void VDispTable::setDataVectors( vector< string > i_ze, vector< string > i_woff, vector< string > i_noise )
{
    if( fData )
    {
        fData->setDataVectors( i_ze, fAz_min, fAz_max, i_woff, i_noise );
        fData->setWidthLengthScalingParameters( fWidthScaleParameter, fLengthScaleParameter );
    }
}

void VDispTable::setQualityCuts( int iNtubes_min, double iSize_min, double iLength_min, double iLoss_max )
{
    fQC_Ntubes_min = iNtubes_min;
    fQC_Size_min = iSize_min;
    fQC_Length_min = iLength_min;
    fQC_Loss_max = iLoss_max;
}

/*!
    apply quality cuts
*/
bool VDispTable::isGoodEvent( Ctpars* c )
{
    if( !c )
    {
        return false;
    }
    
    if( c->ntubes <= fQC_Ntubes_min )
    {
        return false;
    }
    if( c->size <= fQC_Size_min )
    {
        return false;
    }
    if( c->length <= fQC_Length_min )
    {
        return false;
    }
    if( c->loss >= fQC_Loss_max )
    {
        return false;
    }
    
    return true;
}


/*!
   loop over MC data and fill disp tables for given ze, woff
*/
bool VDispTable::fillTable( string iMCFile, float i_ze, float i_woff, int iNentries )
{
    char hname[600];
    
    /////////////////////////////////////////////////////////////////////////////////////
    // reset variables
    /////////////////////////////////////////////////////////////////////////////////////
    float disp = 0.;
    float dispPhi = 0.;
    unsigned int azBin = 0;
    
    // assume same pedvars for all telescopes (is this true in the simulations?)
    float i_meanPedvars = 0.;
    float i_meanPedvarsN = 0.;
    
    // reset histograms
    fData->reset();
    for( unsigned int i = 0; i < h2D_AzDispTable.size(); i++ )
    {
        h2D_AzDispTable[i]->Reset();
        h2D_AzDispPhiTable[i]->Reset();
        h2D_AzDispTableN[i]->Reset();
        h2D_AzDispMissTable[i]->Reset();
        
        h3D_AzDispTable[i]->Reset();
        h3D_AzDispPhiTable[i]->Reset();
        h3D_AzDispTableN[i]->Reset();
        h3D_AzDispMissTable[i]->Reset();
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    // prepare MC data files and trees
    /////////////////////////////////////////////////////////////////////////////////////
    
    // data file with MC input
    TFile iFile( iMCFile.c_str() );
    if( iFile.IsZombie() )
    {
        cout << "VDispTable::fillTable error reading MC file " << iMCFile << endl;
        return false;
    }
    
    cout << "filling tables for " << iMCFile << endl;
    
    // get showerpars tree
    TTree* s = ( TTree* )iFile.Get( "showerpars" );
    if( !s )
    {
        cout << "VDispTable::fillTable error finding tree showerpars" << endl;
        return false;
    }
    Cshowerpars* m = new Cshowerpars( s, true, true );
    if( !m )
    {
        return false;
    }
    
    // get tpars trees (fNTel times)
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        sprintf( hname, "Tel_%d/tpars", u + 1 );
        TTree* t = ( TTree* )iFile.Get( hname );
        if( !t )
        {
            cout << "VDispTable::fillTable error finding tree tpars for telescope " << i + 1 << endl;
            continue;
        }
        Ctpars* c = new Ctpars( t, true, 1 );
        if( !c )
        {
            continue;
        }
        
        if( m->fChain->GetEntries() != c->fChain->GetEntries() )
        {
            cout << "\t error: different number of events in showerpars and tpars trees for telescope " << i + 1 << endl;
            cout << "...ignore this telescope" << endl;
            continue;
        }
        cout << "\t telescope " << i + 1 << ", " << c->fChain->GetEntries() << " entries in data tree" << endl;
        
        /////////////////////////////////////////////////////////////////////////////////////
        // loop over all events for this telescope
        if( iNentries < 0 )
        {
            iNentries = c->fChain->GetEntries();
        }
        if( iNentries > c->fChain->GetEntries() )
        {
            iNentries = c->fChain->GetEntries();
        }
        cout << "\t (looping over " << iNentries << " entries)" << endl;
        for( int n = 0; n < iNentries; n++ )
        {
            c->GetEntry( n );
            m->GetEntry( n );
            
            // apply quality cuts
            if( isGoodEvent( c ) )
            {
                i_meanPedvars += c->meanPedvar_Image;
                i_meanPedvarsN++;
                
                // calculate disp (observe sign convention for MCyoff)
                disp  = sqrt( ( c->cen_y + m->MCyoff ) * ( c->cen_y + m->MCyoff ) + ( c->cen_x - m->MCxoff ) * ( c->cen_x - m->MCxoff ) );
                //disp  = sqrt( ( c->cen_y - m->MCyoff ) * ( c->cen_y - m->MCyoff ) + ( c->cen_x - m->MCxoff ) * ( c->cen_x - m->MCxoff ) );
                disp  = TMath::Abs( disp );
                // dispPhi is always between 0 and 90 deg -> estimation of error of major axis
                dispPhi = TMath::Abs( TMath::ATan2( c->sinphi, c->cosphi ) - TMath::ATan2( c->cen_y + m->MCyoff, c->cen_x - m->MCxoff ) );
                if( dispPhi > TMath::PiOver2() * 3. )
                {
                    dispPhi  = 2.*TMath::Pi() - dispPhi;
                }
                else if( dispPhi > TMath::Pi() )
                {
                    dispPhi -= TMath::Pi();
                }
                else if( dispPhi > TMath::PiOver2() )
                {
                    dispPhi  = TMath::Pi() - dispPhi;
                }
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                // fill tables
                azBin = fData->getAzBin( m->MCaz );
                if( azBin < h2D_AzDispTable.size() )
                {
                    // 2D tables
                    h2D_AzDispTable[azBin]->Fill( c->width / c->length, log10( c->size ), disp );
                    h2D_AzDispTableN[azBin]->Fill( c->width / c->length, log10( c->size ) );
                    h2D_AzDispPhiTable[azBin]->Fill( c->width / c->length, log10( c->size ), dispPhi );
                    // 3D tables
                    h3D_AzDispTable[azBin]->Fill( fData->scaleWidthParameter( c->width ), fData->scaleLengthParameter( c->length ), log10( c->size ), disp );
                    h3D_AzDispTableN[azBin]->Fill( fData->scaleWidthParameter( c->width ), fData->scaleLengthParameter( c->length ), log10( c->size ) );
                    h3D_AzDispPhiTable[azBin]->Fill( fData->scaleWidthParameter( c->width ), fData->scaleLengthParameter( c->length ), log10( c->size ), dispPhi );
                }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
    }
    if( i_meanPedvarsN > 0. )
    {
        i_meanPedvars /= i_meanPedvarsN;
    }
    cout << "\t mean pedvars: " << i_meanPedvars << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // loop a second time over all events to estimate error in reconstruction
    cout << endl << "second loop for error estimation..." << endl;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        sprintf( hname, "Tel_%d/tpars", i + 1 );
        TTree* t = ( TTree* )iFile.Get( hname );
        if( !t )
        {
            cout << "VDispTable::fillTable error finding tree tpars for telescope " << i + 1 << endl;
            continue;
        }
        Ctpars* c = new Ctpars( t, true, 1 );
        if( !c )
        {
            continue;
        }
        
        if( m->fChain->GetEntries() != c->fChain->GetEntries() )
        {
            cout << "\t error: different number of events in showerpars and tpars trees for telescope " << i + 1 << endl;
            cout << "...ignore this telescope" << endl;
            continue;
        }
        cout << "\t telescope " << i + 1 << ", " << c->fChain->GetEntries() << " entries in data tree" << endl;
        
        // loop over all events for this telescope
        double miss = 0.;
        int x_bin = 0;
        int y_bin = 0;
        int z_bin = 0;
        double x_s = 0.;
        double y_s = 0.;
        if( iNentries < 0 )
        {
            iNentries = c->fChain->GetEntries();
        }
        if( iNentries > c->fChain->GetEntries() )
        {
            iNentries = c->fChain->GetEntries();
        }
        for( int n = 0; n < iNentries; n++ )
        {
            c->GetEntry( n );
            m->GetEntry( n );
            
            // quality cuts
            if( isGoodEvent( c ) )
            {
                azBin = fData->getAzBin( m->MCaz );
                if( azBin < h2D_AzDispTable.size() )
                {
                    //////////////////////////////////////////////////////////////////////////////////////
                    // 2D
                    // get disp values from histograms
                    x_bin = h2D_AzDispTable[azBin]->GetXaxis()->FindBin( c->width / c->length );
                    y_bin = h2D_AzDispTable[azBin]->GetYaxis()->FindBin( log10( c->size ) );
                    disp =  h2D_AzDispTable[azBin]->GetBinContent( x_bin, y_bin );
                    // calculate reconstructed shower direction
                    x_s = c->cen_x - disp * c->cosphi;
                    y_s = c->cen_y - disp * c->sinphi;
                    // calculate deviation from true shower direction and fill this into a tree
                    miss = sqrt( ( m->MCxoff - x_s ) * ( m->MCxoff - x_s ) + ( m->MCyoff + y_s ) * ( m->MCyoff + y_s ) );
                    h2D_AzDispMissTable[azBin]->Fill( c->width / c->length, log10( c->size ), miss );
                    //////////////////////////////////////////////////////////////////////////////////////
                    // 3D
                    // get disp values from histograms
                    x_bin = h3D_AzDispTable[azBin]->GetXaxis()->FindBin( fData->scaleWidthParameter( c->width ) );
                    y_bin = h3D_AzDispTable[azBin]->GetXaxis()->FindBin( fData->scaleLengthParameter( c->length ) );
                    z_bin = h3D_AzDispTable[azBin]->GetYaxis()->FindBin( log10( c->size ) );
                    disp =  h3D_AzDispTable[azBin]->GetBinContent( x_bin, y_bin, z_bin );
                    // calculate reconstructed shower direction
                    x_s = c->cen_x - disp * c->cosphi;
                    y_s = c->cen_y - disp * c->sinphi;
                    // calculate deviation from true shower direction and fill this into a tree
                    miss = sqrt( ( m->MCxoff - x_s ) * ( m->MCxoff - x_s ) + ( m->MCyoff + y_s ) * ( m->MCyoff + y_s ) );
                    h3D_AzDispMissTable[azBin]->Fill( fData->scaleWidthParameter( c->width ), fData->scaleLengthParameter( c->length ), log10( c->size ), miss );
                }
            }
        }
    }
    
    // fill data carrier
    for( unsigned int i = 0; i < fAz_min.size(); i++ )
    {
        fData->fill( i_ze, i, fAz_min[i], fAz_max[i], i_woff, i_meanPedvars, ( TH2* )h2D_AzDispTable[i], ( TH2* )h2D_AzDispTableN[i], ( TH2* )h2D_AzDispPhiTable[i], ( TH2* )h2D_AzDispMissTable[i], ( TH3* )h3D_AzDispTable[i], ( TH3* )h3D_AzDispTableN[i], ( TH3* )h3D_AzDispPhiTable[i], ( TH3* )h3D_AzDispMissTable[i] );
    }
    
    return true;
}
