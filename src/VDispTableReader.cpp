/*! \file VDispTableReaderReader.cpp
    \brief read tables for angular reconstruction with disp method



*/

#include "VDispTableReader.h"

ClassImp( VDispTableReader )

VDispTableReader::VDispTableReader()
{
    hisList = 0;
    fData = 0;
    ze = 0.;
    az_bin = 0;
    az_min = 0.;
    az_max = 0.;
    woff = 0.;
    pedvar = 0.;
    
    h2D_DispTable = 0;
    h2D_DispPhiTable = 0;
    h2D_DispTableN = 0;
    h2D_DispMissTable = 0;
    
    h3D_DispTable = 0;
    h3D_DispPhiTable = 0;
    h3D_DispTableN = 0;
    h3D_DispMissTable = 0;
    
    fHisto_ListOfVariables.push_back( "WidthOverLength" );
    fHisto_ListOfVariables.push_back( "ScaleLength" );
    fHisto_ListOfVariables.push_back( "ScaleWidth" );
    fHisto_ListOfVariables.push_back( "Size" );
    
    setHistoBinning( "WidthOverLength", 50, 0., 1. );
    setHistoBinning( "ScaleLength", 50, 0., 1. );
    setHistoBinning( "ScaleWidth", 50, 0., 1. );
    setHistoBinning( "Size", 35, 0., 8. );
    
    setWidthLengthScalingParameters();
}

bool VDispTableReader::isHistoBinningSet( string iVariableName )
{
    if( fHisto_binning.find( iVariableName ) == fHisto_binning.end() )
    {
        return false;
    }
    if( fHisto_min.find( iVariableName ) == fHisto_min.end() )
    {
        return false;
    }
    if( fHisto_max.find( iVariableName ) == fHisto_max.end() )
    {
        return false;
    }
    
    return true;
}

bool VDispTableReader::setHistoBinning( string iVariableName, int iBin, double iMin, double iMax )
{
    // test variable name
    bool iFound = false;
    for( unsigned int i = 0; i < fHisto_ListOfVariables.size(); i++ )
    {
        if( fHisto_ListOfVariables[i] == iVariableName )
        {
            iFound = true;
            break;
        }
    }
    if( !iFound )
    {
        cout << "VDispTableReader::setHistoBinning: invalid variable name; allowed are: " << endl;
        for( unsigned int i = 0; i < fHisto_ListOfVariables.size(); i++ )
        {
            cout << "\t" << fHisto_ListOfVariables[i] << endl;
        }
        return false;
    }
    
    // set histogram variables
    fHisto_binning[iVariableName] = iBin;
    fHisto_min[iVariableName] = iMin;
    fHisto_max[iVariableName] = iMax;
    
    return true;
}



/*!
    plot the 2D histograms
*/
bool VDispTableReader::plot( int iEntry )
{
    if( fData && iEntry < fData->GetEntries() )
    {
        fData->GetEntry( iEntry );
        
        if( h2D_DispTable )
        {
            TCanvas* c = new TCanvas( "cDisp", "disp table", 10, 10, 600, 400 );
            c->Draw();
            
            h2D_DispTable->SetTitle( "" );
            h2D_DispTable->SetStats( 0 );
            h2D_DispTable->SetXTitle( "width/length" );
            h2D_DispTable->SetYTitle( "log_{10} size" );
            h2D_DispTable->Draw( "colz" );
            
            // draw errors
            c = new TCanvas( "cDisperror", "disp table (errors)", 110, 10, 600, 400 );
            c->Draw();
            TH2D* h2D_DispTableError = ( TH2D* )h2D_DispTable->Clone( "h2D_DispTableError" );
            for( int i = 1; i <= h2D_DispTableError->GetNbinsX(); i++ )
            {
                for( int j = 1; j <= h2D_DispTableError->GetNbinsY(); j++ )
                {
                    if( h2D_DispTableError->GetBinContent( i, j ) > 0. )
                    {
                        h2D_DispTableError->SetBinContent( i, j, h2D_DispTableError->GetBinError( i, j ) / h2D_DispTableError->GetBinContent( i, j ) );
                    }
                }
            }
            h2D_DispTableError->SetTitle( "" );
            h2D_DispTableError->SetStats( 0 );
            h2D_DispTableError->SetXTitle( "width/length" );
            h2D_DispTableError->SetYTitle( "log_{10} size" );
            h2D_DispTableError->Draw( "colz" );
        }
        if( h2D_DispTableN )
        {
            TCanvas* c = new TCanvas( "cN", "disp table (entries)", 610, 10, 600, 400 );
            c->Draw();
            
            h2D_DispTableN->SetTitle( "" );
            h2D_DispTableN->SetStats( 0 );
            h2D_DispTableN->SetXTitle( "width/length" );
            h2D_DispTableN->SetYTitle( "log_{10} size" );
            h2D_DispTableN->Draw( "colz" );
        }
        if( h2D_DispPhiTable )
        {
            TCanvas* c = new TCanvas( "cPhi", "dispPhi table", 110, 210, 600, 400 );
            c->Draw();
            
            h2D_DispPhiTable->SetTitle( "" );
            h2D_DispPhiTable->SetStats( 0 );
            h2D_DispPhiTable->SetXTitle( "width/length" );
            h2D_DispPhiTable->SetYTitle( "log_{10} size" );
            h2D_DispPhiTable->Draw( "colz" );
        }
        if( h2D_DispMissTable )
        {
            TCanvas* c = new TCanvas( "cMiss", "dispMiss table", 110, 310, 600, 400 );
            c->Draw();
            
            h2D_DispMissTable->SetTitle( "" );
            h2D_DispMissTable->SetStats( 0 );
            h2D_DispMissTable->SetXTitle( "width/length" );
            h2D_DispMissTable->SetYTitle( "log_{10} size" );
            h2D_DispMissTable->Draw( "colz" );
        }
        else
        {
            cout << "h2D_DispMissTable not found" << endl;
        }
    }
    return true;
}


bool VDispTableReader::initialize( bool iRead )
{
    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // prepare data tree for writing
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    if( !iRead )
    {
        if( !isHistoBinningSet( "WidthOverLength" ) )
        {
            return false;
        }
        if( !isHistoBinningSet( "ScaleLength" ) )
        {
            return false;
        }
        if( !isHistoBinningSet( "ScaleWidth" ) )
        {
            return false;
        }
        if( !isHistoBinningSet( "Size" ) )
        {
            return false;
        }
        
        // list of histograms
        hisList = new TList;
        
        h2D_DispTable = new TH2F( "h2D_DispTable", "disp table", fHisto_binning["WidthOverLength"], fHisto_min["WidthOverLength"], fHisto_max["WidthOverLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h2D_DispTable );
        }
        h2D_DispPhiTable = new TH2F( "h2D_DispPhiTable", "disp phi table", fHisto_binning["WidthOverLength"], fHisto_min["WidthOverLength"], fHisto_max["WidthOverLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h2D_DispPhiTable );
        }
        h2D_DispMissTable = new TH2F( "h2D_DispMissTable", "disp miss table", fHisto_binning["WidthOverLength"], fHisto_min["WidthOverLength"], fHisto_max["WidthOverLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h2D_DispMissTable );
        }
        h2D_DispTableN = new TH2F( "h2D_DispTableN", "disp table (entries)", fHisto_binning["WidthOverLength"], fHisto_min["WidthOverLength"], fHisto_max["WidthOverLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h2D_DispTableN );
        }
        
        h3D_DispTable = new TH3F( "h3D_DispTable", "disp table", fHisto_binning["ScaleWidth"], fHisto_min["ScaleWidth"], fHisto_max["ScaleWidth"], fHisto_binning["ScaleLength"], fHisto_min["ScaleLength"], fHisto_max["ScaleLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h3D_DispTable );
        }
        h3D_DispPhiTable = new TH3F( "h3D_DispPhiTable", "disp phi table", fHisto_binning["ScaleWidth"], fHisto_min["ScaleWidth"], fHisto_max["ScaleWidth"], fHisto_binning["ScaleLength"], fHisto_min["ScaleLength"], fHisto_max["ScaleLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h3D_DispPhiTable );
        }
        h3D_DispMissTable = new TH3F( "h3D_DispMissTable", "disp miss table", fHisto_binning["ScaleWidth"], fHisto_min["ScaleWidth"], fHisto_max["ScaleWidth"], fHisto_binning["ScaleLength"], fHisto_min["ScaleLength"], fHisto_max["ScaleLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h3D_DispMissTable );
        }
        h3D_DispTableN = new TH3F( "h3D_DispTableN", "disp table (entries)", fHisto_binning["ScaleWidth"], fHisto_min["ScaleWidth"], fHisto_max["ScaleWidth"], fHisto_binning["ScaleLength"], fHisto_min["ScaleLength"], fHisto_max["ScaleLength"], fHisto_binning["Size"], fHisto_min["Size"], fHisto_max["Size"] );
        if( hisList )
        {
            hisList->Add( h3D_DispTableN );
        }
        
        fData = new TTree( "disp", "disp tables" );
        fData->Branch( "ze", &ze, "ze/F" );
        fData->Branch( "az_bin", &az_bin, "az_bin/i" );
        fData->Branch( "az_min", &az_min, "az_min/F" );
        fData->Branch( "az_max", &az_max, "az_max/F" );
        fData->Branch( "woff", &woff, "woff/F" );
        fData->Branch( "pedvar", &pedvar, "pedvar/F" );
        if( hisList )
        {
            fData->Branch( hisList, 32000, 0 );
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // prepare data tree for reading
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    else
    {
        if( fData )
        {
            fData->SetBranchAddress( "ze", &ze );
            fData->SetBranchAddress( "az_bin", &az_bin );
            fData->SetBranchAddress( "az_min", &az_min );
            fData->SetBranchAddress( "az_max", &az_max );
            fData->SetBranchAddress( "woff", &woff );
            fData->SetBranchAddress( "pedvar", &pedvar );
            
            fData->SetBranchAddress( "h2D_DispTable", &h2D_DispTable );
            fData->SetBranchAddress( "h2D_DispPhiTable", &h2D_DispPhiTable );
            fData->SetBranchAddress( "h2D_DispMissTable", &h2D_DispMissTable );
            fData->SetBranchAddress( "h2D_DispTableN", &h2D_DispTableN );
            
            fData->SetBranchAddress( "h3D_DispTable", &h3D_DispTable );
            fData->SetBranchAddress( "h3D_DispPhiTable", &h3D_DispPhiTable );
            fData->SetBranchAddress( "h3D_DispMissTable", &h3D_DispMissTable );
            fData->SetBranchAddress( "h3D_DispTableN", &h3D_DispTableN );
            
            // fill data vectors
            f_ze.clear();
            f_az_min.clear();
            f_az_max.clear();
            f_woff.clear();
            f_noise.clear();
            fTreeEntry_Map.clear();
            for( int i = 0; i < fData->GetEntries(); i++ )
            {
                fData->GetEntry( i );
                
                int i_found_ze = -99;
                for( unsigned int j = 0; j < f_ze.size(); j++ )
                {
                    if( TMath::Abs( f_ze[j] - ze ) < 1.e-2 )
                    {
                        i_found_ze = j;
                    }
                }
                if( i_found_ze < 0 )
                {
                    f_ze.push_back( ze );
                    i_found_ze = ( int )f_ze.size() - 1;
                }
                
                int i_found_az = -99;
                for( unsigned int j = 0; j < f_az_min.size(); j++ )
                {
                    if( TMath::Abs( f_az_min[j] - az_min ) < 1.e-2 )
                    {
                        i_found_az = j;
                    }
                }
                if( i_found_az < 0 )
                {
                    f_az_min.push_back( az_min );
                    f_az_max.push_back( az_max );
                    i_found_az = ( int )f_az_min.size() - 1;
                }
                int i_found_woff = -99;
                for( unsigned int j = 0; j < f_woff.size(); j++ )
                {
                    if( TMath::Abs( f_woff[j] - woff ) < 1.e-2 )
                    {
                        i_found_woff = j;
                    }
                }
                if( i_found_woff < 0 )
                {
                    f_woff.push_back( woff );
                    i_found_woff = ( int )f_woff.size() - 1;
                }
                int i_found_noise = -99;
                for( unsigned int j = 0; j < f_noise.size(); j++ )
                {
                    if( TMath::Abs( f_noise[j] - pedvar ) < 5.e-1 )
                    {
                        i_found_noise = j;    // is 0.5 too large??
                    }
                }
                if( i_found_noise < 0 )
                {
                    f_noise.push_back( pedvar );
                    i_found_noise = f_noise.size() - 1;
                }
                
                unsigned int i_ID = getTreeEntryID( i_found_ze, i_found_az, i_found_woff, i_found_noise );
                fTreeEntry_Map[i_ID] = i;
            }
        }
    }
    return true;
}

unsigned int VDispTableReader::getTreeEntryID( int i_found_ze, int i_found_az, int i_found_woff, int i_found_noise )
{
    unsigned int i_ID = ( unsigned int )i_found_noise;
    i_ID += ( unsigned int )i_found_woff * 100;
    i_ID += ( unsigned int )i_found_az * 100 * 100;
    i_ID += ( unsigned int )i_found_ze * 100 * 100 * 100;
    
    return i_ID;
}




void VDispTableReader::reset()
{
    if( h2D_DispTable )
    {
        h2D_DispTable->Reset();
    }
    if( h2D_DispPhiTable )
    {
        h2D_DispPhiTable->Reset();
    }
    if( h2D_DispMissTable )
    {
        h2D_DispMissTable->Reset();
    }
    if( h2D_DispTableN )
    {
        h2D_DispTableN->Reset();
    }
}


bool VDispTableReader::fill( float i_ze, unsigned int i_az, float i_az_min, float i_az_max, float i_woff, float i_meanPedvars, TH2* iH2D, TH2* iH2DN, TH2* iH2DPhi, TH2* iH2DMiss, TH3* iH3D, TH3* iH3DN, TH3* iH3DPhi, TH3* iH3DMiss )
{
    if( !iH2D || !iH2DN || !iH2DPhi || !iH2DMiss )
    {
        return false;
    }
    if( !iH3D || !iH3DN || !iH3DPhi || !iH3DMiss )
    {
        return false;
    }
    
    if( !h2D_DispTable )
    {
        return false;
    }
    if( !h2D_DispPhiTable )
    {
        return false;
    }
    if( !h2D_DispTableN )
    {
        return false;
    }
    if( !h2D_DispMissTable )
    {
        return false;
    }
    if( !h3D_DispTable )
    {
        return false;
    }
    if( !h3D_DispPhiTable )
    {
        return false;
    }
    if( !h3D_DispTableN )
    {
        return false;
    }
    if( !h3D_DispMissTable )
    {
        return false;
    }
    
    if( !fData )
    {
        return false;
    }
    
    ze = i_ze;
    az_bin = i_az;
    az_min = i_az_min;
    az_max = i_az_max;
    woff = i_woff;
    pedvar = i_meanPedvars;
    
    // fill 2D histograms
    for( int x = 0; x <= iH2D->GetNbinsX(); x++ )
    {
        for( int y = 0; y <= iH2D->GetNbinsY(); y++ )
        {
            if( iH2D->GetBinContent( x, y ) > 0. )
            {
                h2D_DispTable->SetBinContent( x, y, iH2D->GetBinContent( x, y ) );
                h2D_DispTable->SetBinError( x, y, iH2D->GetBinError( x, y ) );
                h2D_DispPhiTable->SetBinContent( x, y, iH2DPhi->GetBinContent( x, y ) );
                h2D_DispPhiTable->SetBinError( x, y, iH2DPhi->GetBinError( x, y ) );
                h2D_DispMissTable->SetBinContent( x, y, iH2DMiss->GetBinContent( x, y ) );
                h2D_DispMissTable->SetBinError( x, y, iH2DMiss->GetBinError( x, y ) );
                h2D_DispTableN->SetBinContent( x, y, iH2DN->GetBinContent( x, y ) );
            }
        }
    }
    // fill 3D histograms
    for( int x = 0; x <= iH3D->GetNbinsX(); x++ )
    {
        for( int y = 0; y <= iH3D->GetNbinsY(); y++ )
        {
            for( int z = 0; z <= iH3D->GetNbinsZ(); z++ )
            {
                if( iH3D->GetBinContent( x, y, z ) > 0. )
                {
                    h3D_DispTable->SetBinContent( x, y, z, iH3D->GetBinContent( x, y, z ) );
                    h3D_DispTable->SetBinError( x, y, z, iH3D->GetBinError( x, y, z ) );
                    h3D_DispPhiTable->SetBinContent( x, y, z, iH3DPhi->GetBinContent( x, y, z ) );
                    h3D_DispPhiTable->SetBinError( x, y, z, iH3DPhi->GetBinError( x, y, z ) );
                    h3D_DispMissTable->SetBinContent( x, y, z, iH3DMiss->GetBinContent( x, y, z ) );
                    h3D_DispMissTable->SetBinError( x, y, z, iH3DMiss->GetBinError( x, y, z ) );
                    h3D_DispTableN->SetBinContent( x, y, z, iH3DN->GetBinContent( x, y, z ) );
                }
            }
        }
    }
    
    // fill data tree
    fData->Fill();
    
    // check that vectors are filled correctly
    if( fData->GetEntries() - 1 == ( int )( f_ze.size() ) )
    {
        f_ze.push_back( ze );
        f_az_min.push_back( az_min );
        f_az_max.push_back( az_max );
        f_woff.push_back( woff );
        f_noise.push_back( pedvar );
    }
    
    return true;
}


void VDispTableReader::setDataVectors( vector< string > i_ze, vector< float > i_az_min, vector< float > i_az_max, vector< string > i_woff, vector< string > i_noise )
{
    f_ze.clear();
    for( unsigned int i = 0; i < i_ze.size(); i++ )
    {
        f_ze.push_back( atof( i_ze[i].c_str() ) );
    }
    f_az_min = i_az_min;
    f_az_max = i_az_max;
    f_woff.clear();
    for( unsigned int i = 0; i < i_woff.size(); i++ )
    {
        f_woff.push_back( atof( i_woff[i].c_str() ) );
    }
    f_noise.clear();
    for( unsigned int i = 0; i < i_noise.size(); i++ )
    {
        f_noise.push_back( atof( i_noise[i].c_str() ) );
    }
}


void VDispTableReader::terminate()
{

}


void VDispTableReader::print( bool bDetailed )
{
    if( bDetailed )
    {
        cout << endl;
    }
    else
    {
        cout << "\t";
    }
    cout << "disp table (" << GetName() << ")" << endl;
    cout << "\t tables available for " << f_ze.size() << " zenith angle(s), ";
    cout << f_az_min.size() << " azimuth bin(s), ";
    cout << f_woff.size() << " wobble offset(s), ";
    cout << f_noise.size() << " noise value(s)";
    cout << endl;
    if( bDetailed )
    {
        cout << "============================================================================" << endl;
        cout << "entries in tree: " << fData->GetEntries() << endl;
        for( int i = 0; i < fData->GetEntries(); i++ )
        {
            fData->GetEntry( i );
            
            cout << i << ", ze " << ze << ", az_bin " << az_bin << ", az_min " << az_min;
            cout << ",az_max " << az_max << ",woff " << woff << ",pedvar " << pedvar << endl;
        }
        for( unsigned int i = 0; i < f_ze.size(); i++ )
        {
            cout << "ze vector " << f_ze[i] << endl;
        }
        for( unsigned int i = 0; i < f_az_min.size(); i++ )
        {
            cout << "az vector " << f_az_min[i] << "\t" << f_az_max[i] << endl;
        }
        for( unsigned int i = 0; i < f_woff.size(); i++ )
        {
            cout << "woff vector " << f_woff[i] << endl;
        }
        for( unsigned int i = 0; i < f_noise.size(); i++ )
        {
            cout << "noise vector " << f_noise[i] << endl;
        }
    }
}


unsigned int VDispTableReader::getAzBin( float az )
{
    // be sure that az is bin interval [-180., 180.]
    if( az > 180. )
    {
        az -= 360.;
    }
    // expect first bin to be of type [135,-135]
    // otherwise like [-45,45]
    for( unsigned int i = 0; i < f_az_min.size(); i++ )
    {
        if( i == 0 && ( az > f_az_min[0] || az < f_az_max[0] ) )
        {
            return 0;
        }
        else if( az > f_az_min[i] && az < f_az_max[i] )
        {
            return i;
        }
    }
    return 0;
}

int VDispTableReader::getTreeEntryFinder( unsigned int iID )
{
    if( fTreeEntry_Map.find( iID ) != fTreeEntry_Map.end() )
    {
        return fTreeEntry_Map[iID];
    }
    
    cout << "VDispTableReader::getTreeEntryFinder: warning: no entry found for ID " << iID << endl;
    
    for( map<unsigned int, int>::iterator i_iter = fTreeEntry_Map.begin(); i_iter != fTreeEntry_Map.end(); i_iter++ )
    {
        cout << ( *i_iter ).first << "\t" << ( *i_iter ).second << endl;
    }
    
    return 0;
}


void VDispTableReader::getIndexBoundary( unsigned int* iup, unsigned int* ilow, vector< float >& iV, float x )
{
    if( iV.size() == 0 )
    {
        *iup = 0;
        *ilow = 0;
        return;
    }
    
    if( x <= iV[0] )
    {
        *iup = *ilow = 0;
    }
    else if( x >= iV[iV.size() - 1] )
    {
        *iup = *ilow = iV.size() - 1;
    }
    else
    {
        for( unsigned int i = 0; i < iV.size(); i++ )
        {
            if( x > iV[i] )
            {
                *ilow = ( unsigned int )i;
                *iup = ( unsigned int )( i + 1 );
            }
        }
    }
}

float VDispTableReader::getLowerZe( float iZe )
{
    unsigned int i_ze_bin_up  = 0;
    unsigned int i_ze_bin_low = 0;
    getIndexBoundary( &i_ze_bin_up, &i_ze_bin_low, f_ze, iZe );
    
    return f_ze[i_ze_bin_low];
}

float VDispTableReader::getUpperZe( float iZe )
{
    unsigned int i_ze_bin_up  = 0;
    unsigned int i_ze_bin_low = 0;
    getIndexBoundary( &i_ze_bin_up, &i_ze_bin_low, f_ze, iZe );
    
    return f_ze[i_ze_bin_up];
}

/*!
    return tree entry number for corresponding values

    iZe_Inter = 0:   find entry number for closest Ze
    iZe_Inter = 1:   find entry number for closest smaller Ze
    iZe_Inter = 2:   find entry number for closest larger Ze
*/
int VDispTableReader::getTreeEntryFinder( float iZe, float iAz, float iWoff, float iPedvar, int iZe_Inter )
{
    unsigned int i_ze_bin = 0;
    unsigned int i_az_bin = 0;
    unsigned int i_woff_bin = 0;
    unsigned int i_noise_bin = 0;
    
    float i_diff = 1.e9;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // zenith angle
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    // get closest ze value
    if( iZe_Inter == 0 )
    {
        i_ze_bin = 0;
        for( unsigned int i = 0; i < f_ze.size(); i++ )
        {
            if( TMath::Abs( iZe - f_ze[i] ) < i_diff )
            {
                i_ze_bin = i;
                i_diff = TMath::Abs( iZe - f_ze[i] );
            }
        }
    }
    // get closest lower or upper value
    else
    {
        unsigned int i_ze_bin_up  = 0;
        unsigned int i_ze_bin_low = 0;
        getIndexBoundary( &i_ze_bin_up, &i_ze_bin_low, f_ze, iZe );
        if( iZe_Inter == 1 )
        {
            i_ze_bin = i_ze_bin_low;
        }
        else if( iZe_Inter == 2 )
        {
            i_ze_bin = i_ze_bin_up;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // azimuth angle
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    i_az_bin = getAzBin( iAz );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // wobble offset (ignored?)
    /////////////////////////////////////////////////////////////////////////////////////////////////
    i_woff_bin = 0;
    i_diff = 1.e9;
    for( unsigned int i = 0; i < f_woff.size(); i++ )
    {
        if( TMath::Abs( iWoff - f_woff[i] ) < i_diff )
        {
            i_woff_bin = i;
            i_diff = TMath::Abs( iWoff - f_woff[i] );
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // noise value (get closest value)
    /////////////////////////////////////////////////////////////////////////////////////////////////
    i_diff = 1.e9;
    i_noise_bin = 0;
    for( unsigned int i = 0; i < f_noise.size(); i++ )
    {
        if( TMath::Abs( iPedvar - f_noise[i] ) < i_diff )
        {
            i_noise_bin = i;
            i_diff = TMath::Abs( iPedvar - f_noise[i] );
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // now get tree entry number
    unsigned int iEntryID = getTreeEntryID( i_ze_bin, i_az_bin, i_woff_bin, i_noise_bin );
    
    return getTreeEntryFinder( iEntryID );
}

/*!
    scale parameter to get a more uniform distribution

    (after M.Beilicke 2010; VERITAS internal note)
*/
double VDispTableReader::scaleParameter( double iPara, double iScale )
{
    if( iScale <= 0. )
    {
        return 0.;
    }
    
    return 1. / ( 1. + 0.5 * TMath::Exp( 1 - iPara / iScale ) );
}

