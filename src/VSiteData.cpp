#include "VSiteData.h"

VSiteData::VSiteData()
{
    setDebug();
    reset();
}

void VSiteData::reset()
{
    fSiteName = "";
    fHemisphere = "S";
    fSite_asl = 0.;
    fSite_B_N = 0.;
    fSite_B_S = 0.;
    fSite_B_dB = 0.;
    fSiteRequirementID = 0;
    fReferenceSiteName = "";
    
    fObservationTime_s.clear();
    fCameraOffset_deg.clear();
    fArray.clear();
    
    fSiteFileName.clear();
    fSiteFile_exists.clear();
    fSiteFile_Emin.clear();
    fSiteFile_Emax.clear();
    
    fPlottingColor.clear();
    fPlottingLineStyle.clear();
    fPlottingFillStyle.clear();
    fPlottingMarkerStyle.clear();
    fLegend.clear();
    
    fGraphSensitivity.clear();
    fGraphSensitivityInterPolated = 0;
}

VSiteData::~VSiteData()
{
    /*    for( unsigned int i = 0; i < fGraphSensitivity.size(); i++ )
        {
            if( fGraphSensitivity[i] )
            {
                delete fGraphSensitivity[i];
            }
        }
        if( fGraphSensitivityInterPolated )
        {
            delete fGraphSensitivityInterPolated;
        } */
}

/*

   expect all vectors to be of same size

*/
bool VSiteData::checkIntegrity()
{
    unsigned int iS = fSiteFileName.size();
    
    if( iS != fObservationTime_s.size() )
    {
        return false;
    }
    if( iS != fCameraOffset_deg.size() )
    {
        return false;
    }
    if( iS != fArray.size() )
    {
        return false;
    }
    if( iS != fSiteFile_exists.size() )
    {
        return false;
    }
    if( iS != fSiteFile_Emin.size() )
    {
        return false;
    }
    if( iS != fSiteFile_Emax.size() )
    {
        return false;
    }
    if( iS != fPlottingColor.size() )
    {
        return false;
    }
    if( iS != fPlottingLineStyle.size() )
    {
        return false;
    }
    if( iS != fPlottingFillStyle.size() )
    {
        return false;
    }
    if( iS != fPlottingMarkerStyle.size() )
    {
        return false;
    }
    if( iS != fLegend.size() )
    {
        return false;
    }
    if( iS != fGraphSensitivity.size() )
    {
        return false;
    }
    
    return true;
}

/*

   print data set

*/
void VSiteData::print()
{
    cout << fSiteName << "(" << fHemisphere << ") at ";
    cout << fSite_asl << " m, N: " << fSite_B_N << " muG, S: " << fSite_B_S << " muG" << " (reference site " << fReferenceSiteName << ")" << endl;
    for( unsigned int i = 0; i < fSiteFileName.size(); i++ )
    {
        cout << "\t " << fObservationTime_s[i] / 3600. << "h, array " << fArray[i] << ", offset " << fCameraOffset_deg[i] <<  " deg";
        cout << "energy range [" << fSiteFile_Emin[i] << ", " << fSiteFile_Emax[i] << "] TeV: " << endl;
        cout << "\t " << fSiteFileName[i] << endl;
        cout << "\t (color " << fPlottingColor[i] << ", line style " << fPlottingLineStyle[i] << ")" << endl;
    }
}

/*

    read list of data sets and add them to the fData vector

*/
bool VSiteData::addDataSet( string iDataList, unsigned int iSiteCounter, string iDirectionString )
{
    ifstream is;
    is.open( iDataList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VSiteData::addDataSets error opening data set txt file: " << iDataList << endl;
        return false;
    }
    string is_line;
    string iTemp = "";
    string iTempSite = "";
    
    ///////////////////////
    // get list of sites
    // (sites with the same name are combined)
    vector< string > iListOfSites;
    unsigned int z = 0;
    while( getline( is, is_line ) )
    {
        istringstream is_stream( is_line );
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            // ignore lines with '#' in the beginning
            if( iTemp.size() > 0 && iTemp.substr( 0, 1 ) == "#" )
            {
                continue;
            }
            // short version of list files indicated by "!" in first column
            if( iTemp.size() > 0 && iTemp.substr( 0, 1 ) == "!" )
            {
                ostringstream iSN;
                iSN << "D" << z;
                iTemp = iSN.str();
                z++;
            }
            bool bFound = false;
            for( unsigned int i = 0; i < iListOfSites.size(); i++ )
            {
                if( iTemp == iListOfSites[i] )
                {
                    bFound = true;
                }
            }
            if( !bFound )
            {
                iListOfSites.push_back( iTemp );
            }
        }
    }
    // no more sites
    if( iSiteCounter >= iListOfSites.size() )
    {
        return false;
    }
    cout << "reading from " << iDataList << " : " << iSiteCounter << "=" << iListOfSites[iSiteCounter] << endl;
    
    // reset stream (??!?!?!)
    is.close();
    is.open( iDataList.c_str(), ifstream::in );
    
    // reset counter
    z = 0;
    bool iShortDataListVersion = false;
    while( getline( is, is_line ) )
    {
        istringstream is_stream( is_line );
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            // ignore lines with '#' in the beginning
            if( iTemp == "#" )
            {
                continue;
            }
            // short version of list files indicated by "!" in first column
            if( iTemp == "!" )
            {
                iShortDataListVersion = true;
                ostringstream iSN;
                iSN << "D" << z;
                iTemp = iSN.str();
            }
            if( iTemp != iListOfSites[iSiteCounter] )
            {
                z++;
                continue;
            }
            fSiteName = iTemp;
        }
        // not need for most analysis
        if( !iShortDataListVersion )
        {
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fHemisphere;
                if( fHemisphere != "S" && fHemisphere != "N" )
                {
                    cout << "Warning: unknown hemisphere: " << fHemisphere << endl;
                }
                if( fHemisphere == "S" )
                {
                    fSiteRequirementID = 0;
                }
                else if( fHemisphere == "N" )
                {
                    fSiteRequirementID = 3;
                }
            }
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fReferenceSiteName;
            }
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fSite_asl;
            }
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fSite_B_N;
            }
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fSite_B_S;
            }
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fSite_B_dB;
            }
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            // add direction string
            if( iTemp.find( "NIM" ) != string::npos && iDirectionString.size() > 0 )
            {
                iTemp.insert( iTemp.find( "NIM" ), iDirectionString );
            }
            fSiteFileName.push_back( iTemp );
            fSiteFile_exists.push_back( false );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fArray.push_back( iTemp );
        }
        else
        {
            fArray.push_back( "" );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fObservationTime_s.push_back( atof( iTemp.c_str() ) );
        }
        else
        {
            fObservationTime_s.push_back( 0. );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fCameraOffset_deg.push_back( atof( iTemp.c_str() ) );
        }
        else
        {
            fCameraOffset_deg.push_back( 0. );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fSiteFile_Emin.push_back( atof( iTemp.c_str() ) );
        }
        else
        {
            fSiteFile_Emin.push_back( 1.e-5 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fSiteFile_Emax.push_back( atof( iTemp.c_str() ) );
        }
        else
        {
            fSiteFile_Emax.push_back( 1.e5 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fPlottingColor.push_back( atoi( iTemp.c_str() ) );
        }
        else
        {
            fPlottingColor.push_back( 1 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fPlottingLineStyle.push_back( atoi( iTemp.c_str() ) );
        }
        else
        {
            fPlottingLineStyle.push_back( 1 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fPlottingFillStyle.push_back( atoi( iTemp.c_str() ) );
        }
        else
        {
            fPlottingFillStyle.push_back( 1 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            is_stream >> iTemp;
            fPlottingMarkerStyle.push_back( atoi( iTemp.c_str() ) );
        }
        else
        {
            fPlottingMarkerStyle.push_back( 21 );
        }
        if( !(is_stream>>std::ws).eof() )
        {
            fLegend.push_back( is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ) );
        }
        else
        {
            fLegend.push_back( "" );
        }
        // fill up remaining vectors
        fGraphSensitivity.push_back( 0 );
        fGraphSensitivityInterPolated = 0;
        z++;
    }
    is.close();
    
    if( fDebug )
    {
        cout << "VSiteData::addDataSet: integrity: " << checkIntegrity() << endl;
    }
    
    // add observing time to legend
    char hname[200];
    for( unsigned int i = 0; i < fObservationTime_s.size(); i++ )
    {
        if( fLegend[i] != "NO_LEGEND" )
        {
            if( fObservationTime_s[i] < 18. )
            {
                sprintf( hname, "%s, %s", fLegend[i].c_str(), fArray[i].c_str() );
            }
            else
            {
                // don't plot array layout name into legend
                sprintf( hname, "%s, %s", fLegend[i].c_str(), fArray[i].c_str() );
                sprintf( hname, "%s", fLegend[i].c_str() );
            }
            fLegend[i] = hname;
        }
    }
    
    // set data files and data directories depending on the analysis type
    // (some of the following statements refer to outdated
    //  directory structures)
    for( unsigned int i = 0; i < fSiteFileName.size(); i++ )
    {
        stringstream sobs( stringstream::in | stringstream::out );
        sobs << fixed << TMath::Nint( fObservationTime_s[i] );
        if( fSiteFileName[i].find( "/" ) != string::npos )
        {
            iTemp = "data/" + fSiteFileName[i] + "." + fArray[i] + ".";
            iTemp += sobs.str();
            iTemp += "s.root";
            fSiteFileName[i] = iTemp;
        }
        else if( fSiteFileName[i].find( "DESY" ) != string::npos )
        {
            iTemp = "data/DESY/" + fSiteFileName[i] + "." + fArray[i] + ".";
            iTemp += sobs.str();
            iTemp += "s.root";
            fSiteFileName[i] = iTemp;
        }
        else if( fSiteFileName[i].find( "IFAE" ) != string::npos
                 && fSiteFileName[i].find( "20130901" ) != string::npos )
        {
            iTemp = "data/IFAE_Sept2013/" + fSiteFileName[i] + ".root";
            fSiteFileName[i] = iTemp;
        }
        else if( fSiteFileName[i].find( "IFAE" ) != string::npos && fSiteFileName[i].find( "19122013" ) != string::npos )
        {
            iTemp = "data/IFAE_Dec2013/" + fSiteFileName[i] + ".root";
            fSiteFileName[i] = iTemp;
        }
        else if( fSiteFileName[i].find( "IFAE" ) != string::npos && fSiteFileName[i].find( "20141028" ) != string::npos )
        {
            iTemp = "data/IFAE_Dec2014/" + fSiteFileName[i] + ".root";
            fSiteFileName[i] = iTemp;
        }
        else if( fSiteFileName[i].find( "Fermi" ) != string::npos )
        {
            fSiteFileName[i] = "data/Fermi/Fermi_approx.root";
        }
        else
        {
            fSiteFileName[i] += ".root";
        }
    }
    
    // check that file exists
    for( unsigned int i = 0; i < fSiteFileName.size(); i++ )
    {
        TFile iF( fSiteFileName[i].c_str() );
        if( iF.IsZombie() )
        {
            fSiteFile_exists[i] = false;
        }
        
        fSiteFile_exists[i] = true;
        iF.Close();
    }
    
    return true;
}

/*
 *
 *
 * */
TGraphAsymmErrors* VSiteData::getCombinedSensitivityGraph( bool iInterpolate, string iDirectionString, bool iIntegratedSensitivity )
{
    TGraphAsymmErrors* iGraphSensitivity = new TGraphAsymmErrors( 1 );
    if( fSiteFileName.size() > 0 )
    {
        iGraphSensitivity->SetLineColor( fPlottingColor[0] );
        iGraphSensitivity->SetMarkerColor( fPlottingColor[0] );
        iGraphSensitivity->SetFillColor( fPlottingColor[0] );
        iGraphSensitivity->SetLineStyle( fPlottingLineStyle[0] );
    }
    
    unsigned int z = 0;
    for( unsigned int f = 0; f < fSiteFileName.size(); f++ )
    {
        string iFileName = fSiteFileName[f].c_str();
        if( iDirectionString.size() > 0 )
        {
            iFileName.insert( iFileName.find( "NIM" ), iDirectionString );
        }
        TFile* iFile = new TFile( iFileName.c_str() );
        if( iFile->IsZombie() )
        {
            return 0;
        }
        TH1F* h = 0;
        // offaxis angle < 1.e-2 indicate usage of
        // on-axis histograms
        if( fCameraOffset_deg[f] < 1.e-2 )
        {
            if( iIntegratedSensitivity )
            {
                h = ( TH1F* )iFile->Get( "IntSens" );
            }
            else
            {
                h = ( TH1F* )iFile->Get( "DiffSens" );
            }
        }
        else
        {
            TH2F* h2 = 0;
            if( iIntegratedSensitivity )
            {
                h2 = ( TH2F* )iFile->Get( "IntSens_offaxis" );
            }
            else
            {
                h2 = ( TH2F* )iFile->Get( "DiffSens_offaxis" );
            }
            if( h2 )
            {
                char hname[200];
                if( iIntegratedSensitivity )
                {
                    sprintf( hname, "IntSens_%d", ( int )fCameraOffset_deg[f] * 100 );
                }
                else
                {
                    sprintf( hname, "DiffSens_%d", ( int )fCameraOffset_deg[f] * 100 );
                }
                h = ( TH1F* )h2->ProjectionX( hname, h2->GetYaxis()->FindBin( fCameraOffset_deg[f] ),
                                              h2->GetYaxis()->FindBin( fCameraOffset_deg[f] ) );
            }
        }
        if( !h )
        {
            return 0;
        }
        
        if( fDebug )
        {
            cout << "VSiteData::getCombinedSensitivityGraph: reading for site " << fSiteName << "[";
            cout << fSiteFile_Emin[f] << ", " << fSiteFile_Emax[f] << "] from file: " << endl;
            cout << "\t" << fSiteFileName[f].c_str() << endl;
        }
        
        for( int i = 1; i <= h->GetNbinsX(); i++ )
        {
            if( h->GetXaxis()->GetBinCenter( i ) > log10( fSiteFile_Emin[f] ) && h->GetXaxis()->GetBinCenter( i ) <= log10( fSiteFile_Emax[f] ) )
            {
                if( h->GetBinContent( i ) > 0. )
                {
                    iGraphSensitivity->SetPoint( z, h->GetXaxis()->GetBinCenter( i ), h->GetBinContent( i ) );
                    if( h->GetBinError( i ) < h->GetBinContent( i ) )
                    {
                        iGraphSensitivity->SetPointError( z, h->GetXaxis()->GetBinWidth( i ) / 2., h->GetXaxis()->GetBinWidth( i ) / 2.,
                                                          h->GetBinError( i ), h->GetBinError( i ) );
                    }
                    else
                    {
                        iGraphSensitivity->SetPointError( z, h->GetXaxis()->GetBinWidth( i ) / 2., h->GetXaxis()->GetBinWidth( i ) / 2., 0., 0. );
                    }
                    z++;
                    continue;
                }
                else if( iInterpolate )
                {
                    if( i > 1 && i < h->GetNbinsX() )
                    {
                        if( h->GetBinContent( i - 1 ) > 0. && h->GetBinContent( i + 1 ) > 0. )
                        {
                            iGraphSensitivity->SetPoint( z, h->GetXaxis()->GetBinCenter( i ),
                                                         0.5 * ( h->GetBinContent( i - 1 ) + h->GetBinContent( i + 1 ) ) );
                            iGraphSensitivity->SetPointError( z, h->GetXaxis()->GetBinWidth( i ) / 2., h->GetXaxis()->GetBinWidth( i ) / 2.,
                                                              0.5 * ( h->GetBinError( i - 1 ) + h->GetBinError( i + 1 ) ),
                                                              0.5 * ( h->GetBinError( i - 1 ) + h->GetBinError( i + 1 ) ) );
                            cout << "\t (" << fSiteName << ") interpolate point at energy " << h->GetXaxis()->GetBinCenter( i ) << ": ";
                            cout << h->GetBinError( i - 1 ) << "\t" << h->GetBinError( i + 1 ) << endl;
                            cout << endl;
                            z++;
                            continue;
                        }
                    }
                }
            }
        }
        iFile->Close();
    }
    
    return iGraphSensitivity;
}


