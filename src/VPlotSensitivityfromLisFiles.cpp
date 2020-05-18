/*! \file VPlotSensitivityfromLisFiles
    \brief plot sensitivity parameters from lis files

*/

#include "VPlotSensitivityfromLisFiles.h"

VLisFileData::VLisFileData()
{
    fID = 0;
    fName = "";
    fFileName = "";
}

//////////////////////////////////////////////////////////////////////////////////

VPlotSensitivityfromLisFiles::VPlotSensitivityfromLisFiles()
{
    fDebug = true;
    
    // set variable names and min/max
    fVarName.push_back( "E1" );
    fVarMin[fVarName.back()] =  -2.0;
    fVarMax[fVarName.back()] =   2.5;
    fVarName.push_back( "E2" );
    fVarMin[fVarName.back()] =  -2.0;
    fVarMax[fVarName.back()] =   2.5;
    fVarName.push_back( "Emean" );
    fVarMin[fVarName.back()] =  -2.0;
    fVarMax[fVarName.back()] =   2.5;
    
    fVarName.push_back( "NTel" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  10.;
    
    fVarName.push_back( "NPix" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  20.;
    
    fVarName.push_back( "Amp" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.3;
    
    fVarName.push_back( "BorderTh" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  20.;
    
    fVarName.push_back( "ImageTh" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  20.;
    
    fVarName.push_back( "Eres" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.;
    
    fVarName.push_back( "Hist" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.6;
    
    fVarName.push_back( "AngRes68" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  0.5;
    
    fVarName.push_back( "AngRes80" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  0.5;
    
    fVarName.push_back( "Tr" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.e2;
    
    fVarName.push_back( "S_diff_D" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.e2;
    
    fVarName.push_back( "S_diff" );
    fVarMin[fVarName.back()] =  1.e-15;
    fVarMax[fVarName.back()] =  1.e-3;
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e1;
    
    fVarName.push_back( "S_int_D" );
    fVarMin[fVarName.back()] =  0.;
    fVarMax[fVarName.back()] =  1.e2;
    
    fVarName.push_back( "S_int" );
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e1;
    
    fVarName.push_back( "N_Gamma" );
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e5;
    
    fVarName.push_back( "N_Proton" );
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e6;
    
    fVarName.push_back( "N_Electron" );
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e6;
    
    fVarName.push_back( "N_Nuclei" );
    fVarMin[fVarName.back()] =  1.e-3;
    fVarMax[fVarName.back()] =  1.e6;
}

bool VPlotSensitivityfromLisFiles::addLisFile( string iFile, string iCut )
{
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VPlotSensitivityfromLisFiles::addLisFile: error opening input file: " << iFile << endl;
        return false;
    }
    
    // read file and fill data vector
    string is_line;
    string temp;
    bool bMean = false;
    bool bContinue = false;
    while( getline( is, is_line ) )
    {
        if( is_line.substr( 0, 4 ) == "##+ " )
        {
            if( iCut.size() > 0 )
            {
                if( iCut != is_line.substr( is_line.find( "passing" ), is_line.find( "cuts" ) - is_line.find( "passing" ) - 1 ) )
                {
                    bContinue = false;
                    continue;
                }
                else
                {
                    bContinue = true;
                }
            }
            fData.push_back( new VLisFileData() );
            fData.back()->fID = fData.size();
            fData.back()->fName = is_line.substr( is_line.find( "passing" ), is_line.find( "cuts" ) - is_line.find( "passing" ) - 1 );
            fData.back()->fFileName = iFile;
            cout << "adding new data set " << fData.size() << endl;
            continue;
        }
        else if( !bContinue )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) == "#" )
        {
            if( is_line.find( "<lg E>" ) < is_line.size() )
            {
                bMean = true;
            }
            continue;
        }
        
        if( fData.size() > 0 )
        {
            istringstream is_stream( is_line );
            
            for( unsigned int i = 0; i < fVarName.size(); i++ )
            {
                if( i == 2 && !bMean )
                {
                    fData.back()->fVar[fVarName[i]].push_back( log10( 0.5 * ( TMath::Power( 10., fData.back()->fVar["E1"].back() ) + TMath::Power( 10., fData.back()->fVar["E2"].back() ) ) ) );
                    continue;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fData.back()->fVar[fVarName[i]].push_back( atof( temp.c_str() ) );
                }
            }
        }
    }
    is.close();
    
    return true;
}

void VPlotSensitivityfromLisFiles::listDataSets()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        cout << i << "\t" << fData[i]->fID << "\t" << fData[i]->fName << "\t" << fData[i]->fFileName << endl;
    }
}

bool VPlotSensitivityfromLisFiles::printDataSet( unsigned int iID )
{
    // get ID index
    iID = getID_Index( iID );
    if( iID == 9999 )
    {
        cout << "VPlotSensitivityfromLisFiles::printDataSet error: data set ID not found: " << iID << "\t" << fData.size() << endl;
        return 0;
    }
    
    cout << fData[iID]->fID << "\t" << fData[iID]->fFileName << endl;
    map< string, vector< double > >::iterator fVarIter;
    for( fVarIter = fData[iID]->fVar.begin(); fVarIter != fData[iID]->fVar.end(); ++fVarIter )
    {
        cout << ( *fVarIter ).first << "\t" << ( *fVarIter ).second.size();
        for( unsigned int i = 0; i < ( *fVarIter ).second.size(); i++ )
        {
            cout << "\t" << ( *fVarIter ).second[i];
        }
        cout << endl;
    }
    return true;
}

bool VPlotSensitivityfromLisFiles::checkVarName( string V )
{
    for( unsigned int i = 0; i < fVarName.size(); i++ )
    {
        if( V == fVarName[i] )
        {
            return true;
        }
    }
    return false;
}

unsigned int VPlotSensitivityfromLisFiles::getID_Index( unsigned int iID )
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i]->fID == iID )
        {
            return i;
        }
    }
    
    return 9999;
}

TCanvas* VPlotSensitivityfromLisFiles::plot( string iVName, unsigned int iID, TCanvas* c, Style_t iLineStyle, Style_t iMarkerStyle )
{
    // get ID index
    iID = getID_Index( iID );
    if( iID == 9999 )
    {
        cout << "VPlotSensitivityfromLisFiles::plot error: data set ID not found: " << iID << "\t" << fData.size() << endl;
        return 0;
    }
    
    if( !checkVarName( iVName ) )
    {
        cout << "VPlotSensitivityfromLisFiles::plot error: variable not found: " << iVName << endl;
        return 0;
    }
    
    char hname[800];
    char htitle[800];
    
    if( !c )
    {
        sprintf( hname, "c_%s_%d", iVName.c_str(), iID );
        sprintf( htitle, "%s (%d)", iVName.c_str(), iID );
        c = new TCanvas( hname, htitle, 10, 10, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        
        sprintf( hname, "hnull_%s_%d", iVName.c_str(), iID );
        TH1D* hnull = new TH1D( hname, "", 100, fVarMin["E1"], fVarMax["E2"] );
        hnull->SetMaximum( fVarMax[iVName] );
        hnull->SetMinimum( fVarMin[iVName] );
        hnull->SetXTitle( "log_{10} energy [TeV]" );
        hnull->SetYTitle( iVName.c_str() );
        hnull->SetStats( 0 );
        hnull->Draw();
    }
    
    // fill graph
    
    TGraph* g = new TGraph( 1 );
    setGraphPlottingStyle( g, iID + 1, 2., iMarkerStyle, 1, 0, iLineStyle );
    
    for( unsigned int i = 0; i < fData[iID]->fVar[iVName].size(); i++ )
    {
        g->SetPoint( i, fData[iID]->fVar["Emean"][i], fData[iID]->fVar[iVName][i] );
    }
    
    g->Draw( "pc" );
    
    return c;
}

TCanvas* VPlotSensitivityfromLisFiles::plot_AllIDs( string iVName, Style_t iLineStyle, Style_t iMarkerStyle, TCanvas* c )
{
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    TCanvas* d = plot( iVName, fData[0]->fID, c, iLineStyle, iMarkerStyle );
    if( d )
    {
        for( unsigned int i = 1; i < fData.size(); i++ )
        {
            plot( iVName, fData[i]->fID, d, iLineStyle, iMarkerStyle );
        }
    }
    return d;
}

void VPlotSensitivityfromLisFiles::listVariableNames()
{
    for( unsigned int i = 0; i < fVarName.size(); i++ )
    {
        cout << fVarName[i] << endl;
    }
}

bool VPlotSensitivityfromLisFiles::removeDataSet( unsigned int iID )
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( iID == fData[i]->fID )
        {
            fData.erase( fData.begin() + i );
            return true;
        }
    }
    
    cout << "VPlotSensitivityfromLisFiles::removeDataSet error: data set ID not found: " << iID << "\t" << fData.size() << endl;
    return false;
}

/*
   apply cuts:

   variable v > 0.: remove all data sets with v == value
   vairable v < 0.: remove all data sets with v != value

*/
bool VPlotSensitivityfromLisFiles::applycuts( double amp, double NTel, double NPix )
{
    set< unsigned int > fRemoveID;
    // amplituden cut
    if( TMath::Abs( amp ) > 1.e-2 )
    {
        for( unsigned int i = 0; i < fData.size(); i++ )
        {
            if( fData[i]->fVar["Amp"].size() > 0 )
            {
                // remove all data sets with this amplitude
                if( amp > 0. && TMath::Abs( fData[i]->fVar["Amp"][0] - amp ) < 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
                // remove all data sets with a different amplitude
                if( amp < 0. && TMath::Abs( fData[i]->fVar["Amp"][0] + amp ) > 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
            }
        }
    }
    // pixel cut
    if( TMath::Abs( NPix ) > 1.e-2 )
    {
        for( unsigned int i = 0; i < fData.size(); i++ )
        {
            if( fData[i]->fVar["NPix"].size() > 0 )
            {
                if( NPix > 0. && TMath::Abs( fData[i]->fVar["NPix"][0] - NPix ) < 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
                if( NPix < 0. && TMath::Abs( fData[i]->fVar["NPix"][0] + NPix ) > 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
            }
        }
    }
    // ntel cut
    if( TMath::Abs( NTel ) > 1.e-2 )
    {
        for( unsigned int i = 0; i < fData.size(); i++ )
        {
            if( fData[i]->fVar["NTel"].size() > 0 )
            {
                if( NTel > 0. && TMath::Abs( fData[i]->fVar["NTel"][0] - NTel ) < 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
                if( NTel < 0. && TMath::Abs( fData[i]->fVar["NTel"][0] + NTel ) > 1.e-2 )
                {
                    fRemoveID.insert( fData[i]->fID );
                }
            }
        }
    }
    
    // now remove all marked data sets
    set< unsigned int >::iterator a;
    
    for( a = fRemoveID.begin(); a != fRemoveID.end(); ++a )
    {
        cout << "removing data set with ID " << *a << endl;
        removeDataSet( *a );
    }
    
    return true;
}

void VPlotSensitivityfromLisFiles::compare_cuts()
{
    cout << "ID \t ImageTh \t BorderTh \t Amp \t Npix \t NTel " << endl;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        cout << fData[i]->fID << "\t";
        if( fData[i]->fVar["ImageTh"].size() > 0 )
        {
            cout << fData[i]->fVar["ImageTh"][0];
        }
        cout << "\t";
        if( fData[i]->fVar["BorderTh"].size() > 0 )
        {
            cout << fData[i]->fVar["BorderTh"][0];
        }
        cout << "\t";
        if( fData[i]->fVar["Amp"].size() > 0 )
        {
            cout << fData[i]->fVar["Amp"][0];
        }
        cout << "\t";
        if( fData[i]->fVar["NPix"].size() > 0 )
        {
            cout << fData[i]->fVar["NPix"][0];
        }
        cout << "\t";
        if( fData[i]->fVar["NTel"].size() > 0 )
        {
            cout << fData[i]->fVar["NTel"][0];
        }
        cout << "\t";
        cout << endl;
    }
}
