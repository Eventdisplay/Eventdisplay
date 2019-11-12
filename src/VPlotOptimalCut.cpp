/*
    add variables by hand
*/

#include "VPlotOptimalCut.h"

VPlotOptimalCut::VPlotOptimalCut( string iFile )
{
    fFile = new TFile( iFile.c_str() );
    if( fFile->IsZombie() )
    {
        cout << "VPlotOptimalCut: error finding file " << iFile << endl;
        return;
    }
    fData = ( TTree* )fFile->Get( "topt" );
    if( !fData )
    {
        cout << "VPlotOptimalCut: error finding data tree in " << iFile << endl;
        return;
    }
    
    fListOfVariables.push_back( "MSCW" );
    fListOfVariables.push_back( "MSCL" );
    fListOfVariables.push_back( "EChi2" );
    fListOfVariables.push_back( "EmissionHeight" );
    fListOfVariables.push_back( "EmissionHeightChi2" );
    //    fListOfVariables.push_back( "SizeSecondMax" );
}

void VPlotOptimalCut::listVariables()
{
    for( unsigned int i = 0; i < fListOfVariables.size(); i++ )
    {
        cout << fListOfVariables[i] << endl;
    }
}

void VPlotOptimalCut::plotAll( int iSourceStrength, bool bPrint )
{
    int i_opt = findOptimalCut( iSourceStrength );
    for( unsigned int i = 0; i < fListOfVariables.size(); i++ )
    {
        plotHistograms( fListOfVariables[i], i_opt, bPrint );
    }
}

double VPlotOptimalCut::getCutValue( string iVariable, int iEntryNumber )
{
    if( fData && fData->GetBranchStatus( iVariable.c_str() ) )
    {
        double iValue = 0.;
        fData->SetBranchAddress( iVariable.c_str(), &iValue );
        if( fData->GetEntry( iEntryNumber ) > 0 )
        {
            fData->ResetBranchAddresses();
            return iValue;
        }
    }
    return -9999.;
}



void VPlotOptimalCut::plotHistograms( string iVariable, int i_opt_n, bool bPrint )
{
    if( !fFile )
    {
        return;
    }
    
    char hname[200];
    char htitle[200];
    sprintf( hname, "c%s", iVariable.c_str() );
    sprintf( htitle, "%s", iVariable.c_str() );
    TCanvas* c = new TCanvas( hname, htitle, 10, 10, 500, 400 );
    c->SetGridx( 0 );
    c->SetGridy( 0 );
    c->Draw();
    
    sprintf( hname, "hOn%s", iVariable.c_str() );
    TH1F* hon = ( TH1F* )fFile->Get( hname );
    if( !hon )
    {
        return;
    }
    hon->SetStats( 0 );
    hon->SetTitle( "" );
    hon->GetXaxis()->SetTitle( iVariable.c_str() );
    hon->Draw();
    sprintf( hname, "hOff%s", iVariable.c_str() );
    TH1F* hoff = ( TH1F* )fFile->Get( hname );
    if( !hoff )
    {
        return;
    }
    hoff->Draw( "same" );
    
    double i_opt_min = getCutValue( iVariable + "_min", i_opt_n );
    double i_opt_max = getCutValue( iVariable + "_max", i_opt_n );
    
    cout << "Optimal cut values (" << iVariable << "): " << i_opt_min << ", " << i_opt_max << endl;
    
    // size is on a log scale, also EChi2 and EmissionheightChi2
    if( iVariable.find( "SizeSecond" ) != string::npos || iVariable.find( "Chi2" ) != string::npos )
    {
        if( i_opt_min > 0. )
        {
            i_opt_min = log10( i_opt_min );
        }
        if( i_opt_max > 0. )
        {
            i_opt_max = log10( i_opt_max );
        }
    }
    
    if( i_opt_min > -9998. )
    {
        TLine* iL = new TLine( i_opt_min, hon->GetMinimum(), i_opt_min, hon->GetMaximum() );
        iL->SetLineStyle( 2 );
        iL->SetLineColor( 3 );
        iL->Draw();
    }
    if( i_opt_max > -9998. )
    {
        TLine* iL = new TLine( i_opt_max, hon->GetMinimum(), i_opt_max, hon->GetMaximum() );
        iL->SetLineStyle( 2 );
        iL->SetLineColor( 4 );
        iL->Draw();
    }
    
    // print all canvases
    if( bPrint )
    {
        sprintf( hname, "%s.eps", iVariable.c_str() );
        c->Print( hname );
    }
    
}

/*

    plot observing time vs value of a variable

*/
void VPlotOptimalCut::plotOptimalCuts( string iVariable, int iSourceStrength, bool iMax, double size, double iMaxObs )
{
    if( !fData )
    {
        return;
    }
    
    char hname[200];
    char htitle[200];
    sprintf( hname, "c%s_opt_%d_%d_%d", iVariable.c_str(), iMax, ( int )size, iSourceStrength );
    sprintf( htitle, "optimal cut (%s, %d, %d, %d)", iVariable.c_str(), iMax, ( int )size, iSourceStrength );
    TCanvas* c = new TCanvas( hname, htitle, 610, 10, 500, 400 );
    c->SetGridx( 0 );
    c->SetGridy( 0 );
    c->Draw();
    if( !iMax )
    {
        sprintf( hname, "obs5sigma[%d]:%s_min", iSourceStrength, iVariable.c_str() );
    }
    else
    {
        sprintf( hname, "obs5sigma[%d]:%s_max", iSourceStrength, iVariable.c_str() );
    }
    if( size > 0 )
    {
        sprintf( htitle, "TMath::Abs( SizeSecondMax_min - %f ) < 10.&&obs5sigma[%d]<%f", size, iSourceStrength, iMaxObs );
    }
    else
    {
        sprintf( htitle, "SizeSecondMax_min>0.&&obs5sigma[%d]<%f", iSourceStrength, iMaxObs );
    }
    cout << hname << endl;
    cout << htitle << endl;
    fData->Draw( hname, htitle );
}

int VPlotOptimalCut::findOptimalCut( int iSourceStrength, double size )
{
    if( !fData )
    {
        return -1;
    }
    
    // hardwired: total number of source strengths
    double iOpt5sigma[5];
    if( iSourceStrength > 4 )
    {
        return -1;
    }
    fData->SetBranchAddress( "obs5sigma", iOpt5sigma );
    double iSize = 0.;
    fData->SetBranchAddress( "SizeSecondMax_min", &iSize );
    
    double iopt = 1.e20;
    int iopt_n = -100;
    for( int i = 0; i < fData->GetEntries(); i++ )
    {
        fData->GetEntry( i );
        
        if( size > 0. )
        {
            if( TMath::Abs( size - iSize ) > 10 )
            {
                continue;
            }
        }
        
        if( iOpt5sigma[iSourceStrength] > 0. && iOpt5sigma[iSourceStrength] < iopt )
        {
            iopt = iOpt5sigma[iSourceStrength];
            iopt_n = i;
        }
    }
    if( iopt_n >= 0 )
    {
        fData->Show( iopt_n );
    }
    else
    {
        cout << "no optimum found" << endl;
    }
    
    fData->ResetBranchAddresses();
    
    return iopt_n;
    
}

void VPlotOptimalCut::plotScanParameter( string iVar1, string iVar2 )
{
    // iVar1 = "MSCW_max";
    
    string varY = "sig[4]";  //// plot significance vs variable
    
    const int nVar = 4; //number of variables
    /// fixed parameter cut values (if not scanned)
    double fix_MSCW_max = 0.35;
    double fix_MSCL_max = 0.7;
    double fix_theta2_max = 0.008;
    double fix_SizeSecondMax_min = 700;
    
    /// add fixed parameters
    char cutexp[nVar][600];
    int cVar = 0;
    if( iVar1 != "MSCW_max" && iVar2 != "MSCW_max" )
    {
        sprintf( cutexp[cVar], "MSCW_max == %f", fix_MSCW_max );
        cVar++;
    }
    if( iVar1 != "MSCL_max" && iVar2 != "MSCL_max" )
    {
        sprintf( cutexp[cVar], "MSCL_max == %f", fix_MSCL_max );
        cVar++;
    }
    if( iVar1 != "theta2_max" && iVar2 != "theta2_max" )
    {
        sprintf( cutexp[cVar], "theta2_max == %f", fix_theta2_max );
        cVar++;
    }
    if( iVar1 != "SizeSecondMax_min" && iVar2 != "SizeSecondMax_min" )
    {
        sprintf( cutexp[cVar], "SizeSecondMax_min == %f", fix_SizeSecondMax_min );
        cVar++;
    }
    
    char hname[600]; //TEST
    sprintf( hname, "SizeSecondMax_min == %f", fix_SizeSecondMax_min ); //TEST
    
    vector< string > iDraw;
    vector< string > iCut;
    vector< string > iName;
    vector< string > iTitle;
    iName.push_back( iVar1 );
    iDraw.push_back( varY + ":" + iVar1 );
    iCut.push_back( hname );
    iTitle.push_back( hname );
    
    TCanvas* c = new TCanvas( "c", "c", 100, 100, 450, 450 );
    //c->SetLeftMargin( 0.13 );
    c->Draw();
    
    int i = 0; //TEST
    string plot_option = "";
    
    fData->SetMarkerStyle( 20 );
    fData->SetMarkerSize( 1.0 );
    fData->Draw( iDraw[i].c_str(), iCut[i].c_str(), plot_option.c_str() );
    
}
