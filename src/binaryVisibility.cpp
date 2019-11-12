/*! \file visibility.cc
    \brief plot object and moon elevation vs MJD

    libnova package required

*/


#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"
#include "TROOT.h"


#include "VLibNovaStar.h"
#include "VLibNovaSunAndMoon.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

bool fDebug = false;

struct sRunParameter
{
    string fObservatory_name;
    double fObservatory_latitude;
    double fObservatory_longitude;
    
    double fSunMaxElevation;
    double fMoonMaxElevation;
    double fObservingMinElevation;
    double fObservingMaxElevation;
    double fObservingMinElevation_plotting;
    
    string fObject_name;
    string fObject_printname;
    double fObject_dec;
    double fObject_ra;
    double fObject_orbit;
    double fObject_t0;
    unsigned int fObject_nphases;
};

sRunParameter getRunParameter( string ifile )
{
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error opening runparameter file: " << ifile << endl;
        exit( -1 );
    }
    string is_line;
    string temp;
    
    sRunParameter a;
    a.fObject_orbit = -1;
    a.fObject_t0 = -1;
    a.fObject_nphases = 0;
    a.fObject_printname = "";
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() ==  0 )
        {
            continue;
        }
        istringstream is_stream( is_line );
        is_stream >> temp;
        if( temp != "*" )
        {
            continue;
        }
        is_stream >> temp;
        if( temp == "OBSERVATORY" )
        {
            is_stream >> a.fObservatory_longitude;
            is_stream >> a.fObservatory_latitude;
            is_stream >> a.fObservatory_name;
        }
        else if( temp == "MAX_SUN_ELEVATION" )
        {
            is_stream >> a.fSunMaxElevation;
        }
        else if( temp == "MAX_MOON_ELEVATION" )
        {
            is_stream >> a.fMoonMaxElevation;
        }
        else if( temp == "MINMAX_OBSERVING_ELEVATION" )
        {
            is_stream >> a.fObservingMinElevation;
            is_stream >> a.fObservingMaxElevation;
            is_stream >> a.fObservingMinElevation_plotting;
        }
        else if( temp == "OBJECT" )
        {
            is_stream >> a.fObject_ra;
            is_stream >> a.fObject_dec;
            a.fObject_name = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
        }
        else if( temp == "ORBITALPARAMETER" )
        {
            is_stream >> a.fObject_orbit;
            is_stream >> a.fObject_t0;
            is_stream >> a.fObject_nphases;
        }
        else if( temp == "PRINTNAME" )
        {
            is_stream >> a.fObject_printname;
        }
        else
        {
            cout << "unknown parameter: " << endl;
            cout << is_line << endl;
        }
    }
    is.close();
    
    return a;
}


int main( int argc, char* argv[] )
{

    //////////////////////////////////////////////
    // read input parameters
    if( argc != 5 )
    {
        cout << endl;
        cout << "./binaryVisibility <MJDstart> <MJDstop> <binning [min]> <parameter file>" << endl;
        cout << endl;
        cout << "   <binning [min]> : typically 5 min (shorter=more accurate, but slower)" << endl;
        cout << endl;
        exit( 0 );
    }
    cout << endl;
    cout << "binaryVisibility" << endl;
    cout << "----------------" << endl;
    cout << endl;
    // start date in MJD
    double fMJDStart = atof( argv[1] );
    // stopp date in MJD
    double fMJDStopp = atof( argv[2] );
    if( fMJDStart >= fMJDStopp )
    {
        cout << "invalid MJD range " << fMJDStart << " " << fMJDStopp << endl;
        exit( 0 );
    }
    // binning in time
    double fMJDInt = atof( argv[3] );
    if( fMJDInt <= 0. )
    {
        cout << "binning must be > 0." << endl;
        exit( 0 );
    }
    string fParaFile = "runparameter.dat";
    if( argc == 5 )
    {
        fParaFile = argv[4];
    }
    
    sRunParameter fRunParameter = getRunParameter( fParaFile );
    
    // binning and number of MJD days
    double fMJDBin = 1. / ( 24.*60. / fMJDInt );
    unsigned int fMJDDays = ( unsigned int )( ( fMJDStopp - fMJDStart ) + 0.5 );
    
    cout << "calculating observability for " << fRunParameter.fObject_name;
    cout << " (ra=" << fRunParameter.fObject_ra << ", dec=" << fRunParameter.fObject_dec << ")" << endl;
    if( fRunParameter.fObject_orbit > 0. )
    {
        cout << "\t orbital parameters: " << fRunParameter.fObject_orbit << " " << fixed << fRunParameter.fObject_t0;
        cout << " " << fRunParameter.fObject_nphases << endl;
    }
    cout << "\t require minimum elevation of " << fRunParameter.fObservingMinElevation << " deg" << endl;
    cout << "\t number of days: " << fMJDDays;
    cout << ", binning in MJD: " << fMJDBin << endl;
    
    ///////////////////////////////////////
    // astronomical objects
    
    // sun and moon rise time
    VLibNovaSunAndMoon* fSunAndMoon = new VLibNovaSunAndMoon( fRunParameter.fObservatory_longitude, fRunParameter.fObservatory_latitude );
    
    // object
    VLibNovaStar* fObject = new VLibNovaStar( fRunParameter.fObject_ra, fRunParameter.fObject_dec, fRunParameter.fObservatory_longitude, fRunParameter.fObservatory_latitude );
    
    // time and elevation
    int fNMJD = ( int )( ( fMJDStopp - fMJDStart ) / fMJDBin );
    // vector of length 3: no moon, illumination < 30%, illumination > 30%
    vector< double > fMoonFrac_min;
    vector< double > fMoonFrac_max;
    fMoonFrac_min.push_back( -1. );
    fMoonFrac_max.push_back( 0. );
    fMoonFrac_min.push_back( 0.0 );
    fMoonFrac_max.push_back( 0.3 );
    fMoonFrac_min.push_back( 0.3 );
    fMoonFrac_max.push_back( 1.0 );
    vector< int > fNObj( fMoonFrac_min.size(), 0 );
    vector< double > i_temp;
    vector< vector< double > > fObjMJD( 3, i_temp );
    vector< vector< double > > fObjElevation( 3, i_temp );
    
    vector< double > fMoonMJD;
    vector< double > fMoonElevation;
    vector< double > fMoonIllumination;
    
    /////////////////////////////////////////////
    // loop over time interval
    cout << "\t number of intervals in time: " << fNMJD << endl;
    
    double fMJD;
    double iMoonE = 0.;
    double iSunE = 0.;
    double iObjE = 0.;
    double iMoonF = 0.;
    
    // total time in different phase bins
    vector< vector< double > > fObsTimePhase;
    for( unsigned int i = 0; i < fMoonFrac_min.size(); i++ )
    {
        i_temp.assign( fRunParameter.fObject_nphases, 0. );
        fObsTimePhase.push_back( i_temp );
    }
    double fPhase = 0.;
    
    ///////////////////////////////////////////////////
    // loop over all time bins
    for( int i = 0; i < fNMJD; i++ )
    {
        fMJD = fMJDStart + fMJDBin * i;
        
        fSunAndMoon->setMJD( fMJD );
        iSunE = fSunAndMoon->getSunElevation();
        
        iMoonE = fSunAndMoon->getMoonElevation();
        iMoonF = fSunAndMoon->getMoonDisk();
        
        iObjE = fObject->getElevation( fMJD );
        
        ///////////////////////////
        // orbital phase
        if( fRunParameter.fObject_t0 > 0. && fRunParameter.fObject_orbit > 0. )
        {
            fPhase = ( fMJD + 2400000.5 - fRunParameter.fObject_t0 ) / fRunParameter.fObject_orbit
                     - int( ( fMJD + 2400000.5 - fRunParameter.fObject_t0 ) / fRunParameter.fObject_orbit - 0.5 );
            if( fPhase > 1. )
            {
                fPhase -= 1.;
            }
        }
        else
        {
            fPhase = 0;
        }
        
        ///////////////////////////
        // fill object elevation
        
        // object above minimum elevation and sun below horizon
        if( iObjE > fRunParameter.fObservingMinElevation && iSunE < fRunParameter.fSunMaxElevation && iMoonE < fRunParameter.fMoonMaxElevation )
        {
            if( fDebug )
            {
                cout << "MDJ " << fMJD << "\t Phase " << fPhase << "\t Phase * 10 " << ( int )( fPhase * 10. ) << " Moon " << iMoonF << " " << iMoonE << endl;
            }
            
            // loop over all moon illuminations
            for( unsigned int j = 0; j < fMoonFrac_min.size(); j++ )
            {
                if( iMoonF > fMoonFrac_min[j] && iMoonF < fMoonFrac_max[j] )
                {
                    // moon below horizon: this is always filled into the first bin
                    if( iMoonE < 0. )
                    {
                        fObjMJD[0].push_back( fMJD );
                        fObjElevation[0].push_back( iObjE );
                        fObsTimePhase[0][( int )( fPhase * 10. )] += fMJDInt / 60.;
                    }
                    else
                    {
                        fObjMJD[j].push_back( fMJD );
                        fObjElevation[j].push_back( iObjE );
                        fObsTimePhase[j][( int )( fPhase * 10. )] += fMJDInt / 60.;
                    }
                }
            }
        }
        
        // fill moon elevation
        fMoonMJD.push_back( fMJD );
        fMoonElevation.push_back( iMoonE );
        fMoonIllumination.push_back( iMoonF );
    }
    
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // draw everything
    TApplication app( "app", &argc, argv );
    gROOT->SetStyle( "Plain" );
    
    TCanvas* cCan = new TCanvas( "cCan", "observability", 10, 10, 1000, 500 );
    cCan->Draw();
    
    TH1D* hnull = new TH1D( "hnull", "", 100, fMJDStart, fMJDStopp );
    hnull->SetStats( 0 );
    hnull->SetMaximum( fRunParameter.fObservingMaxElevation );
    hnull->SetMinimum( fRunParameter.fObservingMinElevation_plotting );
    hnull->SetXTitle( "MJD" );
    hnull->SetYTitle( "elevation [deg]" );
    
    hnull->Draw();
    
    // one graph for every moon illumation and every MJD
    vector< TGraph* > i_graph;
    vector< vector< TGraph* > > gObject;
    int z = 0;
    int i_oldMJD = 0;
    int i_newMJD = 0;
    //   for( unsigned int i = 0; i < fMoonFrac.size(); i++ )
    for( unsigned int f = fMoonFrac_min.size(); f > 0; f-- )
    {
        unsigned int i = f - 1;
        cout << "Moon illumination " << fMoonFrac_min[i] << ", color " << i + 2 << endl;
        i_graph.clear();
        z = 0;
        i_oldMJD = 0;
        i_newMJD = 0;
        for( unsigned int j = 0; j < fObjMJD[i].size(); j++ )
        {
            i_newMJD = ( int )( fObjMJD[i][j] );
            if( i_newMJD > i_oldMJD )
            {
                if( z > 0 )
                {
                    i_graph.back()->Draw( "p" );
                }
                
                i_graph.push_back( new TGraph( 1 ) );
                i_graph.back()->SetLineWidth( 3 );
                i_graph.back()->SetLineColor( i + 2 );
                i_graph.back()->SetMarkerColor( i + 2 );
                i_graph.back()->SetMarkerStyle( 7 );
                i_oldMJD = i_newMJD;
                z = 0;
            }
            if( i_graph.back() )
            {
                i_graph.back()->SetPoint( z, fObjMJD[i][j], fObjElevation[i][j] );
                z++;
            }
        }
        gObject.push_back( i_graph );
    }
    
    // moon
    i_oldMJD = 0;
    i_newMJD = 0;
    vector< TGraph* > gMoon;
    z = 0;
    for( unsigned int j = 0; j < fMoonMJD.size(); j++ )
    {
        i_newMJD = ( int )( fMoonMJD[j] );
        
        if( fMoonElevation[j] < 0. || gMoon.size() == 0 )
        {
            if( z > 0 )
            {
                gMoon.back()->Draw( "C" );
            }
            
            gMoon.push_back( new TGraph( 1 ) );
            gMoon.back()->SetLineWidth( 1 );
            gMoon.back()->SetLineStyle( 3 );
            for( unsigned int l = 0; l < fMoonFrac_min.size(); l++ )
            {
                if( fMoonIllumination[j] > fMoonFrac_min[l] && fMoonIllumination[j] < fMoonFrac_max[l] )
                {
                    gMoon.back()->SetLineColor( l + 2 );
                    gMoon.back()->SetMarkerColor( l + 2 );
                    break;
                }
            }
            
            i_oldMJD = i_newMJD;
            z = 0;
        }
        if( gMoon.back() )
        {
            gMoon.back()->SetPoint( z, fMoonMJD[j], fMoonElevation[j] );
            z++;
        }
    }
    
    // draw lines at different phases
    cout << endl;
    char hname[400];
    if( fRunParameter.fObject_orbit > 0. && fRunParameter.fObject_t0 > 0 && fRunParameter.fObject_nphases > 0 )
    {
        double fPhaseBin = 1. / ( double )fRunParameter.fObject_nphases;
        int lastBin = -1;
        unsigned int nloop = fRunParameter.fObject_nphases * ( int )fRunParameter.fObject_orbit;
        for( unsigned int i = 0; i < nloop; i++ )
        {
            fMJD = fMJDStart + i * ( fMJDStopp - fMJDStart ) / ( double )nloop;
            
            fPhase = ( fMJD + 2400000.5 - fRunParameter.fObject_t0 ) / fRunParameter.fObject_orbit - int( ( fMJD + 2400000.5 - fRunParameter.fObject_t0 ) / fRunParameter.fObject_orbit - 0.5 );
            if( fPhase > 1. )
            {
                fPhase -= 1.;
            }
            
            if( ( int )( fPhase / fPhaseBin ) != lastBin && i > 0 && fPhase / fPhaseBin - ( int )( fPhase / fPhaseBin ) < 0.3 )
            {
            
                TLine* iL = new TLine( fMJD, fRunParameter.fObservingMinElevation_plotting, fMJD, fRunParameter.fObservingMaxElevation );
                iL->SetLineStyle( 2 );
                if( ( int )( fPhase / fPhaseBin ) == 0 )
                {
                    iL->SetLineWidth( 3 );
                }
                iL->Draw();
                lastBin = ( int )( fPhase / fPhaseBin );
                
                if( ( ( fMJDStopp - fMJDStart ) / fRunParameter.fObject_orbit * fRunParameter.fObject_nphases ) < 20. )
                {
                    sprintf( hname, "%.1f", fPhase );
                    TText* iT = new TText( fMJD + 0.05, 0.9 * ( fRunParameter.fObservingMaxElevation - fRunParameter.fObservingMinElevation_plotting ), hname );
                    iT->SetTextSize( 0.6 * iT->GetTextSize() );
                    iT->Draw();
                }
            }
        }
    }
    
    // print to file
    if( fRunParameter.fObject_printname.size() > 0 )
    {
        sprintf( hname, "Visibility%s_MJD%d.eps", fRunParameter.fObject_printname.c_str(), ( int )fMJDStart );
        cCan->Print( hname );
        sprintf( hname, "Visibility%s_MJD%d.jpg", fRunParameter.fObject_printname.c_str(), ( int )fMJDStart );
        cCan->Print( hname );
    }
    // available time for observations
    cout << endl;
    cout << "observability: " << endl;
    double fMaxObs = 0.;
    vector< double > fObsTot( fMoonFrac_min.size(), 0. );
    for( unsigned int i = 0; i < fMoonFrac_min.size(); i++ )
    {
        cout << "\t Moon illumination [" << fMoonFrac_min[i] << "," << fMoonFrac_max[i] << "]" << endl;
        for( unsigned int j = 0; j < fObsTimePhase[i].size(); j++ )
        {
            cout << "\t Phase " << setprecision( 2 ) << 1. / ( double )fObsTimePhase[i].size()  * j << " - " << 1. / ( double )fObsTimePhase[i].size()  * ( j + 1 ) << ": ";
            cout << setprecision( 2 ) << fObsTimePhase[i][j] << " h" << endl;
            fObsTot[i] += fObsTimePhase[i][j];
            if( fObsTimePhase[i][j] > fMaxObs )
            {
                fMaxObs = fObsTimePhase[i][j];
            }
        }
    }
    cout << "Total observation time: " << endl;
    for( unsigned int i = 0; i < fMoonFrac_min.size(); i++ )
    {
        cout << "\t Moon illumination [" << fMoonFrac_min[i] << "," << fMoonFrac_max[i] << "]: " << fObsTot[i] << " h" << endl;
    }
    // plot observation time vs phase
    if( fRunParameter.fObject_orbit > 0. && fRunParameter.fObject_t0 > 0. )
    {
        TCanvas* cCObs = new TCanvas( "cObs", "observation time", 850, 10, 400, 400 );
        cCObs->Draw();
        
        vector< TH1D* > hObsTime;
        
        for( unsigned int f = fMoonFrac_min.size(); f > 0; f-- )
        {
            unsigned int i = f - 1;
            sprintf( hname, "hObs_%d", i );
            hObsTime.push_back( new TH1D( hname, "", fObsTimePhase[i].size(), 0., 1. ) );
            hObsTime.back()->SetStats( 0 );
            hObsTime.back()->SetLineColor( i + 2 );
            hObsTime.back()->SetLineWidth( 2 );
            hObsTime.back()->SetXTitle( "orbital phase" );
            hObsTime.back()->SetYTitle( "available observation time/phase interval [h]" );
            hObsTime.back()->GetYaxis()->SetTitleOffset( 1.3 );
            hObsTime.back()->SetMaximum( fMaxObs * 1.2 );
            for( unsigned int j = 0; j < fObsTimePhase[i].size(); j++ )
            {
                hObsTime.back()->SetBinContent( j + 1, fObsTimePhase[i][j] );
            }
            if( f == fMoonFrac_min.size() )
            {
                hObsTime.back()->Draw();
            }
            else
            {
                hObsTime.back()->Draw( "same" );
            }
        }
        // print results to file
        if( fRunParameter.fObject_printname.size() > 0 )
        {
            sprintf( hname, "VisibilityTot%s_MJD%d.eps", fRunParameter.fObject_printname.c_str(), ( int )fMJDStart );
            cCObs->Print( hname );
            sprintf( hname, "VisibilityTot%s_MJD%d.jpg", fRunParameter.fObject_printname.c_str(), ( int )fMJDStart );
            cCObs->Print( hname );
        }
    }
    
    app.Run();
}
