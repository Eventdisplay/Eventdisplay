/*! \file VTS.analyzeMuonRings.cpp
    \brief apply cuts for muon rings and calculate the slope of size vs radius
*/

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

#include "Ctelconfig.h"
#include "Cshowerpars.h"
#include "Ctpars.h"
#include "VGlobalRunParameter.h"
#include "VDB_Connection.h"

#include <iostream>
#include <vector>

using namespace std;

int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "VTS.analyzeMuonRings (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "--------------------------------" << endl;
    if( argc < 6 )
    {
        cout << "VTS.analyzeMuonRings <input file> <create plot (yes/no)> <write to DB (yes/no)> <DB user name> <DB password> [plot dir]" << endl;
        cout << endl;
        exit( 0 );
    }
    cout << endl;
    
    string ifile = argv[1];
    string iplot = argv[2];
    string iDB = argv[3];
    string iDB_user = argv[4];
    string iDB_passwd = argv[5];
    string dir = ".";
    if( argc >= 7 )
    {
        dir = argv[6];
    }
    
    bool bplot = false;
    bool bDB = false;
    if( iplot == "yes" || iplot == "y" )
    {
        bplot = true;
    }
    if( iDB == "yes" || iDB == "y" )
    {
        bDB = true;
    }
    /////// HARDWIRED CUTS /////////
    Int_t CUTmuonValid = 1;
    Int_t CUThoughMuonValid = 1;
    Float_t CUTLOmuonRadius = 0.4;
    Float_t CUTUPmuonRadius = 1.55;
    Float_t CUTUPmuonRSigma = 0.15;
    Float_t CUTUPmuonSize = 10000.;
    
    // arrays for muon size vs radius
    Float_t  x[10000];
    Float_t  y[10000];
    Float_t xe[10000];
    Float_t ye[10000];
    Int_t j = 0;
    
    // open input file
    TFile* f = new TFile( ifile.c_str() );
    if( f->IsZombie() )
    {
        cout << "error opening input file: " << ifile << endl;
        exit( 0 );
    }
    // get telconfig tree
    TTree* tel = ( TTree* )f->Get( "telconfig" );
    if( !tel )
    {
        cout << "error finding tree telconfig" << endl;
        exit( 0 );
    }
    Ctelconfig* ctel = new Ctelconfig( tel );
    if( !ctel )
    {
        exit( 0 );
    }
    ctel->GetEntry( 0 );
    /// number of telescopes
    const unsigned int ntel = ctel->NTel;
    TGraphErrors* gRadiusSize[ntel];
    /// slope, number of muon events passing cuts
    double slope[ntel];
    double slopeErr[ntel];
    int nMuons[ntel];
    bool bMinEvents[ntel];
    int minEvents = 3;
    // get showerpars tree
    TTree* s = ( TTree* )f->Get( "showerpars" );
    if( !s )
    {
        cout << "error finding tree showerpars" << endl;
        exit( 0 );
    }
    Cshowerpars* m = new Cshowerpars( s, false, true );
    if( !m )
    {
        exit( 0 );
    }
    if( m->fChain->GetEntries() < 3 )
    {
        cout << "error: not enough events, exiting" << endl;
        exit( 0 );
    }
    m->GetEntry( 1 );
    int RunNumber = m->runNumber;
    // get tpars trees
    char hname[500];
    for( unsigned int i = 0; i < ntel; i++ )
    {
        nMuons[i] = 0;
        /// read tpars tree
        sprintf( hname, "Tel_%d/tpars", i + 1 );
        TTree* t = ( TTree* )f->Get( hname );
        if( !t )
        {
            cout << "no tree tpars for telescope " << i + 1 << " skipping this telescope" << endl;
            continue;
        }
        Ctpars* c = new Ctpars( t, false, 1 );
        if( !c )
        {
            exit( 0 );
        }
        if( m->fChain->GetEntries() != c->fChain->GetEntries() )
        {
            cout << "error: different number of events in showerpars and tpars trees for telescope " << i + 1 << endl;
            cout << "...ignore this telescope" << endl;
            continue;
        }
        /// check if muon parameters are in the tpars tree
        if( !c->fChain->GetBranchStatus( "houghMuonValid" ) )
        {
            cout << endl << "error: input file does not contain muon parameters, exiting" << endl << endl;
            exit( 0 );
        }
        cout << "telescope " << i + 1 << " entries in tpars tree: " << c->fChain->GetEntries() << endl;
        int Nentries = c->fChain->GetEntries();
        /// initialize arrays
        for( int k = 0; k < 10000; k++ )
        {
            j = 0;
            x[k] = 0;
            y[k] = 0;
            xe[k] = 0;
            ye[k] = 0;
        }
        cout << "\t (looping over entries)" << endl;
        for( int n = 0; n < Nentries; n++ )
        {
            c->GetEntry( n );
            //// apply cuts ////
            if( c->muonValid != CUTmuonValid && c->houghMuonValid != CUThoughMuonValid )
            {
                continue;
            }
            if( ( c->muonValid == CUTmuonValid && c->muonRadius > CUTLOmuonRadius && c->muonRSigma < CUTUPmuonRSigma && ( sqrt( pow( c->muonX0 , 2 ) + pow( c->muonY0 , 2 ) ) + c->muonRadius < CUTUPmuonRadius ) && c->muonSize < CUTUPmuonSize ) || c->houghMuonValid == CUThoughMuonValid )
            {
                if( j > 10000 )
                {
                    cout << "VTS.analyzeMuonRings Warning: Found more than 10000 muon rings. Will ignore the rest of the muons for now" << endl;
                    break;
                }
                x[j] = c->muonIPCorrectedSize;
                y[j] = c->muonRadius;
                xe[j] = 0.0; //no error on the size
                ye[j] = c->muonRSigma;
                j++;
            }
        }
        gRadiusSize[i] = new TGraphErrors( j, x, y, xe, ye );
        nMuons[i] = j;
        cout << "\t events passing cuts: " << j << endl;
    }
    
    //// fit muon size vs radius
    TF1* g2 = new TF1( "g2", "[0]*x" );
    g2->SetParName( 0, "slope" );
    for( unsigned int i = 0; i < ntel; i++ )
    {
        cout << "telescope " << i + 1 << endl;
        slope[i] = 0;
        slopeErr[i] = 0;
        bMinEvents[i] = false;
        if( nMuons[i] > minEvents )
        {
            bMinEvents[i] = true;
            gRadiusSize[i]->Fit( g2 );
            //inverse slope from fit
            slope[i] = ( 1.0 / g2->GetParameter( 0 ) );
            slopeErr[i] = ( 1.0 / g2->GetParameter( 0 ) ) * ( 1.0 / g2->GetParameter( 0 ) ) * ( g2->GetParError( 0 ) );
        }
        cout << "\t size vs. radius: " << slope[i] << " +/- " << slopeErr[i] << endl;
    }
    
    //// plot muon size vs radius
    if( bplot )
    {
        TCanvas* canRadiusSize[ntel];
        char ctitle[100];
        TF1* fitRadiusSize[ntel];
        for( unsigned int i = 0; i < ntel; i++ )
        {
            if( !bMinEvents[i] )
            {
                continue;
            }
            sprintf( ctitle, "muon_tel_%i", i + 1 );
            canRadiusSize[i] = new TCanvas( ctitle, ctitle, 10, 10, 600, 600 );
            canRadiusSize[i]->cd();
            gRadiusSize[i]->GetXaxis()->SetRangeUser( 0.0, 11000.0 );
            gRadiusSize[i]->GetXaxis()->SetLimits( 0.0, 11000.0 );
            gRadiusSize[i]->GetXaxis()->SetTitleOffset( 1.1 );
            gRadiusSize[i]->GetXaxis()->SetTitle( "size (dc)" );
            gRadiusSize[i]->GetYaxis()->SetRangeUser( 0.0, 1.9 );
            gRadiusSize[i]->GetYaxis()->SetLimits( 0.0, 1.9 );
            gRadiusSize[i]->GetYaxis()->SetTitleOffset( 1.3 );
            gRadiusSize[i]->GetYaxis()->SetLabelSize( 0.035 );
            gRadiusSize[i]->GetYaxis()->SetTitle( "radius (deg)" );
            gRadiusSize[i]->SetTitle( "" );
            gRadiusSize[i]->SetLineColor( 1 );
            gRadiusSize[i]->SetMarkerColor( 1 );
            gRadiusSize[i]->SetMarkerStyle( 20 );
            gRadiusSize[i]->SetMarkerSize( 1.0 );
            fitRadiusSize[i] = gRadiusSize[i]->GetFunction( "g2" );
            if( fitRadiusSize[i] )
            {
                fitRadiusSize[i]->SetLineColor( 1 );
                fitRadiusSize[i]->SetLineStyle( 2 );
                fitRadiusSize[i]->SetRange( 0.0, 11000.0 );
            }
            gRadiusSize[i]->Draw( "APZ" );
            sprintf( ctitle, "%s/muon_%i_T%i.pdf", dir.c_str(), RunNumber, i + 1 );
            canRadiusSize[i]->Print( ctitle );
        }
    }
    
    /// write slope and error to DB
    if( bDB )
    {
        VGlobalRunParameter* globalrunpar = new VGlobalRunParameter();
        string DBserver = globalrunpar->getDBServer();
        string DBstring = DBserver + "/VOFFLINE?local";
        cout << "writing to database: " << DBstring << endl;
        
        VDB_Connection my_connection( DBstring.c_str(), iDB_user.c_str(), iDB_passwd.c_str() );
        if( !my_connection.Get_Connection_Status() )
        {
            cout << "error: no connection to database, exiting" << endl;
            exit( 0 );
        }
        /// check if a row exists for this version of evndisp in the table tblCalib_Description
        unsigned int EDversion = VGlobalRunParameter::getEVNDISP_VERSION_UI();
        char c_query[1000];
        sprintf( c_query, "SELECT * FROM tblCalib_Description WHERE package='EVNDISP' AND version=%u", EDversion );
        if( !my_connection.make_query( c_query ) )
        {
            cout << "error: DB" << endl;
            exit( 0 );
        }
        TSQLResult* dbDescription = my_connection.Get_QueryResult();
        bool bNewDescription = false;
        int calib_ID = 0;
        TSQLRow* Description_row = dbDescription->Next();
        if( !Description_row || !Description_row->GetField( 0 ) )
        {
            cout << "\t no existing row for EVNDISP version " << EDversion << " in tblCalib_Description" << endl;
            bNewDescription = true;
        }
        else
        {
            calib_ID = atoi( Description_row->GetField( 0 ) );
        }
        
        /// if missing evndisp version create a row in the table tblCalib_Description
        if( bNewDescription )
        {
            sprintf( c_query, "INSERT INTO tblCalib_Description (package,version,calib_type,val1_str,val2_str) VALUES ('EVNDISP',%u,'MUON_RING','Slope dc/deg','Number of rings')", EDversion );
            if( !my_connection.make_query( c_query ) )
            {
                cout << "error: failed to create new row in tblCalib_Description" << endl;
                exit( 0 );
            }
            cout << "\t created new row for EVNDISP version " << EDversion << " in tblCalib_Description" << endl;
            sprintf( c_query, "SELECT * FROM tblCalib_Description WHERE package='EVNDISP' AND version=%u", EDversion );
            if( !my_connection.make_query( c_query ) )
            {
                cout << "error: DB" << endl;
                exit( 0 );
            }
            TSQLResult* dbDescriptionNew = my_connection.Get_QueryResult();
            TSQLRow* DescriptionNew_row = dbDescriptionNew->Next();
            if( !DescriptionNew_row || !DescriptionNew_row->GetField( 0 ) )
            {
                cout << "\t still no row for EVNDISP version " << EDversion << " in tblCalib_Description" << endl;
                exit( 0 );
            }
            else
            {
                calib_ID = atoi( Description_row->GetField( 0 ) );
            }
        }
        cout << "\t calib_ID = " << calib_ID << endl;
        if( calib_ID == 0 )
        {
            cout << "error: invalid calib_ID, exiting" << endl;
            exit( 0 );
        }
        /// write muon slope results to table tblCalib_Run_Data
        TSQLResult* dbRun[ntel];
        TSQLRow* Run_row[ntel];
        bool bNewData[ntel];
        for( unsigned int i = 0; i < ntel; i++ )
        {
        
            if( !bMinEvents[i] )
            {
                continue;
            }
            bNewData[i] = false;
            sprintf( c_query, "SELECT calib_id,run_id,telescope_id FROM tblCalib_Run_Data WHERE calib_id=%i AND run_id=%i AND telescope_id=%i", calib_ID, RunNumber, i );
            if( !my_connection.make_query( c_query ) )
            {
                cout << "error: failed to connect to tblCalib_Run_Data" << endl;
                exit( 0 );
            }
            dbRun[i] = my_connection.Get_QueryResult();
            Run_row[i] = dbRun[i]->Next();
            if( !Run_row[i] || !Run_row[i]->GetField( 0 ) )
            {
                cout << "\t no existing row for T" << i + 1 << " in tblCalib_Run_Data" << endl;
                bNewData[i] = true;
            }
            if( bNewData[i] )
            {
                sprintf( c_query, "INSERT INTO tblCalib_Run_Data (calib_id,run_id,telescope_id,val1,err1,val2) VALUES (%i, %i, %i, %.1f, %.1f, %i )", calib_ID, RunNumber, i, slope[i], slopeErr[i], nMuons[i] );
                
                if( !my_connection.make_query( c_query ) )
                {
                    cout << "error: failed to create new row in tblCalib_Run_Data for telescope " << i + 1 << endl;
                    exit( 0 );
                }
                cout << "\t created new row in tblCalib_Run_Data for telescope " << i + 1 << endl;
            }
            else
            {
                sprintf( c_query, "UPDATE tblCalib_Run_Data SET val1 = %.1f, err1 = %.1f, val2 = %i WHERE calib_id=%i AND run_id=%i AND telescope_id=%i", slope[i], slopeErr[i], nMuons[i], calib_ID, RunNumber, i );
                if( !my_connection.make_query( c_query ) )
                {
                    cout << "error: failed to update row in tblCalib_Run_Data for telescope " << i + 1 << endl;
                    exit( 0 );
                }
                cout << "\t updated row in tblCalib_Run_Data for telescope " << i + 1 << endl;
            }
        }
        
    }
}
