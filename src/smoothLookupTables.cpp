/*! \file smoothLookupTables
    \brief smooth lookup tables

*/

#include "TClass.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TKey.h"
#include "TSystem.h"
#include "TROOT.h"

#include "VGlobalRunParameter.h"
#include "VHistogramUtilities.h"
#include "VInterpolate2DHistos.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// list with noise levels
map< string, vector< int > > fNoiseLevel;
string fCurrentTelescopeType;

// directory level
unsigned int fDirectoryLevel = 0;
void copyDirectory( TDirectory* source, const char* hx = 0 );

// flag if woff_0500 should be copied
bool fCopy_woff_0500 = false;

TH2F* smooth2DHistogram( TH2F* h, TH2F* hNevents )
{
    if( !h )
    {
        return 0;
    }
    cout << "\t\t smooth histogram using _fit_ algorithm: (min events per bin = 20): ";
    cout << h->GetName() << endl;
    
    VInterpolate2DHistos i_inter;
    
    string iHisName = h->GetName();
    if( iHisName.find( "energy" ) != string::npos )
    {
        return  i_inter.doLogLinearExtrapolation( h,
                "fitexpo",
                hNevents,
                20 );
    }
    else
    {
        return  i_inter.doLogLinearExtrapolation( h,
                "fitpol2",
                hNevents,
                20 );
    }
    return 0;
}

/*
 * get histogram name for event histogram
 *
 */
string getNeventsHistoName( string iName )
{
    vector< string > iHNmae;
    iHNmae.push_back( "mpv" );
    iHNmae.push_back( "median" );
    iHNmae.push_back( "mean" );
    for( unsigned int i = 0; i < iHNmae.size(); i++ )
    {
        if( iName.find( iHNmae[i] ) != string::npos )
        {
            size_t start_pos = iName.find( iHNmae[i] );
            return iName.replace( start_pos, iHNmae[i].length(), "nevents" );
        }
    }
    return "";
}

/*
 * noise directory names are determined in the lookup table code using the
 * mean pedvar level. These can vary by a small avound from simulation to
 * simulation file. We search here therefore for very similar noise levels,
 * and return those directory names if available
 */
string check_for_similar_noise_values( string iTelescopeType, const char* hx )
{
    string iTemp;
    if( !hx )
    {
        cout << "check_for_similar_noise_values warning: empty string " << endl;
        return iTemp;
    }
    iTemp = hx;
    
    if( iTemp.find( "_" ) != string::npos )
    {
        // get noise level from directory name
        int i_noise = atoi( iTemp.substr( iTemp.find( "_" ) + 1, iTemp.size() ).c_str() );
        
        // check if a similar noise level exists for this telescope type
        if( fNoiseLevel.find( iTelescopeType ) != fNoiseLevel.end() )
        {
            for( unsigned int i = 0; i < fNoiseLevel[iTelescopeType].size(); i++ )
            {
                if( TMath::Abs( fNoiseLevel[iTelescopeType][i] - i_noise ) < 10 )
                {
                    char hname[200];
                    sprintf( hname, "NOISE_%05d", fNoiseLevel[iTelescopeType][i] );
                    iTemp = hname;
                    cout << "\t found similar noise level for telescope type " << iTelescopeType;
                    cout << ", save into directory: " << iTemp << "\t" << fNoiseLevel[iTelescopeType][i];
                    cout << " (" << fNoiseLevel[iTelescopeType].size() << ")" << endl;
                    return iTemp;
                }
            }
            // new noise level for existing telescope type
            fNoiseLevel[iTelescopeType].push_back( i_noise );
            cout << "\t new noise level directory for existing telescope type " << iTelescopeType << ": ";
            cout << iTemp << "(" << fNoiseLevel[iTelescopeType].size() << ")" << endl;
        }
        else
        {
            vector< int > i_tempV;
            i_tempV.push_back( i_noise );
            // new noise level for new telescope type
            fNoiseLevel[iTelescopeType] = i_tempV;
            cout << "\t new noise level directory for new telescope type " << iTelescopeType << ": ";
            cout << iTemp << "(" << fNoiseLevel[iTelescopeType].size() << ")" << endl;
        }
    }
    
    return iTemp;
    
    
}


/*
 *
 */
int main( int argc, char* argv[] )
{
    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter fRunPara;
            cout << fRunPara.getEVNDISP_VERSION() << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    
    VGlobalRunParameter* iT = new VGlobalRunParameter();
    cout << endl;
    cout << "smoothLookupTables (" << iT->getEVNDISP_VERSION() << ")" << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    
    // print help
    if( argc < 2 )
    {
        cout << "combine several tables from different files into one single table file" << endl << endl;
        cout << "combineLookupTables <input table file name> <output file name>" << endl;
        cout << endl;
        cout << "(should be used with care; check effect of smoothing first with VPlotLookupTables)" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string fIFile = argv[1];
    string fOFile = argv[2];
    
    //////////////////////////////////////
    // open output lookup table file
    TFile* fROFile = new TFile( fOFile.c_str(), "RECREATE" );
    if( fROFile->IsZombie() )
    {
        cout << "error while opening output lookup file: " << fOFile << endl;
        exit( EXIT_FAILURE );
    }
    
    //////////////////////////////////////////////////////////////////
    // loop over all lookup tables and copy them into the main file
    TFile* fIn = new TFile( fIFile.c_str() );
    if( fIn->IsZombie() )
    {
        cout << "error while opening file: " << fIFile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "now reading file: " << fIFile << endl;
    
    //loop on all entries of this directory
    TKey* key;
    TIter nextkey( fIn->GetListOfKeys() );
    while( ( key = ( TKey* )nextkey() ) )
    {
        const char* classname = key->GetClassName();
        TClass* cl = gROOT->GetClass( classname );
        if( !cl )
        {
            continue;
        }
        if( cl->InheritsFrom( "TDirectory" ) )
        {
            TDirectory* iSource = ( TDirectory* )key->ReadObj();
            if( !iSource )
            {
                continue;
            }
            fROFile->cd();
            const char* hname = iSource->GetName();
            cout << "\t copying directory for telescope type " << hname << endl;
            copyDirectory( iSource, hname );
        }
    }
    fIn->Close();
    
    fROFile->Close();
    cout << endl;
    cout << "finished..." << endl;
}


/*!
 *   from http://root.cern.ch/phpBB2/viewtopic.php?t=2789
 *
 */
void copyDirectory( TDirectory* source, const char* hx )
{
    if( !source )
    {
        cout << "copyDirectory error: zero pointer do directory " << endl;
        return;
    }
    //copy all objects and subdirs of directory source as a subdir of the current directory
    TDirectory* savdir = gDirectory;
    TDirectory* adir = 0;
    
    // 1. case: top directory exists (tel_...)
    if( hx )
    {
        fDirectoryLevel = 0;
        fCurrentTelescopeType = hx;
    }
    
    if( fDirectoryLevel == 1 )
    {
        string noise_dir = check_for_similar_noise_values( fCurrentTelescopeType, source->GetName() );
        adir = ( TDirectory* )savdir->Get( noise_dir.c_str() );
    }
    else
    {
        adir = ( TDirectory* )savdir->Get( source->GetName() );
    }
    fDirectoryLevel++;
    
    if( !adir )
    {
        // 2. case: make top directory
        if( hx )
        {
            adir = savdir->mkdir( hx );
        }
        else
        {
            adir = savdir->mkdir( source->GetName() );
        }
        if( !adir )
        {
            cout << "error while creating directory " << source->GetName() << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    adir->cd();
    //loop on all entries of this directory
    TKey* key;
    TIter nextkey( source->GetListOfKeys() );
    while( ( key = ( TKey* )nextkey() ) )
    {
        const char* classname = key->GetClassName();
        TClass* cl = gROOT->GetClass( classname );
        if( !cl )
        {
            continue;
        }
        if( cl->InheritsFrom( "TDirectory" ) )
        {
            string iName = key->GetName();
            if( iName == "histos1D" )
            {
                continue;
            }
            // do not copy first off-axis directory
            // (for CTA sims: this camera distance is
            // covered by on-axis sims)
            if( fCopy_woff_0500 && iName == "woff_0500" )
            {
                continue;
            }
            source->cd( key->GetName() );
            TDirectory* subdir = gDirectory;
            adir->cd();
            copyDirectory( subdir );
            adir->cd();
        }
        else
        {
            source->cd();
            TObject* obj = key->ReadObj();
            string iName = obj->GetName();
            if( iName.find( "new" ) < string::npos )
            {
                cout << gDirectory->GetPath() << endl;
            }
            adir->cd();
            string iDirName = gDirectory->GetName();
            ///////////////////////////////////
            // smooth histograms
            // get histogram for event counting histogram
            // smooth median, mean, and mpv histograms
            // copy only median and mpv histogram
            if( iName.find( "median" ) != string::npos
                    || iName.find( "Median" ) != string::npos
                    || iName.find( "mpv" ) != string::npos
                    || iName.find( "mean" ) != string::npos
              )
            {
                string iNeventsHistoName = getNeventsHistoName( iName );
                TH2F* iHEvents = ( TH2F* )source->Get( iNeventsHistoName.c_str() );
                if( iHEvents )
                {
                    obj = ( TObject* )smooth2DHistogram( ( TH2F* )obj, iHEvents );
                }
            }
            ///////////////////////////////////
            cout << "\t writing " << iName << " to ";
            cout << adir->GetPath() << endl;
            obj->Write( iName.c_str() );
            delete obj;
        }
    }
    adir->SaveSelf( kTRUE );
    savdir->cd();
}
