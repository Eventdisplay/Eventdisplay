/*! \file combineLookupTables
    \brief combine different lookup tablefiles into a single tablefile

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

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// list with noise levels
struct noiselevel
{
    vector< int > noise;
    vector< int > noiseFileName;
};
map< string, noiselevel > fNoiseLevel;
string fCurrentTelescopeType;

// directory level
unsigned int fDirectoryLevel = 0;
void copyDirectory( TDirectory* source, const char* hx = 0, unsigned int i_noiselevel_filename = 0 );

// flag if woff_0500 should be copied
bool fCopy_woff_0500 = false;

/*
 * extract noise from file name 
 * (is only a second order correction and depends of course on the
 * particular file nameing)
*/
unsigned int extract_noisefrom_filename( string iFN )
{
      unsigned int i_noise = 0;

      if( iFN.find( "noise" ) != string::npos )
      {
          size_t i_beg = iFN.find("noise")+5;
          size_t i_fin = iFN.find("_", i_beg) - i_beg;
          i_noise = atoi( iFN.substr(i_beg, i_fin).c_str() );
      }
      return i_noise;
}

/*
 * noise directory names are determined in the lookup table code using the
 * mean pedvar level. These can vary by a small amount from simulation to
 * simulation file. We search here therefore for very similar noise levels,
 * and return those directory names if available
 */
string check_for_similar_noise_values( string iTelescopeType, const char* hx, unsigned int noise_from_filename )
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
            for( unsigned int i = 0; i < fNoiseLevel[iTelescopeType].noise.size(); i++ )
            {
                if( TMath::Abs( fNoiseLevel[iTelescopeType].noise[i] - i_noise ) < 10 )
                {
                    char hname[200];
                    sprintf( hname, "NOISE_%05d", fNoiseLevel[iTelescopeType].noise[i] );
                    iTemp = hname;
                    cout << "\t found similar noise level for telescope type " << iTelescopeType;
                    cout << ", save into directory: " << iTemp << "\t" << fNoiseLevel[iTelescopeType].noise[i];
                    cout << " (" << fNoiseLevel[iTelescopeType].noise.size() << ")" << endl;
                    return iTemp;
                }
            }
            // use file name (should be second choisce only)
            if( noise_from_filename > 0 )
            {
                  for( unsigned int i = 0; i < fNoiseLevel[iTelescopeType].noiseFileName.size(); i++ )
                  {
                      if( (int)noise_from_filename == fNoiseLevel[iTelescopeType].noiseFileName[i] )
                      {
                          if( i < fNoiseLevel[iTelescopeType].noise.size() )
                          {
                              char hname[200];
                              sprintf( hname, "NOISE_%05d", fNoiseLevel[iTelescopeType].noise[i] );
                              iTemp = hname;
                              cout << "\t found similar noise level for telescope type " << iTelescopeType;
                              cout << ", save into directory: " << iTemp << "\t" << fNoiseLevel[iTelescopeType].noise[i];
                              cout << " (" << fNoiseLevel[iTelescopeType].noise.size() << ")" << endl;
                              cout << "\t using file name for noise level identification: " << iTelescopeType;
                              cout << ", noise levels: " << i_noise << ", " << fNoiseLevel[iTelescopeType].noise[i];
                              cout << " (file noise level " << noise_from_filename << ")" << endl;
                              return iTemp;
                           }
                       }
                   }
            }
            // new noise level for existing telescope type
            fNoiseLevel[iTelescopeType].noise.push_back( i_noise );
            fNoiseLevel[iTelescopeType].noiseFileName.push_back( noise_from_filename );
            cout << "\t new noise level directory for existing telescope type " << iTelescopeType << ": ";
            cout << iTemp << " (" << fNoiseLevel[iTelescopeType].noise.size();
            cout << ", file name " << noise_from_filename << ")";
            cout << endl;
        }
        else
        {
            vector< int > i_tempV;
            i_tempV.push_back( i_noise );
            vector< int > i_tempN;
            i_tempN.push_back( noise_from_filename );
            // new noise level for new telescope type
            fNoiseLevel[iTelescopeType].noise = i_tempV;
            fNoiseLevel[iTelescopeType].noiseFileName = i_tempN;
            cout << "\t new noise level directory for new telescope type " << iTelescopeType << ": ";
            cout << iTemp << " (" << fNoiseLevel[iTelescopeType].noise.size();
            cout << ", file name " << noise_from_filename << ")";
            cout << endl;
        }
    }
    
    return iTemp;
    
    
}

/*
 * read a list of lookup tables files
 *
 */
vector< string > readListOfFiles( string iFile )
{
    vector< string > iList;
    
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "error while reading file list " << iFile << endl;
        cout << "exiting...." << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    
    while( getline( is, is_line ) )
    {
        iList.push_back( is_line );
    }
    
    is.close();
    
    return iList;
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
    cout << "combineLookupTables (" << iT->getEVNDISP_VERSION() << ")" << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    
    // print help
    if( argc < 2 )
    {
        cout << "combine several tables from different files into one single table file" << endl << endl;
        cout << "combineLookupTables <file with list of tables> <output file name> [do not copy woff_0500 directory (default = 0 = false)] " << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string fListOfFiles = argv[1];
    string fOFile       = argv[2];
    if( argc > 3 )
    {
        fCopy_woff_0500 = ( bool )( atoi( argv[3] ) );
    }
    
    /////////////////////////////////////
    // get list of lookup table files
    vector< string > fInFiles = readListOfFiles( fListOfFiles );
    unsigned int nFiles = fInFiles.size();
    if( nFiles == 0 )
    {
        cout << "error: no files in file list" << endl;
        cout << "exiting...." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "combining " << nFiles << " table files into " << fOFile << endl;
    
    //////////////////////////////////////
    // open combined lookup table file
    TFile* fROFile = new TFile( fOFile.c_str(), "RECREATE" );
    if( fROFile->IsZombie() )
    {
        cout << "error while opening combined file: " << fOFile << endl;
        exit( EXIT_FAILURE );
    }
    
    //////////////////////////////////////////////////////////////////
    // loop over all lookup tables and copy them into the main file
    TFile* fIn = 0;
    for( unsigned int f = 0; f < fInFiles.size(); f++ )
    {
        fIn = new TFile( fInFiles[f].c_str() );
        if( fIn->IsZombie() )
        {
            cout << "error while opening file: " << fInFiles[f] << endl;
            continue;
        }
        cout << "now reading file " << f << ": " << fInFiles[f] << endl;

        unsigned int noise_from_filename = extract_noisefrom_filename( fInFiles[f] );
        if( noise_from_filename > 0 )
        {
                cout << "\t noise ID from file name: " << noise_from_filename << endl;
        }
        
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
                copyDirectory( iSource, hname, noise_from_filename );
            }
        }
        fIn->Close();
    }
    
    fROFile->Close();
    cout << endl;
    cout << "finished..." << endl;
}


/*!
 *   from http://root.cern.ch/phpBB2/viewtopic.php?t=2789
 *
 */
void copyDirectory( TDirectory* source, const char* hx, unsigned int noise_from_filename )
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
        string noise_dir = check_for_similar_noise_values( fCurrentTelescopeType, source->GetName(), noise_from_filename );
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
            copyDirectory( subdir, 0, noise_from_filename );
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
            // copy only median and mpv histogram
            if( iName.find( "median" ) != string::npos
                    || iName.find( "Median" ) != string::npos
                    || iName.find( "mpv" ) != string::npos
                    || iName.find( "mean" ) != string::npos
                    || iName.find( "nevents" ) != string::npos
              )
            {
                adir->cd();
                cout << "\t writing " << iName << " to ";
                cout << adir->GetPath() << endl;
                obj->Write( iName.c_str() );
            }
            delete obj;
        }
    }
    adir->SaveSelf( kTRUE );
    savdir->cd();
}
