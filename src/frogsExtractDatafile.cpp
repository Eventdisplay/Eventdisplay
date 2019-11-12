#include <iostream>
#include <stdlib.h>

#include "VEvndispRunParameter.h"
#include "VTableLookupRunParameter.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TKey.h"
#include "TDirectoryFile.h"

using namespace std;

void CopyDir( TDirectory* source ) ;

int main( int argc, char* argv[] )
{
    if( argc != 6 )
    {
        printf( "Error, must give right arguments 'inputfile' 'treename' 'startevent' 'nevents' 'outputfilename', exiting...\n" );
        exit( 1 );
    }
    
    // input file name
    char inputfile[200] = "" ;
    sprintf( inputfile, "%s", argv[1] ) ;
    
    // TTree name to extract
    char treename[200] = "" ;
    sprintf( treename, "%s", argv[2] ) ;
    
    // starting and N events
    int startevent = atoi( argv[3] ) ;
    int nevents    = atoi( argv[4] ) ;
    
    // output file name
    char outputfile[1000] = "" ;
    sprintf( outputfile, "%s", argv[5] ) ;
    
    printf( " inputfile: '%s'\n", inputfile ) ;
    printf( " treename:  '%s'\n", treename ) ;
    printf( " startevent:'%d'\n", startevent ) ;
    printf( " nevents:   '%d'\n", nevents ) ;
    printf( " outputfile:'%s'\n", outputfile ) ;
    
    // open input root file
    TFile* mscwfile = new TFile( inputfile, "READ" ) ;
    
    // open output extracted file
    TFile* extractfile = new TFile( outputfile, "RECREATE" ) ;
    
    // loop over objects in the file
    // see root.cern.ch/drupal/content/how-read-objects-file
    // and ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch11s02.html
    // for more explaination of the below TIter and TKey algorithm
    TIter next( mscwfile->GetListOfKeys() );
    TKey* key ;
    
    // loop over objects in the TFile
    cout << endl;
    while( ( key = ( TKey* )next() ) )
    {
        string nam( key->GetName() ) ;
        printf( "examining: (%s,%s)\n", key->GetClassName(), nam.c_str() ) ;
        //if ( strcmp( key->GetClassName(), "TTree"     ) == 0 ) printf("   is TTree!\n");
        //if ( strcmp( treename           , nam.c_str() ) == 0 ) printf("   name matches target!\n");
        if( ( strcmp( key->GetClassName(), "TTree" ) == 0 ) && ( strcmp( treename           , nam.c_str() ) == 0 ) )
        {
            // then its our extraction tree
            printf( "  extracting %s %s\n", key->GetClassName(), nam.c_str() ) ;
            
            mscwfile->cd() ;
            TTree* mscwtree = ( TTree* ) mscwfile->Get( nam.c_str() ) ;
            char selectionstring[200] = "" ;
            //sprintf( selectionstring, " Entry$ >= %d && Entry$ <%d", startevent, startevent+nevents ) ;
            sprintf( selectionstring, " eventNumber >= %d && eventNumber <%d", startevent, startevent + nevents ) ;
            extractfile->cd() ;
            TTree* newtree = mscwtree->CopyTree( selectionstring ) ;
            newtree->Write() ;
        }
        else
        {
            // else its just a regular object we need to copy
            printf( "  copying %s %s\n", key->GetClassName(), nam.c_str() ) ;
            extractfile->cd() ;
            if( strcmp( key->GetClassName(), "TTree" ) == 0 )
            {
                TTree* oldtree = ( TTree* ) mscwfile->Get( nam.c_str() ) ;
                TTree* newtree = oldtree->CloneTree();
                newtree->Write() ;
            }
            else if( strcmp( key->GetClassName(), "VEvndispRunParameter" ) == 0 )
            {
                VEvndispRunParameter* obj  = ( VEvndispRunParameter* ) mscwfile->Get( nam.c_str() ) ;
                VEvndispRunParameter* obj2 = ( VEvndispRunParameter* ) obj->Clone() ;
                obj2->Write( nam.c_str() ) ;
            }
            else if( strcmp( key->GetClassName(), "VTableLookupRunParameter" ) == 0 )
            {
                VTableLookupRunParameter* obj  = ( VTableLookupRunParameter* ) mscwfile->Get( nam.c_str() ) ;
                VTableLookupRunParameter* obj2 = ( VTableLookupRunParameter* ) obj->Clone() ;
                obj2->Write( nam.c_str() ) ;
            }
            else if( strcmp( key->GetClassName(), "TDirectoryFile" ) == 0 )
            {
                TDirectoryFile* obj  = ( TDirectoryFile* ) mscwfile->Get( nam.c_str() ) ;
                CopyDir( obj ) ;
            }
            else
            {
                printf( "frogsExtractDatafile Error, unsupported copying of object '%s'\n", key->GetClassName() ) ;
                return 1 ;
            }
        }
        
    }
    
    // cleanup
    mscwfile->Close();
    extractfile->Close();
}


void CopyDir( TDirectory* source )
{
    //copy all objects and subdirs of directory source as a subdir of the current directory
    source->ls();
    TDirectory* savdir = gDirectory;
    TDirectory* adir = savdir->mkdir( source->GetName() );
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
        if( cl->InheritsFrom( TDirectory::Class() ) )
        {
            source->cd( key->GetName() );
            TDirectory* subdir = gDirectory;
            adir->cd();
            CopyDir( subdir );
            adir->cd();
        }
        else if( cl->InheritsFrom( TTree::Class() ) )
        {
            TTree* T = ( TTree* )source->Get( key->GetName() );
            adir->cd();
            TTree* newT = T->CloneTree( -1, "fast" );
            newT->Write();
        }
        else
        {
            source->cd();
            TObject* obj = key->ReadObj();
            adir->cd();
            obj->Write();
            delete obj;
        }
    }
    adir->SaveSelf( kTRUE );
    savdir->cd();
}
