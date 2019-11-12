#include <iostream>
#include <fstream>
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
bool isHighestCycle( TKey* tkey, TFile* tfile ) ;

int main( int argc, char* argv[] )
{
    if( argc != 3 )
    {
        printf( "Error, must give right arguments 'mergelistfile' 'outputfilename', exiting...\n" );
        exit( 1 );
    }
    
    // input file name
    char mergelistfile[200] = "" ;
    sprintf( mergelistfile, "%s", argv[1] ) ;
    
    // output file name
    char outputfile[1000] = "" ;
    sprintf( outputfile, "%s", argv[2] ) ;
    
    printf( " mergelistfile: '%s'\n", mergelistfile ) ;
    printf( " outputfile:    '%s'\n", outputfile ) ;
    
    // TTree name to merge
    vector<string> mergetrees ;
    mergetrees.push_back( "data" ) ;
    mergetrees.push_back( "frogspars" ) ;
    
    // read in list of files to merge
    cout << endl;
    vector<string> filestomerge ;
    ifstream fin( mergelistfile ) ;
    string   file_line ;
    while( getline( fin, file_line ) )
    {
        printf( "input merging file: '%s'\n", file_line.c_str() ) ;
        filestomerge.push_back( file_line ) ;
    }
    
    // open output TFile
    TFile* outf = new TFile( outputfile, "RECREATE" ) ;
    
    // check that each file has each TTree in 'mergetrees'
    cout << endl;
    TFile* inpf = 0 ;
    bool mtreeexists[mergetrees.size()] ;
    for( unsigned int i_file = 0 ; i_file < filestomerge.size() ; i_file++ )
    {
    
        // reset our 'which trees were found' array
        for( unsigned int i_mtree = 0 ; i_mtree < mergetrees.size() ; i_mtree++ )
        {
            mtreeexists[i_mtree] = false ;
        }
        
        // loop over each object in the TFile
        inpf = new TFile( filestomerge[i_file].c_str(), "READ" ) ;
        TIter next( inpf->GetListOfKeys() );
        TKey* key ;
        while( ( key = ( TKey* )next() ) )
        {
        
            // check if the object matches one of the TTrees in 'mergetrees'
            string objectname( key->GetName() ) ;
            string objectclass( key->GetClassName() ) ;
            for( unsigned int i_mtree = 0 ; i_mtree < mergetrees.size() ; i_mtree++ )
            {
                // if the object is a TTree, and has a name from 'mergetrees', then record that it exists
                if( objectname.compare( mergetrees[i_mtree] ) == 0 && objectclass.compare( "TTree" ) == 0 )
                {
                    if( ! mtreeexists[i_mtree] )
                    {
                        printf( "check: found merging target TTree '%-10s' in file '%-s'\n", mergetrees[i_mtree].c_str(), filestomerge[i_file].c_str() );
                    }
                    mtreeexists[i_mtree] = true ;
                }
            }
        }
        
        // check that all trees were found
        for( unsigned int i_mtree = 0 ; i_mtree < mergetrees.size() ; i_mtree++ )
        {
        
            // if a tree is missing, throw error and exit
            if( ! mtreeexists[i_mtree] )
            {
                printf( "Error, could not find TTree '%s' in file '%s', which is needed for merging, exiting...\n", mergetrees[i_mtree].c_str(), filestomerge[i_file].c_str() );
                exit( 1 );
            }
        }
        
    }
    
    
    // loop over objects in the file
    // see root.cern.ch/drupal/content/how-read-objects-file
    // and ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch11s02.html
    // for more explaination of the below TIter and TKey algorithm
    inpf = new TFile( filestomerge[0].c_str(), "READ" ) ;
    TIter next( inpf->GetListOfKeys() );
    TKey* key ;
    bool skipflag = false ;
    cout << endl << "Initial copy of not-merge-target objects" << endl;
    while( ( key = ( TKey* )next() ) )
    {
    
        string objectname( key->GetName() ) ;
        string objectclass( key->GetClassName() ) ;
        int    objectcycle = ( int ) key->GetCycle() ;
        printf( "examining %s '%s' (cycle=%d) ...\n", objectclass.c_str(), objectname.c_str(), objectcycle ) ;
        
        // check if the current key is a TTree on the
        skipflag = false ;
        for( unsigned int i_mtree = 0 ; i_mtree < mergetrees.size() ; i_mtree++ )
        {
            if( objectname.compare( mergetrees[i_mtree] ) == 0 && objectclass.compare( "TTree" ) == 0 )
            {
                skipflag = true ;
                break ;
            }
        }
        if( skipflag )
        {
            printf( " -> skipping target TTree '%s' (merging is handled last)\n", objectname.c_str() ) ;
            continue ;
        }
        
        // only keep highest cycle number for each key, skip all older cycles
        if( ! isHighestCycle( key, inpf ) )
        {
            printf( " -> skipping target %s '%s' (%d) because it is an older version of an existing key\n", objectclass.c_str(), objectname.c_str(), objectcycle );
            continue ;
        }
        
        // copy wanted objects
        outf->cd() ;
        if( objectclass.compare( "TTree" ) == 0 )
        {
            printf( " -> copying %s '%s'\n", objectclass.c_str(), objectname.c_str() ) ;
            TTree* oldtree = ( TTree* ) inpf->Get( objectname.c_str() ) ;
            TTree* newtree = oldtree->CloneTree();
            newtree->Write() ;
        }
        else if( objectclass.compare( "VEvndispRunParameter" ) == 0 )
        {
            printf( " -> copying %s '%s'\n", objectclass.c_str(), objectname.c_str() ) ;
            VEvndispRunParameter* obj  = ( VEvndispRunParameter* ) inpf->Get( objectname.c_str() ) ;
            VEvndispRunParameter* obj2 = ( VEvndispRunParameter* ) obj->Clone() ;
            obj2->Write( objectname.c_str() ) ;
        }
        else if( objectclass.compare( "VTableLookupRunParameter" ) == 0 )
        {
            printf( " -> copying %s '%s'\n", objectclass.c_str(), objectname.c_str() ) ;
            VTableLookupRunParameter* obj  = ( VTableLookupRunParameter* ) inpf->Get( objectname.c_str() ) ;
            VTableLookupRunParameter* obj2 = ( VTableLookupRunParameter* ) obj->Clone() ;
            obj2->Write( objectname.c_str() ) ;
        }
        else if( objectclass.compare( "TDirectoryFile" ) == 0 )
        {
            printf( " -> copying %s '%s'\n", objectclass.c_str(), objectname.c_str() ) ;
            TDirectoryFile* obj  = ( TDirectoryFile* ) inpf->Get( objectname.c_str() ) ;
            CopyDir( obj ) ;
        }
        else
        {
            printf( " -> error, can't copy %s '%s'\n", objectclass.c_str(), objectname.c_str() ) ;
        }
        
    }
    inpf->Close() ;
    
    // merge each tree
    cout << endl << "Merging Target Objects" << endl;
    for( unsigned int i_mtree = 0 ; i_mtree < mergetrees.size() ; i_mtree++ )
    {
        cout << "Merging TTree '" << mergetrees[i_mtree] << endl;
        
        // loop over list of files, adding them to tchain
        TChain* dataChain = new TChain() ;
        for( unsigned int i_file = 0 ; i_file < filestomerge.size() ; i_file++ )
        {
            // open the file to merge
            //printf("    looking for ttree '%s' in file %s\n", mergetrees[i_mtree].c_str(), filestomerge[i_file].c_str() ) ;
            //printf("adding file '%s'\n", filestomerge[i_file].c_str() ) ;
            //char mergefile[500] = "" ;
            inpf = new TFile( filestomerge[i_file].c_str(), "READ" ) ;
            TIter next2( inpf->GetListOfKeys() );
            TKey* key2 ;
            
            // loop over objects, to make sure the ttree 'mergetrees[i_mtree]' is there
            while( ( key2 = ( TKey* )next2() ) )
            {
                string objectname( key2->GetName() ) ;
                string objectclass( key2->GetClassName() ) ;
                int    objectcycle = ( int ) key2->GetCycle()       ;
                
                // only keep highest cycle number for each key, skip all older cycles
                if( ! isHighestCycle( key2, inpf ) )
                {
                    printf( " -> skipping target %s '%s' (%d) because it is an older version of an existing key\n", objectclass.c_str(), objectname.c_str(), objectcycle );
                    continue ;
                }
                
                // when we find the target tree, add it to the dataChain
                // and stop looping over objects in this file
                if( objectname.compare( mergetrees[i_mtree] ) == 0 && objectclass.compare( "TTree" ) == 0 )
                {
                    char fulltreename[200] = "" ;
                    //printf("   -> merging tree '%s' from file '%s'\n", mergetrees[i_mtree].c_str(), filestomerge[i_file] ) ;
                    //cout << "extracting from '" << filestomerge[i_file] << endl;
                    sprintf( fulltreename, "%s/%s;%d", filestomerge[i_file].c_str(), mergetrees[i_mtree].c_str(), objectcycle ) ;
                    printf( " -> merging the TTree '%s' to TChain\n", fulltreename ) ;
                    dataChain->Add( fulltreename ) ;
                    break;
                }
                
            }
            inpf->Close() ;
        }
        
        // write TTree to file
        cout << " -> writing merged TTree '" << mergetrees[i_mtree].c_str() << "'" << endl;
        outf->cd() ;
        cout << "ls before:" << endl;
        outf->ls() ;
        TTree* dataTree = ( TTree* )( dataChain->CloneTree() ) ;
        dataTree->Write( mergetrees[i_mtree].c_str() ) ;
        cout << "ls after:" << endl;
        outf->ls() ;
        
    }
    
    outf->Close() ;
    cout << endl;
    cout << "Wrote file '" << outputfile << "', exiting..." << endl;
}


// copyies everything in 'source' to whatever gDirectory is selected
void CopyDir( TDirectory* source )
{
    //copy all objects and subdirs of directory source as a subdir of the current directory
    //source->ls();
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


// loop over all objects in file,
// if object at 'tkey' is the latest version
// of itself, return true, else return false
// warning: will also report true if 'tkey' doesn't exist in 'tfile',
//   so make sure its there before you use it
bool isHighestCycle( TKey* tkey, TFile* tfile )
{
    bool debug = false ;
    string objname( tkey->GetName() ) ;
    string objclass( tkey->GetClassName() ) ;
    int    objcycle = tkey->GetCycle()   ;
    if( debug )
    {
        cout << "isHighestCycle( " << objclass << " '" << objname << "' (" << objcycle << ")" << endl;
    }
    if( debug )
    {
        printf( "  obj:  %-25s '%-20s' %2d\n", objclass.c_str(), objname.c_str(), objcycle ) ;
    }
    TIter next( tfile->GetListOfKeys() );
    TKey* key ;
    
    while( ( key = ( TKey* )next() ) )
    {
        string itername( key->GetName() ) ;
        string iterclass( key->GetClassName() ) ;
        int    itercycle = key->GetCycle()       ;
        if( debug )
        {
            printf( "  iter: %-25s '%-20s' %2d\n", iterclass.c_str(), itername.c_str(), itercycle ) ;
        }
        
        if( objname.compare( itername ) == 0 && objclass.compare( iterclass ) == 0 )
        {
            if( objcycle < itercycle )
            {
                if( debug )
                {
                    cout << "    object is older than this iter" << endl;
                }
                if( debug )
                {
                    cout << "  isHighestCycle()=false!" << endl;
                }
                return false ;
            }
            else
            {
                if( debug )
                {
                    cout << "    object is exact match to or newer than this iter" << endl;
                }
            }
        }
        
    }
    
    if( debug )
    {
        cout << "  isHighestCycle()=true!" << endl;
    }
    return true ;
}


