/*! \class VDLWriter
 *  \brief DL2 (MVA values) writer
 *
 */

#include "VDL2Writer.h"

/*!
 *
 */
VDL2Writer::VDL2Writer( string iConfigFile )
{
    fWriteDataTree = false;
    fWriteDL2Tree = false;
    fWriteDL3Tree = false;
    readConfigFile( iConfigFile );

    fDL2DataTree = new VDL2Tree( "DL2EventTree", "DL2 tree" );
}

/////////////////////////////////////////////////////////////////////////////////////

VDL2Writer::~VDL2Writer()
{
    for( unsigned int i = 0; i < fTMVA.size(); i++ )
    {
        if( fTMVA[i] ) delete fTMVA[i];
    }
}

bool VDL2Writer::writeDataTrees( TFile *iOutputFile, TChain *c )
{
    if( !iOutputFile ) return false;
    iOutputFile->cd();
    if( fWriteDataTree && c )
    {
        cout << "writing data tree to " << iOutputFile->GetName() << endl;
        c->Merge( iOutputFile, 0, "keep" );
    }
    if( fWriteDL2Tree && fDL2DataTree && fDL2DataTree->getDL2Tree() )
    {
        cout << "writing DL2 tree to " << iOutputFile->GetName() << endl;
        fDL2DataTree->getDL2Tree()->Write();
    }
    return true;
}


bool VDL2Writer::readConfigFile( string iConfigFile )
{
    ifstream is;
    is.open( iConfigFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error opening run parameter file ";
        cout << iConfigFile << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    string temp;
    cout << endl;
    cout << "run parameter file (" << iConfigFile << ")" << endl;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            is_stream >> temp;
            cout << is_line << endl;
            if( temp == "DATAFILE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fdatafile.push_back( temp );
                }
            }
            else if( temp == "DATATREE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fWriteDataTree = bool( atoi(temp.c_str() ) );
                }
            }
            else if( temp == "DL2TREE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fWriteDL2Tree = bool( atoi(temp.c_str() ) );
                }
            }
            else if( temp == "DL3TREE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> temp;
                    fWriteDL3Tree = bool( atoi(temp.c_str() ) );
                }
            }
            // TMVA values
            else if( temp == "TMVAPARAMETER" )
            {
                is_stream >> fInstrumentEpoch;
                is_stream >> fTMVA_MVAMethod;
                is_stream >> fTMVA_MVAMethodCounter;
                is_stream >> fTMVAWeightFileIndex_Emin;
                is_stream >> fTMVAWeightFileIndex_Emax;
                is_stream >> fTMVAWeightFileIndex_Zmin;
                is_stream >> fTMVAWeightFileIndex_Zmax;
                is_stream >> fTMVAEnergyStepSize;
                string iWeightFileDirectory;
                is_stream >> iWeightFileDirectory;
                string iWeightFileName;
                is_stream >> iWeightFileName;
                fTMVAWeightFile.push_back( gSystem->ExpandPathName( iWeightFileDirectory.c_str() ) );
                // check if path name is complete
                // note: this method returns FALSE if one **can** access the file
                if( gSystem->AccessPathName( fTMVAWeightFile.back().c_str() ) )
                {
                    VGlobalRunParameter iRunPara;
                    fTMVAWeightFile.back() = iRunPara.getDirectory_EVNDISPAnaData() 
                                            + "/" + fTMVAWeightFile.back();
                    if( gSystem->AccessPathName( fTMVAWeightFile.back().c_str() ) )
                    {
                        cout << "error, weight file directory not found: ";
                        cout << fTMVAWeightFile.back() << endl;
                        cout << "exiting..." << endl;
                        exit( EXIT_FAILURE );
                    }
                }
                fTMVAWeightFile.back() += "/" + iWeightFileName;
            }
        }
    }
    return true;
}

/*
 * get tmva wobble / distance bin
 * (if available)
*/
unsigned int VDL2Writer::get_tmva_distance_bin( float dist )
{
    if( dist_min.size() == 0 ) return 0;
    for( unsigned int b = 0; b < dist_min.size(); b++ )
    {
        if( dist >= dist_min[b] && dist < dist_max[b] )
        {
            return b;
        }
    }
    return 999;
}

/*
 *
 *  event loop 
 *
 */
bool VDL2Writer::fill( CData* d )
{
    // make sure that all data pointers exist
    if( !d )
    {
        cout << "VDL2Writer::fill error: no data tree" << endl;
        return false;
    }
     
    if( !initializeTMVAEvaluators( d ) )
    {
        return false;
    }
    
    //////////////////////////////////////////////////////////////////
    // print some run information
    cout << endl;
    cout << "DL2 Eventlist filling" << endl;
    
    ///////////////////////////////////////////////////////
    // get full data set and loop over all entries
    ///////////////////////////////////////////////////////
    Long64_t d_nentries = d->fChain->GetEntries();
    d_nentries = 1000;
    cout << "\t number of data events in source tree: " << d_nentries << endl;

    unsigned int i_tmva_dist_bin = 0;
    float i_tmva = 0.;
    int iCut_QC = 0;
    
    // loop over all events
    for( Long64_t i = 0; i < d_nentries; i++ )
    {
        d->GetEntry( i );

        if( d->Xoff < -50. || d->Yoff  < -50. || d->ErecS < 0. ) iCut_QC = -1;
        else iCut_QC = 0;

        i_tmva_dist_bin = get_tmva_distance_bin( sqrt( d->Xoff*d->Xoff 
                                                     + d->Yoff*d->Yoff ) );
        if( i_tmva_dist_bin < fTMVA.size() )
        {
           i_tmva = fTMVA[i_tmva_dist_bin]->evaluate();
        }
        else
        {
           i_tmva = -99.;
           iCut_QC = -2.;
        }
        fDL2DataTree->fillEvent( d, 
                                 i_tmva,
                                 iCut_QC );
    }
    if( fDL2DataTree->getDL2Tree() )
    {
        cout << "\t number of events in DL2 tree: " << fDL2DataTree->getDL2Tree()->GetEntries() << endl;
    }
    
    return true;
}

bool VDL2Writer::initializeTMVAEvaluators( CData *d )
{
    bool fZombie = true;
    for( unsigned int i = 0; i < fTMVAWeightFile.size(); i++ )
    {
        fTMVA.push_back( new VTMVA_eval_dist( fTMVA_MVAMethod,
                                              fTMVA_MVAMethodCounter,
                                              fTMVAWeightFile[i],
                                              fTMVAWeightFileIndex_Emin, fTMVAWeightFileIndex_Emax,
                                              fTMVAWeightFileIndex_Zmin, fTMVAWeightFileIndex_Zmax,
                                              fTMVAEnergyStepSize, fInstrumentEpoch,
                                              d ) );
        fZombie = fZombie || fTMVA.back()->fIsZombie;
    }
    return fZombie;
}

/////////////////////////////////////////////////////////
VTMVA_eval_dist::VTMVA_eval_dist( string i_fTMVA_MVAMethod,
                                  unsigned int i_fTMVA_MVAMethodCounter,
                                  string i_fTMVAWeightFile,
                                  unsigned int i_fTMVAWeightFileIndex_Emin, unsigned int i_fTMVAWeightFileIndex_Emax,
                                  unsigned int i_fTMVAWeightFileIndex_Zmin, unsigned int i_fTMVAWeightFileIndex_Zmax,
                                  double i_fTMVAEnergyStepSize, 
                                  string i_fInstrumentEpoch,
                                  CData *d )
{
    fTMVAEvaluator = new VTMVAEvaluator();
    fTMVAEvaluator->setTMVAMethod( i_fTMVA_MVAMethod, i_fTMVA_MVAMethodCounter );
   
    if( !fTMVAEvaluator->initializeWeightFiles( 
            i_fTMVAWeightFile, 
            i_fTMVAWeightFileIndex_Emin, i_fTMVAWeightFileIndex_Emax,
            i_fTMVAWeightFileIndex_Zmin, i_fTMVAWeightFileIndex_Zmax, 
            i_fTMVAEnergyStepSize, i_fInstrumentEpoch ) )
    {
        cout << "initTMVAEvaluator: error while initializing TMVA weight files" << endl;
        cout << "exiting... " << endl;
        exit( EXIT_FAILURE );
    }

    fTMVAEvaluator->initializeDataStrutures( d );

    fIsZombie = fTMVAEvaluator->IsZombie();
}

VTMVA_eval_dist::~VTMVA_eval_dist()
{
   if( fTMVAEvaluator )
   {
       delete fTMVAEvaluator;
   }
}

double VTMVA_eval_dist::evaluate()
{
    if( fIsZombie ) return -1.;
    if( !fTMVAEvaluator ) return -1.;

    fTMVAEvaluator->evaluate();
    return fTMVAEvaluator->getTMVA_EvaluationResult();
}
