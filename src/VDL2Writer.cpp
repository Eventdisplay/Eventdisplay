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
    fTMVAEvaluator = 0;
    fTMVAEnergyStepSize = 0.;
    fTMVA_MVAMethodCounter = 0;
    fTMVAWeightFileIndex_Emin = 0;
    fTMVAWeightFileIndex_Emax = 0;
    fTMVAWeightFileIndex_Zmin = 0;
    fTMVAWeightFileIndex_Zmax = 0;

    readConfigFile( iConfigFile );

    // tree with cut results for each event
    fCut_MVA = 0.;
    fEventTreeCuts = new TTree( "fEventTreeCuts", "event cuts" );
    fEventTreeCuts->Branch( "MVA", &fCut_MVA, "MVA/F" );
}

/////////////////////////////////////////////////////////////////////////////////////

VDL2Writer::~VDL2Writer()
{
    if( fTMVAEvaluator )
    {
        delete fTMVAEvaluator;
    }
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
            if( temp == "SIMULATIONFILE_DATA" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fdatafile;
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
                fTMVAWeightFile = gSystem->ExpandPathName( iWeightFileDirectory.c_str() );
                // check if path name is complete
                // note: this method returns FALSE if one **can** access the file
                if( gSystem->AccessPathName( fTMVAWeightFile.c_str() ) )
                {
                    fTMVAWeightFile = VGlobalRunParameter::getDirectory_EVNDISPAnaData() + "/" + fTMVAWeightFile;
                    if( gSystem->AccessPathName( fTMVAWeightFile.c_str() ) )
                    {
                        cout << "error, weight file directory not found: ";
                        cout << fTMVAWeightFile << endl;
                        cout << "exiting..." << endl;
                        exit( EXIT_FAILURE );
                    }
                }
                fTMVAWeightFile += "/" + iWeightFileName;
            }
        }
    }
    return true;
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
     
    if( !initializeTMVAEvaluator( d ) )
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
    cout << "\t total number of data events: " << d_nentries << endl;
    
    // loop over all events
    for( Long64_t i = 0; i < d_nentries; i++ )
    {
        d->GetEntry( i );

        fTMVAEvaluator->evaluate();
        fillEventDataTree( fTMVAEvaluator->getTMVA_EvaluationResult() );
    }
    
    return true;
}


/* 
 * fill tree with cut numbers and MVA values per event
 *
 * 1    hhEcutTrigger_R (color: 1, marker: 20)
 * 2    hhEcutFiducialArea_R (color: 2, marker: 20)
 * 3    hhEcutStereoQuality_R (color: 3, marker: 20)
 * 4    hhEcutTelType_R (color: 4, marker: 20)
 * 5    hhEcutDirection_R (color: 5, marker: 20)
 * 6    hhEcutEnergyReconstruction_R (color: 6, marker: 20)
 * 7    hhEcutGammaHadron_R (color: 7, marker: 20)
 *
 *  1. Events passing gamma/hadron separation cut and direction cut
 *  fEventTreeCuts->Draw("MVA", "Class==5" );
 *
 *  2. Events passing gamma/hadron separation cut and not direction cut
 *  fEventTreeCuts->Draw("MVA", "Class==0" );
 *
 *  3. Events before gamma/hadron separation cut and before direction cut
 *   fEventTreeCuts->Draw("MVA", "Class==0||Class==7||Class==5", "");
 *
 */
void VDL2Writer::fillEventDataTree( float iMVA )
{
      if( fEventTreeCuts )
      {
          fCut_MVA = iMVA;
          fEventTreeCuts->Fill();
      }
}

bool VDL2Writer::initializeTMVAEvaluator( CData *d )
{
    fTMVAEvaluator = new VTMVAEvaluator();
    fTMVAEvaluator->setTMVAMethod( fTMVA_MVAMethod, fTMVA_MVAMethodCounter );
   
    if( !fTMVAEvaluator->initializeWeightFiles( 
            fTMVAWeightFile, 
            fTMVAWeightFileIndex_Emin, fTMVAWeightFileIndex_Emax,
            fTMVAWeightFileIndex_Zmin, fTMVAWeightFileIndex_Zmax, 
            fTMVAEnergyStepSize, fInstrumentEpoch ) )
    {
        cout << "initTMVAEvaluator: error while initializing TMVA weight files" << endl;
        cout << "exiting... " << endl;
        exit( EXIT_FAILURE );
    }

    fTMVAEvaluator->initializeDataStrutures( d );

    return !fTMVAEvaluator->IsZombie();
}
