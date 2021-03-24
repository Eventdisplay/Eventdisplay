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
    readConfigFile( iConfigFile );

    // tree with cut results for each event

    runNumber = 0;
    eventNumber = 0;
    ArrayPointing_Azimuth = 0.;
    ArrayPointing_Elevation = 0.;
    MCze = 0.;
    MCaz = 0.;
    MCxoff = 0.;
    MCyoff = 0.;
    MCe0 = 0.;
    NImages = 0;
    ImgSel = 0;
    ErecS = 0.;
    dES = 0.;
    EChi2S = 0.;
    Az = 0.;
    Ze = 0.;
    Xoff = 0.;
    Yoff = 0.;
    Xoff_intersect = 0.;
    Yoff_intersect = 0.;
    Xoff_derot = 0.;
    Yoff_derot = 0.;
    Chi2 = 0.;
    Xcore = 0.;
    Ycore = 0.;
    MCxcore = 0.;
    MCycore = 0.;
    DispNImages = 0;
    MSCW = 0.;
    MSCL = 0.;
    SizeSecondMax = 0.;
    EmissionHeight = 0.;
    EmissionHeightChi2 = 0.;
    NTtype = 0;
    fCut_MVA = 0.;

    for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
    {
        R[i] = 0.;
        ES[i] = 0.;
        TtypeID[i] = 0;
        NImages_Ttype[i] = 0;
    }
    // not included in mini trees:
    // size
    // width
    // length
    // dist
    //

    fDL2DataTree = new TTree( "data", "event data (shortened version)" );
    fDL2DataTree->Branch( "MVA", &fCut_MVA, "MVA/F" );
}

/////////////////////////////////////////////////////////////////////////////////////

VDL2Writer::~VDL2Writer()
{
    for( unsigned int i = 0; i < fTMVA.size(); i++ )
    {
        if( fTMVA[i] ) delete fTMVA[i];
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
            if( temp == "OFFSETBIN" )
            {
                is_stream >> temp;
                dist_mean.push_back( temp );
                float temp_f;
                is_stream >> temp_f;
                dist_min.push_back( temp_f );
                is_stream >> temp_f;
                dist_max.push_back( temp_f );
            }
            else if( temp == "SIMULATIONFILE_DATA" )
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
                for( unsigned int i = 0; i < dist_mean.size(); i++ )
                {
                    fTMVAWeightFile.push_back( gSystem->ExpandPathName( (iWeightFileDirectory+dist_mean[i]).c_str() ) );
                    // check if path name is complete
                    // note: this method returns FALSE if one **can** access the file
                    if( gSystem->AccessPathName( fTMVAWeightFile.back().c_str() ) )
                    {
                        fTMVAWeightFile.back() = VGlobalRunParameter::getDirectory_EVNDISPAnaData() 
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
    cout << "\t number of data events in source tree: " << d_nentries << endl;

    unsigned int i_dist_bin = 0;
    float dist = 0.;
    
    // loop over all events
    for( Long64_t i = 0; i < d_nentries; i++ )
    {
        d->GetEntry( i );

        if( d->Xoff < -50. || d->Yoff  < -50. ) continue;

        dist = sqrt( d->Xoff*d->Xoff + d->Yoff*d->Yoff );

        i_dist_bin = 999;
        for( unsigned int b = 0; b < dist_min.size(); b++ )
        {
            if( dist >= dist_min[b] && dist < dist_max[b] )
            {
                i_dist_bin = b;
                break;
            }
        }

        if( i_dist_bin < fTMVA.size() )
        {
           fillEventDataTree( fTMVA[i_dist_bin]->evaluate() );
        }
           
    }
    if( fDL2DataTree )
    {
        cout << "\t number of events in DL2 tree: " << fDL2DataTree->GetEntries() << endl;
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
 *  fDL2DataTree->Draw("MVA", "Class==5" );
 *
 *  2. Events passing gamma/hadron separation cut and not direction cut
 *  fDL2DataTree->Draw("MVA", "Class==0" );
 *
 *  3. Events before gamma/hadron separation cut and before direction cut
 *   fDL2DataTree->Draw("MVA", "Class==0||Class==7||Class==5", "");
 *
 */
void VDL2Writer::fillEventDataTree( float iMVA )
{
      if( fDL2DataTree )
      {
          fCut_MVA = iMVA;
          fDL2DataTree->Fill();
      }
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
