/*! \class VPEReader

     reading class for pe files


*/

#include <VPEReader.h>

VPEReader::VPEReader( string isourcefile, vector< unsigned int > iTeltoAna, VDetectorGeometry* iD, bool iDebug )
{
    fDebug = iDebug;
    setEventStatus( 1 );
    if( fDebug )
    {
        cout << "VPEReader::VPEReader" << endl;
    }
    
    degrad = 45. / atan( 1. );
    
    fSourceFileName = isourcefile;
    
    fTeltoAna = iTeltoAna;
    
    // set up detector geometry
    fDetectorGeometry = iD;
    if( !fDetectorGeometry )
    {
        cout << "VPEReader::VPEReader: error: no valid detector geometry" << endl;
        cout << "...exiting" << endl;
        exit( -1 );
    }
    fNTelescopes = fDetectorGeometry->getNumTelescopes();
    fArrayTrigger = false;
    fNLocalTrigger = 0;
    for( unsigned int i = 0; i < fNTelescopes; i++ )
    {
        fDetectorGeometry->setTelID( i );
        fNChannel.push_back( fDetectorGeometry->getNumChannels() );
        fLocalTrigger.push_back( false );
    }
    
    // pe data is always MC data
    fMC = true;
    
    // open source file and init tree
    init();
}


VPEReader::~VPEReader()
{

}


/*!

    call this once per run (open source file, init tree, ...)
*/
bool VPEReader::init()
{
    if( fDebug )
    {
        cout << " VPEReader::init() " << endl;
    }
    
    // init data vectors
    for( unsigned int i = 0; i < fNTelescopes; i++ )
    {
        valarray< double > i_temp( 0., fNChannel[i] );
        vector< bool > i_tempB( fNChannel[i], true );
        vector< bool > i_tempF( fNChannel[i], false );
        fSums.push_back( i_temp );
        vector< valarray< double > > i_temp_VV;
        for( unsigned int t = 0; t < VDST_MAXTIMINGLEVELS; t++ )
        {
            i_temp_VV.push_back( i_temp );
        }
        fTracePulseTiming.push_back( i_temp_VV );
        fFullHitVec.push_back( i_tempB );
        fFullTrigVec.push_back( i_tempB );
        fHiLo.push_back( i_tempF );
        fNumberofFullTrigger.push_back( 0 );
        fTelAzimuth.push_back( 0. );
        fTelElevation.push_back( 0. );
    }
    
    // open file
    fPE_file = new TFile( fSourceFileName.c_str() );
    if( fPE_file->IsZombie() )
    {
        cout << "Error reading PE source file: " << fSourceFileName << endl;
        cout << "... exiting ... " << endl;
        exit( -1 );
    }
    
    // set up all trees
    cout << "reading PE file: checking trees..." << endl;
    // get single pe tree vector
    
    char hname[2000];
    bool bFirst = true;
    map< unsigned int, unsigned int > iTelGrisu = fDetectorGeometry->getTelIDGrisu();
    for( unsigned int i = 0; i < fNTelescopes; i++ )
    {
        // check if this telescope should be analyzed
        bool bAna = false;
        for( unsigned int t = 0; t < fTeltoAna.size(); t++ ) if( fTeltoAna[t] == i )
            {
                bAna = true;
            }
            
        if( bAna )
        {
            unsigned int iTelID = i;
            if( iTelGrisu.find( i ) != iTelGrisu.end() )
            {
                iTelID = iTelGrisu[i];
            }
            sprintf( hname, "tPE_T%d", iTelID );
            
            TTree* t = ( TTree* )fPE_file->Get( hname );
            if( !t )
            {
                cout << "VPEReader::init(): error reading pe tree " << hname << endl;
                cout << "...exiting" << endl;
                exit( -1 );
            }
            else
            {
                cout << "\t found " << hname << " for telescope " << i + 1 << endl;
            }
            fPE_Tree.push_back( new VPETree( t ) );
            fPE_Tree.back()->initialize( bFirst );
        }
        else
        {
            fPE_Tree.push_back( 0 );
        }
    }
    
    fPE_eventtype = 0;
    fPE_treeEvent = 0;
    fPE_primary = 0;
    fPE_az = 0.;
    fPE_ze = 0.;
    
    return true;
}


bool VPEReader::setTelescopeID( unsigned int iTelID )
{
    if( iTelID < fNTelescopes )
    {
        fTelID = iTelID;
    }
    else
    {
        return false;
    }
    
    return true;
}


bool VPEReader::getNextEvent()
{
    if( fDebug )
    {
        cout << "VPEReader::getNextEvent()" << endl;
    }
    
    for( unsigned int i = 0; i < fSums.size(); i++ )
    {
        fSums[i] = 0.;
        for( unsigned int j = 0; j < fTracePulseTiming[i].size(); j++ )
        {
            fTracePulseTiming[i][j] = 0.;
        }
        fLocalTrigger[i] = false;
    }
    
    bool bLocalTrigger = false;
    fNLocalTrigger = 0;
    fArrayTrigger = false;
    
    // loop over all trees
    for( unsigned int i = 0; i < fPE_Tree.size(); i++ )
    {
        if( !fPE_Tree[i] )
        {
            continue;
        }
        
        if( fPE_treeEvent >= ( unsigned int )fPE_Tree[i]->getEntries() )
        {
            setEventStatus( 999 );
            return false;
        }
        bLocalTrigger = false;
        
        fPE_Tree[i]->getEntry( fPE_treeEvent );
        
        // get MC info from first tree
        
        if( fPE_Tree[i]->isFull() )
        {
            fPE_eventnumber = fPE_Tree[i]->fEventNumber + 1;
            fPE_primary = ( int )fPE_Tree[i]->fPrimaryType;
            fPE_energy = fPE_Tree[i]->fPrimaryEnergy;
            fPE_xcore =  fPE_Tree[i]->fXcore;
            fPE_ycore =  fPE_Tree[i]->fYcore;
            fPE_xcos = fPE_Tree[i]->fXcos;
            fPE_ycos = fPE_Tree[i]->fYcos;
            fPE_az   = atan2( fPE_xcos, fPE_ycos );
            fPE_ze   = 1. - ( fPE_xcos * fPE_xcos + fPE_ycos * fPE_ycos );
            if( fPE_ze < 0. )
            {
                fPE_ze = 0.;
            }
            fPE_ze = acos( sqrt( fPE_ze ) );
            if( TMath::Abs( fPE_ze ) < 0.1 )
            {
                fPE_az += TMath::Pi();
            }
            fPE_Tel_xoff = -1.*fPE_Tree[i]->fXsource;
            fPE_Tel_yoff = -1.*fPE_Tree[i]->fYsource;
            for( unsigned int t = 0; t < fTelElevation.size(); t++ )
            {
                fTelElevation[t] =  90. - fPE_ze * degrad;
            }
            for( unsigned int t = 0; t < fTelAzimuth.size(); t++ )
            {
                fTelAzimuth[t] = fPE_az * degrad;
            }
            // add wobble offset
            double az = 0.;
            double ze = 0.;
            VAstronometry::vlaDtp2s(  -1.*fPE_Tel_xoff * TMath::DegToRad(),
                                      -1.*fPE_Tel_yoff * TMath::DegToRad(),
                                      fPE_az,
                                      0.5*TMath::Pi()-ze,
                                      &az, &ze );

            fPE_ze = 0.5*TMath::Pi()-ze;
            fPE_az = az;
        }
        
        if( fPE_Tree[i]->v_f_time && fPE_Tree[i]->v_f_ID )
        {
            for( unsigned int v = 0; v < fPE_Tree[i]->v_f_ID->size(); v++ )
            {
                if( fPE_Tree[i]->v_f_ID->at( v ) < fSums[i].size() )
                {
                    fSums[i][fPE_Tree[i]->v_f_ID->at( v )]++;
                    bLocalTrigger = true;
                }
            }
        }
        fLocalTrigger[i] = false;
        if( bLocalTrigger )
        {
            fLocalTrigger[i] = true;
            fNLocalTrigger++;
        }
        
        if( fDebug )
        {
            if( fPE_Tree[i]->isFull() )
            {
                cout << "tree " << i << "\t" <<  "entry " << fPE_treeEvent << "\t, event number: " << fPE_Tree[i]->fEventNumber << endl;
                cout << "\t primary energy [TeV]: " << fPE_Tree[i]->fPrimaryEnergy << endl;
                cout << "\t core position [m]: " << fPE_Tree[i]->fXcore << ", " << fPE_Tree[i]->fYcore << endl;
                cout << "\t source direction and offsets: " << fPE_Tree[i]->fXcos << "\t" << fPE_Tree[i]->fYcos << "\t" << fPE_Tree[i]->fXsource << "\t" << fPE_Tree[i]->fYsource << endl;
                if( fPE_Tree[i]->v_f_time && fPE_Tree[i]->v_f_ID )
                {
                    cout << "\t vector sizes: " << fPE_Tree[i]->v_f_time->size() << "\t" << fPE_Tree[i]->v_f_ID->size() << endl;
                    if( fPE_Tree[i]->v_f_time->size() == fPE_Tree[i]->v_f_ID->size() )
                    {
                        for( unsigned int v = 0; v < fPE_Tree[i]->v_f_time->size(); v++ )
                        {
                            cout << "\t\t" << fPE_Tree[i]->v_f_ID->at( v ) << "\t" << fPE_Tree[i]->v_f_time->at( v ) << endl;
                        }
                        cout << endl;
                    }
                }
                cout << "sums: " << endl;
                for( unsigned int j = 0; j < fSums[i].size(); j++ )
                {
                    if( fSums[i][j] > 0. )
                    {
                        cout << "\t" << i << "\t" << j << "\t" << fSums[i][j] << endl;
                    }
                }
                cout << endl;
            }
        }
    }
    
    // fill data vectors
    for( unsigned int i = 0; i < fNTelescopes; i++ )
    {
        for( unsigned int j = 0; j < fNChannel[i]; j++ )
        {
            for( unsigned int k = 0; k < fTracePulseTiming[i].size(); k++ )
            {
                fTracePulseTiming[i][j][k] = 0.;
            }
            fHiLo[i][j] = ( bool )0;
        }
    }
    
    // get local trigger time
    if( fMC )
    {
        fLTtime.clear();
        for( unsigned int i = 0; i < fNTelescopes; i++ )
        {
            fLTtime.push_back( 0 );
            fLDTtime.push_back( 0 );
        }
    }
    
    // array trigger requires simply 2 telescope with non-zero size
    if( fNLocalTrigger > 1 )
    {
        fArrayTrigger = true;
    }
    else
    {
        fArrayTrigger = false;
    }
    // single telescope analysis
    if( fNTelescopes == 1 && fNLocalTrigger == 1 )
    {
        fArrayTrigger = true;
    }
    
    // increment tree event number
    fPE_treeEvent++;
    
    return true;
}


std::pair<bool, uint32_t> VPEReader::getChannelHitIndex( uint32_t hit )
{
    if( hit < fSums[fTelID].size() )
    {
        return std::make_pair( true, hit );
    }
    return std::make_pair( false, ( uint32_t ) 0 );
}


uint32_t VPEReader::getHitID( uint32_t i )
{
    if( i < fSums[fTelID].size() )
    {
        return i;
    }
    return 0;
}


void VPEReader::setTrigger( vector<bool> iImage, vector<bool> iBorder )
{
    if( fFullTrigVec[fTelID].size() != iImage.size() || fFullTrigVec[fTelID].size() != iBorder.size() )
    {
        cout << "VPEReader::setTrigger error: trigger/image/border vectors have different sizes ";
        cout << fFullTrigVec[fTelID].size() << "\t" << iImage.size() << "\t" << iBorder.size() << endl;
    }
    fNumberofFullTrigger[fTelID] = 0;
    for( unsigned int i = 0; i < fFullTrigVec[fTelID].size(); i++ )
    {
        if( iImage[i] || iBorder[i] )
        {
            fFullTrigVec[fTelID][i] = true;
        }
        else
        {
            fFullTrigVec[fTelID][i] = false;
        }
        fNumberofFullTrigger[fTelID]++;
    }
}


vector< double > VPEReader::getTelElevation()
{
    return fTelElevation;
}


vector< double > VPEReader::getTelAzimuth()
{
    return fTelAzimuth;
}
