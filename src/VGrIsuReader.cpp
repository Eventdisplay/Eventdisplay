/*! \class VGrIsuReader
    \brief reading of GrIsu output files

     - fdefaultPed must be the same value as in GrIsu for the trace library mode

    \attention
       - assume that 'P' lines are always before 'S' or 'R' lines

       - fillBackgroundfromTraceLibrary should be checked

*/

#include <VGrIsuReader.h>

/*!
  \param iDetGeo pointer to camera geometry
  \param i_sourcefile source file name
  \param i_sumwindow vector with summation windows for all telescopes (for pedestal calculation)
  \param i_telnumberoffset offset in telescope numbering (different in for grisu versions < 3.0)

  - default pedestal values are hard coded: fdefaultPed, fdefaultPedvars
*/

VGrIsuReader::VGrIsuReader( VDetectorGeometry* iDetGeo, unsigned int intel, string iExtPed, vector< int > i_sumwindow, bool iDebug, int igrisuseed, double ifadcscale )
{
    fDebug = iDebug;
    fDetectorGeo = iDetGeo;
    fExternalPedFile = iExtPed;
    fSumWindow = i_sumwindow;
    setEventStatus( 1 );
    fPedestalMode = false;
    fNPedestalEvents = 0;
    fDataFormat = "grisu";
    fLastWithData = false;
    fSourceFileName = "";
    fTraceFileName = "";
    fExternalPedFile = iExtPed;
    setEventStatus( 1 );
    fMultiGrIsuReader = false;
    fIgnoreCFGVersions = false;
    fNoiseTraceStart = 0;
    bNoiseTraceFilled = false;
    fdefaultPedvars = 3.;                         // should be greater than zero, otherwise VAnalyzer declares pixel dead
    
    degrad = 45. / atan( 1. );
    
    // constants
    fTelescopeID = 0;
    fTelNumberOffset = 1;
    fSampleOffset = 0;
    fMCNdead = 0;
    fMCNdeadboard = 0;
    fEventNumber = 0;
    fEventNumberPreviousEvent = 0;
    fLocalTrigger.assign( 255, false );
    fFADCScale = ifadcscale;
    fNTel =  intel;
    // standard eventtype is IGNORE
    fEventType = 0;
    
    // attention! hard coded numbers for pedestals
    // this must be the same as in the grisu pilot file (in the case of using the tracelibrary)
    fdefaultPed = fDetectorGeo->getDefaultPedestal();
    fdefaultPedUI = ( uint8_t )fdefaultPed;
    
    // initialize data vectors
    initVectors();
    
    // random generator (needed for background generation)
    int iseed = igrisuseed;
    fRandomGen = new TRandom3( iseed );
    cout << "\t VGrIsuReader::VGrIsuReader(): random number generator seed: " << fRandomGen->GetSeed() << endl;
    
    // read the pedestals
    readPeds();
}


VGrIsuReader::VGrIsuReader( VDetectorGeometry* iDetGeo, unsigned int intel, string i_sourcefile, vector< int > i_sumwindow, int i_telnumberoffset, int i_sampleoffset, double ifadcscale, bool iDebug, int igrisuseed, string iExtPed, bool iIgnoreCFGFiles )
{
    fDebug = iDebug;
    if( fDebug )
    {
        cout << "VGrIsuReader::VGrIsuReader" << endl;
    }
    fPedestalMode = false;
    fNPedestalEvents = 0;
    fDetectorGeo = iDetGeo;
    fDataFormat = "grisu";
    fLastWithData = false;
    fSourceFileName = i_sourcefile;
    fTraceFileName = "";
    fExternalPedFile = iExtPed;
    setEventStatus( 1 );
    fMultiGrIsuReader = false;
    fIgnoreCFGVersions = iIgnoreCFGFiles;
    fNoiseTraceStart = 0;
    bNoiseTraceFilled = false;
    
    degrad = 45. / atan( 1. );
    
    // constants
    fTelescopeID = 0;
    fTelNumberOffset = i_telnumberoffset;
    fSampleOffset = i_sampleoffset;
    fSumWindow = i_sumwindow;
    fMCNdead = 0;
    fMCNdeadboard = 0;
    fEventNumber = 0;
    fEventNumberPreviousEvent = 0;
    fLocalTrigger.assign( 255, false );
    fFADCScale = ifadcscale;
    fNTel =  intel;
    // standard eventtype is IGNORE
    fEventType = 0;
    
    // attention! hard coded numbers for pedestals
    // this must be the same as in the grisu pilot file (in the case of using the tracelibrary)
    fdefaultPed = fDetectorGeo->getDefaultPedestal();
    fdefaultPedUI = ( uint8_t )fdefaultPed;
    fdefaultPedvars = 3.;                         // should be greater than zero, otherwise VAnalyzer declares pixel dead
    
    // initialize data vectors
    initVectors();
    
    // open source data file
    openDataFile( false );
    
    // open files with pedestals
    openDataFile( true );
    
    // random generator (needed for background generation)
    int iseed = igrisuseed;
    //   iseed = atoi( fSourceFileName.substr( fSourceFileName.rfind( "/" ) + 4, 6 ).c_str() );
    fRandomGen = new TRandom3( iseed );
    cout << "\t VGrIsuReader::VGrIsuReader(): random number generator seed: " << fRandomGen->GetSeed() << endl;
    
    // read some camera infos
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        readCamera( i );
    }
    
    // get and test grisu version
    fGrisuVersion = getGrisuVersion();
    if( ( fDetectorGeo->getGrIsuVersion() >= 411 || fGrisuVersion >= 411 ) && fGrisuVersion != fDetectorGeo->getGrIsuVersion() && fGrisuVersion != 0 )
    {
        if( fIgnoreCFGVersions )
        {
            cout << "VGrIsuReader warning: different version of GrIsu simulation file and configuration file " << endl;
            cout << "\t GrIsu simulation file: " << fGrisuVersion << endl;
            cout << "\t configuration file: " << fDetectorGeo->getGrIsuVersion() << endl;
        }
        else
        {
            cout << "VGrIsuReader error: incompatible version of GrIsu simulation file and configuration file " << endl;
            cout << "\t GrIsu simulation file: " << fGrisuVersion << endl;
            cout << "\t configuration file: " << fDetectorGeo->getGrIsuVersion() << endl;
            cout << "...exiting" << endl;
            exit( 0 );
        }
    }
    if( fGrisuVersion > 0 )
    {
        cout << "\t GrIsu simulations uses configuration file with version number " << fGrisuVersion << endl;
    }
    
    // read the pedestals
    readPeds();
}


VGrIsuReader::~VGrIsuReader()
{
    // if source file was zipped, zip it again
    if( fBZipped )
    {
        cout << "VGrIsuReader::~VGrIsuReader() info: zipping source file" << endl;
        string i_com = "gzip -v -f  ";
        i_com += fSourceFileName;
        gSystem->Exec( i_com.c_str() );
    }
    else if( fBZipped2 )
    {
        cout << "VGrIsuReader::~VGrIsuReader() info: removing unzipped source file" << endl;
        string i_com = "rm -v -f  ";
        i_com += fSourceFileName;
        gSystem->Exec( i_com.c_str() );
    }
}


void VGrIsuReader::closeDataFile( bool iPeds )
{
    if( iPeds && is_ped )
    {
        is_ped.close();
    }
}


void VGrIsuReader::openDataFile( bool iPeds )
{
    if( !iPeds )
    {
        // check if sourcefile is zipped
        fBZipped = false;
        fBZipped2 = false;
        if( fSourceFileName.size() > 4 &&  fSourceFileName.size() - fSourceFileName.rfind( ".gz" ) == 3 )
        {
            cout << "VGrIsuReader::VGrIsuReader info: found gzipped file, unzipping..." << endl;
            string i_com = "gunzip -v -f  ";
            i_com += fSourceFileName;
            gSystem->Exec( i_com.c_str() );
            fSourceFileName = fSourceFileName.substr( 0, fSourceFileName.size() - 3 );
            fBZipped = true;
        }
        if( fSourceFileName.size() > 4 &&  fSourceFileName.size() - fSourceFileName.rfind( ".bz2" ) == 4 )
        {
            cout << "VGrIsuReader::VGrIsuReader info: found bzipped2 file, unzipping..." << endl;
            string i_com = "bunzip2 -f -c  ";
            i_com += fSourceFileName;
            fSourceFileName = fSourceFileName.substr( 0, fSourceFileName.size() - 4 );
            i_com += " > " + fSourceFileName;
            cout << "\t (" << i_com << ")" << endl;
            gSystem->Exec( i_com.c_str() );
            fBZipped2 = true;
        }
        
        // open the source file
        is.open( fSourceFileName.c_str() );
        if( !is.is_open() )
        {
            cout << "VGrIsuReader::VGrIsuReader() - error opening file (" << fSourceFileName << ")" << endl;
            exit( -1 );
        }
        if( is )
        {
            cout << "VGrIsuReader::VGrIsuReader() opened: " << fSourceFileName << endl;
        }
    }
    else
    {
        // open the pedestal file (shouldn't be zipped)
        if( fExternalPedFile.size() > 0 )
        {
            is_ped.open( fExternalPedFile.c_str() );
        }
        // or use file with simulations data (should be opened second -> unzipped)
        else
        {
            is_ped.open( fSourceFileName.c_str() );
        }
        cout << fExternalPedFile << "\t external file " << endl;
        if( !is_ped.is_open() )
        {
            cout << "VGrIsuReader::VGrIsuReader() - error opening pedestal file (" << fExternalPedFile << ")" << endl;
            exit( -1 );
        }
        if( is_ped )
        {
            cout << "reading pedestals from pedestal lines in file ";
        }
        if( fExternalPedFile.size() > 0 )
        {
            cout << fExternalPedFile << endl;
        }
        else
        {
            cout << fSourceFileName << endl;
        }
    }
}


void VGrIsuReader::resetEvent()
{
    fSelectedHitChan.assign( fNTel, 0 );
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fFullTrigVec[i].assign( fMaxChannels[i], false );
        fNumberofFullTrigger[i] = 0;
    }
}


/*!
     For GrIsu files without pedestal lines (tag "P"), default values for the pedestal and its
     variation are taken (fdefaultPed, fdefaultPedvars)
*/
void VGrIsuReader::readPeds( unsigned int iDummy )
{
    if( fDebug )
    {
        cout << "VGrIsuReader::readPeds() " << iDummy << endl;
    }
    
    char hname[200];
    
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        // reset everything
        valarray<double> v_dtemp( 0., fMaxChannels[i] );
        valarray<double> vp_dtemp( fdefaultPed, fMaxChannels[i] );
        vector< valarray< double > > vv_dtemp;
        fPeds[i] = vp_dtemp;
        fPeds[i] = fdefaultPed;
        fPedvars[i] = v_dtemp;
        v_dtemp  = fdefaultPedvars;
        for( uint16_t s = 0; s < fNumSamples[i] + 1; s++ )
        {
            vv_dtemp.push_back( v_dtemp );
        }
        fVPedvars[i] = vv_dtemp;
        fPedvars[i] = fdefaultPedvars;
        fPedRMS[i] = v_dtemp;
        fPedRMS[i] = fdefaultPedvars;
        vector< TH1F* > i_hpeds;
        for( unsigned int j = 0; j < fMaxChannels[i]; j++ )
        {
            sprintf( hname, "hped_%d_%d_%d", i, fSumWindow[i], j );
            TH1F* iht = ( TH1F* )gDirectory->Get( hname );
            if( iht != 0 )
            {
                iht->Delete();
            }
            i_hpeds.push_back( new TH1F( hname, "", 50 * fSumWindow[i], 0, 50 * fSumWindow[i] ) );
        }
        fhPeds.push_back( i_hpeds );
        
        if( fSumWindow[i] <= 0. )
        {
            cout << "VGrIsuReader::readPeds() error: sumwindow <= 0. for telescope " << i << endl;
            cout << "\t Is the number of telescopes set correctly via -ntelescopes?" << endl;
            exit( 0 );
        }
    }
    
    // case A:
    //    using background traces (read pedestal variations from trace library file)
    //      one tree entry per tube -> entry number = tubenumber
    //      numbering in fPeds is MC numbering
    if( fTraceFileName.size() > 0 )
    {
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            readPedsfromLibrary( i );
        }
    }
    // case B:
    //    pedestals are taken from GrIsu "P" lines
    else
    {
        readPedsfromPlines();
    }
}


/*!
    the tree names with the background traces have the following naming convention:
    Telescope 1: pedTree_0
    Telescope 2: pedTree_1
    ....

*/
void VGrIsuReader::readPedsfromLibrary( unsigned int iTelID )
{
    if( fDebug )
    {
        cout << "VGrIsuReader::readPedsfromLibrary()" << endl;
    }
    // tree names with telescope numbering _0 -> telescope 0, _1 -> telescope 1
    if( iTelID >= fNTel )
    {
        cout << "VGrIsuReader::readPedsfromLibrary: telescope ID out of range: " << iTelID << "\t" << fNTel << endl;
        return;
    }
    unsigned int i = iTelID;
    
    fTraceFile->cd();
    char iTreeName[200];
    sprintf( iTreeName, "pedTree_%d", i );
    TTree* iTree = ( TTree* )gDirectory->Get( iTreeName );
    if( iTree == 0 )
    {
        cout << "VGrIsuReader::readPeds() error, pedestal tree not found: " << iTreeName << endl;
        exit( -1 );
    }
    double iped, ipedvar, isumwindow;
    unsigned int tubeNumber;
    iTree->SetBranchAddress( "tubeNumber", &tubeNumber );
    iTree->SetBranchAddress( "ped", &iped );
    iTree->SetBranchAddress( "pedvar", &ipedvar );
    iTree->SetBranchAddress( "sumwindow", &isumwindow );
    // first set all pedestals to zero (tree contains only non-zero pedestals)
    for( unsigned int j = 0; j < fPeds[i].size(); j++ )
    {
        fPeds[i][j] = 0.;
    }
    for( unsigned int j = 0; j < fPedvars[i].size(); j++ )
    {
        fPedvars[i][j] = 0.;
    }
    // loop over all tubes ( tubeNumber is real data camera file tube numbering)
    for( int j = 0; j < iTree->GetEntries(); j++ )
    {
        // test number of channels in tree
        if( j == ( int )fMaxChannels[i] )
        {
            break;
        }
        iTree->GetEntry( j );
        // in trace library, tube  numbering is according to real data camera file
        //       transform here to MC numbering (with fPixelConvertVecM[][])
        if( tubeNumber < fPeds[i].size() )
        {
            fPeds[i][fPixelConvertVecM[i][tubeNumber]] = iped;
        }
        if( tubeNumber < fPedvars[i].size() )
        {
            fPedvars[i][fPixelConvertVecM[i][tubeNumber]] = ipedvar;
        }
        if( TMath::Abs( fSumWindow[i] - isumwindow ) > 1.e-4 )
        {
            cout << "VGrIsuReader::readPeds() error, incompatible sum window sizes: " << isumwindow << "\t" << fSumWindow[i] << "\t" << i << endl;
            exit( -1 );
        }
    }
    iTree->ResetBranchAddresses();
    return;
}


void VGrIsuReader::readPedsfromPlines()
{
    if( fDebug )
    {
        cout << "VGrIsuReader::readPedsfromPlines" << endl;
    }
    if( TMath::Abs( fFADCScale - 1. ) > 0.01 )
    {
        cout << "VGrIsuReader::readPedsfromPlines: scale fadc traces by " << fFADCScale << endl;
    }
    
    if( is_ped.is_open() )
    {
        is_ped.close();
    }
    // open the pedestal file (shouldn't be zipped)
    if( fExternalPedFile.size() > 0 )
    {
        is_ped.open( fExternalPedFile.c_str() );
    }
    // or use file with simulations data (should be opened second -> unzipped)
    else
    {
        is_ped.open( fSourceFileName.c_str() );
    }
    cout << "\t read noise from: " << fExternalPedFile << endl;
    if( !is_ped )
    {
        cout << "VGrIsuReader::VGrIsuReader() - error opening pedestal file (" << fExternalPedFile << ")" << endl;
        exit( -1 );
    }
    
    string is_line;
    string is_Temp;
    
    unsigned int i_telID = 0;
    unsigned int i_channel = 0;
    vector<double> i_val;
    unsigned int i_val_Size = 0;
    vector<double> i_valMean;
    unsigned int i_valMean_Size = 0;
    vector<double> i_valMean2;
    double meanRMS2 = 0.;
    double mean = 0.;
    double mean2 = 0.;
    double rms = 0.;
    int i_sampleTemp;
    
    // fixed number of pedestal signal, maximum is 1500 samples
    vector<uint8_t> i_pedSample( 1500, 0 );
    for( unsigned int i_telID = 0; i_telID < fNTel; i_telID++ )
    {
        vector< vector<uint8_t> > i_pedChannelSample( fMaxChannels[i_telID], i_pedSample );
        fGrIsuPeds.push_back( i_pedChannelSample );
    }
    i_pedSample.reserve( 1500 );
    
    bool pedFound = false;
    while( getline( is_ped, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            // assume, that pedestals are always before(!) simulation and trace data
            if( is_line.substr( 0, 1 ) == "R" || is_line.substr( 0, 1 ) == "S" )
            {
                if( !pedFound )
                {
                    //////////////////////////////////////////////////////////////////////////////////////////////
                    // what to do, if there are no "P" lines in the GrIsu file?
                    //    setting everything to standard values
                    for( unsigned int i_telID = 0; i_telID < fNTel; i_telID++ )
                    {
                        for( unsigned int j = 0; j < fMaxChannels[i_telID]; j++ )
                        {
                            fPeds[i_telID][j] = fdefaultPed;
                            fPedvars[i_telID][j] = fdefaultPedvars;
                        }
                    }
                }
                return;
                //////////////////////////////////////////////////////////////////////////////////////////////
            }
            if( is_line.substr( 0, 1 ) == "P" )
            {
                pedFound = true;
                meanRMS2 = 0.;
                i_val.clear();
                i_pedSample.clear();
                istringstream is_stream( is_line );
                is_stream >> is_Temp;             // "P"
                is_stream >> is_Temp;             // placeholder
                is_stream >> is_Temp;             // telescope number
                //	    if( i_telID != (unsigned int)(atoi( is_Temp.c_str() ) - fTelNumberOffset) ) continue;
                i_telID = atoi( is_Temp.c_str() ) - fTelNumberOffset;
                if( i_telID >= fNTel )
                {
                    continue;
                }
                is_stream >> is_Temp;             // channel number
                i_channel = ( uint32_t )( atoi( is_Temp.c_str() ) - fTelNumberOffset );
                if( i_channel > fMaxChannels[i_telID] )
                {
                    continue;
                }
                fhPeds[i_telID][i_channel]->Reset();
                is_stream >> is_Temp;
                while( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> is_Temp;
                    i_val.push_back( atof( is_Temp.c_str() ) );
                    i_val.back() = ( int )( fdefaultPed + ( i_val.back() - fdefaultPed ) / fFADCScale );
                    meanRMS2 += ( double )i_val.back() * ( double )i_val.back();
                    i_sampleTemp = atoi( is_Temp.c_str() );
                    i_sampleTemp  = ( int )( fdefaultPed + ( i_sampleTemp - fdefaultPed ) / fFADCScale );
                    i_pedSample.push_back( ( uint8_t )i_sampleTemp );
                }
                fGrIsuPeds[i_telID][i_channel] = i_pedSample;
                fGrIsuPedsN[i_telID][i_channel] = i_pedSample.size();
                
                // calculate pedestal variations for different summation windows
                if( i_channel < fMaxChannels[i_telID] )
                {
                    i_val_Size = i_val.size();
                    if( i_val_Size > 0 )
                    {
                        meanRMS2 /= ( double )i_val_Size;
                    }
                    for( unsigned int w = 1; w < fVPedvars[i_telID].size(); w++ )
                    {
                        // calculate window mean
                        i_valMean.clear();
                        i_valMean2.clear();
                        // sum up in fSumWindow long bunches
                        //		  for( unsigned int s = 0; s < i_val_Size; s += fSumWindow[i_telID] )
                        for( unsigned int s = 0; s < i_val_Size; s += w )
                        {
                            mean = mean2 = rms = 0.;
                            if( s + w > i_val_Size )
                            {
                                break;
                            }
                            for( unsigned int j = s; j < s + w; j++ )
                            {
                                mean += i_val[j];
                            }
                            if( TMath::Abs( mean ) > 1.e-4 )
                            {
                                i_valMean.push_back( mean );
                                if( ( int )w == fSumWindow[i_telID] )
                                {
                                    fhPeds[i_telID][i_channel]->Fill( mean );
                                }
                                i_valMean2.push_back( mean * mean );
                            }
                        }
                        i_valMean_Size = i_valMean.size();
                        // calculate mean and mean2 (E[x**2])
                        mean = 0.;
                        mean2 = 0.;
                        rms = 0.;
                        for( unsigned int s = 0; s < i_valMean_Size; s++ )
                        {
                            mean += i_valMean[s];
                            mean2 += i_valMean2[s];
                        }
                        if( i_valMean_Size > 0 )
                        {
                            mean /= ( double )i_valMean_Size;
                            mean2 /= ( double )i_valMean_Size;
                        }
                        if( w != 0 && ( int )w == fSumWindow[i_telID] )
                        {
                            fPeds[i_telID][i_channel] = mean / ( double )w;
                        }
                        // "normal" rms
                        if( w != 0 )
                        {
                            meanRMS2 = meanRMS2 - mean * mean / w / w;
                        }
                        if( meanRMS2 >= 0. && ( int )w == fSumWindow[i_telID] )
                        {
                            fPedRMS[i_telID][i_channel] = sqrt( meanRMS2 );
                        }
                        rms = mean2 - mean * mean;
                        if( rms > 0 )
                        {
                            rms = sqrt( rms );
                        }
                        fVPedvars[i_telID][w][i_channel] = rms;
                        if( ( int )w == fSumWindow[i_telID] )
                        {
                            fPedvars[i_telID][i_channel] = rms;
                        }
                    }
                }
            }
        }
    }
    if( is_ped.is_open() )
    {
        is_ped.close();
    }
    if( fDebug )
    {
        cout << " ..." << endl;
    }
}


/*!
     this is finetuned to the grisudet output
*/
bool VGrIsuReader::getNextEvent()
{
    if( fDebug )
    {
        cout << "VGrIsuReader::getNextEvent " << fPedestalMode << endl;
    }
    
    bool iSuccess = false;
    if( fPedestalMode )
    {
        iSuccess = getNextPedestalEvent();
    }
    else
    {
        iSuccess = getNextShowerEvent();
    }
    
    if( !iSuccess )
    {
        setEventStatus( 999 );
    }
    
    if( fDebug )
    {
        cout << "VGrIsuReader::getNextEvent, event status: " << getEventStatus() << endl;
    }
    
    return iSuccess;
}


/*!
   make pedestal events from grisudet 'P' lines,
   extract traces by random start values
*/
bool VGrIsuReader::getNextPedestalEvent()
{
    if( fDebug )
    {
        cout << "VGrIsuReader::getNextPedestalEvent()" << endl;
    }
    // increment event number
    fEventNumber++;
    // is maximal number of pedestal events reached?
    if( fEventNumber > fNPedestalEvents )
    {
        return false;
    }
    // set eventtype PEDESTAL
    fEventType = 2;
    unsigned int iPStart = 0;
    // loop over all telescopes
    for( unsigned int t = 0; t < fNTel; t++ )
    {
        // take background generated by grisu ('P' lines)
        for( unsigned int i = 0; i < fMaxChannels[t]; i++ )
        {
            // get random start value
            iPStart = fRandomGen->Integer( fGrIsuPeds[t][i].size() - fNumSamples[t] );
            for( unsigned int j = 0; j < fNumSamples[t]; j++ )
            {
                fSamplesVec[t][i][j] = ( uint8_t )fGrIsuPeds[t][i][iPStart + j];
            }
        }
    }
    // set mc values
    fMC_primary = 1;
    fMC_energy = 0.;
    fMC_X = 0.;
    fMC_Y = 0.;
    fMC_Xcos = 0.;
    fMC_Ycos = 0.;
    fMC_Ze = 0.;
    fMC_Az = 0.;
    // set local trigger
    fLocalTrigger.assign( fNTel, true );
    return true;
}


bool VGrIsuReader::getNextShowerEvent()
{
    if( fDebug )
    {
        cout << "bool VGrIsuReader::getNextShowerEvent()" << endl;
    }
    // set eventtype IGNORE
    fEventType = 1;
    // is there anything left in the file?
    if( is.eof() )
    {
        return false;
    }
    
    unsigned int i_telID = 0;
    uint32_t i_channel = 0;
    uint16_t i_nsample = 0;
    
    cout << fNumSamples[0] << endl;
    
    string is_line = "";
    string is_Temp = "";
    bool i_eventFound = false;
    
    int iRecord = 0;
    int iEvent = 0;
    
    // read triggered pixels from GrIsu file
    sp = is.tellg();
    char iSRecord;
    while( !is.eof() )
    {
        is >> iSRecord;
        is_line = iSRecord;
        // skip pedestal lines
        if( is_line == "P" )
        {
            continue;
        }
        //////////////////////////////////////////////////////
        // MC DATA (SHOWER LINE)
        //////////////////////////////////////////////////////
        else if( is_line == "S" )
        {
            if( fDebug )
            {
                cout << "bool VGrIsuReader::getNextShowerEvent(), S line" << endl;
            }
            if( i_eventFound )
            {
                is.putback( iSRecord );
                break;
            }
            string a_temp;
            is >> a_temp;
            // get shower data
            //	 is >> fEventNumber;
            fEventNumber = ( uint32_t )( atoi( a_temp.c_str() ) );
            fEventNumber += 1;                    // does eventnumbers start at 0, in eventdisplay at 1
            // check event numbers:
            // -- first event should be event 1
            // -- no missing events
            // (everything else points to a problem during the simulation runs)
            /*	 if( (fEventNumberPreviousEvent == 0 && fEventNumber != 1) || (fEventNumber-fEventNumberPreviousEvent) != 1 )
                 {
                     cout << "VGrIsuReader:: problem with event numbers" << endl;
                     cout << "\t current event number: " << fEventNumber << endl;
                     cout << "\t event number of previous event: " << fEventNumberPreviousEvent << endl;
                     return false;
                     } */
            fEventNumberPreviousEvent = fEventNumber;
            is >> fMC_primary;                    // but in grisudet: primary always zero
            is >> fMC_energy;
            is >> fMC_X;
            is >> fMC_Y;
            // direction of telescope plane, shower direction maybe offset
            is >> fMC_Xcos;
            is >> fMC_Ycos;
            is >> is_Temp;
            if( is_Temp == "nan" )
            {
                fMC_Xoffset = 0.;
            }
            else
            {
                fMC_Xoffset = atof( is_Temp.c_str() );
            }
            is >> is_Temp;
            if( is_Temp == "nan" )
            {
                fMC_Yoffset = 0.;
            }
            else
            {
                fMC_Yoffset = atof( is_Temp.c_str() );
            };
            // signs changed with GrIsu 4.12 (from camera to sky)
            if( fGrisuVersion > 412 )
            {
                fMC_Xoffset *= -1.;
                fMC_Yoffset *= -1.;
            }
            // calculate zenith and azimuth angle from direction cosinii
            // (GM 20090728)	 fMC_Az = atan2( fMC_Xcos, fMC_Ycos ) + 180./degrad;
            fMC_Az = atan2( fMC_Xcos, fMC_Ycos );
            
            fMC_Ze = 1. - ( fMC_Xcos * fMC_Xcos + fMC_Ycos * fMC_Ycos );
            if( fMC_Ze < 0. )
            {
                fMC_Ze = 0.;
            }
            fMC_Ze = acos( sqrt( fMC_Ze ) );
            // all telescopes are pointing in the same direction
            for( unsigned int i = 0; i < fNTel; i++ )
            {
                fTelElevation[i] = 90. - fMC_Ze * degrad;
                fTelAzimuth[i] = fMC_Az * degrad;
            }
            // add wobble offset
            double az = 0.;
            double ze = 0.;
            VSkyCoordinatesUtilities::getRotatedShowerDirection( fMC_Ze * degrad, fMC_Az * degrad, fMC_Yoffset, -1.*fMC_Xoffset, ze, az );
            fMC_Ze = ze / degrad;
            fMC_Az = VSkyCoordinatesUtilities::adjustAzimuthToRange( az ) / degrad;
            resetEvent();
            i_eventFound = true;
            if( fDebug )
            {
                cout << "\t Event (S): ";
                cout << fEventNumber << "\t" << fMC_primary << "\t" << fMC_energy << "\t" << fMC_X << "\t" << fMC_Y << "\t";
                cout << fMC_Xcos << "\t" << fMC_Ycos << "\t" << fMC_Xoffset << "\t" << fMC_Yoffset << endl;
            }
        }
        //////////////////////////////////////////////////////
        // TRIGGER RECORD
        //////////////////////////////////////////////////////
        else if( is_line == "W" )
        {
            if( fDebug )
            {
                cout << "bool VGrIsuReader::getNextShowerEvent(), W line (" << fNTel << " telescopes)" << endl;
            }
            float i_T = 0.;
            is >> is_Temp;
            is >> i_T;
            fArrayTrigger = int( i_T + 0.5 );
            // create trigger word
            fLocalTrigger.clear();
            fLocalTrigger.assign( fNTel, false );
            ///////////////////////////
            // no multiple raw data reader defined, information of all telescopes is in the data
            if( !fMultiGrIsuReader )
            {
                // read the local triggers
                if( fNTel > fLocalTrigger.size() )
                {
                    cout << "VGrIsuReader::getNextShowerEvent() warning: too many telescopes for trigger word" << endl;
                }
                for( unsigned int i = 0; i < fNTel; i++ )
                {
                    is >> i_T;
                    if( int( i_T + 0.5 ) )
                    {
                        fLocalTrigger[i] = true;
                    }
                    else
                    {
                        fLocalTrigger[i] = false;
                    }
                }
                // read the local trigger times
                for( unsigned int i = 0; i < fNTel; i++ )
                {
                    is >> fLocalTriggerTime[i];
                }
                // read the local delayed trigger times
                for( unsigned int i = 0; i < fNTel; i++ )
                {
                    is >> fLocalDelayedTriggerTime[i];
                }
            }
            ///////////////////////////
            // multiple raw data reader, one telescope only
            else
            {
                is >> i_T;
                if( int( i_T + 0.5 ) )
                {
                    fLocalTrigger[0] = true;
                }
                else
                {
                    fLocalTrigger[0] = false;
                }
                is >> fLocalTriggerTime[0];
                is >> fLocalDelayedTriggerTime[0];
            }
            ///////////////////////////
            // read the grisu cores
            is >> fXimpactrot;
            is >> fYimpactrot;
            if( fDebug )
            {
                cout << "\t" << fArrayTrigger << "\t";
                for( unsigned int t = 0; t < fLocalTrigger.size(); t++ )
                {
                    cout << fLocalTrigger[t] << "\t";
                }
                for( unsigned int t = 0; t < fLocalTriggerTime.size(); t++ )
                {
                    cout << fLocalTriggerTime[t] << "\t";
                }
                for( unsigned int t = 0; t < fLocalDelayedTriggerTime.size(); t++ )
                {
                    cout << fLocalDelayedTriggerTime[t] << "\t";
                }
                cout << fXimpactrot << "\t" << fYimpactrot << endl;
            }
        }
        //////////////////////////////////////////////////////
        // FADC RECORD
        //////////////////////////////////////////////////////
        else if( is_line == "R" )
        {
            //	 if( fDebug ) cout << "bool VGrIsuReader::getNextShowerEvent(), R line" << endl;
            is >> iEvent;
            if( iEvent + 1 != ( int )fEventNumber )
            {
                cout << "VGrIsuReader::getNextEvent(): error in reading event number ";
                cout << fEventNumber << "\t" << iEvent << "\t" << fSourceFileName << endl;
            }
            // telescope ID
            is >> i_telID;
            i_telID -= fTelNumberOffset;
            if( i_telID >= fNTel )
            {
                //	    cout << "VGrIsuReader::getNextEvent(): error in telescope number: " << endl;
                //	    cout << i_telID << "\t" << fNTel << endl;
                continue;
            }
            // set telescope number in detector geometry
            fDetectorGeo->setTelID( i_telID );
            // channel number
            is >> iRecord;
            iRecord -= fTelNumberOffset;
            i_channel = ( uint32_t )iRecord;
            if( i_channel > fMaxChannels[i_telID] )
            {
                cout << "VGrIsuReader::getNextEvent(): error in channel numbers " << i_channel << " " << iRecord << " " << fTelNumberOffset << " " << " (event " << fEventNumber << ")" << endl;
                // this means something bad happened to the grisu file
                cout << "\t check grisu file" << endl;
                return false;
            }
            // in old grisu version, only trigger channels were in grisudet output file -> trigger vector
            fFullTrigVec[i_telID][( unsigned int )i_channel] = 1;
            // reset hilo switch
            fHiLo[i_telID][i_channel] = false;
            // get hilo settings (new grisuversion only)
            //         if( fDetectorGeo->getGrIsuVersion() >= 413 )
            if( fGrisuVersion >= 413 )
            {
                is >> iRecord;
                fHiLo[i_telID][i_channel] = iRecord;
            }
            // get FADC trace
            i_nsample = fNumSamples[i_telID];
            for( int i = 0; i < ( int )( i_nsample + fSampleOffset ); i++ )
            {
                is >> iRecord;
                iRecord = ( int )( ( iRecord - fdefaultPed ) / fFADCScale + fdefaultPed );
                // overflow
                //	       if( fDetectorGeo->getGrIsuVersion() < 413 && iRecord > (int)fDetectorGeo->getLowGainThreshold()[i_telID] )
                if( fGrisuVersion < 413 && iRecord > ( int )fDetectorGeo->getLowGainThreshold()[i_telID] )
                {
                    fHiLo[i_telID][i_channel] = true;
                }
                
                if( i >= fSampleOffset )
                {
                    ftempSampleUI[i - fSampleOffset] = iRecord;
                    ftempSampleFL[i - fSampleOffset] = ( ( iRecord - fdefaultPed ) / fFADCScale + fdefaultPed );
                }
            }
            for( unsigned int i = 0; i < fSamplesVec[i_telID][i_channel].size(); i++ )
            {
                // apply low gain for large pulses (>255)
                //	    if( fDetectorGeo->getGrIsuVersion() < 413 && fHiLo[i_telID][i_channel] )
                if( fGrisuVersion < 413 && fHiLo[i_telID][i_channel] )
                {
                    double iT = ftempSampleFL[i];
                    iT -= fdefaultPed;
                    //       iT /= fDetectorGeo->getLowGainMultiplier()[i_telID];
                    iT /= fDetectorGeo->getLowGainMultiplier_Trace()[i_telID];
                    iT += fdefaultPed;
                    ftempSampleFL[i] = ( int )( iT + fPeds[i_telID][i_channel] - fdefaultPed );
                }
                // fill sample vector
                if( ftempSampleFL[i] > fDetectorGeo->getLowGainThreshold()[i_telID] )
                {
                    fSamplesVec[i_telID][i_channel][i] = fDetectorGeo->getLowGainThreshold()[i_telID];
                }
                else
                {
                    fSamplesVec[i_telID][i_channel][i] = ( uint8_t )ftempSampleFL[i];
                }
            }
        }
    }
    
    // only make background trace if there was an arraytrigger
    if( fArrayTrigger != 0 || fMultiGrIsuReader )
    {
        // now produce all traces with no trigger and noise  (random electronic noise + NSB) (only hit channels)
        // two possibilities:
        //        - take background from trace library
        if( fTraceFileName.size() > 0 )
        {
            fillBackgroundfromTraceLibrary();
        }
        //        - take background generated by grisu
        else
        {
            fillBackgroundfromPlines();
        }
    }
    
    return true;
}


std::vector< uint8_t >  VGrIsuReader::getSamplesVec()
{
    if( fSelectedHitChan[fTelescopeID] < fMaxChannels[fTelescopeID] )
    {
        return fSamplesVec[fTelescopeID][fSelectedHitChan[fTelescopeID]];
    }
    ftempSample.resize( fNumSamples[fTelescopeID], 0 );
    return ftempSample;
}


std::pair<bool, uint32_t> VGrIsuReader::getChannelHitIndex( uint32_t hit )
{
    if( hit < fMaxChannels[fTelescopeID] )
    {
        return std::make_pair( true, hit );
    }
    return std::make_pair( false, ( uint32_t ) 0 );
    
}


bool VGrIsuReader::getHiLo( uint32_t i )
{
    if( i < fMaxChannels[fTelescopeID] )
    {
        return ( bool )fHiLo[fTelescopeID][i];
    }
    return 0;
}


uint32_t VGrIsuReader::getHitID( uint32_t i )
{
    if( i < fMaxChannels[fTelescopeID] )
    {
        return i;
    }
    return 0;
}


void VGrIsuReader::selectHitChan( uint32_t hit )
{
    if( hit < fMaxChannels[fTelescopeID] )
    {
        fSelectedHitChan[fTelescopeID] = hit;
    }
    else
    {
        fSelectedHitChan[fTelescopeID] = 0;
    }
}


bool VGrIsuReader::setTelescopeID( unsigned int i_tel )
{
    if( i_tel < fNTel )
    {
        fTelescopeID = i_tel;
    }
    else
    {
        return false;
    }
    
    return true;
}


/*!
     the reader needs to following information from the camera files:
        - trigger information (which tubes exist)
    - analysis pixel (tubes can be switched on/off)
    - convertion from real data tube numbering to MC numbering

    \param iTel telescope number
*/
bool VGrIsuReader::readCamera( int iTel )
{
    if( fDebug )
    {
        cout << "bool VGrIsuReader::readCamera( int iTel ) " << iTel << endl;
    }
    
    // set telescope number in detector geometry
    fDetectorGeo->setTelID( iTel );
    
    if( fDetectorGeo->getTrigger().size() < fFullHitVec[iTel].size() )
    {
        cout << "VGrIsuReader::readCameraInfo: no trigger info" << endl;
        return false;
    }
    else
    {
        for( unsigned int i = 0; i < fFullHitVec[iTel].size(); i++ )
        {
            fFullHitVec[iTel][i] = ( bool )fDetectorGeo->getTrigger()[i];
        }
    }
    // fFullAnaVec()[i]: -1: dead channel, 0: channel does not exist, 1 channel o.k.
    if( fDetectorGeo->getAnaPixel().size() < fFullAnaVec[iTel].size() )
    {
        cout << "VGrIsuReader::readCameraInfo: no trigger info" << endl;
        return false;
    }
    else
    {
        for( unsigned int i = 0; i < fFullAnaVec[iTel].size(); i++ )
        {
            fFullAnaVec[iTel][i] = fDetectorGeo->getAnaPixel()[i];
        }
    }
    return true;
}


bool VGrIsuReader::readPixelMix( int iTel )
{
    if( fDebug )
    {
        cout << "bool VGrIsuReader::readPixelMix( iTel ) " << iTel << endl;
    }
    
    // set telescope number in detector geometry
    fDetectorGeo->setTelID( iTel );
    
    // reading vector for convertion of camera pixel numbers from real data numbering to MC numbering
    if( fTraceFileName.size() > 0 )
    {
        if( fDetectorGeo->getRealPixel().size() < fPixelConvertVecM[iTel].size() )
        {
            for( unsigned int i = 0; i < fPixelConvertVecM[iTel].size(); i++ )
            {
                fPixelConvertVecM[iTel][i] = i;
                fPixelConvertVecR[iTel][i] = i;
            }
            //        cout << "VGrIsuReader::readPixelMix: no tube numbers mixing info" << endl;
            //        return false;
        }
        else
        {
            for( unsigned int i = 0; i < fPixelConvertVecM[iTel].size(); i++ )
            {
                fPixelConvertVecM[iTel][i] = fDetectorGeo->getMCPixel()[i];
                fPixelConvertVecR[iTel][i] = fDetectorGeo->getRealPixel()[i];
            }
        }
    }
    return true;
}


bool VGrIsuReader::openTraceLibrary()
{
    if( fDebug )
    {
        cout << "void VGrIsuReader::openTraceLibrary()" << endl;
    }
    if( fTraceFileName.size() > 0 )
    {
        fTraceFile = new TFile( fTraceFileName.c_str() );
        if( fTraceFile->IsZombie() )
        {
            cout << "VGrIsuReader::VGrIsuReader(): error, trace library not found " << fTraceFileName << endl;
            exit( -1 );
        }
        cout << "VGrIsuReader::VGrIsuReader(): reading tracefile: " << fTraceFileName << endl;
        char itree[200];
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            sprintf( itree, "trace_%d", i );
            fTraceTree.push_back( ( TTree* )gDirectory->Get( itree ) );
            if( fTraceTree.back() == 0 )
            {
                cout << "VGrIsuReader::VGrIsuReader(): warning, trace tree not found: " << itree << endl;
            }
            fTraceTree.back()->SetBranchAddress( "fTrace", fTrace );
        }
    }
    return true;
}


void VGrIsuReader::initVectors()
{
    if( fDebug )
    {
        cout << "void VGrIsuReader::initVectors() " << fNTel << "\t" << endl;
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            cout <<  "\t for telescope " << i + 1 << ": ";
            cout << "nchannel: " << fDetectorGeo->getNChannels( i );
            cout << "\t nsample: " << fDetectorGeo->getNSamples( i ) << endl;
        }
    }
    vector<bool> v_temp;
    vector< unsigned int> v_temp_int;
    vector< int> v_temp_int2;
    vector< vector< uint8_t > > v_tempSample;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fMaxChannels.push_back( fDetectorGeo->getNChannels( i ) );
        fNumSamples.push_back( fDetectorGeo->getNSamples( i ) );
        
        v_temp.resize( fMaxChannels.back(), true );
        fFullHitVec.push_back( v_temp );
        v_temp_int2.resize( fMaxChannels.back(), 1 );
        fFullAnaVec.push_back( v_temp_int2 );
        v_temp_int.resize( fMaxChannels.back(), 0 );
        fPixelConvertVecM.push_back( v_temp_int );
        fPixelConvertVecR.push_back( v_temp_int );
        fNumChannelsHit.push_back( fFullHitVec.back().size() );
        v_temp.assign( fMaxChannels.back(), false );
        fFullTrigVec.push_back( v_temp );
        fNumberofFullTrigger.push_back( 0 );
        fHiLo.push_back( v_temp );
        
        fGrIsuPedsN.push_back( v_temp_int );
        
        ftempSampleUI.resize( fNumSamples.back(), 0 );
        ftempSampleFL.resize( fNumSamples.back(), 0 );
        ftempSample.resize( fNumSamples.back(), 0 );
        v_tempSample.resize( fMaxChannels.back(), ftempSample );
        fSamplesVec.push_back( v_tempSample );
        
        valarray<double> v_dtemp( 0., fMaxChannels.back() );
        fPeds.push_back( v_dtemp );
        fPedvars.push_back( v_dtemp );
        fPedRMS.push_back( v_dtemp );
        vector< valarray<double> > vv_dtemp;
        for( uint16_t s = 0; s < fNumSamples[i] + 1; s++ )
        {
            vv_dtemp.push_back( v_dtemp );
        }
        fVPedvars.push_back( vv_dtemp );
        
        // trigger vectors
        fLocalTriggerTime.push_back( 0. );
        fLocalDelayedTriggerTime.push_back( 0. );
        
        // pointing vectors
        fTelElevation.push_back( 0. );
        fTelAzimuth.push_back( 0. );
    }
    fSelectedHitChan.resize( fNTel, 0 );
}


void VGrIsuReader::fillBackgroundfromTraceLibrary()
{
    if( fDebug )
    {
        cout << "void VGrIsuReader::fillBackgroundfromTraceLibrary()" << endl;
    }
    int i_traceValue = 0;
    double i_newPed = 0.;
    int i_newPedN = 0;
    cout << "void VGrIsuReader::fillBackgroundfromTraceLibrary()" << endl;
    
    for( unsigned int h = 0; h < fNTel; h++ )
    {
        // get randomly a event from the trace library
        fTraceTree[h]->GetEntry( fRandomGen->Integer( ( int )fTraceTree[h]->GetEntries() ) );
        // loop over all tubes (MC camera file tube numbering)
        for( unsigned int i = 0; i < fMaxChannels[h]; i++ )
        {
            i_newPed = 0.;
            i_newPedN = 0;
            for( unsigned int j = 0; j < fNumSamples[h]; j++ )
            {
                i_traceValue = 0;
                if( fFullHitVec[h][i] && fFullAnaVec[h][i] == 1 )
                {
                    // hit pixel -> only add noise
                    if( fFullTrigVec[h][i] )
                    {
                        // use pedestals from data file -> first subtract the pedestal from Grisu (fdefaultped),
                        // then add the real one
                        i_traceValue = ( int )( fSamplesVec[h][i][j] - fdefaultPed );
                        i_traceValue += fTrace[fPixelConvertVecR[h][i]][j];
                        if( i_traceValue < 0 )
                        {
                            cout << "VGrIsuReader::fillBackgroundfromTraceLibrary() warning trace value < 0: (a)";
                            cout << h << "\t" << i << "\t" << j << "\t";
                            cout << i_traceValue << "\t" << fFullTrigVec[h][i] << "\t" << fSamplesVec[h][i][j] << "\t";
                            cout << fEventNumber << "\t" << fTraceTree[h]->GetReadEvent() << "\t";
                            cout << fdefaultPed << "\t" << fPeds[h][i] << endl;
                        }
                        fSamplesVec[h][i][j] = ( uint8_t )i_traceValue;
                        if( i_traceValue != 0 )
                        {
                            i_newPed += ( double )i_traceValue;
                            i_newPedN++;
                        }
                    }
                    // silent pixel -> set pedestal and add noise
                    else
                    {
                        i_traceValue = ( int )( fTrace[fPixelConvertVecR[h][i]][j] );
                        if( i_traceValue < 0 )
                        {
                            cout << "VGrIsuReader::fillBackgroundfromTraceLibrary() warning trace value < 0: (b) ";
                            cout << h << "\t" << i << "\t" << j << "\t";
                            cout << i_traceValue << "\t" << fFullTrigVec[h][i] << "\t" << fSamplesVec[h][i][j] << "\t";
                            cout << fEventNumber << "\t" << fTraceTree[h]->GetReadEvent() << "\t";
                            cout << fdefaultPed << "\t" << fPeds[h][i] << endl;
                        }
                        fSamplesVec[h][i][j] = ( uint8_t )i_traceValue;
                        if( i_traceValue != 0 )
                        {
                            i_newPed += ( double )i_traceValue;
                            i_newPedN++;
                        }
                    }
                }
                else
                {
                    fSamplesVec[h][i][j] = 0;
                }
            }
        }
    }
}


void VGrIsuReader::fillBackgroundfromPlines()
{
    bool bAdd = false;
    if( fExternalPedFile.size() > 0 || fMultiGrIsuReader )
    {
        bAdd = true;
    }
    
    for( unsigned int h = 0; h < fNTel; h++ )
    {
        // take background generated by grisu ('P' lines)
        for( unsigned int i = 0; i < fMaxChannels[h]; i++ )
        {
            // select a random start point in long background trace
            unsigned int nstart = fRandomGen->Integer( fGrIsuPeds[h][i].size() - fNumSamples[h] - 1 );
            
            if( fFullHitVec[h][i] )
            {
                // generate a trace from background
                if( !fLocalTrigger[h] && fFullTrigVec[h][i] == 0 )
                {
                    for( unsigned int j = 0; j < fNumSamples[h]; j++ )
                    {
                        fSamplesVec[h][i][j] = ( uint8_t )fGrIsuPeds[h][i][nstart + j];
                        fHiLo[h][i] = false;
                    }
                }
                // add background to existing trace
                else if( bAdd )
                {
                    unsigned int iTemp = 0;
                    for( unsigned int j = 0; j < fNumSamples[h]; j++ )
                    {
                        iTemp = ( unsigned int )fSamplesVec[h][i][j] + ( unsigned int )( fGrIsuPeds[h][i][nstart + j] - fdefaultPed );
                        if( iTemp > fDetectorGeo->getLowGainThreshold()[h] )
                        {
                            fSamplesVec[h][i][j] = fDetectorGeo->getLowGainThreshold()[h];
                        }
                        else
                        {
                            fSamplesVec[h][i][j] += ( uint8_t )( fGrIsuPeds[h][i][nstart + j] - fdefaultPed );
                        }
                    }
                }
            }
        }
    }
}


/*!
   this preliminary, grisu gives no trigger vector in current version
*/
void VGrIsuReader::setTrigger( vector<bool> iImage, vector<bool> iBorder )
{
    if( fFullTrigVec[fTelescopeID].size() != iImage.size() || fFullTrigVec[fTelescopeID].size() != iBorder.size() )
    {
        cout << "VGrIsuReader::setTrigger error: trigger/image/border vectors have different sizes ";
        cout << fFullTrigVec[fTelescopeID].size() << "\t" << iImage.size() << "\t" << iBorder.size() << endl;
    }
    for( unsigned int i = 0; i < fFullTrigVec[fTelescopeID].size(); i++ )
    {
        if( iImage[i] || iBorder[i] )
        {
            fFullTrigVec[fTelescopeID][i] = true;
            fNumberofFullTrigger[fTelescopeID]++;
        }
        else
        {
            fFullTrigVec[fTelescopeID][i] = false;
        }
    }
}


bool VGrIsuReader::hasLocalTrigger( unsigned int iTel )
{
    if( fDebug )
    {
        cout << "VGrIsuReader::hasLocalTrigger" << endl;
    }
    if( iTel < getLocalTrigger().size() )
    {
        return getLocalTrigger()[iTel];
    }
    return false;
}


unsigned int VGrIsuReader::getNTelLocalTrigger()
{
    unsigned int iNTrig = 0;
    for( unsigned int i = 0; i < getLocalTrigger().size(); i++ )
    {
        if( getLocalTrigger()[i] )
        {
            iNTrig++;
        }
    }
    return iNTrig;
}


void VGrIsuReader::setTraceFile( string i_Trf )
{
    if( fDebug )
    {
        cout << "void VGrIsuReader::setTraceFile " << i_Trf << endl;
    }
    fTraceFileName = i_Trf;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        readPixelMix( i );
    }
    if( i_Trf.size() > 0 )
    {
        // open trace library for background generation
        openTraceLibrary();
        // read pedestals again
        readPeds();
    }
}


TH1F* VGrIsuReader::getPedHisto( unsigned int itel, unsigned int ichannel )
{
    if( itel < fhPeds.size() )
    {
        if( ichannel < fhPeds[itel].size() )
        {
            return fhPeds[itel][ichannel];
        }
    }
    return 0;
}


void VGrIsuReader::setRandomDead( int iNC, int iNB )
{
    fMCNdead = iNC;
    fMCNdeadboard = iNB;
    
    if( fMCNdead )
    {
        cout << "VGrIsuReader: setting " << fMCNdead << " channels randomly dead" << endl;
    }
    if( fMCNdeadboard )
    {
        cout << "VGrIsuReader: setting " << fMCNdeadboard << "boards randomly dead" << endl;
    }
    
    // loop over all telescopes and set channel dead in all of them
    for( unsigned int iTel = 0; iTel < fNTel; iTel++ )
    {
        // set dead channel in read camera
        if( fMCNdead > 0 )
        {
            int iN = 0;
            int iPix = 0;
            int iNstep = 0;
            do
            {
                iPix = fRandomGen->Integer( ( int )fFullAnaVec[iTel].size() );
                // check if this one is alread dead
                if( fFullAnaVec[iTel][iPix] == 1 )
                {
                    fFullAnaVec[iTel][iPix] = -1;
                    iN++;
                }
                iNstep++;
                // avoid too large loop
                if( iNstep > 10000 )
                {
                    cout << "VGrIsuReader::readCamera: too many dead channels, break setting of random dead channels" << endl;
                    break;
                }
            }
            while( iN < fMCNdead );
        }
        // add a dead board (iNboard pixels in a row)
        unsigned int iNboard = 10;
        if( fMCNdeadboard > 0 )
        {
            int iN = 0;
            int iPix = 0;
            do
            {
                iPix = fRandomGen->Integer( ( int )fFullAnaVec[iTel].size() );
                iPix = iPix / 10 * 10;
                // don't care about if anything is already dead
                for( unsigned j = 0; j < iNboard; j++ )
                {
                    if( iPix + j < fFullAnaVec[iTel].size() )
                    {
                        fFullAnaVec[iTel][iPix + j] = -1;
                    }
                    fFullAnaVec[iTel][iPix] = -1;
                }
                iN++;
            }
            while( iN < fMCNdeadboard );
        }
    }
}


unsigned int VGrIsuReader::getGrisuVersion()
{
    // set stream to beginning of file
    is.seekg( 0, ios::beg );
    
    string is_line;
    string is_Temp;
    
    unsigned int iV = 0;
    
    sp = is.tellg();
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            if( is_line.substr( 0, 1 ) == "P" || is_line.substr( 0, 1 ) == "R" || is_line.substr( 0, 1 ) == "S" )
            {
                is.seekg( sp );
                return 0;
            }
            else if( is_line.substr( 0, 1 ) == "V" )
            {
                if( is_line.find( "." ) < is_line.size() )
                {
                    iV = atoi( is_line.substr( 1, is_line.find( "." ) ).c_str() ) * 100;
                }
                if( is_line.rfind( "." ) < is_line.size() )
                {
                    iV += atoi( is_line.substr( is_line.find( "." ) + 1, is_line.rfind( "." ) - is_line.find( "." ) - 1 ).c_str() ) * 10;
                }
                if( is_line.rfind( "." ) < is_line.size() )
                {
                    iV += atoi( is_line.substr( is_line.rfind( "." ) + 1, is_line.size() ).c_str() );
                }
                return iV;
            }
        }
    }
    return iV;
}


uint8_t VGrIsuReader::getNoiseSample( unsigned int iTel, uint32_t iHitID, unsigned int iSample, bool iNewTrace )
{
    if( iHitID < getMaxChannels() && iSample < fNumSamples[iTel] )
    {
        if( iNewTrace )
        {
            fNoiseTraceStart = fRandomGen->Integer( fGrIsuPedsN[iTel][iHitID] - fNumSamples[iTel] );
        }
        return fGrIsuPeds[iTel][iHitID][fNoiseTraceStart + iSample];
    }
    
    return 0;
}


vector< uint8_t >& VGrIsuReader::getNoiseVec( unsigned int iTel, uint32_t iHitID, bool iNewTrace )
{
    if( iTel < fGrIsuPeds.size() && iHitID < fGrIsuPeds[iTel].size() )
    {
        if( !bNoiseTraceFilled )
        {
            iNewTrace = true;
            bNoiseTraceFilled = true;
        }
        if( iNewTrace )
        {
            fNoiseTraceStart = fRandomGen->Integer( fGrIsuPeds[iTel][iHitID].size() - fNumSamples[iTel] );
            
            for( unsigned int j = 0; j < fNumSamples[iTel]; j++ )
            {
                fSamplesVec[iTel][iHitID][j] = ( uint8_t )fGrIsuPeds[iTel][iHitID][fNoiseTraceStart + j];
            }
        }
        return fSamplesVec[iTel][iHitID];
    }
    
    return vv8;
}


vector< vector< uint8_t > >& VGrIsuReader::getFullNoiseVec( unsigned int iTel )
{
    if( iTel < fGrIsuPeds.size() )
    {
        return fGrIsuPeds[iTel];
    }
    
    return v8;
}


vector< uint8_t >& VGrIsuReader::getFullNoiseVec( unsigned int iTel, unsigned int iChannel )
{
    if( iTel < fGrIsuPeds.size() && iChannel < fGrIsuPeds[iTel].size() )
    {
        return fGrIsuPeds[iTel][iChannel];
    }
    
    return vv8;
}


void VGrIsuReader::assignGrisuPeds( unsigned int n )
{
    fGrIsuPeds.clear();
    vector< uint8_t > a( n, 0 );
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        vector< vector<uint8_t> > i_pedChannelSample( fMaxChannels[i], a );
        fGrIsuPeds.push_back( i_pedChannelSample );
    }
}


uint8_t VGrIsuReader::getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace )
{
    if( channel < fMaxChannels[fTelescopeID] && sample < fSamplesVec[fTelescopeID][channel].size() )
    {
        return fSamplesVec[fTelescopeID][channel][sample];
    }
    
    iNewNoiseTrace = false;
    return 0;
}

void VGrIsuReader::setDefaultPed( double iD )
{
    fdefaultPed = iD;
    fdefaultPedUI = ( uint8_t )fdefaultPed;
}
