/*! \class VCameraRead
    \brief reading of telescope/camera geometry file (.cam and GrIsu .cfg files)

    all values in these class are static and should not change during a run

    observe that for EVNDISP version > 3 only .cfg are supported


*/

#include "VCameraRead.h"

VCameraRead::VCameraRead()
{
    fDebug = false;
    fNTel = 0;
    fTelID = 0;
    // configuration file type
    fCFGtype = 0;
    // pixel type
    fPixelType = 1;
    // default directory for cfg files
    fConfigDir = "../data/detector_geometry/";
    // default pedestal
    fDefPed = 20.;
    fFADCRange = 256;
    // default number of patches
    fNPatches = 91;
    // default GrIsu version
    fGrIsuVersion = 500;
    // low gain multiplier
    fLowGainIsSet = false;
    // coordinate transformers
    setCoordinateTransformer( 1., 1. );
    // default source type
    fsourcetype = 3;
    // default counting of pixels starts from 1
    fPixelCountingFromZero = false;
}


/*!
     if telescope ID is out of range, ID is set to zero (good idea??)
*/
bool VCameraRead::setTelID( unsigned int iTel )
{
    if( iTel < fNTel )
    {
        fTelID = iTel;
    }
    else
    {
        fTelID = 0;
        cout << " VCameraRead::setTelID() error: telescope ID out of range: requested " << iTel << ", current " << fTelID << ", ntel: " << fNTel << endl;
        //      cout << "                                setting ID to zero" << endl;
        cout << "   ...exiting" << endl;
        exit( EXIT_FAILURE );
        return false;
    }
    return true;
}


bool VCameraRead::initialize( unsigned int iNtel, vector< string > iCamera )
{
    if( fDebug )
    {
        cout << "VCameraRead::initialize " << iNtel << endl;
    }
    if( iNtel > iCamera.size() )
    {
        cout << "VCameraRead::initialize error: number of telescopes larger than camera vector ";
        cout << iNtel << "\t" << iCamera.size() << endl;
        return false;
    }
    fNTel = iNtel;
    resetTelVectors();
    fCameraName = iCamera;
    // get number of channels
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        string iCameraFile = fConfigDir + iCamera[i] + ".cam";
        ifstream inFileStream( iCameraFile.c_str() );
        if( !inFileStream )
        {
            cout << "VCameraRead::initialize camera geometry file not found: " << iCamera[i] << endl;
            cout << iCameraFile << endl;
            return false;
        }
        string i_Line;
        unsigned int zaehler = 0;
        while( getline( inFileStream, i_Line ) )
        {
            zaehler++;
        }
        fCNChannels[i] = zaehler;
        inFileStream.close();
    }
    resetCamVectors();
    return true;
}


/*!
     initialization necessary for .cam files

*/
bool VCameraRead::initialize( unsigned int i_Ntel, unsigned int i_Nchannel )
{
    if( fDebug )
    {
        cout << "VCameraRead::initialize " << i_Ntel << endl;
    }
    fNTel = i_Ntel;
    resetTelVectors();
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fCNChannels[i] = i_Nchannel;
        fCNSamples[i] = 128;                       // actual sample size is set later in VImageBaseAnalyzer (from reader)
        fSample_time_slice[i] = 2.;
    }
    resetCamVectors();
    return true;
}


bool VCameraRead::initialize( unsigned int i_Ntel, vector< unsigned int > i_Nchannel )
{
    if( fDebug )
    {
        cout << "VCameraRead::initialize " << i_Ntel << endl;
    }
    fNTel = i_Ntel;
    resetTelVectors();
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fCNChannels[i] = i_Nchannel[i];
        fCNSamples[i] = 128;                       // actual sample size is set later in VImageBaseAnalyzer (from reader)
        fSample_time_slice[i] = 2.;
    }
    resetCamVectors( false );
    return true;
}


bool VCameraRead::readCameraFile( string iCameraFile, unsigned int i_Nchannel )
{
    if( fDebug )
    {
        cout << "VCameraRead::readCameraFile " << iCameraFile << "\t" << i_Nchannel << endl;
    }
    initialize( 1, i_Nchannel );
    return readCameraFile( iCameraFile );
}


/*! read in camera geometry from .cam files

     \param iCameraFile file name for camera geometry file

     \attention
       - fine tuned to layout of .cam files
       - exit(-1) if camera geometry file not found

     \return true if reading succesful
*/
bool VCameraRead::readCameraFile( string iCameraFile )
{
    if( fDebug )
    {
        cout << "VCameraRead::readCameraFile " << iCameraFile << endl;
    }
    iCameraFile.insert( 0, fConfigDir );
    iCameraFile.append( ".cam" );
    std::ifstream inFileStream( iCameraFile.c_str() );
    if( !inFileStream )
    {
        std::cout << "VCameraRead::readCameraFile() error: camera geometry file not found: " << iCameraFile << std::endl;
        exit( EXIT_FAILURE );
        return false;
    }
    
    if( fCNChannels[fTelID] == 0 )
    {
        cout << "VCameraRead::readCameraFile() error: VCameraRead not initialized" << endl;
        return false;
    }
    
    string i_char;
    int i_ch;
    int i_mix;
    
    std::vector<int> i_temp;
    
    std::string i_Line;
    unsigned int zaehler = 0;
    while( getline( inFileStream, i_Line ) )
    {
        if( zaehler >= fCNChannels[fTelID] )
        {
            cout << "VCameraRead::readCameraFile() error: number of channels invalid" << endl;
            return false;
        }
        if( i_Line.size() > 0 )
        {
            std::istringstream is_stream( i_Line );
            is_stream >> i_char;
            is_stream >> i_ch;
            is_stream >> i_char;
            is_stream >> fXTube[fTelID][i_ch];
            is_stream >> i_char;
            is_stream >> fYTube[fTelID][i_ch];
            is_stream >> i_char;
            is_stream >> fRTube[fTelID][i_ch];
            // reading neighbours
            is_stream >> i_char;
            unsigned int j = 0;
            while( !( is_stream >> std::ws ).eof() && i_char.substr( 0, 1 ) == "N" )
            {
                if( j < fNeighbour[fTelID][i_ch].size() )
                {
                    is_stream >> fNeighbour[fTelID][i_ch][j];
                }
                is_stream >> i_char;
                j++;
            }
            // maybe there is some information about triggers and dead channels
            if( !( is_stream >> std::ws ).eof() && i_char.substr( 0, 4 ) == "TRIG" )
            {
                is_stream >> fTrigTube[fTelID][i_ch];
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> i_char;
                if( i_char.substr( 0, 3 ) == "ANA" )
                {
                    is_stream >> fAnaTube[fTelID][i_ch];
                }
            }
            // get convertion from MC tube numbering to real data tube numbering
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> i_char;
                if( i_char.substr( 0, 3 ) == "MIX" )
                {
                    is_stream >> i_mix;
                    fXim.push_back( ( unsigned int )i_mix );
                }
            }
            zaehler++;
        }
    }
    
    if( fXim.size() > 0 )
    {
        // revert mixing vector
        unsigned int i_max = 0;
        // ?? (GM)
        for( unsigned int i = 0; i < fXim.size(); i++ ) if( fXim[i] > i_max )
            {
                i_max = fXim[i];
            }
        if( i_max > fXim.size() )
        {
            cout << "error in mixing vector size: " << i_max << "\t" << fXim.size() << endl;
            exit( EXIT_FAILURE );
        }
        fMix.assign( fXim.size(), 0 );
        for( unsigned int i = 0; i < fXim.size(); i++ )
        {
            fMix[fXim[i]] = i;
        }
        // preliminary!!!
        fMix[0] = 0;
        //   for( unsigned int i = 0; i < fXim.size(); i++ ) cout << i << "\t" << fXim[i] << "\t" << fMix[i] << endl;
    }
    
    inFileStream.close();
    
    return true;
}


/*!

    - use .cam file for trace library option (preliminary)

    \par iFile name of .cfg GrIsu detector configuration file

*/
bool VCameraRead::readGrisucfg( string iFile, unsigned int iNTel )
{
    if( fDebug )
    {
        cout << "VCameraRead::readGrisucfg " << iFile << endl;
    }
    fNTel = iNTel;
    
    iFile.insert( 0, fConfigDir );
    std::ifstream inFileStream( iFile.c_str() );
    if( !inFileStream )
    {
        cout << "VCameraRead::readGrisucfg() error: grisu cfg file not found: " << iFile << endl;
        exit( EXIT_FAILURE );
        return false;
    }
    
    string iline = "";
    string i_char = "";
    unsigned int i_telID = 0;
    unsigned int i_telID_SIMU = 0;
    unsigned int i_chan = 0;
    unsigned int i_NN = 0;
    unsigned int i_NTelcfg = 0;
    
    while( getline( inFileStream, iline ) )
    {
        // '*' in line
        // to be more stable, allow whitespaces (' ' and '\t' ) before '*'
        std::size_t index = iline.find_first_not_of( " \t" );
        if( index == string::npos )
        {
            continue;    //line has only whitespace
        }
        if( iline[index] != '*' )
        {
            continue;    //first non-whitespace char is not '*'
        }
        
        //        if( iline.substr( 0, 1 ) != "*" && ) continue;
        istringstream i_stream( iline );
        // GrIsu version (4.0.0 = 400, 4.1.1 = 411)
        if( iline.find( "VERSN" ) < iline.size() )
        {
            if( iline.find( "." ) < iline.size() )
            {
                fGrIsuVersion = atoi( iline.substr( iline.find( "VERSN" ) + 5, iline.find( "." ) ).c_str() ) * 100;
            }
            if( iline.rfind( "." ) < iline.size() )
            {
                fGrIsuVersion += atoi( iline.substr( iline.find( "." ) + 1, iline.rfind( "." ) - iline.find( "." ) - 1 ).c_str() ) * 10;
            }
            if( iline.rfind( "." ) < iline.size() )
            {
                fGrIsuVersion += atoi( iline.substr( iline.rfind( "." ) + 1, iline.size() ).c_str() );
            }
            cout << endl;
            cout << "reading detector configuration file (GrIsu version " << fGrIsuVersion << ")" << endl;
        }
        // telescope file type
        if( iline.find( "CFGTYPE" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> fCFGtype;
            cout << "\t (file type " << fCFGtype << ")" << endl;
        }
        // nubmer of telescopes
        if( iline.find( "NBRTL" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_NTelcfg;
            if( i_NTelcfg < fNTel )
            {
                fNTel = i_NTelcfg;
            }
            if( fDebug )
            {
                cout << "VCameraRead: fNTel = " << fNTel << endl;
            }
            resetTelVectors();
        }
        // telescope IDs
        // (used if telescope number in simulation is different then in eventdisplay or cfg file)
        else if( iline.find( "TELID" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            i_stream >> i_telID_SIMU;
            if( i_telID_SIMU > 0 )
            {
                i_telID_SIMU -= 1;
            }
            fTelIDGrisu[i_telID] = i_telID_SIMU;
        }
        // telescope positions
        else if( iline.find( "TLLOC" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fTelXpos[i_telID];
                i_stream >> fTelYpos[i_telID];
                i_stream >> fTelZpos[i_telID];
            }
            if( fDebug )
            {
                cout << "VCameraRead telescope positions: " << fTelXpos[i_telID] << "\t" << fTelYpos[i_telID] << "\t" << fTelZpos[i_telID] << " (id=" << i_telID << ")" << endl;
            }
        }
        // FADC
        // (there is only one FADC record for all telescopes)
        else if( iline.find( "FADCS" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> fDefPed;
            i_stream >> fFADCRange;
            //////////////////////////////////////////////
            // IMPORTANT: IGNORING SAMPLE SETTINGS FROM CFG FILE HERE!!!! (GM)
            i_stream >> fCNSamples[0];
            if( fsourcetype != 1 && fsourcetype != 5 )
            {
                fCNSamples[0] = 128;
            }
            for( unsigned int i = 1; i < fNTel; i++ )
            {
                fCNSamples[i] = fCNSamples[0];
            }
            //////////////////////////////////////////////
            // hi/lo gains
            if( fGrIsuVersion >= 411 )
            {
                i_stream >> i_char;
                i_stream >> i_char;
                if( !i_stream.eof() )
                {
                    i_stream >> fLowGainMultiplier_Trace[0];
                    for( unsigned int i = 1; i < fNTel; i++ )
                    {
                        fLowGainMultiplier_Trace[i] = fLowGainMultiplier_Trace[0];
                    }
                    i_stream >> fLowGainActivator[0];
                    for( unsigned int i = 1; i < fNTel; i++ )
                    {
                        fLowGainActivator[i] = fLowGainActivator[0];
                    }
                    fLowGainIsSet = true;
                }
            }
        }
        // length of a time slice
        else if( iline.find( "SIMUL" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> fSample_time_slice[0];
            for( unsigned int i = 1; i < fNTel; i++ )
            {
                fSample_time_slice[i] = fSample_time_slice[0];
            }
        }
        // low gain multiplier (not a grisu line)
        else if( iline.find( "LOWMULT" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fLowGainMultiplier_Trace[i_telID];
            }
            fLowGainIsSet = true;
        }
        // mirror design
        else if( iline.find( "MIROR" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fTelRad[i_telID];
                i_stream >> fMirFocalLength[i_telID];
            }
            if( fDebug )
            {
                cout << "VCameraRead mirrors: " << i_telID + 1 << " " << fTelRad[i_telID] << " " << fMirFocalLength[i_telID] << endl;
            }
        }
        // number of pixels
        else if( iline.find( "CAMRA" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fCNChannels[i_telID];
                fTelID = i_telID;
            }
            // read max number of neighbours
            i_stream >> i_char;
            i_stream >> i_char;
            if( !i_stream.eof() )
            {
                i_stream >> fMaxNeighbour[i_telID];
            }
            
            if( fDebug )
            {
                cout << "VCameraRead: total number of channels: " << fCNChannels[i_telID] << " (" << i_telID + 1 << ")" << endl;
            }
            if( fCFGtype == 1 )
            {
                for( unsigned int t = 1; t < fCNChannels.size(); t++ )
                {
                    fCNChannels[t] = fCNChannels[0];
                }
            }
        }
        // camera rotation (not a original grisu line)
        else if( iline.find( "CAROT" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fCameraRotation[i_telID];
                fTelID = i_telID;
            }
        }
        // place scale change ( scale factor, offsets, camera rotation)
        else if( iline.find( "CAMPLATE" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            if( i_telID < fNTel )
            {
                i_stream >> fCameraScaleFactor[i_telID];
                i_stream >> fCameraCentreOffset[i_telID];
                i_stream >> fCameraRotation[i_telID];
            }
        }
        // tube stuff
        else if( iline.find( "PMPIX" ) < iline.size() )
        {
            if( fXTubeMM.size() == 0 )
            {
                resetCamVectors();
            }
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            i_stream >> i_chan;
            if( i_chan == 0 )
            {
                fPixelCountingFromZero = true;
            }
            i_chan = adjustPixelCounting( i_chan );
            if( fGrIsuVersion >= 400 )
            {
                i_stream >> i_char;
                i_stream >> i_char;
            }
            if( i_telID < fNTel && i_chan < fCNChannels[i_telID] )
            {
                i_stream >> fXTubeMM[i_telID][i_chan];
                i_stream >> fYTubeMM[i_telID][i_chan];
                i_stream >> fRTubeMM[i_telID][i_chan];
                i_stream >> i_char;
                i_stream >> i_char;     //geom. efficiency; q.e. table number
                i_stream >> i_char;
                i_stream >> i_char;     //single pe signal rel. amplitude fluctuation; signal rise time
                i_stream >> i_char;
                i_stream >> i_char;     //signal fall time; RC coupling
                if( fGrIsuVersion >= 600 )                      // new records:
                {
                    i_stream >> i_char >> i_char >> i_char >> i_char >> i_char; //signal rise time low gain; signal fall time low gain; RC coupling low gain; pulse lookup table number high gain; pulse lookup table number low gain
                }
                i_stream >> fTrigTube[i_telID][i_chan];
                i_stream >> fAnaTube[i_telID][i_chan];
                i_stream >> fTOff[i_telID][i_chan];
                i_stream >> fGain[i_telID][i_chan];
            }
        }
        // neighbour list
        else if( iline.find( "NGHBR" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            vector< unsigned int > i_pixN;
            i_stream >> i_telID;
            if( i_telID > 0 )
            {
                i_telID -= 1;
            }
            i_stream >> i_chan;
            i_chan = adjustPixelCounting( i_chan );
            i_stream >> i_NN;
            if( i_telID < fNTel  && i_chan < fCNChannels[i_telID] )
            {
                fNNeighbour[i_telID][i_chan] = i_NN;
                // hexagonal camera for VTS: edge pixel have a lower number of neighbour pixel
                if( fNNeighbour[i_telID][i_chan] < fMaxNeighbour[i_telID] )
                {
                    fEdgePixel[i_telID][i_chan] = true;
                }
            }
            if( i_NN > fMaxNeighbour[i_telID] )
            {
                cout << "VCameraRead::readGrisucfg() error: maximal number of neighbours wrong " << ( int )i_NN << "\t" << fMaxNeighbour[i_telID] << endl;
                continue;
            }
            if( i_telID < fNTel && i_chan < fCNChannels[i_telID] )
            {
                for( unsigned int j = 0; j < fMaxNeighbour[i_telID]; j++ )
                {
                    if( !i_stream.eof() )
                    {
                        i_stream >> fNeighbour[i_telID][i_chan][j];
                        // grisu starts at 1 with counting, evndisp at 0
                        fNeighbour[i_telID][i_chan][j] = adjustPixelCounting( fNeighbour[i_telID][i_chan][j] );
                    }
                    else
                    {
                        fNeighbour[i_telID][i_chan][j] = -1;
                    }
                }
            }
        }
        // pattern trigger
        else if( iline.find( "PSTON" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> fNPatches;
        }
        else if( iline.find( "PATCH" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            vector< int > i_tPatch( 19, 0 );
            for( int i = 0; i < 19; i++ )
            {
                i_stream >> i_tPatch[i];
            }
            fPatch.push_back( i_tPatch );
        }
        else if( iline.find( "PIXGB" ) < iline.size() )
        {
            if( fXTubeMM.size() == 0 )
            {
                resetCamVectors();
            }
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            fPixelType = ( unsigned int )atoi( i_char.c_str() );
            // square pixel allow 8 neighbours
            if( fPixelType == 3 )
            {
                fMaxNeighbour[i_telID] = 8;
            }
            i_stream >> i_char;
            i_stream >> i_char;
            float r = atof( i_char.c_str() );
            for( unsigned int c = 0; c < fCNChannels[i_telID]; c++ )
            {
                fRTubeMM[i_telID][c] = r;
            }
        }
        else if( iline.find( "PIXFI" ) < iline.size() )
        {
            i_stream >> i_char;
            i_stream >> i_char;
            i_stream >> i_char;
            readPixelFile( i_char );
        }
        
    }
    if( fCFGtype == 1 )
    {
        fillTelescopeVectors();
    }
    
    // convert from mm to deg
    convertMMtoDeg();
    
    // stretch and move camera
    stretchAndMoveCamera();
    
    // rotate the camera
    rotateCamera();
    
    // clean neighbour lists
    cleanNeighbourList();
    
    // set camera centre tube index
    setCameraCentreTubeIndex();
    
    if( fDebug )
    {
        cout << "END: VCameraRead::readGrisucfg " << iFile << endl;
    }
    
    return true;
}

/*
 *  find pixel which is located at the centre of the camera
 *
 *  Note: assume here that the centre of the coordinate system
 *        is the centre of the camera
 */
void VCameraRead::setCameraCentreTubeIndex()
{
    fCameraCentreTubeIndex.clear();
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        if( i < fXTube.size() && i < fYTube.size() && i < fRotXTube.size() && i < fRotYTube.size() )
        {
            double iCentreDistance = 1.e99;
            unsigned int iCentreTube = 0;
            for( unsigned int j = 0; j < fXTube[i].size() && j < fYTube[i].size(); j++ )
            {
                if( sqrt( fXTube[i][j]*fXTube[i][j] + fYTube[i][j]*fYTube[i][j] ) < iCentreDistance )
                {
                    iCentreDistance = sqrt( fXTube[i][j] * fXTube[i][j] + fYTube[i][j] * fYTube[i][j] );
                    iCentreTube = j;
                }
            }
            fCameraCentreTubeIndex.push_back( iCentreTube );
        }
        else
        {
            fCameraCentreTubeIndex.push_back( 999999 );
        }
    }
}

float VCameraRead::getMaximumFOV_deg()
{
    float i_degEdge = 0.;
    float iDist = 0.;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        if( i < fXTube.size() && i < fYTube.size() && i < fRTube.size() )
        {
            for( unsigned int j = 0; j < fXTube[i].size(); j++ )
            {
                iDist = sqrt( fXTube[i][j] * fXTube[i][j] + fYTube[i][j] * fYTube[i][j] ) + fRTube[i][j];
                if( iDist > i_degEdge )
                {
                    i_degEdge = iDist;
                }
            }
        }
    }
    return i_degEdge;
}

void VCameraRead::cleanNeighbourList()
{
    for( unsigned int i = 0; i < fNeighbour.size(); i++ )
    {
        for( unsigned int j = 0; j < fNeighbour[i].size(); j++ )
        {
            if( fNNeighbour[i][j] != fNeighbour[i][j].size() )
            {
                fNeighbour[i][j].erase( fNeighbour[i][j].begin() + fNNeighbour[i][j], fNeighbour[i][j].end() );
            }
        }
    }
}


/*
    read external pixel file
*/
void VCameraRead::readPixelFile( string iFile )
{
    if( fDebug )
    {
        cout << "VCameraRead::readPixelFile " << iFile << endl;
    }
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VCameraRead error: cannot find file with pixel information: " << iFile << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    string itemp;
    unsigned int i_chan = 0;
    unsigned int i_telID = 0;
    unsigned int i_NN = 0;
    
    while( getline( is, is_line ) )
    {
        if( is_line.substr( 0, 1 ) != "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        is_stream >> itemp;
        is_stream >> itemp;
        if( itemp == "PIXLC" )
        {
            is_stream >> i_chan;
            i_chan = adjustPixelCounting( i_chan );
            if( i_telID < fNTel  && i_chan < fCNChannels[i_telID] )
            {
                is_stream >> fXTubeMM[i_telID][i_chan];
                is_stream >> fYTubeMM[i_telID][i_chan];
                is_stream >> i_NN;
                if( i_telID < fNTel  && i_chan < fCNChannels[i_telID] )
                {
                    fNNeighbour[i_telID][i_chan] = i_NN;
                }
                if( i_NN > fMaxNeighbour[i_telID] )
                {
                    cout << "VCameraRead::readGrisucfg() error: maximal number of neighbours wrong ";
                    cout << ( int )i_NN << "\t" << fMaxNeighbour[i_telID] << endl;
                    continue;
                }
                for( unsigned int j = 0; j < fMaxNeighbour[i_telID]; j++ )
                {
                    if( !( is_stream >> std::ws ).eof() )
                    {
                        if( j < fNeighbour[i_telID][i_chan].size() )
                        {
                            is_stream >> fNeighbour[i_telID][i_chan][j];
                        }
                        else
                        {
                            int a = 0;
                            is_stream >> a;
                            fNeighbour[i_telID][i_chan].push_back( a );
                        }
                        fNeighbour[i_telID][i_chan][j] = adjustPixelCounting( fNeighbour[i_telID][i_chan][j] );
                    }
                    else
                    {
                        if( j < fNeighbour[i_telID][i_chan].size() )
                        {
                            fNeighbour[i_telID][i_chan][j] = -1;
                        }
                        else
                        {
                            fNeighbour[i_telID][i_chan].push_back( -1 );
                        }
                    }
                }
            }
        }
    }
}


void VCameraRead::fillTelescopeVectors()
{
    if( fDebug )
    {
        cout << "VCameraRead::fillTelescopeVectors()" << endl;
    }
    for( unsigned int t = 1; t < fTelRad.size(); t++ )
    {
        fTelRad[t] = fTelRad[0];
    }
    for( unsigned int t = 1; t < fMirFocalLength.size(); t++ )
    {
        fMirFocalLength[t] = fMirFocalLength[0];
    }
    for( unsigned int t = 1; t < fNMirrors.size(); t++ )
    {
        fNMirrors[t] = fNMirrors[0];
    }
    for( unsigned int t = 1; t < fMirrorArea.size(); t++ )
    {
        fMirrorArea[t] = fMirrorArea[0];
    }
    for( unsigned int t = 1; t < fCameraFieldofView.size(); t++ )
    {
        fCameraFieldofView[t] = fCameraFieldofView[0];
    }
    for( unsigned int t = 1; t < fCNChannels.size(); t++ )
    {
        fCNChannels[t] = fCNChannels[0];
    }
    for( unsigned int t = 1; t < fCameraScaleFactor.size(); t++ )
    {
        fCameraScaleFactor[t] = fCameraScaleFactor[0];
    }
    for( unsigned int t = 1; t < fCameraCentreOffset.size(); t++ )
    {
        fCameraCentreOffset[t] = fCameraCentreOffset[0];
    }
    for( unsigned int t = 1; t < fCameraRotation.size(); t++ )
    {
        fCameraRotation[t] = fCameraRotation[0];
    }
    for( unsigned int t = 1; t < fXTubeMM.size(); t++ )
    {
        for( unsigned int i = 0; i < fXTubeMM[t].size(); i++ )
        {
            if( i < fXTubeMM[0].size() )
            {
                fXTubeMM[t][i] = fXTubeMM[0][i];
            }
            if( i < fYTubeMM[0].size() )
            {
                fYTubeMM[t][i] = fYTubeMM[0][i];
            }
            if( i < fRTubeMM[0].size() )
            {
                fRTubeMM[t][i] = fRTubeMM[0][i];
            }
            if( i < fTrigTube[0].size() )
            {
                fTrigTube[t][i] = fTrigTube[0][i];
            }
            if( i < fAnaTube[0].size() )
            {
                fAnaTube[t][i] = fAnaTube[0][i];
            }
            if( i < fTOff[0].size() )
            {
                fTOff[t][i] = fTOff[0][i];
            }
            if( i < fGain[0].size() )
            {
                fGain[t][i] = fGain[0][i];
            }
        }
    }
    for( unsigned int t = 1; t < fNNeighbour.size(); t++ )
    {
        for( unsigned int i = 0; i < fNNeighbour[t].size(); i++ )
        {
            if( i < fNNeighbour[0].size() )
            {
                fNNeighbour[t][i] = fNNeighbour[0][i];
            }
            if( t < fEdgePixel.size() && i < fEdgePixel[0].size() )
            {
                fEdgePixel[t][i] = fEdgePixel[0][i];
            }
            for( unsigned int j = 0; j < fNeighbour[0][i].size(); j++ )
            {
                if( j < fNeighbour[t][i].size() )
                {
                    fNeighbour[t][i][j] = fNeighbour[0][i][j];
                }
                else
                {
                    fNeighbour[t][i].push_back( fNeighbour[0][i][j] );
                }
            }
        }
    }
}

//void VCameraRead::print( bool bDetailed ){};
void VCameraRead::print( bool bDetailed )
{
    if( fDebug )
    {
        cout << "VCameraRead::print()" << endl;
    }
    cout << "Number of telescopes: " << fNTel << endl;
    cout << "telescope positions: " << endl;
    for( unsigned int i = 0; i < fTelXpos.size(); i++ )
    {
        cout << "Telescope " << i + 1 << " (type " << fTelType[i] << "): ";
        cout << "X: " << fTelXpos[i] << " [m], ";
        cout << "Y: " << fTelYpos[i] << " [m], ";
        cout << "Z: " << fTelZpos[i] << " [m], ";
        cout << "R: " << fTelRad[i] << " [m], ";
        cout << "Focal length: " << fMirFocalLength[i] << " [m], ";
        cout << "FOV: " << fCameraFieldofView[i] << " [deg], ";
        cout << "NChannel: " << fCNChannels[i];
        // number of samples are determined with the first event, possible not known yet
        //        cout << "NSamples: " << fCNSamples[i] << endl;
        if( fXTube.size() > 0 )
        {
            if( bDetailed )
            {
                cout << "\t Tube geometry: " << fXTube[i].size() << "\t" << fYTube[i].size() << "\t" << fRTube[i].size() << endl;
                for( unsigned int j = 0; j < fXTube[i].size(); j++ )
                {
                    cout << "\t\t" << i << "\t" << j << "\t";
                    cout << setw( 6 ) << setprecision( 4 ) << fXTube[i][j] << "\t";
                    cout << setw( 6 ) << setprecision( 4 ) << fYTube[i][j] << "\t";
                    cout << setw( 6 ) << setprecision( 4 ) << fRTube[i][j] << "\t";
                    cout << fTrigTube[i][j] << "\t" << fAnaTube[i][j] << "\t";
                    cout << fTOff[i][j] << "\t" << fGain[i][j] << "\t";
                    cout << "N ";
                    cout << "(" << fNNeighbour[i][j] << "," << fNeighbour[i][j].size() << ") ";
                    for( unsigned int k = 0; k < fNeighbour[i][j].size(); k++ )
                    {
                        cout << fNeighbour[i][j][k] << "\t";
                    }
                    cout << endl;
                }
            }
            else
            {
                cout << endl;
            }
        }
    }
}


void VCameraRead::convertMMtoDeg()
{
    if( fDebug )
    {
        cout << "VCameraRead::convertMMtoDeg" << endl;
    }
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        for( unsigned int j = 0; j < fXTube[i].size(); j++ )
        {
            // transform coordinates
            fXTubeMM[i][j] *= fCoordinateTransformerX;
            fYTubeMM[i][j] *= fCoordinateTransformerY;
            
            fXTube[i][j] = atan2( ( double )fXTubeMM[i][j] / 1000., ( double )fMirFocalLength[i] ) * 45. / atan( 1. );
            fYTube[i][j] = atan2( ( double )fYTubeMM[i][j] / 1000., ( double )fMirFocalLength[i] ) * 45. / atan( 1. );
            fRTube[i][j] = atan2( ( double )fRTubeMM[i][j] / 1000., ( double )fMirFocalLength[i] ) * 45. / atan( 1. );
        }
    }
}


/*!
 *  stretch and move camera
 */
void VCameraRead::stretchAndMoveCamera()
{
    if( fDebug )
    {
        cout << "VCameraRead::stretchAndMoveCamera" << endl;
    }
    
    cout << "camera plate scaled by";
    // stretch
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        // print only scale factors significantly different to 1
        if( TMath::Abs( fCameraScaleFactor[i] - 1. ) > 1.e-4 )
        {
            cout << " T" << i + 1 << ": " << fCameraScaleFactor[i];
            for( unsigned int j = 0; j < fXTube[i].size(); j++ )
            {
                fXTube[i][j] *= fCameraScaleFactor[i];
                fXTube[i][j] += fCameraCentreOffset[i];
                fYTube[i][j] *= fCameraScaleFactor[i];
                fYTube[i][j] += fCameraCentreOffset[i];
            }
        }
    }
    cout << endl;
}


/*!
 *  rotate camera counterclockwise
 */
void VCameraRead::rotateCamera()
{
    if( fDebug )
    {
        cout << "VCameraRead::rotateCamera " << endl;
    }
    
    cout << "camera rotation (in deg) of ";
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        // print only rotations significantly different to zero
        if( TMath::Abs( fCameraRotation[i] ) > 1.e-3 || fNTel < 5 )
        {
            cout << " T" << i + 1 << ": "  << fCameraRotation[i];
        }
        double iR = fCameraRotation[i] * TMath::DegToRad();
        if( i < fXTube.size() && i < fYTube.size() && i < fRotXTube.size() && i < fRotYTube.size() )
        {
            for( unsigned int j = 0; j < fXTube[i].size(); j++ )
            {
                if( j < fRotXTube[i].size() && j < fRotYTube[i].size() && j < fYTube[i].size() )
                {
                    fRotXTube[i][j] = fXTube[i][j] * cos( iR ) - fYTube[i][j] * sin( iR );
                    fRotYTube[i][j] = fXTube[i][j] * sin( iR ) + fYTube[i][j] * cos( iR );
                }
                else
                {
                    cout << "VCameraRead::rotateCamera() error: invalid vector sizes (expeced " << i << "," << j << "): ";
                    cout << fXTube[i].size() << "\t" << fYTube[i].size() << "\t" << fRotXTube[i].size() << "\t" << fRotYTube[i].size() << endl;
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
            }
        }
        else
        {
            cout << "VCameraRead::rotateCamera() error: invalid vector sizes (expeced " << i << "): ";
            cout << fXTube.size() << "\t" << fYTube.size() << "\t" << fRotXTube.size() << "\t" << fRotYTube.size() << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    cout << endl;
}


/*!
 */
void VCameraRead::resetTelVectors()
{
    fCameraName.assign( fNTel, "camera" );
    fCameraScaleFactor.assign( fNTel, 1. );
    fCameraCentreOffset.assign( fNTel, 0. );
    fCameraRotation.assign( fNTel, 0. );
    fTelIDGrisu.clear();
    fTelType.clear();
    //////////////////////////////////////////////////////////////////////////
    // each telescope is different
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fTelType.push_back( i + 1 );
    }
    // set fTelType to same value for similar telescopes
    //////////////////////////////////////////////////////////////////////////
    fTelXpos.assign( fNTel, 0. );
    fTelYpos.assign( fNTel, 0. );
    fTelZpos.assign( fNTel, 0. );
    fTelRad.assign( fNTel, 7. );
    fCNChannels.assign( fNTel, 0 );
    fCNSamples.assign( fNTel, 128 );               // actual sample size is set later in VImageBaseAnalyzer
    fSample_time_slice.assign( fNTel, 2. );
    fMirFocalLength.assign( fNTel, 12. );
    fNMirrors.assign( fNTel, 0 );
    fMirrorArea.assign( fNTel, 0. );
    fCameraFieldofView.assign( fNTel, 3.5 );
    fLowGainMultiplier_Trace.assign( fNTel, 6.0 );
    fLowGainActivator.assign( fNTel, 255 );
    //maximal number of neighbours is 6 (for circular pixel type)
    //maximal number of neighbours is 8 (for square pixel type) - pSCT
    fMaxNeighbour.assign( fNTel, 6 );
    // set default values for array of four telescopes
    //  later this values are overwritten by the values from the .cfg file
    if( fNTel == 4 )
    {
        fTelXpos[0] = 0.;
        fTelXpos[1] = -69.282;
        fTelXpos[2] = 69.282;
        fTelXpos[3] = 0.;
        fTelYpos[0] = 0.;
        fTelYpos[1] = 40.;
        fTelYpos[2] = 40.;
        fTelYpos[3] = -80.;
        fTelZpos[0] = 0.;
        fTelZpos[1] = 0.;
        fTelZpos[2] = 0.;
        fTelZpos[3] = 0.;
    }
}


/*!
   2d vectors [telescope][pixel]
*/
void VCameraRead::resetCamVectors( bool bMaxN )
{
    vector< float > i_tel;
    vector< int > i_pix;
    vector< unsigned int > i_pixN;
    vector< vector< int > > i_N;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        i_tel.assign( fCNChannels[i], 0. );
        fXTube.push_back( i_tel );
        fYTube.push_back( i_tel );
        fRTube.push_back( i_tel );
        fRotXTube.push_back( i_tel );
        fRotYTube.push_back( i_tel );
        fXTubeMM.push_back( i_tel );
        fYTubeMM.push_back( i_tel );
        fRTubeMM.push_back( i_tel );
        fTOff.push_back( i_tel );
        i_tel.assign( fCNChannels[i], 1. );
        fGain.push_back( i_tel );
        i_pix.assign( fCNChannels[i], 1 );
        i_pixN.assign( fCNChannels[i], 1 );
        fTrigTube.push_back( i_pix );
        fAnaTube.push_back( i_pix );
    }
    resetNeighbourLists( bMaxN );
}


void VCameraRead::resetNeighbourLists( bool bMaxN )
{
    if( fDebug )
    {
        cout << "VCameraRead::resetNeighbourLists() " << endl;
    }
    vector< int > i_pix;
    vector< unsigned int > i_pixN;
    vector< vector< int > > i_N;
    fNeighbour.clear();
    fNNeighbour.clear();
    vector< bool > b_pix;
    fEdgePixel.clear();
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        if( bMaxN )
        {
            i_pix.assign( fMaxNeighbour[i], -1 );
        }
        i_N.assign( fCNChannels[i], i_pix );
        if( bMaxN )
        {
            i_pixN.assign( fCNChannels[i], 1 );
        }
        else
        {
            i_pixN.assign( fCNChannels[i], 0 );
        }
        fNeighbour.push_back( i_N );
        fNNeighbour.push_back( i_pixN );
        b_pix.assign( fCNChannels[i], false );
        fEdgePixel.push_back( b_pix );
    }
}

/*

     use positions and tupe radius to prepare a NN neighbour list

     Note: for VTS, the list of neighbours is given in the configuration file

     (currently used for DST files only)

*/
bool VCameraRead::makeNeighbourList( vector< float > iNeighbourDistanceFactor,
                                     vector< bool > iSquarePixels )
{
    ////////////////////////////////////
    // loop over all telescopes in list
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        /////////////////////////////////////////
        // check if all tubes have the same size
        // (algorithm works only for tubes with the same size)
        double i_TubeRadius_0 = 0.;
        
        if( getTubeRadius_MM( i ).size() > 0 )
        {
            i_TubeRadius_0 = getTubeRadius_MM( i )[0];
            for( unsigned int j = 0; j < getTubeRadius_MM( i ).size(); j++ )
            {
                if( TMath::Abs( getTubeRadius_MM( i )[j] - i_TubeRadius_0 ) > 1.e-2 )
                {
                    cout << "VCameraRead::makeNeighbourList error: found tubes with different sizes in camera " << i + 1 << endl;
                    cout << "(algorithm works only for cameras where all pixels have the same size)" << endl;
                    return false;
                }
            }
        }
        else
        {
            continue;
        }
        //////////////////////////////////////
        // get minimum distance between tubes
        double iTubeDistance_min = 1.e5;
        for( unsigned int j = 0; j < getTubeRadius_MM( i ).size(); j++ )
        {
            for( unsigned int k = 0; k < getTubeRadius_MM( i ).size(); k++ )
            {
                if( j == k )
                {
                    continue;
                }
                double itemp = sqrt( ( getX_MM( i )[j] - getX_MM( i )[k] ) * ( getX_MM( i )[j] - getX_MM( i )[k] )
                                     + ( getY_MM( i )[j] - getY_MM( i )[k] ) * ( getY_MM( i )[j] - getY_MM( i )[k] ) );
                if( itemp < iTubeDistance_min )
                {
                    iTubeDistance_min = itemp;
                }
            }
        }
        iTubeDistance_min *= 0.5;
        
        if( i < iNeighbourDistanceFactor.size() )
        {
            iTubeDistance_min *= iNeighbourDistanceFactor[i];
        }
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Now find all neighbouring pixels
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        double maxDist = 1.1; // Distance in diameters to adjacent pixel
        // In the case of square pixels, increase the distance to go over the gaps in square cameras.
        if( iSquarePixels[i] )
        {
            maxDist = 1.9;
        }
        vector <unsigned int> nnForEdgePixels;
        nnForEdgePixels.assign( getTubeRadius_MM( i ).size(), 0 );
        for( unsigned int j = 0; j < getTubeRadius_MM( i ).size(); j++ )
        {
            for( unsigned int k = 0; k < j; k++ )
            {
                // In filling the list of neighbours for cleaning, we allow also diagonal neighbours for square pixels, hence the sqrt(2)
                double itemp = sqrt( ( getX_MM( i )[j] - getX_MM( i )[k] ) * ( getX_MM( i )[j] - getX_MM( i )[k] )
                                     + ( getY_MM( i )[j] - getY_MM( i )[k] ) * ( getY_MM( i )[j] - getY_MM( i )[k] ) );
                if( itemp < 2.1 * sqrt( 2 ) * iTubeDistance_min ) // Multiply by 2 because iTubeDistance_min is the pixel radius
                {
                    if( getAnaPixel( i )[k] > 0 || getAnaPixel( i )[k] == -822 )
                    {
                        // number of neighbours
                        fNNeighbour[i][j]++;
                        // list of neighbours
                        fNeighbour[i][j].push_back( k );
                    }
                    if( getAnaPixel( i )[j] > 0 || getAnaPixel( i )[j] == -822 )
                    {
                        // number of neighbours
                        fNNeighbour[i][k]++;
                        // list of neighbours
                        fNeighbour[i][k].push_back( j );
                    }
                }
                
                
                // In counting neighbours to determine edge pixels, consider only directly adjacent neighbours (no diagonals)
                // The tolerance is set to 20% of the pixel diameter due to the CHEC-S and ASTRI cameras,
                // where pixels in the same row/column were observed to be shifted in x/y positions due to the camera curvature.
                float maxShift = 0.2 * 2 * iTubeDistance_min; // Multiply by 2 because iTubeDistance_min is the pixel radius
                if( ( iSquarePixels[i] && ( abs( getX_MM( i )[j] - getX_MM( i )[k] ) < maxShift || abs( getY_MM( i )[j] - getY_MM( i )[k] ) < maxShift ) )
                        || !iSquarePixels[i] )
                {
                    if( itemp < 2 * maxDist * iTubeDistance_min )
                    {
                        if( getAnaPixel( i )[k] > 0 || getAnaPixel( i )[k] == -822 )
                        {
                            // number of neighbours
                            nnForEdgePixels[j]++;
                        }
                        if( getAnaPixel( i )[j] > 0 || getAnaPixel( i )[j] == -822 )
                        {
                            // number of neighbours
                            nnForEdgePixels[k]++;
                        }
                    }
                }
            }
        }
        // calculate maximum number of neighbours for this pixel type
        fMaxNeighbour[i] = 0;
        unsigned int maxForEdgePixels = 0;
        for( unsigned int j = 0; j < fNNeighbour[i].size(); j++ )
        {
            if( fNNeighbour[i][j] > fMaxNeighbour[i] )
            {
                fMaxNeighbour[i] = fNNeighbour[i][j];
            }
        }
        for( unsigned int j = 0; j < nnForEdgePixels.size(); j++ )
        {
            if( nnForEdgePixels[j] > maxForEdgePixels )
            {
                maxForEdgePixels = nnForEdgePixels[j];
            }
        }
        
        /////////////////////////////////////////////////////////////
        // determine edge pixels
        /////////////////////////////////////////////////////////////
        
        // The number of adjacent neighbours is 6 or 4 for hexagonal and square cameras respectivley
        if( maxForEdgePixels == 6 || maxForEdgePixels == 4 )
        {
            for( unsigned int j = 0; j < nnForEdgePixels.size(); j++ )
            {
                if( nnForEdgePixels[j] < maxForEdgePixels )
                {
                    fEdgePixel[i][j] = true;
                }
                else
                {
                    fEdgePixel[i][j] = false;
                }
            }
        }
        else
        {
            cout << "\nVCameraRead::makeNeighbourList() error: maxForEdgePixels = " << maxForEdgePixels
                 << " (expeced 6 or 4 for hexagonal or square pixels)" << endl;
            cout << "Tel now - " << fTelType[i] << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    return true;
}


float VCameraRead::getOuterEdgeDistance( unsigned int i )
{
    if( i > getX().size() - 1 )
    {
        cout << "VCameraRead::getOuterEdgeDistance(): channel index out of range, " << i << "\t" << getX().size() << endl;
        return 0.;
    }
    
    double iDist = sqrt( getX()[i] * getX()[i] + getY()[i] * getY()[i] );
    iDist += fRTube[fTelID][i];
    
    return iDist;
}


bool VCameraRead::setLowGainMultiplier_Trace( unsigned int iTel, double ival )
{
    if( iTel < fLowGainMultiplier_Trace.size() )
    {
        fLowGainMultiplier_Trace[iTel] = ival;
    }
    else
    {
        cout << "VCameraRead::setLowGainMultiplier_Trace: invalid low gain multiplier, set value to 1" << endl;
        return false;
    }
    return true;
}

bool VCameraRead::setLowGainThreshold( unsigned int iTel, unsigned int ival )
{
    if( iTel < fLowGainActivator.size() )
    {
        fLowGainActivator[iTel] = ival;
    }
    else
    {
        cout << "VCameraRead::setLowGainThreshold: invalid low gain threshold, set value to 0" << endl;
        return false;
    }
    return true;
}

/*
 * get position of telescope type in list of telescope types
 *
 * (this is a slow routine, don't call it too often)
 *
*/
unsigned int VCameraRead::getTelType_Counter( ULong64_t iTelType )
{
    set< ULong64_t > s;
    
    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
        s.insert( fTelType[i] );
    }
    
    unsigned int z = 0;
    for( set< ULong64_t >::iterator i_s = s.begin(); i_s != s.end(); i_s++ )
    {
        if( *i_s == iTelType )
        {
            return z;
        }
        z++;
    }
    
    // unsuccessfull search
    return 9999;
}


vector<ULong64_t> VCameraRead::getTelType_list()
{
    vector<ULong64_t> t;
    
    set< ULong64_t > s;
    
    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
        s.insert( fTelType[i] );
    }
    
    set< ULong64_t >::iterator it_s;
    for( it_s = s.begin(); it_s != s.end(); it_s++ )
    {
        t.push_back( *it_s );
    }
    
    return t;
}

/*!

     read varios detector parameters from the DB

*/
bool VCameraRead::readDetectorGeometryFromDB( string iDBStartTime, bool iReadRotationsFromDB )
{
    if( fDebug )
    {
        cout << "VCameraRead::readDetectorGeometryFromDB" << endl;
    }
    cout << "VCameraRead::readDetectorGeometryFromDB for " << iDBStartTime;
    if( iReadRotationsFromDB )
    {
        cout << " (read rotations)";
    }
    cout << endl;
    
    if( iDBStartTime.size() < 8 )
    {
        cout << "VCameraRead::readDetectorGeometryFromDB error: no valid SQL data for getting DB detector geometry: " << iDBStartTime << endl;
        return false;
    }
    
    // read camera rotations from DB
    if( iReadRotationsFromDB )
    {
        stringstream iTempS;
        iTempS << getDBServer() << "/VOFFLINE";
        char c_query[800];
        sprintf( c_query, "select telescope_id, version, pmt_rotation from tblPointing_Monitor_Camera_Parameters where start_date <= \"%s\" AND end_date > \"%s\" ", iDBStartTime.substr( 0, 10 ).c_str(), iDBStartTime.substr( 0, 10 ).c_str() );
        
        //std::cout<<"VCameraRead::readDetectorGeometryFromDB "<<std::endl;
        VDB_Connection my_connection( iTempS.str().c_str() , "readonly", "" ) ;
        if( !my_connection.Get_Connection_Status() )
        {
            cout << "VCameraRead: failed to connect to database server" << endl;
            cout << "\t server: " <<  iTempS.str() << endl;
            return false;
        }
        if( !my_connection.make_query( c_query ) )
        {
            return false;
        }
        TSQLResult* db_res = my_connection.Get_QueryResult();
        
        
        int iNRows = db_res->GetRowCount();
        vector< int > iVersion( fCameraRotation.size(), -99 );
        for( int j = 0; j < iNRows; j++ )
        {
            TSQLRow* db_row = db_res->Next();
            if( !db_row )
            {
                continue;
            }
            
            int itelID = -99;
            double iRot = -99.;
            if( db_row->GetField( 0 ) )
            {
                itelID = atoi( db_row->GetField( 0 ) );
            }
            // check entry version
            if( itelID >= 0 && itelID < ( int )iVersion.size() )
            {
                if( db_row->GetField( 1 ) && atoi( db_row->GetField( 1 ) ) > iVersion[itelID] )
                {
                    if( db_row->GetField( 2 ) )
                    {
                        iRot = atof( db_row->GetField( 2 ) ) * TMath::RadToDeg();
                        if( itelID < ( int )fCameraRotation.size() )
                        {
                            fCameraRotation[itelID] = -1.* iRot;
                        }
                    }
                    iVersion[itelID] = atoi( db_row->GetField( 1 ) );
                }
            }
        }
        //       i_DB->Close();
    }
    
    cout << "\t (rotations from DB [deg]: ";
    for( unsigned int i = 0; i < fCameraRotation.size(); i++ )
    {
        cout << " T" << i + 1 << ": " << fCameraRotation[i];
    }
    cout << ")" << endl;
    
    // rotate the camera
    rotateCamera();
    if( fDebug )
    {
        cout << "rotateCamera() Finished" << endl;
    }
    return true;
}

bool VCameraRead::setLengthOfSampleTimeSlice( unsigned int iTelID, float iSample_time_slice )
{
    if( iTelID < fSample_time_slice.size() )
    {
        fSample_time_slice[iTelID] = iSample_time_slice;
        return true;
    }
    return false;
}

vector< unsigned int > VCameraRead::getNumChannelVector()
{
    vector< unsigned int > iN;
    for( unsigned int i = 0; i < fXTube.size(); i++ )
    {
        iN.push_back( fXTube[fTelID].size() );
    }
    return iN;
}

/*
 * adjust pixel IDs
 * (e.g. handle differences between counting from zero or one)
 *
 * default:
 * - expect cfg files to start from 1 (Grisu convention)
 * - change to starting from zero (evndisp convention)
 */
unsigned int VCameraRead::adjustPixelCounting( unsigned int i_chan )
{
    if( fPixelCountingFromZero )
    {
        return i_chan;
    }
    
    return i_chan - 1;
}

