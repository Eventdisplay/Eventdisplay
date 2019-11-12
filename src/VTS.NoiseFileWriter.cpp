/* \file VTS.NoiseFileWriter
 * \brief read a vbf file with pedestal events and
 *        write a GRISU style noise file
 *        to be used for Monte Carlo files without
 *        noise simulations
 *
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <VBFDataReader.h>

using namespace std;

void help()
{
    cout << endl;
    cout << "./VTS.NoiseFileWriter <ntel> <nchannel> <input vbf file> <output grisu files> <length of noise traces [samples]" << endl;
    cout << endl;
    cout << "   (a good length for the noise traces is 20000)" << endl;
    cout << endl;
}

int main( int argc, char* argv[] )
{
    if( argc != 6 )
    {
        help();
        exit( EXIT_SUCCESS );
    }
    unsigned int fNTel   = ( unsigned int )atoi( argv[1] );
    unsigned int fNChannel = ( unsigned int )atoi( argv[2] );
    string fPedestalFile = argv[3];
    string fGrisuFile    = argv[4];
    bool   fDebug        = false;
    unsigned int fTLength = ( unsigned int )atoi( argv[5] );
    
    cout << endl;
    cout << "VTS.NoiseFileWriter " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "---------------------------------------" << endl;
    cout << endl;
    
    //////////////////////////////////////////////////////////////////////////////////
    // read pedestal traces from VBF file
    // and fill them into a vector
    
    vector< vector< vector< uint8_t > > > fNoiseVector;
    
    vector< uint8_t > iTemp;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        vector< vector< uint8_t > > iiTemp;
        for( unsigned int j = 0; j < fNChannel; j++ )
        {
            iiTemp.push_back( iTemp );
        }
        fNoiseVector.push_back( iiTemp );
    }
    
    cout << "reading pedestal events from VBF file " << fPedestalFile << endl;
    
    /////////////////////////////////////////////////////
    unsigned int iNevents = 0;
    // open VBF file with noise traces
    VBFDataReader i_tempReader( fPedestalFile, 2, fNTel, fDebug );
    while( i_tempReader.getNextEvent() )
    {
        // telescopes
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            i_tempReader.setTelescopeID( i );
            
            // channels
            for( unsigned int j = 0; j < i_tempReader.getNumChannelsHit(); j++ )
            {
                unsigned int i_channelHitID = 0;
                try
                {
                    i_channelHitID = i_tempReader.getHitID( j );
                }
                catch( ... )
                {
                    continue;
                }
                if( i_channelHitID < fNChannel )
                {
                    for( unsigned int n = 0; n < i_tempReader.getNumSamples(); n++ )
                    {
                        if( i_tempReader.getSample( i_channelHitID, n ) > 0 )
                        {
                            fNoiseVector[i][j].push_back( i_tempReader.getSample( i_channelHitID, n ) );
                        }
                    }
                    
                }
            }
        }
        iNevents++;
    };
    
    cout << "successfully read noise vector for " << iNevents << " events " << endl;
    if( fNoiseVector.size() > 0 && fNoiseVector[0].size() > 0 )
    {
        cout << "\t length of noise vector: " << fNoiseVector[0][0].size() << endl;
        // check length of traces
        if( fTLength > fNoiseVector[0][0].size() )
        {
            cout << "\t warning: requested noise sample length too long" << endl;
            fTLength = fNoiseVector[0][0].size();
            cout << "\t will write only " << fTLength << " samples to file" << endl;
        }
    }
    else
    {
        cout << "Error: no noise vector filled" << endl;
        exit( EXIT_FAILURE );
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    // write noise files in grisu style
    
    // open output file
    ofstream i_offile;
    i_offile.open( fGrisuFile.c_str(), ios::out | ios::trunc );
    if( !i_offile.is_open() )
    {
        cout << "error opening GRISU NSB file: " << fGrisuFile << endl;
        cout << "exit..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << endl;
    cout << "writing noise file: " << fGrisuFile << endl;
    
    i_offile << "V 4.1.6" << endl;
    
    for( unsigned int i = 0; i < fNoiseVector.size(); i++ )
    {
        for( unsigned int j = 0; j < fNoiseVector[i].size(); j++ )
        {
            i_offile << "P 0 " << i + 1 << " " << j + 1;
            //                 for( unsigned int p = 0; p < fNoiseVector[i][j].size(); p++ )
            for( unsigned int p = 0; p < fTLength; p++ )
            {
                i_offile << " " << ( int )fNoiseVector[i][j][p];
            }
            i_offile << endl;
        }
    }
    i_offile.close();
    
    cout << "closed grisu file..." << endl;
}
