/*! \class VSpecialChannel
    \brief read and administrate list of special channels (e.g. L2 channels)


*/

#include <VSpecialChannel.h>

VSpecialChannel::VSpecialChannel( unsigned int iTelID )
{
    fTelID = iTelID;
    
    fDebug = false;
    
    // zombie until successfull reading of a file
    fIsZombie = true;
    
    fSpecialChannelFile = "";
}

/*!
   open a text file with special channel information
*/
bool VSpecialChannel::readSpecialChannels( int iRun, string ifile, string iDirectory )
{
    if( ifile.size() >  0 )
    {
        // special channel file is in camera geometry directory
        ifile = iDirectory + "/" + ifile;
        // open text file
        ifstream is;
        is.open( ifile.c_str(), ifstream::in );
        if( !is )
        {
            cout << "error reading special channels (L2 channels, etc) for telescope " << getTelID() + 1 << " from " << ifile << endl;
            return false;
        }
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "reading special channels (L2 channels, etc) from: " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << ifile << endl;
        
        string is_line;
        string is_temp;
        
        int iRunMin = 0;
        int iRunMax = 0;
        bool bRunFound = false;
        
        while( getline( is, is_line ) )
        {
            if( is_line.size() <= 0 )
            {
                continue;
            }
            if( is_line.substr( 0, 1 ) != "*" )
            {
                continue;
            }
            
            istringstream is_stream( is_line );
            is_stream >> is_temp;
            
            // check range of runs
            is_stream >> is_temp;
            if( is_temp == "RUN" )
            {
                if( bRunFound )
                {
                    break;
                }
                is_stream >> is_temp;
                iRunMin = atoi( is_temp.c_str() );
                if( iRun < iRunMin )
                {
                    continue;
                }
                is_stream >> is_temp;
                iRunMax = atoi( is_temp.c_str() );
                if( iRun > iRunMax )
                {
                    continue;
                }
                bRunFound = true;
                // is this too optimistic??
                fIsZombie = false;
            }
            if( bRunFound )
                // check telescope numbering
            {
                if( atoi( is_temp.c_str() ) != ( int )( getTelID() + 1 ) )
                {
                    continue;
                }
                
                ///////////////////////////////////////
                // special channel types
                //
                // known types:
                //              L2:     L2 signal to determine crate jitter
                //              HIGHQE: high efficiency PMT gain factor
                //
                ///////////////////////////////////////
                is_stream >> is_temp;
                
                ///////////////////////////////////////
                // L2 channels
                if( is_temp == "L2" )
                {
                    // number of L2 channels
                    is_stream >> is_temp;
                    int i_n = atoi( is_temp.c_str() );
                    getFADCstopTrigChannelID().clear();
                    getFADCstopTrigChannelID().assign( i_n, 1000000 );
                    if( i_n != ( int )getFADCstopTrigChannelID().size() )
                    {
                        cout << "WARNING: number of special channels disagrees with expected: " << i_n << "\t" << getFADCstopTrigChannelID().size() << endl;
                        cout << "\t (all channel numbers are ignored!)" << endl;
                    }
                    else
                    {
                        // now read L2 channels
                        for( unsigned int t = 0; t < getFADCstopTrigChannelID().size(); t++ )
                        {
                            is_stream >> is_temp;
                            getFADCstopTrigChannelID()[t] = ( unsigned int )( atoi( is_temp.c_str() ) );
                        }
                    }
                }
                //////////////////////////////////////////
                // channel status
                else if( is_temp == "STATUS" )
                {
                    if( !(is_stream>>std::ws).eof() )
                    {
                        // get status flag
                        unsigned int iStatusFlag = 1;
                        is_stream >> iStatusFlag;
                        fChannelStatus.clear();
                        while( !(is_stream>>std::ws).eof() )
                        {
                            is_stream >> is_temp;
                            unsigned int iC = ( unsigned int )atoi( is_temp.c_str() );
                            fChannelStatus[iC] = iStatusFlag;
                        }
                    }
                    
                }
                ///////////////////////////////////////
                // HIGHQE channels
                else if( is_temp == "HIGHQE" )
                {
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> is_temp;
                        unsigned int iC = ( unsigned int )atoi( is_temp.c_str() );
                        if( fHIGHQE_gainfactor.find( iC ) == fHIGHQE_gainfactor.end() )
                        {
                            if( !(is_stream>>std::ws).eof() )
                            {
                                is_stream >> is_temp;
                                double iF = atof( is_temp.c_str() );
                                fHIGHQE_gainfactor[iC] = iF;
                            }
                        }
                        else
                        {
                            // what is an appropriate reaction here?
                            cout << "VSpecialChannel::readSpecialChannels warning: found channel " << iC << " twice " << endl;
                        }
                    }
                }
            }
        }
        is.close();
        //////////////////////////////////////////////////////////////
        // print channels to screen
        
        // print L2 channels
        if( getFADCstopTrigChannelID().size() > 0 )
        {
            cout << "Telescope " << getTelID() + 1 << ": ";
            cout << "setting special channels with L2 signals: ";
            for( unsigned int t = 0; t < getFADCstopTrigChannelID().size(); t++ )
            {
                cout << getFADCstopTrigChannelID()[t] << "  ";
            }
            cout << endl;
        }
        // print HIQE_gainfactors
        if( fHIGHQE_gainfactor.size() > 0 )
        {
            cout << "Telescope " << getTelID() + 1 << ": ";
            cout << "setting " << fHIGHQE_gainfactor.size() << " HIGHQE gain factors: ";
            map< unsigned int, double >::iterator it;
            for( it = fHIGHQE_gainfactor.begin(); it != fHIGHQE_gainfactor.end(); it++ )
            {
                cout << ( *it ).first << " (" << ( *it ).second << ") ";
            }
            cout << endl;
        }
        // print status bit
        if( fChannelStatus.size() > 0 )
        {
            cout << "Telescope " << getTelID() + 1 << ": ";
            cout << "setting " << fChannelStatus.size() << " status bits: ";
            map< unsigned int, unsigned int >::iterator it;
            for( it = fChannelStatus.begin(); it != fChannelStatus.end(); it++ )
            {
                cout << ( *it ).first << " (" << ( *it ).second << ") ";
            }
            cout << endl;
        }
    }
    
    return true;
}

/*
 * read throughput correction from file
 */
bool VSpecialChannel::readThroughput( string iEpoch, string ifile, string iDirectory, unsigned int iNChannel )
{
       if( ifile.size() == 0 )
       {
            return true;
       }
       ifile = iDirectory + "/" + ifile;
       ifstream is;
        is.open( ifile.c_str(), ifstream::in );
       if( !is )
       {
            cout << "error reading throughput correction for telescope " << getTelID() + 1 << " from " << ifile << endl;
            return false;
       }
       cout << "Telescope " << getTelID() + 1 << ": ";
       cout << "reading throughput correction from: " << endl;
       cout << "Telescope " << getTelID() + 1 << ": ";
       cout << ifile << endl;

       string is_line;
       string is_temp;
       while( getline( is, is_line ) )
       {
            if( is_line.size() <= 0 )
            {
                    continue;
            }
            if( is_line.substr( 0, 1 ) != "*" )
            {
                    continue;
            }
            
            istringstream is_stream( is_line );
            is_stream >> is_temp;

            // check epoch
            is_stream >> is_temp;
            if( is_temp == "s" )
            {
                is_stream >> is_temp;
                if( is_temp == iEpoch )
                {
                    unsigned int t = 0;
                    double iSFactor = 1.;
                    do
                    {
                        is_stream >> iSFactor;
                        if( t == getTelID() )
                        {
                            break;
                        }
                        t++;     
                    } while( !is_stream.eof() );
                    if( iSFactor > 0. )
                    {
                        for( unsigned int i = 0; i < iNChannel; i++ )
                        {
                            if( fHIGHQE_gainfactor.find( i ) != fHIGHQE_gainfactor.end() )
                            {
                                 if( fHIGHQE_gainfactor[i] > 0. )
                                 {
                                      fHIGHQE_gainfactor[i] /= iSFactor;
                                 }
                            }
                            else
                            {
                                 fHIGHQE_gainfactor[i] = 1./iSFactor;
                            }
                        }
                        cout << "Telescope " << getTelID() + 1 << ": ";
                        cout << "applying throughput correction of " << iSFactor << endl;
                    }
                }
            }
       }
       return true;
}

void VSpecialChannel::reset()
{
    fFADCstopTrigChannelID.clear();
}

double VSpecialChannel::getHIGHQE_gainfactor( unsigned int chanID )
{
    if( fHIGHQE_gainfactor.find( chanID ) != fHIGHQE_gainfactor.end() )
    {
        return fHIGHQE_gainfactor[chanID];
    }
    
    return -1;
}

unsigned int VSpecialChannel::getChannelStatus( unsigned int chanID )
{
    if( fChannelStatus.find( chanID ) != fChannelStatus.end() )
    {
        return fChannelStatus[chanID];
    }
    
    return 1;
}
