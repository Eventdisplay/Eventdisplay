/*! \class VVirtualDataReader
    \brief wrapper classe for data reading


*/

#include <VVirtualDataReader.h>

// VVirtualDataReader::VVirtualDataReader() throw (VFileException)
VVirtualDataReader::VVirtualDataReader()
{

}

double VVirtualDataReader::getSample_double( unsigned channel, unsigned sample, bool iNewNoiseTrace )
{
    if( has16Bit() )
    {
        return ( double )getSample16Bit( channel, sample, iNewNoiseTrace );
    }
    
    return ( double )getSample( channel, sample, iNewNoiseTrace );
}
