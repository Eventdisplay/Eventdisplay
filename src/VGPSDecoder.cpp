#include "VGPSDecoder.h"

void VGPSDecoder::decode( uint32_t word0, uint32_t word1, uint32_t word2 )
throw()
{

    // hopefully code below is endian independant!
    
    uint16_t TimeArray[5];
    TimeArray[0] = word0 & 0x0000FFFF;
    TimeArray[1] = ( word0 & 0xFFFF0000 ) >> 16;
    TimeArray[2] = word1 & 0x0000FFFF;
    TimeArray[3] = ( word1 & 0xFFFF0000 ) >> 16;
    TimeArray[4] = word2 & 0x0000FFFF;
    
    //   uint32_t timeVal[3];
    //   timeVal[0]=word0;
    //   timeVal[1]=word1;
    //   timeVal[2]=word2;
    //   uint16_t *TimeArray = (uint16_t*)timeVal;
    
    // what follows as copied from Scott:
    
    fGPSStatus = ( TimeArray[1] >> 4 & 0xF );
    
    fGPSDays = 100 * ( TimeArray[0] & 0x000F )
               + 10 * ( TimeArray[1] >> 12 & 0x000F )
               + ( TimeArray[1] >> 8 & 0x000F );
               
    fGPSHrs  = 10 * ( TimeArray[1] >> 4 & 0x000F )   + ( TimeArray[1] >> 0 & 0x000F );
    fGPSMins = 10 * ( TimeArray[2] >> 12 & 0x000F )  + ( TimeArray[2] >> 8 & 0x000F );
    fGPSSecs = 10 * ( TimeArray[2] >> 4 & 0x000F )   + ( TimeArray[2] >> 0 & 0x000F );
    
    fGPSSecs += 1E-1 * ( TimeArray[3] >> 12 & 0x000F ) +
                1E-2 * ( TimeArray[3] >> 8 & 0x000F )  +
                1E-3 * ( TimeArray[3] >> 4 & 0x000F )  +
                1E-4 * ( TimeArray[3] >> 0 & 0x000F )  +
                1E-5 * ( TimeArray[4] >> 12 & 0x000F ) +
                1E-6 * ( TimeArray[4] >> 8 & 0x000F )  +
                1E-7 * ( TimeArray[4] >> 4 & 0x000F );
                
    return;
}


void VGPSDecoder::printDataToStream( std::ostream& stream ) throw()
{
    stream << "Status:          " << ( int ) getStatus() << std::endl;
    stream << "Days:            " << ( int ) getDays() << std::endl;
    stream << "Hours:           " << ( int ) getHrs() << std::endl;
    stream << "Minutes:         " << ( int ) getMins() << std::endl;
    stream << "Seconds:         " << ( double ) getSecs() << std::endl;
}
