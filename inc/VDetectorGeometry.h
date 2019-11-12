//! VDetectorGeometry  geometry of all cameras and telescope arrangement (note: most parameters are in VCameraRead)

#ifndef VDETECTORGEOMETRY_H
#define VDETECTORGEOMETRY_H

#include <iostream>
#include <vector>

#include "VCameraRead.h"

using namespace std;

class VDetectorGeometry : public VCameraRead
{
    private:
        vector< unsigned int > fNSamples;
        vector< unsigned int > fNChannels;
        vector< bool >         fSampleWarning;
        
    public:
        VDetectorGeometry() {}
        VDetectorGeometry( unsigned int iNTel, bool iDebug = false );
        VDetectorGeometry( unsigned int iNTel, vector< string > iCamera, string iDir, bool iDebug = false,
                           float iCoordinateTransformerX = 1., float iCoordinateTransformerY = 1., int iSourceType = 3 );
        ~VDetectorGeometry() {}
        vector< unsigned int > getNChannels()
        {
            return fNChannels;
        }
        unsigned int   getNChannels( unsigned int iTelID )
        {
            if( iTelID < fNChannels.size() )
            {
                return fNChannels[iTelID];
            }
            else
            {
                return 0;
            }
        }
        vector< unsigned int > getNSamples()
        {
            return fNSamples;
        }
        void           addDataVector( unsigned int iNTel, vector< unsigned int > iNChannels );
        unsigned int   getNSamples( unsigned int iTelID )
        {
            if( iTelID < fNSamples.size() )
            {
                return fNSamples[iTelID];
            }
            else
            {
                return 0;
            }
        }
        void           setNChannels( unsigned int iTelID, unsigned int iNChannels );
        void           setNSamples( unsigned int iTelID, unsigned int iNSamples, bool iForceSet = false );
};
#endif
