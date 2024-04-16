//! VNoiseFileReader read noise traces from a file

#ifndef VNoiseFileReader_H
#define VNoiseFileReader_H

#include "VGrIsuReader.h"

#include <iostream>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VNoiseFileReader
{
    protected:

        bool fDebug;
        bool fZombie;

        unsigned int fNoiseFileType;              // 0 = grisu style
        string       fNoiseFileName;

        VGrIsuReader* fGrIsuReader;

        // placeholders
        std::valarray<double> v;
        std::vector< std::valarray<double> > vv;
        std::vector< vector< uint8_t > > v8;
        std::vector< uint8_t > vv8;

    public:

        VNoiseFileReader( unsigned int iType = 0, string iFileName = "" );
        ~VNoiseFileReader() {}

        valarray<double>&          getPeds();
        valarray<double>&          getPedvars();
        vector< valarray<double> >& getPedvarsAllSumWindows();
        valarray<double>&          getPedRMS();
        vector< vector< vector< uint8_t > > > getFullNoiseVec();
        vector< vector< uint8_t > >& getFullNoiseVec( unsigned int iTel );
        vector< uint8_t >&         getFullNoiseVec( unsigned int iTel, int iChannel );
        uint16_t                   getNumSamples()
        {
            if( fGrIsuReader )
            {
                return fGrIsuReader->getNumSamples();
            }
            else
            {
                return 0;
            }
        }

        bool init( VDetectorGeometry* iD, unsigned int intel, vector<int> iSW, bool iDebug = false, int iseed = 0, double iFADCorrect = 1. );
        bool isZombie()
        {
            return fZombie;
        }
        uint8_t            getNoiseSample( unsigned int iTel, uint32_t iHitID, unsigned int iSample, bool iNewTrace = true );
        vector< uint8_t >& getNoiseVec( unsigned int iTel, uint32_t iHitID, bool iNewTrace = true );
        void  setDefaultGrisuPed( double iB );
        void  setSumWindow( unsigned int iTelID, int isw );
};
#endif
