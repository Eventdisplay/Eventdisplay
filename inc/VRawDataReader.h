//! VRawDataReader steering class for reading raw data format

#ifndef VRAWDATAREADER_H
#define VRAWDATAREADER_H

#include "VBaseRawDataReader.h"
#include <VRawDataFileRead.h>
#include <VRawEventParser.h>
#include <VRawDataExceptions.h>
#include <VVirtualDataReader.h>
#include <VSimulationDataReader.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VRawDataReader : public VBaseRawDataReader
{

    private:
        //   vector< vector< bool > > fFullTrigVec;  // MS: This might be here by mistake..?
        
    protected:
        vector< VRawDataFileRead* > fFileRead;
        vector< VRawEventData* > fEvent;
        
    public:
        VRawDataReader( string,
                        int isourcetype,
                        unsigned int iNTel,
                        bool iDebug );
        virtual ~VRawDataReader();
        
        bool                        getNextEvent();
};
#endif
