//! VGammaHadronCutsStatistics keep track of efficiency of different cuts

#ifndef VGammaHadronCutsStatistics_H
#define VGammaHadronCutsStatistics_H

#include <TTree.h>

#include <bitset>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VGammaHadronCutsStatistics : public TNamed
{
    private:
    
        TTree* fData;              //!
        
        bitset< 32 >           fCut_bitset;  //!
        unsigned long int      fCut_bitset_ulong;  //!
        vector< string >       fCutName;  //!
        vector< unsigned int > fCutCounter;  //!
        
    public:
    
        // enum for efficiency counting of the different types of cuts
        enum EN_AnaCutsStats { eTot, eMC_XYoff, eXYoff, eStereoQuality, eArrayChi2, eNImages, eMSC_Quality,
                               eErec, eCorePos, eLTrig, eSizeSecondMax, eTelType, eDirection, eIsGamma, eEnergyRec,
                               eError
                             };
                             
                             
        VGammaHadronCutsStatistics();
        ~VGammaHadronCutsStatistics() {};
        
        void         fill();
        unsigned int getCounterValue( unsigned int iCut );
        TTree*       getDataTree()
        {
            return fData;
        }
        void         initialize( string iName );
        void         printCutStatistics();
        void         reset();
        void         setCutCounter( unsigned int iCut, unsigned int iValue );
        void         terminate();
        void         updateCutCounter( unsigned int iCut );
        
        ClassDef( VGammaHadronCutsStatistics, 3 );
};


#endif
