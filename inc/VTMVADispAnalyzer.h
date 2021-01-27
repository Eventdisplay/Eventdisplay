// VTMVADispAnalyzer TMVA based modified disp analysis

#ifndef VTMVADispAnalyzer_H
#define VTMVADispAnalyzer_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TMath.h"

#include "TMVA/Config.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

using namespace std;

class VTMVADispAnalyzer
{
    private:
    
        bool fDebug;
        bool bZombie;
        string fDispType;
        
        vector<ULong64_t> fTelescopeTypeList;
        map< ULong64_t, TMVA::Reader* > fTMVAReader;
        
        float fWidth;
        float fLength;
        float fWoL;
        float fSize;
        float fNtubes;
        float fPedvar;
        float fTGrad;
        float fZe;
        float fAz;
        float fLoss;
        float fDist;
        float fFui;
        float fAsymm;
        float fXcore;
        float fYcore;
        float fcross;
        float fRcore;
        float fEHeight;
        
    public:
    
        VTMVADispAnalyzer( string iFile, vector< ULong64_t > iTelTypeList, string iDispType = "BDTDisp" );
        ~VTMVADispAnalyzer() {}
        
        float evaluate( float iWidth, float iLength, float iSize, float iAsymm, float iLoss,
                        float iTGrad, float icen_x, float icen_y, float xoff_4, float yoff_4,
                        ULong64_t iTelType, float iZe, float iAz, float iRcore,
                        float iEHeight = -1., float iDist = -1., float iFui = -1., float iNtubes = -1 );
        bool isZombie()
        {
            return bZombie;
        }
        void terminate();
        
};

#endif
