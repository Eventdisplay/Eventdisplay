//! VFrogsParameters storage class for Frogs data
#ifndef VFROGSPARAMETERS_H
#define VFROGSPARAMETERS_H

#include "VGlobalRunParameter.h"
#include "TTree.h"

#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using namespace std;

class VFrogsParameters
{
    private:
        bool fDebug;
        TTree* fTreeFrogs;                        //!< output tree
        
    public:
    
        VFrogsParameters();
        ~VFrogsParameters();
        
        void fill()
        {
            if( fTreeFrogs )
            {
                fTreeFrogs->Fill();
            }
        }
        TTree* getTree()
        {
            return fTreeFrogs;
        }
        void initTree( string, string );
        void printParameters();                   //!< write tree parameters to standard output
        void reset();                             //!< reset all tree variable to standard values
        
        unsigned int fNTel;                       //!< number of telescopes
        
        int   frogsEventID;
        int   frogsGSLConStat;
        int   frogsNB_iter;
        int   frogsNImages;
        ULong64_t frogsSelectedImages;
        float frogsXS;
        float frogsXSerr;
        float frogsYS;
        float frogsYSerr;
        float frogsXP;
        float frogsXPerr;
        float frogsYP;
        float frogsYPerr;
        float frogsXPGC;
        float frogsYPGC;
        float frogsEnergy;
        float frogsEnergyerr;
        float frogsLambda;
        float frogsLambdaerr;
        float frogsGoodnessImg;
        int   frogsNpixImg;
        float frogsGoodnessBkg;
        int   frogsNpixBkg;
        
        float frogsXPStart;
        float frogsYPStart;
        float frogsXPED;
        float frogsYPED;
        float frogsXSStart;
        float frogsYSStart;
        
        float frogsTelGoodnessImg0;
        float frogsTelGoodnessImg1;
        float frogsTelGoodnessImg2;
        float frogsTelGoodnessImg3;
        float frogsTelGoodnessBkg0;
        float frogsTelGoodnessBkg1;
        float frogsTelGoodnessBkg2;
        float frogsTelGoodnessBkg3;
        
        Float_t         frogsZe;
        Float_t         frogsAz;
        Float_t		frogsXS_derot;
        Float_t		frogsYS_derot;
        double 		frogsR[VDST_MAXTELESCOPES];
        
};
#endif
