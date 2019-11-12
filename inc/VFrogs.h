#ifndef VFROGS_H_INC
#define VFROGS_H_INC

#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"

#include "frogs.h"

#include "VEvndispData.h"
#include "VFrogsParameters.h"
#include <VDetectorGeometry.h>
#include <VGrIsuAnalyzer.h>

#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <string>
#include <vector>
#include <valarray>
#include <sstream>

using namespace std;
#define VFROGSNEPOCH 10

class VFrogs : public VEvndispData, public VGrIsuAnalyzer
{
    public:
        VFrogs();
        ~VFrogs();
        
        // vectors for readTableFrogs mscw runNumber and Erec
        vector<int>    fTableEventNumber;
        vector<double> fTableEnergy;
        
        vector<int>    fAnasumRunNumber;
        vector<int>    fAnasumEventNumber;
        
        void  doFrogsStuff( int, string );                      //!< do the actual analysis (called for each event)
        int   getFrogsEventID()
        {
            return frogsEventID;
        }
        int   getFrogsGSLConStat()
        {
            return frogsGSLConStat;
        }
        int   getFrogsNB_iter()
        {
            return frogsNB_iter;
        }
        int   getFrogsNImages()
        {
            return frogsNImages;
        }
        ULong64_t   getFrogsSelectedImages()
        {
            return frogsSelectedImages;
        }
        float getFrogsXS()
        {
            return frogsXS;
        }
        float getFrogsXSerr()
        {
            return frogsXSerr;
        }
        float getFrogsYS()
        {
            return frogsYS;
        }
        float getFrogsYSerr()
        {
            return frogsYSerr;
        }
        float getFrogsXP()
        {
            return frogsXP;
        }
        float getFrogsXPerr()
        {
            return frogsXPerr;
        }
        float getFrogsYP()
        {
            return frogsYP;
        }
        float getFrogsYPerr()
        {
            return frogsYPerr;
        }
        float getFrogsXPGC()
        {
            return frogsXPGC;
        }
        float getFrogsYPGC()
        {
            return frogsYPGC;
        }
        float getFrogsEnergy()
        {
            return frogsEnergy;
        }
        float getFrogsEnergyerr()
        {
            return frogsEnergyerr;
        }
        float getFrogsLambda()
        {
            return frogsLambda;
        }
        float getFrogsLambdaerr()
        {
            return frogsLambdaerr;
        }
        float getFrogsGoodnessImg()
        {
            return frogsGoodnessImg;
        }
        int   getFrogsNpixImg()
        {
            return frogsNpixImg;
        }
        float getFrogsGoodnessBkg()
        {
            return frogsGoodnessBkg;
        }
        int   getFrogsNpixBkg()
        {
            return frogsNpixBkg;
        }
        float getFrogsXPStart()
        {
            return frogsXPStart;
        }
        float getFrogsYPStart()
        {
            return frogsYPStart;
        }
        float getFrogsXPED()
        {
            return frogsXPED;
        }
        float getFrogsYPED()
        {
            return frogsYPED;
        }
        float getFrogsXSStart()
        {
            return frogsXSStart;
        }
        float getFrogsYSStart()
        {
            return frogsYSStart;
        }
        float getFrogsTelGoodnessImg( int i )
        {
            return frogsTelGoodnessImg[i];
        }
        float getFrogsTelGoodnessBkg( int i )
        {
            return frogsTelGoodnessBkg[i];
        }
        
        void initAnalysis();
        void initFrogsTree();
        void initOutput();
        void reset();
        void terminate();
        float transformTelescopePosition( int iTel, float i_ze, float i_az, int axis );
        float transformShowerPosition( float i_ze, float i_az, float xcore, float ycore, float zcore, int axis );
        float transformPosition( float i_ze, float i_az, float x, float y, float z, int axis, bool bInv );
        void readTableFrogs();
        double getFrogsStartEnergy( int eventNumber );
        int getFrogsAnasumNumber( int eventNumber, int runNumber );
        void finishFrogs( TFile* f );
        void transformResults();
        TFile* mscwFrogsFile;
        TFile* AnasumFrogsFile;
        
        
    private:
        struct 		     frogs_imgtmplt_in frogs_convert_from_ed( int eventNumber, double inEnergy, string fArrayEpoch );
        VEvndispData*         fData;                    //!< pointer to data class
        
        int frogsRecID;
        string templatelistname;
        string fparamfile;
        
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
        
        bool  fInitialized;                        //!< true after initialization
        int   fStartEnergyLoop;
        
        double frogsTemplateMu0[499];
        double frogsTemplateMu1[499];
        double frogsTemplateMu2[499];
        double frogsTemplateMu3[499];
        
        float frogsTelGoodnessImg[4];
        float frogsTelGoodnessBkg[4];
        
        //const int nepoch = 10 ; // limit on how many epochs we should consider
        double frogsLowerThresh[VFROGSNEPOCH]; // must be same as nepoch!
        double frogsFirstParam[ VFROGSNEPOCH];
        double frogsSecondParam[VFROGSNEPOCH];
        double frogsDCtoPE[VFROGSNEPOCH];
        double frogsPMTNoise[VFROGSNEPOCH]; //pmt electronic noise
        bool   frogsMinimization;
        double frogsDeltaXS;
        double frogsDeltaYS;
        double frogsDeltaXP;
        double frogsDeltaYP;
        double frogsDeltaLog10e;
        double frogsDeltaLambda;
        int    frogsInterpOrder;
        bool   frogsCheating;
        int    frogsNBEventCalib;
        
        Float_t         frogsZe;
        Float_t         frogsAz;
        Float_t		frogsXS_derot;
        Float_t		frogsYS_derot;
        Float_t		frogsR[VDST_MAXTELESCOPES];
        
        unsigned long ffrogsRandomSeed;            	 // random seed for frogs differential evolution.
        
        void processParamFile() ;
};
#endif
