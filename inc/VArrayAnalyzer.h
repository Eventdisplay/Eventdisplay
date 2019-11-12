//! VArrayAnalyzer class for analyzing VERITAS data (full array analysis)
#ifndef VARRAYANALYZER_H
#define VARRAYANALYZER_H

#include "TMath.h"

#include "VEvndispData.h"
#include "VDispAnalyzer.h"
#include "VGrIsuAnalyzer.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VMLPAnalyzer.h"
#include "VShowerParameters.h"
#include "VSkyCoordinatesUtilities.h"
#include "VStarCatalogue.h"

#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

class VArrayAnalyzer : public VEvndispData, public VGrIsuAnalyzer
{
    private:
    
        bool fDebug;
        bool fInitialized;                        //!< true after initialization
        
        vector< VDispAnalyzer* > fDispAnalyzer;
        
        vector< double > fMeanPointingMismatch;   //!< mean pointing mismatch between eventdisplay and vbf (per telescope)
        vector< double > fNMeanPointingMismatch;
        
        // temporary variables needed for array reconstruction
        vector< unsigned int > telID;
        vector< float > x;
        vector< float > y;
        vector< float > w;
        vector< float > l;
        vector< float > m;
        vector< float > phi;
        vector< float > sinphi;
        vector< float > cosphi;
        vector< float > width;
        vector< float > length;
        vector< float > asym;
        vector< float > loss;
        vector< float > dist;
        vector< float > pedvar;
        vector< float > tgrad;
        vector< float > xcore;
        vector< float > ycore;
        vector< float > cen_y;
        vector< float > cen_x;
        vector< float > Xoff;
        vector< float > Yoff;
        vector< float > ltrig;
        vector< float > az;
        vector< float > ze;
        vector<ULong64_t> teltype;
        
        vector< float > xtelnew;
        vector< float > ytelnew;
        vector< float > ztelnew;
        
        // private functions
        
        void calcShowerDirection_and_Core();      //!< calculate shower core and direction
        void checkPointing();                     //!< check for mismatching between different pointing values
        void initializeDispAnalyzer( unsigned int iStereoMethodID );
        void prepareforCoreReconstruction( unsigned int iMeth, float xs, float ys );
        void prepareforDirectionReconstruction( unsigned int iMethIndex, unsigned int iReconstructionMethod );
        bool fillSimulationEvent( bool iArrayTrigger );
        bool fillShowerDirection( unsigned int iMeth, float xoff, float yoff, float stds );
        bool fillShowerCore( unsigned int iMeth, float ximp, float yimp ); //!< fill shower core results into VEvndispData
        double getMeanPointingMismatch( unsigned int iTel );
        string getTMVAFileNameForAngularReconstruction( unsigned int iStereoMethodID, string iBDTFileName = "BDTDisp_BDT_" );
        void initEvent();                         //!< reset vectors, etc. (called for each event)
        int  rcs_method_0( unsigned int );        //!< GrIsu reconstruction method 1(!)
        int  rcs_method_3( unsigned int );
        int  rcs_method_4( unsigned int );
        int  rcs_method_5( unsigned int, unsigned int );
        int  rcs_method_7( unsigned int );
        int  rcs_method_8( unsigned int );
        float recalculateImagePhi( double, double );
        void selectShowerImages( unsigned int );  //!< select shower images to be used in determinate of shower coordinates
        bool testMinAngleBetweenImages( unsigned int iMethod, float iangdiff );
        //!< transform telescope positions into shower coordinates
        void transformTelescopePosition( int iTel, float ize, float iaz, bool i_MC );
        bool updatePointingToStarCatalogue( unsigned int iTelescopeID );
        
    public:
    
        VArrayAnalyzer();
        ~VArrayAnalyzer();
        void doAnalysis();
        void calcTelescopePointing();
        void generateReducedPointingTreeData();
        void updatePointingToArbitraryTime( int iMJD, double iTime ) ;
        void initAnalysis();
        void initOutput();
        void terminate( bool iDebug_IO = false );
        void initTree();
};
#endif
