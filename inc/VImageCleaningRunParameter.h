//! VImageCleaningRunParameter run parameters for image cleaning

#ifndef VImageCleaningRunParameter_H
#define VImageCleaningRunParameter_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VImageCleaningRunParameter
{

    public:
    
        string       fName;
        unsigned int fTelID;
        
        unsigned int fImageCleaningMethod;   // 0: standard two level cleaning; 1: time cluster cleaning, 2: Maxim..., 3: trace correlation method, 4: time two-level
        
        /////////////////////////////////////////////////
        // standard two-level image/border cleaning
        double fimagethresh;              // parameter for image threshold
        double fborderthresh;             // parameter for border threshold
        double fbrightnonimagetresh;      // parameter for bright pixels threshold
        
        bool fUseFixedThresholds;         // use fixed image/border thresholds instead of multiples of pedestal variances
        
        /////////////////////////////////////////////////
        // time cluster cleaning
        double ftimecutpixel;             //  parameter for time cut between pixels
        double ftimecutcluster;           //  parameter for time cut between clusters
        int    fminpixelcluster;          //  parameter for minimum number of pixels in cluster
        int    floops;                    //  parameter for number of loops for border pixel finding
        
        /////////////////////////////////////////////////
        // trace correlation cleaning
        double fCorrelationCleanBoardThresh;  //  parameter for lower border threshold
        double fCorrelationCleanCorrelThresh; //  parameter for trace correlation level (0.6-1.0)
        int    fCorrelationCleanNpixThresh;   //  Maximum Number of pixels to apply correlation cleaning to (eg 10-15)
        
        /////////////////////////////////////////////////
        // cluster cleaning
        int fnmaxcluster;
        double fminsizecluster;
        
        /////////////////////////////////////////////////
        // time two-level
        double ftimediff;                 // parameter for time constraint between next neighbor pixels
        
        /////////////////////////////////////////////////
        // optimized next-neigbour cleaning
        // (note that some parameters in NNcleaninginputcard)
        double fNNOpt_FakeImageProb;
        vector< bool >  fNNOpt_ActiveNN;
        vector< string > fNNOpt_Multiplicities;
        unsigned int fNNOpt_nRings;
        float fNNOpt_CoincWinLimit;
        bool fNNOpt_ifNNoptNoTimeing;
        bool fNNOpt_ifExplicitSampleTimeSlice;
        float fNNOpt_sampleTimeSlice;
        unsigned int fNNOpt_nBinsADC;
        
        VImageCleaningRunParameter( string iName = "ImageCleaningParameter" );
        ~VImageCleaningRunParameter() {};
        
        string getImageCleaningMethod();
        unsigned int getImageCleaningMethodIndex()
        {
            return fImageCleaningMethod;
        }
        string getName()
        {
            return fName;
        }
        bool initialize();
        void print();
        bool setImageCleaningMethod( string iMethod );
        void setTelID( unsigned int iTelID )
        {
            fTelID = iTelID;
        }
        
};

#endif
