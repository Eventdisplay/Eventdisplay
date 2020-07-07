//! VImageAnalyzerData  data class with pixel and image analysis results

#ifndef VImageAnalyzerData_H
#define VImageAnalyzerData_H

#include "TList.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <VImageAnalyzerHistograms.h>
#include <VSpecialChannel.h>
#include <VImageParameter.h>

#include <valarray>
#include <vector>

using namespace std;

class VImageAnalyzerData
{
    private:
    
        unsigned int fTelID;
        unsigned int fNChannels;
        unsigned int fMaxChannels;
        unsigned int fNSamples;
        
        unsigned int fTraceIntegrationMethod;
        
    public:
    
        bool fFillMeanTraces;
        bool fFillPulseSum;
        
        // image analysis data
        VImageParameter* fImageParameter;             //!< image parameter from geometrical analysis
        VImageParameter* fImageParameterLogL;         //!< image parameter from log likelihood analysis
        VImageAnalyzerHistograms* fAnaHistos;               //!< analysis test histograms
        
        valarray<double> fTemplateMu;
        valarray<double> fSums;
        valarray<double> fSums2;
        vector<bool> fLLEst;                      //!< true for dead channels with estimated sum from loglikelihood fit
        unsigned int fNDead;
        vector<unsigned int> fDead;               //!< dead channel (bit coded)
        vector<unsigned int> fMasked;             //! masked channels
        vector<bool> fDeadRecovered;              //!< dead channel recovered for example by smoothing
        vector<unsigned int> fDeadUI;             //!< dead channel (largest bit)
        unsigned int fLowGainNDead;
        vector<unsigned int> fLowGainDead;        //!< dead channel (bit coded); low gain channels
        vector<unsigned int> fLowGainDeadUI;      //!< dead channel (largest bit); low gain channels
        vector<bool> fLowGainDeadRecovered;       //!< dead channel recovered for example by smoothing (low gain)
        
        int fRandomMakeDeadChannelsSeed;
        TRandom3* fRandomMakeDeadChannels;        //!< random generator for setting randomly channels dead
        vector< bool > fHiLo;                     //!< true if low gain is set
        vector< unsigned short int > fZeroSuppressed;  //!< 0=not suppressed; bit 1 = charge only, bit 2 = sample only
        unsigned int fpulsetiming_tzero_index;    //!< position of first 50% value in pulse timing level vector
        unsigned int fpulsetiming_width_index;    //!< position of second 50% value in pulse timing level vector
		unsigned int fpulsetiming_triggertime_index;    //!< trigger time 
        vector< valarray<double> > fPulseTimingUncorrected; //!< pulse timing at certain fraction of pulse maxima (uncorrected values)
        vector< valarray<double> > fPulseTimingCorrected; //!< pulse timing at certain fraction of pulse maxima (corrected values)
        valarray< double >     fPulseTimingAverageTime; //!< average pulse time
        valarray< double >     fPulseTimingAverageTimeCorrected; //!< average pulse time (corrected values)
        valarray<unsigned int> fTCorrectedSumFirst;
        valarray<unsigned int> fTCorrectedSumLast;
        valarray<unsigned int> fTCorrectedSum2First;
        valarray<unsigned int> fTCorrectedSum2Last;
        valarray<unsigned int> fCurrentSummationWindow;
        valarray<unsigned int> fCurrentSummationWindow_2;
        valarray<double> fTraceWidth;             //!< trace width
        valarray<double> fTraceMax;               //!< trace maximum
        valarray<unsigned int> fTraceN255;
        valarray<double> fRawTraceMax;            //!< raw trace maximum
        // cleaning
        vector<bool> fImage;
        vector<int> fImageUser;                   //!< channels which are enabled/disables by the user
        vector<bool> fBorder;
        vector<bool> fTrigger;                    //!< MS: pixels selected by trigger algorithm
        vector<bool> fBrightNonImage;             //!< pixel above non image threshold
        vector<bool> fImageBorderNeighbour;       //!< image and border pixel plus their neighbours
        vector<bool> fBorderBorderNeighbour;      //!< image and border pixel plus their neighbours and neighbours of neighbours
        // time cleaning
        vector<int> fClusterNpix;                 //!< HP: numer of pixels in cluster
        vector<unsigned int> fClusterID;                   //!< HP: cluster ID
        int fMainClusterID;                       //!< HP: main cluster ID
        vector<double> fClusterSize;              //!< HP: size of the cluster
        vector<double> fClusterTime;              //!< HP: weighted mean time of the cluster
        vector<double> fClusterCenx;              //!< HP: cenX of the cluster
        vector<double> fClusterCeny;              //!< HP: cenY of the cluster
        int fncluster_cleaned;                    //!< HP: number of clusters
        int fncluster_uncleaned;                  //!< HP: number of clusters before cluster rejection
        // correlation cleaning
        vector<double> fCorrelationCoefficient;  // correlation coefficient
        // time since run start
        double fTimeSinceRunStart;                //!< time since run start
        double fTimeRunStart;                     //!< time of first event in run
        
        // special channels
        VSpecialChannel* fSpecialChannel;
        
        // FADCstop info
        vector< double > fFADCstopTZero;
        vector< double > fFADCstopSum;
        
        // mean pulse histograms
        TList* hMeanPulses;
        vector< TProfile2D* > hMeanPulseHigh;     //!< high gain mean pulse
        vector< TProfile2D* > hMeanPulseLow;      //!< low gainpulse  mean pulse
        TList* hPulseSum;
        vector< TH1F* > hPulseSumHigh;            //!< integrated charge for high gain channels
        vector< TH1F* > hPulseSumLow;             //!< integrated charge for low gain channels
        
        // dummy vector
        vector< unsigned int > iDummyVectorUI;
        
        VImageAnalyzerData( unsigned int iTelID, unsigned int iShortTree = 0, bool bCalibration = false );
        ~VImageAnalyzerData() {}
        
        void                     fillPulseSum( unsigned int, double, bool );
        TList*                   getMeanPulseHistograms()
        {
            return hMeanPulses;
        }
        TList*                   getIntegratedChargeHistograms()
        {
            return hPulseSum;
        }
        vector< unsigned int >&  getFADCstopTrigChannelID();
        double                   getHIGHQE_gainfactor( unsigned int iChannel );
        int                      getRandomDeadChannel()
        {
            return fRandomMakeDeadChannels->Integer( fDead.size() );
        }
        int                      getRandomDeadChannelSeed()
        {
            return fRandomMakeDeadChannelsSeed;
        }
        valarray<double>&        getTraceAverageTime( bool iCorrected );
		valarray<double>&        getTTrigger();
        valarray<double>&        getTZeros( bool iCorrected );
        valarray<double>&        getTraceWidth( bool iCorrected );
        VSpecialChannel*         getSpecialChannel()
        {
            return fSpecialChannel;
        }
        void                     initialize( unsigned int iChannels, unsigned int iMaxChannels, bool iDebug, int iseed,
						 unsigned int iSamples, unsigned int ipulsetiminglevel, unsigned int iTzeroIndex, unsigned int iWidthIndex,
                                                 unsigned int iTTindex );
        void                     initializeMeanPulseHistograms();
        void                     initializeIntegratedChargeHistograms();
        bool                     readSpecialChannels( int iRunNumber, string iFile, string iDirectory );
        void                     setTraceIntegrationMethod( unsigned iN = 1 )
        {
            fTraceIntegrationMethod = iN;
        }
        void                     setTrace( unsigned int iChannel, vector< double > fT, bool iHiLo, double fPeds );
};
#endif
