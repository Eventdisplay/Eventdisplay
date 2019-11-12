//! VStereoHistograms holds all stereo histograms

#ifndef VSTEREOHISTOGRAMS_H
#define VSTEREOHISTOGRAMS_H

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TList.h"
#include "TTree.h"

#include <iostream>
#include <map>
#include <string>

#include "VUtilities.h"

using namespace std;

class VStereoHistograms
{
    private:
    
        bool bIsOn;
        string fHisSuffix;
        double fBinSize;
        double fBinSizeUC;
        double fBinSizeEnergy;
        double fBinSizeTime;
        double fSkyMapSizeXmin;
        double fSkyMapSizeXmax;
        double fSkyMapSizeYmin;
        double fSkyMapSizeYmax;
        double fTimeMin;
        double fTimeMax;
        
        int    fRunNumber;
        
        bool readHistograms( TList*, string );
        
    public:
    
        TList* hisList;                           //!< list with all histograms
        TList* hListStereoParameterHistograms;    //! list with stereo parameters
        //! list with random forest parameters
        TList* hListRandomForestParameterHistograms;
        TList* hListParameterHistograms;          //!< list with parameter histograms
        map< string, TH1* > hListNameofParameterHistograms;
        TList* hListEnergyHistograms;             //!< list with energy histograms
        TList* hListSkyMaps;                      //!< list with sky maps
        vector< string > hListNameofSkyMaps;      //!< list with histogram names of sky maps
        TList* hListSkyMapsUC;                    //!< list with sky maps (uncorrelated bins)
        
        // data quality monitoring
        TH1D* hTriggerPatternBeforeCuts;
        TH1D* hTriggerPatternAfterCuts;
        TH1D* hImagePatternBeforeCuts;
        TH1D* hImagePatternAfterCuts;
        
        // parameter histograms
        TH1D* htheta2;                            //!< Theta2 Histogram
        TH1D* hmean_width;                        //!< Mean Width Histogram
        TH1D* hmean_length;                       //!< Mean Length Histogram
        TH1D* hmean_dist;                         //!< Mean Distance Histogram
        TH2D* hcore;                              //!< Shower Core map (ground plane)
        TH1D* hmscw;                              //!< MSCW histogram
        TH1D* hmscl;                              //!< MSCL histogram
        TH2D* hmsc;                               //!< MSCW vs MSCL histogram
        TH1D* hZetaTau;                           //!< Theta2 plots toward Zeta Tau
        TH1D* hemiss;                             //!< mean emission height
        TH1D* hemissC2;                           //!< mean emission height Chi2
        TH1D* herecChi2;                          //!< chi2 from energy reconstruction
        
        // ratio of signal to background area (from energy dependent theta2 cut)
        TH1D* hmap_MeanSignalBackgroundAreaRatio; //!< signal to background area ratio
        TH1D* hmap_MeanSignalBackgroundAreaRatioUC; //!< signal to background area ratio
        
        // random forest histograms
        TH1D* hrf;                                //!< random forest classifier
        
        // energy histograms (logarithmic energy axis)
        TH1D* herecCounts;                        //!< reconstructed energy
        TH1D* hDuration1DtimeBinned;              //!< duration of the time bin taking time mask into account
        TH1D* hRealDuration1DtimeBinned;          //!< duration of the time bin taking time mask and dead time fraction into account
        TH2D* herecCounts2DtimeBinned;            //!< reconstructed energy vs observing time (2D)
        TH2D* herecCounts2D_vs_distance;              //!< reconstructed energy vs distance to camera centre
        //time-dependent differential energy spectrum
        TH2D* herecWeights;                       //!< weights vs.  reconstructed energy
        // energy histograms (linear energy axis)
        TH1D* hLinerecCounts;                     //!< reconstructed energy
        TH2D* hLinerecCounts2DtimeBinned;
        TH2D* hLinerecWeights;                    //!< weights vs.  reconstructed energy
        
        // sky maps (uncorrelated)
        TH2D* hmap_stereoUC;                      //!< Sky map (correlated bins)
        TH2D* hmap_alphaUC;                       //!< Background normalisation map (correlated bins)
        TH2D* hmap_alphaNormUC;                   //!< Background normalisation map  (correlated bins)
        TH2D* h_combine_map_alpha_offUC;                   //!< Background normalisation map, off map for on run (correlated bins)
        TH2D* h_combine_map_stereo_onUC;
        TH2D* h_combine_map_stereo_offUC;
        
        // sky maps (correlated)
        TH2D* hmap_stereo;                        //!< Sky map (correlated bins)
        TH2D* hmap_alpha;                         //!< Background normalisation map (correlated bins)
        TH2D* hmap_alphaNorm;                     //!< Background normalisation map  (correlated bins)
        TH2D* h_combine_map_alpha_off;
        TH2D* h_combine_map_stereo_on;
        TH2D* h_combine_map_stereo_off;
        
        // rate lists
        TList* hisRateList;
        TH1D* hrate_1sec;                         //!< Event Rate Histogram (1 second bins)
        TH1D* hrate_10sec;                        //!< Event Rate Histogram (10 second bins)
        TH1D* hrate_1min;                         //!< Event Rate Histogram (1 minute bins)
        
        VStereoHistograms( string i_hsuffix, double ibinsize, double ibinsizeUC, double iEnergyBinSize,
                           double iTimeBinSize, double iTimeMin, double iTimeMax,  bool ion );
        ~VStereoHistograms();
        void defineHistograms();
        void deleteParameterHistograms();
        void deleteSkyPlots();
        int  getRunNumber()
        {
            return fRunNumber;
        }
        void makeRateHistograms( double, double );
        bool readParameterHistograms();
        bool readSkyPlots();
        void scaleDistributions( double );
        void setRunNumber( int iRun )
        {
            fRunNumber = iRun;
        }
        void setSkyMapSize( double xmin, double xmax, double ymin, double ymax );
        void writeObjects( string, string, TObject* );
        void writeHistograms();
};
#endif
