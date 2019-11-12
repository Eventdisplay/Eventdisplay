//! VStereoMaps   model the background

#ifndef VStereoMaps_H
#define VStereoMaps_H

#include "CData.h"

#include "VAnaSumRunParameter.h"
#include "VExclusionRegions.h"
#include "VRadialAcceptance.h"

#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

#include <iostream>

using namespace std;

struct sRE_REGIONS
{
    int noff;                                     //!< number of source regions
    vector< double > xoff;                        //!< x-position of center of off source region
    vector< double > yoff;                        //!< y-position of center of off source region
    vector< double > roff;                        //!< radius of off source region
};

class VStereoMaps
{
    private:
    
        VAnaSumRunParameterDataClass fRunList;
        CData* fData;
        
        // theta2 cut (might be energy dependent)
        double fTheta2Cut_Max;
        
        double fTargetShiftWest;
        double fTargetShiftNorth;
        
        bool   fTMPL_RE_nMaxoffsource;
        
        // regions excluded from sky maps
        vector< VListOfExclusionRegions* > fListOfExclusionRegions;
        // radial acceptances
        VRadialAcceptance* fAcceptance;
        
        bool bUncorrelatedSkyMaps;
        bool fNoSkyPlots;                         //!< do full sky analysis (if false, analyse source region only)
        
        TH2D* hmap_stereo;
        TH2D* hmap_alpha;
        TH1D* hmap_ratio;
        
        TRandom3* fRandom;
        
        int fInitRun;
        
        void makeTwoDStereo_BoxSmooth( double, double, double, double, double );
        
        // theta2 calculation
        unsigned int fTheta2_length;
        vector< double > fTheta2;
        vector< double > fTheta2_weight;
        vector< double > fTheta2_weightREonly;
        vector< double > fTheta2_All;
        
        void initialize_theta2();
        
        // RING BACKGROUND MODEL
        TFile* fRM_file;
        
        bool fill_RingBackgroundModel( double, double, int, bool );
        bool initialize_RingBackgroundModel( bool iIsOn );
        void RM_calculate_norm();
        void RM_getAlpha( bool );
        
        // REFLECTED REGION MODEL:
        vector< vector< sRE_REGIONS > > fRE_off;  //!< off region parameters
        double fRE_roffTemp;                      //!< radius of off source region
        
        bool fill_ReflectedRegionModel( double, double, int, bool );
        bool fill_ReflectedRegionModel( double, double, int, bool, double& i_theta2 );
        void RE_getAlpha( bool iIsOn );
        bool initialize_ReflectedRegionModel();
        void initialize_ReflectedRegionHistograms();
        
        // histograms related to reflected region model
        TH2D* hRE_NRegions;
        TTree* hRE_regions;
        
        bool initialize_Histograms();
        TList* hAuxHisList;                       //!< histograms needed for various calculations
        
        // some variables needed for efficient filling
        double f_RE_binXW;
        double f_RE_binYW;
        int f_RE_xstart;
        int f_RE_xstopp;
        int f_RE_ystart;
        int f_RE_ystopp;
        double f_RE_AreaNorm;
        int f_RE_WW;
        int f_RE_WN;
        
        // etc
        void   cleanup();                         // delete all objects not needed anymore
        bool   defineAcceptance();
        //!< return if event is in on region
        bool   fillOn( double x_sky, double y_sky, int irun, bool ishapecuts, double& i_theta2 );
        //!< return if event is in off region
        bool   fillOff( double x_sky, double y_sky, int irun, bool ishapecuts, double& i_theta2 );
        double phiInt( double );
        
    public:
        TH1D* hAux_theta2On;                      //
        TH1D* hAux_theta2Off;
        TH1D* hAux_theta2Ratio;
        
        VStereoMaps( bool, int, bool );
        ~VStereoMaps();
        
        void              calculateTheta2( bool, double, double );
        bool              fill( bool is_on, double x_sky, double y_sky, double theta2CutMax,
                                int irun, bool ishapecuts, double& i_theta2 );
        void              finalize( bool iIsOn , double OnOff_Alpha = 1.0 );
        VRadialAcceptance*      getAcceptance()
        {
            return fAcceptance;
        }
        unsigned int      getTheta2_length()
        {
            return fTheta2_length;
        }
        vector< double >& getTheta2()
        {
            return fTheta2;
        }
        vector< double >& getTheta2_weigth()
        {
            return fTheta2_weight;
        }
        vector< double >& getTheta2_weigthREonly()
        {
            return fTheta2_weightREonly;
        }
        vector< double >& getTheta2_All()
        {
            return fTheta2_All;
        }
        TList*            getAux_hisList()
        {
            return hAuxHisList;
        }
        void              setData( CData* c )
        {
            fData = c;
        }
        void              setHistograms( TH2D*, TH2D*, TH1D* );
        void              setNoSkyPlots( bool iS )
        {
            fNoSkyPlots = iS;
        }
        void              setRunList( VAnaSumRunParameterDataClass iL );
        void              setTargetShift( double iW, double iN );
        void              setRegionToExclude( vector< VListOfExclusionRegions* > iF );
};
#endif
