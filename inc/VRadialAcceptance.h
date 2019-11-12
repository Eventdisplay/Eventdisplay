//!< VRadialAcceptance radial acceptance for a given point on the sky

#ifndef VACCEPTANCE_H
#define VACCEPTANCE_H

#include <cmath>
#include <iostream>
#include <string>

#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TList.h"

#include "CData.h"
#include "VAnaSumRunParameter.h"
#include "VExclusionRegions.h"
#include "VGammaHadronCuts.h"
#include "VUtilities.h"

using namespace std;

// fit function for radial acceptance
Double_t VRadialAcceptance_fit_acceptance_function( Double_t* x, Double_t* par );

// class definition
class VRadialAcceptance
{
    private:
    
        // if true: only a subset of histograms/functions is written to disk
        bool fproduction_shortIO;
        
        bool fAcceptanceFunctionDefined;
        
        // position of gamma-ray source and size
        double fSourcePosition_X;
        double fSourcePosition_Y;
        double fSourcePosition_Radius;
        double fMaxDistanceAllowed;
        double fCut_CameraFiducialSize_max;
        
        // energy reconstruction method
        unsigned int fEnergyReconstructionMethod;
        
        // regions excluded from background
        vector< VListOfExclusionRegions* > fListOfExclusionRegions;
        
        // calculate acceptance from data
        VGammaHadronCuts*    fCuts;
        VAnaSumRunParameter* fRunPar;
        
        // histogram lists
        TList* hList;
        TList* histListProductionIO;
        TList* hListNormalizeHistograms;
        TList* hListFitHistograms;
        
        // 1D radial acceptance histograms
        TH1F* fRadialAcceptance;                      // average acceptance curve
        TF1*  fRadialAcceptanceFit;                   // fit to average acceptance curve
        vector< TH1F* > fRadialAcceptance_perRun;     // run dependent acceptance curves
        vector< TF1* >  fRadialAcceptanceFit_perRun;  // fit to run dependent acceptance curves
        
        // histograms and variables for normalisation of 1D radial acceptance curves
        double fAzCut_min;
        double fAzCut_max;
        TH1F* hscale;
        vector< TH1F* >  hscaleRun;                //!< run dependent scaling histograms
        vector< TH1F* >  hscaleRunRatio;           //!< run dependent ratio scaling histograms
        vector< TH2F* >  hAreaExcluded2D;          //!< run dependent 2D histograms used for scaling
        
        // debug histograms
        TH1F* hPhiDist;
        TH1F* hPhiDistDeRot;
        
        vector< double > fPhiMin;                   //!< Phi bins (limits)
        vector< double > fPhiMax;                   //!< Phi bins (limits)
        vector< TH1F* >  hAccPhi;                   //!< azimuth angle dependent acceptance curves (azimuth angle in camera coordinates)
        vector< TH1F* >  hAccPhiDerot;              //!< azimuth angle dependent acceptance curves (azimuth angle in derotated camera coordinates)
        unsigned int     fAccZeFitMinBin;          //!< range (in bins) for normalisation of acceptance curves
        unsigned int     fAccZeFitMaxBin;
        
        // acceptance vs yoff vs xoff
        TH2F* hXYAccTotDeRot;
        TH1F* hXYAccTotDeRotPhiDependentSlice ;
        TH1F* hXYAccTotDeRotRadiusDependentSlice000 ;
        vector< TH2F* > hXYAccRun;
        
        // number of raw files used to calculate acceptance
        double fNumberOfRawFiles;
        
        // 2D acceptance stuff
        int    f2DAcceptanceMode ; // USE2DACCEPTANCE
        double f2DBinNormalizationConstant ; // USE2DACCEPTANCE
        
        // for PhiDependentSlices - 1D histogram of a constant-radius band, bin positions depend on phi
        double phi_minphi ; // minimum phi range
        double phi_maxphi ; // maximum phi range
        int    phi_nbins ; // number of phi bins
        double phi_minradius ; // min radius of phi slices
        double phi_maxradius ; // max radius of phi slices
        
        // for RadiusDependentSlices - 1D histogram of constant-phi band, bin positions depend on radius
        double rad_minrad ; // minimum radius
        double rad_maxrad ; // maximum radius
        int    rad_nbins  ; // number of radius bins
        double rad_phiwidth ; // radians +- N, S, E, or W
        // if rad_phiwidth = pi/8, then RadiusDependentSliceN will go from phi=15pi/8 to phi=pi/8
        
        double eventphi ; // phi of event
        double eventradius ; // radius of event
        int eventcount ;
        
        string fExtraHistogramDir ;
        int fExtraHistogramMode ;
        
        vector <TH2F*> hXYAccImgSel ;
        vector <TH2F*> hXYAccImgSelPreDeRot ;
        vector <TH1F*> hXYAccImgSelRadiusDependentSlice000 ; // pie slice centered on   0 deg +- rad_phiwidth
        vector <TH1F*> hXYAccImgSelPhiDependentSlice ; // donut band from phi_minradius to phi_maxradius
        
        vector <TH2F*> hXYAccNImages ;
        vector <TH2F*> hXYAccNImagesPreDeRot ;
        vector <TH1F*> hXYAccNImagesRadiusDependentSlice000 ; // pie slice centered on   0 deg +- rad_phiwidth
        vector <TH1F*> hXYAccNImagesPhiDependentSlice ; // donut band from phi_minradius to phi_maxradius
        
        // get acceptance curves from a file
        TFile* fAccFile;
        
        // reset all variables
        void reset();
        
    public:
    
        // use acceptance curve from simulation
        VRadialAcceptance();
        // set data source and cuts for acceptance curve calculation
        VRadialAcceptance( VGammaHadronCuts* iCuts, VAnaSumRunParameter* irun, double iMaxDistanceAllowed = -99. );
        // use acceptance curve from this file
        VRadialAcceptance( string ifile, int irun = -1 );
        ~VRadialAcceptance();
        
        int    calculateAverageRadialAcceptanceCurveFromRuns( TDirectory* iDirectory );
        // correct run-wise radial acceptances for exclusion regions
        bool   correctRadialAcceptancesForExclusionRegions( TDirectory* iDirectory, unsigned int iRunNumber );
        int    fillAcceptanceFromData( CData* c, int entry, double x_rotJ2000, double y_rotJ2000 );
        double getAcceptance( double x, double y );   //!< return radial acceptance
        double getNumberofRawFiles()
        {
            return fNumberOfRawFiles;
        }
        bool   isExcluded( double, double );                                            //!< region excluded from analysis
        bool   isExcludedfromBackground( double, double );                              //!< region excluded from background analysis
        bool   isExcludedfromSource( double, double );                                  //!< region excluded from source analyis
        void   setAzCut( double iAzMin = -1.e9, double iAzMax = 1.e9 )
        {
            fAzCut_min = iAzMin;    //!< cut on Az (shower directory)
            fAzCut_max = iAzMax;
        }
        void   setEnergyReconstructionMethod( unsigned int iEMethod = 1 )
        {
            fEnergyReconstructionMethod = iEMethod;
        }
        void   setProductionIO( bool  production_shortIO = false )
        {
            fproduction_shortIO = production_shortIO;
        }
        // set source position, radius, and minimal distance between source and background
        void   setSource( double x, double y, double r, double imaxdist = 5. );
        void   setRegionToExcludeAcceptance( vector< VListOfExclusionRegions* > iF );
        bool   terminate( TDirectory* iDirectory );
        
        // for simple 2D acceptance
        // acceptance = closest derotated bin content / NormalizationConstant
        // calling this function updates the value of f2DBinNormalizationConstant
        double calculate2DBinNormalizationConstant( double radius = 0.3 ) ; // radius in degrees // USE2DACCEPTANCE
        
        // write two files <basename>.dat and <basename>.meta
        // the .dat file contains each bin index, bin center, and bin content
        // the .meta file contains the number of bins, and the column names (headers) of the dat file
        void Write2DHistToTextFile( TH2F* hist, string& basename ) ;  // USE2DACCEPTANCE
        
        // write 1D slice histogram to text file
        void Write1DHistToTextFile( TH1F* hist, string& basename ) ;  // USE2DACCEPTANCE
        void Write1DHistToTextFile( TH1F* hist, string& basename, int histtype ) ;  // USE2DACCEPTANCE
        
        // will write many histograms to directory 'dirname'
        // each hist will be saved as two text files (not root files!)
        //
        //void WriteHistsToDirectory( string &dirname ) ;
        
        // sets the mode for the class to use
        // mode=0 : normal operation, radial acceptance used
        // mode=1 : use normalized 2d histogram to compute acceptance
        // doesn't depend on zenith angle
        // see getAcceptance() for more info
        int Set2DAcceptanceMode( int mode = 0 ) ; // USE2DACCEPTANCEndif
        
        // will write several extra histograms to text files in directory 'histdir'
        void SetExtraHistogramDirectory( string histdir ) ;
        // if ehm > 0, then will make a bunch of extra histograms, and
        // write them to text files in a directory set by
        int  SetExtraHistogramMode( int ehm ) ;
        
        
};

#endif
