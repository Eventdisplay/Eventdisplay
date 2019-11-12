//! VDispAnalyzer wrapper class for all modified disp analysis

#ifndef VDispAnalyzer_H
#define VDispAnalyzer_H

#include "TFile.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

#include "VDispTableAnalyzer.h"
#include "VGrIsuAnalyzer.h"
#include "VMLPAnalyzer.h"
#include "VStatistics.h"
#include "VTMVADispAnalyzer.h"

using namespace std;

class VDispAnalyzer
{
    private:
    
        bool                fDebug;
        
        bool                bZombie;
        
        string              fDispMethod;
        
        VMLPAnalyzer*       fMLPAnalyzer;
        VDispTableAnalyzer* fDispTableAnalyzer;
        VTMVADispAnalyzer*  fTMVADispAnalyzer;
        
        float fAxesAngles_min;
        unsigned int fNImages_min;
        float fdistance_max;
        float floss_max;
        bool  fDispErrorWeighting;
        float fDispErrorExponential;
        
        // disp direction reconstruction
        float f_disp;
        float f_dispE;
        float f_dispDiff;
        float f_xs;
        float f_ys;
        float f_angdiff;
        vector< float > fdisp_xs_T;
        vector< float > fdisp_ys_T;
        vector< float > fdisp_xy_weight_T;
        
        // disp direction error
        vector< float > fdisp_error_T;
        
        // disp energy reconstruction
        float fdisp_energy;
        float fdisp_energy_chi;
        float fdisp_energy_dEs;
        float fdisp_energy_median;
        float fdisp_energy_medianAbsoluteError;
        vector< float > fdisp_energy_T;
        unsigned int   fdisp_energy_NT;
        
        // disp core reconstruction
        vector< float > fdisp_core_T;
        
        vector<ULong64_t> fTelescopeTypeList;
        
        void calculateMeanShowerDirection( vector< float > v_x, vector< float > v_y, vector< float > v_weight,
                                           float& xs, float& ys, float& dispdiff, unsigned int iMaxN );
                                           
    public:
    
        VDispAnalyzer();
        ~VDispAnalyzer() {}
        
        void calculateCore( unsigned int i_ntel, float iArrayElevation, float iArrayAzimuth,
                            double* itelX, double* itelY, double* itelZ,
                            ULong64_t* iTelType,
                            double* img_size, double* img_cen_x, double* img_cen_y,
                            double* img_cosphi, double* img_sinphi,
                            double* img_width, double* img_length, double* img_asym,
                            double* img_tgrad, double* img_loss, int* img_ntubes,
                            double* img_weight,
                            double xoff_4, double yoff_4,
                            double* iR,
                            double xcore, double ycore,
                            double xs, double ys,
                            double* img_fui );
                            
        void calculateEnergies( unsigned int i_ntel, float iArrayElevation, float iArrayAzimuth,
                                ULong64_t* iTelType,
                                double* img_size, double* img_cen_x, double* img_cen_y,
                                double* img_cosphi, double* img_sinphi,
                                double* img_width, double* img_length, double* img_asym,
                                double* img_tgrad, double* img_loss, int* img_ntubes,
                                double* img_weight,
                                double xoff_4, double yoff_4,
                                double* iR, double iEHeight,
                                double iMCEnergy = -1.,
                                double* img_fui = 0 );
                                
        void  calculateMeanDirection( float& xs, float& ys,
                                      vector< float > x, vector< float > y,
                                      vector< float > cosphi, vector< float > sinphi,
                                      vector< float > v_disp, vector< float > v_weight,
                                      float& dispdiff,
                                      float x_off4 = -999., float yoff_4 = -999. );
                                      
                                      
        void calculateMeanDirection( unsigned int i_ntel, float iArrayElevation, float iArrayAzimuth,
                                     ULong64_t* iTelType,
                                     double* img_size, double* img_cen_x, double* img_cen_y,
                                     double* img_cosphi, double* img_sinphi,
                                     double* img_width, double* img_length, double* img_asym,
                                     double* img_tgrad, double* img_loss, int* img_ntubes,
                                     double* img_weight,
                                     double xoff_4, double yoff_4,
                                     vector< float > dispErrorT,
                                     double* img_fui );
                                     
        void calculateExpectedDirectionError( unsigned int i_ntel, float iArrayElevation, float iArrayAzimuth,
                                              ULong64_t* iTelType,
                                              double* img_size, double* img_cen_x, double* img_cen_y,
                                              double* img_cosphi, double* img_sinphi,
                                              double* img_width, double* img_length, double* img_asym,
                                              double* img_tgrad, double* img_loss, int* img_ntubes,
                                              double* img_weight,
                                              double xoff_4, double yoff_4,
                                              double* img_fui );
                                              
        float evaluate( float iWidth, float iLength, float iAsymm, float iDist,
                        float iSize, float iPedvar, float itgrad, float iLoss,
                        float icen_x, float icen_y, float xoff_4, float yoff_4, ULong64_t iTelType,
                        float iZe, float iAz, float iRcore = -99., float iFui = -1., bool b2D = true );
        float getAngDiff()
        {
            return f_angdiff;
        }
        float getDisp()
        {
            return f_disp;
        }
        float getDispDiff()
        {
            return f_dispDiff;
        }
        float getDispE()
        {
            return f_dispE;
        }
        float getDispErrorT( unsigned int iTelescopeNumber );
        float getCoreDistance( unsigned int iTelescopeNumber );
        float getEnergy();
        float getEnergyChi2();
        float getEnergydES();
        float getEnergyMedian();
        float getEnergyMedianAbsoluteError();
        float getEnergyT( unsigned int iTelescopeNumber );
        unsigned int   getEnergyNT()
        {
            return fdisp_energy_NT;
        }
        
        float getXcoordinate_disp()
        {
            return f_xs;
        }
        float getXcoordinate_disp( unsigned int i, float x = -999., float cosphi = -999. );
        float getYcoordinate_disp()
        {
            return f_ys;
        }
        float getYcoordinate_disp( unsigned int i, float y = -999., float sinphi = -999. );
        vector< float >& getXYWeight_disp()
        {
            return fdisp_xy_weight_T;
        }
        float getXYWeight_disp( unsigned int i )
        {
            if( i < fdisp_xy_weight_T.size() )
            {
                return fdisp_xy_weight_T[i];
            }
            return -999.;
        }
        bool  initialize( string iFile, string iDispMethod, string iDispType = "BDTDisp" );
        bool  isZombie()
        {
            return bZombie;
        }
        void  setDebug( bool iFDebug = false )
        {
            fDebug = iFDebug;
        }
        void  setDispErrorWeighting( bool iW = false, float iWeight = 5. )
        {
            fDispErrorWeighting = iW;
            fDispErrorExponential = iWeight;
        }
        void  setQualityCuts( unsigned int iNImages_min = 0, float iAxesAngles_min = 0., float imaxdist = 1.e5, float imaxloss = 1. )
        {
            fAxesAngles_min = iAxesAngles_min;
            fNImages_min    = iNImages_min;
            fdistance_max   = imaxdist;
            floss_max       = imaxloss;
        }
        void  setTelescopeTypeList( vector<ULong64_t> iTelescopeTypeList );
        void  setZombie( bool iB = true )
        {
            bZombie = iB;
        }
        void  terminate();
};

#endif
