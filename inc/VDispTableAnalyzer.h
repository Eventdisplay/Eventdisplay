//! VDispTableAnalyzer get DISP values from tables

#ifndef VDispTableAnalyzer_H
#define VDispTableAnalyzer_H

#include "TFile.h"
#include "TMath.h"

#include <bitset>
#include <iostream>
#include <map>
#include <vector>

#include "VDispTableReader.h"

using namespace std;


struct sDispTableEventData
{
    float x;
    float y;
    float cosphi;
    float sinphi;
    float disp;
};

class VDispTableAnalyzer
{
    private:
    
        bool bDebug;
        bool bZombie;
        
        TFile* fFile;
        VDispTableReader* fData;
        
        float f_disp;
        float f_dispE;
        float f_disp_Phi;
        float f_disp_PhiE;
        float f_disp_Miss;
        
        vector< float > x_disp;          // x coordinate of disp calculation
        vector< float > y_disp;          // y coordinate of disp calculation
        
        float calculateDisp( float iWidth, float iLength, float iSize, float iPedvar, float iZe, float iAz, bool b2D = true );
        double interpolate( double w1, double ze1, double w2, double ze2, double ze, bool iCos = false );
        
    public:
    
        VDispTableAnalyzer( string iFile = "" );
        ~VDispTableAnalyzer();
        void  calculateMeanDirection( float& xs, float& ys, vector< float > x, vector< float > y,
                                      vector< float > cosphi, vector< float > sinphi, vector< float > v_disp, vector< float > v_weight );
        float evaluate( float iWidth, float iLength, float iSize, float iPedvar, float iZe, float iAz, bool b2D = true );
        float getDisp()
        {
            return f_disp;
        }
        float getDispE()
        {
            return f_dispE;
        }
        float getDispPhi()
        {
            return f_disp_Phi;
        }
        float getDispPhiE()
        {
            return f_disp_PhiE;
        }
        float getDispMiss()
        {
            return f_disp_Miss;
        }
        float getXcoordinate_disp( unsigned int i )
        {
            if( i < x_disp.size() )
            {
                return x_disp[i];
            }
            else
            {
                return -9999.;
            }
        }
        float getYcoordinate_disp( unsigned int i )
        {
            if( i < y_disp.size() )
            {
                return y_disp[i];
            }
            else
            {
                return -9999.;
            }
        }
        bool  isZombie()
        {
            return bZombie;
        }
        void  terminate();
};
#endif
