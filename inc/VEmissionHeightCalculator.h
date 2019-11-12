//! VEmissionHeightCalculator calculate emission height (and possibly get systematic error in energy reconstruction)

#ifndef VEmissionHeightCalculator_H
#define VEmissionHeightCalculator_H

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "VUtilities.h"

using namespace std;

class VEmissionHeightCalculator
{
    private:
    
        bool fDebug;
        
        unsigned int fNTelPairs;
        double fEmissionHeight;
        double fEmissionHeightChi2;
        vector< double > fEmissionHeightT;
        
        unsigned int fNTel;
        vector< double > fTelX;
        vector< double > fTelY;
        vector< double > fTelZ;
        
        double getTelescopeDistanceSC( unsigned int iTel1, unsigned int iTel2, double az, double z );
        double imageDistance( double c1x, double c2x, double c1y, double c2y );
        
    public:
    
        VEmissionHeightCalculator();
        ~VEmissionHeightCalculator();
        double getEmissionHeight( double* cen_x, double* cen_y, double* size, double az, double el );
        double getMeanEmissionHeight()
        {
            return fEmissionHeight;
        }
        double getMeanEmissionHeightChi2()
        {
            return fEmissionHeightChi2;
        }
        vector< double >& getEmissionHeights()
        {
            return fEmissionHeightT;
        }
        unsigned int getNTelPairs()
        {
            return fNTelPairs;
        }
        void setTelescopePositions( vector< float > x, vector< float > y, vector< float > z );
        void setTelescopePositions( unsigned int ntel, double* x, double* y, double* z );
};
#endif
