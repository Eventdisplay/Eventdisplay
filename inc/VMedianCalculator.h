// VMedianCalculator calculate approximate median of a set

#ifndef VMedianCalculator_H
#define VMedianCalculator_H

#include "TMath.h"

#include <iostream>
#include <vector>

using namespace std;

class VMedianCalculator
{
    private:
        int nDim_exact;          // exact median until this value
        
        vector< double > x;
        int n_counter;
        
        float eta;
        float mean_x;
        float mean_xx;
        float quantiles[3];                        // 0.16, 0.5, 0.84 quantiles
        float prob[3];                             // probabilities (hardwired 0.16, 0.5, 0.84)
        
    public:
        VMedianCalculator();
        ~VMedianCalculator() {}
        
        void   fill( double x );
        double getMean();
        double getMedian( float& medianwidth, int& N );
        double getMedianWidth();
        int    getN()
        {
            return n_counter;
        }
        double getRMS();
        void   reset();
        void   setEta( double iEta = 0.01 )
        {
            eta = iEta;
        }
        void   setNExact( int n = 1000 )
        {
            nDim_exact = n;
        }
};

#endif
