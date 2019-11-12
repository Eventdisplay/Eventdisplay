//! VMLPAnalyzer MLP based disp method

#ifndef VMLPAnalyzer_H
#define VMLPAnalyzer_H

#include "TFile.h"
#include "TMultiLayerPerceptron.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace std;

class VMLPAnalyzer
{
    private:
    
        bool bZombie;
        
        TFile* fFile;
        TMultiLayerPerceptron* fMLP;
        
        Double_t fmlp_var[10];                    //!< input variables for the MLP
        
    public:
    
        VMLPAnalyzer( string iFile = "" );
        ~VMLPAnalyzer() {}
        float evaluate( float& width, float& length, float& asym, float& size, float& dist );
        bool isZombie()
        {
            return bZombie;
        }
        void terminate();
};
#endif
