//! VInterpolate2DHistos interpolate empty bins in 2D histograms

#ifndef VINTERPOLATE2DHISTOS_H
#define VINTERPOLATE2DHISTOS_H

#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <iostream>
#include <string>

using namespace std;

class VInterpolate2DHistos : public TNamed
{
    private:
        TRandom3* fRandom;
        
    public:
    
        VInterpolate2DHistos( int iseed = 0 );
        ~VInterpolate2DHistos() {}
        
        TH2F* doSimpleInterpolation( TH2F*, string, int, int, bool, TH2F* hNevents = 0, int iMinEvents = 0 );
        TH2F* doGaussianInterpolation( TH2F* h, string iname, TH2F* hNevents = 0, int nGausN = 1, double nWidth = 1. );
        TH2F* doLogLinearExtrapolation( TH2F* h, string iname, TH2F* hNevents = 0, int iMinEvents = 20 );
        
        ClassDef( VInterpolate2DHistos, 1 );
};
#endif
