// VSpectralWeight calculate energy dependent spectral weight

#ifndef VSpectralWeight_H
#define VSpectralWeight_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TMath.h"
#include "TObject.h"

using namespace std;

class VSpectralWeight : public TObject
{
    private:
    
        bool fDebug;
        
        double fMCSpectralIndex;
        double fMCMinEnergy_TeV_Lin;
        double fMCMaxEnergy_TeV_Lin;
        double fMCMaxConeAngle_deg;
        double fMCSimulatedEvents;
        double fMCFluxConstant;
        TF1*   fEnergySpectrum;
        
        double fIndex;
        double fSpectralWeightAlpha;
        
        void   calculateMCFluxConstant();
        
        
    public:
    
        VSpectralWeight();
        ~VSpectralWeight() {}
        
        double getSpectralWeight( double iEnergy_TeV_Lin );
        void   print();
        void   setEnergySpectrum( TF1* f )
        {
            fEnergySpectrum = f;
        }
        void   setMCParameter( double iMCSpectralIndex = 2., double iMCEnergy_min_TeV_Lin = 0.03,
                               double iMCEnergy_max_TeV_Lin = 200., double iMCMaxConeAngle_deg = 0., double iMCSimulatedEvents = 0. );
        void   setSpectralIndex( double iG, bool iPrint = false );
        
        ClassDef( VSpectralWeight, 2 );
};

#endif
