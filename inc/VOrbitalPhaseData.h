//! VOrbitalPhase data element / calculation of orbital phases

#ifndef VORBITALPHASE_H
#define VORBITALPHASE_H

#include "VFluxAndLightCurveUtilities.h"

#include "TMath.h"
#include <TNamed.h>

#include <iostream>
#include <string>

using namespace std;

class VOrbitalPhaseData
{
    public:
    
        string fName;
        double fZeroPhase_MJD;
        double fOrbit_days;
        double fOrbit_days_error_low;
        double fOrbit_days_error_high;
        
        VOrbitalPhaseData();
        virtual ~VOrbitalPhaseData() {}
        
        double getOrbitalPhase( double iMJD );
        void print();
        bool isSet();
        
        ClassDef( VOrbitalPhaseData, 1 );
};

#endif
