//! VStar data class for star catalogue

#ifndef VStar_H
#define VStar_H

#include "TMath.h"
#include "TObject.h"

#include "VAstronometry.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VStar : public TObject
{
    public:
    
        unsigned int fStarID;
        string fStarName;
        double fDec2000;
        double fRA2000;
        double fDecCurrentEpoch;
        double fRACurrentEpoch;
        double fRunGalLong1958;
        double fRunGalLat1958;
        double fBrightness_V;
        double fBrightness_B;
        double fMajorDiameter;                        // this is either the source diameter or the positional error
        double fMinorDiameter;
        double fPositionAngle;
        double fMajorDiameter_68;
        double fMinorDiameter_68;
        double fPositionAngle_68;
        double fSignificance;
        double fSpectralIndex;
        double fSpectralIndexError;
        string fSpectrumType;
        double fCutOff_MeV;
        double fCutOffError_MeV;
        vector< double > fFluxEnergyMin;
        vector< double > fFluxEnergyMax;
        vector< double > fFlux;
        vector< double > fFluxError;
        bool   fVariability;
        vector< string > fOtherNames;
        string fType;
        vector< string > fAssociations;
        int fQualityFlag;
        
        VStar();
        ~VStar() {};
        
        double getDistance( VStar* iStar1, VStar* iStar2 = 0 );
        void   printStar();
        
        ClassDef( VStar, 4 );
};


#endif
