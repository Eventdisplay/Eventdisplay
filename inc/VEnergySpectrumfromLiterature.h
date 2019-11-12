//! VEnergySpectrumfromLiterature provide spectra from literature

#ifndef VEnergySpectrumfromLiterature_H
#define VEnergySpectrumfromLiterature_H

#include "VPlotUtilities.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TObject.h"
#include "TRandom.h"
#include "TSystem.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct sData
{
    string           Name;
    string           Observatory;
    string           Source;
    string           Reference;
    string           Comment;
    unsigned int     Type;
    vector< double > Parameter;
    vector< double > ParError;
    double           EnergyRange_min;             // TeV
    double           EnergyRange_max;             // TeV
    vector< double > FluxV_energy;                // TeV
    vector< double > FluxV_DiffFlux;              // cm^-2 s^-1 TeV^-1
    vector< double > FluxV_DiffFluxError;         // cm^-2 s^-1 TeV^-1
};

struct sEnergyFun
{
    string           Name;
    string           Description;
    unsigned int     NumParameters;
};

class VEnergySpectrumfromLiterature : public VPlotUtilities
{
    private:
    
        bool  bIsZombie;
        
        bool  fPlottingLogEnergyAxis;
        float fPlottingMinEnergy;
        float fPlottingMaxEnergy;
        float fPlottingYaxisMin;
        float fPlottingYaxisMax;
        float fPlottingMultiplierIndex;
        
        vector< sEnergyFun > fEnergyFun;
        
        vector< sData > fData;
        
        vector< TH1D* > fRandomErrorHistograms;
        
        bool checkIDRange( unsigned int iID );
        
        unsigned int fIntegral_ID;
        TF1*         fIntegral_TF1;
        double       fIntegral_x[1000];
        double       fIntegral_y[1000];
        
        bool         prepare_integration( unsigned int iID = 0, double iEmin = 1.e-3, double iEmax = 1.e3 );
        void         setFunctions();
        
    public:
    
        VEnergySpectrumfromLiterature( string ifile = "", bool iprint = true );
        ~VEnergySpectrumfromLiterature() {}
        
        sData  getEnergySpectrumDataField( unsigned int iID = 0 );
        TF1*   getEnergySpectrum( unsigned int iID = 0, bool bLogEnergy = true, double iEnergyMin_Lin = -99., double iEnergyMax_Lin = -99. );
        TGraphAsymmErrors* getEnergySpectrumWithErrors( unsigned int iID = 0, bool bLogEnergy = true );
        TGraphErrors* getDifferentialFluxPoints( unsigned int iID = 0, bool bLogEnergy = true );
        double getIntegralFlux( double iEmin = 1., double iEmax = 1.e3,  unsigned int iID = 0 );
        double getPowerLaw_Index( unsigned int iID = 0 );
        double getPowerLaw_FluxConstant_at1TeV( unsigned int iID = 0 );
        TH1D*  getRandomErrorHistograms( unsigned int i )
        {
            if( i < fRandomErrorHistograms.size() )
            {
                return fRandomErrorHistograms[i];
            }
            else
            {
                return 0;
            }
        }
        bool   isValidID( unsigned int iID )
        {
            return checkIDRange( iID );
        }
        bool   isZombie()
        {
            return bIsZombie;
        }
        
        void listValues();
        void listValues( unsigned int iID );
        TCanvas* plot( unsigned int iID = 0, TCanvas* c = 0, bool iPlotY = true );
        TCanvas* plot( string iselection, TCanvas* c = 0 );
        bool readValuesFromFile( string ifile = "$OBS_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues.dat", bool iPrint = true );
        
        void setPlottingLogEnergyAxis( bool iB = true )
        {
            fPlottingLogEnergyAxis = iB;
        }
        void setPlottingMultiplierIndex( float iS = 0. )
        {
            fPlottingMultiplierIndex = iS;
        }
        void setPlottingEnergyRangeLinear( float xmin = 0.03, float xmax = 100. )
        {
            fPlottingMinEnergy = xmin;
            fPlottingMaxEnergy = xmax;
        }
        void setPlottingEnergyRangeLog( float xmin = -1.5, float xmax = 2. )
        {
            fPlottingMinEnergy = TMath::Power( 10., xmin );
            fPlottingMaxEnergy = TMath::Power( 10., xmax );
        }
        void setPlottingYaxis( float iMin = 1.e-17, float iMax = 1.e-7 )
        {
            fPlottingYaxisMin = iMin;
            fPlottingYaxisMax = iMax;
        }
        
        ClassDef( VEnergySpectrumfromLiterature, 5 );
};
#endif
