//! VInstrumentResponseFunctionData data class for response functions (e.g. angular resolution)

#ifndef VInstrumentResponseFunctionData_H
#define VInstrumentResponseFunctionData_H

#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"

#include <iostream>
#include <string>
#include <vector>

#include "CData.h"
#include "VHistogramUtilities.h"
#include "VPlotUtilities.h"
#include "VSpectralWeight.h"

using namespace std;

class VInstrumentResponseFunctionData : public TObject, public VHistogramUtilities
{
    private:
    
        // data tree
        CData*   fData;               //!
        
        // energy binning
        int     fHistogrambinningEnergy_TeV_Log;
        double  fHistogrambinningEnergy_Min_Tev_Log;
        double  fHistogrambinningEnergy_Max_Tev_Log;

        // angular binning
        int     fHistogrambinningAngular_Log;
	double  fHistogrambinningAngular_Min_Log;
	double  fHistogrambinningAngular_Max_Log;
        
        // array centre
        double  fArrayCentre_X;
        double  fArrayCentre_Y;
        
        TList*   calculateResolution( TH2D* iHistogram, TGraphErrors* iResult, string iHistoName, 
                                      double iContainmentProbability );
        double   getResolutionErrorfromToyMC( double i68, double iN );
        int      testResponseFunctionType( string iType );
        
    public:
    
        // list of function types
        vector< string > fListofResponseFunctionTypes;
        
        // basic data
        string  fName;
        string  fType;                // descripes type of response function (e.g. angular resolution or core resolution)
        int     fType_numeric;        //    (same as integer)
        
        // characteristics (all angles in [deg])
        double  fZe;
        int     fAz_bin;
        double  fAz_min;
        double  fAz_max;
        double  fXoff;
        double  fYoff;
        double  fWobble;
        int     fNoise;
        double  fPedvars;
        double  fSpectralIndex;
        
        double  fMCMaxCoreRadius;
        unsigned int fNTel;
        
        unsigned int fEnergyReconstructionMethod;
        
        // list of histogram types
        enum    E_HISTOID { E_DIFF, E_DIFF2, E_LOGDIFF, E_NIMAG, E_DIST, E_ERROR, E_RELA,
                                    E_DIFF_MC, E_DIFF2_MC, E_LOGDIFF_MC };
        TList*                     fHistogramList;
        vector< TH2D* >            f2DHisto;
        vector< TGraphErrors* >    fResolutionGraph;
        vector< double >           fContainmentProbability;
        
        VInstrumentResponseFunctionData();
        ~VInstrumentResponseFunctionData();
        void   fill( double iWeight );
        TList* getListofHistograms()
        {
            return fHistogramList;
        }
        bool   initialize( string iName, string iType, unsigned int iNTel, double iMCMaxCoreRadius = 500. );
        void   setArrayCentre( double iX = 0., double iY = 0. )
        {
            fArrayCentre_X = iX;
            fArrayCentre_Y = iY;
        }
        void   setDataTree( CData* iData )
        {
            fData = iData;
        }
        void   setData( double iZe = -99., int iAz_bin = -99, double iAz_min = -99., double iAz_max = -99.,
                        int iNoise = -99, double iPedvars = -99., double iIndex = -99.,
                        double iXoff = -99., double iYoff = -99. );
        void   setEnergyReconstructionMethod( unsigned int iMethod = 0 )
        {
            fEnergyReconstructionMethod = iMethod;
        }
	void   setHistogramEbinning( int iN = 60, double iMin = -2.0, double iMax = 4.0 )
        {
            fHistogrambinningEnergy_TeV_Log = iN;
            fHistogrambinningEnergy_Min_Tev_Log = iMin;
            fHistogrambinningEnergy_Max_Tev_Log = iMax;
        }
        void   setHistogramLogAngbinning( int iN = 20, double iMin = -4.0, double iMax = 1.0 )
	{
            fHistogrambinningAngular_Log = iN;
            fHistogrambinningAngular_Min_Log = iMin;
            fHistogrambinningAngular_Max_Log = iMax;
        }
        bool   terminate( double iContainmentProbability, double iContainmentProbabilityError );
        
        ClassDef( VInstrumentResponseFunctionData, 9 );
};

#endif
