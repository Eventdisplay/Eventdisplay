//! VCTARequirements plot CTA requirements on top of sensitivity and angular resolution plots

#ifndef VCTARequirements_H
#define VCTARequirements_H

#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"

#include "VCTASensitivityRequirements.h"
#include "VCTASensitivityRequirementsUpdated2017.h"
#include "VPlotUtilities.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class VCTARequirements : public VPlotUtilities
{
    private:
    
        int fSetOfRequirementID;
        string fRequirementTitle;
        
        string fRequirementsDirectory;
        
        string fPlotSubSystemRequirement;
        double fRequirementsGraphLineWidth;
        
        TGraph* fReqDifferentialSensitivity; // vs log10 energy [TeV]
        TGraph* fReqEffectiveArea;           // vs log10 energy [TeV]
        TGraph* fReqAngularResolution;       // vs log10 energy [TeV]
        TGraph* fReqEnergyResolution;        // vs log10 energy [TeV]
        
        TGraph* fReqDifferentialSensitivity2014;  // original requirement from 2014
        TGraph* fReqDifferentialSensitivityxF;
        TGraph* fReqDifferentialSensitivityx2;
        TGraph* fReqDifferentialSensitivityx3;
        TGraph* fReqDifferentialSensitivityx4;
        TGraph* fReqDifferentialSensitivityLST;
        TGraph* fReqDifferentialSensitivityMST;
        TGraph* fReqDifferentialSensitivitySST;
        TGraph* fReqAngularResolutionLST;
        TGraph* fReqAngularResolutionMST;
        TGraph* fReqAngularResolutionSST;
        
        bool  fRequirementsScaling;
        float fRequirementsScalingFactor;
        float fRequirementsScalingEnergy;
        float fRequirementsScalingStretch;
        bool  fRequirementsSystematics;
        
        void plotRequirementsSystematic( TGraph *g );
        void plotRequirements( TGraph* g, bool iLog = false, bool iLine = false, bool iSystematics = false );
        
    public:
    
        VCTARequirements();
        ~VCTARequirements() {}
        
        double  getFOVRequirement( double E_lin_TeV );
        TGraph* getRequiredDifferentalSensitivity()
        {
            return fReqDifferentialSensitivity;
        }
        string  getTitle()
        {
            return fRequirementTitle;
        }
        TGraph* plotRequirement_AngularResolution( TCanvas* c );
        TGraph* plotRequirement_EffectiveArea( TCanvas* c );
        TGraph* plotRequirement_EnergyResolution( TCanvas* c );
        TGraph* plotRequirement_DifferentialSensitivity( TCanvas* c );
        void    setRequirementsGraphLineWidth( double iLineWidth = 1. )
        {
            fRequirementsGraphLineWidth = iLineWidth;
        }
        bool    setRequirement( string iRequirement = "South-50h",
                                float iRequirementsScaling = 1. );     // updated (Dec 2017 requirements)
        void    setPlotSubSystemRequirement( string iR = "" )
        {
            fPlotSubSystemRequirement = iR;
        }
        void   setRequirementsPlotSystematics( bool iPlotSystematics = false )
        {
            fRequirementsSystematics = iPlotSystematics;
        } 
        void    setPlotRequirementsScaling( bool iPlotCTARequirementsScaling = false,
                                            float iRequirementsScalingFactor = 0.3,
                                            float iRequirementsScalingEnergy_TeV = 0.2,
                                            float iRequirementsScalingStretch = 5. )
        {
            fRequirementsScaling = iPlotCTARequirementsScaling;
            fRequirementsScalingFactor = iRequirementsScalingFactor;
            fRequirementsScalingEnergy = iRequirementsScalingEnergy_TeV;
            fRequirementsScalingStretch = iRequirementsScalingStretch;
        }
        void    setRequirementsDirectory( string iReqDir = "" );
};

#endif

