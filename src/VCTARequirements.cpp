/*!  VCTARequirements

    CTA requirements

*/

#include "VCTARequirements.h"

VCTARequirements::VCTARequirements()
{
    fSetOfRequirementID = 0;
    fRequirementTitle = "";

    fReqDifferentialSensitivity = 0;
    fReqEffectiveArea = 0;
    fReqAngularResolution = 0;
    fReqEnergyResolution = 0;

    fReqDifferentialSensitivity2014 = 0; 
    fReqDifferentialSensitivityxF = 0;
    fReqDifferentialSensitivityx2 = 0;
    fReqDifferentialSensitivityx3 = 0;
    fReqDifferentialSensitivityx4 = 0;
    fReqDifferentialSensitivityLST = 0;
    fReqDifferentialSensitivityMST = 0;
    fReqDifferentialSensitivitySST = 0;
    fReqAngularResolutionLST = 0;
    fReqAngularResolutionMST = 0;
    fReqAngularResolutionSST = 0;
    
    setRequirementsGraphLineWidth();
    setRequirementsDirectory();
    setPlotSubSystemRequirement();
    setPlotRequirementsScaling();
    setRequirementsPlotSystematics();
}


TGraph* VCTARequirements::plotRequirement_DifferentialSensitivity( 
        TCanvas* c )
{
    if( !c )
    {
        return 0;
    }
    
    if( !fReqDifferentialSensitivity )
    {
        cout << "VCTARequirements::plotRequirement_DifferentialSensitivity: error, no requirement given" << endl;
        return 0;
    }
    
    c->cd();
    plotRequirements( fReqDifferentialSensitivity, true, true, fRequirementsSystematics );
    if( fPlotSubSystemRequirement == "LST" && fReqDifferentialSensitivityLST )
    {
        plotRequirements( fReqDifferentialSensitivityLST, true, true, fRequirementsSystematics );
    }
    else if( fPlotSubSystemRequirement == "MST" && fReqDifferentialSensitivityMST )
    {
        plotRequirements( fReqDifferentialSensitivityMST, true, true, fRequirementsSystematics );
    }
    else if( fPlotSubSystemRequirement == "SST" && fReqDifferentialSensitivitySST )
    {
        plotRequirements( fReqDifferentialSensitivitySST, true, true, fRequirementsSystematics );
    }
    if( fPlotSubSystemRequirement.size() > 0
            && fReqDifferentialSensitivityxF
            && fReqDifferentialSensitivityx2 && fReqDifferentialSensitivityx3 && fReqDifferentialSensitivityx4 )
    {
        plotRequirements( fReqDifferentialSensitivityxF, true, true, fRequirementsSystematics );
        plotRequirements( fReqDifferentialSensitivityx2, true, true, fRequirementsSystematics );
        plotRequirements( fReqDifferentialSensitivityx3, true, true, fRequirementsSystematics );
        plotRequirements( fReqDifferentialSensitivityx4, true, true, fRequirementsSystematics );
    }
    if( fRequirementsScaling && fReqDifferentialSensitivity2014 )
    {
        plotRequirements( fReqDifferentialSensitivity2014, true, true, fRequirementsSystematics );
    }
    
    return fReqDifferentialSensitivity;
}

TGraph* VCTARequirements::plotRequirement_EnergyResolution( TCanvas* c )
{
    if( !c )
    {
        return 0;
    }
    
    if( !fReqEnergyResolution )
    {
        cout << "VCTARequirements::plotRequirement_EnergyResolution: error, no requirement given" << endl;
        return 0;
    }
    
    c->cd();
    plotRequirements( fReqEnergyResolution, true, true );
    
    return fReqEnergyResolution;
}

TGraph* VCTARequirements::plotRequirement_AngularResolution( TCanvas* c )
{
    if( !c )
    {
        return 0;
    }
    
    if( !fReqAngularResolution )
    {
        cout << "VCTARequirements::plotRequirement_AngularResolution: error, no requirement given" << endl;
        return 0;
    }
    
    c->cd();
    plotRequirements( fReqAngularResolution, false, true );
    
    if( fPlotSubSystemRequirement == "LST" && fReqAngularResolutionLST )
    {
        plotRequirements( fReqAngularResolutionLST, true, true );
    }
    if( fPlotSubSystemRequirement == "MST" && fReqAngularResolutionMST )
    {
        plotRequirements( fReqAngularResolutionMST, true, true );
    }
    if( fPlotSubSystemRequirement == "SST" && fReqAngularResolutionSST )
    {
        plotRequirements( fReqAngularResolutionSST, true, true );
    }
    
    return fReqAngularResolution;
}

TGraph* VCTARequirements::plotRequirement_EffectiveArea( TCanvas* c )
{
    if( !c )
    {
        return 0;
    }
    
    if( !fReqEffectiveArea )
    {
        cout << "VCTARequirements::plotRequirement_EffectiveArea: error, no requirement given" << endl;
        return 0;
    }
    
    c->cd();

    plotRequirements( fReqEffectiveArea, true, true );
    
    return fReqAngularResolution;
}

/*
 * plot requirements systematic
 *
 * hardwired CTA systematics, check code for details
 *
 */
void VCTARequirements::plotRequirementsSystematic( TGraph *g )
{
    if( !g )
    {
        return;
    }
    double x = 0.;
    double y = 0.;
    TGraphAsymmErrors *iGraphSys = new TGraphAsymmErrors( 1 );
    iGraphSys->SetFillStyle( 3244 );
    iGraphSys->SetFillStyle( 1001 );
    iGraphSys->SetFillColor( g->GetLineColor() );
    iGraphSys->SetFillColorAlpha( g->GetLineColor(), 0.10 );
    double fsys = 0.3;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        iGraphSys->SetPoint( i, x, y );
        /*if( x < log10( 0.05 ) )
        {
            fsys = 0.4;
        }
        else if( x < log10( 0.2 ) )
        {
            fsys = 0.4;
        } */
        // fsys changes continously 
        // from 30 to 50% below 200 GeV
        if( x < log10( 0.2 ) )
        {
            fsys = 0.3 - 0.2 * (x + 0.7);
        }
        else
        {
            fsys = 0.3;
        }
        iGraphSys->SetPointEYlow( i, y * fsys );
    }

    iGraphSys->Draw( "3" );
}

void VCTARequirements::plotRequirements( TGraph* g, bool iLog, bool iLine, bool iSystematics )
{
    if( !g )
    {
        return;
    }
    
    // plot requirements
    if( iLine )
    {
        if( iSystematics )
        {
            plotRequirementsSystematic( g );
        }
        g->Draw( "l" );
        return;
    }
    
    // plot requirements 
    double x = 0;
    double y = 0;
    double y_low = 0.;
    
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        
        if( iLog )
        {
            y_low = 0.5 * y;
        }
        else
        {
            y_low = 0.8 * y;
        }
        TArrow* a = new TArrow( x, y, x, y_low, 0.01, "|-|>" );
        setArrowPlottingStyle( a, g->GetMarkerColor(), 1, 2 );
        //		a->Draw();
    }
}

/*

    field of view requirements: A-PERF-0010, A-PERF-0020

*/
double VCTARequirements::getFOVRequirement( double iE_lin_TeV )
{
    if( iE_lin_TeV < 0.1 )
    {
        return 2.5;
    }
    
    // 0.1 - 300 TeV
    
    return 3.;
}

/*
 * update requirements and subsystem requirements (Dec 2017)
 *
 * see Jama or subsystem requirements document for all details
 *
 */
bool VCTARequirements::setRequirement( string iRequirement, float fRequirementsScalingFactor )
{
    
    if( iRequirement.size() == 0 )
    {
        cout << "\t no requirement given" << endl;
        return false;
    }
    cout << "Setting CTA requirement for " << iRequirement << endl;
    
    // check if this is a subsystem requirement
    if( iRequirement.find( "LST" ) != string::npos )
    {
        setPlotSubSystemRequirement( "LST" );
    }
    else if( iRequirement.find( "MST" ) != string::npos )
    {
        setPlotSubSystemRequirement( "MST" );
    }
    else if( iRequirement.find( "SST" ) != string::npos )
    {
        setPlotSubSystemRequirement( "SST" );
    }
    if( fPlotSubSystemRequirement.size() > 0 )
    {
        cout << "\t plotting subsystem requirements for " << fPlotSubSystemRequirement << endl;
    }

    fRequirementTitle = iRequirement;
    /////////////////////////////
    // angular resolution
    string iGFName = fRequirementsDirectory + "/" + iRequirement + "-AngRes.dat";
    ifstream iTMPAngRes( iGFName.c_str() );
    if( iTMPAngRes.good() )
    {
        fReqAngularResolution = new TGraph( iGFName.c_str() );
        if( !fReqAngularResolution->IsZombie() && fReqAngularResolution->GetN() > 0 )
        {
            setGraphPlottingStyle( fReqAngularResolution, 2, fRequirementsGraphLineWidth, 20, 1, 0, 2 );
        }
        else
        {
            fReqAngularResolution = 0;
        }
    }
    else
    {
        fReqAngularResolution = 0;
    }
    
    /////////////////////////////
    // Energy resolution
    iGFName = fRequirementsDirectory + "/" + iRequirement + "-ERes.dat";
    ifstream iTMPERes( iGFName.c_str() );
    if( iTMPERes.good() )
    {
        fReqEnergyResolution = new TGraph( iGFName.c_str() );
        if( !fReqEnergyResolution->IsZombie() && fReqEnergyResolution->GetN() > 0 )
        {
            setGraphPlottingStyle( fReqEnergyResolution, 2, fRequirementsGraphLineWidth, 20, 1, 0, 2 );
        }
        else
        {
            fReqEnergyResolution = 0;
        }
    }
    else
    {
        fReqEnergyResolution = 0;
    }
    
    /////////////////////////////
    // sensitivity requirement
    iGFName = fRequirementsDirectory + "/" + iRequirement + ".dat";
    fReqDifferentialSensitivity = new TGraph( iGFName.c_str() );
    if( !fReqDifferentialSensitivity->IsZombie() && fReqDifferentialSensitivity->GetN() > 0 )
    {
        setGraphPlottingStyle( fReqDifferentialSensitivity, 2, fRequirementsGraphLineWidth, 20, 1., 0, 2 );
    }
    else
    {
        fReqDifferentialSensitivity = 0;
    }
    /// sensitivity scaling
    if( TMath::Abs( fRequirementsScalingFactor - 1. ) > 1.e-3 )
    {
        double x = 0;
        double y = 0;
        for( int k = 0; k < fReqDifferentialSensitivity->GetN(); k++ )
        {
            fReqDifferentialSensitivity->GetPoint( k, x, y );
            fReqDifferentialSensitivity->SetPoint( k, x, y * fRequirementsScalingFactor );
        }
    }
    
    /////////////////////////////
    // effective area requirement
    // (only defined for 30m)
    if( iRequirement.find( "30m" ) != string::npos )
    {
        iGFName = fRequirementsDirectory + "/" + iRequirement + "-EffectiveArea.dat";
        fReqEffectiveArea = new TGraph( iGFName.c_str() );
        if( !fReqEffectiveArea->IsZombie() && fReqEffectiveArea->GetN() > 0 )
        {
            setGraphPlottingStyle( fReqEffectiveArea, 2, fRequirementsGraphLineWidth, 20, 1, 0, 2 );
        }
        else
        {
            fReqEffectiveArea = 0;
        }
    }
    else
    {
        fReqEffectiveArea = 0;
    }
    
    return true;
}

/*
 * set directory with requirements
 *
 */
void VCTARequirements::setRequirementsDirectory( string iReqDir )
{
    if( iReqDir.size() == 0 )
    {
        const char* evn_dir = gSystem->Getenv( "CTA_EVNDISP_AUX_DIR" );
        if( evn_dir )
        {
            fRequirementsDirectory = evn_dir;
            fRequirementsDirectory += "/Requirements/";
        }
        else
        {
            cout << "Error setting requirements directory ($CTA_EVNDISP_AUX_DIR not found)" << endl;
            fRequirementsDirectory = "";
        }
    }
    else
    {
        fRequirementsDirectory = iReqDir;
    }
}
