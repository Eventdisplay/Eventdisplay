/*  VLowGainCalibrationFromPulsShape
    macro to calculate high/low gain multiplier from pulse shapes for high and low gain

    input:
    pulse shape in ascii format
    (e.g. see $VERITAS_EVNDISP_AUX/Calibration/CareSimulations/VERITASLowGainPulseShapesUpgradePMTFromFADC.txt
              $VERITAS_EVNDISP_AUX/Calibration/CareSimulations/VERITASHighGainPulseShapesUpgradePMTFromFADC.txt)

    output:
    LOWGAINMULTIPLIER_SUM line as required in calibrationlist.LowGain.dat

    see VERITAS wiki (EVNDISP manual) for a detailed description:
    Eventdisplay_Manual:_calibration#Low-gain_multiplier_calibration_for_CARE_simulations

    IMPORTANT:

    - do not change any of the default parameters, as they reflect the current standing of the VERITAS analysis

*/


#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TList.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"

class VLowGainCalibrationFromPulsShape
{
    private:
    
        double fFADCSampling;                         // sampling time in [ns]
        
        // time axis is always samples (not [ns]!)
        TGraph* fHighGainPulse;
        TGraph* fLowGainPulse;
        
        //To hold more than one trace.
        TMultiGraph* fHighGainPulseM;
        TMultiGraph* fLowGainPulseM;
        
        TGraph* fLGRatio_window;
        
        // start and stop of high-gain integration window (determined by trace integration algorithm)
        double fHighGainIntegrationWindowStart;
        double fHighGainIntegrationWindowStopp;
        
        double fWindowStartRelT0_sample;
        double fHighGainWindowLength_sample;
        
        // start and stop of low-gain integration window (fixed)
        double fLowGainWindowStart_sample;
        double fLowGainWindowStop_sample;
        
        // trace integration method (1 or 2)
        unsigned int fTraceIntegrationMethod;
        
        double fFADCLOHIGHGAINRATIO;                   // from CARE simulation configuration file (same keyword)
        
        
        double  getTraceTZero( TGraph* g );
        double  getTraceMaxIntegral_windowStart( TGraph* g, double iWindowLength_sample );
        double  integrate( TGraph* g, double iWindowStartRelT0_sample, double iWindowLength_sample, int iIntegrationMethod );
        TGraph* readPulseShape( string iFile, bool iGrisuFormat, double iSampling );
        TMultiGraph* readMultiplePulseShape( string iFile, bool iGrisuFormat, double iSampling );
        
        
    public:
        VLowGainCalibrationFromPulsShape();
        ~VLowGainCalibrationFromPulsShape() {}
        void analyse();
        void calculateChargeFraction();
        bool readPulseShapes( string iHighGainPulse, string iLowGainPulse, bool iGrisuFormat = false, double iSampling_ns = 2.0 );
        bool readMultiplePulseShapes( string iHighGainPulse, string iLowGainPulse, bool iGrisuFormat = false, double iSampling_ns = 2.0 );
        void plot( double y_min = -99., double y_max = -99 );
        void printIntegrationParameters();
        void plotM( TCanvas* c = 0, double x_min = 0, double x_max = 24, double y_min = -1.1, double y_max = 0.1 );
        void setFADCLOHIGHGAINRATIO( double iValue = 0.099 )
        {
            fFADCLOHIGHGAINRATIO = iValue;
        }
        void setTraceIntegrationParameters( double iStartRelT0_sample = -1., int iHighGainWindowLength_sample = 16.,
                                            int iLowGainWindowStart_sample = 4, int iLowGainWindowStop_sample = 18 );
        void setTraceIntegrationMethod( int iTraceIntegrationMethod = 2, bool iSetDefaultParameters = true );
};

/////////////////////////////

VLowGainCalibrationFromPulsShape::VLowGainCalibrationFromPulsShape()
{
    fFADCSampling = 2.;
    fHighGainPulse = 0;
    fLowGainPulse = 0;
    fHighGainPulseM = 0;
    fLowGainPulseM = 0;
    fLGRatio_window = 0;
    
    setTraceIntegrationParameters();
    setTraceIntegrationMethod();
    setFADCLOHIGHGAINRATIO();
    
    fHighGainIntegrationWindowStart = -99.;
    fHighGainIntegrationWindowStopp = -99.;
}

void VLowGainCalibrationFromPulsShape::printIntegrationParameters()
{
    cout <<  endl;
    cout << "**** Trace integration parameters: ***" << endl;
    cout << "Integrationsmethod: " << fTraceIntegrationMethod << endl;
    cout << "Sample length [ns]: " << fFADCSampling << endl;
    cout << "HG window length [samples]: " << fHighGainWindowLength_sample << endl;
    cout << "Window start relative T0 [samples]: " << fWindowStartRelT0_sample << endl;
    cout << "LG window start/stop [samples]: " << fLowGainWindowStart_sample << ", " << fLowGainWindowStop_sample << endl;
    cout << endl;
}

/*

    set trace integration method (and possibly the default parameters)

    note: for double pass method use the DP_1 method here

*/
void VLowGainCalibrationFromPulsShape::setTraceIntegrationMethod( int iTraceIntegrationMethod, bool iSetDefaultParameters )
{
    fTraceIntegrationMethod = iTraceIntegrationMethod;
    
    // set default parameters for the chosen trace integration method
    if( iSetDefaultParameters )
    {
        if( iTraceIntegrationMethod == 1 )
        {
            setTraceIntegrationParameters( -1, 16, 4, 16 );
        }
        else if( iTraceIntegrationMethod == 2 )
        {
            setTraceIntegrationParameters( 0, 16, 4, 16 );
        }
    }
}



/*

     run the low-gain calibration analysis with standard parameters and plot results

*/
void VLowGainCalibrationFromPulsShape::analyse()
{
    readPulseShapes( "VERITASHighGainPulseShapeUpgradePMTFromFADC.txt", "VERITASLowGainPulseShapesUpgradePMTFromFADC.txt" );
    calculateChargeFraction();
    plot();
}

/*


*/
void VLowGainCalibrationFromPulsShape::setTraceIntegrationParameters( double iStartRelT0_sample, int iHighGainWindowLength_sample,
        int iLowGainWindowStart_sample, int iLowGainWindowStop_sample )
{
    fWindowStartRelT0_sample     = iStartRelT0_sample;
    fHighGainWindowLength_sample = iHighGainWindowLength_sample;
    fLowGainWindowStart_sample   = iLowGainWindowStart_sample;
    fLowGainWindowStop_sample    = iLowGainWindowStop_sample;
}

/*

    read a pulse shape from an ASCII file and fill it into a TGraph

*/
TGraph* VLowGainCalibrationFromPulsShape::readPulseShape( string iFile, bool iGrisuFormat, double iSampling )
{
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        return 0;
    }
    
    TGraph* g = new TGraph( 1 );
    
    string iTemp = "";
    string is_line = "";
    double x = 0.;
    double y = 0.;
    int z = 0;
    double iSampleOffset = 0.;
    while( getline( is, is_line ) )
    {
        if( is_line.find( "*" ) != string::npos )
        {
            //			continue;
            break;
        }
        istringstream is_stream( is_line );
        if( !iGrisuFormat )
        {
            if( !is_stream.eof() )
            {
                is_stream >> x;
                if( fFADCSampling > 0. )
                {
                    x /= fFADCSampling;
                }
            }
        }
        else
        {
            x = iSampleOffset;
        }
        if( !is_stream.eof() )
        {
            is_stream >> y;
        }
        if( x > 1000. )
        {
            continue;
        }
        g->SetPoint( z, x, y );
        z++;
        iSampleOffset += iSampling;
    }
    return g;
}

/*

    read multiple pulse shapes from an ASCII file and fill it into a TMultiGraph.
    Assume the traces are separated by a line starting with *.

    Colors are hardcoded. If there is only one trace, it is assumed to be a hg trace -> black.
    If there are several, the first one will be red, the others grey.

*/
TMultiGraph* VLowGainCalibrationFromPulsShape::readMultiplePulseShape( string iFile, bool iGrisuFormat, double iSampling )
{
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        return 0;
    }
    
    TMultiGraph* m = new TMultiGraph();
    TGraph* g = 0;
    
    string iTemp = "";
    string is_line = "";
    double x = 0.;
    double y = 0.;
    int z = 0;
    int col = 2;
    double iSampleOffset = 0.;
    while( getline( is, is_line ) )
    {
        if( is_line.find( "*" ) != string::npos )
        {
            //End of trace.
            g->SetLineColor( col );
            g->SetLineWidth( 2 );
            m->Add( g );
            g = 0;
            col = 15;
            continue;
        }
        istringstream is_stream( is_line );
        if( !iGrisuFormat )
        {
            if( !is_stream.eof() )
            {
                is_stream >> x;
                if( fFADCSampling > 0. )
                {
                    x /= fFADCSampling;
                }
            }
        }
        else
        {
            x = iSampleOffset;
        }
        if( !is_stream.eof() )
        {
            is_stream >> y;
        }
        if( x > 1000. )
        {
            continue;
        }
        if( !g )
        {
            g = new TGraph( 0 );
        }
        
        g->SetPoint( z, x, y );
        z++;
        iSampleOffset += iSampling;
    }
    if( g )
    {
        g->SetLineWidth( 2 );
        g->SetLineColor( col == 2 ? 1 : 15 );
        m->Add( g );
    }
    return m;
}

/*

     read high and low-gain pulse from ascii files

*/
bool VLowGainCalibrationFromPulsShape::readPulseShapes( string iHighGainPulse, string iLowGainPulse, bool iGrisuFormat, double iSampling )
{
    cout << "Reading high gain pulse from " << iHighGainPulse << endl;
    fHighGainPulse = readPulseShape( iHighGainPulse, iGrisuFormat, iSampling );
    if( !fHighGainPulse || fHighGainPulse->GetN() < 1 )
    {
        cout << "error reading high gain pulse shape (" << fHighGainPulse << ")" << endl;
        return false;
    }
    fHighGainPulse->SetName( "HighGainPulseShape" );
    fHighGainPulse->SetTitle( "" );
    cout << "Reading low gain pulse from " << iLowGainPulse << endl;
    fLowGainPulse = readPulseShape( iLowGainPulse, iGrisuFormat, iSampling );
    if( !fLowGainPulse || fLowGainPulse->GetN() < 1 )
    {
        cout << "error reading low gain pulse shape" << endl;
        return false;
    }
    fLowGainPulse->SetName( "LowGainPulseShape" );
    fLowGainPulse->SetLineColor( 2 );
    fLowGainPulse->SetTitle( "" );
    
    return true;
}

/*

     read multiplie high and low-gain pulses from ascii files

*/
bool VLowGainCalibrationFromPulsShape::readMultiplePulseShapes( string iHighGainPulse, string iLowGainPulse, bool iGrisuFormat, double iSampling )
{
    cout << "Reading high gain pulse from " << iHighGainPulse << endl;
    fHighGainPulseM = readMultiplePulseShape( iHighGainPulse, iGrisuFormat, iSampling );
    if( !fHighGainPulseM )
    {
        cout << "error reading high gain pulse shape (" << fHighGainPulse << ")" << endl;
        return false;
    }
    fHighGainPulseM->SetName( "HighGainPulseShape" );
    fHighGainPulseM->SetTitle( "" );
    cout << "Reading low gain pulse from " << iLowGainPulse << endl;
    fLowGainPulseM = readMultiplePulseShape( iLowGainPulse, iGrisuFormat, iSampling );
    if( !fLowGainPulseM )
    {
        cout << "error reading low gain pulse shape" << endl;
        return false;
    }
    fLowGainPulseM->SetName( "LowGainPulseShape" );
    fLowGainPulseM->SetTitle( "" );
    
    return true;
}

/*

    plot all results

*/
void VLowGainCalibrationFromPulsShape::plot( double y_min, double y_max )
{
    TCanvas* c = new TCanvas( "cP", "high and low-gain pulse shapes" );
    c->Draw();
    
    if( !fHighGainPulse )
    {
        return;
    }
    
    fHighGainPulse->Draw( "al" );
    fLowGainPulse->Draw( "l" );
    fHighGainPulse->GetXaxis()->SetTitle( "sample #" );
    
    TLegend* iL = new TLegend( 0.55, 0.42, 0.85, 0.62 );
    iL->AddEntry( fHighGainPulse, "high-gain pulse", "l" );
    iL->AddEntry( fLowGainPulse, "low-gain pulse", "l" );
    iL->Draw();
    
    if( fHighGainIntegrationWindowStart > 0. && fHighGainIntegrationWindowStopp > 0. )
    {
        TLine* iL_start = new TLine( fHighGainIntegrationWindowStart, 0., fHighGainIntegrationWindowStart, -1. );
        iL_start->SetLineStyle( 2 );
        iL_start->Draw();
        TLine* iL_stopp = new TLine( fHighGainIntegrationWindowStopp, 0., fHighGainIntegrationWindowStopp, -1. );
        iL_stopp->SetLineStyle( 2 );
        iL_stopp->Draw();
    }
    
    if( fLGRatio_window && fLGRatio_window->GetN() > 0 )
    {
        TCanvas* d = new TCanvas( "cLGRatio", "LG ratio", 200, 200, 700, 500 );
        d->Draw();
        
        fLGRatio_window->SetTitle( "" );
        fLGRatio_window->SetLineWidth( 2. );
        fLGRatio_window->SetMarkerStyle( 21 );
        fLGRatio_window->Draw( "alp" );
        fLGRatio_window->GetXaxis()->SetTitle( "integration window (samples)" );
        char hname[500];
        fLGRatio_window->GetYaxis()->SetTitleOffset( 1.05 );
        if( y_min > 0. )
        {
            fLGRatio_window->SetMinimum( y_min );
        }
        if( y_max > 0. )
        {
            fLGRatio_window->SetMaximum( y_max );
        }
        sprintf( hname, "SW-dependent LG multiplier (HG %d)",
                 ( int )fHighGainWindowLength_sample );
        fLGRatio_window->GetYaxis()->SetTitle( hname );
    }
}

/*

    plot multigraphs possibly holding several traces each.

*/
void VLowGainCalibrationFromPulsShape::plotM( TCanvas* c, double x_min, double x_max, double y_min, double y_max )
{
    if( !c )
    {
        c = new TCanvas( "cP", "high and low-gain pulse shapes" );
    }
    c->Draw();
    
    if( !fHighGainPulseM || !fLowGainPulseM )
    {
        return;
    }
    
    fHighGainPulseM->Draw( "al" );
    fLowGainPulseM->Draw( "l" );
    fHighGainPulseM->Draw( "l" );
    fHighGainPulseM->GetXaxis()->SetTitle( "sample #" );
    fHighGainPulseM->GetXaxis()->SetLimits( x_min, x_max );
    fHighGainPulseM->GetHistogram()->SetMinimum( y_min );
    fHighGainPulseM->GetHistogram()->SetMaximum( y_max );
    c->SetGrid();
    c->Update();
    
    
    TLegend* iL = new TLegend( 0.55, 0.42, 0.85, 0.62 );
    iL->AddEntry( fHighGainPulseM, "high-gain pulse", "l" );
    iL->AddEntry( fLowGainPulseM, "low-gain pulse", "l" );
    //iL->Draw();
    
}

/*

*/
void VLowGainCalibrationFromPulsShape::calculateChargeFraction()
{
    cout << endl;
    printIntegrationParameters();
    
    // make sure that all data is available
    if( !fHighGainPulse || !fLowGainPulse )
    {
        return;
    }
    
    /////////////////////////////////////////
    // from CARE simulations configuration file
    // CARE uses a hi/lo gain ratio on the amplitudes that describes the characteristics of
    // the system where the lo gain is linear.
    // * FADCLOHIGHGAINRATIO 0 0.099
    cout << "Hi/lo normalization: " << fFADCLOHIGHGAINRATIO << endl;
    cout << endl;
    
    // hilo gain graph and vector
    fLGRatio_window = new TGraph( 1 );
    int z = 0;
    vector< double > iLGRatio_V;
    
    ///////////////////////////////////////////////////////////
    // integrate high gain pulse
    cout << "**** high-gain integration: ***" << endl;
    double iHighGainIntegral = integrate( fHighGainPulse, fWindowStartRelT0_sample, fHighGainWindowLength_sample, fTraceIntegrationMethod );
    
    ///////////////////////////////////////////////////////////
    // loop over all relevant window sizes and integrate low-gain pulses
    for( int i = fLowGainWindowStart_sample; i <= fLowGainWindowStop_sample; i++ )
    {
        cout << "Low-gain integration: " << endl;
        double iLowGainIntegral = integrate( fLowGainPulse, fWindowStartRelT0_sample, i, fTraceIntegrationMethod )  * fFADCLOHIGHGAINRATIO;
        cout << "**** low/high gain ratio ***" << endl;
        cout << "Window length " << i << ": " << iHighGainIntegral << "\t" << iLowGainIntegral;
        if( iLowGainIntegral > 0. && iHighGainIntegral > 0. )
        {
            cout << " ratio: " << iHighGainIntegral / iLowGainIntegral;
            cout << " LG : " << iHighGainIntegral / iLowGainIntegral;
            fLGRatio_window->SetPoint( z, i, iHighGainIntegral / iLowGainIntegral );
            iLGRatio_V.push_back( iHighGainIntegral / iLowGainIntegral );
            z++;
        }
        cout << endl;
        cout << "*******" << endl;
    }
    
    ///////////////////////////////////////////////////////////
    // print results to screen so that they can be copied directly into
    // calibrationlist.LowGainForCare.dat
    for( unsigned int itel = 1; itel <= 4; itel++ )
    {
        cout << "* LOWGAINMULTIPLIER_SUM " << itel << " 0     999999   ";
        cout << "  " << fTraceIntegrationMethod << "  ";
        cout << fHighGainWindowLength_sample << "   ";
        cout << fLowGainWindowStart_sample << "   ";
        cout << fLowGainWindowStop_sample << "   ";
        
        for( unsigned int i = 0; i < iLGRatio_V.size(); i++ )
        {
            cout << "    " << iLGRatio_V[i];
        }
        cout << endl;
    }
    
}

/*
    get trace T0

   (defined as position of half max (rising edge))

*/
double VLowGainCalibrationFromPulsShape::getTraceTZero( TGraph* g )
{
    if( !g )
    {
        return -999.;
    }
    
    
    double x_t0 = 0.;
    
    double x = 0.;
    double y = 0.;
    double ymax_half = 0.;
    /////////////////////////////////
    // get position of trace maximum
    // (note negative trace values)
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        if( y < ymax_half )
        {
            ymax_half = y;
        }
    }
    /////////////////////////////////
    // get position of trace half maximum
    // (note negative trace values)
    ymax_half *= 0.5;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        if( y < ymax_half )
        {
            x_t0 = x;
            break;
        }
    }
    return x_t0;
}

/*
    search along trace for maximum sum
    return window start

    (same as VTraceHandler::getQuickMaximumSum)

*/
double VLowGainCalibrationFromPulsShape::getTraceMaxIntegral_windowStart( TGraph* g, double iWindowLength_sample )
{
    if( !g )
    {
        cout << "VLowGainCalibrationFromPulsShape::getTraceAverageTime error, no graph" << endl;
        return -999.;
    }
    vector< double > fpTrace;
    vector< double > fSample;
    double x = 0.;
    double y = 0.;
    // initialize data vector
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        fSample.push_back( x );
        fpTrace.push_back( -1.*y );
    }
    float xmax = 0.;
    float charge = 0.;
    float t = 0;
    float t_start = fSample[0];
    
    for( unsigned int i = 0; i < fSample.size(); i++ )
    {
        xmax += fpTrace[i];
        if( fSample[i] - t_start  > iWindowLength_sample )
        {
            break;
        }
    }
    unsigned int lolimit = 0;
    unsigned int uplimit = 0;
    for( unsigned int i = 0; i < fSample.size(); i++ )
    {
        // get index for end of integration window
        unsigned int j_window_end = 0;
        for( unsigned int j = i; j < fSample.size(); j++ )
        {
            if( fSample[j] - fSample[i] > iWindowLength_sample )
            {
                j_window_end = j;
                break;
            }
        }
        if( j_window_end == 0 && fSample.size() > 0 )
        {
            j_window_end = fSample.size() - 1;
        }
        // new maximum
        if( charge < xmax )
        {
            charge = xmax;
            lolimit = i;
            uplimit = j_window_end;
        }
        xmax = xmax - fpTrace[i] + fpTrace[j_window_end];
    }
    // extract charge for small window **********************************
    
    
    float tcharge = 0.;
    float arrtime = 0.;
    // arrival times *****************************************************
    cout << "\t maximum integral parameters: ";
    cout << lolimit << "\t" << uplimit << "\t" << fSample[lolimit] << "\t" << fSample[uplimit] << endl;
    for( unsigned int k = lolimit; k < uplimit; k++ )
    {
        tcharge += fpTrace[k] * fSample[k];
    }
    cout << "\t" << tcharge << "\t" << charge << "\t" << tcharge / charge << endl;
    if( charge != 0. )
    {
        arrtime = tcharge / charge;
    }
    
    // get start and stopp of trace integration for high-gain pulse
    string iName = g->GetName();
    if( iName.find( "High" ) != string::npos )
    {
        fHighGainIntegrationWindowStart = fSample[lolimit];
        fHighGainIntegrationWindowStopp = fSample[uplimit];
    }
    
    return fSample[lolimit];
}

/*

     integrate a given pulse

     iIntegrationMethod:

     1: start at Tzero
     2: start at average trace time (from sliding window)

*/

double VLowGainCalibrationFromPulsShape::integrate( TGraph* g, double iWindowStartRelT0_sample,
        double iWindowLength_sample,
        int iIntegrationMethod )
{
    if( !g )
    {
        return -999.;
    }
    
    double x_t0 = 0.;
    // window start in T0
    if( iIntegrationMethod == 1 )
    {
        x_t0 = getTraceTZero( g );
    }
    // average pulse time used for window start
    // sliding window mechanism givining maximum integral
    else if( iIntegrationMethod == 2 )
    {
        x_t0 = getTraceMaxIntegral_windowStart( g, iWindowLength_sample );
    }
    else
    {
        cout << "VLowGainCalibrationFromPulsShape::integrate: unknown trace integration method:";
        cout << fTraceIntegrationMethod << endl;
        return -999.;
    }
    
    //////////////////////////////////
    // determine start and end of integration window relative to typical trace time
    double x_start = x_t0 + iWindowStartRelT0_sample;
    double x_stop  = x_start + iWindowLength_sample;
    
    int N = 10000;
    double i_sum = 0.;
    for( int i = 0; i < iWindowLength_sample * N; i++ )
    {
        i_sum += TMath::Abs( g->Eval( x_start + ( double )i / N ) );
    }
    i_sum /= N;
    
    cout << "\tT0: " << x_t0 << " (T0 offset: " << iWindowStartRelT0_sample << ")" << endl;
    cout << "\tintegration window of length " << iWindowLength_sample << ": ";
    cout << "[" << x_start << ", " << x_stop << "]" << endl;
    cout << "\ttrace sum: " << i_sum << endl;
    
    return i_sum;
}
