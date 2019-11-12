/*
Plot average high/low gain pulses from DST files. Output similar to the one from LowGainCalibrationPulseShapeforCare.C .

Usage example:

	VPlotPulseShapesFromDST data( "$VERITAS_USER_DATA_DIR/analysis/Results/v480_DSTtest/64081.DST.root");
	data.readPulseShapes( 2, 42 ); 	//Tel, channel.
	c = new TCanvas("e", "c", 1200,800);
	data.plot(c,-6,18); //canvas, minx, maxx,miny,maxy.


*/


#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TProfile2D.h"
#include "TList.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TFile.h"

using namespace std;

class VPlotPulseShapesFromDST
{

    public:
    
        TGraphErrors* fHighGainPulse;
        TGraphErrors* fLowGainPulse;
        TMultiGraph* fHighGainPulseM;
        TMultiGraph* fLowGainPulseM;
        
        TGraphErrors* readPulseShape( TProfile2D* hist, int starBin, int StopBin );
        
        VPlotPulseShapesFromDST( TString filename );
        ~VPlotPulseShapesFromDST() {};
        bool readPulseShapes( int tel, int pixel );
        void plot( TCanvas* c = 0, double x_min = -4, double x_max = 20, double y_min = -1.1, double y_max = 0.1 );
        
        TFile* file;
        
};

/*
Constructor + read in DST file.
*/

VPlotPulseShapesFromDST::VPlotPulseShapesFromDST( TString filename )
{
    file = new TFile( filename.Data(), "read" );
    if( !file || file->IsZombie() )
    {
        cout << "VPlotPulseShapesFromDST error: File " << filename << " not found." << endl;
        file = 0;
    }
    
    fHighGainPulse = 0;
    fLowGainPulse = 0;
    fHighGainPulseM = new TMultiGraph();
    fLowGainPulseM = new TMultiGraph();
}

/*
DST file has mean pulses as 2D profile histograms (pulse vs log(collected charge) ).
Get the average pulse using TH2D::ProjectionY().

Note: Bins are currently hardcoded for high gain, unsaturated lg, saturated lg for sw 16.

*/
bool VPlotPulseShapesFromDST::readPulseShapes( int tel, int pixel )
{
    TString name;
    name.Form( "meanPulses/hPulseHigh_%d_%d", tel, pixel );
    TProfile2D* high = ( TProfile2D* )file->Get( name.Data() );
    fHighGainPulse = readPulseShape( high, 25, 30 );
    fHighGainPulse->SetLineColor( kBlack );
    fHighGainPulse->SetLineWidth( 2 );
    fHighGainPulseM->Add( fHighGainPulse );
    
    name.Form( "meanPulses/hPulseLow_%d_%d", tel, pixel );
    TProfile2D* low = ( TProfile2D* )file->Get( name.Data() );
    fLowGainPulse = readPulseShape( low, 30, 34 );
    fLowGainPulse->SetLineColor( kRed );
    fLowGainPulse->SetLineWidth( 2 );
    fLowGainPulseM->Add( fLowGainPulse );
    
    /*	TGraphErrors * temp = readPulseShape( low, 34, 36 );
    	temp->SetLineColor( 15 );
    	temp->SetLineWidth( 2 );
    	fLowGainPulseM->Add( temp );
    */
    TGraphErrors* temp = readPulseShape( low, 37, 45 );
    temp->SetLineColor( 15 );
    temp->SetLineWidth( 2 );
    fLowGainPulseM->Add( temp );
    
    return ( fHighGainPulse && fLowGainPulse );
}

TGraphErrors* VPlotPulseShapesFromDST::readPulseShape( TProfile2D* hist, int startBin, int stopBin )
{
    TH1D* projection = hist->ProjectionY( "_py" , startBin, stopBin );
    projection->Scale( -1.0 / projection->GetMaximum() );
    TGraphErrors* graph = new TGraphErrors( projection );
    return graph;
}

void VPlotPulseShapesFromDST::plot( TCanvas* c , double x_min, double x_max, double y_min, double y_max )
{
    if( !c )
    {
        c = new TCanvas( "cD", "high and low-gain pulse shapes" );
    }
    c->Draw();
    
    if( !fHighGainPulseM )
    {
        return;
    }
    
    fHighGainPulseM->Draw( "alp" );
    fLowGainPulseM->Draw( "lp" );
    fHighGainPulseM->GetXaxis()->SetTitle( "sample #" );
    fHighGainPulseM->GetXaxis()->SetLimits( x_min, x_max );
    fHighGainPulseM->GetXaxis()->SetNdivisions( 12, 5, 0, false );
    
    fHighGainPulseM->GetHistogram()->SetMinimum( y_min );
    fHighGainPulseM->GetHistogram()->SetMaximum( y_max );
    
    c->SetGrid();
    c->Update();
    
}

