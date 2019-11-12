//! VPlotArrayReconstruction plot core and direction reconstruction in sky or camera view

#ifndef VPlotArrayReconstruction_H
#define VPlotArrayReconstruction_H

#include "TCanvas.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TLine.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMarker.h"
#include "TText.h"
#include "TTree.h"

#include <bitset>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "VGlobalRunParameter.h"
#include "VPlotUtilities.h"

using namespace std;

class VPlotArrayReconstruction : public VPlotUtilities, public TObject
{
    private:
    
        bool bZombie;
        int  fEventNumber;
        
        string fName;
        
        TFile* fFile;
        TTree* fshowerpars;
        vector< TTree* > ftpars;
        
        TTree* ftcors;
        
        unsigned short int fnimages[VDST_MAXTELESCOPES];
        vector< bool > fImgSel;
        float fMCEnergy;
        float fchi2[VDST_MAXTELESCOPES];
        float fMCxcore;
        float fMCycore;
        float fxcore[VDST_MAXTELESCOPES];
        float fycore[VDST_MAXTELESCOPES];
        float fMCaz;
        float fMCze;
        float fze[VDST_MAXTELESCOPES];
        float faz[VDST_MAXTELESCOPES];
        float fMCxoff;
        float fMCyoff;
        float fxoff[VDST_MAXTELESCOPES];
        float fyoff[VDST_MAXTELESCOPES];
        vector< bool > fltrig;
        vector< float > fsize;
        vector< float > fcen_x;
        vector< float > fcen_y;
        vector< float > fphi;
        vector< float > flength;
        vector< float > fwidth;
        vector< float > fdist;
        vector< float > floss;
        vector< UShort_t >   fntubes;
        
        double fSizeCut;
        
        vector< float > fphotons;
        
        vector< float > telPos_x;
        vector< float > telPos_y;
        vector< float > telPos_r;
        
        TCanvas* cGround;
        TH2D* hGround;
        
        TCanvas* cCamera;
        TH2D* hCamera;
        
        bool bPlotImageAxes;
        bool bPlotLocalTrigger;
        double fPlotxmin;
        double fPlotxmax;
        double fPlotymin;
        double fPlotymax;
        double fPlotPixelSize;
        double fPlotFOV;
        double fPlotTelescopeScale;
        
        int getEvent( int iEventNumber = 0 );
        
        int get_corsikaIOreader_event( int iEventNumber = 0 );
        int get_eventdisplay_event( int iEventNumber = 0 );
        void read_corsikaIOreader_file();
        void read_eventdisplay_file();
        
    public:
    
        VPlotArrayReconstruction( string iname = "", string ifile = "" );
        ~VPlotArrayReconstruction();
        
        TCanvas* cameraPlot( int iEventNumber = 0 );
        TCanvas* groundPlot( int iEventNumber = 0 );
        void plotEvent( int iEventNumber = 0 );
        bool isZombie()
        {
            return bZombie;
        }
        void nextEvent();
        void printEvent( int iEventNumber );
        void printEvents();
        void printTelescopes();
        void printToFile( string iFileName, string iFileType = ".eps" );
        void setPlotLocalTrigger( bool iB = true )
        {
            bPlotLocalTrigger = iB;
        }
        void setPlotImageAxis( bool iB = true )
        {
            bPlotImageAxes = iB;
        }
        void setPlotGroundCoordinates( double x1 = -550., double x2 = 550., double y1 = -550., double y2 = 550. );
        void setPlotFOV( double iFOV = 8. )
        {
            fPlotFOV = iFOV;
        }
        void setPlotPixelSize( double iP = 0.1 )
        {
            fPlotPixelSize = iP;
        }
        void setPlotTelescopeScale( double iScale = 2. )
        {
            fPlotTelescopeScale = iScale;
        }
        void setSizeCut( double iS = 400 )
        {
            fSizeCut = iS;
        }
        
        ClassDef( VPlotArrayReconstruction, 2 );
};
#endif
