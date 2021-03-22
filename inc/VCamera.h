//! VRCamera    camera plotting routines
#ifndef VRCAMERA_H
#define VRCAMERA_H

#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TGaxis.h>
#include <TGFrame.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TText.h>
#include <TColor.h>

#include "VEvndispData.h"
#ifndef NOVBF
#include "VGPSDecoder.h"
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VCamera
{
    private:
        bool        fDebug;                       //!< debug switch
        
        VEvndispData*      fData;                 //!< pointer to data class
        unsigned int fTelescope;                  //!< this is the camera of telescope fTelescope
        
        enum cameraMode {C_CHARGE, C_TZERO, C_TRIGGER, C_HIT, C_HILO, C_TIMING, C_SUMWINDOW, C_SUMWINDOWSTART,
                         C_PEDMEAN, C_PEDVAR, C_PEDMEANLOW, C_PEDVARLOW, C_GAINS, C_GAINVARS, C_GAINSLOW, C_GAINVARSLOW,
                         C_TOFF, C_TOFFLOW, C_LOWGAIN, C_CALTZERO, C_CALTZEROLOW, C_STATUS, C_STATUSLOW,
                         C_L1, C_HV, C_CURRENTS,
                         C_TRIGGER_EVNDISP, C_TEMPLATE, C_MODEL3D, C_CLUSTERID, C_PE,
                        };
                        
        unsigned int fcameraModus;                //!< camera modus (trigger/charge/hit/timing/etc.)
        bool fBoolAllinOne;                       //!< plot all images in one camera
        TPad* fCanvas;                            //!< canvas for camera
        vector<TLatex* > fTextEvent;              //!< event info text (eventnumber, eventtime, ...) and analyis text
        vector<TLatex* > fTextMC;                 //!< MC parameters
        TLatex* fTextTelescopeN;                  //!< telescope numbering
        TLatex* fTextEventPlotPaper;
        vector<TEllipse*> fgraphTubes;            //!< graphical representation of tubes
        vector<TEllipse*> fgraphTubesEntry;       //!< graphical representation of tubes entry
        vector<TMarker*> fgraphMarker;            //!< a marker to be plot on top of the tubes
        TEllipse* fCameraOuterEdge;
        TEllipse* fCameraFOV;
        TEllipse* fAnaEllipse;                    //!< reconstructed image
        TEllipse* fAnaEllipse1;                   //!< reconstructed image
        TEllipse* fAnaEllipse2;                   //!< reconstructed image
        TEllipse* fAnaEllipseLL;                  //!< reconstructed image (loglik.)
        TMarker* fAnaShowerDir;                   //!< reconstructed shower direction
        TMarker* fAnaShowerCentroid;              //!< image centroid
        TMarker* fMCShowerDir;                    //!< MC shower direction
        TMarker* fCameraCentreDir;                //!< camera centre
        TEllipse* fCameraCentreEllipse;
        double fmaxPlot;                          //!< relative size of camera (0.5 = full canvas)
        double fmaxRad;                           //!< relativ maximum of radius for inner ellipses (1=maximum = radius of fgraphTubes)
        double fScaleMax;
        double fdist_edgeX;                        //!< maximum distance from canvas center to edge of camera
        double fdist_edgeY;                        //!< maximum distance from canvas center to edge of camera
        double fmax_dist_edge;
        unsigned int fEventCounter;               //!< event number
        int fTubeSelected;                        //!< number of mouse pointer selected channel (-1=nothing selected)
        vector< int > fTubeSelectedV;                        //!< number of mouse pointer selected channel (-1=nothing selected)
        
        bool bFixScale;
        
        int fCurrentTimeSlice;                    //!< time slice number for C_TIMING
        
        valarray< double > fPMTData;              //!< array with PMT data values
        
        vector< TEllipse* > fTheta2Circle;
        
        // things plotted while showing the channel numbers
        int fPrintChannel;                        //!< print channel numbers (0=none,1=channels,2=tubes)
        vector<TText*> fTextChannelNumber;        //!< channel number
        vector<TText*> fTextTubeNumber;           //!< tube number
        TGaxis* fCameraXaxis;
        TGaxis* fCameraYaxis;
        TLine* fEllipseLine;
        TLine* fEllipseLine_minus;
        TLine* fEllipseLine_plus;
        TLine* fCenterLine;
        
        // stuff need for color scheme mode
        unsigned int fPlotColor;                  //!< color scheme plot (1) or circles size plot (0)
        TGaxis* fColourAxis;                      //!< axis on the right to color palette
        int fncolors;                             //!< number of colors in current style
        int fndivz;                               //!< number of contours in current color style
        
        valarray<bool> fDeadChan;                 //!< vector with dead channels
        
        bool fFirstTelescopeToDraw;
        int  fTelescopeEllipseColor;
        
        // color definitions
        unsigned int fColorEmpty;
        unsigned int fColorImage;
        unsigned int fColorImageUser;
        unsigned int fColorBorder;
        unsigned int fColorEstimated;
        unsigned int fColorDead;
        unsigned int fColorFADCTrig;
        unsigned int fColorSum;
        unsigned int fColorTrigger;
        unsigned int fColorHit;
        unsigned int fColorTiming;
        unsigned int fColorPedMean;
        unsigned int fColorPedVar;
        unsigned int fFillStylePed;
        unsigned int fFillStyleEmpty;
        unsigned int fFillStyleDead;
        unsigned int fFillStylePos;
        unsigned int fFillStyleNeg;
        unsigned int fFillStyleFADCTrig;
        
        bool fAnaVis;                             //!< for pedestrians....
        bool fPlotPaper;                          //!< nicer plot for papers, no small text, no dead channels
        
        double         convertX( double x, double iOffSet = 0.5 );        //!< convert from camera to canvas coordinates
        double         convertY( double y, double iOffSet = 0.5 );        //!< convert from camera to canvas coordinates
        void           drawMuonResults();         //!< draw muon analysis results Martin
        void           drawAnaResults();          //!< draw analysis results in form of ellipse
        void           drawEventText();           //!< draw basic event info and analysis results
        void           drawStarsInFOV();
                vector< unsigned int > getDrawingMask( unsigned int iSelectionMask, valarray<double> i_data, double iMinValue = -999., bool iLowGain = false );
        //!< get maximum, exclude dead channels
        double         getMax( valarray<double>& );
        //!< get minimum, exclude dead channels
        double         getMin( valarray<double>& );
        //!< get minimum/maximum, exclude dead channels
        void           getMinMax( valarray<double>&, double& imin, double& imax, vector< unsigned int > iPixelDrawMask );
        // plot ellipse line
        void plot_ellipseLine( double iSign = 0. );
        //!< recalculate radii for fgraphTubesEntry (normalize radii)
        valarray<double>& rescaleSums( valarray<double>&, bool );
        //!< fill fgraphTubesEntry values
        void           setPMTColorOnOff( const vector<bool>&, int iColor, int iFillPos, int iFillNeg );
        void           setPMTColorOff( const vector<bool>& );
        //!< fill fgraphTubesEntry values
        void           setPMTColorScheme( vector<unsigned int>, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle, bool i_scale, bool iDrawDead, bool iLowGain = false );
        void           setPMTColorScheme( vector<float>, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle,
                                                                          bool i_scale, bool iDrawDead, bool iLowGain = false );
        void           setPMTColorScheme( valarray<unsigned int>, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle, bool i_scale, bool iDrawDead, bool iLowGain = false );
        void           setPMTColorScheme( valarray<double>, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle,
                                                                          bool i_scale, bool iDrawDead = false, bool iLowGain = false );
        void           setUpCamera();             //!< initialize the camera
        //!< setup camera in charge/timing mode
        void           setPMTColorForChargeTiming();
        
    public:
        VCamera() {}
        //!< standard constructor for telescope iTelescope
        VCamera( unsigned int iTelescope, VEvndispData* fData );
        virtual ~VCamera() {}                       //!< clean up
        void draw( double, int, bool );           //!< draw camera
        int getChannel( int x, int y );           //!< get channel at canvas position x,y after mouse click
        //!< get channel at canvas position x,y after mouse click
        int getChannel( int x, int y, TObject* objSel );
        unsigned int getTelescopeNumber()         //!< get telescope number to this camera (T1 = 0 )
        {
            return fTelescope;
        }
        unsigned int getPlotColor()               //!< color plot modus (1=color, 2=bw)
        {
            return fPlotColor;
        }
        void hideSelectedChannel();               //!< hide currently selected channel
        void showSelectedChannel( int, bool );    //!< plot selected channel bold
        void setCanvas( TPad* iCanvas );          //!< set canvas for current camera display modus
        void setCurrentTimeSlice( int iS )
        {
            fCurrentTimeSlice = iS;
        }
        void setFirstTelescopeToDraw()
        {
            fFirstTelescopeToDraw = true;
        }
        void setFixScale( bool iB )               //!< set a fix scale for color plots
        {
            bFixScale = iB;
        }
        void setMode( unsigned int i_mod );       //!< set camera modus (trigger/charge/hit/timing)
        void setPlotColor( unsigned int i_plot )  //!< set color plot modus
        {
            fPlotColor = i_plot;
        }
        void setPrintChannels( int );             //!< display channel numbers (0=none,1=channel,2=tube)
        void setTubeSelected( int iTube )
        {
            fTubeSelected = iTube;
        }
};
#endif                                            // VRCAMERA_H
