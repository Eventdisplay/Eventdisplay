//! VDisplay GUI class for event display
//     fCanvasCharge->Connect( "ProcessedEvent(Int_t , Int_t , Int_t , TObject* )", "LCamera", this, "Print( Int_t, Int_t, Int_t, TObject* )" );

#ifndef VRDISPLAY_H
#define VRDISPLAY_H

#ifdef __APPLE__
#ifdef __CINT__
#undef __GNUC__
#endif
#endif

#include <TApplication.h>
#include <TCanvas.h>
#include <TError.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGClient.h>
#include <TGComboBox.h>
#include <TGFileDialog.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGMenu.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGraph.h>
#include <TGStatusBar.h>
#include <TGTableLayout.h>
#include <TGTab.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGWindow.h>
#include <TH1D.h>
#include <TQObject.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
#include <RQ_OBJECT.h>

#include <VCamera.h>
#include <VDisplayBirdsEye.h>
#include <VEventLoop.h>
#include <VPETree.h>

#include <bitset>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VDisplay : public TGMainFrame
{
        RQ_OBJECT( "VDisplay" )
        
        //! menu identifier
        enum E_menuIdentifer
        {
            M_FILE_PRINTC, M_FILE_PRINTF, M_FILE_PRINTA, M_FILE_PRINTCA, M_FILE_PRINTB, M_FILE_PRINTTG,
            M_FILE_MOVIE, M_FILE_EXIT, M_FILE_DUMP,
            M_HISTOS_FILL, M_HISTOS_DRAW,
            M_OPT_COL_SCHE, M_OPT_BW_SCHE, M_OPT_SING, M_OPT_FIX, M_OPT_DEAD, M_OPT_IMAGE,
            B_AUTOSTART, B_NGOTO, B_NAUTO, B_TELN, B_NINC, B_VIEW,
            B_NTRI, B_NATRI, B_NAIMA, B_NEXT, B_NCALIB,
            B_SET, B_RESET, B_TEXTCUT,
            B_TALL,
            B_DCOLOR, B_DNONE, B_DCHANNEL, B_DTUBE, B_DCOOR,
            B_CAUTO, B_NSCALE,
            // calibration signals
            P_LPED, P_CPED, P_LGAIN, P_CGAIN, P_LTOFF, P_CTOFF, P_LPAD, P_CPAD,
            P_TPED, P_TGAIN, P_TTOFF, P_TPAD,
            B_NWIN, B_NSUM, B_NIMA, B_NBOR,
            B_FITPMT, B_ADDPMT, B_REMOVEPMT, B_RESETPMT, B_NFADC, B_NVIEW
        };
        // camera display identifier (each of them will appear in the drop down menu for a camera display)
        // note: these have to be copied into VCamera
        enum E_cameraIdent {C_CHARGE, C_TZERO, C_TRIGGER, C_HIT, C_HILO, C_TIMING, C_SUMWINDOW, C_SUMWINDOWSTART,
                            C_PEDMEAN, C_PEDVAR, C_PEDMEANLOW, C_PEDVARLOW, C_GAINS, C_GAINVARS, C_GAINSLOW, C_GAINVARSLOW,
                            C_TOFF, C_TOFFLOW, C_LOWGAIN, C_CALTZERO, C_CALTZEROLOW, C_STATUS, C_STATUSLOW,
                            C_L1, C_HV, C_CURRENTS,
                            C_TRIGGER_EVNDISP, C_TEMPLATE, C_MODEL3D, C_CLUSTERID, C_PE
                           };
        //! FADC/ana tab identifier
        enum E_fadcIDENT {F_FADC, F_ANA};
        
    private:
        bool fDebug;
        VEventLoop* fEventLoop;                   //!< main event loop, steering of data reading
        //  vector<VCamera* > fCamera;                   //!< vector of fNTelescopes cameras
        map<unsigned int, VCamera* > fCamera;     //!< vector of fNTelescopes cameras
        vector< unsigned int > fTelescopesToShow; //!< vector with telescope numbers for plotting
        VDisplayBirdsEye* fBirdsEye;                 //!< drawing the telescopes from above
        
        unsigned int fTelescope;                  //!< telesope to draw (first telescope = 0)
        bool fBoolDrawOne;                        //!< draw only one telescope or all X in window
        bool fBoolDrawAllinOne;                   //!< draw all results into one camera
        
        TGCompositeFrame* fFrameTable;            //!< main frame, containing all widgets
        TGLayoutHints* fL1;
        TGLayoutHints* fL2;
        TGLayoutHints* fL3;
        TGLayoutHints* fL4;
        TGLayoutHints* fL5;
        TGLayoutHints* fL6;
        TGLayoutHints* fL7;
        TGLayoutHints* fL8;
        TGStatusBar* fStatusBar;
        TGMenuBar* fMenuBar;
        TGPopupMenu* fMenuFile;
        TGPopupMenu* fMenuPrint;
        TGPopupMenu* fMenuOpt;
        bool   fBool_M_OPT_COL_SCHE_Checked;
        bool   fBool_M_OPT_BW_SCHE_Checked;
        TGCompositeFrame* fFrameControl;
        TGHorizontalFrame* fHorizontalEvent;
        TGTextButton* fButtonNext;
        TGTextButton* fButtonExit;
        TGTextButton* fButtonDump;
        TGLabel* fLabelGoto;
        TGNumberEntry* fNEntryGoto;
        TGTextButton* fButtonAutorun;
        TGComboBox* fComboTelescopeN;
        TGComboBox* fComboCameraView;
        map< unsigned int, TPad* > fPadsCamera;
        TCanvas* fCanvasCamera;
        unsigned int fCanvasesNx;
        unsigned int fCanvasesNy;
        vector<TGCompositeFrame*> fFrameCamera;
        vector< TPad* > fPadTelescope;
        TRootEmbeddedCanvas* fEmCanvasCamera;
        vector<string> fStringCameraTab;
        TGTab* fTabAna;
        TGCompositeFrame* fFrameFADC;
        TGCompositeFrame* fCompFADC;
        TRootEmbeddedCanvas* fEmFADC;
        TCanvas* fCanvasFADC;
        TRootEmbeddedCanvas* fEmFADCText;
        TCanvas* fCanvasFADCText;
        vector<TText* > fTextFADC;
        TGHorizontalFrame* fGroupFADC;
        TGTextButton* fButtonFADCFit;
        TGTextButton* fButtonFADCset;
        TGTextButton* fButtonFADCunset;
        TGTextButton* fButtonFADCreset;
        TGLabel* fLabelFADCsearch;
        TGNumberEntry* fNEntryFADCsearch;
        TGCompositeFrame* fFrameBird;
        TGCompositeFrame* fCompBird;
        TRootEmbeddedCanvas* fEmBird;
        TCanvas* fCanvasBird;
        
        TGCompositeFrame* fFrameCal;
        TRootEmbeddedCanvas* fEmCal;
        TCanvas* fCanvasCal;
        
        TGCompositeFrame* fFrameTgrad;
        TRootEmbeddedCanvas* fEmTgrad;
        TCanvas* fCanvasPixelHisto;
        
        TGCompositeFrame* fFrameAna;
        TGCompositeFrame* fCompAna;
        TRootEmbeddedCanvas* fEmAna;
        TCanvas* fCanvasAna;
        TGCompositeFrame* fFrameOpt;
        TGCompositeFrame* fCompOpt;
        TGGroupFrame* fGroupOptRun;
        TGLabel* fLabelOptAutoRun;
        TGNumberEntry* fNEntryOAutoRun;
        TGLabel* fLabelOptInc;
        TGNumberEntry* fNEntryOInc;
        
        TGGroupFrame* fGroupTelFrame;
        TGHorizontalFrame* fGroupAnaFrame;
        TGHorizontalFrame* fGroupSetFrame;
        // frame with display radio buttons
        TGButtonGroup* fGroupOptDis;
        TGRadioButton* fRadioB1;
        TGRadioButton* fRadioB2;
        TGRadioButton* fRadioB3;
        TGRadioButton* fRadioB4;
        // frame with telescope radio buttons
        TGButtonGroup* fGroupOptTel;
        TGRadioButton* fRadioTA;
        vector< TGRadioButton* > fRadioTel;
        
        // frame with cut optiones
        TGGroupFrame* fGroupOptCut;
        TGLabel* fLabelOptATri;
        TGNumberEntry* fNEntryOATri;
        TGLabel* fLabelOptAIma;
        TGNumberEntry* fNEntryOAIma;
        TGLabel* fLabelOptCuts;
        TGTextEntry* fTextOptCuts;
        
        TGTextButton* fButtonOptSet;
        TGTextButton* fButtonOptReset;
        TGCompositeFrame* fFrameInfo;
        TGCompositeFrame* fCompInfo;
        TRootEmbeddedCanvas* fEmInfo;
        TCanvas* fCanvasInfo;
        TPaveText* fPaveInfo;
        
        E_cameraIdent fCameraDisplay;             //!< which camera display is active
        bool fBoolDrawImageTraces;                //!< draw individual image traces, not sum signal into FADC canvas
        unsigned int fSelectedChan;               //!< selected FADC channel
        E_fadcIDENT fAnaDisplay;                  //! select FADC/analysis  tab
        TH1D* fHisFADC;                           //!< FADC histogram
        string fHisFADCDrawString;
        TF1* fF1Ped;                              //!< FADC pedestal
        TGraph* fGraphFADC;                       //!< graph to indicate summation window
        TGraph* fGraphFADC_2;                       //!< graph to indicate summation window (second summation window)
        TLine* fLineFADC;                         //!< line to indicate Tzero
        TLine* fTALineFADC;
        bool fBoolFADC;                           //!< draw FADC
        
        int fBWNum;                               //! Number of colours in BW palette
        int fBWPalette[50];                       //! BW (greyscale) palette
        
        bool fCameraMovie;                        //!< plot every new camera view
        string fMovieFileName;                    //!< name of gifs for camera movie + number
        unsigned int fMoviePictNumber;            //!< movie picture number
        
        unsigned int fNumEventIncrement;          //!< increment for nextEvent()
        bool fAutoRunStatus;                      //!< true = autorunmodus is on
        unsigned int fTimingSleep;                //!< pause between each event in autorunmodus (microseconds)
        bool fCameraTiming;                       //!< last tab was timing tab (read again data)
        
        void     bookHistos();                    //!< book histograms
        //!< not all plots make sense
        bool     checkPlotIntentions( unsigned int iCamTab );
        void     defineGui();                     //!< construct all the widgets
        void     drawFADC( bool iFit );           //!< draw FADC histogram
        void     drawCalibrationHistos();         //!< draw calibration histos
        bool     drawTgradGraphs();               //!< draw timing graphs
        bool     drawImageBorderCharge();         //!< draw image border graphs
        bool     drawImageBorderTZero();          //!< draw image border graphs
        void     drawPixelHistos();               //!< draw pixel histos
        void     dumpDeadChannels();              //!< dump dead channels to screen
        void     dumpImageBorderPixels();         //!< pring image/border pixels
        TH1D*    fillFADC( int i_channel, TH1D* i_his,  bool bReset = true );          //!< fill FADC histogram
        void     makeMoviePicture();              //!< make a new camera movie picture
        void     plotFADCFit( int );              //!< plot FADC trace with fit function
        void     printCanvas( TPad* );            //!< print canvas in eps,ps,gif,root
        void     processEvent();
        //!< process input of buttons, menues, etc.
        virtual Bool_t ProcessMessage( Long_t msg, Long_t parm1, Long_t );
        void     resetDisplay();                  //!< reset display
        void     resetRunOptions();               //!< reset run variables to standard values
        void     searchChannel( int );            //!< search channel and plot the FADC trace
        void     setCameraPads( bool iFieldView );//!< define pads position and size in case of more than one telescope
        void     setColorScheme();                //!< switch on/off color scheme for all cameras
        void     setFADCText();                   //!< set information for channel (number, hit status, ped, etc.)
        void     setInfoText();                   //!< set text with general infos on analysis parameters
        void     showWarning();                   //!< plot a TText with warning that functionality is not available
        void     subprocessCheckButton( Long_t ); //!< process check button interaction
        void     subprocessComboBox( Long_t );    //!< process combo box interactions
        void     subprocessButton( Long_t );      //!< process button interaction
        void     subprocessMenu( Long_t );        //!< process menu interaction
        void     subprocessRadioButton( Long_t ); //!< process radio buttion interactions
        void     subprocessTextChanged( Long_t ); //!< process text changed interactions
        void     subprocessTextEnter( Long_t );   //!< process text entered interactions
        
    public:
        VDisplay();
        //!< standard constructor, sourceFile = source data file
        VDisplay( const TGWindow* p, unsigned int w, unsigned int h, VEventLoop* iEventLoop );
        virtual ~VDisplay();                      //!< destructor
        
        //  void setDisplayTraceFit(bool infit){fTraceFit=infit;}
        
        // slots
        void     CloseWindow();                   //!< close application (TGMainFrame method)
        void     selectAnaTab( Int_t );           //!< select tab for FADC/analysis
        //!< select channel for FADC histogram by clicking on channel
        void     selectChannel( Int_t, Int_t, Int_t, TObject* );
        void     updateCamera( Int_t );           //!< update camera view with new event
        
#ifndef __APPLE__
        ClassDef( VDisplay, 2 )
#endif
};
#endif                                            // VRDISPLAY_H
