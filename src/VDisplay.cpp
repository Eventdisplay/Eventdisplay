/*!  \class VDisplay
  \brief  GUI class for event display

  \bug middle mouse button deletes object in canvas
  \bug nextEvent with selection: bad expression in selection -> this event is skipped
  \bug reading of gain filename, seg fault if directory = calibration/gain


*/

#include "VDisplay.h"

VDisplay::VDisplay() : TGMainFrame( 0, 0, 0 )
{
}


//! standard constructor, sourceFile = source data file
/*!
  \param p              client
  \param iEventLoop     pointer to parent eventloop
  \param h              window height
  \param w              window width
*/
VDisplay::VDisplay( const TGWindow* p, unsigned int h, unsigned int w, VEventLoop* iEventLoop ) : TGMainFrame( p, h, w )
{
    fDebug = iEventLoop->getDebugFlag();
    // set the right color palettes
    TStyle* plain  = new TStyle( "Plain", "Plain Style (no colors/fill areas)" );
    plain->SetCanvasBorderMode( 0 );
    plain->SetPadBorderMode( 0 );
    plain->SetPadColor( 0 );
    plain->SetCanvasColor( 0 );
    plain->SetTitleColor( 0 );
    plain->SetStatColor( 0 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetPalette( 1 );
    
    fEventLoop = iEventLoop;
    
    // set telescope number
    fTelescope = fEventLoop->getTelID();          // telescope ID
    
    // setting the camera display classes
    for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
    {
        fTelescopesToShow.push_back( fEventLoop->getTeltoAna()[i] );
        fCamera[fEventLoop->getTeltoAna()[i]] = new VCamera( fEventLoop->getTeltoAna()[i], iEventLoop->getData() );
    }
    // drawing mode in the beging is each telescope in seperate canvas
    fBoolDrawOne = false;
    fBoolDrawAllinOne = false;
    fBool_M_OPT_COL_SCHE_Checked = true;
    fBool_M_OPT_BW_SCHE_Checked = false;
    // birds view
    fBirdsEye = new VDisplayBirdsEye();
    if( fEventLoop->getRunParameter()->fPlotPaper )
    {
        fBirdsEye->setPlotPaper();
    }
    fBirdsEye->setNTel( fEventLoop->getTeltoAna().size() );
    fBirdsEye->setData( ( VEvndispData* )fEventLoop );
    // fadc canvas
    fBoolFADC = false;
    fBoolDrawImageTraces = false;
    
    fSelectedChan = 0;
    fCameraTiming = false;
    fCameraMovie = false;
    fMoviePictNumber = 0;
    fAutoRunStatus = false;
    
    fCanvasesNx = 1;
    fCanvasesNy = 1;
    fBWNum = 50;
    int iOff = 0;
    int iCol = 8800;
    for( int i = 0; i < fBWNum; i++ )
    {
        TColor* color = new TColor( iCol + i, 1 - float( i + iOff ) / ( fBWNum + iOff ),
                                    1 - float( i + iOff ) / ( fBWNum + iOff ), 1 - float( i + iOff ) / ( fBWNum + iOff ), "" );
        color->GetNumber();
        fBWPalette[i] = iCol + i;
    }
    // standard is to display events with telescope trigger only
    if( fEventLoop->getNTel() == 1 )
    {
        fEventLoop->setCutNArrayTrigger( 1 );
    }
    else
    {
        fEventLoop->setCutNArrayTrigger( 2 );
    }
    fEventLoop->setCutNArrayImages( 0 );
    
    // beginning of gui definitions
    defineGui();
    setCameraPads( false );
    updateCamera( fCameraDisplay );
    
    // signals (GM)
    //   to be consistent, should be replace by Associate() mechanismus
    //    but: how to read out mouse position on canvas?
    //  leave it like it is, because signal/slot mechanism is more powerful
    
    fTabAna->Connect( "Selected(Int_t )", "VDisplay", this, "selectAnaTab( Int_t )" );
    fCanvasCamera->Connect( "ProcessedEvent(Int_t , Int_t , Int_t , TObject* )", "VDisplay", this, "selectChannel( Int_t, Int_t, Int_t, TObject* )" );
    fCanvasCamera->SetHighLightColor( 10 );
    fCanvasCamera->SetEditable( false );
    
    char i_wname[300];
    sprintf( i_wname, "%s event analysis and display %s (%s) %s",
             fEventLoop->getRunParameter()->getObservatory().c_str(),
             fEventLoop->getEventDisplayVersion().c_str(),
             fEventLoop->getRunParameter()->getSVN_VERSION().c_str(),
             fEventLoop->getRunParameter()->fRunTitle.c_str() );
    SetWindowName( i_wname );
    MapSubwindows();
    Layout();
    // layout not really good -> window not resizable
    SetMWMHints( kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize | kMWMDecorMinimize | kMWMDecorMenu, kMWMFuncAll |  kMWMFuncResize    | kMWMFuncMaximize | kMWMFuncMinimize, kMWMInputModeless );
    MapWindow();
    // end of gui
    
    resetDisplay();
    setInfoText();
    bookHistos();
}


void VDisplay::CloseWindow()
{
    fEventLoop->shutdown();
    gApplication->Terminate( 0 );
}


/*!
   these histograms are only needed for plotting purpuse (like FADC trace)
*/
void VDisplay::bookHistos()
{
    // FADC histogram settings
    //   number of samples is adjusted later during the filling
    fHisFADC = new TH1D( "fHisFADC", "Channel 0", 64, 0., 64. );
    fHisFADC->SetLineColor( 9 );
    fHisFADC->SetLineWidth( 2 );
    fHisFADC->SetXTitle( "sample number" );
    fHisFADCDrawString = "hist";
    fF1Ped = new TF1( "ped", "[0]", fHisFADC->GetXaxis()->GetXmin(), fHisFADC->GetXaxis()->GetXmax() );
    fF1Ped->SetLineColor( 44 );
    fF1Ped->SetLineWidth( 2 );
    fF1Ped->SetLineStyle( 2 );
    fGraphFADC = new TGraph( 4 );
    fGraphFADC->SetFillColor( 38 );
    fGraphFADC_2 = new TGraph( 4 );
    fLineFADC = new TLine();
    fLineFADC->SetLineColor( 1 );
    fLineFADC->SetLineStyle( 2 );
    fLineFADC->SetLineWidth( 2 );
    fTALineFADC = new TLine();
    fTALineFADC->SetLineColor( 9 );
    fTALineFADC->SetLineStyle( 3 );
    fTALineFADC->SetLineWidth( 2 );
}


/*!
  clean everything

  \attention{ many deletes still miss}
*/
VDisplay::~VDisplay()
{
    Cleanup();
}


/*!
   \param priCanvas  point to canvas to be printed

*/
void VDisplay::printCanvas( TPad* priCanvas )
{

    const Char_t* filetypes[] =
    {
        "PDF files",    "*.pdf",
        "EPS files",    "*.eps",
        "PS files",    "*.ps",
        "GIF files",    "*.gif",
        "ROOT files",    "*.root",
        "ROOT macros",   "*.C",
        "All files",     "*",
        0,               0
    };
    
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    new TGFileDialog( fClient->GetRoot(), this, kFDSave, &fi );
    if( !fi.fFilename )
    {
        return;
    }
    
    priCanvas->Print( fi.fFilename );
    
    priCanvas->Update();
}


/*!
 */
Bool_t VDisplay::ProcessMessage( Long_t msg, Long_t parm1, Long_t )
{
    switch( GET_MSG( msg ) )
    {
        case kC_COMMAND:
            switch( GET_SUBMSG( msg ) )
            {
                case kCM_MENU:
                    subprocessMenu( parm1 );
                    break;
                case kCM_CHECKBUTTON:
                    subprocessCheckButton( parm1 );
                    break;
                case kCM_BUTTON:
                    subprocessButton( parm1 );
                    break;
                case kCM_RADIOBUTTON:
                    subprocessRadioButton( parm1 );
                    break;
                case kCM_COMBOBOX:
                    subprocessComboBox( parm1 );
                    break;
                default:
                    break;
            }
        case kC_TEXTENTRY:
            switch( GET_SUBMSG( msg ) )
            {
                case kTE_ENTER:
                    subprocessTextEnter( parm1 );
                    break;
                case kTE_TEXTCHANGED:
                    subprocessTextChanged( parm1 );
                    break;
                default:
                    break;
            }
        default:
            break;
    }
    return true;
}


/*!
  call this if
  - a new event should be displayed
  - display options changed
  - camera tab changed

  \param i camera tab identifier (e.g. charge/trigger/hit/..., see E_cameraIdent)
*/
void VDisplay::updateCamera( Int_t i )
{
    if( fDebug )
    {
        cout << "VDisplay::updateCamera " << i << "\t" << fCamera.size() << "\t" << fTelescopesToShow.size() << endl;
    }
    
    fCanvasCamera->SetEditable( true );
    // get tab identification
    fCameraDisplay = E_cameraIdent( i );
    if( fBoolDrawAllinOne )
    {
        fCanvasCamera->Clear();
    }
    // check if it makes sense to plot something
    if( !checkPlotIntentions( i ) )
    {
        return;
    }
    
    // ==================================================
    // timing tab
    if( E_cameraIdent( i ) == C_TIMING && !fBoolDrawAllinOne &&
            fEventLoop->getRunParameter()->fsourcetype != 4 && fEventLoop->getRunParameter()->fsourcetype != 7 )
    {
        // don't do anything for first event
        if( fEventLoop->getEventNumber() == 0 )
        {
            return;
        }
        
        // can't do timing tab and all in one camera
        if( fBoolDrawAllinOne )
        {
            return;
        }
        
        fCameraTiming = true;
        unsigned int winsize = 4;
        
        // get maximum value in any tube
        vector< double > i_max;
        for( unsigned t = 0; t < fTelescopesToShow.size(); t++ )
        {
            fEventLoop->getAnalyzer()->setTelID( fTelescopesToShow[t] );
            // (GM) only need maximum sum, don't call doAnalysis() (resets all image parameters)
            // fEventLoop->getAnalyzer()->doAnalysis();
            fEventLoop->getAnalyzer()->calcSums( ( int )fEventLoop->getRunParameter()->fsumfirst[t],
                                                 ( int )( fEventLoop->getRunParameter()->fsumfirst[t] + fEventLoop->getRunParameter()->fsumwindow_1[t] ),
                                                 false, false, 1 );
            i_max.push_back( fEventLoop->getAnalyzer()->getSums().max() );
        }
        // loop over the trace in steps of winsize
        // for long samples only loop over first half
        int i_subS = 0;
        if( fEventLoop->getNSamples() == 64 )
        {
            i_subS = 28;
        }
        for( unsigned int j = 0; j < fEventLoop->getNSamples() - winsize - i_subS; j++ )
        {
            for( unsigned t = 0; t < fTelescopesToShow.size(); t++ )
            {
                fEventLoop->getAnalyzer()->setTelID( fTelescopesToShow[t] );
                fEventLoop->getAnalyzer()->calcSums( ( int )j, ( int )j + winsize, false, false, 1 );
                fEventLoop->getAnalyzer()->gainCorrect();
                if( fEventLoop->getAnalyzer()->getImageCleaner() )
                {
                    fEventLoop->getAnalyzer()->getImageCleaner()->cleanImageFixed( 30.0, 15.0 );
                }
                if( fBoolDrawOne )
                {
                    fCamera[fTelescopesToShow[t]]->setCanvas( fCanvasCamera );
                }
                else
                {
                    fCamera[fTelescopesToShow[t]]->setCanvas( fPadsCamera[fTelescopesToShow[t]] );
                }
                
                fCamera[fTelescopesToShow[t]]->setMode( i );
                fCamera[fTelescopesToShow[t]]->setCurrentTimeSlice( j );
                fCamera[fTelescopesToShow[t]]->draw( i_max[t], fEventLoop->getEventNumber(), fBoolDrawAllinOne );
                fCanvasCamera->Update();
                
                gSystem->ProcessEvents();
                
                if( !checkPlotIntentions( i ) )
                {
                    return;
                }
                
                if( fCameraMovie )
                {
                    makeMoviePicture();
                }
                
                if( !fCameraTiming )
                {
                    break;
                }
                
                fCanvasCamera->Update();
                usleep( fTimingSleep );
            }
            if( !fCameraTiming )
            {
                break;
            }
        }
    }
    // =========================================================
    // all other tabs except timing tab
    else if( E_cameraIdent( i ) != C_TIMING )
    {
        // if previous tab was the timing tab, recalculate image parameters
        if( fCameraTiming )
        {
            fCameraTiming = false;
            for( unsigned int t = 0; t < fTelescopesToShow.size(); t++ )
            {
                fEventLoop->getAnalyzer()->setTelID( fTelescopesToShow[t] );
                fEventLoop->getAnalyzer()->doAnalysis();
            }
            //  updateCamera( fCameraDisplay );
            if( fTabAna->GetCurrent() == 1 )
            {
                drawFADC( false );
            }
            gSystem->ProcessEvents();
            if( !checkPlotIntentions( i ) )
            {
                return;
            }
        }
        // loop over all telescopes we want to see
        for( unsigned int t = 0; t < fTelescopesToShow.size(); t++ )
        {
            fEventLoop->getAnalyzer()->setTelID( fTelescopesToShow[t] );
            if( t == 0 )
            {
                fCamera[fTelescopesToShow[t]]->setFirstTelescopeToDraw();
            }
            if( fEventLoop->getAnalyzer()->getImageParameters() )
            {
                // draw one camera into the camera canvas
                if( fBoolDrawOne || fBoolDrawAllinOne )
                {
                    fCamera[fTelescopesToShow[t]]->setCanvas( fCanvasCamera );
                }
                // draw more than on camera into the camera canvas
                else
                {
                    fCamera[fTelescopesToShow[t]]->setCanvas( fPadsCamera[fTelescopesToShow[t]] );
                }
                fCamera[fTelescopesToShow[t]]->setMode( i );
                // draw the camera
                fCamera[fTelescopesToShow[t]]->draw( 0., fEventLoop->getEventNumber(), fBoolDrawAllinOne );
            }
        }
        // plot the movie
        if( fCameraMovie )
        {
            makeMoviePicture();
        }
        fCanvasCamera->Update();
    }
    // draw the calibration histos depending on camera
    if( fEventLoop->getEventNumber() != 0 &&  fTabAna->GetCurrent() == 2 )
    {
        drawCalibrationHistos();
    }
    
    // draw the timing graphs
    if( fEventLoop->getEventNumber() != 0 && fTabAna->GetCurrent() == 3 )
    {
        drawPixelHistos();
    }
    
    fCanvasCamera->SetEditable( false );
    if( fDebug )
    {
        cout << "END VDisplay::updateCamera " << i << "\t" << fCamera.size() << endl;
    }
}


/*!
   file names of movies are NAME_01.gif, NAME_02.gif, ...
*/
void VDisplay::makeMoviePicture()
{
    string suffix;
    char i_Temp[500];
    fMoviePictNumber++;
    
    if( fMoviePictNumber % fEventLoop->getNTel() == 0 || fCameraDisplay != C_TIMING )
    {
        // can be changed to any graphical output
        sprintf( i_Temp, "_%.5d.gif", fMoviePictNumber );
        suffix = i_Temp;
        // no printout Info in <TCanvas::Print>: ....
        gErrorIgnoreLevel = 1;
        fCanvasCamera->Print( ( fMovieFileName + suffix ).c_str() );
    }
}


/*!
   these are the tabs on the right side of the display

   \param it tab number
*/
void VDisplay::selectAnaTab( Int_t it )
{
    if( fDebug )
    {
        cout << "VDisplay::selectAnaTab" << endl;
    }
    fAnaDisplay = E_fadcIDENT( it );
    // FADC tab
    if( it == 1 )
    {
        drawFADC( false );
        if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
        {
            fButtonFADCset->SetState( kButtonUp );
            fButtonFADCunset->SetState( kButtonUp );
            fButtonFADCreset->SetState( kButtonUp );
            fButtonFADCFit->SetState( kButtonUp );
        }
        else
        {
            fButtonFADCset->SetState( kButtonDisabled );
            fButtonFADCunset->SetState( kButtonDisabled );
            fButtonFADCreset->SetState( kButtonDisabled );
            fButtonFADCFit->SetState( kButtonDisabled );
        }
    }
    // calibration tab
    else if( it == 2 )
    {
        drawCalibrationHistos();
    }
    // birds eye tab
    else if( it == 6 )
    {
        fBirdsEye->setData( fEventLoop );
        fBirdsEye->draw( fCanvasBird );
    }
    // tgrad tab
    else if( it == 3 )
    {
        drawPixelHistos();
    }
}


/*!
  \attention
  fSelectedChan should always be >= 200000 (convention)

  \param bt mouse button
  \param x  mouse position x
  \param y  mouse position y
  \param objSel pointer to selected object

*/
void VDisplay::selectChannel( Int_t bt, Int_t x, Int_t y , TObject* objSel )
{
    if( fEventLoop->getEventNumber() < 1 )
    {
        return;
    }
    if( fBoolDrawAllinOne )
    {
        return;
    }
    if( bt == kButton1Up )
    {
        // get pad and set current telescope number to this pad
        for( unsigned int j = 0; j < fPadsCamera.size(); j++ )
        {
            if( objSel == fPadsCamera[j] && j < fEventLoop->getNTel() )
            {
                fTelescope = j;
                if( fDebug )
                {
                    cout << "VDisplay::selectChannel, selected telescope: " << fTelescope + 1 << "\t" << fCamera.size() << endl;
                }
            }
        }
        // update tgrad pad?
        if( fTabAna->GetCurrent() == 3 )
        {
            //drawPixelHistos();
        }
        // get tube
        unsigned int iSelect = 0;
        fSelectedChan = 0;
        for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
        {
            iSelect = fCamera[fTelescopesToShow[i]]->getChannel( x, y, objSel );
            if( iSelect >= 200000 && iSelect < 200000 + fEventLoop->getAnalyzer()->getImage().size() + 1 )
            {
                fSelectedChan = iSelect;
                fTelescope = fTelescopesToShow[i];
                if( fDebug )
                {
                    cout << "VDisplay::selectChannel, selected tube: " << fSelectedChan - 200000 << " (telescope " << fTelescope + 1 << ")" << endl;
                }
                break;
            }
        }
        // hit a tube
        if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() + 1 )
        {
            fBoolFADC = true;
            if( fTabAna->GetCurrent() == 2 )
            {
                drawCalibrationHistos();
            }
            else if( fTabAna->GetCurrent() == 3 )
            {
                drawPixelHistos();
            }
            else
            {
                if( fTabAna->GetCurrent() != 1 )
                {
                    fTabAna->SetTab( 1 );
                }
                
                drawFADC( false );
                fButtonFADCset->SetState( kButtonUp );
                fButtonFADCunset->SetState( kButtonUp );
                fButtonFADCreset->SetState( kButtonUp );
                fButtonFADCFit->SetState( kButtonUp );
            }
        }
        // hit elsewhere
        else
        {
            fSelectedChan = 0;
            fBoolFADC = false;
            fCanvasFADC->Clear();
            if( fTabAna->GetCurrent() == 1 )
            {
                drawFADC( false );
                fButtonFADCset->SetState( kButtonDisabled );
                fButtonFADCunset->SetState( kButtonDisabled );
                fButtonFADCreset->SetState( kButtonDisabled );
                fButtonFADCFit->SetState( kButtonDisabled );
            }
            else if( fTabAna->GetCurrent() == 2 )
            {
                drawCalibrationHistos();
            }
            else if( fTabAna->GetCurrent() == 3 )
            {
                drawPixelHistos();
            }
        }
    }
    else if( bt == kButton2Up )                   // middle mouse button
    {
        // nothing happens
    }
}


/*!
   \param i_channel channel number
*/
void VDisplay::searchChannel( int i_channel )
{
    fSelectedChan = 200000 + i_channel;
    if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
    {
        fBoolFADC = true;
        // if current tab is the FADC tab, draw FADC trace of selected channel
        if( fTabAna->GetCurrent() == 1 )
        {
            drawFADC( false );
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->showSelectedChannel( i_channel, true );
            }
            fButtonFADCset->SetState( kButtonUp );
            fButtonFADCunset->SetState( kButtonUp );
            fButtonFADCreset->SetState( kButtonUp );
            fButtonFADCFit->SetState( kButtonUp );
        }
    }
}


/*!
 *
 * fill a histogram with the FADCS trace for the given channel

   \param i_channel channel number
   \param i_his pointer to FADC histogram

*/
TH1D* VDisplay::fillFADC( int i_channel, TH1D* i_his, bool bReset )
{
    if( !i_his )
    {
        return 0;
    }
    if( bReset )
    {
        i_his->Reset();
    }
    ///////////////////////////////////////////////////////////////////////////////
    // fill trace into histogram
    ///////////////////////////////////////////////////////////////////////////////
    if( fEventLoop->hasFADCData()
            && fEventLoop->getZeroSuppressed()[i_channel] < 2
            && fEventLoop->getNSamples() > 0 )
    {
        // first set the number of bins of the histogram according to the number of samples
        int i_nbins = ( int )fEventLoop->getNSamples();
        if( fEventLoop->getDigitalFilterMethod() > 0 )
        {
            i_nbins *= fEventLoop->getDigitalFilterUpSample();
        }
        if( i_nbins != i_his->GetNbinsX() )
        {
            i_his->SetBins( i_nbins, 0., double( i_nbins ) );
        }
        pair< bool, uint32_t > i_hitIndexPair = fEventLoop->getReader()->getChannelHitIndex( i_channel );
        if( !i_hitIndexPair.first )
        {
            for( int i = 1; i <= i_his->GetNbinsX(); i++ )
            {
                i_his->SetBinContent( i, 0. );
            }
        }
        else
        {
            // use trace handler
            // (allows to take digital filtering into account)
            VTraceHandler* i_TraceHandler = fEventLoop->getTraceHandler();
            if( !i_TraceHandler )
            {
                return 0;
            }
            i_TraceHandler->setDigitalFilterParameters( fEventLoop->getDigitalFilterMethod(),
                    fEventLoop->getDigitalFilterUpSample(),
                    fEventLoop->getDigitalFilterPoleZero() );
            i_TraceHandler->setTrace( fEventLoop->getReader(),
                                      fEventLoop->getNSamples(),
                                      fEventLoop->getPeds( fEventLoop->getHiLo()[i_channel] )[i_channel],
                                      i_channel, i_hitIndexPair.second,
                                      fEventLoop->getHiLo()[i_channel]*fEventLoop->getLowGainMultiplier_Trace() );
                                      
            vector< double > i_trace = i_TraceHandler->getTrace();
            for( unsigned int t = 0; t < i_trace.size(); t++ )
            {
                if( bReset )
                {
                    i_his->SetBinContent( t + 1, -1.*i_trace[t] );
                }
                else
                {
                    i_his->SetBinContent( t + 1, i_his->GetBinContent( t + 1 ) - i_trace[t] );
                }
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////////
    // DST: fill timing values into histogram
    ///////////////////////////////////////////////////////////////////////////////
    else if( !fEventLoop->hasFADCData()
             && fEventLoop->getZeroSuppressed()[i_channel] != 3
             && fEventLoop->getZeroSuppressed()[i_channel] != 1 )
    {
        // first set the number of bins of the histogram according to the number of timing bins
        // number of bins is: number of pulse time levels + one bin at the beginning and one at the end of the trace
        //const int nbinsX = fEventLoop->getRunParameter()->fpulsetiminglevels.size() + 2 + 1;
        const unsigned int nbinsX_max = 10000;
        double xbins[nbinsX_max];
        double y[nbinsX_max];
        unsigned int nbinsX = 0;
        // determinate number of valid bins
        for( unsigned int i = 0; i < fEventLoop->getRunParameter()->fpulsetiminglevels.size(); i++ )
        {
            if( fEventLoop->getPulseTiming()[i][i_channel] > 1.e-3 )
            {
                nbinsX++;
            }
        }
        nbinsX += 2;
        nbinsX += 1;
        if( nbinsX > nbinsX_max )
        {
            return 0;
        }
        
        // set pedestal value
        for( unsigned int i = 0; i < nbinsX; i++ )
        {
            y[i] = fEventLoop->getPeds()[i_channel];
        }
        
        // first bin is at time zero
        xbins[0] = 0.;
        unsigned int i_maxPV = fEventLoop->getRunParameter()->fpulsetiming_max_index;
        // set pulse times
        unsigned int z = 1;
        for( unsigned int i = 0; i < fEventLoop->getRunParameter()->fpulsetiminglevels.size(); i++ )
        {
            if( fEventLoop->getPulseTiming()[i][i_channel] > 1.e-3 )
            {
                if( i <= i_maxPV )
                {
                    xbins[z] = fEventLoop->getPulseTiming()[i][i_channel];
                }
                else
                {
                    xbins[z] = fEventLoop->getPulseTiming()[fEventLoop->getRunParameter()->fpulsetiminglevels.size() - i - 1][i_channel]
                               + fEventLoop->getPulseTiming()[i][i_channel];
                }
                
                y[i + 1] += fEventLoop->getRunParameter()->fpulsetiminglevels[i] * fEventLoop->getTraceMax()[i_channel];
                z++;
            }
        }
        for( unsigned int i = z; i < nbinsX; i++ )
        {
            xbins[i] = xbins[i - 1] + 2.;
        }
        // check that all bins are fine
        bool bBroken = false;
        for( unsigned int i = 1; i < nbinsX; i++ )
        {
            if( xbins[i] < xbins[i - 1] )
            {
                bBroken = true;
            }
        }
        if( !bBroken )
        {
            i_his->SetBins( nbinsX - 1, xbins );
            for( int i = 1; i <= i_his->GetNbinsX(); i++ )
            {
                i_his->SetBinContent( i, -1.*y[i - 1] );
            }
        }
        else
        {
            i_his->SetBins( 1, 0., 1. );
            i_his->SetBinContent( 1, 0. );
        }
    }
    
    return i_his;
}


/*!

     draw the FADC trace

     for option tracefit, draw as well the result of the fit

*/
void VDisplay::drawFADC( bool iFit )
{
    if( fDebug )
    {
        cout << "VDisplay::drawFADC " << iFit << "\t" << fCamera.size() <<  endl;
    }
    // don't draw anything before first event
    if( fEventLoop->getEventNumber() < 1 )
    {
        return;
    }
    // don't draw FADC trace in PE mode, write channel information only
    if( fEventLoop->getReader()->getDataFormatNum() == 6 )
    {
        // write only channels infos into FADC window
        if( fSelectedChan >= 200000 )
        {
            setFADCText();
        }
        return;
    }
    
    // from here on: draw fHisFADC (at least)
    fHisFADC->Reset();
    char histitle[200];
    fHisFADC->SetLineColor( 9 );
    if( TMath::Abs( fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope )
                    - TMath::Floor( fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) ) ) < 1.e-3 )
    {
        sprintf( histitle, "sample number [%dns]", ( int )fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) );
    }
    else
    {
        sprintf( histitle, "sample number [%.1fns]", fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) );
    }
    fHisFADC->GetXaxis()->SetTitle( histitle );
	fHisFADC->GetXaxis()->SetTitleColor( 1 );
    
    TH1D* iTraceFits;
    iTraceFits = 0;
    
    // make fadc canvas editable
    fCanvasFADC->SetEditable( 1 );
    fCanvasFADC->cd();
    
    fEventLoop->getAnalyzer()->setTelID( fTelescope );
    // reset max pixel number in 'show channel' number entry (might change from telescope to telescope)
    fNEntryFADCsearch->SetLimits( TGNumberFormat::kNELLimitMinMax, 0, fEventLoop->getAnalyzer()->getImage().size() - 1 );
    // plot trace of one channel (click on channel)
    if( fSelectedChan >= 200000 && fEventLoop->getZeroSuppressed()[fSelectedChan - 200000] < 2 )
    {
        sprintf( histitle, "Channel #%d (Telescope %d)", fSelectedChan - 200000, fTelescope + 1 );
        
        // fill histogram with fadc trace
        fHisFADC = fillFADC( fSelectedChan - 200000, fHisFADC );
        fHisFADC->SetLineStyle( 1 );
        fHisFADC->SetLineColor( 1 );
        fHisFADC->SetMarkerColor( 1 );
        fHisFADC->SetLineWidth( 2 );
        if( fSelectedChan - 200000 < fEventLoop->getHiLo().size() )
        {
            // mark overflow channels
            if( fEventLoop->getHiLo()[fSelectedChan - 200000] )
            {
                fHisFADC->SetLineWidth( 4 );
                fHisFADC->SetLineStyle( 5 );
            }
            // set histogram colors according to signal/border/background
            if( fEventLoop->getAnalyzer()->getImage()[fSelectedChan - 200000] )
            {
                fHisFADC->SetLineColor( 2 );
            }
            else if( fEventLoop->getAnalyzer()->getBorder()[fSelectedChan - 200000] )
            {
                fHisFADC->SetLineColor( 8 );
            }
            else if( fEventLoop->getAnalyzer()->getDead( fEventLoop->getHiLo()[fSelectedChan - 200000] )[fSelectedChan - 200000] )
            {
                fHisFADC->SetLineColor( 14 );
            }
            else
            {
                fHisFADC->SetLineColor( 9 );
            }
            
            fF1Ped->SetParameter( 0, -1.* fEventLoop->getAnalyzer()->getPeds( fEventLoop->getHiLo()[fSelectedChan - 200000] )[fSelectedChan - 200000] );
        }
        ///////////////////////////////////////////////////////////
        // draw everything (trace, pedestal)
        setFADCText();
        fHisFADC->SetTitle( histitle );
        if( fEventLoop->getHiLo()[fSelectedChan - 200000] )
        {
            fHisFADC->SetMaximum( -1.*0.1 * fEventLoop->getPed_min( fEventLoop->getHiLo()[fSelectedChan - 200000] ) );
        }
        else
        {
            fHisFADC->SetMaximum( -1.*0.8*fEventLoop->getPed_min( fEventLoop->getHiLo()[fSelectedChan - 200000] ) );
        }
        //fHis <<FADC->SetMaximum( -1111 );
        fHisFADC->SetStats( 0 );
        fCanvasFADC->SetEditable( 1 );
        fCanvasFADC->cd();
        fHisFADC->SetMarkerSize( 2 ) ;
        fHisFADC->Draw( fHisFADCDrawString.c_str() );
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // plot all image signals into one canvas, click beside camera and use switch in 'option' menu
    else if( fBoolDrawImageTraces && !( fSelectedChan >= 200000 ) )
    {
        bool i_FADCdrawn = false;
        sprintf( histitle, "image traces (Telescope %d)", fTelescope + 1 );
        // get maximum of all image traces
        double i_traceMax = 0.;
        for( unsigned int i = 0; i < fEventLoop->getAnalyzer()->getImage().size(); i++ )
        {
            // first get maximum signal in signal traces
            if( fEventLoop->getAnalyzer()->getImage()[i] )
            {
                fHisFADC->Reset();
                if( fEventLoop->getZeroSuppressed()[i] < 2 )
                {
                    fHisFADC = fillFADC( i, fHisFADC );
                }
                if( fHisFADC->GetMinimum() < i_traceMax )
                {
                    i_traceMax = fHisFADC->GetMinimum();
                }
            }
        }
        // now plot all traces (color loop)
        int i_TraceColor = 1;
        int i_TraceStyle = 1;
        for( unsigned int i = 0; i < fEventLoop->getAnalyzer()->getImage().size(); i++ )
        {
            if( fEventLoop->getAnalyzer()->getImage()[i] )
            {
                fHisFADC->Reset();
                if( fEventLoop->getZeroSuppressed()[i] < 2 )
                {
                    fHisFADC = fillFADC( i, fHisFADC );
                }
                fHisFADC->SetMinimum( i_traceMax );
                fHisFADC->SetMaximum( -1111 );
                fHisFADC->SetStats( 0 );
                fHisFADC->SetLineStyle( i_TraceStyle );
                fHisFADC->SetLineWidth( 2 );
                fHisFADC->SetLineColor( i_TraceColor );
                fHisFADC->SetMarkerColor( i_TraceColor );
                fCanvasFADC->SetEditable( 1 );
                fCanvasFADC->cd();
                if( !i_FADCdrawn )
                {
                    fHisFADC->SetTitle( histitle );
                    fHisFADC->DrawCopy( fHisFADCDrawString.c_str() );
                    i_FADCdrawn = true;
                }
                else
                {
                    string iTemp = fHisFADCDrawString + " same";
                    fHisFADC->DrawCopy( iTemp.c_str() );
                }
                i_TraceColor++;
                if( i_TraceColor > 20 )
                {
                    i_TraceColor = 1;
                    if( i_TraceStyle < 7 )
                    {
                        i_TraceStyle++;
                    }
                    else
                    {
                        i_TraceStyle = 1;
                    }
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // plot sum signal (sum of all image pixels), click beside camera for that
    else if( !( fSelectedChan >= 200000 ) )
    {
        sprintf( histitle, "average signal (Telescope %d)", fTelescope + 1 );
        int i_image = 0;
        for( unsigned int i = 0; i < fEventLoop->getAnalyzer()->getImage().size(); i++ )
        {
            if( fEventLoop->getAnalyzer()->getImage()[i] )
            {
                if( fEventLoop->getZeroSuppressed()[i] < 2 )
                {
                    bool bReset = ( i_image == 0 );	//Reset trace histogram iff we are at the first image pixel.
                    fHisFADC = fillFADC( i, fHisFADC, bReset );
                }
                i_image++;
            }
        }
        if( i_image > 0 )
        {
            fHisFADC->Scale( 1. / double( i_image ) );
        }
        fCanvasFADCText->SetEditable( 1 );
        fCanvasFADCText->Clear();
        fCanvasFADCText->Update();
        fCanvasFADCText->SetEditable( 0 );
        fHisFADC->SetLineWidth( 2 );
        fHisFADC->SetLineStyle( 1 );
        fHisFADC->SetLineColor( 1 );
        fHisFADC->SetMarkerColor( 1 );
        fHisFADC->SetTitle( histitle );
        fHisFADC->SetMinimum( -1111 );
        fHisFADC->SetMaximum( -1111 );
        fHisFADC->SetStats( 0 );
        fCanvasFADC->SetEditable( 1 );
        fCanvasFADC->cd();
        fHisFADC->Draw( fHisFADCDrawString.c_str() );
    }
    fCanvasFADC->Update();
    
    if( fSelectedChan >= 200000 && fSelectedChan - 200000 < fEventLoop->getAnalyzer()->getTCorrectedSumFirst().size() )
    {
        // plot box which indicates integration (summation window) window
        if( fEventLoop->hasFADCData() )
        {
            // plot box for second integration window
            fGraphFADC_2->SetPoint( 0, ( double )fEventLoop->getAnalyzer()->getTCorrectedSum2First()[fSelectedChan - 200000], gPad->GetUymin() );
            fGraphFADC_2->SetPoint( 1, ( double )fEventLoop->getAnalyzer()->getTCorrectedSum2First()[fSelectedChan - 200000], gPad->GetUymax() );
            fGraphFADC_2->SetPoint( 2, ( double )fEventLoop->getAnalyzer()->getTCorrectedSum2Last()[fSelectedChan - 200000], gPad->GetUymax() );
            fGraphFADC_2->SetPoint( 3, ( double )fEventLoop->getAnalyzer()->getTCorrectedSum2Last()[fSelectedChan - 200000], gPad->GetUymin() );
            fGraphFADC_2->SetFillColor( 41 );
            fGraphFADC_2->SetFillStyle( 3005 );
            // plot box for first integration window
            fGraphFADC->SetPoint( 0, ( double )fEventLoop->getAnalyzer()->getTCorrectedSumFirst()[fSelectedChan - 200000], gPad->GetUymin() );
            fGraphFADC->SetPoint( 1, ( double )fEventLoop->getAnalyzer()->getTCorrectedSumFirst()[fSelectedChan - 200000], gPad->GetUymax() );
            fGraphFADC->SetPoint( 2, ( double )fEventLoop->getAnalyzer()->getTCorrectedSumLast()[fSelectedChan - 200000], gPad->GetUymax() );
            fGraphFADC->SetPoint( 3, ( double )fEventLoop->getAnalyzer()->getTCorrectedSumLast()[fSelectedChan - 200000], gPad->GetUymin() );
            fGraphFADC->Draw( "f" );
            fGraphFADC_2->Draw( "f" );
        }
        // draw a line to indicate TZero
        if( fEventLoop->getAnalyzer()->getRawTZeros()[fSelectedChan - 200000] > 0. )
        {
            fLineFADC->SetX1( fEventLoop->getAnalyzer()->getRawTZeros()[fSelectedChan - 200000] );
            fLineFADC->SetX2( fEventLoop->getAnalyzer()->getRawTZeros()[fSelectedChan - 200000] );
            fLineFADC->SetY1( gPad->GetUymax() );
            fLineFADC->SetY2( gPad->GetUymin() );
            fLineFADC->Draw( "same" );
        }
        // draw a line to indicate average trace time
        if( fEventLoop->getAnalyzer()->getTraceAverageTime( false ).size() > 0
                && fEventLoop->getAnalyzer()->getTraceAverageTime( false )[fSelectedChan - 200000] > 0. )
        {
            fTALineFADC->SetX1( fEventLoop->getAnalyzer()->getTraceAverageTime( false )[fSelectedChan - 200000] );
            fTALineFADC->SetX2( fEventLoop->getAnalyzer()->getTraceAverageTime( false )[fSelectedChan - 200000] );
            fTALineFADC->SetY1( gPad->GetUymax() );
            fTALineFADC->SetY2( gPad->GetUymin() );
            fTALineFADC->Draw( "same" );
        }
        // trace line indicating pedestals
        fF1Ped->SetRange( fHisFADC->GetXaxis()->GetXmax(), fHisFADC->GetXaxis()->GetXmin() );
        fF1Ped->SetParameter( 0, -1.* fEventLoop->getAnalyzer()->getPeds( fEventLoop->getHiLo()[fSelectedChan - 200000] )[fSelectedChan - 200000] );
        if( !fBoolDrawImageTraces || fSelectedChan >= 200000 )
        {
            string itemp = fHisFADCDrawString + " same";
            fHisFADC->Draw( itemp.c_str() );
        }
        // change color of low gain pedestal is used
        if( fEventLoop->getHiLo()[fSelectedChan - 200000] )
        {
            fF1Ped->SetLineColor( 50 );
        }
        else
        {
            fF1Ped->SetLineColor( 44 );
        }
        // for individual pulses, draw fit function
        fF1Ped->Draw( "same" );
        if( iTraceFits )
        {
            iTraceFits->Draw( "lsame" );
        }
    }
    if( fSelectedChan >= 200000
            && fEventLoop->getZeroSuppressed()[fSelectedChan - 200000] > 1 )
    {
        setFADCText();
        sprintf( histitle, "Channel #%d (Telescope %d)", fSelectedChan - 200000, fTelescope + 1 );
        fHisFADC->SetTitle( histitle );
        fHisFADC->SetStats( 0 );
        fCanvasFADC->SetEditable( 1 );
        fCanvasFADC->cd();
        fHisFADC->Draw();
        string i_textZS = "ZERO SUPPRESSED";
        if( fEventLoop->getZeroSuppressed()[fSelectedChan - 200000] == 2 )
        {
            i_textZS = "SAMPLE SUPPRESSED";
        }
        TText* iT = new TText( 2, 0.5, i_textZS.c_str() );
        iT->Draw();
        fCanvasFADC->Update();
    }
    fCanvasFADC->RedrawAxis();
    
    fCanvasFADC->Update();
}


/*!

 */
void VDisplay::processEvent()
{
    if( fDebug )
    {
        cout << "VDisplay::processEvent() " << endl;
    }
    char c_ev[200];
    
    do
    {
        for( unsigned int i = 0; i < fNumEventIncrement; i++ )
        {
            // stop if there is no next event
            if( !fEventLoop->nextEvent() )
            {
                if( fEventLoop->getReader()->getEventStatus() > 998 )
                {
                    return;
                }
                if( fEventLoop->getReader()->getEventStatus() != fEventLoop->getExpectedEventStatus() )
                {
                    fEventLoop->incrementNumberofIncompleteEvents();
                    i--;
                    continue;
                }
                return;
            }
            // draw everything if event number is incremented by fNumEventIncrement
            if( ( i + 1 ) % fNumEventIncrement == 0 )
            {
                updateCamera( fCameraDisplay );
                if( fTabAna->GetCurrent() == 1 && !fCameraTiming )
                {
                    drawFADC( false );
                }
                else if( fTabAna->GetCurrent() == 3 )
                {
                    drawPixelHistos();
                }
                else if( fTabAna->GetCurrent() == 6 )
                {
                    fBirdsEye->draw( fCanvasBird );
                }
                fEventLoop->incrementNumberofGoodEvents();
            }
            // look for input in auto runstatus
            if( fAutoRunStatus )
            {
                gSystem->ProcessEvents();
            }
            sprintf( c_ev, "now at event %d", fEventLoop->getEventNumber() );
            fStatusBar->SetText( c_ev, 1 );
            if( !fAutoRunStatus )
            {
                fNEntryGoto->SetNumber( fEventLoop->getEventNumber() );
            }
        }
    }
    while( fAutoRunStatus );
}


void VDisplay::resetRunOptions()
{
    fTimingSleep =  10000;                        // in ms
    fNumEventIncrement = 1;
}


void VDisplay::resetDisplay()
{
    resetRunOptions();
    fEventLoop->resetRunOptions();
    
    fSelectedChan = 0;
    fCameraTiming = false;
    fBoolFADC = false;
    
    fCanvasCamera->SetEditable( true );
    //     fCanvasCamera->Clear();
    fCanvasCamera->SetEditable( false );
    
    fCanvasFADC->Clear();
}


/*!
    calculate number of pads from number of cameras
    calculate pad size, etc.

    \par iFieldView true if pads are distributed like telescopes in the field
*/
void VDisplay::setCameraPads( bool iFieldView )
{
    fPadsCamera.clear();
    char i_text[300];
    double iPadXlow, iPadYlow, iPadXup, iPadYup;
    
    // each pad beside another
    if( !iFieldView )
    {
        int iNpadX, iNpadY;
        double iPadXstep, iPadYstep;
        
        // calculate number of pads in each direction
        iNpadX = ( int )sqrt( ( double )fEventLoop->getTeltoAna().size() );
        iNpadY = ( fEventLoop->getTeltoAna().size() + iNpadX - 1 ) / iNpadX;
        if( iNpadX != iNpadY )
        {
            int mini = min( iNpadX, iNpadY ) + 1;
            if( mini * mini >= ( int )fEventLoop->getTeltoAna().size() )
            {
                iNpadX = mini;
                iNpadY = mini;
            }
        }
        
        iPadXstep = 1. / ( double )iNpadX;
        iPadYstep = 1. / ( double )iNpadY;
        
        // clear the canvas camera
        fCanvasCamera->SetEditable( true );
        fCanvasCamera->Clear();
        fCanvasCamera->cd();
        unsigned int z = 0;
        for( int iy = 0; iy < iNpadY; iy++ )
        {
            iPadYlow = 1. - ( iy + 1 ) * iPadYstep;
            for( int ix = 0; ix < iNpadX; ix++ )
            {
                // calculate pad size
                iPadXlow = ix * iPadXstep;
                iPadXup  = ( ix + 1 ) * iPadXstep;
                iPadYup = 1. - iy * iPadYstep;
                // set pad size name
                sprintf( i_text, "%s_tel%d", "pad", ( int )fPadsCamera.size() );
                // create new pad
                if( z < fEventLoop->getTeltoAna().size() )
                {
                    fPadsCamera[fEventLoop->getTeltoAna()[z]] = new TPad( i_text, i_text, iPadXlow, iPadYlow, iPadXup, iPadYup );
                    fPadsCamera[fEventLoop->getTeltoAna()[z]]->Draw();
                    fPadsCamera[fEventLoop->getTeltoAna()[z]]->SetEditable( false );
                }
                z++;
                fCanvasCamera->Update();
            }
        }
    }
    // cameras are distributed like telescopes in the field
    // this works only for showers with ze=0 deg !!!!
    else
    {
        // get largest distance in x and y between telescopes
        // assume that one telescope is at (0.,0.)
        // this the size of one pad and the scale factor
        // for the calculation of the pads position
        double idistmax = 0.;
        double idistmin = 10000.;
        // normalisation factor for convertion [m]->padsize
        double inorm = 1.;
        // setoff to avoid overlapping pads (min x-dist between pads is 69m)
        double iOffset = 0.;
        // loop over all telescopes
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            if( fabs( fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[i]] ) > idistmax )
            {
                idistmax = fabs( fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[i]] );
            }
            if( fabs( fEventLoop->getDetectorGeo()->getTelYpos()[fEventLoop->getTeltoAna()[i]] ) > idistmax )
            {
                idistmax = fabs( fEventLoop->getDetectorGeo()->getTelYpos()[fEventLoop->getTeltoAna()[i]] );
            }
            if( fabs( fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[i]] ) < idistmin && fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[i]] != 0. )
            {
                idistmin = fabs( fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[i]] );
            }
        }
        // in this case, don't know what to do
        if( idistmax == 0. || idistmin == 0. )
        {
            cout << "VDisplay::setCameraPads error: don't get the pad positions right, check telescope positions" << endl;
            return;
        }
        // 3 pads beside each other, maxdist is 1xpadsize
        iOffset = 1. / 3 - idistmin / idistmax * 1. / 3.;
        iOffset /= 2.;
        inorm = 1. / 3. / idistmax;
        unsigned int z = 0;
        // now set the pads for all telescopes
        fCanvasCamera->SetEditable( true );
        fCanvasCamera->Clear();
        fCanvasCamera->cd();
        for( unsigned int j = 0; j < fEventLoop->getTeltoAna().size(); j++ )
        {
            iPadXlow =  fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[j]] * inorm - 1. / 6 + 0.5 + iOffset ;
            iPadYlow =  fEventLoop->getDetectorGeo()->getTelYpos()[fEventLoop->getTeltoAna()[j]] * inorm - 1. / 6 + 0.5 + iOffset;
            iPadXup =  fEventLoop->getDetectorGeo()->getTelXpos()[fEventLoop->getTeltoAna()[j]] * inorm + 1. / 6 + 0.5 - iOffset;
            iPadYup =  fEventLoop->getDetectorGeo()->getTelYpos()[fEventLoop->getTeltoAna()[j]] * inorm + 1. / 6 + 0.5 - iOffset;
            sprintf( i_text, "%s_tel%d", "pad", ( int )fPadsCamera.size() );
            if( z < fEventLoop->getTeltoAna().size() )
            {
                fPadsCamera[fEventLoop->getTeltoAna()[z]] = new TPad( i_text, i_text, iPadXlow, iPadYlow, iPadXup, iPadYup );
                fPadsCamera[fEventLoop->getTeltoAna()[z]]->Draw();
            }
            fCanvasCamera->Update();
        }
        fCanvasCamera->SetEditable( false );
    }
}

/*

   text printed below FADC trace

*/
void VDisplay::setFADCText()
{
    if( fDebug )
    {
        cout << "VDisplay::setFADCText " << fSelectedChan << endl;
    }
    if( fSelectedChan < 200000 )
    {
        return;
    }
    char cTemp[500];
    unsigned int iChannel = unsigned( fSelectedChan - 200000 );
    
    // don't use this for photodiode
    if( iChannel >= fEventLoop->getAnalyzer()->getCurrentSumWindow().size() )
    {
        return;
    }
    
    fCanvasFADCText->SetEditable( true );
    
    // text placement
    float xL = 0.02;
    float yT = 0.92;
    float ystep = 0.04 * 1.5;
    float textsize = 0.035 * 1.5;
    
    fCanvasFADCText->Clear();
    fTextFADC.clear();
    
    fEventLoop->getData()->setTelID( fTelescope );
    // channel number
    TString iFADCtext = "";
    if( fEventLoop->getDBPixelDataReader() && fEventLoop->getDBPixelDataReader()->getDBStatus() )
    {
        int iFADCmodule = fEventLoop->getDBPixelDataReader()->getFADC_module( fTelescope, fSelectedChan - 200000 );
        int iFADCchannel = fEventLoop->getDBPixelDataReader()->getFADC_channel( fTelescope, fSelectedChan - 200000 );
        iFADCtext.Form( ", FADC %d/%d", iFADCmodule, iFADCchannel );
    }
    ostringstream iSTRText;
    iSTRText << "telescope " << fTelescope + 1;
    iSTRText << "channel " << fSelectedChan - 200000 << iFADCtext.Data();
    iSTRText << "(NN " << fEventLoop->getDetectorGeometry()->getNNeighbours()[iChannel] << ": ";
    for( unsigned int n = 0; n < fEventLoop->getDetectorGeometry()->getNNeighbours()[iChannel]; n++ )
    {
        iSTRText << " " << fEventLoop->getDetectorGeometry()->getNeighbours()[iChannel][n];
    }
    iSTRText << ")";
    fTextFADC.push_back( new TText( xL, yT, iSTRText.str().c_str() ) );
    // L1/HV/currents (if available)
    if( fEventLoop->getDBPixelDataReader() && fEventLoop->getDBPixelDataReader()->getDBStatus() )
    {
        sprintf( cTemp, "L1 rate %.2e Hz, HV %d V, I %.2f muA", fEventLoop->getL1Rate( fSelectedChan - 200000 ),
                 ( int )fEventLoop->getHV( fSelectedChan - 200000 ),
                 fEventLoop->getCurrent( fSelectedChan - 200000 ) );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    else
    {
        fTextFADC.push_back( new TText( xL, yT, "" ) );
    }
    // integration window size
    sprintf( cTemp, "integration window: from sample %d to %d", fEventLoop->getAnalyzer()->getTCorrectedSumFirst()[iChannel], fEventLoop->getAnalyzer()->getTCorrectedSumLast()[iChannel] );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // print # of PE if available
    if( fEventLoop->getPE().size() > 0 && iChannel < fEventLoop->getPE().size() )
    {
        sprintf( cTemp, "Pe: %d", ( int )fEventLoop->getPE()[iChannel] );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    // pedestal
    sprintf( cTemp, "pedestal %.2f (low gain: %.2f)", fEventLoop->getAnalyzer()->getPeds()[iChannel], fEventLoop->getAnalyzer()->getPeds( true )[iChannel] );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // pedestal variance
    if( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] > 0 )
    {
        sprintf( cTemp, "ped var %.2f (low gain: %.2f), 1st window %d, LG mult %.2f",
                 fEventLoop->getAnalyzer()->getPedvars( false, fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] )[iChannel],
                 fEventLoop->getAnalyzer()->getPedvars( true, fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] )[iChannel],
                 fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] ,
                 fEventLoop->getLowGainMultiplier_Sum( fEventLoop->getRunParameter()->fTraceIntegrationMethod[ fEventLoop->getAnalyzer()->getTelID() ] ,
                         fEventLoop->getRunParameter()->fsumwindow_1[ fEventLoop->getAnalyzer()->getTelID() ] , fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] ) );
    }
    else if( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] == 0 )
    {
        sprintf( cTemp, "no pedestal variance for 0 summation window" );
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    if( fEventLoop->getAnalyzer()->getSumWindow_2() > 0 )
    {
        sprintf( cTemp, "ped var %.2f (low gain: %.2f), 2nd window %d, LG mult %.2f",
                 fEventLoop->getAnalyzer()->getPedvars( false, fEventLoop->getAnalyzer()->getCurrentSumWindow_2()[iChannel] )[iChannel],
                 fEventLoop->getAnalyzer()->getPedvars( true, fEventLoop->getAnalyzer()->getCurrentSumWindow_2()[iChannel] )[iChannel],
                 fEventLoop->getAnalyzer()->getCurrentSumWindow_2()[iChannel] ,
                 fEventLoop->getLowGainMultiplier_Sum( fEventLoop->getRunParameter()->fTraceIntegrationMethod[ fEventLoop->getAnalyzer()->getTelID() ]  ,
                         fEventLoop->getRunParameter()->fsumwindow_2[ fEventLoop->getAnalyzer()->getTelID() ]  , fEventLoop->getAnalyzer()->getCurrentSumWindow_2()[iChannel] ) );
    }
    else if( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] == 0 )
    {
        sprintf( cTemp, "no pedestal variance for 0 summation window" );
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    if( fEventLoop->getAnalyzer()->isDoublePass() && fEventLoop->getAnalyzer()->getSumWindow_Pass1() > 0 )
    {
        int iSW = fEventLoop->getAnalyzer()->getSumWindow_Pass1();
        sprintf( cTemp, "ped var %.2f (low gain: %.2f), DP1 window %d, LG mult %.2f",
                 fEventLoop->getAnalyzer()->getPedvars( false, iSW )[iChannel],
                 fEventLoop->getAnalyzer()->getPedvars( true, iSW )[iChannel],
                 iSW  ,
                 fEventLoop->getLowGainMultiplier_Sum( fEventLoop->getRunParameter()->fTraceIntegrationMethod_pass1[ fEventLoop->getAnalyzer()->getTelID() ],
                         fEventLoop->getRunParameter()->fsumwindow_pass1[ fEventLoop->getAnalyzer()->getTelID() ] , iSW ) );
    }
    else if( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] == 0 )
    {
        sprintf( cTemp, "no pedestal variance for 0 summation window" );
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // gain
    sprintf( cTemp, "gain %.2f (low gain channel: %.2f)", fEventLoop->getAnalyzer()->getGains()[iChannel], fEventLoop->getAnalyzer()->getGains( true )[iChannel] );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    if( fEventLoop->getAnalyzer()->getHIGHQE_gainfactor( iChannel ) > 0. )
    {
        sprintf( cTemp, "HIGHQE gain correction: %.2f", fEventLoop->getAnalyzer()->getHIGHQE_gainfactor( iChannel ) );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    // time offset
    //    sprintf( cTemp, "toffset %.2f (low gain channel: %.2f, FADC stop: %.2f)", fEventLoop->getAnalyzer()->getTOffsets()[iChannel], fEventLoop->getAnalyzer()->getTOffsets(true)[iChannel], fEventLoop->getAnalyzer()->getFADCStopOffsets()[iChannel]);
    //   fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // L2 FADC stop
    if( !fEventLoop->isDST_MC() )
    {
        sprintf( cTemp, "L2 channels FADC stop: %.2f %.2f %.2f %.2f",
                 fEventLoop->getAnalyzer()->getFADCstopTZero()[0],
                 fEventLoop->getAnalyzer()->getFADCstopTZero()[1],
                 fEventLoop->getAnalyzer()->getFADCstopTZero()[2],
                 fEventLoop->getAnalyzer()->getFADCstopTZero()[3] );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    // low gain multiplier
    sprintf( cTemp, "low gain multiplier (trace): %.2f", fEventLoop->getAnalyzer()->getLowGainMultiplier_Trace() );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    
    // pulse sum
    sprintf( cTemp, "pulse sum %.1f (2ndWi: %.1f) pulse max %.1f (raw max: %.1f) pulse width: %.1f", fEventLoop->getData()->getSums()[iChannel],
             fEventLoop->getData()->getSums2()[iChannel],
             fEventLoop->getData()->getTraceMax()[iChannel],
             fEventLoop->getData()->getTraceRawMax()[iChannel],
             fEventLoop->getData()->getTraceWidth()[iChannel] );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // pulse tzeros
    sprintf( cTemp, "T0 %.1f (average t: %.1f, toffset correction: %.1f, FADC crate trigger: %.1f)",
             fEventLoop->getPulseTime()[iChannel],
             fEventLoop->getTraceAverageTime()[iChannel],
             fEventLoop->getAnalyzer()->getTOffsets()[iChannel],
             fEventLoop->getAnalyzer()->getFADCStopOffsets()[iChannel] );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // pulse timing
    sprintf( cTemp, "pulse timing (raw): " );
    for( unsigned int p = 0; p < fEventLoop->getRunParameter()->fpulsetiminglevels.size(); p++ )
    {
        sprintf( cTemp, "%s %d%% : %.1f ", cTemp, ( int )( fEventLoop->getRunParameter()->fpulsetiminglevels[p] * 100. ), fEventLoop->getPulseTiming( false )[p][iChannel] );
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    // sum / pedvar
    double i_var = 0.;
    if( fEventLoop->getRunParameter()->fsourcetype != 6 &&
            iChannel < fEventLoop->getAnalyzer()->getCurrentSumWindow().size() )
    {
        unsigned int iSumWindow = fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel];
        if( iChannel < fEventLoop->getAnalyzer()->getPedvars( iSumWindow ).size() && iChannel < fEventLoop->getAnalyzer()->getGains().size() )
        {
            if( iChannel < fEventLoop->getAnalyzer()->getPedvars( iSumWindow, true ).size() )
            {
                if( fEventLoop->getAnalyzer()->getPedvars( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel], fEventLoop->getHiLo()[iChannel] )[iChannel]  > 0. )
                {
                    i_var = fEventLoop->getAnalyzer()->getSums()[iChannel] * fEventLoop->getAnalyzer()->getGains()[iChannel]
                            / fEventLoop->getAnalyzer()->getPedvars( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel], fEventLoop->getHiLo()[iChannel] )[iChannel];
                }
            }
        }
    }
    char iClusterString[200] = "";
    if( fEventLoop->getData()->getClusterID().size() > iChannel && fEventLoop->getData()->getClusterID()[iChannel] != 0 )
    {
        sprintf( iClusterString, ", cluster ID: %d",  fEventLoop->getData()->getClusterID()[iChannel] ) ;
    }
    sprintf( cTemp, "pulse sum/pedestal variation (integration window %d): %.2f%s", fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel], i_var, iClusterString );
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    if( fEventLoop->getData()->getImage()[iChannel] )
    {
        fTextFADC.back()->SetTextColor( 2 );
    }
    else if( fEventLoop->getData()->getBorder()[iChannel] )
    {
        fTextFADC.back()->SetTextColor( 8 );
    }
    else
    {
        fTextFADC.back()->SetTextColor( 1 );
    }
    
    // Template Expectation Value
    
    if( fEventLoop->getRunParameter()->ffrogsmode == 1 )
    {
        sprintf( cTemp, "Mu %.2f (%.2f) ImgGood %.2f BkgGood %.2f  (%.2f) Frogs Energy %.2f", fEventLoop->getData()->getTemplateMu()[iChannel], 5.3 * fEventLoop->getData()->getTemplateMu()[iChannel], fEventLoop->getData()->getFrogsParameters()->frogsGoodnessImg, fEventLoop->getData()->getFrogsParameters()->frogsGoodnessBkg, pow( fEventLoop->getAnalyzer()->getSums()[iChannel] / fEventLoop->getAnalyzer()->getPedvars( fEventLoop->getAnalyzer()->getCurrentSumWindow()[iChannel] )[iChannel], 2.0 ) - 1.0, fEventLoop->getData()->getFrogsParameters()->frogsEnergy );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    
    // Model3D Value
    
    if( fEventLoop->getRunParameter()->fUseDisplayModel3D )
    {
        sprintf( cTemp, "Model3D: %.2f", fEventLoop->getData()->getModel3DMu()[iChannel] );
        fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    }
    
    // dead channel text
    
    //fTextFADC.push_back( new TText( xL, yT, "" ) );
    sprintf( cTemp, "  " );
    bitset<15> idead = fEventLoop->getDead()[iChannel];
    if( idead.test( 11 ) )
    {
        sprintf( cTemp, "high gain: %s (%s)", fEventLoop->getDeadChannelText()[11].c_str(), idead.to_string<char, char_traits<char>, allocator<char> >().c_str() );
    }
    else if( idead.test( 9 ) )
    {
        sprintf( cTemp, "high gain: %s (%s)", fEventLoop->getDeadChannelText()[9].c_str(), idead.to_string<char, char_traits<char>, allocator<char> >().c_str() );
    }
    else
    {
        for( unsigned int s = 0; s < idead.size(); s++ )
        {
            if( idead.test( s ) )
            {
                if( s < fEventLoop->getDeadChannelText().size() )
                {
                    sprintf( cTemp, "high gain: %s (%s)", fEventLoop->getDeadChannelText()[s].c_str(), idead.to_string<char, char_traits<char>, allocator<char> >().c_str() );
                }
                break;
            }
        }
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    sprintf( cTemp, "  " );
    idead = fEventLoop->getDead( true )[iChannel];
    if( idead.test( 11 ) )
    {
        sprintf( cTemp, "low gain: %s (%s)", fEventLoop->getDeadChannelText()[11].c_str(), idead.to_string<char, char_traits<char>, allocator<char> >().c_str() );
    }
    else
    {
        for( unsigned int s = 0; s < idead.size(); s++ )
        {
            if( idead.test( s ) )
            {
                if( s < fEventLoop->getDeadChannelText().size() )
                {
                    sprintf( cTemp, "low gain: %s (%s)", fEventLoop->getDeadChannelText()[s].c_str(), idead.to_string<char, char_traits<char>, allocator<char> >().c_str() );
                }
                break;
            }
        }
    }
    fTextFADC.push_back( new TText( xL, yT, cTemp ) );
    
    for( unsigned int i = 0; i < fTextFADC.size(); i++ )
    {
        fTextFADC[i]->SetTextSize( textsize );
        fTextFADC[i]->SetTextFont( 42 );
        fTextFADC[i]->DrawTextNDC( xL, yT, fTextFADC[i]->GetTitle() );
        yT -= ystep;
    }
    fCanvasFADCText->Update();
    fCanvasFADCText->SetEditable( false );
}


/*!
    display basic run infos at 'information' tab

*/
void VDisplay::setInfoText()
{
    fCanvasInfo->SetEditable( true );
    vector<TText*> i_Text;
    char c_Text[300];
    
    float xL = 0.02;
    float yT = 0.96;
    float ystep = 0.03;
    float textsize = 0.030;
    
    fCanvasInfo->Clear();
    
    sprintf( c_Text, "run number: %d", fEventLoop->getRunNumber() );
    i_Text.push_back( new TText( xL, yT, c_Text ) );
    i_Text.back()->SetTextSize( i_Text.back()->GetTextSize() * 1.1 );
    i_Text.push_back( new TText( xL, yT, "" ) );
    sprintf( c_Text, "data file: %s ", fEventLoop->getDataFileName().c_str() );
    i_Text.push_back( new TText( xL, yT, c_Text ) );
    i_Text.push_back( new TText( xL, yT, "" ) );
    i_Text.back()->SetTextSize( i_Text.back()->GetTextSize() * 1.1 );
    
    // run info
    sprintf( c_Text, "Target: %s (ra,dec)=(%.2f,%.2f), wobble E %.2f, N %.2f", fEventLoop->getRunParameter()->fTargetName.c_str(),
             fEventLoop->getRunParameter()->fTargetRA,
             fEventLoop->getRunParameter()->fTargetDec,
             fEventLoop->getRunParameter()->fWobbleEast,
             fEventLoop->getRunParameter()->fWobbleNorth );
    i_Text.push_back( new TText( xL, yT, c_Text ) );
    i_Text.push_back( new TText( xL, yT, "" ) );
    
    // calibration and analysis parameters
    // (do this only if there are less than 10 telescopes)
    if( fEventLoop->getTeltoAna().size() < 10. )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            fEventLoop->getAnalyzer()->setTelID( fEventLoop->getTeltoAna()[i] );
            sprintf( c_Text, "Telescope %d:", fEventLoop->getTeltoAna()[i] + 1 );
            i_Text.push_back( new TText( xL, yT, c_Text ) );
            if( fEventLoop->getReader()->getDataFormat() == "rawdata"  ||  fEventLoop->getReader()->getDataFormat() == "Rawvbf" )
            {
                TString ilpedString;
                if( fEventLoop->getRunParameter()->fPedLowGainFile != "" )
                {
                    ilpedString = fEventLoop->getRunParameter()->fPedLowGainFile;
                }
                else
                {
                    ilpedString.Form( "%d", fEventLoop->getRunParameter()->fPedLowGainFileNumber[fEventLoop->getTeltoAna()[i]] );
                }
                sprintf( c_Text, "    pedestal (low): %d (%s), rel.gain %d time offset %d ",
                         fEventLoop->getRunParameter()->fPedFileNumber[fEventLoop->getTeltoAna()[i]],
                         ilpedString.Data(),
                         fEventLoop->getRunParameter()->fGainFileNumber[fEventLoop->getTeltoAna()[i]],
                         fEventLoop->getRunParameter()->fTOffFileNumber[fEventLoop->getTeltoAna()[i]] );
                i_Text.push_back( new TText( xL, yT, c_Text ) );
            }
            if( fEventLoop->getAnalyzer()->getImageCleaningParameter() )
            {
                sprintf( c_Text, "    image/border threshold %.2f/%.2f", fEventLoop->getAnalyzer()->getImageCleaningParameter()->fimagethresh,
                         fEventLoop->getAnalyzer()->getImageCleaningParameter()->fborderthresh );
                i_Text.push_back( new TText( xL, yT, c_Text ) );
            }
            sprintf( c_Text, "    window size %d (%d), window start (not time corrected) %d",
                     fEventLoop->getAnalyzer()->getSumWindow(), fEventLoop->getAnalyzer()->getSumWindow_2(),
                     fEventLoop->getAnalyzer()->getSumFirst() );
            i_Text.push_back( new TText( xL, yT, c_Text ) );
        }
    }
    
    // draw everything
    for( unsigned int i = 0; i < i_Text.size(); i++ )
    {
        i_Text[i]->SetTextSize( textsize );
        i_Text[i]->SetTextFont( 42 );
        i_Text[i]->DrawTextNDC( xL, yT, i_Text[i]->GetTitle() );
        yT -= ystep;
    }
    fCanvasInfo->Update();
    fCanvasInfo->SetEditable( false );
}


void VDisplay::drawCalibrationHistos()
{
    if( fDebug )
    {
        cout << "VDisplay::drawCalibrationHistos()" << endl;
    }
    
    // is there a channel selected
    unsigned int iChannel = 0;
    if( fSelectedChan < 200000 )
    {
        iChannel = 0;
    }
    else
    {
        iChannel = fSelectedChan - 200000;
    }
    
    double iMeanDistributionValue = -99.;
    
    TH1F* ihis = 0;
    TH1F* ihis2 = 0;
    fCanvasCal->cd();
    // decide what to draw, peds, gains, or toffs
    if( E_cameraIdent( fCameraDisplay ) == C_GAINS || E_cameraIdent( fCameraDisplay ) == C_GAINSLOW
            || E_cameraIdent( fCameraDisplay ) == C_GAINVARS || E_cameraIdent( fCameraDisplay ) == C_GAINVARSLOW )
    {
        if( fSelectedChan >= 200000 )
        {
            ihis = fEventLoop->getCalData( fTelescope )->getHistoGain( fTelescope, iChannel );
            ihis2 = fEventLoop->getCalData( fTelescope )->getHistoGain( fTelescope, iChannel, true );
        }
        else
        {
            ihis = fEventLoop->getCalData( fTelescope )->getGainDist();
            ihis2  = fEventLoop->getCalData( fTelescope )->getGainDist( true );
        }
        if( ihis )
        {
            ihis->SetAxisRange( 0., 4. );
        }
    }
    else if( E_cameraIdent( fCameraDisplay ) == C_TOFF || E_cameraIdent( fCameraDisplay ) == C_TOFFLOW )
    {
        if( fSelectedChan >= 200000 )
        {
            ihis = fEventLoop->getCalData( fTelescope )->getHistoToff( fTelescope, iChannel );
            ihis2 = fEventLoop->getCalData( fTelescope )->getHistoToff( fTelescope, iChannel, true );
        }
        else
        {
            ihis = fEventLoop->getCalData( fTelescope )->getToffsetDist();
            ihis2 = fEventLoop->getCalData( fTelescope )->getToffsetDist( true );
        }
    }
    else if( E_cameraIdent( fCameraDisplay ) == C_PEDMEAN || E_cameraIdent( fCameraDisplay ) == C_PEDMEANLOW )
    {
        if( fSelectedChan >= 200000 )
        {
            ihis = fEventLoop->getCalData( fTelescope )->getHistoPed( fTelescope, iChannel, fEventLoop->getRunParameter()->fsumwindow_1[fTelescope],
                    false, fEventLoop->getTelType( fTelescope ) );
            ihis2 = fEventLoop->getCalData( fTelescope )->getHistoPed( fTelescope, iChannel, fEventLoop->getRunParameter()->fsumwindow_1[fTelescope],
                    true, fEventLoop->getTelType( fTelescope ) );
        }
        else
        {
            ihis = fEventLoop->getCalData( fTelescope )->getPedDist();
            ihis2  = fEventLoop->getCalData( fTelescope )->getPedDist( true );
        }
        if( ihis )
        {
            ihis->SetAxisRange( 0., 250. );
        }
    }
    // pedvar distributions
    else if( E_cameraIdent( fCameraDisplay ) == C_PEDVAR  || E_cameraIdent( fCameraDisplay ) == C_PEDVARLOW )
    {
        if( fSelectedChan >= 200000 )
        {
            ihis = fEventLoop->getCalData( fTelescope )->getHistoPed( fTelescope, iChannel, fEventLoop->getRunParameter()->fsumwindow_1[fTelescope],
                    false, fEventLoop->getTelType( fTelescope ) );
            ihis2 = fEventLoop->getCalData( fTelescope )->getHistoPed( fTelescope, iChannel, fEventLoop->getRunParameter()->fsumwindow_1[fTelescope],
                    true, fEventLoop->getTelType( fTelescope ) );
        }
        else
        {
            ihis = fEventLoop->getCalData( fTelescope )->getPedvarsDist();
            ihis2  = fEventLoop->getCalData( fTelescope )->getPedvarsDist( true );
        }
    }
    // tzero distributions
    else if( E_cameraIdent( fCameraDisplay ) == C_CALTZERO || E_cameraIdent( fCameraDisplay ) == C_CALTZEROLOW )
    {
        if( fSelectedChan >= 200000 )
        {
            ihis = fEventLoop->getCalData( fTelescope )->getHistoAverageTzero( fTelescope, iChannel );
            ihis2 = fEventLoop->getCalData( fTelescope )->getHistoAverageTzero( fTelescope, iChannel, true );
            fEventLoop->setTelID( fTelescope );
            iMeanDistributionValue = fEventLoop->getAverageTZeros()[iChannel];
        }
        else
        {
            ihis = fEventLoop->getCalData( fTelescope )->getAverageTzerosetDist();
            ihis2 = 0;
            //            ihis2 = fEventLoop->getCalData( fTelescope )->getAverageTzerosetDist( true );
            fEventLoop->setTelID( fTelescope );
            iMeanDistributionValue = fEventLoop->getMeanAverageTZero();
        }
    }
    // low gain multiplier
    else if( E_cameraIdent( fCameraDisplay ) == C_LOWGAIN )
    {
        ihis = fEventLoop->getCalData( fTelescope )->getLowGainMultiplierDistribution();
        if( ihis )
        {
            ihis->SetAxisRange( 0., 10. );
        }
    }
    // L1 rates
    else if( E_cameraIdent( fCameraDisplay ) == C_L1 )
    {
        if( fEventLoop->getDBPixelDataReader() )
        {
            ihis = fEventLoop->getDBPixelDataReader()->getL1Histogram( fTelescope, fEventLoop->getEventMJD(),
                    fEventLoop->getEventTime() );
        }
    }
    // HV
    else if( E_cameraIdent( fCameraDisplay ) == C_HV )
    {
        if( fEventLoop->getDBPixelDataReader() )
        {
            ihis = fEventLoop->getDBPixelDataReader()->getHVHistogram( fTelescope, fEventLoop->getEventMJD(),
                    fEventLoop->getEventTime() );
        }
    }
    // Currents
    else if( E_cameraIdent( fCameraDisplay ) == C_CURRENTS )
    {
        if( fEventLoop->getDBPixelDataReader() )
        {
            ihis = fEventLoop->getDBPixelDataReader()->getCurrentsHistogram( fTelescope, fEventLoop->getEventMJD(),
                    fEventLoop->getEventTime() );
        }
    }
    
    ////////////////////////////////////////////
    // plot everything
    if( ihis )
    {
        ihis->SetLineWidth( 2 );
        ihis->Draw();
        if( ihis2 && ihis2->GetEntries() > 0 )
        {
            ihis2->SetLineWidth( 2 );
            ihis2->SetLineColor( 2 );
            ihis2->Draw( "sames" );
        }
        if( iMeanDistributionValue > -98. )
        {
            TLine* iL = new TLine( iMeanDistributionValue, 0., iMeanDistributionValue, ihis->GetMaximum() );
            iL->SetLineStyle( 2 );
            iL->Draw();
        }
    }
    else
    {
        fCanvasCal->Clear();
    }
    
    fCanvasCal->Update();
}

/*

   draw histograms into the corresponding tab

*/
void VDisplay::drawPixelHistos()
{

    fCanvasPixelHisto->Clear();
    fCanvasPixelHisto->Divide( 1, 3 );
    
    TPad* iP = ( TPad* )fCanvasPixelHisto->cd( 3 );
    if( !drawTgradGraphs() )
    {
        iP->Clear();
    }
    
    iP = ( TPad* )fCanvasPixelHisto->cd( 1 );
    if( !drawImageBorderCharge() )
    {
        iP->Clear();
    }
    
    iP = ( TPad* )fCanvasPixelHisto->cd( 2 );
    if( !drawImageBorderTZero() )
    {
        iP->Clear();
    }
    
    fCanvasPixelHisto->Update();
}


/*
     draw a histogram with the timing for each pixel

*/
bool VDisplay::drawImageBorderTZero()
{
    fEventLoop->getData()->setTelID( fTelescope );
    
    // get maximum/minimum  tzero per pixel
    double y_max = 0.;
    double y_min = 0.;
    int z = 0;
    for( unsigned int i = 0; i < fEventLoop->getData()->getPulseTime().size(); i++ )
    {
        if( fEventLoop->getData()->getImage()[i] || fEventLoop->getData()->getBorder()[i] )
        {
            if( fEventLoop->getData()->getTraceAverageTime( true )[i] > y_max )
            {
                y_max = fEventLoop->getData()->getTraceAverageTime( true )[i];
            }
            if( fEventLoop->getData()->getTraceAverageTime( true )[i] < y_min && fEventLoop->getData()->getTraceAverageTime( true )[i] > -90. )
            {
                y_min = fEventLoop->getData()->getTraceAverageTime( true )[i];
            }
        }
        if( fEventLoop->getData()->getSums()[i] > 0. )
        {
            z++;
        }
    }
    if( z < 2 )
    {
        return false;
    }
    double y_diff = 0.1 * ( y_max - y_min );
    y_max += y_diff;
    y_min -= y_diff;
    
    TH1D* hAll = new TH1D( "hAll", "", 100, y_min, y_max );
    hAll->SetStats( 0 );
    hAll->SetLineWidth( 2 );
    hAll->SetYTitle( "# of pixel" );
    hAll->GetXaxis()->SetTitleSize( 0.06 );
    hAll->GetXaxis()->SetLabelSize( 0.05 );
    hAll->GetXaxis()->SetTitleOffset( 0.7 );
    hAll->GetYaxis()->SetTitleSize( 0.06 );
    hAll->GetYaxis()->SetLabelSize( 0.05 );
    hAll->GetYaxis()->SetTitleOffset( 0.7 );
    hAll->GetXaxis()->SetTitleColor( 1 );
    
    TH1D* hImage = new TH1D( "hImage", "", 100, y_min, y_max );
    hImage->SetStats( 0 );
    hImage->SetLineWidth( 2 );
    hImage->SetLineColor( 4 );
    
    if( fEventLoop->getSumWindowStart_T_method() == 2 )
    {
        hAll->SetXTitle( "taverage" );
    }
    else if( fEventLoop->getSumWindowStart_T_method() == 1 )
    {
        hAll->SetXTitle( "tzero" );
    }
    for( unsigned int i = 0; i < fEventLoop->getPulseTime().size(); i++ )
    {
        if( fEventLoop->getPulseTime()[i] > 0. )
        {
            hAll->Fill( fEventLoop->getPulseTime()[i] );
            if( fEventLoop->getData()->getImage()[i] || fEventLoop->getData()->getBorder()[i] )
            {
                hImage->Fill( fEventLoop->getPulseTime() [i] );
            }
        }
    }
    if( hAll->GetEntries() > 0 )
    {
        gPad->SetLogy( 1 );
    }
    hAll->DrawCopy();
    
    if( hImage->GetEntries() > 0 )
    {
        hImage->DrawCopy( "same" );
    }
    
    delete hAll;
    delete hImage;
    
    return true;
}


bool VDisplay::drawImageBorderCharge()
{
    fEventLoop->getData()->setTelID( fTelescope );
    
    // get maximum charge per pixel
    double y_max = 0.;
    int z = 0;
    for( unsigned int i = 0; i < fEventLoop->getData()->getSums().size(); i++ )
    {
        //		if( fEventLoop->getData()->getImage()[i] || fEventLoop->getData()->getBorder()[i] )
        {
            if( fEventLoop->getData()->getSums()[i] > y_max )
            {
                y_max = fEventLoop->getData()->getSums()[i];
            }
        }
        if( fEventLoop->getData()->getSums()[i] > 0. )
        {
            z++;
        }
    }
    if( y_max > 0. )
    {
        y_max = log10( y_max );
    }
    else
    {
        y_max = 1.;
    }
    
    if( z < 2 )
    {
        return false;
    }
    
    TH1D* hAll = new TH1D( "hAll", "", 100, 0., y_max * 1.3 );
    hAll->SetStats( 0 );
    hAll->SetLineWidth( 2 );
    hAll->SetXTitle( "log_{10} charge" );
    hAll->SetYTitle( "# of pixel" );
    hAll->GetXaxis()->SetTitleSize( 0.06 );
    hAll->GetXaxis()->SetLabelSize( 0.05 );
    hAll->GetXaxis()->SetTitleOffset( 0.7 );
    hAll->GetYaxis()->SetTitleSize( 0.06 );
    hAll->GetYaxis()->SetLabelSize( 0.05 );
    hAll->GetYaxis()->SetTitleOffset( 0.7 );
    hAll->GetXaxis()->SetTitleColor( 1 );
    
    TH1D* hImage = new TH1D( "hImage", "", 100, 0., y_max * 1.3 );
    hImage->SetStats( 0 );
    hImage->SetLineWidth( 2 );
    hImage->SetLineColor( 2 );
    TH1D* hBorder = new TH1D( "hBorder", "", 100, 0., y_max * 1.3 );
    hBorder->SetStats( 0 );
    hBorder->SetLineWidth( 2 );
    hBorder->SetLineColor( 3 );
    
    for( unsigned int i = 0; i < fEventLoop->getData()->getSums().size(); i++ )
    {
        if( fEventLoop->getData()->getSums()[i] > 0. )
        {
            hAll->Fill( log10( fEventLoop->getData()->getSums()[i] ) );
        }
        if( fEventLoop->getData()->getImage()[i] )
        {
            hImage->Fill( log10( fEventLoop->getData()->getSums()[i] ) );
        }
        if( fEventLoop->getData()->getBorder()[i] && !fEventLoop->getData()->getImage()[i] )
        {
            hBorder->Fill( log10( fEventLoop->getData()->getSums()[i] ) );
        }
    }
    if( hAll->GetEntries() > 0 )
    {
        gPad->SetLogy( 1 );
    }
    hAll->DrawCopy();
    
    if( hBorder->GetEntries() > 0 )
    {
        hBorder->DrawCopy( "same" );
    }
    if( hImage->GetEntries() > 0 )
    {
        hImage->DrawCopy( "same" );
    }
    
    // draw telescope number
    char tname[200];
    sprintf( tname, "T%d", fTelescope + 1 );
    TText* iT = new TText( 0.5 * hImage->GetXaxis()->GetXmax(), 0.8 * hImage->GetYaxis()->GetXmax(), tname );
    iT->SetNDC();
    iT->SetX( 0.7 );
    iT->SetY( 0.7 );
    iT->SetTextSize( 3.*iT->GetTextSize() );
    iT->Draw();
    
    delete hAll;
    delete hImage;
    delete hBorder;
    
    return true;
}


/*!

    This draws the timing gradient graphs in the tgrad tab

*/
bool VDisplay::drawTgradGraphs()
{
    if( fDebug )
    {
        cout << "VDisplay::drawTgradGraphs() ";
        if( fTabAna )
        {
            cout << fTabAna->GetCurrent();
        }
        cout << " (should be 3)" << endl;
    }
    
    fEventLoop->getData()->setTelID( fTelescope );
    // get time gradient (for DP1)
    TGraphErrors* xgraph = fEventLoop->getData()->getXGraph( false );
    if( !xgraph || xgraph->GetN() < 1 )
    {
        return false;
    }
    
    double y_min = 0.;
    double y_max = 20.;
    double x_min =  1.e10;
    double x_max = -1.e10;
    y_min = 1.e5;
    y_max = -1.e5;
    double x = 0.;
    double y = 0.;
    double yE = 0.;
    for( int i = 0; i < xgraph->GetN(); i++ )
    {
        xgraph->GetPoint( i, x, y );
        yE = xgraph->GetErrorY( i );
        if( y + yE > y_max )
        {
            if( yE < y )
            {
                y_max = y + yE;
            }
            else
            {
                y_max = y * 1.2;
            }
        }
        if( y - yE < y_min )
        {
            if( yE < y )
            {
                y_min = y - yE;
            }
            else
            {
                y_min = y * 0.8;
            }
        }
        if( x > x_max )
        {
            x_max = x;
        }
        if( x < x_min )
        {
            x_min = x;
        }
    }
    double y_add = 0.2 * ( y_max - y_min );
    y_min -= y_add;
    y_max += y_add;
    
    if( TMath::Abs( y_min - y_max ) < 1.e-2 )
    {
        y_min = -0.5;
        y_max = 0.5;
    }
    
    x_max += 1.2;
    x_min -= 1.2;
    
    // fix range to readout window
    y_min = -5.;  // (just to make it look better)
    if( fEventLoop->getNSamples() > 0 )
    {
        y_max = fEventLoop->getNSamples();
    }
    else if( fEventLoop->getImageCleaningParameter( false, fEventLoop->getTelID() )
             && fEventLoop->getImageCleaningParameter( false, fEventLoop->getTelID() )->fNNOpt_nBinsADC > 0 )
    {
        y_max = fEventLoop->getImageCleaningParameter( false, fEventLoop->getTelID() )->fNNOpt_nBinsADC;
    }
    if( fEventLoop->getDigitalFilterMethod() > 0 )
    {
        y_max *= fEventLoop->getDigitalFilterUpSample();
    }
    
    TH2F* h1 = new TH2F( "h1", "", 1, x_min, x_max, 1, y_min, y_max );
    h1->SetStats( 0 );
    h1->SetTitle( "" );
    h1->GetXaxis()->SetTitle( "PMT position on long axis [degrees]" );
    h1->GetXaxis()->SetTitleSize( 0.06 );
    h1->GetXaxis()->SetLabelSize( 0.05 );
    h1->GetXaxis()->SetTitleOffset( 0.8 );
    
    char histitle[200];
    if( TMath::Abs( fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope )
                    - TMath::Floor( fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) ) ) < 1.e-3 )
    {
        sprintf( histitle, "pulse time [samples (%dns)]", ( int )fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) );
    }
    else
    {
        sprintf( histitle, "pulse time [samples (%.1fns)]", fEventLoop->getDetectorGeo()->getLengthOfSampleTimeSlice( fTelescope ) );
    }
    h1->GetYaxis()->SetTitle( histitle );
    h1->GetYaxis()->SetTitleSize( 0.06 );
    h1->GetYaxis()->SetLabelSize( 0.05 );
    h1->GetYaxis()->SetTitleOffset( 0.7 );
    
    h1->GetXaxis()->SetTitle( "PMT position on long axis [deg]" );
    h1->DrawCopy();
    if( xgraph && xgraph->GetN() > 1 )
    {
        xgraph->SetMarkerSize( 0.8 );
        xgraph->SetMarkerColor( 4 );
        xgraph->SetMarkerStyle( 20 );
        xgraph->Draw( "P" );
    }
    
    delete h1;
    
    return true;
}


/*!

 */
void VDisplay::defineGui()
{
    char i_text[200];
    // layout hints
    fL1 = new TGLayoutHints( kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2 );
    fL2 = new TGLayoutHints( kLHintsLeft | kLHintsBottom , 0, 0, 0, 0 );
    fL4 = new TGLayoutHints( kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0 );
    //  fL5 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kFixedWidth, 2, 2, 2, 2);
    fL5 = new TGLayoutHints( kLHintsTop );
    fL6 = new TGLayoutHints( kLHintsTop | kLHintsLeft, 2, 2, 2, 2 );
    fL8 = new TGLayoutHints( kLHintsTop | kLHintsRight, 2, 2, 2, 2 );
    
    fL3 = new TGLayoutHints( kLHintsTop | kLHintsRight , 0, 0, 0, 0 );
    fL7 = new TGLayoutHints( kLHintsBottom | kLHintsRight, 0, 0, 0, 0 );
    
    // main menu bar
    
    fMenuBar = new TGMenuBar( this, 1, 1, kHorizontalFrame );
    fMenuFile = new TGPopupMenu( fClient->GetRoot() );
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry( "&quit", M_FILE_EXIT );
    fMenuFile->Associate( this );
    fMenuBar->AddPopup( "&file", fMenuFile, fL6 );
    fMenuPrint = new TGPopupMenu( fClient->GetRoot() );
    fMenuPrint->AddEntry( "&print camera...", M_FILE_PRINTC );
    fMenuPrint->AddEntry( "&print FADC...", M_FILE_PRINTF );
    fMenuPrint->AddEntry( "&print calibration...", M_FILE_PRINTCA );
    fMenuPrint->AddEntry( "&print analysis...", M_FILE_PRINTA );
    fMenuPrint->AddEntry( "&print histo...", M_FILE_PRINTTG );
    fMenuPrint->AddEntry( "&print shower rec....", M_FILE_PRINTB );
    fMenuPrint->AddEntry( "&make camera movie...", M_FILE_MOVIE );
    fMenuPrint->Associate( this );
    fMenuBar->AddPopup( "&print", fMenuPrint, fL6 );
    fMenuOpt = new TGPopupMenu( fClient->GetRoot() );
    fMenuOpt->AddEntry( "color scale", M_OPT_COL_SCHE );
    fMenuOpt->AddEntry( "B&W scale", M_OPT_BW_SCHE );
    fMenuOpt->AddEntry( "plot all FADC image pixel samples", M_OPT_SING );
    fMenuOpt->AddEntry( "fix scale", M_OPT_FIX );
    fMenuOpt->AddSeparator();
    fMenuOpt->AddEntry( "dump dead channels", M_OPT_DEAD );
    fMenuOpt->AddEntry( "dump image/borderchannels", M_OPT_IMAGE );
    fMenuOpt->Associate( this );
    fMenuBar->AddPopup( "&options", fMenuOpt, fL6 );
    AddFrame( fMenuBar, fL4 );
    fMenuOpt->CheckEntry( M_OPT_COL_SCHE );
    
    // status bar (at the bottom of the window)
    
    fStatusBar = new TGStatusBar( this, 50, 10, kHorizontalFrame );
    Int_t parts[] = {45, 30, 25};
    // three divisions (file name, eventnumber, ??)
    fStatusBar->SetParts( parts, 3 );
    AddFrame( fStatusBar, new TGLayoutHints( kLHintsBottom | kLHintsExpandX, 0, 0, 0, 0 ) );
    
    string iTemp = "reading file ";
    iTemp += fEventLoop->getDataFileName();
    fStatusBar->SetText( iTemp.c_str(), 0 );
    fStatusBar->SetText( "now at event ", 1 );
    
    // main frame
    fFrameTable = new TGCompositeFrame( this, 100, 100, kRaisedFrame );
    AddFrame( fFrameTable, fL1 );
    
    UInt_t ywm = GetHeight();
    UInt_t xwm = GetWidth();
    // camera canvas is quadratic
    int iCanvasCameraSize = ( int )( 570. / 1000. * ( double )xwm );
    int iCanvasFADCSIZEX = ( int )( 400. / 700. * ( double )ywm * 1.03 );
    int iCanvasFADCSIZEY = ( int )( 330. / 700. * ( double )ywm );
    
    fEmCanvasCamera = new TRootEmbeddedCanvas( "mainCanvas", fFrameTable, iCanvasCameraSize, iCanvasCameraSize );
    fCanvasCamera = ( TCanvas* )fEmCanvasCamera->GetCanvas();
    fFrameTable->AddFrame( fEmCanvasCamera, fL2 );
    
    // event control
    fFrameControl = new TGCompositeFrame( fFrameTable, iCanvasCameraSize, 100 );
    
    fHorizontalEvent = new TGHorizontalFrame( fFrameControl, 100, 60 );
    fButtonNext = new TGTextButton( fHorizontalEvent, "&next", B_NEXT );
    fButtonNext->Associate( this );
    fHorizontalEvent->AddFrame( fButtonNext, fL6 );
    fButtonAutorun = new TGTextButton( fHorizontalEvent, "start &auto run", B_AUTOSTART );
    fButtonAutorun->Associate( this );
    fHorizontalEvent->AddFrame( fButtonAutorun, fL6 );
    fLabelGoto  = new TGLabel( fHorizontalEvent, "goto event: " );
    fHorizontalEvent->AddFrame( fLabelGoto, fL6 );
    fNEntryGoto = new TGNumberEntry( fHorizontalEvent, 1, 10, B_NGOTO, ( TGNumberFormat::EStyle ) 0 );
    fNEntryGoto->SetLimits( TGNumberFormat::kNELLimitMin, 0 );
    fNEntryGoto->Associate( this );
    fHorizontalEvent->AddFrame( fNEntryGoto, fL4 );
    fComboCameraView = new TGComboBox( fHorizontalEvent, B_VIEW );
    // the entries must be the same as E_cameraIdent
    fComboCameraView->AddEntry( "charge", 0 );
    fComboCameraView->AddEntry( "tzeros", 1 );
    fComboCameraView->AddEntry( "trigger", 2 );
    fComboCameraView->AddEntry( "hit", 3 );
    fComboCameraView->AddEntry( "hilo", 4 );
    fComboCameraView->AddEntry( "timing", 5 );
    fComboCameraView->AddEntry( "sum window", 6 );
    fComboCameraView->AddEntry( "sum window start", 7 );
    fComboCameraView->AddEntry( "pedestal (high)", 8 );
    fComboCameraView->AddEntry( "pedvar (high)", 9 );
    fComboCameraView->AddEntry( "pedestal (low)", 10 );
    fComboCameraView->AddEntry( "pedvar (low)", 11 );
    fComboCameraView->AddEntry( "gains (high)", 12 );
    fComboCameraView->AddEntry( "gainvar (high)", 13 );
    fComboCameraView->AddEntry( "gains (low)", 14 );
    fComboCameraView->AddEntry( "gainvar (low)", 15 );
    fComboCameraView->AddEntry( "toff (high)", 16 );
    fComboCameraView->AddEntry( "toff (low)", 17 );
    fComboCameraView->AddEntry( "lowgain mult", 18 );
    fComboCameraView->AddEntry( "tzero average (high)", 19 );
    fComboCameraView->AddEntry( "tzero average (low)", 20 );
    fComboCameraView->AddEntry( "status (high)", 21 );
    fComboCameraView->AddEntry( "status (low)", 22 );
    fComboCameraView->AddEntry( "L1 rates", 23 );
    fComboCameraView->AddEntry( "HV", 24 );
    fComboCameraView->AddEntry( "currents", 25 );
    fComboCameraView->AddEntry( "trigger-evndisp", 26 );
    fComboCameraView->AddEntry( "template (frogs)", 27 );
    fComboCameraView->AddEntry( "model3D", 28 );
    fComboCameraView->AddEntry( "clusterID", 29 );
    fComboCameraView->AddEntry( "PE", 30 );
    fComboCameraView->Select( 0 );
    fComboCameraView->Associate( this );
    fComboCameraView->Resize( 110, 20 );
    
    fComboTelescopeN = new TGComboBox( fHorizontalEvent, B_TELN );
    if( fEventLoop->getTeltoAna().size() > 1 )
    {
        fComboTelescopeN->AddEntry( "All telescopes", 0 );
        fComboTelescopeN->AddEntry( "All in one", -1 );
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            sprintf( i_text, "Telescope %d", fEventLoop->getTeltoAna()[i] + 1 );
            fComboTelescopeN->AddEntry( i_text, fEventLoop->getTeltoAna()[i] + 1 );
        }
    }
    else if( fEventLoop->getTeltoAna().size() == 1 )
    {
        sprintf( i_text, "Telescope %d", fEventLoop->getTeltoAna()[0] + 1 );
        fComboTelescopeN->AddEntry( i_text, 0 );
    }
    // field view works currently only with four telescopes (assume VERITAS configuration)
    // (GM) works only for four telescopes
    //  if( fEventLoop->getTeltoAna().size() == 4 ) fComboTelescopeN->AddEntry( "Field view", -2 );
    // all in one camera only if there is more than one telescopes
    // default is plotting all telescopes if there are less than 10 telescopes
    fComboTelescopeN->Select( 0 );
    fComboTelescopeN->Associate( this );
    fComboTelescopeN->Resize( 110, 20 );
    fButtonExit = new TGTextButton( fHorizontalEvent, "&Quit", M_FILE_EXIT );
    fButtonExit->SetTextColor( 2, false );
    fButtonExit->Associate( this );
    fHorizontalEvent->AddFrame( fButtonExit, fL8 );
    fButtonDump = new TGTextButton( fHorizontalEvent, "&Dump", M_FILE_DUMP );
    fButtonDump->SetTextColor( 2, false );
    fButtonDump->Associate( this );
    fHorizontalEvent->AddFrame( fButtonDump, fL8 );
    fHorizontalEvent->AddFrame( fComboTelescopeN, fL6 );
    fHorizontalEvent->AddFrame( fComboCameraView, fL6 );
    fFrameControl->AddFrame( fHorizontalEvent, fL4 );
    fFrameTable->AddFrame( fFrameControl, fL5 );
    
    /*  fEmCanvasCamera = new TRootEmbeddedCanvas( "mainCanvas", fFrameTable, iCanvasCameraSize, iCanvasCameraSize );
      fCanvasCamera = (TCanvas*)fEmCanvasCamera->GetCanvas();
      fFrameTable->AddFrame( fEmCanvasCamera, fL2 ); */
    
    // start with charge display
    fCameraDisplay = C_CHARGE;
    fCanvasesNx = ( unsigned int )sqrt( ( double )fCamera.size() );
    fCanvasesNy = fCamera.size() / fCanvasesNx;
    
    // tab with analysis canvas, button
    fTabAna = new TGTab( fFrameTable, 400, 600 );
    
    // general information tab
    fFrameInfo = fTabAna->AddTab( "info" );
    fCompInfo = new TGCompositeFrame( fFrameInfo, 60, 20, kVerticalFrame );
    fFrameInfo->AddFrame( fCompInfo, fL6 );
    fEmInfo = new TRootEmbeddedCanvas( "canvasInfo", fCompInfo, iCanvasFADCSIZEX, int( iCanvasCameraSize - 10 ) );
    fCompInfo->AddFrame( fEmInfo );
    fCanvasInfo = fEmInfo->GetCanvas();
    fPaveInfo = new TPaveText( 0.01, 0.01, 0.99, 0.99 );
    fPaveInfo->SetFillColor( 10 );
    
    // FADC tab
    fFrameFADC = fTabAna->AddTab( "FADC" );
    fAnaDisplay = F_FADC;
    fCompFADC = new TGCompositeFrame( fFrameFADC, 60, 20, kVerticalFrame );
    fFrameFADC->AddFrame( fCompFADC, fL6 );
    fEmFADC = new TRootEmbeddedCanvas( "canvasFADC", fCompFADC, iCanvasFADCSIZEX, iCanvasFADCSIZEY );
    fCanvasFADC = fEmFADC->GetCanvas();
    fCanvasFADC->SetHighLightColor( 10 );
    fCanvasFADC->SetTickx( 1 );
    fCanvasFADC->SetTicky( 1 );
    fCompFADC->AddFrame( fEmFADC, fL6 );
    fEmFADCText = new TRootEmbeddedCanvas( "canvasFADCText", fCompFADC, iCanvasFADCSIZEX, int( iCanvasFADCSIZEY * 0.7 ) );
    fCanvasFADCText = fEmFADCText->GetCanvas();
    fCanvasFADCText->SetHighLightColor( 10 );
    fCompFADC->AddFrame( fEmFADCText, fL6 );
    fGroupFADC = new TGHorizontalFrame( fCompFADC, 100, 60 );
    fCompFADC->AddFrame( fGroupFADC, fL6 );
    fButtonFADCFit = new TGTextButton( fGroupFADC, "fit", B_FITPMT );
    fGroupFADC->AddFrame( fButtonFADCFit, fL6 );
    fButtonFADCFit->Associate( this );
    fButtonFADCFit->SetState( kButtonDisabled );
    fButtonFADCset = new TGTextButton( fGroupFADC, "add channel", B_ADDPMT );
    fGroupFADC->AddFrame( fButtonFADCset, fL6 );
    fButtonFADCset->Associate( this );
    fButtonFADCset->SetState( kButtonDisabled );
    fButtonFADCunset = new TGTextButton( fGroupFADC, "remove channel", B_REMOVEPMT );
    fGroupFADC->AddFrame( fButtonFADCunset, fL6 );
    fButtonFADCunset->Associate( this );
    fButtonFADCunset->SetState( kButtonDisabled );
    fButtonFADCreset = new TGTextButton( fGroupFADC, "reset", B_RESETPMT );
    fGroupFADC->AddFrame( fButtonFADCreset, fL6 );
    fButtonFADCreset->Associate( this );
    fButtonFADCreset->SetState( kButtonDisabled );
    fLabelFADCsearch = new TGLabel( fGroupFADC, "show channel:" );
    fGroupFADC->AddFrame( fLabelFADCsearch, fL6 );
    fNEntryFADCsearch = new TGNumberEntry( fGroupFADC, 1, 4, B_NFADC, ( TGNumberFormat::EStyle ) 0 );
    fNEntryFADCsearch->SetLimits( TGNumberFormat::kNELLimitMinMax, 0, fEventLoop->getAnalyzer()->getImage().size() - 1 );
    fNEntryFADCsearch->Associate( this );
    fGroupFADC->AddFrame( fNEntryFADCsearch, fL6 );
    // calibration tab
    fFrameCal = fTabAna->AddTab( "calibration" );
    fCompAna = new TGCompositeFrame( fFrameCal, 60, 20, kVerticalFrame );
    fFrameCal->AddFrame( fCompAna, fL6 );
    fEmCal = new TRootEmbeddedCanvas( "canvasCal", fCompAna, iCanvasFADCSIZEX, iCanvasFADCSIZEY );
    fCanvasCal = fEmCal->GetCanvas();
    fCompAna->AddFrame( fEmCal );
    // timing graphs tab
    fFrameTgrad = fTabAna->AddTab( "histos" );
    fCompAna = new TGCompositeFrame( fFrameTgrad, 60, 20, kVerticalFrame );
    fFrameTgrad->AddFrame( fCompAna, fL6 );
    fEmTgrad = new TRootEmbeddedCanvas( "canvasTgrad", fCompAna, iCanvasFADCSIZEX, int( iCanvasCameraSize - 10 ) );
    fCanvasPixelHisto = fEmTgrad->GetCanvas();
    fCanvasPixelHisto->Divide( 1, 3 );
    fCompAna->AddFrame( fEmTgrad );
    // Analysis tab
    fFrameAna = fTabAna->AddTab( "analysis" );
    fCompAna = new TGCompositeFrame( fFrameAna, 60, 20, kVerticalFrame );
    fFrameAna->AddFrame( fCompAna, fL6 );
    fEmAna = new TRootEmbeddedCanvas( "canvasAna", fCompAna,  iCanvasFADCSIZEX, int( iCanvasCameraSize - 10 ) );
    fCanvasAna = fEmAna->GetCanvas();
    fCanvasAna->Divide( 2, 3 );
    fCompAna->AddFrame( fEmAna );
    // option tab
    // frame with run options
    fFrameOpt = fTabAna->AddTab( "options" );
    fCompOpt =  new TGCompositeFrame( fFrameOpt, 60, 20, kVerticalFrame );
    fFrameOpt->AddFrame( fCompOpt, fL6 );
    fGroupOptRun = new TGGroupFrame( fCompOpt, "autorun options", kHorizontalFrame );
    fCompOpt->AddFrame( fGroupOptRun, fL4 );
    fLabelOptAutoRun = new TGLabel( fGroupOptRun, "timing delay (ms):" );
    fGroupOptRun->AddFrame( fLabelOptAutoRun, fL6 );
    fNEntryOAutoRun = new TGNumberEntry( fGroupOptRun, 5000, 7, B_NAUTO, ( TGNumberFormat::EStyle ) 0 );
    fNEntryOAutoRun->SetLimits( TGNumberFormat::kNELLimitMin, 1 );
    fGroupOptRun->AddFrame( fNEntryOAutoRun, fL6 );
    fNEntryOAutoRun->Associate( this );
    fLabelOptInc = new TGLabel( fGroupOptRun, "event increment:" );
    fGroupOptRun->AddFrame( fLabelOptInc, fL6 );
    fNEntryOInc = new TGNumberEntry( fGroupOptRun, 1, 7, B_NINC, ( TGNumberFormat::EStyle ) 0 );
    fNEntryOInc->SetLimits( TGNumberFormat::kNELLimitMin, 1 );
    fNEntryOInc->Associate( this );
    fGroupOptRun->AddFrame( fNEntryOInc, fL6 );
    
    // frame with display options
    fGroupOptDis = new TGButtonGroup( fCompOpt, "coordinate system", kHorizontalFrame );
    fCompOpt->AddFrame( fGroupOptDis, fL4 );
    //  fChButtonColor = new TGCheckButton( fGroupOptDis, "&color scheme" , B_DCOLOR  );
    //  fChButtonColor->Associate( this );
    fRadioB1 = new TGRadioButton( fGroupOptDis, "&no IDs" , B_DNONE );
    fRadioB2 = new TGRadioButton( fGroupOptDis, "&channel IDs ", B_DCHANNEL );
    fRadioB3 = new TGRadioButton( fGroupOptDis, "&tubes IDs", B_DTUBE );
    fRadioB4 = new TGRadioButton( fGroupOptDis, "c&oordinate system", B_DCOOR );
    fRadioB1->Associate( this );
    fRadioB2->Associate( this );
    fRadioB3->Associate( this );
    fRadioB4->Associate( this );
    fRadioB1->SetState( kButtonDown );
    fGroupOptDis->Show();
    
    // frame with telescope selections
    
    fGroupTelFrame = new TGGroupFrame( fCompOpt, "run parameters" );
    fCompOpt->AddFrame( fGroupTelFrame, fL4 );
    
    fGroupOptTel = new TGButtonGroup( fGroupTelFrame, "telescopes", kHorizontalFrame );
    fGroupTelFrame->AddFrame( fGroupOptTel, fL4 );
    
    fRadioTA = new TGRadioButton( fGroupOptTel, "all", B_TALL );
    fRadioTA->Associate( this );
    // don't write more than 8 of this radio buttons (no space for more...)
    if( fEventLoop->getTeltoAna().size() < 8 && fEventLoop->getTeltoAna().size() > 1 )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            sprintf( i_text, "%d ", fEventLoop->getTeltoAna()[i] + 1 );
            fRadioTel.push_back( new TGRadioButton( fGroupOptTel, i_text, 100 + fEventLoop->getTeltoAna()[i] ) );
            fRadioTel.back()->Associate( this );
        }
    }
    fRadioTA->SetState( kButtonDown );
    fGroupOptTel->Show();
    
    // frame with analysis parameters/cuts
    fGroupAnaFrame = new TGHorizontalFrame( fGroupTelFrame, 60, 20 );
    fGroupTelFrame->AddFrame( fGroupAnaFrame, fL4 );
    
    // frame with analysis parameters
    /*	fGroupOptAna = new TGGroupFrame( fGroupAnaFrame, "analysis parameters", kVerticalFrame );
    	fGroupAnaFrame->AddFrame( fGroupOptAna, fL4 );
    	fLabelOptSum = new TGLabel( fGroupOptAna, "summation start:" );
    	fGroupOptAna->AddFrame( fLabelOptSum, fL6 );
    	fNEntryOSum = new TGNumberEntry( fGroupOptAna, fEventLoop->getAnalyzer()->getSumFirst(), 3, B_NSUM, ( TGNumberFormat::EStyle ) 0 );
    	fNEntryOSum->SetLimits( TGNumberFormat::kNELLimitMinMax, 0, 64 );
    	fGroupOptAna->AddFrame( fNEntryOSum, fL6 );
    	fLabelOptWin = new TGLabel( fGroupOptAna, "summation window:" );
    	fGroupOptAna->AddFrame( fLabelOptWin, fL6 );
    	fNEntryOWin = new TGNumberEntry( fGroupOptAna, fEventLoop->getAnalyzer()->getSumWindow(), 3, B_NWIN, ( TGNumberFormat::EStyle ) 0 );
    	fNEntryOWin->SetLimits( TGNumberFormat::kNELLimitMinMax, 0, 64 );
    	fNEntryOWin->SetState( false );               // pedestal calculation depends on summation window !!!
    	fGroupOptAna->AddFrame( fNEntryOWin, fL6 );
    	fLabelOptIma = new TGLabel( fGroupOptAna, "image threshold:" );
    	fGroupOptAna->AddFrame( fLabelOptIma, fL6 );
            if( fEventLoop->getAnalyzer()->getImageCleaningParameter() )
            {
                fNEntryOIma = new TGNumberEntry( fGroupOptAna, fEventLoop->getAnalyzer()->getImageCleaningParameter()->fimagethresh, 5, B_NIMA, ( TGNumberFormat::EStyle ) 2 );
            }
            else
            {
                fNEntryOIma = new TGNumberEntry( fGroupOptAna, 5., 5, B_NIMA, ( TGNumberFormat::EStyle ) 2 );
            }
    	fNEntryOIma->SetLimits( TGNumberFormat::kNELLimitMin, 0. );
    	fGroupOptAna->AddFrame( fNEntryOIma, fL6 );
    	fLabelOptBor = new TGLabel( fGroupOptAna, "border threshold:" );
    	fGroupOptAna->AddFrame( fLabelOptBor, fL6 );
            if( fEventLoop->getAnalyzer()->getImageCleaningParameter() )
            {
                fNEntryOBor = new TGNumberEntry( fGroupOptAna, fEventLoop->getAnalyzer()->getImageCleaningParameter()->fborderthresh, 5, B_NBOR, ( TGNumberFormat::EStyle ) 2 );
            }
            else
            {
                fNEntryOBor = new TGNumberEntry( fGroupOptAna, 2.5, 5, B_NBOR, ( TGNumberFormat::EStyle ) 2 );
            }
    	fNEntryOBor->SetLimits( TGNumberFormat::kNELLimitMin, 0. );
    	fGroupOptAna->AddFrame( fNEntryOBor, fL6 ); */
    
    // frame with cuts
    fGroupOptCut = new TGGroupFrame( fGroupAnaFrame, "cuts", kVerticalFrame );
    fGroupAnaFrame->AddFrame( fGroupOptCut, fL4 );
    fLabelOptATri = new TGLabel( fGroupOptCut, "# telescope trigger:" );
    fGroupOptCut->AddFrame( fLabelOptATri, fL6 );
    fNEntryOATri = new TGNumberEntry( fGroupOptCut, fEventLoop->getCutNArrayTrigger(), 7, B_NATRI, ( TGNumberFormat::EStyle ) 0 );
    fNEntryOATri->SetLimits( TGNumberFormat::kNELLimitMin, 0 );
    fGroupOptCut->AddFrame( fNEntryOATri, fL6 );
    fNEntryOATri->Associate( this );
    fLabelOptAIma = new TGLabel( fGroupOptCut, "# images:" );
    fGroupOptCut->AddFrame( fLabelOptAIma, fL6 );
    fNEntryOAIma = new TGNumberEntry( fGroupOptCut, 1, 7, B_NAIMA, ( TGNumberFormat::EStyle ) 0 );
    fNEntryOAIma->SetLimits( TGNumberFormat::kNELLimitMin, 0 );
    fGroupOptCut->AddFrame( fNEntryOAIma, fL6 );
    fNEntryOAIma->Associate( this );
    fLabelOptCuts = new TGLabel( fGroupOptCut, "analysis cuts: " );
    fGroupOptCut->AddFrame( fLabelOptCuts );
    fTextOptCuts = new TGTextEntry( fGroupOptCut, new TGTextBuffer( 100 ), B_TEXTCUT );
    fTextOptCuts->Resize( 10, fTextOptCuts->GetDefaultHeight() );
    fTextOptCuts->Associate( this );
    fGroupOptCut->AddFrame( fTextOptCuts, fL4 );
    
    fGroupSetFrame = new TGHorizontalFrame( fCompOpt, 60, 20 );
    fCompOpt->AddFrame( fGroupSetFrame, fL6 );
    fButtonOptSet = new TGTextButton( fGroupSetFrame, "&set", B_SET );
    fButtonOptSet->Associate( this );
    fGroupSetFrame->AddFrame( fButtonOptSet, fL6 );
    fButtonOptReset = new TGTextButton( fGroupSetFrame, "&reset", B_RESET );
    fButtonOptReset->Associate( this );
    fGroupSetFrame->AddFrame( fButtonOptReset, fL6 );
    
    // shower reconstruction (array view and reconstructed shower parameters)
    // birds eye tab
    fFrameBird = fTabAna->AddTab( "shower" );
    fCompBird = new TGCompositeFrame( fFrameBird, 60, 20, kVerticalFrame );
    fFrameBird->AddFrame( fCompBird, fL6 );
    // expect canvas to be a square
    fEmBird = new TRootEmbeddedCanvas( "canvasBird", fCompBird, iCanvasFADCSIZEX, iCanvasFADCSIZEX );
    fCanvasBird = fEmBird->GetCanvas();
    fCanvasBird->SetHighLightColor( 10 );
    fCanvasBird->SetTickx( 1 );
    fCanvasBird->SetTicky( 1 );
    fCompBird->AddFrame( fEmBird, fL6 );
    
    fFrameTable->AddFrame( fTabAna, fL3 );
}


/*!
    evaluate menu bar
*/
void VDisplay::subprocessMenu( Long_t parm1 )
{

    switch( parm1 )
    {
        // printing
        case M_FILE_PRINTC:
            printCanvas( fCanvasCamera );
            break;
        case M_FILE_PRINTF:
            printCanvas( fCanvasFADC );
            break;
        case M_FILE_PRINTCA:
            printCanvas( fCanvasCal );
            break;
        case M_FILE_PRINTA:
            printCanvas( fCanvasAna );
            break;
        case M_FILE_PRINTTG:
            printCanvas( fCanvasPixelHisto );
            break;
        case M_FILE_PRINTB:
            printCanvas( fCanvasBird );
            break;
        case M_FILE_MOVIE:
            if( fCameraMovie == true )
            {
                fCameraMovie = false;
                fMenuPrint->UnCheckEntry( M_FILE_MOVIE );
                fMoviePictNumber = 0;
            }
            else
            {
                const Char_t* filetypes[] =
                {
                    "All files",     "*", "GIF files",    "*.gif",
                    "PS files",    "*.ps", "EPS files",    "*.eps",
                    "ROOT files",    "*.root", "ROOT macros",   "*.C",
                    0,               0
                };
                TGFileInfo fi;
                fi.fFileTypes = filetypes;
                new TGFileDialog( fClient->GetRoot(), this, kFDSave, &fi );
                if( !fi.fFilename )
                {
                    break;
                }
                fMenuPrint->CheckEntry( M_FILE_MOVIE );
                fCameraMovie = true;
                fMovieFileName = fi.fFilename;
            }
            break;
        // quit
        case M_FILE_EXIT:
            CloseWindow();
            break;
            
        // plot signal amplitudes with different colors
        case M_OPT_COL_SCHE:
        
            gStyle->SetPalette( 1 );
            gStyle->SetNumberContours( 100 );
            if( fBoolDrawAllinOne )
            {
                break;    // can't draw several images in color scheme into one camera
            }
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                if( fCamera[fTelescopesToShow[i]]->getPlotColor() == 2 )
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 1 );
                    fMenuOpt->UnCheckEntry( M_OPT_BW_SCHE );
                    fMenuOpt->CheckEntry( M_OPT_COL_SCHE );
                }
                else if( fCamera[fTelescopesToShow[i]]->getPlotColor() == 0 )
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 1 );
                    fMenuOpt->CheckEntry( M_OPT_COL_SCHE );
                }
                else
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 0 );
                    fMenuOpt->UnCheckEntry( M_OPT_COL_SCHE );
                }
            }
            updateCamera( fCameraDisplay );
            break;
            
        // plot signal amplitudes with grey scale
        case M_OPT_BW_SCHE:
        
            if( fBoolDrawAllinOne )
            {
                break;    // can't draw several images in color scheme into one camera
            }
            gStyle->SetPalette( fBWNum, fBWPalette );
            gStyle->SetNumberContours( 100 );
            
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                if( fCamera[fTelescopesToShow[i]]->getPlotColor() == 2 )
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 0 );
                    fMenuOpt->UnCheckEntry( M_OPT_BW_SCHE );
                }
                else if( fCamera[fTelescopesToShow[i]]->getPlotColor() == 1 )
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 2 );
                    fMenuOpt->UnCheckEntry( M_OPT_COL_SCHE );
                    fMenuOpt->CheckEntry( M_OPT_BW_SCHE );
                }
                else
                {
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 2 );
                    fMenuOpt->CheckEntry( M_OPT_BW_SCHE );
                }
            }
            updateCamera( fCameraDisplay );
            break;
            
        // plot signals with fixed scale
        
        case M_OPT_FIX:
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                if( fCamera[fTelescopesToShow[i]]->getPlotColor() == 2 )
                {
                    fCamera[fTelescopesToShow[i]]->setFixScale( false );
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 0 );
                    fMenuOpt->UnCheckEntry( M_OPT_FIX );
                }
                else
                {
                    fCamera[fTelescopesToShow[i]]->setFixScale( true );
                    fCamera[fTelescopesToShow[i]]->setPlotColor( 2 );
                    fMenuOpt->CheckEntry( M_OPT_FIX );
                }
            }
            updateCamera( fCameraDisplay );
            break;
            
        // plot all traces of one telescope in FADC canvas
        case M_OPT_SING:
            fBoolDrawImageTraces = !fBoolDrawImageTraces;
            drawFADC( false );
            if( fBoolDrawImageTraces )
            {
                fMenuOpt->CheckEntry( M_OPT_SING );
            }
            else
            {
                fMenuOpt->UnCheckEntry( M_OPT_SING );
            }
            break;
        default:
            break;
            
        // dump dead channel list
        case M_OPT_DEAD:
            dumpDeadChannels();
            break;
        // dump list of image/border pixels
        case M_OPT_IMAGE:
            dumpImageBorderPixels();
            break;
    }
}


void VDisplay::subprocessCheckButton( Long_t parm1 )
{
    switch( parm1 )
    {
        default:
            break;
    }
}


void VDisplay::subprocessButton( Long_t parm1 )
{
    switch( parm1 )
    {
        case B_AUTOSTART:
            if( fAutoRunStatus )
            {
                fAutoRunStatus = false;
                fButtonAutorun->SetText( "&start autorun" );
            }
            else
            {
                fAutoRunStatus = true;
                fButtonAutorun->SetText( "&stop autorun" );
            }
            if( fAutoRunStatus )
            {
                fEventLoop->setNextEventStatus( true );
                processEvent();
                fEventLoop->setNextEventStatus( false );
            }
            break;
        case B_NEXT:
            if( fEventLoop->getNextEventStatus() )
            {
                fEventLoop->setNextEventStatus( false );
                fButtonNext->SetText( "next" );
            }
            else
            {
                fEventLoop->setNextEventStatus( true );
                fButtonNext->SetText( "stop" );
                // skip Monte Carlo events with no trigger
                if( fEventLoop->getReader()->isMC() && fEventLoop->getEventNumber() == 0 )
                {
                    if( fEventLoop->getNTel() == 1 )
                    {
                        fEventLoop->setCutNArrayTrigger( 1 );
                    }
                    else
                    {
                        fEventLoop->setCutNArrayTrigger( 2 );
                    }
                }
                // process event
                processEvent();
                fEventLoop->setNextEventStatus( false );
                fButtonNext->SetText( "next" );
            }
            break;
        case M_FILE_EXIT:
            CloseWindow();
            break;
        case M_FILE_DUMP:
            fEventLoop->dumpTreeData();
            dumpImageBorderPixels();
            break;
        case B_SET:
        {
            fTimingSleep = ( unsigned int )fNEntryOAutoRun->GetNumber();
            fNumEventIncrement = ( unsigned int )fNEntryOInc->GetNumber();
            vector<unsigned int> iTel;
            if( fRadioTA->GetState() )
            {
                for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
                {
                    iTel.push_back( fEventLoop->getTeltoAna()[i] );
                }
            }
            else
            {
                for( unsigned int i = 0; i < fRadioTel.size(); i++ )
                {
                    if( fRadioTel[i]->GetState() )
                    {
                        iTel.push_back( i );
                    }
                }
            }
            // reset all previous settings
            fEventLoop->resetRunOptions();
            if( iTel.size() == 1 )
            {
                fEventLoop->setCutSingleTelescope( true );
            }
            else
            {
                fEventLoop->setCutSingleTelescope( false );
            }
            // set new options
            for( unsigned int i = 0; i < iTel.size(); i++ )
            {
                fEventLoop->setTelID( iTel[i] );
                fEventLoop->setCutString( fTextOptCuts->GetText() );
                fEventLoop->setCutNArrayTrigger( ( int )fNEntryOATri->GetNumber() );
                fEventLoop->setCutNArrayImages( ( int )fNEntryOAIma->GetNumber() );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                setInfoText();
                for( unsigned int i = 0; i < iTel.size(); i++ )
                {
                    fEventLoop->setTelID( iTel[i] );
                    fEventLoop->getAnalyzer()->doAnalysis();
                    fEventLoop->getArrayAnalyzer()->doAnalysis();
                    updateCamera( fCameraDisplay );
                }
                if( fTabAna->GetCurrent() == 1 )
                {
                    drawFADC( false );
                }
            }
        }
        break;
        case B_RESET:
            fEventLoop->resetRunOptions();
            resetRunOptions();
            fNEntryOAutoRun->SetNumber( fTimingSleep );
            fNEntryOInc->SetNumber( fNumEventIncrement );
            fNEntryOATri->SetNumber( fEventLoop->getCutNArrayTrigger() );
            fNEntryOAIma->SetNumber( fEventLoop->getCutNArrayImages() );
            fTextOptCuts->SetText( fEventLoop->getCutString().c_str() );
            fRadioB1->SetState( kButtonDown );
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 0 );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                updateCamera( fCameraDisplay );
                setInfoText();
                fEventLoop->getAnalyzer()->doAnalysis();
                updateCamera( fCameraDisplay );
                if( fTabAna->GetCurrent() == 1 )
                {
                    drawFADC( false );
                }
            }
            break;
        case B_FITPMT:
            if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
            {
                plotFADCFit( fSelectedChan - 200000 );
            }
            break;
        case B_ADDPMT:
            if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
            {
                fEventLoop->getAnalyzer()->setTelID( fTelescope );
                //                fEventLoop->getAnalyzer()->addImageChannel( fSelectedChan - 200000 );
                if( fEventLoop->getAnalyzer()->getImageCleaner() )
                {
                    fEventLoop->getAnalyzer()->getImageCleaner()->addImageChannel( fSelectedChan - 200000 );
                }
                fEventLoop->getAnalyzer()->doAnalysis();
                updateCamera( fCameraDisplay );
                if( fTabAna->GetCurrent() == 1 )
                {
                    drawFADC( false );
                }
            }
            break;
        case B_REMOVEPMT:
            if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
            {
                fEventLoop->getAnalyzer()->setTelID( fTelescope );
                //                fEventLoop->getAnalyzer()->removeImageChannel( fSelectedChan - 200000 );
                if( fEventLoop->getAnalyzer()->getImageCleaner() )
                {
                    fEventLoop->getAnalyzer()->getImageCleaner()->removeImageChannel( fSelectedChan - 200000 );
                }
                fEventLoop->getAnalyzer()->doAnalysis();
                updateCamera( fCameraDisplay );
                if( fTabAna->GetCurrent() == 1 )
                {
                    drawFADC( false );
                }
            }
            break;
        case B_RESETPMT:
            if( fSelectedChan >= 200000 && fSelectedChan < 200000 + fEventLoop->getAnalyzer()->getImage().size() )
            {
                fEventLoop->getAnalyzer()->setTelID( fTelescope );
                //                fEventLoop->getAnalyzer()->resetImageChannel( fSelectedChan - 200000 );
                if( fEventLoop->getAnalyzer()->getImageCleaner() )
                {
                    fEventLoop->getAnalyzer()->getImageCleaner()->resetImageChannel( fSelectedChan - 200000 );
                }
                fEventLoop->getAnalyzer()->doAnalysis();
                updateCamera( fCameraDisplay );
                if( fTabAna->GetCurrent() == 1 )
                {
                    drawFADC( false );
                }
            }
            break;
            
        default:
            break;
    }
}


void VDisplay::subprocessRadioButton( Long_t parm1 )
{
    switch( parm1 )
    {
        case B_DNONE:
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 0 );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                updateCamera( fCameraDisplay );
            }
            break;
        case B_DCHANNEL:
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 1 );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                updateCamera( fCameraDisplay );
            }
            break;
        case B_DTUBE:
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 2 );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                updateCamera( fCameraDisplay );
            }
            break;
        case B_DCOOR:
            for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 3 );
            }
            if( fEventLoop->getEventNumber() > 0 )
            {
                updateCamera( fCameraDisplay );
            }
            break;
        case B_TALL:
            break;
        default:
            break;
    }
    if( parm1 >= 100 )
    {
        fEventLoop->setTelID( ( unsigned int )( parm1 - 100 ) );
        fNEntryOATri->SetNumber( fEventLoop->getCutNArrayTrigger() );
        fNEntryOAIma->SetNumber( fEventLoop->getCutNArrayImages() );
        fTextOptCuts->SetText( fEventLoop->getCutString().c_str() );
    }
}


void VDisplay::subprocessTextEnter( Long_t parm1 )
{
    switch( parm1 )
    {
        case B_NGOTO:
        {
            char c_ev[200];
            sprintf( c_ev, "searching for event %d", int( fNEntryGoto->GetNumber() ) );
            fStatusBar->SetText( c_ev, 1 );
            fEventLoop->gotoEvent( int( fNEntryGoto->GetNumber() ) );
            if( int( fNEntryGoto->GetNumber() ) != 0 )
            {
                sprintf( c_ev, "now at event %d", fEventLoop->getEventNumber() );
                fStatusBar->SetText( c_ev, 1 );
                updateCamera( fCameraDisplay );
                // update analysis tabs if necessary
                if( fTabAna->GetCurrent() == 1 && !fCameraTiming )
                {
                    drawFADC( false );
                }
                else if( fTabAna->GetCurrent() == 6 )
                {
                    fBirdsEye->draw( fCanvasBird );
                }
                else if( fTabAna->GetCurrent() == 3 )
                {
                    drawPixelHistos();
                }
            }
        }
        break;
        case B_NFADC:
            searchChannel( int( fNEntryFADCsearch->GetNumber() ) );
            break;
        default:
            break;
    }
}


void VDisplay::subprocessTextChanged( Long_t parm1 )
{
    switch( parm1 )
    {
        case B_NFADC:
            searchChannel( int( fNEntryFADCsearch->GetNumber() ) );
            break;
        default:
            break;
    }
}


void VDisplay::subprocessComboBox( Long_t parm1 )
{
    fTelescopesToShow.clear();
    if( fDebug )
    {
        cout << "VDisplay::subprocessComboBox " << parm1 << "\t" << fCamera.size() << "\t" << fComboTelescopeN->GetSelected() << endl;
    }
    
    fCameraDisplay = ( E_cameraIdent )fComboCameraView->GetSelected();
    // draw all telescopes (each in one pad, one pad beside another)
    if( fComboTelescopeN->GetSelected() == 0 )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            fTelescopesToShow.push_back( fEventLoop->getTeltoAna()[i] );
        }
        setCameraPads( false );
        fBoolDrawOne = false;
        fBoolDrawAllinOne = false;
        if( fBool_M_OPT_COL_SCHE_Checked )
        {
            fMenuOpt->CheckEntry( M_OPT_COL_SCHE );
        }
        else
        {
            fMenuOpt->UnCheckEntry( M_OPT_COL_SCHE );
        }
        if( fBool_M_OPT_BW_SCHE_Checked )
        {
            fMenuOpt->CheckEntry( M_OPT_BW_SCHE );
        }
        else
        {
            fMenuOpt->UnCheckEntry( M_OPT_BW_SCHE );
        }
        setColorScheme();
        // check radio buttons for channel printings
        for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
        {
            if( fRadioB1->GetState() == 1 )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 0 );
            }
            if( fRadioB2->GetState() == 1 )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 1 );
            }
            if( fRadioB3->GetState() == 1 )
            {
                fCamera[fTelescopesToShow[i]]->setPrintChannels( 2 );
            }
        }
    }
    // draw all cameras into one
    else if( fComboTelescopeN->GetSelected() == -1 )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            fTelescopesToShow.push_back( fEventLoop->getTeltoAna()[i] );
        }
        fBoolDrawAllinOne = true;
    }
    // draw all telescope in pads distributed according to their position in the field
    else if( fComboTelescopeN->GetSelected() == -2 )
    {
        for( unsigned int i = 0; i < fCamera.size(); i++ )
        {
            fTelescopesToShow.push_back( i );
        }
        setCameraPads( true );
        fBoolDrawOne = false;
        fBoolDrawAllinOne = false;
        setColorScheme();
    }
    ///////////////////////////////////
    // draw only one telescope
    else
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            if( fEventLoop->getTeltoAna()[i] == ( unsigned int )( fComboTelescopeN->GetSelected() - 1 ) )
            {
                fTelescopesToShow.push_back( fEventLoop->getTeltoAna()[i] );
                break;
            }
        }
        fBoolDrawOne = true;
        fBoolDrawAllinOne = false;
        fTelescope = fTelescopesToShow.back();
        if( fTabAna->GetCurrent() == 1 )
        {
            drawFADC( false );
            fCamera[fTelescopesToShow.back()]->setTubeSelected( fSelectedChan - 200000 );
        }
        else if( fTabAna->GetCurrent() == 3 )
        {
            drawPixelHistos();
        }
        setColorScheme();
        // check radio buttons for channel printings
        if( fRadioB1->GetState() == 1 )
        {
            fCamera[fTelescopesToShow.back()]->setPrintChannels( 0 );
        }
        if( fRadioB2->GetState() == 1 )
        {
            fCamera[fTelescopesToShow.back()]->setPrintChannels( 1 );
        }
        if( fRadioB3->GetState() == 1 )
        {
            fCamera[fTelescopesToShow.back()]->setPrintChannels( 2 );
        }
    }
    updateCamera( fCameraDisplay );
}


/*!

 */
void VDisplay::showWarning()
{
    gPad->Clear();
    TText i_Text( 0.2, 0.5, "functionality not available" );
    i_Text.Draw();
    gPad->Update();
}


/*!
    plotting of calibration data of all cameras in one pad makes no sense
*/
bool VDisplay::checkPlotIntentions( unsigned int i )
{
    if( fBoolDrawAllinOne )
    {
        if( E_cameraIdent( i ) != C_CHARGE && E_cameraIdent( i ) != C_TRIGGER && E_cameraIdent( i ) != C_TZERO )
        {
            showWarning();
            return false;
        }
    }
    return true;
}


/*!
     assures by change from one telescope view to all telescope view that all cameras appear
     in the same plotting mode

*/
void VDisplay::setColorScheme()
{
    if( fDebug )
    {
        cout << "VDisplay::setColorScheme " << fCamera.size() << endl;
    }
    // if color scheme, show in all pads
    if( fMenuOpt->IsEntryChecked( M_OPT_COL_SCHE ) )
    {
        for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
        {
            fCamera[fTelescopesToShow[i]]->setPlotColor( 1 );
        }
    }
    else if( fMenuOpt->IsEntryChecked( M_OPT_BW_SCHE ) )
    {
        for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
        {
            fCamera[fTelescopesToShow[i]]->setPlotColor( 2 );
        }
    }
    else
    {
        for( unsigned int i = 0; i < fTelescopesToShow.size(); i++ )
        {
            fCamera[fTelescopesToShow[i]]->setPlotColor( 0 );
        }
    }
}


void VDisplay::plotFADCFit( int iChanID )
{
    if( fDebug )
    {
        cout << "VDisplay::plotFADCFit() " << iChanID << endl;
    }
    
    drawFADC( true );
}


void VDisplay::dumpDeadChannels()
{
    if( fDebug )
    {
        cout << "VDisplay::dumpDeadChannels()" << endl;
    }
    
    if( fEventLoop )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            fEventLoop->setTelID( fEventLoop->getTeltoAna()[i] );
            fEventLoop->printDeadChannels();
            fEventLoop->printDeadChannels( true );
        }
    }
}

void VDisplay::dumpImageBorderPixels()
{
    if( fDebug )
    {
        cout << "VDisplay::dumpImageBorderPixels()" << endl;
    }
    
    if( fEventLoop )
    {
        for( unsigned int i = 0; i < fEventLoop->getTeltoAna().size(); i++ )
        {
            fEventLoop->setTelID( fEventLoop->getTeltoAna()[i] );
            
            multimap< double, unsigned int > i_ImagePixel;
            multimap< double, unsigned int >::iterator it_ImagePixel;
            multimap< double, unsigned int > i_BorderPixel;
            multimap< double, unsigned int >::iterator it_BorderPixel;
            for( unsigned int j = 0; j < fEventLoop->getImage().size(); j++ )
            {
                if( j >= fEventLoop->getSums().size() )
                {
                    continue;
                }
                if( j >= fEventLoop->getBorder().size() )
                {
                    continue;
                }
                if( fEventLoop->getImage()[j] )
                {
                    i_ImagePixel.insert( make_pair( fEventLoop->getSums()[j], j ) );
                }
                if( fEventLoop->getBorder()[j] )
                {
                    i_BorderPixel.insert( make_pair( fEventLoop->getSums()[j], j ) );
                }
            }
            if( i_ImagePixel.size() > 0 )
            {
                cout << "Tel " << i + 1 << ": image pixels ";
                for( it_ImagePixel = i_ImagePixel.begin(); it_ImagePixel != i_ImagePixel.end(); it_ImagePixel++ )
                {
                    cout << it_ImagePixel->second << ",";
                }
                cout << endl;
            }
            if( i_BorderPixel.size() > 0 )
            {
                cout << "Tel " << i + 1 << ": border pixels ";
                for( it_BorderPixel = i_BorderPixel.begin(); it_BorderPixel != i_BorderPixel.end(); it_BorderPixel++ )
                {
                    cout << it_BorderPixel->second << ",";
                }
                cout << endl;
            }
        }
    }
}
