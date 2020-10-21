/*! \class VCamera
    \brief camera plotting routines

*/

#include "VCamera.h"

VCamera::VCamera( unsigned int iTel, VEvndispData* iData )
{
	fDebug = iData->getDebugFlag();
	
	fTelescope = iTel;
	fData = iData;
	fData->getDetectorGeometry()->setTelID( fTelescope );
	fCanvas = 0;
	fEventCounter = 0;
	
	fmaxPlot = 0.45;                              /* 2*fmaxPlot = Radius of camera, total pad size is 1.0x1.0 */
	fmaxRad = 0.85;                               /* the tube signal is at its maximum fmaxRad * tube radius */
	fScaleMax = 0.;                               /* maximum value in tube */
	
	fPrintChannel = 0;                            /* no printing of channel numbers */
        fTubeSelected = -1;
	
	fBoolAllinOne = false;                        /* plot tube entries of all cameras into one camera */
	
	fPlotColor = 1;                               /* signals are represented by the size of the filled circle in the PMT */
	fAnaVis = true;
	fPlotPaper = fData->getRunParameter()->fPlotPaper;
	
	bFixScale = false;                            /* set fixed scale for colz plots */
	
	setUpCamera();
	
	fData->initializeDataReader();
	
	// get a nicer color palette and get some informations about colors/contours
	gStyle->SetPalette( 1 );
	gStyle->SetNumberContours( 100 );
	fncolors = gStyle->GetNumberOfColors();
	fndivz   = gStyle->GetNumberContours();
	
	fCurrentTimeSlice = -1;
	
	fFirstTelescopeToDraw = false;
}


/*!
     define colors, pmts, text, etc.
*/
void VCamera::setUpCamera()
{
	// define plot colors and styles
	fColorEmpty = 10;
	fColorImage = 2;
	fColorImageUser = 6;
	fColorBorder = 8;
	fColorEstimated = 6;
	fColorDead = 13;
	if( fPlotPaper )
	{
		fColorDead = 0;
	}
	fColorFADCTrig = 28;
	if( fPlotPaper )
	{
		fColorFADCTrig = 0;
	}
	fColorSum = 4;
	fColorTrigger = 3;
	fColorHit = 2;
	fColorTiming = 4;
	fColorPedMean = 4;
	fColorPedVar = 4;
	fFillStylePed = 1001;
	fFillStyleEmpty = 1001;
	fFillStyleDead = 1001;
	if( fPlotPaper )
	{
		fFillStyleDead = 0;
	}
	fFillStyleFADCTrig = 3144;
	if( fPlotPaper )
	{
		fFillStyleFADCTrig = 0;
	}
	fFillStylePos = 1001;
	fFillStyleNeg = 0;
	fTelescopeEllipseColor = 1;
	
	
	// get maximum distance of tubes from centre, rescale to canvas NDC
	fdist_edgeX = 0.;
	fdist_edgeY = 0.;
	for( unsigned int i = 0; i < fData->getDetectorGeo()->getNumChannels(); i++ )
	{
		if( TMath::Abs( fData->getDetectorGeo()->getX()[i] ) + fData->getDetectorGeo()->getTubeRadius()[i] > fdist_edgeX )
		{
			fdist_edgeX = TMath::Abs( fData->getDetectorGeo()->getX()[i] ) + fData->getDetectorGeo()->getTubeRadius()[i];
		}
		if( TMath::Abs( fData->getDetectorGeo()->getY()[i] ) + fData->getDetectorGeo()->getTubeRadius()[i] > fdist_edgeY )
		{
			fdist_edgeY = TMath::Abs( fData->getDetectorGeo()->getY()[i] ) + fData->getDetectorGeo()->getTubeRadius()[i];
		}
	}
	fmax_dist_edge = fData->getDetectorGeo()->getMaximumFOV_deg();

	// array with pmt data (rescaled data)
	fPMTData.resize( int( fData->getDetectorGeo()->getNumChannels() ), 0. );
	
	// rescale to canvas NDC, shift coordinate system by 0.5
	double x, y, rx, ry;
	char c_number[100];
	/*! Attention:
	    there are maybe more channels defined in the camera layout than there
	    are in the data!
	    Never loop over fgraphTubes/fgraphTubesEntry without testing the
	    radius (GetR1()>0.)
	    This is hopefully only valid for the prototype. */
	
	double iMaxDist = 0.;
	for( int i = 0; i < int( fData->getDetectorGeo()->getNumChannels() ); i++ )
	{
		// calculate new coordinates in canvas system (shifted by 0.5 to center of canvas)
		x = convertX( fData->getDetectorGeo()->getX()[i] );
		y = convertY( fData->getDetectorGeo()->getY()[i] );
		rx = convertX( fData->getDetectorGeo()->getTubeRadius()[i], 0. );
		ry = convertY( fData->getDetectorGeo()->getTubeRadius()[i], 0. );
		if( sqrt( fData->getDetectorGeo()->getX()[i]*fData->getDetectorGeo()->getX()[i]
				  + fData->getDetectorGeo()->getY()[i]*fData->getDetectorGeo()->getY()[i] ) > iMaxDist )
		{
			iMaxDist = sqrt( fData->getDetectorGeo()->getX()[i] * fData->getDetectorGeo()->getX()[i]
							 + fData->getDetectorGeo()->getY()[i] * fData->getDetectorGeo()->getY()[i] )
					   + fData->getDetectorGeo()->getTubeRadius()[i];
		}
		// PMTs (outer shell)
		fgraphTubes.push_back( new TEllipse( x, y, rx, ry ) );
		fgraphTubes.back()->SetFillColor( fColorEmpty );
		fgraphTubes.back()->SetFillStyle( 0 );
		fgraphTubes.back()->SetUniqueID( 200000 + i );
		fgraphTubes.back()->SetLineColor( 15 );
		if( fData->getDetectorGeo()->isEdgePixel()[i] )
		{
		   fgraphTubes.back()->SetLineColor( 9 );
                }
		// PMT values
		fgraphTubesEntry.push_back( new TEllipse( x, y, rx * fmaxRad * 0.5, ry * fmaxRad * 0.5 ) );
		fgraphTubesEntry.back()->SetLineColor( 10 );
		fgraphTubesEntry.back()->SetFillColor( 10 );
		fgraphTubesEntry.back()->SetFillStyle( 0 );
		fgraphTubesEntry.back()->SetUniqueID( 200000 + i );
		// channel numbering
		sprintf( c_number, "%d", i );
		fTextChannelNumber.push_back( new TText( fgraphTubes.back()->GetX1() - fgraphTubes.back()->GetR1() / 1.4,
							 fgraphTubes.back()->GetY1() - fgraphTubes.back()->GetR2() / 2., c_number ) );
		fTextChannelNumber.back()->SetTextFont( 42 );
		fTextChannelNumber.back()->SetTextSize( 0.015 );
		fTextChannelNumber.back()->SetUniqueID( 200000 + i );
		// tube numbering
		sprintf( c_number, "%d", i + 1 );
		fTextTubeNumber.push_back( new TText( fgraphTubes.back()->GetX1() - fgraphTubes.back()->GetR1() / 1.4,
						      fgraphTubes.back()->GetY1() - fgraphTubes.back()->GetR2() / 2., c_number ) );
		fTextTubeNumber.back()->SetTextFont( 42 );
		fTextTubeNumber.back()->SetTextSize( 0.015 );
		fTextTubeNumber.back()->SetUniqueID( 200000 + i );
		fTextTubeNumber.back()->SetTextColor( 2 );
		// tube markers
		fgraphMarker.push_back( new TMarker() );
		fgraphMarker.back()->SetMarkerStyle( 5 );
		fgraphMarker.back()->SetUniqueID( 200000 + i );
	}
	// size of camera
	// (from edge of outermost pixels)
	fCameraOuterEdge = new TEllipse( convertX( 0. ), convertY( 0. ), convertX( iMaxDist, 0. ), convertY( iMaxDist, 0. ) );
	fCameraOuterEdge->SetLineWidth( 2 );
	fCameraOuterEdge->SetFillStyle( 0 );
	// (from FOV entry in detector geometry)
	fCameraFOV = new TEllipse( convertX( 0. ), convertY( 0. ), convertX( fData->getDetectorGeo()->getFieldofView()[getTelescopeNumber()] / 2., 0. ),
							   convertY( fData->getDetectorGeo()->getFieldofView()[getTelescopeNumber()] / 2., 0. ) );
	fCameraFOV->SetLineWidth( 2 );
	fCameraFOV->SetLineColor( 14 );
	fCameraFOV->SetLineStyle( 2 );
	fCameraFOV->SetFillStyle( 0 );
	// ellipse representing the shower image
	fAnaEllipse = new TEllipse();
	fAnaEllipse->SetLineColor( 7 );
	fAnaEllipse->SetFillStyle( 0 );
	fAnaEllipse->SetLineWidth( 2 );
	fAnaEllipse->SetUniqueID( 2999 );
	// circles representing the muon ring radius +/- 1 sigma_radius
	fAnaEllipse1 = new TEllipse();
	fAnaEllipse1->SetLineColor( 5 );
	fAnaEllipse1->SetFillStyle( 0 );
	fAnaEllipse1->SetLineWidth( 1 );
	fAnaEllipse1->SetUniqueID( 2990 );
	fAnaEllipse2 = new TEllipse();
	fAnaEllipse2->SetLineColor( 5 );
	fAnaEllipse2->SetFillStyle( 0 );
	fAnaEllipse2->SetLineWidth( 2 );
	fAnaEllipse2->SetUniqueID( 2991 );
	
	// lines showing the image axis and rotation
	fCenterLine = new TLine( 0., 0., 0., 0. );
	fEllipseLine = new TLine( 0., 0., 0., 0. );
	// ellipse representing the shower image (log likelihood)
	fAnaEllipseLL = new TEllipse();
	fAnaEllipseLL->SetLineColor( 5 );
	fAnaEllipseLL->SetFillStyle( 0 );
	fAnaEllipseLL->SetLineWidth( 2 );
	fAnaEllipseLL->SetUniqueID( 2998 );
	// event text (top left corner) + image parameter line at bottom (last entry in vector)
	float i_TextSize = 0.025 * 0.7;
	float i_TextX = 0.01;
	fTextEvent.push_back( new TLatex( i_TextX, 0.97, "run/event number" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.91, "Max. Channels" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.88, "Num Samples" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.85, "Num Trigger" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.85, "Num Tubes" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.82, "Num LowGain" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.79, "Num Dead" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.01, "Parameters" ) );
	fTextEvent.push_back( new TLatex( i_TextX, 0.01, "ParametersLL" ) );
	// set text font and sizes
	for( unsigned int i = 0; i < fTextEvent.size(); i++ )
	{
		fTextEvent[i]->SetTextFont( 42 );
		fTextEvent[i]->SetTextSize( i_TextSize );
	}
	// MC text (top bottom corner)
	fTextMC.push_back( new TLatex( i_TextX, 0.48, "Primary: " ) );
	fTextMC.push_back( new TLatex( i_TextX, 0.45, "Energy [TeV]: " ) );
	fTextMC.push_back( new TLatex( i_TextX, 0.39, "C_{X}: C_{Y}:" ) );
	fTextMC.push_back( new TLatex( i_TextX, 0.36, "Xcos: " ) );
	fTextMC.push_back( new TLatex( i_TextX, 0.33, "Ycos: " ) );
	fTextMC.push_back( new TLatex( i_TextX, 0.30, "X_{off}: Y_{off}: " ) );
	// set text font and sizes
	for( unsigned int i = 0; i < fTextMC.size(); i++ )
	{
		fTextMC[i]->SetTextFont( 42 );
		fTextMC[i]->SetTextSize( i_TextSize );
	}
	// telescope number (lower right corner)
	fTextTelescopeN = new TLatex( 0.85, 0.39, "TX" );
	fTextTelescopeN->SetTextFont( 42 );
	fTextTelescopeN->SetTextSize( i_TextSize * 3 );
	// event numbers (plot papers )
	fTextEventPlotPaper = new TLatex( 01., 0.85, "EE" );
	fTextEventPlotPaper->SetTextFont( 42 );
	fTextEventPlotPaper->SetTextSize( i_TextSize * 2 );
	// camera scale axis (left+top)
	fCameraXaxis = new TGaxis( convertX( -1.* fdist_edgeX ) , 0.97, convertX( fdist_edgeX ), 0.97, -1.*fdist_edgeX, fdist_edgeX, 510, "+L" );
	fCameraXaxis->SetLabelSize( 0.02 );
	fCameraXaxis->SetLineColor( 42 );
	fCameraXaxis->SetLabelColor( 42 );
	fCameraYaxis = new TGaxis( 0.96, convertY( -1. * fdist_edgeY ) , 0.96, convertY( fdist_edgeY ), -1. * fdist_edgeY, fdist_edgeY, 510, "+L" );
	fCameraYaxis->SetLabelSize( 0.02 );
	fCameraYaxis->SetLineColor( 43 );
	fCameraYaxis->SetLabelColor( 43 );
	// theta2 circles for all in one
	for( int t = 0; t < 2 * ( int )iMaxDist; t++ )
	{
		fTheta2Circle.push_back( new TEllipse( 0., 0., t * 0.5 ) );
	}
}

/*!
 */
void VCamera::draw( double i_max, int iEventNumber, bool iAllinOne )
{
	if( fDebug )
	{
		cout << "VCamera::draw event " << iEventNumber << "\t" << fData->getTelID() << endl;
	}
	if( fCanvas == 0 )
	{
		return;
	}
	fScaleMax = i_max;
	fBoolAllinOne = iAllinOne;

        fData->setTelID( fData->getTelID() );
	
	fCanvas->cd();
	fCanvas->SetEditable( true );

        valarray<double >  i_va_temp( fData->getSums().size() );
	
	if( !iAllinOne )
	{
		fCanvas->Clear();
	}
	
	if( iEventNumber != 0 )                       /* first time, there is no content in the tubes */
	{
		switch( fcameraModus )
		{
			case C_CHARGE:
				if( fPlotColor )
				{
					setPMTColorScheme( fData->getSums(), getDrawingMask( 1, i_va_temp ),  false,  100., 0., "charge [d.c.]", true );
				}
				else
				{
					setPMTColorForChargeTiming();
				}
				break;
			case C_TIMING:
				if( fPlotColor )
				{
					setPMTColorScheme( fData->getSums(), getDrawingMask( 1, i_va_temp ), false,  0., fScaleMax, "charge [d.c]", true );
				}
				else
				{
					setPMTColorForChargeTiming();
				}
				break;
			case C_TRIGGER:
				setPMTColorOnOff( fData->getReader()->getFullTrigVec(), fColorTrigger, fColorTrigger, fFillStylePos );
				break;
			case C_HILO:
				setPMTColorOnOff( fData->getHiLo(), fColorTrigger, fColorTrigger, fFillStylePos );
				break;
			case C_HIT:
				setPMTColorOnOff( fData->getReader()->getFullHitVec(), fColorHit, fColorHit, fFillStylePos );
				break;
			case C_SUMWINDOW:
				setPMTColorScheme( fData->getCurrentSumWindow(), getDrawingMask( 1, i_va_temp ), false, 100., 0., "summation window", false, false );
				break;
			case C_SUMWINDOWSTART:
				setPMTColorScheme( fData->getTCorrectedSumFirst(), getDrawingMask( 1, i_va_temp ), false, 100., 0., "summation window start", false, false );
				break;
			case C_PEDMEAN:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getPeds(), getDrawingMask( 1, i_va_temp ), false, 100., 0., "charge", false );
				}
				else
				{
					setPMTColorScheme( fData->getPeds(), getDrawingMask( 1, i_va_temp ), false, 10., 25., "charge", false );
				}
				break;
			case C_PEDVAR:
			{
				if( fData->getRunParameter()->fsourcetype != 7 || fData->getReader()->hasFADCTrace() )
				{
					valarray< double > i_pedvars( 0., fData->getPedvars().size() );
					for( unsigned int ii = 0; ii < i_pedvars.size(); ii++ )
					{
						i_pedvars[ii] = fData->getPedvars( fData->getCurrentSumWindow()[ii] )[ii];
					}
					if( !bFixScale )
					{
						setPMTColorScheme( i_pedvars, getDrawingMask( 1, i_va_temp ), false, 100., 0., "charge", false );
					}
					else
					{
						setPMTColorScheme( i_pedvars, getDrawingMask( 1, i_va_temp ), false, 5., 20., "charge", false );
					}
				}
			}
			break;
			case C_GAINS:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getGains(),  getDrawingMask( 1, i_va_temp ),false, 100., 0., "gain", false );
				}
				else
				{
					setPMTColorScheme( fData->getGains(), getDrawingMask( 1, i_va_temp ), false, 0.8, 1.2, "gain", false );
				}
				break;
			case C_GAINVARS:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getGainvars(), getDrawingMask( 1, i_va_temp ), false, 100., 0., "gain", false );
				}
				else
				{
					setPMTColorScheme( fData->getGainvars(), getDrawingMask( 1, i_va_temp ), false, 0.8, 1.2, "gain", false );
				}
				break;
			case C_TOFF:
				setPMTColorScheme( fData->getTOffsets(), getDrawingMask( 1, i_va_temp ), false, -2., 2., "time offset [samples]", false );
				break;
			case C_CALTZERO:
			{
				float i_min = 1.e5;
				float i_max = 0.;
				for( unsigned int i = 0; i < fData->getAverageTZeros().size(); i++ )
				{
					if( fData->getAverageTZeros()[i] > i_max )
					{
						i_max = fData->getAverageTZeros()[i];
					}
					if( fData->getAverageTZeros()[i] > 0. && fData->getAverageTZeros()[i] < i_min )
					{
						i_min = fData->getAverageTZeros()[i];
					}
					if( i_min == i_max )
					{
						i_min = i_max - 1.;
					}
				}
				setPMTColorScheme( fData->getAverageTZeros(), getDrawingMask( 1, i_va_temp ), false, i_min, i_max, "time [samples]", false );
			}
			break;
			case C_TZERO:
                                // default: draw image/border
                                // for TTRIGGER: draw all triggered pixel
                                if ( fData->getSumWindowStart_T_method() == 3 )
                                {
                                    setPMTColorScheme( fData->getPulseTime(), 
                                          getDrawingMask( 2, fData->getPulseTime(), -998 ), 
                                          false, 100., 0., "time [samples]", false );
                                }
                                else
                                {
                                    setPMTColorScheme( fData->getPulseTime(), 
                                          getDrawingMask( 4, fData->getPulseTime(), -998 ), 
                                          false, 100., 0., "time [samples]", false );
                                }
				break;
			case C_PEDMEANLOW:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getPeds( true ), getDrawingMask( 1, i_va_temp ), false, 100., 0., "charge", false );
				}
				else
				{
					setPMTColorScheme( fData->getPeds( true ), getDrawingMask( 1, i_va_temp ), false, 10., 25., "charge", false );
				}
				break;
			case C_PEDVARLOW:
			{
				valarray< double > i_pedvars( 0., fData->getPedvars().size() );
				for( unsigned int ii = 0; ii < i_pedvars.size(); ii++ )
				{
					i_pedvars[ii] = fData->getPedvars( fData->getCurrentSumWindow()[ii], true )[ii];
				}
				if( !bFixScale )
				{
					setPMTColorScheme( i_pedvars, getDrawingMask( 1, i_va_temp ), false, 100., 0., "charge", false );
				}
				else
				{
					setPMTColorScheme( i_pedvars, getDrawingMask( 1, i_va_temp ), false, 5., 20., "charge", false );
				}
			}
			break;
			case C_GAINSLOW:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getGains( true ), getDrawingMask( 1, i_va_temp ), false, 100., 0., "gain", false );
				}
				else
				{
					setPMTColorScheme( fData->getGains( true ), getDrawingMask( 1, i_va_temp ), false, 0.8, 1.2, "gain", false );
				}
				break;
			case C_GAINVARSLOW:
				if( !bFixScale )
				{
					setPMTColorScheme( fData->getGainvars( true ), getDrawingMask( 1, i_va_temp ), false, 100., 0., "gain", false );
				}
				else
				{
					setPMTColorScheme( fData->getGainvars( true ), getDrawingMask( 1, i_va_temp ), false, 0.8, 1.2, "gain", false );
				}
				break;
			case C_TOFFLOW:
				setPMTColorScheme( fData->getTOffsets( true ), getDrawingMask( 1, i_va_temp ), false, -2., 2., "time offset [samples]", false );
				break;
			case C_CALTZEROLOW:
				setPMTColorScheme( fData->getAverageTZeros( true ), getDrawingMask( 1, i_va_temp ), false, 100., 0.,  "time (cal tzero) [samples]", false );
				break;
			case C_STATUS:
				setPMTColorScheme( fData->getDeadUI(), getDrawingMask( 3, i_va_temp ), false, 0., 12., "high gain channel status", false, true );
				break;
			case C_STATUSLOW:
				setPMTColorScheme( fData->getDeadUI( true ), getDrawingMask( 3, i_va_temp ), false, 0., 12., "low gain channel status", false, true, true );
				break;
			case C_LOWGAIN:
				setPMTColorScheme( fData->getLowGainMultiplier_Camera(), getDrawingMask( 1, i_va_temp ), false, 100., 0., "multiplier", false );
				break;
			case C_L1:
				setPMTColorScheme( fData->getL1Rates(), getDrawingMask( 0, i_va_temp ), false, 100., 0., "L1 rates", false, false, false );
				break;
			case C_HV:
				setPMTColorScheme( fData->getHV(), getDrawingMask( 0, i_va_temp ), false, 100., 0., "HV [V]", false, false, false );
				break;
			case C_CURRENTS:
				setPMTColorScheme( fData->getCurrents(), getDrawingMask( 0, i_va_temp ), false, 100., 0., "Currents", false, false, false );
				break;
			case C_TRIGGER_EVNDISP:
				setPMTColorOnOff( fData->getTrigger(), fColorTrigger, fColorTrigger, fFillStylePos );
				break;
			case C_TEMPLATE:
				if( fData->getRunParameter()->ffrogsmode == 1 )
				{
				  double minSum = 0;
				  double maxSum = 0;
				  getMinMax( fData->getSums(), minSum, maxSum, getDrawingMask( 1, i_va_temp ) );
				  setPMTColorScheme( fData->getTemplateMu(), getDrawingMask( 1, i_va_temp ), false,  minSum, maxSum, "FROGS signal [d.c.]", false );
				}
				break;
			case C_MODEL3D:
				if( fData->getRunParameter()->fUseDisplayModel3D )
				{
					double minSum = 0;
					double maxSum = 0;
					getMinMax( fData->getSums(), minSum, maxSum, getDrawingMask( 1, i_va_temp ) );
					setPMTColorScheme( fData->getModel3DMu(), getDrawingMask( 1, i_va_temp ), false,  minSum, maxSum, "Model3D signal [d.c.]", false );
					setPMTColorOff( fData->getModel3DClean() );
				}
				break;
			case C_CLUSTERID:
				setPMTColorScheme( fData->getClusterID(), getDrawingMask( 1, i_va_temp ), false, 100, 0, "Cluster ID", true, false );	
				break;
                        case C_PE:
                                setPMTColorScheme( fData->getPE(), getDrawingMask( 1, i_va_temp ), false,  1., fData->getPE().max()+5, "pe", false );
				break;
                                
			default:
				break;
		}
		// draw the tubes
		if( !fBoolAllinOne )
		{
			for( unsigned int i = 0; i < fgraphTubes.size(); i++ )
			{
				if( fgraphTubes[i]->GetR1() > 0. )
				{
					fgraphTubes[i]->Draw();
				}
			}
			// draw the tube entries
			for( unsigned int i = 0; i < fgraphTubesEntry.size(); i++ )
			{
				if( fgraphTubes[i]->GetR1() > 0. )
				{
					fgraphTubesEntry[i]->Draw();
				}
			}
		}
		else
		{
			if( fCameraOuterEdge && fFirstTelescopeToDraw )
			{
				fCameraOuterEdge->Draw();
			}
		}
		if( fCameraFOV && fFirstTelescopeToDraw )
		{
			if( fBoolAllinOne )
			{
				fCameraFOV->SetR1( convertX( fData->getDetectorGeo()->getMaximumFOV_deg(), 0. ) );
				fCameraFOV->SetR2( convertY( fData->getDetectorGeo()->getMaximumFOV_deg(), 0. ) );
			}
			else
			{
				fCameraFOV->SetR1( fData->getDetectorGeo()->getFieldofView()[getTelescopeNumber()] / 2. );
				fCameraFOV->SetR2( fData->getDetectorGeo()->getFieldofView()[getTelescopeNumber()] / 2. );
			}
			fCameraFOV->Draw();
		}
		// mark zero suppressed channels
		if( !fBoolAllinOne )
		{
			for( unsigned int i = 0; i < fgraphTubesEntry.size(); i++ )
			{
                                // only charges are available
				if( fData->getZeroSuppressed()[i] > 1 )
				{
					TMarker* iZeroSupp = new TMarker( fgraphTubes[i]->GetX1(), fgraphTubes[i]->GetY1(), 5 );
                                        if( fData->getZeroSuppressed()[i] > 2 )
                                        {
                                            iZeroSupp->SetMarkerColor( 5 );
                                        }
					iZeroSupp->Draw();
				}
			}
		}
/*		// mark pixel used for LL fit
		if( fData->getImageParameters()->Fitstat >= 0 && !fBoolAllinOne )
		{
			for( unsigned int i = 0; i < fgraphMarker.size(); i++ )
			{
				if( i < fData->getImageBorderNeighbour().size() &&  fData->getImageBorderNeighbour()[i] )
				{
					if( fgraphMarker[i] )
					{
						fgraphMarker[i]->SetMarkerStyle( 2 );
						fgraphMarker[i]->SetX( fgraphTubesEntry[i]->GetX1() );
						fgraphMarker[i]->SetY( fgraphTubesEntry[i]->GetY1() );
						fgraphMarker[i]->Draw();
					}
				}
			}

		} */
                // mark pixels with >1 pe
                // (in PE plot, mark image/border pixels)
                if( fData->getReader()->getPE().size() > 0 )
                {
			for( unsigned int i = 0; i < fgraphMarker.size(); i++ )
			{
                            if( fcameraModus != C_PE )
                            {
                                if( i < fData->getReader()->getPE().size() 
                                        && fData->getReader()->getPE()[i] > 2
                                        && fgraphMarker[i] )
                                {
                                        fgraphMarker[i]->SetMarkerStyle( 5 );
                                        fgraphMarker[i]->SetMarkerColor( 2 );
                                        fgraphMarker[i]->SetX( fgraphTubesEntry[i]->GetX1() );
                                        fgraphMarker[i]->SetY( fgraphTubesEntry[i]->GetY1() );
                                        fgraphMarker[i]->Draw();
                                }
                             }
                            else 
                            {
                                if( i < fData->getImage().size() && i < fData->getBorder().size() 
                                   && (fData->getImage()[i] || fData->getBorder()[i] )
                                   && fgraphMarker[i] )
                                {
                                        fgraphMarker[i]->SetMarkerStyle( 5 );
                                        fgraphMarker[i]->SetMarkerColor( 2 );
                                        fgraphMarker[i]->SetX( fgraphTubesEntry[i]->GetX1() );
                                        fgraphMarker[i]->SetY( fgraphTubesEntry[i]->GetY1() );
                                        fgraphMarker[i]->Draw();
                                }
                            }
                         }
                }



		// draw pixel recovered by the correlation image cleaning
		if( !fBoolAllinOne && fData->getImageCleaningParameter()->getImageCleaningMethod() == "TWOLEVELANDCORRELATION" )
		{
			for( unsigned int i = 0; i < fgraphMarker.size(); i++ )
			{
				if( i < fData->getBorderCorrelationCoefficient().size() &&  fData->getBorderCorrelationCoefficient()[i] > 1.e-2 )
				{
					if( fgraphMarker[i] )
					{
						fgraphMarker[i]->SetX( fgraphTubesEntry[i]->GetX1() );
						fgraphMarker[i]->SetY( fgraphTubesEntry[i]->GetY1() );
						fgraphMarker[i]->Draw();
					}
				}
			}
		}
		// draw tube/channel numbers
		if( fPrintChannel != 0 && fPrintChannel < 3 && !fBoolAllinOne )
		{
			for( unsigned int i = 0; i < fgraphTubes.size(); i++ )
			{
				if( fgraphTubes[i]->GetR1() > 0. && fPrintChannel == 1 )
				{
					fTextChannelNumber[i]->DrawTextNDC( fTextChannelNumber[i]->GetX(),
														fTextChannelNumber[i]->GetY(), fTextChannelNumber[i]->GetTitle() );
				}
				else if( fgraphTubes[i]->GetR1() > 0. && fPrintChannel == 2 )
				{
					fTextTubeNumber[i]->DrawTextNDC( fTextTubeNumber[i]->GetX(), fTextTubeNumber[i]->GetY(), fTextTubeNumber[i]->GetTitle() );
				}
			}
			// draw axis with camera sizes
			fCameraXaxis->Draw();
			fCameraYaxis->Draw();
		}
		else if( fPrintChannel == 3 )
		{
			// draw axis with camera sizes (coordinate system)
			fCameraXaxis->Draw();
			fCameraYaxis->Draw();
			for( unsigned int t = 0; t < fTheta2Circle.size(); t++ )
			{
				fTheta2Circle[t]->SetLineWidth( 1 );
				fTheta2Circle[t]->SetLineStyle( 1 );
				fTheta2Circle[t]->SetLineColor( 9 );
				fTheta2Circle[t]->SetX1( convertX( 0. ) );
				fTheta2Circle[t]->SetY1( convertY( 0. ) );
				fTheta2Circle[t]->SetR1( convertX( 0.5 * ( t + 1 ), 0. ) );
				fTheta2Circle[t]->SetR2( convertY( 0.5 * ( t + 1 ), 0. ) );
				fTheta2Circle[t]->SetFillStyle( 0 );
				fTheta2Circle[t]->Draw();
			}
		}
		// draw the ellipse from image parameters
		// don't draw anything for hit/timing/ped mean/ped var/gains/toffsets
		if( fData->getRunParameter()->fmuonmode == 0
	       && ( fcameraModus == C_CHARGE || fcameraModus == C_TRIGGER || fcameraModus == C_TZERO || fcameraModus == C_TEMPLATE ) )
		{
			drawAnaResults();
		}
		// is this a muon analysis?
		if( fData->getRunParameter()->fmuonmode == 1
	       && ( fcameraModus == C_CHARGE || fcameraModus == C_TRIGGER || fcameraModus == C_TZERO ) )
		{
			drawMuonResults();
		}
		
		// draw bottom line with results from image calculation
		drawEventText();
		
		// draw stars in field of view
		drawStarsInFOV();
		
	}
	//   fCanvas->Update();
	fEventCounter++;
	fCanvas->SetEditable( false );
	if( fTubeSelected > 0 && !fBoolAllinOne )
	{
		showSelectedChannel( fTubeSelected, true );
	}
	fFirstTelescopeToDraw = false;
}


/*!
    calculate radii of filled circles plotted according to the signal in each PMT

    image/border/background pixels are color coded (red/green/blue)

    dead channels are black hashed

    recovered dead channels (loglikelihood) are red hashed
*/
void VCamera::setPMTColorForChargeTiming()
{
	if( fDebug )
	{
		cout << "VCamera::setPMTColorForChargeTiming() " << fData->getTelID() << "\t" << fData->getReader()->getTelescopeID() << endl;
	}
	
	if( fPMTData.size() != fData->getSums().size() )
	{
		fPMTData.resize( fData->getSums().size(), 0. );
	}
	// copy data (will be rescaled to fit into plotted tubes)
	for( unsigned int i = 0; i < fPMTData.size(); i++ )
	{
		fPMTData[i] = fData->getSums()[i];
	}
	// rescale data to maximum values of 1
	fPMTData = rescaleSums( fPMTData , false );
	
	// all colors > 10 are greyish/brown. Start at 1 again
	int iTelescopeColor = ( fTelescope % 10 ) + 1;
	if( iTelescopeColor % 10 == 0 )
	{
		iTelescopeColor += 1;
	}
	
	// loop over all pmts and set tube radii and colors
	for( unsigned int i = 0; i < fPMTData.size(); i++ )
	{
		// image/border/background pixels are
		//   - not dead, or a dead channel recovered by the loglikelihood method
		//   - a hit
		//   - not switched of by the user
		if( ( !( fData->getDead( i, fData->getHiLo()[i] ) && !( fData->getDeadRecovered()[i] && ( fData->getImage()[i] || fData->getBorder()[i] ) ) ) || fData->getLLEst()[i] ) && fData->getImageUser()[i] != -1 )
		{
			// set tube radii
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1()*abs( fPMTData[i] ) * fmaxRad );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2()*abs( fPMTData[i] ) * fmaxRad );
			if( fPMTData[i] >= 0. )
			{
				fgraphTubesEntry[i]->SetFillStyle( fFillStylePos );
			}
			else
			{
				fgraphTubesEntry[i]->SetFillStyle( fFillStyleNeg );
			}
			// image pixels
			if( fData->getImage()[i] )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorImage );
				fgraphTubesEntry[i]->SetFillColor( fColorImage );
				if( fBoolAllinOne )
				{
					fgraphTubesEntry[i]->SetLineColor( iTelescopeColor );
					fgraphTubesEntry[i]->SetFillColor( iTelescopeColor );
				}
			}
			// border pixels
			else if( fData->getBorder()[i] )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorBorder );
				fgraphTubesEntry[i]->SetFillColor( fColorBorder );
				if( fBoolAllinOne )
				{
					fgraphTubesEntry[i]->SetLineColor( iTelescopeColor );
					fgraphTubesEntry[i]->SetFillColor( iTelescopeColor );
				}
			}
			// dead channels estimated by the loglikelihood method
			else if( fData->getDead( i, fData->getHiLo()[i] ) && ( fData->getLLEst()[i] || fData->getDeadRecovered()[i] ) )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorEstimated );
				fgraphTubesEntry[i]->SetFillColor( fColorEstimated );
				fgraphTubesEntry[i]->SetFillStyle( fFillStyleDead );
				if( fBoolAllinOne )
				{
					fgraphTubesEntry[i]->SetLineColor( iTelescopeColor );
					fgraphTubesEntry[i]->SetFillColor( iTelescopeColor );
				}
			}
			// background pixels
			else
			{
				fgraphTubesEntry[i]->SetLineColor( fColorSum );
				fgraphTubesEntry[i]->SetFillColor( fColorSum );
				if( fBoolAllinOne )
				{
					fgraphTubesEntry[i]->SetR1( 0. );
					fgraphTubesEntry[i]->SetR2( 0. );
				}
			}
		}
		// dead pixels
		else if( !fBoolAllinOne )
		{
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1() * fmaxRad );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2() * fmaxRad );
			fgraphTubesEntry[i]->SetLineColor( fColorDead );
			fgraphTubesEntry[i]->SetFillColor( fColorDead );
			fgraphTubesEntry[i]->SetFillStyle( fFillStyleDead );
			if( fData->getImageUser()[i] == -1 )
			{
				fgraphTubesEntry[i]->SetFillColor( fColorImageUser );
			}
		}
		else
		{
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1() * fmaxRad );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2() * fmaxRad );
			fgraphTubesEntry[i]->SetLineColor( 0 );
			fgraphTubesEntry[i]->SetFillColor( 0 );
			fgraphTubesEntry[i]->SetFillStyle( 0 );
		}
		
		// check if pixel is a L2 trigger signal
		for( unsigned int l = 0; l < fData->getFADCstopTrig().size(); l++ )
		{
			if( i == fData->getFADCstopTrig()[l] && !fBoolAllinOne )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillStyle( fFillStyleFADCTrig );
				break;
			}
		}
	}
}


void VCamera::drawEventText()
{
	if( fDebug )
	{
		cout << "VCamera::drawEventText()" << endl;
	}
	char iText[600];
	
	// big T1/T2... in lower right corner
	//   red means telescope has triggered
	if( fData->getReader()->hasLocalTrigger( fTelescope ) )
	{
		fTextTelescopeN->SetTextColor( 2 );
	}
	else
	{
		fTextTelescopeN->SetTextColor( 1 );
	}
	sprintf( iText, "T%d", fTelescope + 1 );
	fTextTelescopeN->SetTitle( iText );
	fTextTelescopeN->SetNDC( true );
	if( !fBoolAllinOne )
	{
		fTextTelescopeN->DrawLatex( 0.85, 0.1, fTextTelescopeN->GetTitle() );
	}
	else
	{
		sprintf( iText, "All" );
		fTextTelescopeN->SetTitle( iText );
		//        fTextTelescopeN->DrawLatex( 0.85, 0.1, fTextTelescopeN->GetTitle() );
	}
	
	// muon text in upper right corner
	if( fData->getRunParameter()->fmuonmode == 1 )
	{
		sprintf( iText, "Detected muon  %d", int( fData->getImageParameters()->muonValid ) );
		fTextEvent[0]->SetTitle( iText );
		sprintf( iText, "Muon radius %.2f +-%.2f", fData->getImageParameters()->muonRadius, fData->getImageParameters()->muonRSigma );
		fTextEvent[1]->SetTitle( iText );
		sprintf( iText, "Muon size %.2f", fData->getImageParameters()->muonSize );
		fTextEvent[2]->SetTitle( iText );
		
		float i_TextX = 0.72;
		float i_TextY = 0.97;
		float i_TextdY = 0.03;
		unsigned int i_toShow = 3;
		// show only event number for fBoolAllinOne
		if( fBoolAllinOne && fTelescope == 0 )
		{
			i_toShow = 3;
		}
		else if( fBoolAllinOne )
		{
			i_toShow = 1;
		}
		for( unsigned int i = 0; i < i_toShow; i++ )
		{
			fTextEvent[i]->SetNDC( true );
			if( !fPlotPaper )
			{
				fTextEvent[i]->DrawLatex( i_TextX, i_TextY -= i_TextdY, fTextEvent[i]->GetTitle() );
			}
		}
	}
	// big letters for plotpaper options
	if( fPlotPaper && fTelescope == fData->getTeltoAna()[0] )
	{
		sprintf( iText, "Run: %d Event: %d", fData->getRunNumber(), int( fData->getReader()->getEventNumber() ) );
		if( fCurrentTimeSlice >= 0 )
		{
			sprintf( iText, "%s FADC %d", iText, fCurrentTimeSlice );
		}
		fTextEventPlotPaper->SetNDC( true );
		fTextEventPlotPaper->SetTitle( iText );
		fTextEventPlotPaper->DrawLatex( 0.02, 0.95, fTextEventPlotPaper->GetTitle() );
	}
	
	// event text in upper left corner
	// no GPS times for MC
	if( fData->getReader()->isMC() )
	{
		sprintf( iText, "Run: %d Event: %d  Type: %d (%d) Trig: %d", fData->getRunNumber(), int( fData->getReader()->getEventNumber() ),
				 int( fData->getReader()->getNewEventType() ),
				 int( fData->getReader()->getATEventType() ),
				 int( fData->getReader()->getLocalTriggerType( fData->getReader()->getTelescopeID() ) ) );
	}
#ifndef NOVBF
	// GPS times are expected only for VBF data files
	else
	{
		VGPSDecoder GPSDecoder;
		GPSDecoder.decode( fData->getReader()->getGPS0(), fData->getReader()->getGPS1(), fData->getReader()->getGPS2(), fData->getReader()->getGPS3(), fData->getReader()->getGPS4() );
		sprintf( iText, "Run: %d Event: %d  Type: %d (%d) GPS: %d %d : %d : %d : %.5f", fData->getRunNumber(), int( fData->getReader()->getEventNumber() ), int( fData->getReader()->getNewEventType() ), int( fData->getReader()->getATEventType() ), int( fData->getReader()->getATGPSYear() + 2000 ), GPSDecoder.getDays(), GPSDecoder.getHrs(), GPSDecoder.getMins(), GPSDecoder.getSecs() );
	}
#endif
	fTextEvent[0]->SetTitle( iText );
	// get local trigger list
	if( fBoolAllinOne )
	{
		sprintf( iText, "local trigger: " );
		for( unsigned int t = 0; t < fData->getNTel(); t++ )
		{
			if( fData->getReader()->hasLocalTrigger( t ) )
			{
				sprintf( iText, "%s %d", iText, t + 1 );
			}
		}
	}
	else
	{
		sprintf( iText, "Max channel %d", int( fData->getReader()->getMaxChannels() ) );
	}
	fTextEvent[1]->SetTitle( iText );
	sprintf( iText, "Num Samples %d", int( fData->getNSamples() ) );
	fTextEvent[2]->SetTitle( iText );
	sprintf( iText, "Num Trigger %d", fData->getReader()->getNumberofFullTrigger() );
	fTextEvent[3]->SetTitle( iText );
	sprintf( iText, "Num Tubes %d", fData->getImageParameters()->ntubes );
	fTextEvent[4]->SetTitle( iText );
	sprintf( iText, "Num LowGain %d", fData->getImageParameters()->nlowgain );
	fTextEvent[5]->SetTitle( iText );
	int i_numtrig = 0;
	int i_numtrig2 = 0;
	for( unsigned int i = 0; i < fData->getDead().size(); i++ ) if( fData->getDead( i, false ) )
		{
			i_numtrig++;
		}
	for( unsigned int i = 0; i < fData->getDead( true ).size(); i++ ) if( fData->getDead( i, true ) )
		{
			i_numtrig2++;
		}
	sprintf( iText, "Num Dead %d/%d", i_numtrig, i_numtrig2 );
	fTextEvent[6]->SetTitle( iText );
	
	if( fData->getImageParameters()->ntubes > 0 )
	{
		sprintf( iText, "MLL: c_x=%.2f,c_y=%.2f,dist=%.2f,length=%.3f,width=%.3f,#alpha=%.2f,size=%.0f,Fitstat=%d",
				 fData->getImageParametersLogL()->cen_x,
				 fData->getImageParametersLogL()->cen_y,
				 fData->getImageParametersLogL()->dist,
				 fData->getImageParametersLogL()->length,
				 fData->getImageParametersLogL()->width,
				 fData->getImageParametersLogL()->alpha,
				 fData->getImageParametersLogL()->size,
				 fData->getImageParametersLogL()->Fitstat );
		fTextEvent[fTextEvent.size() - 1]->SetTitle( iText );
		sprintf( iText, "GEO: c_x=%.2f,c_y=%.2f,dist=%.2f,length=%.3f,width=%.3f,size=%.0f/%.0f,loss=%.2f,lossDead=%.2f,tgrad=%.2f, fui=%.2f",
				 fData->getImageParameters()->cen_x,
				 fData->getImageParameters()->cen_y,
				 fData->getImageParameters()->dist,
				 fData->getImageParameters()->length,
				 fData->getImageParameters()->width,
				 fData->getImageParameters()->size,
				 fData->getImageParameters()->size2,
				 fData->getImageParameters()->loss,
				 fData->getImageParameters()->lossAndDead,
				 fData->getImageParameters()->tgrad_x,
                 fData->getImageParameters()->fui );
		fTextEvent[fTextEvent.size() - 2]->SetTitle( iText );
	}
	
	float i_TextX = 0.01;
	float i_TextY = 1.00;
	float i_TextdY = 0.03;
	unsigned int i_toShow = fTextEvent.size() - 2;
	// show only event number for fBoolAllinOne
	if( fBoolAllinOne && fTelescope == 0 )
	{
		i_toShow = 2;
	}
	else if( fBoolAllinOne )
	{
		i_toShow = 1;
	}
	for( unsigned int i = 0; i < i_toShow; i++ )
	{
		fTextEvent[i]->SetNDC( true );
		if( !fPlotPaper )
		{
			fTextEvent[i]->DrawLatex( i_TextX, i_TextY -= i_TextdY, fTextEvent[i]->GetTitle() );
		}
	}
	// draw results at bottom of canvas
	if( fData->getImageParameters()->ntubes > 0 && !fBoolAllinOne )
	{
		fTextEvent[fTextEvent.size() - 2]->SetNDC( true );
		if( !fPlotPaper )
		{
			fTextEvent[fTextEvent.size() - 2]->DrawLatex( i_TextX, 0.01, fTextEvent[fTextEvent.size() - 2]->GetTitle() );
		}
		if( fData->getRunParameter()->fImageLL )
		{
			fTextEvent[fTextEvent.size() - 1]->SetNDC( true );
			if( !fPlotPaper )
			{
				fTextEvent[fTextEvent.size() - 1]->DrawLatex( i_TextX, 0.04, fTextEvent[fTextEvent.size() - 1]->GetTitle() );
			}
		}
	}
	
	// draw MC infos (if fBoolAllinOne only for first telescope)
	if( fData->getReader()->isMC() && fTextMC.size() > 0 && ( !fBoolAllinOne || fTelescope == 0 ) )
	{
		sprintf( iText, "Primary: %d", fData->getShowerParameters()->MCprimary );
		fTextMC[0]->SetTitle( iText );
		sprintf( iText, "Energy [TeV]: %.2f", fData->getShowerParameters()->MCenergy );
		fTextMC[1]->SetTitle( iText );
		sprintf( iText, "C_{X}: %.2f C_{Y}: %.2f", fData->getShowerParameters()->MCxcore, fData->getShowerParameters()->MCycore );
		fTextMC[2]->SetTitle( iText );
		sprintf( iText, "Xcos: %.3f (Ze: %.2f)", fData->getShowerParameters()->MCxcos, fData->getShowerParameters()->MCze );
		fTextMC[3]->SetTitle( iText );
		sprintf( iText, "Ycos: %.3f (Az: %.2f)", fData->getShowerParameters()->MCycos, fData->getShowerParameters()->MCaz );
		fTextMC[4]->SetTitle( iText );
		sprintf( iText, "X_{off}: %.3f Y_{off}: %.3f", fData->getShowerParameters()->MCTel_Xoff, fData->getShowerParameters()->MCTel_Yoff );
		fTextMC[5]->SetTitle( iText );
		i_TextY = 0.25;
		for( unsigned int i = 0; i < fTextMC.size(); i++ )
		{
			fTextMC[i]->SetNDC( true );
			if( !fPlotPaper )
			{
				fTextMC[i]->DrawLatex( i_TextX, i_TextY -= i_TextdY, fTextMC[i]->GetTitle() );
			}
		}
	}
}

/*
 * return a mask for pixels to be drawn
 *
 * iSelection Mask:
 *  
 *  0 - draw all channels
 *  1 - draw all non-dead channels
 *  2 - draw all values above min value
 *  3 - draw all dead channels
 *  4 - draw all image border channels
 *
 *  Drawing:
 *  0 - don't draw (dead channel color)
 *  1 - draw
 *  2 - don't draw (white)
 *
 */
vector< unsigned int > VCamera::getDrawingMask( unsigned int iSelectionMask, valarray<double> i_data, double iMinValue, bool iLowGain )
{
       if( iSelectionMask == 0 )
       {
           vector< unsigned int > iS( i_data.size(), 1 );
           return iS;
       }

       vector< unsigned int > iS( i_data.size(), 0 );
       if( iSelectionMask == 1 )
       {
            for( unsigned int i = 0; i < i_data.size(); i++ )
            {
                if( fData->getDead( i, fData->getHiLo()[i] ) )
                {
                    iS[i] = 0;
                }
                else if( fData->getImageUser()[i] == -1 )
                {
                   iS[i] = 0;
                } 
                else
                {
                   iS[i] = 1;
                }
             }
       }
       else if( iSelectionMask == 2 )
       {
             for( unsigned int i = 0; i < i_data.size(); i++ )
             {
                  if( i_data[i] > iMinValue )
                  {
                       iS[i] = 1;
                  }
                  else
                  {
                       iS[i] = 2;
                  }
              }
       }
       else if( iSelectionMask == 3 )
       {
            for( unsigned int i = 0; i < i_data.size(); i++ )
            {
                 if( fData->getDead( i, fData->getHiLo()[i] ) )
                 {
                      iS[i] = 0;
                 }
             }
       }
       else if( iSelectionMask == 4 )
       {
            for( unsigned int i = 0; i < i_data.size(); i++ )
            {
                 if( fData->getDead( i, fData->getHiLo()[i] ) )
                 {
                      iS[i] = 0;
                 }
                 if( fData->getImage()[i] || fData->getBorder()[i] )
                 {
                      iS[i] = 1;
                 }
             }
       }
       return iS;
}
             


void VCamera::setPMTColorScheme( vector< float > v_value, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle,
								 bool i_scale, bool i_DrawDead, bool iLowGain )
{
	valarray<double> it( v_value.size() );
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
		it[i] = ( double )( v_value[i] );
	}
	
	setPMTColorScheme( it, iPixelDrawMask, i_select, zmin, zmax, i_axisTitle, i_scale, i_DrawDead, iLowGain );
}

void VCamera::setPMTColorScheme( vector<unsigned int> v_value, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle, bool i_scale, bool i_DrawDead, bool iLowGain )
{
	valarray<double> it( v_value.size() );
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
		it[i] = ( double )( v_value[i] + 0.5 );
	}
	
	setPMTColorScheme( it, iPixelDrawMask, i_select, zmin, zmax, i_axisTitle, i_scale, i_DrawDead, iLowGain );
}


void VCamera::setPMTColorScheme( valarray<unsigned int> v_value, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle, bool i_scale, bool i_DrawDead, bool iLowGain )
{
	valarray<double> it( v_value.size() );
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
		it[i] = ( double )( v_value[i] + 0.5 );
	}
	
	setPMTColorScheme( it, iPixelDrawMask, i_select, zmin, zmax, i_axisTitle, i_scale, i_DrawDead, iLowGain );
}


/*!
   different tubes values are represented by a color scheme

   image pixels have full radius (defined in fmaxRad)
   border pixels have 80% radius
   others have 50% radius

   \par v_value  array with values to plot
   \par i_select timing plot (min/max values predefined)
   \par zmin     minimum value of color scheme axis
   \par zmax     maximum value of color scheme axis (if zmin=100.&&zmax=0. then min/max values are taken from v_value)
   \par i_axisTitle title of color axis
   \par i_scale  scale size of circle draw according to type (image, border, etc)
   \par i_DrawDead  draw contents of dead channels

*/
void VCamera::setPMTColorScheme( valarray<double> v_value, vector< unsigned int > iPixelDrawMask, bool i_select, double zmin, double zmax, string i_axisTitle,
				 bool i_scale, bool i_DrawDead, bool iLowGain )
{
	if( fDebug )
	{
		cout << "VCamera::setPMTColorScheme" << endl;
	}
	// do not draw PMTs for allinone mode
	if( fBoolAllinOne )
	{
		return;
	}
	// set up the color scheme
	// this is some kind of a mess, but it works like that:
	// - get maximum/minum values in v_value (or take it from input parameters zmin/zmax)
	// -
	TBox* iBox = new TBox();
	double wlmin, wlmax;
	
	// channel status takes only 12 values -> only 12 colours
	if( i_DrawDead )
	{
		gStyle->SetNumberContours( 12 );
		fncolors = gStyle->GetNumberOfColors();
		fndivz   = gStyle->GetNumberContours();
	}
	else
	{
		gStyle->SetNumberContours( 100 );
		fncolors = gStyle->GetNumberOfColors();
		fndivz   = gStyle->GetNumberContours();
	}
	//////////////////////////////////////////////////////////////
	// timing plot
	if( i_select )
	{
		double i_min = 64.;                       // maximal number of samples (hard wired number of samples??)
		double i_max = 0.;                        // minimum number of samples
		for( unsigned int i = 0; i < v_value.size(); i++ )
		{
			if( fData->getBorder()[i] || fData->getImage()[i] )
			{
				if( v_value[i] > i_max && v_value[i] < 64. )
				{
					i_max = v_value[i];
				}
				if( v_value[i] < i_min && v_value[i] > 0. )
				{
					i_min = v_value[i];
				}
			}
		}
		// get maxima from input parameters (zmax=0.) or from above
		if( zmax == 0. )
		{
			wlmax = i_max * 1.005;
		}
		else
		{
			wlmax = zmax;
		}
		// get minima from input parameters (zmin=0.) or from above
		if( zmin == 100. )
		{
			wlmin = i_min * 0.995;
		}
		else
		{
			wlmin = zmin;
		}
	}
	//////////////////////////////////////////////////////////////
	else
	{
		// get maximum and minima in v_value
		if( zmax == 0. && zmin == 100. )
		{
			getMinMax( v_value, wlmin, wlmax, iPixelDrawMask );
			wlmax *= 1.02;
                        if( wlmin > 0. ) wlmin *= 0.98;
                        else             wlmin *= 1.02;
		}
		else
		{
			wlmax = zmax;
			wlmin = zmin;
		}
	}
	// if all values equals zero set a range of +-0.1
	if( wlmax < wlmin )
	{
		wlmax = 10.;
		wlmin = 0.;
	}
	if( fabs( wlmin ) < 1.e-4 && fabs( wlmax ) < 1.e-4 )
	{
		wlmin = -0.1;
		wlmax = 0.1;
	}
	double wls  = wlmax - wlmin;
	double scale = double( fndivz ) / wls;
	int color, theColor;
	double w1, w2;
	// coordinates for zaxis
	double x1 = 0.9;
	double x2 = 0.93;
	double ymin = 0.7;
	double ymax = 0.95;
	double y1, y2;
	
	// plot axis on the right side
	for( int i = 0; i < fndivz; i++ )
	{
		w1 = wlmin + wls / fndivz * i;
		if( w1 < wlmin )
		{
			w1 = wlmin;
		}
		w2 = wlmax;
		if( i < fndivz - 1 )
		{
			w2 = wlmin + wls / fndivz * ( i + 1 );
		}
		if( w2 <= wlmin )
		{
			continue;
		}
		y1 = ymin + ( w1 - wlmin ) * ( ymax - ymin ) / wls;
		y2 = ymin + ( w2 - wlmin ) * ( ymax - ymin ) / wls;
		color = int( 0.01 + ( w1 - wlmin ) * scale );
		theColor = int( ( color + 0.99 ) * float( fncolors ) / float( fndivz ) );
		iBox->SetFillColor( gStyle->GetColorPalette( theColor ) );
		iBox->DrawBox( x1, y1, x2, y2 );
	}
	fColourAxis = new TGaxis( x2, ymin, x2, ymax, wlmin, wlmax, 10, "+L" );
	fColourAxis->SetTitle( i_axisTitle.c_str() );
	fColourAxis->SetLabelSize( 0.02 );
	fColourAxis->SetTitleSize( 0.02 );
	fColourAxis->SetTitleOffset( 1.6 );
	fColourAxis->Draw();

	// now assign colours the PMTs
	//  image pixels have full radius
	//  border pixels have 60% radius
	//  other pixels have 20% radius
	double scaler = 1.;
	//////////////////////////////////////////////
	// expect length of v_value to be npixel
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
                if( iPixelDrawMask[i] )
		{
			if( fData->getImage()[i] )
			{
				scaler = 1.;
			}
			else if( fData->getBorder()[i] )
			{
				scaler = 0.6;
			}
			else
			{
				if( i_select )
				{
					scaler = 0.;
				}
				else
				{
					scaler = 0.20;
				}
			}
			if( !i_scale )
			{
				scaler = 1.;
			}
			w1 = v_value[i];
                        if(!i_scale && w1 < wlmin )
                        {
                            scaler = 0.06;
                        }
			
			color = int( 0.01 + ( w1 - wlmin ) * scale );
			theColor = int( ( color + 0.99 ) * float( fncolors ) / float( fndivz ) );
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1() * fmaxRad * scaler );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2() * fmaxRad * scaler );
			fgraphTubesEntry[i]->SetLineColor( gStyle->GetColorPalette( theColor ) );
			fgraphTubesEntry[i]->SetFillColor( gStyle->GetColorPalette( theColor ) );
			fgraphTubesEntry[i]->SetFillStyle( 1001 );
			if( fData->getDeadRecovered( iLowGain )[i] )
			{
				fgraphTubesEntry[i]->SetFillStyle( 3013 );
			}
		}
		// dead channels
                else if( iPixelDrawMask[i] == 0  )
		// else if( !fBoolAllinOne )
		{
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1() * fmaxRad );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2() * fmaxRad );
			fgraphTubesEntry[i]->SetLineColor( fColorDead );
			fgraphTubesEntry[i]->SetFillColor( fColorDead );
			fgraphTubesEntry[i]->SetFillStyle( fFillStyleDead );
			if( fData->getImageUser()[i] == -1 )
			{
				fgraphTubesEntry[i]->SetFillColor( fColorImageUser );
			}
		}
		for( unsigned int l = 0; l < fData->getFADCstopTrig().size(); l++ )
		{
			if( i == fData->getFADCstopTrig()[l] )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillStyle( fFillStyleFADCTrig );
				break;
			}
		}
	}
}



/*!
   PMT are filled according to values in vector of bool

   \par v_value vector of bool with value for each PMT
   \par iColor  line color of PMT
   \par iFillColor fill color of PMT
   \par iFillStyle fill style of PMT

*/
void VCamera::setPMTColorOnOff( const vector<bool>& v_value, int iColor, int iFillColor, int iFillStyle )
{
	if( fDebug )
	{
		cout << "VCamera::setPMTColorOnOff" << endl;
	}
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
		// (preli, GM) don't know why FullTrigVec is sometimes larger than number of channels
		if( i >= fgraphTubesEntry.size() )
		{
			break;
		}
		// PMTs which are not that and a hit
		if( !fData->getDead( i, fData->getHiLo()[i] ) )
		{
			// true: draw filled circle according to the fill color/style
			if( v_value[i] )
			{
				fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1()* fmaxRad );
				fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2()* fmaxRad );
				fgraphTubesEntry[i]->SetLineColor( iColor );
				fgraphTubesEntry[i]->SetFillColor( iFillColor );
				fgraphTubesEntry[i]->SetFillStyle( iFillStyle );
			}
			// false: draw empty circle
			else
			{
				fgraphTubesEntry[i]->SetFillColor( 10 );
				fgraphTubesEntry[i]->SetLineColor( 10 );
			}
		}
		// dead channels
		else if( !fBoolAllinOne )
		{
			fgraphTubesEntry[i]->SetR1( fgraphTubes[i]->GetR1() * fmaxRad );
			fgraphTubesEntry[i]->SetR2( fgraphTubes[i]->GetR2() * fmaxRad );
			fgraphTubesEntry[i]->SetLineColor( fColorDead );
			fgraphTubesEntry[i]->SetFillColor( fColorDead );
			fgraphTubesEntry[i]->SetFillStyle( fFillStyleDead );
		}
		for( unsigned int l = 0; l < fData->getFADCstopTrig().size(); l++ )
		{
			if( i == fData->getFADCstopTrig()[l] )
			{
				fgraphTubesEntry[i]->SetLineColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillColor( fColorFADCTrig );
				fgraphTubesEntry[i]->SetFillStyle( fFillStyleFADCTrig );
				break;
			}
		}
	}
}

void VCamera::setPMTColorOff( const vector<bool>& v_value )
{
	for( unsigned int i = 0; i < v_value.size(); i++ )
	{
		if( i >= fgraphTubesEntry.size() )
		{
			break;
		}
		if( !fData->getDead( i, fData->getHiLo()[i] ) )
		{
			if( ! v_value[i] )   // draw empty circle
			{
				fgraphTubesEntry[i]->SetFillColor( 10 );
				fgraphTubesEntry[i]->SetLineColor( 10 );
			}
		}
	}
}

/*  draw muon ring
     double muonX0;              //!< center of muon ring X-coord
     double muonY0;              //!< center of muon ring Y-coord
     double muonRadius;          //!< radius of muon ring
     double muonRSigma;          //!< std. dev. of radius of muon ring
     double muonSize;            //!< total amount of light in muon ring
     int    muonValid;           //!< 0/1 depending on wether it satisfies criteria

*/
void VCamera::drawMuonResults()
{
	if( fData->getImageParameters()->ntubes > 0 && fData->getImageParameters()->muonRadius > 0. )
	{
		// draw the muon ring
		// transform to local pad coordinates
		fAnaEllipse->SetX1( convertX( fData->getImageParameters()->muonX0 ) );
		fAnaEllipse->SetY1( convertY( fData->getImageParameters()->muonY0 ) );
		fAnaEllipse->SetR1( convertX( fData->getImageParameters()->muonRadius, 0. ) );
		fAnaEllipse->SetR2( convertY( fData->getImageParameters()->muonRadius, 0. ) );
		fAnaEllipse->SetTheta( 0 );
		fAnaEllipse->Draw();
		
		fAnaEllipse1->SetX1( convertX( fData->getImageParameters()->muonX0 ) );
		fAnaEllipse1->SetY1( convertY( fData->getImageParameters()->muonY0 ) );
		fAnaEllipse1->SetR1( convertX( fData->getImageParameters()->muonRadius -
									   fData->getImageParameters()->muonRSigma, 0. ) );
		fAnaEllipse1->SetR2( convertY( fData->getImageParameters()->muonRadius -
									   fData->getImageParameters()->muonRSigma, 0. ) );
		fAnaEllipse1->SetTheta( 0 );
		fAnaEllipse1->Draw();
		
		fAnaEllipse2->SetX1( convertX( fData->getImageParameters()->muonX0 ) );
		fAnaEllipse2->SetY1( convertY( fData->getImageParameters()->muonY0 ) );
		fAnaEllipse2->SetR1( convertX( fData->getImageParameters()->muonRadius +
									   fData->getImageParameters()->muonRSigma, 0. ) );
		fAnaEllipse2->SetR2( convertY( fData->getImageParameters()->muonRadius +
									   fData->getImageParameters()->muonRSigma, 0. ) );
		fAnaEllipse2->SetTheta( 0 );
		fAnaEllipse2->Draw();
		
	}
}


/*

     draw results of image parameterisation

*/
void VCamera::drawAnaResults()
{
	unsigned int iMethod = 0;
	if( fData->getRunParameter()->fPlotAllInOneMethod < fData->getShowerParameters()->fNMethods )
	{
		iMethod = fData->getRunParameter()->fPlotAllInOneMethod;
	}
	// don't plot the ellipse and image line if the image was not used in the array reconstruction
	// (applies for 'all in one' only
	if( fBoolAllinOne && fData->getTelID() < fData->getShowerParameters()->fTelIDImageSelected[iMethod].size() )
	{
		if( !fData->getShowerParameters()->fTelIDImageSelected[iMethod][fData->getTelID()] )
		{
			return;
		}
	}
	
	if( fBoolAllinOne )
	{
		fTelescopeEllipseColor = 0;
		for( unsigned int i = 0; i < fData->getShowerParameters()->fTelIDImageSelected[iMethod].size(); i++ )
		{
			if( fData->getShowerParameters()->fTelIDImageSelected[iMethod][i] )
			{
				fTelescopeEllipseColor++;
			}
			if( fData->getTelID() == i )
			{
				break;
			}
		}
		if( fTelescopeEllipseColor >= 10 )
		{
			fTelescopeEllipseColor++;
		}
	}
	else
	{
		fTelescopeEllipseColor = 1;
	}
	
	if( fData->getImageParameters()->ntubes > 0 && fData->getRunParameter()->fTargetName != "laser" && !fData->getRunParameter()->fPlotRaw )
	{
		// draw the analysis ellipse
		// transform to local pad coordinates
		fAnaEllipse->SetX1( convertX( fData->getImageParameters()->cen_x ) );
		fAnaEllipse->SetY1( convertY( fData->getImageParameters()->cen_y ) );
		fAnaEllipse->SetR1( convertX( fData->getImageParameters()->length, 0. ) );
		fAnaEllipse->SetR2( convertY( fData->getImageParameters()->width, 0. ) );
		fAnaEllipse->SetTheta( fData->getImageParameters()->phi * 180. / TMath::Pi() );
		if( fBoolAllinOne )
		{
			fAnaEllipse->SetLineColor( fTelescopeEllipseColor );
		}
		else
		{
			if( fData->getImageParameters()->Fitstat < 0. )
			{
				fAnaEllipse->SetLineColor( 7 );
			}
			else
			{
				fAnaEllipse->SetLineColor( 5 );
			}
		}
		fAnaEllipse->Draw();
		fAnaEllipse1->SetX1( convertX( fData->getImageParameters()->cen_x ) );
		fAnaEllipse1->SetY1( convertY( fData->getImageParameters()->cen_y ) );
		fAnaEllipse1->SetR1( convertX( 2.*fData->getImageParameters()->length, 0. ) );
		fAnaEllipse1->SetR2( convertY( 2.*fData->getImageParameters()->width, 0. ) );
		fAnaEllipse1->SetTheta( fData->getImageParameters()->phi * 180. / TMath::Pi() );
		if( fBoolAllinOne )
		{
			fAnaEllipse1->SetLineColor( fTelescopeEllipseColor );
		}
		else
		{
			if( fData->getImageParameters()->Fitstat < 0. )
			{
				fAnaEllipse1->SetLineColor( 7 );
			}
			else
			{
				fAnaEllipse1->SetLineColor( 5 );
			}
		}
		fAnaEllipse1->SetLineStyle( 2 );
		fAnaEllipse1->Draw();
		// draw telescope numbers
		if( fBoolAllinOne )
		{
			char hname[50];
			sprintf( hname, "T%d", fData->getTelID() + 1 );
			TText* iT = new TText( convertX( fData->getImageParameters()->cen_x ), convertY( fData->getImageParameters()->cen_y ),
								   hname );
			iT->SetTextFont( 42 );
			iT->SetTextColor( fTelescopeEllipseColor );
			iT->SetTextSize( 0.4 * iT->GetTextSize() );
			iT->Draw();
		}
		// draw the analysis ellipse after loglikelihood recovering of dead channels
		if( fData->getRunParameter()->fImageLL )
		{
			fAnaEllipseLL->SetX1( convertX( fData->getImageParametersLogL()->cen_x ) );
			fAnaEllipseLL->SetY1( convertY( fData->getImageParametersLogL()->cen_y ) );
			fAnaEllipseLL->SetR1( convertX( fData->getImageParametersLogL()->length, 0. ) );
			fAnaEllipseLL->SetR2( convertY( fData->getImageParametersLogL()->width, 0. ) );
			fAnaEllipseLL->SetTheta( fData->getImageParametersLogL()->phi * 180. / TMath::Pi() );
			// draw different line style if fit didn't worked well
			if( fData->getImageParametersLogL()->Fitstat < 3 )
			{
				fAnaEllipseLL->SetLineStyle( 2 );
			}
			else
			{
				fAnaEllipseLL->SetLineStyle( 1 );
			}
			if( fBoolAllinOne )
			{
				fAnaEllipseLL->SetLineColor( fTelescopeEllipseColor );
			}
			else
			{
				fAnaEllipseLL->SetLineColor( 5 );
			}
			fAnaEllipseLL->Draw();
		}
		// draw a line from the center through the long axis of the ellipse
		if( fAnaVis )
		{
			if( !fBoolAllinOne )
			{
				fCenterLine->SetX1( convertX( 0. ) );
				fCenterLine->SetY1( convertY( 0. ) );
				fCenterLine->SetX2( convertX( fData->getImageParameters()->cen_x ) );
				fCenterLine->SetY2( convertY( fData->getImageParameters()->cen_y ) );
				fCenterLine->SetLineColor( 1 );
				fCenterLine->SetLineStyle( 4 );
				fCenterLine->SetLineWidth( 2 );
				fCenterLine->Draw();
			}
			double i_scale1 = 2.;
			double i_scale2 = 2.;
			if( fBoolAllinOne )
			{
				i_scale1 = 30.;
				i_scale2 = 30.;
			}
			double i_x1 = convertX( fData->getImageParameters()->cen_x )
						  + i_scale1 * fAnaEllipse->GetR1() * cos( fAnaEllipse->GetTheta() * TMath::Pi() / 180. );
			double i_x2 = convertX( fData->getImageParameters()->cen_x )
						  - i_scale2 * fAnaEllipse->GetR1() * cos( fAnaEllipse->GetTheta() * TMath::Pi() / 180. );
			double i_y1 = convertY( fData->getImageParameters()->cen_y )
						  + i_scale1 * fAnaEllipse->GetR1() * sin( fAnaEllipse->GetTheta() * TMath::Pi() / 180. );
			double i_y2 = convertY( fData->getImageParameters()->cen_y )
						  - i_scale2 * fAnaEllipse->GetR1() * sin( fAnaEllipse->GetTheta() * TMath::Pi() / 180. );
						  
			fEllipseLine->SetX1( i_x1 );
			fEllipseLine->SetY1( i_y1 );
			fEllipseLine->SetX2( i_x2 );
			fEllipseLine->SetY2( i_y2 );
			if( fBoolAllinOne )
			{
				fEllipseLine->SetLineColor( fTelescopeEllipseColor );
			}
			else
			{
				fEllipseLine->SetLineColor( 1 );
			}
			if( !fBoolAllinOne )
			{
				fEllipseLine->SetLineStyle( 1 );
			}
			else
			{
				fEllipseLine->SetLineStyle( 1 );
			}
			fEllipseLine->SetLineWidth( 1 );
			fEllipseLine->Draw();
			// draw image centroids
			if( !fBoolAllinOne )
			{
				fAnaShowerCentroid = new TMarker( convertX( fData->getImageParameters()->cen_x ),
								  convertY( fData->getImageParameters()->cen_y ), 3 );
				fAnaShowerCentroid->Draw();
			}
			// draw reconstructed shower direction
			// require successfull reconstruction
			if( fData->getShowerParameters()->fShower_Chi2[iMethod] >= 0 )
			{
				if( fData->getDetectorGeo()->getGrIsuVersion() >= 412 )
				{
					fAnaShowerDir = new TMarker( convertX( fData->getShowerParameters()->fShower_Xoffset[iMethod] ),
								     convertY( -1.*fData->getShowerParameters()->fShower_Yoffset[iMethod] ), 29 );
				}
				else
				{
					fAnaShowerDir = new TMarker( convertX( fData->getShowerParameters()->fShower_Xoffset[iMethod] ),
								     convertY( fData->getShowerParameters()->fShower_Yoffset[iMethod] ), 29 );
				}
				fAnaShowerDir->SetMarkerColor( 6 );
				fAnaShowerDir->SetMarkerSize( 2. );
				fAnaShowerDir->SetMarkerStyle( 29 );
				fAnaShowerDir->Draw();
				// draw markers from disp methods
				if( fBoolAllinOne && iMethod < fData->getShowerParameters()->fShower_Xoff_DISP.size() )
				{
					// first draw marker for std method (ID0)
					if( fData->getDetectorGeo()->getGrIsuVersion() >= 412 )
					{
						fAnaShowerDir = new TMarker( convertX( fData->getShowerParameters()->fShower_Xoffset[0] ),
									     convertY( -1.*fData->getShowerParameters()->fShower_Yoffset[0] ), 29 );
					}
					else
					{
						fAnaShowerDir = new TMarker( convertX( fData->getShowerParameters()->fShower_Xoffset[0] ),
									     convertY( fData->getShowerParameters()->fShower_Yoffset[0] ), 29 );
					}
					fAnaShowerDir->SetMarkerColor( 6 );
					fAnaShowerDir->SetMarkerSize( 2. );
					fAnaShowerDir->SetMarkerStyle( 3 );
					fAnaShowerDir->Draw();
				}
				if( iMethod < fData->getShowerParameters()->fShower_Xoff_DISP.size() )
				{
					unsigned int d_min = 0;
					unsigned int d_max = fData->getShowerParameters()->fShower_Xoff_DISP[iMethod].size();
					if( !fBoolAllinOne )
					{
						d_min = fData->getTelID();
						d_max = fData->getTelID() + 1;
					}
                                        int z_color = 1;
					for( unsigned int d = d_min; d < d_max; d++ )
					{
						if( fData->getShowerParameters()->fShower_Xoff_DISP[iMethod][d] > -90. )
						{
							fAnaShowerDir = new TMarker( convertX( fData->getShowerParameters()->fShower_Xoff_DISP[iMethod][d] ),
										     convertY( fData->getShowerParameters()->fShower_Yoff_DISP[iMethod][d] ), 29 );
							if( fBoolAllinOne )
							{
								fAnaShowerDir->SetMarkerColor( z_color );
							}
							else
							{
								if( fData->getShowerParameters()->fShower_Weight_DISP[iMethod][d] > 0.5 )
								{
									fAnaShowerDir->SetMarkerColor( 4 );
								}
								else
								{
									fAnaShowerDir->SetMarkerColor( 3 );
								}
							}
							fAnaShowerDir->SetMarkerSize( 2. );
							fAnaShowerDir->SetMarkerStyle( 3 );
							fAnaShowerDir->Draw();
                                                        z_color++;
                                                        if( z_color >= 10 )
                                                        {
                                                            z_color++;
                                                        }
						}
					}
				}
			}
			// draw MC shower direction
			if( fData->getReader()->isMC() )
			{
				if( fData->getDetectorGeo()->getGrIsuVersion() >= 412 )
				{
					fMCShowerDir  = new TMarker( convertX( fData->getShowerParameters()->MCTel_Xoff ),
												 convertY( -1.*fData->getShowerParameters()->MCTel_Yoff ), 29 );
				}
				else
				{
					fMCShowerDir  = new TMarker( convertX( fData->getShowerParameters()->MCTel_Xoff ),
												 convertY( fData->getShowerParameters()->MCTel_Yoff ), 29 );
				}
				fMCShowerDir->SetMarkerColor( 1 );
				fMCShowerDir->SetMarkerSize( 2. );
				fMCShowerDir->Draw();
			}
			// draw Model3D shower
			if( fData->getRunParameter()->fUseDisplayModel3D )
			{
				fModel3DShowerDir = new TMarker( convertX( fData->getModel3DParameters()->fXoffModel3D ), convertY( -1.*fData->getModel3DParameters()->fYoffModel3D ), 29 );
				fModel3DShowerDir->SetMarkerColor( 3 );
				fModel3DShowerDir->SetMarkerSize( 2. );
				fModel3DShowerDir->Draw();
			}
			//Draw FROGS reconstruction
			if( fData->getRunParameter()->ffrogsmode && fData->getFrogsParameters() )
			{
				fFrogsShowerDir = new TMarker( convertX( fData->getFrogsParameters()->frogsXS ), convertY( 1.*fData->getFrogsParameters()->frogsYS ), 29 );
				fFrogsShowerDir->SetMarkerColor( 7 );
				fFrogsShowerDir->SetMarkerSize( 2. );
				fFrogsShowerDir->Draw();
			}
			// camera center
			fCameraCentreDir = new TMarker( convertX( 0. ), convertY( 0. ), 5 );
			fCameraCentreDir->Draw();
			// draw 0.5 deg circle
			fCameraCentreEllipse = new TEllipse( convertX( 0. ), convertY( 0. ), convertX( 0.5, 0. ), convertY( 0.5, 0. ) );
			fCameraCentreEllipse->SetLineStyle( 3 );
			fCameraCentreEllipse->SetFillStyle( 0 );
			fCameraCentreEllipse->Draw();
		}
	}
	else
	{
		fAnaEllipse->SetR1( 0. );
		fAnaEllipse->SetR2( 0. );
		fAnaEllipse->Draw();
	}
}


/*!
     this function is called after a mouseclick in the camera canvas

     selected channel is highlighted by bold tube, clicking anywhere resets that

     \return unique id of channel clicked on (200000+channel number), if missed return -1

*/

int VCamera::getChannel( int px, int py )
{
	return getChannel( px, py, 0 );
}


int VCamera::getChannel( int px, int py, TObject* objSel )
{
	double x = fCanvas->PadtoX( fCanvas->AbsPixeltoX( px ) );
	double y = fCanvas->PadtoY( fCanvas->AbsPixeltoY( py ) );
	double dist;
	
	fCanvas->SetEditable( true );
	fCanvas->cd();
	if( fTubeSelected >= 0 && fTubeSelected < ( int )fgraphTubes.size() )
	{
		hideSelectedChannel();
	}
	
	for( unsigned int i = 0; i < fgraphTubes.size(); i++ )
	{
		dist = sqrt( ( x - fgraphTubes[i]->GetX1() ) * ( x - fgraphTubes[i]->GetX1() ) + ( y - fgraphTubes[i]->GetY1() ) * ( y - fgraphTubes[i]->GetY1() ) );
		if( dist < fgraphTubes[i]->GetR1() )
		{
			showSelectedChannel( i, false );
			fCanvas->Update();
			fCanvas->SetEditable( false );
			return fgraphTubes[i]->GetUniqueID();
		}
	}
	fCanvas->Update();
	fCanvas->SetEditable( false );
	fCanvas->Update();
	fCanvas->SetEditable( false );
	return -1;
}


/*!
    channel selected by mouse is plotted with stronger line width
*/
void VCamera::showSelectedChannel( int i_channel, bool i_delete )
{
        fTubeSelectedV.clear();
	fCanvas->SetEditable( true );
	fCanvas->cd();
	if( i_delete )
	{
		hideSelectedChannel();
	}
	if( i_channel >= 0 && i_channel < ( int )fgraphTubes.size() )
	{
		fgraphTubes[i_channel]->SetLineWidth( 2 );
		fgraphTubes[i_channel]->SetLineColor( 1 );
		fgraphTubes[i_channel]->Draw();
		fTubeSelected = fgraphTubes[i_channel]->GetUniqueID() - 200000;
                fTubeSelectedV.push_back( fTubeSelected );
                if( fData && fData->getDetectorGeometry() && fData->getDetectorGeometry()->setTelID( fTelescope ) )
                {
                    if( i_channel < (int)fData->getDetectorGeometry()->getNNeighbours().size() )
                    {
                        for( unsigned int n = 0; n < fData->getDetectorGeometry()->getNNeighbours()[i_channel]; n++ )
                        {
                                 int iN = fData->getDetectorGeometry()->getNeighbours()[i_channel][n];
                                 if( iN < (int)fgraphTubes.size() )
                                 {
                                     fgraphTubes[iN]->SetLineColor( 8 );
                                     fgraphTubes[iN]->Draw();
                                     fTubeSelectedV.push_back( iN );
                                 }
                         }
                     }
                }
		fCanvas->Update();
	 }
	fCanvas->SetEditable( false );
}


/*!
    channel deselected by mouse is plotted with thinner line width
*/
void VCamera::hideSelectedChannel()
{
        for( unsigned int i = 0; i < fTubeSelectedV.size(); i++ )
        {
             if( fTubeSelectedV[i] > 0 && fTubeSelectedV[i] < ( int )fgraphTubes.size() )
             {
		fgraphTubes[fTubeSelectedV[i]]->SetLineWidth( 1 );
		fgraphTubes[fTubeSelectedV[i]]->SetLineColor( 15 );
		fgraphTubes[fTubeSelectedV[i]]->Draw();
             }
	}
	fTubeSelected = -1;
        fTubeSelectedV.clear();
}


/*!
     rescale data valarray entries to maximum = 1.

     values < -999999. and > 999999. are ignored

     \param v_value input vector
     \param iOffset only positive values in data value, subtract first minimal value

     \return   vector with rescaled values
*/
valarray<double>& VCamera::rescaleSums( valarray<double>& v_value, bool iOffset )
{
	if( fDebug )
	{
		cout << "VCamera::rescaleSums " << v_value.size() << "\t" << fData->getReader()->getFullHitVec().size() << endl;
	}
	double imax = -999999;
	double imin = 999999.;
	unsigned int i_v_valueSize = v_value.size();
	if( fScaleMax == 0. )
	{
		imax =  -999999.;
		// maximum/minimum
		for( unsigned i = 0; i < i_v_valueSize; i++ )
		{
			if( ( !fData->getDead( i, fData->getHiLo()[i] ) || fData->getLLEst()[i] ) )
			{
				if( abs( v_value[i] ) > imax )
				{
					imax = abs( v_value[i] );
				}
				if( abs( v_value[i] ) < imin )
				{
					imin = abs( v_value[i] );
				}
			}
		}
	}
	else
	{
		imax = fScaleMax;
	}
	
	if( iOffset )
	{
		for( unsigned i = 0; i < i_v_valueSize; i++ )
		{
			v_value[i] -= imin;
		}
		imax -= imin;
	}
	// rescale to maximum = 1.
	if( imax != 0. )
	{
		for( unsigned i = 0; i < i_v_valueSize; i++ )
		{
			v_value[i] /= imax;
		}
	}
	return v_value;
}


void VCamera::setCanvas( TPad* iCanvas )
{
	fCanvas = iCanvas;
}


void VCamera::setMode( unsigned int i_mod )
{
	fcameraModus = i_mod;
}


void VCamera::setPrintChannels( int i_pri )
{
	fPrintChannel = i_pri;
}


double VCamera::getMax( valarray<double>& i_val )
{
	unsigned int iSize = i_val.size();
	if( iSize == 0 )
	{
		return 0.;
	}
	// assume that nothing is smaller than that
	double max = -1.e10;
	//   for( unsigned int i = 0; i < i_val.size(); i++ ) if( !fData->getDead()[i] ) { max = i_val[i]; break; }
	for( unsigned int i = 0; i < iSize; i++ ) if( !fData->getDead( i, fData->getHiLo()[i] ) && i_val[i] > max )
		{
			max = i_val[i];
		}
	return max;
}


double VCamera::getMin( valarray<double>& i_val )
{
	unsigned int iSize = i_val.size();
	if( iSize == 0 )
	{
		return 0.;
	}
	// assume that nothing is bigger than that
	double min = 1.e10;
	//   for( unsigned int i = 0; i < i_val.size(); i++ ) if( !fData->getDead()[i] ) { min = i_val[i]; break; }
	for( unsigned int i = 0; i < iSize; i++ ) if( !fData->getDead( i, fData->getHiLo()[i] ) && i_val[i] < min )
		{
			min = i_val[i];
		}
	return min;
}


void VCamera::getMinMax( valarray<double>& i_val, double& imin, double& imax, vector< unsigned int > iPixelDrawMask )
{
	unsigned int iSize = i_val.size();
	if( iSize == 0 )
	{
		imin = 0.;
		imax = 0.;
	}
	double min = 1.e10;
	double max = -1.e10;
	for( unsigned int i = 0; i < iSize; i++ )
	{
                if( iPixelDrawMask[i] == 1 )
                {
			if( i_val[i] < min )
			{
				min = i_val[i];
			}
			if( i_val[i] > max )
			{
				max = i_val[i];
			}
		}
	}
	imin = min;
	imax = max;
}


double VCamera::convertX( double i_x, double i_off )
{
	double iDist_edge = fdist_edgeX;
	if( fBoolAllinOne )
	{
		iDist_edge = fmax_dist_edge * 0.9;
	}
	if( iDist_edge == 0. )
	{
		return 0.;
	}
	return ( i_x / iDist_edge * fmaxPlot + i_off );
}


double VCamera::convertY( double i_y, double i_off )
{
	double iDist_edge = fdist_edgeY;
	if( fBoolAllinOne )
	{
		iDist_edge = fmax_dist_edge * 0.9;
	}
	if( iDist_edge == 0. )
	{
		return 0.;
	}
	return ( i_y / iDist_edge * fmaxPlot + i_off );
}

/*

   draw a marker for each star in the FOV

*/
void VCamera::drawStarsInFOV()
{
	// check if star catalogue is available
	if( !fData->getStarCatalogue() )
	{
		return;
	}
	
	// get pointing of telescope
	float iTel_dec = 0.;
	float iTel_ra  = 0.;
	if( fTelescope < fData->getPointing().size() )
	{
		iTel_dec = fData->getPointing()[fTelescope]->getTelDec() * TMath::RadToDeg();
		iTel_ra  = fData->getPointing()[fTelescope]->getTelRA() * TMath::RadToDeg();
	}
	else
	{
		iTel_dec = fData->getArrayPointing()->getTelDec() * TMath::RadToDeg();
		iTel_ra  = fData->getArrayPointing()->getTelRA() * TMath::RadToDeg();
	}
	double iScale = 1.;
	if( fTelescope < fData->getDetectorGeometry()->getCameraScaleFactor().size() )
	{
		iScale = fData->getDetectorGeometry()->getCameraScaleFactor()[fTelescope];
	}
	
	vector< VStar* > iStar = fData->getStarCatalogue()->getListOfStarsinFOV();
	
	// draw a marker for each star
	double x_rot = 0.;
	double y_rot = 0.;
	char hname[200];
	for( unsigned int i = 0; i < iStar.size(); i++ )
	{
		if( iStar[i] && iStar[i]->fBrightness_B < fData->getRunParameter()->fMinStarBrightness_B )
		{
			double y = -1. * ( iStar[i]->fDecCurrentEpoch - iTel_dec );
			double x = 0.;
			if( cos( iTel_dec * TMath::DegToRad() ) != 0. )
			{
				x = -1. * ( iStar[i]->fRACurrentEpoch - iTel_ra ) * cos( iTel_dec * TMath::DegToRad() );
			}
			fData->getArrayPointing()->derotateCoords( fData->getEventMJD(), fData->getEventTime(), x, y, x_rot, y_rot );
			
			TMarker* iM = new TMarker( convertX( -1.*x_rot * iScale ), convertY( y_rot * iScale ), 5 );
			iM->SetMarkerColor( 2 );
			iM->Draw();
			sprintf( hname, "BMAG %.1f", iStar[i]->fBrightness_B );
			TText* iT = new TText( convertX( -1.*x_rot ), convertY( y_rot ), hname );
			iT->SetTextAngle( 45. );
			iT->SetTextSize( 0.0175 );
			iT->SetTextColor( 2 );
			iT->Draw();
		}
	}
	
}


