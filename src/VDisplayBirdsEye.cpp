/*! \class VDisplayBirdsEye

    draw telescope and shower image from shower perspective

    all telescope positions, images, etc. are in shower coordinates

*/

#include <VDisplayBirdsEye.h>

VDisplayBirdsEye::VDisplayBirdsEye()
{
    fDebug = false;
    if( fDebug )
    {
        cout << "VDisplayBirdsEye::VDisplayBirdsEye" << endl;
    }
    fPlotPaper = false;
    // half width of field (in [m], draw it from +- fFieldX/Y)
    fFieldX = 300.;
    fFieldY = 300.;
    fFieldCentreX = 0.;
    fFieldCentreY = 0.;
}


/*!
    allow only a maximum of four telescopes (preliminary)

    (this will be overwritten at a later stage with the tru
     number of telescopes)
*/
void VDisplayBirdsEye::setNTel( unsigned int iNTel )
{
    if( iNTel < 5 )
    {
        fNTel = iNTel;
    }
    else
    {
        fNTel = 4;
    }
    
    if( fNTel == 2 )
    {
        fFieldX = 350;
        fFieldY = 350;
    }
}


void VDisplayBirdsEye::setPlotPaper()
{
    fPlotPaper = true;
    fFieldX = 200;
    fFieldY = 200;
}


/*!
    scale x ground coordinate in [m] to canvas coordinates
*/
double VDisplayBirdsEye::scaleX( double iS, double iCentre )
{
    return ( ( iS - iCentre ) / 2. / fFieldX * fXScale );
}


/*!
    scale y ground coordinate in [m] to canvas coordinates
*/
double VDisplayBirdsEye::scaleY( double iS, double iCentre )
{
    return ( ( iS - iCentre ) / 2. / fFieldY * fYScale );
}


double VDisplayBirdsEye::convertX( double iX )
{
    return scaleX( iX, fFieldCentreX ) + 0.5;
}


double VDisplayBirdsEye::convertY( double iY )
{
    return scaleX( iY, fFieldCentreY ) + 0.5;
}


void VDisplayBirdsEye::setParameters( VImageParameter* iPar )
{
    if( fParameters.size() == fData->getTeltoAna().size() )
    {
        fParameters.clear();
    }
    fParameters.push_back( iPar );
}


void VDisplayBirdsEye::draw( TPad* iPad )
{
    if( fDebug )
    {
        cout << "VDisplayBirdsEye::draw( TPad *iPad )" << endl;
    }
    iPad->SetEditable( true );
    iPad->cd();
    gPad->Clear();
    // no scaling, canvas is a square (otherwise distortion)
    fXScale = 1.;
    fYScale = 1.;
    // calculate telescope positions
    setGeometry();
    // draw the telescopes
    drawTelescopes_with_sizeAxis();
    // draw event info and reconstructed shower parameters
    drawEventText();
    iPad->Update();
    //   iPad->SetEditable( false );
}


void VDisplayBirdsEye::drawEventText()
{
    char iText [2000];
    
    if( fPlotPaper )
    {
        return;
    }
    
    unsigned int iM = fData->getRunParameter()->fPlotAllInOneMethod;
    
    fTextRec.clear();
    // run and event number:
    sprintf( iText, "Run/Event: %d/%d", fData->getRunNumber(), fData->getEventNumber() );
    fTextRec.push_back( new TText( 0.02, 0.96, iText ) );
    // MC text
    if( fData->getReader()->isMC() )
    {
        if( fData->getShowerParameters()->MCFirstInteractionHeight > 0 )
        {
            sprintf( iText, "Corsika production: run %d, shower %d", fData->getShowerParameters()->MCCorsikaRunID, fData->getShowerParameters()->MCCorsikaShowerID );
            fTextRec.push_back( new TText( 0.02, 0.18, iText ) );
            sprintf( iText, "MC: First interaction at %.3f km a.s.l. = %.3f g/cm^2",
                     fData->getShowerParameters()->MCFirstInteractionHeight / 1000.,
                     fData->getShowerParameters()->MCFirstInteractionDepth );
            fTextRec.push_back( new TText( 0.02, 0.15, iText ) );
        }
        sprintf( iText, "MC: E=%.2f TeV, Ze=%.0f deg, Az=%.0f deg",
                 fData->getShowerParameters()->MCenergy, fData->getShowerParameters()->MCze, fData->getShowerParameters()->MCaz );
        fTextRec.push_back( new TText( 0.02, 0.12, iText ) );
        sprintf( iText, "MC: Xoff=%.2f deg, Yoff=%.2f deg, Xcore=%.0f m, Ycore=%.0f m",
                 fData->getShowerParameters()->MCTel_Xoff, fData->getShowerParameters()->MCTel_Yoff,
                 fData->getShowerParameters()->MCxcore, fData->getShowerParameters()->MCycore );
    }
    // shower reconstruction text
    stringstream i_stext;
    i_stext << fData->getShowerParameters()->fShowerNumImages[iM];
    i_stext << " tel in reco (ID" << iM << "):";
    for( unsigned int i = 0; i < fData->getNTel(); i++ )
    {
        if( fData->getShowerParameters()->fTelIDImageSelected_list[iM][i] )
        {
	    i_stext << " " << ( int )( i + 1 );
        }
    }
    fTextRec.push_back( new TText( 0.02, 0.09, i_stext.str().c_str() ) );
    // triggered events
    // (any trigger condition)
    i_stext.str("");
    i_stext << fData->getShowerParameters()->fNTrig << " tel triggered: ";
    for( unsigned int i = 0; i < fData->getShowerParameters()->fNTrig; i++ )
    {
        i_stext << " " << fData->getShowerParameters()->fTrig_list[i] + 1;
    }
    fTextRec.push_back( new TText( 0.02, 0.06, i_stext.str().c_str() ) );
    
    sprintf( iText, "Ze, Az: [%.1f, %.1f] deg, X/Yoff: [%.2f, %.2f] deg, X/Ycore: [%.0f, %.0f ] m",
             fData->getShowerParameters()->fShowerZe[iM], fData->getShowerParameters()->fShowerAz[iM],
             fData->getShowerParameters()->fShower_Xoffset[iM], fData->getShowerParameters()->fShower_Yoffset[iM],
             fData->getShowerParameters()->fShowerXcore[iM], fData->getShowerParameters()->fShowerYcore[iM] );
    fTextRec.push_back( new TText( 0.02, 0.03, iText ) );
    
    for( unsigned int i = 0; i < fTextRec.size(); i++ )
    {
        fTextRec[i]->SetTextSize( 0.03 );
        fTextRec[i]->SetTextFont( 42 );
        fTextRec[i]->Draw();
    }
}


void VDisplayBirdsEye::drawTelescopes_with_sizeAxis()
{
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        // telescopes (circles)
        fElTel[i]->SetX1( convertX( fTelPosX[i] ) );
        fElTel[i]->SetY1( convertY( fTelPosY[i] ) );
        fElTel[i]->SetR1( scaleX( fTelRad[i] ) );
        fElTel[i]->SetR2( scaleY( fTelRad[i] ) );
        
        fElTelImage[i]->SetX1( convertX( fTelPosX[i] ) );
        fElTelImage[i]->SetY1( convertY( fTelPosY[i] ) );
        fElTelImage[i]->SetR1( scaleX( fTelRad[i] ) );
        fElTelImage[i]->SetR2( scaleY( fTelRad[i] ) );
    }
    
    // get min/max size or nphotons
    double iMin = 1.e20;
    double iMax = 0.;
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        fData->setTelID( fData->getTeltoAna()[i] );
        float iS = fData->getImageParameters()->size;
        if( iS > 0 && log10( iS ) < iMin )
        {
            iMin = log10( iS );
        }
        if( iS > 0 && log10( iS ) > iMax )
        {
            iMax = log10( iS );
        }
    }
    setColorAxisDataVector_minmax( iMin * 0.99, iMax * 1.01 );
    
    // fill histogram
    unsigned int i_tel_with_images = 0;
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        fData->setTelID( fData->getTeltoAna()[i] );
        if( fData->getImageParameters()->size > 0. )
        {
            i_tel_with_images++;
            // plot telescopes
            fElTelImage[i]->SetFillStyle( 1001 );
            fElTelImage[i]->SetR1( scaleX( fTelRad[i] + 3.*fTelRad[i] * ( log10( fData->getImageParameters()->size ) - iMin ) / ( iMax - iMin ) ) );
            fElTelImage[i]->SetR2( scaleX( fTelRad[i] + 3.*fTelRad[i] * ( log10( fData->getImageParameters()->size ) - iMin ) / ( iMax - iMin ) ) );
            fElTelImage[i]->SetFillColor( getColorAxisColor( log10( fData->getImageParameters()->size ) ) );
            fElTelImage[i]->SetLineColor( getColorAxisColor( log10( fData->getImageParameters()->size ) ) );
            fElTelImage[i]->Draw();
            
            fElTel[i]->SetFillStyle( 4001 );
            fElTel[i]->SetFillColor( getColorAxisColor( log10( fData->getImageParameters()->size ) ) );
            fElTel[i]->SetLineColor( 1 );
            fElTel[i]->Draw();
        }
        else
        {
            fElTel[i]->SetFillStyle( 4001 );
            fElTel[i]->SetLineColor( 1 );
            fElTel[i]->Draw();
        }
        if( fData->getReader()->hasLocalTrigger( fData->getTeltoAna()[i] ) )
        {
            TMarker* iLT = new TMarker( fElTel[i]->GetX1(), fElTel[i]->GetY1(), 2 );
            if( fData->getImageParameters()->size > 0. )
            {
                iLT->SetMarkerColor( 2 );
            }
            iLT->Draw();
        }
        // telescope identifier (T1,T2,etc.); draw only if less than 10 telescopes
        if( fData->getTeltoAna().size() < 10 )
        {
            fTextTel[i]->Draw();
        }
        else
        {
            if( fData->getReader()->hasLocalTrigger( fData->getTeltoAna()[i] ) )
            {
                fTextTel[i]->Draw();
            }
        }
    }
    // plot shower axis
    if( i_tel_with_images )
    {
        TGaxis* g = getColorAxisAxis( 0.92, 0.94, 0.70, 0.98, "log_{10} size", 10 );
        if( g && iMax > 0. )
        {
            g->Draw();
        }
    }
    
    drawImageLines_and_Corepositions();
    
}


/*!
      draw most of the graphical stuff
      - telescopes as circles
      - image lines (lines connect center of camera with centroids)
      - MC shower core
      - reconstructed shower cores
*/
void VDisplayBirdsEye::drawTelescopes()
{
    if( fDebug )
    {
        cout << "VDisplayBirdsEye::drawTelescopes()" << endl;
    }
    
    ///////////////////////////////
    // now loop over all telescopes
    ///////////////////////////////
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        // telescopes (circles)
        fElTel[i]->SetX1( convertX( fTelPosX[i] ) );
        fElTel[i]->SetY1( convertY( fTelPosY[i] ) );
        fElTel[i]->SetR1( scaleX( fTelRad[i] ) );
        fElTel[i]->SetR2( scaleY( fTelRad[i] ) );
        if( i + 1 < 10 )
        {
            fElTel[i]->SetFillColor( i + 1 );
        }
        else
        {
            fElTel[i]->SetFillColor( 1 );
        }
        // triggered telescopes are bold/red lined
        if( fData->getReader()->hasLocalTrigger( fData->getTeltoAna()[i] ) )
        {
            fElTel[i]->SetLineColor( 2 );
            fElTel[i]->SetLineWidth( 2 );
        }
        else
        {
            fElTel[i]->SetLineColor( 1 );
            fElTel[i]->SetLineWidth( 1 );
        }
        fElTel[i]->Draw();
        // telescope identifier (T1,T2,etc.); draw only if less than 10 telescopes
        if( fData->getTeltoAna().size() < 10 )
        {
            fTextTel[i]->Draw();
        }
    }
    drawImageLines_and_Corepositions();
}

void VDisplayBirdsEye::drawImageLines_and_Corepositions()
{
    // reconstruction method to draw
    unsigned int iM = fData->getRunParameter()->fPlotAllInOneMethod;
    // image line coordinates
    float i_x1 = 0.;
    float i_y1 = 0.;
    float i_x2 = 0.;
    float i_y2 = 0.;
    float i_x1_SC = 0.;
    float i_y1_SC = 0.;
    float i_x2_SC = 0.;
    float i_y2_SC = 0.;
    float i_z = 0.;
    float i_xcos = 0.;
    float i_ycos = 0.;
    // telescope are transformed into shower coordinates using MC real shower direction
    // telescopes are transformed into telescope pointing coordinates
    if( fData )
    {
        double degrad = 45. / atan( 1. );
        double i_y = fData->getShowerParameters()->fShower_Yoffset[iM];
        double i_x = fData->getShowerParameters()->fShower_Xoffset[iM];
        if( fData->getTeltoAna().size() > 0 && fData->getPointing()[fData->getTeltoAna()[0]] )
        {
            i_xcos = sin( ( 90. - fData->getPointing()[fData->getTeltoAna()[0]]->getTelElevation() + i_x ) / degrad ) *
                     sin( ( fData->getPointing()[fData->getTeltoAna()[0]]->getTelAzimuth() + i_y - 180. ) / degrad );
            i_ycos = sin( ( 90. - fData->getPointing()[fData->getTeltoAna()[0]]->getTelElevation() + i_x ) / degrad ) *
                     cos( ( fData->getPointing()[fData->getTeltoAna()[0]]->getTelAzimuth() + i_y - 180. ) / degrad );
        }
        else
        {
            cout << "VDisplayBirdsEye::drawImageLines_and_Corepositions() error: no telescope coordinates found" << endl;
            return;
        }
    }
    
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        // image lines
        if( fData->getShowerParameters()->fTelIDImageSelected[iM].size() > 0
                && fData->getShowerParameters()->fTelIDImageSelected[iM][fData->getTeltoAna()[i]] )
        {
            // image line coordinates in ground coordinates
            // image lines connects the reconstructed shower position with the centroid position
            // (assuming telescopes are pointing into shower direction)
            double i_y = -1.*fData->getShowerParameters()->fShower_Yoffset[iM];
            double i_x = -1.*fData->getShowerParameters()->fShower_Xoffset[iM];
            if( fData->getDetectorGeo()->getGrIsuVersion() >= 412 )
            {
                i_y *= -1.;
            }
            
            if( i_y < 99998. && i_x < 99998. )
            {
                double i_cen_x = ( fParameters[i]->cen_x + i_x );
                double i_cen_y = ( fParameters[i]->cen_y + i_y );
                
                i_x1 = convertX( fTelPosX[i] + 2.*fFieldX * cos( atan2( fMCSign * i_cen_y, i_cen_x ) ) );
                i_y1 = convertY( fTelPosY[i] + 2.*fFieldY * sin( atan2( fMCSign * i_cen_y, i_cen_x ) ) );
                i_x2 = convertX( fTelPosX[i] - 2.*fFieldX * cos( atan2( fMCSign * i_cen_y, i_cen_x ) ) );
                i_y2 = convertY( fTelPosY[i] - 2.*fFieldY * sin( atan2( fMCSign * i_cen_y, i_cen_x ) ) );
                
                double xC = convertX( fData->getShowerParameters()->MCxcore_SC );
                double yC = convertY( fData->getShowerParameters()->MCycore_SC );
                if( fData->getShowerParameters()->fShowerXcore_SC[iM] > -99998. )
                {
                    xC = convertX( fData->getShowerParameters()->fShowerXcore_SC[iM] );
                    yC = convertY( fData->getShowerParameters()->fShowerYcore_SC[iM] );
                }
                
                if( ( i_x1 - xC ) * ( i_x1 - xC ) + ( i_y1 - yC ) * ( i_y1 - yC ) > ( i_x2 - xC ) * ( i_x2 - xC ) + ( i_y2 - yC ) * ( i_y2 - yC ) )
                {
                    i_x1 = convertX( fTelPosX[i] );
                    i_y1 = convertY( fTelPosY[i] );
                }
                else
                {
                    i_x2 = convertX( fTelPosX[i] );
                    i_y2 = convertY( fTelPosY[i] );
                }
                
                // set image line coordinates
                fLiImage[i]->SetX1( i_x1 );
                fLiImage[i]->SetX2( i_x2 );
                fLiImage[i]->SetY1( i_y1 );
                fLiImage[i]->SetY2( i_y2 );
                fLiImage[i]->Draw();
            }
            
        }
        else
        {
            // don't draw anything
        }
        
    }
    // draw MC core position as a black cross
    if( fData->getReader()->isMC() )
    {
        fMarkerMCCore = new TMarker( convertX( fData->getShowerParameters()->MCxcore_SC ), convertY( fData->getShowerParameters()->MCycore_SC ), 28 );
        fMarkerMCCore->SetMarkerColor( 1 );
        fMarkerMCCore->SetMarkerSize( 1.5 );
        fMarkerMCCore->Draw();
    }
    // draw reconstructed shower cores with different colors (method 0 = color 2, method 1 = color 3, method 2 = color 4, ...)
    // draw the cores only for succesfull reconstruction (Chi2>=0.)
    fMarkerCore.clear();
    fMarkerCore.push_back( new TMarker( convertX( fData->getShowerParameters()->fShowerXcore_SC[iM] ), convertY( fData->getShowerParameters()->fShowerYcore_SC[iM] ), 3 ) );
    fMarkerCore.back()->SetMarkerColor( 1 );
    fMarkerCore.back()->SetMarkerSize( 3 );
    if( fData->getShowerParameters()->fShower_Chi2[iM] >= 0. )
    {
        fMarkerCore.back()->Draw();
    }
    //Draw FROGS reconstruction
    if( fData->getRunParameter()->ffrogsmode  && fData->getFrogsParameters() )
    {
        fMarkerFrogsCore = new TMarker( convertX( fData->getFrogsParameters()->frogsXP ), convertY( fData->getFrogsParameters()->frogsYP ), 29 );
        fMarkerFrogsCore->SetMarkerColor( 7 );
        fMarkerFrogsCore->SetMarkerSize( 2. );
        fMarkerFrogsCore->Draw();
    }
    // draw coordinate system
    // two arrow somewhere in the lower left corner of the canvas
    // the red arrow indicates the x-axis, the black the y-axis
    // (GM) not sure if the coordinate systems are correct
    float i_startx = 55.;
    i_startx = 0.;
    float i_stopx = i_startx + 35.;
    float i_starty = 55.;
    i_starty = 0.;
    float i_stopy = i_starty + 35.;
    if( fData )
    {
        tel_impact( i_xcos, i_ycos, i_startx, i_starty, 0., &i_x1_SC, &i_y1_SC, &i_z, false );
        tel_impact( i_xcos, i_ycos, i_stopx, i_starty, 0., &i_x2_SC, &i_y2_SC, &i_z, false );
        fAxis_SC_X = new TArrow( convertX( i_x1_SC ), convertY( i_y1_SC ), convertX( i_x2_SC ), convertY( i_y2_SC ), 0.03, "|>" );
        fAxis_SC_X->SetFillColor( 2 );
        fAxis_SC_X->SetLineColor( 2 );
        //      fAxis_SC_X->Draw();
        tel_impact( i_xcos, i_ycos, i_startx, i_starty, 0., &i_x1_SC, &i_y1_SC, &i_z, false );
        tel_impact( i_xcos, i_ycos, i_startx, i_stopy, 0., &i_x2_SC, &i_y2_SC, &i_z, false );
        fAxis_SC_Y = new TArrow( convertX( i_x1_SC ), convertY( i_y1_SC ), convertX( i_x2_SC ), convertY( i_y2_SC ), 0.03, "|>" );
        fAxis_SC_Y->SetFillColor( 1 );
        //      fAxis_SC_Y->Draw();
    }
}


void VDisplayBirdsEye::setData( VEvndispData* iData )
{
    if( fDebug )
    {
        cout << " VDisplayBirdsEye::setData( VEvndispData *iData )" << endl;
    }
    fData = iData;
    if( fData->getReader()->isGrisuMC() && fData->getDetectorGeo()->getGrIsuVersion() < 412 )
    {
        fMCSign =  -1.;
    }
    else
    {
        fMCSign = -1.;
    }
}


/*!
     calculate telesope positions and image orientation in schower coordinates
*/
void VDisplayBirdsEye::setGeometry()
{
    if( fDebug )
    {
        cout << "VDisplayBirdsEye::setGeometry()" << endl;
    }
    // reset all data vectors
    fNTel = fData->getTeltoAna().size();
    fTelPosX.clear();
    fTelPosY.clear();
    fTelRad.clear();
    fElTel.clear();
    fElTelImage.clear();
    fLiImage.clear();
    fTextTel.clear();
    
    double i_meanX = 0.;
    double i_meanY = 0.;
    double i_meanN = 0.;
    // calculate centre of array and  size of field to display
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        // get telescopes in shower coordinates
        fData->setTelID( fData->getTeltoAna()[i] );
        fTelPosX.push_back( fData->getImageParameters( fData->getRunParameter()->fImageLL )->Tel_x_SC );
        fTelPosY.push_back( fData->getImageParameters( fData->getRunParameter()->fImageLL )->Tel_y_SC );
        i_meanX += fTelPosX.back();
        i_meanY += fTelPosY.back();
        i_meanN++;
    }
    if( i_meanN > 0. )
    {
        fFieldCentreX = i_meanX / i_meanN;
        fFieldCentreY = i_meanY / i_meanN;
    }
    fFieldX = 0.;
    fFieldY = 0.;
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        if( fTelPosX[i] > fFieldX - fFieldCentreX )
        {
            fFieldX = fTelPosX[i] - fFieldCentreX;
        }
        if( fTelPosY[i] > fFieldY - fFieldCentreY )
        {
            fFieldY = fTelPosY[i] - fFieldCentreY;
        }
    }
    fFieldX = max( fFieldX, fFieldY ) + 350.;
    fFieldY = fFieldX;
    
    char i_text[200];
    for( unsigned int i = 0; i < fData->getTeltoAna().size(); i++ )
    {
        fData->setTelID( fData->getTeltoAna()[i] );
        // set telescope radius
        fTelRad.push_back( fData->getDetectorGeo()->getTelRadius()[fData->getTeltoAna()[i]] );
        // define telescopes at circles in canvas
        fElTel.push_back( new TEllipse( convertX( fTelPosX[i] ), convertY( fTelPosY[i] ), scaleX( fTelRad[i] ), scaleY( fTelRad[i] ) ) );
        // define shower image
        fElTelImage.push_back( new TEllipse( convertX( fTelPosX[i] ), convertY( fTelPosY[i] ), scaleX( fTelRad[i] ), scaleY( fTelRad[i] ) ) );
        // get image parameters
        fParameters.push_back( fData->getImageParameters( fData->getRunParameter()->fImageLL ) );
        // set telescope label
        if( fData->getTeltoAna().size() < 10 )
        {
            sprintf( i_text, "T%d", fData->getTeltoAna()[i] + 1 );
        }
        else
        {
            sprintf( i_text, "T%d", fData->getTeltoAna()[i] + 1 );
        }
        fTextTel.push_back( new TText( convertX( fTelPosX[i] ) + 0.02, convertY( fTelPosY[i] ), i_text ) );
        fTextTel.back()->SetTextFont( 42 );
        if( fData->getTeltoAna().size() < 10 )
        {
            fTextTel.back()->SetTextSize( 0.5 * fTextTel.back()->GetTextSize() );
        }
        else
        {
            fTextTel.back()->SetTextSize( 0.4 * fTextTel.back()->GetTextSize() );
            fTextTel.back()->SetTextAngle( 45. );
        }
        fTextTel.back()->SetTextColor( 13 );
        // set line representing the main image shower axis
        fLiImage.push_back( new TLine( 0., 0., 1., 1. ) );
        fLiImage.back()->SetLineStyle( 3 );
        fLiImage.back()->SetLineWidth( 1 );
        fLiImage.back()->SetLineColor( 50 );
    }
}
