//! VDisplayBirdsEye   draw telescope and shower image from shower perspective (eventdisplay display)

#ifndef VDisplayBirdsEye_H
#define VDisplayBirdsEye_H

#include <iostream>
#include <string>
#include <vector>

#include "TArrow.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "TText.h"

#include "VEvndispData.h"
#include "VGrIsuAnalyzer.h"
#include "VImageParameter.h"
#include "VPlotUtilities.h"

class VDisplayBirdsEye : public VGrIsuAnalyzer, public VPlotUtilities
{
    private:
        bool fDebug;
        double fMCSign;                           //!< +-1 to correct for sign errors in MC
        unsigned int fNTel;
        double fFieldX;
        double fFieldY;
        double fFieldCentreX;
        double fFieldCentreY;
        double fXScale;
        double fYScale;
        vector<double> fTelPosX;
        vector<double> fTelPosY;
        vector<double> fTelRad;
        
        vector< TEllipse* > fElTel;
        vector< TEllipse* > fElTelImage;
        vector< TMarker* > fMarkerCore;
        TMarker* fMarkerMCCore;
        TArrow* fAxis_SC_X;
        TArrow* fAxis_SC_Y;
        vector< TLine* > fLiImage;
        vector< TText* > fTextTel;
        vector< TText* > fTextRec;
        
        VEvndispData* fData;
        vector< VImageParameter* > fParameters;
        
        bool fPlotPaper;
        
        double convertX( double );
        double convertY( double );
        double scaleX( double a, double iCentre = 0. );
        double scaleY( double a, double iCentre = 0. );
        void drawEventText();
        void drawImagesinCamera();
        void drawImageLines_and_Corepositions();
        void drawTelescopes();
        void drawTelescopes_with_sizeAxis();
        void setGeometry();
        
    public:
    
        VDisplayBirdsEye();
        ~VDisplayBirdsEye() {}
        void draw( TPad* );
        void setData( VEvndispData* );
        void setNTel( unsigned int iNTel );
        void setParameters( VImageParameter* iParameters );
        void setPlotPaper();
};
#endif
