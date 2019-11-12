//! VModel3DData data class for Model3D

#ifndef VMODEL3DDATA_H
#define VMODEL3DDATA_H

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

#include <valarray>
#include <vector>

using namespace std;

class VModel3DData
{

    private:
    
    public:
    
        ///////////////// DATA ////////////////////
        bool fDebug3D; // set for debug mode
        
        bool fGoodPointing3D; // if pointing is same for all telescopes
        bool fGoodLnLTable3D; // if LnL table is read correctly
        bool fGoodEvent3D; // if event passes cuts for Model3D analysis
        
        unsigned int fNdim3D; // dimensional space (x,y,z)
        
        double fexNoise3D; // HARD-WIRED
        vector<double> fDCPE;  // for each telescope
        
        unsigned int fNTel3D; // total number of telescopes
        vector<unsigned int> fNpix3D; // total number of pixels for each telescope
        
        unsigned int fNtotPixel3D; // total number of selected pixels in full array
        
        //// telescope configuration parameters ////
        // for each telescope
        vector< vector< double > > flocTel3D; // location of telescope on ground
        vector<double> fMarea3D;              // mirror area (m^2)
        // for each telescope and pixel
        vector< vector< double > > fomegapix3D; // solid angle for each telescope and pixel
        
        //// pointing for each telescope (sub-array pointing not implemented) ////
        vector<double> fTze3D; // zenith
        vector<double> fTaz3D; // azimuth
        vector<double> fT3D;   // "T" vector of telescope pointing
        // for each telescope and pixel
        vector< vector< double > > fcosptheta3D;   // cos of angle between pixel center and telescope axis
        vector< vector< double > > fpX3D; // "p" vector (X in ground coordinates)
        vector< vector< double > > fpY3D; // "p" vector (Y in ground coordinates)
        vector< vector< double > > fpZ3D; // "p" vector (Z in ground coordinates)
        
        //// general coordinate frame conversion ////
        // for three dimensional space
        vector<double> fxg3D;   // X ground coordinate (unit base vectors)
        vector<double> fyg3D;   // Y ground coordinate (unit base vectors)
        vector<double> fzg3D;   // Z ground coordinate (unit base vectors)
        // for each telescope
        vector<double> fxsg3D;  // X sky unit base vector in ground coordinate frame
        vector<double> fysg3D;  // Y sky unit base vector in ground coordinate frame
        vector<double> fzsg3D;  // Z sky unit base vector in ground coordinate frame
        
        //// Model3D calculated signals ////
        // for each telescope and pixel
        vector< vector< double > > fMu3D; // Model3D calculated signal
        vector<double> fMuTel3D;          // Model3D total signal in each telescope
        
        ////// pixel signals (sums), cleaning, and pedvars ////
        // for each telescope and pixel
        vector< vector< bool > > fClean3D;         // pixels passing cleaning
        vector< vector< float > > fMeasuredSum3D; // measured signal
        vector< vector< float > > fPedvar3D;      // measured pedvars
        
        ////// model start parameters ////
        double fStartSel3D;    // model: elevation of shower direction (deg)
        double fStartSaz3D;    // model: azimuth of shower direction (deg)
        double fStartXcore3D;  // model: X shower core in ground coordinates
        double fStartYcore3D;  // model: Y shower core in ground coordinates
        double fStartSmax3D;   // model: height of shower maximum (along axis)
        double fStartsigmaL3D; // model: longitudinal (3D-length)
        double fStartsigmaT3D; // model: transverse (3D-width)
        double fStartNc3D;     // model: number of Cherenkov photons in shower
        
        //// model fit parameters ////
        double fSel3D;         // model: elevation of shower direction (deg)
        double fSaz3D;         // model: azimuth of shower direction (deg)
        double fXcore3D;       // model: X shower core in ground coordinates
        double fYcore3D;       // model: Y shower core in ground coordinates
        double fSmax3D;        // model: height of shower maximum (along axis)
        double fsigmaL3D;      // model: longitudinal (3D-length)
        double fsigmaT3D;      // model: transverse (3D-width)
        double fNc3D;          // model: number of Cherenkov photons in shower
        
        //// slant depth and reduced Width ////
        double fDepth3D;       // model: slant depth of shower maximum
        double fRWidth3D;      // model: reduced 3D-width
        double fErrRWidth3D;   // model: error in reduced 3D-width
        
        //// angular distance between Model3D and Hillas direction
        double fOmega3D;       // model and hillas: direction difference
        
        //// error in fit parameters ////
        double fErrorSel3D;    // model: elevation of shower direction (deg)
        double fErrorSaz3D;    // model: azimuth of shower direction (deg)
        double fErrorXcore3D;  // model: X shower core in ground coordinates
        double fErrorYcore3D;  // model: Y shower core in ground coordinates
        double fErrorSmax3D;   // model: height of shower maximum (along axis)
        double fErrorsigmaL3D; // model: longitudinal (3D-length)
        double fErrorsigmaT3D; // model: transverse (3D-width)
        double fErrorNc3D;     // model: number of Cherenkov photons in shower
        
        ////// derived shower direction from model ////
        double fXoffModel3D; // model shower direction
        double fYoffModel3D; // model shower direction
        double fXoffDeRot3D; // model shower direction (derotated)
        double fYoffDeRot3D; // model shower direction (derotated)
        
        ////// goodness-of-fit ////
        bool fConverged3D;              // fit converged
        unsigned int fNDF3D;            // number of degrees of freedom
        vector<unsigned int> fNDFTel3D; // NDF for each telescope
        double fStartGOF3D;             // starting goodness of fit
        double fGOF3D;                  // best-fit goodness of fit
        
        //// FUNCTIONS ////
        //// initialize ////
        void initModel3D( unsigned int iNTel3D, vector<unsigned int>& iNpix3D );
        void initEventModel3D(); // reset (called for each event)
        //// vector geometry ////
        void norm3D( double& ax, double& ay, double& az ); // normalize the vector
        void cross3D( double ax, double ay, double az, double bx, double by, double bz, double& cx, double& cy, double& cz ); // compute the cross product
        double dot3D( double ax, double ay, double az, double bx, double by, double bz ); // compute the dot product
        
        VModel3DData();
        ~VModel3DData();
        
};
#endif
