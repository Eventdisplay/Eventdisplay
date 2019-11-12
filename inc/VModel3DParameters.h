//! VModel3DParameters storage class for Frogs data
#ifndef VMODEL3DPARAMETERS_H
#define VMODEL3DPARAMETERS_H

#include "VGlobalRunParameter.h"
#include "TTree.h"

#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using namespace std;

class VModel3DParameters
{
    private:
        TTree* fTreeModel3D;                         //!< output tree
        
    public:
        VModel3DParameters();
        ~VModel3DParameters();
        
        void fill()
        {
            if( fTreeModel3D )
            {
                fTreeModel3D->Fill();
            }
        }
        TTree* getTree()
        {
            return fTreeModel3D;
        }
        void initTree( string, string );
        void reset();        //!< reset all tree variable to standard values
        
        int eventNumber;
        // Model3D parameters, JG //
        float fSel3D;    // elevation of shower direction (deg)
        float fSaz3D;    // azimuth of shower direction (deg)
        float fXcore3D;  // shower core in ground coordinates
        float fYcore3D;  // shower core in ground coordinates
        float fSmax3D;   // height of shower maximum (along the shower axis)
        float fsigmaL3D; // longitudinal (3D-length)
        float fsigmaT3D; // transverse (3D-width)
        float fNc3D;     // total number of Cherenkov photons emitted by the shower
        float fXoffModel3D;  // model sky direction
        float fYoffModel3D;  // model sky direction
        float fXoffDeRot3D; // model shower direction (derotated)
        float fYoffDeRot3D; // model shower direction (derotated)
        
        float fGoodness3D;   // model goodness of fit
        float fDepth3D;      // model: slant depth of shower maximum
        float fRWidth3D;     // model: reduced 3D-width
        float fErrRWidth3D;  // model: error in reduced 3D-width
        float fOmega3D;      //model and hillas: direction difference
        bool fConverged3D;   // model: fit converged
        /// debug ///
        float fStartGoodness3D;   // model goodness of fit
        float fStartSel3D;    // Start elevation of shower direction (deg)
        float fStartSaz3D;    // Start azimuth of shower direction (deg)
        float fStartXcore3D;  // Start shower core in ground coordinates
        float fStartYcore3D;  // Start shower core in ground coordinates
        float fStartSmax3D;   // Start height of shower maximum (along the shower axis)
        float fStartsigmaL3D; // Start longitudinal (3D-length)
        float fStartsigmaT3D; // Start transverse (3D-width)
        float fStartNc3D;     // Start total number of Cherenkov photons emitted by the shower
        float fErrorSel3D;    // Start elevation of shower direction (deg)
        float fErrorSaz3D;    // Start azimuth of shower direction (deg)
        float fErrorXcore3D;  // Start shower core in ground coordinates
        float fErrorYcore3D;  // Start shower core in ground coordinates
        float fErrorSmax3D;   // Start height of shower maximum (along the shower axis)
        float fErrorsigmaL3D; // Start longitudinal (3D-length)
        float fErrorsigmaT3D; // Start transverse (3D-width)
        float fErrorNc3D;     // Start total number of Cherenkov photons emitted
        
};
#endif
