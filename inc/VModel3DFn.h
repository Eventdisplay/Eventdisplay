//! VModel3DFn

#ifndef VMODEL3DFN_H
#define VMODEL3DFN_H

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

#include <valarray>
#include <vector>

#include "TMath.h"

using namespace std;

class VModel3DFn
{
    private:
        //// model fit parameters ////
        double fSel;      // model: elevation of shower direction (deg)
        double fSaz;      // model: azimuth of shower direction (deg)
        double fXcore;    // model: X shower core in ground coordinates
        double fYcore;    // model: Y shower core in ground coordinates
        double fSmax;     // model: height of shower maximum along the shower axis
        double fsigmaL;   // model: longitudinal (3D-length)
        double fsigmaT;   // model: transverse (3D-width)
        double fNc;       // model: number of Cherenkov photons in shower
        
        unsigned int fNdim3D;
        unsigned int fNTel3D;
        
        //// geometry for model fit ////
        vector< vector< double > > flocTel3D; // location of telescopes
        
        vector<double> fs;   // "s" unit vector along the shower axis
        vector<double> fxBo; // "xB" vector of telescope optical center and shower barycenter B at (0,0,0) ground
        vector< vector< double > > fxB; // "xB" vector of telescope optical center and shower barycenter B for each telescope
        
        double freq( double x ) const;
        
        double dot( double ax, double ay, double az, double bx, double by, double bz ) const
        {
            return ( ax * bx ) + ( ay * by ) + ( az * bz );
        }
        
    public:
    
        VModel3DFn();
        ~VModel3DFn();
        
        void setModel3DFn( const vector< vector< double > >& locTel3D, const vector<double>& param );
        
        void calcModel3DFn( unsigned int iTel, unsigned int iPix, const double ipX, const double ipY, const double ipZ, double& val, vector<double>& grad ) const;
        
};
#endif
