//! VPointingCorrectionsTreeReader reading class for pointing corrections

#ifndef VPointingCorrectionsTreeReader_H
#define VPointingCorrectionsTreeReader_H

#include "TChain.h"
#include "TTree.h"

#include <iostream>

using namespace std;

class VPointingCorrectionsTreeReader
{
    private:
      float fPointingErrorX;
      float fPointingErrorY;

      TChain* fTree;

    public:
       VPointingCorrectionsTreeReader( TChain *t = 0 );
       bool initialize();
       int getEntry( Long64_t );
       Long64_t getEntries();

       float getCorrected_cen_x( float cen_x );
       float getCorrected_cen_y( float cen_y );
       float getCorrected_phi( float cen_x, float cen_y, float d, float s, float sdevxy );

};

#endif
