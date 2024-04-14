//! VDST   data summarizer class (stores pixel sums and times in a tree)
#ifndef VDST_H
#define VDST_H

#include "TFile.h"
#include "TTree.h"

#include <bitset>
#include <iostream>

#include "VGlobalRunParameter.h"
#include "VImageBaseAnalyzer.h"
#include "VImageCleaning.h"
#include "VDSTTree.h"

///////////////////////////////////////////////////////////////////////////////////
// MAXIMUM NUMBER OF TELESCOPES AND CHANNELS IS DEFINED IN EVNDISP_definition.h
///////////////////////////////////////////////////////////////////////////////////

using namespace std;

class VDST : public VImageBaseAnalyzer, public VDSTTree
{
    private:
        bool fBLaser;
        VImageCleaning* fVImageCleaning;
        TFile* fDSTfile;

        bool fDSTini;

        bool writeCalibrationData();

    public:
        VDST( bool iMode, bool iMC );
        ~VDST();
        void fill();                              //!< fill dst tree
        TFile* getDSTFile()
        {
            return fDSTfile;
        }
        void initialize();
        void terminate();                         //!< write dst tree do disk
};
#endif
