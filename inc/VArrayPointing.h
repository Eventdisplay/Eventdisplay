//! VArrayPointing relevant pointing information for array (e.g. reference pointing)

#ifndef VArrayPointing_H
#define VArrayPointing_H

#include "TMath.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>

#include "VSkyCoordinates.h"

using namespace std;

class VArrayPointing : public VSkyCoordinates
{
    private:
    
        TTree* fPointingTree;
        
        // smaller version of fPointingTree
        // contains pointing info at 1/second, instead of 1/event
        // pointing is from interpolated
        TTree* fPntReduced;
        
        void initializePointingTree();
        
    public:
    
        VArrayPointing( bool bInitTree = true );
        ~VArrayPointing() {}
        void fillPointingTree( bool bIsMC );
        void fillPntReduced();
        void terminate( bool iDebug_IO = false, bool bIsMC = false );
        
};

#endif
