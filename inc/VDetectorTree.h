//! VDetectorTree tree with basic telescope parameters (positions, ...)

#ifndef VDETECTORTREE
#define VDETECTORTREE

#include "TTree.h"

#include "VDetectorGeometry.h"

#include <iostream>

using namespace std;

class VDetectorTree
{
    private:
        bool fDebug;
        
        TTree* fTreeDet;
        
    public:
    
        unsigned int fNTel;
        float fTelxpos;
        float fTelypos;
        float fTelzpos;
        
        VDetectorTree();
        ~VDetectorTree();
        bool fillDetectorTree( VDetectorGeometry* iDet );
        bool readDetectorTree( VDetectorGeometry* iDet, TTree* iTree, bool iCTA );
        TTree* getTree()
        {
            return fTreeDet;
        }
};
#endif
