//! VPETreeReader write for pe data to root files

#ifndef VPETree_H
#define VPETree_H

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

class VPETree
{
    private:
    
        bool fReadFullTree;
        
        TTree* fTree;
        
    public:
    
        unsigned int fEventNumber;
        unsigned int fPrimaryType;
        float        fPrimaryEnergy;
        float        fXcore;
        float        fYcore;
        float        fXcos;
        float        fYcos;
        float        fXsource;
        float        fYsource;
        std::vector< unsigned int >* v_f_ID;
        std::vector< float >* v_f_time;
        TBranch* b_v_f_ID;
        TBranch* b_v_f_time;
        
        VPETree( TTree* t );
        bool initialize( bool bFull );
        int getEntries();
        int getEntry( int );
        TTree* getDataTree()
        {
            return fTree;
        }
        bool   isFull()
        {
            return fReadFullTree;
        }
};
#endif
