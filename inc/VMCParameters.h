//! VMCParameters.h MC data storage class

#ifndef VMCParameters_H
#define VMCParameters_H

#include "TTree.h"

#include <iostream>

using namespace std;

class VMCParameters
{
    private:
    
        bool fDebug;
        
        TTree* fTree;
        
    public:
    
        unsigned int runNumber;
        unsigned int eventNumber;
        short int MCprimary;
        float MCenergy;                           //!< energy in [TeV]
        float MCxcore;
        float MCycore;
        float MCzcore;
        float MCxcos;
        float MCycos;
        float MCaz;
        float MCze;
        float MCTel_Xoff;                         //!< source offset in MC in deg (grisudet telescope coordinate system)
        float MCTel_Yoff;                         //!< source offset in MC in deg (grisudet telescope coordinate system)
        float MCFirstInteractionHeight;
        float MCFirstInteractionDepth;
        int MCCorsikaRunID;
        int MCCorsikaShowerID;
        short int ArrayTrigger;
        
        VMCParameters( bool iDebug = false );
        ~VMCParameters();
        
        void fill()
        {
            if( fTree )
            {
                fTree->Fill();
            }
        }
        void initTree();
        TTree* getTree()
        {
            return fTree;
        }
        void reset();
};
#endif
