//! DL2 Tree definition
#ifndef VDL2Tree_H
#define VDL2Tree_H

#include <string>
#include "TTree.h"
#include "CData.h"


class VDL2Tree
{
    private:

    TTree *fDL2EventTree;
    bool   fDL2WriteFullEventTree;
    UInt_t fDL2_runNumber;
    UInt_t fDL2_eventNumber;
    float fDL2_MCaz;
    float fDL2_MCel;
    float fDL2_MCe0;
    float fDL2_MCxoff;
    float fDL2_MCyoff;
    float fDL2_ArrayPointing_Elevation;
    float fDL2_ArrayPointing_Azimuth;
    float fDL2_az;
    float fDL2_el;
    float fDL2_xoff;
    float fDL2_yoff;
    float fDL2_erec;
    UChar_t fDL2_nimages;
    UChar_t fDL2_Cut_Class;
    float fDL2_Cut_MVA;
    Int_t fDL2_Cut_QC;

    void reset_tree_variables();

    public:

    VDL2Tree( string dl2_tree_name = "DL2EventTree",
              string dl2_tree_title = "DL2 tree",
              bool isMCtree = false );
   ~VDL2Tree();
    void fillEvent( CData *c, float mva, int qc );
    TTree *getDL2Tree() { return fDL2EventTree; }

};

#endif
