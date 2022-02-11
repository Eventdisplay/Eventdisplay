/* \class VDL2Tree 
 * \brief VDL2Tree data class
 *
 */

#include "VDL2Tree.h"

VDL2Tree::VDL2Tree( string dl2_tree_name, 
                    string dl2_tree_title,
                    bool isMCtree )
{
    reset_tree_variables();

    fDL2EventTree = new TTree( dl2_tree_name.c_str(), dl2_tree_title.c_str() );
    fDL2EventTree->Branch( "runNumber", &fDL2_runNumber, "runNumber/i" );
    // save disk space: don't write eventNumber
    // fDL2EventTree->Branch( "eventNumber", &fDL2_eventNumber, "eventNumber/i" );
    if( isMCtree )
    {
        fDL2EventTree->Branch( "MCaz", &fDL2_MCaz, "MCaz/F" );
        fDL2EventTree->Branch( "MCel", &fDL2_MCel, "MCel/F" );
        fDL2EventTree->Branch( "MCe0", &fDL2_MCe0, "MCe0/F" );
    }
    // save disk space: don't write Xoff, Yoff
    // fDL2EventTree->Branch( "MCxoff", &fDL2_MCxoff, "MCxoff/F" );
    // fDL2EventTree->Branch( "MCyoff", &fDL2_MCyoff, "MCyoff/F" );
    fDL2EventTree->Branch( "ArrayPointing_Azimuth", &fDL2_ArrayPointing_Azimuth, "ArrayPointing_Azimuth/F" );
    fDL2EventTree->Branch( "ArrayPointing_Elevation", &fDL2_ArrayPointing_Elevation, "ArrayPointing_Elevation/F" );
    fDL2EventTree->Branch( "az", &fDL2_az, "az/F" );
    fDL2EventTree->Branch( "el", &fDL2_el, "el/F" );
    // save disk space: don't write Xoff, Yoff
    // fDL2EventTree->Branch( "xoff", &fDL2_xoff, "xoff/F" );
    // fDL2EventTree->Branch( "yoff", &fDL2_yoff, "yoff/F" );
    fDL2EventTree->Branch( "erec", &fDL2_erec, "erec/F" );
    fDL2EventTree->Branch( "nimages", &fDL2_nimages, "nimages/b" );
    fDL2EventTree->Branch( "CutClass", &fDL2_Cut_Class, "Class/b" );
    fDL2EventTree->Branch( "QC", &fDL2_Cut_QC, "QC/I" );
    fDL2EventTree->Branch( "MVA", &fDL2_Cut_MVA, "MVA/F" );
}

void VDL2Tree::fillEvent( CData *c, 
                          float mva,
                          int qc )
{
    if( fDL2EventTree && c )
    {
          fDL2_runNumber = (UInt_t)c->runNumber;
          fDL2_eventNumber = (UInt_t)c->eventNumber;
          fDL2_MCaz = c->MCaz;
          fDL2_MCel = 90. - c->MCze;
          fDL2_MCxoff = c->MCxoff;
          fDL2_MCyoff = c->MCyoff;
          fDL2_MCe0 = c->MCe0;
          fDL2_ArrayPointing_Azimuth = c->ArrayPointing_Azimuth;
          fDL2_ArrayPointing_Elevation = c->ArrayPointing_Elevation;
          fDL2_az = c->Az;
          fDL2_el = 90. - c->Ze;
          fDL2_xoff = c->Xoff;
          fDL2_yoff = c->Yoff;
          fDL2_erec = c->ErecS;
          fDL2_nimages = (UChar_t)c->NImages;

          fDL2_Cut_Class = 0;
          fDL2_Cut_QC = qc;
          fDL2_Cut_MVA = mva;
          fDL2EventTree->Fill();
    }
}

void VDL2Tree::reset_tree_variables()
{
    fDL2_runNumber = 0.;
    fDL2_eventNumber = 0;
    fDL2_MCaz = 0.;
    fDL2_MCel = 0.;
    fDL2_MCe0 = 0.;
    fDL2_MCxoff = 0.;
    fDL2_MCyoff = 0.;
    fDL2_ArrayPointing_Elevation = 0.;
    fDL2_ArrayPointing_Azimuth = 0.;
    fDL2_az = 0.;
    fDL2_el = 0.;
    fDL2_xoff = 0.;
    fDL2_yoff = 0.;
    fDL2_erec = 0.;
    fDL2_nimages = 0;
    fDL2_Cut_Class = 0;
    fDL2_Cut_QC = 0;
    fDL2_Cut_MVA = 0.;
}
