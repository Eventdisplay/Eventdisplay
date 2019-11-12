/*! \class VGammaHadronCutsStatistics

    keeps track of why events are cut away

    fills a tree with same number of events as data tree (use AddFriend...)

    Note: VGammaHadronCutsStatistics::newEvents() should be called for all events

*/

#include "VGammaHadronCutsStatistics.h"

VGammaHadronCutsStatistics::VGammaHadronCutsStatistics()
{
    fData = 0;
    
    fData = 0;
    
    reset();
}

void VGammaHadronCutsStatistics::initialize()
{
    fCutName.push_back( "Tot               " );
    fCutName.push_back( "MC_XYoff          " );
    fCutName.push_back( "XYoff             " );
    fCutName.push_back( "StereoQuality     " );
    fCutName.push_back( "  QC ArrayChi2    " );
    fCutName.push_back( "  QC NImages      " );
    fCutName.push_back( "  QC MSC_Quality  " );
    fCutName.push_back( "  QC Erec         " );
    fCutName.push_back( "  QC CorePos      " );
    fCutName.push_back( "  QC LTrig        " );
    fCutName.push_back( "SizeSecondMax     " );
    fCutName.push_back( "TelType           " );
    fCutName.push_back( "Direction         " );
    fCutName.push_back( "IsGamma           " );
    fCutName.push_back( "EnergyRec         " );
    fCutName.push_back( "Unkown cut (problem?) " );
    
    fData = new TTree( "GammaHadronCutsStats", "cut statistics for gamma/hadron cuts" );
    fData->Branch( "cut", &fCut_bitset_ulong, "cut/l" );
}

void VGammaHadronCutsStatistics::reset()
{
    fCutCounter.clear();
    
    fCut_bitset.reset();
    fCut_bitset_ulong = 0;
    
    // the vector will have exactly the name of EN_AnaCutsStats
    for( unsigned int i = 0; i < fCutName.size(); i++ )
    {
        fCutCounter.push_back( 0 );
    }
}

unsigned int VGammaHadronCutsStatistics::getCounterValue( unsigned int iCut )
{
    if( iCut < fCutCounter.size() )
    {
        return fCutCounter[iCut];
    }
    
    return 0;
}

void VGammaHadronCutsStatistics::fill()
{
    // fill data tree
    fCut_bitset_ulong = fCut_bitset.to_ulong();
    fData->Fill();
    
    // reset bit counter
    fCut_bitset.reset();
    fCut_bitset_ulong = 0;
}

/*
   make sure that last event is filled
*/
void VGammaHadronCutsStatistics::terminate()
{
    if( getCounterValue( eTot ) != fData->GetEntries() )
    {
        fCut_bitset_ulong = fCut_bitset.to_ulong();
        fData->Fill();
    }
}

void VGammaHadronCutsStatistics::updateCutCounter( unsigned int iCut )
{
    if( iCut < fCutCounter.size() )
    {
        fCutCounter[iCut]++;
        fCut_bitset.set( iCut, true );
    }
}

void VGammaHadronCutsStatistics::setCutCounter( unsigned int iCut, unsigned int iValue )
{
    if( iCut < fCutCounter.size() )
    {
        fCutCounter[iCut] = iValue;
    }
}

void VGammaHadronCutsStatistics::printCutStatistics()
{
    cout << endl;
    cout << "\t cut statistics: " << endl;
    for( unsigned int i = 0; i < fCutName.size(); i++ )
    {
        cout << "\t\t" << fCutName[i] << "\t";
        if( i < fCutCounter.size() )
        {
            cout << fCutCounter[i] << endl;
        }
        else
        {
            cout << "no cut defined" << endl;
        }
    }
    cout << endl;
}
