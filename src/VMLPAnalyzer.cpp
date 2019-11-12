/*! \class VMLPAnalyzer
    \brief MLP based disp method


*/

#include "VMLPAnalyzer.h"

VMLPAnalyzer::VMLPAnalyzer( string iFile )
{
    bZombie = false;
    fFile = 0;
    fMLP = 0;
    
    fFile = new TFile( iFile.c_str() );
    if( fFile->IsZombie() )
    {
        bZombie = true;
        return;
    }
    fMLP = ( TMultiLayerPerceptron* )fFile->Get( "TMultiLayerPerceptron" );
    if( !fMLP )
    {
        cout << "VMLPAnalyzer::VMLPAnalyzer Error: cannot find MLP object in " << iFile << endl;
        bZombie = true;
        return;
    }
    cout << "\t opened MLP data for direction reconstruction from " << fFile->GetName() << endl;
}


void VMLPAnalyzer::terminate()
{
    if( fFile )
    {
        fFile->Close();
    }
}


float VMLPAnalyzer::evaluate( float& width, float& length, float& asym, float& size, float& dist )
{
    if( size < 0. )
    {
        return -99999.;
    }
    if( length < 0. )
    {
        return -99999.;
    }
    if( !fMLP )
    {
        return -99999.;
    }
    if( bZombie )
    {
        return -99999.;
    }
    
    fmlp_var[0] = width / length;
    fmlp_var[1] = asym;
    fmlp_var[2] = dist;
    fmlp_var[3] = log10( size );
    
    return fMLP->Evaluate( 0, fmlp_var );
}
