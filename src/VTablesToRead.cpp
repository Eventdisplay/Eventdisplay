/* \class VTablesToRead
   \brief service class for lookup tables


*/

#include "VTablesToRead.h"

VTablesToRead::VTablesToRead( unsigned int n_tabletypes, int iNTel )
{
    fNTel = ( unsigned int )iNTel;
    fNLT_types = n_tabletypes;
    
    for( unsigned int t = 0; t < fNLT_types; t++ )
    {
        vector< TH2F* > i_hnull;
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            i_hnull.push_back( 0 );
        }
        hMedian[t] = i_hnull;
        hSigma[t]  = i_hnull;
        
        value_T[t]       = new double[fNTel];
        value_T_sigma[t] = new double[fNTel];
        
        value[t] = -99.;
        value_Chi2[t] = -99.;
        value_dE[t] = -99.;
    }
    
    
    reset();
    
}


void VTablesToRead::reset()
{
    for( unsigned int t = 0; t < fNLT_types; t++ )
    {
        value[t] = -99.;
        value_Chi2[t] = -99.;
        value_dE[t] = -99.;
        
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            value_T[t][i] = -99.;
            value_T_sigma[t][i] = -99.;
        }
    }
}
