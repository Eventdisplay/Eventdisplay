#define Cmodel3Dpars_cxx
#include "Cmodel3Dpars.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Cmodel3Dpars::Loop()
{
    if( fChain == 0 )
    {
        return;
    }
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for( Long64_t jentry = 0; jentry < nentries; jentry++ )
    {
        Long64_t ientry = LoadTree( jentry );
        if( ientry < 0 )
        {
            break;
        }
        nb = fChain->GetEntry( jentry );
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
    }
}
