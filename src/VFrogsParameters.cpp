/*! \class VFrogsParameters
    \brief  VFrogsParameters storage class for shower Frogs reconstruction

*/

#include "VFrogsParameters.h"

VFrogsParameters::VFrogsParameters()
{
    fDebug = false;
    fTreeFrogs = 0;
    reset();
    
}

void VFrogsParameters::initTree( string iName, string iTitle )
{

    printf( "FROGPUT In VFrogsParameters %s %s\n", iName.c_str(), iTitle.c_str() );
    
    fTreeFrogs = new TTree( iName.c_str(), iTitle.c_str() );
    fTreeFrogs->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fTreeFrogs->SetAutoSave( 1000 );
    
    //  Branches:
    fTreeFrogs->Branch( "frogsEventID", &frogsEventID, "frogsEventID/I" );
    fTreeFrogs->Branch( "frogsGSLConStat", &frogsGSLConStat, "frogsGSLConStat/I" );
    fTreeFrogs->Branch( "frogsNB_iter", &frogsNB_iter, "frogsNB_iter/I" );
    fTreeFrogs->Branch( "frogsNImages", &frogsNImages, "frogsNImages/I" );
    fTreeFrogs->Branch( "frogsSelectedImages", &frogsSelectedImages, "frogsSelectedImages/l" );
    fTreeFrogs->Branch( "frogsXS", &frogsXS, "frogsXS/F" );
    fTreeFrogs->Branch( "frogsXSerr", &frogsXSerr, "frogsXSerr/F" );
    fTreeFrogs->Branch( "frogsYS", &frogsYS, "frogsYS/F" );
    fTreeFrogs->Branch( "frogsYSerr", &frogsYSerr, "frogsYSerr/F" );
    fTreeFrogs->Branch( "frogsXP", &frogsXP, "frogsXP/F" );
    fTreeFrogs->Branch( "frogsXPerr", &frogsXPerr, "frogsXPerr/F" );
    fTreeFrogs->Branch( "frogsYP", &frogsYP, "frogsYP/F" );
    fTreeFrogs->Branch( "frogsYPerr", &frogsYPerr, "frogsYPerr/F" );
    fTreeFrogs->Branch( "frogsXPGC", &frogsXPGC, "frogsXPGC/F" );
    fTreeFrogs->Branch( "frogsYPGC", &frogsYPGC, "frogsYPGC/F" );
    fTreeFrogs->Branch( "frogsEnergy", &frogsEnergy, "frogsEnergy/F" );
    fTreeFrogs->Branch( "frogsEnergyerr", &frogsEnergyerr, "frogsEnergyerr/F" );
    fTreeFrogs->Branch( "frogsLambda", &frogsLambda, "frogsLambda/F" );
    fTreeFrogs->Branch( "frogsLambdaerr", &frogsLambdaerr, "frogsLambdaerr/F" );
    fTreeFrogs->Branch( "frogsGoodnessImg", &frogsGoodnessImg, "frogsGoodnessImg/F" );
    fTreeFrogs->Branch( "frogsNpixImg", &frogsNpixImg, "frogsNpixImg/I" );
    fTreeFrogs->Branch( "frogsGoodnessBkg", &frogsGoodnessBkg, "frogsGoodnessBkg/F" );
    fTreeFrogs->Branch( "frogsNpixBkg", &frogsNpixBkg, "frogsNpixBkg/I" );
    fTreeFrogs->Branch( "frogsTelGoodnessImg0", &frogsTelGoodnessImg0, "frogsTelGoodnessImg0/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessImg1", &frogsTelGoodnessImg1, "frogsTelGoodnessImg1/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessImg2", &frogsTelGoodnessImg2, "frogsTelGoodnessImg2/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessImg3", &frogsTelGoodnessImg3, "frogsTelGoodnessImg3/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessBkg0", &frogsTelGoodnessBkg0, "frogsTelGoodnessBkg0/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessBkg1", &frogsTelGoodnessBkg1, "frogsTelGoodnessBkg1/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessBkg2", &frogsTelGoodnessBkg2, "frogsTelGoodnessBkg2/F" );
    fTreeFrogs->Branch( "frogsTelGoodnessBkg3", &frogsTelGoodnessBkg3, "frogsTelGoodnessBkg3/F" );
    
    fTreeFrogs->Branch( "frogsXPStart", &frogsXPStart, "frogsXPStart/F" );
    fTreeFrogs->Branch( "frogsYPStart", &frogsYPStart, "frogsYPStart/F" );
    fTreeFrogs->Branch( "frogsXPED", &frogsXPED, "frogsXPED/F" );
    fTreeFrogs->Branch( "frogsYPED", &frogsYPED, "frogsYPED/F" );
    fTreeFrogs->Branch( "frogsXSStart", &frogsXSStart, "frogsXSStart/F" );
    fTreeFrogs->Branch( "frogsYSStart", &frogsYSStart, "frogsYSStart/F" );
    
    fTreeFrogs->Branch( "frogsZe", &frogsZe, "frogsZe/F" );
    fTreeFrogs->Branch( "frogsAz", &frogsAz, "frogsAz/F" );
    fTreeFrogs->Branch( "frogsXS_derot", &frogsXS_derot, "frogsXS_derot/F" );
    fTreeFrogs->Branch( "frogsYS_derot", &frogsYS_derot, "frogsYS_derot/F" );
    
    TString temp = TString::Format( "frogsR[%d]/D", fNTel );
    fTreeFrogs->Branch( "frogsR" , &frogsR, temp.Data() );
}



void VFrogsParameters::reset()
{

    //  0 all the values
    frogsEventID     = 0;
    frogsGSLConStat  = 0;
    frogsNB_iter     = 0;
    frogsNImages     = 0;
    frogsSelectedImages = 0;
    frogsXS          = 0.0;
    frogsXSerr       = 0.0;
    frogsYS          = 0.0;
    frogsYSerr       = 0.0;
    frogsXP          = 0.0;
    frogsXPerr       = 0.0;
    frogsYP          = 0.0;
    frogsYPerr       = 0.0;
    frogsEnergy      = 0.0;
    frogsEnergyerr   = 0.0;
    frogsLambda      = 0.0;
    frogsLambdaerr   = 0.0;
    frogsGoodnessImg = 0.0;
    frogsNpixImg     = 0;
    frogsGoodnessBkg  = 0.0;
    frogsNpixBkg      = 0;
    frogsTelGoodnessImg0 = 0.;
    frogsTelGoodnessImg1 = 0.;
    frogsTelGoodnessImg2 = 0.;
    frogsTelGoodnessImg3 = 0.;
    frogsTelGoodnessBkg0 = 0.;
    frogsTelGoodnessBkg1 = 0.;
    frogsTelGoodnessBkg2 = 0.;
    frogsTelGoodnessBkg3 = 0.;
    frogsXPStart = 0.;
    frogsYPStart = 0.;
    frogsXPED    = 0.;
    frogsYPED    = 0.;
    frogsXSStart = 0.;
    frogsYSStart = 0.;
    frogsXS_derot = 0.;
    frogsYS_derot = 0.;
    frogsZe = 0;
    frogsAz = 0;
    for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
    {
        frogsR[i] = 0;
    }
}

