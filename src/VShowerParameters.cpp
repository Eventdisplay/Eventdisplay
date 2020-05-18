/*! \class VShowerParameters
    \brief  VShowerParameters storage class for shower data

    variables called XXX_SC are in shower coordinates

    some variables are the same for each event (data format, source file, etc.),
    not really necessary to store them in this tree (but convenient)

*/

#include "VShowerParameters.h"

VShowerParameters::VShowerParameters( int iNTel, unsigned int iShortTree, unsigned int iNMethod )
{
    fDebug = false;
    if( fDebug )
    {
        cout << "VShowerParameters::VShowerParameters( int iNTel ) " << iNTel << endl;
    }
    fMC = false;
    fTreeSC = 0;
    fNTel = iNTel;
    fNMethods = iNMethod;
    fShortTree = iShortTree;
    reset();
    
    vector< bool > it( fNTel, 0 );
    vector< float > i_f( fNTel, -99. );
    for( unsigned int i = 0; i < fNMethods; i++ )
    {
        fTelIDImageSelected.push_back( it );
        fTelIDImageSelected_bitcode[i] = 0;
        fShower_Xoff_DISP.push_back( i_f );
        fShower_Yoff_DISP.push_back( i_f );
        fShower_Weight_DISP.push_back( i_f );
    }
}


/*!
   initialise this data tree
   \param iName tree name
   \param iTitle tree title
   \param iMC  source data is Monte Carlo, include MC parameter in tree
*/
void VShowerParameters::initTree( string iName, string iTitle, bool iMC )
{
    char i_des[600];
    fTreeSC = new TTree( iName.c_str(), iTitle.c_str() );
    fTreeSC->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fTreeSC->SetAutoSave( 1024 * 1024 );          // autosave when 100 Mbytes written
    
    fTreeSC->Branch( "runNumber", &runNumber, "runNumber/I" );
    fTreeSC->Branch( "eventNumber",  &eventNumber,  "eventNumber/I" );
    fTreeSC->Branch( "eventStatus", &eventStatus, "eventStatus/i" );
    fTreeSC->Branch( "MJD",  &MJD,  "MJD/I" );
    fTreeSC->Branch( "Time",  &time,  "time/D" );
    fTreeSC->Branch( "NTel", &fNTelescopes, "NTel/i" );
    if( fShortTree < 1 )
    {
        fTreeSC->Branch( "TelElevation", fTelElevation, "TelElevation[NTel]/F" );
        fTreeSC->Branch( "TelAzimuth", fTelAzimuth, "TelAzimuth[NTel]/F" );
    }
    
    if( !iMC && fShortTree < 1 )
    {
        fTreeSC->Branch( "TelDec", fTelDec, "TelDec[NTel]/F" );
        fTreeSC->Branch( "TelRA", fTelRA, "TelRA[NTel]/F" );
    }
    if( fShortTree < 1 && !iMC )
    {
        fTreeSC->Branch( "TelElevationVBF", fTelElevationVBF, "TelElevationVBF[NTel]/F" );
        fTreeSC->Branch( "TelAzimuthVBF", fTelAzimuthVBF, "TelAzimuthVBF[NTel]/F" );
        fTreeSC->Branch( "TelPointingMismatch", fTelPointingMismatch, "TelPointingMismatch[NTel]/F" );
        fTreeSC->Branch( "PointingErrorX", fTelPointingErrorX, "PointingErrorX[NTel]/F" );
        fTreeSC->Branch( "PointingErrorY", fTelPointingErrorY, "PointingErrorY[NTel]/F" );
        
        sprintf( i_des, "Tel_x_SC[%d]/F", fNTel );
        fTreeSC->Branch( "Tel_x_SC", fTel_x_SC, i_des );
        sprintf( i_des, "Tel_y_SC[%d]/F", fNTel );
        fTreeSC->Branch( "Tel_y_SC", fTel_y_SC, i_des );
        sprintf( i_des, "Tel_z_SC[%d]/F", fNTel );
        fTreeSC->Branch( "Tel_z_SC", fTel_z_SC, i_des );
    }
    fTreeSC->Branch( "ArrayPointingElevation", &fArrayPointing_Elevation, "ArrayPointingElevation/F" );
    fTreeSC->Branch( "ArrayPointingAzimuth", &fArrayPointing_Azimuth, "ArrayPointingAzimuth/F" );
    fTreeSC->Branch( "ArrayPointing_deRotationAngle_deg", &fArrayPointing_deRotationAngle_deg, "ArrayPointing_deRotationAngle_deg/F" );
    
    fTreeSC->Branch( "WobbleN", &fWobbleNorth, "WobbleN/F" );
    fTreeSC->Branch( "WobbleE", &fWobbleEast, "WobbleE/F" );
    
    fTreeSC->Branch( "NTrig", &fNTrig, "NTrig/i" );
    fTreeSC->Branch( "LTrig", &fLTrig, "LTrig/l" );
    fTreeSC->Branch( "Trig_list", fTrig_list, "Trig_list[NTrig]/s" );
    fTreeSC->Branch( "Trig_type", fTrig_type, "Trig_type[NTrig]/s" );
    ///////////////////////////////////
    // reconstructed shower parameters
    
    // number of different methods
    sprintf( i_des, "NMethods/i" );
    fTreeSC->Branch( "NMethods", &fNMethods, i_des );
    if( fShortTree < 1 )
    {
        sprintf( i_des, "MethodID[NMethods]/s" );
        fTreeSC->Branch( "MethodID", fMethodID, i_des );
    }
    
    sprintf( i_des, "NImages[NMethods]/s" );
    fTreeSC->Branch( "NImages", fShowerNumImages, i_des );
    
    sprintf( i_des, "ImgSel[NMethods]/l" );
    fTreeSC->Branch( "ImgSel", fTelIDImageSelected_bitcode , i_des );
    sprintf( i_des, "ImgSel_list[NMethods][%d]/b", VDST_MAXTELESCOPES );
    fTreeSC->Branch( "ImgSel_list", fTelIDImageSelected_list, i_des );
    
    sprintf( i_des, "img2_ang[NMethods]/F" );
    fTreeSC->Branch( "img2_ang", fiangdiff, i_des );
    
    sprintf( i_des, "Ze[NMethods]/F" );
    fTreeSC->Branch( "Ze", fShowerZe, i_des );
    sprintf( i_des, "Az[NMethods]/F" );
    fTreeSC->Branch( "Az", fShowerAz, i_des );
    sprintf( i_des, "Xoff[NMethods]/F" );
    fTreeSC->Branch( "Xoff", fShower_Xoffset, i_des );
    sprintf( i_des, "Yoff[NMethods]/F" );
    fTreeSC->Branch( "Yoff", fShower_Yoffset, i_des );
    sprintf( i_des, "XoffDeRot[NMethods]/F" );
    fTreeSC->Branch( "XoffDeRot", fShower_XoffsetDeRot, i_des );
    sprintf( i_des, "YoffDeRot[NMethods]/F" );
    fTreeSC->Branch( "YoffDeRot", fShower_YoffsetDeRot, i_des );
    if( fShortTree < 1 )
    {
        sprintf( i_des, "stds[NMethods]/F" );
        fTreeSC->Branch( "stds", fShower_stdS, i_des );
    }
    if( !iMC && fShortTree < 1 )
    {
        sprintf( i_des, "dec[NMethods]/F" );
        fTreeSC->Branch( "dec", fDec, i_des );
        sprintf( i_des, "ra[NMethods]/F" );
        fTreeSC->Branch( "ra", fRA, i_des );
    }
    ///////////////////////////////////////////////////////////////
    // (debug only)
    // should be commented for general analysis
    if( fShortTree < 1 )
    {
        sprintf( i_des, "NPairXY/i" );
        fTreeSC->Branch( "NPairXY", &fShower_NPair, i_des );
        sprintf( i_des, "PairXs[NPairXY]/F" );
        fTreeSC->Branch( "PairXs", fShower_PairXS, i_des );
        sprintf( i_des, "PairYs[NPairXY]/F" );
        fTreeSC->Branch( "PairYs", fShower_PairYS, i_des );
        sprintf( i_des, "PairXd[NPairXY]/F" );
        fTreeSC->Branch( "PairXd", fShower_PairXD, i_des );
        sprintf( i_des, "PairYd[NPairXY]/F" );
        fTreeSC->Branch( "PairYd", fShower_PairYD, i_des );
        sprintf( i_des, "PairAngDiff[NPairXY]/F" );
        fTreeSC->Branch( "PairAngDiff", fShower_PairAngDiff, i_des );
        sprintf( i_des, "PairDispWeight[NPairXY]/F" );
        fTreeSC->Branch( "PairDispWeight", fShower_PairDispWeight, i_des );
    }
    // (end debug only)
    ///////////////////////////////////////////////////////////////
    sprintf( i_des, "Xcore[NMethods]/F" );
    fTreeSC->Branch( "Xcore", fShowerXcore, i_des );
    sprintf( i_des, "Ycore[NMethods]/F" );
    fTreeSC->Branch( "Ycore", fShowerYcore, i_des );
    if( fShortTree < 1 )
    {
        sprintf( i_des, "Xcore_SC[NMethods]/F" );
        fTreeSC->Branch( "Xcore_SC", fShowerXcore_SC, i_des );
        sprintf( i_des, "Ycore_SC[NMethods]/F" );
        fTreeSC->Branch( "Ycore_SC", fShowerYcore_SC, i_des );
        sprintf( i_des, "stdp[NMethods]/F" );
        fTreeSC->Branch( "stdp", fShower_stdP, i_des );
    }
    sprintf( i_des, "Chi2[NMethods]/F" );
    fTreeSC->Branch( "Chi2", fShower_Chi2, i_des );
    
    //DispDiff
    fTreeSC->Branch( "DispDiff", &fDispDiff, "DispDiff[NMethods]/F" );
    
    // MC parameters
    if( iMC )
    {
        fTreeSC->Branch( "MCprim", &MCprimary, "MCprimary/I" );
        fTreeSC->Branch( "MCe0", &MCenergy, "MCenergy/F" );
        fTreeSC->Branch( "MCxcore", &MCxcore, "MCxcore/F" );
        fTreeSC->Branch( "MCycore", &MCycore, "MCycore/F" );
        if( fShortTree < 1 )
        {
            fTreeSC->Branch( "MCxcos", &MCxcos, "MCxcos/F" );
            fTreeSC->Branch( "MCycos", &MCycos, "MCycos/F" );
        }
        fTreeSC->Branch( "MCze", &MCze, "MCze/F" );
        fTreeSC->Branch( "MCaz", &MCaz, "MCaz/F" );
        fTreeSC->Branch( "MCxoff", &MCTel_Xoff, "MCxoff/F" );
        fTreeSC->Branch( "MCyoff", &MCTel_Yoff, "MCyoff/F" );
        
        if( fShortTree < 1 )
        {
            fTreeSC->Branch( "MCxcore_SC", &MCxcore_SC, "MCxcore_SC/F" );
            fTreeSC->Branch( "MCycore_SC", &MCycore_SC, "MCycore_SC/F" );
            fTreeSC->Branch( "MCzcore_SC", &MCzcore_SC, "MCzcore_SC/F" );
        }
        fTreeSC->Branch( "MCCorsikaShowerID", &MCCorsikaShowerID, "MCCorsikaShowerID/I" );
        fTreeSC->Branch( "MCCorsikaRunID", &MCCorsikaRunID, "MCCorsikaRunID/I" );
        fTreeSC->Branch( "MCFirstInteractionHeight", &MCFirstInteractionHeight, "MCFirstInteractionHeight/F" );
        fTreeSC->Branch( "MCFirstInteractionDepth", &MCFirstInteractionDepth, "MCFirstInteractionDepth/F" );
        
    }
}


void VShowerParameters::reset()
{
    reset( VDST_MAXTELESCOPES );
}


void VShowerParameters::reset( unsigned int iNTel )
{
    runNumber = 0;
    eventNumber = 0;
    eventStatus = 0;
    time = 0.;
    MJD = 0;
    
    fNTelescopes = 1;
    fNTrig = 0;
    fLTrig = 0;
    
    fArrayPointing_Elevation = 0.;
    fArrayPointing_Azimuth = 0.;
    fArrayPointing_deRotationAngle_deg = 0.;
    fWobbleNorth     = 0.0;
    fWobbleEast      = 0.0;
    
    // reconstructed shower parameters
    fNumImages = 0;
    for( unsigned int i = 0; i < fNMethods; i++ )
    {
        fMethodID[i] = 0;
        fShowerNumImages[i] = 0;
        fTelIDImageSelected_bitcode[i] = 0;
        fShowerZe[i] = -90.;
        fShowerAz[i] = 0.;
        fDec[i] = -90.;
        fRA[i] = 0.;
        fShower_stdS[i] = 0.;
        fShower_Xoffset[i] = 0.;
        fShower_Yoffset[i] = 0.;
        fShower_XoffsetDeRot[i] = 0.;
        fShower_YoffsetDeRot[i] = 0.;
        fShowerXcore[i] = -9999.;
        fShowerYcore[i] = -9999.;
        fShowerXcore_SC[i] = -9999.;
        fShowerYcore_SC[i] = -9999.;
        fShower_stdP[i] = 0.;
        fShower_Chi2[i] = -1.;
        fiangdiff[i] = 0.0;
        fDispDiff[i] = -9999.;
        
        fShower_NPair = 0;
        fShower_PairXS[i] = 0.;
        fShower_PairYS[i] = 0.;
        fShower_PairXD[i] = 0.;
        fShower_PairYD[i] = 0.;
        fShower_PairAngDiff[i] = 0.;
        fShower_PairDispWeight[i] = 0.;
    }
    for( unsigned int j = 0; j < iNTel; j++ )
    {
        for( unsigned int i = 0; i < fNMethods; i++ )
        {
            fTelIDImageSelected_list[i][j] = 0;
        }
        fTrig_list[j] = 0;
        fTrig_type[j] = 0;
        fTelElevation[j] = 0.;
        fTelAzimuth[j] = 0.;
        fTelElevationVBF[j] = 0.;
        fTelAzimuthVBF[j] = 0.;
        fTelPointingMismatch[j] = 0.;
        fTelDec[j] = 0.;
        fTelRA[j] = 0.;
        fTelPointingErrorX[j] = 0.;
        fTelPointingErrorY[j] = 0.;
        fTel_x_SC[j] = 0.;
        fTel_y_SC[j] = 0.;
        fTel_z_SC[j] = 0.;
    }
    // MC parameters
    MCprimary = 0;
    MCenergy = 0.;
    MCxcore = 0.;
    MCycore = 0.;
    MCzcore = 0.;
    MCxcos = 0.;
    MCycos = 0.;
    MCxcore_SC = 0.;
    MCycore_SC = 0.;
    MCzcore_SC = 0.;
    MCaz = 0.;
    MCze = 0.;
    MCTel_Xoff = 0.;
    MCTel_Yoff = 0.;
    
    MCCorsikaShowerID = 0;
    MCFirstInteractionHeight = 0;
    MCFirstInteractionDepth = 0;
    MCCorsikaRunID = 0;
    
    resetDispVectors();
}


void VShowerParameters::printParameters()
{
    cout << "Shower parameters: " << endl;
    cout << runNumber << "\t" << eventNumber << " (status " << eventStatus << ")\t" <<  MJD;
    cout << time << "\t" << fNTelescopes << endl;
    cout << "Trigger: " << endl;
    cout << fNTrig << "\t" << fLTrig << endl;
    
    cout << "reconstructed parameters: " << endl;
    for( unsigned int i = 0; i < fNMethods; i++ )
    {
        cout << "Method " << i << endl;
        cout << "\t" << fShowerZe[i] << "\t" << fShowerAz[i] << "\t";
        cout << fShower_Xoffset[i] << "\t" << fShower_Yoffset[i] << "\t";
        cout << fShowerXcore[i] << "\t" << fShowerYcore[i] << endl;
        cout << fShowerXcore_SC[i] << "\t" << fShowerYcore_SC[i] << endl;
        cout << fiangdiff[i] << endl;
    }
    if( fMC )
    {
        cout << "\t" << MCenergy << "\t" << MCprimary;
        cout << "\t" << MCxcos << "\t" << MCycos;
        cout << "\t" << MCxcore << "\t" << "\t" << MCycore;
        cout << "\t" << MCxcore_SC << "\t" << MCycore_SC << endl;
    }
}


void VShowerParameters::resetDispVectors()
{
    for( unsigned int i = 0; i < fShower_Yoff_DISP.size(); i++ )
    {
        for( unsigned int j = 0; j < fShower_Yoff_DISP[i].size(); j++ )
        {
            fShower_Yoff_DISP[i][j] = -99.;
        }
    }
    for( unsigned int i = 0; i < fShower_Xoff_DISP.size(); i++ )
    {
        for( unsigned int j = 0; j < fShower_Xoff_DISP[i].size(); j++ )
        {
            fShower_Xoff_DISP[i][j] = -99.;
        }
    }
    for( unsigned int i = 0; i < fShower_Weight_DISP.size(); i++ )
    {
        for( unsigned int j = 0; j < fShower_Weight_DISP[i].size(); j++ )
        {
            fShower_Weight_DISP[i][j] = -99.;
        }
    }
}


void VShowerParameters::addDISPPoint( unsigned int iTelID, unsigned int iMethod, float x, float y, float w )
{
    if( iMethod < fShower_Xoff_DISP.size() && iTelID < fShower_Xoff_DISP[iMethod].size() )
    {
        fShower_Xoff_DISP[iMethod][iTelID] = x;
        fShower_Yoff_DISP[iMethod][iTelID] = y;
        fShower_Weight_DISP[iMethod][iTelID] = w;
    }
}
