/*! \file VExclusionRegions
    \brief handling of exclusion regions due to e.g. bright stars

    includes data classes for bright star exclusion regions,
    exlusion regions, coordinate transformations for exclusion regions

*/

#include "VExclusionRegions.h"


//==================================================================================
// star exclusion regions settings
//==================================================================================

VBrightStarExclusionSettings::VBrightStarExclusionSettings()
{
    fStarMinBrightness      = -100.;
    fStarBand               = "B";
    fStarExlusionRadius_DEG = -1.;
}

VBrightStarExclusionSettings::~VBrightStarExclusionSettings() {}

/*
 *  operator used to sort bright star settings according to their magnitude
 */
bool VBrightStarExclusionSettings::operator<( const VBrightStarExclusionSettings& s1 ) const
{
    return fStarMinBrightness < s1.fStarMinBrightness;
}


//==================================================================================
// list of exclusion regions
//==================================================================================

VListOfExclusionRegions::VListOfExclusionRegions()
{
    fExcludeFromBackground_West  = -9999.;
    fExcludeFromBackground_North = -9999.;
    fExcludeFromBackground_CameraCentre_x = -9999.;
    fExcludeFromBackground_CameraCentre_y = -9999.;
    fExcludeFromBackground_Radius1 = 0.;
    fExcludeFromBackground_Radius2 = 0.;
    fExcludeFromBackground_RotAngle = 0.;
    fExcludeFromBackground_StarID = -1;
    fExcludeFromBackground_StarName = "";
    fExcludeFromBackground_DecJ2000 = -9999.;
    fExcludeFromBackground_RAJ2000  = -9999.;
    fExcludeFromBackground_B_Band   = true;
    fExcludeFromBackground_StarBrightness_V = -9999.;
    fExcludeFromBackground_StarBrightness_B = -9999.;
}

VListOfExclusionRegions::~VListOfExclusionRegions() {}

bool VListOfExclusionRegions::operator<( const VListOfExclusionRegions& s1 ) const
{
    return fExcludeFromBackground_StarBrightness_B < s1.fExcludeFromBackground_StarBrightness_B;
}

/*

    check if these coordinates are inside an exclusion region

    returns true if direction is inside the region

*/
bool VListOfExclusionRegions::isInsideExclusionRegion( double x_cam_deg, double y_cam_deg, double r_offset_deg )
{

    if( TMath::Power( ( ( x_cam_deg - fExcludeFromBackground_CameraCentre_x ) * TMath::Cos( fExcludeFromBackground_RotAngle * TMath::DegToRad() )
                        + ( y_cam_deg - fExcludeFromBackground_CameraCentre_y ) * TMath::Sin( fExcludeFromBackground_RotAngle * TMath::DegToRad() ) )
                      / ( fExcludeFromBackground_Radius1 + r_offset_deg ), 2. )
            + TMath::Power( ( ( x_cam_deg - fExcludeFromBackground_CameraCentre_x ) * TMath::Sin( fExcludeFromBackground_RotAngle * TMath::DegToRad() )
                              - ( y_cam_deg - fExcludeFromBackground_CameraCentre_y ) * TMath::Cos( fExcludeFromBackground_RotAngle * TMath::DegToRad() ) )
                            / ( fExcludeFromBackground_Radius2 + r_offset_deg ), 2. )
            < 1. )
    {
        return true;
    }
    
    return false;
}


//==================================================================================
// exlusion regions handler
//==================================================================================

VExclusionRegions::VExclusionRegions()
{
    setStarCatalogue( "Hipparcos_MAG8_1997.dat" );
}

VExclusionRegions::~VExclusionRegions()
{
    for( unsigned int i = 0; i < fExclusionRegions.size(); i++ )
    {
        if( fExclusionRegions[i] )
        {
            delete fExclusionRegions[i];
        }
    }
    fExclusionRegions.clear();
    for( unsigned int i = 0; i < fBrightStarSettings.size(); i++ )
    {
        if( fBrightStarSettings[i] )
        {
            delete fBrightStarSettings[i];
        }
    }
}

/*

   add a new bright star setting

*/
void VExclusionRegions::addBrightStarSettings( double iStarMinBrightness, double iStarExlusionRadius_DEG, string iStarBand )
{
    fBrightStarSettings.push_back( new VBrightStarExclusionSettings() );
    fBrightStarSettings.back()->fStarMinBrightness = iStarMinBrightness;
    if( iStarExlusionRadius_DEG > 0. )
    {
        fBrightStarSettings.back()->fStarExlusionRadius_DEG = iStarExlusionRadius_DEG;
    }
    if( iStarBand.size() > 0 )
    {
        fBrightStarSettings.back()->fStarBand = iStarBand;
    }
    sort( fBrightStarSettings.begin(), fBrightStarSettings.end() );
}

void VExclusionRegions::printBrightStarSettings()
{
    for( unsigned int i = 0; i < fBrightStarSettings.size(); i++ )
    {
        if( fBrightStarSettings[i] )
        {
            cout << "\tbright stars (magnitude brighter than ";
            cout << fBrightStarSettings[i]->fStarMinBrightness;
            cout << ", exclusion radius ";
            cout << fBrightStarSettings[i]->fStarExlusionRadius_DEG << " deg, ";
            cout << fBrightStarSettings[i]->fStarBand;
            cout << "-band) excluded" << endl;
        }
    }
}

/*

     read star catalogue and set exclusion regions according to
     bright star settings

*/
bool VExclusionRegions::setExclusionRegionsFromStarCatalogue( VStarCatalogue* iStarCatalogue )
{
    if( !iStarCatalogue )
    {
        cout << "VExclusionRegions::setExclusionRegionsFromStarCatalogue() error: star catalogue not found: ";
        cout << fStarCatalogue << endl;
        return false;
    }
    // remove double entries in the star catalogue
    iStarCatalogue->purge();
    if( iStarCatalogue->getListOfStarsinFOV().size() < 1 )
    {
        cout << "VExclusionRegions::setExclusionRegionsFromStarCatalogue(): ";
        cout << "no bright stars found in the field of view (";
        cout << fStarCatalogue << ")" << endl;
        return false;
    }
    
    //////////////////////////////////////////////////////////
    // print out star brightness requirements for exclusion region search
    printBrightStarSettings();
    
    //////////////////////////////////////////////////////////
    // set up list of exclusion regions from star catalogue
    double i_brightness = 100.;
    for( unsigned int i = 0; i < iStarCatalogue->getListOfStarsinFOV().size(); i++ )
    {
        // get list of stars in the relevant FOV
        if( !iStarCatalogue->getListOfStarsinFOV()[i] )
        {
            continue;
        }
        
        // compare this star with bright star settings
        for( unsigned int b = 0; b < fBrightStarSettings.size(); b++ )
        {
            if( !fBrightStarSettings[b] )
            {
                continue;
            }
            
            // get magnitude in correct band
            if( fBrightStarSettings[b]->fStarBand == "V" )
            {
                i_brightness = iStarCatalogue->getListOfStarsinFOV()[i]->fBrightness_V;
            }
            else if( fBrightStarSettings[b]->fStarBand == "B" )
            {
                i_brightness = iStarCatalogue->getListOfStarsinFOV()[i]->fBrightness_B;
            }
            else
            {
                cout << "VExclusionRegions::setExclusionRegions(): ";
                cout << "unknown band" << endl;
                continue;
            }
            
            // check brightness
            if( i_brightness < fBrightStarSettings[b]->fStarMinBrightness )
            {
                // check size of exclusion region
                if( fBrightStarSettings[b]->fStarExlusionRadius_DEG <= 0. )
                {
                    continue;
                }
                
                // check if a new exclusion region is required
                VListOfExclusionRegions* i_ExclusionRegion = 0;
                
                // check if this star is already in the list of exclusion regions
                for( unsigned int e = 0; e < fExclusionRegions.size(); e++ )
                {
                    if( fExclusionRegions[e]->fExcludeFromBackground_StarID >= 0
                            && fExclusionRegions[e]->fExcludeFromBackground_StarID == ( int )iStarCatalogue->getListOfStarsinFOV()[i]->fStarID )
                    {
                        i_ExclusionRegion = fExclusionRegions[e];
                    }
                }
                if( !i_ExclusionRegion )
                {
                    fExclusionRegions.push_back( new VListOfExclusionRegions() );
                    i_ExclusionRegion = fExclusionRegions.back();
                }
                // fill the exclusion region
                i_ExclusionRegion->fExcludeFromBackground_RAJ2000  = iStarCatalogue->getListOfStarsinFOV()[i]->fRA2000;
                i_ExclusionRegion->fExcludeFromBackground_DecJ2000 = iStarCatalogue->getListOfStarsinFOV()[i]->fDec2000;
                // exclusion radius might be already set to be larger
                if( fBrightStarSettings[b]->fStarExlusionRadius_DEG
                        > i_ExclusionRegion->fExcludeFromBackground_Radius1 )
                {
                    i_ExclusionRegion->fExcludeFromBackground_Radius1 = fBrightStarSettings[b]->fStarExlusionRadius_DEG;
                }
                if( fBrightStarSettings[b]->fStarExlusionRadius_DEG
                        > i_ExclusionRegion->fExcludeFromBackground_Radius2 )
                {
                    i_ExclusionRegion->fExcludeFromBackground_Radius2 = fBrightStarSettings[b]->fStarExlusionRadius_DEG;
                }
                // (will be filled later)
                i_ExclusionRegion->fExcludeFromBackground_North = 0.;
                i_ExclusionRegion->fExcludeFromBackground_West = 0.;
                i_ExclusionRegion->fExcludeFromBackground_CameraCentre_x = 0.;
                i_ExclusionRegion->fExcludeFromBackground_CameraCentre_y = 0.;
                i_ExclusionRegion->fExcludeFromBackground_StarID = ( int )iStarCatalogue->getListOfStarsinFOV()[i]->fStarID;
                i_ExclusionRegion->fExcludeFromBackground_StarName = iStarCatalogue->getListOfStarsinFOV()[i]->fStarName;
                i_ExclusionRegion->fExcludeFromBackground_StarBrightness_V = iStarCatalogue->getListOfStarsinFOV()[i]->fBrightness_V;
                i_ExclusionRegion->fExcludeFromBackground_StarBrightness_B = iStarCatalogue->getListOfStarsinFOV()[i]->fBrightness_B;
                if( fBrightStarSettings[b]->fStarBand == "B" )
                {
                    i_ExclusionRegion->fExcludeFromBackground_B_Band = true;
                }
                else
                {
                    i_ExclusionRegion->fExcludeFromBackground_B_Band = false;
                }
            }
        }
    }
    iStarCatalogue->purge();
    
    return true;
}

/*
 *  initialize exclusion regions
 *
 *  if star catalogue given: set exclusion regions due to bright stars
 *
 */
bool VExclusionRegions::initializeExclusionRegions( VStarCatalogue* iStarCatalogue,
        double iSkyMapCentre_ra_deg,
        double iSkyMapCentre_dec_deg,
        double iCameraCentre_ra_deg,
        double iCameraCentre_dec_deg )
{
    if( iStarCatalogue )
    {
        setExclusionRegionsFromStarCatalogue( iStarCatalogue );
    }
    
    // sort exclusion regions in magnitude
    sort( fExclusionRegions.begin(), fExclusionRegions.end() );
    
    // the following doesn't work if different sources are analyzed with RA <> 360 deg
    /////////////////////////////////////////////////////////
    // calculate camera coordinates for regions to exclude
    // (these are calculated relative to the sky map centre)
    for( unsigned int k = 0 ; k < fExclusionRegions.size(); k++ )
    {
        if( fExclusionRegions[k]->fExcludeFromBackground_DecJ2000 > -90. )
        {
            // set coordinates of exclusion regions
            // in sky map coordinates:
            fExclusionRegions[k]->fExcludeFromBackground_West  = VSkyCoordinatesUtilities::getTargetShiftWest( iSkyMapCentre_ra_deg, iSkyMapCentre_dec_deg,
                    fExclusionRegions[k]->fExcludeFromBackground_RAJ2000,
                    fExclusionRegions[k]->fExcludeFromBackground_DecJ2000 );
            fExclusionRegions[k]->fExcludeFromBackground_North = VSkyCoordinatesUtilities::getTargetShiftNorth( iSkyMapCentre_ra_deg, iSkyMapCentre_dec_deg,
                    fExclusionRegions[k]->fExcludeFromBackground_RAJ2000,
                    fExclusionRegions[k]->fExcludeFromBackground_DecJ2000 );
            if( TMath::Abs( fExclusionRegions[k]->fExcludeFromBackground_North ) < 1.e-4 )
            {
                fExclusionRegions[k]->fExcludeFromBackground_North = 0.;
            }
            if( TMath::Abs( fExclusionRegions[k]->fExcludeFromBackground_West ) < 1.e-4 )
            {
                fExclusionRegions[k]->fExcludeFromBackground_West  = 0.;
            }
        }
        else
        {
            // calculate ra/dec of exclusion regions
            VSkyCoordinatesUtilities::getWobbledDirection( fExclusionRegions[k]->fExcludeFromBackground_North,
                    -1.*fExclusionRegions[k]->fExcludeFromBackground_West,
                    iSkyMapCentre_dec_deg,
                    iSkyMapCentre_ra_deg,
                    fExclusionRegions[k]->fExcludeFromBackground_DecJ2000,
                    fExclusionRegions[k]->fExcludeFromBackground_RAJ2000 );
                    
        }
        // in coordinates relative to the camera centre
        fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_x  = VSkyCoordinatesUtilities::getTargetShiftWest( iCameraCentre_ra_deg, iCameraCentre_dec_deg,
                fExclusionRegions[k]->fExcludeFromBackground_RAJ2000,
                fExclusionRegions[k]->fExcludeFromBackground_DecJ2000 );
        fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_y = VSkyCoordinatesUtilities::getTargetShiftNorth( iCameraCentre_ra_deg, iCameraCentre_dec_deg,
                fExclusionRegions[k]->fExcludeFromBackground_RAJ2000,
                fExclusionRegions[k]->fExcludeFromBackground_DecJ2000 );
        if( TMath::Abs( fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_x ) < 1.e-4 )
        {
            fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_x = 0.;
        }
        if( TMath::Abs( fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_y ) < 1.e-4 )
        {
            fExclusionRegions[k]->fExcludeFromBackground_CameraCentre_y = 0.;
        }
    }
    
    return true;
}

/*
    add an exclusion region

*/
void VExclusionRegions::addExclusionRegion( double iNorth, double iWest,
        double iRa, double iDec,
        double iR1, double iR2, double iAngle,
        string iName )
{
    fExclusionRegions.push_back( new VListOfExclusionRegions() );
    fExclusionRegions.back()->fExcludeFromBackground_North = iNorth;
    fExclusionRegions.back()->fExcludeFromBackground_West = iWest;
    fExclusionRegions.back()->fExcludeFromBackground_RAJ2000 = iRa;
    fExclusionRegions.back()->fExcludeFromBackground_DecJ2000 = iDec;
    fExclusionRegions.back()->fExcludeFromBackground_Radius1 = iR1;
    fExclusionRegions.back()->fExcludeFromBackground_Radius2 = iR2;
    fExclusionRegions.back()->fExcludeFromBackground_RotAngle = iAngle;
    fExclusionRegions.back()->fExcludeFromBackground_StarName = iName;
    
}

void VExclusionRegions::printExclusionRegions()
{
    if( fExclusionRegions.size() == 0 )
    {
        return;
    }
    ios::fmtflags f( cout.flags() );
    cout << "\t region(s) excluded from background estimation (all values in [deg]): " << endl;
    for( unsigned int l = 0; l < fExclusionRegions.size(); l++ )
    {
        if( !fExclusionRegions[l] )
        {
            continue;
        }
        streamsize ss = std::cout.precision();
        cout << "\t";
        cout << "     region " << setw( 2 ) << l + 1 << ":" << std::fixed;
        cout << " N " << setw( 6 ) << setprecision( 3 ) << fExclusionRegions[l]->fExcludeFromBackground_North;
        cout << " W " << setw( 6 ) << setprecision( 3 ) << fExclusionRegions[l]->fExcludeFromBackground_West;
        cout << " Y " << setw( 6 ) << setprecision( 3 ) << fExclusionRegions[l]->fExcludeFromBackground_CameraCentre_y;
        cout << " X " << setw( 6 ) << setprecision( 3 ) << fExclusionRegions[l]->fExcludeFromBackground_CameraCentre_x;
        cout << " RAJ2000 " << setw( 6 ) << setprecision( 4 ) << fExclusionRegions[l]->fExcludeFromBackground_RAJ2000;
        cout << " DECJ2000 " << setw( 6 ) << setprecision( 4 ) << fExclusionRegions[l]->fExcludeFromBackground_DecJ2000;
        cout << " R1 " << setw( 4 ) << setprecision( 2 ) << fExclusionRegions[l]->fExcludeFromBackground_Radius1;
        cout << " R2 " << setw( 4 ) << setprecision( 2 ) << fExclusionRegions[l]->fExcludeFromBackground_Radius2;
        cout << " theta " << setw( 4 ) << setprecision( 2 ) << fExclusionRegions[l]->fExcludeFromBackground_RotAngle;
        cout << " ID " << fExclusionRegions[l]->fExcludeFromBackground_StarID;
        cout << " Bmag = " << fExclusionRegions[l]->fExcludeFromBackground_StarBrightness_B;
        if( fExclusionRegions[l]->fExcludeFromBackground_StarName.size() > 0 )
        {
            cout << " (" << fExclusionRegions[l]->fExcludeFromBackground_StarName << ")";
        }
        cout << endl;
        cout.precision( ss );
    }
    cout.flags( f );
}

bool VExclusionRegions::writeExclusionRegionTree( int iOnRun )
{
    char hname[200];
    if( iOnRun < 0 )
    {
        sprintf( hname, "user defined list of regions excluded from background calculation" );
    }
    else
    {
        sprintf( hname, "list of regions excluded from background calculation (run %d)", iOnRun );
    }
    TTree tEx( "tExcludedRegions", hname );
    float x = 0.;
    float y = 0.;
    float x_cam = 0.;
    float y_cam = 0.;
    float r1 = 0.;
    float r2 = 0.;
    float theta = 0;
    float decJ2000 = 0.;
    float raJ2000 = 0.;
    int id = 0;
    float iStarBrightness_V = 0;
    float iStarBrightness_B = 0;
    tEx.Branch( "x", &x, "x/F" );
    tEx.Branch( "y", &y, "y/F" );
    tEx.Branch( "x_cam", &x_cam, "x_cam/F" );
    tEx.Branch( "y_cam", &y_cam, "y_cam/F" );
    tEx.Branch( "r1", &r1, "r1/F" );
    tEx.Branch( "r2", &r2, "r2/F" );
    tEx.Branch( "theta", &theta, "theta/F" );
    tEx.Branch( "decj2000", &decJ2000, "decJ2000/F" );
    tEx.Branch( "raj2000", &raJ2000, "raJ2000/F" );
    tEx.Branch( "star_id", &id, "star_id/I" );
    tEx.Branch( "Vmag", &iStarBrightness_V, "Vmag/F" );
    tEx.Branch( "Bmag", &iStarBrightness_B, "Bmag/F" );
    
    for( unsigned int i = 0; i < fExclusionRegions.size(); i++ )
    {
        x = fExclusionRegions[i]->fExcludeFromBackground_West;
        y = fExclusionRegions[i]->fExcludeFromBackground_North;
        x_cam = fExclusionRegions[i]->fExcludeFromBackground_CameraCentre_x;
        y_cam = fExclusionRegions[i]->fExcludeFromBackground_CameraCentre_y;
        r1 = fExclusionRegions[i]->fExcludeFromBackground_Radius1;
        r2 = fExclusionRegions[i]->fExcludeFromBackground_Radius2;
        theta = fExclusionRegions[i]->fExcludeFromBackground_RotAngle;
        decJ2000 = fExclusionRegions[i]->fExcludeFromBackground_DecJ2000;
        raJ2000 = fExclusionRegions[i]->fExcludeFromBackground_RAJ2000;
        id = fExclusionRegions[i]->fExcludeFromBackground_StarID;
        iStarBrightness_V = fExclusionRegions[i]->fExcludeFromBackground_StarBrightness_V;
        iStarBrightness_B = fExclusionRegions[i]->fExcludeFromBackground_StarBrightness_B;
        
        tEx.Fill();
    }
    
    tEx.Write();
    
    return true;
}

bool VExclusionRegions::readExclusionRegionTree( TFile* f, int iRunOn )
{
    if( !f )
    {
        return false;
    }
    char hname[200];
    if( iRunOn < 0 )
    {
        sprintf( hname, "total_1/stereo/tExcludedRegions" );
    }
    else
    {
        sprintf( hname, "run_%d/stereo/tExcludedRegions", iRunOn );
    }
    
    TTree* tEx = ( TTree* )f->Get( hname );
    if( !tEx )
    {
        tEx = ( TTree* )f->Get( "tExcludedRegions" );
        if( !tEx )
        {
            return false;
        }
    }
    
    float x = 0.;
    float y = 0.;
    float x_cam = 0.;
    float y_cam = 0.;
    float r1 = 0.;
    float r2 = -99.;
    float theta = 0.;
    bool  bOldStyleExclusionRegions = false;
    float decJ2000 = 0.;
    float raJ2000 = 0.;
    int id = 0;
    float iV = 0.;
    float iB = 0.;
    tEx->SetBranchAddress( "x", &x );
    tEx->SetBranchAddress( "y", &y );
    if( tEx->GetBranchStatus( "x_cam" ) )
    {
        tEx->SetBranchAddress( "x_cam", &x_cam );
        tEx->SetBranchAddress( "y_cam", &y_cam );
    }
    // keep backwards compatibility to circular exclusion regions
    if( tEx->GetBranchStatus( "r" ) )
    {
        tEx->SetBranchAddress( "r", &r1 );
        theta = 0.;
        bOldStyleExclusionRegions = true;
    }
    else
    {
        tEx->SetBranchAddress( "r1", &r1 );
        tEx->SetBranchAddress( "r2", &r2 );
        tEx->SetBranchAddress( "theta", &theta );
    }
    tEx->SetBranchAddress( "star_id", &id );
    if( tEx->GetBranch( "decj2000" ) )
    {
        tEx->SetBranchAddress( "decj2000", &decJ2000 );
    }
    if( tEx->GetBranch( "raj2000" ) )
    {
        tEx->SetBranchAddress( "raj2000", &raJ2000 );
    }
    if( tEx->GetBranch( "Vmag" ) )
    {
        tEx->SetBranchAddress( "Vmag", &iV );
    }
    if( tEx->GetBranch( "Bmag" ) )
    {
        tEx->SetBranchAddress( "Bmag", &iB );
    }
    
    for( unsigned int i = 0; i < tEx->GetEntries(); i++ )
    {
        tEx->GetEntry( i );
        fExclusionRegions.push_back( new VListOfExclusionRegions() );
        fExclusionRegions.back()->fExcludeFromBackground_West = x;
        fExclusionRegions.back()->fExcludeFromBackground_North = y;
        fExclusionRegions.back()->fExcludeFromBackground_CameraCentre_x = x_cam;
        fExclusionRegions.back()->fExcludeFromBackground_CameraCentre_y = y_cam;
        fExclusionRegions.back()->fExcludeFromBackground_Radius1 = r1;
        if( bOldStyleExclusionRegions )
        {
            fExclusionRegions.back()->fExcludeFromBackground_Radius2 = r1;
        }
        else
        {
            fExclusionRegions.back()->fExcludeFromBackground_Radius2 = r2;
        }
        fExclusionRegions.back()->fExcludeFromBackground_RotAngle = theta;
        fExclusionRegions.back()->fExcludeFromBackground_StarID = id;
        fExclusionRegions.back()->fExcludeFromBackground_DecJ2000 = decJ2000;
        fExclusionRegions.back()->fExcludeFromBackground_RAJ2000  = raJ2000;
        fExclusionRegions.back()->fExcludeFromBackground_StarBrightness_V = iV;
        fExclusionRegions.back()->fExcludeFromBackground_StarBrightness_B = iB;
    }
    
    delete tEx;
    
    return true;
}

double VExclusionRegions::getLargestStarExlusionRadius()
{
    double iExRadius = 0.;
    
    for( unsigned int i = 0; i < fBrightStarSettings.size(); i++ )
    {
        if( fBrightStarSettings[i] && fBrightStarSettings[i]->fStarExlusionRadius_DEG > iExRadius )
        {
            iExRadius = fBrightStarSettings[i]->fStarExlusionRadius_DEG;
        }
    }
    
    return iExRadius;
}
