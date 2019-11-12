//-*-mode:c++; mode:font-lock;-*-

/**
 * \file CorrectionParameters.cpp
 * \ingroup SEphem
 * \brief This is a one-line description of this cpp file.
 *
 * Here is a tedious verbose multi-line description of all
 * the details of the code, more than you would
 * ever want to read. Generally, all the important documentation
 * goes in the .cpp files.
 *
 * Original Author: Stephen Fegan
 * $Date: 2010/03/08 07:39:54 $
 * $Tag$
 *
 **/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<ctime>

#include"CorrectionParameters.h"

using namespace SEphem;

const double CorrectionParameters::sc_lim_az_cw =     270.0 * ANGLE_RADPERDEG;
const double CorrectionParameters::sc_lim_az_cc =    -270.0 * ANGLE_RADPERDEG;
const double CorrectionParameters::sc_inversion_tol = 0.0001 * ANGLE_RADPERDEG;
const int    CorrectionParameters::sc_inversion_it_max = 20;

bool CorrectionParameters::
doAzElCorrections( double& az_driveangle, double& el_driveangle,
                   const double& tel_az_driveangle, bool do_corrections ) const
{
    // In:  az_driveangle     - azimuth of source [0 to 2pi]
    //      el_driveangle     - elevation of source [0 to pi/2]
    //      tel_az_driveangle - scope azimuth drive angle [-3pi/2 to 3pi/2]
    //      do_corrections    - enable corrections
    //
    // Out: az_driveangle     - target azimuth drive angle [-3pi/2 to 3pi/2]
    //      el_driveangle     - target elevation drive angle
    //
    // Return: true -  target can be pointed at
    //         false - target cannot be pointed at due to mis-alignments
    
#if 0
    Debug::stream()
            << "DO: In "
            << Angle::toDeg( el_driveangle ) << ' '
            << Angle::toDeg( az_driveangle ) << std::endl;
#endif
            
    double az = az_driveangle;
    double el = el_driveangle;
    
    double az1;
    double az2;
    
    if( ( enable_offsets ) && ( do_corrections ) && ( enable_corrections ) )
    {
        // ----------------------------------------------------------------------
        // 1: Compensate for flexure in the mount
        // ----------------------------------------------------------------------
        
        // I can see this perhaps been more appropriate as correction 4,
        // but lets see how we get on with it here. Also there is an
        // elevation effect on azimuth, and an effect /of/ azimuth on
        // both elevation and azimuth but I expect them to be small.
        
        // zero at zenith
        el = el - ( flex_el_A * cos( el ) + flex_el_B * sin( 2 * el ) );
        
        // Transform to Cartesians for next corrections
        double ud = sin( el );
        double ns = cos( el ) * cos( az );
        double ew = cos( el ) * sin( az );
        
        // ----------------------------------------------------------------------
        // 2: Rotate to account for non-level AZ plane
        // ----------------------------------------------------------------------
        
        Angle::rotateCartesians( -az_ew, ew, ud );
        Angle::rotateCartesians( -az_ns, ns, ud );
        
        // ----------------------------------------------------------------------
        // 3: Account for focal plane mis-alignment perp to the EL drive axis
        // 4: EL drive axis mis-alignment in the UD-EW plane
        // ----------------------------------------------------------------------
        
        // I cannot think of a way of inverting 3/4 without doing it
        // iteratively in the forwards direction.
        
        // Initial guess: just do the transformations in
        // reverse.. should get us close to the correct point.
        az = atan2( ew, ns );
        el = atan2( ud, sqrt( ns * ns + ew * ew ) );
        
        ud = 0.0;
        ns = cos( -fp_az );
        ew = sin( -fp_az );
        
        Angle::rotateCartesians( el, ns, ud );
        Angle::rotateCartesians( -el_udew, ew, ud );
        Angle::rotateCartesians( az, ns, ew );
        
        // xaz and xel are estimates of the drive az and el which would
        // give the real az and el if they were uncorrected using the
        // undoCorrections routines below
        double xaz = atan2( ew, ns );
        double xel = atan2( ud, sqrt( ew * ew + ns * ns ) );
        
        SphericalCoords here = SphericalCoords::makeLatLongRad( el, az );
        SphericalCoords there;
        
        int i = 0;
        do
        {
            // ---------- This code copied from undoCorrections below -----------
            ud = 0;
            ns = cos( fp_az );
            ew = sin( fp_az );
            Angle::rotateCartesians( xel, ns, ud );
            Angle::rotateCartesians( el_udew, ew, ud );
            Angle::rotateCartesians( xaz, ns, ew );
            // ------------------------------------------------------------------
            
            there = SphericalCoords::makeLatLongRad( atan2( ud, sqrt( ew * ew + ns * ns ) ),
                    atan2( ew, ns ) );
                    
#if 0
            Debug::stream()
                    << i << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << Angle::toDeg( el ) << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << Angle::toDeg( az ) << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << there.latitudeDeg() << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << there.longitudeDeg() << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << Angle::toDeg( xel ) << ' '
                    << std::fixed << std::setw( 7 ) << std::setprecision( 4 )
                    << Angle::toDeg( xaz ) << std::endl;
#endif
                    
            xaz = Angle::makeRad( xaz - ( there.longitudeRad() - az ) );
            xel = Angle::makeRad( xel - ( there.latitudeRad() - el ) ).radPM180();
            i++;
        }
        while( ( here.separation( there ) > sc_inversion_tol ) &&
                ( i < sc_inversion_it_max ) );
                
        if( here.separation( there ) > sc_inversion_tol )
        {
            return false;
        }
        
        az = xaz;
        el = xel;
        
        // ----------------------------------------------------------------------
        // 5: Scale for gear ratio
        // ----------------------------------------------------------------------
        
        // CW & CCW cases must be handled seperately in scaling
        az1 = fmod( fmod( az, Angle::sc_twoPi ) + Angle::sc_twoPi, Angle::sc_twoPi );
        az2 = az1 - Angle::sc_twoPi;
        
        az1 = az / az_ratio;
        az2 = az / az_ratio;
        el = el / el_ratio;
    }
    else
    {
        az1 = fmod( fmod( az, Angle::sc_twoPi ) + Angle::sc_twoPi, Angle::sc_twoPi );
        az2 = az1 - Angle::sc_twoPi;
    }
    
    // --------------------------------------------------------------------------
    // 6: Subtract Offsets
    // --------------------------------------------------------------------------
    
    if( ( enable_offsets ) && ( do_corrections ) )
    {
        az1 = az1 - az_offset;
        az2 = az2 - az_offset;
        el = el - el_offset;
    }
    
    az1 = fmod( fmod( az1, Angle::sc_twoPi ) + Angle::sc_twoPi, Angle::sc_twoPi );
    az2 = az1 - Angle::sc_twoPi;
    
    // never happen
    if( ( az1 > sc_lim_az_cw ) && ( az2 < sc_lim_az_cc ) )
    {
        return false;
    }
    else if( az1 > sc_lim_az_cw )
    {
        az = az2;
    }
    else if( az2 < sc_lim_az_cc )
    {
        az = az1;
    }
    else if( fabs( az1 - tel_az_driveangle ) <= fabs( az2 - tel_az_driveangle ) )
    {
        az = az1;
    }
    else
    {
        az = az2;
    }
    
#if 0
    if( el > Angle::sc_halfPi )      -- THIS IS WRONG WE JUST HAVE TO ALLOW THE
    {
        -- ELEVATION TO GO GREATER THAN 90 AND LESS
        el = Angle::sc_Pi - el;
        -- THAN ZERO OR PREFERABLY USE THE INBUILT
        az = az + Angle::sc_Pi;
        -- OFFSETS
    }
#endif
    
    el_driveangle = el;
    az_driveangle = az;
    
#if 0
    Debug::stream()
            << "DO: Out: "
            << Angle::toDeg( el_driveangle ) << ' '
            << Angle::toDeg( az_driveangle ) << std::endl;
#endif
            
    return true;
}


void CorrectionParameters::
undoAzElCorrections( double& az_driveangle, double& el_driveangle,
                     bool do_corrections ) const
{
    if( ( !enable_offsets ) || ( !do_corrections ) )
    {
        SphericalCoords sc =
            SphericalCoords::makeLatLongRad( el_driveangle, az_driveangle );
        az_driveangle = sc.phiRad();
        el_driveangle = sc.latitudeRad();
        return;
    }
    
#if 0
    Debug::stream()
            << "UNDO: In: "
            << Angle::toDeg( el_driveangle ) << ' '
            << Angle::toDeg( az_driveangle ) << std::endl;
#endif
            
    // --------------------------------------------------------------------------
    // 1: Add offsets
    // --------------------------------------------------------------------------
    
    double az = az_driveangle + az_offset;
    double el = el_driveangle + el_offset;
    
#if 0
    if( el > Angle::sc_halfPi )      -- THIS IS WRONG WE JUST HAVE TO ALLOW THE
    {
        -- ELEVATION TO GO GREATER THAN 90 AND LESS
        el = Angle::sc_Pi - el;
        -- THAN ZERO OR PREFERABLY USE THE INBUILT
        az = az + Angle::sc_Pi;
        -- OFFSETS
    }
#endif
    
    if( enable_corrections )
    {
        // ----------------------------------------------------------------------
        // 2: Remove gear ratio scaling
        // ----------------------------------------------------------------------
        
        az = fmod( az * az_ratio + Angle::sc_twoPi, Angle::sc_twoPi );
        el = el * el_ratio;
        
        // ----------------------------------------------------------------------
        // 3: EL drive axis mis-alignment in the UD-EW plane
        // 4: Focal plane mis-alignment perp to the EL drive axis
        // ----------------------------------------------------------------------
        
        // Pretend scope is at 0,0 -- find direction of the focal plane
        // mis-aligment and rotate it up to the El
        double ud = 0;
        double ns = cos( fp_az );
        double ew = sin( fp_az );
        Angle::rotateCartesians( el, ns, ud );
        
        // Rotate coords down in the UD-EW plane.
        Angle::rotateCartesians( el_udew, ew, ud );
        
        // Rotate to correct AZ.
        Angle::rotateCartesians( az, ns, ew );
        
        // ----------------------------------------------------------------------
        // 5: Rotate to account for non-level AZ plane
        // ----------------------------------------------------------------------
        
        Angle::rotateCartesians( az_ns, ns, ud );
        Angle::rotateCartesians( az_ew, ew, ud );
        
        // Back to Spherical coordinates
        az = atan2( ew, ns );
        el = atan2( ud, sqrt( ew * ew + ns * ns ) );
        
        // ----------------------------------------------------------------------
        // 6: Compensate for flexure in the mount
        // ----------------------------------------------------------------------
        
        // I can see this perhaps been more appropriate as correction 3,
        // but lets see how we get on with it here. Also there is an
        // elevation effect on azimuth, and an effect /of/ azimuth on
        // both elevation and azimuth but I expect them to be small.
        
        double real_el = el;
        double last_el;
        
        int i = 0;
        do
        {
            last_el = el;
            el = real_el + ( flex_el_A * cos( el ) + flex_el_B * sin( 2 * el ) );
            i++;
        }
        while( ( fabs( last_el - el ) > sc_inversion_tol ) &&
                ( i < sc_inversion_it_max ) );
    }
    
    SphericalCoords sc = SphericalCoords::makeLatLongRad( el, az );
    az_driveangle = sc.phiRad();
    el_driveangle = sc.latitudeRad();
    
#if 0
    Debug::stream()
            << "UNDO: Out: "
            << Angle::toDeg( el_driveangle ) << ' '
            << Angle::toDeg( az_driveangle ) << std::endl;
#endif
}


bool CorrectionParameters::save( const char* filename )
const
{
    std::ofstream stream( filename );
    if( !stream )
    {
        return false;
    }
    stream << enable_offsets << std::endl
           << enable_corrections << std::endl
           << az_ratio << std::endl
           << el_ratio << std::endl
           << az_offset << std::endl
           << el_offset << std::endl
           << az_ns << std::endl
           << az_ew << std::endl
           << el_udew << std::endl
           << fp_az << std::endl
           << flex_el_A << std::endl
           << flex_el_B << std::endl
           << enable_vff << std::endl
           << el_pos_vff_s << std::endl
           << el_pos_vff_t << std::endl
           << el_neg_vff_s << std::endl
           << el_neg_vff_t << std::endl
           << az_pos_vff_s << std::endl
           << az_pos_vff_t << std::endl
           << az_neg_vff_s << std::endl
           << az_neg_vff_t << std::endl;
    return true;
}


bool CorrectionParameters::load( const char* filename )
{
    std::ifstream stream( filename );
    if( !stream )
    {
        return false;
    }
    stream >> enable_offsets
           >> enable_corrections
           >> az_ratio
           >> el_ratio
           >> az_offset
           >> el_offset
           >> az_ns
           >> az_ew
           >> el_udew
           >> fp_az
           >> flex_el_A
           >> flex_el_B
           >> enable_vff
           >> el_pos_vff_s
           >> el_pos_vff_t
           >> el_neg_vff_s
           >> el_neg_vff_t
           >> az_pos_vff_s
           >> az_pos_vff_t
           >> az_neg_vff_s
           >> az_neg_vff_t;
    if( !stream )
    {
        return false;
    }
    return true;
}


std::string CorrectionParameters::loadFilename( unsigned scope_num )
{
    std::ostringstream stream;
    stream << "corrections_t" << scope_num + 1 << ".dat";
    return stream.str();
}


std::string CorrectionParameters::saveFilename( unsigned scope_num )
{
    time_t t = time( 0 );
    struct tm* tm;
    tm = gmtime( &t );
    std::ostringstream stream;
    stream << "corrections_t" << scope_num + 1 << '_'
           << std::setw( 2 ) << std::setfill( '0' ) << tm->tm_year % 100
           << std::setw( 2 ) << std::setfill( '0' ) << tm->tm_mon + 1
           << std::setw( 2 ) << std::setfill( '0' ) << tm->tm_mday << ".dat";
    return stream.str();
}



#ifdef TESTMAIN
int main( int argc, char** argv )
{
    CorrectionParameters cp;
    cp.fp_az = Angle::frDeg( 0.35 );
    
    double tar_el = 0;
    double tar_az = 0;
    for( tar_el = Angle::frDeg( 85 );
            tar_el <= Angle::frDeg( 90.0001 );
            tar_el += Angle::frDeg( 0.100000001 ) )
    {
        double tel_el = tar_el;
        double tel_az = tar_az;
        bool good = cp.doAzElCorrections( tel_az, tel_el, tar_az, true );
        
        double unc_el = tel_el;
        double unc_az = tel_az;
        cp.undoAzElCorrections( unc_az, unc_el, true );
        
        Debug::stream()
                << std::fixed << std::showpos
                << good << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( tar_el ) << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( tar_az ) << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( tel_el ) << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( tel_az ) << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( unc_el ) << ' '
                << std::setw( 9 ) << std::setprecision( 5 )
                << Angle::toDeg( unc_az ) << std::endl;
    }
}
#endif
