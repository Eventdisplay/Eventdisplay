//-*-mode:c++; mode:font-lock;-*-

/**
 * \file SphericalCoords.h
 * \ingroup SEphem
 * \brief This is a one-line description of this cpp file.
 *
 * Here is a tedious verbose multi-line description of all
 * the details of the code, more than you would
 * ever want to read. Generally, all the important documentation
 * goes in the .cpp files.
 *
 * Original Author: Stephen Fegan
 *
 **/

#ifndef SEPHEM_SPHERICALCOORDS_H
#define SEPHEM_SPHERICALCOORDS_H

#include<cassert>
#include<cmath>

#include"Angle.h"

namespace SEphem
{

    class SphericalCoords
    {
        public:
            SphericalCoords(): m_theta(), m_phi() {}
            SphericalCoords( const Angle& theta, const Angle& phi ):
                m_theta( theta ), m_phi( phi )
            {
                rationalize();
            }
            
            // ----------------------------------------------------------------------
            // Getters
            // ----------------------------------------------------------------------
            
            const Angle& theta() const
            {
                return m_theta;
            }
            const Angle& phi() const
            {
                return m_phi;
            }
            
            double thetaRad() const
            {
                return m_theta.rad();
            }
            double thetaDeg() const
            {
                return m_theta.deg();
            }
            double thetaRot() const
            {
                return m_theta.rot();
            }
            double phiRad() const
            {
                return m_phi.rad();
            }
            double phiDeg() const
            {
                return m_phi.deg();
            }
            double phiRot() const
            {
                return m_phi.rot();
            }
            double phiHrs() const
            {
                return m_phi.hrs();
            }
            
            Angle latitude() const
            {
                return m_theta.coAngle();
            }
            Angle longitude() const
            {
                return phi();
            }
            
            double latitudeRad() const
            {
                return m_theta.coAngleRadPM180();
            }
            double latitudeDeg() const
            {
                return m_theta.coAngleDegPM180();
            }
            double latitudeRot() const
            {
                return m_theta.coAngleRotPM180();
            }
            double longitudeRad() const
            {
                return phiRad();
            }
            double longitudeDeg() const
            {
                return phiDeg();
            }
            double longitudeRot() const
            {
                return phiRot();
            }
            
            bool isPole() const
            {
                return( ( m_theta == 0 ) || ( m_theta == Angle::sc_Pi ) );
            }
            
            double z() const
            {
                return cos( m_theta );
            }
            double x() const
            {
                return sin( m_theta ) * cos( m_phi );
            }
            double y() const
            {
                return sin( m_theta ) * sin( m_phi );
            }
            void cartesian( double& x, double& y, double& z ) const;
            
            Angle separation( const SphericalCoords& c ) const;
            Angle directionTo( const SphericalCoords& c ) const;
            Angle directionFrom( const SphericalCoords& c ) const;
            Angle compassDirectionTo( const SphericalCoords& c ) const;
            Angle compassDirectionFrom( const SphericalCoords& c ) const;
            void separationAndDirectionTo( const SphericalCoords& c,
                                           Angle& s, Angle& d ) const;
                                           
            // ----------------------------------------------------------------------
            // Setters
            // ----------------------------------------------------------------------
            
            void set( const Angle& theta, const Angle& phi )
            {
                m_theta = theta;
                m_phi = phi;
                rationalize();
            }
            void setRad( double theta, double phi )
            {
                m_theta = theta;
                m_phi = phi;
                rationalize();
            }
            void setDeg( double theta, double phi )
            {
                m_theta = Angle::frDeg( theta );
                m_phi = Angle::frDeg( phi );
                rationalize();
            }
            void setTheta( const Angle& theta )
            {
                set( theta, m_phi );
            }
            void setPhi( const Angle& phi )
            {
                set( m_theta, phi );
            }
            void setThetaPhi( const Angle& theta, const Angle& phi )
            {
                set( theta, phi );
            }
            void setThetaPhiRad( double theta, double phi )
            {
                setRad( theta, phi );
            }
            void setThetaPhiDeg( double theta, double phi )
            {
                setDeg( theta, phi );
            }
            void setLatLong( const Angle& latitude, const Angle& longitude )
            {
                set( Angle::sc_halfPi - latitude.rad(), longitude );
            }
            void setLatLongRad( double latitude, double longitude )
            {
                setRad( Angle::sc_halfPi - latitude, longitude );
            }
            void setLatLongDeg( double latitude, double longitude )
            {
                setDeg( 90 - latitude, longitude );
            }
            
            // ----------------------------------------------------------------------
            // Operations
            // ----------------------------------------------------------------------
            
            void rotate( const Angle& phi, const Angle& theta, const Angle& psi )
            {
                rotateRad( phi, theta, psi );
            }
            void rotate( const SphericalCoords& thetaPhi, const Angle& psi )
            {
                rotateRad( thetaPhi.phi(), thetaPhi.theta(), psi );
            }
            
            void rotateRad( double phi, double theta, double psi );
            void rotateDeg( double phi, double theta, double psi )
            {
                rotateRad( Angle::frDeg( phi ), Angle::frDeg( theta ), Angle::frDeg( psi ) );
            }
            void rotateRot( double phi, double theta, double psi )
            {
                rotateRad( Angle::frRot( phi ), Angle::frRot( theta ), Angle::frRot( psi ) );
            }
            
            // ----------------------------------------------------------------------
            // Wobble
            // ----------------------------------------------------------------------
            
            inline void wobble( const double theta_rad, const double phi_rad );
            inline bool wobble_inv( const double theta_rad, const double phi_rad,
                                    const double tol_rad = 1e-11 );
                                    
            // ----------------------------------------------------------------------
            // Static functions and constants
            // ----------------------------------------------------------------------
            
            static SphericalCoords make( Angle theta, Angle phi );
            static SphericalCoords makeRad( double theta, double phi );
            static SphericalCoords makeDeg( double theta, double phi );
            static SphericalCoords makeRot( double theta, double phi );
            
            static SphericalCoords makeLatLong( Angle latitude, Angle longitude );
            static SphericalCoords makeLatLongRad( double latitude, double longitude );
            static SphericalCoords makeLatLongDeg( double latitude, double longitude );
            static SphericalCoords makeLatLongRot( double latitude, double longitude );
            
        private:
            void rationalize();
            
            Angle m_theta;
            Angle m_phi;
    };                                            // SphericalCoords
    
    inline void SphericalCoords::rationalize()
    {
        if( m_theta.rad() > Angle::sc_Pi )
        {
            m_theta.setRad( m_theta.rad180() );
            m_phi.rotate( Angle::sc_Pi );
        }
    }
    
    inline void SphericalCoords::rotateRad( double phi, double theta, double psi )
    {
        m_phi.rotateRad( psi );
        const double sin_m_theta = sin( m_theta );
        const double x = sin_m_theta * cos( m_phi );
        const double y = sin_m_theta * sin( m_phi );
        const double z = cos( m_theta );
        const double cos_theta = cos( theta );
        const double sin_theta = sin( theta );
        const double xx = x * cos_theta + z * sin_theta;
        const double zz = z * cos_theta - x * sin_theta;
        // ATAN2 BETTEAR THAN ACOS FOR SMALL ANGLES
        set( atan2( sqrt( y * y + xx * xx ), zz ), atan2( y, xx ) );
        m_phi.rotateRad( phi );
    }
    
    inline Angle
    SphericalCoords::separation( const SphericalCoords& c ) const
    {
        // -------------------------------------------------------------------
        // ASTRONOMICAL ALGORITHMS RECOMMENDS NOT USING USUAL FORMULA
        // SINCE IT IS INACCURATE AT SMALL ANGLES, INSTEAD USE CARTESIANS
        // -------------------------------------------------------------------
#if HELL_FREEZES_OVER
        double cosang = cos( c.m_theta ) * cos( m_theta ) +
                        sin( c.m_theta ) * sin( m_theta ) * cos( c.m_phi - m_phi );
        return acos( cosang );
#else
        double sindphi = sin( c.m_phi - m_phi );
        double cosdphi = cos( c.m_phi - m_phi );
        double sinth1 = sin( c.m_theta );
        double costh1 = cos( c.m_theta );
        double sinth2 = sin( m_theta );
        double costh2 = cos( m_theta );
        double x = sinth2 * costh1 - costh2 * sinth1 * cosdphi;
        double y = sinth1 * sindphi;
        double z = costh2 * costh1 + sinth2 * sinth1 * cosdphi;
        return atan2( sqrt( x * x + y * y ), z );
#endif
    }
    
    inline Angle
    SphericalCoords::directionTo( const SphericalCoords& c ) const
    {
        double sindphi;
        double cosdphi;
        if( isPole() )
        {
            sindphi = sin( c.m_phi ), cosdphi = cos( c.m_phi );
        }
        else
        {
            sindphi = sin( c.m_phi - m_phi ), cosdphi = cos( c.m_phi - m_phi );
        }
        double sinth1 = sin( c.m_theta );
        double costh1 = cos( c.m_theta );
        double sinth2 = sin( m_theta );
        double costh2 = cos( m_theta );
        return atan2( -sinth1 * costh2 * cosdphi + costh1 * sinth2, sinth1 * sindphi );
    }
    
    inline void
    SphericalCoords::separationAndDirectionTo( const SphericalCoords& c,
            Angle& s, Angle& d ) const
    {
        double sindphi;
        double cosdphi;
        if( isPole() )
        {
            sindphi = sin( c.m_phi ), cosdphi = cos( c.m_phi );
        }
        else
        {
            sindphi = sin( c.m_phi - m_phi ), cosdphi = cos( c.m_phi - m_phi );
        }
        double sinth1 = sin( c.m_theta );
        double costh1 = cos( c.m_theta );
        double sinth2 = sin( m_theta );
        double costh2 = cos( m_theta );
        double x = sinth2 * costh1 - costh2 * sinth1 * cosdphi;
        double y = sinth1 * sindphi;
        double z = costh2 * costh1 + sinth2 * sinth1 * cosdphi;
        s = atan2( sqrt( x * x + y * y ), z );
        d = atan2( x, y );
    }
    
    inline Angle
    SphericalCoords::directionFrom( const SphericalCoords& c ) const
    {
        return c.directionTo( *this );
    }
    
    inline Angle
    SphericalCoords::compassDirectionTo( const SphericalCoords& c ) const
    {
        return directionTo( c ).coAngleRad();
    }
    
    inline Angle
    SphericalCoords::compassDirectionFrom( const SphericalCoords& c ) const
    {
        return directionFrom( c ).coAngleRad();
    }
    
    inline void
    SphericalCoords::cartesian( double& x, double& y, double& z ) const
    {
        z = cos( m_theta );
        double rho = sqrt( 1 - z * z );
        x = rho * cos( m_phi );
        y = rho * sin( m_phi );
    }
    
    inline void
    SphericalCoords::wobble( const double theta_rad, const double phi_rad )
    {
        SphericalCoords c( theta_rad, Angle::sc_Pi - phi_rad );
        c.rotate( phi(), theta(), 0 );
        *this = c;
    }
    
    inline bool
    SphericalCoords::wobble_inv( const double theta_rad, const double phi_rad,
                                 const double tol_rad )
    {
        SphericalCoords c_test( theta_rad, -phi_rad );
        c_test.rotateRad( phiRad(), thetaRad(), 0 );
        
        SphericalCoords c_test_wob( theta_rad, Angle::sc_Pi - phi_rad );
        c_test_wob.rotateRad( c_test.phiRad(), c_test.thetaRad(), 0 );
        
        unsigned nrounds = 50;
        while( separation( c_test_wob ).rad() > tol_rad )
        {
            if( nrounds-- == 0 )
            {
                return false;
            }
            
            c_test_wob.rotateRad( 0, -thetaRad(), -phiRad() );
            c_test_wob.rotateRad( c_test.phiRad(), c_test.thetaRad(), M_PI );
            
            c_test = c_test_wob;
            
            c_test_wob.setRad( theta_rad, Angle::sc_Pi - phi_rad );
            c_test_wob.rotateRad( c_test.phiRad(), c_test.thetaRad(), 0 );
            
#if 0
            std::cout << phi().hmsString( 4 ) << ' '
                      << latitude().dmsString( 4 ) << ' '
                      << c_test.phi().hmsString( 4 ) << ' '
                      << c_test.latitude().dmsString( 4 ) << ' '
                      << c_test_wob.phi().hmsString( 4 ) << ' '
                      << c_test_wob.latitude().dmsString( 4 ) << '\n';
#endif
        }
        
        *this = c_test;
        return true;
    }
    
    inline SphericalCoords SphericalCoords::make( Angle theta, Angle phi )
    {
        return SphericalCoords( theta, phi );
    }
    
    inline SphericalCoords
    SphericalCoords::makeRad( double theta, double phi )
    {
        return SphericalCoords( Angle::makeRad( theta ), Angle::makeRad( phi ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeDeg( double theta, double phi )
    {
        return SphericalCoords( Angle::makeDeg( theta ), Angle::makeDeg( phi ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeRot( double theta, double phi )
    {
        return SphericalCoords( Angle::makeRot( theta ), Angle::makeRot( phi ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeLatLong( Angle latitude, Angle longitude )
    {
        return SphericalCoords( Angle::makeCoAngleRad( latitude ),
                                Angle::makeRad( longitude ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeLatLongRad( double latitude, double longitude )
    {
        return SphericalCoords( Angle::makeCoAngleRad( latitude ),
                                Angle::makeRad( longitude ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeLatLongDeg( double latitude, double longitude )
    {
        return SphericalCoords( Angle::makeCoAngleDeg( latitude ),
                                Angle::makeDeg( longitude ) );
    }
    
    inline SphericalCoords
    SphericalCoords::makeLatLongRot( double latitude, double longitude )
    {
        return SphericalCoords( Angle::makeCoAngleRot( latitude ),
                                Angle::makeRot( longitude ) );
    }
    
}                                                 // namespace SEphem
#endif                                            // SEPHEM_SPHERICALCOORDS_H
