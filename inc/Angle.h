//-*-mode:c++; mode:font-lock;-*-

/**
 * \file Angle.h
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

#ifndef SEPHEM_ANGLE_H
#define SEPHEM_ANGLE_H

#include<cmath>
#include<string>

#define ANGLE_DEGPERRAD  57.2957795130823208768
#define ANGLE_ROTPERRAD  0.15915494309189533576
#define ANGLE_HRSPERRAD  3.81971863420548805845
#define ANGLE_RADPERDEG  0.01745329251994329576
#define ANGLE_RADPERROT  6.28318530717958647692
#define ANGLE_RADPERHRS  0.26179938779914943653

namespace SEphem
{

    class Angle
    {
        public:
        
            // Constructor, Destructor and Assignment
            
            Angle(): m_angle() {}
            Angle( double a ): m_angle( a )
            {
                rationalize();
            }
            Angle( const Angle& a ): m_angle( a.m_angle ) { }
            Angle& operator=( const Angle& a )
            {
                m_angle = a.m_angle;
                return *this;
            }
            
            // ----------------------------------------------------------------------
            // Getters
            // ----------------------------------------------------------------------
            
            operator double() const
            {
                return m_angle;
            }
            
            double rad() const
            {
                return m_angle;
            }
            double radians() const
            {
                return rad();
            }
            double deg() const
            {
                return toDeg( rad() );
            }
            double degrees() const
            {
                return deg();
            }
            double rot() const
            {
                return toRot( rad() );
            }
            double rotations() const
            {
                return rot();
            }
            double hrs() const
            {
                return toHrs( rad() );
            }
            
            double radPM() const
            {
                return toRadPM( m_angle );
            }
            double radiansPM() const
            {
                return radPM();
            }
            double degPM() const
            {
                return toDeg( radPM() );
            }
            double degreesPM() const
            {
                return degPM();
            }
            double rotPM() const
            {
                return toRot( radPM() );
            }
            double rotationsPM() const
            {
                return radPM();
            }
            double hrsPM() const
            {
                return toHrs( radPM() );
            }
            
            double rad180() const
            {
                return toRad180( m_angle );
            }
            double radians180() const
            {
                return rad180();
            }
            double deg180() const
            {
                return toDeg( rad180() );
            }
            double degrees180() const
            {
                return deg180();
            }
            double rot180() const
            {
                return toRot( rad180() );
            }
            double rotations180() const
            {
                return rad180();
            }
            double hrs180() const
            {
                return toHrs( rad180() );
            }
            
            double radPM180() const
            {
                return toRadPM180( m_angle );
            }
            double radiansPM180() const
            {
                return radPM180();
            }
            double degPM180() const
            {
                return toDeg( radPM180() );
            }
            double degreesPM180() const
            {
                return degPM180();
            }
            double rotPM180() const
            {
                return toRot( radPM180() );
            }
            double rotationsPM180() const
            {
                return radPM180();
            }
            double hrsPM180() const
            {
                return toHrs( radPM180() );
            }
            
            Angle coAngle() const
            {
                return coAngleRad();
            }
            
            double coAngleRad() const
            {
                return toCoAngle( m_angle );
            }
            double coAngleDeg() const
            {
                return toDeg( coAngleRad() );
            }
            double coAngleRot() const
            {
                return toRot( coAngleRad() );
            }
            double coAngleHrs() const
            {
                return toHrs( coAngleRad() );
            }
            
            double coAngleRadPM() const
            {
                return toRadPM( toCoAngle( m_angle ) );
            }
            double coAngleDegPM() const
            {
                return toDeg( coAngleRadPM() );
            }
            double coAngleRotPM() const
            {
                return toRot( coAngleRadPM() );
            }
            double coAngleHrsPM() const
            {
                return toHrs( coAngleRadPM() );
            }
            
            double coAngleRad180() const
            {
                return toRad180( toCoAngle( m_angle ) );
            }
            double coAngleDeg180() const
            {
                return toDeg( coAngleRad180() );
            }
            double coAngleRot180() const
            {
                return toRot( coAngleRad180() );
            }
            double coAngleHrs180() const
            {
                return toHrs( coAngleRad180() );
            }
            
            double coAngleRadPM180() const
            {
                return toRadPM180( toCoAngle( m_angle ) );
            }
            double coAngleDegPM180() const
            {
                return toDeg( coAngleRadPM180() );
            }
            double coAngleRotPM180() const
            {
                return toRot( coAngleRadPM180() );
            }
            double coAngleHrsPM180() const
            {
                return toHrs( coAngleRadPM180() );
            }
            
            void hms( unsigned& h, unsigned& m, unsigned& s, unsigned& f,
                      unsigned sec_digits = 1 ) const;
            void dmsPM180( bool& negative,
                           unsigned& d, unsigned& m, unsigned& s, unsigned& f,
                           unsigned sec_digits = 1 ) const;
            void dmsPM360( bool& negative,
                           unsigned& d, unsigned& m, unsigned& s, unsigned& f,
                           unsigned sec_digits = 1 ) const;
                           
            std::string hmsString( unsigned sec_digits = 1, bool hmsSep = false ) const;
            std::string hmsPMString( unsigned sec_digits = 1, bool hmsSep = false ) const;
            std::string dmsString( unsigned sec_digits = 1, bool dmsSep = false ) const
            {
                return dmsPM180String( sec_digits, dmsSep );
            }
            std::string dmsPM180String( unsigned sec_digits = 1, bool dmsSep = false )
            const;
            std::string dmsPM360String( unsigned sec_digits = 1, bool dmsSep = false )
            const;
            std::string degString( unsigned dec_digits = 1 ) const;
            std::string degPMString( unsigned dec_digits = 1 ) const;
            std::string deg180String( unsigned dec_digits = 1 ) const;
            std::string degPM180String( unsigned dec_digits = 1 ) const;
            
            Angle separation( const Angle& a ) const
            {
                Angle b = a - m_angle;
                return b;
            }
            double separationRad( const Angle& a ) const
            {
                return toRad180( m_angle - a.rad() );
            }
            double separationDeg( const Angle& a ) const
            {
                return toDeg( separationRad( a ) );
            }
            double separationRot( const Angle& a ) const
            {
                return toRot( separationRad( a ) );
            }
            double separationHrs( const Angle& a ) const
            {
                return toHrs( separationRad( a ) );
            }
            
            double x() const
            {
                return cos( m_angle );
            }
            double y() const
            {
                return sin( m_angle );
            }
            void cartesian( double& x, double& y ) const
            {
                x = Angle::x();
                y = Angle::y();
            }
            
            unsigned long bar( unsigned radix, bool round = false ) const
            {
                return toBAR( rad(), radix, round );
            }
            
            // ----------------------------------------------------------------------
            // Setters
            // ----------------------------------------------------------------------
            
            void setRad( double a )
            {
                m_angle = a;
                rationalize();
            }
            void setRadians( double a )
            {
                setRad( a );
            }
            void setDeg( double a )
            {
                setRad( frDeg( a ) );
            }
            void setDegrees( double a )
            {
                setDeg( a );
            }
            void setRot( double a )
            {
                setRad( frRot( a ) );
            }
            void setRotations( double a )
            {
                setRot( a );
            }
            void setHrs( double a )
            {
                setRad( frHrs( a ) );
            }
            
            void setCoAngleRad( double a )
            {
                m_angle = frCoAngle( a );
                rationalize();
            }
            void setCoAngleDeg( double a )
            {
                setCoAngleRad( frDeg( a ) );
            }
            void setCoAngleRot( double a )
            {
                setCoAngleRad( frRot( a ) );
            }
            void setCoAngleHrs( double a )
            {
                setCoAngleRad( frHrs( a ) );
            }
            
            bool setFromHMSString( std::string str );
            bool setFromDMSString( std::string str );
            inline bool setCoAngleFromDMSString( std::string str );
            
            void setBAR( unsigned long a, unsigned radix, bool round = false )
            {
                setRad( frBAR( a, radix, round ) );
            }
            
            // ----------------------------------------------------------------------
            // Operations
            // ----------------------------------------------------------------------
            
            void rotate( Angle a )
            {
                m_angle += a.m_angle;
                if( m_angle > sc_twoPi )
                {
                    m_angle -= sc_twoPi;
                }
            }
            void rotateRad( double a )
            {
                m_angle += a;
                rationalize();
            }
            void rotateDeg( double a )
            {
                rotateRad( frDeg( a ) );
            }
            void rotateHrs( double a )
            {
                rotateRad( frHrs( a ) );
            }
            
            Angle& operator+= ( const Angle a )
            {
                m_angle += a.m_angle;
                rationalize();
                return *this;
            }
            Angle& operator-= ( const Angle a )
            {
                m_angle -= a.m_angle;
                rationalize();
                return *this;
            }
            
            Angle& operator+= ( double a )
            {
                m_angle += a;
                rationalize();
                return *this;
            }
            Angle& operator-= ( double a )
            {
                m_angle -= a;
                rationalize();
                return *this;
            }
            
            // ----------------------------------------------------------------------
            // Static functions and constants
            // ----------------------------------------------------------------------
            
            static double frDeg( double a )
            {
                return a * sc_radPerDeg;
            }
            static double frRot( double a )
            {
                return a * sc_radPerRot;
            }
            static double frHrs( double a )
            {
                return a * sc_radPerHrs;
            }
            static double fromDeg( double a )
            {
                return frDeg( a );
            }
            static double fromRot( double a )
            {
                return frRot( a );
            }
            static double fromHrs( double a )
            {
                return frHrs( a );
            }
            static double toDeg( double a )
            {
                return a * sc_degPerRad;
            }
            static double toRot( double a )
            {
                return a * sc_rotPerRad;
            }
            static double toHrs( double a )
            {
                return a * sc_hrsPerRad;
            }
            
            static unsigned long toBAR( double a, unsigned radix, bool round = false );
            static double frBAR( unsigned long a, unsigned radix, bool round = false );
            
            static double frCoAngle( double a )
            {
                return a <= sc_halfPi ? sc_halfPi - a : sc_5HalfPi - a;
            }
            static double fromCoAngle( double a )
            {
                return frCoAngle( a );
            }
            static double toCoAngle( double a )
            {
                return frCoAngle( a );
            }
            
            inline static double toRadPM( double a );
            inline static double toRad180( double a );
            inline static double toRadPM180( double a );
            
#if 0
            static double frHMSString( const std::string& a );
            static double frDMSString( const std::string& a );
            static std::string toHMSString( double a );
            static std::string toDMSString( double a );
#endif
            
            static Angle makeRad( double a )
            {
                return Angle( a );
            }
            static Angle makeDeg( double a )
            {
                return Angle( frDeg( a ) );
            }
            static Angle makeRot( double a )
            {
                return Angle( frRot( a ) );
            }
            static Angle makeHrs( double a )
            {
                return Angle( frHrs( a ) );
            }
            
            static Angle makeCoAngleRad( double a )
            {
                return Angle( toCoAngle( a ) );
            }
            static Angle makeCoAngleDeg( double a )
            {
                return makeCoAngleRad( frDeg( a ) );
            }
            static Angle makeCoAngleRot( double a )
            {
                return makeCoAngleRad( frRot( a ) );
            }
            static Angle makeCoAngleHrs( double a )
            {
                return makeCoAngleRad( frHrs( a ) );
            }
            
            static Angle makeBAR( unsigned long a, unsigned radix, bool round = false )
            {
                return Angle( frBAR( a, radix, round ) );
            }
            
            static void rotateCartesians( const Angle& a, double& x, double& y );
            
            static const double sc_halfPi;
            static const double sc_Pi;
            static const double sc_3HalfPi;
            static const double sc_twoPi;
            static const double sc_5HalfPi;
            static const double sc_degPerRad;
            static const double sc_rotPerRad;
            static const double sc_hrsPerRad;
            static const double sc_radPerDeg;
            static const double sc_radPerRot;
            static const double sc_radPerHrs;
            
        private:
            double m_angle;
            
            inline void rationalize();
    };                                            // class Angle
    
    inline void Angle::rationalize()
    {
        if( m_angle < 0 )
        {
            m_angle = fmod( m_angle, sc_twoPi ) + sc_twoPi;
        }
        else if( m_angle >= sc_twoPi )
        {
            m_angle = fmod( m_angle, sc_twoPi );
        }
    }
    
    inline double Angle::toRadPM( double a )
    {
        if( a <= sc_Pi )
        {
            return a;
        }
        else
        {
            return a - sc_twoPi;
        }
    }
    
    inline double Angle::toRad180( double a )
    {
        if( a <= sc_Pi )
        {
            return a;
        }
        else
        {
            return sc_twoPi - a;
        }
    }
    
    inline double Angle::toRadPM180( double a )
    {
        if( a <= sc_halfPi )
        {
            return a;
        }
        else if( a <= sc_3HalfPi )
        {
            return sc_Pi - a;
        }
        else
        {
            return a - sc_twoPi;
        }
    }
    
    inline void Angle::rotateCartesians( const Angle& a, double& x, double& y )
    {
        double c = cos( a.rad() );
        double s = sin( a.rad() );
        double t = x * c - y * s;
        y = y * c + x * s;
        x = t;
    }
    
    inline bool Angle::setCoAngleFromDMSString( std::string str )
    {
        Angle a;
        if( a.setFromDMSString( str ) )
        {
            setCoAngleRad( a.rad() );
        }
        else
        {
            return false;
        }
        return true;
    }
    
}                                                 // namespace SEphem
#endif                                            // SEPHEM_ANGLE_H
