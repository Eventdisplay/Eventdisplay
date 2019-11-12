//-*-mode:c++; mode:font-lock;-*-

/**
 * \file Angle.cpp
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

#include<string>
#include<sstream>
#include<iostream>
#include<iomanip>

#include"Angle.h"

const double SEphem::Angle::sc_halfPi = 1.57079632679489661923;
const double SEphem::Angle::sc_Pi = 3.14159265358979323846;
const double SEphem::Angle::sc_3HalfPi = 4.71238898038468985769;
const double SEphem::Angle::sc_twoPi = 6.28318530717958647692;
const double SEphem::Angle::sc_5HalfPi = 7.85398163397448309615;
const double SEphem::Angle::sc_degPerRad = ANGLE_DEGPERRAD;
const double SEphem::Angle::sc_rotPerRad = ANGLE_ROTPERRAD;
const double SEphem::Angle::sc_hrsPerRad = ANGLE_HRSPERRAD;
const double SEphem::Angle::sc_radPerDeg = ANGLE_RADPERDEG;
const double SEphem::Angle::sc_radPerRot = ANGLE_RADPERROT;
const double SEphem::Angle::sc_radPerHrs = ANGLE_RADPERHRS;

static const unsigned divisor10[] =
{
    1, 10, 100, 1000, 10000, 100000,
    1000000, 10000000, 100000000
};

void SEphem::Angle::hms( unsigned& h, unsigned& m, unsigned& s, unsigned& f,
                         unsigned sec_digits ) const
{
    unsigned multiplier = divisor10[sec_digits];
    unsigned iangle = unsigned( floor( hrs() * 60 * 60 * multiplier + 0.5 ) );
    
    h = iangle / ( 60 * 60 * multiplier );
    m = ( iangle / ( 60 * multiplier ) ) % 60;
    s = ( iangle / multiplier ) % 60;
    f = iangle % multiplier;
}


void SEphem::Angle::dmsPM180( bool& negative, unsigned& d,
                              unsigned& m, unsigned& s, unsigned& f,
                              unsigned sec_digits ) const
{
    unsigned multiplier = divisor10[sec_digits];
    
    double df = degPM180();
    unsigned iangle = unsigned( floor( fabs( df ) * 60 * 60 * multiplier + 0.5 ) );
    
    negative = ( df < 0 );
    d = iangle / ( 60 * 60 * multiplier );
    m = ( iangle / ( 60 * multiplier ) ) % 60;
    s = ( iangle / multiplier ) % 60;
    f = iangle % multiplier;
}


void SEphem::Angle::dmsPM360( bool& negative, unsigned& d,
                              unsigned& m, unsigned& s, unsigned& f,
                              unsigned sec_digits ) const
{
    unsigned multiplier = divisor10[sec_digits];
    
    double df = degPM();
    unsigned iangle = unsigned( floor( fabs( df ) * 60 * 60 * multiplier + 0.5 ) );
    
    negative = ( df < 0 );
    d = iangle / ( 60 * 60 * multiplier );
    m = ( iangle / ( 60 * multiplier ) ) % 60;
    s = ( iangle / multiplier ) % 60;
    f = iangle % multiplier;
}


std::string SEphem::Angle::hmsString( unsigned sec_digits, bool hmsSep )
const
{
    unsigned hours;
    unsigned mins;
    unsigned secs;
    unsigned fsec;
    hms( hours, mins, secs, fsec, sec_digits );
    
    std::ostringstream stream;
    stream << std::setfill( '0' )
           << std::setw( 2 ) << std::setprecision( 2 ) << hours << ( hmsSep ? 'h' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << mins << ( hmsSep ? 'm' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << secs;
    if( sec_digits )stream << '.' << std::setw( sec_digits )
                               << std::setprecision( sec_digits ) << fsec;
    if( hmsSep )
    {
        stream << 's';
    }
    
    return stream.str();
}


std::string SEphem::Angle::hmsPMString( unsigned sec_digits, bool hmsSep )
const
{
    unsigned hours;
    unsigned mins;
    unsigned secs;
    unsigned fsec;
    bool negative = false;
    
    if( rad() <= Angle::sc_Pi )
    {
        hms( hours, mins, secs, fsec, sec_digits );
        negative = false;
    }
    else
    {
        Angle::makeRad( -rad() ).hms( hours, mins, secs, fsec, sec_digits );
        negative = true;
    }
    
    std::ostringstream stream;
    stream << ( negative ? '-' : '+' ) << std::setfill( '0' )
           << std::setw( 2 ) << std::setprecision( 2 ) << hours << ( hmsSep ? 'h' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << mins << ( hmsSep ? 'm' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << secs;
    if( sec_digits )stream << '.' << std::setw( sec_digits )
                               << std::setprecision( sec_digits ) << fsec;
    if( hmsSep )
    {
        stream << 's';
    }
    
    return stream.str();
}


std::string SEphem::Angle::dmsPM180String( unsigned sec_digits, bool dmsSep )
const
{
    bool negative;
    unsigned degs;
    unsigned mins;
    unsigned secs;
    unsigned fsec;
    dmsPM180( negative, degs, mins, secs, fsec, sec_digits );
    
    std::ostringstream stream;
    stream << ( negative ? '-' : '+' )
           << std::setfill( '0' )
           << std::setw( 2 ) << std::setprecision( 2 ) << degs << ( dmsSep ? 'd' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << mins << ( dmsSep ? 'm' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << secs;
    if( sec_digits )stream << '.' << std::setw( sec_digits )
                               << std::setprecision( sec_digits ) << fsec;
    if( dmsSep )
    {
        stream << 's';
    }
    
    return stream.str();
}


std::string SEphem::Angle::dmsPM360String( unsigned sec_digits, bool dmsSep )
const
{
    bool negative;
    unsigned degs;
    unsigned mins;
    unsigned secs;
    unsigned fsec;
    dmsPM360( negative, degs, mins, secs, fsec, sec_digits );
    
    std::ostringstream stream;
    stream << ( negative ? '-' : '+' )
           << std::setfill( '0' )
           << std::setw( 3 ) << std::setprecision( 3 ) << degs << ( dmsSep ? 'd' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << mins << ( dmsSep ? 'm' : ':' )
           << std::setw( 2 ) << std::setprecision( 2 ) << secs;
    if( sec_digits )stream << '.' << std::setw( sec_digits )
                               << std::setprecision( sec_digits ) << fsec;
    if( dmsSep )
    {
        stream << 's';
    }
    
    return stream.str();
}


std::string
SEphem::Angle::degString( unsigned dec_digits ) const
{
    std::ostringstream stream;
    int ideg = int( floor( deg() * divisor10[dec_digits] + 0.5 ) ) / divisor10[dec_digits];
    int fdeg = int( floor( deg() * divisor10[dec_digits] + 0.5 ) ) % divisor10[dec_digits];
    
    stream << '+' << ideg;
    if( dec_digits )
        stream << '.' << std::setw( dec_digits ) << std::setprecision( dec_digits )
               << std::setfill( '0' ) << fdeg;
    return stream.str();
}


std::string
SEphem::Angle::degPMString( unsigned dec_digits ) const
{
    std::ostringstream stream;
    double degPM = Angle::degPM();
    int ideg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) / divisor10[dec_digits];
    int fdeg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) % divisor10[dec_digits];
        
    stream << ( degPM >= 0 ? '+' : '-' ) << ideg;
    if( dec_digits )
        stream << '.' << std::setw( dec_digits ) << std::setprecision( dec_digits )
               << std::setfill( '0' ) << fdeg;
    return stream.str();
}


std::string
SEphem::Angle::degPM180String( unsigned dec_digits ) const
{
    std::ostringstream stream;
    double degPM = degPM180();
    int ideg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) / divisor10[dec_digits];
    int fdeg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) % divisor10[dec_digits];
        
    stream << ( degPM >= 0 ? '+' : '-' ) << ideg;
    if( dec_digits )
        stream << '.' << std::setw( dec_digits ) << std::setprecision( dec_digits )
               << std::setfill( '0' ) << fdeg;
    return stream.str();
}


std::string
SEphem::Angle::deg180String( unsigned dec_digits ) const
{
    std::ostringstream stream;
    double degPM = deg180();
    int ideg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) / divisor10[dec_digits];
    int fdeg =
        int( floor( fabs( degPM ) * divisor10[dec_digits] + 0.5 ) ) % divisor10[dec_digits];
        
    stream << ( degPM >= 0 ? '+' : '-' ) << ideg;
    if( dec_digits )
        stream << '.' << std::setw( dec_digits ) << std::setprecision( dec_digits )
               << std::setfill( '0' ) << fdeg;
    return stream.str();
}


bool SEphem::Angle::setFromHMSString( std::string str )
{
    unsigned hours = 0;
    unsigned mins = 0;
    unsigned secs = 0;
    unsigned fracs = 0;
    unsigned frac10s = 1;
    unsigned i = 0;
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            hours = hours * 10 + ( str[i++] - '0' );
        }
        else if( ( str[i] == ':' ) || ( str[i] == 'h' ) || ( str[i] == ' ' ) )
        {
            i++;
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            mins = mins * 10 + ( str[i++] - '0' );
        }
        else if( ( str[i] == ':' ) || ( str[i] == 'm' ) || ( str[i] == ' ' ) )
        {
            i++;
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            secs = secs * 10 + ( str[i++] - '0' );
        }
        else if( str[i] == '.' )
        {
            i++;
            break;
        }
        else if( str[i] == 's' )
        {
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            fracs = fracs * 10 + ( str[i++] - '0' );
            frac10s = frac10s * 10;
        }
        else if( str[i] == 's' )
        {
            break;
        }
        else
        {
            return false;
        }
    }
    
    setHrs( double( hours ) + double( mins ) / 60 + double( secs ) / ( 60 * 60 ) +
            double( fracs ) / double( frac10s ) / ( 60 * 60 ) );
            
    return true;
}


bool SEphem::Angle::setFromDMSString( std::string str )
{
    unsigned degs = 0;
    unsigned mins = 0;
    unsigned secs = 0;
    unsigned fracs = 0;
    unsigned frac10s = 1;
    unsigned i = 0;
    
    bool negative = false;
    if( str[i] == '-' )
    {
        negative = true;
        i++;
    }
    else if( str[i] == '+' )
    {
        negative = false;
        i++;
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            degs = degs * 10 + ( str[i++] - '0' );
        }
        else if( ( str[i] == ':' ) || ( str[i] == 'd' ) || ( str[i] == ' ' ) )
        {
            i++;
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            mins = mins * 10 + ( str[i++] - '0' );
        }
        else if( ( str[i] == ':' ) || ( str[i] == 'm' ) || ( str[i] == ' ' ) )
        {
            i++;
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            secs = secs * 10 + ( str[i++] - '0' );
        }
        else if( str[i] == '.' )
        {
            i++;
            break;
        }
        else if( str[i] == 's' )
        {
            break;
        }
        else
        {
            return false;
        }
    }
    
    while( i < str.length() )
    {
        if( ( str[i] >= '0' ) && ( str[i] <= '9' ) )
        {
            fracs = fracs * 10 + ( str[i++] - '0' );
            frac10s = frac10s * 10;
        }
        else if( str[i] == 's' )
        {
            i++;
            break;
        }
        else
        {
            return false;
        }
    }
    
    if( i < str.length() )
    {
        return false;
    }
    
    setDeg( ( negative ? -1 : 1 ) *
            ( double( degs ) + double( mins ) / 60 + double( secs ) / ( 60 * 60 ) +
              double( fracs ) / double( frac10s ) / ( 60 * 60 ) ) );
              
    return true;
}


unsigned long SEphem::Angle::toBAR( double a, unsigned radix, bool round )
{
    unsigned long maxbar = ( 0x00000001 << radix );
    return ( unsigned long )( floor( fmod( toRot( a ) * double( maxbar ) + ( round ? 0.5 : 0.0 ),
                                           maxbar ) ) );
}


double SEphem::Angle::frBAR( unsigned long a, unsigned radix, bool round )
{
    unsigned long maxbar = ( 0x00000001 << radix );
    return frRot( ( double( a % maxbar ) + ( round ? 0.5 : 0.0 ) ) / double( maxbar ) );
}
