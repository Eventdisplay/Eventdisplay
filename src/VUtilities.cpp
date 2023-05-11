/*! \file VUtilities.cpp

 */

#include "VUtilities.h"

//! Make a lowercase copy of s:
string VUtilities::lowerCase( string& s )
{
    char* buf = new char[s.length()];
    s.copy( buf, s.length() );
    for( unsigned int i = 0; i < s.length(); i++ )
    {
        buf[i] = tolower( buf[i] );
    }
    string r( buf, s.length() );
    delete [] buf;
    return r;
}


//! Make an uppercase copy of s:
string VUtilities::upperCase( string& s )
{
    char* buf = new char[s.length()];
    s.copy( buf, s.length() );
    for( unsigned int i = 0; i < s.length(); i++ )
    {
        buf[i] = toupper( buf[i] );
    }
    string r( buf, s.length() );
    delete [] buf;
    return r;
}


/*!
    \brief delete pointers in an STL sequence container
*/
template<class Seq> void VUtilities::purge( Seq& c )
{
    typename Seq::iterator i;
    for( i = c.begin(); i != c.end(); i++ )
    {
        delete *i;
        *i = 0;
    }
}

/*

    test if a certain file exist

*/
string VUtilities::testFileLocation( string iFile, string iDirectory, bool bEVNDISPDATA, bool bExitIfFailure )
{
    TString iPath( iFile.c_str() );
    gSystem->ExpandPathName( iPath );
    ifstream is;
    is.open( iPath.Data(), ifstream::in );
    if( !is )
    {
        string itemp;
        if( bEVNDISPDATA )
        {
            itemp = VGlobalRunParameter::getDirectory_EVNDISPAnaData();
            itemp += "/" + iDirectory + "/" + iFile;
        }
        else
        {
            itemp = "./" + iDirectory + "/" + iFile;
        }
        ifstream is2;
        is2.open( itemp.c_str(), ifstream::in );
        if( !is2 )
        {
            cout << "testFileLocation: Error opening file: " << iFile << endl;
            cout << "(not found in current directory and in " << iDirectory << ")" << endl;
            iFile = "";
            if( bExitIfFailure )
            {
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            return iFile;
        }
        iFile = itemp;
    }
    return iFile;
}

/*

    remove all white spaces from a string

*/
string VUtilities::removeSpaces( string stringIn )
{
    string::size_type pos = 0;
    bool spacesLeft = true;
    
    while( spacesLeft )
    {
        pos = stringIn.find( " " );
        if( pos != string::npos )
        {
            stringIn.erase( pos, 1 );
        }
        else
        {
            spacesLeft = false;
        }
    }
    
    return stringIn;
}

/*

   remove all leading spaces in a string

*/
string VUtilities::remove_leading_spaces( string stringIn )
{
    string::size_type pos = 0;
    bool spacesLeft = true;
    
    while( spacesLeft )
    {
        pos = stringIn.find( " " );
        if( pos != 0 )
        {
            pos = stringIn.find( "\t" );
        }
        
        if( pos == 0 && pos != string::npos )
        {
            stringIn.erase( pos, 1 );
        }
        else
        {
            spacesLeft = false;
        }
    }
    
    return stringIn;
}

/*

   remove leading and trailing spaces

*/
string VUtilities::trim_spaces( string str, string whitespace )
{
    const size_t strBegin = str.find_first_not_of( whitespace );
    if( strBegin == std::string::npos )
    {
        return "";    // no content
    }
    
    const size_t strEnd = str.find_last_not_of( whitespace );
    const size_t strRange = strEnd - strBegin + 1;
    
    return str.substr( strBegin, strRange );
}


/*

   search and replace a certain letter in a string

*/
string VUtilities::search_and_replace( string i1, string iO, string iN )
{
    size_t j;
    for( ; ( j = i1.find( iO ) ) != string::npos ; )
    {
        i1.replace( j, iO.length(), iN );
    }
    return i1;
}


/* ------------------- line_point_distance --------------------- */
/**
 *  Distance between a straight line and a point in space
 *
 *  @param  x1,y1,z1  reference point on the line
 *  @param  cx,cy,cz  direction cosines of the line
 *  @param  x,y,z     point in space
 *
 *  @return distance
 *
 *
 * OBSERVE:  ground coordinates should be in CORSIKA coordinate system
 *
 *
 */
double VUtilities::line_point_distance( double x1, double y1, double z1, double ze, double az, double x, double y, double z )
{
    double alt = 90. - ze;
    az = 180. - az;
    
    double cx = -1.*cos( alt * TMath::DegToRad() ) * cos( az * TMath::DegToRad() );
    double cy = -1.*cos( alt * TMath::DegToRad() ) * sin( az * TMath::DegToRad() );
    double cz = sin( alt * TMath::DegToRad() );
    
    double a1 = ( y - y1 ) * cz - ( z - z1 ) * cy;
    double a2 = ( z - z1 ) * cx - ( x - x1 ) * cz;
    double a3 = ( x - x1 ) * cy - ( y - y1 ) * cx;
    double a  = a1 * a1 + a2 * a2 + a3 * a3;
    double b = cx * cx + cy * cy + cz * cz;
    
    if( a < 0. || b <= 0. )
    {
        return -1;
    }
    
    return sqrt( a / b );
}

unsigned int VUtilities::count_number_of_textblocks( string str )
{
    str = trim_spaces( str );
    
    unsigned int z = 0;
    string iTemp;
    istringstream is_stream( str );
    while( !( is_stream >> std::ws ).eof() )
    {
        is_stream >> iTemp;
        z++;
    }
    
    return z;
}
