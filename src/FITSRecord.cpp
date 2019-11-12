// -*- mode: c++; c-basic-offset: 4 -*-
//
// A class for writing FITS binary tables in a friendly manner
//
// This code is copied from   https://evlio.readthedocs.io
//
// \author Karl Kosack <karl.kosack@cea.fr>
//

#include <iostream>
#include <iomanip>
#include <fitsio.h>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <cstring>
#include <map>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>

#include "FITSRecord.h"

#define DEBUG 0

using namespace std;


/**
 * Construct a FITSRecord from an already open FITS filepointer
 * (must already be at the correct HDU)
 */
FITSRecord::FITSRecord( fitsfile* fptr )
    : _fptr( fptr ), _own_fptr( false ), _rowcount( 0 ), _verbose( 0 ),
      _extension_version( 0 )
{
    if( _fptr == NULL )
    {
        throw FITSRecordError( "FITS filepointer cannot be NULL!" );
    }
}



/**
 * Construct a FITSRecord, with automatic handling of opening the FITS file.
 *
 * \param fits_filename: name of input data
 * \param template_name: name of FITS .tpl file to use as model for the file
 * \param extension_name: which HDU to move to within the template
 * \param extension_version: version number of HDU, 0 for automatic
 */
FITSRecord::FITSRecord( std::string fits_filename,
                        std::string template_filename,
                        std::string extension_name, int extension_version )
    : 	 _fptr( NULL ), _own_fptr( true ), _rowcount( 0 ), _verbose( 0 ),
         //		 _extension_version( 0 ), _is_writable( true ), _is_finished( true )
         _extension_version( 0 ), _is_finished( true ), _is_writable( true )
{


    //string templatename = locateTemplate( template_filename );
    // this segfaults out when it tries to get environment variable FITSRECORD_PATH
    string templatename = template_filename ;
    
    if( DEBUG )
    {
        cout << "DEBUG: " << fits_filename << "  " << template_filename << endl;
    }
    
    int status = 0;
    string fullname = fits_filename + "(" + templatename + ")";
    
    // first try to open the file
    if( _verbose >= 2 )
    {
        cout << "Looking for: " << fits_filename << endl;
    }
    fits_open_file( &_fptr, fits_filename.c_str(), READWRITE, &status );
    
    if( status )
    {
        if( _verbose >= 2 )
        {
            cout << "REOPENING " << fits_filename << endl;
        }
        fitsfile* newfptr;
        fits_reopen_file( _fptr, &newfptr, &status );
    }
    
    // try to create the file
    if( status )
    {
    
        // the file may not exist, so try to create it:
        
        if( _verbose >= 2 )
        {
            cout << "CREATING " << fullname.c_str() << endl;
        }
        
        status = 0;
        fits_create_file( &_fptr, fullname.c_str(), &status );
        
    }
    
    // finally give up
    if( status )
    {
        throw FITSRecordError( "FITSRecord create file: " + getFITSError( status ) );
    }
    
    
    if( extension_name != "" )
    {
    
        if( DEBUG )
        {
            cout << "MOVING TO EXTENSION: " << extension_name << endl;
        }
        
        fits_movnam_hdu( _fptr, BINARY_TBL,
                         const_cast<char*>( extension_name.c_str() ),
                         extension_version, &status );
    }
    else
    {
    
        int type = 0;
        fits_movabs_hdu( _fptr, 1, &type, &status );
        
    }
    
    if( status )
    {
        throw FITSRecordError( "FITSRecord hdu: " + getFITSError( status ) );
    }
    
    // check the version number:
    int real_extension_version;
    char comment[100];
    fits_read_key( _fptr,  TINT, "EXTVER",  &real_extension_version,
                   comment, &status ) ;
    if( !status )
    {
        _extension_version = real_extension_version;
    }
    
    writeDate();
    
}

void
FITSRecord::writeDate()
{
    int status = 0;
    fits_write_date( _fptr, &status );
}


/**
 * Returns a structure containing all the information about the column
 * as read from the FITS file.
 */
FITSRecord::ColInfo
FITSRecord::
lookupColInfo( string column )
{

    FITSRecord::ColInfo col;
    int colnum = -1;
    int status = 0;
    
    char name[100];
    strncpy( name, column.c_str(), 100 ) ;
    
    fits_get_colnum( _fptr, CASEINSEN, name , &colnum, &status );
    
    if( status )
    {
        throw FITSRecordError( "Couldn't find column " + column
                               + " in the table (check the template or code for typos)" );
    }
    
    col.name = name;
    col.num = colnum;
    
    fits_get_coltype( _fptr, col.num,
                      &col.typecode, &col.repeat, &col.width,
                      &status );
                      
    return col;
}


/**
 * Helper method to add a column value to the map. Looks up the
 * column number, etc and checks the template definitions, warning
 * if the user-specified type doesn't match the template. Note
 * that this is ok since CFITSIO will automatically convert
 * between types, however truncation or loss of precision may
 * occur.
 */
void
FITSRecord::
addToMap( string column, Val* colval )
{

    // look up the column number of the column (and its datatype).
    
    ColInfo colinfo = lookupColInfo( column );
    
    if( _verbose >= 1 )
        cout << " adding Column: " << column
             << " with colnum " << colinfo.num
             << ", type=" << colinfo.typecode
             << ", and repeat=" << colinfo.repeat
             << endl;
             
             
    // check that the data types are the same (ignoring the case when
    // ULONG and LONG types are used, which are always stored as LONG
    // with an offset by FITS convention)
    if( _verbose >= 1 && ( colval->getType() != colinfo.typecode )
            && !( colval->getType() == TSTRING  && colinfo.typecode == TBYTE )
            && !( colval->getType() == TBYTE  && colinfo.typecode == TSTRING )
      )
    {
        cout << __FILE__ << ": WARNING: FITS column type for " << column
             << " (" << colval->getType() << ") "
             << "doesn't agree with template definition "
             << " (" << colinfo.typecode << ") "
             << "Type conversion may occur." << endl;
    }
    
    
    // check that the data typeshes are the same:
    if( !colval->isVariableLength() && colval->getSize() != colinfo.repeat )
    {
        cout <<  __FILE__ << ": WARNING: FITS column length for " << column
             << " (" << colval->getSize() << ") "
             << "doesn't agree with template definition "
             << " (" << colinfo.repeat << ") Updating column length..." << endl;
             
        int status = 0;
        fits_modify_vector_len(	_fptr, colinfo.num,
                                colval->getSize(),
                                &status );
        if( status )
            throw FITSRecordError( string( __func__ ) + ": modify vector length: "
                                   + getFITSError( status ) );
                                   
    }
    
    
    _colmap[colinfo.num] = colval;
    _namemap[colinfo.num] = column;
    
    
}



FITSRecord::~FITSRecord()
{

    // clean up the colmap
    
    map<int, Val*>::iterator it;
    for( it = _colmap.begin(); it != _colmap.end(); ++it )
    {
        if( it->second )
        {
            delete it->second;
        }
    }
    
    // if we own the fptr, close the file and finish
    
    if( _own_fptr )
    {
    
        int status = 0;
        
        if( _is_writable )
        {
            finishWriting();
        }
        
        fits_close_file( _fptr, &status );
        if( status )
            throw FITSRecordError( string( __func__ ) + " "
                                   + getFITSError( status ) );
                                   
        if( _verbose >= 1 )
        {
            cout << "File closed. _fptr=" << _fptr << endl;
        }
    }
    
    
}

/**
 * Call this when all rows have been written to update the NAXIS2
 * keyword in the header with the proper value (number of rows
 * written)
 */
void
FITSRecord::
finishWriting()
{

    int status = 0;
    
    fits_update_key( _fptr, TULONG, "NAXIS2",
                     &_rowcount, ( char* )"number of events", &status );
    if( status )
        cout << "ERROR: KEYWORD UPDATE FAILED FOR NAXIS2."
             << " File may be corrupt." << endl;
             
    // write the Checksum too:
    
    fits_write_chksum( _fptr, &status );
    if( status )
        cout << "ERROR: CHECKSUM couldn't be written."
             << " File may be corrupt." << endl;
             
}

/**
 * Helper function to print the FITSRecord to a stream. You don't need
 * to call this directly, rather just use the << operator.
 *
 * \code
 * FITSRecord rec;
 * ...
 * cout << rec << endl;
 * \endcode
 */
ostream&
FITSRecord::
print( ostream& stream )
{

    stream << setw( 6 ) << _rowcount << ": ";
    
    map<int, Val*>::iterator it;
    
    for( it = _colmap.begin(); it != _colmap.end(); ++it )
    {
    
        int colnum = it->first;
        Val* colval = it->second;
        
        stream << _namemap[colnum] << "="
               << *colval << "; ";
               
    }
    
    return stream;
    
}


/**
 * Try to find the specified template file.  Looks in FITSRECORD_PATH
 * environment variable for a list of possibilities if the file isn't
 * found in the local directory.
 */
string
FITSRecord::
locateTemplate( string template_filename )
{


    vector<string> directories;
    string delimiters = ":";
    string str = "";
    
    try
    {
        str = getenv( "FITSRECORD_PATH" );
    }
    catch( ... )
    {
    
    }
    
    directories.push_back( "." );
    
    string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    string::size_type pos     = str.find_first_of( delimiters, lastPos );
    
    while( string::npos != pos || string::npos != lastPos )
    {
        directories.push_back( str.substr( lastPos, pos - lastPos ) );
        lastPos = str.find_first_not_of( delimiters, pos );
        pos = str.find_first_of( delimiters, lastPos );
    }
    
    // try each directory in the path
    for( vector<string>::iterator dir = directories.begin();
            dir != directories.end(); dir++ )
    {
    
        string templatefile = *dir + "/" + template_filename;
        
        if( DEBUG )
            cout << "DEBUG: searching " << *dir
                 << " for " << template_filename << endl;
                 
        // try to open the specified file:
        ifstream infile( templatefile.c_str() );
        if( infile )
        {
            infile.close();
            return templatefile;
        }
        
    }
    
    throw FITSRecordError( "Template '" + template_filename
                           + "' not found in FITSRECORD_PATH" );
                           
}

string
FITSRecord::
getFITSError( int status )
{

    string strerr = "";
    
    if( status )
    {
        char text[80];
        fits_get_errstatus( status, text );
        strerr = text;
        
        // get any other error messages on the stack:
        while( fits_read_errmsg( text ) )
        {
            strerr += "\n" + string( text );
        }
        
    }
    
    return " [" + strerr + "] ";
    
}



/**
 * Write one row of the table using the values of all currently mapped
 * variables.  Any non-mapped variable that is defined in the FITS
 * template will get written as void.
 *
 * \throws a FITSRecordError if the writing fails
 */
void
FITSRecord::
write() throw( FITSRecordError )
{

    if( !_is_writable )
    {
        throw FITSRecordError( "can't call write() on a read-only FITSRecord" );
    }
    
    int status = 0;
    
    map<int, Val*>::iterator it;
    
    for( it = _colmap.begin(); it != _colmap.end(); ++it )
    {
    
        int colnum = it->first;
        Val* colval = it->second;
        
        fits_write_col( _fptr,
                        colval->getType(),
                        colnum,
                        _rowcount + 1,
                        1,	// first element
                        colval->getSize(), // n elements
                        colval->getVoidPointerToValue() ,
                        &status );
        if( status )
            throw FITSRecordError( "write failed for column "
                                   + _namemap[colnum] + getFITSError( status ) );
                                   
                                   
        if( DEBUG )
        {
            long nrows;
            fits_get_num_rows( _fptr, &nrows, &status );
            cout << "DEBUG: nrows=" << nrows << endl;
        }
        
    }
    
    _rowcount++;
    
}

/**
 * If you know how many events will be written, you can call this
 * function first to allocate space for all the rows (speeding up the
 * writing). Otherwise, rows are added dynamically
 */
void
FITSRecord::
allocateRows( size_t num_rows )
{

    int status = 0;
    if( _verbose >= 1 )
    {
        cout << "Reserving " << num_rows << " rows..." << endl;
    }
    fits_insert_rows( _fptr, _rowcount, num_rows, &status ) ;
    
    if( status )
    {
        throw FITSRecordError( "allocateRows: " + getFITSError( status ) );
    }
    
    
    
}


/**
 * Delete rows from end of table starting with start_row (in case you
 * allocated too many with allocateRows())
 */
void
FITSRecord::
truncateRows( size_t start_row )
{

    int status = 0;
    
    if( start_row >= _rowcount )
    {
        cout << "truncateRows: there are no rows to truncate after " << start_row << endl;
        return;
    }
    
    if( _verbose >= 1 )
    {
        cout << "Deleting the last " << _rowcount - start_row << " rows..." << endl;
    }
    
    fits_delete_rows( _fptr, start_row, _rowcount - start_row, &status );
    
    if( status )
    {
        throw FITSRecordError( "truncateRows: " + getFITSError( status ) );
    }
    
    
}


/**
 * Returns true if the column name is found in the file or
 * template. This is useful for checking what template is being used
 * (e.g. does it have simulation parameters?)
 */
bool
FITSRecord::
hasColumn( std::string colname )
{

    int status = 0;
    int colnum;
    
    char name[100]; // needed since doesn't work properly with
    // c_str(), which returns a const char*
    strncpy( name, colname.c_str(), 100 );
    
    fits_get_colnum( _fptr, CASEINSEN, name, &colnum, &status );
    
    if( status == 0 )
    {
        return true;
    }
    
    return false;
}




bool FITSRecord::printFITSTypes()
{

    double d = 0;
    cout << "double -> " << getFITSType( d ) << endl;
    
    float f = 0;
    cout << "float  -> " << getFITSType( f ) << endl;
    
    int i = 0;
    cout << "int    -> " << getFITSType( i ) << endl;
    
    unsigned long ul = 0;
    cout << "ulong  -> " << getFITSType( ul ) << endl;
    
    long l = 0;
    cout << "long    -> " << getFITSType( l ) << endl;
    
    long long ll = 0;
    cout << "longlong-> " << getFITSType( ll ) << endl;
    
    bool b = 0;
    cout << "bool    -> " << getFITSType( b ) << endl;
    
    char c = 0;
    cout << "char    -> " << getFITSType( c ) << endl;
    
    unsigned char uc = 0;
    cout << "uchar   -> " << getFITSType( uc ) << endl;
    
    
    unsigned short us = 0;
    cout << "ushort  -> " << getFITSType( us ) << endl;
    
    string s;
    cout << "string  -> " << getFITSType( s ) << endl;
    
    return true;
    
}


FITSRecord::FITSRecord( std::string url )
    : _fptr( NULL ), _own_fptr( true ), _rowcount( 0 ), _verbose( 0 ),
      _extension_version( 0 ),  _is_finished( true ), _is_writable( false )
{

    int status = 0;
    
    if( _verbose >= 1 )
    {
        cout << "OPENING URL: " << url << endl;
    }
    fits_open_file( &_fptr, url.c_str(), READONLY, &status );
    
    if( status )
    {
        throw FITSRecordError( "FITSRecord open file  " + url
                               + ": " + getFITSError( status ) );
    }
    
    // check the version number:
    // int real_extension_version;
    // char comment[100];
    // fits_read_key( _fptr,  TINT, "EXTVER",  &real_extension_version,
    // 		   comment, &status) ;
    // if (!status)
    // 	_extension_version = real_extension_version;
    
    if( _verbose > 1 )
    {
        cout << "EXTENSION VERSION: " <<  _extension_version << endl;
    }
    
    _is_finished = false; // now ready to read
    
}


/**
 * Reads a single row of the input table into the data structure
 * defined by mapColumnToVar().
 *
 * \returns true if a data row was read, false if no data is available
 * (or file is finished)
 */
bool
FITSRecord::
read()
{

    int status = 0;
    int nulval = 0;
    int anynull = 0;
    
    if( _is_finished )
    {
        return false;
    }
    
    map<int, Val*>::iterator it;
    for( it = _colmap.begin(); it != _colmap.end(); ++it )
    {
    
        int colnum = it->first;
        Val* colval = it->second;
        
        fits_read_col( _fptr,
                       colval->getType(),
                       colnum,
                       _rowcount + 1, // first row
                       1, // first elelement
                       colval->getSize(), //num elements
                       &nulval,
                       colval->getVoidPointerToValue(),
                       &anynull,
                       &status );
                       
        if( status )
        {
            _is_finished = true;
            if( status != BAD_ROW_NUM ) // normal end of file
            {
                throw FITSRecordError( "read:" + getFITSError( status ) );
            }
            break;
        }
        
    }
    
    _rowcount++;
    return true;
    
}


/**
 * Allow printing of FITSRecords
 */
std::ostream& operator<<( std::ostream& stream, FITSRecord& rec )
{
    return rec.print( stream );
}

std::ostream& operator<<( std::ostream& stream, Val& colval )
{
    return colval.stream( stream );
}
