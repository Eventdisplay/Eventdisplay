/*
 * exception.h -- base-class for all VERITAS exceptions
 * by Filip Pizlo, 2002, 2003
 */

#ifndef V_EXCEPTION_H
#define V_EXCEPTION_H

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <exception>
#include <typeinfo>
#include <sstream>
#include <vector>
#include <new>

class VException:
    public std::exception
{
    private:
        std::string desc, msg;
        
        mutable std::string total;
        
        std::string throw_file;
        unsigned throw_line;
        
        std::vector< std::string > comments;
        
        std::ostringstream last_comment;
        
        void makeTotal() const
        {
            std::ostringstream ret;
            
            ret << desc.c_str();
            
            if( msg.length() != 0 )
            {
                ret << ": " << msg.c_str();
            }
            
            if( throw_file.length() != 0 &&
                    throw_line != 0 )
            {
                ret << ": thrown in " << throw_file << " on line " << throw_line;
            }
            
            for( std::vector< std::string >::const_iterator
                    i = comments.begin();
                    i != comments.end();
                    ++i )
            {
                ret << ": " << *i;
            }
            
            std::string last_cmnt_str = last_comment.str();
            if( last_cmnt_str.length() > 0 )
            {
                ret << ": " << last_cmnt_str;
            }
            
            total = ret.str();
        }
        
    protected:
    
        VException( const std::string& exception_type,
                    const std::string& thrown_by ):
            desc( exception_type )
        {
            msg.append( "In " );
            msg.append( thrown_by );
        }
        
        void setType( std::string _desc )
        {
            desc = _desc;
            makeTotal();
        }
        
        void setType( const char* _desc )
        {
            if( _desc != NULL )
            {
                setType( std::string( _desc ) );
            }
        }
        
        void setMessage( std::string _msg )
        {
            msg = _msg;
            makeTotal();
        }
        
        void setMessage( const char* _msg )
        {
            if( _msg != NULL )
            {
                setMessage( std::string( _msg ) );
            }
        }
        
        void setStrings( const char* msg, const char* arg )
        {
            setType( msg );
            if( arg != NULL )
            {
                setMessage( arg );
            }
        }
        
    public:
    
        VException()
        {
        }
        
        VException( const VException& other ):
            std::exception(),
            desc( other.desc ),
            msg( other.msg ),
            throw_file( other.throw_file ),
            throw_line( other.throw_line ),
            comments( other.comments )
        {
            last_comment.str( other.last_comment.str() );
            endComment();
        }
        
        virtual ~VException()
        throw()
        {
        }
        
        virtual const char* what() const
        throw()
        {
            makeTotal();
            return total.c_str();
        }
        
        std::ostream& stream() throw()
        {
            return last_comment;
        }
        
        void endComment() throw()
        {
            std::string last_cmnt_str = last_comment.str();
            last_comment.str( "" );
            if( last_cmnt_str.length() > 0 )
            {
                comments.push_back( last_cmnt_str );
            }
        }
        
        void addComment( const char* comment )
        {
            endComment();
            
            comments.push_back( std::string( comment ) );
        }
        
        void setThrowLocation( const char* filename,
                               unsigned line_no )
        {
            throw_file = filename;
            throw_line = line_no;
        }
        
        const char* getType() const
        throw()
        {
            return desc.c_str();
        }
        
        const char* getMessage() const
        throw()
        {
            return msg.c_str();
        }
        
        const char* getThrowFile() const
        throw()
        {
            return throw_file.c_str();
        }
        
        unsigned getThrowLine() const
        throw()
        {
            return throw_line;
        }
        
        const std::vector< std::string > getComments() const
        {
            std::vector< std::string > ret( comments );
            std::string last_cmnt_str = last_comment.str();
            if( last_cmnt_str.length() > 0 )
            {
                ret.push_back( last_cmnt_str );
            }
            return ret;
        }
        
        virtual void printTypeSpecificExceptionDetails( std::ostream& output )
        const
        {
        }
        
};

inline std::ostream& operator<<( std::ostream& stream, VException& x )
{
    stream << "Exception: " << x.getType() << ": " << x.getMessage() << std::endl;
    stream << "Code location: " << x.getThrowFile() << ":" << x.getThrowLine() << std::endl;
    x.printTypeSpecificExceptionDetails( stream );
    std::vector< std::string > comments = x.getComments();
    for( std::vector< std::string >::iterator
            i = comments.begin();
            i != comments.end();
            ++i )
    {
        stream << *i << std::endl;
    }
    
    return stream;
}


class VSystemException: public VException
{
    private:
    
        int sys_errno;
        
        void composeError( const char* msg )
        {
            std::ostringstream buf;
            buf << "System error: " << strerror( sys_errno );
            setType( buf.str().c_str() );
            setMessage( msg );
        }
        
    public:
    
        VSystemException( const char* msg = NULL ):
            sys_errno( errno )
        {
            composeError( msg );
        }
        
        VSystemException( int _errno, const char* msg = NULL ):
            sys_errno( _errno )
        {
            composeError( msg );
        }
        
        int getErrno()
        {
            return sys_errno;
        }
        
};

class VInvalidArgumentException: public VException
{
    public:
        VInvalidArgumentException( const char* msg = NULL )
        {
            setStrings( "Invalid Argument", msg );
        }
};

class VAssertionFailedException: public VException
{
    public:
        VAssertionFailedException( const char* exp,
                                   const char* file,
                                   unsigned line )
        {
            setType( "Assertion failed" );
            std::ostringstream buf;
            buf << "For \"" << exp << "\"";
            setMessage( buf.str().c_str() );
            setThrowLocation( file, line );
        }
};

// helpful stuff

#define V_THROW(e_class,msg) { \
        e_class __e_to_throw(msg); \
        __e_to_throw.setThrowLocation(__FILE__,__LINE__); \
        throw __e_to_throw; \
    }

#define V_ASSERT(exp) ((void)((exp)?0: \
                              (throw VAssertionFailedException(#exp,__FILE__,__LINE__), 0)))
#endif
