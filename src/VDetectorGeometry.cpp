/*! \class VDetectorGeometry
    \brief description of telescopes, camera and FADC system


*/

#include "VDetectorGeometry.h"

VDetectorGeometry::VDetectorGeometry( unsigned int iNTel, vector< string > iCamera, string iDir,
                                      bool iDebug, float iCoordinateTransformerX, float iCoordinateTransformerY,
                                      int iSourceType )
{
    fDebug = iDebug;
    if( fDebug )
    {
        cout << "VDetectorGeometry::VDetectorGeometry" << endl;
    }
    if( iNTel > iCamera.size() )
    {
        cout << "VDetectorGeometry::VDetectorGeometry  error: number of telescopes larger than camera vector ";
        cout << iNTel << "\t" << iCamera.size() << endl;
        exit( EXIT_FAILURE );
    }
    setCoordinateTransformer( iCoordinateTransformerX, iCoordinateTransformerY );
    setSourceType( iSourceType );
    
    // set directory with all configuration files
    setConfigDir( iDir );
    
    // detector configuration file from GrIsu (.cfg)
    if( iCamera[0].find( ".cfg" ) < iCamera[0].size() || iCamera[0].find( ".txt" ) < iCamera[0].size() )
    {
        if( fDebug )
        {
            cout << "VDetectorGeometry::VDetectorGeometry: .cfg file detected" << endl;
        }
        readGrisucfg( iCamera[0], iNTel );
    }
    // camera configuration from .cam file (default telescope positions)
    else
    {
        if( !initialize( iNTel, iCamera ) )
        {
            cout << "DetectorGeometry::VDetectorGeometry error in camera reader" << endl;
            exit( EXIT_FAILURE );
        }
        for( unsigned int i = 0; i < iNTel; i++ )
        {
            setTelID( i );
            readCameraFile( iCamera[i] );
        }
    }
    
    for( unsigned int i = 0; i < iNTel; i++ )
    {
        // set channels and samples
        setTelID( i );
        fNChannels.push_back( getNumChannels() );
        fNSamples.push_back( getNumSamples() );
        fSampleWarning.push_back( true );
    }
    
    //    if( fDebug ) print();
    if( fDebug )
    {
        cout << "END: VDetectorGeometry::VDetectorGeometry" << endl;
    }
}


VDetectorGeometry::VDetectorGeometry( unsigned int iNTel, bool iDebug )
{
    fDebug = iDebug;
    if( fDebug )
    {
        cout << "VDetectorGeometry::VDetectorGeometry (ntel=" << iNTel << ")" << endl;
    }
}


void VDetectorGeometry::setNSamples( unsigned int iTelID, unsigned int iNSamples, bool iForceSet )
{
    if( iTelID < fNSamples.size() )
    {
        // ignore if samples length from cfg is shorter than request sample length (except if iForceSet is set)
        if( fNSamples[iTelID] < iNSamples && !iForceSet )
        {
            if( fSampleWarning[iTelID] )
            {
                fSampleWarning[iTelID] = false;
            }
        }
        else
        {
            fNSamples[iTelID] = iNSamples;
        }
    }
}


void VDetectorGeometry::addDataVector( unsigned int iNTel, vector< unsigned int > iNChannels )
{
    for( unsigned int i = 0; i < iNChannels.size(); i++ )
    {
        fNChannels.push_back( iNChannels[i] );
        fNSamples.push_back( 0 );
        fSampleWarning.push_back( true );
    }
    initialize( iNTel, iNChannels );
}

void VDetectorGeometry::setNChannels( unsigned int iTelID, unsigned int iNChannels )
{
    if( iTelID < fNChannels.size() )
    {
        fNChannels[iTelID] = iNChannels;
    }
}
