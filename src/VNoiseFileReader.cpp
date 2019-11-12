/*!  \class VNoiseFileReader
     \brief read noise traces from a file



*/

#include "VNoiseFileReader.h"

VNoiseFileReader::VNoiseFileReader( unsigned int iType, string iFileName )
{
    fNoiseFileType = iType;
    fNoiseFileName = iFileName;
    
    fDebug = false;
    fZombie = false;
    
    fGrIsuReader = 0;
    
    if( fNoiseFileType == 0 )
    {
        cout << "VNoiseFileReader: noise file type is grisu style" << endl;
    }
    else
    {
        fZombie = true;
        cout << "VNoiseFileReader: unknown noise file type: " << fNoiseFileType << endl;
    }
}


bool VNoiseFileReader::init( VDetectorGeometry* iD, unsigned int intel, vector<int> iSW, bool iDebug, int iseed, double iFADCorrect )
{
    if( iD == 0 )
    {
        fZombie = true;
    }
    if( iSW.size() == 0 )
    {
        fZombie = true;
    }
    if( fZombie )
    {
        return false;
    }
    
    // initialize grisu reader
    if( fNoiseFileType == 0 )
    {
        fGrIsuReader = new VGrIsuReader( iD, intel,  fNoiseFileName, iSW, iDebug, iseed, iFADCorrect );
    }
    
    return true;
}


vector< uint8_t >& VNoiseFileReader::getNoiseVec( unsigned int iTel, uint32_t iHitID,  bool iNewTrace )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getNoiseVec( iTel, iHitID, iNewTrace );
    }
    
    return vv8;
}


uint8_t VNoiseFileReader::getNoiseSample( unsigned int iTel, uint32_t iHitID, unsigned int iSample, bool iNewTrace )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getNoiseSample( iTel, iHitID, iSample, iNewTrace );
    }
    
    return 0;
}


void VNoiseFileReader::setDefaultGrisuPed( double iB )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        cout << "\t VNoiseFileReader::setDefaultGrisuPed to " << iB << endl;
        fGrIsuReader->setDefaultPed( iB );
    }
}


void VNoiseFileReader::setSumWindow( unsigned int iTelID, int isw )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        fGrIsuReader->setSumWindow( iTelID, isw );
    }
}


valarray<double>& VNoiseFileReader::getPeds()
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getPeds();
    }
    
    return v;
}


valarray<double>& VNoiseFileReader::getPedvars()
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getPedvars();
    }
    
    return v;
}


vector< valarray<double> >& VNoiseFileReader::getPedvarsAllSumWindows()
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getPedvarsAllSumWindows();
    }
    
    return vv;
}


valarray<double>& VNoiseFileReader::getPedRMS()
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getPedRMS();
    }
    
    return v;
}


vector< vector< vector< uint8_t > > > VNoiseFileReader::getFullNoiseVec()
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getFullNoiseVec();
    }
    
    vector< vector< vector< uint8_t > > > a;
    
    return a;
}


vector< vector< uint8_t > >& VNoiseFileReader::getFullNoiseVec( unsigned int iTel )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getFullNoiseVec( iTel );
    }
    
    return v8;
}


vector< uint8_t >& VNoiseFileReader::getFullNoiseVec( unsigned int iTel, int iChannel )
{
    if( fNoiseFileType == 0 && fGrIsuReader )
    {
        return fGrIsuReader->getFullNoiseVec( iTel, iChannel );
    }
    
    return vv8;
}
