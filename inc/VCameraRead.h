//! VCameraRead reading of camera geometry file (.cam and GrIsu .cfg files)
#ifndef VCAMERAREAD_H
#define VCAMERAREAD_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TMath.h"
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>

#include "VGlobalRunParameter.h"
#include "VDB_Connection.h"

using namespace std;

class VCameraRead : public VGlobalRunParameter
{
    private:
    
    protected:
        bool fDebug;
        unsigned int fGrIsuVersion;               //!< GrIsu Version
        unsigned int fCFGtype;                    //!< cfg file type (0=std, 1 = mirrors and pixels for 1 telelescope only)
        unsigned int fTelID;                      //!< telescope ID
        
        int          fsourcetype;
        //!< telescope ID for multiple data readers
        map< unsigned int, unsigned int > fTelIDGrisu;
        unsigned int fNTel;                       //!< number of telescopes
        string fConfigDir;                        //!< directory with geometry config files
        // telescope type
        vector< ULong64_t > fTelType;
        // telescope positions
        vector< float > fTelXpos;                 //!< telescope X position in [m]
        vector< float > fTelYpos;                 //!< telescope Y position in [m]
        vector< float > fTelZpos;                 //!< telescope Z position in [m]
        vector< float > fTelRad;                  //!< telescope mirror radius in [m]
        // mirror design
        vector< float > fMirFocalLength;          //!< mirror focal length in [m]
        vector< unsigned int > fNMirrors;         //!< number of mirros
        vector< float > fMirrorArea;              //!< mirror area
        // camera data
        vector< string > fCameraName;             //!< camera name
        vector< float > fCameraScaleFactor;       //!< camera stretch factor
        vector< float > fCameraCentreOffset;      //!< camera centre offset [deg]
        vector< float > fCameraRotation;          //!< camera rotation (clockwise)
        vector< float > fCameraFieldofView;       //!< camera field of view [deg]
        unsigned int fPixelType;                  //!< pixel type: 1 = circle, 2 = hex, 3 = square
        vector< unsigned int > fCNChannels;       //!< number of channels
        vector< unsigned int > fCNSamples;        //!< number of FADC samples
        vector< float > fSample_time_slice;       //!< length of time slice
        // 2D camera data
        vector< vector<float> > fXTube;           //!< x-position of tube in [deg] (in camera coordinates)
        vector< vector<float> > fYTube;           //!< y-position of tube in [deg] (in camera coordiantes)
        vector< vector<float> > fRTube;           //!< tube radius in [deg]
        vector< vector<float> > fRotXTube;        //!< x-position of tube in [deg] (rotated, x-axis is now horizontal for telescope at elevation 90)
        vector< vector<float> > fRotYTube;        //!< y-position of tube in [deg] (rotated, x-axis is now horizontal for telescope at elevation 90)
        vector< vector<float> > fXTubeMM;         //!< x-position of tube in [mm]
        vector< vector<float> > fYTubeMM;         //!< y-position of tube in [mm]
        vector< vector<float> > fRTubeMM;         //!< tube radius in [mm]
        vector< vector<int> > fTrigTube;          //!< pixel takes part in trigger
        vector< vector<int> > fAnaTube;           //!< pixel takes part in analysis
        vector< vector< bool > > fEdgePixel;      //!< this pixel is a edge pixel
        vector< vector< unsigned int > > fNNeighbour;  //!< number of neighbours
        vector< vector< vector<int> > > fNeighbour; //!< neighbour identifier
        vector< unsigned int > fMaxNeighbour;               //!< maximal number of neighbours
        vector< unsigned int > fCameraCentreTubeIndex;  //!< vector index of centre tube (vector length: [ntel]
        // pattern trigger
        int fNPatches;                            //!< number of trigger patches
        vector< vector<int> > fPatch;             //!< pattern trigger patches
        
        vector<unsigned int> fMix;                //!< pixel number in real data .cam file
        vector<unsigned int> fXim;                //!< pixel number in real data .cam file
        // calibration data
        double fDefPed;                           //!< default pedestal
        int fFADCRange;                           //!< FADC range
        vector< vector<float> > fGain;            //!< gain
        vector< vector<float> > fTOff;            //!< toffsets
        // electronics
        bool                  fLowGainIsSet;       //!< low gain multiplier is set in cfg file
        vector< double >      fLowGainMultiplier_Trace;  //!< low gain multiplier (usually 6.)
        vector< unsigned int> fLowGainActivator;  //!< threshold for low gain activation (usually 255)
        
        void                 cleanNeighbourList();
        void                 convertMMtoDeg();    //!< convert camera vectors from mm to deg
        void                 fillTelescopeVectors();
        void                 readPixelFile( string );
        void                 resetCamVectors( bool bMaxN = true );
        void                 resetNeighbourLists( bool bMaxN = true );
        void                 resetTelVectors();
        
        // coordinate transformer
        float    fCoordinateTransformerX;
        float    fCoordinateTransformerY;
        
    public:
        VCameraRead();
        ~VCameraRead() {};
        vector<int>&         getAnaPixel()        //!< get vector with analysis pixels (not switched off by MC/user)
        {
            return fAnaTube[fTelID];
        }
        //!< get vector with analysis pixels (not switched off by MC/user)
        vector<int>&         getAnaPixel( unsigned int iTel )
        {
            return fAnaTube[iTel];
        }
        unsigned int getCameraCentreTubeIndex()
        {
            if( fTelID < fCameraCentreTubeIndex.size() )
            {
                return fCameraCentreTubeIndex[fTelID];
            }
            else
            {
                return 9999;
            }
        }
        string               getCameraName()      //!< get camera name
        {
            return fCameraName[fTelID];
        }
        double               getDefaultPedestal() //!< get default pedestal
        {
            return fDefPed;
        }
        int                  getFADCRange()
        {
            return fFADCRange;
        }
        vector<float>&       getFocalLength()     //! get focal length [m]
        {
            return fMirFocalLength;
        }
        vector< unsigned int>& getNMirrors()
        {
            return fNMirrors;
        }
        vector< float >&     getMirrorArea()
        {
            return fMirrorArea;
        }
        vector<float>&       getFieldofView()
        {
            return fCameraFieldofView;
        }
        unsigned int         getGrIsuVersion()
        {
            return fGrIsuVersion;
        }
        bool                 isLowGainSet()
        {
            return fLowGainIsSet;
        }
        vector<float>&       getCameraScaleFactor()
        {
            return fCameraScaleFactor;
        }
        vector<float>&       getCameraCentreOffset()
        {
            return fCameraCentreOffset;
        }
        vector<float>&       getCameraRotation()
        {
            return fCameraRotation;
        }
        vector<unsigned int> getLowGainThreshold()
        {
            return fLowGainActivator;
        }
        vector< double >     getLowGainMultiplier_Trace()
        {
            return fLowGainMultiplier_Trace;
        }
        vector<unsigned int>& getMCPixel()        //!< get pixel number v[real]=mc (pixel numbering different in MC and reality)
        {
            return fMix;
        }
        unsigned int         getMaxNeighbour()
        {
            return fMaxNeighbour[fTelID];
        }
        float                getMaximumFOV_deg();
        vector<vector<int> >& getNeighbours()     //!< neighbour identifier
        {
            return fNeighbour[fTelID];
        }
        vector<unsigned int>& getNNeighbours()    //!< number of neighbours
        {
            return fNNeighbour[fTelID];
        }
        int                  getNPatches()        //!< return number of patches
        {
            return fNPatches;
        }
        vector<vector<int> > getPatch()           //!< return patch
        {
            return fPatch;
        }
        vector< unsigned int > getNumChannelVector();
        unsigned int         getNumChannels()     //!< get number of channels in this camera
        {
            return fXTube[fTelID].size();
        }
        unsigned int         getNumSamples()      //!< get number of FADC samples
        {
            return fCNSamples[fTelID];
        }
        float                getLengthOfSampleTimeSlice( unsigned int iTelID )
        {
            if( iTelID < fSample_time_slice.size() )
            {
                return fSample_time_slice[iTelID];
            }
            else
            {
                return 0;
            }
        }
        unsigned int         getNumTelescopes()   //!< get number of telescopes
        {
            return fNTel;
        }
        unsigned int         getNTel()            //!< get number of telescopes
        {
            return fNTel;
        }
        unsigned int         getPixelType()       //!< get pixel type
        {
            return fPixelType;
        }
        //!< get output edge distance of this tube (for current telID)
        float                getOuterEdgeDistance( unsigned int fTubeID );
        vector<float>&       getTubeRadius()      //!< get tube radius
        {
            return fRTube[fTelID];
        }
        //!< get tube radius
        vector<float>&       getTubeRadius( unsigned int iTelID )
        {
            return fRTube[iTelID];
        }
        //!< get tube radius
        vector<float>&       getTubeRadius_MM( unsigned int iTelID )
        {
            return fRTubeMM[iTelID];
        }
        vector<unsigned int> getRealPixel()       //!< get pixel number v[mc]=real
        {
            return fXim;
        }
        unsigned int         getTelID()           //!< return telescope ID
        {
            return fTelID;
        }
        map< unsigned int, unsigned int > getTelIDGrisu()
        {
            return fTelIDGrisu;
        }
        map< unsigned int, unsigned int > getTelID_matrix()
        {
            return getTelIDGrisu();
        }
        vector<float>        getTelRadius()       //!< return vector of radius of telescope dishes
        {
            return fTelRad;
        }
        vector<ULong64_t>&     getTelType()
        {
            return fTelType;
        }
        vector<ULong64_t>      getTelType_list();
        unsigned int         getTelType_Counter( ULong64_t iTelType );
        vector<float>&       getTelXpos()         //!< return vector with x positions of telescopes
        {
            return fTelXpos;
        }
        vector<float>&       getTelYpos()         //!< return vector with y positions of telescopes
        {
            return fTelYpos;
        }
        vector<float>&       getTelZpos()         //!< return vector with z positions of telescopes
        {
            return fTelZpos;
        }
        vector<int>          getTrigger()         //!< get trigger vector
        {
            return fTrigTube[fTelID];
        }
        vector<float>&       getX()               //!< get x-position of tube for current telescope
        {
            return fRotXTube[fTelID];
        }
        //!< get x-position of tube
        vector<float>&       getX( unsigned int iTel )
        {
            return fRotXTube[iTel];
        }
        vector<float>&       getX_MM( unsigned int iTel )
        {
            return fXTubeMM[iTel];
        }
        vector<float>&       getY()               //!< get y-position of tube for current telescope
        {
            return fRotYTube[fTelID];
        }
        //!< get y-position of tube
        vector<float>&       getY( unsigned int iTel )
        {
            return fRotYTube[iTel];
        }
        vector<float>&       getY_MM( unsigned int iTel )
        {
            return fYTubeMM[iTel];
        }
        vector<float>&       getXUnrotated()      //!< get x-position of tube for current telescope
        {
            return fXTube[fTelID];
        }
        //!< get x-position of tube for current telescope
        vector<float>&       getXUnrotated( unsigned int iTel )
        {
            return fXTube[iTel];
        }
        vector<float>&       getYUnrotated()      //!< get y-position of tube for current telescope
        {
            return fYTube[fTelID];
        }
        //!< get x-position of tube for current telescope
        vector<float>&       getYUnrotated( unsigned int iTel )
        {
            return fYTube[iTel];
        }
        
        bool                 initialize( unsigned int i_Ntel, vector< string > iCamera );
        //!< set number of telescopes and channels (only necessary for readCameraFile() )
        bool                 initialize( unsigned int iNtel, unsigned int iNchannel );
        //!< set number of telescopes and channels (only necessary for readCameraFile() )
        bool                 initialize( unsigned int iNtel, vector< unsigned int > iNchannel );
        vector<bool>&        isEdgePixel()
        {
            return fEdgePixel[fTelID];
        }
        bool                 makeNeighbourList( vector< float > iNeighbourSearchScaleFactor,
                                                vector< unsigned int > iNeighbourMultiplicity,
                                                vector< bool > iSquarePixels );
        void                 print( bool bDetailed = true );             //!< print all data vectors to stdout
        //!< read in camera geometry
        bool                 readCameraFile( string );
        //!< read in camera geometry
        bool                 readCameraFile( string, unsigned int );
        //!< read telescope geometry from grisu cfg file
        bool                 readDetectorGeometryFromDB( string iDBStartTime, bool iReadRotationsFromDB = true );
        bool                 readGrisucfg( string iFile, unsigned int fNTel );
        void                 setCameraCentreTubeIndex();
        void                 setConfigDir( string iDir )
        {
            fConfigDir = iDir;
        }
        void                 setFADCRange( int iFADCRange )
        {
            fFADCRange = iFADCRange;
        }
        bool                 setLengthOfSampleTimeSlice( unsigned int iTelID, float iSample_time_slice );
        void                 setTelID_matrix( map< unsigned int, unsigned int > m )
        {
            fTelIDGrisu = m;
        }
        bool                 setLowGainMultiplier_Trace( unsigned int iTel, double ival );
        bool                 setLowGainThreshold( unsigned int iTel, unsigned int ival );
        void                 setCoordinateTransformer( float iX = 1., float iY = 1. )
        {
            fCoordinateTransformerX = iX;
            fCoordinateTransformerY = iY;
        }
        //!< set telescope ID for getters
        void                 setSourceType( int iST )
        {
            fsourcetype = iST;
        }
        bool                 setTelID( unsigned int iID );
        void                 rotateCamera();
        void                 stretchAndMoveCamera();
};
#endif
