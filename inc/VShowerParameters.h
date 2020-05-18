//! VShowerParameters storage class for shower data

#ifndef VSHOWERPARAMETERS_H
#define VSHOWERPARAMETERS_H

#include "VGlobalRunParameter.h"
#include "TTree.h"

#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

using namespace std;

class VShowerParameters
{
    private:
        bool fDebug;
        bool fMC;
        TTree* fTreeSC;                           //!< output tree
        unsigned int fNTel;                       //!< number of telescopes
        
        unsigned int fShortTree;
        
        void resetDispVectors();
        
    public:
    
        int runNumber;
        int eventNumber;
        int MJD;
        double time;
        
        unsigned int eventStatus;
        
        // geometry and run parameters
        unsigned int fNTelescopes;                //!< number of telescopes
        
        // pointing information (array pointing)
        float fArrayPointing_Elevation;
        float fArrayPointing_Azimuth;
        float fArrayPointing_deRotationAngle_deg;
        // pointing information per telescope
        // calculated from time and source position
        float fTelElevation[VDST_MAXTELESCOPES];
        float fTelAzimuth[VDST_MAXTELESCOPES];
        // from vbf
        float fTelElevationVBF[VDST_MAXTELESCOPES];
        float fTelAzimuthVBF[VDST_MAXTELESCOPES];
        // difference between calculation and vbf
        float fTelPointingMismatch[VDST_MAXTELESCOPES];
        float fTelDec[VDST_MAXTELESCOPES];
        float fTelRA[VDST_MAXTELESCOPES];
        float fTelPointingErrorX[VDST_MAXTELESCOPES];
        float fTelPointingErrorY[VDST_MAXTELESCOPES];
        
        unsigned int fNMethods;
        // array reconstruction method
        uint16_t fMethodID[VDST_MAXRECMETHODS];
        
        // trigger information
        unsigned int fNTrig;                      //!< number of telescopes with local triggers
        ULong64_t fLTrig;                         //!< local trigger vector, bit coded
        //!< list of telescopes with local triggers (length of array is fNTrig)
        unsigned short int fTrig_list[VDST_MAXTELESCOPES];
        unsigned short int fTrig_type[VDST_MAXTELESCOPES];
        // reconstruction parameters
        unsigned int fNumImages;                  //!< total number of images
        // C. Duke 19Oct06
        //!< images selected bitcoded
        ULong64_t fTelIDImageSelected_bitcode[VDST_MAXRECMETHODS];
        //!< list of telescopes with images selected (length of array is fNMethods)
        uint8_t fTelIDImageSelected_list[VDST_MAXRECMETHODS][VDST_MAXTELESCOPES];
        float fTel_x_SC[VDST_MAXTELESCOPES];      //!< array of telescope x location in array plane
        float fTel_y_SC[VDST_MAXTELESCOPES];      //!< array of telescope x location in array plane
        float fTel_z_SC[VDST_MAXTELESCOPES];      //!< array of telescope x location in array plane
        
        float fiangdiff[VDST_MAXRECMETHODS];
        
        //!< true, if image of telescope is used for the reconstruction (one for each method)
        vector< vector< bool > > fTelIDImageSelected;
        //!< number of images used for shower reconstruction
        uint16_t fShowerNumImages[VDST_MAXRECMETHODS];
        float fShowerZe[VDST_MAXRECMETHODS];      //!< shower zenith angle in [deg]
        float fShowerAz[VDST_MAXRECMETHODS];      //!< shower azimuth angle in [deg]
        
        //------------------------------
        float fWobbleNorth;
        float fWobbleEast;
        //----------------------------
        float fDec[VDST_MAXRECMETHODS];           //!< declination of shower [deg]
        float fRA[VDST_MAXRECMETHODS];            //!< right ascension of shower [deg]
        float fShower_Xoffset[VDST_MAXRECMETHODS];//!< offset of direction from camera center [deg]
        float fShower_Yoffset[VDST_MAXRECMETHODS];//!< offset of direction from camera center [deg]
        vector< vector< float > > fShower_Xoff_DISP;
        vector< vector< float > > fShower_Yoff_DISP;
        vector< vector< float > > fShower_Weight_DISP;
        
        // (debug use only)
        unsigned int fShower_NPair;               //! number of pairs
        float fShower_PairXS[VDST_MAXRECMETHODS];  //! pairwise X (observe: debugging only) (geo)
        float fShower_PairYS[VDST_MAXRECMETHODS];  //! pairwise Y (observe: debugging only) (geo)
        float fShower_PairXD[VDST_MAXRECMETHODS];  //! pairwise X (observe: debugging only) (disp)
        float fShower_PairYD[VDST_MAXRECMETHODS];  //! pairwise Y (observe: debugging only) (disp)
        float fShower_PairAngDiff[VDST_MAXRECMETHODS];  //! pairwise angdiff (observe: debugging only)
        float fShower_PairDispWeight[VDST_MAXRECMETHODS];  //! pairwise disp weight (observe: debugging only)
        // (end debug use only)
        
        //!< offset of direction from camera center [deg] (derotated)
        float fShower_XoffsetDeRot[VDST_MAXRECMETHODS];
        //!< offset of direction from camera center [deg] (derotated);
        float fShower_YoffsetDeRot[VDST_MAXRECMETHODS];
        float fShower_stdS[VDST_MAXRECMETHODS];   //!< std (radius) about source point
        float fShowerXcore[VDST_MAXRECMETHODS];   //!< shower core in ground coordinates
        float fShowerYcore[VDST_MAXRECMETHODS];   //!< shower core in ground coordinates
        float fShowerXcore_SC[VDST_MAXRECMETHODS];//!< shower core in shower coordinates
        float fShowerYcore_SC[VDST_MAXRECMETHODS];//!< shower core in shower coordinates
        float fShower_stdP[VDST_MAXRECMETHODS];   //!<  std (radius) about impact point
        float fShower_Chi2[VDST_MAXRECMETHODS];   //!<  chi2 value where appropriate, < 0. for no reconstruction (-99. or angle between lines for two-images events)
        
        float fDispDiff[VDST_MAXRECMETHODS];      //!< difference in disp event direction between telescopes
        
        int MCprimary;
        float MCenergy;                           //!< MC energy in [TeV]
        float MCxcore;                            //!< MC core position in ground coordinates (x)
        float MCycore;                            //!< MC core position in ground coordinates (y)
        float MCzcore;                            //!< MC core position in ground coordinates (z)
        float MCxcos;
        float MCycos;
        float MCaz;
        float MCze;
        float MCTel_Xoff;                         //!< MC source offset in MC in deg (grisudet telescope coordinate system)
        float MCTel_Yoff;                         //!< MC source offset in MC in deg (grisudet telescope coordinate system)
        
        float MCxcore_SC;                         //!< MC core position in shower coordinates
        float MCycore_SC;                         //!< MC core position in shower coordinates
        float MCzcore_SC;                         //!< MC core position in shower coordinates
        int	       MCCorsikaRunID;
        int	       MCCorsikaShowerID;
        float          MCFirstInteractionHeight;
        float	       MCFirstInteractionDepth;
        
        VShowerParameters( int iNTel = 4, unsigned int iShortTree = 0, unsigned int iNMethods = 1 );
        ~VShowerParameters();
        
        void addDISPPoint( unsigned int iTelID, unsigned int iMethod, float x, float y, float idispw = -999. );
        void fill()
        {
            if( fTreeSC )
            {
                fTreeSC->Fill();
            }
        }
        unsigned int getMaxNTelescopes() const
        {
            return VDST_MAXTELESCOPES;
        }
        unsigned int getNArrayReconstructionMethods() const
        {
            return fNMethods;
        }
        unsigned int getNArrayMaxReconstructionMethods() const
        {
            return VDST_MAXRECMETHODS;
        }
        TTree* getTree()
        {
            return fTreeSC;
        }
        void initTree( string, string, bool );
        bool isMC()                               //!< is data Monte Carlo?
        {
            return fMC;
        }
        void printParameters();                   //!< write tree parameters to standard output
        void reset();                             //!< reset all tree variable to standard values
        void reset( unsigned int iNTel );         //!< reset all tree variable to standard values
        void setMC()
        {
            fMC = true;
        }
        void setNArrayReconstructionMethods( unsigned int iNMethods )
        {
            fNMethods = iNMethods;
        }
};
#endif
