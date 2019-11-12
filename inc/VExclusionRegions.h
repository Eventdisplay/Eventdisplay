//! VExclusionRegions  handling of exclusion regions due to e.g. bright stars

#ifndef VExclusionRegions_H
#define VExclusionRegions_H

#include "VSkyCoordinatesUtilities.h"
#include "VStarCatalogue.h"

#include "TFile.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

///////////////////////////////////////////////
// setting of bright stars

class VBrightStarExclusionSettings
{
    public:
        double fStarMinBrightness;
        string fStarBand;
        double fStarExlusionRadius_DEG;
        
        VBrightStarExclusionSettings();
        virtual ~VBrightStarExclusionSettings();
        
        bool operator<( const VBrightStarExclusionSettings& s1 ) const;
};

///////////////////////////////////////////////
// list of exclusion regions

class VListOfExclusionRegions
{
    public:
    
        double fExcludeFromBackground_North;    //[deg] relative to sky map centre
        double fExcludeFromBackground_West;     //[deg] relative to sky map centre
        double fExcludeFromBackground_CameraCentre_x;   // [deg] relative to camera centre
        double fExcludeFromBackground_CameraCentre_y;   // [deg] relative to camera centre
        double fExcludeFromBackground_DecJ2000; //[deg]
        double fExcludeFromBackground_RAJ2000;  //[deg]
        double fExcludeFromBackground_Radius1;   //[deg]
        double fExcludeFromBackground_Radius2;   //[deg]
        double fExcludeFromBackground_RotAngle;  //[deg]
        int    fExcludeFromBackground_StarID;
        string fExcludeFromBackground_StarName;
        bool   fExcludeFromBackground_B_Band;
        double fExcludeFromBackground_StarBrightness_V;
        double fExcludeFromBackground_StarBrightness_B;
        
        bool operator<( const VListOfExclusionRegions& s1 ) const;
        
        VListOfExclusionRegions();
        virtual ~VListOfExclusionRegions();
        bool isInsideExclusionRegion( double x_cam_deg, double y_cam_deg,
                                      double r_offset_deg = 0. );
};

//////////////////////////////////////////////
// exlusion regions handler

class VExclusionRegions
{
    private:
    
        bool setExclusionRegionsFromStarCatalogue( VStarCatalogue* iStarCatalogue );
        
    public:
    
        // list of exclusion regions
        vector< VListOfExclusionRegions* > fExclusionRegions;
        
        // star exclusion regions
        string fStarCatalogue;
        vector< VBrightStarExclusionSettings* > fBrightStarSettings;
        
        
        VExclusionRegions();
        virtual ~VExclusionRegions();
        
        void addBrightStarSettings( double iStarMinBrightness,
                                    double iStarExlusionRadius_DEG = -1.,
                                    string iStarBand = "" );
        void addBrightStarSettings( vector< VBrightStarExclusionSettings* > iBS )
        {
            fBrightStarSettings = iBS;
        }
        void addExclusionRegion( double iNorth, double iWest,
                                 double iRa, double iDec,
                                 double iR1, double iR2, double iAngle,
                                 string iName = "" );
        vector< VBrightStarExclusionSettings* > getBrightStarSettings()
        {
            return fBrightStarSettings;
        }
        vector< VListOfExclusionRegions* > getListOfExclusionRegions()
        {
            return fExclusionRegions;
        }
        double getLargestStarExlusionRadius();
        bool initializeExclusionRegions( VStarCatalogue* iCatalogue,
                                         double iSkyMapCentre_ra_deg,
                                         double iSkyMapCentre_dec_deg,
                                         double iCameraCentre_ra_deg,
                                         double iCameraCentre_dec_deg );
        void printBrightStarSettings();
        void printExclusionRegions();
        bool readExclusionRegionTree( TFile* f, int iRunOn );
        void setStarCatalogue( string iCatalogue )
        {
            fStarCatalogue = iCatalogue;
        }
        bool writeExclusionRegionTree( int iOnRun = -1 );
        
        ClassDef( VExclusionRegions, 1 );
        
};

#endif
