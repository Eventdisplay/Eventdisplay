//! VStarCatalogue star catalogue

#ifndef VSTARCATALOGUE_H
#define VSTARCATALOGUE_H

#include "TFile.h"
#include "TMath.h"
#include "TObject.h"
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TSystem.h>
#include "TTree.h"

#include "VAstronometry.h"
#include "VGlobalRunParameter.h"
#include "VSkyCoordinatesUtilities.h"
#include "VStar.h"
#include "VUtilities.h"
#include "VDB_Connection.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class VStarCatalogue : public TObject, public VGlobalRunParameter
{
    private:
        bool   fDebug;
        
        string fCatalogue;
        unsigned int fCatalogueVersion;
        
        vector< VStar* > fStars;
        vector< VStar* > fStarsinFOV;
        
        // telescope pointing
        unsigned int fTel_telescopeID;
        double       fTel_deRotationAngle_deg;
        double       fTel_ra;
        double       fTel_dec;
        double       fTel_camerascale;
        
        bool readCatalogue();
        VStar* readCommaSeparatedLine_Fermi( string, int, VStar* );
        VStar* readCommaSeparatedLine_Fermi_Catalogue( string, int, VStar* );
        VStar* readCommaSeparatedLine_Fermi2nd_Catalogue( string, int, VStar* );
        VStar* readCommaSeparatedLine_FAVA( string, int, VStar* );
        
    public:
    
        VStarCatalogue();
        ~VStarCatalogue();
        bool          init( double MJD );
        bool          init( double MJD, string iCatalogue );
        
        //    void          getStar( unsigned int ID, double &dec, double &ra, double &brightness );
        double        getDistanceToClosestStar( double x_cam_deg, double y_cam_deg );
        unsigned int  getNStar()
        {
            return fStars.size();
        }
        VStar*        getStar( unsigned int ID );
        double        getStarBrightness( unsigned int ID, string iBand = "B" );
        string        getStarCatalogueName()
        {
            return fCatalogue;
        }
        unsigned int  getStarID( unsigned int ID );
        string        getStarName( unsigned int ID );
        double        getStar_l( unsigned int ID );
        double        getStar_b( unsigned int ID );
        double        getStarMajorDiameter( unsigned int ID );
        double        getStarSpectralIndex( unsigned int ID );
        double        getStarDec2000( unsigned int ID );
        double        getStarRA2000( unsigned int ID );
        double        getStarDecCurrentEpoch( unsigned int ID );
        double        getStarRACurrentEpoch( unsigned int ID );
        vector< double > getStarFluxEnergyMin( unsigned int ID );
        vector< double > getStarFluxEnergyMax( unsigned int ID );
        vector< double > getStarFlux( unsigned int ID );
        vector< double > getStarFluxError( unsigned int ID );
        vector< string > getStarOtherNames( unsigned int ID );
        string        getStarType( unsigned int ID );
        vector< string > getStarAssociations( unsigned int ID );
        vector< VStar* > getListOfStars()
        {
            return fStars;
        }
        vector< VStar* > getListOfStarsinFOV()
        {
            return fStarsinFOV;
        }
        void          printCatalogue( unsigned int i_nRows = 0, double iMinBrightness = 999999., string iBand = "B" );
        void          printStarsInFOV();
        void          printStarsInFOV( double iMinBrightness, string iBand = "B" );
        void          purge();
        bool          readVERITASsourcesfromDB( string );
        unsigned int  setFOV( double ra_deg, double dec_deg, double FOV_x, double FOV_y, bool bJ2000 = true, double iBrightness = 9999., string iBand = "B" );
        unsigned int  setFOV( string ra_hour, string dec, double FOV_x, double FOV_y, bool bJ2000 = true );
        void          setTelescopePointing( unsigned int iTelID = 0, double iDerotationAngle = 0.,
                                            double ra_deg = -99., double dec_deg = -99., double iCameraScale = 1. );
        bool          writeCatalogueToRootFile( string iRootFile );
        
        bool          checkTextBlocks( string iL, unsigned int iV );
        
        ClassDef( VStarCatalogue, 8 );
};
#endif
