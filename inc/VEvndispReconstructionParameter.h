//! VEvndispReconstructionParameter event cuts for array analysis

#ifndef VARRAYANALYSISCUTS_H
#define VARRAYANALYSISCUTS_H

#include "VEvndispRunParameter.h"
#include "VImageParameter.h"
#include "VStarCatalogue.h"
#include "VUtilities.h"

#include <TNamed.h>
#include <TSystem.h>

#include <bitset>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

class VEvndispReconstructionCut : public TNamed
{
    private:
    
        bool   testNtubes( int ntubes );
        
    public:
    
        bool   fCutSet;
        int    fCut_ntubes_min;
        int    fCut_ntubes_max;
        
        double fCut_double_min;
        double fCut_double_max;
        int    fCut_int_min;
        int    fCut_int_max;
        
        
        VEvndispReconstructionCut();
        ~VEvndispReconstructionCut() {}
        bool  isCutSet()
        {
            return fCutSet;
        }
        void  print( unsigned int telType, string iCutName );
        bool  test( int iValue, int ntubes );
        bool  test( double iValue, int ntubes );
        void  setCutValues( bool iCutSet = false );
        void  setCutValues( int iMin, int iMax, int ntubes_min = -99999, int ntubes_max = -99999 );
        void  setCutValues( double iMin, double iMax, int ntubes_min = -99999, int ntubes_max = -99999 );
        
        ClassDef( VEvndispReconstructionCut, 1 );
};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

class VEvndispReconstructionParameterData : public TNamed
{
    public:
    
        bool          fDebug;
        
        unsigned int  fMethodID;
        unsigned int  fNTel_type;
        
        int           fNImages_min;
        int           fNImages_max;
        double        fAxesAngles_min;
        
        // [telescope type] (not telescope number!!)
        vector< unsigned int >                           fTelescopeType;
        vector< map< string, VEvndispReconstructionCut* > > fTelescopeTypeCut;
        
        vector< unsigned int > fL2TriggerType;
        vector< bool >         fLocalUseImage;
        
        string            fDISP_MLPFileName;
        vector< double >  fDISP_TMVAZenithBin;
        vector< string >  fDISP_TMVAFileNameVector;
        vector< unsigned int > fDISP_MAXTelescopes;
        vector< unsigned int > fDISP_MAXMethodID;
        
        bool              fUseEventdisplayPointing;
        
        VEvndispReconstructionParameterData() {}
        VEvndispReconstructionParameterData( unsigned int iNtel_type, set< ULong64_t > fTel_type );
        ~VEvndispReconstructionParameterData() {}
        unsigned int      getTelescopeType_counter( unsigned int iTelType );
        bool              isValidKeyword( string iKeyword );
        bool              isValidTelescopeType( unsigned int iTelType );
        void              setDebug( bool iB = false )
        {
            fDebug = iB;
        }
        void              setMethodID( unsigned int iMethodID );
        bool              test( unsigned int iTelType, string iVarName, double iVarD, int iNtubes );
        bool              test( unsigned int iTelType, string iVarName, int iVarI, int iNtubes );
        bool              testUserImage( unsigned int iTelType );
        bool              testL2TriggerType( unsigned int iTel, unsigned int iTelType, unsigned short int iLocalTriggerType );
        
        ClassDef( VEvndispReconstructionParameterData, 4 );
        
};

class VEvndispReconstructionParameter : public TNamed
{
    private:
    
        bool   fDebug;
        
        VEvndispRunParameter*  fRunPara;
        
        unsigned int fNTel_type;
        vector< ULong64_t > fTel_type_perTelescope;         // list of telescope types; one entry per telescope (length #n telescopes)
        set< ULong64_t > fTel_type;
        vector< ULong64_t> fTel_typeUniqueVector;           // list of telescope types; length of vector: #of telescope types;
        vector< VEvndispReconstructionParameterData* > fReconstructionParameterData;  // one element per rec cut / method
        
        void addNewMethod( unsigned int iMethodID );
        bool fillImageCleaningParameter( vector< string > iTemp,
                                         int t_temp,
                                         vector< VImageCleaningRunParameter* >& iImageCleaningParameters );
                                         
        bool readKeyWord_FADCANALYSIS( vector< string > iTemp, int t_temp );
        bool readKeyWord_FADCDOUBLEPASS( vector< string > iTemp, int t_temp );
        bool readKeyWord_FADCSUMMATIONWINDOW( vector< string > iTemp, int t_temp );
        bool readKeyWord_FADCSUMMATIONSTART( vector< string > iTemp, int t_temp );
        bool readKeyWord_CLEANING( vector< string > iTemp, int t_temp );
        bool readKeyWord_BRIGHTSTARS( vector< string > iTemp );
        bool readKeyWord_LLEDGEFIT( vector< string > iTemp, int t_temp );
        bool readKeyWord_FORCELL( vector< string > iTemp );
        bool readKeyWord_CreateIPRdatabase( vector< string > iTemp );
        bool readKeyWord_IPRdatabaseFile( vector< string > iTemp );
        bool readKeyWord_ReadIPRfromDST( vector< string > iTemp );
        bool readKeyWord_ReadIPRfromDatabase( vector< string > iTemp );
        bool readKeyWord_IPRdatabase( vector< string > iTemp );
        bool readKeyWord_WriteGraphsToFile( vector< string > iTemp );
        bool readKeyWord_GraphsFile( vector< string > iTemp );
        bool readKeyWord_RECMETHOD( vector< string > iTemp, int t_temp, ULong64_t t_type );
        bool readKeyWord_NEIGHBOURS( vector< string > iTemp, int t_temp );
        bool readKeyWord_SQUARE( vector< string > iTemp, int t_temp );
        
        void reset();
        
    public:
    
        VEvndispReconstructionParameter();
        VEvndispReconstructionParameter( vector< ULong64_t > itel_type, VEvndispRunParameter* iRunPara );
        ~VEvndispReconstructionParameter() {}
        
        bool             applyArrayAnalysisCuts( unsigned int iMeth, unsigned int iTel, unsigned int iTelType,
                VImageParameter* iImageParameter, unsigned short int iLocalTriggerType,
                VStarCatalogue* iStar = 0 );
        unsigned int     getNReconstructionCuts()
        {
            return fReconstructionParameterData.size();
        }
        VEvndispReconstructionParameterData* getReconstructionParameterData( unsigned int iMethod );
        int              getTelescopeType_counter( ULong64_t t );
        int              getTelescopeType_counter_from_MirrorArea( ULong64_t t );
        int              getTelescopeType_counter_from_MirrorArea_and_PixelSize( ULong64_t t );
        vector <int >    getTelescopeType_counterVector( ULong64_t t );
        vector< int >    getTelescopeType_counter_from_MirrorAreaVector( ULong64_t t );
        vector< int >    getTelescopeType_counter_from_MirrorArea_and_PixelSizeVector( ULong64_t t );
        void             print_arrayAnalysisCuts();
        unsigned int     read_arrayAnalysisCuts( string ifile );
        void             setDebug( bool iD = false )
        {
            fDebug = iD;
        }
        
        ClassDef( VEvndispReconstructionParameter, 29 );
};
#endif
