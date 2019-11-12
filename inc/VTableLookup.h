//! VTableLookup calculation of mean scaled variables and energies using MC filled tables

#ifndef VTABLELOOKUP
#define VTABLELOOKUP

#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"

#include "VStatistics.h"
#include "VTableLookupDataHandler.h"
#include "VTableLookupRunParameter.h"
#include "VTablesToRead.h"
#include "VTableCalculator.h"

#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

using namespace std;

class VTableCalculatorData
{
    public:
    
        bool   fEnergy;
        string fDirectoryName;      // lookup table is stored in this directory in the table file
        string fFillVariable;       // 1D variable name
        
        // n-dimensional vectors [tel_type][NSB][ze][woff][az]
        vector< vector< vector< vector< vector< VTableCalculator* > > > > > fTable;
        
        double fValueNormalizationRange_min;
        double fValueNormalizationRange_max;
        
        VTableCalculatorData();
        ~VTableCalculatorData() {}
        bool    assertTableVector( unsigned int wobble_bin );
        void    print();
        void    terminate( TFile* iFile );
        
};

class VTableLookup
{
    private:
        // list of lookup table types
        // (E_MSCW, E_MSCL, E_EREC are expected to be the first three types)
        enum E_TLvalue { E_MSCW, E_MSCL, E_EREC, E_TGRA, E_FRGO };
        
        VTableLookupRunParameter* fTLRunParameter; // lookup table run parameter
        
        int fNTel;
        
        //////////////////////////////////////////////////////////////////////
        // event data
        VTableLookupDataHandler* fData;          // event data tree
        int fNumberOfIgnoredEvents;              // number of events ignored due to evendisp error code
        
        //////////////////////////////////////////////////////////////////////
        // lookup tables
        
        // root file with lookup tables and pointers to directories
        TFile* fLookupTableFile;
        
        // lookup table parameter space
        
        // definition of azimuth bins (same for all tables)
        vector< double > fTableAzLowEdge;
        vector< double > fTableAzUpEdge;
        
        // telescope types
        vector< ULong64_t >                     fTableTelTypes;
        // NSB level [tel_type]
        vector< vector< double > >              fTableNoiseLevel;
        // zenith angle [tel_type][NSB][ze]
        vector< vector< vector< double > > >    fTableZe;
        // wobble offsets [tel_type][NSB][ze][woff]
        vector< vector< vector< vector< double > > > > fTableDirectionOffset;
        
        // lookup tables
        map< unsigned int, VTableCalculatorData* > fTableData;
        
        // used for calculations
        VTableCalculator*       fTableCalculator;
        
        // big and ugly
        VTablesToRead* s_NupZupWup;
        VTablesToRead* s_NupZupWlow;
        VTablesToRead* s_NupZup;
        VTablesToRead* s_NupZlowWup;
        VTablesToRead* s_NupZlowWlow;
        VTablesToRead* s_NupZlow;
        VTablesToRead* s_Nup;
        VTablesToRead* s_NlowZupWup;
        VTablesToRead* s_NlowZupWlow;
        VTablesToRead* s_NlowZup;
        VTablesToRead* s_NlowZlowWup;
        VTablesToRead* s_NlowZlowWlow;
        VTablesToRead* s_NlowZlow;
        VTablesToRead* s_Nlow;
        VTablesToRead* s_N;
        
        //////////////////////////////////////////////////////////////////////
        // intermediate values which change from event to event
        vector< bool > fTelToAnalyze;         // telescopes with data in this event
        double         fMeanNoiseLevel;       // pedeval variations (average over all telescopes)
        vector< double > fNoiseLevel;         // pedestal variances per telescope
        unsigned int   fNNoiseLevelWarnings;
        
        //////////////////////////////////////////////////////////////////////
        // private functions
        
        void             calculateMSFromTables( VTablesToRead* s );
        bool             cut( bool bWrite = false );  // apply cuts on successfull reconstruction to input data
        void             fillLookupTable();
        int              getAzBin( double az );
        void             getIndexBoundary( unsigned int* ib, unsigned int* il, vector< double >& iV, double x );
        vector< string > getSortedListOfDirectories( TDirectory* );
        void             getTables( unsigned int inoise, unsigned int ize, unsigned int iwoff, unsigned int iaz, unsigned int tel, VTablesToRead* s );
        unsigned int     getTelTypeCounter( unsigned int iTel, bool iStopIfError = false );
        unsigned int     getWobbleBin( double w );
        void             initializeLookupTableDataVector();
        void             interpolate( VTablesToRead* s1, double w1, VTablesToRead* s2, double w2, VTablesToRead* s, double w, bool iCos = false );
        void             readLookupTable();
        void             readNoiseLevel( bool bWriteToRunPara = true ); // read noise level from pedvar histograms of data files
        bool             sanityCheckLookupTableFile( bool iPrint = false );
        bool             setInputFiles( vector< string > iInputFiles ); // set input files from evndisp
        void             setMCTableFiles_forTableReading( string, string ); // set MC table file names (reading tables)
        void             setMCTableFiles_forTableWriting( string, double, int, map< ULong64_t, double >, string, string ); // set MC table file names
        void             setOutputFile( string outputfile, string writeoption, string tablefile ); // set file for output tree with mscw, energy, etc.
        void             setSpectralIndex( double iS );
        
    public:
        VTableLookup( VTableLookupRunParameter* iTLRunParameter );
        ~VTableLookup() {}
        
        double getMaxTotalTime()
        {
            return fData->getMaxTotalTime();
        }
        Long64_t getNEntries()
        {
            if( fData )
            {
                return fData->getNEntries();
            }
            else
            {
                return 0;
            }
        }
        bool   initialize();
        void   loop();                              // loop over all events
        void   terminate();
        
};
#endif
