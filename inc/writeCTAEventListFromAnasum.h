#include "stdint.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TMinuit.h"

#include "VEvndispRunParameter.h"
#include "VSkyCoordinatesUtilities.h"
#include "VStereoAnalysis.h"
#include "VAnaSumRunParameter.h"
#include "VDBRunInfo.h"
#include "VPointingDB.h"
#include "VTimeMask.h"
#include "VDeadTime.h"
#include "VInstrumentResponseFunctionReader.h"

#include "FITSRecord.h"
#include "fitsio.h"

#define KNRM "\x1B[0m"
#define KRED "\x1B[31m"
#define KYEL "\x1B[33m"

const double toRad = M_PI / 180.0 ;
const double toDeg = 1 / toRad ;
const unsigned int NTEL = 4 ;

struct Ntimespan
{
    double beg ; // MET
    double end ; // MET
} ;

struct azbinspec
{
    int azbin ;
    double begaz ;
    double endaz ;
    double centaz ;
    double degwid ;
} ;

// 68%-containment-radius * R68_TO_1SIGMA = Gaussian-Function-Sigma-Parameter
const double R68_TO_1SIGMA = 0.6624305 ;

// psf calculation structs,
// used in calcPSFvalues
struct inputpsf
{
    double emin      ;
    double emax      ;
    double offsetmin ;
    double offsetmax ;
    double psffunc   ; // PSF_KING or PSF_GAUSS
    double contrad68 ;
    double contrad80 ;
} ;
struct outputpsf
{
    double scale   ;
    double sigma_1 ;
    double sigma_2 ;
    double sigma_3 ;
    double ampl_2  ;
    double ampl_3  ;
} ;
// psf calculation types
// gauss function (only sigma and scale(sigma))
// king function (scale, 3x sigmas, and 2x amplitudes)
enum { PSF_GAUSS, PSF_KING } ;

// king parameter reading mode
// KING_GUESS  : Scrape the r68/r80 containment radii and attempt to fit a king
//               function to them.  Warning, tends to fit each parameter space
//               point 50% of the time...
// KING_NATIVE : Get natively-fit king sigma/gamma parameters that were fitted to the
//               same histogram that the r68/r80 values are calculated from,
//               much more reliable, but may not be in old effective area files
enum { KING_GUESS, KING_NATIVE } ;

// effective area calculation structs,
// used in calcEFFAREAvalues()
struct inputeffarea
{
    double emin      ;
    double emax      ;
    double offsetmin ;
    double offsetmax ;
} ;
struct outputeffarea
{
    double effarea      ;
    double effarea_reco ;
} ;

// background calculation structs
struct inputbackground
{
    double camx_min ; // deg
    double camx_max ; // deg
    double camy_min ; // deg
    double camy_max ; // deg
    double emin     ; // TeV
    double emax     ; // TeV
} ;
struct outputbackground
{
    double bgd      ; // 1/s/MeV/sr
    double bgd_reco ; // 1/s/MeV/sr
} ;

// structs for calcAvgPointing()
struct pointinginput
{
    TChain* ch ;
    double beg ;
    double end ;
} ;
struct pointingoutput
{
    double ra  ;
    double dec ;
    double alt ;
    double az  ;
} ;

inline bool file_exists( const std::string& name ) ;
string      exec_shell( char* cmd ) ;
std::string get_env_var( std::string const& key ) ;

void   calcBACKGROUNDvalues( const inputbackground inp, outputbackground& out ) ;
int    calcKingParameters( double rad68, double rad80, double& sigma, double& gamma ) ;            // using newton method
int    calcKingParametersViaSecant( double rad68, double rad80, double& sigma, double& gamma ) ;   // using secant method
int    calcKingParametersRootFit( double rad68, double rad80, double& sigma, double& gamma ) ;     // using TMinuit
void   getHumanDateTimes( VEvndispRunParameter* evrunpar, char dateobs[], char timeobs[],
                          char dateend[], char timeend[] ) ;
void   getRunTimes( VEvndispRunParameter* evrunpar, double met_zero, double& tstart,
                    double& tstop, double& telapse, double& startMJD, double& stopMJD ) ;
double parseTimeArg( char timechar[] ) ;
void   calcAvgPointing( const pointinginput inp, pointingoutput& out ) ;
void   calcGTIsFromTimeMask( VTimeMask* vtm, const double tInterval, const double maxSeconds,
                             const double startMJD, const double stopMJD,
                             vector<double>& gti_beg, vector<double>& gti_end ) ;
void   telmask_clear( bool telmask[] ) ;
void   telmask_set( bool telmask[], int tel, bool state ) ;
bool   readImgSel( int ImgSel, int tel ) ;
void   reduceTimeSpans( const double t1_beg,  const double t1_end,
                        const double t2_beg,  const double t2_end,
                        double&       red_beg, double&       red_end ) ;
void   printChunkList( const double tmin, const double tmax, vector<Ntimespan>& chunklist ) ;
int    checkChunkList( vector<Ntimespan>& chunklist, double time ) ;
void   print_fits_error( int status ) ;
void   checkArgBounds( char flag[], int narg, int iarg ) ;
void   printDeadTimeMask( const int ichunk, const double beg, const double end, const vector<bool> timemask ) ;
string getEffAreaFileName( char* veritas_data_dir, unsigned int run_edver, char* run_epoch, int run_atmo, char* run_cut, float run_wobble, char* sim_type, char* eff_particle, int eff_zenith, int eff_noise ) ;
double getPedVarFromRunSummary( TTree* runsum, int runid ) ;
double getWobbleRadiusFromRunSummary( TTree* runsum, int runid ) ;
double get_pedvar_from_eff( char* fname ) ;
double get_ang_dist( double ang1, double ang2 ) ;
string command_line_regex( string input, string regcmd ) ;
vector<double> get_unique_values_in_branch( string rootfile, string treename, string branchname ) ;

class VMagnetZenith
{
    public :
        string      fEffFile ;
        vector<int> zeniths ;
        VMagnetZenith( string efile ) ;
        string print() ;
        int magnetize( double zen ) ;
} ;

class VMagnetAzimuth
{
    public :
        string         fEffFile ;
        vector<double> azBins ;
        vector<double> azMins ;
        vector<double> azMaxs ;
        vector<double> azCens ;
        VMagnetAzimuth( string efile ) ;
        string print() ;
        int    magnetize( double az ) ;
} ;

struct noiseAndPedvarStruct
{
    int   noise  ;
    float pedvar ;
} ;

class VMagnetNoise
{
    public :
        string         fEffFile  ;
        int            fAzbin    ;
        int            fZenith   ;
        vector<double> pedvars   ;
        vector<int>    noises    ;
        VMagnetNoise( string efile, int azbin, int zenith ) ;
        void printAll() ;
        string printNoises() ;
        int magnetize( double pedvar ) ; // returns the closest noise!
} ;

class VMagnetOffset
{
    public :
        string fEffArea ;
        double binwidth ;
        vector<double> offsets ;
        vector<double> offsetupperedges ;
        vector<double> offsetloweredges ;
        VMagnetOffset( string efile ) ;
        string print()  ;
        void   printAll() ;
        double magnetize( double offset ) ;
        string gernotFormat( double offset ) ;
        
} ;

struct archiveStruct
{
    int    zenith   ;
    int    azbin    ;
    int    noise    ;
    double specind  ;
    double offset   ;
} ;

// convert an archiveStruct struct to a root-compatible string
// (for giving root objects a name
string getArchiveROOTString( archiveStruct archline )
{
    char tmp[1000] = "" ;
    sprintf( tmp, "ze%02d_azbin%02d_offset%04.2f_noise%04d_specind%3.1f",
             archline.zenith, archline.azbin, archline.offset, archline.noise, archline.specind ) ;
    string str = tmp ;
    return str ;
} ;

// for printing a single archive as a justified string
string getArchiveString( archiveStruct archline )
{
    char tmp[1000] = "" ;
    sprintf( tmp, "ze=%02d azbin=%02d offset=%04.2f noise=%04d specind=%3.1f",
             archline.zenith, archline.azbin, archline.offset, archline.noise, archline.specind ) ;
    string str = tmp ;
    return str ;
} ;


// stuct for storing psf information from a single point in the effective area parameter space
struct psfDataStruct
{
    vector<double> energy ; // energy (log10 TeV)
    vector<double> rad68  ; // 60% containment radius (deg)
    vector<double> rad80  ; // 80% containment radius (deg)
    vector<double> nativeking_sigma ; // natively-fit king function sigma
    vector<double> nativeking_gamma ; // natively-fit king function gamma
} ;

const int DATA_IN_ARCHIVE_BUT_NOT_LOADED = -1 ;
const int DATA_NOT_IN_ARCHIVE = -2 ;
class VArchiveEffArea
{
    public:
        string                fEffFile ;
        
        vector<archiveStruct> fArchive ;
        vector<bool>          fLoaded  ;
        vector<psfDataStruct> fPsf     ;
        vector<VInstrumentResponseFunctionReader*> irf ;
        vector<TH2D*>         fenmigmatrix ;
        vector<TH1D*>         fPsf68King ;
        vector<TH1D*>         fPsf68Nominal ;
        vector<int>           fEffEntry ;
        
        VArchiveEffArea( string efile ) ;
        void addArchive( archiveStruct archline ) ;
        int  checkIfLoaded( archiveStruct archline ) ;
        bool checkIfInArchive( archiveStruct archline ) ;
        void onlyOneParameterPointAllowed( archiveStruct archline, vector<int> matent, CEffArea* c ) ;
        void loadAllArchives() ;
        void printAllArchiveFilenames() ;
} ;

// Effective Area Table Format
// _BASIC : 2D array of effarea vs Energy & Offset, effarea vs MC & Reco energies
enum { EFFAREATABLE_BASIC } ;

// Point Spread Function Table Format
// _BASIC : 2D array of gaussian fit parameters vs Energy and Offset
enum { PSFTABLE_BASIC, PSFTABLE_GAUSS } ;

// Energy Dispersion Table Format
enum { ENERGYDISPERSIONTABLE_BASIC } ;

// Background Format
enum { BACKGROUNDTABLE_BASIC, BACKGROUNDTABLE_3D } ;

// background flavors
// BCK_FLAT   : Background() = constant
// BCK_OFFSET : Background(offset) = peak * Exp(-offset^2/sigma^2)
enum { BCK_FLAT, BCK_OFFSET } ;


// for using root libraries to fit a king function when given
// the 68% and 80% containment radii (in deg)
#define KINGFITTER_NPAR 2 // only 2 parameters to fit in a king function, sigma and gamma

// because of how TMinuit does its fitting, you can only provide your
// fitting datapoints via global arrays, and your fit/minimization functions
// as separate functions (and not, as I wanted, as class-member functions).
// My ROOT-induced despair knows no bounds.
Float_t KingFitXArray[   KINGFITTER_NPAR] ;
Float_t KingFitYArray[   KINGFITTER_NPAR] ;
Float_t KingFitYArrayErr[KINGFITTER_NPAR] ;
double kingfunc( Float_t r, Double_t* par ) ;
void   kingfcn( Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag ) ;

// simple class for fitting the 68 and 80% containment radii with
// a king function, built on top of TMinuit
class KingFitter
{
    public :
        // run these functions in the order listed here!
        KingFitter( double r68, double r80 ) ;  // 1
        void   fit() ;                          // 2
        double getSigma()
        {
            return fSigma ;
        } ; // 3   units:deg
        double getGamma()
        {
            return fGamma ;
        } ; // 3   units:none
    private :
        double fRad68 ;
        double fRad80 ;
        double fSigma ;
        double fGamma ;
} ;

struct psfBlockPartStruct
{
    double rad68 ;
    double rad80 ;
    double sigma ;
    double gamma ;
    int    i_offset ;
    int    i_en ;
    int    archindex ;
    bool   goodfit ; // true if fit converged
    bool   wasinterpolated ; // true if value was copied from another point
    double nativeking_sigma ;
    double nativeking_gamma ;
} ;
string formatPsfBlockStructString( psfBlockPartStruct psf );

class VPsfBlock
{
    public:
        vector<vector<psfBlockPartStruct> > fPsfBlock ;
        VPsfBlock( int nOffRows, int nEnRows ) ;
        void formatBlock_goodfit();
        void formatBlock_wasinterpolated();
        void formatBlock_zerorad();
        void formatBlock_sanefit();
        int  scanBlockForAlternateFit( int i_off, int i_en ) ;
        
} ;

float similarity_zenith( int target,   int compare ) ;
float similarity_noise( int target,   int compare ) ;
float similarity_specind( float target, float compare ) ;
float similarity_wobble( float target, float compare ) ;
float similarity_azbin( int target,   int compare ) ;

