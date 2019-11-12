
#include <string>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

// Model3D best-fit (shower) //
Float_t f2Sel3D;	// JG: elevation of shower direction (deg)
Float_t f2Saz3D;	// JG: azimuth of shower direction (deg)
Float_t f2Xcore3D;	// JG: shower core in ground coordinates
Float_t f2Ycore3D;	// JG: shower core in ground coordinates
Float_t f2Smax3D;	// JG: height of shower maximum (along the shower axis)
Float_t f2sigmaL3D;	// JG: longitudinal (3D-length)
Float_t f2sigmaT3D;	// JG: transverse (3D-width)
Float_t f2Nc3D;		// JG: total number of Cherenkov photons emitted by the shower
Float_t f2Goodness3D;	// NAH, Goodness i.e how much a shower is gamma like
Float_t f2Xoff3D;	// NAH, Offset 3D
Float_t f2Yoff3D;	// NAH, Offset 3D
Float_t f2MCe0;
Float_t f2RWidth3D;
Float_t f2Depth3D;

Float_t f2EnergySel3D;
Float_t f2EnergySaz3D;
Float_t f2EnergyXcore3D;
Float_t f2EnergyYcore3D;
Float_t f2EnergySmax3D;
Float_t f2EnergysigmaL3D;
Float_t f2EnergysigmaT3D;
Float_t f2EnergyNc3D;
Float_t f2EnergyGoodness3D;
Float_t f2EnergyXoff3D;
Float_t f2EnergyYoff3D;
Float_t f2EnergyRWidth3D;
Float_t f2EnergyDepth3D;

Float_t f2ReturnSel3D;
Float_t f2ReturnSaz3D;
Float_t f2ReturnsigmaL3D;
Float_t f2ReturnsigmaT3D;
Float_t f2ReturnSmax3D;
Float_t f2ReturnNc3D;
Float_t f2ReturnRWidth3D;
Float_t f2ReturnDepth3D;
Float_t f2ReturnGoodness3D;
Float_t f2ReturnXoff3D;
Float_t f2ReturnYoff3D;
Float_t f2ReturnXcore3D;
Float_t f2ReturnYcore3D;

Float_t f2PlasementSel3D;
Float_t f2PlasementSaz3D;
Float_t f2PlasementsigmaL3D;
Float_t f2PlasementsigmaT3D;
Float_t f2PlasementSmax3D;
Float_t f2PlasementNc3D;
Float_t f2PlasementRWidth3D;
Float_t f2PlasementDepth3D;
Float_t f2PlasementXcore3D;
Float_t f2PlasementYcore3D;
Float_t f2PlasementGoodness3D;
Float_t f2PlasementXoff3D;
Float_t f2PlasementYoff3D;

Float_t f2ErrorSel3D;
Float_t f2ErrorSaz3D;
Float_t f2ErrorXcore3D;
Float_t f2ErrorYcore3D;
Float_t f2ErrorSmax3D;
Float_t f2ErrorsigmaL3D;
Float_t f2ErrorsigmaT3D;
Float_t f2ErrorNc3D;
Float_t f2ErrorRWidth3D;

Float_t fSel3D;
Float_t fsigmaL3D;
Float_t fsigmaT3D;
Float_t fNc3D;
Float_t fRWidth3D;
Float_t fSaz3D;
Float_t fSmax3D;
Float_t fDepth3D;
Float_t fXcore3D;
Float_t fYcore3D;
Float_t fXoff3D;
Float_t fYoff3D;
Float_t fGoodness3D;
Float_t fMCe0;

Float_t fErrorSel3D;
Float_t fErrorsigmaL3D;
Float_t fErrorsigmaT3D;
Float_t fErrorNc3D;
Float_t fErrorRWidth3D;
Float_t fErrorSaz3D;
Float_t fErrorSmax3D;
Float_t fErrorXcore3D;
Float_t fErrorYcore3D;

using namespace std;

char NotIncluded[20000] = {"These were NOT included:"};
char IncludedSize[20000] = {"These where the sizes:"};

bool compare( const pair<float_t, float_t>& i, const pair<float_t, float_t>& j )
{
    return i.first < j.first;
}

bool compareInt( const pair<Int_t, float_t>& i, const pair<Int_t, float_t>& j )
{
    return i.first < j.first;
}

int main( int argc, char* argv[] )
{
    //freopen("Iwattoknow.txt","w",stdout);
    string Sim;
    string TelS;
    string SumWin;
    string Oljud;
    string InOut;
    if( argc == 6 )
    {
        Sim    = argv[1];
        TelS   = argv[2];
        SumWin = argv[3];
        Oljud  = argv[4];
        InOut  = argv[5];
    }
    else
    {
        cout << "usage: ./create_energy3d_referencetables [Sim] [TelS] [SumWin] [Noice] [INOUT]" << endl;
        cout << "./create_energy3d_referencetables CARE 6 21 50 $YOURDIR/YOURSUBDIR" << endl;
        return 0;
    }
    
    char* SIMTYPE = new char[Sim.size() + 1 ];
    SIMTYPE[ Sim.size() ] = 0;
    memcpy( SIMTYPE, Sim.c_str(), Sim.size() );
    
    char* TELSETUP = new char[TelS.size() + 1 ];
    TELSETUP[ TelS.size() ] = 0;
    memcpy( TELSETUP, TelS.c_str(), TelS.size() );
    
    char* SUMWIN = new char[SumWin.size() + 1 ];
    SUMWIN[ SumWin.size() ] = 0;
    memcpy( SUMWIN, SumWin.c_str(), SumWin.size() );
    
    char* NOICE = new char[Oljud.size() + 1 ];
    NOICE[ Oljud.size() ] = 0;
    memcpy( NOICE, Oljud.c_str(), Oljud.size() );
    
    char* INOUT = new char[InOut.size() + 1 ];
    INOUT[ InOut.size() ] = 0;
    memcpy( INOUT, InOut.c_str(), InOut.size() );
    
    const char* Deg[10] 	= 	{"00", "20", "30", "35", "40", "45", "50", "55", "60", "65"};
    int deg[10] = 			{90, 70, 60, 55, 50, 45, 40, 35, 30, 25};
    
    //Rw	Sel	 saz	    Xcor	Ycor      Smax   sigL   sigT   Nc
    float_t Prospars[10][9] = { 	{0.003, 0.0001,  0.0000001, 0.00003,    0.0001,   0.002, 0.02,  0.005, 0.0005},//00ok
        {0.004, 0.00015, 0.00025,   0.0000015,  0.000015, 0.04,  0.005, 0.005, 0.0015},//20ok
        {0.005, 0.0002,  0.0002,    0.0000015,  0.000015, 0.045, 0.006, 0.005, 0.0015},//30ok
        {0.005, 0.00022, 0.00019,   0.0000014,  0.000015, 0.045, 0.010, 0.007, 0.0015},//35
        {0.006, 0.00023, 0.00018,   0.0000012,  0.000015, 0.045, 0.016, 0.008, 0.0015},//40
        {0.007, 0.00025, 0.00017,   0.000001,   0.000015, 0.046, 0.02,  0.01,  0.0015},//45ok
        {0.008, 0.00028, 0.00016,   0.0000007,  0.000020, 0.046, 0.05,  0.01,  0.0015},//50ok
        {0.010, 0.0003,  0.00015,   0.0000004,  0.000025, 0.047, 0.08,  0.01,  0.0015},//55ok
        {0.012, 0.00032, 0.00012,   0.00000025, 0.00005,  0.048, 0.16,  0.015, 0.0015},//60ok
        {0.015, 0.00035, 0.0001,    0.00000015, 0.0001,   0.05,  0.2,   0.02,  0.0015}
    };//65ok
    
    int nentries;
    
    char Outname[150];
    sprintf( Outname, "%s/energy3d_referencetables_V%s_ATM%s_NOISE%s_%s.root", INOUT, TELSETUP, SUMWIN, NOICE, SIMTYPE );
    
    //NAH, creating the outputfile
    TFile* newfile =  new TFile( Outname, "recreate" );
    
    //NAH, create the tree that is wanted
    TTree* ModelPars3D = new TTree( "ModelPars3D", "ModelPars3D" );
    
    //NAH, setting branches and subbranches, also known as leaves
    ModelPars3D->Branch( "Sel3D", 		&f2Sel3D, 	"fSel3D/F" );
    ModelPars3D->Branch( "sigmaL3D", 	&f2sigmaL3D, 	"fsigmaL3D/F" );
    ModelPars3D->Branch( "sigmaT3D", 	&f2sigmaT3D, 	"fsigmaT3D/F" );
    ModelPars3D->Branch( "Nc3D", 		&f2Nc3D, 	"fNc3D/F" );
    ModelPars3D->Branch( "RWidth3D", 	&f2RWidth3D, 	"fRWidth3D/F" );
    ModelPars3D->Branch( "Saz3D",		&f2Saz3D, 	"fSaz3D/F" );
    ModelPars3D->Branch( "Smax3D",		&f2Smax3D, 	"fSmax3D/F" );
    ModelPars3D->Branch( "Depth3D", 	&f2Depth3D, 	"fDepth3D/F" );
    ModelPars3D->Branch( "Xcore3D", 	&f2Xcore3D, 	"fXcore3D/F" );
    ModelPars3D->Branch( "Ycore3D", 	&f2Ycore3D, 	"fYcore3D/F" );
    ModelPars3D->Branch( "Xoff3D", 		&f2Xoff3D, 	"fXoff3D/F" );
    ModelPars3D->Branch( "Yoff3D", 		&f2Yoff3D, 	"fYoff3D/F" );
    ModelPars3D->Branch( "Goodness3D", 	&f2Goodness3D, 	"fGoodness3D/F" );
    ModelPars3D->Branch( "MCe0", 		&f2MCe0, 	"fMCe0/F" );
    
    ModelPars3D->Branch( "EnergySel3D", 	&f2EnergySel3D, 	"fEnergySel3D/F" );
    ModelPars3D->Branch( "EnergysigmaL3D", 	&f2EnergysigmaL3D, 	"fEnergysigmaL3D/F" );
    ModelPars3D->Branch( "EnergysigmaT3D", 	&f2EnergysigmaT3D, 	"fEnergysigmaT3D/F" );
    ModelPars3D->Branch( "EnergyNc3D", 	&f2EnergyNc3D, 		"fEnergyNc3D/F" );
    ModelPars3D->Branch( "EnergyRWidth3D", 	&f2EnergyRWidth3D, 	"fEnergyRWidth3D/F" );
    ModelPars3D->Branch( "EnergySmax3D", 	&f2EnergySmax3D, 	"fEnergySmax3D/F" );
    ModelPars3D->Branch( "EnergySaz3D", 	&f2EnergySaz3D, 	"fEnergySaz3D/F" );
    ModelPars3D->Branch( "EnergyDepth3D", 	&f2EnergyDepth3D, 	"fEnergyDepth3D/F" );
    ModelPars3D->Branch( "EnergyXcore3D", 	&f2EnergyXcore3D, 	"fEnergyXcore3D/F" );
    ModelPars3D->Branch( "EnergyYcore3D", 	&f2EnergyYcore3D, 	"fEnergyYcore3D/F" );
    ModelPars3D->Branch( "EnergyXoff3D", 	&f2EnergyXoff3D, 	"fEnergyXoff3D/F" );
    ModelPars3D->Branch( "EnergyYoff3D", 	&f2EnergyYoff3D, 	"fEnergyYoff3D/F" );
    ModelPars3D->Branch( "EnergyGoodness3D", &f2EnergyGoodness3D, 	"fEnergyGoodness3D/F" );
    
    ModelPars3D->Branch( "ReturnSel3D", 	&f2ReturnSel3D, 	"fReturnSel3D/F" );
    ModelPars3D->Branch( "ReturnsigmaL3D", 	&f2ReturnsigmaL3D, 	"fReturnsigmaL3D/F" );
    ModelPars3D->Branch( "ReturnsigmaT3D", 	&f2ReturnsigmaT3D, 	"fReturnsigmaT3D/F" );
    ModelPars3D->Branch( "ReturnNc3D", 	&f2ReturnNc3D, 		"fReturnNc3D/F" );
    ModelPars3D->Branch( "ReturnRWidth3D", 	&f2ReturnRWidth3D, 	"fReturnRWidth3D/F" );
    ModelPars3D->Branch( "ReturnSaz3D", 	&f2ReturnSaz3D, 	"fReturnSaz3D/F" );
    ModelPars3D->Branch( "ReturnSmax3D", 	&f2ReturnSmax3D, 	"fReturnSmax3D/F" );
    ModelPars3D->Branch( "ReturnDepth3D",	&f2ReturnDepth3D, 	"fReturnDepth3D/F" );
    ModelPars3D->Branch( "ReturnXcore3D", 	&f2ReturnXcore3D, 	"fReturnXcore3D/F" );
    ModelPars3D->Branch( "ReturnYcore3D", 	&f2ReturnYcore3D, 	"fReturnYcore3D/F" );
    ModelPars3D->Branch( "ReturnXoff3D",	&f2ReturnXoff3D, 	"fReturnXoff3D/F" );
    ModelPars3D->Branch( "ReturnYoff3D", 	&f2ReturnYoff3D, 	"fReturnYoff3D/F" );
    ModelPars3D->Branch( "ReturnGoodness3D", &f2ReturnGoodness3D, 	"fReturnGoodness3D/F" );
    
    ModelPars3D->Branch( "PlasementSel3D", 		&f2PlasementSel3D, 	"fPlasementSel3D/F" );
    ModelPars3D->Branch( "PlasementsigmaL3D", 	&f2PlasementsigmaL3D, 	"fPlasementsigmaL3D/F" );
    ModelPars3D->Branch( "PlasementsigmaT3D", 	&f2PlasementsigmaT3D, 	"fPlasementsigmaT3D/F" );
    ModelPars3D->Branch( "PlasementNc3D", 		&f2PlasementNc3D, 	"fPlasementNc3D/F" );
    ModelPars3D->Branch( "PlasementRWidth3D", 	&f2PlasementRWidth3D, 	"fPlasementRWidth3D/F" );
    ModelPars3D->Branch( "PlasementSaz3D", 		&f2PlasementSaz3D, 	"fPlasementSaz3D/F" );
    ModelPars3D->Branch( "PlasementSmax3D", 	&f2PlasementSmax3D, 	"fPlasementSmax3D/F" );
    ModelPars3D->Branch( "PlasementDepth3D", 	&f2PlasementDepth3D, 	"fPlasementDepth3D/F" );
    ModelPars3D->Branch( "PlasementXcore3D", 	&f2PlasementXcore3D, 	"fPlasementXcore3D/F" );
    ModelPars3D->Branch( "PlasementYcore3D", 	&f2PlasementYcore3D, 	"fPlasementYcore3D/F" );
    ModelPars3D->Branch( "PlasementXoff3D", 	&f2PlasementXoff3D, 	"fPlasementXoff3D/F" );
    ModelPars3D->Branch( "PlasementYoff3D", 	&f2PlasementYoff3D, 	"fPlasementYoff3D/F" );
    ModelPars3D->Branch( "PlasementGoodness3D", 	&f2PlasementGoodness3D, "fPlasementGoodness3D/F" );
    
    ModelPars3D->Branch( "ErrorSel3D",	&f2ErrorSel3D, 		"fErrorSel3D/F" );
    ModelPars3D->Branch( "ErrorsigmaL3D", 	&f2ErrorsigmaL3D, 	"fErrorsigmaL3D/F" );
    ModelPars3D->Branch( "ErrorsigmaT3D", 	&f2ErrorsigmaT3D, 	"fErrorsigmaT3D/F" );
    ModelPars3D->Branch( "ErrorNc3D", 	&f2ErrorNc3D, 		"fErrorNc3D/F" );
    ModelPars3D->Branch( "ErrorRWidth3D", 	&f2ErrorRWidth3D, 	"fErrorRWidth3D/F" );
    ModelPars3D->Branch( "ErrorSaz3D", 	&f2ErrorSaz3D, 		"fErrorSaz3D/F" );
    ModelPars3D->Branch( "ErrorSmax3D", 	&f2ErrorSmax3D, 	"fErrorSmax3D/F" );
    ModelPars3D->Branch( "ErrorXcore3D",	&f2ErrorXcore3D, 	"fErrorXcore3D/F" );
    ModelPars3D->Branch( "ErrorYcore3D", 	&f2ErrorYcore3D, 	"fErrorYcore3D/F" );
    
    
    int addnentries = 0;
    int abba = 0;
    char SizeName[20000];
    
    vector < pair< float_t, float_t > > Sel3DForSorting;
    vector < pair< float_t, float_t > > sigmaL3DForSorting;
    vector < pair< float_t, float_t > > sigmaT3DForSorting;
    vector < pair< float_t, float_t > > Nc3DForSorting;
    vector < pair< float_t, float_t > > RWidth3DForSorting;
    vector < pair< float_t, float_t > > Saz3DForSorting;
    vector < pair< float_t, float_t > > Smax3DForSorting;
    vector < pair< float_t, float_t > > Depth3DForSorting;
    vector < pair< float_t, float_t > > Xcore3DForSorting;
    vector < pair< float_t, float_t > > Ycore3DForSorting;
    vector < pair< float_t, float_t > > Xoff3DForSorting;
    vector < pair< float_t, float_t > > Yoff3DForSorting;
    vector < pair< float_t, float_t > > Goodness3DForSorting;
    vector < float_t > MCE;
    
    vector < pair < float_t, float_t > > Sel3Dreturn;
    vector < pair < float_t, float_t > > sigmaL3Dreturn;
    vector < pair < float_t, float_t > > sigmaT3Dreturn;
    vector < pair < float_t, float_t > > Nc3Dreturn;
    vector < pair < float_t, float_t > > RWidth3Dreturn;
    vector < pair < float_t, float_t > > Saz3Dreturn;
    vector < pair < float_t, float_t > > Smax3Dreturn;
    vector < pair < float_t, float_t > > Depth3Dreturn;
    vector < pair < float_t, float_t > > Xcore3Dreturn;
    vector < pair < float_t, float_t > > Ycore3Dreturn;
    vector < pair < float_t, float_t > > Xoff3Dreturn;
    vector < pair < float_t, float_t > > Yoff3Dreturn;
    vector < pair < float_t, float_t > > Goodness3Dreturn;
    
    vector < pair < float_t, float_t > > ErrorSel3DForSorting;
    vector < pair < float_t, float_t > > ErrorsigmaL3DForSorting;
    vector < pair < float_t, float_t > > ErrorsigmaT3DForSorting;
    vector < pair < float_t, float_t > > ErrorNc3DForSorting;
    vector < pair < float_t, float_t > > ErrorRWidth3DForSorting;
    vector < pair < float_t, float_t > > ErrorSaz3DForSorting;
    vector < pair < float_t, float_t > > ErrorSmax3DForSorting;
    vector < pair < float_t, float_t > > ErrorXcore3DForSorting;
    vector < pair < float_t, float_t > > ErrorYcore3DForSorting;
    
    vector < int > Averagingvector;
    
    char SIM[5];
    if( strcmp( SIMTYPE, "CARE" ) == 0 )
    {
        sprintf( SIM, "1200" );
    }
    else
    {
        sprintf( SIM, "6500" );
    }
    
    
    for( int ideg = 0; ideg < 10; ideg += 1 ) 	//10 0-9
    {
        cout << " WintSum: " << SUMWIN << " Deg: " << Deg[ ideg ] << " Noice: " << NOICE << endl;
        
        char Inname[150];
        sprintf( Inname, "%s/ze%sdeg_offset0.5deg_NSB%sMHz_9%s%s.root", INOUT, Deg[ ideg ], NOICE, TELSETUP, SIM );
        
        bool anka = true;
        TFile Zom( Inname );
        char NotWorkName[20000];
        
        if( Zom.IsZombie() || Zom.TestBit( TFile::kRecovered ) )
        {
            sprintf( NotWorkName, "%s /%s/%s/%s", NotIncluded, SUMWIN, Deg[ ideg ], NOICE );
            strcpy( NotIncluded, NotWorkName );
            cout << "Because file IsZombie " << Inname << " won't be used"  << endl;
        }
        else
        {
            TFile* f1 = new TFile( Inname );
            anka = ( TTree* )f1->Get( "showerpars" );
            if( anka == 0 )
            {
                sprintf( NotWorkName, "%s %s/%s/%s", NotIncluded, SUMWIN, Deg[ ideg ], NOICE );
                strcpy( NotIncluded, NotWorkName );
                sprintf( SizeName, "%s %s/%s/%s/%d", IncludedSize, SUMWIN, Deg[ ideg ], NOICE, 0 );
                strcpy( IncludedSize, SizeName );
                cout << "Because of missing branches " << Inname << " won't be used" << endl;
            }
            else
            {
                cout << "Name of outputfile: " << Outname << endl;
                cout << "Have opened: " << Inname << endl;
                abba = 1;
                
                /// data ////
                TTree* s1 = ( TTree* )f1->Get( "model3Dpars" );
                
                TTree* s2 = ( TTree* )f1->Get( "showerpars" );
                cout << "Trees have been loaded." << endl;
                
                /// set branches /////////
                s1->SetBranchAddress( "Sel3D", 		&fSel3D );
                s1->SetBranchAddress( "sigmaL3D", 	&fsigmaL3D );
                s1->SetBranchAddress( "sigmaT3D", 	&fsigmaT3D );
                s1->SetBranchAddress( "Nc3D", 		&fNc3D );
                s1->SetBranchAddress( "RWidth3D",	&fRWidth3D );
                s1->SetBranchAddress( "Saz3D",		&fSaz3D );
                s1->SetBranchAddress( "Smax3D", 	&fSmax3D );
                s1->SetBranchAddress( "Depth3D", 	&fDepth3D );
                s1->SetBranchAddress( "Xcore3D", 	&fXcore3D );
                s1->SetBranchAddress( "Ycore3D", 	&fYcore3D );
                s1->SetBranchAddress( "Xoff3D", 	&fXoff3D );
                s1->SetBranchAddress( "Yoff3D", 	&fYoff3D );
                s1->SetBranchAddress( "Goodness3D", 	&fGoodness3D );
                s2->SetBranchAddress( "MCe0", 		&fMCe0 );
                
                s1->SetBranchAddress( "ErrorSel3D", 	&fErrorSel3D );
                s1->SetBranchAddress( "ErrorsigmaL3D", 	&fErrorsigmaL3D );
                s1->SetBranchAddress( "ErrorsigmaT3D", 	&fErrorsigmaT3D );
                s1->SetBranchAddress( "ErrorNc3D", 	&fErrorNc3D );
                s1->SetBranchAddress( "ErrRWidth3D", 	&fErrorRWidth3D );
                s1->SetBranchAddress( "ErrorSaz3D",	&fErrorSaz3D );
                s1->SetBranchAddress( "ErrorSmax3D", 	&fErrorSmax3D );
                s1->SetBranchAddress( "ErrorXcore3D", 	&fErrorXcore3D );
                s1->SetBranchAddress( "ErrorYcore3D", 	&fErrorYcore3D );
                
                nentries = s1->GetEntries();
                //cout << "nentries are: " << nentries <<  " and addnenties are " << addnentries << endl;
                
                int inputlength = 0;
                for( int newcount2 = 0; newcount2 < nentries; newcount2 += 1 )
                {
                    s1->GetEntry( newcount2 );	//Bad values will be selected against here
                    s2->GetEntry( newcount2 );
                    
                    if(
                        fMCe0 		> 0.03 &&  abs( sqrt( fXcore3D * fXcore3D + fYcore3D * fYcore3D ) ) <= 500/*178*/ &&
                        fNc3D	 	> 0.0 &&
                        fSel3D	 	> 0.0 && fSel3D > deg[ ideg ] - 2.5 && fSel3D < deg[ ideg ] + 2.5 &&
                        fSmax3D 	> 0.0 &&
                        fsigmaL3D 	> 0.0 &&
                        fsigmaT3D 	> 0.0 &&
                        fRWidth3D 	> 0.0 &&
                        fErrorNc3D 	> 0.0 &&
                        fErrorSel3D 	> 0.0 &&
                        fErrorSaz3D 	> 0.0 &&
                        fErrorXcore3D 	> 0.0 &&
                        fErrorYcore3D 	> 0.0 &&
                        fErrorSmax3D 	> 0.0 &&
                        fErrorsigmaL3D	> 0.0 &&
                        fErrorsigmaT3D	> 0.0 &&
                        fErrorRWidth3D	> 0.0 &&
                        abs( fXoff3D )	< 100 &&
                        
                        abs( fErrorRWidth3D / fRWidth3D )	<= Prospars[ ideg ][0] &&
                        abs( fErrorSel3D / fSel3D )		<= Prospars[ ideg ][1] &&
                        abs( fErrorSaz3D / fSaz3D )		<= Prospars[ ideg ][2] &&
                        abs( fErrorXcore3D / fXcore3D )	<= Prospars[ ideg ][3] &&
                        abs( fErrorYcore3D / fYcore3D )	<= Prospars[ ideg ][4] &&
                        abs( fErrorSmax3D / fSmax3D )	<= Prospars[ ideg ][5] &&
                        abs( fErrorsigmaL3D / fsigmaL3D )	<= Prospars[ ideg ][6] &&
                        abs( fErrorsigmaT3D / fsigmaT3D )	<= Prospars[ ideg ][7] &&
                        abs( fErrorNc3D / fNc3D )		<= Prospars[ ideg ][8] &&
                        
                        !TMath::IsNaN( fMCe0 ) &&
                        !TMath::IsNaN( fSel3D ) &&
                        !TMath::IsNaN( fsigmaL3D ) &&
                        !TMath::IsNaN( fsigmaT3D ) &&
                        !TMath::IsNaN( fRWidth3D ) &&
                        !TMath::IsNaN( fSaz3D ) &&
                        !TMath::IsNaN( fSmax3D ) &&
                        !TMath::IsNaN( fDepth3D ) &&
                        !TMath::IsNaN( fXcore3D ) &&
                        !TMath::IsNaN( fYcore3D ) &&
                        !TMath::IsNaN( fXoff3D ) &&
                        !TMath::IsNaN( fYoff3D ) &&
                        !TMath::IsNaN( fGoodness3D ) &&
                        !TMath::IsNaN( fErrorsigmaL3D ) &&
                        !TMath::IsNaN( fErrorsigmaT3D ) &&
                        !TMath::IsNaN( fErrorRWidth3D ) &&
                        !TMath::IsNaN( fErrorSaz3D ) &&
                        !TMath::IsNaN( fErrorSmax3D ) &&
                        !TMath::IsNaN( fErrorXcore3D ) &&
                        !TMath::IsNaN( fErrorYcore3D ) &&
                        !TMath::IsNaN( fErrorSel3D ) )
                    {
                    
                        //HERE the creation of average sets begins
                        if( inputlength == 0 )
                        {
                            Sel3DForSorting.push_back(	pair < float_t, int > (	fSel3D, 	addnentries + inputlength ) ); // 0-90 with  res 0.1, units deg
                            sigmaL3DForSorting.push_back(	pair < float_t, int > (	fsigmaL3D, 	addnentries + inputlength ) ); // normal 3,  res 0.1-0.5, units km
                            sigmaT3DForSorting.push_back(	pair < float_t, int > (	fsigmaT3D, 	addnentries + inputlength ) ); // normal 15, res 0.1-0.5, units m
                            Nc3DForSorting.push_back(	pair < float_t, int > (	fNc3D, 		addnentries + inputlength ) ); // normal 15, res 0.1-0.5,ln(x) e^x = #
                            RWidth3DForSorting.push_back(	pair < float_t, int > (	fRWidth3D, 	addnentries + inputlength ) ); // normal 3,  res 0.1-0.5, units m
                            Saz3DForSorting.push_back(	pair < float_t, int > (	fSaz3D, 	addnentries + inputlength ) ); // 0-360 with res 0.1, units deg
                            Smax3DForSorting.push_back(	pair < float_t, int > (	fSmax3D, 	addnentries + inputlength ) ); // normal 3,  res 0.1-0.5, units km
                            Depth3DForSorting.push_back(	pair < float_t, int > (	fDepth3D, 	addnentries + inputlength ) ); // normal 250-300, res ?, units optical depth
                            Xcore3DForSorting.push_back(	pair < float_t, int > (	fXcore3D, 	addnentries + inputlength ) ); // -178 - +178,    res 0.1-0.5, units m
                            Ycore3DForSorting.push_back(	pair < float_t, int > (	fYcore3D, 	addnentries + inputlength ) ); // -178 - +178,    res 0.1-0.5, units m
                            Xoff3DForSorting.push_back(	pair < float_t, int > (	fXoff3D, 	addnentries + inputlength ) ); // normal 3, res 0.1, units deg
                            Yoff3DForSorting.push_back(	pair < float_t, int > (	fYoff3D, 	addnentries + inputlength ) ); // normal 3, res 0.1, units deg
                            Goodness3DForSorting.push_back(	pair < float_t, int > (	fGoodness3D, 	addnentries + inputlength ) ); // ?
                            
                            ErrorSel3DForSorting.push_back(	pair < float_t, float_t > (	fSel3D, 	fErrorSel3D ) );
                            ErrorsigmaL3DForSorting.push_back(	pair < float_t, float_t > (	fsigmaL3D, 	fErrorsigmaL3D ) );
                            ErrorsigmaT3DForSorting.push_back(	pair < float_t, float_t > (	fsigmaT3D, 	fErrorsigmaT3D ) );
                            ErrorNc3DForSorting.push_back(	pair < float_t, float_t > (	fNc3D, 		fErrorNc3D ) );
                            ErrorRWidth3DForSorting.push_back(	pair < float_t, float_t > (	fRWidth3D, 	fErrorRWidth3D ) );
                            ErrorSaz3DForSorting.push_back(	pair < float_t, float_t > (	fSaz3D, 	fErrorSaz3D ) );
                            ErrorSmax3DForSorting.push_back(	pair < float_t, float_t > (	fSmax3D, 	fErrorSmax3D ) );
                            ErrorXcore3DForSorting.push_back(	pair < float_t, float_t > (	fXcore3D, 	fErrorXcore3D ) );
                            ErrorYcore3DForSorting.push_back(	pair < float_t, float_t > (	fYcore3D, 	fErrorYcore3D ) );
                            
                            Averagingvector.push_back( 1 );
                            
                            inputlength++;
                            MCE.push_back( fMCe0 );
                            
                        }
                        else
                        {
                            int istherecounter = 0;
                            for( int reloop = addnentries; reloop < addnentries + inputlength; reloop++ )
                            {
                                /* inside of the range, average */
                                if( istherecounter == 0 &&
                                
                                        fsigmaL3D 	> sigmaL3DForSorting.at( reloop ).first -	0.25	&& //0.2-0.21
                                        fsigmaL3D 	< sigmaL3DForSorting.at( reloop ).first +	0.25	&& //0.2-0-21
                                        fsigmaT3D 	> sigmaT3DForSorting.at( reloop ).first -	0.25	&& //0.4-1
                                        fsigmaT3D 	< sigmaT3DForSorting.at( reloop ).first +	0.25	&& //0.4-1
                                        fSmax3D 	> Smax3DForSorting.at( reloop ).first -	0.25	&& //0.21-0.8
                                        fSmax3D 	< Smax3DForSorting.at( reloop ).first +	0.25	&& //0.21-0.8
                                        fNc3D > Nc3DForSorting.at( reloop ).first -			0.25	&& //0.06
                                        fNc3D < Nc3DForSorting.at( reloop ).first +		0.25	&& //0.06
                                        sqrt( fXcore3D * fXcore3D + fYcore3D * fYcore3D ) > sqrt( Xcore3DForSorting.at( reloop ).first * Xcore3DForSorting.at( reloop ).first + Ycore3DForSorting.at( reloop ).first * Ycore3DForSorting.at( reloop ).first ) - 0.5 && //0.5
                                        sqrt( fXcore3D * fXcore3D + fYcore3D * fYcore3D ) < sqrt( Xcore3DForSorting.at( reloop ).first * Xcore3DForSorting.at( reloop ).first + Ycore3DForSorting.at( reloop ).first * Ycore3DForSorting.at( reloop ).first ) + 0.5 && //0.5
                                        fYcore3D / ( Ycore3DForSorting.at( reloop ).first + 0.000000001 ) >= 0.0 &&
                                        fXcore3D / ( Xcore3DForSorting.at( reloop ).first + 0.000000001 ) >= 0.0
                                  )
                                {
                                
                                    //HERE I need to average the two values
                                    Sel3DForSorting.at( reloop ).first	= ( Sel3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSel3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    sigmaL3DForSorting.at( reloop ).first 	= ( sigmaL3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fsigmaL3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    sigmaT3DForSorting.at( reloop ).first 	= ( sigmaT3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fsigmaT3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Nc3DForSorting.at( reloop ).first 	= ( Nc3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fNc3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    RWidth3DForSorting.at( reloop ).first 	= ( RWidth3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fRWidth3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Saz3DForSorting.at( reloop ).first 	= ( Saz3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSaz3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Smax3DForSorting.at( reloop ).first 	= ( Smax3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSmax3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Depth3DForSorting.at( reloop ).first 	= ( Depth3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fDepth3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Xcore3DForSorting.at( reloop ).first 	= ( Xcore3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fXcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Ycore3DForSorting.at( reloop ).first 	= ( Ycore3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fYcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Xoff3DForSorting.at( reloop ).first 	= ( Xoff3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fXoff3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Yoff3DForSorting.at( reloop ).first 	= ( Yoff3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fYcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    Goodness3DForSorting.at( reloop ).first 	= ( Goodness3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fGoodness3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    
                                    ErrorSel3DForSorting.at( reloop ).first 	 = ( ErrorSel3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSel3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorsigmaL3DForSorting.at( reloop ).first = ( ErrorsigmaL3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fsigmaL3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorsigmaT3DForSorting.at( reloop ).first = ( ErrorsigmaT3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fsigmaT3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorNc3DForSorting.at( reloop ).first 	 = ( ErrorNc3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fErrorNc3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorRWidth3DForSorting.at( reloop ).first = ( ErrorRWidth3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fRWidth3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorSaz3DForSorting.at( reloop ).first 	 = ( ErrorSaz3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSaz3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorSmax3DForSorting.at( reloop ).first 	 = ( ErrorSmax3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fSmax3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorXcore3DForSorting.at( reloop ).first  = ( ErrorXcore3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fXcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorYcore3DForSorting.at( reloop ).first  = ( ErrorYcore3DForSorting.at( reloop ).first * Averagingvector.at( reloop ) + fYcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    
                                    ErrorSel3DForSorting.at( reloop ).second 	  = ( ErrorSel3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorSel3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorsigmaL3DForSorting.at( reloop ).second = ( ErrorsigmaL3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorsigmaL3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorsigmaT3DForSorting.at( reloop ).second = ( ErrorsigmaT3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorsigmaT3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorNc3DForSorting.at( reloop ).second 	  = ( ErrorNc3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorNc3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorRWidth3DForSorting.at( reloop ).second = ( ErrorRWidth3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorRWidth3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorSaz3DForSorting.at( reloop ).second 	  = ( ErrorSaz3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorSaz3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorSmax3DForSorting.at( reloop ).second   = ( ErrorSmax3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorSmax3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorXcore3DForSorting.at( reloop ).second  = ( ErrorXcore3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorXcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    ErrorYcore3DForSorting.at( reloop ).second  = ( ErrorYcore3DForSorting.at( reloop ).second * Averagingvector.at( reloop ) + fErrorYcore3D ) / ( Averagingvector.at( reloop ) + 1 );
                                    
                                    Averagingvector.at( reloop ) = Averagingvector.at( reloop ) + 1;
                                    istherecounter = reloop;
                                }
                            }
                            if( istherecounter == 0 )  /* outside of the range, just add*/
                            {
                                Sel3DForSorting.push_back(	pair < float_t, int > (	fSel3D, 	addnentries + inputlength ) );
                                sigmaL3DForSorting.push_back(	pair < float_t, int > (	fsigmaL3D, 	addnentries + inputlength ) );
                                sigmaT3DForSorting.push_back(	pair < float_t, int > (	fsigmaT3D, 	addnentries + inputlength ) );
                                Nc3DForSorting.push_back(	pair < float_t, int > (	fNc3D, 		addnentries + inputlength ) );
                                
                                RWidth3DForSorting.push_back(	pair < float_t, int > (	fRWidth3D, 	addnentries + inputlength ) );
                                Saz3DForSorting.push_back(	pair < float_t, int > (	fSaz3D, 	addnentries + inputlength ) );
                                Smax3DForSorting.push_back(	pair < float_t, int > (	fSmax3D, 	addnentries + inputlength ) );
                                Depth3DForSorting.push_back(	pair < float_t, int > (	fDepth3D, 	addnentries + inputlength ) );
                                Xcore3DForSorting.push_back(	pair < float_t, int > (	fXcore3D, 	addnentries + inputlength ) );
                                Ycore3DForSorting.push_back(	pair < float_t, int > (	fYcore3D, 	addnentries + inputlength ) );
                                Xoff3DForSorting.push_back(	pair < float_t, int > (	fXoff3D, 	addnentries + inputlength ) );
                                Yoff3DForSorting.push_back(	pair < float_t, int > (	fYoff3D, 	addnentries + inputlength ) );
                                Goodness3DForSorting.push_back(	pair < float_t, int > (	fGoodness3D, 	addnentries + inputlength ) );
                                
                                ErrorSel3DForSorting.push_back(	pair < float_t, float_t > (	fSel3D, 	fErrorSel3D ) );
                                ErrorsigmaL3DForSorting.push_back(	pair < float_t, float_t > (	fsigmaL3D, 	fErrorsigmaL3D ) );
                                ErrorsigmaT3DForSorting.push_back(	pair < float_t, float_t > (	fsigmaT3D, 	fErrorsigmaT3D ) );
                                ErrorNc3DForSorting.push_back(	pair < float_t, float_t > (	fNc3D, 		fErrorNc3D ) );
                                ErrorRWidth3DForSorting.push_back(	pair < float_t, float_t > (	fRWidth3D, 	fErrorRWidth3D ) );
                                ErrorSaz3DForSorting.push_back(	pair < float_t, float_t > (	fSaz3D, 	fErrorSaz3D ) );
                                ErrorSmax3DForSorting.push_back(	pair < float_t, float_t > (	fSmax3D, 	fErrorSmax3D ) );
                                ErrorXcore3DForSorting.push_back(	pair < float_t, float_t > (	fXcore3D, 	fErrorXcore3D ) );
                                ErrorYcore3DForSorting.push_back(	pair < float_t, float_t > (	fYcore3D, 	fErrorYcore3D ) );
                                
                                Averagingvector.push_back( 1 );
                                
                                MCE.push_back( fMCe0 );
                                inputlength++;
                            }
                        }
                    }
                }
                addnentries += inputlength;
                cout << inputlength << " of " << nentries << endl;
                sprintf( SizeName, "%s %s/%s/%d", IncludedSize, SUMWIN, NOICE, addnentries );
                strcpy( IncludedSize, SizeName );
            }
        }
    }
    
    //make a sorting loop
    //This should sort my vectors of pairs according to the first coloumn and have the energy follow that coloumn
    sort( Sel3DForSorting.begin(), 		Sel3DForSorting.end(), 		compare );	//NAH, sorting acc to the energy
    sort( sigmaL3DForSorting.begin(),	sigmaL3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( sigmaT3DForSorting.begin(),	sigmaT3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( Nc3DForSorting.begin(), 		Nc3DForSorting.end(), 		compare );	//NAH, sorting acc to the energy
    sort( RWidth3DForSorting.begin(),	RWidth3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( Smax3DForSorting.begin(),		Smax3DForSorting.end(),		compare );	//NAH, sorting acc to the energy
    sort( Saz3DForSorting.begin(), 		Saz3DForSorting.end(),		compare );	//NAH, sorting acc to the energy
    sort( Depth3DForSorting.begin(),	Depth3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( Xcore3DForSorting.begin(),	Xcore3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    sort( Ycore3DForSorting.begin(),	Ycore3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    sort( Xoff3DForSorting.begin(),		Xoff3DForSorting.end(),		compare );	//NAH, sorting acc to the energy
    sort( Yoff3DForSorting.begin(),		Yoff3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( Goodness3DForSorting.begin(),	Goodness3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    
    sort( ErrorSel3DForSorting.begin(), 	ErrorSel3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( ErrorsigmaL3DForSorting.begin(),	ErrorsigmaL3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( ErrorsigmaT3DForSorting.begin(),	ErrorsigmaT3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( ErrorNc3DForSorting.begin(), 	ErrorNc3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( ErrorRWidth3DForSorting.begin(),	ErrorRWidth3DForSorting.end(), 	compare );	//NAH, sorting acc to the energy
    sort( ErrorSaz3DForSorting.begin(), 	ErrorSaz3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    sort( ErrorSmax3DForSorting.begin(),	ErrorSmax3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    sort( ErrorXcore3DForSorting.begin(),	ErrorXcore3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    sort( ErrorYcore3DForSorting.begin(),	ErrorYcore3DForSorting.end(),	compare );	//NAH, sorting acc to the energy
    
    for( int sortback = 0; sortback < addnentries; sortback++ )
    {
        Sel3Dreturn.push_back(	pair < int, int > (	Sel3DForSorting.at(	sortback ).second, 	sortback ) );
        sigmaL3Dreturn.push_back(	pair < int, int > (	sigmaL3DForSorting.at(	sortback ).second, 	sortback ) );
        sigmaT3Dreturn.push_back(	pair < int, int > (	sigmaT3DForSorting.at(	sortback ).second, 	sortback ) );
        Nc3Dreturn.push_back(	pair < int, int > (	Nc3DForSorting.at(	sortback ).second,	sortback ) );
        RWidth3Dreturn.push_back(	pair < int, int > (	RWidth3DForSorting.at(	sortback ).second,	sortback ) );
        Saz3Dreturn.push_back(	pair < int, int > (	Saz3DForSorting.at(	sortback ).second,	sortback ) );
        Smax3Dreturn.push_back(	pair < int, int > (	Smax3DForSorting.at(	sortback ).second, 	sortback ) );
        Depth3Dreturn.push_back(	pair < int, int > (	Depth3DForSorting.at(	sortback ).second,	sortback ) );
        Xcore3Dreturn.push_back(	pair < int, int > (	Xcore3DForSorting.at(	sortback ).second, 	sortback ) );
        Ycore3Dreturn.push_back(	pair < int, int > (	Ycore3DForSorting.at(	sortback ).second,	sortback ) );
        Xoff3Dreturn.push_back(	pair < int, int > (	Xoff3DForSorting.at(	sortback ).second, 	sortback ) );
        Yoff3Dreturn.push_back(	pair < int, int > (	Yoff3DForSorting.at(	sortback ).second,	sortback ) );
        Goodness3Dreturn.push_back(	pair < int, int > (	Goodness3DForSorting.at( sortback ).second, 	sortback ) );
    }
    
    sort( Sel3Dreturn.begin(),	Sel3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( sigmaL3Dreturn.begin(),	sigmaL3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( sigmaT3Dreturn.begin(),	sigmaT3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Nc3Dreturn.begin(),	Nc3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( RWidth3Dreturn.begin(),	RWidth3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Saz3Dreturn.begin(),	Saz3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Smax3Dreturn.begin(),	Smax3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Depth3Dreturn.begin(),	Depth3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Xcore3Dreturn.begin(),	Xcore3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Ycore3Dreturn.begin(),	Ycore3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Xoff3Dreturn.begin(),	Xoff3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Yoff3Dreturn.begin(),	Yoff3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    sort( Goodness3Dreturn.begin(),	Goodness3Dreturn.end(),	compareInt );	//NAH, sorting acc to the energy
    cout << " " << endl;
    
    for( int gothroughagain = 0; gothroughagain < addnentries; gothroughagain++ )
    {
    
        f2Sel3D = 	Sel3DForSorting.at( gothroughagain ).first;
        f2sigmaL3D = 	sigmaL3DForSorting.at( gothroughagain ).first;
        f2sigmaT3D = 	sigmaT3DForSorting.at( gothroughagain ).first;
        f2Nc3D = 	Nc3DForSorting.at( gothroughagain ).first;
        f2RWidth3D = 	RWidth3DForSorting.at( gothroughagain ).first;
        f2Saz3D = 	Saz3DForSorting.at( gothroughagain ).first;
        f2Smax3D = 	Smax3DForSorting.at( gothroughagain ).first;
        f2Depth3D = 	Depth3DForSorting.at( gothroughagain ).first;
        f2Xcore3D = 	Xcore3DForSorting.at( gothroughagain ).first;
        f2Ycore3D = 	Ycore3DForSorting.at( gothroughagain ).first;
        f2Xoff3D = 	Xoff3DForSorting.at( gothroughagain ).first;
        f2Yoff3D = 	Yoff3DForSorting.at( gothroughagain ).first;
        f2Goodness3D = 	Goodness3DForSorting.at( gothroughagain ).first;
        f2MCe0 = 	MCE[			 gothroughagain];
        
        f2EnergySel3D = 	Sel3DForSorting.at(	gothroughagain ).second;
        f2EnergysigmaL3D = 	sigmaL3DForSorting.at(	gothroughagain ).second;
        f2EnergysigmaT3D = 	sigmaT3DForSorting.at(	gothroughagain ).second;
        f2EnergyNc3D = 		Nc3DForSorting.at(	gothroughagain ).second;
        f2EnergyRWidth3D = 	RWidth3DForSorting.at(	gothroughagain ).second;
        f2EnergySaz3D = 	Saz3DForSorting.at(	gothroughagain ).second;
        f2EnergySmax3D = 	Smax3DForSorting.at(	gothroughagain ).second;
        f2EnergyDepth3D = 	Depth3DForSorting.at(	gothroughagain ).second;
        f2EnergyXcore3D = 	Xcore3DForSorting.at(	gothroughagain ).second;
        f2EnergyYcore3D = 	Ycore3DForSorting.at(	gothroughagain ).second;
        f2EnergyXoff3D = 	Xoff3DForSorting.at(	gothroughagain ).second;
        f2EnergyYoff3D = 	Yoff3DForSorting.at(	gothroughagain ).second;
        f2EnergyGoodness3D = 	Goodness3DForSorting.at( gothroughagain ).second;
        
        f2ReturnSel3D = 	Sel3Dreturn.at( gothroughagain ).first; //NOTE I use to have them as .second, but af2ter testing appered I should have them as .first
        f2ReturnsigmaL3D = 	sigmaL3Dreturn.at(	gothroughagain ).first;
        f2ReturnsigmaT3D = 	sigmaT3Dreturn.at(	gothroughagain ).first;
        f2ReturnNc3D = 		Nc3Dreturn.at( gothroughagain ).first;
        f2ReturnRWidth3D = 	RWidth3Dreturn.at(	gothroughagain ).first;
        f2ReturnSaz3D = 	Saz3Dreturn.at( gothroughagain ).first;
        f2ReturnSmax3D = 	Smax3Dreturn.at(	gothroughagain ).first;
        f2ReturnDepth3D = 	Depth3Dreturn.at(	gothroughagain ).first;
        f2ReturnXcore3D = 	Xcore3Dreturn.at(	gothroughagain ).first;
        f2ReturnYcore3D = 	Ycore3Dreturn.at(	gothroughagain ).first;
        f2ReturnXoff3D = 	Xoff3Dreturn.at(	gothroughagain ).first;
        f2ReturnYoff3D = 	Yoff3Dreturn.at(	gothroughagain ).first;
        f2ReturnGoodness3D = 	Goodness3Dreturn.at(	gothroughagain ).first;
        
        
        f2PlasementSel3D = 	Sel3Dreturn.at( gothroughagain ).second;
        f2PlasementsigmaL3D = 	sigmaL3Dreturn.at(	gothroughagain ).second;
        f2PlasementsigmaT3D = 	sigmaT3Dreturn.at(	gothroughagain ).second;
        f2PlasementNc3D = 	Nc3Dreturn.at( gothroughagain ).second;
        f2PlasementRWidth3D = 	RWidth3Dreturn.at(	gothroughagain ).second;
        f2PlasementSaz3D = 	Saz3Dreturn.at( gothroughagain ).second;
        f2PlasementSmax3D = 	Smax3Dreturn.at(	gothroughagain ).second;
        f2PlasementDepth3D = 	Depth3Dreturn.at(	gothroughagain ).second;
        f2PlasementXcore3D = 	Xcore3Dreturn.at(	gothroughagain ).second; //NOTE I use to have them as .second, but af2ter testing appered I should have them as .second
        f2PlasementYcore3D = 	Ycore3Dreturn.at(	gothroughagain ).second;
        f2PlasementXoff3D = 	Xoff3Dreturn.at(	gothroughagain ).second;
        f2PlasementYoff3D = 	Yoff3Dreturn.at(	gothroughagain ).second;
        f2PlasementGoodness3D = Goodness3Dreturn.at(	gothroughagain ).second;
        
        f2ErrorSel3D = 		ErrorSel3DForSorting.at(	gothroughagain ).second;
        f2ErrorsigmaL3D = 	ErrorsigmaL3DForSorting.at(	gothroughagain ).second;
        f2ErrorsigmaT3D = 	ErrorsigmaT3DForSorting.at(	gothroughagain ).second;
        f2ErrorNc3D = 		ErrorNc3DForSorting.at(	gothroughagain ).second;
        f2ErrorRWidth3D = 	ErrorRWidth3DForSorting.at(	gothroughagain ).second;
        f2ErrorSaz3D = 		ErrorSaz3DForSorting.at(	gothroughagain ).second;
        f2ErrorSmax3D = 	ErrorSmax3DForSorting.at(	gothroughagain ).second;
        f2ErrorXcore3D = 	ErrorXcore3DForSorting.at(	gothroughagain ).second;
        f2ErrorYcore3D = 	ErrorYcore3DForSorting.at(	gothroughagain ).second;
        
        ModelPars3D->Fill();
    }
    
    delete[] SIMTYPE;
    delete[] TELSETUP;
    delete[] SUMWIN;
    delete[] NOICE;
    delete[] INOUT;
    
    if( abba == 1 )
    {
        cout << "Writing and closing file" << endl;
        newfile->Write();	//NAH, write the outputfile
        newfile->Close();	//NAH, properly close the outputfile
    }
    
    cout << NotIncluded << endl;
    cout << endl;
    cout << IncludedSize << endl;
}
