
#include <string>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"

//For the InputOutputfile
Float_t InOutSel3D;		// elevation of shower direction (deg)
Float_t InOutsigmaL3D;		// longitudinal (3D-length)
Float_t InOutsigmaT3D;		// transverse (3D-width)
Float_t InOutNc3D;     		// total number of Cherenkov photons emitted by the shower
Float_t InOutRWidth3D;
Float_t InOutSaz3D;		// azimuth of shower direction (deg)
Float_t InOutSmax3D;		// height of shower maximum (along the shower axis)
Float_t InOutDepth3D;
Float_t InOutXcore3D;		// shower core in ground coordinates
Float_t InOutYcore3D;		// shower core in ground coordinates
Float_t InOutGoodness3D;	// Goodness i.e how much a shower is gamma like
Float_t InOutXoff3D;		// Xcore offset 3D
Float_t InOutYoff3D;		// Ycore offset 3D
Float_t InOutRedSigT3D;		// RedSigT3D
Float_t InOutRedSigL3D;		// RedSigL3D
Float_t InOutMCe0;

Float_t InOutSel3DErr;		// error values for the same
Float_t InOutsigmaL3DErr;
Float_t InOutsigmaT3DErr;
Float_t InOutNc3DErr;
Float_t InOutRWidth3DErr;
Float_t InOutSaz3DErr;
Float_t InOutSmax3DErr;
Float_t InOutXcore3DErr;
Float_t InOutYcore3DErr;

Double_t InOutpedvars1;		// pedvars to get right noise
Double_t InOutpedvars2;
Double_t InOutpedvars3;
Double_t InOutpedvars4;

Double_t TelElevation;		// telescope elevation to get correct biascorrection


//For the Referencefiles
Float_t RefSel3D;	// elevation of shower direction (deg)
Float_t RefsigmaL3D;	// longitudinal (3D-length)
Float_t RefsigmaT3D; 	// transverse (3D-width)
Float_t RefNc3D;     	// total number of Cherenkov photons emitted by the shower
Float_t RefRWidth3D;	// redused width 3D
Float_t RefSaz3D;	// azimuth of shower direction (deg)
Float_t RefSmax3D;	// height of shower maximum (along the shower axis)
Float_t RefDepth3D;	// shower depth
Float_t RefXcore3D;	// shower core in ground coordinates
Float_t RefYcore3D;	// shower core in ground coordinates
Float_t RefGoodness3D;	// Goodness i.e how much a shower is gamma like
Float_t RefXoff3D;	// Xcore offset 3D
Float_t RefYoff3D;	// Ycore offset 3D
Float_t RefEnergy3D;	// Energy3D

Float_t RefMCe0;

Float_t RefEnergySel3D;	// position of elevation of shower direction and so on in the reference dataset
Float_t RefEnergysigmaL3D;
Float_t RefEnergysigmaT3D;
Float_t RefEnergyNc3D;
Float_t RefEnergyRWidth3D;
Float_t RefEnergySaz3D;
Float_t RefEnergySmax3D;
Float_t RefEnergyDepth3D;

Float_t RefReturnSel3D;		// To be able to jump between the different sorted values for all relevant reference parameters
Float_t RefReturnsigmaL3D;
Float_t RefReturnsigmaT3D;
Float_t RefReturnNc3D;
Float_t RefReturnRWidth3D;
Float_t RefReturnSaz3D;
Float_t RefReturnSmax3D;
Float_t RefReturnDepth3D;

Float_t RefPlasementSel3D;		// again to be able to jump between the relevatn reference parameters
Float_t RefPlasementsigmaL3D;
Float_t RefPlasementsigmaT3D;
Float_t RefPlasementNc3D;
Float_t RefPlasementRWidth3D;
Float_t RefPlasementSaz3D;
Float_t RefPlasementSmax3D;
Float_t RefPlasementDepth3D;

Float_t RefErrorSel3D;		// Error values for the reference values
Float_t RefErrorsigmaL3D;
Float_t RefErrorsigmaT3D;
Float_t RefErrorNc3D;
Float_t RefErrorRWidth3D;
Float_t RefErrorSaz3D;
Float_t RefErrorSmax3D;
Float_t RefErrorDepth3D;
Float_t RefErrorXcore3D;
Float_t RefErrorYcore3D;

Float_t EL;	// higher and lower values to do interpolation to get better energy estimation
Float_t EH;
Float_t selL;
Float_t selH;
Float_t sigLL;
Float_t sigTL;
Float_t NcL;
Float_t RWL;
Float_t SML;
Float_t DL;
Float_t CEH;
Float_t sigLH;
Float_t sigTH;
Float_t NcH;
Float_t RWH;
Float_t SMH;
Float_t DH;

// things that does not need to be in the main body of code
int SIMend;
const char* Noice[10];
const char* Deg[10] = {"00", "20", "30", "35", "40", "45", "50", "55", "60", "65"};	// These are Zenith angels
int intDeg[10] = {0, 20, 30, 35, 40, 45, 50, 55, 60, 65};					// Same angles for another purpouse
int Oldnoiceloop = 42;	// Only needed if run locally

using namespace std;

bool sortingfloat( const pair<Float_t, Float_t>& i, const pair<Float_t, Float_t>& j )
{
    return i.first < j.first;
}


int main( int argc, char* argv[] )
{

    // Getting the input information //
    string Inputfilenamelist;
    string indirname;
    string telver;
    string wintsum;
    string REFSIMULATIONS;
    string BIASCORRLIST;
    
    if( argc > 8 )
    {
        cout << "Usage: ./Energy3d [3D model analysed .root file list (should be runnumbers)] [Directory where the 3D-modeled files are and where the output will go] [list of tel vertions] [list of winter summer]  " << endl;
        cout << "Example for real data: ./Energy3d listdir/RunList.txt $VERITAS_USER_DATA_DIR/analysis/CARE/v500-dev08/CARE_June1425/ 6 21 CARE energy3d_biascorrection.txt" << endl;
        cout << " OR " << endl;
        cout << "Example for sim data ./Energy3d listdir/RunList.txt $VERITAS_USER_DATA_DIR/analysis/energy3d/simanalysis/GRISU/ 5 22 GRISU energy3d_biascorrection.txt" << endl;
        return 0;
    }
    else
    {
        Inputfilenamelist = argv[1];
        indirname	  = argv[2];
        telver		  = argv[3];
        wintsum		  = argv[4];
        REFSIMULATIONS 	  = argv[5];
        BIASCORRLIST	  = argv[6];
    }
    
    int DOBIASCORRECTION = 1;	// If 0, no biascorrection is done. If 1, biascorrection is done
    
    // The Inputfile name //
    vector< string > InNameLista;
    string line;
    int listlength = 0;
    ifstream inputa( Inputfilenamelist.c_str() );
    for( string line; getline( inputa, line ); )
    {
        InNameLista.push_back( line.c_str() ) ;
        listlength++;
    }
    
    // Getting the biascorrection values from the inputlist
    vector< string > INBIASCORR;
    string rader;
    int biascorrlength = 0;
    ifstream inputbias( BIASCORRLIST.c_str() );
    for( string rader; getline( inputbias, rader ); )
    {
        INBIASCORR.push_back( rader.c_str() ) ;
        biascorrlength++;
    }
    
    // The inputdirectory //
    char* InDir = new char[ indirname.size() + 1 ];
    InDir[ indirname.size() ] = 0;
    memcpy( InDir, indirname.c_str(), indirname.size() );
    
    // Telescope and wint/sum information
    const char* VERSION = telver.c_str();
    const char* WinterSummer = wintsum.c_str();
    
    // Refsimulation information //
    if( REFSIMULATIONS == "CARE" )
    {
        Noice[0] = "50";
        Noice[1] = "80";
        Noice[2] = "120";
        Noice[3] = "170";
        Noice[4] = "230";
        Noice[5] = "290";
        Noice[6] = "370";
        Noice[7] = "450";
        Noice[8] = "0";
        Noice[9] = "0";
        SIMend = 1200;
    }
    else  	//GRISU
    {
        Noice[0] = "075";
        Noice[1] = "100";
        Noice[2] = "150";
        Noice[3] = "200";
        Noice[4] = "250";
        Noice[5] = "325";
        Noice[6] = "425";
        Noice[7] = "550";
        Noice[8] = "750";
        Noice[9] = "1000";
        SIMend = 6500;
    }
    
    // Getting the data from the input root file //
    for( int runnumberLength = 0; runnumberLength < listlength; runnumberLength++ )
    {
        const char* theinfile = InNameLista[ runnumberLength ].c_str();
        char RunnumberFile[150];
        
        sprintf( RunnumberFile, "%s/%s.root", InDir, theinfile );
        cout << "Opening " << RunnumberFile << endl;
        
        TFile* InTfile = new TFile( RunnumberFile , "Update" );
        
        TTree* InTTrees = ( TTree* )InTfile->Get( "showerpars" );
        TTree* InTTreem = ( TTree* )InTfile->Get( "model3Dpars" );
        TTree* InTTreep = ( TTree* )InTfile->Get( "pointingData" );
        TTree* t1 = ( TTree* )InTfile->Get( "Tel_1/calib_1" );
        TTree* t2 = ( TTree* )InTfile->Get( "Tel_2/calib_2" );
        TTree* t3 = ( TTree* )InTfile->Get( "Tel_3/calib_3" );
        TTree* t4 = ( TTree* )InTfile->Get( "Tel_4/calib_4" );
        
        // Setting addresses to the indata //
        InTTreem->SetBranchAddress( "Sel3D", 		&InOutSel3D );
        InTTreem->SetBranchAddress( "sigmaL3D", 	&InOutsigmaL3D );
        InTTreem->SetBranchAddress( "sigmaT3D", 	&InOutsigmaT3D );
        InTTreem->SetBranchAddress( "Nc3D", 		&InOutNc3D );
        InTTreem->SetBranchAddress( "RWidth3D", 	&InOutRWidth3D );
        InTTreem->SetBranchAddress( "Saz3D",		&InOutSaz3D );
        InTTreem->SetBranchAddress( "Smax3D", 		&InOutSmax3D );
        InTTreem->SetBranchAddress( "Depth3D", 		&InOutDepth3D );
        InTTreem->SetBranchAddress( "Xcore3D", 		&InOutXcore3D );
        InTTreem->SetBranchAddress( "Ycore3D", 		&InOutYcore3D );
        InTTreem->SetBranchAddress( "Xoff3D", 		&InOutXoff3D );
        InTTreem->SetBranchAddress( "Yoff3D", 		&InOutYoff3D );
        InTTreem->SetBranchAddress( "Goodness3D", 	&InOutGoodness3D );
        
        InTTrees->SetBranchAddress( "MCe0", 		&InOutMCe0 );
        t1->SetBranchAddress( "pedvar", 		&InOutpedvars1 );
        t2->SetBranchAddress( "pedvar", 		&InOutpedvars2 );
        t3->SetBranchAddress( "pedvar", 		&InOutpedvars3 );
        t4->SetBranchAddress( "pedvar", 		&InOutpedvars4 );
        
        InTTreem->SetBranchAddress( "ErrorSel3D", 	&InOutSel3DErr );
        InTTreem->SetBranchAddress( "ErrorsigmaL3D", 	&InOutsigmaL3DErr );
        InTTreem->SetBranchAddress( "ErrorsigmaT3D", 	&InOutsigmaT3DErr );
        InTTreem->SetBranchAddress( "ErrorNc3D", 	&InOutNc3DErr );
        InTTreem->SetBranchAddress( "ErrRWidth3D", 	&InOutRWidth3DErr );
        InTTreem->SetBranchAddress( "ErrorSaz3D",	&InOutSaz3DErr );
        InTTreem->SetBranchAddress( "ErrorSmax3D", 	&InOutSmax3DErr );
        InTTreem->SetBranchAddress( "ErrorXcore3D", 	&InOutXcore3DErr );
        InTTreem->SetBranchAddress( "ErrorYcore3D", 	&InOutYcore3DErr );
        
        InTTreep->SetBranchAddress( "TelElevation", 	&TelElevation );
        
        int InEntries = InTTrees->GetEntries();
        int t1Entries = t1->GetEntries();
        int t2Entries = t2->GetEntries();
        int t3Entries = t3->GetEntries();
        int t4Entries = t4->GetEntries();
        
        // Preparation for filling vectors with inputdata //
        vector < pair< Float_t, Float_t > > ForSort3D;
        vector < vector < Float_t > > Pars3DForSort;
        vector < Float_t > ForSorting3DParameters;
        
        // For faster search, getting the Nc3D, then sorting Nc3D according to its size
        for( int InSort = 0; InSort < InEntries; InSort++ )
        {
            InTTreem->GetEntry( InSort );
            ForSort3D.push_back( pair < Float_t, int > ( InOutNc3D, InSort ) );
        }
        sort( ForSort3D.begin(), ForSort3D.end(), sortingfloat );	//I am using this later, it is correct about speed increase, but it might be faster if I get the parameters according to Nc3D as I say I should...
        
        //float_t AverageElevation = 0;
        
        //float_t AverageArray[10] = {0,0,0,0,0,0,0,0,0,0};
        //float_t AverageWeightArray[10] = {1,1,1,1,1,1,1,1,1,1};
        //vector < const char > AveragingLocation;
        vector < int > AveragingLocation2;
        float_t Zenithangle;
        // Getting input according to sorted Nc3D //
        for( int InS = 0; InS < InEntries; InS++ )
        {
            //InTTreem->GetEntry( ForSort3D[InS].second );	//Whith this I should change all Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ] to Pars3DForSort[ EnergySearchIn ] or alternatives to EnergySearchIn
            //InTTrees->GetEntry( ForSort3D[InS].second );
            //InTTreep->GetEntry( ForSort3D[InS].second );
            InTTreem->GetEntry( InS );
            InTTrees->GetEntry( InS );
            InTTreep->GetEntry( InS );
            ForSorting3DParameters.push_back(	InOutSel3D );		//0
            ForSorting3DParameters.push_back(	InOutsigmaL3D );		//1
            ForSorting3DParameters.push_back(	InOutsigmaT3D ); 	//2
            ForSorting3DParameters.push_back(	InOutNc3D );		//3
            ForSorting3DParameters.push_back(	InOutRWidth3D );		//4
            ForSorting3DParameters.push_back(	InOutSaz3D );		//5
            ForSorting3DParameters.push_back(	InOutSmax3D );		//6
            ForSorting3DParameters.push_back(	InOutDepth3D );		//7
            ForSorting3DParameters.push_back(	InOutXcore3D );		//8
            ForSorting3DParameters.push_back(	InOutYcore3D );		//9
            ForSorting3DParameters.push_back(	InOutGoodness3D );	//10
            ForSorting3DParameters.push_back(	InS );			//11
            ForSorting3DParameters.push_back(	InOutSel3DErr );		//12
            ForSorting3DParameters.push_back(	InOutsigmaL3DErr );	//13
            ForSorting3DParameters.push_back(	InOutsigmaT3DErr ); 	//14
            ForSorting3DParameters.push_back(	InOutNc3DErr );		//15
            ForSorting3DParameters.push_back(	InOutRWidth3DErr );	//16
            ForSorting3DParameters.push_back(	InOutSaz3DErr );		//17
            ForSorting3DParameters.push_back(	InOutSmax3DErr );	//18
            ForSorting3DParameters.push_back(	InOutXcore3DErr );	//19
            ForSorting3DParameters.push_back(	InOutYcore3DErr );	//20
            ForSorting3DParameters.push_back(	InOutMCe0 );		//21
            Zenithangle = 90 - TelElevation;
            ForSorting3DParameters.push_back(	Zenithangle );		//22
            if( Zenithangle <= 10.0 ) 	//00
            {
                //AverageArray[0] += TelElevation;
                //AverageWeightArray[0] += 1;
                //AveragingLocation.push_back("00");
                AveragingLocation2.push_back( 0 );
            }
            else if( Zenithangle > 10.0 && Zenithangle <= 25.0 )  //20
            {
                //AverageArray[1] += TelElevation;
                //AverageWeightArray[1] += 1;
                //AveragingLocation.push_back("20");
                AveragingLocation2.push_back( 1 );
            }
            else if( Zenithangle > 25.0 && Zenithangle <= 32.5 )  //30
            {
                //AverageArray[2] += TelElevation;
                //AverageWeightArray[2] += 1;
                //AveragingLocation.push_back("30");
                AveragingLocation2.push_back( 2 );
            }
            else if( Zenithangle > 32.5 && Zenithangle <= 37.5 )  //35
            {
                //AverageArray[3] += TelElevation;
                //AverageWeightArray[3] += 1;
                //AveragingLocation.push_back("35");
                AveragingLocation2.push_back( 3 );
            }
            else if( Zenithangle > 37.5 && Zenithangle <= 42.5 )  //40
            {
                //AverageArray[4] += TelElevation;
                //AverageWeightArray[4] += 1;
                //AveragingLocation.push_back("40");
                AveragingLocation2.push_back( 4 );
            }
            else if( Zenithangle > 42.5 && Zenithangle <= 47.5 )  //45
            {
                //AverageArray[5] += TelElevation;
                //AverageWeightArray[5] += 1;
                //AveragingLocation.push_back("45");
                AveragingLocation2.push_back( 5 );
            }
            else if( Zenithangle > 47.5 && Zenithangle <= 52.5 )  //50
            {
                //AverageArray[6] += TelElevation;
                //AverageWeightArray[6] += 1;
                //AveragingLocation.push_back("50");
                AveragingLocation2.push_back( 6 );
            }
            else if( Zenithangle > 52.5 && Zenithangle <= 57.5 )  //55
            {
                //AverageArray[7] += TelElevation;
                //AverageWeightArray[7] += 1;
                //AveragingLocation.push_back("55");
                AveragingLocation2.push_back( 7 );
            }
            else if( Zenithangle > 57.5 && Zenithangle <= 62.5 )  //60
            {
                //AverageArray[8] += TelElevation;
                //AverageWeightArray[8] += 1;
                //AveragingLocation.push_back("60");
                AveragingLocation2.push_back( 8 );
            }
            else  // if ( TelElevation > 62.5 ) //65
            {
                //AverageArray[9] += TelElevation;
                //AverageWeightArray[9] += 1;
                //AveragingLocation.push_back("65");
                AveragingLocation2.push_back( 9 );
            }
            Pars3DForSort.push_back( ForSorting3DParameters );
            
            //AverageElevation += TelElevation/InEntries;
            ForSorting3DParameters.clear();
        }
        /*for ( averaging = 0; averaging < 10; averaging ++){
        	AverageArray[averaging]=AverageArray[averaging]/AverageWeightArray[averaging];c
        }*/
        
        cout << "Indata is loaded and sorted, with " << InEntries << " enrties. " << endl;
        
        // Getting the Pedvars to get the corect Noice level (NOTE, this is here so that the program can be run locally as well as on batch, this info could technically be given to the program) //
        double Totalpedvar = 0;
        
        for( int ped1 = 0; ped1 < t1Entries; ped1++ )
        {
            t1->GetEntry( ped1 );
            Totalpedvar += InOutpedvars1;
        }
        for( int ped2 = 0; ped2 < t2Entries; ped2++ )
        {
            t2->GetEntry( ped2 );
            Totalpedvar += InOutpedvars2;
        }
        for( int ped3 = 0; ped3 < t3Entries; ped3++ )
        {
            t3->GetEntry( ped3 );
            Totalpedvar += InOutpedvars3;
        }
        for( int ped4 = 0; ped4 < t4Entries; ped4++ )
        {
            t4->GetEntry( ped4 );
            Totalpedvar += InOutpedvars4;
        }
        double meanpedvar = Totalpedvar / ( t1Entries + t2Entries + t3Entries + t4Entries );
        int noiceloop;
        int Uselength;	// neede to make the biassearch vector function correctly //
        
        // Setting up the name of the reference file needed //
        char RefName[150];
        if( REFSIMULATIONS == "CARE" )
        {
            if( meanpedvar <= 4.195 )
            {
                noiceloop = 0;	//50
            }
            else if( meanpedvar > 4.195 && meanpedvar <= 5.495 )
            {
                noiceloop = 1;	//80
            }
            else if( meanpedvar > 5.495 && meanpedvar <= 6.83 )
            {
                noiceloop = 2;	//120
            }
            else if( meanpedvar > 6.83 && meanpedvar <= 7.89 )
            {
                noiceloop = 3;	//170
            }
            else if( meanpedvar > 7.89 && meanpedvar <= 9.07 )
            {
                noiceloop = 4;	//230
            }
            else if( meanpedvar > 9.07 && meanpedvar <= 10.15 )
            {
                noiceloop = 5;	//290
            }
            else if( meanpedvar > 10.15 && meanpedvar <= 11.35 )
            {
                noiceloop = 6;	//370
            }
            else   //V
            {
                noiceloop = 7;	//450
            }
            sprintf( RefName, "%s/energy3d_referencetables_V%s_ATM%s_NOISE%s_CARE.root", InDir, VERSION, WinterSummer,  Noice[noiceloop] );
            Uselength = 8;
        }
        else  	//GRISU
        {
            if( meanpedvar <= 3.705 )
            {
                noiceloop = 0;	//075
            }
            else if( meanpedvar > 3.705 && meanpedvar <= 4.32 )
            {
                noiceloop = 1;	//100
            }
            else if( meanpedvar > 4.32 && meanpedvar <= 5.0725 )
            {
                noiceloop = 2;	//150
            }
            else if( meanpedvar > 5.0725 && meanpedvar <= 5.68 )
            {
                noiceloop = 3;	//200
            }
            else if( meanpedvar > 5.68 && meanpedvar <= 6.405 )
            {
                noiceloop = 4;	//250
            }
            else if( meanpedvar > 6.405 && meanpedvar <= 7.278 )
            {
                noiceloop = 5;	//325
            }
            else if( meanpedvar > 7.278 && meanpedvar <= 8.27 )
            {
                noiceloop = 6;	//425
            }
            else if( meanpedvar > 8.27 && meanpedvar <= 9.50 )
            {
                noiceloop = 7;	//550
            }
            else if( meanpedvar > 9.50 && meanpedvar <= 11.00 )
            {
                noiceloop = 8;	//750
            }
            else
            {
                noiceloop = 9;	//1000
            }
            sprintf( RefName, "%s/energy3d_referencetables_V%s_ATM%s_NOISE%s_GRISU.root", InDir, VERSION, WinterSummer,  Noice[noiceloop] );
            Uselength = 10;
        }
        cout << "Reference file being used is: " << RefName << endl;
        
        // Getting Reference data //
        vector< Float_t > RefRow;
        vector < vector < Float_t > > RefVect2D;
        
        TFile* RefTfile = new TFile( RefName );
        TTree* PointerToRefTtreeM = ( TTree* )RefTfile->Get( "ModelPars3D" );
        int NRef = PointerToRefTtreeM->GetEntries();
        
        // Cheking if a new reference data set is needed, if needed the set is gotten, if not it goes on after this part (NOTE, this is needed if one wants to do a list of runs locally and sequentially) //
        if( Oldnoiceloop != noiceloop )
        {
            PointerToRefTtreeM->SetBranchAddress( "Sel3D", 		&RefSel3D );			//These are the values from ForSorting.first, thus should be the final "Sel"values
            PointerToRefTtreeM->SetBranchAddress( "sigmaL3D",	&RefsigmaL3D );
            PointerToRefTtreeM->SetBranchAddress( "sigmaT3D",	&RefsigmaT3D );
            PointerToRefTtreeM->SetBranchAddress( "Nc3D", 		&RefNc3D );
            PointerToRefTtreeM->SetBranchAddress( "RWidth3D", 	&RefRWidth3D );
            PointerToRefTtreeM->SetBranchAddress( "Saz3D",	 	&RefSaz3D );
            PointerToRefTtreeM->SetBranchAddress( "Smax3D", 	&RefSmax3D );
            PointerToRefTtreeM->SetBranchAddress( "Depth3D", 	&RefDepth3D );
            PointerToRefTtreeM->SetBranchAddress( "Xcore3D", 	&RefXcore3D );
            PointerToRefTtreeM->SetBranchAddress( "Ycore3D", 	&RefYcore3D );
            PointerToRefTtreeM->SetBranchAddress( "Goodness3D", 	&RefGoodness3D );
            PointerToRefTtreeM->SetBranchAddress( "Xoff3D", 	&RefXoff3D );
            PointerToRefTtreeM->SetBranchAddress( "Yoff3D", 	&RefYoff3D );
            PointerToRefTtreeM->SetBranchAddress( "MCe0", 		&RefMCe0 );
            
            PointerToRefTtreeM->SetBranchAddress( "EnergySel3D", 	&RefEnergySel3D );		//These are the values from ForSorting.second, thus should go into the return
            PointerToRefTtreeM->SetBranchAddress( "EnergysigmaL3D", &RefEnergysigmaL3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergysigmaT3D", &RefEnergysigmaT3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergyNc3D", 	&RefEnergyNc3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergyRWidth3D", &RefEnergyRWidth3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergySaz3D",	&RefEnergySaz3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergySmax3D", 	&RefEnergySmax3D );
            PointerToRefTtreeM->SetBranchAddress( "EnergyDepth3D",  &RefEnergyDepth3D );
            
            PointerToRefTtreeM->SetBranchAddress( "ReturnSel3D", 	&RefReturnSel3D );		//These are the values from return.first, thus should go into the Energy
            PointerToRefTtreeM->SetBranchAddress( "ReturnsigmaL3D", &RefReturnsigmaL3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnsigmaT3D", &RefReturnsigmaT3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnNc3D", 	&RefReturnNc3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnRWidth3D", &RefReturnRWidth3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnSaz3D",	&RefReturnSaz3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnSmax3D", 	&RefReturnSmax3D );
            PointerToRefTtreeM->SetBranchAddress( "ReturnDepth3D",  &RefReturnDepth3D );
            
            PointerToRefTtreeM->SetBranchAddress( "PlasementSel3D",		&RefPlasementSel3D );		//These are the values from return.second, thus should go into the get the Sel
            PointerToRefTtreeM->SetBranchAddress( "PlasementsigmaL3D",	&RefPlasementsigmaL3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementsigmaT3D", 	&RefPlasementsigmaT3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementNc3D", 		&RefPlasementNc3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementRWidth3D", 	&RefPlasementRWidth3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementSaz3D",		&RefPlasementSaz3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementSmax3D", 	&RefPlasementSmax3D );
            PointerToRefTtreeM->SetBranchAddress( "PlasementDepth3D", 	&RefPlasementDepth3D );
            
            PointerToRefTtreeM->SetBranchAddress( "ErrorRWidth3D", 		&RefErrorRWidth3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorSmax3D",	 	&RefErrorSmax3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorsigmaT3D", 		&RefErrorsigmaT3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorsigmaL3D", 		&RefErrorsigmaL3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorXcore3D",	 	&RefErrorXcore3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorYcore3D",	 	&RefErrorYcore3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorSel3D",	 	&RefErrorSel3D );
            PointerToRefTtreeM->SetBranchAddress( "ErrorSaz3D",	 	&RefErrorSaz3D );
            
            // Fills 2Dvector with Reference values //
            for( int i_ref = 0; i_ref < NRef; i_ref++ )
            {
                PointerToRefTtreeM->GetEntry( i_ref );
                RefRow.push_back( RefSel3D );		//0	//All of are sorted according to themselvs, eccept for the Energy
                RefRow.push_back( RefsigmaL3D );	//1
                RefRow.push_back( RefsigmaT3D ); 	//2
                RefRow.push_back( RefNc3D ); 		//3
                RefRow.push_back( RefRWidth3D );	//4
                RefRow.push_back( RefSaz3D );		//5
                RefRow.push_back( RefSmax3D );		//6
                RefRow.push_back( RefDepth3D );		//7
                RefRow.push_back( RefXcore3D );		//8
                RefRow.push_back( RefYcore3D );		//9
                RefRow.push_back( RefXoff3D );		//10
                RefRow.push_back( RefYoff3D );		//11
                RefRow.push_back( RefGoodness3D ); 	//12
                RefRow.push_back( RefMCe0 );		//13	//This one is not sorted, because it is the energy to which 16-25 is pointing, the so called original positions.
                
                RefRow.push_back( RefEnergySel3D );	//14	//THESE points to the original position, "parallel" to the Energy
                RefRow.push_back( RefEnergysigmaL3D );	//15
                RefRow.push_back( RefEnergysigmaT3D ); 	//16
                RefRow.push_back( RefEnergyNc3D ); 	//17
                RefRow.push_back( RefEnergyRWidth3D ); 	//18
                RefRow.push_back( RefEnergySaz3D );	//19
                RefRow.push_back( RefEnergySmax3D );	//20
                RefRow.push_back( RefEnergyDepth3D ); 	//21
                
                RefRow.push_back( RefReturnSel3D );	//22	//THESE points to the original position of the parameters
                RefRow.push_back( RefReturnsigmaL3D );	//23
                RefRow.push_back( RefReturnsigmaT3D ); 	//24
                RefRow.push_back( RefReturnNc3D ); 	//25
                RefRow.push_back( RefReturnRWidth3D ); 	//26
                RefRow.push_back( RefReturnSaz3D );	//27
                RefRow.push_back( RefReturnSmax3D );	//28
                RefRow.push_back( RefReturnDepth3D ); 	//29
                
                RefRow.push_back( RefPlasementSel3D );		//30 //THESE points to the original position of the Energy belonging to the resorted parameters
                RefRow.push_back( RefPlasementsigmaL3D );	//31
                RefRow.push_back( RefPlasementsigmaT3D ); 	//32
                RefRow.push_back( RefPlasementNc3D ); 		//33
                RefRow.push_back( RefPlasementRWidth3D );	//34
                RefRow.push_back( RefPlasementSaz3D );		//35
                RefRow.push_back( RefPlasementSmax3D );		//36
                RefRow.push_back( RefPlasementDepth3D );	//37
                
                RefRow.push_back( RefErrorSel3D );	//38
                RefRow.push_back( RefErrorsigmaL3D );	//39
                RefRow.push_back( RefErrorsigmaT3D );	//40
                RefRow.push_back( RefErrorNc3D ); 	//41
                RefRow.push_back( RefErrorRWidth3D );	//42
                RefRow.push_back( RefErrorSaz3D );	//43
                RefRow.push_back( RefErrorSmax3D );	//44
                RefRow.push_back( RefErrorXcore3D );	//45
                RefRow.push_back( RefErrorYcore3D );	//46
                
                
                RefVect2D.push_back( RefRow );
                RefRow.clear();
            }
        }
        Oldnoiceloop = noiceloop;
        
        //Bias correction data //
        // temporary variables to store scanned values //
        int SimTelV;
        int SimNoise;
        int SimEl;
        int SimWintSum;
        int SimCorrType;
        float FirstCorr;
        float SecondCorr;
        float ThirdCorr;
        
        // storage vectors for storing values for later searching //
        vector< int > VSimTelV;
        vector< int > VSimNoise;
        vector< int > VSimEl;
        vector< int > VSimWintSum;
        vector< int > VSimCorrType;
        vector< float > VFirstCorr;
        vector< float > VSecondCorr;
        vector< float > VThirdCorr;
        vector< string > VSimCorr;
        
        string SC;
        // loop over lines in file and store the name of the corection and the three corresponding values //
        for( int bcl = 0; bcl < biascorrlength; bcl++ )
        {
            char* inbiascorrsecond = new char [INBIASCORR[bcl].length() + 1];
            strcpy( inbiascorrsecond, INBIASCORR[bcl].c_str() );
            sscanf( inbiascorrsecond, "%d_%d_%d_%d_%d %f %f %f" , &SimTelV, &SimNoise, &SimEl, &SimWintSum, &SimCorrType, &FirstCorr, &SecondCorr, &ThirdCorr ) ;
            VFirstCorr.push_back( FirstCorr );
            VSecondCorr.push_back( SecondCorr );
            VThirdCorr.push_back( ThirdCorr );
            char SimCorr[100];
            if( SimEl == 0 )
            {
                sprintf( SimCorr, "%d_%d_%d%d_%d_%d", SimTelV, SimNoise, SimEl, SimEl, SimWintSum, SimCorrType );
            }
            else
            {
                sprintf( SimCorr, "%d_%d_%d_%d_%d", SimTelV, SimNoise, SimEl, SimWintSum, SimCorrType );
            }
            SC = SimCorr;
            VSimCorr.push_back( SC );
        }
        
        char WantedBiascorr[100];
        float_t BiasCorr[ Uselength ][10];
        
        for( int somevalue = 0; somevalue < 10; somevalue++ ) 	// Setting up a mane to search for //
        {
            sprintf( WantedBiascorr, "9%s%d_%s_%s_%s_1", VERSION, SIMend, Noice[noiceloop], Deg[ somevalue ], WinterSummer );
            string WBias = WantedBiascorr;
            
            
            for( int findCorr = 0 ; findCorr < biascorrlength - 3; findCorr += 3 ) 	// Getting the usefull correction values
            {
                if( VSimCorr[ findCorr ].compare( WBias ) == 0 )
                {
                    BiasCorr[ somevalue ][ 0 ] = intDeg[ somevalue ];
                    BiasCorr[ somevalue ][ 1 ] = VFirstCorr[  findCorr ];
                    BiasCorr[ somevalue ][ 2 ] = VSecondCorr[ findCorr ];
                    BiasCorr[ somevalue ][ 3 ] = VThirdCorr[  findCorr ];
                    BiasCorr[ somevalue ][ 4 ] = VFirstCorr[  findCorr + 1 ];
                    BiasCorr[ somevalue ][ 5 ] = VSecondCorr[ findCorr + 1 ];
                    BiasCorr[ somevalue ][ 6 ] = VThirdCorr[  findCorr + 1 ];
                    BiasCorr[ somevalue ][ 7 ] = VFirstCorr[  findCorr + 2 ];
                    BiasCorr[ somevalue ][ 8 ] = VSecondCorr[ findCorr + 2 ];
                    BiasCorr[ somevalue ][ 9 ] = VThirdCorr[  findCorr + 2 ];
                }
            }
        }
        
        // Setting up the biascorrection, well tested using landau and second degree polynomial//
        TF1* BIASCORlow  = new TF1( "BIASCORlow", "landau", 0.05, 0.5 );
        TF1* BIASCORmid  = new TF1( "BIASCORmid", "pol2", 0.5, 5.0 );
        TF1* BIASCORhigh = new TF1( "BIASCORhigh", "pol2", 0.5, 5.0 );
        
        
        const int NumberOfParameters = 8;
        int TheEnergyIsAt = 13;
        int RefEnergyStart = TheEnergyIsAt + 1;
        int RefReturnStart = RefEnergyStart + NumberOfParameters;
        int RefPlasementStart = RefReturnStart + NumberOfParameters;
        
        cout << "Reference dataset for noice " << Noice[noiceloop] << " have been loaded with " << NRef << " entires. " << endl;
        
        
        /////////PRE-work is done, on to the analysis////////
        
        // Preparing for Energy3d search //
        double div[ 10 ] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};	// Will give output every 10 percent of data analysis, just to keep check //
        int numdiv = 0;
        
        Float_t* ENER = new Float_t[ InEntries ];
        Float_t* GBU  = new Float_t[ InEntries ];
        Float_t* Pros = new Float_t[ InEntries ];
        
        Float_t* Pros0 = new Float_t[ InEntries ];
        Float_t* Pros1 = new Float_t[ InEntries ];
        Float_t* Pros2 = new Float_t[ InEntries ];
        Float_t* Pros3 = new Float_t[ InEntries ];
        Float_t* Pros4 = new Float_t[ InEntries ];
        Float_t* Pros5 = new Float_t[ InEntries ];
        Float_t* Pros6 = new Float_t[ InEntries ];
        
        Float_t* RefErrSel3D	= new Float_t[ InEntries ];
        Float_t* RefErrSaz3D	= new Float_t[ InEntries ];
        Float_t* RefErrNc3D 	= new Float_t[ InEntries ];
        Float_t* RefErrsigL3D 	= new Float_t[ InEntries ];
        Float_t* RefErrsigT3D 	= new Float_t[ InEntries ];
        Float_t* RefErrRW3D 	= new Float_t[ InEntries ];
        Float_t* RefErrXcor3D 	= new Float_t[ InEntries ];
        Float_t* RefErrYcor3D 	= new Float_t[ InEntries ];
        Float_t* RefErrSmax3D 	= new Float_t[ InEntries ];
        
        Float_t XandY;
        Float_t XandYlow;
        Float_t XandYhigh;
        
        Float_t Mean_Energy3D = 0;
        Float_t Mean_GOODBADUGLY = 0;
        
        Float_t Mean_ProsXandY = 0;
        Float_t Mean_ProsXcor = 0;
        Float_t Mean_ProsYcor = 0;
        Float_t Mean_ProsSel = 0;
        //Float_t Mean_ProsSaz = 0;
        Float_t Mean_ProssigL = 0;
        Float_t Mean_ProssigT = 0;
        Float_t Mean_ProsSmax = 0;
        Float_t Mean_ProsDepth = 0;
        Float_t Mean_ProsNc = 0;
        Float_t Mean_ProsRW = 0;
        
        Float_t Mean_RefErrSel3D = 0;
        Float_t	Mean_RefErrsigL3D = 0;
        Float_t	Mean_RefErrsigT3D = 0;
        Float_t	Mean_RefErrNc3D = 0;
        Float_t	Mean_RefErrRW = 0;
        Float_t	Mean_RefErrSaz3D = 0;
        Float_t	Mean_RefErrSmax3D = 0;
        Float_t	Mean_RefErrXcor3D = 0;
        Float_t	Mean_RefErrYcor3D = 0;
        
        int InPutValueSearch = 3;
        
        int InPutRef = InPutValueSearch + 12;
        int HighSearch;
        int LowSearch;
        int ChangeValueslower;
        int ChangeValueshigher;
        
        Float_t DatavalErrlower;
        Float_t DatavalErrhigher;
        
        int Zerocounter = 0;
        int completeanalysiscounter = 0;
        int fullcounter = 0;
        int closecounter = 0;
        
        int Emisscounter = 0;
        
        // Setting up for biascorrection values //
        //int Bval;
        float_t biaslowone;
        float_t biaslowtwo;
        float_t biaslowthree;
        
        float_t biasmidone;
        float_t biasmidtwo;
        float_t biasmidthree;
        
        float_t biashione;
        float_t biashitwo;
        float_t biashithree;
        
        vector< Float_t > VEnergy3D;
        vector< Float_t > VSel;
        //vector< Float_t > VSaz;
        vector< Float_t > VsigT;
        vector< Float_t > VsigL;
        vector< Float_t > VSmax;
        vector< Float_t > VDepth;
        vector< Float_t > VXcor;
        vector< Float_t > VYcor;
        vector< Float_t > VNc;
        vector< Float_t > VRw;
        
        vector< Float_t > VProsXandY;
        vector< Float_t > VProsSel;
        //vector< Float_t > VProsSaz;
        vector< Float_t > VProssigT;
        vector< Float_t > VProssigL;
        vector< Float_t > VProsSmax;
        vector< Float_t > VProsDepth;
        vector< Float_t > VProsXcor;
        vector< Float_t > VProsYcor;
        vector< Float_t > VProsNc;
        vector< Float_t > VProsRw;
        
        vector< Float_t > VRefErrSel3D;
        vector< Float_t > VRefErrSaz3D;
        vector< Float_t > VRefErrSmax3D;
        vector< Float_t > VRefErrNc3D;
        vector< Float_t > VRefErrsigL3D;
        vector< Float_t > VRefErrsigT3D;
        vector< Float_t > VRefErrRW;
        vector< Float_t > VRefErrXcor3D;
        vector< Float_t > VRefErrYcor3D;
        
        vector< Float_t > VGOODBADUGLY;
        
        // Loop over all events //
        for( int EnergySearchIn = 0 ; EnergySearchIn < InEntries; EnergySearchIn++ )
        {
        
            biaslowone =   BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 1 ];
            biaslowtwo =   BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 2 ];
            biaslowthree = BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 3 ];
            
            biasmidone =   BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 4 ];
            biasmidtwo =   BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 5 ];
            biasmidthree = BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 6 ];
            
            biashione =    BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 7 ];
            biashitwo =    BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 8 ];
            biashithree =  BiasCorr[ AveragingLocation2[ForSort3D.at( EnergySearchIn ).second] ][ 9 ];
            
            BIASCORlow->SetParameters( biaslowone, biaslowtwo, biaslowthree );
            BIASCORmid->SetParameters( biasmidone, biasmidtwo, biasmidthree );
            BIASCORhigh->SetParameters( biashione,  biashitwo,  biashithree );
            
            cout << " once " << biaslowone << endl;
            int singel = 0;
            Float_t Currentbest 	= 10000000;
            
            Float_t Currentbestlow 	= 10000000;
            Float_t CurrentEnergylow = 0;
            Float_t CurrentSellow 	= 10000000;
            Float_t CurrentsigLlow 	= 10000000;
            Float_t CurrentsigTlow 	= 10000000;
            Float_t CurrentSmaxlow 	= 10000000;
            Float_t CurrentNclow 	= 10000000;
            Float_t CurrentRwlow 	= 10000000;
            Float_t CurrentDepthlow	= 10000000;
            
            Float_t Currentbesthigh	= 10000000;
            Float_t CurrentEnergyhigh = 0;
            Float_t CurrentSelhigh	= 10000000;
            Float_t CurrentsigLhigh	= 10000000;
            Float_t CurrentsigThigh	= 10000000;
            Float_t CurrentSmaxhigh	= 10000000;
            Float_t CurrentNchigh	= 10000000;
            Float_t CurrentRwhigh 	= 10000000;
            Float_t CurrentDepthhigh = 10000000;
            
            Float_t Energy3Dlow = 0;
            Float_t Sellow = 0;
            //	Float_t Sazlow = 0;
            Float_t sigTlow = 0;
            Float_t sigLlow = 0;
            Float_t Smaxlow = 0;
            Float_t Depthlow = 0;
            Float_t Xcorlow = 0;
            Float_t Ycorlow = 0;
            Float_t Nclow = 0;
            Float_t Rwlow = 0;
            
            Float_t ProsXandYlow = 0;
            Float_t ProsSellow = 0;
            //	Float_t ProsSazlow = 0;
            Float_t ProssigTlow = 0;
            Float_t ProssigLlow = 0;
            Float_t ProsSmaxlow = 0;
            Float_t ProsDepthlow = 0;
            Float_t ProsXcorlow = 0;
            Float_t ProsYcorlow = 0;
            Float_t ProsNclow = 0;
            Float_t ProsRwlow = 0;
            
            Float_t RefErrSel3Dlow = 0;
            Float_t RefErrSaz3Dlow = 0;
            Float_t RefErrSmax3Dlow = 0;
            Float_t RefErrNc3Dlow = 0;
            Float_t RefErrsigL3Dlow = 0;
            Float_t RefErrsigT3Dlow = 0;
            Float_t RefErrRWlow = 0;
            Float_t RefErrXcor3Dlow = 0;
            Float_t RefErrYcor3Dlow = 0;
            
            Float_t Energy3Dhigh = 0;
            Float_t Selhigh = 0;
            //	Float_t Sazhigh = 0;
            Float_t sigThigh = 0;
            Float_t sigLhigh = 0;
            Float_t Smaxhigh = 0;
            Float_t Depthhigh = 0;
            Float_t Xcorhigh = 0;
            Float_t Ycorhigh = 0;
            Float_t Nchigh = 0;
            Float_t Rwhigh = 0;
            
            Float_t ProsXandYhigh = 0;
            Float_t ProsSelhigh = 0;
            //	Float_t ProsSazhigh = 0;
            Float_t ProssigThigh = 0;
            Float_t ProssigLhigh = 0;
            Float_t ProsSmaxhigh = 0;
            Float_t ProsDepthhigh = 0;
            Float_t ProsXcorhigh = 0;
            Float_t ProsYcorhigh = 0;
            Float_t ProsNchigh = 0;
            Float_t ProsRwhigh = 0;
            
            Float_t RefErrSel3Dhigh = 0;
            Float_t RefErrSaz3Dhigh = 0;
            Float_t RefErrSmax3Dhigh = 0;
            Float_t RefErrNc3Dhigh = 0;
            Float_t RefErrsigL3Dhigh = 0;
            Float_t RefErrsigT3Dhigh = 0;
            Float_t RefErrRWhigh = 0;
            Float_t RefErrXcor3Dhigh = 0;
            Float_t RefErrYcor3Dhigh = 0;
            
            // Ignoring values that can't propperly be given a propper energy due to failed 3D-model reconstruction //
            if(
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] <= 0.0 ||
                
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 12 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 13 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 14 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 15 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 16 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 17 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 18 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 19 ] <= 0.0 ||
                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 20 ] <= 0.0 ||
                
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 5 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ) ||
                TMath::IsNaN( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ InPutRef ] ) )
            {
            
                Zerocounter++;
                
                VGOODBADUGLY.push_back( 0 );
                
                VProsSel.push_back( 0 );
                VProssigL.push_back( 0 );
                VProssigT.push_back( 0 );
                VProsNc.push_back( 0 );
                VProsRw.push_back( 0 );
                //VProsSaz.push_back(0);
                VProsSmax.push_back( 0 );
                VProsDepth.push_back( 0 );
                VProsXcor.push_back( 0 );
                VProsYcor.push_back( 0 );
                
                VEnergy3D.push_back( -100 );
                VSel.push_back( 1 );
                //VSaz.push_back(1);
                VsigL.push_back( 1 );
                VsigT.push_back( 1 );
                VSmax.push_back( 1 );
                VDepth.push_back( 1 );
                VXcor.push_back( 1 );
                VYcor.push_back( 1 );
                VNc.push_back( 1 );
                VRw.push_back( 1 );
                VProsXandY.push_back( 1 );
                
                VRefErrSel3D.push_back( 0 );
                VRefErrsigL3D.push_back( 0 );
                VRefErrsigT3D.push_back( 0 );
                VRefErrNc3D.push_back( 0 );
                VRefErrRW.push_back( 0 );
                VRefErrSaz3D.push_back( 0 );
                VRefErrSmax3D.push_back( 0 );
                VRefErrXcor3D.push_back( 0 );
                VRefErrYcor3D.push_back( 0 );
                
                Emisscounter = 1;
            }
            else    // Here the real adding of Energy to the events will be done, for the events that have been successfully reconstructed with model3d //
            {
            
                // Changing biascorrecton values is needed //
                /*Bval = 0;
                while ( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] > BiasCorr[ Bval ][ 0 ] && Bval < Uselength ){
                	Bval++;
                }*/
                
                /*biaslowone =   BiasCorr[Bval-1][1]+(BiasCorr[Bval][1]-BiasCorr[Bval-1][1])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biaslowtwo =   BiasCorr[Bval-1][2]+(BiasCorr[Bval][2]-BiasCorr[Bval-1][2])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biaslowthree = BiasCorr[Bval-1][3]+(BiasCorr[Bval][3]-BiasCorr[Bval-1][3])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                
                biasmidone =   BiasCorr[Bval-1][4]+(BiasCorr[Bval][4]-BiasCorr[Bval-1][4])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biasmidtwo =   BiasCorr[Bval-1][5]+(BiasCorr[Bval][5]-BiasCorr[Bval-1][5])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biasmidthree = BiasCorr[Bval-1][6]+(BiasCorr[Bval][6]-BiasCorr[Bval-1][6])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                
                biashione =    BiasCorr[Bval-1][7]+(BiasCorr[Bval][7]-BiasCorr[Bval-1][7])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biashitwo =    BiasCorr[Bval-1][8]+(BiasCorr[Bval][8]-BiasCorr[Bval-1][8])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));
                biashithree =  BiasCorr[Bval-1][9]+(BiasCorr[Bval][9]-BiasCorr[Bval-1][9])*((Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-BiasCorr[Bval-1][0])/(BiasCorr[Bval][0]-BiasCorr[Bval-1][0]));*/
                
                /*biaslowone =   BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 1 ];
                biaslowtwo =   BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 2 ];
                biaslowthree = BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 3 ];
                
                biasmidone =   BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 4 ];
                biasmidtwo =   BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 5 ];
                biasmidthree = BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 6 ];
                
                biashione =    BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 7 ];
                biashitwo =    BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 8 ];
                biashithree =  BiasCorr[ AveragingLocation2[ForSort3D.at(EnergySearchIn).second] ][ 9 ];
                
                BIASCORlow->SetParameters(  biaslowone, biaslowtwo, biaslowthree);
                BIASCORmid->SetParameters(  biasmidone, biasmidtwo, biasmidthree);
                BIASCORhigh->SetParameters( biashione,  biashitwo,  biashithree);*/
                
                // Setting up parameters do find best searchrange //
                completeanalysiscounter++;
                ChangeValueslower = 0;
                ChangeValueshigher = 0;
                Float_t ChangeValueMid = 0;
                Float_t DataVal = Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ InPutValueSearch ];
                Float_t DataErr = abs( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ InPutRef ] );
                DatavalErrlower = DataVal - DataErr;
                DatavalErrhigher = DataVal + DataErr;
                Float_t RefHigh = RefVect2D[ NRef - 1 ][InPutValueSearch];
                Float_t RefLow = RefVect2D[ 0 ][InPutValueSearch];
                
                ////// FINDING the best searchrange////
                if( DatavalErrlower < RefLow && RefHigh < DatavalErrhigher && ( DataVal < RefLow || DataVal > RefHigh ) ) // Data+-Err outside range and Data outside range, overlapping bad
                {
                    ChangeValueslower = 0;
                    ChangeValueshigher = NRef - 1;
                }
                else if( DatavalErrlower < RefLow && RefHigh < DatavalErrhigher && DataVal > RefLow && DataVal < RefHigh )  // Data+-Err outside range, Data inside range, overlapping good
                {
                    ChangeValueslower = 0;
                    ChangeValueshigher = NRef - 1;
                    while( RefVect2D[ ChangeValueMid ][InPutValueSearch] <= DataVal && ChangeValueslower + ChangeValueshigher <= NRef - 1 )
                    {
                        ChangeValueMid++;
                    }
                }
                else if( DatavalErrlower < RefLow && DatavalErrhigher >= RefLow && RefHigh >= DatavalErrhigher )  // Data-Err outside range, Data+Err inside range, overlapping
                {
                    ChangeValueslower  = 0;
                    while( RefVect2D[ ChangeValueshigher ][InPutValueSearch] <= DatavalErrhigher && ChangeValueslower + ChangeValueshigher < NRef - 1 )
                    {
                        ChangeValueshigher++;
                    }
                    while( RefVect2D[ ChangeValueMid ][InPutValueSearch] < DataVal && ChangeValueslower + ChangeValueshigher <= NRef - 1 )
                    {
                        ChangeValueMid++;
                    }
                }
                else if( RefLow < DatavalErrlower && DatavalErrlower <= RefHigh && RefHigh <= DatavalErrhigher )  // Data-Err inside range, Data+Err outside range, overlapping
                {
                    while( RefVect2D[ ChangeValueslower ][InPutValueSearch] <= DatavalErrlower && ChangeValueslower < NRef - 1 )
                    {
                        ChangeValueslower++;
                    }
                    ChangeValueshigher = NRef - 1;
                }
                else if( DatavalErrhigher < RefLow )  // Data+-Err outside range, not overlapping, lower
                {
                    ChangeValueslower = 0;
                    ChangeValueshigher = 10000;
                }
                else if( DatavalErrlower > RefHigh )  // Data+-Err outside range, not overlapping, higher
                {
                    ChangeValueslower = NRef - 10001;
                    ChangeValueshigher = NRef - 1;
                }
                else
                {
                    while( RefVect2D[ ChangeValueslower ][InPutValueSearch] <= DatavalErrlower && ChangeValueslower < NRef - 1 ) // Data+-Err inside range, overlapping, best
                    {
                        ChangeValueslower++;
                    }
                    while( RefVect2D[ ChangeValueslower + ChangeValueshigher ][InPutValueSearch] <= DatavalErrhigher && ChangeValueslower + ChangeValueshigher < NRef - 1 )
                    {
                        ChangeValueshigher++;
                    }
                    while( RefVect2D[ ChangeValueMid ][InPutValueSearch] <= DataVal && ChangeValueslower + ChangeValueshigher <= NRef - 1 )
                    {
                        ChangeValueMid++;
                    }
                    ChangeValueshigher = ChangeValueslower + ChangeValueshigher;
                }
                
                //////// Change searchwidth if to small (<1000), arbetrary but tested value, works well//////
                if( ChangeValueshigher - ChangeValueslower <= 1000 && ChangeValueMid <= NRef - 1 && ChangeValueMid >= 0 ) // searchwidth<1000, inside range: make 1000
                {
                    if( ChangeValueMid + 500 >= NRef - 1 )
                    {
                        LowSearch = NRef - 1001;
                        HighSearch = NRef - 1;
                    }
                    else if( ChangeValueMid - 500 <= 0 )
                    {
                        LowSearch = 0;
                        HighSearch = 1000;
                    }
                    else
                    {
                        HighSearch = ChangeValueMid + 500;
                        LowSearch = ChangeValueMid - 500;
                    }
                }
                else   // Any Width > 1000
                {
                    LowSearch = ChangeValueslower;
                    HighSearch = ChangeValueshigher;
                }
                
                XandYlow = 0;
                XandYhigh = 0;
                XandY = 0;
                int test21 = 0;
                int test22 = 0;
                Currentbest 		= 10000000;
                
                Currentbestlow 		= 10000000;
                CurrentEnergylow 	= 0;
                CurrentSellow 		= 10000000;
                CurrentsigLlow		= 10000000;
                CurrentsigTlow 		= 10000000;
                CurrentSmaxlow		= 10000000;
                CurrentNclow		= 10000000;
                CurrentRwlow		= 10000000;
                CurrentDepthlow 	= 10000000;
                Currentbestlow 		= 10000000;
                
                Currentbesthigh		= 10000000;
                CurrentEnergyhigh 	= 0;
                CurrentSelhigh 		= 10000000;
                CurrentsigLhigh 	= 10000000;
                CurrentsigThigh		= 10000000;
                CurrentSmaxhigh		= 10000000;
                CurrentNchigh		= 10000000;
                CurrentRwhigh		= 10000000;
                CurrentDepthhigh 	= 10000000;
                
                // Searching the range //
                for( int ES = LowSearch ; ES <= HighSearch; ES++ )
                {
                    // FINDING A Reference-set that has a few values one one side of the indata that should indicate that the ref energy is higher than the indata //
                    if(	Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] )
                    {
                        test21++;
                        XandY = ( 1.0 / 7 ) * (
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] -
                                           Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][0] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] ) );
                                           
                        // If closer than the previous value choose it //
                        if( XandY < Currentbest )
                        {
                            Currentbest 		= XandY;
                            CurrentEnergylow 	= RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ];
                            CurrentSellow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0];
                            CurrentsigLlow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1];
                            CurrentsigTlow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2];
                            CurrentNclow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3];
                            CurrentRwlow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4];
                            CurrentSmaxlow 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6];
                            CurrentDepthlow 	= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7];
                            
                            ProsSellow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] );
                            ProssigLlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] );
                            ProssigTlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] );
                            ProsNclow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                              Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] );
                            ProsRwlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                              Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] );
                            ProsSmaxlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] );
                            ProsDepthlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] );
                            ProsXcorlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] );
                            ProsYcorlow  	= (	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] );
                                                  
                            RefErrSel3Dlow  	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 38 ] );
                            RefErrsigL3Dlow 	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 39 ] );
                            RefErrsigT3Dlow 	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 40 ] );
                            RefErrNc3Dlow  		= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 41 ] );
                            RefErrRWlow  		= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 42 ] );
                            RefErrSaz3Dlow  	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 43 ] );
                            RefErrSmax3Dlow 	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 44 ] );
                            RefErrXcor3Dlow 	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 45 ] );
                            RefErrYcor3Dlow 	= (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 46 ] );
                            
                            Energy3Dlow  	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ] );
                            
                            Sellow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] );
                            sigLlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] );
                            sigTlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] );
                            Nclow  		= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] );
                            Rwlow  		= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] );
                            Smaxlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] );
                            Depthlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] );
                            Xcorlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] );
                            Ycorlow  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] );
                            
                            ProsXandYlow  = XandY;
                        }
                    }
                    
                    // FINDING A Reference-set that has a few values one the other side of the indata that should indicate that the ref energy is lower than the indata //
                    if(	Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] < RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] &&
                            //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6] &&
                            Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] > RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] )
                    {
                        test22++;
                        XandYhigh = ( 1.0 / 7 ) * (
                                        //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+0]][0]-
                                        //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+0]][0]) +
                                        //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1]-
                                        //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1]) +
                                        //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2]-
                                        //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2]) +
                                        ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                               Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] ) +
                                        ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                               Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] ) +
                                        //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6]-
                                        //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6]) +
                                        ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                               Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] ) );
                                               
                        // If closer than the previous value choose it //
                        if( XandYhigh < Currentbesthigh )
                        {
                            Currentbesthigh 	= XandYhigh;
                            CurrentEnergyhigh 	= RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ];
                            CurrentSelhigh 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0];
                            CurrentsigLhigh		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1];
                            CurrentsigThigh		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2];
                            CurrentNchigh		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3];
                            CurrentRwhigh 		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4];
                            CurrentSmaxhigh		= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6];
                            CurrentDepthhigh 	= RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7];
                            
                            ProsSelhigh  	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] );
                            ProssigLhigh 	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] );
                            ProssigThigh 	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] );
                            ProsNchigh  	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] );
                            ProsRwhigh  	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] );
                            ProsSmaxhigh 	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] );
                            ProsDepthhigh	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] );
                            ProsXcorhigh  	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] );
                            ProsYcorhigh  	= ( ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] );
                                                  
                            RefErrSel3Dhigh  = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 38 ] );
                            RefErrsigL3Dhigh = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 39 ] );
                            RefErrsigT3Dhigh = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 40 ] );
                            RefErrNc3Dhigh   = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 41 ] );
                            RefErrRWhigh     = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 42 ] );
                            RefErrSaz3Dhigh  = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 43 ] );
                            RefErrSmax3Dhigh = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 44 ] );
                            RefErrXcor3Dhigh = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 45 ] );
                            RefErrYcor3Dhigh = (	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 46 ] );
                            
                            Energy3Dhigh  	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ] );
                            
                            Selhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] );
                            sigLhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] );
                            sigThigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] );
                            Nchigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] );
                            Rwhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] );
                            Smaxhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] );
                            Depthhigh 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] );
                            Xcorhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] );
                            Ycorhigh  	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] );
                            
                            ProsXandYhigh  = XandYhigh;
                        }
                    }
                }
                
                // If two references sets are found mark as good and go on //
                if( test21 >= 1 && test22 >= 1 )
                {
                    fullcounter++;
                    VGOODBADUGLY.push_back( 1 );
                }
                else  	//If two were not found, search for the clossest value and mark as good
                {
                    closecounter++;
                    XandY = 42;
                    singel = 1;
                    Currentbest = 100000000;
                    VGOODBADUGLY.push_back( 1 );
                    for( int ES = LowSearch ; ES <= HighSearch; ES++ )
                    {
                        XandY = ( 1.0 / 7 ) * (
                                    //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+0]][0]-
                                    //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+0]][0]) +
                                    //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1]-
                                    //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+1]][1]) +
                                    //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2]-
                                    //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+2]][2]) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] ) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] ) +
                                    //(abs(RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6]-
                                    //Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6])/RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart+InPutValueSearch]][RefPlasementStart+6]][6]) +
                                    ( abs( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                           Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7] ) / RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] ) );
                                           
                        // If closer than the previous value, choose it //
                        if( XandY < Currentbest )
                        {
                            Currentbest = XandY;
                            if( singel == 1 )
                            {
                                VProsSel.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] -
                                                      Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] );
                                VProssigL.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] );
                                VProssigT.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] );
                                VProsNc.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                                      Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] );
                                VProsRw.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                                      Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] );
                                VProsSmax.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] );
                                VProsDepth.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] );
                                VProsXcor.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] );
                                VProsYcor.push_back(	( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] -
                                                          Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] );
                                                          
                                VEnergy3D.push_back(	RefVect2D[ RefVect2D[RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ] );
                                VSel.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] );
                                VsigL.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] );
                                VsigT.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] );
                                VNc.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] );
                                VRw.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] );
                                VSmax.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] );
                                VDepth.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] );
                                VXcor.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] );
                                VYcor.push_back(	RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] );
                                
                                VRefErrSel3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 38 ] );
                                VRefErrsigL3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 39 ] );
                                VRefErrsigT3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 40 ] );
                                VRefErrNc3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 41 ] );
                                VRefErrRW.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 42 ] );
                                VRefErrSaz3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 43 ] );
                                VRefErrSmax3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 44 ] );
                                VRefErrXcor3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 45 ] );
                                VRefErrYcor3D.push_back(	RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 46 ] );
                                
                                VProsXandY.push_back( XandY );
                                singel = 0;
                            }
                            else   // If no value is there set it as the first, but liekly bad, value //
                            {
                                VProsSel[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 0 ] ;
                                VProssigL[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 1 ] ;
                                VProssigT[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 2 ] ;
                                VProsNc[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] -
                                                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 3 ] ;
                                VProsRw[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] -
                                                Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 4 ] ;
                                VProsSmax[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 6 ] ;
                                VProsDepth[0] = ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] -
                                                  Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 7 ] ;
                                VProsXcor[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 8 ] ;
                                VProsYcor[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] -
                                                    Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] ) / Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 9 ] ;
                                                    
                                VEnergy3D[0] 	= ( RefVect2D[ RefVect2D[RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ TheEnergyIsAt ] );
                                VSel[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 0]][0] );
                                VsigL[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 1]][1] );
                                VsigT[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 2]][2] );
                                VNc[0] 		= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 3]][3] );
                                VRw[0] 		= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 4]][4] );
                                VSmax[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 6]][6] );
                                VDepth[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 7]][7] );
                                VXcor[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 8]][8] );
                                VYcor[0] 	= ( RefVect2D[RefVect2D[RefVect2D[ES][RefEnergyStart + InPutValueSearch]][RefPlasementStart + 9]][9] );
                                
                                VRefErrSel3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 38 ] );
                                VRefErrsigL3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 39 ] );
                                VRefErrsigT3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 40 ] );
                                VRefErrNc3D[0]  	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 41 ] );
                                VRefErrRW[0] 		= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 42 ] );
                                VRefErrSaz3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 43 ] );
                                VRefErrSmax3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 44 ] );
                                VRefErrXcor3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 45 ] );
                                VRefErrYcor3D[0] 	= ( RefVect2D[ RefVect2D[ RefVect2D[ ES ][ RefEnergyStart + InPutValueSearch]][ RefReturnStart]][ 46 ] );
                                
                                VProsXandY[0] =  XandY;
                            }
                        }
                    }
                }
            }
            
            Mean_Energy3D = 0;
            Mean_GOODBADUGLY = 0;
            
            Mean_ProsXandY = 0;
            Mean_ProsSel = 0;
            Mean_ProssigL = 0;
            Mean_ProssigT = 0;
            Mean_ProsSmax = 0;
            Mean_ProsDepth = 0;
            Mean_ProsRW = 0;
            Mean_ProsNc = 0;
            Mean_ProsXcor = 0;
            Mean_ProsYcor = 0;
            
            Mean_RefErrSel3D = 0;
            Mean_RefErrsigL3D = 0;
            Mean_RefErrsigT3D = 0;
            Mean_RefErrNc3D = 0;
            Mean_RefErrRW = 0;
            Mean_RefErrSaz3D = 0;
            Mean_RefErrSmax3D = 0;
            Mean_RefErrXcor3D = 0;
            Mean_RefErrYcor3D = 0;
            int CurrentSize = int( VEnergy3D.size() );
            Emisscounter = 0;
            
            // Testing if there are two closest values, if so interpolate to get better energy value //
            if( CurrentEnergylow > 0 && CurrentEnergyhigh > 0 )
            {
                // Change low and high value if neede //
                if( CurrentEnergyhigh < CurrentEnergylow )
                {
                    EL = CurrentEnergyhigh;
                    EH = CurrentEnergylow;
                    selL = CurrentSelhigh;
                    selH = CurrentSellow;
                    sigLL = CurrentsigLhigh;
                    sigLH = CurrentsigLlow;
                    sigTL = CurrentsigThigh;
                    sigTH = CurrentsigTlow;
                    NcL = CurrentNchigh;
                    NcH = CurrentNclow;
                    RWL = CurrentRwhigh;
                    RWH = CurrentRwlow;
                    SML = CurrentSmaxhigh;
                    SMH = CurrentSmaxlow;
                    DL = CurrentDepthhigh;
                    DH = CurrentDepthlow;
                }
                else  	//Did not need to change high and low values
                {
                    EL = CurrentEnergylow;
                    EH = CurrentEnergyhigh;
                    selL = CurrentSellow;
                    selH = CurrentSelhigh;
                    sigLL = CurrentsigLlow;
                    sigLH = CurrentsigLhigh;
                    sigTL = CurrentsigTlow;
                    sigTH = CurrentsigThigh;
                    NcL = CurrentNclow;
                    NcH = CurrentNchigh;
                    RWL = CurrentRwlow;
                    RWH = CurrentRwhigh;
                    SML = CurrentSmaxlow;
                    SMH = CurrentSmaxhigh;
                    DL = CurrentDepthlow;
                    DH = CurrentDepthhigh;
                }
                
                
                Float_t Addon = 0;
                Float_t Addondivision = 0;
                /*	if ( ( selH < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0] && selL > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]) || (selH > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0] && selL < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0])) {
                		Addon = Addon + ( EH-EL)*( (Pars3DForSort[ForSort3D.at(EnergySearchIn).second][0]-selL)/(selH-selL) );
                		Addondivision++;
                	}
                	if ( ( sigLH < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][1] && sigLL > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][1]) || (sigLH > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][1] && sigLL < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][1])) {
                		Addon = Addon + ( EH-EL)*( (Pars3DForSort[ForSort3D.at(EnergySearchIn).second][1]-sigLL)/(sigLH-sigLL) );
                		Addondivision++;
                	}
                	if ( ( sigTH < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][2] && sigTL > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][2]) || (sigTH > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][2] && sigTL < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][2])) {
                		Addon = Addon + ( EH-EL)*( (Pars3DForSort[ForSort3D.at(EnergySearchIn).second][2]-sigTL)/(sigTH-sigTL) );
                		Addondivision++;
                	}*/
                if( ( NcH < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][3] && NcL > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][3] ) || ( NcH > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][3] && NcL < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][3] ) )
                {
                    Addon = Addon + ( EH - EL ) * ( ( Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][3] - NcL ) / ( NcH - NcL ) );
                    Addondivision++;
                }
                if( ( RWH < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][4] && RWL > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][4] ) || ( RWH > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][4] && RWL < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][4] ) )
                {
                    Addon = Addon + ( EH - EL ) * ( ( Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][4] - RWL ) / ( RWH - RWL ) );
                    Addondivision++;
                }
                /*if ( ( SMH < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][6] && SML > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][6]) || (SMH > Pars3DForSort[ForSort3D.at(EnergySearchIn).second][6] && SML < Pars3DForSort[ForSort3D.at(EnergySearchIn).second][6])) {
                	Addon = Addon + ( EH-EL)*( (Pars3DForSort[ForSort3D.at(EnergySearchIn).second][6]-SML)/(SMH-SML) );
                	Addondivision++;
                }*/
                if( ( DH < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][7] && DL > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][7] ) || ( DH > Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][7] && DL < Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][7] ) )
                {
                    Addon = Addon + ( EH - EL ) * ( ( Pars3DForSort[ForSort3D.at( EnergySearchIn ).second][7] - DL ) / ( DH - DL ) );
                    Addondivision++;
                }
                
                Mean_Energy3D = EL + Addon / Addondivision;
                Mean_GOODBADUGLY = 	VGOODBADUGLY[0];
                
                Mean_ProsXandY =	( ProsXandYlow  	+ ProsXandYhigh ) / 2;	//Should these be treated the same way as the energy?
                Mean_ProsSel =	( ProsSellow  	+ ProsSelhigh	) / 2;
                Mean_ProssigL =	( ProssigLlow  	+ ProssigLhigh	) / 2;
                Mean_ProssigT =	( ProssigTlow  	+ ProssigThigh ) / 2;
                Mean_ProsSmax =	( ProsSmaxlow  	+ ProsSmaxhigh ) / 2;
                Mean_ProsDepth =	( ProsDepthlow  	+ ProsDepthhigh ) / 2;
                Mean_ProsRW =	( ProsRwlow	+ ProsRwhigh	) / 2;
                Mean_ProsNc =	( ProsNclow  	+ ProsNchigh	) / 2;
                Mean_ProsXcor =	( ProsXcorlow  	+ ProsXcorhigh ) / 2;
                Mean_ProsYcor =	( ProsYcorlow  	+ ProsYcorhigh	) / 2;
                
                Mean_RefErrSel3D =	( RefErrSel3Dlow  	+ RefErrSel3Dhigh	) / 2;
                Mean_RefErrsigL3D =	( RefErrsigL3Dlow  	+ RefErrsigL3Dhigh	) / 2;
                Mean_RefErrsigT3D =	( RefErrsigT3Dlow  	+ RefErrsigT3Dhigh	) / 2;
                Mean_RefErrNc3D =	( RefErrNc3Dlow  	+ RefErrNc3Dhigh	) / 2;
                Mean_RefErrRW =	( RefErrRWlow  		+ RefErrRWhigh	) / 2;
                Mean_RefErrSaz3D =	( RefErrSaz3Dlow  	+ RefErrSaz3Dhigh	) / 2;
                Mean_RefErrSmax3D =	( RefErrSmax3Dlow  	+ RefErrSmax3Dhigh	) / 2;
                Mean_RefErrXcor3D =	( RefErrXcor3Dlow  	+ RefErrXcor3Dhigh	) / 2;
                Mean_RefErrYcor3D =	( RefErrYcor3Dlow  	+ RefErrYcor3Dhigh	) / 2;
                
                Emisscounter++;
            }
            else  	//If only closest values were found, choose them
            {
                Mean_Energy3D = 	VEnergy3D[0];
                Mean_GOODBADUGLY = 	VGOODBADUGLY[0];
                
                Mean_ProsXandY =	VProsXandY[0];
                Mean_ProsSel = 		VProsSel[0];
                Mean_ProssigL =		VProssigL[0];
                Mean_ProssigT =		VProssigT[0];
                Mean_ProsSmax =		VProsSmax[0];
                Mean_ProsDepth =	VProsDepth[0];
                Mean_ProsRW =		VProsRw[0];
                Mean_ProsNc =		VProsNc[0];
                Mean_ProsXcor =		VProsXcor[0];
                Mean_ProsYcor =		VProsYcor[0];
                
                Mean_RefErrSel3D 	= VRefErrSel3D[0];
                Mean_RefErrsigL3D 	= VRefErrsigL3D[0];
                Mean_RefErrsigT3D 	= VRefErrsigT3D[0];
                Mean_RefErrNc3D 	= VRefErrNc3D[0];
                Mean_RefErrRW 		= VRefErrRW[0];
                Mean_RefErrSaz3D 	= VRefErrSaz3D[0];
                Mean_RefErrSmax3D 	= VRefErrSmax3D[0];
                Mean_RefErrXcor3D 	= VRefErrXcor3D[0];
                Mean_RefErrYcor3D 	= VRefErrYcor3D[0];
                
                Emisscounter = 1;
            }
            CurrentSize = Emisscounter;
            Mean_Energy3D = Mean_Energy3D / CurrentSize;
            
            // Clear all values from previous event to make sure nothing gets moved to another event //
            VEnergy3D.clear();
            VSel.clear();
            //VSaz.clear();
            VsigL.clear();
            VsigT.clear();
            VSmax.clear();
            VDepth.clear();
            VNc.clear();
            VXcor.clear();
            VYcor.clear();
            VRw.clear();
            
            VProsXandY.clear();
            VProsSel.clear();
            //VProsSaz.clear();
            VProssigL.clear();
            VProssigT.clear();
            VProsSmax.clear();
            VProsDepth.clear();
            VProsNc.clear();
            VProsXcor.clear();
            VProsYcor.clear();
            VProsRw.clear();
            
            VRefErrSel3D.clear();
            VRefErrsigL3D.clear();
            VRefErrsigT3D.clear();
            VRefErrNc3D.clear();
            VRefErrRW.clear();
            VRefErrSaz3D.clear();
            VRefErrSmax3D.clear();
            VRefErrXcor3D.clear();
            VRefErrYcor3D.clear();
            VGOODBADUGLY.clear();
            
            if( DOBIASCORRECTION == 1 )
            {
                if( Mean_Energy3D >= 0.1 && Mean_Energy3D <= 0.5 )
                {
                    Mean_Energy3D = Mean_Energy3D * ( 1.0 / ( 1 + BIASCORlow->Eval( Mean_Energy3D ) ) );
                }
                else if( Mean_Energy3D > 0.5 && Mean_Energy3D <= 5 )
                {
                    Mean_Energy3D = Mean_Energy3D * ( 1.0 / ( 1 + BIASCORmid->Eval( Mean_Energy3D ) ) );
                }
                else if( Mean_Energy3D > 5 && Mean_Energy3D <= 50 )
                {
                    Mean_Energy3D = Mean_Energy3D * ( 1.0 / ( 1 + BIASCORhigh->Eval( Mean_Energy3D ) ) );
                }
            }
            
            cout << "again " << biaslowone << endl;
            
            /*if ( VEnergy3D[i_in2]>=0.1 && VEnergy3D[i_in2]<=0.5 ){
            	VEnergy3D[i_in2] = VEnergy3D[i_in2]*(1.0/(1+BIASCORlow->Eval(VEnergy3D[i_in2])));
            }
            else if (VEnergy3D[i_in2]>0.5 && VEnergy3D[i_in2]<=5 ){
            	VEnergy3D[i_in2]=VEnergy3D[i_in2]*(1.0/(1+BIASCORmid->Eval(VEnergy3D[i_in2])));
            	}
            else if (VEnergy3D[i_in2]>5 && VEnergy3D[i_in2]<=50 ){
            	VEnergy3D[i_in2]=VEnergy3D[i_in2]*(1.0/(1+BIASCORhigh->Eval(VEnergy3D[i_in2])));
            }*/
            // Putting all values into the right place according to original input //
            ENER[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_Energy3D;
            GBU[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ]	= Mean_GOODBADUGLY;
            Pros[  int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsXandY;
            Pros0[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsSel;
            Pros1[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProssigL;
            Pros2[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProssigT;
            Pros3[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsNc;
            Pros4[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsRW;
            Pros5[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsSmax;
            Pros6[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_ProsDepth;
            
            RefErrSel3D[  int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrSel3D;
            RefErrsigL3D[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrsigL3D;
            RefErrsigT3D[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrsigT3D;
            RefErrNc3D[   int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrNc3D;
            RefErrRW3D[   int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrRW;
            RefErrSaz3D[  int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrSaz3D;
            RefErrSmax3D[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrSmax3D;
            RefErrXcor3D[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrXcor3D;
            RefErrYcor3D[ int( Pars3DForSort[ ForSort3D.at( EnergySearchIn ).second ][ 11 ] ) ] 	= Mean_RefErrYcor3D;
            
            //cout << EnergySearchIn << endl;
            if( EnergySearchIn == ( int )( InEntries * div[ numdiv ] - 1 ) ) 	//To get an idea of how far along the program is, mentioned before
            {
                cout << 100 * div[ numdiv ] << " persent is done, " << EnergySearchIn << " of " << InEntries << endl;
                numdiv++;
            }
        }
        
        // Add branches to the Inputfile. //
        Float_t Energy3D;
        Float_t PRO;
        Float_t PRO0;
        Float_t PRO1;
        Float_t PRO2;
        Float_t PRO3;
        Float_t PRO4;
        Float_t PRO5;
        Float_t PRO6;
        
        Float_t BRefErrRW3D;
        Float_t BRefErrSmax3D;
        Float_t BRefErrsigL3D;
        Float_t BRefErrsigT3D;
        Float_t BRefErrSel3D;
        Float_t BRefErrNc3D;
        Float_t BRefErrSaz3D;
        Float_t BRefErrXcor3D;
        Float_t BRefErrYcor3D;
        
        Float_t GoodBadUgly;
        
        TBranch* NewBranchGBU 	= InTTreem->Branch( "GOODBADUGLY", 	&GoodBadUgly, 	"GOODBADUGLY/F" );
        TBranch* NewBranch 	= InTTreem->Branch( "energy3d", 	&Energy3D, 	"Energy3D/F" );
        TBranch* NewBranch2 	= InTTreem->Branch( "PROS_Total3d", 	&PRO, 		"PROS/F" );
        TBranch* NewBranchp0 	= InTTreem->Branch( "PROS_Sel3d", 	&PRO0, 		"PROS1/F" );
        TBranch* NewBranchp1 	= InTTreem->Branch( "PROS_sigL3d", 	&PRO1, 		"PROS1/F" );
        TBranch* NewBranchp2 	= InTTreem->Branch( "PROS_sigT3d", 	&PRO2, 		"PROS2/F" );
        TBranch* NewBranchp3 	= InTTreem->Branch( "PROS_Nc3d", 	&PRO3, 		"PROS3/F" );
        TBranch* NewBranchp4 	= InTTreem->Branch( "PROS_RWidth3d", 	&PRO4, 		"PROS4/F" );
        TBranch* NewBranchp5 	= InTTreem->Branch( "PROS_Smax3d", 	&PRO5, 		"PROS5/F" );
        TBranch* NewBranchp6 	= InTTreem->Branch( "PROS_Depth3d", 	&PRO6,		"PROS6/F" );
        
        TBranch* NewRefErrBranch0 = InTTreem->Branch( "RefErr_Sel3d",    &BRefErrSel3D,  "ErrTSel3D/F" );
        TBranch* NewRefErrBranch1 = InTTreem->Branch( "RefErr_Saz3d",    &BRefErrSaz3D,  "ErrTSaz3D/F" );
        TBranch* NewRefErrBranch2 = InTTreem->Branch( "RefErr_Smax3d",   &BRefErrSmax3D, "ErrTSmax3D/F" );
        TBranch* NewRefErrBranch3 = InTTreem->Branch( "RefErr_sigmaL3d", &BRefErrsigL3D, "ErrTsigmaL3D/F" );
        TBranch* NewRefErrBranch4 = InTTreem->Branch( "RefErr_sigmaT3d", &BRefErrsigT3D, "ErrTsigmaT3D/F" );
        TBranch* NewRefErrBranch5 = InTTreem->Branch( "RefErr_RWidth3d", &BRefErrRW3D,   "ErrTRWidth3D/F" );
        TBranch* NewRefErrBranch6 = InTTreem->Branch( "RefErr_Nc3d",     &BRefErrNc3D,   "ErrTNc3D/F" );
        TBranch* NewRefErrBranch7 = InTTreem->Branch( "RefErr_Xcore3d",  &BRefErrXcor3D, "ErrTXcore3D/F" );
        TBranch* NewRefErrBranch8 = InTTreem->Branch( "RefErr_Ycore3d",  &BRefErrYcor3D, "ErrTYcore3D/F" );
        
        // Filling the branches //
        for( int FillEnergy3DBranch = 0; FillEnergy3DBranch < InEntries; FillEnergy3DBranch++ )
        {
            Energy3D 	= ENER[ FillEnergy3DBranch ];
            GoodBadUgly 	= GBU[ FillEnergy3DBranch ];
            PRO 		= Pros[ FillEnergy3DBranch ];
            
            PRO0 = Pros0[ FillEnergy3DBranch ];
            PRO1 = Pros1[ FillEnergy3DBranch ];
            PRO2 = Pros2[ FillEnergy3DBranch ];
            PRO3 = Pros3[ FillEnergy3DBranch ];
            PRO4 = Pros4[ FillEnergy3DBranch ];
            PRO5 = Pros5[ FillEnergy3DBranch ];
            PRO6 = Pros6[ FillEnergy3DBranch ];
            
            BRefErrSel3D 	= RefErrSel3D[  FillEnergy3DBranch ];
            BRefErrSaz3D 	= RefErrSaz3D[  FillEnergy3DBranch ];
            BRefErrSmax3D 	= RefErrSmax3D[ FillEnergy3DBranch ];
            BRefErrsigL3D 	= RefErrsigL3D[ FillEnergy3DBranch ];
            BRefErrsigT3D 	= RefErrsigT3D[ FillEnergy3DBranch ];
            BRefErrRW3D 	= RefErrRW3D[   FillEnergy3DBranch ];
            BRefErrNc3D 	= RefErrNc3D[   FillEnergy3DBranch ];
            BRefErrXcor3D 	= RefErrXcor3D[ FillEnergy3DBranch ];
            BRefErrYcor3D 	= RefErrYcor3D[ FillEnergy3DBranch ];
            
            NewBranchGBU->Fill();
            
            NewBranch->Fill();
            NewBranch2->Fill();
            
            NewBranchp0->Fill();
            NewBranchp1->Fill();
            NewBranchp2->Fill();
            NewBranchp3->Fill();
            NewBranchp4->Fill();
            NewBranchp5->Fill();
            NewBranchp6->Fill();
            
            NewRefErrBranch0->Fill();
            NewRefErrBranch1->Fill();
            NewRefErrBranch2->Fill();
            NewRefErrBranch3->Fill();
            NewRefErrBranch4->Fill();
            NewRefErrBranch5->Fill();
            NewRefErrBranch6->Fill();
            NewRefErrBranch7->Fill();
            NewRefErrBranch8->Fill();
        }
        
        // Save only the new version of the tree //
        InTfile->Write( "", TObject::kOverwrite );
        delete RefTfile;
        delete InTfile;
        RefVect2D.clear();
        
        // Extra output if wanted for different checks //
        /*	cout << RunnumberFile << " is done." << endl;
        	cout << "# of Zeros or NaN sets (can't be Energized): " << Zerocounter << endl;
        	cout << "# of events going to full analysis: " << completeanalysiscounter << endl;
        	cout << "# of fully done means: " << fullcounter << endl;
        	cout << "# of closest fit: "<< closecounter << ", fully done means and closes fit should add up to full analysis." << endl;*/
    }
}
