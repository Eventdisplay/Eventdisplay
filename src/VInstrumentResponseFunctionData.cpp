/*! \class VInstrumentResponseFunctionData
    \brief data class for instrumental response functions



*/


#include "VInstrumentResponseFunctionData.h"

VInstrumentResponseFunctionData::VInstrumentResponseFunctionData()
{
    fType = "";
    fType_numeric = 0;
    fName = "";
    fNTel = 0;
    fMCMaxCoreRadius = 0.;
    
    fData = 0;
    
    fListofResponseFunctionTypes.push_back( "angular_resolution" );    // fType_numeric == 0
    fListofResponseFunctionTypes.push_back( "core_resolution" );       // fType_numeric == 1
    fListofResponseFunctionTypes.push_back( "energy_resolution" );     // fType_numeric == 2
    
    setEnergyReconstructionMethod();
    
    fHistogramList = 0;
    setHistogrambinning();
    setArrayCentre();
}

VInstrumentResponseFunctionData::~VInstrumentResponseFunctionData()
{
    if( fHistogramList )
    {
        fHistogramList->Delete();
    }
}

bool VInstrumentResponseFunctionData::initialize( string iName, string iType, unsigned int iNTel, double iMCMaxCoreRadius )
{
    fType_numeric = testResponseFunctionType( iType );
    if( fType_numeric < 0 )
    {
        return false;
    }
    fType = iType;
    
    fName = iName;
    fNTel = iNTel;
    fMCMaxCoreRadius = iMCMaxCoreRadius;
    if( fMCMaxCoreRadius < 1.e-2 )
    {
        fMCMaxCoreRadius = 500.;
    }
    
    fHistogramList = new TList();
    
    // histograms
    vector< string > iHisName;
    vector< string > iHisXaxisName;
    vector< string > iHisYaxisName;
    vector< int >    iHisNbinsX;
    vector< double > iHisXmin;
    vector< double > iHisXmax;
    vector< int >    iHisNbinsY;
    vector< double > iHisYmin;
    vector< double > iHisYmax;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // NOTE:   look at E_HISTOID before changing anything here
    //
    //         the sequence of the histogram definition must be
    //         exactly the same as the one of E_HISTOID
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // angular resolution plots
    if( fType == "angular_resolution" )
    {
        // angular difference vs. reconstructed energy
        iHisName.push_back( "AngE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "angular diff. (R,MC) [deg]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 9000 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // (angular difference)^2 vs. reconstructed energy
        iHisName.push_back( "AngE0_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "(angular diff.)^{2} (R,MC) [deg^{2}]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 9000 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // log(angular difference) vs reconstructed energy
        iHisName.push_back( "AngE0Log_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "log(angular diff. (R,MC) [deg])" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -4. );
        iHisYmax.push_back( 1. );
        // angular resolution vs number of images per telescope
        iHisName.push_back( "AngNImages" + fName );
        iHisXaxisName.push_back( "number of images" );
        iHisYaxisName.push_back( "angular diff. (R,MC) [deg]" );
        iHisNbinsX.push_back( fNTel );
        iHisXmin.push_back( 0.5 );
        iHisXmax.push_back( 0.5 + fNTel );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // angular resolution vs core distance
        iHisName.push_back( "AngCoreDistance" + fName );
        iHisXaxisName.push_back( "distance to array center [m]" );
        iHisYaxisName.push_back( "angular diff. (R,MC) [deg]" );
        iHisNbinsX.push_back( 50 );
        iHisXmin.push_back( 0. );
        iHisXmax.push_back( fMCMaxCoreRadius );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // angular error vs. energy
        iHisName.push_back( "AngErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "angular error (R,MC) [deg]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( -5. );
        iHisYmax.push_back( 5. );
        // not defined here
        iHisName.push_back( "AngRelativeErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "relative angular error (R,MC)" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 1 );
        iHisYmin.push_back( -4.5 );
        iHisYmax.push_back( 4.5 );
        
        // angular difference vs. true energy
        iHisName.push_back( "AngEMC_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "angular diff. (R,MC) [deg]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 9000 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // (angular difference)^2 vs. true energy
        iHisName.push_back( "AngEMC_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "(angular diff.)^{2} (R,MC) [deg^{2}]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 9000 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4.5 );
        // log(angular difference) vs true energy
        iHisName.push_back( "AngEMCLog_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "log(angular diff. (R,MC) [deg])" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -4. );
        iHisYmax.push_back( 1. );
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // core resolution plots
    else if( fType == "core_resolution" )
    {
        // core position difference vs. reconstructed energy
        iHisName.push_back( "CoreE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 300 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300. );
        // (core position difference)^2 vs. reconstructed energy
        iHisName.push_back( "CoreE0_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "(core position diff.)^{2} (R,MC) [m^{2}]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 900 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300.*300. );
        // log(angular difference) vs reconstructed energy
        iHisName.push_back( "CoreE0Log_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "log(core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -1. );
        iHisYmax.push_back( 3. );
        // core position resolution vs number of images per telescope
        iHisName.push_back( "CoreNImages" + fName );
        iHisXaxisName.push_back( "number of images" );
        iHisYaxisName.push_back( "core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( fNTel );
        iHisXmin.push_back( 0.5 );
        iHisXmax.push_back( 0.5 + fNTel );
        iHisNbinsY.push_back( 300 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300. );
        // core position resolution vs core distance
        iHisName.push_back( "CoreCoreDistance" + fName );
        iHisXaxisName.push_back( "distance to array center [m]" );
        iHisYaxisName.push_back( "core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( 50 );
        iHisXmin.push_back( 0. );
        iHisXmax.push_back( fMCMaxCoreRadius );
        iHisNbinsY.push_back( 300 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300. );
        // core position error vs energy
        iHisName.push_back( "CoreErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "core position error (R,MC) [m]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 600 );
        iHisYmin.push_back( -300. );
        iHisYmax.push_back( 300. );
        // not defined here
        iHisName.push_back( "CoreRelativeErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "core position angular error (R,MC)" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 1 );
        iHisYmin.push_back( -5. );
        iHisYmax.push_back( 5. );
        // core position difference vs. true energy
        iHisName.push_back( "CoreEMC_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 300 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300. );
        // (core position difference)^2 vs. true energy
        iHisName.push_back( "CoreEMC_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "(core position diff.)^{2} (R,MC) [m^{2}]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 900 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 300.*300. );
        // log(angular difference) vs true energy
        iHisName.push_back( "CoreEMCLog_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "log(core position diff. (R,MC) [m]" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -1. );
        iHisYmax.push_back( 3. );
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // energy resolution plots
    else if( fType == "energy_resolution" )
    {
        // energy difference vs. reconstructed energy
        iHisName.push_back( "EnergE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "log_{10} E_{rec} - log_{10} E_{MC}" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 4500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // (energy difference)^2 vs. reconstructed energy
        iHisName.push_back( "EnergE0_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "(log_{10} E_{rec} - log_{10} E_{MC})^2" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 4500 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4. );
        // log(angular difference) vs true energy (probably useless)
        iHisName.push_back( "EnergyE0Log_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "log(log_{10} E_{rec} - log_{10} E_{MC})" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -1. );
        iHisYmax.push_back( 3. );
        // energy resolution vs number of images per telescope
        iHisName.push_back( "EnergNImages" + fName );
        iHisXaxisName.push_back( "number of images" );
        iHisYaxisName.push_back( "log_{10} E_{rec} - log_{10} E_{MC}" );
        iHisNbinsX.push_back( fNTel );
        iHisXmin.push_back( 0.5 );
        iHisXmax.push_back( 0.5 + fNTel );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // energy resolution vs core distance
        iHisName.push_back( "EnergCoreDistance" + fName );
        iHisXaxisName.push_back( "distance to array center [m]" );
        iHisYaxisName.push_back( "log_{10} E_{rec} - log_{10} E_{MC}" );
        iHisNbinsX.push_back( 50 );
        iHisXmin.push_back( 0. );
        iHisXmax.push_back( fMCMaxCoreRadius );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // energy reconstruction error vs. energy (used for energy systematics)
        iHisName.push_back( "EnergErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "log_{10} E_{rec} - log_{10} E_{MC}" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // not defined here
        iHisName.push_back( "EnergyRelativeErrorE0_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{rec} [TeV]" );
        iHisYaxisName.push_back( "#Delta energy resolution" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 2500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // energy difference vs. true energy
        iHisName.push_back( "EnergEMC_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "log_{10} E_{rec} - log_{10} E_{MC}" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 4500 );
        iHisYmin.push_back( -2. );
        iHisYmax.push_back( 2. );
        // (energy difference)^2 vs. true energy
        iHisName.push_back( "EnergEMC_2_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "(log_{10} E_{rec} - log_{10} E_{MC})^2" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 4500 );
        iHisYmin.push_back( 0. );
        iHisYmax.push_back( 4. );
        // log(angular difference) vs true energy (probably useless)
        iHisName.push_back( "EnergyEMCLog_" + fName );
        iHisXaxisName.push_back( "log_{10} energy_{MC} [TeV]" );
        iHisYaxisName.push_back( "log(log_{10} E_{rec} - log_{10} E_{MC})" );
        iHisNbinsX.push_back( fHistogrambinningEnergy_TeV_Log );
        iHisXmin.push_back( fHistogrambinningEnergy_Min_Tev_Log );
        iHisXmax.push_back( fHistogrambinningEnergy_Max_Tev_Log );
        iHisNbinsY.push_back( 100 );
        iHisYmin.push_back( -1. );
        iHisYmax.push_back( 3. );
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // create histograms
    /////////////////////////////////////////////////////////////////////////////////////////////////
    char iname[1000] ;
    for( unsigned int i = 0; i < iHisName.size(); i++ )
    {
        // 2D histo
        f2DHisto.push_back( new TH2D( ( "h" + iHisName[i] ).c_str(), "", iHisNbinsX[i], iHisXmin[i],
                                      iHisXmax[i], iHisNbinsY[i],
                                      iHisYmin[i], iHisYmax[i] ) );
        f2DHisto.back()->SetXTitle( iHisXaxisName[i].c_str() );
        f2DHisto.back()->SetYTitle( iHisYaxisName[i].c_str() );
        fHistogramList->Add( f2DHisto.back() );
        
        // corresponding resolution graph
        fResolutionGraph.push_back( new TGraphErrors( 1 ) );
        fResolutionGraph.back()->SetName( ( "g" + iHisName[i] ).c_str() );
        fResolutionGraph.back()->SetTitle();
        fHistogramList->Add( fResolutionGraph.back() );
        
        // most of these histograms wont use any king function fitting, but creating the extra
        // histograms keeps things more organized
        sprintf( iname, "gKingSigma%s_i%d", iHisName[i].c_str(), i ) ;
        fResolutionKingSigmaGraph.push_back( new TGraphErrors( 1 ) );
        fResolutionKingSigmaGraph.back()->SetName( iname ) ;
        fResolutionKingSigmaGraph.back()->SetTitle( "Fitted King Function Sigma Parameter vs Energy" ) ;
        fResolutionKingSigmaGraph.back()->GetXaxis()->SetTitle( iHisXaxisName[i].c_str() ) ;
        fResolutionKingSigmaGraph.back()->GetYaxis()->SetTitle( "Sigma (deg)" ) ;
        fHistogramList->Add( fResolutionKingSigmaGraph.back() );
        
        sprintf( iname, "gKingGamma%s_i%d", iHisName[i].c_str(), i ) ;
        fResolutionKingGammaGraph.push_back( new TGraphErrors( 1 ) );
        fResolutionKingGammaGraph.back()->SetName( iname ) ;
        fResolutionKingGammaGraph.back()->SetTitle() ;
        fResolutionKingGammaGraph.back()->GetXaxis()->SetTitle( iHisXaxisName[i].c_str() ) ;
        fResolutionKingGammaGraph.back()->GetYaxis()->SetTitle( "Fitted King Function Gamma (unitless) vs Energy" ) ;
        fHistogramList->Add( fResolutionKingGammaGraph.back() );
        
        // containment probability
        fContainmentProbability.push_back( 0. );
    }
    
    return true;
}

int VInstrumentResponseFunctionData::testResponseFunctionType( string iType )
{
    for( unsigned int i = 0; i < fListofResponseFunctionTypes.size(); i++ )
    {
        if( iType == fListofResponseFunctionTypes[i] )
        {
            return ( int )i;
        }
    }
    
    cout << "VInstrumentResponseFunctionData::testResponseFunctionType() error: type not found: " << iType << endl;
    cout << "\t available reponse function types are: " << endl;
    for( unsigned int i = 0; i < fListofResponseFunctionTypes.size(); i++ )
    {
        cout << "\t" << fListofResponseFunctionTypes[i] << endl;
    }
    
    return -99;
}

/*
    expect that all quality checks happended before
*/
void VInstrumentResponseFunctionData::fill( double iWeight )
{
    if( !fData )
    {
        return;
    }
    
    // simple quality check (FOV shouldn't be larger than 50 deg)
    if( fData->getXoff() < -50. || fData->getYoff() < -50. )
    {
        return;
    }
    
    // get reconstructed energy
    double iErec_lin = -99.e6;
    iErec_lin = fData->getEnergy_TeV() ;
    if( iErec_lin < 0. )
    {
        return;
    }
    
    double iDiff = -99.e6;
    double iError = -99.e6;
    double iErrorRelative = -99.e6;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // angular resolution
    if( fType_numeric == 0 )
    {
        // angular difference
        iDiff = sqrt( ( fData->getXoff() - fData->MCxoff ) * ( fData->getXoff() - fData->MCxoff ) +
                      ( fData->getYoff() - fData->MCyoff ) * ( fData->getYoff() - fData->MCyoff ) );
        // error
        iError = sqrt( fData->getXoff() * fData->getXoff() + fData->getYoff() * fData->getYoff() ) -
                 sqrt( fData->MCxoff * fData->MCxoff + fData->MCyoff * fData->MCyoff );
        // relative error (not sure if it is useful)
        iErrorRelative = -99.e6;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // core resolution
    else if( fType_numeric == 1 )
    {
        // core difference
        iDiff = sqrt( ( fData->getXcore_M() - fData->MCxcore ) * ( fData->getXcore_M() - fData->MCxcore ) +
                      ( fData->getYcore_M() - fData->MCycore ) * ( fData->getYcore_M() - fData->MCycore ) );
        // core error
        iError = sqrt( fData->getXcore_M() * fData->getXcore_M() + fData->getYcore_M() * fData->getYcore_M() ) - //was xcore both times
                 sqrt( fData->MCxcore * fData->MCxcore + fData->MCycore * fData->MCycore );
        // relative error
        if( sqrt( fData->MCxcore * fData->MCxcore + fData->MCycore * fData->MCycore ) > 0. )
        {
            iErrorRelative = iError / sqrt( fData->MCxcore * fData->MCxcore + fData->MCycore * fData->MCycore );
        }
        else
        {
            iErrorRelative = -99.e6;
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // energy resolution
    else if( fType_numeric == 2 )
    {
        if( fData->MCe0 > 0. )
        {
            iDiff = TMath::Abs( 1. - iErec_lin / fData->MCe0 );
        }
        iError = iDiff;
        if( fData->MCe0 > 0. )
        {
            iErrorRelative = ( iErec_lin - fData->MCe0 ) / fData->MCe0;
        }
        else
        {
            iErrorRelative = -99.e6;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // fill histograms
    
    // difference vs energy
    if( E_DIFF < f2DHisto.size() && f2DHisto[E_DIFF] )
    {
        f2DHisto[E_DIFF]->Fill( log10( iErec_lin ), iDiff, iWeight );
    }
    // difference vs true energy
    if( E_DIFF_MC < f2DHisto.size() && f2DHisto[E_DIFF_MC] )
    {
        f2DHisto[E_DIFF_MC]->Fill( log10( fData->MCe0 ), iDiff, iWeight );
    }
    // squared difference vs energy
    if( E_DIFF2 < f2DHisto.size() && f2DHisto[E_DIFF2] )
    {
        f2DHisto[E_DIFF2]->Fill( log10( iErec_lin ), iDiff * iDiff, iWeight );
    }
    // squared difference vs reconstructed energy
    if( E_DIFF2_MC < f2DHisto.size() && f2DHisto[E_DIFF2_MC] )
    {
        f2DHisto[E_DIFF2_MC]->Fill( log10( fData->MCe0 ), iDiff * iDiff, iWeight );
    }
    // log10 difference vs energy
    if( E_LOGDIFF < f2DHisto.size() && f2DHisto[E_LOGDIFF] && iDiff > 0. )
    {
        f2DHisto[E_LOGDIFF]->Fill( log10( iErec_lin ), log10(iDiff), iWeight );
    }
    // log10 difference vs true energy
    if( E_LOGDIFF_MC < f2DHisto.size() && f2DHisto[E_LOGDIFF_MC] && iDiff > 0. )
    {
        f2DHisto[E_LOGDIFF_MC]->Fill( log10( fData->MCe0 ), log10(iDiff), iWeight );
    }
    
    // difference vs number of images
    if( E_NIMAG < f2DHisto.size() && f2DHisto[E_NIMAG] )
    {
        f2DHisto[E_NIMAG]->Fill( fData->getNImages(), iDiff, iWeight );
    }
    
    // difference vs core distance
    if( E_DIST < f2DHisto.size() && f2DHisto[E_DIST] )
    {
        f2DHisto[E_DIST]->Fill( sqrt( ( fData->MCxcore - fArrayCentre_X ) * ( fData->MCxcore - fArrayCentre_X ) +
                                      ( fData->MCycore - fArrayCentre_Y ) * ( fData->MCycore - fArrayCentre_Y ) ), iDiff, iWeight );
    }
    
    // error vs energy
    if( E_ERROR < f2DHisto.size() && f2DHisto[E_ERROR] )
    {
        f2DHisto[E_ERROR]->Fill( log10( iErec_lin ), iError, iWeight );
    }
    
    // relative error vs energy
    if( E_RELA < f2DHisto.size() && f2DHisto[E_RELA] && iErrorRelative > -98.e6 )
    {
        f2DHisto[E_RELA]->Fill( log10( iErec_lin ), iErrorRelative, iWeight );
    }
}

/*
 * finalize calculation of instrument response functions
 *
 * calculate containment graphs
*/
bool VInstrumentResponseFunctionData::terminate( double iContainmentProbability )
{
    for( unsigned int i = 0; i < f2DHisto.size(); i++ )
    {
        // calculate XX% values (default is 68%)
        fContainmentProbability[i] = iContainmentProbability;
        if( i != E_RELA )
        {
            calculateResolution( f2DHisto[i], fResolutionGraph[i], f2DHisto[i]->GetName(), iContainmentProbability,
                                 fResolutionKingSigmaGraph[i], fResolutionKingGammaGraph[i] );
        }
        // for relative plots get mean and spread from each bin in the histogram
        else
        {
            get_Profile_from_TH2D( f2DHisto[i], fResolutionGraph[i], "meanS" );
        }
    }
    
    return true;
}


/*!
  Integrate the psf from 0 -> r,
  return the containment fraction, from 0.0 to 1.0
  par is a double[2] list of sigma, gamma
  NOTE: THIS IS NOT THE PSF,
  its the integral of the psf (the containment fraction)
*/
Double_t kingfunc( Double_t* r, Double_t* par )
{
    // from mathematica snippet in KingPsf.nb
    return 1 - pow( 2, -1 + par[1] ) * pow( par[0], -2 + 2 * par[1] ) * pow( par[1] / ( pow( r[0], 2 ) + 2 * par[1] * pow( par[0], 2 ) ), -1 + par[1] ) ;
}


/*!

    calculate threshold (usually 68%) reconstruction accuracy from 2D histogram

    iHistogram is a 2D histogram (e.g. angular difference vs energy)

*/
TList*  VInstrumentResponseFunctionData::calculateResolution( TH2D* iHistogram, TGraphErrors* iResult,
        string iHistoName, double iContainmentProbability,
        TGraphErrors* iResultKingSigma, TGraphErrors* iResultKingGamma )
{
    if( !iHistogram || !iResult )
    {
        return 0;
    }
    
    TH1D* iTemp = 0;
    TH1D* iTempCumulative = 0 ;
    TList* hList = new TList();
    
    char iname[800];
    
    // set number of points in graph
    iResult->Set( iHistogram->GetNbinsX() );
    
    // temporary vectors
    vector< double > vEnergy;
    vector< double > vRes;
    vector< double > vResE;
    vector< double > vKingEnergy;
    vector< double > vKingSigma;
    vector< double > vKingSigmaError;
    vector< double > vKingGamma;
    vector< double > vKingGammaError;
    
    // whether or not to do a king function fit of the current resolution
    bool doKingFit = false ;
    
    // only do king function fitting for
    //    angular resolution mode
    //    difference-vs-energy hist
    //    containment radii < 70 % (so theoretically, this should be True only once for the 0.68 containment radii)
    doKingFit = ( ( fType_numeric == 0 ) && ( f2DHisto[E_DIFF] == iHistogram ) && ( iContainmentProbability < 0.7 ) ) ;
    
    double i_energy = 0.;
    
    //////////////////////////////////////////////////////////////////////////////
    // loop over all energy bins and project each bin into a TH1D
    for( int i = 1; i <= iHistogram->GetNbinsX(); i++ )
    {
    
        // define temporary histogram and fill with projection
        if( iHistoName.size() > 0 )
        {
            sprintf( iname, "%s_%d", iHistoName.c_str(), i );
        }
        else
        {
            sprintf( iname, "iH_%d", i );
        }
        iTemp = iHistogram->ProjectionY( iname, i, i );
        sprintf( iname, "log_{10} E_{0} = %.2f", iHistogram->GetXaxis()->GetBinCenter( i ) );
        iTemp->SetTitle( iname );
        
        i_energy = iHistogram->GetXaxis()->GetBinCenter( i );
        
        //////////////////////////////////////////////////////////
        // calculate containment
        double iTotSum = 0.;
        // get total number of events in histogram
        for( int j = 1; j <= iTemp->GetNbinsX(); j++ )
        {
            iTotSum += iTemp->GetBinContent( j );
        }
        if( iTotSum  > 0. )
        {
            double iTempSum = 0.;
            int iTempNBins = 0;
            for( int j = 1; j <= iTemp->GetNbinsX(); j++ )
            {
                iTempSum += iTemp->GetBinContent( j );
                if( iTemp->GetBinContent( j ) > 0 )
                {
                    iTempNBins++;
                }
                if( iTempSum / iTotSum  > iContainmentProbability )
                {
                    // require at least 20 events to calculate
                    // a good containment radius and RMS
                    if( iTemp->GetEntries() > 20. && iTempNBins > 4 )
                    {
                        vEnergy.push_back( i_energy );
                        vRes.push_back( iTemp->GetBinCenter( j ) );
                        vResE.push_back( iTemp->GetRMS() / sqrt( iTemp->GetEntries() ) );
                    }
                    break;
                }
            }
        }
        
        if( doKingFit )
        {
            TDirectory *current_dir = gDirectory;
            if( !current_dir->FindObject( "kingfunc" ) )
            {
                if( current_dir->mkdir( "kingfunc" ) )
                {
                     current_dir->cd( "kingfunc" );
                }
            }
            else
            {
                current_dir->cd( "kingfunc" );
            }
        
            /////////////////////////////////////////////////////////////
            // fit the cumulative iTemp histogram with a king function //
            /////////////////////////////////////////////////////////////
            double totalEvents      = iTemp->GetEntries() ;
            char kingfuncname[1000] ;
            //double binarea = 0.0 ; // BINSCALING
            
            // only proceed if we have enough events in iTemp for a stable fit
            if( totalEvents > 20.0 )
            {
            
                // setup cumulative iTemp
                iTempCumulative = ( TH1D* ) iTemp->Clone() ;
                if( iHistoName.size() > 0 )
                {
                    sprintf( iname, "%s_%d_cumulative", iHistoName.c_str(), i );
                }
                else
                {
                    sprintf( iname, "iH_%d_cumulative", i );
                }
                iTempCumulative->SetName( iname ) ;
                sprintf( iname, "%s Cumulative (For King Function Fitting)", iTemp->GetTitle() ) ;
                iTempCumulative->SetTitle( iname ) ;
                
                // fill bin contents of iTempCumulative
                double cumulativeEvents = 0.0 ;
                for( int j = 1 ; j <= iTempCumulative->GetNbinsX() ; j++ )
                {
                    cumulativeEvents += iTemp->GetBinContent( j ) ;
                    // binarea = 2 * TMath::Pi() * pow( iTempCumulative->GetXaxis()->GetBinUpEdge(j) , 2 ) ; // BINSCALING
                    iTempCumulative->SetBinContent( j, cumulativeEvents ) ;
                    //iTempCumulative->SetBinContent( j, cumulativeEvents / binarea ) ; // BINSCALING
                }
                
                // convert bin contents from '# of events' to '% of all events'
                iTempCumulative->Scale( 1 / totalEvents ) ;
                // binarea = 2 * TMath::Pi() * pow( iTempCumulative->GetXaxis()->GetBinUpEdge( iTempCumulative->GetNbinsX()-1 ) , 2 ) ; // BINSCALING
                // iTempCumulative->Scale( binarea / totalEvents ) // BINSCALING
                
                // set up a king function for fitting to iTempCumulative
                sprintf( kingfuncname, "%s_kingfunc", iTempCumulative->GetName() ) ;
                //sprintf( kingfuncname, "%s_kingfunc2", iTempCumulative->GetName() ) ;
                TF1* fitfunc = new TF1( kingfuncname, kingfunc, 0.0, 4.0, 2 ) ;
                //TF1 * fitfunc = new TF1( kingfuncname, kingfunc2, 0.0, 4.0, 2 ) ;
                fitfunc->SetParName( 0, "Sigma" ) ;
                fitfunc->SetParName( 1, "Gamma" ) ;
                fitfunc->SetParameter( 0,  0.08 ) ;   // basic ballpark guess for sigma
                fitfunc->SetParameter( 1,  1.9 ) ;    // basic ballpark guess for gamma
                fitfunc->SetParLimits( 1,  0.1   , 50.0 ) ; // gamma > 50.0 only happens when the fit fails
                
                // Fit the king function
                // TH1::Fit() :
                // M : Use TMinuit IMPROVE command to escape local fitting minima
                // E : Better errors estimation using Minos technique
                // V : Verbose text logging mode
                // Q: Quite - don't printout fit results
                // 0: do not plot
                iTempCumulative->Fit( kingfuncname, "MEQ0" );
                
                // save the fitted parameters
                double sigma       = fitfunc->GetParameter( 0 ) ;
                double sigma_error = fitfunc->GetParError( 0 ) ;
                double gamma       = fitfunc->GetParameter( 1 ) ;
                double gamma_error = fitfunc->GetParError( 1 ) ;
                vKingEnergy.push_back( i_energy ) ;
                vKingSigma.push_back( sigma ) ;
                vKingSigmaError.push_back( sigma_error ) ;
                vKingGamma.push_back( gamma ) ;
                vKingGammaError.push_back( gamma_error ) ;
                
                iTemp->Write();
                iTempCumulative->Write();
                fitfunc->Write();
                
            }
            current_dir->cd();
        }
        
        
        // save or delete temporary histograms
        if( iHistoName.size() > 0 )
        {
            if( iTemp )
            {
                hList->Add( iTemp ) ;
            }
            if( iTempCumulative )
            {
                hList->Add( iTempCumulative ) ;
            }
        }
        else
        {
            if( iTemp )
            {
                delete iTemp           ;
            }
            if( iTempCumulative )
            {
                delete iTempCumulative ;
            }
        }
    }
    
    // fill graph
    iResult->Set( ( int )vEnergy.size() );
    for( unsigned i = 0; i < vEnergy.size(); i++ )
    {
        iResult->SetPoint( i, vEnergy[i], vRes[i] );
        iResult->SetPointError( i, 0., vResE[i] );
    }
    
    if( doKingFit )
    {
    
        // fill king function result graphs
        iResultKingSigma->Set( ( int )vKingEnergy.size() ) ;
        iResultKingGamma->Set( ( int )vKingEnergy.size() ) ;
        
        for( unsigned i = 0 ; i < vKingEnergy.size() ; i++ )
        {
            iResultKingSigma->SetPoint( i, vKingEnergy[i], vKingSigma[     i] ) ;
            iResultKingSigma->SetPointError( i, 0.0           , vKingSigmaError[i] ) ;
            iResultKingGamma->SetPoint( i, vKingEnergy[i], vKingGamma[     i] ) ;
            iResultKingGamma->SetPointError( i, 0.0           , vKingGammaError[i] ) ;
        }
        
        iResultKingSigma->Write();
        iResultKingGamma->Write();
    }
    
    return hList;
}




double VInstrumentResponseFunctionData::getResolutionErrorfromToyMC( double i68, double iN )
{
    if( i68 < 0. || iN < 1. )
    {
        return 0.;
    }
    
    // number of times to run the experiment
    const int nRun = 100;
    
    // histogram with results from each experiment
    TH1D h68( "h68", "h68", 1000, 0., 1.5 );
    
    // histogram with angular differences
    TH1D hDiff( "hDiff", "", 1000, 0., 1.0 );
    
    // normal distribution
    TF1 f( "f", "gaus(0)", 0., 5. );
    f.SetParameter( 0, 1. );
    f.SetParameter( 1, 0. );
    // normalized to 2D distribution, see Minuit table 7.1
    f.SetParameter( 2, i68 / sqrt( 2.41 ) );
    
    double x = 0;
    double y = 0.;
    double q[] = { 0 };
    int nq = 1;
    double d[] = { 0.68 };
    
    for( int j = 0; j < nRun; j++ )
    {
        hDiff.Reset();
        for( int i = 0; i < iN; i++ )
        {
            x = f.GetRandom();
            y = f.GetRandom();
            
            hDiff.Fill( sqrt( x * x + y * y ) );
        }
        if( hDiff.GetEntries() > 0 )
        {
            hDiff.GetQuantiles( nq, q, d );
            h68.Fill( q[0] );
        }
    }
    
    return h68.GetRMS();
}


void VInstrumentResponseFunctionData::setData( double iZe, int iAz_bin, double iAz_min, double iAz_max, int iNoise, double iPedvars, double iIndex, double iXoff, double iYoff )
{
    fZe = iZe;
    fAz_bin = iAz_bin;
    fAz_min = iAz_min;
    fAz_max = iAz_max;
    fXoff = iXoff;
    fYoff = iYoff;
    fWobble = sqrt( fXoff * fXoff + fYoff * fYoff );
    fNoise = iNoise;
    fPedvars = iPedvars;
    fSpectralIndex = iIndex;
}

