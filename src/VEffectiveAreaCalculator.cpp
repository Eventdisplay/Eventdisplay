/*! \class VEffectiveAreaCalculator
 *  \brief calculate effective areas and energy spectra
 *
 *  How to add new histograms:
 *  - add a new enum depending on the histogram type
 *  (E_HIS1D, E_HIS1P, E_HIS2D)
 *  - add a call to newEffectiveAreaHistogram()
 *  - add an entry in the enum to string function getEffectiveAreaNamefromEnumInt()
 *  - add the corresponding filling of the histogram
 *
 */

#include "VEffectiveAreaCalculator.h"

/*!
 *  CALLED FOR CALCULATION OF EFFECTIVE AREAS
 *
 *  this constructor is called for FILLING of the effective area tree
 *  (calculation of effective areas)
 *
 */
VEffectiveAreaCalculator::VEffectiveAreaCalculator( VInstrumentResponseFunctionRunParameter* iRunPara,
        vector< VGammaHadronCuts* > icuts,
        TFile* iFile )
{
    fRunPara = iRunPara;
    if( !fRunPara )
    {
        cout << "VEffectiveAreaCalculator: no run parameters given" << endl;
        cout << "...exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    reset();
    fOutputFile = iFile;
    if( fOutputFile )
    {
        fOutputFile->cd();
    }
    
    // no effective area file present
    bNOFILE = true;
    
    // number of bins for histograms
    nbins = fRunPara->fEnergyAxisBins_log10;
    nbins_MC = fRunPara->fEnergyAxisBins_log10;
    
    // bin definition for 2D histograms (allows coarser binning in energy)
    fBiasBin       = fRunPara->fBiasBin;
    fhistoNEbins   = fRunPara->fhistoNEbins;
    fLogAngularBin = fRunPara->fLogAngularBin;
    fResponseMatricesEbinning = fRunPara->fResponseMatricesEbinning;
    
    // Important: changing this means probably that the values used in
    // VEffectiveAreaCalculatorMCHistograms have to be changed as well
    // this should not be changed
    fEnergyAxis_minimum_defaultValue = -2.;
    fEnergyAxis_maximum_defaultValue =  4.;
    fLogAngular_minimum_defaultValue = -4.;
    fLogAngular_maximum_defaultValue =  1.;
    
    // cuts
    fCuts = icuts;
    fZe.push_back( fRunPara->fze );
    fAreaRadius.push_back( fRunPara->fCoreScatterRadius );
    fScatterMode.push_back( fRunPara->fCoreScatterMode );
    for( unsigned int i = 0; i < fZe.size(); i++ )
    {
        fXWobble.push_back( 0. );
        fYWobble.push_back( 0. );
        fNoise.push_back( 0 );
        fPedVar.push_back( 0. );
    }
    setIgnoreEnergyReconstructionCuts( fRunPara->fIgnoreEnergyReconstructionQuality );
    setIsotropicArrivalDirections( fRunPara->fIsotropicArrivalDirections );
    setTelescopeTypeCuts( fRunPara->fTelescopeTypeCuts );
    setWobbleOffset( fRunPara->fXoff, fRunPara->fYoff );
    setNoiseLevel( fRunPara->fnoise, fRunPara->fpedvar );
    
    // spectral weighting class
    fSpectralWeight = new VSpectralWeight();
    setMonteCarloEnergyRange( fRunPara->fMCEnergy_min, fRunPara->fMCEnergy_max, TMath::Abs( fRunPara->fMCEnergy_index ) );
    
    // define output tree (all histograms are written to this tree)
    hisTreeList = new TList();
    // same list, but histograms only
    // (used to reset histograms)
    hisTreeListofHistograms = new TList();
    // list for temporary histograms
    hisVList    = new TList();
    
    // Gaussian function for approximating response matrix
    fGauss = new TF1( "fGauss", "gaus", -2.5, 2.5 );
    
    // these histograms are filled into the output tree
    char hname[400];
    char htitle[400];
    
    // define range and binning of energy axis of effective area plots
    cout << "histogram parameters (bins, log10(Emin), log10(Emax)): " << nbins;
    cout << ", " << fEnergyAxis_minimum_defaultValue << ", " << fEnergyAxis_maximum_defaultValue << endl;
    sprintf( htitle, "title" );
    
    newEffectiveAreaHistogram( "1D", E_Emc,
                               "energy spectrum",
                               "energy_{MC} [TeV]",
                               "entries",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EmcUW,
                               "energy spectrum (unweighted)",
                               "energy_{MC} [TeV]",
                               "entries",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_Ecut,
                               "energy spectrum, after cuts",
                               "energy_{MC} [TeV]",
                               "entries",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EcutUW,
                               "unweighted energy spectrum, after cuts",
                               "energy_{MC} [TeV]",
                               "entries (unweighted)",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EcutNoTh2,
                               "energy spectrum, no direction cut, after cuts",
                               "energy_{MC} [TeV]",
                               "entries (no theta2 cut)",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_Ecut500,
                               "energy spectrum, after cuts",
                               "energy_{MC} [TeV]",
                               "entries",
                               500, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EcutRec,
                               "energy spectrum, no direction cutRecs",
                               "energy_{rec} [TeV]",
                               "entries",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EcutRecUW,
                               "unweighted energy spectrum, after cutsRecs",
                               "energy_{rec} [TeV]",
                               "entries (unweighted)",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    newEffectiveAreaHistogram( "1D", E_EcutRecNoTh2,
                               "energy spectrum, no direction cut, no direction cut, after cuts",
                               "energy_{rec} [TeV]",
                               "entries (no theta2 cut)",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue );
    // weighted rate
    // (use CTA binning, 5 bins per decade)
    newEffectiveAreaHistogram( "1D", E_WeightedRate,
                               "weighted rates",
                               "energy_{rec} [TeV]",
                               "entries",
                               30, -2.9, 3.1 );
    // weighted rate
    // (finner binning, primarily used for VTS analysis)
    newEffectiveAreaHistogram( "1D", E_WeightedRate005,
                               "weighted rates",
                               "energy_{rec} [TeV]",
                               "entries",
                               120, -2.9, 3.1 );
                               
    // individual cuts
    vector< enum E_HIS1D > iCutName;
    iCutName.push_back( E_EcutTrigger );
    iCutName.push_back( E_EcutFiducialArea );
    iCutName.push_back( E_EcutStereoQuality );
    iCutName.push_back( E_EcutTelType );
    iCutName.push_back( E_EcutDirection );
    iCutName.push_back( E_EcutEnergyReconstruction );
    iCutName.push_back( E_EcutGammaHadron );
    for( unsigned int i = 0; i < iCutName.size(); i++ )
    {
        newEffectiveAreaHistogram( "1D", iCutName[i],
                                   "energy spectrum, cut selection",
                                   "energy_{MC} [TeV]",
                                   "entries",
                                   nbins, fEnergyAxis_minimum_defaultValue,
                                   fEnergyAxis_maximum_defaultValue );
    }
    
    sprintf( hname, "gEffAreaMC" );
    gEffAreaMC = new TGraphAsymmErrors( 1 );
    gEffAreaMC->SetName( hname );
    gEffAreaMC->SetTitle( "effective area vs E_{MC}" );
    hisTreeList->Add( gEffAreaMC );
    
    sprintf( hname, "gEffAreaRec" );
    gEffAreaRec = new TGraphAsymmErrors( 1 );
    gEffAreaRec->SetName( hname );
    gEffAreaRec->SetTitle( "effective area vs E_{rec}" );
    hisTreeList->Add( gEffAreaRec );
    
    sprintf( hname, "gEffAreaNoTh2MC" );
    gEffAreaNoTh2MC = new TGraphAsymmErrors( 1 );
    gEffAreaNoTh2MC->SetName( hname );
    gEffAreaNoTh2MC->SetTitle( "effective area vs E_{MC} (no direction cut)" );
    hisTreeList->Add( gEffAreaNoTh2MC );
    
    sprintf( hname, "gEffAreaNoTh2Rec" );
    gEffAreaNoTh2Rec = new TGraphAsymmErrors( 1 );
    gEffAreaNoTh2Rec->SetName( hname );
    gEffAreaNoTh2Rec->SetTitle( "effective area vs E_{rec} (no direction cut)" );
    hisTreeList->Add( gEffAreaNoTh2Rec );
    
    newEffectiveAreaHistogram( "1P", E_EmcSWeight,
                               "eMCeweight",
                               "log_{10} energy_{MC} [TeV]",
                               "spectral weight",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               -1, 0., 1.e12 );
    newEffectiveAreaHistogram( "1P", E_EsysMCRelative,
                               "energy reconstruction",
                               "log_{10} energy_{MC} [TeV]",
                               "energy bias (E_{rec}-E_{MC})/E_{MC}",
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               -1, -1000., 1000., "s" );
                               
    newEffectiveAreaHistogram( "2D", E_EsysMCRelativeRMS,
                               "energy reconstruction",
                               "energy_{MC} [TeV]",
                               "energy bias (E_{rec}-E_{MC})/E_{MC}",
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               1000, -5., 5., "" );
    newEffectiveAreaHistogram( "2D", E_EsysMCRelative2D,
                               "energy reconstruction",
                               "energy_{MC} [TeV]",
                               "energy bias E_{rec}/E_{MC}",
                               300, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               fBiasBin, 0., 3., "" );
    newEffectiveAreaHistogram( "2D", E_EsysMCRelative2DNoDirectionCut,
                               "energy reconstruction, after gamma-selection cuts",
                               "energy_{MC} [TeV]",
                               "energy bias E_{rec}/E_{MC}",
                               300, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               fBiasBin, 0., 3., "" );
    newEffectiveAreaHistogram( "2D", E_Esys2D,
                               "energy reconstruction",
                               "energy_{MC} [TeV]",
                               "log_{10} E_{rec} - log_{10} E_{MC}",
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               100, -0.98, 2.02, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrix,
                               "migration matrix",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrixFine,
                               "migration matrix, fine binning",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               fResponseMatricesEbinning, -2.3, 2.7,
                               fResponseMatricesEbinning, -2.3, 2.7, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrixQC,
                               "migration matrix, after quality cuts",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               nbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrixFineQC,
                               "migration matrix, fine binning",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               fResponseMatricesEbinning, -2., 2.7,
                               fResponseMatricesEbinning, -2.3, 2.7, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrixNoDirectionCut,
                               "migration matrix",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue,
                               fhistoNEbins, fEnergyAxis_minimum_defaultValue,
                               fEnergyAxis_maximum_defaultValue, "" );
    newEffectiveAreaHistogram( "2D", E_ResponseMatrixFineNoDirectionCut,
                               "migration matrix, fine binning",
                               "energy_{rec} [TeV]",
                               "energy_{MC} [TeV]",
                               fResponseMatricesEbinning, -2., 2.7,
                               fResponseMatricesEbinning, -2.3, 2.7, "" );
                               
    // log angular difference histogram (vs true energy)
    sprintf( hname, "hAngularLogDiffEmc_2D" );
    hAngularLogDiffEmc_2D = new TH2D( hname, "log angular difference histogram (vs true energy)",
                                      25, -1.9, 3.5,
                                      100., -4., 1. );
    hAngularLogDiffEmc_2D->SetXTitle( "energy_{MC} [TeV]" );
    hisTreeList->Add( hAngularLogDiffEmc_2D );
    hisTreeListofHistograms->Add( hAngularLogDiffEmc_2D );
    
    // angular resolution graphs
    for( unsigned int i = 0; i < fRunPara->fAzMin.size(); i++ )
    {
        fGraph_AngularResolution68p.push_back( 0 );
        fGraph_AngularResolution80p.push_back( 0 );
        fGraph_AngularResolutionKingSigma.push_back( 0 );
        fGraph_AngularResolutionKingGamma.push_back( 0 );
        hVAngularLogDiffEmc_2D.push_back( 0 );
    }
    
    
    // tree with results from effective area calculation
    fEffArea = new TTree( "fEffArea", "effective area values" );
    fEffArea->Branch( "ze", &ze, "ze/F" );
    fEffArea->Branch( "az", &fAzBin, "az/I" );
    fEffArea->Branch( "azMin", &fMinAz, "azMin/F" );
    fEffArea->Branch( "azMax", &fMaxAz, "azMax/F" );
    if( !fRunPara->fEffArea_short_writing )
    {
        fEffArea->Branch( "Xoff", &fXoff, "Xoff/F" );
        fEffArea->Branch( "Yoff", &fYoff, "Yoff/F" );
    }
    fEffArea->Branch( "Woff", &fWoff, "Woff/F" );
    fEffArea->Branch( "noise", &fTNoise, "noise/I" );
    if( !fRunPara->fEffArea_short_writing )
    {
        fEffArea->Branch( "noisePE", &fTNoisePE, "noisePE/F" );
    }
    fEffArea->Branch( "pedvar", &fTPedvar, "pedvar/F" );
    fEffArea->Branch( "index", &fSpectralIndex, "index/F" );
    fEffArea->Branch( "nbins", &nbins, "nbins/I" );
    fEffArea->Branch( "e0", e0, "e0[nbins]/F" );  // log10( energy ) in [TeV]
    fEffArea->Branch( "eff", eff, "eff[nbins]/F" ); // effective area vs MC energy
    fEffArea->Branch( "eff_error", eff_error, "eff_error[nbins]/F" );
    fEffArea->Branch( "esys_rel", esys_rel, "esys_rel[nbins]/F" );
    fEffArea->Branch( "effNoTh2", effNoTh2, "effNoTh2[nbins]/F" );
    fEffArea->Branch( "effNoTh2_error", effNoTh2_error, "effNoTh2_error[nbins]/F" );
    fEffArea->Branch( "Rec_eff", Rec_eff, "Rec_eff[nbins]/F" ); // effective area vs reconstructed energy (approximation)
    fEffArea->Branch( "Rec_eff_error", Rec_eff_error, "Rec_eff_error[nbins]/F" );
    fEffArea->Branch( "Rec_effNoTh2", Rec_effNoTh2, "Rec_effNoTh2[nbins]/F" );
    fEffArea->Branch( "Rec_effNoTh2_error", Rec_effNoTh2_error, "Rec_effNoTh2_error[nbins]/F" );
    fEffArea->Branch( "Rec_angRes_p68", Rec_angRes_p68, "Rec_angRes_p68[nbins]/F" );
    fEffArea->Branch( "Rec_angRes_p80", Rec_angRes_p80, "Rec_angRes_p80[nbins]/F" );
    fEffArea->Branch( "Rec_angRes_kingSigma", Rec_angRes_kingSigma, "Rec_angRes_kingSigma[nbins]/F" );
    fEffArea->Branch( "Rec_angRes_kingGamma", Rec_angRes_kingGamma, "Rec_angRes_kingGamma[nbins]/F" );
    if( !fRunPara->fEffArea_short_writing )
    {
        fEffArea->Branch( hisTreeList, 64000, 1 );
    }
    // For reconstructing the response matrices
    fEffArea->Branch( "nbins_MC_Res", &nbins_MC_Res, "nbins_MC_Res/I" );
    fEffArea->Branch( "e_MC_Res", e_MC_Res, "e_MC_Res[nbins_MC_Res]/F" );
    fEffArea->Branch( "e_Rec_Res", e_Rec_Res, "e_Rec_Res[nbins_MC_Res]/F" );
    fEffArea->Branch( "e_Rec_Res_Err", e_Rec_Res_Err, "e_Rec_Res_Err[nbins_MC_Res]/F" );
    fEffArea->SetMarkerStyle( 20 );
    
    fAcceptance_AfterCuts_tree = new TTree( "Acceptance_AfterCuts" , "Info to construct background map" );
    fAcceptance_AfterCuts_tree->Branch( "Xoff_aC", &fXoff_aC, "Xoff_aC/F" );
    fAcceptance_AfterCuts_tree->Branch( "Yoff_aC", &fYoff_aC, "Yoff_aC/F" );
    fAcceptance_AfterCuts_tree->Branch( "Xoff_derot_aC", &fXoff_derot_aC, "Xoff_derot_aC/F" );
    fAcceptance_AfterCuts_tree->Branch( "Yoff_derot_aC", &fYoff_derot_aC, "Yoff_derot_aC/F" );
    fAcceptance_AfterCuts_tree->Branch( "Erec", &fErec, "Erec/F" );
    fAcceptance_AfterCuts_tree->Branch( "EMC", &fEMC, "EMC/F" );
    fAcceptance_AfterCuts_tree->Branch( "CRweight", &fCRweight, "CRweight/F" );
    
    fsolid_angle_norm_done = false;
    fsolid_angle_norm = 1.;
    
    ////////////////////////////////////
    // tree with DL2 event information (DL2 event) for each event
    fDL2_extendedTrees = false;
    // note adaptions to save disk space
    fDL2_runNumber = 0.;
    fDL2_eventNumber = 0;
    fDL2_MCaz = 0.;
    fDL2_MCel = 0.;
    fDL2_MCe0 = 0.;
    fDL2_MCxoff = 0.;
    fDL2_MCyoff = 0.;
    fDL2_ArrayPointing_Elevation = 0.;
    fDL2_ArrayPointing_Azimuth = 0.;
    fDL2_az = 0.;
    fDL2_el = 0.;
    fDL2_xoff = 0.;
    fDL2_yoff = 0.;
    fDL2_erec = 0.;
    fDL2_nimages = 0;
    fDL2_Cut_Class = 0;
    fDL2_Cut_MVA = 0.;
    if( fRunPara->fWriteEventdatatrees != "FALSE" )
    {
        fDL2EventTree = new TTree( "DL2EventTree", "DL2 tree" );
        fDL2EventTree->Branch( "runNumber", &fDL2_runNumber, "runNumber/i" );
        fDL2EventTree->Branch( "MCaz", &fDL2_MCaz, "MCaz/F" );
        fDL2EventTree->Branch( "MCel", &fDL2_MCel, "MCel/F" );
        fDL2EventTree->Branch( "MCe0", &fDL2_MCe0, "MCe0/F" );
        fDL2EventTree->Branch( "ArrayPointing_Azimuth", &fDL2_ArrayPointing_Azimuth, "ArrayPointing_Azimuth/F" );
        fDL2EventTree->Branch( "ArrayPointing_Elevation", &fDL2_ArrayPointing_Elevation, "ArrayPointing_Elevation/F" );
        fDL2EventTree->Branch( "az", &fDL2_az, "az/F" );
        fDL2EventTree->Branch( "el", &fDL2_el, "el/F" );
        // save disk space: don't write Xoff, Yoff (switched on for debugging)
        fDL2EventTree->Branch( "xoff", &fDL2_xoff, "xoff/F" );
        fDL2EventTree->Branch( "yoff", &fDL2_yoff, "yoff/F" );
        fDL2EventTree->Branch( "erec", &fDL2_erec, "erec/F" );
        fDL2EventTree->Branch( "nimages", &fDL2_nimages, "nimages/b" );
        fDL2EventTree->Branch( "CutClass", &fDL2_Cut_Class, "Class/b" );
        fDL2EventTree->Branch( "MVA", &fDL2_Cut_MVA, "MVA/F" );
        
        // extended trees
        if( fDL2_extendedTrees )
        {
            fDL2EventTree->Branch( "eventNumber", &fDL2_eventNumber, "eventNumber/i" );
            fDL2EventTree->Branch( "MCxoff", &fDL2_MCxoff, "MCxoff/F" );
            fDL2EventTree->Branch( "MCyoff", &fDL2_MCyoff, "MCyoff/F" );
            fDL2EventTree->Branch( "Xcore", &fDL2_Xcore, "Xcore/F" );
            fDL2EventTree->Branch( "Ycore", &fDL2_Ycore, "Ycore/F" );
            fDL2EventTree->Branch( "Xoff_intersect", &fDL2_Xoff_intersect, "Xoff_intersect/F" );
            fDL2EventTree->Branch( "Yoff_intersect", &fDL2_Yoff_intersect, "Yoff_intersect/F" );
            fDL2EventTree->Branch( "img2_ang", &fDL2_img2_ang, "img2_ang/F" );
            fDL2EventTree->Branch( "SizeSecondMax", &fDL2_SizeSecondMax, "SizeSecondMax/F" );
            fDL2EventTree->Branch( "NTelPairs", &fDL2_NTelPairs, "NTelPairs/i" );
            fDL2EventTree->Branch( "MSCW", &fDL2_MSCW, "MSCW/F" );
            fDL2EventTree->Branch( "MSCL", &fDL2_MSCL, "MSCL/F" );
            fDL2EventTree->Branch( "EmissionHeight", &fDL2_EmissionHeight, "EmissionHeight/F" );
            fDL2EventTree->Branch( "EmissionHeightChi2", &fDL2_EmissionHeightChi2, "EmissionHeightChi2/F" );
            fDL2EventTree->Branch( "size", fDL2_size, "size[nimages]/F" );
            fDL2EventTree->Branch( "dist", fDL2_dist, "dist[nimages]/F" );
            fDL2EventTree->Branch( "loss", fDL2_loss, "loss[nimages]/F" );
            fDL2EventTree->Branch( "fui", fDL2_fui, "fui[nimages]/F" );
            fDL2EventTree->Branch( "cross", fDL2_cross, "cross[nimages]/F" );
            fDL2EventTree->Branch( "asym", fDL2_asym, "asym[nimages]/F" );
            fDL2EventTree->Branch( "tgrad_x", fDL2_tgrad_x, "tgrad_x[nimages]/F" );
            fDL2EventTree->Branch( "R", &fDL2_R, "R[nimages]/F" );
            fDL2EventTree->Branch( "DispDiff", &fDL2_DispDiff, "DispDiff/F" );
            fDL2EventTree->Branch( "dESabs", &fDL2_dESabs, "dESabs/F" );
            fDL2EventTree->Branch( "NTrig", &fDL2_NTrig, "NTrig/i" );
            fDL2EventTree->Branch( "meanPedvar_Image", &fDL2_meanPedvar_Image, "meanPedvar_Image/F" );
            fDL2EventTree->Branch( "ES", fDL2_ES, "ES[nimages]/F" );
        }
    }
    else
    {
        fDL2EventTree = 0;
    }
    fDL2WriteFullEventTree = false;
    if( fRunPara->fWriteEventdatatrees == "FULLTREES" )
    {
        fDL2WriteFullEventTree = true;
    }
}

vector< TH1D* > VEffectiveAreaCalculator::initializeHistogramsVectorH1D( TH1D* h, string iName, unsigned int i )
{
    char hname[400];
    vector< TH1D* > iT_TH1D;
    
    for( unsigned int j = 0; j < fVMinAz.size(); j++ )
    {
        sprintf( hname, "%s_%d_%d", iName.c_str(), i, j );
        if( h )
        {
            iT_TH1D.push_back( ( TH1D* )h->Clone( hname ) );
            if( hisVList )
            {
                hisVList->Add( iT_TH1D.back() );
            }
        }
        else
        {
            iT_TH1D.push_back( 0 );
        }
    }
    return iT_TH1D;
}

vector< TH2D* > VEffectiveAreaCalculator::initializeHistogramsVectorH2D( TH2D* h, string iName, unsigned int i )
{
    char hname[400];
    vector< TH2D* > iT_TH2D;
    
    for( unsigned int j = 0; j < fVMinAz.size(); j++ )
    {
        sprintf( hname, "%s_%d_%d", iName.c_str(), i, j );
        if( h )
        {
            iT_TH2D.push_back( ( TH2D* )h->Clone( hname ) );
            if( hisVList )
            {
                hisVList->Add( iT_TH2D.back() );
            }
        }
        else
        {
            iT_TH2D.push_back( 0 );
        }
    }
    return iT_TH2D;
}

vector< TProfile* > VEffectiveAreaCalculator::initializeHistogramsVectorHProfile( TProfile* h, string iName, unsigned int i )
{
    char hname[400];
    vector< TProfile* > iT_TProfile;
    
    for( unsigned int j = 0; j < fVMinAz.size(); j++ )
    {
        sprintf( hname, "%s_%u_%u", iName.c_str(), i, j );
        if( h )
        {
            iT_TProfile.push_back( ( TProfile* )h->Clone( hname ) );
            if( hisVList )
            {
                hisVList->Add( iT_TProfile.back() );
            }
        }
        else
        {
            iT_TProfile.push_back( 0 );
        }
    }
    return iT_TProfile;
}



/*

   set up vectors with histograms for azimuth and spectral index bins

*/
void VEffectiveAreaCalculator::initializeHistograms( vector< double > iAzMin, vector< double > iAzMax, vector< double > iSpectralIndex )
{
    fVMinAz = iAzMin;
    fVMaxAz = iAzMax;
    fVSpectralIndex = iSpectralIndex;
    
    // temporary histograms for effective area calculation
    for( unsigned int i = 0; i < fVSpectralIndex.size(); i++ )
    {
        map< int, TH1D* >::iterator h_HIS1D_iterator;
        for( h_HIS1D_iterator = h_HIS1D.begin();
                h_HIS1D_iterator !=  h_HIS1D.end();
                h_HIS1D_iterator++ )
        {
            hV_HIS1D[h_HIS1D_iterator->first].push_back(
                initializeHistogramsVectorH1D(
                    h_HIS1D_iterator->second,
                    getEffectiveAreaNamefromEnumInt( h_HIS1D_iterator->first, "1D" ) + "V",
                    i ) );
        }
        map< int, TProfile* >::iterator h_HIS1P_iterator;
        for( h_HIS1P_iterator = h_HIS1P.begin();
                h_HIS1P_iterator !=  h_HIS1P.end();
                h_HIS1P_iterator++ )
        {
            hV_HIS1P[h_HIS1P_iterator->first].push_back(
                initializeHistogramsVectorHProfile(
                    h_HIS1P_iterator->second,
                    getEffectiveAreaNamefromEnumInt( h_HIS1P_iterator->first, "1P" ) + "V",
                    i ) );
        }
        map< int, TH2D* >::iterator h_HIS2D_iterator;
        for( h_HIS2D_iterator = h_HIS2D.begin();
                h_HIS2D_iterator !=  h_HIS2D.end();
                h_HIS2D_iterator++ )
        {
            hV_HIS2D[h_HIS2D_iterator->first].push_back(
                initializeHistogramsVectorH2D(
                    h_HIS2D_iterator->second,
                    getEffectiveAreaNamefromEnumInt( h_HIS2D_iterator->first, "2D" ) + "V",
                    i ) );
        }
    }
}


/*!
 *
 *  CALLED TO USE (READ) EFFECTIVE AREAS
 *
 *  this constructor is called to GET the effective area from a file and calculate energy spectra
 *
 *  called from anasum
 *
 */
VEffectiveAreaCalculator::VEffectiveAreaCalculator( string iInputFile, double azmin, double azmax, double ipedvar,
        double iSpectralIndex, vector< double > iMCZe,
        int iSmoothIter, double iSmoothThreshold, int iEffectiveAreaVsEnergyMC, bool iLikelihoodAnalysis, bool iIsOn )
{
    reset();
    
    bLikelihoodAnalysis = iLikelihoodAnalysis;
    bIsOn = iIsOn;
    fRunPara = 0;
    
    // MC intervals
    fMCZe = iMCZe;
    
    // no effective area file present
    bNOFILE = false;
    // no input file available, return always 1 for all corrections
    if( iInputFile == "NOFILE" )
    {
        bNOFILE = true;
        return;
    }
    // effective area smoothing
    fSmoothIter      = iSmoothIter;
    fSmoothThreshold = iSmoothThreshold;
    
    // no weighting
    fSpectralWeight = 0;
    
    // effective areas vs E_MC or E_rec
    fEffectiveAreaVsEnergyMC = iEffectiveAreaVsEnergyMC;
    
    // mean effective area
    gMeanEffectiveArea = new TGraphAsymmErrors( 1 );
    gMeanEffectiveArea->SetName( "gMeanEffectiveArea" );
    gMeanEffectiveArea->SetMarkerStyle( 20 );
    gMeanEffectiveArea->SetMarkerColor( fEffectiveAreaVsEnergyMC + 1 );
    gMeanEffectiveArea->SetLineColor( fEffectiveAreaVsEnergyMC + 1 );
    
    gMeanEffectiveAreaMC = new TGraphAsymmErrors( 1 );
    gMeanEffectiveAreaMC->SetName( "gMeanEffectiveAreaMC" );
    gMeanEffectiveAreaMC->SetMarkerStyle( 21 );
    gMeanEffectiveAreaMC->SetMarkerColor( fEffectiveAreaVsEnergyMC + 2 );
    gMeanEffectiveAreaMC->SetLineColor( fEffectiveAreaVsEnergyMC + 2 );
    
    // mean effective area for Time BINS
    
    gTimeBinnedMeanEffectiveArea = new TGraph2DErrors( 1 );
    gTimeBinnedMeanEffectiveArea->SetName( "gTimeBinnedMeanEffectiveArea" );
    gTimeBinnedMeanEffectiveArea->SetMarkerStyle( 20 );
    gTimeBinnedMeanEffectiveArea->SetMarkerColor( fEffectiveAreaVsEnergyMC + 1 );
    gTimeBinnedMeanEffectiveArea->SetLineColor( fEffectiveAreaVsEnergyMC + 1 );
    
    for( int i = 0; i < gTimeBinnedMeanEffectiveArea->GetN(); i++ )
    {
        gTimeBinnedMeanEffectiveArea->SetPoint( i, 0., 0., 0. );
        gTimeBinnedMeanEffectiveArea->SetPointError( i, 0., 0., 0. );
    }
    
    
    // current directory
    fGDirectory = gDirectory;
    
    // test if input file with effective areas exists
    if( iInputFile.find( "IGNOREEFFECTIVEAREA" ) == string::npos )
    {
        iInputFile = VUtilities::testFileLocation( iInputFile, "EffectiveAreas", true );
    }
    if( iInputFile.size() == 0 )
    {
        cout << "VEffectiveAreaCalculator::VEffectiveAreaCalculator error: input file string of zero length" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else if( iInputFile.find( "IGNOREEFFECTIVEAREA" ) != string::npos )
    {
        cout << "ignoring effective areas: ";
        cout << "all energy spectra will be invalid" << endl;
        bNOFILE = true;
    }
    else
    {
        TFile fIn( iInputFile.c_str() );
        if( fIn.IsZombie() )
        {
            cout << "VEffectiveAreaCalculator::VEffectiveAreaCalculator error opening file with effective areas: " << iInputFile << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        cout << "\t reading effective areas from " << fIn.GetName() << endl;
        
        if( !initializeEffectiveAreasFromHistograms( ( TTree* )gDirectory->Get( "fEffArea" ), ( TH1D* )gDirectory->Get( "hEmc" ), azmin, azmax, iSpectralIndex, ipedvar ) )
        {
            cout << "VEffectiveAreaCalculator ERROR: no effective areas found" << endl;
            cout << "all energy spectra will be invalid" << endl;
            bNOFILE = true;
        }
        fIn.Close();
        
        if( fGDirectory )
        {
            fGDirectory->cd();
        }
    }
}

/*
 * get mean value between to azimuth angle
 * taking into account the Eventdisplay azimuth
 * definition
 *
 * handles mix of neg/positive azimuth values
 *
 */
float VEffectiveAreaCalculator::getAzMean( float azmin, float azmax )
{
    // mean azimuth angle
    float iAzMean = 0.;
    if( azmin > 120. && azmax < -120. )
    {
        azmax += 360.;
    }
    else if( azmin < -150. && azmax > 120. )
    {
        azmin += 360.;
    }
    
    iAzMean = 0.5 * ( azmin + azmax );
    if( iAzMean > 180. )
    {
        iAzMean -= 360.;
    }
    
    return iAzMean;
}

/*
 * interpolate between two vectors
 *
 * (with or without 1/cos)
 */
vector< float > VEffectiveAreaCalculator::interpolate_effectiveArea( double iV, double iVLower, double iVupper,
        vector< float > iElower, vector< float > iEupper, bool iCos )
{
    vector< float > i_temp;
    if( iElower.size() == iEupper.size() )
    {
        i_temp.assign( iElower.size(), 0. );
        for( unsigned int i = 0; i < iElower.size(); i++ )
        {
            i_temp[i] = VStatistics::interpolate( iElower[i], iVLower, iEupper[i], iVupper, iV, iCos, 0.5, -90. );
        }
        return i_temp;
    }
    
    return i_temp;
}

/*
 *
 *  CALLED TO USE (READ) EFFECTIVE AREAS
 *
 *  called from anasum
 *
 */
bool VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms( TTree* iEffArea, TH1D* i_hEMC,
        double azmin, double azmax,
        double iSpectralIndex, double ipedvar )
{
    if( !iEffArea )
    {
        return false;
    }
    if( !i_hEMC )
    {
        cout << "----- Warning -----" << endl;
        cout << "  no MC histogram found to determine energy binning " << endl;
        cout << "   assume default binning: ";
        cout << "   " << fEnergyAxis_minimum_defaultValue << " " << fEnergyAxis_maximum_defaultValue;
        cout << " " << nbins;
        cout << endl;
        cout << "---- End of Warning ----" << endl;
        i_hEMC = new TH1D( "hEmc", "", nbins, fEnergyAxis_minimum_defaultValue, fEnergyAxis_maximum_defaultValue );
    }
    if( iEffArea->GetEntries() == 0 )
    {
        cout << "VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms: empty effective area tree" << endl;
        return false;
    }
    
    // mean azimuth angle
    float iAzMean = getAzMean( azmin, azmax );
    
    ////////////////////////////
    // define input tree
    float TazMin = 0.;
    float TazMax = 0.;
    float Tpedvar = 1.;
    iEffArea->SetBranchAddress( "azMin", &TazMin );
    iEffArea->SetBranchAddress( "azMax", &TazMax );
    iEffArea->SetBranchAddress( "pedvar", &Tpedvar );
    iEffArea->SetBranchAddress( "index", &fSpectralIndex );
    iEffArea->SetBranchAddress( "ze", &ze );
    iEffArea->SetBranchAddress( "Woff", &fWoff );
    // MC binning required for systematic error estimation
    iEffArea->SetBranchAddress( "e0", e0 );
    iEffArea->SetBranchAddress( "nbins", &nbins );
    // effective areas vs true energy
    if( fEffectiveAreaVsEnergyMC == 0 )
    {
        iEffArea->SetBranchAddress( "eff", eff );
    }
    // effective area vs reconstruction energy
    // this method should be used for the correction unfolding method
    else if( fEffectiveAreaVsEnergyMC == 1 )
    {
        iEffArea->SetBranchAddress( "Rec_eff", eff );
    }
    // Binned likelihood analysis requires
    // MC effective areas and response matrix
    // Getting MC eff
    if( iEffArea->GetBranchStatus( "nbins_MC_Res" ) )
    {
        iEffArea->SetBranchAddress( "nbins_MC_Res", &nbins_MC_Res );
        // Response Matrix
        iEffArea->SetBranchAddress( "e_MC_Res" , e_MC_Res );
        iEffArea->SetBranchAddress( "e_Rec_Res" , e_Rec_Res );
        iEffArea->SetBranchAddress( "e_Rec_Res_Err" , e_Rec_Res_Err );
    }
    else
    {
        nbins_MC_Res = 0;
    }
    // MC effective areas
    iEffArea->SetBranchAddress( "eff", eff_MC );
    // bias in energy reconstruction
    iEffArea->SetBranchAddress( "esys_rel", esys_rel );
    
    ////////////////////////////////////////////////////////////////////////////////////
    // prepare the energy vectors
    // (binning should be the same for all entries in the effective area tree)
    ////////////////////////////////////////////////////////////////////////////////////
    
    iEffArea->GetEntry( 0 );
    
    if( !i_hEMC )
    {
        cout << "VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms error: no effective area histogram found" << endl;
        return false;
    }
    
    fEff_E0.clear();
    for( int b = 1; b <= i_hEMC->GetNbinsX(); b++ )
    {
        fEff_E0.push_back( i_hEMC->GetBinCenter( b ) );
    }
    
    fVTimeBinnedMeanEffectiveArea.assign( i_hEMC->GetNbinsX(), 0. );
    // temporary vectors filled into the effective area maps later
    // effective areas (energy axis depend on fEffectiveAreaVsEnergyMC)
    vector< float > i_temp_Eff( i_hEMC->GetNbinsX(), 0. );
    vector< float > i_temp_Eff_MC( i_hEMC->GetNbinsX(), 0. );
    // bias in energy reconstruction (energy axis is in E_true)
    vector< float > i_temp_Esys( i_hEMC->GetNbinsX(), 0. );
    
    // temp vectors for binned likelihood analysis
    // used to fill maps
    vector< float > i_e_MC_Res( nbins_MC_Res, 0. );
    vector< float > i_e_Rec_Res( nbins_MC_Res, 0. );
    vector< float > i_e_Rec_Res_Err( nbins_MC_Res, 0. );
    
    cout << "\t selecting effective areas for mean az " << iAzMean << " deg, spectral index ";
    cout << iSpectralIndex << ", noise level " << ipedvar << endl;
    cout << "\t\ttotal number of curves: " << iEffArea->GetEntries();
    cout << ", total number of bins on energy axis: " << fEff_E0.size() << endl;
    
    fNTimeBinnedMeanEffectiveArea = 0;
    fNTimeBinnedMeanEffectiveAreaMC = 0;
    if( gTimeBinnedMeanEffectiveArea )
    {
        for( int i = 0; i < gTimeBinnedMeanEffectiveArea->GetN(); i++ )
        {
            gTimeBinnedMeanEffectiveArea->SetPoint( i, 0., 0., 0. );
            gTimeBinnedMeanEffectiveArea->SetPointError( i, 0., 0., 0. );
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    // fill ze, woff, noise, index vectors
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    bool i_ze_F = false;
    unsigned int i_index_ze = 0;
    bool i_woff_F = false;
    unsigned int i_index_woff = 0;
    bool i_noise_F = false;
    unsigned int i_index_noise = 0;
    bool i_index_F = false;
    unsigned int i_index_index = 0;
    int iIndexAz = 0;
    double iInvMean = 0.;
    double iInvMax = 1.e5;
    
    bool fLotsOfPrintOuts = false;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // loop over all entries in effective area tree
    // (not sure if this is really necessary, in the end a few entries are only needed)
    for( int i = 0; i < iEffArea->GetEntries(); i++ )
    {
        iEffArea->GetEntry( i );
        
        // check binning of effective areas
        if( i_hEMC && ( int )i_hEMC->GetNbinsX() != ( int )fEff_E0.size() )
        {
            cout << "VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms error: effective area curve with different binning";
            cout << " " << i_hEMC->GetNbinsX() << " " << fEff_E0.size() << endl;
            cout << "abort reading effective area tree..." << endl;
            return false;
        }
        ///////////////////////////////////////////////////
        // check the azimuth range
        ///////////////////////////////////////////////////
        
        // expect a couple of az bins and then a last bin with the average azimuth bin; typically [-1000., 1000.]
        // this is the sign to read effective areas
        if( fabs( TazMin ) > 5.e2 || fabs( TazMax ) > 5.e2 )
        {
            iInvMax = 1.e5;
            iEffArea->GetEntry( iIndexAz );
            
            ///////////////////////////////////////////////////
            // zenith angle
            ///////////////////////////////////////////////////
            i_ze_F = false;
            for( unsigned z = 0; z < fZe.size(); z++ )
            {
                // this has to be a relatively large value due
                // to allowed wobble offsets up to 2.0 deg
                if( fabs( fZe[z] - ze ) < 2.0 )
                {
                    i_index_ze = z;
                    i_ze_F = true;
                    break;
                }
            }
            if( !i_ze_F )
            {
                for( unsigned int w = 0; w < fMCZe.size(); w++ )
                {
                    if( fabs( fMCZe[w] - ze ) < 2.0 )
                    {
                        fZe.push_back( fMCZe[w] );
                        break;
                    }
                }
            }
            ///////////////////////////////////////////////////
            // wobble offset
            ///////////////////////////////////////////////////
            if( i_index_ze < fEff_WobbleOffsets.size() )
            {
                i_woff_F = false;
                for( unsigned z = 0; z < fEff_WobbleOffsets[i_index_ze].size(); z++ )
                {
                    if( fabs( fEff_WobbleOffsets[i_index_ze][z] - fWoff ) < 0.05 )
                    {
                        i_index_woff = z;
                        i_woff_F = true;
                        break;
                    }
                }
                if( !i_woff_F )
                {
                    fEff_WobbleOffsets[i_index_ze].push_back( fWoff );
                    i_index_woff = fEff_WobbleOffsets[i_index_ze].size() - 1;
                }
            }
            else
            {
                vector< double > itemp;
                itemp.push_back( fWoff );
                fEff_WobbleOffsets.push_back( itemp );
                if( i_index_ze < fEff_WobbleOffsets.size() )
                {
                    i_index_woff = fEff_WobbleOffsets[i_index_ze].size() - 1;
                }
                else
                {
                    cout << "Error in ";
                    cout << "VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms:";
                    cout << " setting of wobble vector" << endl;
                    cout << i_index_ze << " >= " << fEff_WobbleOffsets.size() << endl;
                    exit( EXIT_FAILURE );
                }
                
            }
            ///////////////////////////////////////////////////
            // noise level
            ///////////////////////////////////////////////////
            if( i_index_ze < fEff_Noise.size() )
            {
                if( i_index_woff < fEff_Noise[i_index_ze].size() )
                {
                    i_noise_F = false;
                    for( unsigned w = 0; w < fEff_Noise[i_index_ze][i_index_woff].size(); w++ )
                    {
                        if( fabs( fEff_Noise[i_index_ze][i_index_woff][w] - Tpedvar ) < 0.005 )
                        {
                            i_index_noise = w;
                            i_noise_F = true;
                            break;
                        }
                    }
                    if( !i_noise_F )
                    {
                        fEff_Noise[i_index_ze][i_index_woff].push_back( Tpedvar );
                        i_index_noise = fEff_Noise[i_index_ze][i_index_woff].size() - 1;
                    }
                }
                else
                {
                    vector< double > itemp;
                    itemp.push_back( Tpedvar );
                    fEff_Noise[i_index_ze].push_back( itemp );
                    i_index_noise = fEff_Noise[i_index_ze][i_index_woff].size() - 1;
                }
            }
            else
            {
                vector< double > itemp;
                itemp.push_back( Tpedvar );
                vector< vector< double > > iitemp;
                iitemp.push_back( itemp );
                fEff_Noise.push_back( iitemp );
                i_index_noise = fEff_Noise[i_index_ze][i_index_woff].size() - 1;
            }
            ///////////////////////////////////////////////////
            // spectral index
            ///////////////////////////////////////////////////
            if( i_index_ze < fEff_SpectralIndex.size() )
            {
                if( i_index_woff < fEff_SpectralIndex[i_index_ze].size() )
                {
                    if( i_index_noise < fEff_SpectralIndex[i_index_ze][i_index_woff].size() )
                    {
                        i_index_F = false;
                        for( unsigned s = 0; s < fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size(); s++ )
                        {
                            if( fabs( fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise][s] - fSpectralIndex ) < 0.005 )
                            {
                                i_index_index = s;
                                i_index_F = true;
                                break;
                            }
                        }
                        if( !i_index_F )
                        {
                            fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].push_back( fSpectralIndex );
                            i_index_index = fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size() - 1;
                        }
                    }
                    else
                    {
                        vector< double > itemp;
                        itemp.push_back( fSpectralIndex );
                        fEff_SpectralIndex[i_index_ze][i_index_woff].push_back( itemp );
                        i_index_index = fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size() - 1;
                    }
                }
                else
                {
                    vector< double > itemp;
                    itemp.push_back( fSpectralIndex );
                    vector< vector< double > > iitemp;
                    iitemp.push_back( itemp );
                    fEff_SpectralIndex[i_index_ze].push_back( iitemp );
                    i_index_index = fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size() - 1;
                }
            }
            else
            {
                vector< double > itemp;
                itemp.push_back( fSpectralIndex );
                vector< vector< double > > iitemp;
                iitemp.push_back( itemp );
                vector< vector< vector< double > > > iiitemp;
                iiitemp.push_back( iitemp );
                fEff_SpectralIndex.push_back( iiitemp );
                i_index_index = fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size() - 1;
            }
            unsigned int i_ID = i_index_index + 100 * ( i_index_noise + 100 * ( i_index_woff + 100 * i_index_ze ) );
            if( fLotsOfPrintOuts )
            {
                cout << i_index_ze << " " << i_index_woff << " " << i_index_noise << " " << i_index_index << "\t" << i_ID << endl;
            }
            
            ///////////////////////////////////////////////////
            // read effective area and load into maps
            ///////////////////////////////////////////////////
            if( nbins <= 0 )
            {
                cout << "WARNING : incomplete effective areas for id " << i_ID << endl;
                cout << i_index_index << "\t" << i_index_noise << "\t" << i_index_woff << "\t" << i_index_ze << endl;
                cout << "in bool VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms(";
                cout << "TTree *iEffArea, double azmin, double azmax, double iSpectralIndex )" << endl;
                cout << "Missing effective area:";
                if( i_index_ze < fZe.size() )
                {
                    cout << " ze = " << fZe[i_index_ze] << " [deg],";
                }
                if( i_index_ze < fEff_WobbleOffsets.size() && i_index_woff < fEff_WobbleOffsets[i_index_ze].size() )
                {
                    cout << " woff = " << fEff_WobbleOffsets[i_index_ze][i_index_woff] << " [deg]";
                }
                if( i_index_ze < fEff_Noise.size() &&
                        i_index_woff < fEff_Noise[i_index_ze].size() &&
                        i_index_noise < fEff_Noise[i_index_ze][i_index_woff].size() )
                {
                    cout << " noise = " << fEff_Noise[i_index_ze][i_index_woff][i_index_noise];
                }
                if( i_index_ze < fEff_SpectralIndex.size() && i_index_woff < fEff_SpectralIndex[i_index_ze].size() &&
                        i_index_noise < fEff_SpectralIndex[i_index_ze][i_index_woff].size() &&
                        i_index_index < fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise].size() )
                {
                    cout << " spectral index = " << fEff_SpectralIndex[i_index_ze][i_index_woff][i_index_noise][i_index_index];
                }
                cout << endl;
                cout << "please check if this falls into your parameter space" << endl;
            }
            // effective areas vs energy (decision on vs MC, rec energy etc. has been made earlier)
            for( unsigned int e = 0; e < fEff_E0.size(); e++ )
            {
                i_temp_Eff[e] = 0.;
                i_temp_Eff_MC[e] = 0.;
                i_temp_Esys[e] = 0.;
                for( int j = 0; j < nbins; j++ )
                {
                    if( TMath::Abs( e0[j] - fEff_E0[e] ) < 1.e-5 )
                    {
                        i_temp_Eff[e] = eff[j];
                        i_temp_Eff_MC[e] = eff_MC[j];
                        i_temp_Esys[e]  = esys_rel[j];
                    }
                }
            }
            fEffArea_map[i_ID]        = i_temp_Eff;
            fEffAreaMC_map[i_ID]      = i_temp_Eff_MC;
            fEff_EsysMCRelative[i_ID] = i_temp_Esys;
            
            // Setting Values for the generating the migration matrix
            i_e_MC_Res.resize( nbins_MC_Res );
            i_e_Rec_Res.resize( nbins_MC_Res );
            i_e_Rec_Res_Err.resize( nbins_MC_Res );
            
            for( int j = 0; j < nbins_MC_Res; j++ )
            {
                i_e_MC_Res[j] = e_MC_Res[j];
                i_e_Rec_Res[j] = e_Rec_Res[j];
            }
            // Assigning Key to Map
            fe_MC_Res_map[i_ID] = i_e_MC_Res;
            fe_Rec_Res_map[i_ID] = i_e_Rec_Res;
            fe_Rec_Res_Err_map[i_ID] = i_e_Rec_Res_Err;
            
            fEntry_map[i_ID] = i;
            
            // this is neeeded only if there are no azimuth dependent effective areas
            iIndexAz++;
            
            continue;
        }
        ///////////////////////////////////////////////////
        // test az bin
        ///////////////////////////////////////////////////
        // expect bin like [135,-135]
        iInvMean = getAzMean( TazMin, TazMax );
        if( TazMax < 0. && TazMin > 0. )
        {
            double iT =  fabs( iAzMean - iInvMean );
            if( iT > 180. )
            {
                iT = fabs( iT - 360. );
            }
            if( iT < iInvMax )
            {
                iInvMax = iT;
                iIndexAz = i;
            }
            if( iAzMean < 0. && iAzMean > TazMax )
            {
                continue;
            }
            else if( iAzMean > 0. && iAzMean < TazMin )
            {
                continue;
            }
        }
        else
        {
            if( fabs( iAzMean - iInvMean ) < iInvMax )
            {
                iInvMax = fabs( iAzMean - iInvMean );
                iIndexAz = i;
            }
        }
    }
    ///////////////////////////////////////////////////
    if( fZe.size() == 0 )
    {
        cout << "VEffectiveAreaCalculator::initializeEffectiveAreasFromHistograms error:";
        cout << " no effective areas found in effective area tree" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "\t use histograms for effective areas (number of zenith angle intervals: " << fZe.size() << "; ";
    for( unsigned int i = 0; i < fZe.size() - 1; i++ )
    {
        cout << fZe[i] << ", ";
    }
    if( fZe.size() > 0 )
    {
        cout << fZe[fZe.size() - 1];
    }
    cout << " deg)" << endl;
    if( fEffectiveAreaVsEnergyMC == 0 )
    {
        cout << "\t (effective area vs MC energy)" << endl;
    }
    else
    {
        cout << "\t (effective area vs reconstructed energy)" << endl;
    }
    
    if( fSmoothIter > 0 )
    {
        smoothEffectiveAreas( fEffArea_map );
    }
    
    ///////////////////////////////////////////////////
    if( fLotsOfPrintOuts )
    {
        cout << "ze size " << fZe.size() << endl;
        for( unsigned int z = 0; z < fZe.size(); z++ )
        {
            cout << "ze " << z << " " << fZe[z] << endl;
            if( z < fEff_WobbleOffsets.size() )
            {
                cout << "\t woff size " << fEff_WobbleOffsets[z].size() << endl;
                for( unsigned int w = 0; w < fEff_WobbleOffsets[z].size(); w++ )
                {
                    cout << "\t woff " << w << " " << fEff_WobbleOffsets[z][w] << endl;
                    if( z < fEff_Noise.size() && w < fEff_Noise[z].size() )
                    {
                        cout << "\t\t noise size " << fEff_Noise[z][w].size() << endl;
                        for( unsigned n = 0; n < fEff_Noise[z][w].size(); n++ )
                        {
                            cout << "\t\t noise " << n << " " << fEff_Noise[z][w][n] << endl;
                            if( z < fEff_SpectralIndex.size() && w < fEff_SpectralIndex[z].size() && n < fEff_SpectralIndex[z][w].size() )
                            {
                                cout << "\t\t\t index size " << fEff_SpectralIndex[z][w][n].size() << endl;
                                for( unsigned s = 0; s < fEff_SpectralIndex[z][w][n].size(); s++ )
                                {
                                    cout << "\t\t\t index " << s << " " << fEff_SpectralIndex[z][w][n][s] << endl;
                                    unsigned int i_ID = s + 100 * ( n + 100 * ( w + 100 * z ) );
                                    if( fEffArea_map.find( i_ID ) != fEffArea_map.end() )
                                    {
                                        cout << "\t\t\t map " << fEffArea_map[i_ID].size() << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ///////////////////////////////////////////////////
    
    return true;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

VEffectiveAreaCalculator::~VEffectiveAreaCalculator()
{
    if( gMeanEffectiveArea )
    {
        delete gMeanEffectiveArea;
    }
    if( gTimeBinnedMeanEffectiveArea )
    {
        delete gTimeBinnedMeanEffectiveArea;
    }
    if( gMeanSystematicErrorGraph )
    {
        delete gMeanSystematicErrorGraph;
    }
    if( gMeanEffectiveAreaMC )
    {
        delete gMeanEffectiveAreaMC;
    }
    if( hMeanResponseMatrix )
    {
        delete hMeanResponseMatrix;
    }
}


void VEffectiveAreaCalculator::reset()
{
    gMeanEffectiveArea = 0;
    gTimeBinnedMeanEffectiveArea = 0;
    gMeanEffectiveAreaMC = 0;
    hMeanResponseMatrix = 0;
    
    fMC_ScatterArea = 0.;
    
    bNOFILE = true;
    fGDirectory = 0;
    
    fSpectralIndex = 2.0;
    
    // Important: changing this means probably that the values used in
    // VEffectiveAreaCalculatorMCHistograms have to be changed as well
    fEnergyAxis_minimum_defaultValue = -2.0;
    fEnergyAxis_maximum_defaultValue = 4.0;
    
    fTNoise = 0;
    fTNoisePE = 0.;
    fTPedvar = 0.;
    
    fAzBin = 0;
    fMinAz = -1.e3;
    fMaxAz = 1.e3;
    fEffectiveAreaVsEnergyMC = 1;
    
    gEffAreaMC = 0;
    gEffAreaRec = 0;
    
    gEffAreaNoTh2MC = 0;
    gEffAreaNoTh2Rec = 0;
    
    hAngularLogDiffEmc_2D = 0;
    fEffArea = 0;
    hisTreeList = 0;
    hisTreeListofHistograms = 0;
    hisVList = 0;
    
    fEffArea = 0;
    ze = 0.;
    nbins = 60;
    for( int i = 0; i < VMAXBINS; i++ )
    {
        e0[i] = 0.;
        eff[i] = 0.;
        eff_error[i] = 0.;
        esys_rel[i] = 0.;
        effNoTh2[i] = 0.;
        effNoTh2_error[i] = 0.;
        Rec_eff[i] = 0.;
        Rec_eff_error[i] = 0.;
        Rec_effNoTh2[i] = 0.;
        Rec_effNoTh2_error[i] = 0.;
        Rec_angRes_p68[i] = 0.;
        Rec_angRes_p80[i] = 0.;
        Rec_angRes_kingSigma[i] = 0.;
        Rec_angRes_kingGamma[i] = 0.;
    }
    
    fEffectiveAreas_meanZe = 0.;
    fEffectiveAreas_meanWoff = 0.;
    fEffectiveAreas_meanPedVar = 0.;
    fEffectiveAreas_meanIndex = 0.;
    fEffectiveAreas_meanN = 0.;
    
    gMeanSystematicErrorGraph = 0;
    
    fXoff_aC = -99;
    fYoff_aC = -99;
    fXoff_derot_aC = -99;
    fYoff_derot_aC = -99;
    fErec = -99;
    fEMC = -99;
    fCRweight = -99;
    
    
}

/*
 * solid angle from MC simulations (cone)
 *
 */
double VEffectiveAreaCalculator::getMCSolidAngleNormalization()
{
    double iSolAngleNorm = 1.;
    if( fCuts.size() > 0 && fCuts[0] && fRunPara )
    {
        if( fRunPara->fViewcone_max > 0.
                && fCuts[0]->fCut_CameraFiducialSize_MC_max > 0.
                && fCuts[0]->fCut_CameraFiducialSize_MC_max < 1000.
                && fCuts[0]->fCut_CameraFiducialSize_MC_max < fRunPara->fViewcone_max )
        {
            // solid angle of simulated showers
            double iSN_mc = ( 1. - cos( fRunPara->fViewcone_max * TMath::DegToRad() ) );
            if( fRunPara->fViewcone_min > 0. )
            {
                iSN_mc -= ( 1. - cos( fRunPara->fViewcone_min * TMath::DegToRad() ) );
            }
            // solid angle of angular bin
            double iSN_cu = ( 1. - cos( fCuts[0]->fCut_CameraFiducialSize_MC_max * TMath::DegToRad() ) );
            if( fCuts[0]->fCut_CameraFiducialSize_MC_min > 0. )
            {
                iSN_cu -= ( 1. - cos( fCuts[0]->fCut_CameraFiducialSize_MC_min * TMath::DegToRad() ) );
            }
            
            if( iSN_mc > 0. )
            {
                iSolAngleNorm = iSN_cu / iSN_mc;
            }
        }
    }
    
    return iSolAngleNorm;
}

/*

   read spectra of MC events from mscw file (filled in eventdisplay)

*/
bool VEffectiveAreaCalculator::getMonteCarloSpectra( VEffectiveAreaCalculatorMCHistograms* iMC_histo )
{
    // get solid angle normalization
    // (important for calculation of effective ares vs wobble offsets using diffuse gamma-ray simulations)
    double iSolAngleNorm = getMCSolidAngleNormalization();
    cout << "VEffectiveAreaCalculator::getMonteCarloSpectra: solid angle normalization factor: " << iSolAngleNorm << endl;
    
    char hname[800];
    // loop over az bins
    // (MC az [-180., 180.])
    for( unsigned int i_az = 0; i_az < fVMinAz.size(); i_az++ )
    {
        // loop over all spectral index
        for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
        {
            if( s < hV_HIS1D[E_Emc].size() && i_az < hV_HIS1D[E_Emc][s].size() )
            {
                sprintf( hname, "hVEmc_%u_%u", s, i_az );
                if( iMC_histo->getHistogram_Emc( i_az, s ) )
                {
                    hV_HIS1D[E_Emc][s][i_az] = ( TH1D* )iMC_histo->getHistogram_Emc( i_az, s )->Clone( hname );
                    if( hV_HIS1D[E_Emc][s][i_az] )
                    {
                        hV_HIS1D[E_Emc][s][i_az]->Scale( iSolAngleNorm );
                    }
                    if( hV_HIS1D[E_Emc][s][i_az] && fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
                    {
                        hV_HIS1D[E_Emc][s][i_az]->Scale( ( 1.0 - fRunPara->fIgnoreFractionOfEvents ) );
                    }
                }
                else
                {
                    hV_HIS1D[E_Emc][s][i_az] = 0;
                }
            }
            if( s < hV_HIS1D[E_EmcUW].size() && i_az < hV_HIS1D[E_EmcUW][s].size() )
            {
                sprintf( hname, "hVEmcUW_%u_%u", s, i_az );
                if( iMC_histo->getHistogram_EmcUnweighted( s ) )
                {
                    hV_HIS1D[E_EmcUW][s][i_az] = ( TH1D* )iMC_histo->getHistogram_EmcUnweighted( s )->Clone( hname );
                }
                else
                {
                    hV_HIS1D[E_EmcUW][s][i_az] = 0;
                }
            }
            // profiles with spectral weights
            if( s < hV_HIS1P[E_EmcSWeight].size() && i_az < hV_HIS1P[E_EmcSWeight][s].size() )
            {
                sprintf( hname, "hVEmcSWeight_%u_%u", s, i_az );
                if( iMC_histo->getHistogram_EmcWeight( i_az, s ) )
                {
                    hV_HIS1P[E_EmcSWeight][s][i_az] = ( TProfile* )iMC_histo->getHistogram_EmcWeight( i_az, s )->Clone( hname );
                    if( hV_HIS1P[E_EmcSWeight][s][i_az] )
                    {
                        hV_HIS1P[E_EmcSWeight][s][i_az]->Scale( iSolAngleNorm );
                    }
                    if( hV_HIS1P[E_EmcSWeight][s][i_az] && fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
                    {
                        hV_HIS1P[E_EmcSWeight][s][i_az]->Scale( ( 1.0 - fRunPara->fIgnoreFractionOfEvents ) );
                    }
                }
                else
                {
                    hV_HIS1P[E_EmcSWeight][s][i_az] = 0;
                }
            }
        }
    }
    
    return true;
}

/*
 *
 *  CALLED FOR CALCULATION OF EFFECTIVE AREAS
 *
 */
bool VEffectiveAreaCalculator::fill( CData* d, VEffectiveAreaCalculatorMCHistograms* iMC_histo, unsigned int iMethod )
{
    bool bDebugCuts = false;          // lots of debug output
    
    // make sure that vectors are initialized
    unsigned int ize = 0;      // should always be zero
    if( ize >= fZe.size() )
    {
        cout << "VEffectiveAreaCalculator::fill error: vectors are not correctly initialized " << fZe.size() << "\t" << ize << endl;
        return false;
    }
    resetHistogramsVectors();
    
    // do not require successfull energy reconstruction
    if( fIgnoreEnergyReconstruction )
    {
        iMethod = 100;
    }
    
    //////////////////////////////////////////////////////////////////
    // total Monte Carlo core scatter area (depends on CORSIKA shower core scatter mode)
    fMC_ScatterArea = 0.;
    if( fScatterMode[ize] == "VIEWCONE" )
    {
        fMC_ScatterArea = fAreaRadius[ize] * fAreaRadius[ize] * TMath::Pi();
    }
    else if( fScatterMode[ize] == "FLAT" )
    {
        fMC_ScatterArea = fAreaRadius[ize] * fAreaRadius[ize] * TMath::Pi() * cos( fZe[ize] * TMath::DegToRad() );
    }
    else
    {
        cout << "VEffectiveAreaCalculator::fill ERROR: unknown CORSIKA scatter mode: " << fScatterMode[ize] << endl;
        return false;
    }
    // reset unique event counter
    //	fUniqueEventCounter.clear();
    int iSuccessfullEventStatistics = 0;
    
    //////////////////////////////////////////////////////////////////
    // print some run information
    cout << endl;
    cout << "calculating effective areas: " << endl;
    cout << "\t zenith angle: " << fZe[ize];
    cout << ", wobble offset (x,y): " << fXWobble[ize] << ", " << fYWobble[ize];
    cout << " (" << sqrt( fXWobble[ize]*fXWobble[ize] + fYWobble[ize]*fYWobble[ize] ) << " deg)" << endl;
    cout << "\t noise level: " << fNoise[ize] << " (pedvar: " << fPedVar[ize] << ")" << endl;
    cout << "\t area (" << fScatterMode[ize] << ") [m^2]: " << fMC_ScatterArea;
    cout << " (scatter radius " << fAreaRadius[ize] << " [m])" << endl;
    cout << "\t energy reconstruction method: " << iMethod << endl;
    if( fIsotropicArrivalDirections )
    {
        cout << "\t assuming isotropic arrival directions" << endl;
    }
    if( fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
    {
        cout << "\t ignore first " << fRunPara->fIgnoreFractionOfEvents * 100. << " % of events" << endl;
    }
    
    cout << endl;
    if( fSpectralWeight )
    {
        fSpectralWeight->print();
    }
    else
    {
        cout << "(no specral weight given)" << endl;
    }
    cout << endl;
    
    // make sure that all data pointers exist
    if( !d || !iMC_histo )
    {
        cout << "VEffectiveAreaCalculator::fill error: no data tree or MC histograms: " << endl;
        cout << "Data tree: " << d << ", MC histograms: " << iMC_histo << endl;
        return false;
    }
    
    // spectral weight
    double i_weight = 1.;
    // reconstructed energy (TeV, log10)
    double eRec = 0.;
    double eRecLin = 0.;
    // MC energy (TeV, log10)
    double eMC = 0.;
    
    ////////////////////////////////////////////////////////////////////////////
    // get MC histograms
    if( !getMonteCarloSpectra( iMC_histo ) )
    {
        cout << "VEffectiveAreaCalculator::fill error while getting MC spectra" << endl;
        return false;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // reset cut statistics
    for( unsigned int i = 0; i < fCuts.size(); i++ )
    {
        fCuts[i]->resetCutStatistics();
    }
    
    ///////////////////////////////////////////////////////
    // get full data set and loop over all entries
    ///////////////////////////////////////////////////////
    Long64_t d_nentries = d->fChain->GetEntries();
    Long64_t i_start = 0;
    if( fRunPara && fRunPara->fIgnoreFractionOfEvents > 0. )
    {
        i_start = ( Long64_t )( fRunPara->fIgnoreFractionOfEvents * d_nentries );
    }
    cout << "\t total number of data events: " << d_nentries << " (start at event " << i_start << ")" << endl;
    
    //--- for the CR normalisation filling Acceptance tree total number of simulated is needed
    //-- WARNING if the rule for the azimuth bin changes in VInstrumentResponseFunctionRunParameter the following line must be adapted!!!!
    unsigned int number_of_az_bin = fRunPara->fAzMin.size();
    // if no azimuth bin, all events are in bin 0. if azimuth bin, all event are in the last bin
    int az_bin_index = 0;
    if( number_of_az_bin > 0 )
    {
        az_bin_index = ( int ) number_of_az_bin - 1;
    }
    // loop over all events
    for( Long64_t i = i_start; i < d_nentries; i++ )
    {
        d->GetEntry( i );
        
        // update cut statistics
        VGammaHadronCuts* iAnaCuts = getGammaHadronCuts( d );
        if( !iAnaCuts )
        {
            continue;
        }
        iAnaCuts->newEvent();
        
        if( bDebugCuts )
        {
            cout << "============================== " << endl;
            cout << "EVENT entry number " << i << endl;
        }
        
        // apply MC cuts
        if( bDebugCuts )
        {
            cout << "#0 CUT MC " << iAnaCuts->applyMCXYoffCut( d->MCxoff, d->MCyoff, false ) << endl;
        }
        
        if( !iAnaCuts->applyMCXYoffCut( d->MCxoff, d->MCyoff, true ) )
        {
            fillDL2EventDataTree( d, ( UChar_t )VGammaHadronCutsStatistics::eMC_XYoff, -1. );
            continue;
        }
        
        // log of MC energy
        eMC = log10( d->MCe0 );
        
        // fill trigger cuts
        fillEcutSub( eMC, E_EcutTrigger );
        
        ////////////////////////////////
        // apply general quality and gamma/hadron separation cuts
        
        // apply reconstruction cuts
        if( bDebugCuts )
        {
            cout << "#1 CUT applyInsideFiducialAreaCut ";
            cout << iAnaCuts->applyInsideFiducialAreaCut();
            cout << "\t" << iAnaCuts->applyStereoQualityCuts( iMethod, false, i, true ) << endl;
        }
        
        // apply fiducial area cuts
        if( !iAnaCuts->applyInsideFiducialAreaCut( true ) )
        {
            fillDL2EventDataTree( d, 2, -1. );
            continue;
        }
        fillEcutSub( eMC, E_EcutFiducialArea );
        
        // apply reconstruction quality cuts
        if( !iAnaCuts->applyStereoQualityCuts( iMethod, true, i , true ) )
        {
            fillDL2EventDataTree( d, 3, -1. );
            continue;
        }
        fillEcutSub( eMC, E_EcutStereoQuality );
        
        // apply telescope type cut (e.g. for CTA simulations)
        if( fTelescopeTypeCutsSet )
        {
            if( bDebugCuts )
            {
                cout << "#2 Cut NTELType " << iAnaCuts->applyTelTypeTest( false ) << endl;
            }
            if( !iAnaCuts->applyTelTypeTest( true ) )
            {
                fillDL2EventDataTree( d, 4, -1. );
                continue;
            }
        }
        fillEcutSub( eMC, E_EcutTelType );
        
        
        //////////////////////////////////////
        // apply direction cut
        //
        // bDirectionCut = false: if direction is inside
        // theta_min and theta_max
        //
        // point source cut; use MC shower direction as reference direction
        bool bDirectionCut = false;
        if( !fIsotropicArrivalDirections )
        {
            if( !iAnaCuts->applyDirectionCuts( true ) )
            {
                bDirectionCut = true;
            }
        }
        // background cut; use (0,0) as reference direction
        // (command line option -d)
        else
        {
            if( !iAnaCuts->applyDirectionCuts( true, 0., 0. ) )
            {
                bDirectionCut = true;
            }
        }
        if( !bDirectionCut )
        {
            fillEcutSub( eMC, E_EcutDirection );
        }
        
        //////////////////////////////////////
        // apply energy reconstruction quality cut
        if( !fIgnoreEnergyReconstruction )
        {
            if( bDebugCuts )
            {
                cout << "#4 EnergyReconstructionQualityCuts ";
                cout << iAnaCuts->applyEnergyReconstructionQualityCuts( iMethod ) << endl;
            }
            if( !iAnaCuts->applyEnergyReconstructionQualityCuts( iMethod, true ) )
            {
                fillDL2EventDataTree( d, 6, -1. );
                continue;
            }
        }
        if( !bDirectionCut )
        {
            fillEcutSub( eMC, E_EcutEnergyReconstruction );
        }
        
        // skip event if no energy has been reconstructed
        // get energy according to reconstruction method
        if( fIgnoreEnergyReconstruction )
        {
            eRec = log10( d->MCe0 );
            eRecLin = d->MCe0;
        }
        else if( d->getEnergy_TeV() > 0. )
        {
            eRec = d->getEnergy_Log10();
            eRecLin = d->getEnergy_TeV();
        }
        else
        {
            fillDL2EventDataTree( d, 6, -1. );
            continue;
        }
        
        /////////////////////////////////////////////////////////
        // fill response matrix after quality and direction cuts
        
        if( !bDirectionCut )
        {
            // loop over all az bins
            for( unsigned int i_az = 0; i_az < fVMinAz.size(); i_az++ )
            {
                if( !testAzimuthInterval( d, fZe[ize], fVMinAz[i_az], fVMaxAz[i_az] ) )
                {
                    continue;
                }
                // loop over all spectral index
                for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
                {
                    fillHistogram( E_2D, E_ResponseMatrixQC, s, i_az, eRec, eMC );
                    fillHistogram( E_2D, E_ResponseMatrixFineQC, s, i_az, eRec, eMC, i_weight );
                }
            }
        }
        
        //////////////////////////////////////
        // apply gamma hadron cuts
        if( bDebugCuts )
        {
            cout << "#3 CUT ISGAMMA " << iAnaCuts->isGamma( i ) << endl;
        }
        if( !iAnaCuts->isGamma( i, true ) )
        {
            if( ( fIsotropicArrivalDirections && !bDirectionCut ) || !fIsotropicArrivalDirections )
            {
                fillDL2EventDataTree( d, 7, iAnaCuts->getTMVA_EvaluationResult() );
            }
            continue;
        }
        if( !bDirectionCut )
        {
            fillEcutSub( eMC, E_EcutGammaHadron );
            fillDL2EventDataTree( d, 5, iAnaCuts->getTMVA_EvaluationResult() );
        }
        // remaining events
        else
        {
            if( !fIsotropicArrivalDirections )
            {
                fillDL2EventDataTree( d, 0, iAnaCuts->getTMVA_EvaluationResult() );
            }
        }
        
        // unique event counter
        // (make sure that map doesn't get too big)
        if( !bDirectionCut && iSuccessfullEventStatistics >= 0 )
        {
            iSuccessfullEventStatistics++;
        }
        
        // loop over all az bins
        for( unsigned int i_az = 0; i_az < fVMinAz.size(); i_az++ )
        {
            if( !testAzimuthInterval( d, fZe[ize], fVMinAz[i_az], fVMaxAz[i_az] ) )
            {
                fillDL2EventDataTree( d, 11, -1 );
                continue;
            }
            
            // fill tree with acceptance information after cuts (needed to construct background model in ctools)
            // NOTE: This tree is currently allways filled with the eventdisplay reconstruction results.
            if( !bDirectionCut && fRunPara->fgetXoff_Yoff_afterCut )
            {
                fXoff_aC = d->Xoff;
                fYoff_aC = d->Yoff;
                fXoff_derot_aC = d->Xoff_derot;
                fYoff_derot_aC = d->Yoff_derot;
                fErec = eRecLin;
                fEMC  = d->MCe0;
                fCRweight = getCRWeight( d->MCe0, hV_HIS1D[E_Emc][0][az_bin_index], true ); //So that the acceptance can be normalised to the CR spectrum.
                // when running on gamma, this should return 1.
                fAcceptance_AfterCuts_tree->Fill();
            }
            
            
            // loop over all spectral index
            for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
            {
                // weight by spectral index
                if( fSpectralWeight )
                {
                    fSpectralWeight->setSpectralIndex( fVSpectralIndex[s] );
                    i_weight = fSpectralWeight->getSpectralWeight( d->MCe0 );
                }
                else
                {
                    i_weight = 0.;
                }
                
                ////////////////////////////////////////////
                // fill effective areas before direction cut
                fillHistogram( E_1D, E_EcutNoTh2, s, i_az, eMC, i_weight );
                fillHistogram( E_1D, E_EcutRecNoTh2, s, i_az, eRec, i_weight );
                // fill response matrix (migration matrix) before
                // direction cut
                fillHistogram( E_2D, E_ResponseMatrixNoDirectionCut, s, i_az, eRec, eMC, i_weight );
                fillHistogram( E_2D, E_ResponseMatrixFineNoDirectionCut, s, i_az, eRec, eMC, i_weight );
                fillHistogram( E_2D, E_EsysMCRelative2DNoDirectionCut, s, i_az, eMC, eRecLin / d->MCe0 );
                
                /////////////////////////
                // apply direction cut
                if( bDirectionCut )
                {
                    continue;
                }
                
                ///////////////////////////////////////////////////////////
                // from here on: after gamma/hadron and after direction cut
                
                // fill true MC energy (hVEmc is in true MC energies)
                fillHistogram( E_1D, E_Ecut, s, i_az, eMC, i_weight );
                fillHistogram( E_1D, E_EcutUW, s, i_az, eMC, 1. );
                fillHistogram( E_1D, E_Ecut500, s, i_az, eMC, i_weight );
                fillHistogram( E_1D, E_EcutRec, s, i_az, eRec, i_weight );
                fillHistogram( E_1D, E_EcutRecUW, s, i_az, eRec, 1. );
                fillHistogram( E_1P, E_EsysMCRelative, s, i_az, eMC, ( eRecLin - d->MCe0 ) / d->MCe0 );
                fillHistogram( E_2D, E_EsysMCRelativeRMS, s, i_az, eMC, ( eRecLin - d->MCe0 ) / d->MCe0 );
                
                fillHistogram( E_2D, E_EsysMCRelative2D, s, i_az, eMC, eRecLin / d->MCe0 );
                fillHistogram( E_2D, E_Esys2D, s, i_az, eMC,  eRec - eMC );
                fillHistogram( E_2D, E_ResponseMatrix, s, i_az, eRec, eMC );
                fillHistogram( E_2D, E_ResponseMatrixFine, s, i_az, eRec, eMC, i_weight );
                // events weighted by CR spectra
                fillHistogram( E_1D, E_WeightedRate, s, i_az, eRec,
                               getCRWeight( d->MCe0, hV_HIS1D[E_Emc][s][i_az],
                                            false, hV_HIS1D[E_WeightedRate][s][i_az] ) );
                fillHistogram( E_1D, E_WeightedRate005, s, i_az, eRec,
                               getCRWeight( d->MCe0, hV_HIS1D[E_Emc][s][i_az],
                                            false, hV_HIS1D[E_WeightedRate005][s][i_az] ) );
            }
        }
        // don't do anything between here and the end of the loop! Never!
    }  // end of loop
    /////////////////////////////////////////////////////////////////////////////
    
    
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // calculate effective areas and fill output trees
    //
    /////////////////////////////////////////////////////////////////////////////
    
    ze = fZe[ize];
    fTNoise = fNoise[ize];
    // WARNING: hardwired values - not used to my knowledge anywhere? (GM)
    fTNoisePE = ( double )( fNoise[ize] ) / 0.15 * 1.e9;
    fTPedvar = fPedVar[ize];
    fXoff = fXWobble[ize];
    fYoff = fYWobble[ize];
    fWoff = sqrt( fXoff * fXoff + fYoff * fYoff );
    
    cout << "calculating effective areas..." << endl;
    
    // loop over all spectral index
    for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
    {
        fSpectralIndex = fVSpectralIndex[s];
        // loop over all az bins
        for( unsigned int i_az = 0; i_az < fVMinAz.size(); i_az++ )
        {
            fAzBin = ( int )i_az;
            fMinAz = fVMinAz[i_az];
            fMaxAz = fVMaxAz[i_az];
            
            resetEffAreaArray( e0 );
            for( int b = 0; b < hV_HIS1D[E_Emc][s][i_az]->GetNbinsX(); b++ )
            {
                e0[b] = hV_HIS1D[E_Emc][s][i_az]->GetBinCenter( b + 1 );
            }
            
            // bayesdivide works only for weights == 1
            // errors might be wrong, since histograms are filled with weights != 1
            if( !binomialDivide( gEffAreaMC, hV_HIS1D[E_Ecut][s][i_az], hV_HIS1D[E_Emc][s][i_az],
                                 eff, eff_error ) )
            {
                cout << "VEffectiveAreaCalculator::fill: error calculating effective area vs MC energy" << endl;
                cout << "s : " << s << " , az: " << i_az << endl;
            }
            if( !binomialDivide( gEffAreaRec, hV_HIS1D[E_EcutRec][s][i_az], hV_HIS1D[E_Emc][s][i_az],
                                 Rec_eff, Rec_eff_error ) )
            {
                cout << "VEffectiveAreaCalculator::fill: error calculating effective area vs rec energy" << endl;
                cout << "s : " << s << " , az: " << i_az << endl;
            }
            if( !binomialDivide( gEffAreaNoTh2MC, hV_HIS1D[E_EcutNoTh2][s][i_az], hV_HIS1D[E_Emc][s][i_az],
                                 effNoTh2, effNoTh2_error ) )
            {
                cout << "VEffectiveAreaCalculator::fill: error calculating effective area vs MC energy" << endl;
                cout << "s : " << s << " , az: " << i_az << endl;
            }
            if( !binomialDivide( gEffAreaNoTh2Rec, hV_HIS1D[E_EcutRecNoTh2][s][i_az], hV_HIS1D[E_Emc][s][i_az],
                                 Rec_effNoTh2, Rec_effNoTh2_error ) )
            {
                cout << "VEffectiveAreaCalculator::fill: error calculating effective area vs rec energy" << endl;
                cout << "s : " << s << " , az: " << i_az << endl;
            }
            // normalize all response matrices
            VHistogramUtilities::normalizeTH2D_y( hV_HIS2D[E_ResponseMatrix][s][i_az] );
            VHistogramUtilities::normalizeTH2D_y( hV_HIS2D[E_ResponseMatrixQC][s][i_az] );
            VHistogramUtilities::normalizeTH2D_y( hV_HIS2D[E_ResponseMatrixNoDirectionCut][s][i_az] );
            
            resetEffAreaArray( esys_rel );
            resetEffAreaArray( Rec_angRes_p68 );
            resetEffAreaArray( Rec_angRes_p80 );
            resetEffAreaArray( Rec_angRes_kingSigma );
            resetEffAreaArray( Rec_angRes_kingGamma );
            
            if( hV_HIS1P.find( E_EsysMCRelative ) != hV_HIS1P.end()
                    && hV_HIS1P[E_EsysMCRelative][s][i_az] )
            {
                for( int b = 0; b < hV_HIS1P[E_EsysMCRelative][s][i_az]->GetNbinsX(); b++ )
                {
                    esys_rel[b] = hV_HIS1P[E_EsysMCRelative][s][i_az]->GetBinContent( b + 1 );
                }
            }
            
            // copy all histograms
            resetHistograms( ize );
            
            map< int, TH1D* >::iterator h_HIS1D_iterator;
            for( h_HIS1D_iterator = h_HIS1D.begin();
                    h_HIS1D_iterator !=  h_HIS1D.end();
                    ++h_HIS1D_iterator )
            {
                copyHistograms( h_HIS1D_iterator->second,
                                hV_HIS1D[h_HIS1D_iterator->first][s][i_az],
                                false );
            }
            map< int, TProfile* >::iterator h_HIS1P_iterator;
            for( h_HIS1P_iterator = h_HIS1P.begin();
                    h_HIS1P_iterator !=  h_HIS1P.end();
                    ++h_HIS1P_iterator )
            {
                copyProfileHistograms( h_HIS1P_iterator->second,
                                       hV_HIS1P[h_HIS1P_iterator->first][s][i_az] );
            }
            map< int, TH2D* >::iterator h_HIS2D_iterator;
            for( h_HIS2D_iterator = h_HIS2D.begin();
                    h_HIS2D_iterator !=  h_HIS2D.end();
                    ++h_HIS2D_iterator )
            {
                copyHistograms( h_HIS2D_iterator->second,
                                hV_HIS2D[h_HIS2D_iterator->first][s][i_az],
                                true );
            }
            
            copyHistograms( hAngularLogDiffEmc_2D, hVAngularLogDiffEmc_2D[i_az], true );
            
            // fill angular resolution vs energy
            fillAngularResolution( i_az, false );
            fillAngularResolution( i_az, true );
            
            // Memory Issue when trying to store all the response matrices in the
            // final effectivea area file.
            // Approximating the response matrix to have a gaussian profile.
            // This approximation breaks down as we approach the energy threshold
            // (mainly due to low statistics) consistent and accurate above energy thresholds
            //
            
            /*nbins_MC_Res = hVResponseMatrixFineQC[s][i_az]->GetYaxis()->GetNbins();
            for ( int i_ybin = 0 ; i_ybin < nbins_MC_Res ; i_ybin++)
            {
                fGauss->SetParameter(0,0);
                fGauss->SetParameter(1,0);
                fGauss->SetParameter(2,0);
            
                // Getting a slice
                TH1D *i_slice = hVResponseMatrixFineQC[s][i_az]->ProjectionX("i_slice_Project", i_ybin,i_ybin);
                // Fitting quietly
                if( i_slice->GetEntries() > 0 )
                {
                    i_slice->Fit("fGauss","0q");
            
                    e_MC_Res[i_ybin] = hVResponseMatrixFineQC[s][i_az]->GetYaxis()->GetBinCenter(i_ybin);
                    e_Rec_Res[i_ybin] = fGauss->GetParameter(1);
                    e_Rec_Res_Err[i_ybin] = fGauss->GetParameter(2);
                }
                delete i_slice;
            } */
            
            fEffArea->Fill();
        }
    }
    
    for( unsigned int i = 0; i < fCuts.size(); i++ )
    {
        fCuts[i]->printCutStatistics();
    }
    
    // print out uniqueness of events
    /*    cout << "event statistics: " << endl;
        if( iSuccessfullEventStatistics > 0 )
        {
           map< unsigned int, unsigned short int>::iterator it;
           for( it = fUniqueEventCounter.begin(); it != fUniqueEventCounter.end(); it++ )
           {
    	  if( (*it).second > 1 )
    	  {
    	     cout << "\t multiple event after cuts at " << (double)((*it).first)/1.e5 << " TeV, used " << (*it).second << " times" << endl;
    	  }
           }
        }
        else iSuccessfullEventStatistics *= -1; */
    if( iSuccessfullEventStatistics < 0 )
    {
        iSuccessfullEventStatistics *= -1;
    }
    cout << "\t total number of events after cuts: " << iSuccessfullEventStatistics << endl;
    
    return true;
}


/*!
 *
 *  CALLED TO USE EFFECTIVE AREAS
 *
 * interpolate between effective area from different zenith angles with cos weights
 *
 *   return value is 1/effective area
 *
 *
 */
double VEffectiveAreaCalculator::getEffectiveArea( double erec, double ze, double woff, double iPedVar, double iSpectralIndex,
        bool bAddtoMeanEffectiveArea )
{
    if( bNOFILE )
    {
        return 1.;
    }
    
    // read effective areas from histograms
    return getEffectiveAreasFromHistograms( erec, ze, woff, iPedVar, iSpectralIndex, bAddtoMeanEffectiveArea );
}



/*
      this function always returns a vector of size 2
*/
vector< unsigned int > VEffectiveAreaCalculator::getUpperLowBins( vector< double > i_values, double d )
{
    vector< unsigned int > i_temp( 2, 0 );
    
    double i_min = 1.e10;
    unsigned int i_min_index = 0;
    double i_max = -1.e10;
    unsigned int i_max_index = 0;
    
    for( unsigned int i = 0; i < i_values.size(); i++ )
    {
        if( i_values[i] < i_min )
        {
            i_min = i_values[i];
            i_min_index = i;
        }
        if( i_values[i] > i_max )
        {
            i_max = i_values[i];
            i_max_index = i;
        }
    }
    if( d < i_min )
    {
        i_temp[0] = i_min_index;
        i_temp[1] = i_min_index;
        return i_temp;
    }
    if( d > i_max )
    {
        i_temp[0] = i_max_index;
        i_temp[1] = i_max_index;
        return i_temp;
    }
    
    // look for closest values
    double i_diff_upper = 1.e10;
    unsigned int i_diff_upper_index = 0;
    double i_diff_lower = 1.e10;
    unsigned int i_diff_lower_index = 0;
    for( unsigned int i = 0; i < i_values.size(); i++ )
    {
        if( i_values[i] - d > 0. && i_values[i] - d < i_diff_upper )
        {
            i_diff_upper = i_values[i] - d;
            i_diff_upper_index = i;
        }
        if( d - i_values[i] > 0. && d - i_values[i] < i_diff_lower )
        {
            i_diff_lower = d - i_values[i];
            i_diff_lower_index = i;
        }
    }
    i_temp[0] = i_diff_lower_index;
    i_temp[1] = i_diff_upper_index;
    
    return i_temp;
}


/*!
 *  CALLED TO USE EFFECTIVE AREAS
 *
 *  return effective area value for given ze, woff, iPedVar, ...
 *
 *
 */
double VEffectiveAreaCalculator::getEffectiveAreasFromHistograms(
    double erec, double ze, double woff,
    double iPedVar, double iSpectralIndex,
    bool bAddtoMeanEffectiveArea )
{
    vector< float > i_eff_temp( fEff_E0.size(), 0. );
    vector< float > i_eff_MC_temp;
    
    // All the resoponse matrix stuff needs to be defined
    // here even if its note being used.
    vector< float > i_e_MC_Res_temp;
    vector< float > i_e_Rec_Res_temp;
    vector< float > i_e_Rec_Res_Err_temp;
    
    vector< vector < float > > i_ze_eff_MC_temp;
    
    // Response Matrix
    vector< vector < float > > i_ze_e_MC_Res_temp;
    vector< vector < float > > i_ze_e_Rec_Res_temp;
    vector< vector < float > > i_ze_e_Rec_Res_Err_temp;
    // log10 of energy
    if( erec <= 0. )
    {
        return 0.;
    }
    float lerec = log10( erec );
    
    // calculate mean values of
    // phase space parameters
    fEffectiveAreas_meanZe     += ze;
    fEffectiveAreas_meanWoff   += woff;
    fEffectiveAreas_meanPedVar += iPedVar;
    fEffectiveAreas_meanIndex   = iSpectralIndex;
    fEffectiveAreas_meanN++;
    
    ////////////////////////////////////////////////////////
    // get upper and lower zenith angle bins
    ////////////////////////////////////////////////////////
    vector< unsigned int > i_ze_bins = getUpperLowBins( fZe, ze );
    vector< vector< float > > i_ze_eff_temp( 2, i_eff_temp );
    
    // Assigning and shaping vectors
    if( bLikelihoodAnalysis && bIsOn )
    {
        // Temp Response Matrix Vectors
        // Typically nbins_REC != nbins_MC hence requiring
        // the vectors being assigned this way
        i_eff_MC_temp.assign( nbins_MC, 0. );
        
        // Temp EffectiveAreas Vector
        i_ze_eff_MC_temp.resize( 2 );
        i_ze_eff_MC_temp[0].resize( i_eff_MC_temp.size() );
        i_ze_eff_MC_temp[1].resize( i_eff_MC_temp.size() );
        
        i_e_MC_Res_temp.assign( nbins_MC_Res, 0. );
        i_e_Rec_Res_temp.assign( nbins_MC_Res, 0. );
        i_e_Rec_Res_Err_temp.assign( nbins_MC_Res, 0. );
        
        i_ze_e_MC_Res_temp.resize( 2 );
        i_ze_e_Rec_Res_temp.resize( 2 );
        i_ze_e_Rec_Res_Err_temp.resize( 2 );
        
        i_ze_e_MC_Res_temp[0].resize( i_e_MC_Res_temp.size() );
        i_ze_e_Rec_Res_temp[0].resize( i_e_Rec_Res_temp.size() );
        i_ze_e_Rec_Res_Err_temp[0].resize( i_e_Rec_Res_Err_temp.size() );
        
        i_ze_e_MC_Res_temp[1].resize( i_e_MC_Res_temp.size() );
        i_ze_e_Rec_Res_temp[1].resize( i_e_Rec_Res_temp.size() );
        i_ze_e_Rec_Res_Err_temp[1].resize( i_e_Rec_Res_Err_temp.size() );
        
    }
    // following: step-by-step interpolation for
    // ze, az, NSB in effective areas
    // loop over all zenith bins
    for( unsigned int i = 0; i < i_ze_bins.size(); i++ )
    {
        ////////////////////////////////////////////////////////
        // get lower and upper wobble offset bins
        ////////////////////////////////////////////////////////
        if( i_ze_bins[i] < fEff_WobbleOffsets.size() )
        {
            vector< unsigned int > i_woff_bins = getUpperLowBins( fEff_WobbleOffsets[i_ze_bins[i]], woff );
            vector< vector< float > > i_woff_eff_temp( 2, i_eff_temp );
            vector< vector< float > > i_woff_eff_MC_temp;
            vector< vector< float > > i_woff_e_MC_Res_temp;
            vector< vector< float > > i_woff_e_Rec_Res_temp;
            vector< vector< float > > i_woff_e_Rec_Res_Err_temp;
            
            if( bLikelihoodAnalysis && bIsOn )
            {
                // Temp EffectiveAreas Vector
                i_woff_e_MC_Res_temp.resize( 2 );
                i_woff_e_Rec_Res_temp.resize( 2 );
                i_woff_e_Rec_Res_Err_temp.resize( 2 );
                i_woff_eff_MC_temp.resize( 2 );
                
                i_woff_e_MC_Res_temp[0].resize( i_e_MC_Res_temp.size() );
                i_woff_e_Rec_Res_temp[0].resize( i_e_Rec_Res_temp.size() );
                i_woff_e_Rec_Res_Err_temp[0].resize( i_e_Rec_Res_Err_temp.size() );
                
                i_woff_eff_MC_temp[0].resize( i_eff_MC_temp.size() );
                i_woff_eff_MC_temp[1].resize( i_eff_MC_temp.size() );
                
                i_woff_e_MC_Res_temp[1].resize( i_e_MC_Res_temp.size() );
                i_woff_e_Rec_Res_temp[1].resize( i_e_Rec_Res_temp.size() );
                i_woff_e_Rec_Res_Err_temp[1].resize( i_e_Rec_Res_Err_temp.size() );
            }
            // wobble offset bins
            for( unsigned int w = 0; w < i_woff_bins.size(); w++ )
            {
                ////////////////////////////////////////////////////////
                // get lower and upper noise bins
                ////////////////////////////////////////////////////////
                if( i_ze_bins[i] < fEff_Noise.size() && i_woff_bins[w] < fEff_Noise[i_ze_bins[i]].size() )
                {
                    vector< unsigned int > i_noise_bins = getUpperLowBins( fEff_Noise[i_ze_bins[i]][i_woff_bins[w]], iPedVar );
                    vector< vector< float > > i_noise_eff_temp( 2, i_eff_temp );
                    vector< vector< float > > i_noise_eff_MC_temp;
                    vector< vector< float > > i_noise_e_MC_Res_temp;
                    vector< vector< float > > i_noise_e_Rec_Res_temp;
                    vector< vector< float > > i_noise_e_Rec_Res_Err_temp;
                    
                    if( bLikelihoodAnalysis && bIsOn )
                    {
                        // Temp EffectiveAreas Vector
                        i_noise_eff_MC_temp.resize( 2 );
                        i_noise_e_MC_Res_temp.resize( 2 );
                        i_noise_e_Rec_Res_temp.resize( 2 );
                        i_noise_e_Rec_Res_Err_temp.resize( 2 );
                        
                        i_noise_e_MC_Res_temp[0].resize( i_e_MC_Res_temp.size() );
                        i_noise_e_Rec_Res_temp[0].resize( i_e_Rec_Res_temp.size() );
                        i_noise_e_Rec_Res_Err_temp[0].resize( i_e_Rec_Res_Err_temp.size() );
                        
                        i_noise_eff_MC_temp[0].resize( i_eff_MC_temp.size() );
                        i_noise_eff_MC_temp[1].resize( i_eff_MC_temp.size() );
                        
                        i_noise_e_MC_Res_temp[1].resize( i_e_MC_Res_temp.size() );
                        i_noise_e_Rec_Res_temp[1].resize( i_e_Rec_Res_temp.size() );
                        i_noise_e_Rec_Res_Err_temp[1].resize( i_e_Rec_Res_Err_temp.size() );
                    }
                    // pedvar / NSB interpolation
                    for( unsigned int n = 0; n < i_noise_bins.size(); n++ )
                    {
                        ////////////////////////////////////////////////////////
                        // get lower and upper spectral index bins
                        ////////////////////////////////////////////////////////
                        if( i_ze_bins[i] < fEff_SpectralIndex.size()
                                && i_woff_bins[w] < fEff_SpectralIndex[i_ze_bins[i]].size()
                                && i_noise_bins[n] < fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]].size() )
                        {
                            vector< unsigned int > i_index_bins = getUpperLowBins( fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]],
                                                                  iSpectralIndex );
                            unsigned int i_ID_0 = i_index_bins[0] + 100 * ( i_noise_bins[n] + 100 * ( i_woff_bins[w] + 100 * i_ze_bins[i] ) );
                            unsigned int i_ID_1 = i_index_bins[1] + 100 * ( i_noise_bins[n] + 100 * ( i_woff_bins[w] + 100 * i_ze_bins[i] ) );
                            i_noise_eff_temp[n] = interpolate_effectiveArea( iSpectralIndex,
                                                  fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                                  fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                                  fEffArea_map[i_ID_0],
                                                  fEffArea_map[i_ID_1], false );
                                                  
                            if( bLikelihoodAnalysis && bIsOn )
                            {
                            
                                // cout << "Map Entry: " << fEntry_map[i_ID_0] << " " <<  fEntry_map[i_ID_1] << endl;
                                // Interpolated over index
                                i_noise_eff_MC_temp[n] = interpolate_effectiveArea( iSpectralIndex,
                                                         fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                                         fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                                         fEffAreaMC_map[i_ID_0],
                                                         fEffAreaMC_map[i_ID_1], false );
                                                         
                                i_noise_e_MC_Res_temp[n] = interpolate_effectiveArea( iSpectralIndex,
                                                           fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                                           fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                                           fe_MC_Res_map[i_ID_0],
                                                           fe_MC_Res_map[i_ID_1], false );
                                                           
                                i_noise_e_Rec_Res_temp[n] = interpolate_effectiveArea( iSpectralIndex,
                                                            fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                                            fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                                            fe_Rec_Res_map[i_ID_0],
                                                            fe_Rec_Res_map[i_ID_1], false );
                                                            
                                i_noise_e_Rec_Res_Err_temp[n] = interpolate_effectiveArea( iSpectralIndex,
                                                                fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                                                fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                                                fe_Rec_Res_Err_map[i_ID_0],
                                                                fe_Rec_Res_Err_map[i_ID_1], false );
                            }
                        }
                        else
                        {
                            cout << "VEffectiveAreaCalculator::getEffectiveAreasFromHistograms error: spectral index index out of range: ";
                            cout << i_ze_bins[i] << " " << fEff_SpectralIndex.size();
                            if( i_ze_bins[i] < fEff_SpectralIndex.size() )
                            {
                                cout << " " << i_woff_bins[w] << " " << fEff_SpectralIndex[i_ze_bins[i]].size();
                            }
                            if( i_woff_bins[w] < fEff_SpectralIndex[i_ze_bins[i]].size() )
                            {
                                cout << " " << i_noise_bins[n] <<  fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]].size();
                            }
                            cout << endl;
                            return -1.;
                        }
                        ////////////////////////////////////////////////////////
                    }
                    i_woff_eff_temp[w] = interpolate_effectiveArea( iPedVar,
                                         fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                         fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                         i_noise_eff_temp[0],
                                         i_noise_eff_temp[1], false );
                    if( bLikelihoodAnalysis && bIsOn )
                    {
                    
                        // Interpolating noise
                        i_woff_eff_MC_temp[w] = interpolate_effectiveArea( iPedVar,
                                                fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                                fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                                i_noise_eff_MC_temp[0],
                                                i_noise_eff_MC_temp[1], false );
                                                
                        i_woff_e_MC_Res_temp[w] = interpolate_effectiveArea( iPedVar,
                                                  fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                                  fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                                  i_noise_e_MC_Res_temp[0],
                                                  i_noise_e_MC_Res_temp[1], false );
                                                  
                        i_woff_e_Rec_Res_temp[w] = interpolate_effectiveArea( iPedVar,
                                                   fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                                   fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                                   i_noise_e_Rec_Res_temp[0],
                                                   i_noise_e_Rec_Res_temp[1], false );
                                                   
                        i_woff_e_Rec_Res_Err_temp[w] = interpolate_effectiveArea( iPedVar,
                                                       fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                                       fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                                       i_noise_e_Rec_Res_Err_temp[0],
                                                       i_noise_e_Rec_Res_Err_temp[1], false );
                    }
                }
                else
                {
                    cout << "VEffectiveAreaCalculator::getEffectiveAreasFromHistograms error: noise index out of range: ";
                    cout << i_ze_bins[i] << " " << fEff_Noise.size();
                    if( i_ze_bins[i] < fEff_Noise.size() )
                    {
                        cout << " " << i_woff_bins[w] << " " << fEff_Noise[i_ze_bins[i]].size() << endl;
                    }
                    cout << endl;
                    return -1.;
                }
            }
            i_ze_eff_temp[i] = interpolate_effectiveArea( woff,
                               fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                               fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                               i_woff_eff_temp[0],
                               i_woff_eff_temp[1], false );
            if( bLikelihoodAnalysis && bIsOn )
            {
                // Interpolating Wobble
                i_ze_eff_MC_temp[i] = interpolate_effectiveArea( woff,
                                      fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                                      fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                                      i_woff_eff_MC_temp[0],
                                      i_woff_eff_MC_temp[1], false );
                                      
                i_ze_e_MC_Res_temp[i] = interpolate_effectiveArea( woff,
                                        fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                                        fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                                        i_woff_e_MC_Res_temp[0],
                                        i_woff_e_MC_Res_temp[1], false );
                                        
                i_ze_e_Rec_Res_temp[i] = interpolate_effectiveArea( woff,
                                         fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                                         fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                                         i_woff_e_Rec_Res_temp[0],
                                         i_woff_e_Rec_Res_temp[1], false );
                                         
                i_ze_e_Rec_Res_Err_temp[i] = interpolate_effectiveArea( woff,
                                             fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                                             fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                                             i_woff_e_Rec_Res_Err_temp[0],
                                             i_woff_e_Rec_Res_Err_temp[1], false );
            }
            
        }
        else
        {
            cout << "VEffectiveAreaCalculator::getEffectiveAreasFromHistograms error: woff index out of range: ";
            cout << i_ze_bins[i] << " " << fEff_WobbleOffsets.size() << endl;
            return -1.;
        }
    }
    i_eff_temp = interpolate_effectiveArea( ze, fZe[i_ze_bins[0]], fZe[i_ze_bins[1]], i_ze_eff_temp[0], i_ze_eff_temp[1], true );
    
    if( bLikelihoodAnalysis && bIsOn )
    {
    
        // Interpolating Zenith
        i_eff_MC_temp = interpolate_effectiveArea( ze, fZe[i_ze_bins[0]], fZe[i_ze_bins[1]],
                        i_ze_eff_MC_temp[0], i_ze_eff_MC_temp[1], true );
        i_e_MC_Res_temp = interpolate_effectiveArea( ze, fZe[i_ze_bins[0]], fZe[i_ze_bins[1]],
                          i_ze_e_MC_Res_temp[0], i_ze_e_MC_Res_temp[1], true );
        i_e_Rec_Res_temp = interpolate_effectiveArea( ze, fZe[i_ze_bins[0]], fZe[i_ze_bins[1]],
                           i_ze_e_Rec_Res_temp[0], i_ze_e_Rec_Res_temp[1], true );
        i_e_Rec_Res_Err_temp = interpolate_effectiveArea( ze, fZe[i_ze_bins[0]], fZe[i_ze_bins[1]],
                               i_ze_e_Rec_Res_Err_temp[0], i_ze_e_Rec_Res_Err_temp[1], true );
                               
    }
    if( fEff_E0.size() == 0 )
    {
        return -1.;
    }
    
    // mean effective area calculation
    if( bAddtoMeanEffectiveArea && fVTimeBinnedMeanEffectiveArea.size() == i_eff_temp.size() )
    {
        for( unsigned int i = 0; i < i_eff_temp.size(); i++ )
        {
            fVTimeBinnedMeanEffectiveArea[i] += i_eff_temp[i];
        }
        fNTimeBinnedMeanEffectiveArea++;
    }
    
    if( bLikelihoodAnalysis && bIsOn )
    {
        // Adding to mean effective area (MC)
        if( bAddtoMeanEffectiveArea && fVTimeBinnedMeanEffectiveAreaMC.size() == i_eff_MC_temp.size() )
        {
            for( unsigned int i = 0; i < fEff_E0.size() ; i++ )
            {
                if( i_eff_MC_temp[i] > 1.e-9 )
                {
                    fVTimeBinnedMeanEffectiveAreaMC[i] += i_eff_MC_temp[i];
                }
            }
            fNTimeBinnedMeanEffectiveAreaMC++;
        }
        // adding to mean response matrix
        addMeanResponseMatrix( i_e_MC_Res_temp, i_e_Rec_Res_temp, i_e_Rec_Res_Err_temp );
    }
    /////////////////////////////////////////
    // effective area for a specific energy
    // (as requested by VStereoAnalysis)
    /////////////////////////////////////////
    float i_eff_e = 1.;
    unsigned int ie0_low = 0;
    unsigned int ie0_up = 0;
    
    if( lerec <= fEff_E0[0] )
    {
        ie0_low = ie0_up = 0;
    }
    else if( lerec > fEff_E0[fEff_E0.size() - 1] )
    {
        ie0_low = ie0_up = fEff_E0.size() - 1;
    }
    else
    {
        for( unsigned int j = 1; j < fEff_E0.size() - 1; j++ )
        {
            if( lerec > fEff_E0[j] )
            {
                ie0_low = j;
                ie0_up = j + 1;
            }
        }
    }
    
    ///////////////////////////////////
    // final result on effective area:
    // linear interpolate between two
    // adjacent energy bins
    ///////////////////////////////////
    i_eff_e = VStatistics::interpolate( i_eff_temp[ie0_low], fEff_E0[ie0_low],
                                        i_eff_temp[ie0_up], fEff_E0[ie0_up],
                                        lerec, false );
                                        
    if( i_eff_e > 0. )
    {
        return 1. / i_eff_e;
    }
    
    return -1.;
}

// reset the sum of effective areas

void VEffectiveAreaCalculator::resetTimeBin()
{
    for( unsigned int i = 0; i < fVTimeBinnedMeanEffectiveArea.size(); i++ )
    {
        fVTimeBinnedMeanEffectiveArea[i] = 0;
    }
    fNTimeBinnedMeanEffectiveArea = 0;
    for( unsigned int i = 0; i < fVTimeBinnedMeanEffectiveAreaMC.size(); i++ )
    {
        fVTimeBinnedMeanEffectiveAreaMC[i] = 0;
    }
    fNTimeBinnedMeanEffectiveAreaMC = 0;
}


/*
 *  CALLED FOR CALCULATION OF EFFECTIVE AREAS
 */
void VEffectiveAreaCalculator::resetHistograms( unsigned int ize )
{
    if( ize >= fZe.size() )
    {
        return;
    }
    
    char htitle[400];
    if( hisTreeListofHistograms )
    {
        TIter next( hisTreeListofHistograms );
        while( TH1* obj = ( TH1* )next() )
        {
            obj->Reset();
            string iTemp = obj->GetTitle();
            if( iTemp.find( "(" ) != string::npos )
            {
                iTemp = iTemp.substr( 0, iTemp.find( "(" ) - 1 );
            }
            sprintf( htitle, "%s (%1.f deg)", iTemp.c_str(), fZe[ize] );
            obj->SetTitle( htitle );
        }
    }
    
    // set titles of effective area graphs
    sprintf( htitle, "effective area vs E_{MC} (%.1f deg)", fZe[ize] );
    gEffAreaMC->SetTitle( htitle );
    
    sprintf( htitle, "effective area vs E_{rec} (%.1f deg)", fZe[ize] );
    gEffAreaRec->SetTitle( htitle );
    
    sprintf( htitle, "effective area vs E_{MC} (%.1f deg)", fZe[ize] );
    gEffAreaNoTh2MC->SetTitle( htitle );
    
    sprintf( htitle, "effective area vs E_{rec} (%.1f deg)", fZe[ize] );
    gEffAreaNoTh2Rec->SetTitle( htitle );
}


void VEffectiveAreaCalculator::setWobbleOffset( double x, double y )
{
    if( fXWobble.size() == 0 )
    {
        fXWobble.push_back( x );
    }
    else
    {
        fXWobble[0] = x;
    }
    if( fYWobble.size() == 0 )
    {
        fYWobble.push_back( y );
    }
    else
    {
        fYWobble[0] = y;
    }
    
}

void VEffectiveAreaCalculator::setNoiseLevel( int iN, double iP )
{
    if( fNoise.size() == 0 )
    {
        fNoise.push_back( iN );
    }
    else
    {
        fNoise[0] = iN;
    }
    
    if( fPedVar.size() == 0 )
    {
        fPedVar.push_back( iP );
    }
    else
    {
        fPedVar[0] = iP;
    }
}

/*!
    set azimuth cut for effective area fillings
*/
void VEffectiveAreaCalculator::setAzimuthCut( int iAzBin, double iAzMin, double iAzMax )
{
    fAzBin = iAzBin;
    fMinAz = iAzMin;
    fMaxAz = iAzMax;
}

void VEffectiveAreaCalculator::resetEffAreaArray( float* v )
{
    if( v )
    {
        for( unsigned int i = 0; i < VMAXBINS; i++ )
        {
            v[i] = 0.;
        }
    }
}


bool VEffectiveAreaCalculator::binomialDivide( TGraphAsymmErrors* g, TH1D* hrec, TH1D* hmc,
        float* i_eff, float* i_eff_error )
{
    if( !g )
    {
        cout << "VEffectiveAreaCalculator::binomialDivide error: no return graph given" << endl;
        return false;
    }
    if( !hrec )
    {
        cout << "VEffectiveAreaCalculator::binomialDivide error: no histogram with reconstructed events given" << endl;
        return false;
    }
    if( !hmc )
    {
        cout << "VEffectiveAreaCalculator::binomialDivide error: no histogram with simulated events given" << endl;
        return false;
    }
    resetEffAreaArray( i_eff );
    resetEffAreaArray( i_eff_error );
    
    int z = 0;
    double pj = 0.;
    double pr = 0.;
    double pm = 0.;
    double sj = 0.;
    
    for( int b = 0; b < hmc->GetNbinsX(); b++ )
    {
        if( hmc->GetBinContent( b + 1 ) > 0 && hrec->GetBinContent( b + 1 ) > 0 )
        {
            pj = hrec->GetBinContent( b + 1 ) / hmc->GetBinContent( b + 1 );
            // error calculation for effective areas
            //  this far from being straightforward!
            //  none of the methods works consistently, therefore the simplest (normal) solution
            //  is used.
            //  note: approach is only correct for unweighted histograms
            //
            // error calculation assuming binomial distributions
            // (see Blobel/Lohrmann; chapter 11.2 (Akzeptanzkorrekturen)
            pr = hrec->GetBinError( b + 1 );
            pm = hmc->GetBinError( b + 1 );
            if( pj != 1. )
            {
                sj = TMath::Abs( ( ( 1. - 2.*pj ) * pr * pr + pj * pj * pm * pm )
                                 / ( hmc->GetBinContent( b + 1 ) * hmc->GetBinContent( b + 1 ) ) );
            }
            else
            {
                sj = 0.;
            }
            sj = sqrt( sj ) * fMC_ScatterArea;
            pj *= fMC_ScatterArea;
            // fill effective area graphs
            g->SetPoint( z, hmc->GetBinCenter( b + 1 ), pj );
            g->SetPointError( z, 0., 0., sj, sj );
            if( i_eff )
            {
                i_eff[b] = pj;
            }
            if( i_eff_error )
            {
                i_eff_error[b] = sj;
            }
            z++;
        }
        else
        {
            if( i_eff )
            {
                i_eff[b] = 0.;
            }
            if( i_eff_error )
            {
                i_eff_error[b] = 0.;
            }
        }
    }
    g->Set( z );
    
    return true;
}

/*

    running mean smoothing of effective area curve
    (not default)

    needs careful testing before using
*/
void VEffectiveAreaCalculator::smoothEffectiveAreas( map< unsigned int, vector< float > > m )
{
    cout << "\t smooth effective areas, parameters: " << fSmoothIter << ", " << fSmoothThreshold << endl;
    
    typedef map< unsigned int, vector< float > >::const_iterator CI;
    
    for( CI p = m.begin(); p != m.end(); ++ p )
    {
        vector< float > itemp = p->second;
        
        for( int l = 0; l < fSmoothIter; l++ )
        {
            for( unsigned int k = 1; k < itemp.size() - 1; k++ )
            {
                if( itemp[k - 1] <= 0 || itemp[k] <= 0. || itemp[k + 1] <= 0. )
                {
                    continue;
                }
                
                // upwards outlier
                if( itemp[k] / itemp[k - 1] > ( 1. + fSmoothThreshold ) && itemp[k] / itemp[k + 1] > ( 1. + fSmoothThreshold ) )
                {
                    itemp[k] = 0.5 * ( itemp[k - 1] + itemp[k + 1] );
                }
                
                // downwards outlier
                if( itemp[k] / itemp[k - 1] < ( 1. - fSmoothThreshold ) && itemp[k] / itemp[k + 1] < ( 1. - fSmoothThreshold ) )
                {
                    itemp[k] = 0.5 * ( itemp[k - 1] + itemp[k + 1] );
                }
            }
        }
    }
}


void VEffectiveAreaCalculator::copyProfileHistograms( TProfile* h1,  TProfile* h2 )
{
    if( !h1 || !h2 )
    {
        return;
    }
    
    string iEOption = h1->GetErrorOption();
    
    for( int b = 0; b <= h2->GetNbinsX(); b++ )
    {
        h1->SetBinContent( b, h2->GetBinContent( b ) * h2->GetBinEntries( b ) );
        
        if( TMath::Abs( h2->GetBinError( b ) ) < 1.e-4 )
        {
            h1->SetBinError( b, 0. );
        }
        else
        {
            if( h2->GetBinEntries( b ) > 0. )
            {
                double iE = h2->GetBinError( b );
                if( iEOption != "s" )
                {
                    iE = h2->GetBinError( b ) * sqrt( h2->GetBinEntries( b ) );
                }
                h1->SetBinError( b,  sqrt( h2->GetBinEntries( b ) * ( h2->GetBinContent( b ) *  h2->GetBinContent( b ) + iE * iE ) ) );
            }
            else
            {
                h1->SetBinError( b, 0. );
            }
        }
        h1->SetBinEntries( b, h2->GetBinEntries( b ) );
    }
    h1->SetEntries( h2->GetEntries() );
}


void VEffectiveAreaCalculator::copyHistograms( TH1* h1,  TH1* h2, bool b2D )
{
    if( !h1 || !h2 )
    {
        return;
    }
    
    if( !b2D )
    {
        for( int b = 0; b <= h2->GetNbinsX(); b++ )
        {
            h1->SetBinContent( b, h2->GetBinContent( b ) );
            h1->SetBinError( b, h2->GetBinError( b ) );
        }
    }
    else
    {
        for( int b = 0; b <= h2->GetNbinsX(); b++ )
        {
            for( int b2 = 0; b2 <= h2->GetNbinsY(); b2++ )
            {
                h1->SetBinContent( b, b2, h2->GetBinContent( b, b2 ) );
                h1->SetBinError( b, b2, h2->GetBinError( b, b2 ) );
            }
        }
    }
    h1->SetEntries( h2->GetEntries() );
}


/*

    return mean effective area for given run

    values are filled ....

*/
TGraphAsymmErrors* VEffectiveAreaCalculator::getMeanEffectiveArea()
{
    if( gMeanEffectiveArea && fEffArea_time.size() > 0 )
    {
        gMeanEffectiveArea->Set( 0 );
        double z = 0;
        double x = 0.;
        double y = 0.;
        
        // set energies
        // energy dimension of fEffArea_time is alway of same size
        for( unsigned int i = 0; i < fEff_E0.size(); i++ )
        {
            gMeanEffectiveArea->SetPoint( i, fEff_E0[i], 0. );
        }
        
        for( unsigned int i = 0; i < fEffArea_time.size(); i++ )
        {
            for( unsigned int j = 0; j < fEffArea_time[i].size(); j++ )
            {
                gMeanEffectiveArea->GetPoint( j, x, y );
                gMeanEffectiveArea->SetPoint( j, x, y + fEffArea_time[i][j] );
            }
        }
        
        z = ( double )fEffArea_time.size();
        for( int i = 0; i < gMeanEffectiveArea->GetN(); i++ )
        {
            if( z > 0. )
            {
                gMeanEffectiveArea->GetPoint( i, x, y );
                gMeanEffectiveArea->SetPoint( i, x, y / z );
            }
        }
        return gMeanEffectiveArea;
    }
    
    return 0;
}


/*

    return mean effective area (MC) for given run

    values are filled ....

*/
TGraphAsymmErrors* VEffectiveAreaCalculator::getMeanEffectiveAreaMC()
{
    cout << "\t\tVEffectiveAreaCalculator::getMeanEffectiveAreaMC() Getting Effective Areas" << endl;
    cout << "Size of MC Effective Area Vector: " << fEffAreaMC_time.size() << endl;
    cout << "Size of REC Effective Area Vector: " << fEffArea_time.size() << endl;
    if( gMeanEffectiveAreaMC && fEffAreaMC_time.size() > 0 )
    {
        gMeanEffectiveArea->Set( 0 );
        double z = 0;
        double x = 0.;
        double y = 0.;
        
        // set energies
        // energy dimension of fEffArea_time is alway of same size
        for( unsigned int i = 0; i < fEff_E0.size(); i++ )
        {
            gMeanEffectiveAreaMC->SetPoint( i, fEff_E0[i], 0. );
        }
        
        for( unsigned int i = 0; i < fEffAreaMC_time.size(); i++ )
        {
            for( unsigned int j = 0; j < fEffAreaMC_time[i].size(); j++ )
            {
                gMeanEffectiveAreaMC->GetPoint( j, x, y );
                gMeanEffectiveAreaMC->SetPoint( j, x, y + fEffAreaMC_time[i][j] );
            }
        }
        
        z = ( double )fEffAreaMC_time.size();
        for( int i = 0; i < gMeanEffectiveAreaMC->GetN(); i++ )
        {
            if( z > 0. )
            {
                gMeanEffectiveAreaMC->GetPoint( i, x, y );
                gMeanEffectiveAreaMC->SetPoint( i, x, y / z );
            }
        }
        return gMeanEffectiveAreaMC;
    }
    
    return 0;
}
/*

    return mean effective area for a given time bin


*/
TGraph2DErrors* VEffectiveAreaCalculator::getTimeBinnedMeanEffectiveArea()
{
    if( gTimeBinnedMeanEffectiveArea )
    {
        int z = 0;
        gTimeBinnedMeanEffectiveArea->Set( 0 );
        for( unsigned int i = 0; i < fEffArea_time.size(); i++ )
        {
            for( unsigned int j = 0; j < fEffArea_time[i].size(); j++ )
            {
                if( fEffArea_time[i][j] > 0. )
                {
                    gTimeBinnedMeanEffectiveArea->SetPoint( z , fEff_E0[j], fEffArea_time[i][j], timebins[i] );
                    gTimeBinnedMeanEffectiveArea->SetPointError( z , 0 , 0 , 0 );
                    z++;
                }
            }
        }
        return gTimeBinnedMeanEffectiveArea;
    }
    
    return 0;
}


/*

   set a new time bin and calculate mean effective area

*/
void VEffectiveAreaCalculator::setTimeBinnedMeanEffectiveArea( double i_time )
{
    timebins.push_back( i_time );
    
    vector < double > inter;
    
    for( unsigned int i = 0; i < fVTimeBinnedMeanEffectiveArea.size(); i++ )
    {
        if( fNTimeBinnedMeanEffectiveArea > 0. )
        {
            inter.push_back( fVTimeBinnedMeanEffectiveArea[i] / fNTimeBinnedMeanEffectiveArea );
        }
        else
        {
            inter.push_back( 0. );
        }
    }
    
    fEffArea_time.push_back( inter );
    if( bLikelihoodAnalysis && bIsOn )
    {
        vector < double > interMC;
        for( unsigned int i = 0; i < fVTimeBinnedMeanEffectiveAreaMC.size(); i++ )
        {
            if( fNTimeBinnedMeanEffectiveAreaMC > 0. )
            {
                interMC.push_back( fVTimeBinnedMeanEffectiveAreaMC[i] / fNTimeBinnedMeanEffectiveAreaMC );
            }
            else
            {
                interMC.push_back( 0. );
            }
        }
        
        fEffAreaMC_time.push_back( interMC );
    }
    
    resetTimeBin();
}


void VEffectiveAreaCalculator::resetHistogramsVectors( )
{
    if( hisVList )
    {
        return;
    }
    
    // loop over all temporary histograms and reset them
    TIter next( hisVList );
    while( TH1* obj = ( TH1* )next() )
    {
        obj->Reset();
    }
}

/*
 * calculate mean systematic error (bias) in energy reconstruction
 *
 * interpolate betwen ze, az, noise, woff bins
 */
TGraphErrors* VEffectiveAreaCalculator::getMeanSystematicErrorHistogram()
{
    if( fEffectiveAreas_meanN <= 0. )
    {
        return 0;
    }
    fEffectiveAreas_meanZe /= fEffectiveAreas_meanN;
    fEffectiveAreas_meanWoff /= fEffectiveAreas_meanN;
    fEffectiveAreas_meanPedVar /= fEffectiveAreas_meanN;
    
    // new graph with mean systematic (bias) error in energy reconstruction
    gMeanSystematicErrorGraph = new TGraphErrors( 1 );
    gMeanSystematicErrorGraph->SetName( "gMeanEnergySystematicError" );
    gMeanSystematicErrorGraph->SetTitle( "" );
    gMeanSystematicErrorGraph->SetMarkerStyle( 20 );
    vector< float > hX;
    vector< float > i_eff_temp;
    
    ////////////////////////////////////////////////////////
    // get upper and lower zenith angle bins
    ////////////////////////////////////////////////////////
    vector< unsigned int > i_ze_bins = getUpperLowBins( fZe, fEffectiveAreas_meanZe );
    vector< vector< float > > i_ze_eff_temp( 2, hX );
    for( unsigned int i = 0; i < i_ze_bins.size(); i++ )
    {
        ////////////////////////////////////////////////////////
        // get lower and upper wobble offset bins
        ////////////////////////////////////////////////////////
        if( i_ze_bins[i] >= fEff_WobbleOffsets.size() )
        {
            continue;
        }
        vector< unsigned int > i_woff_bins = getUpperLowBins( fEff_WobbleOffsets[i_ze_bins[i]], fEffectiveAreas_meanWoff );
        vector< vector< float > > i_woff_eff_temp( 2, hX );
        for( unsigned int w = 0; w < i_woff_bins.size(); w++ )
        {
            ////////////////////////////////////////////////////////
            // get lower and upper noise bins
            ////////////////////////////////////////////////////////
            if( i_ze_bins[i] >= fEff_Noise.size() || i_woff_bins[w] >= fEff_Noise[i_ze_bins[i]].size() )
            {
                continue;
            }
            vector< unsigned int > i_noise_bins = getUpperLowBins( fEff_Noise[i_ze_bins[i]][i_woff_bins[w]], fEffectiveAreas_meanPedVar );
            vector< vector< float > > i_noise_eff_temp( 2, hX );
            for( unsigned int n = 0; n < i_noise_bins.size(); n++ )
            {
                ////////////////////////////////////////////////////////
                // get lower and upper spectral index bins
                ////////////////////////////////////////////////////////
                if( i_ze_bins[i] >= fEff_SpectralIndex.size() || i_woff_bins[w] >= fEff_SpectralIndex[i_ze_bins[i]].size()
                        || i_noise_bins[n] >= fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]].size() )
                {
                    continue;
                }
                vector< unsigned int > i_index_bins = getUpperLowBins( fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]],
                                                      fEffectiveAreas_meanIndex );
                unsigned int i_ID_0 = i_index_bins[0] + 100 * ( i_noise_bins[n] + 100 * ( i_woff_bins[w] + 100 * i_ze_bins[i] ) );
                unsigned int i_ID_1 = i_index_bins[1] + 100 * ( i_noise_bins[n] + 100 * ( i_woff_bins[w] + 100 * i_ze_bins[i] ) );
                i_noise_eff_temp[n] = interpolate_effectiveArea(
                                          fEffectiveAreas_meanIndex,
                                          fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[0]],
                                          fEff_SpectralIndex[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[n]][i_index_bins[1]],
                                          fEff_EsysMCRelative[i_ID_0], fEff_EsysMCRelative[i_ID_1], false );
                // don't interpolate errors, assume they are more or less constant
            }
            i_woff_eff_temp[w] = interpolate_effectiveArea( fEffectiveAreas_meanPedVar,
                                 fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[0]],
                                 fEff_Noise[i_ze_bins[i]][i_woff_bins[w]][i_noise_bins[1]],
                                 i_noise_eff_temp[0], i_noise_eff_temp[1], false );
        }
        i_ze_eff_temp[i] = interpolate_effectiveArea( fEffectiveAreas_meanWoff,
                           fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[0]],
                           fEff_WobbleOffsets[i_ze_bins[i]][i_woff_bins[1]],
                           i_woff_eff_temp[0], i_woff_eff_temp[1], false );
    }
    i_eff_temp = interpolate_effectiveArea( fEffectiveAreas_meanZe,
                                            fZe[i_ze_bins[0]],
                                            fZe[i_ze_bins[1]],
                                            i_ze_eff_temp[0], i_ze_eff_temp[1], true );
    if( fEff_E0.size() == i_eff_temp.size() )
    {
        unsigned int z = 0;
        for( unsigned int i = 0; i < i_eff_temp.size(); i++ )
        {
            // don't expect these values to be exactly zero
            if( TMath::Abs( i_eff_temp[i] ) > 1.e-8 )
            {
                gMeanSystematicErrorGraph->SetPoint( z, fEff_E0[i], i_eff_temp[i] );
                gMeanSystematicErrorGraph->SetPointError( z, 0., 0. );
                z++;
            }
        }
    }
    
    return gMeanSystematicErrorGraph;
}


bool VEffectiveAreaCalculator::setMonteCarloEnergyRange( double iMin, double iMax, double iMCIndex )
{
    if( fSpectralWeight )
    {
        fSpectralWeight->setMCParameter( iMCIndex, iMin, iMax );
        return true;
    }
    
    cout << "VEffectiveAreaCalculator::setMonteCarloEnergyRange error: not spectral weighter defined" << endl;
    return false;
}
/*
 * weight to correct the MC spectrum to the CR spectrum
 *
 * note that this requires the MC spectrum to be a power law, as it is integrated by hand
 *
 * weight derived from dN_CR (E) = w(E) * dN_MC (E)
 *
 */
double VEffectiveAreaCalculator::getCRWeight( double iEMC_TeV_lin, TH1* h, bool for_back_map, TH1* hF )
{
    if( !h || !fRunPara )
    {
        return 1.;
    }
    
    if( !fRunPara->fCREnergySpectrum )
    {
        return 1.;
    }
    
    // energy bin upper and lower bounds
    double iE_min = 0.;
    double iE_max = 0.;
    if( hF )
    {
        iE_min = hF->GetXaxis()->GetBinLowEdge( hF->GetXaxis()->FindBin( log10( iEMC_TeV_lin ) ) );
        iE_max = hF->GetXaxis()->GetBinUpEdge( hF->GetXaxis()->FindBin( log10( iEMC_TeV_lin ) ) );
        iE_min = TMath::Power( 10., iE_min );
        iE_max = TMath::Power( 10., iE_max );
    }
    
    // normalization of MC spectrum
    double c_ig = 0.;
    if( fRunPara->fIgnoreFractionOfEvents > 0. )
    {
        c_ig = fRunPara->fIgnoreFractionOfEvents;
    }
    double c_mc = ( 1.0 - c_ig ) * h->GetEntries()
                  * ( -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. )
                  / ( TMath::Power( fRunPara->fMCEnergy_max, -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. )
                      - TMath::Power( fRunPara->fMCEnergy_min, -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. ) );
                      
                      
    // number of MC events in this energy bin
    double n_mc = 1.;
    if( hF )
    {
        n_mc = c_mc / ( -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. )
               * ( TMath::Power( iE_max , -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. )
                   - TMath::Power( iE_min, -1.*TMath::Abs( fRunPara->fMCEnergy_index ) + 1. ) );
    }
    else
    {
        n_mc = c_mc * TMath::Power( iEMC_TeV_lin, -1.*TMath::Abs( fRunPara->fMCEnergy_index ) );
    }
    
    // number of expected CR events /min/sr (for gammas: /min)
    // (convert to cm2 and min)
    double n_cr = fMC_ScatterArea * 1.e4 * 60.;
    if( hF )
    {
        n_cr *= fRunPara->fCREnergySpectrum->Integral( iE_min, iE_max );
    }
    else
    {
        n_cr *= fRunPara->fCREnergySpectrum->Eval( log10( iEMC_TeV_lin ) );
    }
    
    // (DL3) for the acceptance map construction, the weight is in #/s ()
    if( for_back_map )
    {
        // need to normalize this number considering the cone in which the particle are simulated
        // (the offset bin is not consider so the binning can be changed later).
        if( !fsolid_angle_norm_done )
        {
            Calculate_Bck_solid_angle_norm();
        }
        n_cr = fsolid_angle_norm * n_cr / 60.;
        if( n_mc != 0. )
        {
            return n_cr / n_mc;
        }
        else
        {
            return 0.;
        }
    }
    else
    {
        // getMCSolidAngleNormalization() return a ratio of solid angle (only for gamma? not sure this thing returning something else than 1 here, ever...)
        if( getMCSolidAngleNormalization() > 0. )
        {
            n_cr /= getMCSolidAngleNormalization();    // do we want to leave this line here?
        }
        // the solid angle normalization is done in VSensitivityCalculator
        // return #/min/sr
        if( n_mc != 0. )
        {
            return n_cr / n_mc;
        }
        return 0.;
    }
}

// Calculate the ratio between the solid angle of the cone where the MC have been done and the solid angle of the offset bin considered
void VEffectiveAreaCalculator::Calculate_Bck_solid_angle_norm()
{
    if( fRunPara->fViewcone_max > 0. )
    {
        // solid angle in which the particule have been simulated
        double SolidAngle_MCScatterAngle  =  2 * TMath::Pi() * ( 1. - cos( fRunPara->fViewcone_max * TMath::DegToRad() ) );
        
        fsolid_angle_norm = SolidAngle_MCScatterAngle;
        fsolid_angle_norm_done = true;
    }
    
    return;
}

/*

    check if azimuth is in correct bin

*/
bool VEffectiveAreaCalculator::testAzimuthInterval( CData* d, double iZe, double iMinAz, double iMaxAz )
{
    if( !d )
    {
        return false;
    }
    
    // check at what azimuth bin we are
    if( iZe > 3. )
    {
        // confine MC az to -180., 180.
        if( d->MCaz > 180. )
        {
            d->MCaz -= 360.;
        }
        // expect bin like [135,-135]
        if( iMinAz > iMaxAz )
        {
            if( d->MCaz < iMinAz && d->MCaz > iMaxAz )
            {
                return false;
            }
        }
        // expect bin like [-135,-45.]
        else
        {
            if( d->MCaz < iMinAz || d->MCaz > iMaxAz )
            {
                return false;
            }
        }
    }
    return true;
}

/*

   copy angular resolution values to tree variable

*/
void VEffectiveAreaCalculator::fillAngularResolution( unsigned int i_az, bool iContainment_80p )
{
    if( iContainment_80p && i_az < fGraph_AngularResolution68p.size() && fGraph_AngularResolution68p[i_az] )
    {
        // get first and last energy bin
        double i_emin = 1.e5;
        double i_emax = 1.e-5;
        double x = 0.;
        double y = 0.;
        for( int i = 0; i < fGraph_AngularResolution80p[i_az]->GetN(); i++ )
        {
            fGraph_AngularResolution80p[i_az]->GetPoint( i, x, y );
            if( x < i_emin )
            {
                i_emin = x;
            }
            if( x > i_emax )
            {
                i_emax = x;
            }
        }
        for( int i = 0; i < nbins; i++ )
        {
            if( e0[i] > i_emin && e0[i] < i_emax )
            {
                Rec_angRes_p80[i]  = fGraph_AngularResolution80p[i_az]->Eval( e0[i] );
            }
        }
    }
    else if( !iContainment_80p && i_az < fGraph_AngularResolution68p.size() && fGraph_AngularResolution68p[i_az] )
    {
        double i_emin = 0.;
        double i_emax = 0.;
        double x = 0.;
        double y = 0.;
        for( int i = 0; i < fGraph_AngularResolution68p[i_az]->GetN(); i++ )
        {
            fGraph_AngularResolution68p[i_az]->GetPoint( i, x, y );
            if( x < i_emin )
            {
                i_emin = x;
            }
            if( x > i_emax )
            {
                i_emax = x;
            }
        }
        fGraph_AngularResolution68p[i_az]->GetPoint( 0, i_emin, y );
        fGraph_AngularResolution68p[i_az]->GetPoint( fGraph_AngularResolution68p[i_az]->GetN(), i_emax, y );
        for( int i = 0; i < nbins; i++ )
        {
            if( e0[i] > i_emin && e0[i] < i_emax )
            {
                Rec_angRes_p68[i] = fGraph_AngularResolution68p[i_az]->Eval( e0[i] );
            }
        }
    }
    
    // if on 68% containment, also fill king gamma/sigma parameters
    // (so they only get filled once, not twice (once for each containment radius) )
    if( !iContainment_80p
            && i_az < fGraph_AngularResolutionKingSigma.size()
            && fGraph_AngularResolutionKingSigma[i_az] )
    {
    
        // fill sigma parameter
        double i_emin = 0.;
        double i_emax = 0.;
        double x = 0.;
        double y = 0.;
        for( int i = 0; i < fGraph_AngularResolutionKingSigma[i_az]->GetN(); i++ )
        {
            fGraph_AngularResolutionKingSigma[i_az]->GetPoint( i, x, y );
            if( x < i_emin )
            {
                i_emin = x;
            }
            if( x > i_emax )
            {
                i_emax = x;
            }
        }
        fGraph_AngularResolutionKingSigma[i_az]->GetPoint( 0, i_emin, y );
        fGraph_AngularResolutionKingSigma[i_az]->GetPoint( fGraph_AngularResolutionKingSigma[i_az]->GetN(), i_emax, y );
        for( int i = 0; i < nbins; i++ )
        {
            if( e0[i] > i_emin && e0[i] < i_emax )
            {
                Rec_angRes_kingSigma[i] = fGraph_AngularResolutionKingSigma[i_az]->Eval( e0[i] );
            }
        }
        
        // fill gamma parameter
        i_emin = 0.;
        i_emax = 0.;
        x = 0.;
        y = 0.;
        for( int i = 0; i < fGraph_AngularResolutionKingGamma[i_az]->GetN(); i++ )
        {
            fGraph_AngularResolutionKingGamma[i_az]->GetPoint( i, x, y );
            if( x < i_emin )
            {
                i_emin = x;
            }
            if( x > i_emax )
            {
                i_emax = x;
            }
        }
        fGraph_AngularResolutionKingGamma[i_az]->GetPoint( 0, i_emin, y );
        fGraph_AngularResolutionKingGamma[i_az]->GetPoint( fGraph_AngularResolutionKingGamma[i_az]->GetN(), i_emax, y );
        for( int i = 0; i < nbins; i++ )
        {
            if( e0[i] > i_emin && e0[i] < i_emax )
            {
                Rec_angRes_kingGamma[i] = fGraph_AngularResolutionKingGamma[i_az]->Eval( e0[i] );
            }
        }
        
    }
}

void VEffectiveAreaCalculator::setAngularResolution2D( unsigned int i_az, vector< TH2D* > h )
{
    if( i_az < hVAngularLogDiffEmc_2D.size() && h.size() == 4 )
    {
        hVAngularLogDiffEmc_2D[i_az] = h[3];
    }
}

void VEffectiveAreaCalculator::setAngularResolutionGraph( unsigned int i_az, TGraphErrors* g, bool iAngContainment_80p )
{
    if( i_az < fGraph_AngularResolution68p.size() && g )
    {
        if( iAngContainment_80p )
        {
            fGraph_AngularResolution80p[i_az] = g;
        }
        else
        {
            fGraph_AngularResolution68p[i_az] = g;
        }
    }
}

void VEffectiveAreaCalculator::setAngularResolutionKingSigmaGraph( unsigned int i_az, TGraphErrors* g )
{
    if( i_az < fGraph_AngularResolutionKingSigma.size() && g )
    {
        fGraph_AngularResolutionKingSigma[i_az] = g;
    }
}

void VEffectiveAreaCalculator::setAngularResolutionKingGammaGraph( unsigned int i_az, TGraphErrors* g )
{
    if( i_az < fGraph_AngularResolutionKingGamma.size() && g )
    {
        fGraph_AngularResolutionKingGamma[i_az] = g;
    }
}

void VEffectiveAreaCalculator::fillEcutSub( double iE, enum E_HIS1D iCutIndex )
{
    if( hV_HIS1D.find( iCutIndex ) != hV_HIS1D.end() )
    {
        for( unsigned int i_az = 0; i_az < fVMinAz.size(); i_az++ )
        {
            for( unsigned int s = 0; s < fVSpectralIndex.size(); s++ )
            {
                if( s < hV_HIS1D[iCutIndex].size()
                        && i_az < hV_HIS1D[iCutIndex][s].size() )
                {
                    hV_HIS1D[iCutIndex][s][i_az]->Fill( iE, 1. );
                }
            }
        }
    }
}

/*
 * fill tree with cut numbers and MVA values per event
 * (same cut numbering as in VEffectiveAreaCalculator::fillEcutSub()
 *
 * 1    hhEcutTrigger_R (color: 1, marker: 20)
 * 2    hhEcutFiducialArea_R (color: 2, marker: 20)
 * 3    hhEcutStereoQuality_R (color: 3, marker: 20)
 * 4    hhEcutTelType_R (color: 4, marker: 20)
 * 5    hhEcutDirection_R (color: 5, marker: 20)
 * 6    hhEcutEnergyReconstruction_R (color: 6, marker: 20)
 * 7    hhEcutGammaHadron_R (color: 7, marker: 20)
 * 11   --> removed by testAzimuthInterval cut
 *
 *  1. Events passing gamma/hadron separation cut and direction cut
 *  fDL2EventTree->Draw("MVA", "Class==5" );
 *
 *  2. Events passing gamma/hadron separation cut and not direction cut
 *  fDL2EventTree->Draw("MVA", "Class==0" );
 *
 *  3. Events before gamma/hadron separation cut and before direction cut
 *   fDL2EventTree->Draw("MVA", "Class==0||Class==7||Class==5", "");
 *
 */
void VEffectiveAreaCalculator::fillDL2EventDataTree( CData* c, UChar_t iCutClass, float iMVA )
{
    // apply minimal quality cuts:
    // - require successful direction and energy reconstruction
    // - only applied for DL2 cases, not when data tree is also
    //   copied (otherwise event numbers don't fit)
    if( !fDL2WriteFullEventTree )
    {
        if( c->Xoff < -90.
                || c->Yoff < -90.
                || c->Ze < 0.
                || c->ErecS < 0. )
        {
            return;
        }
    }
    
    if( fDL2EventTree && c )
    {
        fDL2_runNumber = ( UInt_t )c->runNumber;
        fDL2_eventNumber = ( UInt_t )c->eventNumber;
        fDL2_MCaz = c->MCaz;
        fDL2_MCel = 90. - c->MCze;
        fDL2_MCxoff = c->MCxoff;
        fDL2_MCyoff = c->MCyoff;
        fDL2_MCe0 = c->MCe0;
        fDL2_ArrayPointing_Azimuth = c->ArrayPointing_Azimuth;
        fDL2_ArrayPointing_Elevation = c->ArrayPointing_Elevation;
        fDL2_az = c->Az;
        fDL2_el = 90. - c->Ze;
        fDL2_xoff = c->Xoff;
        fDL2_yoff = c->Yoff;
        fDL2_erec = c->ErecS;
        fDL2_nimages = ( UChar_t )c->NImages;
        fDL2_Cut_Class = iCutClass;
        if( fDL2_extendedTrees )
        {
            fDL2_Xcore = c->Xcore;
            fDL2_Ycore = c->Ycore;
            fDL2_Xoff_intersect = c->Xoff_intersect;
            fDL2_Yoff_intersect = c->Yoff_intersect;
            fDL2_img2_ang = c->img2_ang;
            fDL2_SizeSecondMax = c->SizeSecondMax;
            fDL2_NTelPairs = ( UChar_t )c->NTelPairs;
            fDL2_MSCW = c->MSCW;
            fDL2_MSCL = c->MSCL;
            fDL2_EmissionHeight = c->EmissionHeight;
            fDL2_EmissionHeightChi2 = c->EmissionHeightChi2;
            fDL2_DispDiff = c->DispDiff;
            fDL2_dESabs = c->getEnergyDelta();
            fDL2_NTrig = ( UChar_t )c->NTrig;
            fDL2_meanPedvar_Image = c->meanPedvar_Image;
            for( int i = 0; i < c->NImages; i++ )
            {
                fDL2_size[i] = c->size[i];
                fDL2_dist[i] = c->dist[i];
                fDL2_loss[i] = c->loss[i];
                fDL2_fui[i] = c->fui[i];
                fDL2_cross[i] = c->cross[i];
                fDL2_asym[i] = c->asym[i];
                fDL2_tgrad_x[i] = c->tgrad_x[i];
                fDL2_R[i] = c->R[i];
                fDL2_ES[i] = c->ES[i];
            }
        }
        fDL2_Cut_MVA = iMVA;
        fDL2EventTree->Fill();
    }
}


void VEffectiveAreaCalculator::addMeanResponseMatrix( vector <float> i_emc, vector <float> i_erec , vector <float> i_erec_err )
{

    // Getting binning
    float i_binw = i_emc[1] - i_emc[0] ;
    float* i_bins = new float[i_emc.size() + 1];
    vector <float> i_binc( i_emc.size(), 0 );
    
    i_bins[0] = i_emc[0] - i_binw / 2.;
    
    for( unsigned int i = 0; i < i_emc.size() ; i++ )
    {
        i_binc[i] = i_emc[i];
        i_bins[i + 1] = i_bins[i] + i_binw;
    }
    
    // cout << i_emc.size() << " bins: " << i_emc[0] << " - " << i_emc[i_emc.size() -1] << endl;
    int i_nbins = i_binc.size();
    
    TH2F* i_hist = new TH2F( "i_tmp", "i_tmp", i_nbins , i_bins, i_nbins, i_bins );
    TF1* i_gaussian = new TF1( "i_gaussian", "gaus", -2, 2.5 );
    
    // Filling Histograms
    for( int  i = 0; i < i_hist->GetYaxis()->GetNbins(); i++ )
    {
    
        for( unsigned int j = 0; j < i_emc.size(); j++ )
        {
            // Assuming a gaussian shape
            i_gaussian->SetParameter( 0, 1 );
            i_gaussian->SetParameter( 1, i_erec[j] );
            i_gaussian->SetParameter( 2, i_erec_err[j] );
            
            // Filling bins
            if( fabs( i_hist->GetYaxis()->GetBinCenter( i + 1 ) -  i_emc[j] ) < 1.e-4 && fabs( i_emc[j] ) > 1.e-9 )
            {
                double j_tot = 0;
                
                for( int  k = 0; k < i_hist->GetXaxis()->GetNbins(); k++ )
                {
                    j_tot += i_gaussian->Eval( i_hist->GetXaxis()->GetBinCenter( k + 1 ) );
                    i_hist->SetBinContent( k + 1, i + 1,  i_gaussian->Eval( i_hist->GetXaxis()->GetBinCenter( k + 1 ) ) );
                    
                }
                for( int  k = 0; k < i_hist->GetXaxis()->GetNbins(); k++ )
                {
                    double i_entry = i_hist->GetBinContent( k + 1, i + 1 );
                    // j_tot += i_gaussian->Eval( i_hist->GetXaxis()->GetBinCenter(k+1) );
                    i_hist->SetBinContent( k + 1, i + 1,  i_entry / j_tot );
                    
                }
            }
        }
        
        
        
    }
    
    // Normalization not required until final stage
    // Hist scaled by NEntries (normalization will div by (n+const))
    // VHistogramUtilities::normalizeTH2D_y(i_hist);
    
    // Checking if mean histogram exists
    VHistogramUtilities::normalizeTH2D_y( i_hist );
    if( !hMeanResponseMatrix )
    {
        cout << "\t\tVEffectiveAreaCalculator::addMeanResponseMatrix Creating new histogram" << endl;
        
        hMeanResponseMatrix = new TH2D( "hMeanResponseMatrix", "hMeanResponseMatrix", i_nbins , i_bins, i_nbins, i_bins );
        hMeanResponseMatrix->Add( i_hist );
        
    }
    
    
    else
    {
        // cout << "VEffectiveAreaCalculator::addMeanResponseMatrix Adding Hist" << endl;
        hMeanResponseMatrix->Add( i_hist );
        VHistogramUtilities::normalizeTH2D_y( hMeanResponseMatrix );
    }
    
    delete i_hist;
    delete[] i_bins;
    // fNMeanResponseMatrix++ ;
}

/*
 * conversion from enum count for histograms to
 * histogram names
 *
 * (no better solution found)
*/
string VEffectiveAreaCalculator::getEffectiveAreaNamefromEnumInt( int iHisID, string iType )
{
    if( iType == "1D" )
    {
        switch( iHisID )
        {
            case E_Emc:
                return "hEmc";
            case E_Ecut:
                return "hEcut";
            case E_EcutUW:
                return "hEcutUW";
            case E_EcutNoTh2:
                return "hEcutNoTh2";
            case E_Ecut500:
                return "hEcut500";
            case E_EcutRec:
                return "hEcutRec";
            case E_EcutRecUW:
                return "hEcutRecUW";
            case E_EcutRecNoTh2:
                return "hEcutRecNoTh2";
            case E_WeightedRate:
                return "hWeightedRate";
            case E_WeightedRate005:
                return "hWeightedRate005";
            case E_EcutTrigger:
                return "hEcutTrigger";
            case E_EcutFiducialArea:
                return "hEcutFiducialArea";
            case E_EcutStereoQuality:
                return "hEcutStereoQuality";
            case E_EcutTelType:
                return "hEcutTelType";
            case E_EcutDirection:
                return "hEcutDirection";
            case E_EcutEnergyReconstruction:
                return "hEcutEnergyReconstruction";
            case E_EcutGammaHadron:
                return "hEcutGammaHadron";
            case E_EmcUW:
                return "hEmcUW";
        }
    }
    else if( iType == "1P" )
    {
        switch( iHisID )
        {
            case E_EmcSWeight:
                return "hEmcSWeight";
            case E_EsysMCRelative:
                return "hEsysMCRelative";
        }
    }
    else if( iType == "2D" )
    {
        switch( iHisID )
        {
            case E_EsysMCRelativeRMS:
                return "hEsysMCRelativeRMS";
            case E_EsysMCRelative2D:
                return "hEsysMCRelative2D";
            case E_EsysMCRelative2DNoDirectionCut:
                return "hEsysMCRelative2DNoDirectionCut";
            case E_Esys2D:
                return "hEsys2D";
            case E_ResponseMatrix:
                return "hResponseMatrix";
            case E_ResponseMatrixFine:
                return "hResponseMatrixFine";
            case E_ResponseMatrixQC:
                return "hResponseMatrixQC";
            case E_ResponseMatrixFineQC:
                return "hResponseMatrixFineQC";
            case E_ResponseMatrixNoDirectionCut:
                return "hResponseMatrixNoDirectionCut";
            case E_ResponseMatrixFineNoDirectionCut:
                return "hResponseMatrixFineNoDirectionCut";
        }
    }
    
    return "";
}

bool VEffectiveAreaCalculator::newEffectiveAreaHistogram(
    string iType,
    int iHisN, string iHisTitle,
    string iTitleX, string iTitleY,
    int i_nbins_x, double i_xmin, double i_xmax,
    int i_nbins_y, double i_ymin, double i_ymax,
    string iPOpt )
{
    string iHName = getEffectiveAreaNamefromEnumInt( iHisN, iType );
    if( iHName.size() > 0 )
    {
        if( iType == "1D" )
        {
            h_HIS1D[iHisN] = new TH1D( iHName.c_str(),
                                       iHisTitle.c_str(), i_nbins_x,
                                       i_xmin, i_xmax );
            h_HIS1D[iHisN]->Sumw2();
            h_HIS1D[iHisN]->SetXTitle( iTitleX.c_str() );
            h_HIS1D[iHisN]->SetYTitle( iTitleY.c_str() );
            hisTreeList->Add( h_HIS1D[iHisN] );
            hisTreeListofHistograms->Add( h_HIS1D[iHisN] );
        }
        else if( iType == "1P" )
        {
            h_HIS1P[iHisN] = new TProfile( iHName.c_str(),
                                           iHisTitle.c_str(), i_nbins_x,
                                           i_xmin, i_xmax,
                                           i_ymin, i_ymax, iPOpt.c_str() );
            h_HIS1P[iHisN]->Sumw2();
            h_HIS1P[iHisN]->SetXTitle( iTitleX.c_str() );
            h_HIS1P[iHisN]->SetYTitle( iTitleY.c_str() );
            hisTreeList->Add( h_HIS1P[iHisN] );
            hisTreeListofHistograms->Add( h_HIS1P[iHisN] );
        }
        else if( iType == "2D" )
        {
            h_HIS2D[iHisN] = new TH2D( iHName.c_str(),
                                       iHisTitle.c_str(),
                                       i_nbins_x, i_xmin, i_xmax,
                                       i_nbins_y, i_ymin, i_ymax );
            h_HIS2D[iHisN]->SetXTitle( iTitleX.c_str() );
            h_HIS2D[iHisN]->SetYTitle( iTitleY.c_str() );
            hisTreeList->Add( h_HIS2D[iHisN] );
            hisTreeListofHistograms->Add( h_HIS2D[iHisN] );
        }
        return true;
    }
    return false;
}

void VEffectiveAreaCalculator::fillHistogram( int iHisType, int iHisN,
        unsigned int s, unsigned i_az,
        double i_x, double i_w, double i_w2D )
{
    if( iHisType == E_1D
            && hV_HIS1D.find( iHisN ) != hV_HIS1D.end()
            && hV_HIS1D[iHisN][s][i_az] )
    {
        hV_HIS1D[iHisN][s][i_az]->Fill( i_x, i_w );
    }
    else if( iHisType == E_1P
             && hV_HIS1P.find( iHisN ) != hV_HIS1P.end()
             && hV_HIS1P[iHisN][s][i_az] )
    {
        hV_HIS1P[iHisN][s][i_az]->Fill( i_x, i_w );
    }
    else if( iHisType == E_2D
             && hV_HIS2D.find( iHisN ) != hV_HIS2D.end()
             && hV_HIS2D[iHisN][s][i_az] )
    {
        hV_HIS2D[iHisN][s][i_az]->Fill( i_x, i_w, i_w2D );
    }
}


VGammaHadronCuts* VEffectiveAreaCalculator::getGammaHadronCuts( CData* c )
{
    if( !c )
    {
        return 0;
    }
    
    if( fCuts.size() == 1 )
    {
        return fCuts[0];
    }
    for( unsigned int i = 0; i < fCuts.size(); i++ )
    {
        if( fCuts[i] && fCuts[i]->useThisCut( c ) )
        {
            return fCuts[i];
        }
    }
    return 0;
}
