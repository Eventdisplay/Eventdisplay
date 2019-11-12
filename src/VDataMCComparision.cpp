/* \class VDataMCComparision
 *
 *
 * Note: extreme care needs to be taken when reading telescope parameters from the data tree:
 * - some variable arrays have lengths NImages
 * - some variable arrays have lengths NTel
 *
 */

#include "VDataMCComparision.h"

VDataMCComparisionHistogramData::VDataMCComparisionHistogramData( string iName, string iHistogramType, unsigned int iTelescopeID )
{
    fVarName = iName;
    fHistogramType = iHistogramType;
    fTelescopeID = iTelescopeID;
    fHis1D = 0;
    fHis2D_Erec = 0;
    fHis2D_ntubes = 0;
    fHis2D_size = 0;
    fHis2D_sizeHG = 0;
    fHis2D_sizeLG = 0;
    
}

TH2D* VDataMCComparisionHistogramData::newHistogram( string iName, string iYTitle, string iXTitle, int iNbins, double ix_min, double ix_max, double iy_min, double iy_max )
{
    char hname[200];
    if( fTelescopeID > 0 )
    {
        sprintf( hname, "h%s%s_%d_%s", fVarName.c_str(), iName.c_str(), fTelescopeID, fHistogramType.c_str() );
    }
    else
    {
        sprintf( hname, "h%s%s_%s", fVarName.c_str(), iName.c_str(), fHistogramType.c_str() );
    }
    
    TH2D* h = new TH2D( hname, "",  6, iy_min, iy_max, iNbins, ix_min, ix_max );
    h->SetXTitle( iYTitle.c_str() );
    h->SetYTitle( iXTitle.c_str() );
    h->Sumw2();
    h->SetLineColor( 2 );
    h->SetMarkerColor( 2 );
    
    return h;
}


bool VDataMCComparisionHistogramData::initHistogram( string iXTitle, int iNbins, double ix_min, double ix_max )
{
    char hname[200];
    
    if( fTelescopeID > 0 )
    {
        sprintf( hname, "h%s_%d_%s", fVarName.c_str(), fTelescopeID, fHistogramType.c_str() );
    }
    else
    {
        sprintf( hname, "h%s_%s", fVarName.c_str(), fHistogramType.c_str() );
    }
    
    fHis1D = new TH1D( hname, "", iNbins, ix_min, ix_max );
    fHis1D->SetXTitle( iXTitle.c_str() );
    fHis1D->Sumw2();
    
    // variable vs energy
    fHis2D_Erec = newHistogram( "Erec", "log_{10} energy_{rec} [TeV]", iXTitle, iNbins, ix_min, ix_max, -1., 1. );
    
    // variable vs ntubes
    fHis2D_ntubes = newHistogram( "ntubes", "ntubes", iXTitle, iNbins, ix_min, ix_max, 0., 60. );
    
    // variable vs size
    fHis2D_size = newHistogram( "size", "size", iXTitle, iNbins, ix_min, ix_max, 2., 5. );
    
    // variable vs sizeHG
    fHis2D_sizeHG = newHistogram( "sizeHG", "sizeHG", iXTitle, iNbins, ix_min, ix_max, 2., 5. );
    
    // variable vs sizeLG
    fHis2D_sizeLG = newHistogram( "sizeLG", "sizeLG", iXTitle, iNbins, ix_min, ix_max, 2., 5. );
    
    if( fHistogramType == "OFF" )
    {
        fHis1D->SetLineColor( 2 );
        fHis1D->SetMarkerColor( 2 );
    }
    
    return true;
}

/*
 * fill comparison histograms as function of energy, ntubes, size, sizeHG, sizeLG
 *
 */
void VDataMCComparisionHistogramData::fill( double iV, double iWeight, double iLogEnergy_TeV, int i_ntubes, double i_size, double i_fracLow )
{
    if( fHis1D )
    {
        fHis1D->Fill( iV, iWeight );
    }
    if( fHis2D_Erec && iLogEnergy_TeV > -98. )
    {
        fHis2D_Erec->Fill( iLogEnergy_TeV, iV, iWeight );
    }
    if( fHis2D_ntubes && i_ntubes > 0 )
    {
        fHis2D_ntubes->Fill( i_ntubes, iV, iWeight );
    }
    if( fHis2D_size && i_size > 2. )
    {
        fHis2D_size->Fill( log10( i_size ), iV, iWeight );
    }
    if( fHis2D_sizeHG && i_size - i_size * i_fracLow > 0. )
    {
        fHis2D_sizeHG->Fill( log10( i_size - i_size * i_fracLow ), iV, iWeight );
    }
    if( fHis2D_sizeLG && i_size * i_fracLow > 0. )
    {
        fHis2D_sizeLG->Fill( log10( i_size * i_fracLow ), iV, iWeight );
    }
}

////////////////////////////////////////////////////////////////////////////////////////

VDataMCComparision::VDataMCComparision( string iname, int intel, bool iMVAValues )
{
    fNTel = intel;
    fName = iname;
    if( fName != "ON" && fName != "OFF" && fName != "SIMS" && fName != "DIFF" )
    {
        cout << "error: unknown data type: " << fName << endl;
        cout << " allowed data types: ON, OFF, SIMS, DIFF" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    fData = 0;
    fCuts = 0;
    fCalculateMVAValues = iMVAValues;
    
    // spectral weighting (will be set later correctly, as it is run from MC run header)
    fSpectralWeight = new VSpectralWeight();
    //	fSpectralWeight->setMCParameter( 2.0, 0.05, 20. );
    fSpectralWeight->setMCParameter( 2.0, 0.003, 200. );
    // index MC events are weighted to
    fSpectralWeight->setSpectralIndex( 2.5 );
    
    // weighting in azimuth
    hAzWeight = 0;
    
    // setting all variables
    hisList = 0;
    
    // 2D histograms
    hXYcore = 0;
    hAzYcore = 0;
    hYt2 = 0;
    
    fAzRange = false;
    fAzMin = 0.;
    fAzMax = 0.;
    fZeMin = 0.;
    fZeMax = 90.;
    
    fWobbleNorth = 0.;
    fWobbleEast  = 0.;
    fWobbleFromDataTree = false;
    
    setShowerMaximZe_deg();
    
    defineHistograms();
}

/*

  needed only for the calculation of MVA value (not a default)

*/
void VDataMCComparision::initialGammaHadronCuts()
{
    fCuts = new VGammaHadronCuts();
    fCuts->initialize();
    fCuts->resetCutValues();
    // HARDWIRED CUT FILE
    if( !fCuts->readCuts( "$VERITAS_EVNDISP_AUX_DIR/GammaHadronCutFiles/ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-TMVA-BDT.dat" ) )
    {
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fCuts->initializeCuts();
    fCuts->setDataTree( fData );
    fCuts->printCutSummary();
}


void VDataMCComparision::resetTelescopeCoordinates()
{
    fTel_x.clear();
    fTel_y.clear();
    fTel_z.clear();
}

bool VDataMCComparision::setTelescopeCoordinates( double x, double y, double z )
{
    fTel_x.push_back( x );
    fTel_y.push_back( y );
    fTel_z.push_back( z );
    cout << "\t setting telescope coordinates for telescope " << fTel_x.size() << "\t" << x << "\t" << y << "\t" << z << endl;
    return true;
}

void VDataMCComparision::defineHistograms()
{
    hisList = new TList();
    char hname[200];
    
    double core_max = 350.;
    
    fHistoArray[ETHETA2] = new VDataMCComparisionHistogramData( "theta2", fName, 0 );
    fHistoArray[ETHETA2]->initHistogram( "#theta^{2} [deg^{2}]", 100, 0., 0.3 );
    
    fHistoArray[ELTHETA2] = new VDataMCComparisionHistogramData( "ltheta2", fName, 0 );
    fHistoArray[ELTHETA2]->initHistogram( "log_{10} #theta^{2} [deg^{2}]", 30, -5., 2. );
    
    fHistoArray[EMSCW] = new VDataMCComparisionHistogramData( "MSCW", fName, 0 );
    fHistoArray[EMSCW]->initHistogram( "mean reduced scaled width", 500, -5., 10. );
    
    fHistoArray[EMSCL] = new VDataMCComparisionHistogramData( "MSCL", fName, 0 );
    fHistoArray[EMSCL]->initHistogram( "mean reduced scaled length" , 500, -5., 10. );
    
    fHistoArray[EMWR] = new VDataMCComparisionHistogramData( "MWR", fName, 0 );
    fHistoArray[EMWR]->initHistogram( "mean scaled width", 500, -2., 2. );
    
    fHistoArray[EMLR] = new VDataMCComparisionHistogramData( "MLR", fName, 0 );
    fHistoArray[EMLR]->initHistogram( "mean scaled length", 500, -2., 2. );
    
    fHistoArray[EXCORE] = new VDataMCComparisionHistogramData( "Xcore", fName, 0 );
    fHistoArray[EXCORE]->initHistogram( "core position x [m]", 200, -1.*core_max, core_max );
    
    fHistoArray[EYCORE] = new VDataMCComparisionHistogramData( "Ycore", fName, 0 );
    fHistoArray[EYCORE]->initHistogram( "core position y [m]", 200, -1.*core_max, core_max );
    
    fHistoArray[EAEL] = new VDataMCComparisionHistogramData( "ArrayEl", fName, 0 );
    fHistoArray[EAEL]->initHistogram( "elevation [deg]", 180, 0., 90. );
    
    fHistoArray[EAAZ] = new VDataMCComparisionHistogramData( "ArrayAz", fName, 0 );
    fHistoArray[EAAZ]->initHistogram( "azimuth [deg]", 180, 0., 360. );
    
    fHistoArray[ENIMAGES] = new VDataMCComparisionHistogramData( "NImages", fName, 0 );
    fHistoArray[ENIMAGES]->initHistogram( "number of images", 5, 0., 5. );
    
    fHistoArray[EIMGSEL] = new VDataMCComparisionHistogramData( "ImgSel", fName, 0 );
    fHistoArray[EIMGSEL]->initHistogram( "telescope combinations", 16, 0., 16. );
    
    fHistoArray[EEMISSIONHEIGHT] = new VDataMCComparisionHistogramData( "EmissionHeight", fName, 0 );
    fHistoArray[EEMISSIONHEIGHT]->initHistogram( "estimated height of maximum emission [km]", 50, 0., 25. );
    
    fHistoArray[EEREC] = new VDataMCComparisionHistogramData( "Erec", fName, 0 );
    fHistoArray[EEREC]->initHistogram( "log_{10} energy_{MC}", 50, -2., 2. );
    
    fHistoArray[EPEDVAR] = new VDataMCComparisionHistogramData( "pedvar", fName, 0 );
    fHistoArray[EPEDVAR]->initHistogram( "mean pedvar", 100., 5., 10. );
    
    fHistoArray[EMVA] = new VDataMCComparisionHistogramData( "MVA", fName, 0 );
    fHistoArray[EMVA]->initHistogram( "MVA value #Tau", 100, -1., 1. );
    
    //these histograms are only filled for Model3D analysis
    
    fHistoArray[ESIGMAT3D] = new VDataMCComparisionHistogramData( "sigmaT3D", fName, 0 );
    fHistoArray[ESIGMAT3D]->initHistogram( "3D width [m]", 100, 0., 50. );
    
    fHistoArray[ENC3D] = new VDataMCComparisionHistogramData( "Nc3D", fName, 0 );
    fHistoArray[ENC3D]->initHistogram( "ln(Nc)", 100, 0., 20. );
    
    fHistoArray[EDEPTH3D] = new VDataMCComparisionHistogramData( "Depth3D", fName, 0 );
    fHistoArray[EDEPTH3D]->initHistogram( "depth of shower (3D) [g cm^{-2}]", 100, 0., 1000. );
    
    fHistoArray[ERWIDTH3D] = new VDataMCComparisionHistogramData( "RWidth3D", fName, 0 );
    fHistoArray[ERWIDTH3D]->initHistogram( "reduced 3D width [m]", 100, 0, 10. );
    
    fHistoArray[EERRRWIDTH3D] = new VDataMCComparisionHistogramData( "ErrRWidth3D", fName, 0 );
    fHistoArray[EERRRWIDTH3D]->initHistogram( "error in reduced 3D width [m]", 100, 0., 0.2 );
    
    // additional 2D histograms
    sprintf( hname, "hXYcore_%s", fName.c_str() );
    hXYcore = new TH2D( hname, "", 75, -1.*core_max, core_max, 75, -core_max, core_max );
    hXYcore->SetXTitle( "core position X [m]" );
    hXYcore->SetYTitle( "core position Y [m]" );
    hisList->Add( hXYcore );
    
    sprintf( hname, "hAzYcore_%s", fName.c_str() );
    hAzYcore = new TH2D( hname, "", 360, -180., 180., 300, -1.*core_max, core_max );
    hAzYcore->SetXTitle( "core position X [m]" );
    hAzYcore->SetYTitle( "core position Y [m]" );
    hisList->Add( hAzYcore );
    
    sprintf( hname, "hYt2_%s", fName.c_str() );
    hYt2 = new TH2D( hname, "", 100, -6., 0., 800, -200., 200. );
    hYt2->SetXTitle( "log_{10} #Theta^{2}" );
    hYt2->SetYTitle( "core position Y [m]" );
    hisList->Add( hYt2 );
    
    // telescope numbering starts at 1!
    for( int i = 1; i <= fNTel; i++ )
    {
        fHistoSingleTel[ELENGTH].push_back( new VDataMCComparisionHistogramData( "length", fName, i ) );
        fHistoSingleTel[ELENGTH].back()->initHistogram( "length [deg]",  100, 0., 0.5 );
        
        fHistoSingleTel[EWIDTH].push_back( new VDataMCComparisionHistogramData( "width", fName, i ) );
        fHistoSingleTel[EWIDTH].back()->initHistogram( "width [deg]",  100, 0., 0.25 );
        
        fHistoSingleTel[EDIST].push_back( new VDataMCComparisionHistogramData( "dist", fName, i ) );
        fHistoSingleTel[EDIST].back()->initHistogram( "dist [deg]",  100, 0., 2.0 );
        
        fHistoSingleTel[EALPHA].push_back( new VDataMCComparisionHistogramData( "alpha", fName, i ) );
        fHistoSingleTel[EALPHA].back()->initHistogram( "alpha [deg]",  100, 0., 25 );
        
        fHistoSingleTel[ENTUBES].push_back( new VDataMCComparisionHistogramData( "ntubes", fName, i ) );
        fHistoSingleTel[ENTUBES].back()->initHistogram( "ntubes",  100, 0., 200. );
        
        fHistoSingleTel[ENLOWGAIN].push_back( new VDataMCComparisionHistogramData( "nlowgain", fName, i ) );
        fHistoSingleTel[ENLOWGAIN].back()->initHistogram( "nlowgain",  100, 0., 100. );
        
        fHistoSingleTel[ESIZE].push_back( new VDataMCComparisionHistogramData( "size", fName, i ) );
        fHistoSingleTel[ESIZE].back()->initHistogram( "log_{10} size [d.c.]",  100, 1., 8. );
        
        fHistoSingleTel[ESIZEHG].push_back( new VDataMCComparisionHistogramData( "sizeHG", fName, i ) );
        fHistoSingleTel[ESIZEHG].back()->initHistogram( "log_{10} sizeHG [d.c.]",  100, 1., 8. );
        
        fHistoSingleTel[ESIZELG].push_back( new VDataMCComparisionHistogramData( "sizeLG", fName, i ) );
        fHistoSingleTel[ESIZELG].back()->initHistogram( "log_{10} sizeLG [d.c.]",  100, 1., 8. );
        
        fHistoSingleTel[EFRACLOW].push_back( new VDataMCComparisionHistogramData( "fraclow", fName, i ) );
        fHistoSingleTel[EFRACLOW].back()->initHistogram( "fraclow",  100, 0., 1. );
        
        fHistoSingleTel[EMAX1].push_back( new VDataMCComparisionHistogramData( "max1", fName, i ) );
        fHistoSingleTel[EMAX1].back()->initHistogram( "log_{10} size_{max1} [d.c.]",  100, 1., 4. );
        
        fHistoSingleTel[EMAX2].push_back( new VDataMCComparisionHistogramData( "max2", fName, i ) );
        fHistoSingleTel[EMAX2].back()->initHistogram( "log_{10} size_{max2} [d.c.]",  100, 1., 4. );
        
        fHistoSingleTel[EMAX3].push_back( new VDataMCComparisionHistogramData( "max3", fName, i ) );
        fHistoSingleTel[EMAX3].back()->initHistogram( "log_{10} size_{max3} [d.c.]",  100, 1., 4. );
        
        fHistoSingleTel[ELOSS].push_back( new VDataMCComparisionHistogramData( "loss", fName, i ) );
        fHistoSingleTel[ELOSS].back()->initHistogram( "loss",  100, 0., 1. );
        
        fHistoSingleTel[ELOS].push_back( new VDataMCComparisionHistogramData( "los", fName, i ) );
        fHistoSingleTel[ELOS].back()->initHistogram( "length/size [deg] x 1000",  100, 0., 1. );
        
        fHistoSingleTel[EASYM].push_back( new VDataMCComparisionHistogramData( "asym", fName, i ) );
        fHistoSingleTel[EASYM].back()->initHistogram( "asymmetry",  100, -2., 2. );
        
        fHistoSingleTel[ECENX].push_back( new VDataMCComparisionHistogramData( "cen_x", fName, i ) );
        fHistoSingleTel[ECENX].back()->initHistogram( "image centroid (x) [deg]",  100, -2., 2. );
        
        fHistoSingleTel[ECENY].push_back( new VDataMCComparisionHistogramData( "cen_y", fName, i ) );
        fHistoSingleTel[ECENY].back()->initHistogram( "image centroid (y) [deg]",  100, -2., 2. );
        
        fHistoSingleTel[ETGRADX].push_back( new VDataMCComparisionHistogramData( "tgrad_x", fName, i ) );
        fHistoSingleTel[ETGRADX].back()->initHistogram( "time gradient (x)",  100, -30., 30. );
        
        fHistoSingleTel[EMSCWT].push_back( new VDataMCComparisionHistogramData( "mscwt", fName, i ) );
        fHistoSingleTel[EMSCWT].back()->initHistogram( "expected width",  100, 0., 2. );
        
        fHistoSingleTel[EMWRT].push_back( new VDataMCComparisionHistogramData( "mwrt", fName, i ) );
        fHistoSingleTel[EMWRT].back()->initHistogram( "width/expected width",  100, 0., 2. );
        
        fHistoSingleTel[EMSCLT].push_back( new VDataMCComparisionHistogramData( "msclt", fName, i ) );
        fHistoSingleTel[EMSCLT].back()->initHistogram( "expected length",  100, 0., 2. );
        
        fHistoSingleTel[EMLTT].push_back( new VDataMCComparisionHistogramData( "mltt", fName, i ) );
        fHistoSingleTel[EMLTT].back()->initHistogram( "length/expected length",  100, 0., 2. );
        
        fHistoSingleTel[ETELDIST].push_back( new VDataMCComparisionHistogramData( "r", fName, i ) );
        fHistoSingleTel[ETELDIST].back()->initHistogram( "distance telescope - core [m]",  100, 0., 400. );
        
        fHistoSingleTel[ERECRAT].push_back( new VDataMCComparisionHistogramData( "erecratio", fName, i ) );
        fHistoSingleTel[ERECRAT].back()->initHistogram( "E_{rec,T}/E_{rec}",  100, 0., 5. );
        
        fHistoSingleTel[EPEDVART].push_back( new VDataMCComparisionHistogramData( "pedvarT", fName, i ) );
        fHistoSingleTel[EPEDVART].back()->initHistogram( "pedvar (T)",  100, 5., 10. );
        
        sprintf( hname, "hcen_xy%d_%s", i, fName.c_str() );
        hcen_xy.push_back( new TH2D( hname, "", 100, -2., 2., 100, -2., 2. ) );
        hcen_xy.back()->SetXTitle( "image centroid (x) [deg]" );
        hcen_xy.back()->SetYTitle( "image centroid (y) [deg]" );
        hTel2D.push_back( hcen_xy.back() );
        hisList->Add( hcen_xy.back() );
        
        sprintf( hname, "hdistR%d_%s", i, fName.c_str() );
        hdistR.push_back( new TH2D( hname, "", 20, 0., 300., 20, 0., 2. ) );
        sprintf( hname, "distance to T%d [m]", i );
        hdistR.back()->SetXTitle( hname );
        hdistR.back()->SetYTitle( "local distance [deg]" );
        hisList->Add( hdistR.back() );
        hTel2D.push_back( hdistR.back() );
    }
    
    // add all histograms to list of histograms (for 1D and 2D)
    for( map <E_varname, vector< VDataMCComparisionHistogramData* > >::iterator it = fHistoSingleTel.begin(); it != fHistoSingleTel.end(); it++ )
    {
        for( unsigned int i = 0; i < it->second.size(); i++ )
        {
            hisList->Add( it->second[i]->fHis1D );
            hisList->Add( it->second[i]->fHis2D_Erec );
            hisList->Add( it->second[i]->fHis2D_ntubes );
            hisList->Add( it->second[i]->fHis2D_size );
            hisList->Add( it->second[i]->fHis2D_sizeHG );
            hisList->Add( it->second[i]->fHis2D_sizeLG );
            hTel.push_back( it->second[i]->fHis1D );
            hTel2D.push_back( it->second[i]->fHis2D_Erec );
            hTel2D.push_back( it->second[i]->fHis2D_ntubes );
            hTel2D.push_back( it->second[i]->fHis2D_size );
            hTel2D.push_back( it->second[i]->fHis2D_sizeHG );
            hTel2D.push_back( it->second[i]->fHis2D_sizeLG );
        }
    }
    for( map <E_varname, VDataMCComparisionHistogramData* >::iterator it = fHistoArray.begin(); it != fHistoArray.end(); it++ )
    {
        hisList->Add( it->second->fHis1D );
        hisList->Add( it->second->fHis2D_Erec );
        hisList->Add( it->second->fHis2D_ntubes );
        hisList->Add( it->second->fHis2D_size );
        hisList->Add( it->second->fHis2D_sizeHG );
        hisList->Add( it->second->fHis2D_sizeLG );
        hTel.push_back( it->second->fHis1D );
        hTel2D.push_back( it->second->fHis2D_Erec );
        hTel2D.push_back( it->second->fHis2D_ntubes );
        hTel2D.push_back( it->second->fHis2D_size );
        hTel2D.push_back( it->second->fHis2D_sizeHG );
        hTel2D.push_back( it->second->fHis2D_sizeLG );
    }
    
}

/*
 * norm applies on off runs (mult.)
 */

bool VDataMCComparision::setOnOffHistograms( VDataMCComparision* on, VDataMCComparision* off, double norm )
{
    if( !on || !off )
    {
        cout << "on or off not defined" << endl;
        return false;
    }
    norm *= -1.;
    
    hXYcore->Add( on->hXYcore, off->hXYcore, 1., norm );
    hAzYcore->Add( on->hAzYcore, off->hAzYcore, 1., norm );
    hYt2->Add( on->hYt2, off->hYt2, 1., norm );
    
    // 1D histograms
    for( unsigned int j = 0; j < hTel.size(); j++ )
    {
        hTel[j]->Add( on->hTel[j], off->hTel[j], 1., norm );
    }
    // 2d histograms
    for( unsigned int j = 0; j < hTel2D.size(); j++ )
    {
        hTel2D[j]->Add( on->hTel2D[j], off->hTel2D[j], 1., norm );
    }
    
    
    return true;
}

void VDataMCComparision::setEntries( TH1D* iH )
{
    double ie = 0.;
    for( int i = 1; i <= iH->GetNbinsX(); i++ )
    {
        if( iH->GetBinContent( i ) > 0. )
        {
            ie += iH->GetBinContent( i );
        }
    }
    
    iH->SetEntries( ie );
}

void VDataMCComparision::setEntries( TH2D* iH )
{
    double ie = 0.;
    for( int i = 1; i <= iH->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= iH->GetNbinsY(); j++ )
        {
            if( iH->GetBinContent( i, j ) > 0. )
            {
                ie += iH->GetBinContent( i, j );
            }
        }
    }
    iH->SetEntries( ie );
}


/*
 * fill all histograms for both simulations and data
 *
 * stereo and single telescope parameters
 *
 * *** extreme care is needed: some telescope parameters have length of NTel, some NImages ***
 *
 */
bool VDataMCComparision::fillHistograms( string ifile, int iSingleTelescopeCuts )
{
    // quality cuts
    double fCoreMax_QC = 350; //350.;       // cut on core distance
    int    fNImages_min = 2;         // minimum number of images per event
    // stereo cuts
    double theta2_cut = 0.035; // 0.35;
    if( iSingleTelescopeCuts > 0 || iSingleTelescopeCuts == -2 )
    {
        theta2_cut = 0.2;
    }
    if( iSingleTelescopeCuts == -3 )
    {
        theta2_cut = 0.035; //0.035;
    }
    
    double theta2_min = 0.;
    double msw_max = 2.0;
    double msw_min = -1.2;
    double msl_max = 0.5;
    double msl_min = -1.2;
    
    /*
    double msw_max = 1.15;
    double msw_min = 0.05;
    double msl_max = 0.40;
    double msl_min = 0.05;
    */

    double size_min = 200.;
    double size2ndmax_min = 0.;
    
    // single telescope cuts
    int    ntubes_min = 4;
    double alpha_max = 8.;
    double length_min = 0.12;
    double length_max = 0.38;
    double width_max = 0.12;
    double width_min = 0.05;
    double los_max = 0.0003;
    double dist_min = -10.;
    double dist_max = 1.e10;
    
    // chain data files
    TChain* iC = new TChain( "data" );
    if( !iC->Add( ifile.c_str() ) )
    {
        cout << "error while reading data chain" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( iC->GetFile()->Get( "data" ) )
    {
        iC->SetTitle( iC->GetFile()->Get( "data" )->GetTitle() );
        
        // get MC data
        VMonteCarloRunHeader* iMC_H = ( VMonteCarloRunHeader* )iC->GetFile()->Get( "MC_runheader" );
        if( iMC_H )
        {
            if( fSpectralWeight )
            {
                fSpectralWeight->setMCParameter( -1.*iMC_H->spectral_index, iMC_H->E_range[0], iMC_H->E_range[1] );
                fSpectralWeight->print();
                
                // Deals with CARE sims without full information in the run header
                if( fSpectralWeight->getSpectralWeight( 1.0 ) == 0 )
                {
                    cout << "WARNING: Probably some error with the MC runheader" << endl;
                    cout << "WARNING: Setting  default values" << endl;
                    fSpectralWeight->setMCParameter( 2.0, 0.03, 200. );//Hard coded guess values
                    fSpectralWeight->print();
                }
            }
        }
    }
    
    // set this false for stereo cuts
    int fSingleTelescopeCuts = iSingleTelescopeCuts;
    
    if( fName == "SIMS" )
    {
        cout << "\t reading simulations..." << endl;
    }
    fData = new CData( iC, fName == "SIMS" );
    
    int nentries =  fData->fChain->GetEntries();
    cout << "\t entries: " << nentries << " (" << fNTel << " telescopes)" << endl;
    cout << "\t quality cuts: " << endl;
    cout << "\t\t maximum core distance [m]: " << fCoreMax_QC << endl;
    cout << "\t\t minimum number of images per event: " << fNImages_min << endl;
    if( fAzRange )
    {
        cout << "\t\t azimuth cut: [" << fAzMin << ", " << fAzMax << "]" << endl;
    }
    cout << "\t\t zenith cut: [" << fZeMin << ", " << fZeMax << "]" << endl;
    cout << "\t cuts: ";
    
    if( fSingleTelescopeCuts == -1 )
    {
        cout << " stereo cuts (hardwired)" << endl;
        cout << "\t " << msw_min << " < MSCW < " << msw_max << ", " << msl_min << " < MSCL < " << msl_max;
        cout << ", theta2 < " << theta2_cut << " deg2" << endl;
    }
    else if( fSingleTelescopeCuts == -2 )
    {
        cout << "NO CUTS (quality cuts only)" << endl;
    }
    else if( fSingleTelescopeCuts == -3 )
    {
        cout << " Theta2 cut (<" << theta2_cut << " deg2),  ";
        cout << " Size2ndMax cut (>" << size2ndmax_min << ")" << endl;
    }
    else
    {
        cout << " single telescope cuts for Telescope (hardwired): " << fSingleTelescopeCuts << endl;
        cout << "\t ntubes >" << ntubes_min << endl;
        cout << "\t alpha < " << alpha_max << endl;
        cout << "\t size < " << size_min << endl;
        cout << "\t " << length_min << " < length < " << length_max << endl;
        cout << "\t " << width_min << " < width < " << width_max << endl;
        cout << "\t los < " << los_max << endl;
        cout << "\t " << dist_min << " < dist < " << dist_max << endl;
    }
    
    double rdist1 = 0.;
    double theta2 = 0.;
    
    int iOldRun = 0;
    
    double weight = 1.;
    
    // code the different input sources
    // 0 = SIMS
    // 1 = ON
    // 2 = OFF
    int fInput = 0;
    if( fName == "SIMS" )
    {
        fInput = 0;
    }
    else if( fName == "ON" )
    {
        fInput = 1;
    }
    else if( fName == "OFF" )
    {
        fInput = 2;
    }
    
    ////////////////////////////////////////////////////////////
    // VGammaHadronCuts is needed for calulate of MVA values
    if( fCalculateMVAValues )
    {
        initialGammaHadronCuts();
    }
    
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // now loop over all events
    for( int i = 0; i < nentries; i++ )
    {
        fData->GetEntry( i );
        
        // define astro object (each time a new run starts)
        if( iOldRun != fData->runNumber && fWobbleFromDataTree )
        {
        
            TFile* iCurrentFile = iC->GetFile();
            if( iCurrentFile )
            {
                VEvndispRunParameter* iRunPar = ( VEvndispRunParameter* )iCurrentFile->Get( "runparameterV2" );
                if( iRunPar )
                {
                    fWobbleNorth = iRunPar->fWobbleNorth;
                    fWobbleEast  = iRunPar->fWobbleEast;
                    delete iRunPar;
                }
                else
                {
                    cout << "VDataMCComparision::fillHistograms: error reading runparameter for run ";
                    cout << fData->runNumber << endl;
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
            }
            else
            {
                cout << "VDataMCComparision::fillHistograms: error reading file for wobbles ";
                cout << iCurrentFile->GetName() << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            
            cout << "\t now at run " << fData->runNumber << " (" << fData->eventNumber << "):";
            cout << " N " << fWobbleNorth << ", E " << fWobbleEast << endl;
            
            iOldRun = fData->runNumber;
        }
        
        /////////////////////////////////////////////////
        /////////////////////////////////////////////////
        // quality cuts
        
        // nimage cut
        if( fData->NImages < fNImages_min )
        {
            continue;
        }
        // successfull energy reconstruction
        if( fData->EChi2S < 0 || fData->ErecS <= 0 )
        {
            continue;
        }
        // successfull lookup table reconstruction
        if( fData->MSCW < -50. || fData->MSCL < -50. )
        {
            continue;
        }
        // core inside a certain range
        if( sqrt( fData->Xcore * fData->Xcore + fData->Ycore * fData->Ycore ) > fCoreMax_QC )
        {
            continue;
        }
        // cut on image size
        if( fData->SizeSecondMax < size2ndmax_min )
        {
            continue;
        }
        
        // azimuth cuts (should be adjusted to match the azimuth distribution of the data runs)
        if( fAzRange && ( fData->Az < fAzMin || fData->Az > fAzMax ) )
        {
            continue;
        }
        // zenith cut (as we can use only one MC file)
        if( fData->Ze < fZeMin || fData->Ze > fZeMax )
        {
            continue;
        }
        
        /////////////////////////////////////////////////
        // direction cut
        theta2 = fData->theta2;
        // get correct theta2 for wobble runs
        // MC data
        if( fInput == 0 )
        {
            theta2 = ( fData->Yoff - fData->MCyoff ) * ( fData->Yoff - fData->MCyoff )
                     + ( fData->Xoff - fData->MCxoff ) * ( fData->Xoff - fData->MCxoff );
        }
        // data runs (on and off runs)
        else
        {
            float on_x = fData->Xoff_derot + fWobbleEast;
            float on_y = fData->Yoff_derot + fWobbleNorth;
            // off data
            if( fInput == 2 )
            {
                // loop over 5 off regions and use the smallest theta2
                theta2 = 1.e5;
                float off_theta = 90. * TMath::DegToRad();
                for( unsigned int th = 0; th < 5; th++ )
                {
                    float off_x = cos( off_theta ) * fWobbleEast - sin( off_theta ) * fWobbleNorth + fData->Xoff_derot;
                    float off_y = sin( off_theta ) * fWobbleEast + cos( off_theta ) * fWobbleNorth + fData->Yoff_derot;
                    if( off_x * off_x + off_y * off_y < theta2 )
                    {
                        theta2 = off_x * off_x + off_y * off_y;
                    }
                    off_theta += 45. * TMath::DegToRad();
                }
            }
            else
            {
                theta2 = on_x * on_x + on_y * on_y;
            }
        }
        
        // MC energy and spectral weight is filled for simulations only
        weight = 1.;
        if( fInput == 0 && fData->MCe0 > 0 )
        {
            if( fHistoArray[EEREC] )
            {
                fHistoArray[EEREC]->fill( log10( fData->MCe0 ), 1., -99., -99, -99., -99. );
            }
            if( fSpectralWeight )
            {
                weight = fSpectralWeight->getSpectralWeight( fData->MCe0 );
            }
            // weight in Az
            if( hAzWeight )
            {
                weight *= hAzWeight->GetBinContent( hAzWeight->FindBin( fData->Az ) );
            }
        }
        
        
        /////////////////////////////////////////////////////////
        //   ---    apply theta2 cuts ---
        /////////////////////////////////////////////////////////
        if( fSingleTelescopeCuts < -1 )
        {
            if( theta2 <= theta2_min || theta2 > theta2_cut )
            {
                continue;
            }
        }
        
        /////////////////////////////////////////////////////////
        //   ---    apply single telescope cuts  ----
        /////////////////////////////////////////////////////////
        //   apply single telescope cuts to the named telescope
        //   fSingleTelescopeCuts = telescope ID (counting starting at 1)
        if( fSingleTelescopeCuts > 0 )
        {
            // get telescope and check if telescope is available
            int iTelIndex = -99;
            for( int i = 0; i < fData->NImages; i++ )
            {
                if( fSingleTelescopeCuts - 1 == ( int )fData->ImgSel_list[i] )
                {
                    iTelIndex = i;
                    break;
                }
            }
            // required telescope did not participate
            if( iTelIndex < 0 )
            {
                continue;
            }
            // wobble cut
            if( fWobbleNorth != 0. || fWobbleEast != 0. )
            {
                if( theta2 <= theta2_min || theta2 > theta2_cut )
                {
                    continue;
                }
            }
            else
            {
                if( fData->alpha[fSingleTelescopeCuts - 1] > alpha_max )
                {
                    continue;
                }
            }
            if( fData->dist[iTelIndex] < dist_min || fData->dist[iTelIndex] > dist_max )
            {
                continue;
            }
            if( fData->size[iTelIndex]
                    && fData->length[iTelIndex] / fData->size[iTelIndex] > los_max )
            {
                continue;
            }
            
            if( fData->size[iTelIndex]
                    && fData->size[iTelIndex] < size_min )
            {
                continue;
            }
            
            if( fData->length[iTelIndex] < length_min ||  fData->length[iTelIndex] > length_max )
            {
                continue;
            }
            if( fData->width[iTelIndex] < width_min || fData->width[iTelIndex] > width_max )
            {
                continue;
            }
        }
        
        /////////////////////////////////////////////
        // fill histograms with all cuts applied
        
        if( fSingleTelescopeCuts != -1
                || ( theta2 >= 0. && fData->MSCW < msw_max && fData->MSCW > msw_min && fData->MSCL > msl_min && fData->MSCL < msl_max ) )
        {
            if( fHistoArray[ETHETA2] )
            {
                fHistoArray[ETHETA2]->fill( theta2, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[ELTHETA2] )
            {
                fHistoArray[ELTHETA2]->fill( log10( theta2 ), weight, log10( fData->ErecS ) );
            }
            hYt2->Fill( log10( theta2 ), fData->Ycore, weight );
        }
        if( fSingleTelescopeCuts != -1 || ( theta2 >= theta2_min && theta2 < theta2_cut &&  fData->MSCL < msl_max && fData->MSCL > msl_min ) )
        {
            if( fHistoArray[EMSCW] )
            {
                fHistoArray[EMSCW]->fill( fData->MSCW, weight, log10( fData->ErecS ) );
            }
            for( int j = 0; j < fData->NImages; j++ )
            {
                // telescope ID
                UInt_t iTelImage = fData->ImgSel_list[j];
                
                if( fData->MSCWT[iTelImage] > 0. && fHistoSingleTel[EMWRT].size() > 0
                        && iTelImage < fHistoSingleTel[EMWRT].size() && fHistoSingleTel[EMWRT][iTelImage] )
                {
                    fHistoSingleTel[EMWRT][iTelImage]->fill( fData->width[j] / fData->MSCWT[iTelImage],
                            weight, log10( fData->ErecS ),
                            fData->ntubes[j], fData->size[j],
                            fData->fraclow[iTelImage] );
                }
                if( fHistoSingleTel[EMSCWT].size() > 0 && iTelImage < fHistoSingleTel[EMWRT].size()
                        && fHistoSingleTel[EMSCWT][iTelImage] )
                {
                    fHistoSingleTel[EMSCWT][iTelImage]->fill( fData->MSCWT[iTelImage],
                            weight, log10( fData->ErecS ),
                            fData->ntubes[j], fData->size[j],
                            fData->fraclow[iTelImage] );
                }
            }
            if( fHistoArray[EMWR] )
            {
                fHistoArray[EMWR]->fill( fData->MWR, weight, log10( fData->ErecS ) );
            }
            // model3D variables
            if( fHistoArray[ESIGMAT3D] )
            {
                fHistoArray[ESIGMAT3D]->fill( fData->sigmaT3D, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[ENC3D] )
            {
                fHistoArray[ENC3D]->fill( fData->Nc3D, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[EDEPTH3D] )
            {
                fHistoArray[EDEPTH3D]->fill( fData->Depth3D, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[ERWIDTH3D] )
            {
                fHistoArray[ERWIDTH3D]->fill( fData->RWidth3D, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[EERRRWIDTH3D] )
            {
                fHistoArray[EERRRWIDTH3D]->fill( fData->ErrRWidth3D, weight, log10( fData->ErecS ) );
            }
            
        }
        if( fSingleTelescopeCuts != -1
                || ( theta2 >= theta2_min && theta2 < theta2_cut  &&  fData->MSCW < msw_max && fData->MSCW > msw_min ) )
        {
            if( fHistoArray[EMSCL] )
            {
                fHistoArray[EMSCL]->fill( fData->MSCL, weight, log10( fData->ErecS ) );
            }
            // Again: this is extremely bad as some index have a range of
            // NTel and others a range of NImages
            for( int j = 0; j < fData->NImages; j++ )
            {
                // telescope ID
                UInt_t iTelImage = fData->ImgSel_list[j];
                
                if( fData->MSCLT[j] > 0. && fHistoSingleTel[EMLTT].size() > 0
                        && iTelImage < fHistoSingleTel[EMLTT].size() && fHistoSingleTel[EMLTT][iTelImage] )
                {
                    fHistoSingleTel[EMLTT][iTelImage]->fill( fData->length[j] / fData->MSCLT[iTelImage],
                            weight, log10( fData->ErecS ),
                            fData->ntubes[j], fData->size[j],
                            fData->fraclow[iTelImage] );
                }
                if( fHistoSingleTel[EMSCLT].size() > 0 && iTelImage < fHistoSingleTel[EMSCLT].size()
                        && fHistoSingleTel[EMSCLT][iTelImage] )
                {
                    fHistoSingleTel[EMSCLT][j]->fill( fData->MSCLT[iTelImage], weight, log10( fData->ErecS ),
                                                      fData->ntubes[j], fData->size[j],
                                                      fData->fraclow[iTelImage] );
                }
            }
            if( fHistoArray[EMLR] )
            {
                fHistoArray[EMLR]->fill( fData->MLR, weight, log10( fData->ErecS ) );
            }
            if( fData->EmissionHeight > 0. && fHistoArray[EEMISSIONHEIGHT] )
            {
                fHistoArray[EEMISSIONHEIGHT]->fill( getCorrectedEmissionHeight( fData->EmissionHeight, 90. - fData->ArrayPointing_Elevation ),
                                                    weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[ENIMAGES] )
            {
                fHistoArray[ENIMAGES]->fill( fData->NImages, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[EPEDVAR] )
            {
                fHistoArray[EPEDVAR]->fill( fData->meanPedvar_Image, weight, log10( fData->ErecS ) );
            }
            if( fHistoArray[EIMGSEL] )
            {
                fHistoArray[EIMGSEL]->fill( fData->ImgSel, weight, log10( fData->ErecS ) );
            }
            
            if( fCuts )
            {
                fCuts->newEvent();
                fCuts->applyTMVACut( i );
                if( fCuts->getTMVA_EvaluationResult() > -90. && fHistoArray[EMVA] )
                {
                    fHistoArray[EMVA]->fill( fCuts->getTMVA_EvaluationResult(), weight, log10( fData->ErecS ) );
                }
            }
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        // stereo parameters
        if( fSingleTelescopeCuts != -1
                || ( theta2 >= theta2_min && theta2 < theta2_cut
                     && fData->MSCW < msw_max && fData->MSCW > msw_min && fData->MSCL < msl_max && fData->MSCL > msl_min ) )
        {
            if( fHistoArray[EXCORE] )
            {
                fHistoArray[EXCORE]->fill( fData->Xcore, weight, log10( fData->ErecS ) );
            }
            if( fInput == 0 )
            {
                if( fHistoArray[EYCORE] )
                {
                    fHistoArray[EYCORE]->fill( -1.*fData->Ycore, weight, log10( fData->ErecS ) );
                }
                hXYcore->Fill( fData->Xcore, -1.*fData->Ycore, weight );
                hAzYcore->Fill( fData->Az, -1.*fData->Ycore , weight );
            }
            else
            {
                if( fHistoArray[EYCORE] )
                {
                    fHistoArray[EYCORE]->fill( fData->Ycore, weight, log10( fData->ErecS ) );
                }
                hXYcore->Fill( fData->Xcore, fData->Ycore, weight );
                hAzYcore->Fill( fData->Az, fData->Ycore , weight );
            }
            if( fHistoArray[EAEL] )
            {
                fHistoArray[EAEL]->fill( fData->ArrayPointing_Elevation, weight );
            }
            if( fHistoArray[EAAZ] )
            {
                fHistoArray[EAAZ]->fill( fData->ArrayPointing_Azimuth, weight );
            }
            
            
            ///////////////////////////////////////////////////////////////
            // core position relative to each telescope
            for( int j = 0; j < fData->NImages; j++ )
            {
                // telescope ID
                UInt_t iTelImage = fData->ImgSel_list[j];
                
                if( fData->ntubes[j] > ntubes_min )
                {
                    rdist1 = VUtilities::line_point_distance( fData->Ycore, -1.*fData->Xcore, 0., fData->Ze, fData->Az,
                             fTel_y[iTelImage], -1.*fTel_x[iTelImage], fTel_z[iTelImage] );
                    // camera distance is calculated relative to centre of camera (should be: wobble offset position)
                    hdistR[j]->Fill( rdist1, sqrt( fData->cen_x[iTelImage]*fData->cen_x[iTelImage]
                                                   + fData->cen_y[iTelImage]*fData->cen_y[iTelImage] ), weight );
                }
            }
            
            ///////////////////////////////////////////////////////////////
            // single telescope distributions
            //  (loop over all images -> use selected images only)
            for( int t = 0; t < fData->NImages; t++ )
            {
                // get telescope number (0...3)
                int j = fData->ImgSel_list[t];
                
                // telescope wise quality cuts
                if( fData->ntubes[t] > ntubes_min && fData->size[t] > 0. && fData->size[t] > 0. && fData->ErecS > 0. )
                {
                    if( fHistoSingleTel[ELENGTH][j] )
                    {
                        fHistoSingleTel[ELENGTH][j]->fill( fData->length[t], weight, log10( fData->ErecS ),
                                                           fData->ntubes[t], fData->size[t],
                                                           fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EWIDTH][j] )
                    {
                        fHistoSingleTel[EWIDTH][j]->fill( fData->width[t], weight, log10( fData->ErecS ),
                                                          fData->ntubes[t], fData->size[t],
                                                          fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ENTUBES][j] )
                    {
                        fHistoSingleTel[ENTUBES][j]->fill( fData->ntubes[t], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                           fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ENLOWGAIN][j] && fData->nlowgain[j] )
                    {
                        fHistoSingleTel[ENLOWGAIN][j]->fill( fData->nlowgain[j], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                             fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EDIST][j] )
                    {
                        fHistoSingleTel[EDIST][j]->fill( fData->dist[t], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ESIZE][j] )
                    {
                        fHistoSingleTel[ESIZE][j]->fill( log10( fData->size[t] ), weight, log10( fData->ErecS ),
                                                         fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ESIZEHG][j] )
                    {
                        fHistoSingleTel[ESIZEHG][j]->fill( log10( fData->size[t] - fData->size[t]*fData->fraclow[j] ),
                                                           weight, log10( fData->ErecS ), fData->ntubes[t],
                                                           fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ESIZELG][j] && fData->fraclow[j] > 0. )
                    {
                        fHistoSingleTel[ESIZELG][j]->fill( log10( fData->size[t]*fData->fraclow[j] ), weight,
                                                           log10( fData->ErecS ), fData->ntubes[t],
                                                           fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EFRACLOW][j] )
                    {
                        fHistoSingleTel[EFRACLOW][j]->fill( fData->fraclow[j], weight, log10( fData->ErecS ),
                                                            fData->ntubes[t],
                                                            fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EMAX1][j] && fData->max1[j] > 0. )
                    {
                        fHistoSingleTel[EMAX1][j]->fill( log10( fData->max1[j] ), weight, log10( fData->ErecS ),
                                                         fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EMAX2][j] && fData->max2[j] > 0. )
                    {
                        fHistoSingleTel[EMAX2][j]->fill( log10( fData->max2[j] ), weight, log10( fData->ErecS ),
                                                         fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EMAX3][j] && fData->max3[j] > 0. )
                    {
                        fHistoSingleTel[EMAX3][j]->fill( log10( fData->max3[j] ), weight, log10( fData->ErecS ),
                                                         fData->ntubes[t], fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EALPHA][j] )
                    {
                        fHistoSingleTel[EALPHA][j]->fill( fData->alpha[j], weight, log10( fData->ErecS ),
                                                          fData->ntubes[t], fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ELOS][j] && fData->size[t] > 0. )
                    {
                        fHistoSingleTel[ELOS][j]->fill( fData->length[t] / fData->size[t] * 1000., weight,
                                                        log10( fData->ErecS ), fData->ntubes[t],
                                                        fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ELOSS][j] )
                    {
                        fHistoSingleTel[ELOSS][j]->fill( fData->loss[t], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EASYM][j] )
                    {
                        fHistoSingleTel[EASYM][j]->fill( fData->asym[j], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ECENX][j] )
                    {
                        fHistoSingleTel[ECENX][j]->fill( fData->cen_x[j], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ECENY][j] )
                    {
                        fHistoSingleTel[ECENY][j]->fill( fData->cen_y[j], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                         fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ETGRADX][j] )
                    {
                        fHistoSingleTel[ETGRADX][j]->fill( fData->tgrad_x[t], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                           fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ETELDIST][j] )
                    {
                        fHistoSingleTel[ETELDIST][j]->fill( fData->R[t], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                            fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ERECRAT][j] && fData->ErecS > 0. && fData->ES[t] > 0. )
                    {
                        fHistoSingleTel[ERECRAT][j]->fill( fData->ES[t] / fData->ErecS, weight, log10( fData->ErecS ), fData->ntubes[t],
                                                           fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[EPEDVART][j] )
                    {
                        fHistoSingleTel[EPEDVART][j]->fill( fData->meanPedvar_ImageT[j], weight, log10( fData->ErecS ), fData->ntubes[t],
                                                            fData->size[t], fData->fraclow[j] );
                    }
                    if( fHistoSingleTel[ERECRAT][j] && fData->ErecS > 0. )
                        if( hcen_xy[j] )
                        {
                            hcen_xy[j]->Fill( fData->cen_x[j], fData->cen_y[j], weight );
                        }
                }
            }
        }
    }
    return true;
}

void VDataMCComparision::scaleHistograms( string ifile )
{
    TFile f( ifile.c_str() );
    if( f.IsZombie() )
    {
        cout << "histogram scaling failed, no input file: " << ifile << endl;
        return;
    }
    TH1D* hEmc = ( TH1D* )f.Get( "hE0mc" );
    if( !hEmc )
    {
        cout << "histogram scaling failed, no mc histogram" << endl;
        return;
    }
    double iEntries = ( double )hEmc->GetEntries();
    
    TIter next( hisList );
    while( TH1* h = ( TH1* )next() )
    {
        if( iEntries > 0 )
        {
            h->Scale( 1. / iEntries );
        }
    }
}

bool VDataMCComparision::writeHistograms( string iOutFile )
{
    TFile iOut( iOutFile.c_str(), "UPDATE" );
    if( iOut.IsZombie() )
    {
        cout << "VDataMCComparision::writeHistograms: error opening output file: ";
        cout << iOutFile << endl;
        return false;
    }
    // iOut->cd();
    cout << "\t writing results (" << fName << ") to " << iOut.GetName() << endl;
    
    TIter next( hisList );
    while( TObject* h = next() )
    {
        h->Write();
    }
    if( hAzWeight )
    {
        hAzWeight->Write();
    }
    cout << "end writing histograms to disk" << endl;
    
    return true;
}

void VDataMCComparision::setZeRange( double iMin, double iMax )
{
    fZeMin = iMin;
    fZeMax = iMax;
}

void VDataMCComparision::setAzRange( double iMin, double iMax )
{
    fAzMin = 0.;
    fAzMax = 0.;
    fAzRange = false;
    
    if( TMath::Abs( iMin ) > 5. && TMath::Abs( iMax ) > 5. )
    {
        fAzMin = iMin;
        fAzMax = iMax;
        fAzRange = true;
    }
}

/*
 *    return a normalised histogram of the Az distribution of all events
 *
 *    typically filled for On data only
 *
 */
TH1D*  VDataMCComparision::getAzimuthWeightingHistogram( string ifile )
{
    char hname[100];
    sprintf( hname, "hAzWeight_%s", fName.c_str() );
    hAzWeight = new TH1D( hname, "", 360., 0., 360. );
    hAzWeight->SetXTitle( "azimuth [deg]" );
    
    
    // chain data files
    TChain* iC = new TChain( "data" );
    if( !iC->Add( ifile.c_str() ) )
    {
        cout << "error while reading data chain" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    CData* tData = new CData( iC, fName == "SIMS" );
    int nentries =  tData->fChain->GetEntries();
    cout << "Filling Az distributions for entries: " << nentries << endl;
    cout << "(requires loop over all ON files)" << endl;
    
    // now loop over all events
    for( int i = 0; i < nentries; i++ )
    {
        tData->GetEntry( i );
        
        // nimage cut
        if( tData->NImages < 2 )
        {
            continue;
        }
        // successfull energy reconstruction
        if( tData->EChi2S < 0 || tData->ErecS <= 0 )
        {
            continue;
        }
        // azimuth cuts (should be adjusted to match the azimuth distribution of the data runs)
        if( fAzRange && ( tData->Az < fAzMin || tData->Az > fAzMax ) )
        {
            continue;
        }
        // zenith cut (as we can use only one MC file)
        if( tData->Ze < fZeMin || tData->Ze > fZeMax )
        {
            continue;
        }
        
        hAzWeight->Fill( tData->Az );
    }
    
    if( hAzWeight->GetEntries() > 0. )
    {
        hAzWeight->Scale( 1. / hAzWeight->GetEntries() );
    }
    
    delete tData;
    delete iC;
    
    return hAzWeight;
    
}

/*
 * correct emission height for different zenith angle
 * of observations
 *
 */
double VDataMCComparision::getCorrectedEmissionHeight( double iEM, double iZe )
{
    double iCorr = 1.;
    if( cos( fShowerMaxZe_deg * TMath::DegToRad() ) != 0. )
    {
        iCorr = cos( iZe * TMath::DegToRad() ) / cos( fShowerMaxZe_deg * TMath::DegToRad() );
    }
    return iEM * iCorr;
}
