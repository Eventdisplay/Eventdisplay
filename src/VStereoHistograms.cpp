/*! \class VStereoHistograms
 *  \brief  holds all histograms for science analysis (stereo mode): sky maps, energy spectra, etc.
 *
 *
 */

#include "VStereoHistograms.h"

VStereoHistograms::VStereoHistograms( string i_hsuffix, double ibinsize, double ibinsizeUC, double iEnergyBinSize, double iTimeBinSize, double iTimeMin, double iTimeMax, bool iIsOn )
{
    bIsOn = iIsOn;
    fHisSuffix = i_hsuffix;
    fBinSize = ibinsize;
    fBinSizeUC = ibinsizeUC;
    fBinSizeEnergy = iEnergyBinSize;
    fBinSizeTime = iTimeBinSize;
    fTimeMin = iTimeMin;
    fTimeMax = iTimeMax;
    
    fRunNumber = -99;
    
    // this is probably overwritten later
    fSkyMapSizeXmin = -2.;
    fSkyMapSizeXmax =  2.;
    fSkyMapSizeYmin = -2.;
    fSkyMapSizeYmax =  2.;
    
    // histogram lists
    hisList = new TList();
    hListParameterHistograms = new TList();
    hListStereoParameterHistograms = new TList();
    hListRandomForestParameterHistograms = new TList();
    hListEnergyHistograms = new TList();
    hListSkyMaps = new TList();
    hListSkyMapsUC = new TList();
    hisRateList = new TList();
}

VStereoHistograms::~VStereoHistograms()
{
    if( hisList )
    {
        hisList->Delete();
    }
    if( hListParameterHistograms )
    {
        hListParameterHistograms->Delete();
    }
    if( hListStereoParameterHistograms )
    {
        hListStereoParameterHistograms->Delete();
    }
    if( hListRandomForestParameterHistograms )
    {
        hListRandomForestParameterHistograms->Delete();
    }
    if( hListEnergyHistograms )
    {
        hListEnergyHistograms->Delete();
    }
    if( hListSkyMaps )
    {
        hListSkyMaps->Delete();
    }
    if( hListSkyMapsUC )
    {
        hListSkyMapsUC->Delete();
    }
    if( hisRateList )
    {
        hisRateList->Delete();
    }
    
}

/*
 * set sky map bin size [deg]
 *
 */
void VStereoHistograms::setSkyMapSize( double xmin, double xmax, double ymin, double ymax )
{
    fSkyMapSizeXmin = xmin;
    fSkyMapSizeXmax = xmax;
    fSkyMapSizeYmin = ymin;
    fSkyMapSizeYmax = ymax;
}


void VStereoHistograms::defineHistograms()
{
    char i_key[800];
    char i_name[800];
    
    // points required in VStereoAnalysis::combineHistograms()
    h_combine_map_alpha_off = 0;
    h_combine_map_alpha_offUC = 0;
    h_combine_map_stereo_on = 0;
    h_combine_map_stereo_onUC = 0;
    h_combine_map_stereo_off = 0;
    h_combine_map_stereo_offUC = 0;
    
    ///////////////////////////////////
    // limits of 2D sky plots (uneven to get a bin exactly on the source position)
    double xmin = fSkyMapSizeXmin - fBinSize / 2.;
    double xmax = fSkyMapSizeXmax - fBinSize / 2.;
    double xminUC =  fSkyMapSizeXmin - fBinSizeUC / 2.;
    double xmaxUC =  fSkyMapSizeXmax - fBinSizeUC / 2.;
    
    double ymin = fSkyMapSizeYmin - fBinSize / 2.;
    double ymax = fSkyMapSizeYmax - fBinSize / 2.;
    double yminUC =  fSkyMapSizeYmin - fBinSizeUC / 2.;
    double ymaxUC =  fSkyMapSizeYmax - fBinSizeUC / 2.;
    
    if( fBinSize <= 0 )
    {
        cout << "VStereoHistograms::defineHistograms error: invalid binsize for stereo maps: " << fBinSize << endl;
        exit( EXIT_FAILURE );
    }
    int xbin   = ( int )( ( xmax - xmin ) / fBinSize + 0.5 );
    int xbinUC = ( int )( ( xmax - xmin ) / fBinSizeUC + 0.5 );
    int ybin   = ( int )( ( ymax - ymin ) / fBinSize + 0.5 );
    int ybinUC = ( int )( ( ymax - ymin ) / fBinSizeUC + 0.5 );
    
    ///////////////////////////////////
    // limits of energy histograms
    // log energy axis
    double i_emin = -1.9;
    double i_emax =  2.6;   // (clear that this is very high, but ensure that bin numbers can be combined)
    int    i_ebin = int( ( i_emax - i_emin ) / fBinSizeEnergy + 0.5 );
    // linear energy axis
    double i_Linemin = 0.05;
    double i_Linemax = 50.;
    int    i_Linebin = int( ( i_Linemax - i_Linemin ) / 0.01 );
    // 2D time binning
    //fill vector of time bins used for intra-run fluxes. Bins have width specified by user, except for the last bin which does not exceed the run duration.
    int i_tbin = 0;
    vector<double> tbins;
    if( fTimeMax > fTimeMin )
    {
        for( int i_b = 0; i_b * fBinSizeTime < fTimeMax - fTimeMin; i_b++ )
        {
            tbins.push_back( i_b * fBinSizeTime );
        }
        //last time bin is merged to second-to-last if the last bin is less than 10% of the nominal time bin length.
        if( tbins.size() > 1 && ( fTimeMax - fTimeMin ) - tbins.back() < 0.1 * fBinSizeTime )
        {
            tbins.back() = fTimeMax - fTimeMin;
        }
        else
        {
            tbins.push_back( fTimeMax - fTimeMin );
        }
        if( tbins.size() > 0 )
        {
            i_tbin = tbins.size() - 1;
        }
    }
    
    // distance to camera centre axis
    double i_t_offmin = 0.;
    double i_t_offmax = 4.;
    int i_t_offbin = TMath::Nint( ( i_t_offmax - i_t_offmin ) / 0.25 );
    
    sprintf( i_key, "htheta2_%s", fHisSuffix.c_str() );
    sprintf( i_name, "#theta^{2} Histogram (%s)", fHisSuffix.c_str() );
    htheta2 = new TH1D( i_key, i_name, 2000, 0., 1.0 );
    htheta2->SetXTitle( "#theta^{2} [deg^{2}]" );
    htheta2->SetYTitle( "No. of Events" );
    hisList->Add( htheta2 );
    hListParameterHistograms->Add( htheta2 );
    hListStereoParameterHistograms->Add( htheta2 );
    hListNameofParameterHistograms["htheta2"] = htheta2;
    
    sprintf( i_key, "hemiss_%s", fHisSuffix.c_str() );
    sprintf( i_name, "emission height Histogram (%s)", fHisSuffix.c_str() );
    hemiss = new TH1D( i_key, i_name, 500, 0., 500. );
    hemiss->SetXTitle( "mean emission height [km]" );
    hemiss->SetYTitle( "No. of Events" );
    hisList->Add( hemiss );
    hListParameterHistograms->Add( hemiss );
    hListStereoParameterHistograms->Add( hemiss );
    hListNameofParameterHistograms["hemiss"] = hemiss;
    
    sprintf( i_key, "hemissC2_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Chi2 of mean emission height Histogram (%s)", fHisSuffix.c_str() );
    hemissC2 = new TH1D( i_key, i_name, 500, 0., 5000. );
    hemissC2->SetXTitle( "Chi2 of mean emission height" );
    hemissC2->SetYTitle( "No. of Events" );
    hisList->Add( hemissC2 );
    hListParameterHistograms->Add( hemissC2 );
    hListStereoParameterHistograms->Add( hemissC2 );
    hListNameofParameterHistograms["hemissC2"] = hemissC2;
    
    sprintf( i_key, "hEChi2_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Chi2 of energy reconstruction Histogram (%s)", fHisSuffix.c_str() );
    herecChi2 = new TH1D( i_key, i_name, 500, 0., 500. );
    herecChi2->SetXTitle( "Chi2 of enery reconstruction" );
    herecChi2->SetYTitle( "No. of Events" );
    hisList->Add( herecChi2 );
    hListParameterHistograms->Add( herecChi2 );
    hListStereoParameterHistograms->Add( herecChi2 );
    hListNameofParameterHistograms["herecChi2"] = herecChi2;
    
    sprintf( i_key, "hmean_width_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Mean Width Histogram (%s)", fHisSuffix.c_str() );
    hmean_width = new TH1D( i_key, i_name, 100, 0., 0.5 );
    hmean_width->SetXTitle( "Mean Width (#circ)" );
    hmean_width->SetYTitle( "No. of Events" );
    hisList->Add( hmean_width );
    hListParameterHistograms->Add( hmean_width );
    hListStereoParameterHistograms->Add( hmean_width );
    hListNameofParameterHistograms["hmean_width"] = hmean_width;
    
    sprintf( i_key, "hmean_length_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Mean Length Histogram (%s)", fHisSuffix.c_str() );
    hmean_length = new TH1D( i_key, i_name, 100, 0., 1.0 );
    hmean_length->SetXTitle( "Mean Length (#circ)" );
    hmean_length->SetYTitle( "No. of Events" );
    hisList->Add( hmean_length );
    hListParameterHistograms->Add( hmean_length );
    hListStereoParameterHistograms->Add( hmean_length );
    hListNameofParameterHistograms["hmean_length"] = hmean_length;
    
    sprintf( i_key, "hmean_dist_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Mean Dist Histogram (%s)", fHisSuffix.c_str() );
    hmean_dist = new TH1D( i_key, i_name, 100, 0., 3.5 );
    hmean_dist->SetXTitle( "Mean Dist (#circ)" );
    hmean_dist->SetYTitle( "No. of Events" );
    hisList->Add( hmean_dist );
    hListParameterHistograms->Add( hmean_dist );
    hListStereoParameterHistograms->Add( hmean_dist );
    hListNameofParameterHistograms["hmean_dist"] = hmean_dist;
    
    //! setup the sky map histogram
    sprintf( i_key, "hmap_stereo_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Stereo Sky Map (%s)", fHisSuffix.c_str() );
    hmap_stereo = new TH2D( i_key, i_name, xbin, xmin, xmax, xbin, ymin, ymax );
    hmap_stereo->SetXTitle( "X-position on Sky (#circ)" );
    hmap_stereo->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_stereo );
    hListSkyMaps->Add( hmap_stereo );
    hListNameofSkyMaps.push_back( i_key );
    
    //! setup the sky map histogram for the background normalisation
    // this histogram must have the same binning definition as hmap_stereo !!!!
    sprintf( i_key, "hmap_alpha_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Alpha Sky Map (%s)", fHisSuffix.c_str() );
    hmap_alpha = new TH2D( i_key, i_name, xbin, xmin, xmax, ybin, ymin, ymax );
    hmap_alpha->SetXTitle( "X-position on Sky (#circ)" );
    hmap_alpha->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_alpha );
    hListSkyMaps->Add( hmap_alpha );
    hListNameofSkyMaps.push_back( i_key );
    
    sprintf( i_key, "hmap_alphaNorm_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Alpha Sky Map, norm (%s)", fHisSuffix.c_str() );
    hmap_alphaNorm = new TH2D( i_key, i_name, xbin, xmin, xmax, ybin, ymin, ymax );
    hmap_alphaNorm->SetXTitle( "X-position on Sky (#circ)" );
    hmap_alphaNorm->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_alphaNorm );
    hListSkyMaps->Add( hmap_alphaNorm );
    hListNameofSkyMaps.push_back( i_key );
    
    sprintf( i_key, "hmap_MeanSignalBackgroundAreaRatio_%s", fHisSuffix.c_str() );
    sprintf( i_name, "mean signal to background ratio (%s)", fHisSuffix.c_str() );
    hmap_MeanSignalBackgroundAreaRatio = new TH1D( i_key, i_name, 500, 0., 1. );
    hmap_MeanSignalBackgroundAreaRatio->SetXTitle( "mean signal to background ratio" );
    hisList->Add( hmap_MeanSignalBackgroundAreaRatio );
    hListSkyMaps->Add( hmap_MeanSignalBackgroundAreaRatio );
    hListNameofSkyMaps.push_back( i_key );
    
    
    //! setup the sky map histogram
    sprintf( i_key, "hmap_stereoUC_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Stereo Sky Map, uncorrelated (%s)", fHisSuffix.c_str() );
    hmap_stereoUC = new TH2D( i_key, i_name, xbinUC, xminUC, xmaxUC, ybinUC, yminUC, ymaxUC );
    hmap_stereoUC->SetXTitle( "X-position on Sky (#circ)" );
    hmap_stereoUC->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_stereoUC );
    hListSkyMapsUC->Add( hmap_stereoUC );
    hListNameofSkyMaps.push_back( i_key );
    
    //! setup the sky map histogram for the background normalisation
    // this histogram must have the same binning definition as hmap_stereo !!!!
    sprintf( i_key, "hmap_alphaUC_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Alpha Sky Map, uncorrelated (%s)", fHisSuffix.c_str() );
    hmap_alphaUC = new TH2D( i_key, i_name, xbinUC, xminUC, xmaxUC, ybinUC, yminUC, ymaxUC );
    hmap_alphaUC->SetXTitle( "X-position on Sky (#circ)" );
    hmap_alphaUC->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_alphaUC );
    hListSkyMapsUC->Add( hmap_alphaUC );
    hListNameofSkyMaps.push_back( i_key );
    
    sprintf( i_key, "hmap_MeanSignalBackgroundAreaRatioUC_%s", fHisSuffix.c_str() );
    sprintf( i_name, "mean signal to background ratio, uncorrelated (%s)", fHisSuffix.c_str() );
    hmap_MeanSignalBackgroundAreaRatioUC = new TH1D( i_key, i_name, 500, 0., 1. );
    hmap_MeanSignalBackgroundAreaRatioUC->SetXTitle( "mean signal to background ratio" );
    hisList->Add( hmap_MeanSignalBackgroundAreaRatioUC );
    hListSkyMaps->Add( hmap_MeanSignalBackgroundAreaRatioUC );
    hListNameofSkyMaps.push_back( i_key );
    
    sprintf( i_key, "hmap_alphaNormUC_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Alpha Sky Map, uncorrelated, norm (%s)", fHisSuffix.c_str() );
    hmap_alphaNormUC = new TH2D( i_key, i_name, xbinUC, xminUC, xmaxUC, ybinUC, yminUC, ymaxUC );
    hmap_alphaNormUC->SetXTitle( "X-position on Sky (#circ)" );
    hmap_alphaNormUC->SetYTitle( "Y-position on Sky (#circ)" );
    hisList->Add( hmap_alphaNormUC );
    hListSkyMapsUC->Add( hmap_alphaNormUC );
    hListNameofSkyMaps.push_back( i_key );
    
    //! setup the ground plane core map histogram
    sprintf( i_key, "hcore_%s", fHisSuffix.c_str() );
    sprintf( i_name, "Core Position Map (%s)", fHisSuffix.c_str() );
    hcore = new TH2D( i_key, i_name, 100, -200, 200, 100, -200, 200 );
    hcore->SetXTitle( "X-position of Core on Ground (m)" );
    hcore->SetYTitle( "Y-position of Core on Ground (m)" );
    hisList->Add( hcore );
    hListParameterHistograms->Add( hcore );
    hListStereoParameterHistograms->Add( hcore );
    hListNameofParameterHistograms["hcore"] = hcore;
    
    sprintf( i_key, "hmsc_%s", fHisSuffix.c_str() );
    sprintf( i_name, "MSC Histogram (%s)", fHisSuffix.c_str() );
    hmsc = new TH2D( i_key, i_name, 200, -5., 10., 200, -5., 10. );
    hmsc->SetXTitle( "mean scaled width" );
    hmsc->SetYTitle( "mean scaled length" );
    hisList->Add( hmsc );
    hListParameterHistograms->Add( hmsc );
    hListStereoParameterHistograms->Add( hmsc );
    hListNameofParameterHistograms["hmsc"] = hmsc;
    
    //! setup mscw histogram
    sprintf( i_key, "hmscw_%s", fHisSuffix.c_str() );
    sprintf( i_name, "MSCW Histogram (%s)", fHisSuffix.c_str() );
    hmscw = new TH1D( i_key, i_name, 200, -5., 10. );
    hmscw->SetXTitle( "mean scaled width" );
    hmscw->SetYTitle( "No. of Events" );
    hisList->Add( hmscw );
    hListParameterHistograms->Add( hmscw );
    hListStereoParameterHistograms->Add( hmscw );
    hListNameofParameterHistograms["hmscw"] = hmscw;
    
    //! setup mscl histogram
    sprintf( i_key, "hmscl_%s", fHisSuffix.c_str() );
    sprintf( i_name, "MSCL Histogram(%s)", fHisSuffix.c_str() );
    hmscl = new TH1D( i_key, i_name, 200, -5., 10. );
    hmscl->SetXTitle( "mean scaled length" );
    hmscl->SetYTitle( "No. of Events" );
    hisList->Add( hmscl );
    hListParameterHistograms->Add( hmscl );
    hListStereoParameterHistograms->Add( hmscl );
    hListNameofParameterHistograms["hmscl"] = hmscl;
    
    //! setup rf histogram
    sprintf( i_key, "hrf_%s", fHisSuffix.c_str() );
    sprintf( i_name, "rf Histogram (%s)", fHisSuffix.c_str() );
    hrf = new TH1D( i_key, i_name, 200, 0., 1. );
    hrf->SetXTitle( "signal classification (rf)" );
    hrf->SetYTitle( "No. of Events" );
    hisList->Add( hrf );
    hListParameterHistograms->Add( hrf );
    hListRandomForestParameterHistograms->Add( hrf );
    hListNameofParameterHistograms["hrf"] = hrf;
    
    sprintf( i_key, "herecCounts2D_vs_distance_%s", fHisSuffix.c_str() );
    sprintf( i_name, "counting histogram (energy vs distance to camera centre [deg]) (%s)", fHisSuffix.c_str() );
    herecCounts2D_vs_distance = new TH2D( i_key, i_name, i_ebin, i_emin, i_emax, i_t_offbin, i_t_offmin, i_t_offmax );
    herecCounts2D_vs_distance->SetXTitle( "log_{10} energy [TeV]" );
    herecCounts2D_vs_distance->SetZTitle( "No. of events" );
    herecCounts2D_vs_distance->SetYTitle( "distance to camera centre [deg]" );
    hisList->Add( herecCounts2D_vs_distance );
    hListParameterHistograms->Add( herecCounts2D_vs_distance );
    hListEnergyHistograms->Add( herecCounts2D_vs_distance );
    hListNameofParameterHistograms["herecCounts2D_vs_distance"] = herecCounts2D_vs_distance;
    
    
    sprintf( i_key, "herecCounts_%s", fHisSuffix.c_str() );
    sprintf( i_name, "counting histogram (energy) (%s)", fHisSuffix.c_str() );
    herecCounts = new TH1D( i_key, i_name, i_ebin, i_emin, i_emax );
    herecCounts->SetXTitle( "log_{10} energy [TeV]" );
    herecCounts->SetYTitle( "No. of events" );
    hisList->Add( herecCounts );
    hListParameterHistograms->Add( herecCounts );
    hListEnergyHistograms->Add( herecCounts );
    hListNameofParameterHistograms["herecCounts"] = herecCounts;
    
    sprintf( i_key, "herecWeights_%s", fHisSuffix.c_str() );
    sprintf( i_name, "effective area vs. raw energy (%s)", fHisSuffix.c_str() );
    herecWeights = new TH2D( i_key, i_name, i_ebin, i_emin, i_emax, 140, 1., 7. );
    herecWeights->SetXTitle( "log_{10} energy [TeV]" );
    herecWeights->SetYTitle( "log_{10} 1/eweights" );
    hisList->Add( herecWeights );
    hListParameterHistograms->Add( herecWeights );
    hListEnergyHistograms->Add( herecWeights );
    hListNameofParameterHistograms["herecWeights"] = herecWeights;
    
    sprintf( i_key, "hLinerecCounts_%s", fHisSuffix.c_str() );
    sprintf( i_name, "counting histogram (energy) (%s)", fHisSuffix.c_str() );
    hLinerecCounts = new TH1D( i_key, i_name, i_Linebin, i_Linemin, i_Linemax );
    hLinerecCounts->SetXTitle( "energy [TeV]" );
    hLinerecCounts->SetYTitle( "No. of events" );
    hisList->Add( hLinerecCounts );
    hListParameterHistograms->Add( hLinerecCounts );
    hListEnergyHistograms->Add( hLinerecCounts );
    hListNameofParameterHistograms["hLinerecCounts"] = hLinerecCounts;
    
    sprintf( i_key, "hLinerecWeights_%s", fHisSuffix.c_str() );
    sprintf( i_name, "effective area vs. raw energy (%s)", fHisSuffix.c_str() );
    hLinerecWeights = new TH2D( i_key, i_name, i_Linebin, i_Linemin, i_Linemax, 140, 1., 7. );
    hLinerecWeights->SetXTitle( "energy [TeV]" );
    hLinerecWeights->SetYTitle( "log_{10} 1/eweights" );
    hisList->Add( hLinerecWeights );
    hListParameterHistograms->Add( hLinerecWeights );
    hListEnergyHistograms->Add( hLinerecWeights );
    hListNameofParameterHistograms["hLinerecWeights"] = hLinerecWeights;
    
    // rate histograms
    
    sprintf( i_key, "hrate_1sec_%s", fHisSuffix.c_str() );
    hrate_1sec = new TH1D( i_key, "One Second Count rate", 100, i_emin, i_emax );
    hrate_1sec->SetXTitle( "Time [MJD]" );
    hrate_1sec->SetYTitle( "Rate (Hz)" );
    hisRateList->Add( hrate_1sec );
    
    sprintf( i_key, "hrate_10sec_%s", fHisSuffix.c_str() );
    hrate_10sec = new TH1D( i_key, "Ten Second Count rate", 50, i_emin, i_emax );
    hrate_10sec->SetXTitle( "Time [MJD]" );
    hrate_10sec->SetYTitle( "Rate (per 10 seconds)" );
    hisRateList->Add( hrate_10sec );
    
    sprintf( i_key, "hrate_1min_%s", fHisSuffix.c_str() );
    hrate_1min = new TH1D( i_key, "One Minute Count rate", 50, i_emin, i_emax );
    hrate_1min->SetXTitle( "Time [MJD]" );
    hrate_1min->SetYTitle( "Rate (per minute)" );
    hisRateList->Add( hrate_1min );
    
    sprintf( i_key, "hTriggerPatternBeforeCuts_%s", fHisSuffix.c_str() );
    hTriggerPatternBeforeCuts = new TH1D( i_key, "Trigger pattern before cuts", 16, 0., 16. );
    hTriggerPatternBeforeCuts->SetYTitle( "number of events" );
    hisRateList->Add( hTriggerPatternBeforeCuts );
    
    sprintf( i_key, "hTriggerPatternAfterCuts_%s", fHisSuffix.c_str() );
    hTriggerPatternAfterCuts = new TH1D( i_key, "Trigger pattern After cuts", 16, 0., 16. );
    hTriggerPatternAfterCuts->SetYTitle( "number of events" );
    hisRateList->Add( hTriggerPatternAfterCuts );
    
    sprintf( i_key, "hImagePatternBeforeCuts_%s", fHisSuffix.c_str() );
    hImagePatternBeforeCuts = new TH1D( i_key, "image pattern before cuts", 16, 0., 16. );
    hImagePatternBeforeCuts->SetYTitle( "number of events" );
    hisRateList->Add( hImagePatternBeforeCuts );
    
    sprintf( i_key, "hImagePatternAfterCuts_%s", fHisSuffix.c_str() );
    hImagePatternAfterCuts = new TH1D( i_key, "image pattern After cuts", 16, 0., 16. );
    hImagePatternAfterCuts->SetYTitle( "number of events" );
    hisRateList->Add( hImagePatternAfterCuts );
    
    /////
    //// ***************** IMPORTANT ***************************
    ////
    ////      the following histograms are not added up for the
    ////      combined analysis, please do not define any histogram
    ////      after them
    
    ////////////////////////
    // time binned histograms
    if( i_tbin > 0 && tbins.size() > 0 )
    {
        sprintf( i_key, "herecCounts2DtimeBinned_%s", fHisSuffix.c_str() );
        sprintf( i_name, "counting histogram (energy) (%s)", fHisSuffix.c_str() );
        herecCounts2DtimeBinned = new TH2D( i_key, i_name, i_ebin, i_emin, i_emax, i_tbin, &( tbins[0] ) );
        herecCounts2DtimeBinned->SetXTitle( "log_{10} energy [TeV]" );
        herecCounts2DtimeBinned->SetZTitle( "No. of events" );
        herecCounts2DtimeBinned->SetYTitle( "time [s]" );
        hisList->Add( herecCounts2DtimeBinned );
        hListParameterHistograms->Add( herecCounts2DtimeBinned );
        hListEnergyHistograms->Add( herecCounts2DtimeBinned );
        hListNameofParameterHistograms["herecCounts2DtimeBinned"] = herecCounts2DtimeBinned;
        
        sprintf( i_key, "hDuration1DtimeBinned_%s", fHisSuffix.c_str() );
        sprintf( i_name, "Bin Duration histogram (%s)", fHisSuffix.c_str() );
        hDuration1DtimeBinned = new TH1D( i_key, i_name, i_tbin, &( tbins[0] ) );
        hDuration1DtimeBinned->SetZTitle( "Real Duration [s]" );
        hDuration1DtimeBinned->SetXTitle( "time [s]" );
        hisList->Add( hDuration1DtimeBinned );
        hListParameterHistograms->Add( hDuration1DtimeBinned );
        hListEnergyHistograms->Add( hDuration1DtimeBinned );
        hListNameofParameterHistograms["hDuration1DtimeBinned"] = hDuration1DtimeBinned;
        
        sprintf( i_key, "hRealDuration1DtimeBinned_%s", fHisSuffix.c_str() );
        sprintf( i_name, "Real Duration histogram (%s)", fHisSuffix.c_str() );
        hRealDuration1DtimeBinned = new TH1D( i_key, i_name, i_tbin, &( tbins[0] ) );
        hRealDuration1DtimeBinned->SetZTitle( "Real Duration [s]" );
        hRealDuration1DtimeBinned->SetXTitle( "time [s]" );
        hisList->Add( hRealDuration1DtimeBinned );
        hListParameterHistograms->Add( hRealDuration1DtimeBinned );
        hListEnergyHistograms->Add( hRealDuration1DtimeBinned );
        hListNameofParameterHistograms["hRealDuration1DtimeBinned"] = hRealDuration1DtimeBinned;
        
        sprintf( i_key, "hLinerecCounts2DtimeBinned_%s", fHisSuffix.c_str() );
        sprintf( i_name, "counting histogram (energy) (%s)", fHisSuffix.c_str() );
        hLinerecCounts2DtimeBinned = new TH2D( i_key, i_name, i_Linebin, i_Linemin, i_Linemax, i_tbin, &( tbins[0] ) );
        hLinerecCounts2DtimeBinned->SetXTitle( "energy [TeV]" );
        hLinerecCounts2DtimeBinned->SetZTitle( "No. of events" );
        hLinerecCounts2DtimeBinned->SetYTitle( "time [s]" );
        hisList->Add( hLinerecCounts2DtimeBinned );
        hListParameterHistograms->Add( hLinerecCounts2DtimeBinned );
        hListEnergyHistograms->Add( hLinerecCounts2DtimeBinned );
        hListNameofParameterHistograms["hLinerecCounts2DtimeBinned"] = hLinerecCounts2DtimeBinned;
    }
    
    
    
    
    
    // set xtitles for trigger pattern histograms
    vector< string > i_xtitle;
    i_xtitle.push_back( "Tel.1" );
    i_xtitle.push_back( "Tel.2" );
    i_xtitle.push_back( "Tel.1+2" );
    i_xtitle.push_back( "Tel.3" );
    i_xtitle.push_back( "Tel.1+3" );
    i_xtitle.push_back( "Tel.2+3" );
    i_xtitle.push_back( "Tel.1+2+3" );
    i_xtitle.push_back( "Tel.4" );
    i_xtitle.push_back( "Tel.1+4" );
    i_xtitle.push_back( "Tel.2+4" );
    i_xtitle.push_back( "Tel.1+2+4" );
    i_xtitle.push_back( "Tel.3+4" );
    i_xtitle.push_back( "Tel.1+3+4" );
    i_xtitle.push_back( "Tel.2+3+4" );
    i_xtitle.push_back( "Tel.1+2+3+4" );
    for( unsigned int i = 0; i < i_xtitle.size(); i++ )
    {
        hTriggerPatternBeforeCuts->GetXaxis()->SetBinLabel( i + 2, i_xtitle[i].c_str() );
        hTriggerPatternAfterCuts->GetXaxis()->SetBinLabel( i + 2, i_xtitle[i].c_str() );
        hImagePatternBeforeCuts->GetXaxis()->SetBinLabel( i + 2, i_xtitle[i].c_str() );
        hImagePatternAfterCuts->GetXaxis()->SetBinLabel( i + 2, i_xtitle[i].c_str() );
    }
    
    string icname;
    TIter next( hisList );
    while( TH1* obj = ( TH1* )next() )
    {
        icname = obj->ClassName();
        if( icname != "TProfile" )
        {
            obj->Sumw2();
        }
        obj->SetLineWidth( 2 );
        if( !bIsOn )
        {
            obj->SetLineColor( 2 );
            obj->SetMarkerColor( 2 );
        }
        obj->GetYaxis()->SetTitleOffset( 1.2 );
    }
    TIter nextR( hisRateList );
    while( TH1* obj = ( TH1* )nextR() )
    {
        obj->Sumw2();
    }
}


void VStereoHistograms::scaleDistributions( double is )
{
    TIter next( hListParameterHistograms );
    while( TH1* obj = ( TH1* )next() )
    {
        obj->Scale( is );
    }
}


void VStereoHistograms::writeHistograms()
{
    TDirectory* iDir = gDirectory;
    
    TDirectory* wDir = 0;
    
    // write all sky plots into sky histogram directory
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "skyHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "skyHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    hListSkyMaps->Write();
    hListSkyMapsUC->Write();
    // delete sky plots from memory
    deleteSkyPlots();
    
    // write all stereo parameter histograms
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "stereoParameterHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "stereoParameterHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    hListStereoParameterHistograms->Write();
    
    // write all energy histograms
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "energyHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "energyHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    hListEnergyHistograms->Write();
    
    // write all random forest histograms
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "randomForestHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "randomForestHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    hListRandomForestParameterHistograms->Write();
    
    // write all rate histograms
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "rawRateHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "rawRateHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    hisRateList->Write();
    hisRateList->Delete();
    
    // delete all parameter histograms from heap
    deleteParameterHistograms();
    
    iDir->cd();
}


void VStereoHistograms::deleteSkyPlots()
{
    hListSkyMaps->Delete();
    hListSkyMapsUC->Delete();
}


void VStereoHistograms::deleteParameterHistograms()
{
    hListParameterHistograms->Delete();
}


bool VStereoHistograms::readParameterHistograms()
{
    TDirectory* iDir = gDirectory;
    
    if( !readHistograms( hListStereoParameterHistograms, "stereoParameterHistograms" ) )
    {
        return false;
    }
    iDir->cd();
    if( !readHistograms( hListRandomForestParameterHistograms, "randomForestHistograms" ) )
    {
        return false;
    }
    iDir->cd();
    if( !readHistograms( hListEnergyHistograms, "energyHistograms" ) )
    {
        return false;
    }
    iDir->cd();
    if( !readHistograms( hisRateList, "rawRateHistograms" ) )
    {
        return false;
    }
    iDir->cd();
    
    return true;
}


bool VStereoHistograms::readHistograms( TList* iL, string iDir )
{
    TDirectory* wDir = 0;
    
    // write all sky plots into sky histogram directory
    wDir = ( TDirectory* )gDirectory->Get( iDir.c_str() );
    if( !wDir )
    {
        string iT = "stereo/" + iDir;
        wDir = ( TDirectory* )gDirectory->Get( iT.c_str() );
    }
    if( !wDir )
    {
        cout << "VStereoHistograms::readHistograms() directory not found " << iDir << endl;
        return false;
    }
    wDir->cd();
    string iTemp;
    TIter next( wDir->GetListOfKeys() );
    map< string, TH1* >::iterator iter;
    
    while( TKey* key = ( TKey* )next() )
    {
        TNamed* ho = ( TNamed* )key->ReadObj();
        iTemp = ho->GetName();
        
        // get only on/off histograms
        if( iTemp.find( "_diff" ) < iTemp.size() )
        {
            ho->Delete();
            continue;
        }
        
        bool bRead = false;
        if( bIsOn && iTemp.find( "_on" ) < iTemp.size() )
        {
            bRead = true;
        }
        if( !bIsOn && iTemp.find( "_off" ) < iTemp.size() )
        {
            bRead = true;
        }
        
        if( bRead )
        {
            for( iter = hListNameofParameterHistograms.begin(); iter != hListNameofParameterHistograms.end(); iter++ )
            {
                if( iTemp.find( iter->first ) < iTemp.size() )
                {
                    iter->second = ( TH1* )ho;
                }
            }
            hListParameterHistograms->Add( ho );
            if( iL )
            {
                iL->Add( ho );
            }
            
        }
        else
        {
            ho->Delete();
        }
    }
    return true;
}

/*
 * access all run directories and read sky maps from disk
 *
 */
bool VStereoHistograms::readSkyPlots()
{
    TDirectory* wDir = 0;
    
    // sky histogram directory
    wDir = ( TDirectory* )gDirectory->Get( "skyHistograms" );
    if( !wDir )
    {
        wDir = ( TDirectory* )gDirectory->Get( "stereo/skyHistograms" );
    }
    if( !wDir )
    {
        cout << "VStereoHistograms::readSkyPlots() directory not found " << endl;
        return false;
    }
    wDir->cd();
    
    h_combine_map_alpha_off = 0;
    h_combine_map_alpha_offUC = 0;
    h_combine_map_stereo_on = 0;
    h_combine_map_stereo_onUC = 0;
    h_combine_map_stereo_off = 0;
    h_combine_map_stereo_offUC = 0;
    
    //////////////////////////////////////////
    // loop over all histograms in directory
    string iTemp;
    TIter next( wDir->GetListOfKeys() );
    while( TKey* key = ( TKey* )next() )
    {
        iTemp = key->GetName();
        // skip difference and significance histograms
        if( iTemp.find( "_diff" ) < iTemp.size() || iTemp.find( "_sig" ) < iTemp.size() )
        {
            continue;
        }
        
        bool bFill = false;
        if( bIsOn && iTemp.find( "_on" ) < iTemp.size() )
        {
            bFill = true;
        }
        if( !bIsOn && iTemp.find( "_off" ) < iTemp.size() )
        {
            bFill = true;
        }
        
        // read object from disk
        TNamed* ho = ( TNamed* )key->ReadObj();
        
        // histogram needed for VStereoAnalysis::combineHistograms()
        if( iTemp.find( "hmap_alphaNorm_off" ) < iTemp.size() )
        {
            h_combine_map_alpha_off = ( TH2D* )ho;
        }
        if( iTemp.find( "hmap_stereo_on" ) < iTemp.size() )
        {
            h_combine_map_stereo_on = ( TH2D* )ho;
        }
        if( iTemp.find( "hmap_stereo_off" ) < iTemp.size() )
        {
            h_combine_map_stereo_off = ( TH2D* )ho;
        }
        if( iTemp.find( "hmap_alphaNormUC_off" ) < iTemp.size() )
        {
            h_combine_map_alpha_offUC = ( TH2D* )ho;
        }
        if( iTemp.find( "hmap_stereoUC_on" ) < iTemp.size() )
        {
            h_combine_map_stereo_onUC = ( TH2D* )ho;
        }
        if( iTemp.find( "hmap_stereoUC_off" ) < iTemp.size() )
        {
            h_combine_map_stereo_offUC = ( TH2D* )ho;
        }
        
        if( !bFill )
        {
            continue;
        }
        
        /////////////////////////////
        // all other sky histograms
        
        // uncorrelated histograms
        if( iTemp.find( "UC" ) < iTemp.size() )
        {
            hListSkyMapsUC->Add( ho );
            if( iTemp.find( "hmap_stereoUC" ) < iTemp.size() )
            {
                hmap_stereoUC = ( TH2D* )ho;
            }
            if( iTemp.find( "hmap_alphaNormUC" ) < iTemp.size() )
            {
                hmap_alphaNormUC = ( TH2D* )ho;
            }
            if( iTemp.find( "hmap_alphaUC" ) < iTemp.size() )
            {
                hmap_alphaUC = ( TH2D* )ho;
            }
        }
        // correlated histograms
        else
        {
            hListSkyMaps->Add( ho );
            if( iTemp.find( "hmap_stereo" ) < iTemp.size() )
            {
                hmap_stereo = ( TH2D* )ho;
            }
            else if( iTemp.find( "hmap_alphaNorm" ) < iTemp.size() )
            {
                hmap_alphaNorm = ( TH2D* )ho;
            }
            else if( iTemp.find( "hmap_alpha" ) < iTemp.size() )
            {
                hmap_alpha = ( TH2D* )ho;
            }
        }
    }
    return true;
}


void VStereoHistograms::makeRateHistograms( double iStart, double iStopp )
{
    if( !hisRateList )
    {
        return;
    }
    
    if( !hrate_1sec || !hrate_10sec || !hrate_1min )
    {
        return;
    }
    
    // reset bins for rate histograms
    hrate_1sec->Reset();
    hrate_10sec->Reset();
    hrate_1min->Reset();
    
    // 1 bin per second
    int i_nbins = ( int )( ( iStopp - iStart ) * 24. * 60. * 60. );
    hrate_1sec->SetBins( i_nbins, iStart, iStopp );
    i_nbins = ( int )( ( iStopp - iStart ) / 10. * 24. * 60. * 60. );
    hrate_10sec->SetBins( i_nbins, iStart, iStopp );
    i_nbins = ( int )( ( iStopp - iStart ) / 60. * 24. * 60. * 60. );
    hrate_1min->SetBins( i_nbins, iStart, iStopp );
}


/*!

    copying objects over to anasum file (e.g. from effective area file or radial acceptance file)

*/
void VStereoHistograms::writeObjects( string iFile, string iDirectory, TObject* g )
{
    TDirectory* iDir = gDirectory;
    
    TDirectory* wDir = 0;
    
    iDir->cd();
    if( iDirectory.size() > 0 )
    {
        wDir = ( TDirectory* )iDir->Get( iDirectory.c_str() );
        if( !wDir )
        {
            iDir->mkdir( iDirectory.c_str() )->cd();
        }
        else
        {
            wDir->cd();
        }
    }
    wDir = gDirectory;
    
    if( iFile.find( "IGNOREEFFECTIVEAREA" ) == string::npos )
    {
        iFile = VUtilities::testFileLocation( iFile, iDirectory, true );
    }
    if( iFile.size() == 0 )
    {
        iDir->cd();
        return;
    }
    // open input file
    TFile fIn( iFile.c_str(), "READ" );
    if( fIn.IsZombie() )
    {
        iDir->cd();
        return;
    }
    if( iDirectory == "EffectiveAreas" )
    {
        if( g )
        {
            wDir->cd();
            g->Write();
        }
        /////////////////////////////////////////////////
    }
    else if( iDirectory == "RadialAcceptances" )
    {
        string itemp;
        
        TIter next( fIn.GetListOfKeys() );
        while( TKey* key = ( TKey* )next() )
        {
            TNamed* ho = ( TNamed* )key->ReadObj();
            if( ho )
            {
                itemp = ho->GetName();
                if( itemp.find( "AccZe" ) < itemp.size() ||
                        itemp.find( "RadialAcceptance" ) < itemp.size() )
                {
                    wDir->cd();
                    ho->Write();
                    fIn.cd();
                }
            }
        }
    }
    else
    {
        TIter next( fIn.GetListOfKeys() );
        while( TKey* key = ( TKey* )next() )
        {
            TNamed* ho = ( TNamed* )key->ReadObj();
            if( ho )
            {
                wDir->cd();
                ho->Write();
                fIn.cd();
            }
        }
    }
    fIn.Close();
    iDir->cd();
}
