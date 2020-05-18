/*! \file VTableCalculator.cpp
    \brief  general lookup table class (e.g. mscw, mscl, energy)


*/

#include "VTableCalculator.h"

VTableCalculator::VTableCalculator( int intel, bool iEnergy, bool iPE )
{
    setDebug();
    
    // initialize variables for table value normalization
    setNormalizeTableValues();
    
    setConstants( iPE );
    
    if( intel == 0 )
    {
        return;
    }
    
    fEnergy = iEnergy;
    fUseMedianEnergy = 0;
    
    setEventSelectionCut();
    
    for( int i = 0; i < intel; i++ )
    {
        hVMedian.push_back( 0 );
    }
    hMedian = 0;
    hMean = 0;
    
    fWriteTables = false;
    
    fInterPolWidth = 1;
    fInterPolIter = 3;
    
    // lookup table 1D histogram bin limits
    // generally [0,1]
    fBinning1DXlow = 1.e-5;
    fBinning1DXhigh = 1. + 1.e-5;
    // bin limits in log10(GeV) for energy calculation
    if( fEnergy )
    {
        fBinning1DXlow =  1.;
        fBinning1DXhigh = 5.;
    }
    
    fReadHistogramsFromFile = false;
    
}


VTableCalculator::VTableCalculator( string fpara, string hname_add,
                                    bool iWriteTables, TDirectory* iDir,
                                    bool iEnergy, bool iPE, int iUseMedianEnergy )
{
    setDebug();
    
    // initialize variables for table value normalization
    setNormalizeTableValues();
    
    // use 1D histograms to calculate medians (more precise, but needs much more memory and is slower)
    fWrite1DHistograms = true;
    // use median approximation (faster)
    fFillMedianApproximations = false;
    
    setConstants( iPE );
    // using lookup tables to calculate energies
    fEnergy = iEnergy;
    fUseMedianEnergy = iUseMedianEnergy;
    fReadHistogramsFromFile = false;
    
    setEventSelectionCut();
    
    fHName_Add = hname_add;
    
    fName = fpara;
    
    fInterPolWidth = 1;
    fInterPolIter = 3;
    
    setBinning();
    
    fOutDir = iDir;
    if( !fOutDir )
    {
        cout << "VTableCalculator: error data directory in root file does not exist " << fOutDir << "\t" << fpara << endl;
        exit( EXIT_FAILURE );
    }
    if( !fOutDir->cd() )
    {
        cout << "VTableCalculator: error accessing data directory in root file " << fOutDir << "\t" << fpara << endl;
        exit( EXIT_FAILURE );
    }
    fWriteTables = iWriteTables;
    
    int i = 0;
    int j = 0;
    char hname[1000];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // table writing
    if( fWriteTables )
    {
        if( !fOutDir->IsWritable() )
        {
            cout << "VTableCalculator: error data directory in root file not writable" << endl;
            exit( EXIT_FAILURE );
        }
        
        /* HSTOGRAM BOOKING */
        char htitle[1000];
        // median of variable
        sprintf( hname, "%s_median_%s", fpara.c_str(), fHName_Add.c_str() );
        sprintf( htitle, "%s vs. dist. vs. log10 size (median)", fpara.c_str() );
        hMedian = new TH2F( hname, htitle, NumSize, amp_offset, amp_offset + NumSize * amp_delta, NumDist, 0., dist_delta * NumDist );
        hMedian->SetXTitle( "log_{10} size" );
        hMedian->SetYTitle( "distance [m]" );
        if( !fEnergy )
        {
            sprintf( htitle, "%s (median) [deg]", fpara.c_str() );
        }
        else
        {
            sprintf( htitle, "%s (median) [TeV]", fpara.c_str() );
        }
        hMedian->SetZTitle( htitle );
        // mean and rms
        sprintf( hname, "%s_mean_%s", fpara.c_str(), fHName_Add.c_str() );
        sprintf( htitle, "%s vs. dist. vs. log10 size (mean)", fpara.c_str() );
        hMean = new TProfile2D( hname, htitle,
                                NumSize, amp_offset, amp_offset + NumSize * amp_delta,
                                NumDist, 0., dist_delta * NumDist,
                                fBinning1DXlow, fBinning1DXhigh );
        hMean->SetXTitle( "log_{10} size" );
        hMean->SetYTitle( "distance [m]" );
        if( !fEnergy )
        {
            sprintf( htitle, "%s (mean) [deg]", fpara.c_str() );
        }
        else
        {
            sprintf( htitle, "%s (mean) [TeV]", fpara.c_str() );
        }
        hMean->SetZTitle( htitle );
        // 1d histograms for variable distribution
        for( i = 0; i < NumSize; i++ )
        {
            vector< TH1F* > iH1;
            vector< VMedianCalculator* > iM1;
            for( j = 0; j < NumDist; j++ )
            {
                iH1.push_back( 0 );
                iM1.push_back( 0 );
            }
            Oh.push_back( iH1 );
            OMedian.push_back( iM1 );
        }
    }
    /////////////////////////////////////////
    // table reading
    /////////////////////////////////////////
    else
    {
        fReadHistogramsFromFile = false;
        
        if( fUseMedianEnergy == 1 )
        {
            sprintf( hname, "%s_median_%s", fpara.c_str(), fHName_Add.c_str() );
        }
        else if( fUseMedianEnergy == 2 )
        {
            if( fEnergy )
            {
                sprintf( hname, "%s_mpv_%s", fpara.c_str(), fHName_Add.c_str() );
            }
            else
            {
                sprintf( hname, "%s_median_%s", fpara.c_str(), fHName_Add.c_str() );
            }
        }
        else
        {
            if( fEnergy )
            {
                sprintf( hname, "%s_mean_%s", fpara.c_str(), fHName_Add.c_str() );
            }
            else
            {
                sprintf( hname, "%s_median_%s", fpara.c_str(), fHName_Add.c_str() );
            }
        }
        hMedianName = hname;
    }
    
}

/*
 * bin size settings for 1D histograms
 *
 * adaptive binning for energy calculation
 * all other values always between [0,1]
 *
 */

void VTableCalculator::setBinning()
{
    fBinning1DXlow = 0.;
    fBinning1DXhigh = 1.;
    if( fEnergy )
    {
        // lowest energy bin: 5 GeV
        fBinning1DXlow =  0.005;
        // highest energy bin: 300 TeV
        fBinning1DXhigh = 300.;
        HistBins = int( fBinning1DXhigh / 0.005 );
        ///////////////////////////////////////////////////////////////////////////////////
        // set adaptive binning to make sure target on energy resolution is always met
        // unit is TeV
        fBinning1Dxbins = new float[4000];
        // 1 - 100 GeV: bins of 1 GeV
        fBinning1DxbinsN = 0;
        for( unsigned int i = 0; i < 95 ; i++ )
        {
            fBinning1Dxbins[fBinning1DxbinsN] = fBinning1DXlow + i * 1.e-3;
            fBinning1DxbinsN++;
        }
        // 100 GeV - 1 TeV: bins of 0.005 (5 GeV)
        for( unsigned int i = 0; i < 180; i++ )
        {
            fBinning1Dxbins[fBinning1DxbinsN] = 0.1 + i * 5.e-3;
            fBinning1DxbinsN++;
        }
        // 1 TeV to 10 TeV: bins of 0.025 (25 GeV)
        for( unsigned int i = 0; i < 360; i++ )
        {
            fBinning1Dxbins[fBinning1DxbinsN] = 1. + i * 25.e-3;
            fBinning1DxbinsN++;
        }
        // 10 TeV to 300 TeV: bins of 0.1 (100 GeV)
        for( unsigned int i = 0; i < 2900; i++ )
        {
            fBinning1Dxbins[fBinning1DxbinsN] = 10. + i * 100.e-3;
            fBinning1DxbinsN++;
        }
        // counting needs to be n-1
        fBinning1DxbinsN = fBinning1DxbinsN - 1;
    }
    else if( HistBins > 0 )
    {
        float i_width = ( fBinning1DXhigh - fBinning1DXlow ) / ( float )HistBins;
        fBinning1Dxbins = new float[HistBins + 1];
        for( int i = 0; i <= HistBins; i++ )
        {
            fBinning1Dxbins[i] = fBinning1DXlow + i * i_width;
        }
        fBinning1DxbinsN = HistBins;
    }
    
}

bool VTableCalculator::createMedianApprox( int i, int j )
{
    if( i >= 0 && j >= 0 && i < ( int )OMedian.size() && j < ( int )OMedian[i].size() && !OMedian[i][j] )
    {
        OMedian[i][j] = new VMedianCalculator();
    }
    else
    {
        return false;
    }
    return true;
}

bool VTableCalculator::create1DHistogram( int i, int j, double w_first_event )
{
    if( i >= 0 && j >= 0 && i < ( int )Oh.size() && j < ( int )Oh[i].size() && !Oh[i][j] )
    {
        if( !fOutDir->cd() )
        {
            return false;
        }
        char hisname[200];
        char histitle[200];
        int id = i * 1000 + j;
        
        sprintf( hisname , "h%d", id );
        double is1 = hMedian->GetXaxis()->GetBinLowEdge( i + 1 );
        double is2 = hMedian->GetXaxis()->GetBinLowEdge( i + 1 ) + hMedian->GetXaxis()->GetBinWidth( i + 1 );
        double id1 = hMedian->GetYaxis()->GetBinLowEdge( j + 1 );
        double id2 = hMedian->GetYaxis()->GetBinLowEdge( j + 1 ) + hMedian->GetYaxis()->GetBinWidth( j + 1 );
        sprintf( histitle, "%.2f < log10 size < %.2f, %.1f < r < %.1f (%s)", is1, is2, id1, id2, fHName_Add.c_str() );
        
        Oh[i][j] = new TH1F( hisname, histitle, fBinning1DxbinsN, fBinning1Dxbins );
        Oh[i][j]->SetXTitle( fName.c_str() );
        // allow automatic rebinning
#ifdef ROOT6
        Oh[i][j]->GetXaxis()->SetCanExtend( true );
#else
        Oh[i][j]->SetBit( TH1::kCanRebin );
#endif
    }
    else
    {
        return false;
    }
    return true;
}


void VTableCalculator::setConstants( bool iPE )
{
    NumSize = 55;
    amp_offset = 1.5;
    amp_delta = 0.1;
    NumDist = 80;
    dist_delta = 15.;
    HistBins = 500;
    xlow = 0.;
    xhigh = 1.;
    fMinShowerPerBin = 10;
    
    // binning is different for MC with values in PE
    if( iPE )
    {
        NumSize = 45;
        amp_offset = 1.0;
        amp_delta = 0.15;
        
        NumDist = 75;
        dist_delta = 20.;
        HistBins = 500;
    }
}


void VTableCalculator::terminate( TDirectory* iOut, char* xtitle )
{
    if( iOut != 0 )
    {
        fOutDir = iOut;
    }
    
    if( !fOutDir->cd() )
    {
        cout << "Error: unable to reach writing directory ( VTableCalculator::terminate())" << endl;
        return;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // table writing
    /////////////////////////////////////////////////////////////////////////////////////////////
    if( fWriteTables )
    {
        TDirectory* iDir1D = 0;
        // make output directory for 1D histograms
        if( fOutDir )
        {
            iDir1D = fOutDir->mkdir( "histos1D" );
        }
        
        ///////////////////////////////////
        // 2D histograms
        // number of events
        char hname[1000];
        char htitle[1000];
        sprintf( hname, "%s_nevents_%s", fName.c_str(), fHName_Add.c_str() );
        sprintf( htitle, "%s vs. dist. vs. log10 size (# of events)", fName.c_str() );
        TH2F* hNevents = new TH2F( hname, htitle, NumSize, amp_offset, amp_offset + NumSize * amp_delta, NumDist, 0., dist_delta * NumDist );
        hNevents->SetXTitle( "log_{10} size" );
        hNevents->SetYTitle( "distance [m]" );
        hNevents->SetZTitle( "# of events/bin" );
        // most probable of variable
        sprintf( hname, "%s_mpv_%s", fName.c_str(), fHName_Add.c_str() );
        sprintf( htitle, "%s vs. dist. vs. log10 size (mpv)", fName.c_str() );
        TH2F* hMPV = new TH2F( hname, htitle, NumSize, amp_offset, amp_offset + NumSize * amp_delta, NumDist, 0., dist_delta * NumDist );
        hMPV->SetXTitle( "log_{10} size" );
        hMPV->SetYTitle( "distance [m]" );
        if( !fEnergy )
        {
            sprintf( htitle, "%s (mpv) [deg]", fName.c_str() );
        }
        else
        {
            sprintf( htitle, "%s (mpv) [TeV]", fName.c_str() );
        }
        hMPV->SetZTitle( htitle );
        // sigma of median (16-84% (2sigma for Gauss))
        sprintf( hname, "%s_sigma_%s", fName.c_str(), fHName_Add.c_str() );
        sprintf( htitle, "%s vs. dist. vs. log10 size (sigma)", fName.c_str() );
        TH2F* hSigma = new TH2F( hname, htitle, NumSize, amp_offset, amp_offset + NumSize * amp_delta, NumDist, 0., dist_delta * NumDist );
        hSigma->SetXTitle( "log_{10} size" );
        hSigma->SetYTitle( "distance [m]" );
        if( !fEnergy )
        {
            sprintf( htitle, "%s (2xsigma) [deg]", fName.c_str() );
        }
        else
        {
            sprintf( htitle, "%s (2xsigma) [TeV]", fName.c_str() );
        }
        hSigma->SetZTitle( htitle );
        
        ///////////////////////////////////
        // EVALUATION OF HISTOGRAMS
        if( hNevents->GetEntries() > 0 )
        {
            cout << "\t tables: evaluating " << fName << " histograms ";
        }
        
        float med = 0.;
        float sigma = 0.;
        int   nevents = 0;
        
        double i_a[] = { 0.16, 0.5, 0.84 };
        double i_b[] = { 0.0,  0.0, 0.0  };
        
        // loop over all size bin and distance bins
        for( int i = 0; i < NumSize; i++ )
        {
            for( int j = 0; j < NumDist; j++ )
            {
                // use 1D histograms for median calculation
                if( fWrite1DHistograms && Oh[i][j] &&  Oh[i][j]->GetEntries() > fMinShowerPerBin )
                {
                    Oh[i][j]->GetQuantiles( 3, i_b, i_a );
                    med     = i_b[1];
                    sigma   = i_b[2] - i_b[0];
                    nevents = Oh[i][j]->GetEntries();
                }
                // use approx median calculation
                else if( fFillMedianApproximations && OMedian[i][j] && OMedian[i][j]->getN() > fMinShowerPerBin )
                {
                    med     = OMedian[i][j]->getMedian( sigma, nevents );
                }
                else
                {
                    med = 0.;
                    sigma = 0.;
                    nevents = 0;
                }
                hMedian->SetBinContent( i + 1, j + 1, med );
                hMedian->SetBinError( i + 1, j + 1, sigma );
                hSigma->SetBinContent( i + 1, j + 1, sigma );
                if( nevents > fMinShowerPerBin )
                {
                    hNevents->SetBinContent( i + 1, j + 1, nevents );
                }
                else if( fEnergy )
                {
                    hMean->SetBinContent( i + 1, j + 1, 0. );
                }
                if( fEnergy && fWrite1DHistograms )
                {
                    fillMPV( hMPV, i + 1, j + 1, Oh[i][j], med, sigma );
                    hMPV->SetBinError( i + 1, j + 1, sigma );
                }
                // write 1D histograms to file
                if( fWrite1DHistograms && Oh[i][j] )
                {
                    if( fOutDir && Oh[i][j]->GetEntries() > 0 )
                    {
                        fOutDir->cd();
                        iDir1D->cd();
                        Oh[i][j]->Write();
                    }
                    delete Oh[i][j];
                }
                else if( fFillMedianApproximations && OMedian[i][j] )
                {
                    delete OMedian[i][j];
                }
            }
        }
        // write 2D histograms to file
        if( fOutDir && hNevents->GetEntries() > 0 )
        {
            fOutDir->cd();
            if( xtitle && hMedian )
            {
                hMedian->SetTitle( xtitle );
            }
            if( xtitle && hMPV )
            {
                hMPV->SetTitle( xtitle );
            }
            if( hNevents && hMedian )
            {
                hMedian->SetEntries( hNevents->GetEntries() );
            }
            if( hNevents && hSigma )
            {
                hSigma->SetEntries( hNevents->GetEntries() );
            }
            // reduce size of the 2D histograms
            TH2F* h = 0;
            string n;
            if( hMedian )
            {
                n = hMedian->GetName();
                h = VHistogramUtilities::reduce2DHistogramSize( hMedian, n + "_new" );
                if( h )
                {
                    cout << "(" << hMedian->GetEntries() << " entries)";
                    delete hMedian;
                    if( h )
                    {
                        h->SetName( n.c_str() );
                        h->Write();
                    }
                    delete h;
                }
            }
            if( hMPV )
            {
                n = hMPV->GetName();
                h = VHistogramUtilities::reduce2DHistogramSize( hMPV, n + "_new" );
                if( h )
                {
                    delete hMPV;
                    if( h )
                    {
                        h->SetName( n.c_str() );
                        h->Write();
                    }
                    delete h;
                }
            }
            if( hSigma )
            {
                n = hSigma->GetName();
                h = VHistogramUtilities::reduce2DHistogramSize( hSigma, n + "_new" );
                if( h )
                {
                    delete hSigma;
                    if( h )
                    {
                        h->SetName( n.c_str() );
                        h->Write();
                    }
                    delete h;
                }
            }
            if( hMean )
            {
                n = hMean->GetName();
                h = VHistogramUtilities::reduce2DHistogramSize( hMean, n + "_new" );
                if( h )
                {
                    delete hMean;
                    if( h )
                    {
                        h->SetName( n.c_str() );
                        h->Write();
                    }
                    delete h;
                }
            }
            if( hNevents )
            {
                n = hNevents->GetName();
                h = VHistogramUtilities::reduce2DHistogramSize( hNevents, n + "_new" );
                if( h )
                {
                    delete hNevents;
                    if( h )
                    {
                        h->SetName( n.c_str() );
                        h->Write();
                    }
                    delete h;
                }
            }
            cout << endl;
        }
    }
}


/*!
     main calculation routine for lookup tables

     used for table filling and reading

     \param ntel   number of telescopes
     \param r      core distance from each telescopes [ntel]
     \param s      size per telescopes [ntel]
     \param l      loss per telescope
     \param chi2   weight used while filling the lookup table
                   (e.g. a spectral weighting)

     Table filling:
     --------------

     \param chi2   external weight (table filling only)
     \param dE     not used
     \param s_sigma  no used

     Table reading:
     --------------

     Width/length calculation:

     \param w       width/length per telescope [ntel]
     \param mt      expected width/length per telescope [ntel]
     \param chi2    always 0 (meaningless)
     \param dE      always 0 (meaningless)
     \param s_sigma expected sigma for width/length per telescope [ntel]

     return value is mean scaled width/length

     Energy calculation:

     \param w       MCenergy [ntel] (is telescope independent, but easier to implement like that)
     \param mt      lookup table energy per telescope [ntel]
     \param chi2    scatter of energies per telescope
     \param dE      (mean) energy resolution

     return value is energy (linear scale)

*/
double VTableCalculator::calc( int ntel, double* r, double* s, double* l, double *d, double* w,
                               double* mt, double& chi2, double& dE, double* s_sigma )
{
    int tel = 0;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // normalize input parameters values to interval [0,1]
    // (only when this is set via setNormalizeTableValues()
    // should not be used for energy calculation
    // (need to copy the value, otherwise the normalized value will be filled
    //  into the output tree)
    vector< double > w_fill( ntel, 0. );
    for( int i = 0; i < ntel; i++ )
    {
        if( w )
        {
            if( fValueNormalizationRange_min > -9990. && fValueNormalizationRange_max > -9990.
                    && ( fValueNormalizationRange_max - fValueNormalizationRange_min ) != 0. )
            {
                w_fill[i] = ( w[i] - fValueNormalizationRange_min )
                            / ( fValueNormalizationRange_max - fValueNormalizationRange_min );
            }
            else
            {
                w_fill[i] = w[i];
            }
        }
        // no input values given
        else
        {
            w_fill[i] = -9999.;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // fill the tables
    //
    // ntel and size/width etc arrays are expected to be values from same telescope types
    // therefore: ntel = number of telescopes of same type
    // (fill one mscw table for each telescope type)
    ///////////////////////////////////////////////////////////////////////////////////////
    if( fWriteTables )
    {
        // don't allow zero or negative weights
        if( chi2 <= 0. )
        {
            return -99.;
        }
        
        // loop over all telescopes
        double i_logs = 0.;
        int ir = 0;
        int is = 0.;
        int i_Oh_size = 0;
        if( fWrite1DHistograms )
        {
            i_Oh_size = ( int )Oh.size();
        }
        else
        {
            i_Oh_size = ( int )OMedian.size();
        }
        // loop over all telescopes
        for( tel = 0; tel < ntel; tel++ )
        {
            if( s[tel] > 0. && r[tel] >= 0. 
                  && l[tel] < fEventSelectionCut_lossCutMax
                  && d[tel] < fEventSelectionCut_distanceCutMax
                  && w_fill[tel] > fBinning1DXlow && w_fill[tel] < fBinning1DXhigh )
            {
                // check limits (to avoid under/overflows)
                ir = hMedian->GetYaxis()->FindFixBin( r[tel] ) - 1;
                if( ir == NumDist )
                {
                    continue;
                }
                i_logs = log10( s[tel] );
                is = hMedian->GetXaxis()->FindFixBin( i_logs ) - 1;
                if( is == NumSize )
                {
                    continue;
                }
                // reject showers in the first size bin
                if( ir >= 0 && is >= 0 && is < i_Oh_size )
                {
                    if( fWrite1DHistograms && ir < ( int )Oh[is].size() )
                    {
                        if( !Oh[is][ir] && !create1DHistogram( is, ir, w_fill[tel] ) )
                        {
                            continue;
                        }
                        // fill width/length/energy into a 1D and 2D histogram
                        // (chi2 is here an external weight (from e.g. spectral weighting))
                        Oh[is][ir]->Fill( w_fill[tel], chi2 );
                        if( fEnergy )
                        {
                            if( w_fill[tel] < Oh[is][ir]->GetXaxis()->GetXmin() || w_fill[tel] > Oh[is][ir]->GetXaxis()->GetXmax() )
                            {
                                cout << "Energy table filling: value (" << w_fill[tel] << ") out of range:";
                                cout << Oh[is][ir]->GetXaxis()->GetXmin() << "\t" << Oh[is][ir]->GetXaxis()->GetXmax() << endl;
                            }
                        }
                    }
                    if( fFillMedianApproximations && ir < ( int )OMedian.size() )
                    {
                        if( !OMedian[is][ir] )
                        {
                            // weight is ignored in the approx median calculation
                            if( !createMedianApprox( is, ir ) )
                            {
                                continue;
                            }
                        }
                        OMedian[is][ir]->fill( w_fill[tel] );
                    }
                    hMean->Fill( i_logs, r[tel], w_fill[tel], chi2 );
                }
            }
        }
        return -99.;
    }
    /////////////////////////////////////////////////////////
    // END OF writing/filling lookup tables
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    // read table
    // compute mean or mean scaled width/length value
    //
    else
    {
    
        // tables are accessed for the first time: get the from the file
        if( !fReadHistogramsFromFile )
        {
            cout << "read tables from " << fOutDir->GetPath() << endl;
            hMedian = ( TH2F* )fOutDir->Get( hMedianName.c_str() );
            fReadHistogramsFromFile = true;
            
            if( !hMedian )
            {
                cout << "VTableCalculator error: table histograms not found in " << gDirectory->GetName() << endl;
                exit( EXIT_FAILURE );
            }
        }
        
        chi2 = 0.;
        dE   = 0.;
        double med = 0.;
        double sigma = 0.;
        double value = 0.;
        double weight = 0.;
        // energy per telescope
        vector< double > energy_tel;
        vector< double > sigma2_tel;
        vector< double > sigma2_tel_noRadiusWeigth;
        vector< double > sigma_tel;
        vector< bool > good_image;
        
        // reset everything
        for( tel = 0; tel < ntel; tel++ )
        {
            mt[tel] = -99.;
            if( s_sigma )
            {
                s_sigma[tel] = -99.;
            }
        }
        
        ////////////////////////////////////////////////////
        // loop over all telescopes
        ////////////////////////////////////////////////////
        for( tel = 0; tel < ntel; tel++ )
        {
            if( r[tel] >= 0. && s[tel] > 0 
            && l[tel] < fEventSelectionCut_lossCutMax 
            && d[tel] < fEventSelectionCut_distanceCutMax )
            {
                // get expected value and sigma of expected value
                if( hMedian )
                {
                    med   = interpolate( hMedian, log10( s[tel] ), r[tel], false );
                    sigma = interpolate( hMedian, log10( s[tel] ), r[tel], true );
                }
                else if( hVMedian.size() == ( unsigned int )ntel && hVMedian[tel] )
                {
                    med   = interpolate( hVMedian[tel], log10( s[tel] ), r[tel], false );
                    sigma = interpolate( hVMedian[tel], log10( s[tel] ), r[tel], true );
                    if( fDebug && fEnergy )
                    {
                        cout << "\t  double VTableCalculator::calc() getting energy from table for tel " << tel;
                        cout << ", size " << s[tel];
                        cout << " , distance " << r[tel];
                        cout << " :" << med << "\t" << sigma << endl;
                        cout << endl;
                    }
                }
                else
                {
                    med = 0.;
                    sigma = 0.;
                    if( fDebug )
                    {
                        cout << "\t  double VTableCalculator::calc() med equal zero for tel " << tel;
                        cout << ", size " << s[tel] << ", distance " << r[tel];
                        cout << " (no lookup table)" << endl;
                    }
                }
                // accept only values > 0
                // (expeced width/length should be > 0)
                // (log10 energy should be > 0, good reason why we work in GeV here)
                if( med > 0. )
                {
                    mt[tel] = med;
                    if( s_sigma )
                    {
                        s_sigma[tel] = sigma;
                    }
                }
                else
                {
                    mt[tel] = -99.;
                    if( s_sigma )
                    {
                        s_sigma[tel] = -99.;
                    }
                    if( fDebug )
                    {
                        cout << "\t  double VTableCalculator::calc() med equal zero from tables for tel " << tel;
                        cout << ", size " << s[tel] << ", distance " << r[tel] << endl;
                    }
                }
                // weighted mean
                if( med > 0. )
                {
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////
                    // mean scaled calculation
                    if( !fEnergy && sigma > 0. )
                    {
                        // handle showers with (width==0.) correctly
                        if( w_fill[tel] > 0. 
                           && l[tel] < fEventSelectionCut_lossCutMax 
                           && d[tel] < fEventSelectionCut_distanceCutMax )
                        {
                            value  += ( w_fill[tel] - med ) / sigma * ( med * med ) / ( sigma * sigma );
                            weight += ( med * med ) / ( sigma * sigma );
                        }
                    }
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////
                    // energy calculation
                    else if( fEnergy && sigma > 0. )
                    {
                        // store energy per telescope
                        energy_tel.push_back( med );
                        // store expected width of the distribution
                        sigma_tel.push_back( 1. / sigma );
                        // use relative error as weighting (otherwise: significant bias towards lower energies)
                        sigma2_tel.push_back( med / ( sigma * sigma ) );
                        sigma2_tel_noRadiusWeigth.push_back( 1. / ( sigma * sigma ) );
                        float scaleDistance = 400.; // default
                        // add addional weight for events inside or outside the light pool
                        if( r[tel] < scaleDistance )
                        {
                            sigma2_tel.back() = sigma2_tel.back() * 100.;
                        }
                        else
                        {
                            sigma2_tel.back() = sigma2_tel.back() * 100.*exp( -1.*( r[tel] - scaleDistance ) / 200. );
                        }
                        if( l[tel] < fEventSelectionCut_lossCutMax 
                           && d[tel] < fEventSelectionCut_distanceCutMax )
                        {
                            good_image.push_back( true );
                        }
                        else
                        {
                            good_image.push_back( false );
                        }
                    }
                    else
                    {
                        chi2 = -99;
                        dE = -99.;
                        return -99.;
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////////
        // mean scaled value
        // (MSCW/MSCL)
        if( !fEnergy )
        {
            if( weight > 0 )
            {
                return value / weight;
            }
            else
            {
                return -99.;
            }
        }
        ///////////////////////////////////////////////////////////////
        // Energy calculation only
        ///////////////////////////////////////////////////////////////
        // calculate mean energy
        if( energy_tel.size() > 0 && energy_tel.size() == good_image.size() )
        {
            // Occasionally one energy is significantly off and distorts the mean.
            // therefore: get rid of N sigma outliers
            // use robust statistics (median and mean absolute error)
            // Note: applied only to larger events > 4 telescopes
            double median = TMath::Median( energy_tel.size(), &energy_tel[0] );
            double meanAbsoluteError = VStatistics::getMeanAbsoluteError( energy_tel );
            weight = 0.;
            for( unsigned int j = 0; j < energy_tel.size(); j++ )
            {
                if( energy_tel.size() < 5 || TMath::Abs( energy_tel[j] - median ) < meanAbsoluteError * 5 )
                {
                    if( good_image[j] )
                    {
                        value  += energy_tel[j] * sigma2_tel[j];
                        weight += sigma2_tel[j];
                    }
                }
            }
            if( weight > 0. )
            {
                value /= weight;
            }
            
            // loop over number if images with valid entries
            if( energy_tel.size() > 1 )
            {
                chi2 = 0.;
                double z1 = 0;
                for( unsigned int j = 0; j < energy_tel.size(); j++ )
                {
                    if( sigma2_tel_noRadiusWeigth[j] != 0. 
                    && good_image[j] )
                    {
                        chi2 += ( value - energy_tel[j] ) * ( value - energy_tel[j] ) * sigma2_tel_noRadiusWeigth[j];
                        z1++;
                    }
                }
                if( z1 > 1 )
                {
                    chi2 /= ( z1 - 1. );
                }
                else
                {
                    chi2  = -99.;
                }
                dE   = 0.;
                z1 = 0.;
                for( unsigned int j = 0; j < sigma_tel.size(); j++ )
                {
                    if( sigma_tel[j] > 0.
                    && good_image[j] )
                    {
                        dE += 1. / ( energy_tel[j] * sigma_tel[j] );
                        z1++;
                    }
                }
                if( z1 > 0. )
                {
                    dE = dE / z1;
                }
                else
                {
                    dE = 0.;
                }
            }
            // energy reconstructed from single image
            else if( energy_tel.size() == 1 && value > 0. )
            {
                chi2 = 0.;
                dE = 0.;
            }
            // expect here no good energy
            else
            {
                chi2 = -99;
                dE = -99.;
            }
            return value;
        }
        else
        {
            return -99.;
        }
    }
    
    // should never reach this point
    return -99.;
}

void VTableCalculator::setInterpolationConstants( int iwidth, int iinter )
{
    fInterPolWidth = iwidth;
    fInterPolIter = iinter;
}


double VTableCalculator::getWeightMeanBinContent( TH2F* h, int ix0, int iy0, double x, double y )
{
    if( !h )
    {
        return 0.;
    }
    
    if( h->GetBinContent( ix0, iy0 ) == 0. )
    {
        return 0.;
    }
    
    double ibc = 0.;
    
    // get bin centers
    double i_bc_x0 = h->GetXaxis()->GetBinCenter( ix0 );
    double i_bc_y0 = h->GetYaxis()->GetBinCenter( iy0 );
    int ix1 = ix0;
    int iy1 = iy0;
    if( x < i_bc_x0 && ix0 > 1 )
    {
        ix1 = ix0 - 1;
    }
    else if( ix0 < h->GetNbinsX() - 1 )
    {
        ix1 = ix0 + 1;
    }
    if( y < i_bc_y0 && iy0 > 1 )
    {
        iy1 = iy0 - 1;
    }
    else if( iy0 < h->GetNbinsY() - 1 )
    {
        iy1 = iy0 + 1;
    }
    double i_bc_x1 = h->GetXaxis()->GetBinCenter( ix1 );
    double i_bc_y1 = h->GetYaxis()->GetBinCenter( iy1 );
    
    double weight = 0.;
    double dist = 0.;
    
    // first bin (x,y is inside this bin)
    dist = sqrt( ( i_bc_x0 - x ) * ( i_bc_x0 - x ) + ( i_bc_y0 - y ) * ( i_bc_y0 - y ) );
    // return bin content if x,y is very close to bin center
    if( fabs( dist ) < 1.e-5 )
    {
        return h->GetBinContent( ix0, iy0 );
    }
    dist = 1. / dist;
    
    weight += dist;
    ibc += h->GetBinContent( ix0, iy0 ) * dist;
    
    // second bin
    if( h->GetBinContent( ix1, iy0 ) > 0. )
    {
        dist = sqrt( ( i_bc_x1 - x ) * ( i_bc_x1 - x ) + ( i_bc_y0 - y ) * ( i_bc_y0 - y ) );
        if( dist > 0. )
        {
            weight += 1. / dist;
            ibc += h->GetBinContent( ix1, iy0 ) / dist;
        }
    }
    
    // third bin
    if( h->GetBinContent( ix1, iy1 ) > 0. && dist > 0. )
    {
        dist = sqrt( ( i_bc_x1 - x ) * ( i_bc_x1 - x ) + ( i_bc_y1 - y ) * ( i_bc_y1 - y ) );
        if( dist > 0. )
        {
            weight += 1. / dist;
            ibc += h->GetBinContent( ix1, iy1 ) / dist;
        }
    }
    
    // fourth bin
    if( h->GetBinContent( ix0, iy1 ) > 0. )
    {
        dist = sqrt( ( i_bc_x0 - x ) * ( i_bc_x0 - x ) + ( i_bc_y1 - y ) * ( i_bc_y1 - y ) );
        if( dist > 0. )
        {
            weight += 1. / dist;
            ibc += h->GetBinContent( ix0, iy1 ) / dist;
        }
    }
    
    if( weight > 0. )
    {
        return ibc / weight;
    }
    
    return 0.;
}


void VTableCalculator::setVHistograms( vector< TH2F* >& hM )
{
    hVMedian = hM;
    
    fReadHistogramsFromFile = true;
}


TH2F* VTableCalculator::getHistoMedian()
{
    if( !fReadHistogramsFromFile )
    {
        fReadHistogramsFromFile = readHistograms();
    }
    
    return hMedian;
}


bool VTableCalculator::readHistograms()
{
    if( fOutDir )
    {
        hMedian = ( TH2F* )fOutDir->Get( hMedianName.c_str() );
        
        if( hMedian )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    hMedian = 0;
    return false;
}

/*
 * interpolate in x and y for a 2D histogram
 *
 * (either using the value or the error of the value)
 */
double VTableCalculator::interpolate( TH2F* h, double x, double y, bool iError )
{
    if( !h )
    {
        return -999.;
    }
    
    int i_x = h->GetXaxis()->FindFixBin( x );
    int i_y = h->GetYaxis()->FindFixBin( y );
    // handle under and overflows ( bin nBinsX+1 is needed)
    if( i_x == 0 || i_y == 0 || i_x == h->GetNbinsX() || i_y == h->GetNbinsY() )
    {
        if( iError )
        {
            return h->GetBinError( i_x, i_y );
        }
        else
        {
            return h->GetBinContent( i_x, i_y );
        }
    }
    if( x < h->GetXaxis()->GetBinCenter( i_x ) )
    {
        i_x--;
    }
    if( y < h->GetYaxis()->GetBinCenter( i_y ) )
    {
        i_y--;
    }
    
    double e1 = 0.;
    double e2 = 0.;
    double v = 0.;
    
    // first interpolate on distance axis, then on size axis
    if( !iError )
    {
        e1 = VStatistics::interpolate( h->GetBinContent( i_x, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinContent( i_x, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        e2 = VStatistics::interpolate( h->GetBinContent( i_x + 1, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinContent( i_x + 1, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        v = VStatistics::interpolate( e1, h->GetXaxis()->GetBinCenter( i_x ),
                                      e2, h->GetXaxis()->GetBinCenter( i_x + 1 ),
                                      x, false, 0.5, 1.e-5 );
    }
    else
    {
        e1 = VStatistics::interpolate( h->GetBinError( i_x, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinError( i_x, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        e2 = VStatistics::interpolate( h->GetBinError( i_x + 1, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinError( i_x + 1, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
                                       
        v = VStatistics::interpolate( e1, h->GetXaxis()->GetBinCenter( i_x ), e2, h->GetXaxis()->GetBinCenter( i_x + 1 ), x, false, 0.5, 1.e-5 );
    }
    // final check on consistency of results
    // (don't expect to reconstruct anything below 1 GeV)
    if( e1 > 1.e-3 && e2 < 1.e-3 )
    {
        return e1;
    }
    if( e1 < 1.e-3 && e2 > 1.e-3 )
    {
        return e2;
    }
    
    return v;
}

/*

     search most probable value of energy distribution for a give size/radius bin

     (these distributions are at low energies often very skewed)

*/
void VTableCalculator::fillMPV( TH2F* h, int i, int j, TH1F* h1D, double iMedianValue, double iSigmaValue )
{
    if( !h || !h1D )
    {
        return;
    }
    
    // only fit well filled histograms -> otherwise will median
    if( h1D->GetEntries() <= 50 || iMedianValue <= 0. )
    {
        h->SetBinContent( i, j, iMedianValue );
        return;
    }
    // don't do anything if difference between mean and median is <15%
    if( iMedianValue > 0. && TMath::Abs( ( iMedianValue - h1D->GetMean() ) / iMedianValue ) < 0.15 )
    {
        h->SetBinContent( i, j, iMedianValue );
        return;
    }
    
    /////////////////////////////////////////
    // try a Landau fit
    TF1 iLandau( "iLandau", "TMath::Landau(x,[0],[1],0)*[2]", iMedianValue / 3., iMedianValue * 3. );
    iLandau.SetParameters( iMedianValue, iSigmaValue, h1D->GetEntries() );
    // do not allow the most probable value to be more than x3 off the median
    iLandau.SetParLimits( 0, iMedianValue / 3., iMedianValue * 3. );
    h1D->Fit( &iLandau, "QMNR" );
    // require >10% fit probability to use Landau most probable value
    if( TMath::Prob( iLandau.GetChisquare(), iLandau.GetNDF() ) > 0.1
            && iLandau.GetParameter( 0 ) > 0. && iMedianValue / iLandau.GetParameter( 0 ) < 2.5 )
    {
        h->SetBinContent( i, j, iLandau.GetParameter( 0 ) );
        
        if( fDebug )
        {
            cout << "\t\t Landau interpolation for energy tables: " << i << "\t" << j << "\t";
            cout << TMath::Prob( iLandau.GetChisquare(), iLandau.GetNDF() ) << ", median " << iMedianValue;
            cout << ", fit: " << iLandau.GetParameter( 0 ) << "\t" << iLandau.GetParameter( 1 );
            cout << "\t" << iMedianValue / iLandau.GetParameter( 0 );
            cout << "\t" << iLandau.GetParameter( 0 ) /  iLandau.GetParameter( 1 ) << endl;
        }
    }
    else
    {
        h->SetBinContent( i, j, iMedianValue );
    }
    
}

/*
 *  set min/max values for variable normalization
 *
 *  default values are -9999. for both - means no normalization
 */
void VTableCalculator::setNormalizeTableValues( double i_value_min, double i_value_max )
{
    fValueNormalizationRange_min = i_value_min;
    fValueNormalizationRange_max = i_value_max;
}

