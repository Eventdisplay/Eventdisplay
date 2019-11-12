#include "VLowGainCalibrator.h"

VLowGainCalibrator::VLowGainCalibrator( int run, int sw, bool isInnerHigh, bool sw2monitor, TString dir, TString outdir )
{

    fIsOk = false;
    if( outdir == "" )
    {
        outdir = dir;
    }
    
    cout << "Reading from " << dir << endl;
    cout << "Writing to " << outdir << endl;
    
    fRun = run;
    fWindow = sw;
    if( isInnerHigh )
    {
        fChanMon_start = 250;
        fChanMon_stop = 499;
        fChan_start = 0;
        fChan_stop = 250;
    }
    else
    {
        fChanMon_start = 0;
        fChanMon_stop = 250;
        fChan_start = 250;
        fChan_stop = 499;
    }
    
    
    for( int tel = 0; tel < fNTel; tel++ )
    {
        TString iName = TString::Format( "hist_medians_Tel%d_run_%d", tel + 1, run );
        fMonitorChargeHist[tel] = new TH1D( iName.Data(), "monitor charge;qmon (median) ;entries", 1000, 0, 1000 );
    }
    
    TString iName = TString::Format( "%s/%d.DST.root", dir.Data(), run );
    fInfile = new TFile( iName.Data(), "read" );
    if( !fInfile || fInfile->IsZombie() )
    {
        cout << "ERROR: Input file " << iName << " does not exist, exiting immediately..." << endl;
        return ;
    }
    fDsttree = ( TTree* )fInfile->Get( "dst" );
    if( !fDsttree )
    {
        cout << "ERROR: Input file " << iName << " does not contain a tree named dst, exiting immediately..." << endl;
        return ;
    }
    
    fDsttree->SetBranchAddress( "sumfirst", &sumfirst );
    
    if( sw2monitor )
    {
        fDsttree->SetBranchAddress( "sum", &sumHiLo );
        fDsttree->SetBranchAddress( "sum2", &sumMonitor );
    }
    else
    {
        fDsttree->SetBranchAddress( "sum", &sumMonitor );
        fDsttree->SetBranchAddress( "sum2", &sumHiLo );
    }
    fDsttree->SetBranchAddress( "HiLo", &HiLo );
    fDsttree->SetBranchAddress( "sumwindow", &sumwindow );
    fDsttree->SetBranchAddress( "eventNumber", &eventNumber );
    fDsttree->SetBranchAddress( "dead", &dead );
    fDsttree->SetBranchAddress( "RawMax", &RawMax );
    fDsttree->SetBranchAddress( "pulsetiming", &pulsetiming );
    
    for( int tel = 0; tel < fNTel; tel++ )
    {
        iName.Form( "%s/Tel_%d/", outdir.Data(), tel + 1 );
        gSystem->mkdir( iName.Data(), "r" );
        iName.Form( "%s/Tel_%d/%d.lmult.root", outdir.Data(), tel + 1, run );
        fOutfile[tel] = new TFile( iName.Data(), "recreate" );
        if( !fOutfile[tel] || fOutfile[tel]->IsZombie() )
        {
            cout << "ERROR: Output file " << iName << " could not be opened, exiting immediately..." << endl;
            return;
        }
        iName.Form( "slopes_%d_%d", tel + 1, sw );
        fOuttree[tel] = new TTree( iName.Data(), iName.Data() );
        
        fOuttree[tel]->Branch( "channel",	&fTree_Channel,		"channel/I" );
        fOuttree[tel]->Branch( "m",		&fTree_m, 		"m[2]/D" );
        fOuttree[tel]->Branch( "mErr", 		&fTree_mErr, 		"mErr[2]/D" );
        fOuttree[tel]->Branch( "chi2", 		&fTree_chi2, 		"chi2[2]/D" );
        fOuttree[tel]->Branch( "ndf", 		&fTree_ndf, 		"ndf[2]/I" );
        fOuttree[tel]->Branch( "status", 	&fTree_status, 		"status[2]/I" );
        fOuttree[tel]->Branch( "run", 		&fTree_run,		"run/I" );
        fOuttree[tel]->Branch( "tel", 		&fTree_tel,		"tel/I" );
        
        iName.Form( "%s/Tel_%d/%d.debug.root", outdir.Data(), tel + 1, run );
        fDebugfile[tel] = new TFile( iName.Data(), "recreate" );
        if( !fDebugfile[tel] || fDebugfile[tel]->IsZombie() )
        {
            cout << "ERROR: Output file " << iName << " could not be opened, exiting immediately..." << endl;
            return ;
        }
        
        iName.Form( "debug_%d_%d", tel + 1, sw );
        fDebugtree[tel] = new TTree( iName.Data(), iName.Data() );
        fDebugtree[tel]->Branch( "channel",	&fTree_Channel,		"channel/I" );
        fDebugtree[tel]->Branch( "eventNumber",	&fTree_eventNumber,	"eventNumber/I" );
        fDebugtree[tel]->Branch( "level",	&fTree_level,		"level/I" );
        fDebugtree[tel]->Branch( "hilo"	,	&fTree_hilo,		"hilo/I" );
        fDebugtree[tel]->Branch( "QMon",		&fTree_QMon,		"QMon/D" );
        fDebugtree[tel]->Branch( "QMonMean",	&fTree_QMonMean,	"QMonMean/D" );
        fDebugtree[tel]->Branch( "Q",		&fTree_Q,		"Q/D" );
        fDebugtree[tel]->Branch( "TZero",	&fTree_TZero,		"TZero/D" );
        fDebugtree[tel]->Branch( "RawMax",	&fTree_RawMax,		"RawMax/I" );
        fDebugtree[tel]->Branch( "run", 		&fTree_run,		"run/I" );
        fDebugtree[tel]->Branch( "tel", 		&fTree_tel,		"tel/I" );
        
        
    }
    setLowGainMultiplierUsedInDST();
    setMonitorChargeOptions() ;
    setFitOptions();
    
    fTree_run = fRun;
    fDEBUG = false;
    fIsOk = true;
    
    fMaxDSTEvent = 1000000;
    fMinDSTEvent = 0;
    
}


VLowGainCalibrator::~VLowGainCalibrator()
{
    for( int tel = 0; tel < fNTel; tel++ )
    {
        if( fMonitorChargeHist[tel] )
        {
            fMonitorChargeHist[tel]->Delete();
        }
    }
    if( fDsttree )
    {
        fDsttree->Delete();
    }
}


void VLowGainCalibrator::setMonitorChargeOptions( int nLive_min, double sum_min, bool useMedian, double width )
{
    fNLiveMonitor_min = nLive_min;
    fUseMedian = useMedian;
    fSumMonitor_min = sum_min;
    fLightLevelWidth = width;
    
}

void VLowGainCalibrator::setFitOptions( int n_min, double pure_min, double sat_max, double nevents_min, double prob_min, double b_max )
{
    fFitPure_min = pure_min;
    fFitNPoints_min = n_min;
    fFitProb_min = prob_min ;
    fFitB_max = b_max;
    fFitRSat_max = sat_max;
    fFitNEvents_min = nevents_min;
    
}

/*
     fit function (to be wrapped in a TF1 object)
     assume each light level for the monitoring charge is gaussian distributed
     and above a linear background (pol1 = par[1] + par[2]*x)
     par[0] is number of peaks.
*/
Double_t fpeaks( double* x, double* par )
{
    double result = par[1] + par[2] * x[0];
    for( int p = 0; p < par[0]; p++ )
    {
        double norm  = par[3 * p + 3];
        double mean  = par[3 * p + 4];
        double sigma = par[3 * p + 5];
        result += norm * TMath::Gaus( x[0], mean, sigma );
    }
    return result;
}

double VLowGainCalibrator::calcMeanMonitorCharge( int tel, int ientry )
{
    if( ientry > -1 )
    {
        fDsttree->GetEntry( ientry );
    }
    int nLive = 0;
    double mean = 0;
    for( int iChan = fChanMon_start; iChan < fChanMon_stop; iChan++ )
    {
        if( dead[tel][iChan] || HiLo[tel][iChan] )
        {
            continue;
        }
        if( isNewPixel( tel, iChan ) )
        {
            continue;
        }
        if( sumMonitor[tel][iChan] < fSumMonitor_min )
        {
            continue;
        }
        nLive++;
        mean += sumMonitor[tel][iChan];
    }
    
    if( nLive < 100 )
    {
        return -9999;
    }
    else
    {
        return mean / nLive;
    }
    
}

double VLowGainCalibrator::calcMedianMonitorCharge( int tel, int ientry )
{
    if( ientry > -1 )
    {
        fDsttree->GetEntry( ientry );
    }
    int nLive = 0;
    vector<double> Qmon;
    for( int iChan = fChanMon_start; iChan < fChanMon_stop; iChan++ )
    {
        if( dead[tel][iChan] || HiLo[tel][iChan] )
        {
            continue;
        }
        if( isNewPixel( tel, iChan ) )
        {
            continue;
        }
        if( sumMonitor[tel][iChan] < fSumMonitor_min )
        {
            continue;
        }
        nLive++;
        Qmon.push_back( sumMonitor[tel][iChan] );
    }
    
    if( nLive < 100 )
    {
        return -9999;
    }
    
    std::sort( Qmon.begin(), Qmon.end() );
    double median;
    if( nLive  % 2 )
    {
        median = ( Qmon.at( nLive / 2 - 1 ) + Qmon.at( nLive / 2 ) ) / 2.0;
    }
    else
    {
        median = Qmon.at( nLive / 2 ) ;
    }
    return median;
    
}
double VLowGainCalibrator::calcMonitorCharge( int tel, int ientry )
{

    if( fUseMedian )
    {
        return calcMedianMonitorCharge( tel, ientry );
    }
    else
    {
        return calcMeanMonitorCharge( tel, ientry );
    }
    
}


bool VLowGainCalibrator::makeMonitorChargeHists( )
{

    for( int tel = 0; tel < fNTel; tel++ )
    {
        fMonitorChargeHist[tel]->Reset();
    }
    
    for( unsigned int iEntry = fMinDSTEvent; iEntry < fDsttree->GetEntries() && iEntry < fMaxDSTEvent; iEntry++ )
    {
        fDsttree->GetEntry( iEntry );
        
        for( int tel = 0; tel < 4; tel++ )
        {
        
            double median = calcMonitorCharge( tel );
            
            fMonitorChargeHist[tel]->Fill( median );
            
        }
        
    }
    return true;
    
}

void VLowGainCalibrator::resetLightLevels()
{
    for( int i = 0; i < fNTel; i++ )
    {
        resetLightLevels( i );
    }
}

void VLowGainCalibrator::resetLightLevels( int tel )
{
    fNLightLevels[tel] = 0;
    fLightLevelMean[tel].clear();
    fLightLevelMeanError[tel].clear();
    fLightLevelSigma[tel].clear();
    fLightLevelSigmaError[tel].clear();
}


/*

    main function to identify different light levels from the monitoring charge histograms
    and perform some tests on the goodness of the light levels (loop through all telescopes)

*/
bool VLowGainCalibrator::findLightLevels( bool iDraw )
{
    for( int tel = 0; tel < fNTel; tel++ )
    {
        findLightLevels( tel, 2 , iDraw );
        
        int test = checkLightLevels( tel );
        if( test == 0 )
        {
            cout << "       This is very likely a bad HiLo run for telescope " << tel + 1 << "." << endl << endl;
            resetLightLevels( tel );
        }
        else if( test == 2 )
        {
            cout << "         Rerun the light level finder for telescope " << tel + 1 << " with higher peak significance." << endl;
            findLightLevels( tel, 3, iDraw );
            
            int finaltest = checkLightLevels( tel );
            if( finaltest != 1 )
            {
                cout << "       Resetting light levels for telescope " << tel + 1 << " (likely a bad HiLo run)." << endl << endl;
                resetLightLevels( tel );
            }
        }
    }
    return true;
}


/*

     fill the different light levels from the monitoring charge histograms

*/
void VLowGainCalibrator::findLightLevels( int tel, int iPeakSignificance , bool iDraw )
{
    resetLightLevels( tel );
    // check first that monitoring charge histogram exists
    if( ( int )fMonitorChargeHist[tel]->GetEntries() == 0 )
    {
        cout << "ERROR: Could not find histogram with monitoring charges for telescope " << tel + 1 << "." << endl;
        cout << "       Have you run VLowGainCalibrator::makeMonitorChargeHists() before?" << endl;
        return;
    }
    
    TCanvas* c1 = NULL;
    if( iDraw )
    {
        c1 = new TCanvas( "c1", "c1", 10, 10, 800, 600 );
        c1->Divide( 1, 2 );
        c1->cd( 1 );
        fMonitorChargeHist[tel]->Draw();
    }
    
    // exclude the zero light level for finding the light levels
    // (might be greater than zero in case the median is used)
    fMonitorChargeHist[tel]->GetXaxis()->SetRangeUser( 20, 1000 );
    
    // find peaks (there should be 7 if the flasher is ok)
    int iLightLevel = 7;
    TSpectrum* s = new TSpectrum( 2 * iLightLevel );
    if( iDraw )
    {
        s->Search( fMonitorChargeHist[tel], iPeakSignificance, "", 0.05 );
        c1->cd( 2 );
    }
    else
    {
        s->Search( fMonitorChargeHist[tel], iPeakSignificance, "goff", 0.05 );
    }
    
    int iNfound = s->GetNPeaks();
    if( fDEBUG )
    {
        s->Print();
    }
    else
    {
        cout << "Found " << iNfound << " possible light levels in tel " << tel + 1 << endl;
    }
    
    
    bool rebin = false;
    int cnt = 0;
    TH1D* ihist;
    while( iNfound > iLightLevel )
    {
        ihist = ( TH1D* )fMonitorChargeHist[tel]->Clone( "ihist2" );
        ihist->Rebin( cnt + 1 );
        
        if( iDraw )
        {
            ihist->Draw();
            s->Search( ihist, iPeakSignificance, "", 0.05 );
        }
        else
        {
            s->Search( ihist, iPeakSignificance, "goff", 0.05 );
        }
        
        iNfound = s->GetNPeaks();
        
        cout << "Rebin histogram (" << cnt + 1 << ")" << endl;
        if( fDEBUG )
        {
            s->Print();
        }
        
        rebin = true;
        
        // make sure that at some point the while loop is exited
        if( cnt == 3 )
        {
            cout << endl;
            cout << "WARNING: Even after three rebinning attempts more then " << iNfound << " light levels are found for tel " << tel + 1  << "." << endl;
            cout << "         Please manually inspect the histogram for goodness of this calibration run." << endl << endl;
        }
        else if( cnt >= 10 || cnt >= ( fMonitorChargeHist[tel]->GetEntries() / fMonitorChargeHist[tel]->GetNbinsX() ) )
        {
            cout << "ERROR: Too many peaks found even after " << cnt << " rebinning attempts for tel " << tel + 1  << "." << endl;
            cout << "       Set number of peaks to zero! " << endl;
            cout << "       exiting... " << endl;
            iNfound = 0;
            return;
        }
        cnt++;
    }
    TH1D* h2;
    if( rebin )
    {
        h2 = ( TH1D* )ihist->Clone( "h2" );
    }
    else
    {
        h2 = ( TH1D* )fMonitorChargeHist[tel]->Clone( "h2" );
    }
    if( ( int )h2->GetEntries() == 0 )
    {
        cout << "ERROR: Could not find histogram to be used for fitting the light level for tel " << tel + 1  << "." << endl;
        return;
    }
    
    
    // estimate background using a linear fitting method
    TF1* iline = new TF1( "iline", "pol1", 20, 1000 );
    fMonitorChargeHist[tel]->Fit( "iline", "qn" );
    fMonitorChargeHist[tel]->GetXaxis()->SetRangeUser( 0, 1000 );
    
    double par[300];
    par[1] =  0;
    par[2] =  0;
    
    // eliminate peaks with light level zero
    // and too small peaks at the background level
    int npeaks = 0;
    for( int i = 0; i < iNfound; i++ )
    {
        double ix = s->GetPositionX()[i];
        int ibin  = h2->GetXaxis()->FindBin( ix );
        double iy = h2->GetBinContent( ibin );
        
        if( ibin <= 1 )
        {
            continue;
        }
        if( iy - TMath::Sqrt( iy ) < iline->Eval( ix ) )
        {
            continue;
        }
        
        // fill paramter estimates for the gaussian fits
        par[3 * npeaks + 3] = iy;
        par[3 * npeaks + 4] = ix;
        par[3 * npeaks + 5] = 3;
        npeaks++;
    }
    par[0] = npeaks;
    cout << "Finally used " << npeaks << " peaks for fitting the light level parameters" << endl;
    
    if( npeaks > iLightLevel )
    {
        cout << "ERROR: Found more than " << iLightLevel << " light levels (this should not happen). " << endl
             << "       Please check if distinct light levels exist." << endl;
        return;
    }
    else
    {
        cout << "Now fitting: Be patient" << endl;
        
        TF1* fit = new TF1( "fit", fpeaks, 0, 1000, 3 + 3 * npeaks );
        TVirtualFitter::Fitter( h2, 3 + 3 * npeaks );
        fit->SetParameters( par );
        fit->SetParLimits( 0, npeaks, npeaks );
        fit->SetNpx( 1000 );
        
        TString opt = "Q";
        if( fDEBUG )
        {
            opt = "";
        }
        
        if( iDraw )
        {
            h2->Fit( "fit", opt );
            c1->Update();
            TString iName = TString::Format( "light_t%d_%d", tel + 1, fWindow );
            c1->SetName( iName.Data() );
            fDebugfile[tel]->cd();
            c1->Write();
            c1->Close();
        }
        else
        {
            h2->Fit( "fit", "N" );
        }
        
        TString iStatus( gMinuit->fCstatu );
        if( iStatus.Contains( "CONVERGED" ) )
        {
        
            if( fit->GetParameter( 1 ) > 10 )
            {
                cout << endl;
                cout << "WARNING: The estimated linear background seems high for tel " << tel + 1  << "." << endl;
                cout << "         Please manually inspect the histogram for goodness of this calibration run." << endl << endl;
            }
            
            fNLightLevels[tel] = npeaks;
            if( fDEBUG )
            {
                cout << "#light levels in tel " << tel + 1 << " : " << npeaks << endl;
            }
            
            for( unsigned int i = 0; i < fNLightLevels[tel]; i++ )
            {
                fLightLevelMean[tel].push_back( fit->GetParameter( 3 * i + 4 ) );
                fLightLevelMeanError[tel].push_back( fit->GetParError( 3 * i + 4 ) );
                if( fDEBUG )
                {
                    cout << "mean:   " << fLightLevelMean[tel][i] << " +/- " << fLightLevelMeanError[tel][i] << endl;
                }
                
                fLightLevelSigma[tel].push_back( TMath::Abs( fit->GetParameter( 3 * i + 5 ) ) );
                fLightLevelSigmaError[tel].push_back( fit->GetParError( 3 * i + 5 ) );
                if( fDEBUG )
                {
                    cout << "sigma:  " << fLightLevelSigma[tel][i] << " +/- " << fLightLevelSigmaError[tel][i] << endl;
                }
            }
        }
        else
        {
            cout << "ERROR: Fit of the peaks did not converge for telescope " << tel + 1 << "!" << endl;
            cout << "Fit status: " << iStatus << endl;
            resetLightLevels( tel );
        }
        if( fit )
        {
            //fit->Delete();
        }
    }
    h2->Delete();
    iline->Delete();
    
    delete s;
}


/*

   Routine to test the goodness of the found light levels

     Return values include the following checks:
     - Do light levels exist?
       If not, return value == 0
     - Are the peaks distinct/well seperated?
       If not, return value == 2

   If everything seems ok so far, return value should be 1.

*/
int VLowGainCalibrator::checkLightLevels( int tel )
{

    // check if light levels have been set
    if( fNLightLevels[tel] == 0 )
    {
        cout << endl << "ERROR: No light levels found for telescope " << tel + 1 << "." << endl;
        return 0;
    }
    
    //  check if light levels are distinct
    //  default: they should not overlap at the 2 sigma level
    bool twopeaks = false;
    for( unsigned int i = 0; i < fNLightLevels[tel] - 1; i++ )
    {
        for( unsigned int k = ( i + 1 ); k < fNLightLevels[tel]; k++ )
        {
            if( fLightLevelMean[tel][i] < fLightLevelMean[tel][k] )
            {
                if( fLightLevelMean[tel][i] + fLightLevelWidth * fLightLevelSigma[tel][i] < fLightLevelMean[tel][k] - fLightLevelWidth * fLightLevelSigma[tel][k] )
                {
                    continue;
                }
                else
                {
                    cout << "WARNING: Found two light levels close to each other. " << endl;
                    if( fDEBUG ) cout << fLightLevelMean[tel][i] << " + " << fLightLevelWidth << "*" << fLightLevelSigma[tel][i]
                                          << " and " << fLightLevelMean[tel][k] << " - " << fLightLevelWidth << "*" << fLightLevelSigma[tel][k] << endl;
                    twopeaks = true;
                }
            }
            else if( fLightLevelMean[tel][i] > fLightLevelMean[tel][k] )
            {
                if( fLightLevelMean[tel][i] - fLightLevelWidth * fLightLevelSigma[tel][i] > fLightLevelMean[tel][k] + fLightLevelWidth * fLightLevelSigma[tel][k] )
                {
                    continue;
                }
                else
                {
                    cout << "WARNING: Found two light levels close to each other. " << endl;
                    if( fDEBUG ) cout << fLightLevelMean[tel][i] << " + " << fLightLevelWidth << "*" << fLightLevelSigma[tel][i]
                                          << " and " << fLightLevelMean[tel][k] << " - " << fLightLevelWidth << "*" <<  fLightLevelSigma[tel][k] << endl;
                    twopeaks = true;
                }
            }
            // don't compare the same entries/values (should not happen, but who knows...)
            else if( TMath::Abs( fLightLevelMean[tel][i] - fLightLevelMean[tel][k] ) < 1e-3 )
            {
                continue;
            }
            else
            {
                cout << "Warning: Wow, you should have never reached this point... " << endl;
            }
        }
    }
    if( twopeaks )
    {
        return 2;
    }
    else
    {
        return 1;
    }
    
}




bool VLowGainCalibrator::calculateMeanCharges()
{
    //initialise vectors
    for( int tel = 0; tel < fNTel; tel++ )
    {
        for( int iChan = fChan_start; iChan < fChan_stop; iChan++ )
        {
            for( int hilo = 0; hilo < 2; hilo++ )
            {
                fN[tel][iChan][hilo].assign( fNLightLevels[tel], 0 );
                fNSat[tel][iChan][hilo].assign( fNLightLevels[tel], 0 );
                fY[tel][iChan][hilo].assign( fNLightLevels[tel], 0 );
                fY2[tel][iChan][hilo].assign( fNLightLevels[tel], 0 );
                
            }//hilo
        }//chan
    }//tel
    
    for( unsigned int ientry = fMinDSTEvent; ientry < fDsttree->GetEntries() && ientry < fMaxDSTEvent ; ientry++ )
    {
    
        fDsttree->GetEntry( ientry );
        
        for( int tel = 0; tel < fNTel; tel++ )
        {
        
            fTree_tel = tel + 1;
            
            //if( fNLightLevels[tel] == 0 ) continue;	//don't bother if there are no light levels.
            
            double qmon = calcMonitorCharge( tel );
            int level = -1;
            for( unsigned int ilevel = 0; ilevel < fNLightLevels[tel]; ilevel++ )
            {
                if( qmon < fLightLevelMean[tel][ilevel] + fLightLevelWidth * fLightLevelSigma[tel][ilevel] &&  qmon > fLightLevelMean[tel][ilevel] - fLightLevelWidth * fLightLevelSigma[tel][ilevel] )
                {
                    level = ilevel;
                    break;
                }
            } //levels
            
            //first calculate mean/sdev of start time for hi/lo channels.
            
            double start[2];
            start[0] = 0;
            start[1] = 0;
            
            double start2[2];
            start2[0] = 0;
            start2[1] = 0;
            
            double N[2];
            N[0] = 0;
            N[1] = 0;
            
            for( int iChan = fChan_start; iChan < fChan_stop; iChan++ )
            {
            
                if( dead[tel][iChan] )
                {
                    continue;
                }
                
                N	[ HiLo[tel][iChan] ] ++;
                
                start	[ HiLo[tel][iChan] ] += sumfirst[tel][iChan];
                start2	[ HiLo[tel][iChan] ] += sumfirst[tel][iChan] * sumfirst[tel][iChan];
                
            }//chan
            
            for( int hilo = 0; hilo < 2; hilo++ )
            {
                if( N[hilo] == 0 )
                {
                    continue ;
                }
                start[hilo] /= N[hilo];
                start2[hilo] /= N[hilo];
                start2[hilo] = sqrt( start2[hilo] - start[hilo] * start[hilo] );
            }//hilo
            
            for( int iChan = fChan_start; iChan < fChan_stop; iChan++ )
            {
            
                if( dead[tel][iChan] )
                {
                    continue;
                }
                //todo test this			if( fabs( sumfirst[tel][iChan] - start[ HiLo[tel][iChan] ] ) > 2*start2[ HiLo[tel][iChan] ] ) continue;
                
                if( level > -1 )
                {
                    fN   [tel][iChan][ HiLo[tel][iChan] ][level]++;
                    fNSat[tel][iChan][ HiLo[tel][iChan] ][level] += ( RawMax[tel][iChan] == 255 ? 1 : 0 );
                    fY   [tel][iChan][ HiLo[tel][iChan] ][level] += sumHiLo[tel][iChan] / ( HiLo[tel][iChan] ? fLMult[tel] : 1.0 );
                    fY2  [tel][iChan][ HiLo[tel][iChan] ][level] += TMath::Power( sumHiLo[tel][iChan] / ( HiLo[tel][iChan] ? fLMult[tel] : 1.0 ) , 2 ) ;
                }
                if( isDebugChannel( iChan ) )
                {
                    fTree_eventNumber = eventNumber;
                    fTree_Channel = iChan;
                    fTree_level = level;
                    fTree_hilo = HiLo[tel][iChan];
                    fTree_Q = sumHiLo[tel][iChan] / ( HiLo[tel][iChan] ? fLMult[tel] : 1.0 ) ;
                    fTree_QMon = qmon;
                    fTree_RawMax = RawMax[tel][iChan];
                    fTree_TZero = pulsetiming[tel][0][iChan];
                    if( level > -1 )
                    {
                        fTree_QMonMean = fLightLevelMean[tel][level];
                    }
                    else
                    {
                        fTree_QMonMean = -1;
                    }
                    fDebugtree[tel]->Fill();
                }
            }//chan
            
        }//tel
        
    }//ientry
    
    //now fix normalization etc.
    
    for( int tel = 0; tel < fNTel; tel++ )
    {
        for( int iChan = fChan_start; iChan < fChan_stop; iChan++ )
        {
            for( int hilo = 0; hilo < 2; hilo++ )
            {
                for( unsigned int ilevel = 0; ilevel < fNLightLevels[tel]; ilevel++ )
                {
                    if( fN[tel][iChan][hilo][ilevel] == 0 )
                    {
                        continue;
                    }
                    fY[tel][iChan][hilo][ilevel] /=  fN[tel][iChan][hilo][ilevel];
                    fY2[tel][iChan][hilo][ilevel] /= fN[tel][iChan][hilo][ilevel];
                }//level
            }//hilo
        }//chan
    }//tel
    return true;
}

bool VLowGainCalibrator::doTheFit()
{

    for( int tel = 0; tel < fNTel; tel++ )
    {
    
        int nGood[3] = { 0, 0, 0 };
        int nBad[3] = { 0, 0, 0 };
        
        
        fTree_tel = tel + 1;
        
        for( int iChan = fChan_start; iChan < fChan_stop; iChan++ )
        {
        
            fTree_Channel = iChan;
            
            for( int hilo = 0; hilo < 2; hilo++ )
            {
                TString name = TString::Format( "graph_t%d_c%d_%d_%d", tel + 1, iChan, fWindow, hilo );
                TGraphErrors* t = new TGraphErrors( 0 );
                t->SetName( name.Data() );
                t->SetTitle( name.Data() );
                TF1* f = new TF1( "f", "[0]*x", 0, 1000 );
                int N = 0;
                for( unsigned int i = 0; i < fNLightLevels[tel]; i++ )
                {
                    if( fN[tel][iChan][hilo][i] > fFitNEvents_min 									//min number events
                            && ( double )fN[tel][iChan][hilo][i] / ( fN[tel][iChan][0][i] + fN[tel][iChan][1][i] ) > fFitPure_min	//min purity
                            && ( double )fNSat[tel][iChan][hilo][i] / fN[tel][iChan][hilo][i] < fFitRSat_max				//max saturation
                      )
                    {
                    
                        t->SetPoint( N, fLightLevelMean[tel][i], fY[tel][iChan][hilo][i] );
                        double eX, eY;
                        eX = fLightLevelSigma[tel][i];
                        eY = sqrt( fY2[tel][iChan][hilo][i] - fY[tel][iChan][hilo][i] * fY[tel][iChan][hilo][i] ) / sqrt( fN[tel][iChan][hilo][i] );
                        t->SetPointError( N, eX, eY );
                        N++;
                    }
                }// light levels
                if( N < fFitNPoints_min )
                {
                    if( fDEBUG )
                    {
                        cout << "Warning: Less than " << fFitNPoints_min << " point" << ( fFitNPoints_min == 1 ? "" : "s" ) << " for Tel " << tel + 1 << ", channel " << iChan  << " " << hilo ;
                        cout << ", not fitting." << endl;
                    }
                    fTree_status[hilo] = NO_POINTS;
                    fTree_m[hilo] = -1;
                    fTree_mErr[hilo] = 99999;
                    fTree_chi2[hilo] = 99999;
                    fTree_ndf[hilo]	= -1;
                    
                }
                else
                {
                
                    TString opt = "";
                    if( !fDEBUG )
                    {
                        opt = "Q";
                    }
                    t->Fit( f, opt.Data() );
                    
                    if( fDEBUG || isDebugChannel( iChan ) )
                    {
                        fDebugfile[tel]->cd();
                        t->Write();
                        
                    }
                    
                    fTree_m[hilo] = f->GetParameter( 0 );
                    fTree_mErr[hilo] = f->GetParError( 0 );
                    //double b = f->GetParameter(1);
                    //double berr = f->GetParError(1);
                    
                    fTree_chi2[hilo] = f->GetChisquare();
                    fTree_ndf[hilo]	= f->GetNDF();
                    double p = TMath::Prob( fTree_chi2[hilo], fTree_ndf[hilo] );
                    
                    if( p < fFitProb_min )
                    {
                        fTree_status[hilo] = BAD_CHI2;
                    }
                    //else if( TMath::Abs( b ) / berr  > fFitB_max ) fTree_status[hilo] = NOT_PROPORTIONAL;
                    else
                    {
                        fTree_status[hilo] = GOOD;
                    }
                    
                    
                }
                f->Delete();
                t->Delete();
                
            }//hilo
            
            if( fTree_status[0] == GOOD )
            {
                nGood[0]++;
            }
            else
            {
                nBad[0]++;
            }
            
            if( fTree_status[1] == GOOD )
            {
                nGood[1]++;
            }
            else
            {
                nBad[1]++;
            }
            
            if( fTree_status[0] == GOOD && fTree_status[1] == GOOD )
            {
                nGood[2]++;
            }
            else
            {
                nBad[2]++;
            }
            
            fOuttree[tel]->Fill();
        }//chan
        
        TString out = TString::Format( "RUN %d T%d good/bad: %3d/%3d high,  %3d/%3d low,  %3d/%3d total." , fRun, tel + 1, nGood[0], nBad[0],  nGood[1], nBad[1],  nGood[2], nBad[2] );
        cout << out << endl;
        
    }//tel
    
    return true;
}

bool VLowGainCalibrator::terminate( )
{
    for( int tel = 0; tel < fNTel; tel++ )
    {
    
        if( fOutfile[tel] )
        {
            fOutfile[tel]->cd();
            if( fOuttree[tel] )
            {
                fOuttree[tel]->Write();
            }
        }
        
        if( fDebugfile[tel] )
        {
            fDebugfile[tel]->cd();
            if( fDebugtree[tel] && fDebugtree[tel]->GetEntries() > 0 )
            {
                fDebugtree[tel]->Write();
            }
            if( fOuttree[tel] )
            {
                fOuttree[tel]->Write();
            }
        }
        
        if( fOutfile[tel] )
        {
            fOutfile[tel]->Close();
        }
        if( fDebugfile[tel] )
        {
            fDebugfile[tel]->Close();
        }
    }
    return true;
}


bool VLowGainCalibrator::isNewPixel( int tel, int iChan )
{
    //catches runs that had test pixels (with the new PMTs) in T3. The gain/gain vs HV was different and we don't want to use them to determine the monitor charge.
    if( tel != 2 )
    {
        return false;
    }
    if( fRun < 54239 || fRun > 63195 )
    {
        return false;
    }
    if( iChan == 280 || iChan == 281 || iChan == 340 || iChan == 341 || iChan == 342 || iChan == 406 || iChan == 408 )
    {
        return true;
    }
    if( iChan == 291 || iChan == 292 || iChan == 253 || iChan == 254 || iChan == 419 || iChan == 420 )
    {
        return true;
    }
    return false;
}




