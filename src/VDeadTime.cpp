/*! \class VDeadTime
    \brief dead time calculator

    two methods are applied:

    i)  using the distribution of time differences
    ii) using the ten MHz clock and busy counter (scalar)

*/

#include "VDeadTime.h"

VDeadTime::VDeadTime( bool iIsOn )
{
    bIsOn = iIsOn;
    
    fDeadTimeMiss = 0.;
    fDeadTimeMS = 0.;
    fDeadTimeFrac = 0.;
    
    fDeadTimeFrac_status = 0;
    setDeadTimeConsistencyCheck_allowedDifference();
    
    fTFitMin = 0.004;
    fTFitMax = 0.015;
    
    hTimeDiff = 0;
    hTimeDiff2D = 0;
    hTimeDiffLog = 0;
    hFTimeDiff = 0;
    hgDeadTime = 0;
    hNEventTime = 0;
    
    // scalar
    hScalarClock = 0;
    hScalarBusy = 0;
    hScalarDeadTimeFraction = 0;
    
    fScalarDeadTimeFrac = 0.;
    fScalarDeadTimeChi2 = 0.;
    
    hisList = 0;
    reset();
}


void VDeadTime::reset()
{
    ft0 = 0.;
    fRunStart = 0.;
    
    if( hTimeDiff )
    {
        hTimeDiff->Reset();
    }
    if( hTimeDiffLog )
    {
        hTimeDiffLog->Reset();
    }
    if( hTimeDiff2D )
    {
        hTimeDiff2D->Reset();
    }
    if( hScalarClock )
    {
        hScalarClock->Reset();
    }
    if( hScalarBusy )
    {
        hScalarBusy->Reset();
    }
    if( hScalarDeadTimeFraction )
    {
        hScalarDeadTimeFraction->Reset();
    }
    
    fScalarClockLastEvent = 0;
    fScalarBusyLastEvent = 0;
    fFirstScalarEvent = true;
    fScalarTimeOfFirstEvent = 0;
}


void VDeadTime::defineHistograms( float iRunDuration, bool iNoWarning )
{
    hisList = new TList();
    
    char hname[200];
    
    // histograms for calculating dead vs
    if( bIsOn )
    {
        sprintf( hname, "hTimeDiff_on" );
    }
    else
    {
        sprintf( hname, "hTimeDiff_off" );
    }
    hTimeDiff = new TH1D( hname, "time difference between events (lin)", 50000, 0., 0.2 );
    hTimeDiff->SetXTitle( "time difference [s]" );
    hTimeDiff->SetYTitle( "number of entries" );
    hTimeDiff->SetFillColor( 8 );
    hisList->Add( hTimeDiff );
    
    if( bIsOn )
    {
        sprintf( hname, "hTimeDiffLog_on" );
    }
    else
    {
        sprintf( hname, "hTimeDiffLog_off" );
    }
    hTimeDiffLog = new TH1D( hname, "time difference between events (log)", 500, -4., -0.5 );
    hTimeDiffLog->SetXTitle( "log_{10} time difference [s]" );
    hTimeDiffLog->SetYTitle( "number of entries" );
    hTimeDiffLog->SetLineWidth( 2 );
    hTimeDiffLog->SetFillColor( 8 );
    hisList->Add( hTimeDiffLog );
    
    if( bIsOn )
    {
        sprintf( hname, "hTimeDiff2D_on" );
    }
    else
    {
        sprintf( hname, "hTimeDiff2D_off" );
    }
    int ibins = ( int )( iRunDuration / 60.*1.2 / 2. );
    if( ibins < 1 )
    {
        ibins = 1;
    }
    hTimeDiff2D = new TH2D( hname, "time difference between events vs run time", ibins, 0., iRunDuration / 60.*1.2, 4000, 0., 0.05 );
    hTimeDiff2D->SetXTitle( "event time [min]" );
    hTimeDiff2D->SetYTitle( "time difference [s]" );
    hTimeDiff2D->SetZTitle( "number of entries" );
    hTimeDiff2D->SetFillColor( 8 );
    hisList->Add( hTimeDiff2D );
    
    if( bIsOn )
    {
        sprintf( hname, "hgDeadTime_on" );
    }
    else
    {
        sprintf( hname, "hgDeadTime_off" );
    }
    hgDeadTime = new TGraphErrors( hTimeDiff2D->GetNbinsX() );
    hgDeadTime->SetName( hname );
    hgDeadTime->SetTitle( "" );
    hgDeadTime->SetMarkerStyle( 20 );
    hgDeadTime->SetMarkerSize( 2 );
    hisList->Add( hgDeadTime );
    
    if( bIsOn )
    {
        sprintf( hname, "hNEventTime_on" );
    }
    else
    {
        sprintf( hname, "hNEventTime_off" );
    }
    hNEventTime = new TH1D( hname, "rate", ibins, 0., iRunDuration / 60.*1.2 );
    hNEventTime->SetXTitle( "event time [min]" );
    hNEventTime->SetYTitle( "number of entries" );
    hNEventTime->SetFillColor( 8 );
    hisList->Add( hNEventTime );
    
    if( bIsOn )
    {
        sprintf( hname, "hTimeDiffLog_on" );
    }
    
    if( bIsOn )
    {
        sprintf( hname, "hFTimeDiff_on" );
    }
    else
    {
        sprintf( hname, "hFTimeDiff_off" );
    }
    //    hFTimeDiff = new TF1( hname, "expo", fTFitMin, fTFitMax );
    hFTimeDiff = new TF1( hname, "expo", 1.e-4, 0.2 );
    hFTimeDiff->SetLineColor( 2 );
    hisList->Add( hFTimeDiff );
    
    // scalar related histograms
    if( iRunDuration > 0. )
    {
        if( bIsOn )
        {
            sprintf( hname, "hScalarClock_on" );
        }
        else
        {
            sprintf( hname, "hScalarClock_off" );
        }
        hScalarClock = new TH1D( hname, "scalar clock", int( ( iRunDuration * 1.1 ) / 10. ), 0., iRunDuration * 1.1 );
        hScalarClock->Sumw2();
        hScalarClock->SetXTitle( "event time [s]" );
        hScalarClock->SetYTitle( "10MHz counts" );
        hisList->Add( hScalarClock );
        
        if( bIsOn )
        {
            sprintf( hname, "hScalarBusy_on" );
        }
        else
        {
            sprintf( hname, "hScalarBusy_off" );
        }
        hScalarBusy = new TH1D( hname, "scalar busy counter", int( ( iRunDuration * 1.1 ) / 10. ), 0., iRunDuration * 1.1 );
        hScalarBusy->SetXTitle( "event time [s]" );
        hScalarBusy->SetYTitle( "10MHz counts (busy)" );
        hScalarBusy->SetLineColor( 2 );
        hScalarBusy->SetMarkerColor( 2 );
        hScalarBusy->Sumw2();
        hisList->Add( hScalarBusy );
        
        if( bIsOn )
        {
            sprintf( hname, "hScalarDeadTimeFraction_on" );
        }
        else
        {
            sprintf( hname, "hScalarDeadTimeFraction_off" );
        }
        hScalarDeadTimeFraction = new TH1D( hname, "scalar dead time fraction", int( ( iRunDuration * 1.1 ) / 10. ), 0., iRunDuration * 1.1 );
        hScalarDeadTimeFraction->SetXTitle( "event time [s]" );
        hScalarDeadTimeFraction->SetYTitle( "dead time fraction" );
        hScalarDeadTimeFraction->Sumw2();
        hisList->Add( hScalarDeadTimeFraction );
    }
    else if( !iNoWarning )
    {
        cout << "Warning: scalars not used for dead time calculation, run duration zero (" << iRunDuration << ")" << endl;
    }
    
}

double VDeadTime::fillDeadTime( double time, unsigned int* tenMHzClock )
{
    double iDiff = fillTimeDifferenceHistograms( time );
    
    if( tenMHzClock )
    {
        fillTenMHzClockArray( time, tenMHzClock );
    }
    
    return iDiff;
}

void VDeadTime::fillTenMHzClockArray( double time, unsigned int* tenMHzClock )
{
    if( !tenMHzClock || !hScalarClock )
    {
        return;
    }
    // check if this is the first scalar event
    if( fFirstScalarEvent )
    {
        fFirstScalarEvent = false;
        fScalarClockLastEvent = tenMHzClock[0];
        fScalarBusyLastEvent  = tenMHzClock[1];
        fScalarTimeOfFirstEvent = time;
        return;
    }
    
    // check time bins and reset variables if needed
    int nbin = hScalarClock->FindBin( time - fScalarTimeOfFirstEvent );
    
    // fill scalar clock
    if( tenMHzClock[0] >= fScalarClockLastEvent )
    {
        hScalarClock->SetBinContent( nbin, hScalarClock->GetBinContent( nbin ) + double( tenMHzClock[0] - fScalarClockLastEvent ) );
    }
    else
    {
        hScalarClock->SetBinContent( nbin, hScalarClock->GetBinContent( nbin )
                                     + TMath::Power( 2., 32. ) - 1. - double( fScalarClockLastEvent )
                                     +  double( tenMHzClock[0] ) );
    }
    // fill busy clock
    if( tenMHzClock[1] >= fScalarBusyLastEvent )
    {
        hScalarBusy->SetBinContent( nbin, hScalarBusy->GetBinContent( nbin ) + double( tenMHzClock[1] - fScalarBusyLastEvent ) );
    }
    else
    {
        hScalarBusy->SetBinContent( nbin, hScalarBusy->GetBinContent( nbin )
                                    + TMath::Power( 2., 32. ) - 1. - double( fScalarBusyLastEvent )
                                    +  double( tenMHzClock[1] ) );
    }
    
    fScalarClockLastEvent = tenMHzClock[0];
    fScalarBusyLastEvent  = tenMHzClock[1];
    
}


double VDeadTime::fillTimeDifferenceHistograms( double time )
{
    double tdiff = time - ft0;
    
    if( tdiff < 100. )
    {
        if( hTimeDiff )
        {
            hTimeDiff->Fill( tdiff );
        }
        if( hTimeDiffLog && tdiff > 0. )
        {
            hTimeDiffLog->Fill( log10( tdiff ) );
        }
        if( fRunStart > 0. && hTimeDiff2D && hNEventTime )
        {
            hTimeDiff2D->Fill( ( time - fRunStart ) / 60., tdiff );
            hNEventTime->Fill( ( time - fRunStart ) / 60. );
        }
    }
    ft0 = time;
    
    if( fRunStart == 0. )
    {
        fRunStart = time;
    }
    
    return tdiff;
}


/*

   calculate dead time and check consistency of result

*/
double VDeadTime::calculateDeadTime()
{
    fScalarDeadTimeFrac = calculateDeadTimeFromScalars();
    
    fDeadTimeFrac = calculateDeadTimeFromTimeDifferences();
    
    // check for consistency between the two methods
    if( fScalarDeadTimeFrac > 0. && TMath::Abs( fScalarDeadTimeFrac - fDeadTimeFrac ) < fDeadTimeConsistencyCheck_allowedDifference )
    {
        fDeadTimeFrac_status = 1;
        return fScalarDeadTimeFrac;
    }
    
    // scalar method failed
    fDeadTimeFrac_status = 0;
    
    return fDeadTimeFrac;
}

/*
    calculate time-dependent dead time

    function returns average deadtime over whole run

*/
double VDeadTime::calculateDeadTimeFromScalars()
{
    if( !hScalarClock || !hScalarBusy || !hScalarDeadTimeFraction )
    {
        return 0.;
    }
    
    double iClock = 0.;
    double iBusy = 0.;
    double iChi2 = 0.;
    double iNDF = 0.;
    // calculate mean dead time fraction and errors on dead time histogram
    for( int i = 1; i <= hScalarClock->GetNbinsX(); i++ )
    {
        iClock += hScalarClock->GetBinContent( i );
        iBusy +=  hScalarBusy->GetBinContent( i );
        
        if( hScalarClock->GetBinContent( i ) > 0. )
        {
            hScalarClock->SetBinError( i, sqrt( hScalarClock->GetBinContent( i ) ) );
        }
        if( hScalarBusy->GetBinContent( i ) > 0. )
        {
            hScalarBusy->SetBinError( i, sqrt( hScalarBusy->GetBinContent( i ) ) );
        }
        
    }
    
    // calculate time-dependent dead time (use binomial errors)
    hScalarDeadTimeFraction->Divide( hScalarBusy, hScalarClock, 1., 1., "B" );
    
    // make sure that dead time is never > 1
    // (this should never happen)
    for( int i = 1; i <= hScalarDeadTimeFraction->GetNbinsX(); i++ )
    {
        if( hScalarDeadTimeFraction->GetBinContent( i ) > 1. )
        {
            hScalarDeadTimeFraction->SetBinContent( i, 1. );
        }
    }
    
    // calculate chi2
    if( iClock > 0. )
    {
        for( int i = 1; i <= hScalarDeadTimeFraction->GetNbinsX(); i++ )
        {
            if( hScalarDeadTimeFraction->GetBinError( i ) > 0. && iClock > 0. )
            {
                iChi2 += ( hScalarDeadTimeFraction->GetBinContent( i ) - iBusy / iClock ) * ( hScalarDeadTimeFraction->GetBinContent( i ) - iBusy / iClock )
                         / hScalarDeadTimeFraction->GetBinError( i ) / hScalarDeadTimeFraction->GetBinError( i );
                iNDF++;
            }
        }
    }
    if( iNDF > 0. )
    {
        fScalarDeadTimeChi2 = iChi2 / iNDF;
    }
    else
    {
        fScalarDeadTimeChi2 = 0.;
    }
    
    
    if( iClock > 0. )
    {
        return iBusy / iClock;
    }
    
    return 0.;
}

double VDeadTime::calculateDeadTimeFromTimeDifferences()
{
    ///////////////////////////////////////////////////
    // dead time calculation from time differences
    if( !hTimeDiff || !hFTimeDiff )
    {
        return 0.;
    }
    if( hTimeDiff->GetEntries() <= 0. )
    {
        return 0.;
    }
    
    // assume exponential distributions of dead times
    TF1 fFit( "fFit", "expo", fTFitMin, fTFitMax );
    hTimeDiff->Fit( &fFit, "Q0R" );
    hFTimeDiff->SetParameter( 0, fFit.GetParameter( 0 ) );
    hFTimeDiff->SetParameter( 1, fFit.GetParameter( 1 ) );
    
    double ix = hTimeDiff->GetBinCenter( 1 );
    double nmiss = 0.;
    
    // go left in delta t histogram from 0.01, first 0 bin defines dead time
    // (fails for very short runs)
    for( int i = hTimeDiff->FindBin( hTimeDiff->GetMean() ); i > 2; i-- )
    {
        // require  at least three zero bins
        if( hTimeDiff->GetBinContent( i ) == 0 && hTimeDiff->GetBinContent( i - 1 ) == 0 && hTimeDiff->GetBinContent( i - 2 ) == 0 )
        {
            ix = hTimeDiff->GetBinCenter( i );
            break;
        }
    }
    
    fDeadTimeMiss = nmiss;
    fDeadTimeMS = ix * 1000.;
    fDeadTimeFrac = 1. - TMath::Power( TMath::E(), ix * hFTimeDiff->GetParameter( 1 ) );
    if( fDeadTimeFrac < 1.e-5 )
    {
        fDeadTimeFrac = 0.;
    }
    
    // get dead time dependent on time diff
    
    int i_np = 0;
    for( int i = 1; i <= hTimeDiff2D->GetNbinsX(); i++ )
    {
        TH1D* h = hTimeDiff2D->ProjectionY( "h_tempXDGS", i, i, "e" );
        if( h && h->GetEntries() > 1000 )
        {
            h->Fit( &fFit, "Q0R" );
            ix = h->GetBinCenter( 1 );
            for( int b = h->FindBin( 0.005 ); b > 0; b-- )
            {
                if( h->GetBinContent( b ) == 0 )
                {
                    ix = h->GetBinCenter( b );
                    if( h->GetBinContent( b ) < 5 )
                    {
                        ix = h->GetBinCenter( b + 1 );
                    }
                    break;
                }
            }
            hgDeadTime->SetPoint( i - 1, hTimeDiff2D->GetXaxis()->GetBinCenter( i ), ( 1. - TMath::Power( TMath::E(), ix * fFit.GetParameter( 1 ) ) ) * 100. );
            double iDE = 0.;
            iDE += ix * TMath::Power( TMath::E(), ix * fFit.GetParameter( 1 ) ) * fFit.GetParError( 1 ) * ix * TMath::Power( TMath::E(), ix * fFit.GetParameter( 1 ) ) * fFit.GetParError( 1 );
            iDE += fFit.GetParameter( 1 ) * TMath::Power( TMath::E(), ix * fFit.GetParameter( 1 ) ) * h->GetBinWidth( 1 ) * fFit.GetParameter( 1 ) * TMath::Power( TMath::E(), ix * fFit.GetParameter( 1 ) ) * h->GetBinWidth( 1 ) / 4.;
            iDE = sqrt( iDE );
            hgDeadTime->SetPointError( i - 1, 0., iDE * 100. );
            
            i_np++;
            
            delete h;
        }
    }
    hgDeadTime->Set( i_np );
    
    return fDeadTimeFrac;
}

bool VDeadTime::checkStatus()
{
    if( fDeadTimeFrac_status == 0 )
    {
        cout << "Warning: no consistent results between time difference and scalar-based dead time calculation; falling back to time differences";
        return false;
    }
    return true;
}


void VDeadTime::printDeadTime()
{
    cout << "\t Tdiff: dead time [ms] " << fDeadTimeMS << ", fraction of missing events: " << fDeadTimeFrac * 100.;
    cout << "% (" << fDeadTimeMiss << ", ";
    if( hTimeDiff )
    {
        cout << hTimeDiff->GetEntries() << ")";
    }
    else
    {
        cout << ")";
    }
    cout << endl;
    checkStatus();
    
    if( hScalarClock && hScalarClock->GetEntries() > 0 )
    {
        cout << "\t Scalar: dead time (fraction of missing events): " << fScalarDeadTimeFrac * 100. << "%";
        cout << " (Chi2 = " << fScalarDeadTimeChi2 << ")" << endl;
    }
}


TList* VDeadTime::getDeadTimeHistograms()
{
    if( hisList )
    {
        return hisList;
    }
    
    return 0;
}


/*

   read time dependent dead time fraction

   double iT_run_s: time into the run
*/
double VDeadTime::getDeadTimeFraction( double iT_run_s, bool iTimeDiff, bool iCheckForConsistentDeadTime )
{
    double iDeadTime = 0.;
    // method of time difference: no time dependence (dead times might change during a run)
	if( iTimeDiff )
    {
        iDeadTime = fDeadTimeFrac;
    }
    // dead time fraction from scalars
    else
    {
        if( hScalarDeadTimeFraction )
        {
            int nbin = hScalarDeadTimeFraction->FindBin( iT_run_s );
            if( nbin > 0 && nbin <= hScalarDeadTimeFraction->GetNbinsX() )
            {
                iDeadTime = hScalarDeadTimeFraction->GetBinContent( nbin );
            }
        }
        else
        {
                 iDeadTime = fScalarDeadTimeFrac;
        }
	}
        // scalar dead time does not provide consistent values for a very small
        // number of runs
        if( iCheckForConsistentDeadTime && iDeadTime > 0.98 )
        {
            return fDeadTimeFrac;
        }
	return iDeadTime;
}

/*

   get mean dead time fraction assuming that each mask entry corresponds to 1 s

*/
double VDeadTime::getDeadTimeFraction( vector< bool > iMask, bool iTimeDiff, bool iCheckForConsistentDeadTime )
{
    double iN = 0.;
    double iD = 0.;
    for( unsigned int i = 0; i < iMask.size(); i++ )
    {
        // only get dead time if mask is open!
        if( iMask[i] )
        {
            iD += getDeadTimeFraction( ( double )i + 0.5, iTimeDiff, iCheckForConsistentDeadTime );
            iN++;
        }
    }
    if( iN > 0. )
    {
        return iD / iN;
    }
    
    return 0.;
}


void VDeadTime::writeHistograms( bool iDebug_IO )
{
    TDirectory* iDir = gDirectory;
    
    TDirectory* wDir = 0;
    
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "deadTimeHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "deadTimeHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    
    if( hisList )
    {
        int i_nbytes = hisList->Write();
        if( iDebug_IO )
        {
            cout << "WRITEDEBUG: dead time histograms (nbytes " << i_nbytes << ")" << endl;
        }
    }
    
    // remove all objects created with new in this class
    hisList->Delete();
    delete hisList;
    hisList = 0;
    
    iDir->cd();
}


bool VDeadTime::readHistograms( TDirectoryFile* iDir )
{
    if( !iDir )
    {
        return false;
    }
    
    if( hisList )
    {
        hisList->Delete();
    }
    else
    {
        hisList = new TList();
    }
    
    cout << "\t reading dead time histograms from file " << iDir->GetPath() << endl;
    
    // time difference histograms
    
    hTimeDiff = ( TH1D* )iDir->Get( "hTimeDiff_on" );
    if( hTimeDiff )
    {
        hisList->Add( hTimeDiff );
    }
    hTimeDiffLog = ( TH1D* )iDir->Get( "hTimeDiffLog_on" );
    if( hTimeDiffLog )
    {
        hisList->Add( hTimeDiffLog );
    }
    hTimeDiff2D = ( TH2D* )iDir->Get( "hTimeDiff2D_on" );
    if( hTimeDiff2D )
    {
        hisList->Add( hTimeDiff2D );
    }
    hgDeadTime = ( TGraphErrors* )iDir->Get( "hgDeadTime_on" );
    if( hgDeadTime )
    {
        hisList->Add( hgDeadTime );
    }
    hNEventTime = ( TH1D* )iDir->Get( "hNEventTime_on" );
    if( hNEventTime )
    {
        hisList->Add( hNEventTime );
    }
    hFTimeDiff = ( TF1* )iDir->Get( "hFTimeDiff_on" );
    if( hFTimeDiff )
    {
        hisList->Add( hFTimeDiff );
    }
    
    // scalar histograms
    hScalarClock = ( TH1D* )iDir->Get( "hScalarClock_on" );
    if( hScalarClock )
    {
        hisList->Add( hScalarClock );
    }
    hScalarBusy = ( TH1D* )iDir->Get( "hScalarBusy_on" );
    if( hScalarBusy )
    {
        hisList->Add( hScalarBusy );
    }
    hScalarDeadTimeFraction = ( TH1D* )iDir->Get( "hScalarDeadTimeFraction_on" );
    if( hScalarDeadTimeFraction )
    {
        hisList->Add( hScalarDeadTimeFraction );
    }
    
    return true;
}

