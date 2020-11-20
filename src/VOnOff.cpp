/*! \class VOnOff
 *  \brief do signal - background histogramming for sky maps and 1D histograms (e.g. energy spectra, mscw histograms, ...)
 *
 */

#include "VOnOff.h"

VOnOff::VOnOff()
{
    fDebug = false;
    
    hList = new TList();
    hQList = new TList();
    hSList = new TList();                         // sky histograms
    hPList = new TList();                         // stereo parameter histograms
    
    hListStereoParameterHistograms = new TList();
    hListRandomForestHistograms = new TList();
    hListEnergyHistograms = new TList();
    hListSkyHistograms = new TList();
    
    h1Dsig = 0;
    hmap_stereo_sig = 0;
    hTheta2_diff = 0;
    
    fMaxSigma = 0.;
    fMaxSigmaX = 0.;
    fMaxSigmaY = 0.;
}


VOnOff::~VOnOff()
{
    // remove all objects in hList
    if( hList )
    {
        TIter next( hList );
        while (TObject *obj = next())
        {
            if( obj && obj->TestBit( kCanDelete ) )
            {
               obj->Delete();
            }
        }
        delete hList;
    } 
    if( hList )
    {
        delete hList;
    }
    if( hQList )
    {
        delete hQList;
    }
    if( hSList )
    {
        delete hSList;
    }
    if( hPList )
    {
        delete hPList;
    }
    if( hListStereoParameterHistograms )
    {
        delete hListStereoParameterHistograms;
    }
    if( hListRandomForestHistograms )
    {
        delete hListRandomForestHistograms;
    }
    if( hListEnergyHistograms )
    {
        delete hListEnergyHistograms;
    }
    if( hListSkyHistograms )
    {
        delete hListSkyHistograms;
    }
}


void VOnOff::setTitles( TH1* his, string iname, string ititle, string ytitle )
{
    string itemp;
    // set names and titles
    itemp = his->GetName();
    itemp.replace( itemp.rfind( "on" ), 2, iname );
    his->SetName( itemp.c_str() );
    itemp = his->GetTitle();
    if( itemp.rfind( "(on)" ) < itemp.size() )
    {
        itemp.replace( itemp.rfind( "(on)" ), 4, ititle );
    }
    else
    {
        itemp += ititle;
    }
    his->SetTitle( itemp.c_str() );
    if( ytitle.size() > 0 )
    {
        his->GetYaxis()->SetTitle( ytitle.c_str() );
    }
}

/*
 *
 *
 */
void VOnOff::createHistograms( TList* ion, TList* il )
{
    string itemp;
    string ihname;
    TIter next_on( ion );
    while( TH1* hon = ( TH1* )next_on() )
    {
        itemp = hon->ClassName();
        ihname = hon->GetName();
        if( itemp == "TH1D" )
        {
            il->Add( new TH1D( *( ( TH1D* )hon ) ) );
        }
        else if( itemp == "TH2D" )
        {
            il->Add( new TH2D( *( ( TH2D* )hon ) ) );
        }
        else if( itemp == "TProfile" )
        {
            il->Add( new TProfile( *( ( TProfile* )hon ) ) );
        }
        else
        {
            cout << "VOnOff::createHistograms: error, unknown class name " << hon->GetName() << "\t" << itemp << endl;
        }
    }
}


void VOnOff::doOnOffforParameterHistograms( TList* iponlist, TList* ipofflist, double i_norm_alpha, bool isCombined )
{
    if( fDebug )
    {
        cout << "VOnOff::doOnOff" << endl;
    }
    string itemp;
    
    hPList->Clear();
    
    ///////////////////////////////////////////////////////////////
    // fill parameter histograms
    
    // create diff histograms
    createHistograms( iponlist, hPList );
    
    TIter nextp( hPList );
    TIter n_onp( iponlist );
    TIter n_offp( ipofflist );
    while( TH1* hTemp = ( TH1* )nextp() )
    {
        hTemp->Reset();
        setTitles( hTemp, "diff", " (ON-OFF)", "" );
        // get on/off histograms
        TH1* hon  = ( TH1* )n_onp();
        TH1* hoff = ( TH1* )n_offp();
        
        // calculate difference
        itemp = hon->GetName();
        // htheta2 histogram (note: calculated from one reflected region only!)
        if( itemp.find( "htheta2" ) == 0 )
        {
            hTemp->Add( hon, hoff, 1., -1. );
            hTheta2_diff = ( TH1D* )hTemp;
        }
        // energy histogram with x-axis in logE or linE
        else if( itemp.find( "herec" ) == 0 || itemp.find( "hLinerec" ) == 0 )
        {
            hTemp->Add( hon, hoff, 1., -1.*i_norm_alpha );
        }
        // all other histograms
        else
        {
            if( !isCombined )
            {
                hTemp->Add( hon, hoff, 1., -1.*i_norm_alpha );
                hoff->Scale( i_norm_alpha );
            }
            else
            {
                hTemp->Add( hon, hoff, 1., -1. );
            }
        }
        
        // fill the corresponding lists
        if( itemp.find( "herec" ) == 0 || itemp.find( "hLinerec" ) == 0 )
        {
            hListEnergyHistograms->Add( hTemp );
        }
        else if( itemp.find( "hrf" ) == 0 )
        {
            hListRandomForestHistograms->Add( hTemp );
        }
        else
        {
            hListStereoParameterHistograms->Add( hTemp );
        }
        hList->Add( hTemp );
    }
    hList->Add( hPList );
}

/*

    on - alpha * off for sky histograms

*/
void VOnOff::doOnOffforSkyHistograms( TList* ionlist, TList* iofflist, TH2D* ialpha )
{
    hSList->Clear();
    
    // create diff histograms
    createHistograms( ionlist, hSList );
    
    string itemp;
    TIter next( hSList );
    TIter n_on( ionlist );
    TIter n_off( iofflist );
    while( TH1* hTemp = ( TH1* )next() )
    {
        hTemp->Reset();
        setTitles( hTemp, "diff", " (ON-OFF)", "" );
        // get on/off histograms
        TH1* hon  = ( TH1* )n_on();
        TH1* hoff = ( TH1* )n_off();
        
        // calculate difference
        itemp = hTemp->ClassName();
        
        // do only 2D histograms
        if( itemp == "TH2D" )
        {
            for( int i = 1; i <= hTemp->GetNbinsX(); i++ )
            {
                int i_a = ialpha->GetXaxis()->FindBin( hTemp->GetXaxis()->GetBinCenter( i ) );
                for( int j = 1; j <= hTemp->GetNbinsY(); j++ )
                {
                    int j_a = ialpha->GetYaxis()->FindBin( hTemp->GetYaxis()->GetBinCenter( j ) );
                    // normalisation factor must be nonzero
                    if( ialpha->GetBinContent( i_a, j_a ) != 0. && hon->GetBinContent( i, j ) > 0. )
                    {
                        hTemp->SetBinContent( i, j, hon->GetBinContent( i, j ) - hoff->GetBinContent( i, j )*ialpha->GetBinContent( i_a, j_a ) );
                        float iE = hon->GetBinContent( i, j ) + ialpha->GetBinContent( i_a, j_a ) * ialpha->GetBinContent( i_a, j_a ) * hoff->GetBinContent( i, j );
                        if( iE > 0. )
                        {
                            hTemp->SetBinError( i, j, sqrt( iE ) );
                        }
                    }
                    else
                    {
                        hTemp->SetBinContent( i, j, -9999. );
                    }
                }
            }
        }
        else if( itemp != "TH1D" )
        {
            cout << "VOnOff::doOnOffforSkyHistograms: error, unknown histogram type " << hTemp->GetName() << "\t" << itemp << endl;
        }
        // add histogram to the corresponding list
        // (we don't want to use alpha histograms)
        itemp = hTemp->GetName();
        if( itemp.find( "alpha" ) == string::npos )
        {
            hList->Add( hTemp );
            hListSkyHistograms->Add( hTemp );
        }
    }
    
    if( fDebug )
    {
        cout << "\t VOnOff::doOnOff" << endl;
    }
}


void VOnOff::doQfactors( TList* ionlist, TList* iofflist, double i_norm )
{
    if( fDebug )
    {
        cout << "VOnOff::doQfactors" << endl;
    }
    hQList->Clear();
    
    double i_sumon;
    double i_sumoff;
    int i_nbins;
    
    TIter q_on( ionlist );
    TIter q_off( iofflist );
    string itemp;
    TH1D* hTemp = 0;
    
    while( TH1* hon = ( TH1* )q_on() )
    {
        TH1* hoff = ( TH1* )q_off();
        itemp = hon->ClassName();
        if( itemp == "TH1D" )
        {
            // qfactor on low bound cut
            hTemp = new TH1D( *( ( TH1D* )hon ) );
            hTemp->Reset();
            hQList->Add( hTemp );
            hList->Add( hTemp );
            setTitles( hTemp, "qlo", " (Q-Factor for Low-Bound Cut)", "significance" );
            
            i_sumon  = 0;
            i_sumoff = 0;
            i_nbins = hTemp->GetNbinsX();
            for( int i = i_nbins; i > 0; i-- )
            {
                i_sumon  += hon->GetBinContent( i );
                i_sumoff += hoff->GetBinContent( i );
                if( i_sumon + i_sumoff > 0 )
                {
                    hTemp->SetBinContent( i, VStatistics::calcSignificance( i_sumon, i_sumoff, i_norm ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
                hTemp->SetBinError( i, 0. );
            }
            
            // qfactor on high bound cut
            hTemp = new TH1D( *( ( TH1D* )hon ) );
            hTemp->Reset();
            hQList->Add( hTemp );
            hList->Add( hTemp );
            setTitles( hTemp, "qhi", " (Q-Factor for High-Bound Cut)", "significance" );
            
            i_sumon  = 0;
            i_sumoff = 0;
            i_nbins = hTemp->GetNbinsX();
            for( int i = 1; i <= i_nbins; i++ )
            {
                i_sumon  += hon->GetBinContent( i );
                i_sumoff += hoff->GetBinContent( i );
                if( i_sumon + i_sumoff > 0 )
                {
                    hTemp->SetBinContent( i, VStatistics::calcSignificance( i_sumon, i_sumoff, i_norm ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
                hTemp->SetBinError( i, 0. );
            }
        }
    }
    if( fDebug )
    {
        cout << "\t VOnOff::doQfactors" << endl;
    }
}

/*
 * fill Li & Ma significance into a 2D map
 *
 * calculated from on, off, and alpha maps
*/
TH2D* VOnOff::do2DSignificance( TH2D* ion, TH2D* ioff, TH2D* ialpha, string ititle )
{
    if( fDebug )
    {
        cout << "VOnOff::do2DSignificance" << endl;
    }
    hmap_stereo_sig = new TH2D( *ion );
    hmap_stereo_sig->Reset();
    hList->Add( hmap_stereo_sig );
    
    double i_sigmax = 0.;
    double i_sigon = 0.;
    double i_sigoff = 0.;
    double i_sigalpha = 0.;
    double i_sigx = 0.;
    double i_sigy = 0.;
    
    for( int j = 1; j <= hmap_stereo_sig->GetNbinsX(); j++ )
    {
        int j_a = ialpha->GetXaxis()->FindBin( hmap_stereo_sig->GetXaxis()->GetBinCenter( j ) );
        for( int k = 1; k <= hmap_stereo_sig->GetNbinsY(); k++ )
        {
            int k_a = ialpha->GetYaxis()->FindBin( hmap_stereo_sig->GetYaxis()->GetBinCenter( k ) );
            double on  = ion->GetBinContent( j, k );
            double off = ioff->GetBinContent( j, k );
            if( on > 0. && ialpha->GetBinContent( j_a, k_a ) > 0. )
            {
                hmap_stereo_sig->SetBinContent( j, k, VStatistics::calcSignificance( on, off, ialpha->GetBinContent( j_a, k_a ) ) );
            }
            else
            {
                hmap_stereo_sig->SetBinContent( j, k, -9999. );
            }
            // get maximum in 2D sky map
            if( hmap_stereo_sig->GetBinContent( j, k ) > i_sigmax )
            {
                i_sigmax = hmap_stereo_sig->GetBinContent( j, k );
                i_sigon = on;
                i_sigoff = off;
                i_sigalpha = ialpha->GetBinContent( j_a, k_a );
                i_sigx = hmap_stereo_sig->GetXaxis()->GetBinCenter( j );
                i_sigy = hmap_stereo_sig->GetYaxis()->GetBinCenter( k );
            }
        }
    }
    setTitles( hmap_stereo_sig, "sig", " (Significance ON-OFF)", "" );
    
    // print maximum significance to screen
    if( i_sigmax > 0. && ititle.size() < 1 )
    {
        cout << "\t " << setprecision( 4 ) << i_sigmax << " (On: " << i_sigon << ", Off: " << i_sigoff;
        cout << ", Alpha: " << i_sigalpha << ", Off/Alpha: " << i_sigoff* i_sigalpha << ")";
        cout << " at (x=" << i_sigx << ", y=" << i_sigy << ")" << endl;
        
        fMaxSigma = i_sigmax;
        fMaxSigmaX = i_sigx;
        fMaxSigmaY = i_sigy;
    }
    else if( ititle.size() < 1 )
    {
        cout << "\t Max stereo two-d significance is zero" << endl;
        fMaxSigma = 0.;
        fMaxSigmaX = 0.;
        fMaxSigmaY = 0.;
    }
    
    if( fDebug )
    {
        cout << "\t VOnOff::do2DSignificance" << endl;
    }
    return hmap_stereo_sig;
}


void VOnOff::writeHistograms( TH2D* hSig, TH2D* hSigUC, TH2D* hDiff, TH2D* hDiffUC )
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
    if( hListSkyHistograms )
    {
        hListSkyHistograms->Write();
    }
    if( hSig )
    {
        hSig->Write();
    }
    if( hSigUC )
    {
        hSigUC->Write();
    }
    if( hDiff )
    {
        hDiff->Write();
    }
    if( hDiffUC )
    {
        hDiffUC->Write();
    }
    
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
    if( hListStereoParameterHistograms )
    {
        hListStereoParameterHistograms->Write();
    }
    
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
    if( hListEnergyHistograms )
    {
        hListEnergyHistograms->Write();
    }
    
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
    if( hListRandomForestHistograms )
    {
        hListRandomForestHistograms->Write();
    }
    
    // write all quality factor histograms
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "qualityFactorHistograms" );
    if( !wDir )
    {
        iDir->mkdir( "qualityFactorHistograms" )->cd();
    }
    else
    {
        wDir->cd();
    }
    if( hQList )
    {
        hQList->Write();
    }
    
    iDir->cd();
}


/*!
 *   fill 1 dimensional distribution of significances from 2D sky maps
 *
 */
void VOnOff::fill1DSignificanceHistogram( double rmax )
{
    h1Dsig = new TH1D( "hsig1D", "1D significance distribution", 500, -10., 60. );
    h1Dsig->Sumw2();
    h1Dsig->SetXTitle( "significance" );
    h1Dsig->SetYTitle( "entries" );
    hList->Add( h1Dsig );
    
    if( !hmap_stereo_sig )
    {
        return;
    }
    
    double ir = 0.;
    for( int j = 0; j < hmap_stereo_sig->GetNbinsX(); j++ )
    {
        for( int k = 0; k < hmap_stereo_sig->GetNbinsY(); k++ )
        {
            ir  = hmap_stereo_sig->GetXaxis()->GetBinCenter( j ) * hmap_stereo_sig->GetXaxis()->GetBinCenter( j );
            ir += hmap_stereo_sig->GetYaxis()->GetBinCenter( k ) * hmap_stereo_sig->GetYaxis()->GetBinCenter( k );
            ir  = sqrt( ir );
            if( ir < rmax && hmap_stereo_sig->GetBinContent( j, k ) != 0. )
            {
                h1Dsig->Fill( hmap_stereo_sig->GetBinContent( j, k ) );
            }
        }
    }
}


TList* VOnOff::getEnergyHistograms()
{
    return hListEnergyHistograms;
}

double VOnOff::getMaxSigma()
{
    if( TMath::Abs( fMaxSigma ) < 1.e-10 )
    {
        return 0.;
    }
    
    return fMaxSigma;
}

double VOnOff::getMaxSigmaX()
{
    if( TMath::Abs( fMaxSigmaX ) < 1.e-10 )
    {
        return 0.;
    }
    
    return fMaxSigmaX;
}

double VOnOff::getMaxSigmaY()
{
    if( TMath::Abs( fMaxSigmaY ) < 1.e-10 )
    {
        return 0.;
    }
    
    return fMaxSigmaY;
}
