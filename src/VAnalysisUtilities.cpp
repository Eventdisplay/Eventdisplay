/*! \class VAnalysisUtilities
    \brief basic utilities

*/

#include "VAnalysisUtilities.h"

VAnalysisUtilities::VAnalysisUtilities()
{
    fDebug = 0;
    
    bZombie = false;
    
    fAnasumDataFile = 0;
    
    fSkyMapCentreDecJ2000 = 0.;
    fSkyMapCentreRAJ2000 = 0.;
    fTargetShiftWest = 0.;
    fTargetShiftNorth = 0.;
    fTargetDec = 0.;
    fTargetRA = 0.;
    fTargetDecJ2000 = 0.;
    fTargetRAJ2000 = 0.;
    
    fRunList_MJD_min = 0.;
    fRunList_MJD_max = 0.;
    
    fAnasumDataFile = 0;
    fEVNDISPVersion = "noVersionSet";
    
    setRunListMJDRange();
    setPhaseFoldingValues();
    
    setRunListCutMJDRange();
    setRunListCutPhaseRange();
}


bool VAnalysisUtilities::openFile( string iname, int irun, bool iStereo, bool iPrint )
{
    fAnasumDataFile = new TFile( iname.c_str() );
    if( fAnasumDataFile->IsZombie() )
    {
        cout << "error opening anasum file " << iname << endl;
        bZombie = true;
        return false;
    }
    VGlobalRunParameter* iPar = ( VGlobalRunParameter* )fAnasumDataFile->Get( "anasumRunParameter" );
    if( iPar )
    {
        fEVNDISPVersion = iPar->getEVNDISP_VERSION();
    }
    
    
    if( iStereo )
    {
        char dname[200];
        bool b_cd = false;
        if( irun == 0 )
        {
            b_cd = fAnasumDataFile->cd( "total/stereo" );
        }
        else if( irun < 0 )
        {
            sprintf( dname, "total_%d/stereo", -1 * irun );
            b_cd = fAnasumDataFile->cd( dname );
        }
        else if( irun == 1 )
        {
            sprintf( dname, "total_%d/stereo", irun );
            b_cd = fAnasumDataFile->cd( dname );
        }
        else
        {
            sprintf( dname, "run_%d/stereo", irun );
            b_cd = fAnasumDataFile->cd( dname );
        }
        if( !b_cd )
        {
            cout << "directory not found for run " << irun << " in " << iname << endl;
            bZombie = true;
            return false;
        }
        TTree* t = ( TTree* )gDirectory->Get( "tRunSummary" );
        if( t )
        {
            readTargetCoordinatesFromtRunSummary( t, irun );
        }
    }
    else
    {
        char hdir[200];
        if( irun == 0 )
        {
            fAnasumDataFile->cd( "total" );
        }
        else if( irun < 0 )
        {
            sprintf( hdir, "total_%d", -1 * irun );
            if( !fAnasumDataFile->cd( hdir ) )
            {
                return 0;
            }
        }
        else
        {
            sprintf( hdir, "run_%d", irun );
            if( !fAnasumDataFile->cd( hdir ) )
            {
                return 0;
            }
        }
    }
    if( iPrint )
    {
        cout << "file open: " << iname << " (";
        if( irun > 0 )
        {
            cout << irun;
        }
        else
        {
            cout << "combined runs";
        }
        cout << ")" << endl;
    }
    
    return true;
}

bool VAnalysisUtilities::readTargetCoordinatesFromtRunSummary( TTree* t, int ion )
{
    if( !t )
    {
        return false;
    }
    
    int iRun;
    t->SetBranchAddress( "runOn", &iRun );
    t->SetBranchAddress( "TargetRA", &fTargetRA );
    t->SetBranchAddress( "TargetDec", &fTargetDec );
    t->SetBranchAddress( "TargetRAJ2000", &fTargetRAJ2000 );
    t->SetBranchAddress( "TargetDecJ2000", &fTargetDecJ2000 );
    if( t->GetBranch( "SkyMapCentreRAJ2000" ) )
    {
        t->SetBranchAddress( "SkyMapCentreRAJ2000", &fSkyMapCentreRAJ2000 );
    }
    else
    {
        fSkyMapCentreRAJ2000 = 0.;
    }
    if( t->GetBranch( "SkyMapCentreDecJ2000" ) )
    {
        t->SetBranchAddress( "SkyMapCentreDecJ2000", &fSkyMapCentreDecJ2000 );
    }
    else
    {
        fSkyMapCentreDecJ2000 = -99.;
    }
    if( t->GetBranch( "TargetShiftWest" ) )
    {
        t->SetBranchAddress( "TargetShiftWest", &fTargetShiftWest );
    }
    else
    {
        fTargetShiftWest = 0.;
    }
    if( t->GetBranch( "TargetShiftNorth" ) )
    {
        t->SetBranchAddress( "TargetShiftNorth", &fTargetShiftNorth );
    }
    else
    {
        fTargetShiftNorth = 0.;
    }
    
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        // order matters!
        if( fSkyMapCentreDecJ2000 < -90. || ( TMath::Abs( fSkyMapCentreRAJ2000 ) < 1.e-8 && TMath::Abs( fSkyMapCentreDecJ2000 ) < 1.e-8 ) )
        {
            if( TMath::Abs( fTargetRAJ2000 ) < 1.e-8  && TMath::Abs( fTargetDecJ2000 ) < 1.e-8 )
            {
                fSkyMapCentreDecJ2000 = fTargetDec;
            }
            else
            {
                fSkyMapCentreDecJ2000 = fTargetDecJ2000;
            }
        }
        if( TMath::Abs( fSkyMapCentreRAJ2000 ) < 1.e-8 )
        {
            if( TMath::Abs( fTargetRAJ2000 ) < 1.e-8  && TMath::Abs( fTargetDecJ2000 ) < 1.e-8 )
            {
                fSkyMapCentreRAJ2000 = fTargetRA;
            }
            else
            {
                fSkyMapCentreRAJ2000 = fTargetRAJ2000;
            }
        }
        if( iRun == ion || ion <= 0 )
        {
            break;
        }
    }
    
    return true;
}

bool VAnalysisUtilities::closeFile()
{
    if( fAnasumDataFile )
    {
        fAnasumDataFile->Close();
    }
    else
    {
        return false;
    }
    
    return true;
}


bool VAnalysisUtilities::readRunList()
{
    vector< int > irunlist;
    
    return readRunList( irunlist );
}

CRunSummary* VAnalysisUtilities::getRunSummaryTree( int iTot )
{
    if( !fAnasumDataFile )
    {
        return 0;
    }
    
    char hname[200];
    if( iTot > 0 )
    {
        sprintf( hname, "total_%d/stereo", iTot );
    }
    else if( iTot < 0 )
    {
        sprintf( hname, "total_%d/stereo", -1 * iTot );
    }
    else
    {
        sprintf( hname, "total/stereo" );
    }
    fAnasumDataFile->cd( hname );
    
    TTree* t = ( TTree* )gDirectory->Get( "tRunSummary" );
    if( !t )
    {
        cout << "Error: no tRunSummary tree found in " << fAnasumDataFile->GetName() << endl;
        return 0;
    }
    CRunSummary* c = new CRunSummary( t );
    
    return c;
}

/*
 * calculate the cumulative significance vs times
 *
*/
TGraph* VAnalysisUtilities::calcCumulativeSig( int iTot, bool iCalculateAlphaProperly )
{

    CRunSummary* ctRunSum = getRunSummaryTree( iTot );
    int nentries;
    if( !ctRunSum )
    {
        cout << "VAnalysisUtilities::calcCumulativeSig Error: Run Summary tree not found." << endl;
        nentries = 0;
    }
    else
    {
        nentries = ctRunSum->fChain->GetEntries();
        if( nentries < 1 )
        {
            cout << "VAnalysisUtilities::calcCumulativeSig Warning: No entries in run summary tree." << endl;
        }
    }
    
    
    double cum_Non = 0;
    double cum_Noff = 0;
    double cum_alpha = 0;
    double cum_alphaNorm = 0;
    double cum_t = 0;
    int cum_N = 0;
    TGraph* gCumSig = new TGraph( 0 );
    
    for( int i = 0; i < nentries; i++ )
    {
        ctRunSum->GetEntry( i );
        // exclude runs with 0 rate and no error on rate
        if( ctRunSum->runOn > 0 && !( TMath::Abs( ctRunSum->Rate ) < 1.e-5 && TMath::Abs( ctRunSum->RateE ) < 1.e-5 ) )
        {
        
            // cumulative significance vs time
            cum_Non += ctRunSum->NOn;
            cum_Noff += ctRunSum->NOff;
            cum_t += ctRunSum->tOn;
            cum_N++;
            
            if( iCalculateAlphaProperly )
            {
                //calculate average alpha weighted by N_i / (1+alpha_i). See Stephane's calculations.
                cum_alpha += ctRunSum->OffNorm * ( ctRunSum->NOn + ctRunSum->NOff ) / ( 1.0 + ctRunSum->OffNorm );
                cum_alphaNorm += ( ctRunSum->NOn + ctRunSum->NOff ) / ( 1.0 + ctRunSum->OffNorm );
            }
            else
            {
                //calculate average alpha weighted by Noff_i (old way, do not use).
                cum_alpha += ctRunSum->OffNorm * ctRunSum->NOff ;
                cum_alphaNorm += ctRunSum->NOff;
            }
            //protect against division by 0 (no on/off events).
            double alpha =  cum_alphaNorm != 0 ?  cum_alpha / cum_alphaNorm : 0;
            
            gCumSig->SetPoint( cum_N, cum_t / 60, VStatistics::calcSignificance( cum_Non, cum_Noff, alpha ) );
            
        }
    }
    gCumSig->SetTitle( "Cumulative Significance vs Time;Observing time [min], not dead time corrected;significance [#sigma]" );
    gCumSig->SetMarkerStyle( 3 );
    gCumSig->SetLineWidth( 2 );
    gCumSig->SetMarkerSize( 2 );
    return gCumSig;
}



bool VAnalysisUtilities::readRunList( vector< int > irunlist, int iTot )
{
    VRunList i_temp;
    
    // reset run list
    fRunList.clear();
    
    CRunSummary* c = getRunSummaryTree( iTot );
    if( !c )
    {
        return false;
    }
    
    fRunList_MJD_min = 1.e14;
    fRunList_MJD_max = 0.;
    
    bool bFoundRun = false;
    
    for( int i = 0; i < c->fChain->GetEntries(); i++ )
    {
        c->GetEntry( i );
        
        if( c->runOn > 0 )
        {
            // check if this run is in runlist
            bFoundRun = false;
            if( irunlist.size() > 0 )
            {
                for( unsigned int r = 0; r < irunlist.size(); r++ )
                {
                    if( c->runOn == irunlist[r] )
                    {
                        bFoundRun = true;
                        break;
                    }
                }
            }
            else
            {
                bFoundRun = true;
            }
            
            if( !bFoundRun )
            {
                continue;
            }
            
            i_temp.runnumber = c->runOn;
            i_temp.MJD = c->MJDOn;
            if( i_temp.runnumber > 0 &&  i_temp.MJD > fRunList_MJD_max )
            {
                fRunList_MJD_max = i_temp.MJD;
            }
            if( i_temp.runnumber > 0 &&  i_temp.MJD < fRunList_MJD_min )
            {
                fRunList_MJD_min = i_temp.MJD;
            }
            i_temp.tOn = c->tOn;
            i_temp.deadTimeFraction = 1. - c->DeadTimeFracOn;
            i_temp.NOn = c->NOn;
            i_temp.NOff = c->NOff;
            i_temp.OffNorm = c->OffNorm;
            i_temp.elevationOn = c->elevationOn;
            i_temp.elevationOff = c->elevationOff;
            i_temp.TargetRAJ2000 = c->TargetRAJ2000;
            i_temp.TargetDecJ2000 = c->TargetDecJ2000;
            i_temp.pedvarsOn = c->pedvarsOn;
            i_temp.alpha = c->OffNorm;
            
            if( fPhase_MJD0 > 0. && fPhase_Period_days > 0. )
            {
                i_temp.phase = ( i_temp.MJD - fPhase_MJD0 ) / fPhase_Period_days;
                i_temp.phase = i_temp.phase - TMath::Floor( i_temp.phase );
            }
            
            // apply run list cuts
            bool iCut = true;
            for( unsigned int i = 0; i < fRunListCut_MJD_min.size(); i++ )
            {
                if( fRunListCut_MJD_min[i] > 0 && i_temp.MJD > fRunListCut_MJD_min[i]
                        && fRunListCut_MJD_max[i] > 0 && i_temp.MJD < fRunListCut_MJD_max[i] )
                {
                    iCut = false;
                }
                else if( fRunListCut_MJD_min[i] < 1. && fRunListCut_MJD_max[i] < 1. )
                {
                    iCut = false;
                }
            }
            if( iCut )
            {
                continue;
            }
            
            iCut = true;
            for( unsigned int i = 0; i < fRunListCut_Phase_min.size(); i++ )
            {
                if( fRunListCut_Phase_min[i] > 0 && i_temp.phase > fRunListCut_Phase_min[i]
                        && fRunListCut_Phase_max[i] > 0 && i_temp.phase < fRunListCut_Phase_max[i] )
                {
                    iCut = false;
                }
                else if( fRunListCut_Phase_min[i] < 0. && fRunListCut_Phase_max[i] < 0. )
                {
                    iCut = false;
                }
            }
            if( iCut )
            {
                continue;
            }
            
            
            fRunList.push_back( i_temp );
        }
    }
    
    return true;
}

/*!

   get a histogram/function/graph/etc from an anasum output file

*/
TObject* VAnalysisUtilities::getHistogram( string hisname, int runnumber, string dirname, double iSlizeY )
{
    if( !fAnasumDataFile )
    {
        return 0;
    }
    
    char dx[600];
    if( runnumber > 1 )
    {
        sprintf( dx, "run_%d/stereo/%s", runnumber, dirname.c_str() );
    }
    else
    {
        if( runnumber == 0 )
        {
            sprintf( dx, "total/stereo/%s", dirname.c_str() );
        }
        else if( runnumber == 1 )
        {
            sprintf( dx, "total_%d/stereo/%s", runnumber, dirname.c_str() );
        }
        else
        {
            sprintf( dx, "total_%d/stereo/%s", -1 * runnumber, dirname.c_str() );
        }
    }
    
    fAnasumDataFile->cd( dx );
    TDirectory* iDir = gDirectory;
    if( !iDir )
    {
        return 0;
    }
    
    TObject* h = ( TObject* )iDir->Get( hisname.c_str() );
    
    if( h && iSlizeY < -9998. )
    {
        return h->Clone();
    }
    else if( h )
    {
        string iClassName = h->ClassName();
        if( iClassName.find( "TH2" ) != string::npos )
        {
            TH2* i_h2 = ( TH2* )h;
            string iN = hisname + "px";
            TH1* i_h = ( TH1* )i_h2->ProjectionX( iN.c_str(), i_h2->GetYaxis()->FindBin( iSlizeY ), i_h2->GetYaxis()->FindBin( iSlizeY ) );
            return i_h->Clone();
        }
    }
    
    
    return 0;
}

TChain* VAnalysisUtilities::getTreeWithSelectedEvents( string iFile, bool iOn )
{
    if( !fAnasumDataFile )
    {
        return 0;
    }
    
    char dname[200];
    char hname[2000];
    if( iOn )
    {
        sprintf( dname, "data_on" );
    }
    else
    {
        sprintf( dname, "data_off" );
    }
    
    TTree* t = 0;
    TChain* c = new TChain( dname );
    TDirectory* iDir = gDirectory;
    
    // get some numbers from the run summary tree
    fAnasumDataFile->cd( "total_1/stereo" );
    
    t = ( TTree* )iDir->Get( "tRunSummary" );
    if( t )
    {
        int iRun;
        t->SetBranchAddress( "runOn", &iRun );
        for( int i = 0; i < t->GetEntries(); i++ )
        {
            t->GetEntry( i );
            
            if( iRun != -1 )
            {
                sprintf( hname, "%s/run_%d/stereo/%s", iFile.c_str(), iRun, dname );
                c->Add( hname );
            }
        }
    }
    
    return c;
}



double VAnalysisUtilities::getNormalisationFactor( int iRun )
{

    if( !fAnasumDataFile )
    {
        return -1;
    }
    
    TDirectory* iDir = gDirectory;
    TTree* t = ( TTree* )iDir->Get( "tRunSummary" );
    
    if( !t )
    {
        return -1.;
    }
    
    CRunSummary* c = new CRunSummary( t );
    
    int nentries = c->fChain->GetEntries();
    
    for( int i = 0; i < nentries; i++ )
    {
        c->GetEntry( i );
        
        double ivalue = c->OffNorm;
        if( iRun == c->runOn )
        {
            closeFile();
            return ivalue;
        }
    }
    return -1.;
    
}

void VAnalysisUtilities::printEnergyThresholds()
{
    for( unsigned int i = 0; i < fRunList.size(); i++ )
    {
        cout << "RUN " << fRunList[i].runnumber;
        cout << "  Threshold [GeV]: " <<  fRunList[i].energyThreshold * 1.e3 << endl;
    }
}

void VAnalysisUtilities::printRunList()
{
    for( unsigned int i = 0; i < getRunList().size(); i++ )
    {
        getRunList()[i].print();
    }
}

void VAnalysisUtilities::printRunList( bool csv )
{
    if( csv )
        for( unsigned int i = 0; i < getRunList().size(); i++ )
        {
            getRunList()[i].print( csv );
        }
    else
        for( unsigned int i = 0; i < getRunList().size(); i++ )
        {
            getRunList()[i].print();
        }
}

void VAnalysisUtilities::setRunListCutMJDRange( double iMJDMin, double iMJDMax )
{
    fRunListCut_MJD_min.clear();
    fRunListCut_MJD_min.push_back( iMJDMin );
    
    fRunListCut_MJD_max.clear();
    fRunListCut_MJD_max.push_back( iMJDMax );
}

void VAnalysisUtilities::setRunListCutMJDRangeVector( vector< double > iMJDMinV, vector< double > iMJDMaxV )
{
    fRunListCut_MJD_min = iMJDMinV;
    fRunListCut_MJD_max = iMJDMaxV;
}

void VAnalysisUtilities::setRunListCutPhaseRange( double iPhaseMin, double iPhaseMax )
{
    fRunListCut_Phase_min.clear();
    fRunListCut_Phase_min.push_back( iPhaseMin );
    
    fRunListCut_Phase_max.clear();
    fRunListCut_Phase_max.push_back( iPhaseMax );
}

void VAnalysisUtilities::setRunListCutPhaseRangeVector( vector< double > iPhaseMinV, vector< double > iPhaseMaxV )
{
    fRunListCut_Phase_min = iPhaseMinV;
    fRunListCut_Phase_max = iPhaseMaxV;
}

vector< int > VAnalysisUtilities::getRunListVector()
{
    vector< int > r;
    
    CRunSummary* c = getRunSummaryTree( 1 );
    int nentries = c->fChain->GetEntries();
    
    for( int i = 0; i < nentries; i++ )
    {
        c->GetEntry( i );
        
        if( c->runOn > 0 )
        {
            r.push_back( c->runOn );
        }
    }
    return r;
}
