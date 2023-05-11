/*! \class VRunSummary
 *  \brief data class for rate/significance results from each run
 *
 *
 */

#include "VRunSummary.h"

VRunSummary::VRunSummary()
{
    fRunSummaryTree = new TTree( "tRunSummary", "anasum results" );
    
    setBranches();
}


bool VRunSummary::setBranches()
{
    if( !fRunSummaryTree )
    {
        return false;
    }
    
    fRunSummaryTree->Branch( "runOn", &runOn, "runOn/I" );
    fRunSummaryTree->Branch( "MJDOn", &MJDOn, "MJDOn/D" );
    fRunSummaryTree->Branch( "MJDOn_runStart", &MJDOn_runStart, "MJDOn_runStart/D" );
    fRunSummaryTree->Branch( "MJDOn_runStopp", &MJDOn_runStopp, "MJDOn_runStopp/D" );
    fRunSummaryTree->Branch( "RunDurationOn", &RunDurationOn, "RunDurationOn/D" );
    
    fRunSummaryTree->Branch( "runOff", &runOff, "runOff/I" );
    fRunSummaryTree->Branch( "MJDOff", &MJDOff, "MJDOff/D" );
    fRunSummaryTree->Branch( "MJDOff_runStart", &MJDOff_runStart, "MJDOff_runStart/D" );
    fRunSummaryTree->Branch( "MJDOff_runStopp", &MJDOff_runStopp, "MJDOff_runStopp/D" );
    fRunSummaryTree->Branch( "RunDurationOff", &RunDurationOff, "RunDurationOff/D" );
    
    fRunSummaryTree->Branch( "TargetName", &fTargetName, "TargetName/C" );
    fRunSummaryTree->Branch( "TargetRA", &fTargetRA, "TargetRA/D" );
    fRunSummaryTree->Branch( "TargetDec", &fTargetDec, "TargetDec/D" );
    fRunSummaryTree->Branch( "TargetRAJ2000", &fTargetRAJ2000, "TargetRAJ2000/D" );
    fRunSummaryTree->Branch( "TargetDecJ2000", &fTargetDecJ2000, "TargetDecJ2000/D" );
    fRunSummaryTree->Branch( "SkyMapCentreRAJ2000", &fSkyMapCentreRAJ2000, "SkyMapCentreRAJ2000/D" );
    fRunSummaryTree->Branch( "SkyMapCentreDecJ2000", &fSkyMapCentreDecJ2000, "SkyMapCentreDecJ2000/D" );
    fRunSummaryTree->Branch( "TargetShiftRAJ2000", &fTargetShiftRAJ2000, "TargetShiftRAJ2000/D" );
    fRunSummaryTree->Branch( "TargetShiftDecJ2000", &fTargetShiftDecJ2000, "TargetShiftDecJ2000/D" );
    fRunSummaryTree->Branch( "TargetShiftWest", &fTargetShiftWest, "TargetShiftWest/D" );
    fRunSummaryTree->Branch( "TargetShiftNorth", &fTargetShiftNorth, "TargetShiftNorth/D" );
    fRunSummaryTree->Branch( "WobbleNorth", &fWobbleNorth, "WobbleNorth/D" );
    fRunSummaryTree->Branch( "WobbleWest", &fWobbleWest, "WobbleWest/D" );
    fRunSummaryTree->Branch( "NTel", &fNTel, "NTel/i" );
    fRunSummaryTree->Branch( "TelList", &fTelList, "TelList/C" );
    fRunSummaryTree->Branch( "tOn", &tOn, "tOn/D" );
    fRunSummaryTree->Branch( "tOff", &tOff, "tOff/D" );
    fRunSummaryTree->Branch( "elevationOn", &elevationOn, "elevationOn/D" );
    fRunSummaryTree->Branch( "azimuthOn", &azimuthOn, "azimuthOn/D" );
    fRunSummaryTree->Branch( "elevationOff", &elevationOff, "elevationOff/D" );
    fRunSummaryTree->Branch( "azimuthOff", &azimuthOff, "azimuthOff/D" );
    fRunSummaryTree->Branch( "Theta2Max", &fTheta2Max, "Theta2Max/D" );
    fRunSummaryTree->Branch( "RawRateOn", &RawRateOn, "RawRateOn/D" );
    fRunSummaryTree->Branch( "RawRateOff", &RawRateOff, "RawRateOff/D" );
    fRunSummaryTree->Branch( "pedvarsOn", &pedvarsOn, "pedvarsOn/D" );
    fRunSummaryTree->Branch( "pedvarsOff", &pedvarsOff, "pedvarsOff/D" );
    fRunSummaryTree->Branch( "NOn", &NOn, "NOn/D" );
    fRunSummaryTree->Branch( "NOff", &NOff, "NOff/D" );
    fRunSummaryTree->Branch( "NOffNorm", &NOffNorm, "NOffNorm/D" );
    fRunSummaryTree->Branch( "OffNorm", &OffNorm, "OffNorm/D" );
    fRunSummaryTree->Branch( "Signi", &Signi, "Signi/D" );
    fRunSummaryTree->Branch( "Rate", &Rate, "Rate/D" );
    fRunSummaryTree->Branch( "RateE", &RateE, "RateE/D" );
    fRunSummaryTree->Branch( "RateOff", &RateOff, "RateOff/D" );
    fRunSummaryTree->Branch( "RateOffE", &RateOffE, "RateOffE/D" );
    fRunSummaryTree->Branch( "DeadTimeFracOn", &DeadTimeFracOn, "DeadTimeFracOn/D" );
    fRunSummaryTree->Branch( "DeadTimeFracOff", &DeadTimeFracOff, "DeadTimeFracOff/D" );
    fRunSummaryTree->Branch( "MaxSigni", &MaxSigni, "MaxSigni/D" );
    fRunSummaryTree->Branch( "MaxSigniX", &MaxSigniX, "MaxSigniX/D" );
    fRunSummaryTree->Branch( "MaxSigniY", &MaxSigniY, "MaxSigniY/D" );
    
    init();
    
    return true;
}


void VRunSummary::init()
{
    runOn = 0;
    runOff = 0;
    MJDOn = 0.;
    MJDOn_runStart = 0.;
    MJDOn_runStopp = 0.;
    RunDurationOn = 0.;
    MJDOff = 0.;
    MJDOff_runStart = 0.;
    MJDOff_runStopp = 0.;
    RunDurationOff = 0.;
    sprintf( fTargetName, "NOTSET" );
    sprintf( fTelList, "NOTSET" );
    fTargetDec = 0.;
    fTargetRA = 0.;
    fTargetDecJ2000 = 0.;
    fTargetRAJ2000 = 0.;
    fSkyMapCentreRAJ2000 = 0.;
    fSkyMapCentreDecJ2000 = 0.;
    fTargetShiftRAJ2000 = 0.;
    fTargetShiftDecJ2000 = 0.;
    fTargetShiftWest = 0.;
    fTargetShiftNorth = 0.;
    fWobbleNorth = 0.;
    fWobbleWest = 0.;
    fNTel = 0;
    tOn = 0.;
    tOff = 0.;
    elevationOn = 0.;
    elevationOff = 0.;
    azimuthOn  = 0.;
    azimuthOff = 0.;
    fTheta2Max = 0.;
    RawRateOn = 0.;
    RawRateOff = 0.;
    pedvarsOn = 0.;
    pedvarsOff = 0.;
    NOn = 0.;
    NOff = 0.;
    NOffNorm = 0.;
    OffNorm = 0.;
    Signi = 0.;
    Rate = 0.;
    RateE = 0.;
    RateOff = 0.;
    RateOffE = 0.;
    DeadTimeFracOn = 0.;
    DeadTimeFracOff = 0.;
    MaxSigni = 0.;
    MaxSigniX = 0.;
    MaxSigniY = 0.;
    
    fRunMJD.clear();
    fTotalExposureOn = 0.;
    fTotalExposureOff = 0.;
    f_exposureOn.clear();
    f_exposureOff.clear();
    fMeanElevationOn = 0.;
    fMeanElevationOff = 0.;
    fNMeanElevation = 0.;
    fMeanAzimuthOn = 0.;
    fMeanAzimuthOff = 0.;
    fMeanDeadTimeOn = 0.;
    fMeanDeadTimeOff = 0.;
    fMeanRawRateOn = 0.;
    fMeanRawRateOff = 0.;
    fMeanPedVarsOn = 0.;
    fMeanPedVarsOff = 0.;
    fTotTargetRA = 0.;
    fTotTargetDec = 0.;
    fTotTargetRAJ2000 = 0.;
    fTotTargetDecJ2000 = 0.;
    
}


void VRunSummary::fill()
{
    initTree();
    fRunSummaryTree->Fill();
}


void VRunSummary::write()
{
    if( fRunSummaryTree )
    {
        fRunSummaryTree->Write();
    }
}


void VRunSummary::print()
{
    initTree();
    char itemp[200];
    
    cout << endl << endl;
    cout << "RUN SUMMARY: " << endl << endl;
    
    for( int i = 0; i < fRunSummaryTree->GetEntries(); i++ )
    {
        fRunSummaryTree->GetEntry( i );
        
        if( runOn == -1 )
        {
            cout << endl;
            cout << "-----------------------------------------------------------------------------(6)" << endl;
            cout << "ALL RUNS (" << fRunSummaryTree->GetEntries() - 1 << " runs)";
        }
        else
        {
            //cout << "RUN " << runOn;
            printf( "RUN %5d", runOn ) ;
            if( runOff != runOn )
            {
                cout << " (" << runOff << ")       ";
            }
        }
        sprintf( itemp, " (%+4.1fN,%+4.1fW) ", fWobbleNorth, fWobbleWest );
        if( fWobbleNorth != 0. || fWobbleWest != 0. )
        {
            cout << itemp;
        }
        if( elevationOn > 2. )
        {
            //cout << " at " << (int)elevationOn << "," << (int)azimuthOn << " deg El.,Az, ";
            printf( " at %2d,%4d deg El., Az, ", ( int )elevationOn, ( int )azimuthOn ) ;
        }
        else
        {
            //cout << " at " << (int)elevationOff << "," << (int)azimuthOff << " deg El.,Az, ";
            printf( " at %2d,%4d deg El., Az, ", ( int )elevationOff, ( int )azimuthOff ) ;
        }
        sprintf( itemp, "%5.2f", tOn / 60. );
        cout << itemp << " min, ";
        sprintf( itemp, "%4d", ( int )NOn );
        cout << "Non: " << itemp;
        sprintf( itemp, "%5.2f", NOffNorm );
        cout << ", Noff: " <<  itemp;
        sprintf( itemp, "%4d", ( int )NOff );
        cout << " (" << itemp << ", norm ";
        sprintf( itemp, "%5.3f", OffNorm );
        cout << itemp << ")";
        sprintf( itemp, "%5.1f", Signi );
        cout << ", " <<  itemp << " sigma";
        sprintf( itemp, "%7.3f", Rate );
        cout << ", Rates: " << itemp;
        sprintf( itemp, "%7.3f", RateE );
        cout << " +/- " << itemp << " gamma/min";
        sprintf( itemp, "%7.3f", RateOff );
        cout << " (background: " << itemp << " events/min)" << endl;
    }
    
}

/*

   called for summary runs only

*/
bool VRunSummary::fill( string iDataDirectory,
                        string i_inputfile_total_directory,
                        vector< VAnaSumRunParameterDataClass > iRunList )
{
    char i_temp[2000];
    
    // current directory
    TDirectory* iCurrentDirectory = gDirectory;
    // copy relevant entries
    TChain* i_runSumChain = new TChain( "total_1/stereo/tRunSummary" );
    for( unsigned int i = 0; i < iRunList.size(); i++ )
    {
        sprintf( i_temp, "%s/%d.anasum.root", iDataDirectory.c_str(), iRunList[i].fRunOn );
        i_runSumChain->Add( i_temp );
    }
    fRunSummaryTree = i_runSumChain->CopyTree( "runOn>0" );
    fRunSummaryTree->SetDirectory( iCurrentDirectory );
    fRunSummaryTree->AutoSave();
    
    // reset variables
    fRunMJD.clear();
    
    fTotalExposureOn = 0.;
    fTotalExposureOff = 0.;
    f_exposureOn.clear();
    f_exposureOff.clear();
    fMeanElevationOn = 0.;
    fMeanElevationOff = 0.;
    fNMeanElevation = 0.;
    fMeanAzimuthOn = 0.;
    fMeanAzimuthOff = 0.;
    fMeanDeadTimeOn = 0.;
    fMeanDeadTimeOff = 0.;
    fMeanRawRateOn = 0.;
    fMeanRawRateOff = 0.;
    fMeanPedVarsOn = 0.;
    fMeanPedVarsOff = 0.;
    
    CRunSummary i_runSum( i_runSumChain );
    double iTargetRA = 0.;
    double iTargetDec = 0.;
    double iTargetRAJ2000 = 0.;
    double iTargetDecJ2000 = 0.;
    
    for( int n = 0; n < i_runSum.fChain->GetEntries(); n++ )
    {
        i_runSum.GetEntry( n );
        for( unsigned int i = 0; i < iRunList.size(); i++ )
        {
            if( i_runSum.runOn == iRunList[i].fRunOn )
            {
                runOn = iRunList[i].fRunOn;
                runOff = iRunList[i].fRunOff;
                
                // add values to run list
                fRunMJD[runOn] = i_runSum.MJDOn;
                fRunMJD[runOff] = i_runSum.MJDOff;
                fTotalExposureOn += i_runSum.tOn;
                fTotalExposureOff += i_runSum.tOff;
                f_exposureOn[runOn] = i_runSum.tOn;
                f_exposureOff[runOff] = i_runSum.tOff;
                fMeanElevationOn += i_runSum.elevationOn;
                fMeanElevationOff += i_runSum.elevationOff;
                fMeanAzimuthOn  = VSkyCoordinatesUtilities::addToMeanAzimuth( fMeanAzimuthOn,  i_runSum.azimuthOn );
                fMeanAzimuthOff = VSkyCoordinatesUtilities::addToMeanAzimuth( fMeanAzimuthOff, i_runSum.azimuthOff );
                fMeanDeadTimeOn += i_runSum.DeadTimeFracOn * tOn;
                fMeanDeadTimeOff += i_runSum.DeadTimeFracOff * tOff;
                fMeanRawRateOn += i_runSum.RawRateOn;
                fMeanRawRateOff += i_runSum.RawRateOff;
                fMeanPedVarsOn += i_runSum.pedvarsOn;
                fMeanPedVarsOff += i_runSum.pedvarsOff;
                fNMeanElevation++;
                
                if( fNMeanElevation == 0. )
                {
                    iTargetRA = i_runSum.TargetRA;
                    iTargetDec = i_runSum.TargetDec;
                    iTargetRAJ2000 = i_runSum.TargetRAJ2000;
                    iTargetDecJ2000 = i_runSum.TargetDecJ2000;
                }
                else if( iTargetRA != fTargetRA || iTargetDec != fTargetDec )
                {
                    iTargetRA = 0.;
                    iTargetDec = 0.;
                    iTargetRAJ2000 = 0.;
                    iTargetDecJ2000 = 0.;
                }
                break;
            }
        }
    }
    // set target coordinates
    fTotTargetRA = iTargetRA;
    fTotTargetDec = iTargetDec;
    fTotTargetRAJ2000 = iTargetRAJ2000;
    fTotTargetDecJ2000 = iTargetDecJ2000;
    
    iCurrentDirectory->cd();
    
    return true;
}


bool VRunSummary::initTree()
{
    if( !fRunSummaryTree )
    {
        return false;
    }
    
    fRunSummaryTree->SetBranchAddress( "runOn", &runOn );
    fRunSummaryTree->SetBranchAddress( "runOff", &runOff );
    fRunSummaryTree->SetBranchAddress( "MJDOn", &MJDOn );
    fRunSummaryTree->SetBranchAddress( "MJDOn_runStart", &MJDOn_runStart );
    fRunSummaryTree->SetBranchAddress( "MJDOn_runStopp", &MJDOn_runStopp );
    fRunSummaryTree->SetBranchAddress( "RunDurationOn", &RunDurationOn );
    fRunSummaryTree->SetBranchAddress( "MJDOff", &MJDOff );
    fRunSummaryTree->SetBranchAddress( "MJDOff_runStart", &MJDOff_runStart );
    fRunSummaryTree->SetBranchAddress( "MJDOff_runStopp", &MJDOff_runStopp );
    fRunSummaryTree->SetBranchAddress( "RunDurationOff", &RunDurationOff );
    fRunSummaryTree->SetBranchAddress( "TargetName", &fTargetName );
    fRunSummaryTree->SetBranchAddress( "TargetRA", &fTargetRA );
    fRunSummaryTree->SetBranchAddress( "TargetDec", &fTargetDec );
    fRunSummaryTree->SetBranchAddress( "TargetRAJ2000", &fTargetRAJ2000 );
    fRunSummaryTree->SetBranchAddress( "TargetDecJ2000", &fTargetDecJ2000 );
    fRunSummaryTree->SetBranchAddress( "SkyMapCentreRAJ2000", &fSkyMapCentreRAJ2000 );
    fRunSummaryTree->SetBranchAddress( "SkyMapCentreDecJ2000", &fSkyMapCentreDecJ2000 );
    fRunSummaryTree->SetBranchAddress( "TargetShiftRAJ2000", &fTargetShiftRAJ2000 );
    fRunSummaryTree->SetBranchAddress( "TargetShiftDecJ2000", &fTargetShiftDecJ2000 );
    fRunSummaryTree->SetBranchAddress( "TargetShiftWest", &fTargetShiftWest );
    fRunSummaryTree->SetBranchAddress( "TargetShiftNorth", &fTargetShiftNorth );
    fRunSummaryTree->SetBranchAddress( "WobbleNorth", &fWobbleNorth );
    fRunSummaryTree->SetBranchAddress( "WobbleWest", &fWobbleWest );
    fRunSummaryTree->SetBranchAddress( "NTel", &fNTel );
    fRunSummaryTree->SetBranchAddress( "TelList", &fTelList );
    fRunSummaryTree->SetBranchAddress( "tOn", &tOn );
    fRunSummaryTree->SetBranchAddress( "tOff", &tOff );
    fRunSummaryTree->SetBranchAddress( "elevationOn", &elevationOn );
    fRunSummaryTree->SetBranchAddress( "elevationOff", &elevationOff );
    fRunSummaryTree->SetBranchAddress( "azimuthOn", &azimuthOn );
    fRunSummaryTree->SetBranchAddress( "azimuthOff", &azimuthOff );
    fRunSummaryTree->SetBranchAddress( "Theta2Max", &fTheta2Max );
    fRunSummaryTree->SetBranchAddress( "RawRateOn", &RawRateOn );
    fRunSummaryTree->SetBranchAddress( "RawRateOff", &RawRateOff );
    fRunSummaryTree->SetBranchAddress( "pedvarsOn", &pedvarsOn );
    fRunSummaryTree->SetBranchAddress( "pedvarsOff", &pedvarsOff );
    fRunSummaryTree->SetBranchAddress( "NOn", &NOn );
    fRunSummaryTree->SetBranchAddress( "NOff", &NOff );
    fRunSummaryTree->SetBranchAddress( "NOffNorm", &NOffNorm );
    fRunSummaryTree->SetBranchAddress( "OffNorm", &OffNorm );
    fRunSummaryTree->SetBranchAddress( "Signi", &Signi );
    fRunSummaryTree->SetBranchAddress( "Rate", &Rate );
    fRunSummaryTree->SetBranchAddress( "RateE", &RateE );
    fRunSummaryTree->SetBranchAddress( "RateOff", &RateOff );
    fRunSummaryTree->SetBranchAddress( "RateOffE", &RateOffE );
    fRunSummaryTree->SetBranchAddress( "DeadTimeFracOn", &DeadTimeFracOn );
    fRunSummaryTree->SetBranchAddress( "DeadTimeFracOff", &DeadTimeFracOff );
    fRunSummaryTree->SetBranchAddress( "MaxSigni", &MaxSigni );
    fRunSummaryTree->SetBranchAddress( "MaxSigniX", &MaxSigniX );
    fRunSummaryTree->SetBranchAddress( "MaxSigniY", &MaxSigniY );
    
    return true;
}
