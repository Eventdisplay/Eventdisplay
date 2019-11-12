/*! \class VPedestalCalculator
    \brief pedestal calculation in time slices



*/

#include "VPedestalCalculator.h"

VPedestalCalculator::VPedestalCalculator()
{
    fDebug = getDebugFlag();
    
    // default parameters (can be adjusted later in initialize()
    fLengthofTimeSlice = 180.;                     // in [s]
    fSumWindow = 24;
    fNPixel = 500;
    fSumFirst = 0;
    
    bCalibrationRun = false;
}


bool VPedestalCalculator::initialize()
{
    return initialize( bCalibrationRun, fNPixel, fLengthofTimeSlice, fSumFirst, fSumWindow );
}


bool VPedestalCalculator::initialize( bool ibCalibrationRun, unsigned int iNPixel, double iLengthofTimeSlice,
                                      int iSumFirst, int iSumWindow, double iRunStartTime, double iRunStoppTime )
{
    if( fDebug )
    {
        cout << "VPedestalCalculator::initialize " << iNPixel << " " << iLengthofTimeSlice << " " << iSumFirst << " " << iSumWindow << endl;
    }
    bCalibrationRun = ibCalibrationRun;
    fLengthofTimeSlice = adjustTimeSliceLength( iLengthofTimeSlice, iRunStartTime, iRunStoppTime );
    fSumFirst = iSumFirst;
    fSumWindow = iSumWindow;
    fNPixel = iNPixel;
    cout << endl;
    cout << "VPedestalCalculator::initialize, adjusted length of time slice from " << iLengthofTimeSlice << "s to " << fLengthofTimeSlice << "s";
    cout << " (length of time slice: " << fLengthofTimeSlice << " [s], sum window starts at ";
    cout << fSumFirst << " with lengths up to " << fSumWindow << " samples)" << endl;
    // test if camera is not too big
    if( getDetectorGeo()->getNChannels()[getTelID()] > fNPixel )
    {
        cout << "=================================" << endl;
        cout << "VPedestalCalculator::initialize error: Telescope: " << getTelID() + 1 << endl;
        cout << "VPedestalCalculator::initialize error: maximum number of allowed channel is " << fNPixel << endl;
        cout << "VPedestalCalculator::initialize error: actual number of allowed channel is " << getDetectorGeo()->getNChannels()[getTelID()] << endl;
        cout << "VPedestalCalculator::initialize info: ignoring pedestal variance calculation" << endl;
        getRunParameter()->fPedestalsInTimeSlices = false;
        return false;
    }
    // test if summation window is not to big
    if( ( unsigned int )fSumWindow > VDST_MAXSUMWINDOW )
    {
        cout << "=================================" << endl;
        cout << "VPedestalCalculator::initialize error: Telescope: " << getTelID() + 1 << endl;
        cout << "VPedestalCalculator::initialize error: maximum size of summation window is: " << VDST_MAXSUMWINDOW << endl;
        cout << "VPedestalCalculator::initialize error: actual size of summation window is: " << fSumWindow << endl;
        cout << "VPedestalCalculator::initialize info: ignoring pedestal variance calculation" << endl;
        getRunParameter()->fPedestalsInTimeSlices = false;
        return false;
    }
    // initialize data readers
    if( !initializeDataReader() )
    {
        cout << "VPedestalCalculator::initialize, error: cannot initialize data readers" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // reset all variables
    reset();
    
    // set up the trees
    TDirectory* iDir = gDirectory;
    
    char hname[200];
    char htitle[200];
    
    vector< float > iped_cal;
    vector< vector< float > > iped_cal2;
    
    for( unsigned int p = 0; p < fNPixel; p++ )
    {
        iped_cal.clear();
        for( int w = 0; w < fSumWindow; w++ )
        {
            iped_cal.push_back( 0. );
        }
        v_temp_pedEntries.push_back( iped_cal );
        v_temp_ped.push_back( iped_cal );
        v_temp_pedvar.push_back( iped_cal );
    }
    
    // vectors
    vector< int > iv_vint;
    vector< float > iv_vdouble;
    vector< double > iv_vdoubleX;
    vector< vector< unsigned int > > iv_ivuint;
    vector< vector< vector< float > > > iv_vvdouble;
    
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        // define data vectors
        v_MJD.push_back( iv_vint );
        v_time.push_back( iv_vdoubleX );
        v_pedEntries.push_back( iv_vvdouble );
        v_ped.push_back( iv_vvdouble );
        v_pedvar.push_back( iv_vvdouble );
        
        setTelID( i );
        unsigned int t = getTeltoAna()[i];
        
        // define the output tree
        if( !bCalibrationRun )
        {
        
            getAnaDirectories()[t]->cd();
            
            sprintf( hname, "tPedVar" );
            sprintf( htitle, "pedestal variations (telescope %d)", t + 1 );
            fTree.push_back( new TTree( hname, htitle ) );
            if( !fTree.back() )
            {
                return false;
            }
            fTree[i]->SetAutoSave( 100000000 );
            fTree[i]->Branch( "runNumber", &runNumber, "runNumber/I" );
            fTree[i]->Branch( "time",  &time,  "time/D" );
            fTree[i]->Branch( "MJD",  &MJD,  "MJD/I" );
            sprintf( hname, "x[%d]/D", fNPixel );
            fTree[i]->Branch( "x", x, hname );
            sprintf( hname, "y[%d]/D", fNPixel );
            fTree[i]->Branch( "y", y, hname );
            sprintf( hname, "xRot[%d]/D", fNPixel );
            fTree[i]->Branch( "xRot", xRot, hname );
            sprintf( hname, "yRot[%d]/D", fNPixel );
            fTree[i]->Branch( "yRot", yRot, hname );
        }
        // fill x and y positions of pixels (unrotated, don't change during run)
        getDetectorGeo()->setTelID( t );
        for( unsigned int p = 0; p < getDetectorGeo()->getNChannels()[getTelID()]; p++ )
        {
            x[p] = getDetectorGeo()->getXUnrotated()[p];
            y[p] = getDetectorGeo()->getYUnrotated()[p];
        }
        // get run number
        runNumber = getRunNumber();
        
        // initialise the pedvars variables
        iped_cal2.clear();
        for( unsigned int p = 0; p < fNPixel; p++ )
        {
            iped_cal.clear();
            for( int w = 0; w < fSumWindow; w++ )
            {
                iped_cal.push_back( 0. );
            }
            iped_cal2.push_back( iped_cal );
        }
        fpedcal_n.push_back( iped_cal2 );
        fpedcal_mean.push_back( iped_cal2 );
        fpedcal_mean2.push_back( iped_cal2 );
        
        // define the time vector
        fTimeVec.push_back( 0 );
    }
    iDir->cd();
    
    if( fDebug )
    {
        cout << "END: VPedestalCalculator::initialize " << endl;
    }
    
    return true;
}


void VPedestalCalculator::fillTimeSlice( unsigned int telID )
{
    // loop over all channels
    for( unsigned int p = 0; p < fpedcal_mean[telID].size(); p++ )
    {
        // loop over all summation windows
        unsigned int iTempSW = fpedcal_mean[telID][p].size();
        if( getRunParameter()->fCalibrationSumWindow < ( int )iTempSW )
        {
            iTempSW = getRunParameter()->fCalibrationSumWindow;
        }
        for( unsigned int w = 0; w < iTempSW; w++ )
        {
            // get pedestal values
            if( fpedcal_n[telID][p][w] > 0. )
            {
                v_temp_pedEntries[p][w] = fpedcal_n[telID][p][w];
                v_temp_ped[p][w]        = fpedcal_mean[telID][p][w] / fpedcal_n[telID][p][w] / ( double )( w + 1 );
                // root use 1/n and not 1./(n-1) in rms calculation: for consistency use the same here (shouldn't matter)
                //	      v_temp_pedvar[p][w]     = sqrt( 1./(fpedcal_n[telID][p][w]-1.) * ( fpedcal_mean2[telID][p][w] - fpedcal_mean[telID][p][w]*fpedcal_mean[telID][p][w]/fpedcal_n[telID][p][w] ) );
                v_temp_pedvar[p][w]     = sqrt( 1. / ( fpedcal_n[telID][p][w] )
                                                * TMath::Abs( fpedcal_mean2[telID][p][w]
                                                        - fpedcal_mean[telID][p][w] * fpedcal_mean[telID][p][w] / fpedcal_n[telID][p][w] ) );
            }
            else
            {
                v_temp_pedEntries[p][w] = 0.;
                v_temp_ped[p][w]        = 0.;
                v_temp_pedvar[p][w]     = 0.;
            }
            fpedcal_n[telID][p][w] = 0.;
            fpedcal_mean[telID][p][w] = 0.;
            fpedcal_mean2[telID][p][w] = 0.;
        }
        // deroate the pixel coordinates
        if( getTelID() < getPointing().size() && getPointing()[getTelID()] )
        {
            getPointing()[getTelID()]->derotateCoords( VSkyCoordinatesUtilities::getUTC( getEventMJD(), getEventTime() ), x[p], y[p], xRot[p], yRot[p] );
        }
    }
    
    // fill the tree
    if( telID < fTree.size() && fTree[telID] )
    {
        fTree[telID]->Fill();
    }
    
    // fill the data vectors
    v_MJD[telID].push_back( getEventMJD() );
    v_time[telID].push_back( getEventTime() );
    v_pedEntries[telID].push_back( v_temp_pedEntries );
    v_ped[telID].push_back( v_temp_ped );
    v_pedvar[telID].push_back( v_temp_pedvar );
}


void VPedestalCalculator::doAnalysis( bool iLowGain )
{
    double t = getEventTime();
    // get right index for tel id
    unsigned int telID = 99;
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getTeltoAna()[i] == getTelID() )
        {
            telID = i;
            break;
        }
    }
    
    // temporary vectors
    if( telID < fTimeVec.size() )
    {
        if( fTimeVec[telID] == 0 )
        {
            fTimeVec[telID] = t;
        }
        ///////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////
        // end of a time slice
        else if( t - fTimeVec[telID] > fLengthofTimeSlice )
        {
            time = t;
            fillTimeSlice( telID );
            fTimeVec[telID] = t;
        }  // if( t - fTimeVec[telID] > fLengthofTimeSlice )
        ///////////////////////////////////////////////////////
        
        double i_tr_sum = 0.;
        // calculate the sums (don't use calcsums because it overwrites getSums() )
        // and fill the histograms
        for( unsigned int i = 0; i < getNChannels(); i++ )
        {
            if( i < getReader()->getNumChannelsHit() )
            {
                unsigned int chanID = 0;
                try
                {
                    chanID = fReader->getHitID( i );
                    
                    // don't use low gain channels
                    if( chanID >= getHiLo().size() || ( iLowGain && !getHiLo()[chanID] ) )
                    {
                        continue;
                    }
                    // check for dead channels
                    if( !getDead()[chanID] && chanID < fpedcal_mean[telID].size() )
                    {
                        fReader->selectHitChan( i );
                        
                        // set digital filter analysis (normally not used for VTS analysis)
                        fTraceHandler->setDigitalFilterParameters( getDigitalFilterMethod(),
                                getDigitalFilterUpSample(),
                                getDigitalFilterPoleZero() );
                                
                        fTraceHandler->setTrace( fReader, getNSamples(), getPeds()[chanID], chanID, i,
                                                 getLowGainMultiplier_Trace()*getHiLo()[chanID] );
                                                 
                        //////////////////////////
                        // loop over all summation windows
                        unsigned int iTempSW = fpedcal_mean[telID][chanID].size();
                        if( getRunParameter()->fCalibrationSumWindow < ( int )iTempSW )
                        {
                            iTempSW = getRunParameter()->fCalibrationSumWindow;
                        }
                        // w = 0 --> summation window 1
                        for( unsigned int w = 0; w < iTempSW; w++ )
                        {
                            // calculate trace sum
                            i_tr_sum = fTraceHandler->getTraceSum( fSumFirst, fSumFirst + ( w + 1 ), true, 1 );
                            if( i_tr_sum > 0. && i_tr_sum < 50.*( w + 1 ) )
                            {
                                if( chanID < fpedcal_n[telID].size() && w < fpedcal_n[telID][chanID].size() )
                                {
                                    fpedcal_n[telID][chanID][w]++;
                                    fpedcal_mean[telID][chanID][w] += i_tr_sum;
                                    fpedcal_mean2[telID][chanID][w] += i_tr_sum * i_tr_sum;
                                }
                                else
                                {
                                    cout << "ERROR (VPedestalCalculator): ";
                                    cout << telID << " " << fpedcal_n.size() << endl;
                                    cout << "\t" << chanID << " " << fpedcal_n[telID].size() << endl;
                                    cout << "\t" << w << " " << fpedcal_n[telID][chanID].size() << endl;
                                }
                            }
                        }
                    }
                }
                catch( ... )
                {
                    if( getDebugLevel() == 0 )
                    {
                        cout << "VPedestalCalculator::doAnalysis(), index out of range (fReader->getHitID) ";
                        cout << i << "(Telescope " << getTelID() + 1 << ", event " << getEventNumber() << ")" << endl;
                        setDebugLevel( 1 );
                    }
                    continue;
                }
            }
        }                                         // for (unsigned int i = 0; i < nhits; i++)
    }                                             // if( telID < fTree.size() && fTree[telID] && telID < fTimeVec.size() )
    
}


void VPedestalCalculator::terminate( bool iWrite, bool iDebug_IO )
{
    if( iWrite )
    {
        cout << endl;
        if( getOutputFile()->IsOpen() && getOutputFile()->cd() )
        {
            TDirectory* iDir = gDirectory;
            
            for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
            {
                if( getTeltoAna()[i] < getAnaDirectories().size() && getAnaDirectories()[getTeltoAna()[i]]->cd() )
                {
                    if( i < fTree.size() && fTree[i] )
                    {
                        cout << "\t number of time slices for Telescope " << getTeltoAna()[i] + 1 << ": " << fTree[i]->GetEntries() << endl;
                        int i_nbytes = fTree[i]->Write();
                        if( iDebug_IO )
                        {
                            cout << "WRITEDEBUG: pedestal trees (nbytes " << i_nbytes << "):";
                            if( getOutputFile() )
                            {
                                cout << getOutputFile()->Get( fTree[i]->GetName() );
                            }
                            cout << endl;
                        }
                    }
                }
            }
            iDir->cd();
        }
    }
    // fill last timing bin
    else
    {
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            // this is not for pointing checks, don't care about pixel rotation
            fillTimeSlice( i );
        }
    }
    
}


void VPedestalCalculator::reset()
{
    runNumber = 0;
    MJD = 0;
    time = 0.;
    for( unsigned int i = 0; i < fNPixel; i++ )
    {
        x[i] = 0.;
        y[i] = 0.;
        xRot[i] = 0.;
        yRot[i] = 0.;
    }
}

/*
 * stretch time slices to match the run lengths
 *
 * (avoid small time bins at the end of runs)
 */
double VPedestalCalculator::adjustTimeSliceLength( double iLengthofTimeSlice, double iRunStartTime, double iRunStoppTime )
{
    if( iRunStartTime < 0. || iRunStoppTime < 0. || iLengthofTimeSlice <= 0. )
    {
        return iLengthofTimeSlice;
    }
    
    double iRunLength = iRunStoppTime - iRunStartTime;
    
    double iNSlices = ( int )( ( iRunStoppTime - iRunStartTime ) / iLengthofTimeSlice );
    if( iNSlices > 0. )
    {
        double iR = iRunLength - iNSlices * iLengthofTimeSlice;
        return iLengthofTimeSlice + iR / iNSlices;
    }
    
    return iLengthofTimeSlice;
}
