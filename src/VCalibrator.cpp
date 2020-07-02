/*! \class VCalibrator
    \brief calibration class, calculation of pedestals, gains, ...

*/

#include "VCalibrator.h"

VCalibrator::VCalibrator()
{
    setCalibrated( false );
    fRaw = true;
    
    opfgain = 0;
    opftoff = 0;
    
    fCalibrationfileVersion = 1;
    
    fPedSingleOutFile = 0;
    fPedPerTelescopeTypeMinCnt = 1.E5;  // minimal counter for IPR measurements
}

/*
 * initialize all pedestals and IPR histograms
 *
 * called once per telescope and per run
 *
 */
bool VCalibrator::initializePedestalHistograms( ULong64_t iTelType, bool iLowGain,
                                                vector< double > minSumPerSumWindow, 
                                                vector< double > maxSumPerSumWindow )
{
    if( getDebugFlag() )
    {
        cout << "void VCalibrator::initializePedestalHistograms()" << endl;
    }
    
    string ioutfile;
    char ic[800];
    if( getDebugFlag() )
    {
        cout << "\t void VCalibrator::initializePedestalHistograms(): creating histograms and files:";
        cout << " tel " << getTelID() + 1 << ", type " << iTelType << endl;
    }
    // define output files
    if( !iLowGain )
    {
        ioutfile = fPedFileNameC[getTelID()] + ".root";
    }
    else
    {
        ioutfile = fLowGainPedFileNameC[getTelID()] + ".root";
    }
    
    // all pedestals are written to one output file (e.g. for CTA DST analysis)
    if( getRunParameter()->fPedestalSingleRootFile )
    {
        if( !fPedSingleOutFile )
        {
            fPedSingleOutFile = new TFile( ioutfile.c_str(), "RECREATE" );
            cout << "opened a single pedestal file for all telescope types: " << fPedSingleOutFile->GetName() << endl;
        }
    }
    else if( getTelID() <  getDetectorGeometry()->getTelType().size()
             && fPedOutFile.find( iTelType ) != fPedOutFile.end() )
    {
        fPedOutFile[iTelType] = new TFile( ioutfile.c_str(), "RECREATE" );
        if( fPedOutFile[iTelType]->IsZombie() )
        {
            cout << "VCalibrator::calculatePedestal error in creating pedestal file: ";
            cout << ioutfile << endl;
            exit( EXIT_FAILURE );
        }
    }
    else
    {
        cout << "VCalibrator::calculatePedestal error in creating pedestal file" << endl;
        cout << " (no single or array of output files defined)" << endl;
        exit( EXIT_FAILURE );
    }
    /////////////////////////////////////////////////////////////////
    // init histograms etc (only done before first event)
    if( fReader->getMaxChannels() )
    {
        unsigned int z = 0;
        // loop over all sumwindows
        for( unsigned int i = 0; i < hped_vec[iTelType].size(); i++ )
        {
                double min = 0.;
                if( i < minSumPerSumWindow.size() )
                {
                    min = minSumPerSumWindow[i] * 0.9;
                }
                double max = 1.;
                if( i < maxSumPerSumWindow.size() )
                {
                    max = maxSumPerSumWindow[i] * 1.1;
                }
                unsigned int nBins = static_cast<int>( max - min );
                if (nBins > 2000) nBins = 2000;
            /////////////////////////////////////////
            // histograms for pedestal calculation
            
            // loop over all channels
            for( unsigned int j = 0; j < hped_vec[iTelType][i].size(); j++ )
            {
                if( !hped_vec[iTelType][i][j] )
                {
                    char pedkey[100];
                    sprintf( pedkey, "hped_%d_%d_%d", ( int )iTelType, i + 1, j );
                    // If no trace, readout window is constant and pedestal does not grow with summation window
                    // if (getNSamples() - getSumFirst() > 0) max *= ( ( double )i + 1. );
                    sprintf( ic, "ped distribution (tel type %d, channel %d, sumwindow %d)", ( int )iTelType, j, i + 1 );
                    if( fPedestalsHistoClonesArray.find( iTelType ) != fPedestalsHistoClonesArray.end() )
                    {
                        hped_vec[iTelType][i][j] = ( TH1F* )fPedestalsHistoClonesArray[iTelType]->ConstructedAt( z );
                        hped_vec[iTelType][i][j]->SetName( pedkey );
                        hped_vec[iTelType][i][j]->SetTitle( ic );
                        hped_vec[iTelType][i][j]->SetBins( nBins, min, max );
                        z++;
                    }
                }
            }
            
            /////////////////////////////////////////
            // histograms for IPR graph calculation
                min = 0.;
                // max /= iMeanGain;
            
            // loop over all channels
            for( unsigned int j = 0; j < hpedPerTelescopeType[iTelType][i].size(); j++ )
            {
                if( !hpedPerTelescopeType[iTelType][i][j] )
                {
                    char pedkey[100];
                    sprintf( pedkey, "hpedPerTelescopeType_%d_%d_%d", ( int )iTelType, i + 1, j );
                    // If no trace, readout window is constant and pedestal does not grow with summation window
                    // if (getNSamples() - getSumFirst() > 0) max *= ( ( double )i + 1. );
                    sprintf( ic, "IPR distribution (tel type %d, channel %d, sumwindow %d)", ( int )iTelType, j, i + 1 );
                    if( fPedestalsHistoClonesArray.find( iTelType ) != fPedestalsHistoClonesArray.end() )
                    {
                        hpedPerTelescopeType[iTelType][i][j] = ( TH1F* )fPedestalsHistoClonesArray[iTelType]->ConstructedAt( z );
                        hpedPerTelescopeType[iTelType][i][j]->SetName( pedkey );
                        hpedPerTelescopeType[iTelType][i][j]->SetTitle( ic );
                        hpedPerTelescopeType[iTelType][i][j]->SetBins( nBins, min, max );
                        z++;
                    }
                }
            }
        }
    }
    setCalibrated();
    
    return true;
}


/*

   pedestal calculation and filling of charge histograms for IPR calculation

 */
void VCalibrator::calculatePedestals( bool iLowGain )
{
    if( getDebugFlag() )
    {
        cout << "void VCalibrator::calculatePedestals()" << endl;
    }
    
    /////////////////////////////////////////////////////////////
    // get telescope type of current telescope
    // all histograms in the following are filled per telescope type
    ULong64_t iTelType = 0;
    if( getTelID() < getDetectorGeometry()->getTelType().size() )
    {
        iTelType = getDetectorGeometry()->getTelType()[getTelID()];
    }
    else
    {
        cout << "VCalibrator::calculatePedestals: error: invalid telescope type : ";
        cout << getTelID() << "\t" << getDetectorGeometry()->getTelType().size() << endl;
        exit( EXIT_FAILURE );
    }
    //////////////////////////////////////////////////////////////////////////////
    // first call of this function (first event?): define histograms and output files for pedestals
    if( !getCalibrated() )
    {
        vector<double> maxSumPerSumWindow;
        vector<double> minSumPerSumWindow;
        // loop over all sumwindows
        for( unsigned int i = 0; i < hped_vec[iTelType].size(); i++ )
        {
            // calculate trace sums  (with calcSums(...iMakingPeds=true))
            // (always use trace integration method 1 here)
            calcSums( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + ( i + 1 ), true, iLowGain, 1 );
            
            // calculate the min/max sum for all channels
            double maxSum = 0.;
            double minSum = 1.e10;
            for( unsigned int j = 0; j < hped_vec[iTelType][i].size(); j++ )
            {
                // exclude low gain channels from pedestal calculation
                if( j < getHiLo().size() )
                {
                    if( iLowGain && !getHiLo()[j] )
                    {
                        continue;
                    }
                    else if( !iLowGain && getHiLo()[j] )
                    {
                        continue;
                    }
                }
                // calculate pedestals only for channels with valid hitbit
                //   important if trying to generate pedestals for zero supressed
                //   runs without injected pedestal events
                if( j < getSums().size() && getSums()[j] > 0. )
                {
                    if( fRunPar->fRunIsZeroSuppressed )
                    {
                        if( fReader->getChannelHitIndex( j ).first )
                        {
                            if (maxSum < getSums()[j])
                            {
                                maxSum = getSums()[j];
                            }
                            if( getSums()[j] > 0. && getSums()[j] < minSum )
                            {
                                minSum = getSums()[j];
                            }
                        }
                    }
                    else
                    {
                        if( maxSum < getSums()[j] )
                        {
                            maxSum = getSums()[j];
                        }
                        if( getSums()[j] > 0. && getSums()[j] < minSum )
                        {
                            minSum = getSums()[j];
                        }
                    }
                }
            }
            maxSumPerSumWindow.push_back( maxSum );
            minSumPerSumWindow.push_back( minSum );
        }
        
        initializePedestalHistograms( iTelType, iLowGain, minSumPerSumWindow, maxSumPerSumWindow );
    }

    // reset and fill vectors
    resetAnaData();
    fillHiLo();
    
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    // fill pedestal sum for current telescope type
    fNumberPedestalEvents[iTelType]++;
    
    // loop over all sumwindows
    for( unsigned int i = 0; i < hped_vec[iTelType].size(); i++ )
    {
        // calculate trace sums  (with calcSums(...iMakingPeds=true))
        // (always use trace integration method 1 here)
        calcSums( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + ( i + 1 ), true, iLowGain, 1 );
        
        // fill pedestal histograms for all channels
        for( unsigned int j = 0; j < hped_vec[iTelType][i].size(); j++ )
        {
            // exclude low gain channels from pedestal calculation
            if( j < getHiLo().size() )
            {
                if( iLowGain && !getHiLo()[j] )
                {
                    continue;
                }
                else if( !iLowGain && getHiLo()[j] )
                {
                    continue;
                }
            }
            // calculate pedestals only for channels with valid hitbit
            //   important if trying to generate pedestals for zero supressed
            //   runs without injected pedestal events
            if( j < getSums().size() && getSums()[j] > 0. )
            {
                if( fRunPar->fRunIsZeroSuppressed )
                {
                    if( fReader->getChannelHitIndex( j ).first )
                    {
                        hped_vec[iTelType][i][j]->Fill( getSums()[j] );
                    }
                }
                else
                {
                    hped_vec[iTelType][i][j]->Fill( getSums()[j] );
                }
            }
        }
        
        //////////////////////////////////////////////////////////
        // calculate trace sums for IPR graph (with calcSums(...iMakingPeds=false)) for charge extraction (not pedestal)
        // (should always be trace integration method 2)
        fTraceHandler->setIPRmeasure( !fRunPar->fCombineChannelsForPedestalCalculation );
        calcSums(
            getSumFirst(),
            getSumFirst() + ( i + 1 ),
            false,
            iLowGain,
            getTraceIntegrationMethod()
        );
        // loop over all channels
        for( unsigned int j = 0; j < hped_vec[iTelType][i].size(); j++ )
        {
            // exclude low gain channels from IPR calculation
            if( j < getHiLo().size() )
            {
                if( iLowGain && !getHiLo()[j] )
                {
                    continue;
                }
                else if( !iLowGain && getHiLo()[j] )
                {
                    continue;
                }
            }
            // calculate IPR charges only for channels with valid hitbit
            //   important if trying to generate pedestals for zero supressed
            //   runs without injected pedestal events
            if( j < getSums().size() && getSums()[j] > 0. )
            {
                if( fRunPar->fRunIsZeroSuppressed )
                {
                    if( fReader->getChannelHitIndex( j ).first )
                    {
                        hpedPerTelescopeType[iTelType][i][j]->Fill( getSums()[j] );
                    }
                }
                else
                {
                    hpedPerTelescopeType[iTelType][i][j]->Fill( getSums()[j] );
                }
            }
        }
        fTraceHandler->setIPRmeasure( false );
    }
}

/*

   write pedestals to disk

   this might be an ascii file and/or a root file

*/
void VCalibrator::writePeds( bool iLowGain, VPedestalCalculator* iPedestalCalculator, bool iWriteAsciiFile )
{
    if( getDebugFlag() )
    {
        cout << "void VCalibrator::writePeds()" << endl;
    }
    
    string ioutfile;
    // boolean to keep track which files have been written
    // (this is historically a mess)
    map< ULong64_t, bool > iFileWritten;
    for( unsigned int i = 0; i < getDetectorGeometry()->getTelType_list().size(); i++ )
    {
        iFileWritten[getDetectorGeometry()->getTelType_list()[i]] = false;
    }
    
    /////////////////////////////////////////////////////////////////
    // loop over all telescopes
    // (note that pedestal might be calculated per telescope type,
    //  but we write the results to disk per telescope)
    for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
    {
        unsigned int t = getTeltoAna()[tel];
        ULong64_t telType = getDetectorGeometry()->getTelType()[t];
        
        // make one pedestal file per telescope
        if( !iLowGain )
        {
            ioutfile = fPedFileNameC[t];
        }
        else
        {
            ioutfile = fLowGainPedFileNameC[t];
        }
        //////////////////////////////////
        // write histograms for first telescope of a certain teltype only
        if( !iFileWritten[telType] )
        {
        
            cout << "Telescope (type) " << telType << endl;
            cout << "\t total number of pedestal events: " << fNumberPedestalEvents[telType] << endl;
            cout << "\t writing ";
            if( iLowGain )
            {
                cout << " low ";
            }
            else
            {
                cout << " high ";
            }
            cout << " gain pedestals to " << endl;
            if( iWriteAsciiFile )
            {
                cout << "\t\t" << ioutfile << endl;
            }
            if( getPedestalRootFile( telType ) )
            {
                cout << "\t\t" << getPedestalRootFile( telType )->GetName() << endl;
            }
            ////////////////////////////////////////////////////////////
            // write pedestal and pedestal variances to an ascii file
            if( iWriteAsciiFile )
            {
                // loop over all channels
                ofstream os( ioutfile.c_str() );
                if( !os )
                {
                    cout << "VCalibrator::writePeds(): ERROR, unable to write pedestals to " << ioutfile << " (" << iLowGain << ")" << endl;
                    exit( EXIT_FAILURE );
                }
                for( unsigned int i = 0; i < hped_vec[telType][0].size(); i++ )
                {
                    // get pedestal and pedestal variances from pedestal histograms
                    // (require at least 100 entries in pedestal events)
                    os << t << " " << i << " ";
                    if( hped_vec[telType][fRunPar->fCalibrationSumWindow - 1][i]
                            && hped_vec[telType][fRunPar->fCalibrationSumWindow - 1][i]->GetEntries() > 100 )
                    {
                        os << hped_vec[telType][fRunPar->fCalibrationSumWindow - 1][i]->GetMean() / ( double )fRunPar->fCalibrationSumWindow << " ";
                    }
                    else
                    {
                        cout << "VCalibrator::writePeds(): WARNING, less than 100 events";
                        if( hped_vec[telType][fRunPar->fCalibrationSumWindow - 1][i] )
                        {
                            cout << "(" << hped_vec[telType][fRunPar->fCalibrationSumWindow - 1][i]->GetEntries() << ")";
                        }
                        cout << ", setting pedestal to 0 for telescope (type) ";
                        cout << telType << ", channel " << i << endl;
                        os << 0. << " ";
                    }
                    // loop over all window sizes
                    for( unsigned int j = 0; j < hped_vec[telType].size(); j++ )
                    {
                        if( hped_vec[telType][j][i] && hped_vec[telType][j][i]->GetEntries() > 100 )
                        {
                            os << hped_vec[telType][j][i]->GetRMS() << " ";
                        }
                        else
                        {
                            os << 0. << " ";
                        }
                    }
                    os << endl;
                }
                os.close();
            } // end writing ascii file
            
            ///////////////////////////////////////////////////////////////////////
            // write histograms to file
            if( !getPedestalRootFile( telType ) || !getPedestalRootFile( telType )->cd() )
            {
                cout << "VCalibrator::writePeds(): error accessing pedestal output file for telescope type " << telType << endl;
                cout << "...exiting" << endl;
                exit( EXIT_FAILURE );
            }
            
            // fill and write pedestal tree
            fillPedestalTree( tel, iPedestalCalculator );
            
            // write 1D histograms to directory calibration_TEL
            std::ostringstream iSname;
            iSname << "distributions_" << telType;
            TDirectory* i_dist = getPedestalRootFile( telType )->mkdir( iSname.str().c_str() );
            if( i_dist->cd() )
            {
                i_dist->cd();
                for( unsigned int i = 0; i < hped_vec[telType].size(); i++ )
                {
                    for( unsigned int j = 0; j < hped_vec[telType][i].size(); j++ )
                    {
                        if( hped_vec[telType][i][j] )
                        {
                            hped_vec[telType][i][j]->Write();
                        }
                    }
                }
                for( unsigned int i = 0; i < hpedPerTelescopeType[telType].size(); i++ )
                {
                    for( unsigned int j = 0; j < hpedPerTelescopeType[telType][i].size(); j++ )
                    {
                        if( hpedPerTelescopeType[telType][i][j] )
                        {
                            hpedPerTelescopeType[telType][i][j]->Write();
                        }
                    }
                }
            }
            iFileWritten[telType] = true;
        }
    }   // end loop over all telescopes
    // delete all histograms
    map< ULong64_t, TClonesArray* >::iterator i_PedestalsHistoClonesArray_iter;
    for( i_PedestalsHistoClonesArray_iter = fPedestalsHistoClonesArray.begin();
            i_PedestalsHistoClonesArray_iter != fPedestalsHistoClonesArray.end(); i_PedestalsHistoClonesArray_iter++ )
    {
        if( i_PedestalsHistoClonesArray_iter->second )
        {
            i_PedestalsHistoClonesArray_iter->second->Delete();
        }
    }
    // close all pedestal files
    cout << "Closing all pedestal files " << endl;
    map< ULong64_t, TFile* >::iterator iPedOutFile_iter;
    for( iPedOutFile_iter = fPedOutFile.begin(); iPedOutFile_iter != fPedOutFile.end(); iPedOutFile_iter++ )
    {
        if( iPedOutFile_iter->second && iPedOutFile_iter->second->IsOpen() )
        {
            cout << "\t closing " << iPedOutFile_iter->second->GetName() << endl;
            iPedOutFile_iter->second->Close();
        }
    }
    
    cout << "all files closed " << endl;
    if( fPedSingleOutFile && fPedSingleOutFile->IsOpen() )
    {
        fPedSingleOutFile->Close();
    }
    if( getDebugFlag() )
    {
        cout << "void VCalibrator::writePeds() END" << endl;
    }
    
    ///////////////////////////////////////////////////////////////////
    // calculate IPR graphs already at this point and write them to
    // the dst calibration file
    if( getRunParameter()->fCombineChannelsForPedestalCalculation
            && getRunParameter()->fPedestalSingleRootFile )
    {
        // write one IPR graph per telescope type and summation window
        
        cout << "Calculating IPR graphs " << endl;
        
        // keep track of telescope type
        map< ULong64_t, bool > iTelDone;
        for( unsigned int i = 0; i < getDetectorGeometry()->getTelType_list().size(); i++ )
        {
            iTelDone[getDetectorGeometry()->getTelType_list()[i]] = false;
        }
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            setTelID( i );
            
            ULong64_t iTelType = getTelType( i );
            
            if( !iTelDone[iTelType] )
            {
                if( iTelType >= hped_vec.size() )
                {
                    if( hped_vec.find( iTelType ) == hped_vec.end() )
                    {
                        continue;
                    }
                }
                for( unsigned int j = 0; j < hped_vec[iTelType].size(); j++ )
                {
                    // summation window
                    int i_sw = j + 1;
                    
                    calculateIPRGraphs( fPedSingleOutFile->GetName(), i_sw, iTelType, i );
                    iTelDone[getTelType( i )] = true;
                }
            }
        }
        writeIPRgraphs( fPedSingleOutFile->GetName() );
    }
    
}

TFile* VCalibrator::getPedestalRootFile( ULong64_t iT )
{
    if( fRunPar->fPedestalSingleRootFile )
    {
        return fPedSingleOutFile;
    }
    
    if( fPedOutFile.find( iT ) != fPedOutFile.end() )
    {
        return fPedOutFile[iT];
    }
    
    return 0;
}

/*

    fill results full pedestal calculation into a root tree

*/
bool VCalibrator::fillPedestalTree( unsigned int tel, VPedestalCalculator* iPedestalInTimeSlices )
{
    // get telescope ID
    unsigned int t = 0;
    if( tel < getTeltoAna().size() )
    {
        t = getTeltoAna()[tel];
    }
    else
    {
        return false;
    }
    
    // get telescope type
    ULong64_t iTelType = 0;
    if( t < getDetectorGeometry()->getTelType().size() )
    {
        iTelType = getDetectorGeometry()->getTelType()[t];
    }
    else
    {
        cout << "VCalibrator::fillPedestalTree: error: invalid telescope type : ";
        cout << t << "\t" << getDetectorGeometry()->getTelType().size() << endl;
        return false;
    }
    if( hped_vec.find( iTelType ) == hped_vec.end() )
    {
        cout << "VCalibrator::fillPedestalTree warning, telescope number out of range ";
        cout << t << "\t" << hped_vec.size() << endl;
        return false;
    }
    cout << "\t filling pedestal tree for telescope " << t << " (telescope type " << iTelType << ")" << endl;
    
    char iname[800];
    char ititle[800];
    
    UInt_t ichannel = 0;
    UInt_t insumw = 0;
    Float_t isumw[VDST_MAXSUMWINDOW];
    Float_t iped = 0.;
    Float_t ipedv[VDST_MAXSUMWINDOW];
    for( int i = 0; i < VDST_MAXSUMWINDOW; i++ )
    {
        isumw[i] = 0.;
        ipedv[i] = 0.;
    }
    UInt_t inevents = 0;
    // pedestals in time slices
    UInt_t TSnSlices = 0;
    if( iPedestalInTimeSlices && t < iPedestalInTimeSlices->v_MJD.size() )
    {
        TSnSlices = iPedestalInTimeSlices->v_MJD[t].size();
        
        cout << "\t number of time slices for pedestals in telescope " << t + 1 << ": " << TSnSlices << endl;
    }
    
    UInt_t TSMJD[VDST_PEDTIMESLICES];
    Double_t TStime[VDST_PEDTIMESLICES];
    UInt_t TSnevents[VDST_PEDTIMESLICES];
    Float_t TSpedmean[VDST_PEDTIMESLICES];
    for( int i = 0; i < VDST_PEDTIMESLICES; i++ )
    {
        TSMJD[i] = 0;
        TStime[i] = 0.;
        TSnevents[i] = 0;
        TSpedmean[i] = 0.;
    }
    // TSpedvar[Summationwindow-1][Time slice]
    vector< Float_t* > TSpedvar;
    
    std::ostringstream iSname;
    iSname << "tPeds_" << iTelType;
    std::ostringstream iStitle;
    iStitle << "pedestals, telescope type " << iTelType;
    TTree* iPedTree = new TTree( iSname.str().c_str(), iStitle.str().c_str() );
    iPedTree->Branch( "channel", &ichannel, "channel/i" );
    iPedTree->Branch( "nsumwindows", &insumw, "nsumwindow/i" );
    iPedTree->Branch( "sumwindow", isumw, "sumwindow[nsumwindow]/F" );
    iPedTree->Branch( "pedmean", &iped, "pedmean/F" );
    iPedTree->Branch( "pedvars", ipedv, "pedvars[nsumwindow]/F" );
    iPedTree->Branch( "nevents", &inevents, "nevents/i" );
    // fill this part only if time dependent pedestals are calculated
    if( iPedestalInTimeSlices )
    {
        iPedTree->Branch( "TSnSlices", &TSnSlices, "TSnSlices/i" );
        iPedTree->Branch( "TSMJD", TSMJD, "TSMJD[TSnSlices]/i" );
        iPedTree->Branch( "TStime", TStime, "TStime[TSnSlices]/D" );
        iPedTree->Branch( "TSnevents", TSnevents, "TSnevents[TSnSlices]/i" );
        iPedTree->Branch( "TSpedmean", TSpedmean, "TSpedmean[TSnSlices]/F" );
        for( unsigned int i = 0; i < hped_vec[iTelType].size(); i++ )
        {
            TSpedvar.push_back( new Float_t[VDST_PEDTIMESLICES] );
            sprintf( iname, "TSpedvar_sw%d[TSnSlices]/F", i + 1 );
            sprintf( ititle, "TSpedvar_sw%d", i + 1 );
            iPedTree->Branch( ititle, TSpedvar[i], iname );
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    // tree filling
    
    //////////////////////////////////////////
    // loop over all channels
    for( unsigned int i = 0; i < hped_vec[iTelType][0].size(); i++ )
    {
        ichannel = ( Int_t )i;
        
        // get pedestal and pedestal variances from pedestal histograms
        if( fRunPar->fCalibrationSumWindow > 0 && hasFADCData() )
        {
            iped = hped_vec[iTelType][fRunPar->fCalibrationSumWindow - 1][i]->GetMean() / ( double )fRunPar->fCalibrationSumWindow;
        }
        else if( !hasFADCData() )
        {
            if( hped_vec[iTelType][fRunPar->fCalibrationSumWindow - 1][i] )
            {
                iped = hped_vec[iTelType][fRunPar->fCalibrationSumWindow - 1][i]->GetMean();
            }
        }
        else
        {
            iped = 0.;
        }
        inevents = ( Int_t )hped_vec[iTelType][fRunPar->fCalibrationSumWindow - 1][i]->GetEntries();
        
        // loop over all summation window sizes
        insumw = ( UInt_t )hped_vec[iTelType].size();
        for( unsigned int j = 0; j < hped_vec[iTelType].size(); j++ )
        {
            isumw[j] = ( Float_t )j + 1;
            ipedv[j] = hped_vec[iTelType][j][i]->GetRMS();
        }
        
        /////////////////////////////////////////////////
        // time dependent pedestals
        if( iPedestalInTimeSlices )
        {
            // loop over all time slices
            TSnSlices = ( UInt_t )iPedestalInTimeSlices->v_MJD[tel].size();
            for( unsigned int ts = 0; ts < iPedestalInTimeSlices->v_MJD[tel].size(); ts++ )
            {
                TSMJD[ts] = iPedestalInTimeSlices->v_MJD[tel][ts];
                TStime[ts] = iPedestalInTimeSlices->v_time[tel][ts];
                if( i < iPedestalInTimeSlices->v_ped[tel][ts].size() && iPedestalInTimeSlices->v_ped[tel][ts][i].size() > 0 )
                {
                    unsigned int iSW_mean = iPedestalInTimeSlices->v_ped[tel][ts][i].size() - 1;
                    if( iSW_mean > ( unsigned int )getRunParameter()->fCalibrationSumWindow && getRunParameter()->fCalibrationSumWindow > 0 )
                    {
                        iSW_mean = getRunParameter()->fCalibrationSumWindow - 1;
                    }
                    TSnevents[ts] = ( UInt_t )iPedestalInTimeSlices->v_pedEntries[tel][ts][i][iSW_mean];
                    TSpedmean[ts] = iPedestalInTimeSlices->v_ped[tel][ts][i][iSW_mean];
                    // check summation windows
                    unsigned int iSW_temp = iPedestalInTimeSlices->v_pedvar[tel][ts][i].size();
                    if( iSW_temp > TSpedvar.size() )
                    {
                        iSW_temp = TSpedvar.size();
                    }
                    for( unsigned int w = 0; w < iSW_temp; w++ )
                    {
                        TSpedvar[w][ts] = iPedestalInTimeSlices->v_pedvar[tel][ts][i][w];
                    }
                }
            }
        }
        
        iPedTree->Fill();
    }
    
    iPedTree->Write();
    
    return true;
}

/*

   calculate average Tzeros and average trace times from data events

   use trace integration method 2 to find weighted position of pulse maximum

*/
void VCalibrator::calculateAverageTZero( bool iLowGain )
{
    if( getDebugFlag() )
    {
        cout << "VCalibrator::calculateAverageTZero() (";
        if( iLowGain )
        {
            cout << "low gain)" << endl;
        }
        else
        {
            cout << "high gain)" << endl;
        }
    }
    
    // use same histogram definitions as in toff calculation
    if( fReader->getMaxChannels() > 0 && htzero.size() == 0 )
    {
        /////////////////////
        // create histograms
        int telID = 0;
        for( unsigned int t = 0; t < getNTel(); t++ )
        {
            telID = t + 1;
            vector< TH1F* > ih_tzero;
            vector< TH1F* > ih_taverage;
            for( unsigned int i = 0; i < getNChannels(); i++ )
            {
                char toffkey[100];
                char ic[100];
                double imin = 0.;
                double imax = ( double )getNSamples();
                sprintf( toffkey, "htzero_%d_%d", telID, i );
                sprintf( ic, "TZero distribution (tel %d, channel %d)", telID, i );
                ih_tzero.push_back( new TH1F( toffkey, ic, 150, imin, imax ) );
                sprintf( toffkey, "htaverage_%d_%d", telID, i );
                sprintf( ic, "TAverageArrTime distribution (tel %d, channel %d)", telID, i );
                ih_taverage.push_back( new TH1F( toffkey, ic, 150, imin, imax ) );
            }
            htzero.push_back( ih_tzero );
            htaverage.push_back( ih_taverage );
        }
        ////////////////////
        // root output file
        string ioutfile;
        for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
        {
            unsigned int i = getTeltoAna()[tel];
            if( !iLowGain )
            {
                ioutfile = fTZeroFileNameC[i] + ".root";
            }
            else
            {
                ioutfile = fLowGainTZeroFileNameC[i] + ".root";
            }
            
            fTZeroOutFile.push_back( new TFile( ioutfile.c_str(), "RECREATE" ) );
            if( fTZeroOutFile.back()->IsZombie() )
            {
                cout << "Error: calculateAverageTZero() error, can't open output file: " << fTZeroOutFile.back()->GetName() << endl;
                exit( EXIT_FAILURE );
            }
            cout << "opened file for average tzeros: " << fTZeroOutFile.back()->GetName() << endl;
            setSpecialChannels();
        }
        
        cout << "calculate average pulse arrival times with summation window " << fRunPar->fCalibrationSumWindow;
        cout << " (start at " << fRunPar->fCalibrationSumFirst << ", window size for max integral finder: ";
        cout << fRunPar->fCalibrationSumWindowAverageTime << ")" << endl;
        cout << "     (require at least " << fRunPar->fCalibrationIntSumMin << " [dc] per channel/event)" << endl;
        
        setCalibrated();
    }
    fillHiLo();
    
    findDeadChans( iLowGain, ( fNumberTZeroEvents[fTelID] == 0 ) );
    
    ////////////////////////
    // calculate average arrival times (trace integration method 2)
    calcSums( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + fRunPar->fCalibrationSumWindowAverageTime, false, false, 2 );
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        if( iLowGain && !getHiLo()[i] )
        {
            continue;
        }
        if( !iLowGain && getHiLo()[i] )
        {
            continue;
        }
        
        // require a min sum for tzero filling
        if( getSums()[i] > fRunPar->fCalibrationIntSumMin && !getDead()[i] && !getMasked()[i] )
        {
            if( getTraceAverageTime( false )[i] > 0. && i < htaverage[fTelID].size() && htaverage[fTelID][i] )
            {
                htaverage[fTelID][i]->Fill( getTraceAverageTime( false )[i] );
            }
        }
    }
    ////////////////////////
    // fill tzeros
    
    // calculate tzeros and sums
    fRaw = false;
    calcTZerosSums( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + fRunPar->fCalibrationSumWindow, 1 );
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        if( iLowGain && !getHiLo()[i] )
        {
            continue;
        }
        if( !iLowGain && getHiLo()[i] )
        {
            continue;
        }
        
        // require a min sum for tzero filling
        if( getSums()[i] > fRunPar->fCalibrationIntSumMin && !getDead()[i] && !getMasked()[i] )
        {
            if( getTZeros()[i] > 0. && i < htzero[fTelID].size() && htzero[fTelID][i] )
            {
                htzero[fTelID][i]->Fill( getTZeros()[i] );
            }
        }
    }
    fNumberTZeroEvents[fTelID]++;
    
}

/*

   write average Tzeros and average trace times from data events to disk

   (histogram and summary tree)

*/
void VCalibrator::writeAverageTZeros( bool iLowGain )
{
    if( getDebugFlag() )
    {
        cout << "VCalibrator::writeAverageTZeros() (";
        if( iLowGain )
        {
            cout << "low gain)" << endl;
        }
        else
        {
            cout << "high gain)" << endl;
        }
    }
    
    for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
    {
        unsigned int t = getTeltoAna()[tel];
        setTelID( t );
        
        if( tel >= fTZeroOutFile.size() || !fTZeroOutFile[tel] )
        {
            cout << "VCalibrator::writeAverageTZeros() error writing average tzeros for telescope " << t + 1 << endl;
            cout << "   (tel " << tel << ", " << fTZeroOutFile.size() << ")" << endl;
            continue;
        }
        cout << "\t writing average tzeros to ";
        cout << fTZeroOutFile[tel]->GetName() << endl;
        cout << "\t calculated from " << fNumberTZeroEvents[t] << " events" << endl;
        
        fTZeroOutFile[tel]->cd();
        if( t < htzero.size() && htaverage.size() )
        {
            for( unsigned int i = 0; i < htzero[t].size(); i++ )
            {
                if( htzero[t][i] )
                {
                    htzero[t][i]->Write();
                }
            }
            for( unsigned int i = 0; i < htaverage[t].size(); i++ )
            {
                if( htaverage[t][i] )
                {
                    htaverage[t][i]->Write();
                }
            }
            TTree* i_tree = fillCalibrationSummaryTree( t, "TZero", htzero[t] );
            if( i_tree )
            {
                i_tree->Write();
            }
            TTree* i_treeAV = fillCalibrationSummaryTree( t, "TAverage", htaverage[t] );
            if( i_treeAV )
            {
                i_treeAV->Write();
            }
        }
        fTZeroOutFile[tel]->Close();
    }
}


/*
 * calculate / fill relative gain and timing (toffsets) histograms
 *
 * use a laser/flasher run
 *
 */
void VCalibrator::calculateGainsAndTOffsets( bool iLowGain )
{
    if( getDebugFlag() )
    {
        cout << "VCalibrator::calculateGainsAndTOffsets()" << endl;
    }
    ////////////////////////////////////////////////
    // initialize output files and histograms
    if( fReader->getMaxChannels() > 0 && hgain.size() == 0 )
    {
        findDeadChans( iLowGain );
        
        // root output files are gain/toff file.root
        if( !iLowGain )
        {
            opfgain = new TFile( ( fGainFileNameC[getTelID()] + ".root" ).c_str(), "RECREATE" );
        }
        else
        {
            opfgain = new TFile( ( fLowGainGainFileNameC[getTelID()] + ".root" ).c_str(), "RECREATE" );
        }
        if( opfgain->IsZombie() )
        {
            cout << "calculateGainsAndTOffsets() error, can't open output file: " << opfgain->GetName() << endl;
            exit( EXIT_FAILURE );
        }
        cout << "recreated gain file: " << opfgain->GetName() << endl;
        cout << "calculate gains and toffsets with summation window " << fRunPar->fCalibrationSumWindow;
        cout << " (start at " << fRunPar->fCalibrationSumFirst << ")" << endl;
        
        // create histograms
        for( unsigned int i = 0; i < getNChannels(); i++ )
        {
            char gainkey[100];
            char ic[100];
            sprintf( gainkey, "hgain_%d", i );
            double imin = 0.0;
            double imax = 5.0;
            sprintf( ic, "gain distribution (tel %d, channel %d)", getTelID() + 1, i );
            hgain.push_back( new TH1F( gainkey, ic, 150, imin, imax ) );
            
            char pulsekey[100];
            sprintf( pulsekey, "hpulse_%d", i );
            imin = 0;
            imax = getNSamples();
            hpulse.push_back( new TProfile( pulsekey, "Mean pulse", ( int )( imax - imin ), imin, imax, -100., 10000. ) );
            
            char spekey[100];
            sprintf( spekey, "htcpulse_%d", i );
            htcpulse.push_back( new TProfile( spekey, "time corrected mean pulse", 100, imin, imax ) );
        }
        if( !iLowGain )
        {
            opftoff = new TFile( ( fToffFileNameC[getTelID()] + ".root" ).c_str(), "RECREATE" );
        }
        else
        {
            opftoff = new TFile( ( fLowGainToffFileNameC[getTelID()] + ".root" ).c_str(), "RECREATE" );
        }
        if( opftoff->IsZombie() )
        {
            cout << "calculateGainsAndTOffsets() error, can't open output file: " << opftoff->GetName() << endl;
            exit( EXIT_FAILURE );
        }
        for( unsigned int i = 0; i < getNChannels(); i++ )
        {
            char toff_vs_sumkey[100];
            char ic[100];
            sprintf( toff_vs_sumkey, "htoff_vs_sum_%d", i );
            double imin = 5;
            double imax = 405;
            htoff_vs_sum.push_back( new TProfile( toff_vs_sumkey, "TOff vs Sum", 20, imin, imax ) );
            
            char toffkey[100];
            sprintf( toffkey, "htoff_%d", i );
            imin = -10;
            imax = 10;
            sprintf( ic, "TOffset distribution (tel %d, channel %d)", getTelID() + 1, i );
            htoff.push_back( new TH1F( toffkey, ic, 150, imin, imax ) );
        }
        
        if( getRunParameter()->fWriteExtraCalibTree )
        {
            opfgain->cd();
            
            //Extra calib output.
            fExtra_sum = new vector<double>( getNChannels(), -99 );
            fExtra_ped = new vector<double>( getNChannels(), -99 );
            fExtra_pedVar = new vector<double>( getNChannels(), -99 );
            fExtra_tzero = new vector<double>( getNChannels(), -99 );
            fExtra_HiLo = new vector<short>( getNChannels(), -99 );
            fExtra_sumfirst = new vector<short>( getNChannels(), -99 );
            fExtra_sumwindow = new vector<short>( getNChannels(), -99 );
            fExtra_dead = new vector<short>( getNChannels(), -99 );
            fExtra_use = new vector<short>( getNChannels(), -99 );
            fExtra_QMon = 0;
            fExtra_nMon = 0;
            fExtra_TZeroMon = 0;
            fExtra_nHiLo = 0;
            fExtra_eventNumber = 0;
            fExtra_nPix = getNChannels();
            TString title = TString::Format( "charges_%d", fTelID + 1 );
            tExtra_ChargeTree = new TTree( title.Data(), "extra calib output (charges/monitor charge per event)" );
            tExtra_ChargeTree->Branch( "eventNumber", &fExtra_eventNumber );
            tExtra_ChargeTree->Branch( "QMon", &fExtra_QMon );
            tExtra_ChargeTree->Branch( "TZeroMon", &fExtra_TZeroMon );
            tExtra_ChargeTree->Branch( "nMon", &fExtra_nMon );
            tExtra_ChargeTree->Branch( "nHiLo", &fExtra_nHiLo );
            tExtra_ChargeTree->Branch( "nPix", &fExtra_nPix );
            tExtra_ChargeTree->Branch( "Q", &fExtra_sum->at( 0 ), "sum[nPix]/D" );
            tExtra_ChargeTree->Branch( "ped", &fExtra_ped->at( 0 ), "ped[nPix]/D" );
            tExtra_ChargeTree->Branch( "pedVar", &fExtra_pedVar->at( 0 ), "pedVar[nPix]/D" );
            tExtra_ChargeTree->Branch( "sumwindow", &fExtra_sumwindow->at( 0 ), "sumwindow[nPix]/S" );
            tExtra_ChargeTree->Branch( "sumfirst", &fExtra_sumfirst->at( 0 ), "sumfirst[nPix]/S" );
            tExtra_ChargeTree->Branch( "dead", &fExtra_dead->at( 0 ), "dead[nPix]/S" );
            tExtra_ChargeTree->Branch( "use", &fExtra_use->at( 0 ), "use[nPix]/S" );
            tExtra_ChargeTree->Branch( "HiLo", &fExtra_HiLo->at( 0 ), "HiLo[nPix]/S" );
            tExtra_ChargeTree->Branch( "tzero", &fExtra_tzero->at( 0 ), "tzero[nPix]/D" );
            
        }
        
        setSpecialChannels();
        
        setCalibrated();
    }
    
    fillHiLo();
    
    ////////////////////////////////////////////////
    //! calculate sums and tzeros
    calcSums( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + fRunPar->fCalibrationSumWindow, false );
    calcTZeros( fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumFirst + fRunPar->fCalibrationSumWindow );
    
    if( fRunPar->fL2TimeCorrect )
    {
        FADCStopCorrect();
    }
    
    //////////////////////////////////////////
    // test if this is a laser event (by sum of total charge in all channels )
    bool i_laser = false;
    double i_total = 0;
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        i_total += getSums()[i];
    }
    // check if image sum is larger than a minimum value
    // (given in command line with -lasermin=$LASERMIN
    if( i_total > fRunPar->fLaserSumMin )
    {
        i_laser = true;
    }
    if( i_laser )
    {
        fNumberGainEvents[fTelID]++;
        // write pulse to disk, pulse histograms for each event in one directory
        char i_name[200];
        char i_title[200];
        TDirectory* i_curDir = ( TDirectory* )opfgain;
        TDirectory* i_pulseDir = 0;
        TH1D* i_pulse = 0;
        if( getRunParameter()->fwriteLaserPulseN > 0 )
        {
            i_curDir->cd();
            sprintf( i_name, "zEvent_%d", getEventNumber() );
            i_pulseDir = gDirectory->mkdir( i_name, i_name );
        }
        
        // need number of saturated PMTs for laser calibration
        int n_lowgain = 0;
        for( unsigned int i = 0; i < fReader->getNumChannelsHit(); i++ )
        {
            fReader->selectHitChan( i );
            if( getHiLo()[i] != 0 && !getDead()[i] )
            {
                n_lowgain += 1;
            }
        }
        // don't do anything if there are not enough channels in low or high gain
        if( !iLowGain && n_lowgain >  5 )
        {
            return;
        }
        // for LG calibration: request a minimum number of 70% of the pixels in LG
        if( iLowGain && n_lowgain < ( int )( 0.7 * ( double )getNChannels() ) )
        {
            return;
        }
        // check for number of saturated channels
        int n_saturated = 0;
        for( unsigned int i = 0; i < fReader->getNumChannelsHit(); i++ )
        {
            if( getTraceRawMax()[i] == 255 )
            {
                n_saturated++;
            }
        }
        // reject event if more than 10 channels are saturated
        if( n_saturated > 10 )
        {
            return;
        }
        
        int counttzero = 0;
        int countsums = 0;
        double m_tzero = 0;
        double m_sums = 0;
        // calculate total sum
        for( unsigned int i = 0; i < getNChannels(); i++ )
        {
            if( iLowGain && !getHiLo()[i] )
            {
                continue;
            }
            if( !iLowGain && getHiLo()[i] )
            {
                continue;
            }
            
            if( ( getSums()[i] > fRunPar->fCalibrationIntSumMin || fRunPar->fLaserSumMin < 0. )  && !getDead()[i] && !getMasked()[i] )
            {
                countsums++;
                m_sums  += getSums()[i];
                if( getTZeros()[i] >= 0. )
                {
                    counttzero++;
                    m_tzero += getTZeros()[i];
                }
            }
        }
        if( countsums > 0 )
        {
            m_sums /= countsums;
        }
        if( counttzero > 0 )
        {
            m_tzero /= counttzero;
        }
        
        if( getRunParameter()->fWriteExtraCalibTree )
        {
            fExtra_nMon = countsums;
            fExtra_QMon = m_sums;
            fExtra_TZeroMon = m_tzero;
            fExtra_nHiLo = n_lowgain;
            fExtra_eventNumber = fEventNumber;
        }
        
        
        // subtract pedestals
        // set number of entries in hpulse to number of pulses added up
        double tcorr = 0.;
        for( unsigned int i = 0; i < getNChannels(); i++ )
        {
            if( getRunParameter()->fWriteExtraCalibTree )
            {
                //Extra calib output.
                fExtra_sum->at( i ) = getSums()[i];
                fExtra_tzero->at( i ) = getTZeros()[i];
                fExtra_ped->at( i ) = getPeds()[i];
                fExtra_HiLo->at( i ) = getHiLo()[i];
                fExtra_sumfirst->at( i ) = getRunParameter()->fCalibrationSumFirst;
                fExtra_sumwindow->at( i ) = getSumWindow();
                fExtra_dead->at( i ) = getDead()[i];
                fExtra_use->at( i ) = 0; //will be set to 1 later is sum>limit etc.
            }
            
            
            int i_entries = ( int )hpulse[i]->GetEntries();
            if( getRunParameter()->fwriteLaserPulseN > 0 )
            {
                i_pulseDir->cd();
                sprintf( i_name, "h_%d", i );
                sprintf( i_title, "event %d, channel %d", getEventNumber(), i );
                i_pulse = new TH1D( i_name, i_title, hpulse[i]->GetNbinsX(), hpulse[i]->GetXaxis()->GetXmin(), hpulse[i]->GetXaxis()->GetXmax() );
            }
            
            if( m_sums > 0.1 || fRunPar->fLaserSumMin < 0. )
            {
                // fill average traces
                if( getRunParameter()->fwriteAverageLaserPulse )
                {
                    double this_content = 0.;
                    unsigned int chanID = 0;
                    int this_bin = 0;
                    for( unsigned int k = 0; k < fReader->getNumChannelsHit(); k++ )
                    {
                        try
                        {
                            chanID = fReader->getHitID( k );
                            if( chanID == i )
                            {
                                for( unsigned int j = 0; j < ( unsigned int )hpulse[i]->GetNbinsX(); j++ )
                                {
                                    this_bin = ( int )( j + fCalData[getTeltoAnaID()]->fFADCStopOffsets[i] + 1 );
                                    if( this_bin > 0 && this_bin <= hpulse[i]->GetNbinsX() )
                                    {
                                        this_content = fReader->getSample_double( chanID, this_bin, ( this_bin == 0 ) ) - getPeds( iLowGain )[i];
                                        if( getRunParameter()->fwriteLaserPulseN > 0 )
                                        {
                                            i_pulse->SetBinContent( this_bin, this_content );
                                        }
                                        hpulse[i]->Fill( this_bin, this_content );
                                        // time corrected pulse
                                        if( getTZeros()[i] > -99. )
                                        {
                                            tcorr = getTZeros()[i] - m_tzero;
                                        }
                                        else
                                        {
                                            tcorr = 0.;
                                        }
                                        htcpulse[i]->Fill( ( double )this_bin - tcorr, this_content );
                                    }
                                }
                            }
                        }
                        catch( ... )
                        {
                            cout << "VCalibrator::calculateGainsAndTOffsets() index out of range " << k << endl;
                        }
                    }
                }
                // check status of hi/lo gain bits
                if( iLowGain && !getHiLo()[i] )
                {
                    continue;
                }
                if( !iLowGain && getHiLo()[i] )
                {
                    continue;
                }
                /////////////////////////////////////////////
                // fill gain and toffset histograms
                if( ( getSums()[i] > fRunPar->fCalibrationIntSumMin || fRunPar->fLaserSumMin < 0. ) && !getDead()[i] && !getMasked()[i] )
                {
                    if( getRunParameter()->fWriteExtraCalibTree )
                    {
                        fExtra_use->at( i ) = 1;
                    }
                    
                    hgain[i]->Fill( ( float )getSums()[i] / m_sums );
                    if( getTZeros()[i] >= 0. )
                    {
                    
                        htoff[i]->Fill( ( float )getTZeros()[i] - m_tzero );
                        htoff_vs_sum[i]->Fill( ( float )getSums()[i], ( float )getTZeros()[i] - m_tzero );
                    }
                }
                
                if( getRunParameter()->fwriteLaserPulseN > 0 )
                {
                    i_pulse->Write();
                }
            }
            if( i_entries != hpulse[i]->GetEntries() )
            {
                hpulse[i]->SetEntries( i_entries + 1 );
            }
            if( i_entries != htcpulse[i]->GetEntries() )
            {
                htcpulse[i]->SetEntries( i_entries + 1 );
            }
        }
        if( getRunParameter()->fwriteLaserPulseN > 0 )
        {
            getRunParameter()->fwriteLaserPulseN--;
        }
        i_curDir->cd();
        
        if( getRunParameter()->fWriteExtraCalibTree )
        {
            tExtra_ChargeTree->Fill();
        }
    }
}

/*
 * write gain and toffsets results to output file
 *
 */
void VCalibrator::writeGains( bool iLowGain )
{
    string iFile;
    for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
    {
        int t = getTeltoAna()[tel];
        setTelID( t );
        if( !iLowGain )
        {
            iFile = fGainFileNameC[t];
        }
        else
        {
            iFile = fLowGainGainFileNameC[t];
        }
        
        cout << "Telescope " << t + 1 << endl;
        cout << "\t total number of laser events: " << fNumberGainEvents[t] << endl;
        if( !opfgain )
        {
            cout << "VCalibrator::writeGains error: no output file for gain values" << endl;
            continue;
        }
        cout << "\t writing gains to " << iFile << " and " << opfgain->GetName() << endl;
        ofstream os( iFile.c_str() );
        if( !os )
        {
            cout << "VCalibrator::writeGains(): ERROR unable to write to " << iFile << endl;
            exit( EXIT_FAILURE );
        }
        else
        {
            for( unsigned int i = 0; i < getNChannels(); i++ )
            {
                os   << i << " " << hgain[i]->GetMean() << " " << hgain[i]->GetRMS() << endl;
            }
        }
        os.close();
        
        opfgain->cd();
        for( unsigned int i = 0; i < hgain.size(); i++ )
        {
            hgain[i]->Write();
        }
        
        if( getRunParameter()->fwriteAverageLaserPulse )
        {
            for( unsigned int i = 0; i < hpulse.size(); i++ )
            {
                hpulse[i]->Write();
                htcpulse[i]->Write();
            }
        }
        if( getRunParameter()->fWriteExtraCalibTree )
        {
            tExtra_ChargeTree->Write();
        }
        
        TTree* iTG = fillCalibrationSummaryTree( t, "gain", hgain );
        if( iTG )
        {
            iTG->Write();
        }
        opfgain->Close();
    }
}

/*

   generic function to write calibration tree for e.g. gains, toffs, average tzeros

*/
TTree* VCalibrator::fillCalibrationSummaryTree( unsigned int itel, string iName, vector<TH1F* > h )
{
    setTelID( itel );
    
    char iname[200];
    char ititle[200];
    
    int ichannel = 0;
    float i_mean = 0.;
    float i_median = 0.;
    float i_rms = 0.;
    
    sprintf( iname, "t%s_%d", iName.c_str(), itel + 1 );
    sprintf( ititle, "%s, telescope %d", iName.c_str(), itel + 1 );
    TTree* t = new TTree( iname, ititle );
    t->Branch( "channel", &ichannel, "channel/I" );
    sprintf( iname, "%s", iName.c_str() );
    sprintf( ititle, "%s/F", iName.c_str() );
    t->Branch( iname, &i_mean, ititle );
    sprintf( iname, "%smed", iName.c_str() );
    sprintf( ititle, "%smed/F", iName.c_str() );
    t->Branch( iname, &i_median, ititle );
    sprintf( iname, "%svar", iName.c_str() );
    sprintf( ititle, "%svar/F", iName.c_str() );
    t->Branch( iname, &i_rms, ititle );
    
    double i_a[] = { 0.5 };
    double i_b[] = { 0.0 };
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        ichannel = ( int )i;
        if( i < h.size() && h[i] && h[i]->GetEntries() > 0 )
        {
            i_mean = h[i]->GetMean();
            i_rms  = h[i]->GetRMS();
            h[i]->GetQuantiles( 1, i_b, i_a );
            i_median = i_b[0];
        }
        else
        {
            i_mean = 0.;
            i_rms  = 0.;
            i_median = 0.;
        }
        t->Fill();
    }
    
    return t;
}


void VCalibrator::writeTOffsets( bool iLowGain )
{
    string iFile;
    for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
    {
        int t = getTeltoAna()[tel];
        setTelID( t );
        if( !iLowGain )
        {
            iFile = fToffFileNameC[t];
        }
        else
        {
            iFile = fLowGainToffFileNameC[t];
        }
        if( !opftoff )
        {
            cout << "VCalibrator::writeTOffsets error writing time offsets to " << iFile << endl;
            continue;
        }
        cout << "\t writing time offsets to " << iFile << " and " << opftoff->GetName() << endl;
        
        ofstream os( iFile.c_str() );
        
        if( !os )
        {
            cout << "VCalibrator::writeTOffsets() ERROR: unable to write to " << iFile << endl;
            exit( EXIT_FAILURE );
        }
        else for( unsigned int i = 0; i < getNChannels(); i++ )
            {
                os << i << " " << htoff[i]->GetMean() << " " << htoff[i]->GetRMS() << endl;
            }
        os.close();
        
        opftoff->cd();
        for( unsigned int i = 0; i < htoff.size(); i++ )
        {
            htoff[i]->Write();
        }
        for( unsigned int i = 0; i < htoff_vs_sum.size(); i++ )
        {
            htoff_vs_sum[i]->SetErrorOption( "S" );
            htoff_vs_sum[i]->Write();
        }
        TTree* iTT = fillCalibrationSummaryTree( t, "toff", htoff );
        if( iTT )
        {
            iTT->Write();
        }
        opftoff->Close();
    }
}


void VCalibrator::terminate( VPedestalCalculator* iP )
{
    if( fRunPar->frunmode == 1 || fRunPar->frunmode == 6 )
    {
        writePeds( fRunPar->frunmode == 6, iP, !fRunPar->fPedestalSingleRootFile );
    }
    else if( fRunPar->frunmode == 2 || fRunPar->frunmode == 5 )
    {
        writeGains( fRunPar->frunmode == 5 );
        writeTOffsets( fRunPar->frunmode == 5 );
    }
    else if( fRunPar->frunmode == 7 || fRunPar->frunmode == 8 )
    {
        writeAverageTZeros( fRunPar->frunmode == 8 );
    }
}


/*!

   read calibration data from the different files or from the database

 */
void VCalibrator::readCalibrationData()
{
    if( fDebug )
    {
        cout << "VCalibrator::readCalibrationData " << endl;
    }
    
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        setTelID( getTeltoAna()[i] );
        
        // read high gain gains
        if( getRunParameter()->frunmode != 2 && getRunParameter()->frunmode != 5 )
        {
            readGains( false );
            // read low gain gains
            if( fLowGainGainFileNameC[getTeltoAna()[i]].size() > 0 )
            {
                readGains( true );
                setLowGainGains();
            }
            // set low gain gains equal to high gain gains
            else
            {
                getGains( true ) = getGains( false );
            }
        }
        
        // read average tzeros
        if( getRunParameter()->frunmode != 2 && getRunParameter()->frunmode != 5
                && getRunParameter()->frunmode != 7 && getRunParameter()->frunmode != 8 )
        {
            readAverageTZeros( false );
            if( fLowGainTZeroFileNameC[getTeltoAna()[i]].size() > 0 )
            {
                readAverageTZeros( true );
            }
            else
            {
                getAverageTZeros( true ) = getAverageTZeros( false );
            }
        }
        
        // read high gain pedestals
        readPeds( fPedFileNameC[getTelID()], false, getSumWindow() );
        
        // use for simulations a fixed value for low-gain pedestal
        if( getRunParameter()->fsimu_lowgain_pedestal_DefaultPed > 0. )
        {
            getPedsLowGain() = getRunParameter()->fsimu_lowgain_pedestal_DefaultPed;
            getPedvarsLowGain() = getPedvars();
            getPedvarsAllSumWindows( true ) = getPedvarsAllSumWindows( false );
            getmeanPedvarsAllSumWindow( true ) = getmeanPedvarsAllSumWindow( false );
            getmeanRMSPedvarsAllSumWindow( true ) = getmeanRMSPedvarsAllSumWindow( false );
        }
        // read low gain pedestals
        //old version (low gain run number)
        else if( fLowGainPedFileNameC[getTeltoAna()[i]].size() > 0 )
        {
            readPeds( fLowGainPedFileNameC[getTeltoAna()[i]], true, getSumWindow() );
        }
        //new version (combined file)
        else if( fNewLowGainPedFileNameC[getTeltoAna()[i]].size() > 0 )
        {
            readPeds( fNewLowGainPedFileNameC[getTeltoAna()[i]], true, getSumWindow() );
        }
        // no low gain peds available -> set low gain peds to high gain peds (not sure if this is a good idea)
        else
        {
            getPedsLowGain() = getPeds();
            getPedvarsLowGain() = getPedvars();
            getPedvarsAllSumWindows( true ) = getPedvarsAllSumWindows( false );
            getmeanPedvarsAllSumWindow( true ) = getmeanPedvarsAllSumWindow( false );
            getmeanRMSPedvarsAllSumWindow( true ) = getmeanRMSPedvarsAllSumWindow( false );
        }
        
        // read low gain multipliers
        readLowGainMultiplier( );
        
        if( getRunParameter()->frunmode != 2 && getRunParameter()->frunmode != 5 )
        {
            // read high gain channel time offsets
            readTOffsets( false );
            // read low gain channel time offsets
            if( fLowGainToffFileNameC.size() > 0 )
            {
                readTOffsets( true );
                setLowGainTOff();
            }
            else
            {
                getTOffsets( true ) = getTOffsets( false );
            }
        }
        
        // read pixel status
        readPixelstatus();
        
        // read summation window dependend values for second summation window
        readPeds( fPedFileNameC[getTeltoAna()[i]], false, getSumWindow_2() );
        
        // use for simulations a fixed value for low-gain pedestal
        if( getRunParameter()->fsimu_lowgain_pedestal_DefaultPed > 0. )
        {
            for( unsigned int g = 0; g < getPedsLowGain().size(); g++ )
            {
                if( fabs( getPedsLowGain()[g] ) < 0.5 || fabs( getPedvarsLowGain()[g] ) < 0.5 )
                {
                    getPedsLowGain()[g] = getRunParameter()->fsimu_lowgain_pedestal_DefaultPed;
                }
            }
        }
        
        // read low gain pedestals
        //old version (low gain run number)
        else if( fLowGainPedFileNameC[getTeltoAna()[i]].size() > 0 )
        {
            readPeds( fLowGainPedFileNameC[getTeltoAna()[i]], true, getSumWindow() );
        }
        //new version (combined file)
        else if( fNewLowGainPedFileNameC[getTeltoAna()[i]].size() > 0 )
        {
            readPeds( fNewLowGainPedFileNameC[getTeltoAna()[i]], true, getSumWindow() );
        }
        // (preli) fill high gain value if low value is not available
        else
        {
            for( unsigned int g = 0; g < getPedsLowGain().size(); g++ )
            {
                if( fabs( getPedsLowGain()[g] ) < 0.5 || fabs( getPedvarsLowGain()[g] ) < 0.5 )
                {
                    getPedsLowGain()[g] = getPeds()[g];
                }
            }
        }
        if( fLowGainPedFileNameC[getTeltoAna()[i]].size() > 0 )
        {
            getCalibrationData()->recoverLowGainPedestals();
        }
        setCalibrated();
    }                                             // end loop over all telescopes
    
}

/*

    reading pedestal from root file

*/
bool VCalibrator::readPeds_from_rootfile( string iFile, bool iLowGain, unsigned int i_SumWindow )
{
    ifstream infile_root;
    infile_root.open( ( iFile + ".root" ).c_str(), ifstream::in );
    
    //  do not read from here for low gain channels
    if( infile_root && usePedestalsInTimeSlices( iLowGain ) )
    {
        infile_root.close();
        iFile += ".root";
        TFile iFPed( iFile.c_str() );
        char hname[200];
        sprintf( hname, "tPeds_%d", fTelID + 1 );
        TTree* tPed = ( TTree* )gDirectory->Get( hname );
        if( !tPed )
        {
            cout << "VCalibrator::readPeds_from_rootfile WARNING:";
            cout << "found root file with pedestals but no pedestal tree (telescope ";
            cout << getTelID() + 1 << "): " << iFile << endl;
            return false;
        }
        cout << "Telescope " << getTelID() + 1;
        if( !iLowGain )
        {
            cout << ": reading pedestals for high gain channels";
        }
        else
        {
            cout << ": reading pedestals for low gain channels";
        }
        cout << " and sumwindow " << i_SumWindow << " from :" << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << iFile << endl;
        
        UInt_t ichannel = 0;
        UInt_t insumw = 0;
        Float_t isumw[VDST_MAXSUMWINDOW];
        for( unsigned int i = 0; i < VDST_MAXSUMWINDOW; i++ )
        {
            isumw[i] = 0.;
        }
        Float_t iped = 0.;
        Float_t ipedv[VDST_MAXSUMWINDOW];
        for( unsigned int i = 0; i < VDST_MAXSUMWINDOW; i++ )
        {
            ipedv[i] = 0.;
        }
        UInt_t iTSnSlices = 0;
        UInt_t iTSMJD[VDST_PEDTIMESLICES];
        for( unsigned int i = 0; i < VDST_PEDTIMESLICES; i++ )
        {
            iTSMJD[i] = 0;
        }
        Double_t iTStime[VDST_PEDTIMESLICES];
        for( unsigned int i = 0; i < VDST_PEDTIMESLICES; i++ )
        {
            iTStime[i] = 0.;
        }
        Float_t iTSpedmean[VDST_PEDTIMESLICES];
        for( unsigned int i = 0; i < VDST_PEDTIMESLICES; i++ )
        {
            iTSpedmean[i] = 0.;
        }
        vector< Float_t* > iTSpedvars;
        for( unsigned int i = 0; i < getNSamples(); i++ )
        {
            iTSpedvars.push_back( new Float_t[VDST_PEDTIMESLICES] );
        }
        
        tPed->SetBranchAddress( "channel", &ichannel );
        tPed->SetBranchAddress( "nsumwindows", &insumw );
        tPed->SetBranchAddress( "sumwindow", isumw );
        tPed->SetBranchAddress( "pedmean", &iped );
        tPed->SetBranchAddress( "pedvars", ipedv );
        tPed->SetBranchAddress( "TSnSlices", &iTSnSlices );
        tPed->SetBranchAddress( "TSMJD", iTSMJD );
        tPed->SetBranchAddress( "TStime", iTStime );
        tPed->SetBranchAddress( "TSpedmean", iTSpedmean );
        char i_branchname[100];
        for( unsigned int i = 0; i < getNSamples(); i++ )
        {
            sprintf( i_branchname, "TSpedvar_sw%d", i + 1 );
            if( tPed->GetBranchStatus( i_branchname ) )
            {
                tPed->SetBranchAddress( i_branchname, iTSpedvars[i] );
            }
        }
        // initialize time dependant pedestal vectors
        if( tPed->GetEntries() > 0 )
        {
            valarray< double > i_ts_temp_peds;
            i_ts_temp_peds.resize( getPeds( false, -99. ).size(), 0. );
            valarray< valarray< double > > i_ts_temp_pedvar;
            i_ts_temp_pedvar.resize( getNSamples() + 1, i_ts_temp_peds );
            
            // for first event only, expect that all channels have same number of time slices
            tPed->GetEntry( 0 );
            
            if( iTSnSlices > 0 )
            {
                // prepare time vectors
                getCalData()->getMJDTS_vector().resize( iTSnSlices, 0 );
                getCalData()->getTimeTS_vector().resize( iTSnSlices, 0. );
                
                // prepare pedestal vector
                getCalData()->getPedsTS_vector( iLowGain ).resize( iTSnSlices, i_ts_temp_peds );
                
                // prepare pedestal variation vector
                getCalData()->getPedvarsVTS_vector( iLowGain ).resize( iTSnSlices, i_ts_temp_pedvar );
            }
        }
        /////////////////////////////////////////////////////////////////////////////////
        // loop over tree entries and read pedestals and pedestal variations from tree
        /////////////////////////////////////////////////////////////////////////////////
        
        // loop over all channels
        for( int i = 0; i < tPed->GetEntries(); i++ )
        {
            tPed->GetEntry( i );
            
            // fix summation window
            if( insumw > getLargestSumWindow() )
            {
                insumw = getLargestSumWindow();
            }
            
            if( ichannel >= getPeds( iLowGain, -99. ).size() )
            {
                continue;
            }
            
            // read pedestal and fill pedestal distributions
            getPeds( iLowGain, -99. )[ichannel] = iped;
            if( !iLowGain )
            {
                getPedDist()->Fill( iped );
            }
            else
            {
                getPedLowGainDist()->Fill( iped );
            }
            // read pedestal variation as function of summation window
            for( unsigned int sw = 0; sw < insumw; sw++ )
            {
                unsigned int icurrent_sw = ( unsigned int )TMath::Nint( isumw[sw] );
                if( ( isumw[sw] == i_SumWindow || ( iLowGain && isumw[sw] == 23 ) ) && ichannel < getPedvars( iLowGain, i_SumWindow, -99. ).size() )
                {
                    getPedvars( iLowGain, i_SumWindow, -99. )[ichannel] = ipedv[sw];
                    if( isumw[sw] == i_SumWindow )
                    {
                        if( !iLowGain )
                        {
                            getPedvarsDist()->Fill( ipedv[sw] );
                        }
                        else
                        {
                            getPedvarsLowGainDist()->Fill( ipedv[sw] );
                        }
                    }
                }
                if( icurrent_sw < getPedvarsAllSumWindows( iLowGain ).size() && ichannel < getPedvarsAllSumWindows( iLowGain )[sw].size() )
                {
                    getPedvarsAllSumWindows( iLowGain )[icurrent_sw][ichannel] = ipedv[sw];
                }
                else
                {
                    cout << "VCalibrator::readPeds error:";
                    cout << " incompatible summation window from root file: ";
                    cout << sw << "\t" << isumw[sw] << "\t" << getPedvarsAllSumWindows( iLowGain ).size() << endl;
                }
                if( sw == 1 && ichannel < getPedrms( iLowGain ).size() )
                {
                    getPedrms( iLowGain )[ichannel] = ipedv[sw];
                }
            }
            if( insumw > 0 && isumw[insumw - 1] < i_SumWindow )
            {
                cout << "VCalibrator::readPeds error:";
                cout << "no pedestal found for this sumwindow : ";
                cout << i_SumWindow << " (" << isumw[insumw - 1] << "), tel, channel " << fTelID + 1 << ", " << ichannel << endl;
                if( iLowGain )
                {
                    cout << "VCalibrator::readPeds: using pedestals for window " << isumw[insumw - 1] << " for low gain channels" << endl;
                }
            }
            // read date for pedestals in time slices
            if( iTSnSlices > 0 )
            {
                // fill vectors
                for( unsigned int ts = 0; ts < iTSnSlices; ts++ )
                {
                    // read time for first entry only
                    if( i == 0 )
                    {
                        getCalData()->getMJDTS_vector()[ts] = iTSMJD[ts];
                        getCalData()->getTimeTS_vector()[ts] = iTStime[ts];
                    }
                    // pedestals and pedvars
                    getCalData()->getPedsTS_vector( iLowGain )[ts][i] = iTSpedmean[ts];
                    for( unsigned int sw = 0; sw < insumw; sw++ )
                    {
                        if( sw < getCalData()->getPedvarsVTS_vector( iLowGain )[ts].size() )
                        {
                            getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw][i] = iTSpedvars[sw][ts];
                        }
                    }
                }                         // for( unsigned int ts = 0; ts < iTSnSlices; ts++ )
            }                             // if( iTSnSlices > 0 )
        }                                 // for( int i = 0; i < tPed->GetEntries(); i++ )
        tPed->ResetBranchAddresses();
        
        /////////////////////////////////////////////////////////////////////////////////
        // calculate mean and RMS pedvars for time slices
        /////////////////////////////////////////////////////////////////////////////////
        valarray< double > its_mean( 0., getNSamples() + 1 );
        valarray< double > its_rms( 0., getNSamples() + 1 );
        getCalData()->getMeanPedvarsVTS_vector( iLowGain ).resize( getCalData()->getPedvarsVTS_vector( iLowGain ).size(), its_mean );
        getCalData()->getMeanRMSPedvarsVTS_vector( iLowGain ).resize( getCalData()->getPedvarsVTS_vector( iLowGain ).size(), its_rms );
        
        for( unsigned int ts = 0; ts < getCalData()->getPedvarsVTS_vector( iLowGain ).size(); ts++ )
        {
            its_mean.resize( getCalData()->getPedvarsVTS_vector( iLowGain )[ts].size(), 0. );
            its_rms.resize( getCalData()->getPedvarsVTS_vector( iLowGain )[ts].size(), 0. );
            
            for( unsigned int sw = 0; sw < getCalData()->getPedvarsVTS_vector( iLowGain )[ts].size(); sw++ )
            {
                double its_n = 0.;
                double its_sum = 0.;
                double its_sum2 = 0.;
                for( unsigned int ts_c = 0; ts_c < getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw].size(); ts_c++ )
                {
                    if( getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw][ts_c] >= 1.0 )
                    {
                        its_n++;
                        its_sum  += getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw][ts_c];
                        its_sum2 += getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw][ts_c] * getCalData()->getPedvarsVTS_vector( iLowGain )[ts][sw][ts_c];
                    }
                }
                if( its_n > 1. )
                {
                    its_mean[sw] = its_sum / its_n;
                    its_rms[sw]  = sqrt( 1. / ( its_n - 1. ) * ( its_sum2 - 1. / its_n * its_sum * its_sum ) );
                }
                else
                {
                    its_mean[sw] = 0.;
                    its_rms[sw]  = 0.;
                }
            }
            getCalData()->getMeanPedvarsVTS_vector( iLowGain )[ts] = its_mean;
            getCalData()->getMeanRMSPedvarsVTS_vector( iLowGain )[ts] = its_rms;
            // everything seemed to worked find
            return true;
        }
        
    }
    
    return false;
}


/*
   read pedestals from text file
*/
bool VCalibrator::readPeds_from_textfile( string iFile, bool iLowGain, unsigned int i_SumWindow )
{
    ifstream infile;
    infile.open( iFile.c_str(), ifstream::in );
    if( !infile )
    {
        cout << "VCalibrator::readPeds_from_textfile error: unable to open pedestal file " << iFile << endl;
        cout << "\t exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    cout << "Telescope " << getTelID() + 1;
    if( !iLowGain )
    {
        cout << ": reading pedestals (txt) for high gain channels";
    }
    else
    {
        cout << ": reading pedestals (txt) for low gain channels";
    }
    cout << " and sumwindow " << i_SumWindow << " from: " << endl;
    cout << "Telescope " << getTelID() + 1 << ": ";
    cout << iFile << endl;
    
    // resetting all vectors
    getPeds( iLowGain, -99. ) = 1.0;
    getPedvars( iLowGain, i_SumWindow, -99. ) = 1.0;
    getPedrms( iLowGain ) = 1.0;
    unsigned int tel = 0;
    unsigned int ch = 0;
    float mean = 0.;
    float rms = 0.;
    int i_testCounter = 0;
    // loop over pedestal file (text format)
    string i_Line;
    while( getline( infile, i_Line ) )
    {
        if( i_Line.size() > 0 )
        {
            std::istringstream is_stream( i_Line );
            // telescope number
            is_stream >> tel;
            if( tel == fTelID )
            {
                i_testCounter++;
                // channel number
                is_stream >> ch;
                // pedestal
                is_stream >> mean;
                if( ch < getPeds( iLowGain, -99. ).size() )
                {
                    getPeds( iLowGain, -99. )[ch] = mean;
                    if( !iLowGain )
                    {
                        getPedDist()->Fill( mean );
                    }
                    else
                    {
                        getPedLowGainDist()->Fill( mean );
                    }
                }
                unsigned int count = 0;
                // pedestal variances
                do
                {
                    count += 1;
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> rms;
                    }
                    // ther might be more than NSamples values in the file; stop when reached NSamples
                    if( ( count == i_SumWindow || ( iLowGain && count == getNSamples() ) )  && ch < getPedvars( iLowGain, i_SumWindow, -99. ).size() )
                    {
                        getPedvars( iLowGain, i_SumWindow, -99. )[ch] = rms;
                        if( count == i_SumWindow )
                        {
                            if( !iLowGain )
                            {
                                getPedvarsDist()->Fill( rms );
                            }
                            else
                            {
                                getPedvarsLowGainDist()->Fill( rms );
                            }
                            
                        }
                    }
                    if( count <  getPedvarsAllSumWindows( iLowGain ).size() && ch < getPedvarsAllSumWindows( iLowGain )[count].size() )
                    {
                        getPedvarsAllSumWindows( iLowGain )[count][ch] = rms;
                    }
                    if( count == 1 && ch < getPedrms( iLowGain ).size() )
                    {
                        getPedrms( iLowGain )[ch] = rms;
                    }
                }
                while( !(is_stream>>std::ws).eof() );
                if( count < i_SumWindow )
                {
                    cout << "VCalibrator::readPeds_from_textfile error:";
                    cout << "no pedestal found for this sumwindow : " << i_SumWindow << " (" << count << "), tel, channel " << tel + 1 << ", " << ch << endl;
                    if( iLowGain )
                    {
                        cout << "VCalibrator::readPeds_from_textfile: using pedestals for window " << count << " for low gain channels" << endl;
                    }
                }
            }
        }
    }
    if( i_testCounter == 1 )
    {
        cout << "VCalibrator::readPeds: problem with reading pedestal file" << endl;
        cout << "\t probably old pedestal file format. Rerun pedestal calculation" << endl;
        cout << "\t exiting....." << endl;
        exit( EXIT_FAILURE );
    }
    
    return true;
}

/*

   read low gain pedestals from combined file

*/
bool VCalibrator::readPeds_from_combinedfile( string iFile, bool iLowGain, unsigned int i_SumWindow )
{

    ifstream infile;
    infile.open( iFile.c_str(), ifstream::in );
    if( !infile )
    {
        cout << "VCalibrator::readPeds_from_combinedfile error: unable to open pedestal file " << iFile << endl;
        return false;
    }
    
    cout << "Telescope " << getTelID() + 1;
    if( !iLowGain )
    {
        cout << ": reading pedestals (combined) for high gain channels";
    }
    else
    {
        cout << ": reading pedestals (combined) for low gain channels";
    }
    cout << " and sumwindow " << i_SumWindow << " from: " << endl;
    cout << "Telescope " << getTelID() + 1 << ": ";
    cout << iFile << endl;
    
    // resetting all vectors
    getPeds( iLowGain, -99. ) = 1.0;
    getPedvars( iLowGain, i_SumWindow, -99. ) = 1.0;
    getPedrms( iLowGain ) = 1.0;
    float mean = 0.;
    float rms = 0.;
    vector<double> vars;
    int i_testCounter = 0;
    
    //read out DB to get current FADC to pixel relations.
    vector<vector<int> > fadc_modules;
    vector<vector<int> > fadc_channels; //per telescope and channel
    
    // loop over long pedestal file (text format)
    string i_Line;
    
    std::map< std::pair<int, int>, vector<int> > runs;
    std::map< std::pair<int, int>, vector<double> > peds;
    std::map< std::pair<int, int>, vector<vector<double> > > pedvars;
    
    while( getline( infile, i_Line ) )
    {
        int i_module, i_chan, i_run, i;
        double i_ped, i_pedvar;
        
        if( i_Line.size() > 0 )
        {
            vector<double> i_pedvars;
            std::istringstream is_stream( i_Line );
            is_stream >> i_module >> i_chan >> i_run >> i >> i >> i_ped;
            while( is_stream >> i_pedvar )
            {
                i_pedvars.push_back( i_pedvar ) ;
            }
            if( i_ped == 0 )   //don't use runs where no pedestal was measured.
            {
                continue;
            }
            std::pair<int, int> itemp( i_module, i_chan );
            runs[itemp].push_back( i_run );
            peds[itemp].push_back( i_ped );
            pedvars[itemp].push_back( i_pedvars );
        }
    }
    
    //now find closest runs for each tel/pixel and fill
    
    for( unsigned int iChan = 0; iChan < getPedvars( iLowGain, i_SumWindow, -99. ).size(); iChan++ )
    {
        int iFADCModule;
        int iFADCChannel;
        iFADCModule = getDBPixelDataReader()->getFADC_module( getTelID(), iChan );
        iFADCChannel = getDBPixelDataReader()->getFADC_channel( getTelID(), iChan );
        
        //special catch for swapped channels at start of season 2013/14. This should be fixed in the data base.
        //see elog http://veritash.sao.arizona.edu:8081/VERITAS-Operations/11544
        if( getRunParameter()->frunnumber >= 69474 && getRunParameter()->frunnumber <= 69641 )
        {
            if( getTelID() == 3 && iChan >= 220 && iChan < 230 )
            {
                iFADCModule = 26;
            }
            else if( getTelID() == 3 && iChan >= 230 && iChan < 240 )
            {
                iFADCModule = 18 ;
            }
        }
        
        std::pair<int, int> temp( iFADCModule , iFADCChannel );
        vector<int> i_runs = runs[temp];
        vector<double> i_peds = peds[temp];
        vector< vector<double> > i_pedvars = pedvars[temp];
        
        
        mean = 0;
        rms = 0;
        vars = std::vector<double>( getPedvarsAllSumWindows( iLowGain ).size() , 0 );
        
        
        if( i_runs.size() == 0 )
        {
            i_testCounter++;
        }
        else
        {
            unsigned int iLower = 0;
            unsigned int iUpper = i_runs.size() - 1;
            for( unsigned int iT = 0; iT < i_runs.size(); iT++ )
            {
                if( i_runs.at( iT ) / 100 < getRunNumber() && i_runs.at( iT ) > i_runs.at( iLower ) )
                {
                    iLower = iT;
                }
                if( i_runs.at( iT ) / 100 > getRunNumber() && i_runs.at( iT ) < i_runs.at( iUpper ) )
                {
                    iUpper = iT;
                }
            }
            
            //int theIndex = iLower;
            //if( TMath::Abs( i_runs.at( iUpper )/100 - getRunNumber() ) < TMath::Abs( i_runs.at( iLower )/100 ) - getRunNumber() )  theIndex = iUpper ;
            int theIndex = iUpper;
            mean = i_peds.at( theIndex );
            vars = std::vector<double>( getPedvarsAllSumWindows( iLowGain ).size() , 0.0 );
            
            for( unsigned int iS = 1; iS < getPedvarsAllSumWindows( iLowGain ).size() && iS < i_pedvars.at( theIndex ).size(); iS++ )
            {
                vars.at( iS ) = i_pedvars.at( theIndex ).at( iS );
            }
            
            rms = vars.size() > 1 ? vars.at( 1 ) : 0.0;
            
        }
        
        getPeds( iLowGain, -99 )[iChan] = mean;
        getPedvars( iLowGain, i_SumWindow, -99. )[iChan] = vars.size() > i_SumWindow ?  vars.at( i_SumWindow ) : 0.0 ;
        getPedrms( iLowGain )[iChan] = rms;
        for( unsigned int iS = 1; iS < getPedvarsAllSumWindows( iLowGain ).size(); iS++ )
        {
            getPedvarsAllSumWindows( iLowGain )[iS][iChan] = vars.at( iS );
        }
        
        if( iLowGain )
        {
            getPedLowGainDist()->Fill( mean ) ;
            getPedvarsLowGainDist()->Fill( vars.size() > i_SumWindow ?  vars.at( i_SumWindow ) : 0.0 );
            
        }
        else
        {
            getPedDist()->Fill( mean );
            getPedvarsDist()->Fill( vars.size() > i_SumWindow ?  vars.at( i_SumWindow ) : 0.0 );
            
        }
        
    }//channel
    
    if( i_testCounter > 40 )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "VCalibrator::readPeds_from_combinedfile(): Warning: did not find valid pedestals for ";
        cout << i_testCounter << " channels in " << iFile << ", trying traditional pedestal text file." << endl;
        return false;
    }
    else
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "VCalibrator::readPeds_from_combinedfile(): found valid pedestals for ";
        cout << getPedvars( iLowGain, i_SumWindow, -99. ).size() - i_testCounter << " channels in " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << iFile << endl;
        return true;
    }
    
}

bool VCalibrator::readPeds_from_grisufile( bool iLowGain, unsigned int i_SumWindow )
{
    cout << "Telescope " << getTelID() + 1 << ": reading peds from 'P' lines (sumwindow " << i_SumWindow;
    if( iLowGain )
    {
        cout << ", low gain)";
    }
    else
    {
        cout << ", high gain)";
    }
    cout << endl;
    
    setPedsFromPLine();
    fReader->setSumWindow( getTelID(), i_SumWindow );
    getPeds( iLowGain, -99. ) = fReader->getPeds();
    getPedvars( iLowGain, i_SumWindow, -99. ) = fReader->getPedvars();
    getPedvarsAllSumWindows( iLowGain ) = fReader->getPedvarsAllSumWindows();
    getPedrms( iLowGain ) = fReader->getPedRMS();
    
    // fill pedestal distributions
    if( !iLowGain )
    {
        for( unsigned int i = 0; i < getPeds( iLowGain, -99. ).size(); i++ )
        {
            getPedDist()->Fill( getPeds( iLowGain, -99. )[i] );
        }
        for( unsigned int i = 0; i < getPedvars( iLowGain, i_SumWindow, -99. ).size(); i++ )
        {
            getPedvarsDist()->Fill( getPedvars( iLowGain, i_SumWindow, -99. )[i] );
        }
    }
    else
    {
        for( unsigned int i = 0; i < getPeds( iLowGain, -99. ).size(); i++ )
        {
            getPedLowGainDist()->Fill( getPeds( iLowGain, -99. )[i] );
        }
        for( unsigned int i = 0; i < getPedvars( iLowGain, i_SumWindow, -99. ).size(); i++ )
        {
            getPedvarsLowGainDist()->Fill( getPedvars( iLowGain, i_SumWindow, -99. )[i] );
        }
    }
    return true;
}

/*!

    read pedestals from root file, text file or get if from the Monte Carlo file

*/
bool VCalibrator::readPeds( string i_pedfile, bool iLowGain, unsigned int i_SumWindow )
{
    if( getDebugFlag() )
    {
        cout << "VCalibrator::readPeds() " << getTelID() + 1 << "\t" << i_pedfile << endl;
    }
    
    // set input file
    //////////////////////////////////////////////
    // reset histograms with pedestal distributions
    if( !iLowGain )
    {
        getPedvarsDist()->Reset();
        getPedDist()->Reset();
    }
    else
    {
        getPedvarsLowGainDist()->Reset();
        getPedLowGainDist()->Reset();
    }
    if( i_pedfile.size() > 0 && getRunParameter()->fsimu_pedestalfile.size() == 0 )
    {
        // read pedestals from root file
        if( !readPeds_from_rootfile( i_pedfile, iLowGain, i_SumWindow ) )
        {
            // for low gain: read peds from a combined file
            if( !getRunParameter()->fuseDB || !iLowGain || !readPeds_from_combinedfile( i_pedfile, iLowGain, i_SumWindow ) )
            {
                // in rare cases: read peds from a text file (no sum window support!)
                readPeds_from_textfile( i_pedfile, iLowGain, i_SumWindow );
            }
        }
    }
    ////////////////////////////////////////////////////////////////
    // getting pedestals directly from MC (grisu) file
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    else if( fReader->getDataFormat() == "grisu" || ( getRunParameter()->fsourcetype == 2 && getRunParameter()->fsimu_pedestalfile.size() > 0 ) )
    {
        readPeds_from_grisufile( iLowGain, i_SumWindow );
    }
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    
    // mean pedvar calculation
    if( getPedvarsAllSumWindows( iLowGain ).size() == getmeanPedvarsAllSumWindow( iLowGain ).size() )
    {
        double i_n = 0.;
        double i_sum = 0.;
        double i_sum2 = 0.;
        for( unsigned int i = 0; i < getPedvarsAllSumWindows( iLowGain ).size(); i++ )
        {
            i_n = 0.;
            i_sum = 0.;
            i_sum2 = 0.;
            for( unsigned int s = 0; s < getPedvarsAllSumWindows( iLowGain )[i].size(); s++ )
            {
                if( getPedvarsAllSumWindows( iLowGain )[i][s] >= 1.0 )
                {
                    i_sum  += getPedvarsAllSumWindows( iLowGain )[i][s];
                    i_sum2 += getPedvarsAllSumWindows( iLowGain )[i][s] * getPedvarsAllSumWindows( iLowGain )[i][s];
                    i_n += 1.;
                }
            }
            if( i_n > 1. )
            {
                getmeanPedvarsAllSumWindow( iLowGain )[i] = i_sum / i_n;
                getmeanRMSPedvarsAllSumWindow( iLowGain )[i] = sqrt( 1 / ( i_n - 1 ) * ( i_sum2 - 1. / i_n * i_sum * i_sum ) );
            }
            else
            {
                getmeanPedvarsAllSumWindow( iLowGain )[i] = 0.;
                getmeanRMSPedvarsAllSumWindow( iLowGain )[i] = 0.;
            }
        }
    }
    
    if( iLowGain )
    {
        setLowGainPedestals();
    }
    
    return true;
}


string VCalibrator::getCalibrationFileName( int iTel, int irun, string iSuffix, string name )
{
    if( irun <= 0 )
    {
        return "";
    }
    stringstream iFileStr;
    
    if( iSuffix == "newlped" )
    {
        iFileStr << getRunParameter()->getDirectory_EVNDISPCalibrationData();
        iFileStr << name ;
        return iFileStr.str();
    }
    iFileStr << getRunParameter()->getDirectory_EVNDISPCalibrationData_perRun();
    iFileStr << "Tel_" << iTel + 1 << "/";
    iFileStr << irun << ".";
    iFileStr << iSuffix;
    
    return iFileStr.str();
}


void VCalibrator::readfromVOFFLINE_DB( int gain_or_toff, string& iFile, vector< unsigned int >& Vchannel, vector< double >& Vmean, vector< double >& Vrms )
{
    int LOW_GAIN = 0;
    if( gain_or_toff != 1 && gain_or_toff != 2 )
    {
        std::cout << "ERROR VCalibrator::readfromVOFFLINE_DB: gain_or_toff must be 1 or 2" << std::endl;
        return;
    }
    
    iFile += "_DB";
    if( getRunParameter()->freadCalibfromDB_versionquery > 0 )
    {
        iFile += "_version";
        char cversion[10];
        sprintf( cversion, "%d", getRunParameter()->freadCalibfromDB_versionquery );
        string sversion = cversion;
        iFile += sversion;
    }
    // to avoid double file
    time_t rawtime;
    struct tm* timeinfo;
    time( &rawtime );
    timeinfo = localtime( &rawtime );
    char today_now[1000]  ;
    sprintf( today_now, "_%04d%02d%02d_%02dh%02dm%02ds", 1900 + timeinfo->tm_year, 1 + timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec );
    string NOW = today_now;
    iFile += NOW;
    
    // fill the data into vectors
    if( !getRunParameter()->freadCalibfromDB_save_file )
    {
        iFile = "";
    }
    
    TString DB_server = getRunParameter()->getDBServer();
    VDB_CalibrationInfo fDB_calibinfo( getRunParameter()->fGainFileNumber[getTelID()], getTelID() + 1, iFile,
                                       gain_or_toff, getRunParameter()->freadCalibfromDB_versionquery, LOW_GAIN, DB_server );
                                       
                                       
    //--- reading the VOFFLINE DB (if something goes wrong in the reading the analysis is interrupted)
    if( fDB_calibinfo.readVOFFLINE() )
    {
        Vmean = fDB_calibinfo.getVectorMean();
        Vrms = fDB_calibinfo.getVectorVariance();
        Vchannel = fDB_calibinfo.getVectorChannelList();
        return;
    }
    else
    {
    
        std::cout << "ERROR: calibration information could not be read from the VOFFLine DB " << std::endl;
        std::cout << "exiting..." << std::endl;
        std::cout << "            please launch EVNDISP without the option: readcalibdb, if you want to be independent from the VOFFLINE DB" << std::endl;
        exit( EXIT_FAILURE );
    }
    
    
}


/*!

 */
void VCalibrator::readGains( bool iLowGain )
{
    if( fDebug )
    {
        cout << "VCalibrator::readGains() " << fReader->getDataFormat();
        cout << "\t low gain: " << iLowGain << endl;
        cout << endl;
    }
    
    string iFile = fGainFileNameC[getTelID()];
    if( iLowGain )
    {
        iFile = fLowGainGainFileNameC[getTelID()];
    }
    int nZero = 0;	//counts number of channels with 0 gain.
    
    // don't read gains for runmode = 2
    if( iFile.size() > 0 && getRunParameter()->frunmode != 2 )
    {
        setGains( 1.0, iLowGain );
        setGainvars( 1.0, iLowGain );
        bool use_default = false;
        vector< unsigned int > VchannelList;
        vector < double > Vmean;
        vector < double > Vvar;
        cout << "Telescope " << getTelID() + 1 << ":";
        cout << " reading relative gains";
        if( iLowGain )
        {
            cout << " for low gain channels ";
        }
        else
        {
            cout << " for high gain channels ";
        }
        cout << "from: " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        // read gain from DB
        if( !iLowGain && getRunParameter()->freadCalibfromDB && getRunParameter()->fcalibrationfile.size() == 0 )
        {
            cout << "VOFFLINE DB" << endl;
            readfromVOFFLINE_DB( 1, iFile, VchannelList, Vmean, Vvar );
        }
        // read gains from external file
        else
        {
            cout << iFile << endl;
            
            ifstream infile;
            infile.open( iFile.c_str(), ifstream::in );
            if( !infile )
            {
                cout << "VCalibrator::readGains() warning: Input file " << iFile << " cannot be opened.\n" ;
                use_default = true;
                if( getRunParameter()->fDBRunType == "flasher" || getRunParameter()->fDBRunType == "laser" )
                {
                    cout << "VCalibrator::readGains() info: all gains set to 1.0\n" ;
                    setGains_DefaultValue( true, iLowGain );
                }
                else
                {
                    if( !getRunParameter()->fNoCalibNoPb )
                    {
                        std::cout << "Error: No gain information available " << std::endl;
                        std::cout << "exiting... " << std::endl;
                        std::cout << "          please set option: -nocalibnoproblem , when launching EVNDISP if you want to continue anyway (all gains will be set to 1) " << std::endl;
                        exit( EXIT_FAILURE );
                        
                    }
                    setGains_DefaultValue( false, iLowGain );
                    cout << "VCalibrator::readGains() info: all gains set to 1.0\n" ;
                }
                
                
                
            }
            
            char buffer[100];
            setGains( 1.0, iLowGain );
            setGainvars( 1.0, iLowGain );
            if( !use_default )
            {
                while( !infile.eof() )
                {
                    infile.getline( buffer, 100 );
                    unsigned int ch;
                    float mean, rms;
                    sscanf( buffer, "%u %f %f", &ch, &mean, &rms );
                    VchannelList.push_back( ch );
                    Vmean.push_back( mean );
                    Vvar.push_back( rms );
                    if( mean == 0 )
                    {
                        nZero++;
                    }
                }
            }
        }
        
        //this is mostly for the nextday analysis. If flasher analysis failed (all gains 0) we want to have all gains 1 instead of eventdisplay exiting with 499 dead pixels.
        //turn on with -nextdaygainhack. Not recommended for regular analysis.
        //this is a bit risky as we can miss 'bad' channels (gains won't be tested...).
        if( getRunParameter()->fNextDayGainHack && nZero > 200 )
        {
            cout << "VCalibrator::readGains(): Warning: " << nZero << " channels in Tel " <<  getTelID() + 1 << " had 0 gain. All rel. gains will be set to 1." << endl;
            setGains( 1., iLowGain );
            setGainvars( 1., iLowGain );
            
        }
        
        else if( !use_default )
        {
            if( VchannelList.size() == Vmean.size() && VchannelList.size() == Vvar.size() )
            {
                for( unsigned int i = 0; i < VchannelList.size(); i++ )
                {
                    if( VchannelList[i] < getGains().size() )
                    {
                        setGains( VchannelList[i], Vmean[i], iLowGain );
                        if( Vmean[i] > 0. )
                        {
                            getGainDist( iLowGain )->Fill( Vmean[i] );
                            getGainVarsDist( iLowGain )->Fill( Vvar[i] );
                        }
                        setGainvars( VchannelList[i], Vvar[i], iLowGain );
                    }
                    else
                    {
                        cout << "VCalibrator::readGains(): channel out of range: " << VchannelList[i] << "\t" << getGains().size() << endl;
                    }
                }
            }
        }
    }
    else if( getRunParameter()->fNoCalibNoPb || getRunParameter()->fIsMC )
    {
        if( getRunParameter()->fPrintSmallArray )
        {
            std::cout << "VCalibrator::readGains() info: Gains are set to 1, Gains are not tested to find dead channels " << std::endl;
        }
        setGains( 1., iLowGain );
        setGainvars( 1., iLowGain );
    }
    else
    {
        std::cout << "Error: No gain information available " << std::endl;
        std::cout << "exiting... " << std::endl;
        std::cout << "          please set option: -nocalibnoproblem , when launching EVNDISP if you want to continue anyway (all gains will be set to 1) " << std::endl;
        exit( EXIT_FAILURE );
    }
    
    // apply additional gain corrections
    getGains() /= getRunParameter()->fGainCorrection[getTelID()];
    
    
    if( getRunParameter()->freadCalibfromDB && !getRunParameter()->freadCalibfromDB_save_file )
    {
        char rm_calib_info_file[800];
        sprintf( rm_calib_info_file, "rm -rf  %s", iFile.c_str() );
        int syst_ret = system( rm_calib_info_file );
        if( syst_ret == -1 )
        {
            std::cout << "VCalibrator::readGains error removing calibration file " << iFile << endl;
        }
    }
    else if( getRunParameter()->freadCalibfromDB && getRunParameter()->freadCalibfromDB_save_file )
    {
        std::cout << "calibration information are stored in  " << iFile << std::endl;
    }
    
    
}

/*
 * read average pulse time (tzero or average pulse time) from root file
 *
 */
bool VCalibrator::readAverageTZeros( bool iLowGain )
{
    if( fDebug )
    {
        cout << "VCalibrator::readAverageTZeros (low gain " << iLowGain << ", method " << getSumWindowStart_T_method() << "): ";
        cout << fTZeroFileNameC[getTelID()] << endl;
    }
    
    string iFile;
    if( getTelID() < fTZeroFileNameC.size() )
    {
        iFile = fTZeroFileNameC[getTelID()];
    }
    if( iLowGain && getTelID() < fLowGainTZeroFileNameC.size() )
    {
        iFile = fLowGainTZeroFileNameC[getTelID()];
    }
    iFile += ".root";
    
    setAverageTZero( 0., iLowGain );
    setAverageTZerovars( 5., iLowGain );
    
    if( iFile.size() > 5 )
    {
        TFile iTZ( iFile.c_str(), "READ" );
        if( iTZ.IsZombie() )
        {
            cout << " VCalibrator::readAverageTZeros error reading average tzero file: " << endl;
            cout << "\t" << iFile << endl;
            return false;
        }
        char hname[200];
        if( getSumWindowStart_T_method() == 2 )
        {
            sprintf( hname, "tTAverage_%d", getTelID() + 1 );
        }
        else
        {
            sprintf( hname, "tTZero_%d", getTelID() + 1 );
        }
        TTree* t = ( TTree* )iTZ.Get( hname );
        if( !t )
        {
            cout << " VCalibrator::readAverageTZeros error finding data tree in average tzero file: " << endl;
            cout << "\t" << iFile << endl;
            if( !getRunParameter()->fNoCalibNoPb )
            {
                cout << "Exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            return false;
        }
        cout << "Telescope " << getTelID() + 1;
        if( !iLowGain )
        {
            cout << ": reading average pulse timing for high gain channels";
        }
        else
        {
            cout << ": reading average pulse timing for low gain channels";
        }
        cout << " from : " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << iFile << endl;
        
        int ichannel = 0;
        float itzero = 0.;
        float itzero_med = 0.;
        float itzero_var = 0.;
        t->SetBranchAddress( "channel", &ichannel );
        if( getSumWindowStart_T_method() == 2 )
        {
            t->SetBranchAddress( "TAverage", &itzero );
            t->SetBranchAddress( "TAveragemed", &itzero_med );
            t->SetBranchAddress( "TAveragevar", &itzero_var );
        }
        else
        {
            t->SetBranchAddress( "TZero", &itzero );
            t->SetBranchAddress( "TZeromed", &itzero_med );
            t->SetBranchAddress( "TZerovar", &itzero_var );
        }
        
        double n_T0 = 0.;
        double n_mean = 0.;
        for( int i = 0; i < t->GetEntries(); i++ )
        {
            t->GetEntry( i );
            
            // choose median of distribution as typical average
            // (distribution is assymetric)
            setAverageTZero( i, itzero_med, iLowGain );
            setAverageTZerovars( i, itzero_var, iLowGain );
            if( getAverageTZeroDist( iLowGain ) )
            {
                getAverageTZeroDist( iLowGain )->Fill( itzero_med );
            }
            if( itzero_med > 0. )
            {
                double d = getDetectorGeometry()->getOuterEdgeDistance( i );
                if( getTelID() < getDetectorGeometry()->getFieldofView().size() &&
                        d < 0.5 * getDetectorGeometry()->getFieldofView()[getTelID()] * getRunParameter()->faverageTZeroFiducialRadius )
                {
                    n_mean += itzero_med;
                    n_T0++;
                }
            }
        }
        if( n_T0 > 0. )
        {
            n_mean /= n_T0;
        }
        setMeanAverageTZero( n_mean, iLowGain );
        cout << "Telescope " << getTelID() + 1 << ": average pulse timing for this telescope is ";
        if( iLowGain )
        {
            cout << "(low gain)";
        }
        cout << getMeanAverageTZero( iLowGain );
        cout << " (calculated from " << n_T0 << " pixels";
        cout << ", pulse timing method " << getSumWindowStart_T_method();
        if( getTelID() < getDetectorGeometry()->getFieldofView().size() && getRunParameter()->faverageTZeroFiducialRadius > 0. &&
                getRunParameter()->faverageTZeroFiducialRadius < 1. )
        {
            cout << ", use inner " << getRunParameter()->faverageTZeroFiducialRadius * 100. << "% of camera";
        }
        cout << ")" << endl;
        return true;
    }
    
    return false;
    
}


void VCalibrator::readTOffsets( bool iLowGain )
{
    if( fDebug )
    {
        cout << "VCalibrator::readTOffsets" << endl;
    }
    
    string iFile = fToffFileNameC[getTelID()];
    if( iLowGain )
    {
        iFile = fLowGainToffFileNameC[getTelID()];
    }
    
    
    if( iFile.size() > 0 && getRunParameter()->frunmode != 2 )
    {
    
        vector< unsigned int > VchannelList;
        vector < double > Vmean;
        vector < double > Vvar;
        bool use_default = false;
        
        cout << "Telescope " << getTelID() + 1 << ":";
        cout << " reading time offsets";
        if( iLowGain )
        {
            cout << " for low gain channels ";
        }
        else
        {
            cout << " for high gain channels ";
        }
        cout << "from: " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        // read toffs from DB
        if( !iLowGain && getRunParameter()->freadCalibfromDB )
        {
            cout << "VOFFLINE DB" << endl;
            readfromVOFFLINE_DB( 2, iFile, VchannelList, Vmean, Vvar );
        }
        // read toffs from external file
        else
        {
            cout << iFile << endl;
            
            ifstream infile;
            infile.open( iFile.c_str(), ifstream::in );
            if( !infile )
            {
                cout << "VCalibrator::readTOffsets() warning: input file " << iFile << " cannot be opened." << endl;
                
                if( getRunParameter()->fDBRunType == "flasher" || getRunParameter()->fDBRunType == "laser" )
                {
                    cout << "Telescope " << getTelID() + 1 << ": ";
                    cout << "VCalibrator::readTOffsets() info: all tOffsets set to 0." << endl;
                    use_default = true;
                }
                else
                {
                    if( !getRunParameter()->fNoCalibNoPb )
                    {
                        std::cout << "Error: No toff information available " << std::endl;
                        std::cout << "exiting... " << std::endl;
                        std::cout << "          please set option: -nocalibnoproblem , when launching EVNDISP if you want to continue anyway (all TOffsets will be set to 0) " << std::endl;
                        exit( EXIT_FAILURE );
                        
                    }
                    cout << "Telescope " << getTelID() + 1 << ": ";
                    cout << "VCalibrator::readTOffsets() info: all tOffsets set to 0.\n" ;
                    use_default = true;
                }
                
                
                
            }
            
            char buffer[100];
            setTOffsets( 0., iLowGain );
            setTOffsetvars( 5., iLowGain );
            if( !use_default )
            {
                while( !infile.eof() )
                {
                    infile.getline( buffer, 100 );
                    unsigned int ch;
                    float mean, rms;
                    sscanf( buffer, "%u %f %f", &ch, &mean, &rms );
                    VchannelList.push_back( ch );
                    Vmean.push_back( mean );
                    Vvar.push_back( rms );
                }
            }
        }
        
        if( !use_default )
        {
            if( VchannelList.size() == Vmean.size() && VchannelList.size() == Vvar.size() )
            {
                for( unsigned int i = 0; i < VchannelList.size(); i++ )
                {
                    if( VchannelList[i] < getTOffsets().size() )
                    {
                        setTOffsets( VchannelList[i], Vmean[i], iLowGain );
                        if( Vmean[i] != 0. )
                        {
                            getToffsetDist( iLowGain )->Fill( Vmean[i] );
                            getToffsetVarsDist( iLowGain )->Fill( Vvar[i] );
                        }
                        if( VchannelList[i] < getTOffsetvars( iLowGain ).size() )
                        {
                            setTOffsetvars( VchannelList[i], Vvar[i], iLowGain );
                        }
                    }
                    else
                    {
                        cout << "Telescope " << getTelID() + 1 << ": ";
                        cout << "VCalibrator::readTOffsets(): channel out of range: " << VchannelList[i] << "\t" << getTOffsets().size() << endl;
                    }
                }
            }
        }
    }
    else if( getRunParameter()->fNoCalibNoPb || iLowGain || getRunParameter()->fIsMC )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        std::cout << "VCalibrator::readTOffsets() info: TOffsets are set to 0";
        if( iLowGain )
        {
            std::cout << " (low gain)";
        }
        cout << std::endl;
        if( !iLowGain )
        {
            cout << "Telescope " << getTelID() + 1 << ": ";
            std::cout << "VCalibrator::readTOffsets() info: TOffsets are not tested to find dead channels " << std::endl;
        }
        setTOffsets( 0., iLowGain );
        setTOffsetvars( 0.1, iLowGain );
    }
    else
    {
    
        std::cout << "Error: No toff information available " << std::endl;
        std::cout << "exiting... " << std::endl;
        std::cout << "          please set option: -nocalibnoproblem , when launching EVNDISP if you want to continue anyway (all TOffsets will be set to 0)" << std::endl;
        exit( EXIT_FAILURE );
        
    }
    
    
    // check if toffs are saved to disk
    if( getRunParameter()->freadCalibfromDB && !getRunParameter()->freadCalibfromDB_save_file )
    {
        char rm_calib_info_file[800];
        sprintf( rm_calib_info_file, "rm -rf  %s", iFile.c_str() );
        int syst_ret = system( rm_calib_info_file );
        if( syst_ret == -1 )
        {
            std::cout << "VCalibrator::readTOffsets error removing calibration file " << iFile << endl;
        }
    }
    else if( getRunParameter()->freadCalibfromDB && getRunParameter()->freadCalibfromDB_save_file )
    {
        std::cout << "calibration information are stored in  " << iFile << std::endl;
    }
    
}


void VCalibrator::initialize()
{
    if( fDebug )
    {
        cout << "VCalibrator::initialize()" << endl;
    }
    
    // set the data readers
    initializeDataReader();
    
    // get all relevation run numbers (e.g. from tzero, ped files)
    // for calibration settings
    getCalibrationRunNumbers();
    
    ////////////////////////////////////////////////////////////////
    // create the calibration data structures
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        setTelID( i );
        
        fCalData.push_back( new VCalibrationData( i,
                            fPedFileNameC[i], fGainFileNameC[i], fToffFileNameC[i], fLowGainPedFileNameC[i],
                            "", "", "", fTZeroFileNameC[i], fLowGainTZeroFileNameC[i],
                            getRunParameter()->getObservatory() ) );
        fCalData.back()->setSumWindows( getSumWindow( i ) );
        fNumberGainEvents.push_back( 0 );
        fNumberTZeroEvents.push_back( 0 );
        
        fCalData.back()->initialize( getNChannels(), getNSamples(), usePedestalsInTimeSlices( false ), usePedestalsInTimeSlices( true ),
                                     ( fReader->getDataFormat() == "grisu"
                                       || getRunParameter()->frunmode == 1 || getRunParameter()->frunmode == 6
                                       || ( getRunParameter()->fsourcetype == 2 && getRunParameter()->fsimu_pedestalfile.size() > 0 ) ),
                                     getRunParameter()->freadCalibfromDB,
                                     isDST_MC(),
                                     getDebugFlag(), getRunParameter()->frunmode, isTeltoAna( i ) );
        fCalData.back()->setReader( fReader );
    }
    
    //////////////////////////////////////////////////////////////////////////
    // define histograms and output files for pedestal calculation
    // (only done for pedestal calculation (rumode = 1 or runmode = 6)
    vector<ULong64_t> i_TelTypeList = getDetectorGeometry()->getTelType_list();
    if( fRunPar->frunmode == 1 || fRunPar->frunmode == 6 )
    {
        for( unsigned int t = 0; t < i_TelTypeList.size(); t++ )
        {
            // set telescope ID
            unsigned int iTelID = 99999;
            for( unsigned int q = 0; q < getDetectorGeometry()->getTelType().size(); q++ )
            {
                if( getDetectorGeometry()->getTelType()[q] == i_TelTypeList[t] )
                {
                    iTelID = q;
                    break;
                }
            }
            if( iTelID == 99999 )
            {
                cout << "VCalibrator::initialize: error: no telID found for telescope type " << i_TelTypeList[t] << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            setTelID( iTelID );
            
            unsigned int i_nclones = 0;
            cout << "VCalibrator::initialize: setting ped histos and files for telescope type " << i_TelTypeList[t] << endl;
            fPedSingleOutFile = 0;
            if( !getRunParameter()->fPedestalSingleRootFile )
            {
                fPedOutFile[i_TelTypeList[t]] = 0;
            }
            vector< vector< TH1F* > > iped_vec;
            vector< vector< TH1F* > > iped_ttype_vec;
            //vector< TH1F* > iped_ttype_vec;
            for( int i = 0; i < getRunParameter()->fCalibrationSumWindow; i++ )
            {
                vector<TH1F* > ihped;
                for( unsigned int j = 0; j < getNChannels(); j++ )
                {
                    ihped.push_back( 0 );
                }
                iped_vec.push_back( ihped );
                i_nclones += ihped.size();
                
                vector<TH1F* > ihped_ttype;
                for( unsigned int j = 0; j < getNChannels(); j++ )
                {
                    ihped_ttype.push_back( 0 );
                }
                iped_ttype_vec.push_back( ihped_ttype );
                i_nclones += ihped_ttype.size();
            }
            
            hped_vec[i_TelTypeList[t]] = iped_vec;
            hpedPerTelescopeType[i_TelTypeList[t]] = iped_ttype_vec;
            // number of pedestal events
            fNumberPedestalEvents[i_TelTypeList[t]] = 0;
            // create a TClonesArray for all pedestals histograms
            fPedestalsHistoClonesArray[i_TelTypeList[t]] = new TClonesArray( "TH1F", i_nclones );
        }
        
    }
    ////////////////////////////////////////////////////////////////
    // read the calibration files
    if( fRunPar->fsourcetype != 6 && fRunPar->fsourcetype != 7 && fRunPar->fsourcetype != 4 )
    {
        // don't read data for pedestal calculation step
        if( fRunPar->frunmode != 1 && fRunPar->frunmode != 6 )
        {
            readCalibrationData();
        }
    }
    // PE mode: set gains to 1
    else
    {
        for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
        {
            int t = getTeltoAna()[tel];
            setTelID( t );
            readGains( false );
            readGains( true );
        }
        // read peds from DST file (for DST MC source file)
        
        if( fRunPar->fsourcetype == 7 || fRunPar->fsourcetype == 4 )
        {
            readCalibrationDatafromDSTFiles( fRunPar->fsourcefile, fRunPar->frunmode == 1 );
        }
    }
    // if needed: write IPR graphs to disk
    if( getRunParameter()->ifCreateIPRdatabase == true && getRunParameter()->ifReadIPRfromDatabase == false )
    {
        writeIPRgraphs();
    }
    
    // initialize dead  channel finder
    initializeDeadChannelFinder();
    
    // read IPR graphs used in time next-neighbour cleaning
    if( getImageCleaningParameter()
            && getImageCleaningParameter()->getImageCleaningMethod() == "TIMENEXTNEIGHBOUR" )
    {
        if( fRunPar->fsourcetype != 6
                && fRunPar->fsourcetype != 7
                && fRunPar->fsourcetype != 4
                && fRunPar->frunmode != 1
                && fRunPar->frunmode != 6 )
        {
            calculateIPRGraphs();
        }
    }
}

void VCalibrator::setCalibrationFileNames()
{
    if( fDebug )
    {
        cout << "VCalibrator::setCalibrationFileNames()" << endl;
    }
    
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( i < fPedFileNameC.size() )
        {
            fPedFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fPedFileNumber[i], "ped" );
        }
        else
        {
            fPedFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fPedFileNumber[i], "ped" ) );
        }
        if( i < fGainFileNameC.size() )
        {
            fGainFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fGainFileNumber[i], "gain" );
        }
        else
        {
            fGainFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fGainFileNumber[i], "gain" ) );
        }
        if( i < fToffFileNameC.size() )
        {
            fToffFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fTOffFileNumber[i], "toff" );
        }
        else
        {
            fToffFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fTOffFileNumber[i], "toff" ) );
        }
        if( i < fTZeroFileNameC.size() )
        {
            fTZeroFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fTZeroFileNumber[i], "tzero" );
        }
        else
        {
            fTZeroFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fTZeroFileNumber[i], "tzero" ) );
        }
        if( i < fPixFileNameC.size() )
        {
            fPixFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fPixFileNumber[i], "pix" );
        }
        else
        {
            fPixFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fPixFileNumber[i], "pix" ) );
        }
        if( i >= fBlockTel.size() )
        {
            fBlockTel.push_back( false );
        }
        if( i < fLowGainPedFileNameC.size() )
        {
            fLowGainPedFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fPedLowGainFileNumber[i], "lped" );
        }
        else
        {
            fLowGainPedFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fPedLowGainFileNumber[i], "lped" ) );
        }
        
        if( i < fNewLowGainPedFileNameC.size() )
        {
            fNewLowGainPedFileNameC[i] = getCalibrationFileName( i, 10, "newlped" , getRunParameter()->fPedLowGainFile );
        }
        else
        {
            fNewLowGainPedFileNameC.push_back( getCalibrationFileName( i, 10, "newlped", getRunParameter()->fPedLowGainFile ) );
        }
        
        if( i < fLowGainMultiplierNameC.size() )
        {
            fLowGainMultiplierNameC[i] = getCalibrationFileName( i, getRunParameter()->fLowGainMultiplierFileNumber[i], "lmult" );
        }
        else
        {
            fLowGainMultiplierNameC.push_back( getCalibrationFileName( i, getRunParameter()->fLowGainMultiplierFileNumber[i], "lmult" ) );
        }
        if( i < fLowGainGainFileNameC.size() )
        {
            fLowGainGainFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fGainLowGainFileNumber[i], "lgain" );
        }
        else
        {
            fLowGainGainFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fGainLowGainFileNumber[i], "lgain" ) );
        }
        if( i < fLowGainToffFileNameC.size() )
        {
            fLowGainToffFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fTOffLowGainFileNumber[i], "ltoff" );
        }
        else
        {
            fLowGainToffFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fTOffLowGainFileNumber[i], "ltoff" ) );
        }
        if( i < fLowGainTZeroFileNameC.size() )
        {
            fLowGainTZeroFileNameC[i] = getCalibrationFileName( i, getRunParameter()->fTZeroLowGainFileNumber[i], "tzero" );
        }
        else
        {
            fLowGainTZeroFileNameC.push_back( getCalibrationFileName( i, getRunParameter()->fTZeroLowGainFileNumber[i], "tzero" ) );
        }
    }
}

/*

    getting calibration run numbers

    ped, gains, toffs (high and low gain)
*/
void VCalibrator::getCalibrationRunNumbers()
{
    if( fDebug )
    {
        cout << "VCalibrator::getCalibrationRunNumbers()" << endl;
    }
    
    // initialize vectors
    setCalibrationFileNames();
    
    // reading peds etc from grisu files
    if( getRunParameter()->fsimu_pedestalfile.size() > 0 )
    {
        cout << "VCalibrator::getCalibrationRunNumbers() info: taking calibration from grisu files" << endl;
        return;
    }
    
    // get calibration run numbers from calibration text file
    int iCaliLines = 0;
    if( getRunParameter()->fcalibrationfile.size() > 0 )
    {
        iCaliLines = getCalibrationRunNumbers_fromCalibFile();
        if( iCaliLines > 0 )
        {
            setCalibrationFileNames();
        }
    }
    int iLowGainCaliLines = 0;
    if( getRunParameter()->fLowGainCalibrationFile.size() > 0 && getRunParameter()->frunmode != 6 && !getRunParameter()->fIsMC )
    {
        iLowGainCaliLines = readLowGainCalibrationValues_fromCalibFile( "LOWGAINPED" ) + readLowGainCalibrationValues_fromCalibFile( "LOWGAINPEDFILE" ) ;
        if( iLowGainCaliLines > 0 )
        {
            setCalibrationFileNames();
        }
    }
    
    // take pedestals from grisu output file ('P'-lines), gains=1, and toff = 0.
    if( iCaliLines == 0 && ( fReader->getDataFormat() == "grisu" || getRunParameter()->fsimu_pedestalfile.size() > 0 ) )
    {
        cout << "VCalibrator::getCalibrationRunNumbers() info: taking calibration from grisu files" << endl;
        return;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // no information found for this runnumber
    if( iCaliLines == 0 )
    {
        // for dst, take dst.ped
        if( fReader->getDataFormat() == "DST" || fReader->getDataFormat() == "PE" )
        {
            for( unsigned int i = 0; i < getNTel(); i++ )
            {
                fPedFileNameC[i] = "dst";
            }
        }
        // plot raw: don't need any calibration data
        else if( getRunParameter()->fPlotRaw )
        {
            for( unsigned int i = 0; i < getNTel(); i++ )
            {
                fPedFileNameC[i] = "";
                fLowGainPedFileNameC[i] = "";
            }
        }
        // calibration run numbers are taken from DB
        else if( getRunParameter()->fsimu_pedestalfile.size() == 0 &&
                 getRunParameter()->fuseDB && getRunParameter()->frunmode != 2 &&
                 getRunParameter()->frunmode != 5 &&
                 getRunParameter()->frunmode != 1 &&
                 getRunParameter()->frunmode != 6 )
        {
            cout << "VCalibrator::getCalibrationRunNumbers() info: using laser file number from DB" << endl;
        }
        // for laser files: take data file run number for pedestals and gains
        else if( getRunParameter()->fDBRunType == "flasher" )
        {
            cout << "VCalibrator::getCalibrationRunNumbers() info: laser/flasher files, using run number for pedestals and gains" << endl;
        }
        // this is an calibration run
        else if( getRunParameter()->fsimu_pedestalfile.size() == 0 && getRunParameter()->fcalibrationrun )
        {
            cout << "VCalibrator::getCalibrationRunNumbers() info: calibration run, using run number for pedestals, tzeros and gains" << endl;
        }
        // this is a MC run
        else if( getRunParameter()->fsimu_pedestalfile.size() == 0 && getRunParameter()->fIsMC )
        {
            cout << "VCalibrator::getCalibrationRunNumbers() info: MC run, using run number for pedestals, tzeros and gains" << endl;
        }
        // no calibration data: problem
        else
        {
            cout << "VCalibrator::getCalibrationRunNumbers() error: no calibration information found for run ";
            cout << getRunNumber();
            if( getRunParameter()->fcalibrationfile.size() > 0 )
            {
                cout << " in file " << getRunParameter()->fcalibrationfile <<  endl;
            }
            else
            {
                cout << ": no file given " << endl;
            }
            cout << "(run type " << getRunParameter()->fDBRunType << ")" << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // append suffixes
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fPixFileNameC[i].size() > 0 )
        {
            fPixFileNameC[i] += ".pix";
        }
    }
    cout << endl;
    
    return;
}

/*
 *  read hilo multiplier from a calibration file for a telescope telescope
 *
 *  the default location for these files is:
 *  $VERITAS_EVNDISP_AUX_DIR/Calibration/calibrationlist.LowGain.dat
 *
 */
int VCalibrator::readLowGainCalibrationValues_fromCalibFile( string iVariable, unsigned int iTelescopeSelect )
{
    int iLinesFound = 0;
    // open file with calibration information
    string is_Temp;
    is_Temp = fRunPar->getDirectory_EVNDISPCalibrationData();
    is_Temp += getRunParameter()->fLowGainCalibrationFile;
    ifstream is;
    is.open( is_Temp.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VCalibrator::readLowGainCalibrationValues_fromCalibFile: error, calibration data file not found: ";
        cout <<  is_Temp << endl;
        exit( EXIT_FAILURE );
    }
    if( iTelescopeSelect != 9999 )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "reading low-gain values (" << iVariable << ") ";
        cout << "from: " << endl;
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << is_Temp << endl;
    }
    
    string is_line;
    string is_tempstring;
    int iRunMin = 0;
    int iRunMax = 0;
    int iTel = -1;
    
    string iValueString;
    for( unsigned int i = 0; i < fBlockTel.size(); i++ )
    {
        fBlockTel[i] = false;
    }
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            // only line with '*' as first characters are valid lines
            is_stream >> is_Temp;
            if( is_Temp != "*" )
            {
                continue;
            }
            
            if( (is_stream>>std::ws).eof() )
            {
                continue;
            }
            is_stream >> is_Temp;
            if( is_Temp == iVariable )
            {
                if( (is_stream>>std::ws).eof() )
                {
                    continue;
                }
                
                if( iVariable == "LOWGAINPEDFILE" )
                {
                    is_stream >> is_tempstring;
                    getRunParameter()->fPedLowGainFile = is_tempstring ;
                    iLinesFound++;
                    continue;
                }
                
                is_stream >> iTel;
                iTel--;              // internal counting starts at 0
                if( (is_stream>>std::ws).eof() )
                {
                    continue;
                }
                is_stream >> iRunMin;
                if( (is_stream>>std::ws).eof() )
                {
                    continue;
                }
                is_stream >> iRunMax;
                // check run number range
                if( getRunNumber() >= iRunMin && getRunNumber() <= iRunMax )
                {
                    if( iVariable == "LOWGAINPED" )
                    {
                        cout << "reading low-gain parameters for run range " << iRunMin << ", " << iRunMax << endl;
                    }
                    if( (is_stream>>std::ws).eof() )
                    {
                        continue;
                    }
                    is_stream >> iValueString;
                    iLinesFound++;
                    if( iTel < 0 )
                    {
                        for( unsigned int i = 0; i < getNTel(); i++ )
                        {
                            if( !fBlockTel[i] )
                            {
                                if( iVariable == "LOWGAINPED" && atoi( iValueString.c_str() ) > 0 )
                                {
                                    getRunParameter()->fPedLowGainFileNumber[i] = atoi( iValueString.c_str() );
                                }
                                else if( iVariable == "LOWGAINMULTIPLIER_TRACE" && atof( iValueString.c_str() ) > 0. )
                                {
                                    setLowGainMultiplier_Trace( atof( iValueString.c_str() ) );
                                }
                                else if( iVariable == "LOWGAINMULTIPLIER_SUM" )
                                {
                                    cout << "VCalibrator::readLowGainCalibrationValues_fromCalibFile() warning: ";
                                    cout << "sumwindow-dependent low-gain multipliers need to be set for each telescope " << endl;
                                    cout << "(ignored this entry)" << endl;
                                }
                            }
                        }
                    }
                    else if( iTel < ( int )getNTel() )
                    {
                        if( iVariable == "LOWGAINPED" && atoi( iValueString.c_str() ) > 0 )
                        {
                            getRunParameter()->fPedLowGainFileNumber[iTel] = atoi( iValueString.c_str() );
                            fBlockTel[iTel] = true;
                        }
                        else if( iVariable == "LOWGAINMULTIPLIER_TRACE" && atof( iValueString.c_str() ) > 0. && iTel == ( int )iTelescopeSelect )
                        {
                            cout << "Telescope " << getTelID() + 1 << ": ";
                            cout << "VCalibrator::readLowGainMultiplier():  LOWGAINMULTIPLIER_TRACE " << iTel + 1 << ": ";
                            cout <<  atof( iValueString.c_str() ) << endl;
                            setLowGainMultiplier_Trace( atof( iValueString.c_str() ) );
                        }
                        else if( iVariable == "LOWGAINMULTIPLIER_SUM" && atoi( iValueString.c_str() ) > 0. && iTel == ( int )iTelescopeSelect )
                        {
                            // low-gain multiplier depend on the trace integration method;
                            unsigned int iTraceIntegrationMethod = ( unsigned int )atoi( iValueString.c_str() );
                            is_stream >> iValueString;
                            int isum = atoi( iValueString.c_str() );
                            std::pair<int, int> iPair( iTraceIntegrationMethod, isum );
                            getLowGainDefaultSumWindows().push_back( iPair );
                            int jmin = 0;
                            int jmax = 0;
                            int jsum = 0;
                            double lmult;
                            if( !( is_stream >> jmin >> jmax ) )
                            {
                                continue;
                            }
                            jsum = jmin;
                            while( is_stream >> lmult )
                            {
                                setLowGainMultiplier_Sum( iTraceIntegrationMethod, isum, jsum, lmult );
                                jsum++;
                            }
                            if( jsum - 1 != jmax )
                            {
                                cout << "Telescope " << getTelID() + 1 << ": ";
                                cout << "VCalibrator::readLowGainCalibrationValues_fromCalibFile: warning,";
                                cout << "you claim to have low gain multpliers from sw " << jmin << " to " << jmax ;
                                cout << ", but the last value read in corresponds to sw " << jsum - 1 << ". Please check your config file, line " << endl;
                                cout << is_line << endl;
                            }
                            
                        }
                    }
                }
            }
        }
    }
    
    return iLinesFound;
}



int VCalibrator::getCalibrationRunNumbers_fromCalibFile()
{
    // open file with calibration information
    string is_Temp;
    is_Temp = fRunPar->getDirectory_EVNDISPCalibrationData();
    is_Temp += getRunParameter()->fcalibrationfile;
    ifstream is;
    is.open( is_Temp.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VCalibrator::getCalibrationRunNumbers error, calibration data file not found: ";
        cout <<  is_Temp << endl;
        exit( EXIT_FAILURE );
    }
    
    string is_line;
    string iPed, iPad, iGain, iToff, iPix, iLowGainPeds, iLowGainMultiplier, iLowGainGains, iLowGainToff;
    int iTel = 0;
    int iCaliLines = 0;
    bool bReset = false;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            // only line with '*' as first characters are valid lines
            is_stream >> is_Temp;
            if( is_Temp != "*" && is_Temp != "*V1" && is_Temp != "*V2" && is_Temp != "*V3" && is_Temp != "*V4" )
            {
                continue;
            }
            
            // get version number
            if( is_Temp == "*V2" )
            {
                fCalibrationfileVersion = 2;
            }
            else if( is_Temp == "*V3" )
            {
                fCalibrationfileVersion = 3;
            }
            else if( is_Temp == "*V4" )
            {
                fCalibrationfileVersion = 4;
            }
            
            // get and check run number
            is_stream >> is_Temp;
            if( atoi( is_Temp.c_str() ) != getRunNumber() )
            {
                continue;
            }
            iCaliLines++;
            
            cout << "reading calibration information (version " << fCalibrationfileVersion << ")" << endl;
            cout << is_line << endl;
            
            // get telescope numbers
            is_stream >> is_Temp;
            iTel = atoi( is_Temp.c_str() ) - 1;
            
            // get ped file number
            is_stream >> iPed;
            
            // get gain file number
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iGain;
            }
            
            // get toff file number
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iToff;
            }
            
            // get pixel status number
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iPix;
            }
            
            // get pad file number
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iPad;
            }
            
            // get low gain pedestal file number
            if( fCalibrationfileVersion > 1 )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iLowGainPeds;
                }
                if( bReset )
                {
                    getRunParameter()->fPedLowGainFileNumber.clear();
                    for( unsigned int i = 0; i < getNTel(); i++ )
                    {
                        getRunParameter()->fPedLowGainFileNumber.push_back( 0 );
                    }
                    bReset = true;
                }
            }
            else
            {
                iLowGainPeds = "-2";
            }
            
            // low gain gains and time offsets
            if( fCalibrationfileVersion > 3 )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iLowGainGains;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iLowGainToff;
                }
            }
            
            // get low gain multiplier file name
            if( fCalibrationfileVersion > 2 )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iLowGainMultiplier;
                }
            }
            
            // iTel < 0: valid for all telescopes
            if( iTel < 0 )
            {
                for( unsigned int i = 0; i < getNTel(); i++ )
                {
                    if( !fBlockTel[i] )
                    {
                        if( iPed != "-1" )
                        {
                            getRunParameter()->fPedFileNumber[i] = atoi( iPed.c_str() );
                        }
                        if( iGain != "-1" )
                        {
                            getRunParameter()->fGainFileNumber[i] = atoi( iGain.c_str() );
                        }
                        if( iToff != "-1" )
                        {
                            getRunParameter()->fTOffFileNumber[i] = atoi( iToff.c_str() );
                        }
                        if( iPix != "-1" )
                        {
                            getRunParameter()->fPixFileNumber[i] = atoi( iPix.c_str() );
                        }
                        if( fCalibrationfileVersion > 1 && iLowGainPeds != "-1" )
                        {
                            getRunParameter()->fPedLowGainFileNumber[i] = atoi( iLowGainPeds.c_str() );
                        }
                        else if( iLowGainPeds == "-2"
                                 && i < getRunParameter()->fPedLowGainFileNumber.size() && getRunParameter()->fPedLowGainFileNumber[i] > 0 )
                        {
                            ostringstream s_temp;
                            s_temp << getRunParameter()->fPedLowGainFileNumber[i];
                            fLowGainPedFileNameC[i] = s_temp.str();
                        }
                        else if( iLowGainPeds == "-1" )
                        {
                            getRunParameter()->fPedLowGainFileNumber[i] = atoi( iLowGainPeds.c_str() );
                        }
                        if( fCalibrationfileVersion > 2 && iLowGainMultiplier != "-1" )
                        {
                            fLowGainMultiplierNameC[i] = iLowGainMultiplier;
                        }
                        if( fCalibrationfileVersion > 3 && iLowGainGains != "-1" )
                        {
                            fLowGainGainFileNameC[i] = iLowGainGains;
                        }
                        if( fCalibrationfileVersion > 3 && iLowGainToff != "-1" )
                        {
                            getRunParameter()->fTOffLowGainFileNumber[i] = atoi( iLowGainToff.c_str() );
                        }
                    }
                }
            }
            else if( iTel < ( int )getNTel() )
            {
                if( iPed != "-1" )
                {
                    getRunParameter()->fPedFileNumber[iTel] = atoi( iPed.c_str() );
                }
                if( iGain != "-1" )
                {
                    getRunParameter()->fGainFileNumber[iTel] = atoi( iGain.c_str() );
                }
                if( iToff != "-1" )
                {
                    getRunParameter()->fTOffFileNumber[iTel] = atoi( iToff.c_str() );
                }
                if( iPix != "-1" )
                {
                    getRunParameter()->fPixFileNumber[iTel] = atoi( iPix.c_str() );
                }
                else
                {
                    fPixFileNameC[iTel] = "";
                }
                if( fCalibrationfileVersion > 1 && iLowGainPeds != "-1" )
                {
                    getRunParameter()->fPedLowGainFileNumber[iTel] = atoi( iLowGainPeds.c_str() );
                }
                if( fCalibrationfileVersion > 2 && iLowGainMultiplier != "-1" )
                {
                    fLowGainMultiplierNameC[iTel] = iLowGainMultiplier;
                }
                if( fCalibrationfileVersion > 3 && iLowGainGains != "-1" )
                {
                    getRunParameter()->fGainLowGainFileNumber[iTel] = atoi( iLowGainGains.c_str() );
                }
                if( fCalibrationfileVersion > 3 && iLowGainToff != "-1" )
                {
                    getRunParameter()->fTOffLowGainFileNumber[iTel] = atoi( iLowGainToff.c_str() );
                }
                fBlockTel[iTel] = true;
            }
        }
    }
    is.close();
    
    return iCaliLines;
}

void VCalibrator::readPixelstatus()
{
    if( fDebug )
    {
        cout << "VCalibrator::readPixelstatus()" << endl;
    }
    
    if( fPixFileNameC[getTelID()].size() > 0 )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "reading pixel status from ";
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << fPixFileNameC[getTelID()] << endl;
        // open pedestal file
        ifstream is;
        is.open( fPixFileNameC[getTelID()].c_str(), ifstream::in );
        if( !is )
        {
            cout << "VCalibrator::readPixelstatus error: unable to open pixel file " << fPixFileNameC[getTelID()] << endl;
            cout << "\t exiting...." << endl;
            exit( EXIT_FAILURE );
        }
        string is_line;
        string iTemp;
        string iTemp2;
        while( getline( is, is_line ) )
        {
            if( is_line.size() > 0 )
            {
                istringstream is_stream( is_line );
                is_stream >> iTemp;
                is_stream >> iTemp2;
                if( atoi( iTemp.c_str() ) < ( int )getChannelStatus().size() )
                {
                    getChannelStatus()[atoi( iTemp.c_str() )] = atoi( iTemp2.c_str() );
                }
                else
                {
                    cout << "VCalibrator::readPixelstatus() warning: index out of range " << iTemp << "\t";
                    cout << getChannelStatus().size() << endl;
                }
            }
        }
        is.close();
    }
}

/*!

     get low gain multiplier (might be summation window dependent)

*/
bool VCalibrator::readLowGainMultiplier( )
{

    //////////////////////////////////////////////////////////////////
    // read low gain multiplier from low gain calibration file
    // (one value per telescope & sumwindow)
    if( getRunParameter()->fLowGainCalibrationFile.size() > 0 && getRunParameter()->frunmode != 6 )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "reading low gain multiplier from " << getRunParameter()->fLowGainCalibrationFile << endl;
        readLowGainCalibrationValues_fromCalibFile( "LOWGAINMULTIPLIER_TRACE", getTelID() );
        readLowGainCalibrationValues_fromCalibFile( "LOWGAINMULTIPLIER_SUM", getTelID() );
        getDetectorGeometry()->setLowGainMultiplier_Trace( getTelID(), getCalData()->getLowGainMultiplier_Trace() );
    }
    //////////////////////////////////////////////////////////////////
    // fill low gain multipliers set in .cfg file
    // (one value per telescope)
    else if( getDetectorGeometry()->isLowGainSet() )
    {
        cout << "Telescope " << getTelID() + 1 << ": ";
        cout << "reading low gain multiplier from cfg file: " << getDetectorGeometry()->getLowGainMultiplier_Trace()[getTelID()] << endl;
        setLowGainMultiplier_Trace( getDetectorGeometry()->getLowGainMultiplier_Trace()[getTelID()] );
    }
    return true;
}


/*

    read calibration data from DST file

*/
bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile, bool iPedOnly )
{
    if( iDSTfile.size() == 0 )
    {
        cout << "VCalibrator::readCalibrationDatafromDSTFiles: no DST file given" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    TFile iF( iDSTfile.c_str() );
    if( iF.IsZombie() )
    {
        cout << "VCalibrator::readCalibrationDatafromDSTFiles error while opening DST tree in " << iDSTfile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    TTree* t = ( TTree* )iF.Get( "calibration" );
    if( !t )
    {
        cout << "VCalibrator::readCalibrationDatafromDSTFiles warning: failed reading calibration tree from file " << iDSTfile << endl;
        return false;
    }
    cout << "reading calibration data from dst file" << endl;
    
    int fTelID = 0;
    unsigned int nPixel = 0;
    unsigned int fnum_sumwindow = 1;
    unsigned int fsumwindow[VDST_MAXSUMWINDOW];
    float fPed_high[VDST_MAXCHANNELS];
    float* fPedvar_high = 0;
    if( !iPedOnly )
    {
        fPedvar_high = new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
    }
    float fPed_low[VDST_MAXCHANNELS];
    float* fPedvar_low = 0;
    if( !iPedOnly )
    {
        fPedvar_low = new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
    }
    float fConv_high[VDST_MAXCHANNELS];
    float fConv_low[VDST_MAXCHANNELS];
    float ftzero[VDST_MAXCHANNELS];
    
    for( unsigned int i = 0; i < VDST_MAXCHANNELS; i++ )
    {
        ftzero[i] = -99.;
    }
    
    t->SetBranchAddress( "TelID", &fTelID );
    t->SetBranchAddress( "NPixel", &nPixel );
    if( t->GetBranchStatus( "num_sumwindow" ) )
    {
        t->SetBranchAddress( "num_sumwindow", &fnum_sumwindow );
        t->SetBranchAddress( "sumwindow", fsumwindow );
    }
    else
    {
        fnum_sumwindow = 0;
    }
    t->SetBranchAddress( "ped_high", fPed_high );
    // check if this is a old or new style pedvar tree
    bool iPedVarTreeTypeNew = true;
    if( !iPedOnly && t->GetLeaf( "pedvar_high" ) )
    {
        string i_temp = ( string )t->GetLeaf( "pedvar_high" )->GetTitle();
        if( i_temp.find_first_of( "[" ) == i_temp.find_last_of( "[" ) )
        {
            iPedVarTreeTypeNew = true;
        }
        else
        {
            iPedVarTreeTypeNew = false;
        }
    }
    if( !iPedVarTreeTypeNew )
    {
        cout << "VCalibrator::readCalibrationDatafromDSTFiles: old style DST file: warning: pedvars will not be set correctly" << endl;
        cout << "   (you might be able to proceed if you don't do FADC trace analysis and use the pedvars for the image cleaning" << endl;
    }
    t->SetBranchAddress( "pedvar_high", fPedvar_high );
    t->SetBranchAddress( "ped_low", fPed_low );
    if( !iPedOnly )
    {
        t->SetBranchAddress( "pedvar_low", fPedvar_low );
    }
    t->SetBranchAddress( "conv_high", fConv_high );
    t->SetBranchAddress( "conv_low", fConv_low );
    if( t->GetBranch( "tzero" ) )
    {
        t->SetBranchAddress( "tzero", ftzero );
    }
    
    // reset histograms with pedestal distributions
    if( getPedvarsDist() )
    {
        getPedvarsDist()->Reset();
    }
    if( getPedDist() )
    {
        getPedDist()->Reset();
    }
    if( getPedvarsLowGainDist() )
    {
        getPedvarsLowGainDist()->Reset();
    }
    if( getPedLowGainDist() )
    {
        getPedLowGainDist()->Reset();
    }
    
    // histogram for average tzero calculation (per telescope type)
    vector<ULong64_t> iTelType = getDetectorGeometry()->getTelType_list();
    vector< TH1D* > iH_averageTZero;
    char hname[200];
    for( unsigned int i = 0; i < iTelType.size(); i++ )
    {
        sprintf( hname, "iH_averageTZero_%d", i );
        if( getNSamples() > 0 )
        {
            iH_averageTZero.push_back( new TH1D( hname, "", 3 * ( int )getNSamples(), 0., ( float )getNSamples() ) );
        }
        else
        {
            iH_averageTZero.push_back( new TH1D( hname, "", 100., 0., 100. ) );
        }
    }
    
    if( getNTel() != ( unsigned int )t->GetEntries() )
    {
        cout << "VCalibrator::readCalibrationDatafromDSTFiles error: mismatch in number of telescopes: " << endl;
        cout << "\t expected: " << getNTel() << endl;
        cout << "\t found: " << t->GetEntries() << endl;
        exit( EXIT_FAILURE );
    }
    // loop over all telescopes
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        setTelID( i );
        
        // no calibration data available for this telescope
        if( nPixel == 0 )
        {
            continue;
        }
        
        // pedestals
        if( nPixel == getPeds( false ).size() )
        {
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                getPeds( false )[p] = fPed_high[p];
                if( getPedDist() )
                {
                    getPedDist()->Fill( fPed_high[p] );
                }
            }
        }
        else
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile )";
            cout << " error: index out of range (peds, high gain): ";
            cout << nPixel << "\t" << getPeds( false ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        if( nPixel == getPeds( true ).size() )
        {
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                getPeds( true )[p] = fPed_low[p];
                if( getPedLowGainDist() )
                {
                    getPedLowGainDist()->Fill( fPed_low[p] );
                }
            }
        }
        else
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile )";
            cout << "error: index out of range (peds, low gain): ";
            cout << nPixel << "\t" << getPeds( true ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        // pedvars
        if( nPixel == getPedvars( false ).size() && !iPedOnly )
        {
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( p < getPedvars( false ).size() )
                {
                    for( unsigned int s = 0; s < fnum_sumwindow; s++ )
                    {
                        if( iPedVarTreeTypeNew )
                        {
                            getPedvars( false, fsumwindow[s] )[p] = fPedvar_high[p * VDST_MAXSUMWINDOW + s];
                            if( fsumwindow[s] == 1 )
                            {
                                getPedrms( false ) = fPedvar_high[p * VDST_MAXSUMWINDOW + s];
                            }
                        }
                        if( fsumwindow[s] == getSumWindow() && getPedvarsDist() )
                        {
                            if( iPedVarTreeTypeNew )
                            {
                                getPedvarsDist()->Fill( fPedvar_high[p * VDST_MAXSUMWINDOW + s] );
                            }
                        }
                    }
                }
            }
        }
        else if( nPixel != getPedvars( false ).size() )
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile ) error: ";
            cout << "index out of range (pedvars, high gain): ";
            cout << nPixel << "\t" << getPedvars( false ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        if( nPixel == getPedvars( true ).size() && !iPedOnly )
        {
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( p < getPedvars( true ).size() )
                {
                    for( unsigned int s = 0; s < fnum_sumwindow; s++ )
                    {
                        if( iPedVarTreeTypeNew )
                        {
                            getPedvars( true, fsumwindow[s] )[p] = fPedvar_low[p * VDST_MAXSUMWINDOW + s];
                            if( fsumwindow[s] == 1 )
                            {
                                getPedrms( true ) = fPedvar_low[p * VDST_MAXSUMWINDOW + s];
                            }
                        }
                        if( fsumwindow[s] == getSumWindow() && getPedvarsDist( true ) )
                        {
                            if( iPedVarTreeTypeNew )
                            {
                                getPedvarsDist( true )->Fill( fPedvar_low[p * VDST_MAXSUMWINDOW + s] );
                            }
                        }
                    }
                }
            }
        }
        else if( nPixel != getPedvars( true ).size() )
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile ) error: ";
            cout << "index out of range (pedvars, low gain): ";
            cout << nPixel << "\t" << getPedvars( true ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        ////////////////
        // relative gains
        if( nPixel == getGains( false ).size() )
        {
            // mean conversion factor
            float i_z = 0;
            float i_meanC = 0;
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( fConv_high[p] > 0. )
                {
                    i_z++;
                    i_meanC += fConv_high[p];
                }
                // fill calibration data
                if( p < getFADCtoPhe( false ).size() )
                {
                    getFADCtoPhe( false )[p] = fConv_high[p];
                    getFADCtoPhe( false )[p] = 1.;   // TMPTMP: IPR graphs are still in dc, not in pe
                }
                if( p < getFADCtoPhe( true ).size() )
                {
                    getFADCtoPhe( true )[p] = fConv_low[p];
                    getFADCtoPhe( true )[p] = 1.;    // TMPTMP: IPR graphs are still in dc, not in pe
                }
            }
            if( i_z > 0. && i_meanC > 0. )
            {
                i_meanC /= i_z;
            }
            else
            {
                i_meanC  = 1.;
            }
            // calculate relative gain
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( fConv_high[p] > 0. && !getRunParameter()->fIgnoreDSTGains )
                {
                    getGains( false )[p] = 1. / fConv_high[p];
                    if( getGainDist( false ) )
                    {
                        getGainDist( false )->Fill( i_meanC / fConv_high[p] );
                    }
                }
                else
                {
                    getGains( false )[p] = 1.;
                }
            }
        }
        else
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile ) error: ";
            cout << "index out of range (gains, high gain): ";
            cout << nPixel << "\t" << getGains( false ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        if( nPixel == getGains( true ).size() )
        {
            // use high-gain conversion and understand differences between low and high gain as multiplication factor
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( getGains( false )[p] > 0. && !getRunParameter()->fIgnoreDSTGains )
                {
                    getGains( true )[p] = getGains( false )[p];
                    if( getGainDist( true ) )
                    {
                        getGainDist( true )->Fill( getGains( true )[p] );
                    }
                }
                else
                {
                    getGains( true )[p] = 1.;
                }
                // assumes fConv_low / fConv_high is the same for all pixels.
                // Ok as long as we don't have low gain in cta.
                if( fConv_low[p] > 0. && fConv_high[p] > 0. )
                {
                    setLowGainMultiplier_Trace( fConv_low[p] / fConv_high[p] );
                }
                else
                {
                    setLowGainMultiplier_Trace( 1. );
                }
            }
        }
        else
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile ) error: ";
            cout << "index out of range (gains, low gain): ";
            cout << nPixel << "\t" << getGains( true ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
        ////////////////
        // tzeros
        if( nPixel == getAverageTZeros( false ).size() )
        {
            double i_mean = 0.;
            double i_nn = 0.;
            // get teltype counter
            unsigned int iTelTypeC = 0;
            for( unsigned int t = 0; t < iTelType.size(); t++ )
            {
                if( iTelType[t] == getDetectorGeometry()->getTelType()[i] )
                {
                    iTelTypeC = t;
                    break;
                }
            }
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                getAverageTZeros( false )[p] = ftzero[p];
                // check if this pixel is inside the averateTZeroFiducialRadius
                if( getRunParameter()->faverageTZeroFiducialRadius < 1. )
                {
                    double d = getDetectorGeometry()->getOuterEdgeDistance( p );
                    if( d > 0.5 * getDetectorGeometry()->getFieldofView()[i] * getRunParameter()->faverageTZeroFiducialRadius )
                    {
                        continue;
                    }
                }
                if( getAverageTZeros( false )[p] > -98. )
                {
                    i_mean += ftzero[p];
                    i_nn++;
                    setAverageTZero( p, ftzero[p], false );
                    if( iTelTypeC < iH_averageTZero.size() && iH_averageTZero[iTelTypeC] )
                    {
                        iH_averageTZero[iTelTypeC]->Fill( ftzero[p] );
                    }
                }
            }
        }
        else
        {
            cout << "bool VCalibrator::readCalibrationDatafromDSTFiles( string iDSTfile ) error: ";
            cout << " index out of range (tzero, high gain): ";
            cout << nPixel << "\t" << getAverageTZeros( false ).size();
            cout << " (telescope " << getTelID() + 1 << ")" << endl;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // get average tzero per telescope type
    // use median, as outliers are expected
    cout << "\t calculate average tzeros" << endl;
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        setTelID( i );
        
        // no calibration data available for this telescope
        if( nPixel == 0 )
        {
            continue;
        }
        // get teltype counter
        unsigned int iTelTypeC = 0;
        for( unsigned int t = 0; t < iTelType.size(); t++ )
        {
            if( iTelType[t] == getDetectorGeometry()->getTelType()[i] )
            {
                iTelTypeC = t;
                break;
            }
        }
        if( iTelTypeC < iH_averageTZero.size() && iH_averageTZero[iTelTypeC] )
        {
            if( iH_averageTZero[iTelTypeC]->GetEntries() > 0. )
            {
                double i_a[] = { 0.5 };
                double i_b[] = { 0.0 };
                iH_averageTZero[iTelTypeC]->GetQuantiles( 1, i_b, i_a );
                setMeanAverageTZero( i_b[0], false );
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // read IPR graph from dst root file (for direct usage or creation of database )
    if( getRunParameter()->ifReadIPRfromDSTFile == true )
    {
        cout << "\t reading IPR graphs for NN image cleaning from DST file" << endl;
        for( int i = 0; i < t->GetEntries(); i++ )
        {
            setTelID( i );
            readIPRGraph_from_DSTFile( iDSTfile, getSumWindow(), getTelType( i ) );
        }
    }
    else if( getRunParameter()->ifReadIPRfromDatabase == true || getRunParameter()->ifCreateIPRdatabase == true )
    {
        cout << "\t reading IPR graphs for NN image cleaning" << endl;
        for( int i = 0; i < t->GetEntries(); i++ )
        {
            setTelID( i );
            calculateIPRGraphs( iDSTfile, getSumWindow(), getTelType( i ), i );
        }
    }
    else if( getRunParameter()->fIPRdatabase.Length() > 0 )
    {
        cout << "\t IPR graphs for NN cleaening will be read from external database: ";
        cout << getRunParameter()->fIPRdatabase << endl;
    }
    
    // cleanup (not necessary, this is done by TFile::Close()
    //	delete [] fPedvar_high;
    //	delete [] fPedvar_low;
    
    iF.Close();
    return true;
}

/*
 *
 * write IPR graphs to disk (one per telescope type)
 *
 */
bool VCalibrator::writeIPRgraphs( string iFile )
{
    TFile* fgraphs = 0;
    if( iFile.size() == 0 )
    {
        fgraphs = new TFile( getRunParameter()->fIPRdatabaseFile, "RECREATE" );
    }
    else
    {
        fgraphs = new TFile( iFile.c_str(), "UPDATE" );
    }
    if( fgraphs->IsZombie() )
    {
        cout << "VCalibrator::writeIPRgraphs error opening IPR output file: " << endl;
        cout << "\t" << fgraphs->GetName() << endl;
        return false;
    }
    
    // tree with conditions for these IRPs
    TTree* header = new TTree( "IPRheader", "IPR parameters" );
    unsigned int sumwin = 0;            // [FADC slices]
    unsigned int Nsamples = 0;          // [FADC slices]
    float ROwin = 0.;                   // [ns]
    float FADCslice = 0.;               // [ns]
    unsigned int SignalExtractor = 0;   // [ Selected Extractor ]
    unsigned int UpSample = 0;          // [Upsample in digital filter]
    header->Branch( "SignalExtractor", &SignalExtractor, "SignalExtractor/i" );
    header->Branch( "SummationWindow", &sumwin, "SummationWindow/i" );
    header->Branch( "Nsamples", &Nsamples, "Nsamples/i" );
    header->Branch( "FADCtimeSlice", &FADCslice, "FADCtimeSlice/F" );
    header->Branch( "ReadoutWindow", &ROwin, "ReadoutWindow/F" );
    header->Branch( "UpSample", &UpSample, "UpSample/i" );
    
    // one graph per telescope ID
    map< ULong64_t, bool > iTelDone;
    for( unsigned int i = 0; i < getDetectorGeometry()->getTelType_list().size(); i++ )
    {
        iTelDone[getDetectorGeometry()->getTelType_list()[i]] = false;
    }
    
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        setTelID( i );
        if( iTelDone[getTelType( i )] )
        {
            continue;
        }
        
        if( hped_vec.find( getTelType( i ) ) == hped_vec.end() )
        {
            continue;
        }
        
        // loop over all summation windows
        for( unsigned int j = 0; j < hped_vec[getTelType( i )].size(); j++ )
        {
            // summation window
            int i_sw = j + 1;
            
            TGraphErrors* g = getIPRGraph( i_sw, false );
            if( !g )
            {
                continue;
            }
            SignalExtractor = getRunParameter()->fTraceIntegrationMethod.at( getTelID() );
            sumwin          = i_sw;
            Nsamples        = getNSamples();
            FADCslice       = getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() ); //[ns]
            ROwin           = getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() ) * ( float )getNSamples(); // [ns]
            if( getDigitalFilterMethod() > 0 )
            {
                UpSample = getDigitalFilterUpSample();
            }
            
            header->Fill();
            g->Write();
        }
        
        iTelDone[getTelType( i )] = true;
    }
    header->Write();
    fgraphs->Close();
    fgraphs->Delete();
    
    return true;
}

unsigned int VCalibrator::getNumberOfEventsUsedInCalibration( int iTelID, int iType )
{
    if( iType == 1 || iType == 6 )
    {
        return getNumberOfEventsUsedInCalibration( fNumberPedestalEvents, iTelID );
    }
    else if( iType == 2 || iType == 5 )
    {
        return getNumberOfEventsUsedInCalibration( fNumberGainEvents, iTelID );
    }
    else if( iType == 7 || iType == 8 )
    {
        return getNumberOfEventsUsedInCalibration( fNumberTZeroEvents, iTelID );
    }
    return 0;
}

unsigned int VCalibrator::getNumberOfEventsUsedInCalibration( vector< int > iE, int iTelID )
{
    if( iTelID > 0 && iTelID < ( int )iE.size() )
    {
        return iE[iTelID];
    }
    else
    {
        int iTot = 0;
        int iN = 0;
        for( unsigned int tel = 0; tel < getTeltoAna().size(); tel++ )
        {
            if( tel < iE.size() )
            {
                iTot += iE[tel];
                iN++;
            }
        }
        if( iN > 0 )
        {
            return iTot / iN;
        }
    }
    
    return 0;
}

unsigned int VCalibrator::getNumberOfEventsUsedInCalibration( map< ULong64_t, int > iE, int iTelID )
{
    ULong64_t iTelType = 0;
    if( iTelID < ( int )getDetectorGeometry()->getTelType().size() )
    {
        iTelType  = getDetectorGeometry()->getTelType()[iTelID];
    }
    else
    {
        return 0;
    }
    
    if( iTelID > 0 && iE.find( iTelType ) != iE.end() )
    {
        return iE[iTelID];
    }
    else
    {
        int iTot = 0;
        int iN = 0;
        map< ULong64_t, int >::iterator iN_iter;
        for( iN_iter = iE.begin(); iN_iter != iE.end(); iN_iter++ )
        {
            iTot += iN_iter->second;
            iN++;
        }
        if( iN > 0 )
        {
            return iTot / iN;
        }
    }
    
    return 0;
}

/*
 * calculate IPR graphs and write them to disk
 *
 * (loop over all telescopes)
 *
*/
bool VCalibrator::calculateIPRGraphs()
{
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        setTelID( getTeltoAna()[i] );
        
        // first find dead channels
        // (ignore low gains here)
        findDeadChans( false );
        // calculate IPR graphs
        calculateIPRGraphs( fPedFileNameC[getTeltoAna()[i]], getSumWindow(), getTelType( getTeltoAna()[i] ), i );
    }
    return true;
}

/*
 * calculate IPR graphs and write them to disk
 *
 * (this is done per telescope type)
 *
 */
bool VCalibrator::calculateIPRGraphs( string iPedFileName, unsigned int iSummationWindow, ULong64_t iTelType, unsigned int i_tel )
{
    TDirectory* iG_CurrentDirectory = gDirectory;
    
    // get an IPR graph
    TGraphErrors* g = getIPRGraph( iSummationWindow, true );
    if( !g )
    {
        cout << "VCalibrator::calculateIPRGraphs info: no IPR graph found for telescope type " << iTelType << endl;
        return false;
    }
    // check suffix of ped file
    if( iPedFileName.find( ".root" ) == string::npos )
    {
        iPedFileName += ".root";
    }
    // open pedestal files
    TFile iPedFile( iPedFileName.c_str() );
    if( iPedFile.IsZombie() )
    {
        cout << "VCalibrator::calculateIPRGraphs error reading IPR graphs from ";
        cout << iPedFileName << endl;
        return false;
    }
    // histograms with IPR distributions are either
    // i) in a directory called distributions_TelType
    // ii) in the current directory
    cout << "Telescope type " << iTelType << ": ";
    cout << "reading IPR histograms for summation window " << iSummationWindow;
    cout << " from ";
    cout << iPedFileName << endl;
    
    stringstream i_Directory( stringstream::in | stringstream::out );
    i_Directory << "distributions_" << iTelType;
    if( iPedFile.Get( i_Directory.str().c_str() ) )
    {
        iPedFile.cd( i_Directory.str().c_str() );
    }
    
    // get charge distribution for first channel as reference histogram
    stringstream i_hIPRname( stringstream::in | stringstream::out );
    i_hIPRname << "hpedPerTelescopeType_" << iTelType << "_" << iSummationWindow << "_" << 1;
    TH1F* href = ( TH1F* )gDirectory->Get( i_hIPRname.str().c_str() );
    
    if( !href )
    {
        cout << " Error: cannot find IPR histogram " << i_hIPRname.str().c_str();
        cout << " in file:" << iPedFileName << " ... exiting " << endl;
        return false;
    }
    
    ///////////////////////////
    // summary histogram
    TH1F* hIPR = 0;
    if( getRunParameter()->fIgnoreDSTGains )
    {
        // work in dc
        hIPR = new TH1F( "", "", int( 1.5 * href->GetNbinsX() + 0.5 ), href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax() );
    }
    else
    {
        // work in pe
        hIPR = new TH1F( "", "", 1000, 0., 100. );
    }
    hIPR->Reset();
    
    //////////////////////////////////////////////////
    // average over all channels in one telescope
    float i_gainCorrect = 1.;
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        if( getDetectorGeometry()->getAnaPixel()[i] > 0
                && i < getDead().size() && !getDead()[i] )
        {
            stringstream i_Hname( stringstream::in | stringstream::out );
            i_Hname << "hpedPerTelescopeType_" << iTelType << "_" << iSummationWindow << "_" << i;
            TH1F* h = ( TH1F* )gDirectory->Get( i_Hname.str().c_str() );
            
            if( h )
            {
                float ped = 0;
                // default: pedestals are subtracted here
                // (for combined channel analysis: charges are filled already pedestal subtracted)
                // apply relative gain correction to integrated charges
                if( getRunParameter()->fCombineChannelsForPedestalCalculation == 0 )
                {
                    ped = getPeds()[i];
                }
                // special treatment for ASTRI telescopes
                else if( getRunParameter()->fCombineChannelsForPedestalCalculation == 2 )
                {
                    stringstream i_Pname( stringstream::in | stringstream::out );
                    i_Pname << "hped_" << iTelType << "_" << iSummationWindow << "_" << i;
                    TH1F* hP = ( TH1F* )gDirectory->Get( i_Pname.str().c_str() );
                    if( hP )
                    {
                        ped = hP->GetMean();
                    }
                }
                for( int j = 1; j <= h->GetNbinsX(); j++ )
                {
                    if( h->GetBinContent( j ) > 0. && getGains()[i] > 0 )
                    {
                        i_gainCorrect = getGains()[i];
                        if( getHIGHQE_gainfactor( i ) > 0. )
                        {
                            i_gainCorrect *= getHIGHQE_gainfactor( i );
                        }
                        if( i_gainCorrect > 0. )
                        {
                            if( hasFADCData() )
                            {
                                hIPR->Fill( ( h->GetBinCenter( j ) - ped * iSummationWindow ) / i_gainCorrect, h->GetBinContent( j ) );
                            }
                            else
                            {
                                hIPR->Fill( ( h->GetBinCenter( j ) - ped ) / i_gainCorrect, h->GetBinContent( j ) );
                            }
                        }
                    }
                }
            }
        }
    }
    
    int z = 0;
    float norm = hIPR->Integral( 1, hIPR->GetNbinsX() );
    if( norm < fPedPerTelescopeTypeMinCnt )  //statistical limit for number of counts
    {
        cout << "Telescope " << iTelType << ": ";
        cout << "VCalibrator::calculateIPRGraphs WARNING: too few statistics to measure IPR curves ";
        cout << "(total counts available: " << norm << ", ";
        cout << "current limit " << fPedPerTelescopeTypeMinCnt << ")" << endl;
    }
    if( norm == 0 )
    {
        cout << "VCalibrator::calculateIPRGraphs ERROR: no counts in IPR histogram !" << endl;
        return false;
    }
    // convert to Rate
    float nsToSec = 1E-9;
    float Tsearch = getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() );
    if( getSearchWindowLast() < getNSamples() )
    {
        Tsearch *= ( getSearchWindowLast() - getSumFirst() ); // [ns]
    }
    else
    {
        Tsearch *= ( getNSamples() - getSumFirst() ); // [ns]
    }
    // digital filter: take upsampling into account
    if( getDigitalFilterMethod() > 0 && getDigitalFilterUpSample() > 0 )
    {
        if( getSearchWindowLast() < getNSamples() * getDigitalFilterUpSample() )
        {
            Tsearch  = getSearchWindowLast() - getSumFirst();
        }
        else
        {
            Tsearch  = getNSamples() * getDigitalFilterUpSample() - getSumFirst();
        }
        Tsearch *= getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() ) / getDigitalFilterUpSample();
    }
    float convToHz = 1.;
    if( nsToSec > 0. && Tsearch > 0. )
    {
        convToHz /= ( nsToSec * Tsearch );
    }
    else if( getRunParameter()->fImageCleaningParameters[i_tel]->fNNOpt_ifExplicitSampleTimeSlice
             && getRunParameter()->fImageCleaningParameters[i_tel]->fNNOpt_sampleTimeSlice > 0
             && getRunParameter()->fImageCleaningParameters[i_tel]->fNNOpt_nBinsADC > 0 )
    {
        // simple peak sensing: sim_telarray uses the maximum bin only
        // The values for sampleTimeSlice and nBinsADC are set in the cleaning parameter file
        // For example, for the currect (Apr 17) ASTRI simulation, it is set in sim_telarray as
        // fadc_mhz = 500 % MHz ==> sampleTimeSlice = 2 ns
        // fadc_sum_bins = nBinsADC = 25 % Number of ADC time intervals actually summed up.
        convToHz /= ( nsToSec
                      * getRunParameter()->fImageCleaningParameters[i_tel]->fNNOpt_sampleTimeSlice
                      * getRunParameter()->fImageCleaningParameters[i_tel]->fNNOpt_nBinsADC );
    }
    
    for( int i = 1; i <= hIPR->GetNbinsX(); i++ )
    {
        if( hIPR->GetBinContent( i ) > 5 )
        {
            double val = convToHz * hIPR->Integral( i, hIPR->GetNbinsX() ) / norm;
            double valerr = convToHz * sqrt( hIPR->Integral( i, hIPR->GetNbinsX() ) ) / norm;
            double charge_pe = hIPR->GetXaxis()->GetBinCenter( i ) * getTelescopeAverageFADCtoPhe();
            double charge_pe_bin_width = 0.5 * hIPR->GetXaxis()->GetBinWidth( i ) * getTelescopeAverageFADCtoPhe();
            
            g->SetPoint( z, charge_pe, val );
            g->SetPointError( z, charge_pe_bin_width, valerr );
            z++;
        }
    }
    
    g->SetMinimum( 1 );
    if( getDigitalFilterMethod() > 0 && getDigitalFilterUpSample() )
    {
        g->SetTitle( Form( "Rate vs Threshold. W_{RO}=%2.1f ns, W_{int}=%2.1f ns", Tsearch,
                           getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() )
                           / getDigitalFilterUpSample() * iSummationWindow ) );
    }
    else
    {
        g->SetTitle( Form( "Rate vs Threshold. W_{RO}=%2.1f ns, W_{int}=%2.1f ns", Tsearch,
                           getDetectorGeometry()->getLengthOfSampleTimeSlice( getTelID() ) * iSummationWindow ) );
    }
    if( getRunParameter()->fIgnoreDSTGains )
    {
        g->GetXaxis()->SetTitle( "Threshold [FADC counts]" );
    }
    else
    {
        g->GetXaxis()->SetTitle( "Threshold [p.e.]" );
    }
    g->GetYaxis()->SetTitle( "Rate above Threshold [Hz]" );
    g->SetName( Form( "IPRcharge_TelType%d_SW%d", ( int )iTelType, iSummationWindow ) );
    hIPR->Delete();
    
    iPedFile.Close();
    
    iG_CurrentDirectory->cd();
    
    return true;
}

bool VCalibrator::readIPRGraph_from_DSTFile( string iDSTFile, unsigned int iSummationWindow, ULong64_t iTelType )
{
    TDirectory* iG_CurrentDirectory = gDirectory;
    
    // IPR graph
    TGraphErrors* gIPR = 0;
    
    TFile* iFileIn = new TFile( iDSTFile.c_str() );
    if( iFileIn->IsZombie() )
    {
        cout << "VCalibrator::readIPRGraph_from_DSTFile: error opening DST file: " << endl;
        cout << "\t" << iDSTFile << endl;
        return false;
    }
    std::ostringstream iSname;
    iSname << "IPRcharge_TelType";
    iSname << iTelType;
    iSname << "_SW" << iSummationWindow;
    
    cout << "\t IPR graph: " << iSname.str() << endl;
    
    gIPR = ( TGraphErrors* )iFileIn->Get( iSname.str().c_str() );
    if( gIPR )
    {
        if( iTelType == 201511619 )
        {
            setIPRGraph( 0, ( TGraphErrors* )gIPR->Clone() );
        }
        setIPRGraph( iSummationWindow, ( TGraphErrors* )gIPR->Clone() );
    }
    else
    {
        cout << "VCalibrator::readIPRGraph_from_DSTFile: warning, could not read IPR graph ";
        cout << iSname.str();
        cout << " from DST file: " << endl;
        cout << "\t" << iDSTFile << endl;
    }
    
    iFileIn->Close();
    
    iG_CurrentDirectory->cd();
    
    return true;
}

