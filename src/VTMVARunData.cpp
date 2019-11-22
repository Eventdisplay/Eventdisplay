/*! \class VTMVARunData
    \brief run parameter data class for TMVA optimization

*/

#include "VTMVARunData.h"

VTMVARunData::VTMVARunData()
{
    setDebug();
    
    fName = "noname";
    
    fTrainGammaHadronSeparation = true;
    fTrainReconstructionQuality = false;  // in development: please ignore
    
    fCheckValidityOfInputVariables = true;
    
    fOutputDirectoryName = "";
    fOutputFileName = "";
    
    fQualityCuts = "1";
    fQualityCutsBkg = "1";
    fQualityCutsSignal = "1";
    fMCxyoffCut = "1";
    fMCxyoffCutSignalOnly = false;
    fAzimuthCut = "1";
    fPrepareTrainingOptions = "SplitMode=random:!V";
    
    fSignalWeight = 1.;
    fBackgroundWeight = 1.;
    fMinSignalEvents = 50;
    fMinBackgroundEvents = 0;
    
    fNTtype = -1;
    
    fReconstructionQualityTarget = "ErecS/MCe0";
    fReconstructionQualityTargetName = "EQuality";
}

/*!

    open data files and make data trees available

*/
bool VTMVARunData::openDataFiles()
{
    if( fDebug )
    {
        cout << "VTMVARunData::openDataFiles()" << endl;
    }
    // open signal trees
    fSignalTree.clear();
    
    for( unsigned int i = 0; i < fSignalFileName.size(); i++ )
    {
        fSignalTree.push_back( new TChain( "data" ) );
        Int_t i_nfiles = fSignalTree.back()->Add( fSignalFileName[i].c_str() );
        if( i_nfiles == 0 )
        {
            cout << "VTMVARunData::openDataFiles() error: no data tree in signal file: " << fSignalFileName[i] << endl;
            cout << "aborting..." << endl;
            return false;
        }
        // get for first tree the number of telescope types
        // (note this is the source of the the Error in <TLeafI::GetLen>: Leaf counter is greater than maximum! messages
        if( fNTtype < 0 )
        {
            Int_t NTtype = 0;
            fSignalTree.back()->SetBranchAddress( "NTtype", &NTtype );
            if( fSignalTree.back()->GetEntries() > 0 )
            {
                fSignalTree.back()->GetEntry( 0 );
            }
            fNTtype = NTtype;
        }
    }
    // open background trees
    fBackgroundTree.clear();
    for( unsigned int i = 0; i < fBackgroundFileName.size(); i++ )
    {
        fBackgroundTree.push_back( new TChain( "data" ) );
        Int_t i_nfiles = fBackgroundTree.back()->Add( fBackgroundFileName[i].c_str() );
        if( i_nfiles == 0 )
        {
            cout << "VTMVARunData::openDataFiles() error: no data tree in background file: " << fBackgroundFileName[i] << endl;
            cout << "aborting..." << endl;
            return false;
        }
    }
    
    if( fDebug )
    {
        cout << "VTMVARunData::openDataFiles() open output files " << endl;
    }
    
    ///////////////////////////////////////////////////////////////////
    // check how many events there are in signal and background trees (after cuts)
    // (note that no cuts are applied here)
    if( fMinSignalEvents > 0 && fMinBackgroundEvents > 0 )
    {
        TEntryList* i_j_SignalList = 0;
        TEntryList* i_j_BackgroundList = 0;
        bool iEnoughEvents = true;
        
        
        // loop over all energy and zenith bins
        for( unsigned int i = 0; i < fEnergyCutData.size(); i++ )
        {
            for( unsigned int j = 0; j < fZenithCutData.size(); j++ )
            {
                if( fMinSignalEvents > 0 )
                {
                    for( unsigned int k = 0; k < fSignalTree.size(); k++ )
                    {
                        if( fSignalTree[k] )
                        {
                            fSignalTree[k]->Draw( ">>+signalList", fQualityCuts && fQualityCutsSignal && fMCxyoffCut && fAzimuthCut && fEnergyCutData[i]->fEnergyCut && fZenithCutData[j]->fZenithCut, "entrylist" );
                            i_j_SignalList = ( TEntryList* )gDirectory->Get( "signalList" );
                        }
                    }
                    
                    if( i_j_SignalList )
                    {
                        cout << "number of signal events in energy bin " << i << " zenith bin " << j << ": ";
                        cout << i_j_SignalList->GetN() << "\t required > " << fMinSignalEvents << endl;
                        cout << "  (cuts are " << fQualityCuts << "&&" << fQualityCutsSignal << "&&" << fMCxyoffCut << "&&" << fAzimuthCut;
                        cout << "&&" << fEnergyCutData[i]->fEnergyCut  << "&&" << fZenithCutData[j]->fZenithCut << ")" << endl;
                        if( i_j_SignalList->GetN() < fMinSignalEvents )
                        {
                            iEnoughEvents = false;
                        }
                        i_j_SignalList->Reset();
                    }
                }
                // background events
                if( fMinBackgroundEvents > 0 )
                {
                    for( unsigned int k = 0; k < fBackgroundTree.size(); k++ )
                    {
                        if( fMCxyoffCutSignalOnly )
                        {
                            fBackgroundTree[k]->Draw( ">>+BackgroundList", fQualityCuts && fQualityCutsBkg
                                                      && fMCxyoffCut && fAzimuthCut && fEnergyCutData[i]->fEnergyCut
                                                      && fZenithCutData[j]->fZenithCut, "entrylist" );
                            i_j_BackgroundList = ( TEntryList* )gDirectory->Get( "BackgroundList" );
                        }
                        else if( fBackgroundTree[k] )
                        {
                            fBackgroundTree[k]->Draw( ">>+BackgroundList", fQualityCuts && fQualityCutsBkg
                                                      && fMCxyoffCut && fAzimuthCut && fEnergyCutData[i]->fEnergyCut
                                                      && fZenithCutData[j]->fZenithCut, "entrylist" );
                            i_j_BackgroundList = ( TEntryList* )gDirectory->Get( "BackgroundList" );
                        }
                    }
                    if( i_j_BackgroundList )
                    {
                        cout << "number of background events in energy bin " << i <<  " zenith bin " << j << ": ";
                        cout << i_j_BackgroundList->GetN();
                        cout << "\t required > " << fMinBackgroundEvents << endl;
                        cout << "  (cuts are " << fQualityCuts << "&&" << fQualityCutsBkg << "&&" << fMCxyoffCut << "&&" << fAzimuthCut;
                        cout << "&&" << fEnergyCutData[i]->fEnergyCut  << "&&" << fZenithCutData[j]->fZenithCut << ")" << endl;
                        if( i_j_BackgroundList->GetN() < fMinBackgroundEvents )
                        {
                            iEnoughEvents = false;
                        }
                        i_j_BackgroundList->Reset();
                    }
                }
            }
        }
        if( !iEnoughEvents )
        {
            cout << endl;
            cout << "ERROR: not enough signal or/and background events" << endl;
            cout << "please adjust energy intervals " << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    ///////////////////////////////////////////////////////////////////
    // open output file
    cout << "output file name size " << fOutputFileName.size() << endl;
    if( fOutputFileName.size() > 0 && fOutputDirectoryName.size() > 0 )
    {
        for( unsigned int i = 0; i < fEnergyCutData.size(); i++ )
        {
            vector< TFile* > output_zenith;
            for( unsigned int j = 0; j < fZenithCutData.size(); j++ )
            {
                stringstream iTempS;
                stringstream iTempS2;
                gSystem->mkdir( fOutputDirectoryName.c_str() );
                if( fEnergyCutData.size() > 1 && fZenithCutData.size() > 1 )
                {
                    iTempS << fOutputDirectoryName << "/" << fOutputFileName << "_" << i << "_" << j << ".root";    // append a _# at the file name
                    iTempS2 << fOutputFileName << "_" << i << "_" << j;
                }
                else if( fEnergyCutData.size() > 1 && fZenithCutData.size() <= 1 )
                {
                    iTempS << fOutputDirectoryName << "/" << fOutputFileName << "_" << i << ".root";    // append a _# at the file name
                    iTempS2 << fOutputFileName << "_" << i;
                }
                else if( fZenithCutData.size() > 1 &&  fEnergyCutData.size() <= 1 )
                {
                    iTempS << fOutputDirectoryName << "/" << fOutputFileName << "_0_" << j << ".root";    // append a _# at the file name
                    iTempS2 << fOutputFileName << "_0_" << i;
                }
                else
                {
                    iTempS << fOutputDirectoryName << "/" << fOutputFileName << ".root";
                    iTempS2 << fOutputFileName;
                }
                output_zenith.push_back( new TFile( iTempS.str().c_str(), "RECREATE" ) );
                if( output_zenith.back()->IsZombie() )
                {
                    cout << "VTMVARunData::openDataFiles() error creating output file " << output_zenith.back()->GetName() << endl;
                    cout << "aborting..." << endl;
                    return false;
                }
                output_zenith.back()->SetTitle( iTempS2.str().c_str() );
                if( fEnergyCutData[i] )
                {
                    fEnergyCutData[i]->Write();
                }
                if( fZenithCutData[j] )
                {
                    fZenithCutData[j]->Write();
                }
            }
            fOutputFile.push_back( output_zenith );
        }
    }
    cout << "output file size " << fOutputFile.size()*fOutputFile[0].size() << endl;
    
    if( fDebug )
    {
        cout << "VTMVARunData::openDataFiles() END" << endl;
    }
    
    return true;
}

/*!
    print run information to screen
*/
void VTMVARunData::print()
{
    cout << endl;
    if( fTrainGammaHadronSeparation )
    {
        cout << "Training gamma/hadron separation" << endl;
    }
    if( fTrainReconstructionQuality )
    {
        cout << "Training reconstruction quality" << endl;
    }
    cout << "MVA Methods and options: " << endl;
    for( unsigned int i = 0; i < fMVAMethod.size(); i++ )
    {
        cout << "METHOD: " << fMVAMethod[i];
        if( i < fMVAMethod_Options.size() )
        {
            cout << "  Options: " << fMVAMethod_Options[i];
        }
        cout << endl;
    }
    cout << endl << "list of variables: " << endl;
    for( unsigned int i = 0; i < fTrainingVariable.size(); i++ )
    {
        cout << "\t" << fTrainingVariable[i];
        if( i < fTrainingVariableType.size() )
        {
            cout << "\t\t" << fTrainingVariableType[i];
        }
        if( i < fTrainingVariable_CutRangeMin.size() )
        {
            cout << "\t" << fTrainingVariable_CutRangeMin[i];
        }
        if( i < fTrainingVariable_CutRangeMax.size() )
        {
            cout << "\t" << fTrainingVariable_CutRangeMax[i];
        }
        if( i < fTrainingVariable_VarProp.size() )
        {
            cout << "\t" << fTrainingVariable_VarProp[i];
        }
        cout << endl;
    }
    if( fSpectatorVariable.size() > 0 )
    {
        cout << endl;
        cout << "list of spectator variables: " << endl;
        for( unsigned int i = 0; i < fSpectatorVariable.size(); i++ )
        {
            cout << "\t spectator: " << fSpectatorVariable[i] << endl;
        }
    }
    cout << endl;
    cout << "pre-training selection cuts: " << fQualityCuts << endl;
    cout << "background specific pre-training selection cuts: " << fQualityCutsBkg << endl;
    cout << "signal specific pre-training selection cuts: " << fQualityCutsSignal << endl;
    cout << "cut on MC arrival directions: " << fMCxyoffCut;
    if( fMCxyoffCutSignalOnly )
    {
        cout << " (signal only)";
    }
    cout << endl;
    cout << "azimuth cut: " << fAzimuthCut << endl;
    cout << endl;
    cout << "prepare training options: " << fPrepareTrainingOptions << endl;
    cout << "energy bin(s) [log10(TeV)] (" << fEnergyCutData.size() << "): ";
    for( unsigned int i = 0; i < fEnergyCutData.size(); i++ )
    {
        cout << "[" << fEnergyCutData[i]->fEnergyCut_Log10TeV_min << ", " << fEnergyCutData[i]->fEnergyCut_Log10TeV_max << "]";
    }
    cout << endl;
    cout << "zenith bin(s) [deg] (" << fZenithCutData.size() << "): ";
    for( unsigned int j = 0; j < fZenithCutData.size(); j++ )
    {
        cout << "[" << fZenithCutData[j]->fZenithCut_min << ", " << fZenithCutData[j]->fZenithCut_max << "]";
    }
    cout << endl;
    // all bins should use same energy reconstruction method
    if( fEnergyCutData.size() > 0 && fEnergyCutData[0] )
    {
        cout << "energy reconstruction method " << fEnergyCutData[0]->fEnergyReconstructionMethod << endl;
    }
    cout << "signal data file(s): " << endl;
    for( unsigned int i = 0; i < fSignalFileName.size(); i++ )
    {
        cout << "\t" << fSignalFileName[i] << endl;
    }
    if( !fTrainReconstructionQuality )
    {
        cout << "background data file(s): " << endl;
        for( unsigned int i = 0; i < fBackgroundFileName.size(); i++ )
        {
            cout << "\t" << fBackgroundFileName[i] << endl;
        }
    }
    cout << "output file: " << fOutputFileName << " (" << fOutputDirectoryName << ")" << endl;
    cout << endl;
    cout << endl;
}

/*!

    read configuration file (ascii format)

*/
bool VTMVARunData::readConfigurationFile( char* iC )
{
    cout << "reading TMVA optimizer configuration from " << iC << endl;
    
    ifstream is;
    is.open( iC, ifstream::in );
    if( !is )
    {
        cout << "VTMVARunData::readConfigurationFile error configuration file not found: " << iC << endl;
        return false;
    }
    
    string is_line;
    string temp;
    fMVAMethod.clear();
    fMVAMethod_Options.clear();
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            if( is_stream.eof() )
            {
                continue;
            }
            
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            if( is_stream.eof() )
            {
                continue;
            }
            
            is_stream >> temp;
            ///////////////////////////////////////////////////////////////////////////////////////////
            // MVA method and options
            ///////////////////////////////////////////////////////////////////////////////////////////
            if( temp == "MVA_METHOD" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> temp;
                    fMVAMethod.push_back( temp );
                }
                if( !is_stream.eof() )
                {
                    is_stream >> temp;
                    fMVAMethod_Options.push_back( temp );
                }
                else
                {
                    fMVAMethod_Options.push_back( "" );
                }
            }
            // Box cuts: kept for backwards compatibility
            if( temp == "OPTIMIZATION_METHOD" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> temp;
                    fMVAMethod.push_back( "BOXCUTS" );
                    fMVAMethod_Options.push_back( temp );
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable OPTIMIZATION_METHOD" << endl;
                    return false;
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////
            // training variables
            if( temp == "VARIABLE" )
            {
                if( !is_stream.eof() )
                {
                    char iV = 'F';
                    if( !is_stream.eof() )
                    {
                        is_stream >> iV;
                    }
                    fTrainingVariableType.push_back( iV );
                    float iR = -1.;
                    if( !is_stream.eof() )
                    {
                        is_stream >> iR;
                    }
                    fTrainingVariable_CutRangeMin.push_back( iR );
                    iR = -1.;
                    if( !is_stream.eof() )
                    {
                        is_stream >> iR;
                    }
                    fTrainingVariable_CutRangeMax.push_back( iR );
                    temp = "NotEnforced";
                    if( !is_stream.eof() )
                    {
                        is_stream >> temp;
                    }
                    fTrainingVariable_VarProp.push_back( temp );
                    fTrainingVariable.push_back( is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ) );
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable VARIABLE" << endl;
                    return false;
                }
            }
            // spectator variables
            if( temp == "SPECTATOR" )
            {
                if( !is_stream.eof() )
                {
                    fSpectatorVariable.push_back( is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ) );
                }
            }
            // preselection cut
            if( temp == "SELECTION_CUTS" )
            {
                if( !is_stream.eof() )
                {
                    fQualityCuts = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable SELECTION_CUTS" << endl;
                    return false;
                }
            }
            // preselection cuts specific to background
            if( temp == "SELECTION_CUTS_BKG" )
            {
                if( !is_stream.eof() )
                {
                    fQualityCutsBkg = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable SELECTION_CUTS_BKG" << endl;
                    return false;
                }
            }
            // preselection cuts specific to signal
            if( temp == "SELECTION_CUTS_SIGNAL" )
            {
                if( !is_stream.eof() )
                {
                    fQualityCutsSignal = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable SELECTION_CUTS_SIG" << endl;
                    return false;
                }
            }
            // MC arrival direction cut
            if( temp == "MCXYOFF" )
            {
                if( !is_stream.eof() )
                {
                    fMCxyoffCut = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable MCXYOFF" << endl;
                    return false;
                }
            }
            if( temp == "MCXYCUTSignalOnly" )
            {
                if( !is_stream.eof() )
                {
                    fMCxyoffCutSignalOnly = ( atoi )( is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str() );
                }
            }
            // cut on azimuth direction
            if( temp == "AZIMUTH" )
            {
                if( !is_stream.eof() )
                {
                    fAzimuthCut = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable AZIMUTH" << endl;
                    return false;
                }
            }
            // prepare training options
            if( temp == "PREPARE_TRAINING_OPTIONS" )
            {
                if( !is_stream.eof() )
                {
                    fPrepareTrainingOptions = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                    fPrepareTrainingOptions = VUtilities::removeSpaces( fPrepareTrainingOptions );
                    // remove all spaces
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable PREPARE_TRAINING_OPTIONS" << endl;
                    return false;
                }
            }
            // check event validity
            if( temp == "CHECKEVENTVALIDITY" )
            {
                if( !is_stream.eof() )
                {
                    int iT = 0;
                    is_stream >> iT;
                    fCheckValidityOfInputVariables = ( bool )iT;
                }
            }
            // signal weight
            if( temp == "SIGNALWEIGHT" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> fSignalWeight;
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable SIGNALWEIGHT" << endl;
                    return false;
                }
            }
            // signal files
            if( temp == "SIGNALFILE" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> temp;
                    fSignalFileName.push_back( temp );
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable SIGNALFILE" << endl;
                    return false;
                }
            }
            // background weight
            if( temp == "BACKGROUNDWEIGHT" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> fBackgroundWeight;
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable BACKGROUNDWEIGHT" << endl;
                    return false;
                }
            }
            // background files
            if( temp == "BACKGROUNDFILE" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> temp;
                    fBackgroundFileName.push_back( temp );
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable BACKGROUNDFILE" << endl;
                    return false;
                }
            }
            // output file
            if( temp == "OUTPUTFILE" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> fOutputDirectoryName;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> fOutputFileName;
                }
                else
                {
                    cout << "VTMVARunData::readConfigurationFile error while reading input for variable OUTPUTFILE" << endl;
                    return false;
                }
            }
            // energy bins
            if( temp == "ENERGYBINS" )
            {
                vector< double > iEnergyCut_Log10TeV_min;
                vector< double > iEnergyCut_Log10TeV_max;
                vector< TCut > iEnergyCut;
                
                // energy reconstruction method (should be 1, unless you know it better)
                unsigned int iEMethod;
                if( !is_stream.eof() )
                {
                    is_stream >> iEMethod;
                }
                
                // read in energy bin
                while( !is_stream.eof() )
                {
                    double iT = 0.;
                    is_stream >> iT;
                    iEnergyCut_Log10TeV_min.push_back( iT );
                }
                // sort
                sort( iEnergyCut_Log10TeV_min.begin(), iEnergyCut_Log10TeV_min.end() );
                // check sanity
                if( iEnergyCut_Log10TeV_min.size() < 2 )
                {
                    cout << "VTMVARunData::readConfigurationFile error: need at least two energy bins " << iEnergyCut_Log10TeV_min.size() << endl;
                    return false;
                }
                // fill maximum bins
                for( unsigned int i = 1; i < iEnergyCut_Log10TeV_min.size(); i++ )
                {
                    iEnergyCut_Log10TeV_max.push_back( iEnergyCut_Log10TeV_min[i] );
                }
                // remove last minimum
                iEnergyCut_Log10TeV_min.pop_back();
                // fill cuts
                for( unsigned int i = 0; i < iEnergyCut_Log10TeV_min.size(); i++ )
                {
                    ostringstream iCut;
                    if( iEMethod == 0 )
                    {
                        iCut << "Erec>0.&&"  << iEnergyCut_Log10TeV_min[i]  <<  "<log10(Erec)&&log10(Erec)<" << iEnergyCut_Log10TeV_max[i];
                    }
                    else
                    {
                        iCut << "ErecS>0.&&" <<  iEnergyCut_Log10TeV_min[i] <<  "<log10(ErecS)&&log10(ErecS)<" << iEnergyCut_Log10TeV_max[i];
                    }
                    iEnergyCut.push_back( iCut.str().c_str() );
                }
                // filling everything into the energy data structure
                fEnergyCutData.clear();
                for( unsigned int i = 0; i < iEnergyCut_Log10TeV_min.size(); i++ )
                {
                    fEnergyCutData.push_back( new VTMVARunDataEnergyCut() );
                    fEnergyCutData.back()->SetName( "fDataEnergyCut" );
                    fEnergyCutData.back()->fEnergyCutBin = 0;
                    fEnergyCutData.back()->fEnergyCut_Log10TeV_min = iEnergyCut_Log10TeV_min[i];
                    fEnergyCutData.back()->fEnergyCut_Log10TeV_max = iEnergyCut_Log10TeV_max[i];
                    fEnergyCutData.back()->fEnergyCut = iEnergyCut[i];
                    fEnergyCutData.back()->fEnergyReconstructionMethod = iEMethod;
                }
            }
            // zenith angle bins (in [deg])
            if( temp == "ZENBINS" || temp == "ZENITHBINS" )
            {
                vector< double > iZenithCut_min;
                vector< double > iZenithCut_max;
                vector< TCut > iZenithCut;
                
                // read in zenith angle bin
                while( !is_stream.eof() )
                {
                    double iT = 0.;
                    is_stream >> iT;
                    iZenithCut_min.push_back( iT );
                }
                // sort
                sort( iZenithCut_min.begin(), iZenithCut_min.end() );
                // check sanity
                if( iZenithCut_min.size() < 2 )
                {
                    cout << "VTMVARunData::readConfigurationFile error: need at least one zenith bin " << iZenithCut_min.size() << endl;
                    return false;
                }
                // fill maximum bins
                for( unsigned int i = 1; i < iZenithCut_min.size(); i++ )
                {
                    iZenithCut_max.push_back( iZenithCut_min[i] );
                }
                // remove last minimum
                iZenithCut_min.pop_back();
                // fill cuts
                for( unsigned int i = 0; i < iZenithCut_min.size(); i++ )
                {
                    ostringstream iCut;
                    iCut << "Ze>0.&&"  << iZenithCut_min[i]  <<  "<Ze&&Ze<" << iZenithCut_max[i];
                    
                    iZenithCut.push_back( iCut.str().c_str() );
                }
                // filling everything into the zenith data structure
                fZenithCutData.clear();
                for( unsigned int i = 0; i < iZenithCut_min.size(); i++ )
                {
                    fZenithCutData.push_back( new VTMVARunDataZenithCut() );
                    fZenithCutData.back()->SetName( "fDataZenithCut" );
                    fZenithCutData.back()->fZenithCutBin = 0;
                    fZenithCutData.back()->fZenithCut_min = iZenithCut_min[i];
                    fZenithCutData.back()->fZenithCut_max = iZenithCut_max[i];
                    fZenithCutData.back()->fZenithCut = iZenithCut[i];
                }
            }
            // minimum number of events
            if( temp == "MINEVENTS" )
            {
                if( !is_stream.eof() )
                {
                    is_stream >> fMinSignalEvents;
                }
                if( !is_stream.eof() )
                {
                    is_stream >> fMinBackgroundEvents;
                }
            }
        }
    }
    
    return true;
}

