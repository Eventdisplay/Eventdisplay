/*! \class VTMVAEvaluator
    \brief use a TMVA weight file for energy dependent gamma/hadron separation

*/

#include "VTMVAEvaluator.h"

VTMVAEvaluator::VTMVAEvaluator()
{
    fIsZombie = false;
    
    setDebug();
    
    fData = 0;
    fTMVAEvaluatorResults = 0;
    
    reset();
}

VTMVAEvaluator::~VTMVAEvaluator()
{
    if( fTMVACutValueFile )
    {
        delete fTMVACutValueFile;
    }
}

void VTMVAEvaluator::reset()
{
    ///////////////////////////////////////////////
    // variables used for gamma/hadron separation
    fNImages = 0.;
    fMSCW = 0.;
    fMSCL = 0.;
    fMWR = 0.;
    fMLR = 0.;
    fEmissionHeight = 0.;
    fEmissionHeightChi2_log10 = 0.;
    fEnergyReconstructionMethod = 1;
    fEChi2S = 0.;
    fEChi2S_log10 = 0.;
    fEChi2S_gt0 = 0.;
    fEChi2S_gt0_bool = 0.;
    fdES = 0.;
    fSizeSecondMax_log10 = 0;
    fTheta2 = 0.;
    fCoreDist = 0.;
    // disp below
    fDispDiff_log10 = 0.;
    fDispDiff_gt0 = 0.;
    fDispDiff_gt0_bool = 0.;
    fDummy = 0.;
    for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
    {
        fImages_Ttype[i] = 0.;
    }
    
    ///////////////////////////////////////////////
    // set default optimization values
    setTMVACutValue();
    setSignalEfficiency();
    setIgnoreTheta2Cut();
    setSpectralIndexForEnergyWeighting();
    setParticleNumberFile();
    setPlotEfficiencyPlotsPerBin();
    setPrintPlotting();
    setSensitivityOptimizationParameters();
    setSensitivityOptimizationFixedSignalEfficiency();
    setSensitivityOptimizationSourceStrength();
    setOptimizeAngularContainment();
    setTMVAAngularContainmentThetaFixedMinRadius();
    setTMVAMethod();
    // default: don't expect that the theta2 cut is performed here
    setTMVAThetaCutVariable( false );
    setTMVAErrorFraction();
    setTMVAAngularContainmentRadiusMax();
    fTMVA_EvaluationResult = -99.;
    fTMVACutValueNoVec = -99.;
    
    setSmoothAndInterpolateMVAValues();
    
    fWeightFileIndex_Emin = 0;
    fWeightFileIndex_Emax = 0;
    fWeightFileIndex_Zmin = 0;
    fWeightFileIndex_Zmax = 0;
    
}

/*

    get list of training variables from TMVA XML file

*/
vector< string > VTMVAEvaluator::getTrainingVariables( string iXMLFile, vector< bool >& iSpectator )
{
    vector< string > iVar;
    
    if( fDebug )
    {
        cout << endl;
        cout << "reading list of variables from TMVA XML file: " << iXMLFile << endl;
    }
    // open TMVA XML file
    // NOTE: extreme dependendence on the structure of the TMVA XML file
    ifstream is;
    is.open( iXMLFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VTMVAEvaluator::getTrainingVariable error: cannot open TMVA weight file: " << iXMLFile << endl;
        return iVar;
    }
    string is_line;
    string iTemp;
    
    int nVar = 0;
    
    while( getline( is, is_line ) )
    {
    
        // number of variables
        if( is_line.find( "NVar=\"" ) != string::npos )
        {
            nVar = atoi( is_line.substr( is_line.find( "NVar=\"" ) + 6, is_line.size() - is_line.find( "\">" ) - 1 ).c_str() );
            if( fDebug )
            {
                cout << "\t reading TMVA XML file: number of variables is " << nVar << endl;
            }
        }
        // AGAIN, NOTE: extreme dependendence on the structure of the TMVA XML file
        if( is_line.find( "Expression=\"" ) != string::npos )
        {
            iVar.push_back( is_line.substr( is_line.find( "Expression=\"" ) + 12, is_line.find( "Label=" ) -
                                            is_line.find( "Expression=\"" ) - 14 ) );
            if( is_line.find( "SpecIndex" ) != string::npos )
            {
                iSpectator.push_back( true );
            }
            else
            {
                iSpectator.push_back( false );
            }
            if( fDebug )
            {
                cout << "\t reading TMVA XML file: new variable: " << iVar.back() << endl;
            }
        }
    }
    is.close();
    
    return iVar;
}

/*

    initialize TMVA readers

*/
bool VTMVAEvaluator::initializeWeightFiles( string iWeightFileName,
        unsigned int iWeightFileIndex_Emin, unsigned int iWeightFileIndex_Emax,
        unsigned int iWeightFileIndex_Zmin, unsigned int iWeightFileIndex_Zmax,
        double iEnergyStepSize, string iInstrumentEpoch,
        string iOptimizationType,
        string iCutID )
{
    //////////////////////////////
    // sanity checks
    if( iWeightFileName.size() == 0 )
    {
        cout << "VTMVAEvaluator::initializeWeightFiles error: no file name" << endl;
        fIsZombie = true;
        return false;
    }
    if( iWeightFileIndex_Emin > iWeightFileIndex_Emax )
    {
        cout << "VTMVAEvaluator::initializeWeightFiles: min energy bin larger than maximum: ";
        cout << iWeightFileIndex_Emin << " > " << iWeightFileIndex_Emax << endl;
        fIsZombie = true;
        return false;
    }
    if( iWeightFileIndex_Zmin > iWeightFileIndex_Zmax )
    {
        cout << "VTMVAEvaluator::initializeWeightFiles: min zenith bin larger than maximum: ";
        cout << iWeightFileIndex_Zmin << " > " << iWeightFileIndex_Zmax << endl;
        fIsZombie = true;
        return false;
    }
    fWeightFileIndex_Emin = iWeightFileIndex_Emin;
    fWeightFileIndex_Emax = iWeightFileIndex_Emax;
    fWeightFileIndex_Zmin = iWeightFileIndex_Zmin;
    fWeightFileIndex_Zmax = iWeightFileIndex_Zmax;
    
    //////////////////////////////
    // reset data vector
    fTMVAData.clear();
    
    /////////////////////////////////////////
    // number of energy and zenith bins bins
    unsigned int iNbinE = iWeightFileIndex_Emax - iWeightFileIndex_Emin + 1;
    unsigned int iNbinZ = iWeightFileIndex_Zmax - iWeightFileIndex_Zmin + 1;
    
    cout << "VTMVAEvaluator::initializeWeightFiles: reading energy and zenith bins from TMVA root files ";
    cout << "(nbinE: " << iNbinE << ", nbinZ: " << iNbinZ << ")" << endl;
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // read energy and zenith binning from root files and check that all neccessary objects are in the file
    unsigned int iMinMissingBin = 0;
    char hname[800];
    for( unsigned int i = 0; i < iNbinE; i++ )
    {
        unsigned int jMinMissingBin = 0;
        for( unsigned int j = 0; j < iNbinZ; j++ )
        {
            bool bGoodRun = true;
            string iFullFileName = setFullMVAFileName( iWeightFileName, 
                                                   iWeightFileIndex_Emin, i,
                                                   iWeightFileIndex_Zmin, j,
                                                   fTMVAMethodName, fTMVAMethodCounter,
                                                   iInstrumentEpoch,
                                                   ".root" );
            string iFullFileNameXML = setFullMVAFileName( iWeightFileName, 
                                                   iWeightFileIndex_Emin, i,
                                                   iWeightFileIndex_Zmin, j,
                                                   fTMVAMethodName, fTMVAMethodCounter,
                                                   iInstrumentEpoch,
                                                   ".weights.xml" );

            TFile* iF = 0;
            if( iFullFileName.size() > 0 )
            {
                 iF = new TFile( iFullFileName.c_str() );
            }
            VTMVARunDataEnergyCut* iEnergyData = 0;
            VTMVARunDataZenithCut* iZenithData = 0;
            if( !iF || iF->IsZombie() )
            {
                bGoodRun = false;
            }
            else
            {
                // energy and zenith bins
                iEnergyData = ( VTMVARunDataEnergyCut* )iF->Get( "fDataEnergyCut" );
                iZenithData = ( VTMVARunDataZenithCut* )iF->Get( "fDataZenithCut" );
                if( !iEnergyData )
                {
                    cout << "No energy cut data: setting goodrun to false" << endl;
                    bGoodRun = false;
                }
                // backwards compatibility
                if( !iZenithData && iInstrumentEpoch != "noepoch" && iInstrumentEpoch != "CTA" )
                {
                    cout << "No zenith cut data: ";
                    cout << " setting goodrun to false" << endl;
                    bGoodRun = false;
                }
                // signal efficiency
                // (note that there are files with different directory structure around)
                sprintf( hname, "Method_%s/%s_%d/MVA_%s_%d_effS", fTMVAMethodName.c_str(),
                         fTMVAMethodName.c_str(), fTMVAMethodCounter,
                         fTMVAMethodName.c_str(), fTMVAMethodCounter );
                         
                if( !iF->Get( hname ) )
                {
                    sprintf( hname, "Method_%s_%d/%s_%d/MVA_%s_%d_effS",
                             fTMVAMethodName.c_str(), fTMVAMethodCounter,
                             fTMVAMethodName.c_str(), fTMVAMethodCounter,
                             fTMVAMethodName.c_str(), fTMVAMethodCounter );
                    if( !iF->Get( hname ) )
                    {
                        cout << "No signal efficiency histogram found (" << hname << ")" << endl;
                        bGoodRun = false;
                    }
                }
            }
            // allow that first files are missing (this happens when there are no training events in the first energy bins)
            if( !bGoodRun )
            {
                if( i == iMinMissingBin || j == jMinMissingBin )
                {
                    cout << "VTMVAEvaluator::initializeWeightFiles() warning: TMVA root file not found or incomplete file (";
                    cout << "ebin " << i << ", zebin " << j << ") " << endl;
                    cout << iFullFileName << endl;
                    if( i == iMinMissingBin )
                    {
                        cout << "  assume this is a low-energy empty bin (bin number " << i << ";";
                        cout << " number of missing bins: " << iMinMissingBin + 1 << ")" << endl;
                        iMinMissingBin++;
                    }
                    if( j == jMinMissingBin )
                    {
                        cout << "  assume this is a zenith empty bin (bin number " << j << ";";
                        cout << " number of missing bins: " << jMinMissingBin + 1 << ")" << endl;
                    }
                    continue;
                }
                else if( i == ( iWeightFileIndex_Emax ) || j == ( iWeightFileIndex_Zmax ) )
                {
                    cout << "VTMVAEvaluator::initializeWeightFiles() warning: TMVA root file not found " << iFullFileName << endl;
                    if( i == ( iWeightFileIndex_Emax ) )
                    {
                        cout << "  assume this is a high-energy empty bin (bin number " << i << ")" << endl;
                        iNbinE--;
                        iWeightFileIndex_Emax--;
                    }
                    if( j == ( iWeightFileIndex_Zmax ) )
                    {
                        cout << "  assume this is a high-zenith empty bin (bin number " << j << ")" << endl;
                        iNbinZ--;
                        iWeightFileIndex_Zmax--;
                    }
                    continue;
                }
                else
                {
                    cout << "VTMVAEvaluator::initializeWeightFiles: warning: problem while initializing energies from TMVA root file ";
                    cout << iFullFileName << endl;
                    cout << "(this might be not a problem if the sensitive energy range of the given array is relatively small)" << endl;
                    continue;
                }
                // continue;
            }
            if( !iEnergyData )
            {
                cout << "VTMVAEvaluator::initializeWeightFiles: warning: problem while reading energies from TMVA root file ";
                cout << iFullFileName << endl;
                fIsZombie = true;
                return false;
            }
            // from here on: expect a good TMVA file
            // initialize one value per energy/zenith bin
            
            // set energy binning:
            //    - one VTMVAEvaluatorData per energy bin
            //    - bins are set for the energy interval read from the root file:
            //      [iEnergyData->fEnergyCut_Log10TeV_min, iEnergyData->fEnergyCut_Log10TeV_max]
            //    - sub-bins given by iEnergyStepSize;
            double e = iEnergyData->fEnergyCut_Log10TeV_min;
            do
            {
                // central data element for this energy bin
                fTMVAData.push_back( new VTMVAEvaluatorData() );
                fTMVAData.back()->fEnergyCut_bin = i;
                fTMVAData.back()->fZenithCut_bin = j + iWeightFileIndex_Zmin;
                // find e_min and e_max
                fTMVAData.back()->fEnergyCut_Log10TeV_min = e;
                if( iEnergyStepSize > 0. )
                {
                    fTMVAData.back()->fEnergyCut_Log10TeV_max = e + iEnergyStepSize;
                }
                else
                {
                    fTMVAData.back()->fEnergyCut_Log10TeV_max = iEnergyData->fEnergyCut_Log10TeV_max;
                }
                e = fTMVAData.back()->fEnergyCut_Log10TeV_max;
                
                // calculate spectral weighted mean energy
                fTMVAData.back()->fSpectralWeightedMeanEnergy_Log10TeV =
                    VMathsandFunctions::getSpectralWeightedMeanEnergy( fTMVAData.back()->fEnergyCut_Log10TeV_min,
                            fTMVAData.back()->fEnergyCut_Log10TeV_max,
                            fSpectralIndexForEnergyWeighting );
                // zenith angle range
                if( iZenithData )
                {
                    fTMVAData.back()->fZenithCut_min = iZenithData->fZenithCut_min;
                    fTMVAData.back()->fZenithCut_max = iZenithData->fZenithCut_max;
                }
                else
                {
                    fTMVAData.back()->fZenithCut_min = 0.;
                    fTMVAData.back()->fZenithCut_max = 89.;
                }
                if( fTMVAData.back()->fZenithCut_max > 89. )
                {
                    fTMVAData.back()->fZenithCut_max = 89.;
                }
                
                fTMVAData.back()->fSignalEfficiency = getSignalEfficiency( iWeightFileIndex_Emin + i,
                                                      fTMVAData.back()->fEnergyCut_Log10TeV_min,
                                                      fTMVAData.back()->fEnergyCut_Log10TeV_max,
                                                      iWeightFileIndex_Zmin + j,
                                                      fTMVAData.back()->fZenithCut_min,
                                                      fTMVAData.back()->fZenithCut_max,
                                                      fTMVAData.size(),
                                                      iEnergyStepSize );
                fTMVAData.back()->fTMVACutValue = getTMVACutValue( iWeightFileIndex_Emin + i,
                                                  fTMVAData.back()->fEnergyCut_Log10TeV_min,
                                                  fTMVAData.back()->fEnergyCut_Log10TeV_max,
                                                  iWeightFileIndex_Zmin + j,
                                                  fTMVAData.back()->fZenithCut_min,
                                                  fTMVAData.back()->fZenithCut_max,
                                                  fTMVAData.size(),
                                                  iEnergyStepSize );
                fTMVAData.back()->fBackgroundEfficiency = -99.;
                fTMVAData.back()->fTMVAOptimumCutValueFound = false;
                fTMVAData.back()->fSourceStrengthAtOptimum_CU = 0.;
                fTMVAData.back()->fAngularContainmentRadius = -99.;
                fTMVAData.back()->fAngularContainmentFraction = -99.;
                sprintf( hname, "bin %d, %.2f < log10(E) < %.2f, %.2f < Ze < %.2f)",
                         ( int )( fTMVAData.size() - 1 ),
                         fTMVAData.back()->fEnergyCut_Log10TeV_min,
                         fTMVAData.back()->fEnergyCut_Log10TeV_max,
                         fTMVAData.back()->fZenithCut_min, 
                         fTMVAData.back()->fZenithCut_max );
                fTMVAData.back()->SetTitle( hname );
                
                sprintf( hname, "MVA%u%u", i, j );
                fTMVAData.back()->fTMVAMethodTag = hname;
                if( iNbinZ > 1 )
                {
                    sprintf( hname, "%u_%u", i, j );
                }
                else
                {
                    sprintf( hname, "%u", i );
                }
                
                fTMVAData.back()->fTMVAMethodTag_2 = hname;
                fTMVAData.back()->fTMVAFileName = iFullFileName;
                fTMVAData.back()->fTMVAFileNameXML = iFullFileNameXML;

                TFile *iTMVAFile = new TFile( fTMVAData.back()->fTMVAFileName.c_str() );
                fTMVAData.back()->hSignalEfficiency = getEfficiencyHistogram( "effS", iTMVAFile, fTMVAData.back()->fTMVAMethodTag_2 );
                fTMVAData.back()->hBackgroundEfficiency = getEfficiencyHistogram( "effB", iTMVAFile, fTMVAData.back()->fTMVAMethodTag_2 );

                if( iEnergyStepSize < 0. )
                {
                    break;
                }
            }
            while( e < ( iEnergyData->fEnergyCut_Log10TeV_max - 0.0001 ) );
            
            iF->Close();
        }//end loop on zenith bins
    }//end loop on energy bins
    
    // after this stage, there should be no energy/zenith bins (both of them are combined)
    
    if( fTMVAData.size() == 0 )
    {
        fIsZombie = true;
        return false;
    }
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        fTMVAData[i]->print();
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // create and initialize TMVA readers
    // loop over all  energy bins: open one weight (XML) file per energy bin
    //looping over spectral energy and zenith angle bins
    for( unsigned int b = 0; b < fTMVAData.size(); b++ )
    {
        fTMVAData[b]->fTMVAReader = new TMVA::Reader();
        if( fDebug )
        {
            cout << "INITIALIZE TMVA file: " << fTMVAData[b]->fTMVAFileName << endl;
        }
        //////////////////////////////////////////
        // set TMVA cut value
        // (optimization later)
        
        // fixed signal efficiency
        if( fTMVACutValueNoVec < -1. && fSignalEfficiencyNoVec > 0. )
        {
            getValuesFromEfficiencyHistograms( b );
        }
        // fixed TMVA cut value
        else if( fTMVACutValueNoVec > -1. )
        {
            fTMVAData[b]->fSignalEfficiency = -99.;
            getValuesFromEfficiencyHistograms( b );
        }
        // no optimization took place
        fTMVAData[b]->fTMVAOptimumCutValueFound = false;
        
        // weight file for this energy bin
        
        if( fDebug )
        {
            cout << "reading TMVA XML weight file: " << fTMVAData[b]->fTMVAFileNameXML << endl;
        }
        
        // get list of training variables
        vector< bool > iVariableIsASpectator;
        vector< string > iTrainingVariables = getTrainingVariables( fTMVAData[b]->fTMVAFileNameXML, iVariableIsASpectator );
        
        // note that the following list of variables must be the same as during training
        for( unsigned int t = 0; t < iTrainingVariables.size(); t++ )
        {
            if( iTrainingVariables[t] == "MSCW" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "MSCW", &fMSCW );
            }
            else if( iTrainingVariables[t] == "MSCL" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "MSCL", &fMSCL );
            }
            else if( iTrainingVariables[t] == "EmissionHeight" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "EmissionHeight", &fEmissionHeight );
            }
            else if( iTrainingVariables[t] == "log10(EmissionHeightChi2)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10(EmissionHeightChi2)", &fEmissionHeightChi2_log10 );
            }
            else if( iTrainingVariables[t] == "NImages" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "NImages", &fNImages );
            }
            else if( iTrainingVariables[t] == "dE" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "dE", &fdES );
            }
            else if( iTrainingVariables[t] == "EChi2" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "EChi2", &fEChi2S );
                fEnergyReconstructionMethod = 0;
            }
            else if( iTrainingVariables[t] == "log10(EChi2)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10(EChi2)", &fEChi2S_log10 );
            }
            else if( iTrainingVariables[t] == "dES" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "dES", &fdES );
            }
            else if( iTrainingVariables[t] == "log10(SizeSecondMax)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10(SizeSecondMax)", &fSizeSecondMax_log10 );
            }
            else if( iTrainingVariables[t] == "EChi2S" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "EChi2S", &fEChi2S );
                fEnergyReconstructionMethod = 1;
            }
            else if( iTrainingVariables[t] == "log10((EChi2S&lt;0)+(EChi2S&gt;0)*EChi2S)" && !iVariableIsASpectator[t] )
            {
                 fTMVAData[b]->fTMVAReader->AddVariable( "log10((EChi2S<0)+(EChi2S>0)*EChi2S)", &fEChi2S_gt0 );
                 fEnergyReconstructionMethod = 1;
            }
            else if( iTrainingVariables[t] == "(EChi2S&lt;=0)" && !iVariableIsASpectator[t] )
            {
                 fTMVAData[b]->fTMVAReader->AddVariable( "(EChi2S<=0)", &fEChi2S_gt0_bool );
            }
            else if( iTrainingVariables[t] == "log10(EChi2S)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10(EChi2S)", &fEChi2S_log10 );
            }
            else if( iTrainingVariables[t] == "(Xoff*Xoff+Yoff*Yoff)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "(Xoff*Xoff+Yoff*Yoff)", &fTheta2 );
                setTMVAThetaCutVariable( true );
            }
            else if( iTrainingVariables[t] == "sqrt(Xcore*Xcore+Ycore*Ycore)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "sqrt(Xcore*Xcore+Ycore*Ycore)", &fCoreDist );
            }
            // disp below
            else if( iTrainingVariables[t] == "log10(DispDiff)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10(DispDiff)", &fDispDiff_log10 );
            }
            else if( iTrainingVariables[t] == "log10((DispDiff&lt;=0)+(DispDiff&gt;0.)*DispDiff)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "log10((DispDiff<=0)+(DispDiff>0.)*DispDiff)", &fDispDiff_gt0 );
            }
            else if( iTrainingVariables[t] == "(DispDiff&lt;=0)" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "(DispDiff<=0)", &fDispDiff_gt0_bool );
            }
            // Note: assume not more then 3 different telescope types
            else if( iTrainingVariables[t] == "NImages_Ttype[0]" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "NImages_Ttype[0]", &fImages_Ttype[0] );
            }
            else if( iTrainingVariables[t] == "NImages_Ttype[1]" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "NImages_Ttype[1]", &fImages_Ttype[1] );
            }
            else if( iTrainingVariables[t] == "NImages_Ttype[2]" && !iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddVariable( "NImages_Ttype[2]", &fImages_Ttype[2] );
            }
            else if( iVariableIsASpectator[t] )
            {
                fTMVAData[b]->fTMVAReader->AddSpectator( iTrainingVariables[t].c_str(), &fDummy );
            }
        }
        if( fDebug )
        {
            cout << "Following " << iTrainingVariables.size();
            cout << " variables have been found and are used for TMVA separation: " << endl;
            for( unsigned int t = 0; t < iTrainingVariables.size(); t++ )
            {
                cout << "\t" << iTrainingVariables[t];
                if( iVariableIsASpectator[t] )
                {
                    cout << " (spectator)";
                }
                cout << endl;
            }
        }
        if( !fTMVAData[b]->fTMVAReader->BookMVA( fTMVAData[b]->fTMVAMethodTag_2.c_str(), fTMVAData[b]->fTMVAFileNameXML.c_str() ) )
        {
            cout << "VTMVAEvaluator::initializeWeightFiles: error while initializing TMVA reader from weight file ";
            cout << fTMVAData[b]->fTMVAFileNameXML << endl;
            fIsZombie = true;
            return false;
        }
        /////////////////////////////////////////////////////////
        // get optimal signal efficiency (from maximum signal/noise ratio)
        /////////////////////////////////////////////////////////
        if( fParticleNumberFileName.size() > 0 )
        {
            cout << endl;
            cout << "======================= optimize sensitivity " << b << " =======================" << endl;
            if( !optimizeSensitivity( b, iOptimizationType, iInstrumentEpoch ) )
            {
                cout << "VTMVAEvaluator::initializeWeightFiles: error while calculating optimized sensitivity" << endl;
                return false;
            }
            cout << "======================= end optimize sensitivity =======================" << endl;
            cout << endl;
        }
        
    }
    
    // smooth and Interpolate
    if( fParticleNumberFileName.size() > 0 && fSmoothAndInterpolateMVAValues )
    {
        smoothAndInterpolateMVAValue( 
                 iWeightFileIndex_Emin, iWeightFileIndex_Emax, 
                 iWeightFileIndex_Zmin, iWeightFileIndex_Zmax, 
                 iEnergyStepSize );
    }
    
    // print some info to screen
    cout << "VTMVAEvaluator: Initialized " << fTMVAData.size() << " MVA readers " << endl;
    
    fillTMVAEvaluatorResults( iCutID );
    
    return true;
}

/*
 *  copy TMVA data vectors
 *
 * VTMVAEvaluatorResults are written to output files and used in
 * anasum, effective areas code and sensitivity calculation
 */
void VTMVAEvaluator::fillTMVAEvaluatorResults( string iCutID )
{
    if( !fTMVAEvaluatorResults )
    {
        fTMVAEvaluatorResults = new VTMVAEvaluatorResults;
        fTMVAEvaluatorResults->SetName( ("TMVAEvaluatorResults"+iCutID).c_str() );
    }
    if( fTMVAEvaluatorResults )
    {
        fTMVAEvaluatorResults->fTMVAData = fTMVAData;
    }
}

TH1D* VTMVAEvaluator::getEfficiencyHistogram( string iName, TFile* iF, string iMethodTag_2 )
{
    if( !iF )
    {
        return 0;
    }
    
    char hname[800];
    sprintf( hname, "Method_%s/%s_%s/MVA_%s_%s_%s", fTMVAMethodName.c_str(),
             fTMVAMethodName.c_str(), iMethodTag_2.c_str(),
             fTMVAMethodName.c_str(), iMethodTag_2.c_str(), iName.c_str() );
             
    // read signal efficiency histogram
    TH1D* eff = ( TH1D* )iF->Get( hname );
    if( !eff )
    {
        sprintf( hname, "Method_%s/%s_%d/MVA_%s_%d_%s", fTMVAMethodName.c_str(),
                 fTMVAMethodName.c_str(), fTMVAMethodCounter,
                 fTMVAMethodName.c_str(), fTMVAMethodCounter, iName.c_str() );
        eff = ( TH1D* )iF->Get( hname );
        if( !eff )
        {
            sprintf( hname, "Method_%s_%d/%s_%d/MVA_%s_%d_%s", 
                     fTMVAMethodName.c_str(), fTMVAMethodCounter,
                     fTMVAMethodName.c_str(), fTMVAMethodCounter,
                     fTMVAMethodName.c_str(), fTMVAMethodCounter, iName.c_str() );
            eff = ( TH1D* )iF->Get( hname );
            if( !eff )
            {
                cout << "VTMVAEvaluator::getEfficiencyHistogram() error finding efficiency histogram " << hname;
                cout << " from " << iF->GetName() << endl;
                return 0;
            }
        }
    }
    return eff;
}

/*

   get TMVA cut values
   (e.g. signal efficiency for a given MVA cut or
         MVA cut for a given signal efficiency

*/

bool VTMVAEvaluator::getValuesFromEfficiencyHistograms( unsigned int b )
{
    if( b >= fTMVAData.size() )
    {
        return false;
    }
    // make sure that default values are set
    if( fTMVAData[b]->fTMVACutValue > -1. )
    {
        fTMVAData[b]->fSignalEfficiency = fTMVAData[b]->fBackgroundEfficiency = -99.;
    }
    else if( fTMVAData[b]->fSignalEfficiency > 0. )
    {
        fTMVAData[b]->fTMVACutValue = fTMVAData[b]->fBackgroundEfficiency = -99.;
    }
    else if( fTMVAData[b]->fBackgroundEfficiency > 0. )
    {
        fTMVAData[b]->fTMVACutValue = fTMVAData[b]->fSignalEfficiency = -99.;
    }
    
    // check file name for consistency
    if( fTMVAData[b]->fTMVAFileName.size() == 0 )
    {
        return false;
    }
    
    TFile* iTMVAFile = new TFile( fTMVAData[b]->fTMVAFileName.c_str() );
    if( iTMVAFile->IsZombie() )
    {
        cout << "VTMVAEvaluator::getValuesFromEfficiencyHistograms() ";
        cout << "error reading TMVA root file: ";
        cout << fTMVAData[b]->fTMVAFileName << endl;
        return false;
    }
    TH1D* effS = getEfficiencyHistogram( "effS", iTMVAFile, fTMVAData[b]->fTMVAMethodTag_2 );
    TH1D* effB = getEfficiencyHistogram( "effB", iTMVAFile, fTMVAData[b]->fTMVAMethodTag_2 );
    if( !effS || !effB )
    {
        return false;
    }
    
    if( fDebug )
    {
        cout << "VTMVAEvaluator::getValuesFromEfficiencyHistograms: evaluating " << iTMVAFile->GetName() << endl;
    }
    // get MVA cut for a given signal efficiency
    if( fTMVAData[b]->fSignalEfficiency > 0. )
    {
        fTMVAData[b]->fTMVACutValue = effS->GetBinCenter( effS->FindLastBinAbove( fTMVAData[b]->fSignalEfficiency ) );
        fTMVAData[b]->fBackgroundEfficiency = effB->GetBinContent( effB->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        
        cout << "TMVA CUT VALUE FOR SIGNAL EFFICIENCY " << fTMVAData[b]->fSignalEfficiency << ": " << fTMVAData[b]->fTMVACutValue;
        cout << " (bin " << effS->FindLastBinAbove( fTMVAData[b]->fSignalEfficiency ) << ")" << endl;
    }
    // get signal efficiency from histogram using TMVA cut value from graph
    else if( fTMVAData[b]->fZenithCut_bin < fTMVACutValueGraph.size() )
    {
        fTMVAData[b]->fTMVACutValue = fTMVACutValueGraph[fTMVAData[b]->fZenithCut_bin]->Eval(
                            fTMVAData[b]->fSpectralWeightedMeanEnergy_Log10TeV );
        if( fTMVAData[b]->fZenithCut_bin < fTMVASignalEfficencyGraph.size()
            && fTMVASignalEfficencyGraph[fTMVAData[b]->fZenithCut_bin] )
        {
              fTMVAData[b]->fSignalEfficiency = fTMVASignalEfficencyGraph[fTMVAData[b]->fZenithCut_bin]->Eval(
                            fTMVAData[b]->fSpectralWeightedMeanEnergy_Log10TeV );
        }
        else
        {
            fTMVAData[b]->fSignalEfficiency     = effS->GetBinContent( effS->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        }
        if( fTMVAData[b]->fZenithCut_bin < fTMVABackgroundEfficencyGraph.size()
            && fTMVABackgroundEfficencyGraph[fTMVAData[b]->fZenithCut_bin] )
        {
              fTMVAData[b]->fBackgroundEfficiency = fTMVABackgroundEfficencyGraph[fTMVAData[b]->fZenithCut_bin]->Eval(
                            fTMVAData[b]->fSpectralWeightedMeanEnergy_Log10TeV );
        }
        else
        {
            fTMVAData[b]->fBackgroundEfficiency = effB->GetBinContent( effB->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        }

        if( fDebug )
        {
            cout << "Signal efficiency for TMVA cut value " << fTMVAData[b]->fTMVACutValue << ": " << fTMVAData[b]->fSignalEfficiency;
            cout << " (bin " << effS->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) << ")" << endl;
        }
    }
    // get signal efficiency from histogram
    else if( fTMVAData[b]->fTMVACutValue > -1. )
    {
        fTMVAData[b]->fSignalEfficiency     = effS->GetBinContent( effS->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        fTMVAData[b]->fBackgroundEfficiency = effB->GetBinContent( effB->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        
        if( fDebug )
        {
            cout << "Signal efficiency for TMVA cut value " << fTMVAData[b]->fTMVACutValue << ": " << fTMVAData[b]->fSignalEfficiency;
            cout << " (bin " << effS->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) << ")" << endl;
        }
    }
    // get MVA cut for a given background efficiency
    else if( fTMVAData[b]->fBackgroundEfficiency > 0. )
    {
        fTMVAData[b]->fTMVACutValue = effB->GetBinCenter( effB->FindLastBinAbove( fTMVAData[b]->fBackgroundEfficiency ) );
        fTMVAData[b]->fSignalEfficiency = effS->GetBinContent( effS->GetXaxis()->FindBin( fTMVAData[b]->fTMVACutValue ) );
        
        cout << "TMVA CUT VALUE FOR SIGNAL EFFICIENCY " << fTMVAData[b]->fBackgroundEfficiency << ": " << fTMVAData[b]->fTMVACutValue;
        cout << " (bin " << effB->FindLastBinAbove( fTMVAData[b]->fBackgroundEfficiency ) << ")" << endl;
    }
    
    iTMVAFile->Close();
    
    return true;
}

/*!

    evaluate this event using the MVA and return passed/not passed

*/
bool VTMVAEvaluator::evaluate()
{
    if( fDebug )
    {
        cout << "VTMVAEvaluator::evaluate (" << fData << ")" << endl;
    }
    // copy event data
    if( fData )
    {
        fNImages        = ( float )fData->getNImages();
        fMSCW           = fData->MSCW;
        fMSCL           = fData->MSCL;
        fMWR            = fData->MWR;
        fMLR            = fData->MLR;
        fEmissionHeight = fData->EmissionHeight;
        if( fData->EmissionHeightChi2 > 0. )
        {
            fEmissionHeightChi2_log10 = TMath::Log10( fData->EmissionHeightChi2 );
        }
        else
        {
            fEmissionHeightChi2_log10 = -10.;    // !!! not clear what the best value is
        }
        // fill according of energy reconstruction method
        fEChi2S          = fData->getEnergyChi2();
        if( fEChi2S > 0. )
        {
            fEChi2S_log10 = TMath::Log10( fEChi2S );
            fEChi2S_gt0 = TMath::Log10( fEChi2S );
            fEChi2S_gt0_bool = 0.;
        }
        else
        {
            fEChi2S_log10 = 0.;    // !!! not clear what the best value is
            fEChi2S_gt0 = 1.;
            fEChi2S_gt0_bool = 1.;
        }
        fdES = fData->getEnergyDelta();
        
        fSizeSecondMax_log10 = fData->SizeSecondMax;
        if( fSizeSecondMax_log10 > 0. )
        {
            fSizeSecondMax_log10 = TMath::Log10( fSizeSecondMax_log10 );
        }
        else
        {
            fSizeSecondMax_log10 = 0.;    // !!! not clear what the best value is
        }
        // for no theta2: set this value extremely small
        if( fTMVAIgnoreTheta2Cut )
        {
            fTheta2 = 1.e-30;
        }
        else
        {
            fTheta2 = fData->getXoff() * fData->getXoff() + fData->getYoff() * fData->getYoff();
        }
        fCoreDist = sqrt( fData->getXcore_M() * fData->getXcore_M() + fData->getYcore_M() * fData->getYcore_M() );
        if( fData->DispDiff > 0. )
        {
            fDispDiff_log10 = log10( fData->DispDiff );
            fDispDiff_gt0  = log10( fData->DispDiff );
            fDispDiff_gt0_bool = 0.;
        }
        else
        {
            fDispDiff_log10 = 100.;
            fDispDiff_gt0 = 1.;
            fDispDiff_gt0_bool = 1.;
        }
        if( fData->NTtype < VDST_MAXTELESCOPES )
        {
            for( int i = 0; i < fData->NTtype; i++ )
            {
                fImages_Ttype[i] = ( float )fData->NImages_Ttype[i];
            }
        }
    }
    else
    {
        return false;
    }
    
    // valid energy is required for TMVA evaluation
    if( fData->getEnergy_TeV() < 0 )
    {
        return false;
    }
    
    // find correct bin (e.g., depending on energy or zenith)
    unsigned int iDataBin = getDataBin( fData->getEnergy_Log10(), fData->getZe() );
    
    fTMVA_EvaluationResult = -99.;
    
    if( iDataBin < fTMVAData.size() )
    {
        if( fDebug )
        {
            cout << "VTMVAEvaluator::evaluate: data bin " << iDataBin;
            cout << ", MVA Method Tag " << fTMVAData[iDataBin]->fTMVAMethodTag;
            cout << ", MVA Cut value " << fTMVAData[iDataBin]->fTMVACutValue;
            cout << endl;
        }
        
        // evaluate MVA for this event
        fTMVA_EvaluationResult = fTMVAData[iDataBin]->fTMVAReader->EvaluateMVA( fTMVAData[iDataBin]->fTMVAMethodTag_2 );

        // evaluate interpolate MVA for this event
        // fTMVA_EvaluationResult = evaluateInterPolateMVA( fData->getEnergy_Log10(), fData->getZe(), iDataBin );

        // apply MVA cut
        // -> can be either
        //    i) a 1D Graph (vector of graphs in ze)
        //    ii) a map for values in energy and zenith
        if( fTMVACutValueGraph.size() > 0
        && fTMVAData[iDataBin]->fZenithCut_bin < fTMVACutValueGraph.size()
        && fTMVACutValueGraph[fTMVAData[iDataBin]->fZenithCut_bin] )
        {
            // fTMVACutValueGraph is a vector of the size of ze bins
            if( fTMVA_EvaluationResult < fTMVACutValueGraph[fTMVAData[iDataBin]->fZenithCut_bin]->Eval( fData->getEnergy_Log10() ) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        else if( fTMVA_EvaluationResult < fTMVAData[iDataBin]->fTMVACutValue )
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    return false;
}

/*
 * calculate MVA for a given event
 *
 * interpolate over different bins and calculate
 * the weighted average
 *
 */
double VTMVAEvaluator::evaluateInterPolateMVA( double iErec_log10TeV, double iZe, unsigned int iDataBin )
{
    vector< double > iW = getDataBinWeights( iErec_log10TeV, iZe );
    double iMVA = 0.;
    double iMVA_tot = 0.;

    // efficiency check: don't interpolate if there is only one energy bin
    // for the TMVAs
    if( fWeightFileIndex_Emax - fWeightFileIndex_Emin == 1 )
    {
        return fTMVAData[iDataBin]->fTMVAReader->EvaluateMVA( fTMVAData[iDataBin]->fTMVAMethodTag_2 );
    }

    for( unsigned int w = 0; w < iW.size(); w++ )
    {
        // ignore very small weights
        if( w < fTMVAData.size() && iW[w] > 0.001
        && fTMVAData[w]->fTMVAReader )
        {
             double t = fTMVAData[w]->fTMVAReader->EvaluateMVA( fTMVAData[w]->fTMVAMethodTag_2 );
             iMVA += iW[w] * t;
             iMVA_tot += iW[w];
        }
    }
    if( iMVA_tot > 0. )
    {
        return iMVA / iMVA_tot;
    }

    return -99.;
}

/*
 * for given E and Ze, get weights depending on the 
 * distance in energy
 *
 * no interpolation in zenith angle
 *
 * return vector always same legnth is fTMVAData
 *
 */
vector< double > VTMVAEvaluator::getDataBinWeights( double iErec_log10TeV, unsigned int iZeBin )
{
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
           if( fTMVAData[i]->fZenithCut_bin == iZeBin )
           {
                return getDataBinWeights( iErec_log10TeV,
                     0.5 * ( fTMVAData[i]->fZenithCut_min
                           + fTMVAData[i]->fZenithCut_max ) );
           }
    }
    vector< double > iW( fTMVAData.size(), 0. );
    return iW;
}

/*
 * for given E and Ze, get weights depending on the 
 * distance in energy
 *
 * no interpolation in zenith angle
 *
 * return vector always same legnth is fTMVAData
 *
 */
vector< double > VTMVAEvaluator::getDataBinWeights( double iErec_log10TeV, double iZe )
{
    vector< double > iW( fTMVAData.size(), 0. );

    double i_tot = 0.;

    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        // get zenith bin for the current zenith (read from fData)
        if( ( iZe > fTMVAData[i]->fZenithCut_min && iZe <= fTMVAData[i]->fZenithCut_max ) || iZe < -998. )
        {
            // mean energy of this energy bin (possibly spectral weighted)
            if( iErec_log10TeV < fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV )
            {
                if( TMath::Abs(fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV-fTMVAData[i]->fEnergyCut_Log10TeV_min) > 0. )
                {
                    iW[i] = 1. - tanh( 2.* TMath::Abs( fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - iErec_log10TeV )
                         / TMath::Abs(fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV-fTMVAData[i]->fEnergyCut_Log10TeV_min)  );
                }
                else
                {
                    iW[i] = 1.;
                }
            }
            else
            {
                if( TMath::Abs(fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV+fTMVAData[i]->fEnergyCut_Log10TeV_max) > 0. )
                {
                    iW[i] = 1. - tanh( 2. * TMath::Abs( fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - iErec_log10TeV )
                          / TMath::Abs(fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV+fTMVAData[i]->fEnergyCut_Log10TeV_max) );
                }
                else
                {
                    iW[i] = 1.;
                }
            }
            i_tot += iW[i];
        }
     }
     if( i_tot > 0. )
     {
         for( unsigned int i = 0; i < iW.size(); i++ )
         {
             iW[i] = iW[i] / i_tot;
         }
     }
     return iW;
}

/*
 *   get bin number for current event
 *
 *   --> ZeBin given
 */
unsigned int VTMVAEvaluator::getDataBin( double iErec_log10TeV, unsigned int iZeBin )
{
      for( unsigned int i = 0; i < fTMVAData.size(); i++ )
      {
           if( fTMVAData[i]->fZenithCut_bin == iZeBin )
           {
                return getDataBin( iErec_log10TeV,
                     0.5 * ( fTMVAData[i]->fZenithCut_min
                           + fTMVAData[i]->fZenithCut_max ) );
           }
      }

      // ze bin not found
      return 9999;
}


/*
 *   get bin number for current event
 */
unsigned int VTMVAEvaluator::getDataBin( double iErec_log10TeV, double iZe )
{
    double       i_Diff_Energy = 1.e10;           // difference between energy of current event and mean bin energy
    double       iMeanEnergy = 0.;
    double       iMeanEnergy_min = 1.e10;
    
    unsigned int iBin = 9999;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        //   get zenith bin for the current zenith (read from fData)
        if( ( iZe > fTMVAData[i]->fZenithCut_min && iZe <= fTMVAData[i]->fZenithCut_max ) || iZe < -998. )
        {
            // mean energy of this energy bin (possibly spectral weighted)
            iMeanEnergy = VMathsandFunctions::getMeanEnergyInBin( 2, fTMVAData[i]->fEnergyCut_Log10TeV_min,
                          fTMVAData[i]->fEnergyCut_Log10TeV_max,
                          fSpectralIndexForEnergyWeighting );
            // check which energy bin is closest (on log scale)
            if( TMath::Abs( iMeanEnergy - iErec_log10TeV ) < i_Diff_Energy )
            {
                i_Diff_Energy = TMath::Abs( iMeanEnergy - iErec_log10TeV );
                iBin = i;
                iMeanEnergy_min = iMeanEnergy;
            }
        }
    }
    if( fDebug && iBin < fTMVAData.size() )
    {
        cout << "VTMVAEvaluator::getDataBin: " << iBin << endl;
        fTMVAData[iBin]->print();
        cout << "\t mean energy " << iMeanEnergy_min;
        cout << ", log10 energy " << iErec_log10TeV << "\t" << i_Diff_Energy ;
        cout << "\t" << fSpectralIndexForEnergyWeighting << endl;
    }
    
    return iBin;
}

bool VTMVAEvaluator::initializeDataStrutures( CData* iC )
{
    fData = iC;
    
    if( !fData )
    {
        fIsZombie = true;
        return false;
    }
    
    return true;
}

/*

   get energy dependent theta2 cut

*/
double VTMVAEvaluator::getOptimalTheta2Cut( double iEnergy_log10TeV, double iZe )
{
    if( fTMVAData.size() == 0 )
    {
        cout << "VTMVAEvaluator::getOptimalTheta2Cut error: empty data vector" << endl;
        return -99.;
    }
    
    unsigned int iDataBin = getDataBin( iEnergy_log10TeV, iZe );
    
    if( iDataBin < fTMVAData.size() )
    {
        return ( fTMVAData[iDataBin]->fAngularContainmentRadius * fTMVAData[iDataBin]->fAngularContainmentRadius );
    }
    
    return 0.;
}

/*
   return a graph with all the theta2 cuts

   (is a memory leak...)

*/
TGraph* VTMVAEvaluator::getOptimalTheta2Cut_Graph()
{
    if( fTMVAData.size() == 0 )
    {
        cout << "VTMVAEvaluator::getOptimalTheta2Cut_Graph error: empty data vector" << endl;
        return 0;
    }
    
    // sort
    map< double, double > i_AngContaint;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] && fTMVAData[i]->fAngularContainmentRadius > 1.e-10 )
        {
            i_AngContaint[0.5 * ( fTMVAData[i]->fEnergyCut_Log10TeV_min + fTMVAData[i]->fEnergyCut_Log10TeV_max )] = fTMVAData[i]->fAngularContainmentRadius;
        }
    }
    
    // fill the graph - energy is at the spectral weighted energy
    TGraph* g = new TGraph( 1 );
    unsigned int z = 0;
    for( map<double, double>::iterator it = i_AngContaint.begin(); it != i_AngContaint.end(); ++it )
    {
        g->SetPoint( z, it->first, it->second );
        z++;
    }
    
    return g;
}

/*
 * plot signal and background efficiencies
 *
 */
TGraphAsymmErrors* VTMVAEvaluator::plotSignalAndBackgroundEfficiencies( 
                          bool iLogY, double iYmin, 
                          double iMVA_min, double iMVA_max )
{
    if( fTMVAData.size() == 0 )
    {
        cout << "TMVAEvaluator::plotSignalAndBackgroundEfficiencies error: signal efficiency vector with size 0" << endl;
        return 0;
    }
    
    // fill graphs
    TGraphAsymmErrors* igSignal = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* igSignalOpt = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* igBck = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* igBckOpt = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* igCVa = new TGraphAsymmErrors( 1 );
    TGraphAsymmErrors* igCVaOpt = new TGraphAsymmErrors( 1 );
    
    unsigned int z_opt = 0;
    unsigned int z_noOpt = 0;
    
    double iMinBck = 1.;
    
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( !fTMVAData[i] )
        {
            continue;
        }
        
        if( fTMVAData[i]->fSignalEfficiency > 0. && fTMVAData[i]->fBackgroundEfficiency > 0. )
        {
            if( fTMVAData[i]->fTMVAOptimumCutValueFound )
            {
                igSignal->SetPoint( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fSignalEfficiency );
                igSignal->SetPointEXlow( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igSignal->SetPointEXhigh( z_opt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                igBck->SetPoint( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fBackgroundEfficiency );
                igBck->SetPointEXlow( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igBck->SetPointEXhigh( z_opt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                igCVa->SetPoint( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fTMVACutValue );
                igCVa->SetPointEXlow( z_opt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igCVa->SetPointEXhigh( z_opt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                z_opt++;
            }
            else if( fTMVAData[i]->fTMVACutValue > -90. )
            {
                igSignalOpt->SetPoint( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fSignalEfficiency );
                igSignalOpt->SetPointEXlow( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igSignalOpt->SetPointEXhigh( z_noOpt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                igBckOpt->SetPoint( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fBackgroundEfficiency );
                igBckOpt->SetPointEXlow( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igBckOpt->SetPointEXhigh( z_noOpt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                igCVaOpt->SetPoint( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fTMVACutValue );
                igCVaOpt->SetPointEXlow( z_noOpt, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                igCVaOpt->SetPointEXhigh( z_noOpt, fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                
                z_noOpt++;
            }
        }
        else
        {
            cout << "VTMVAEvaluator::plotSignalAndBackgroundEfficiencies: ";
            cout << "signal / background efficiency histograms not found in " << endl;
            cout << fTMVAData[i]->fTMVAFileName << endl;
            cout << "Signal and background efficiency: ";
            cout << fTMVAData[i]->fSignalEfficiency << "\t" << fTMVAData[i]->fBackgroundEfficiency << endl;
        }
        if( fTMVAData[i]->fBackgroundEfficiency < iMinBck )
        {
            iMinBck = fTMVAData[i]->fBackgroundEfficiency;
        }
    }
    
    // plot everything
    TCanvas* iCanvas = new TCanvas( "cSignalAndBackgroundEfficiencies", "signal and background efficiencies", 10, 10, 400, 400 );
    iCanvas->SetGridx( 0 );
    iCanvas->SetGridy( 0 );
    iCanvas->SetLeftMargin( 0.13 );
    if( iLogY )
    {
        iCanvas->SetLogy();
    }
    else
    {
        iCanvas->SetLogy( 0 );
    }
    iCanvas->Draw();
    
    double iE_min =  1.e99;
    double iE_max = -1.e99;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] && fTMVAData[i]->fEnergyCut_Log10TeV_min < iE_min )
        {
            iE_min = fTMVAData[i]->fEnergyCut_Log10TeV_min;
        }
        if( fTMVAData[i] && fTMVAData[i]->fEnergyCut_Log10TeV_max > iE_max )
        {
            iE_max = fTMVAData[i]->fEnergyCut_Log10TeV_max;
        }
    }
    
    TH1D* hnull = new TH1D( "hnullcSignalAndBackgroundEfficiencies", "", 100, iE_min, iE_max );
    hnull->SetStats( 0 );
    hnull->SetXTitle( "energy [TeV]" );
    hnull->SetYTitle( "signal/background efficiency" );
    hnull->SetMinimum( iYmin );
    hnull->SetMaximum( 1. );
    plot_nullHistogram( iCanvas, hnull, false, false, 1.5, iE_min, iE_max );
    
    setGraphPlottingStyle( igSignal, 1, 1., 20 );
    setGraphPlottingStyle( igSignalOpt, 1, 1., 24 );
    if( igBck )
    {
        setGraphPlottingStyle( igBck, 2, 1., 21 );
    }
    if( igBckOpt )
    {
        setGraphPlottingStyle( igBckOpt, 2, 1., 25 );
    }
    
    igSignal->Draw( "pl" );
    if( z_noOpt > 0 )
    {
        igSignalOpt->Draw( "pl" );
    }
    if( igBck )
    {
        igBck->Draw( "pl" );
    }
    if( igBckOpt && z_noOpt > 0 )
    {
        igBckOpt->Draw( "pl" );
    }
    if( fPrintPlotting )
    {
        iCanvas->Print( "MVA-SignalBackgroundEfficiency.pdf" );
    }

    // plot MVA cut value
    if( igCVa )
    {
        TCanvas* iCVACanvas = new TCanvas( "iCVACanvas", "MVA cut value", 500, 10, 400, 400 );
        iCVACanvas->SetGridx( 0 );
        iCVACanvas->SetGridy( 0 );
        
        TH1D* hnull = new TH1D( "hnullcMVACuts", "", 100, iE_min, iE_max );
        hnull->SetStats( 0 );
        hnull->SetXTitle( "energy [TeV]" );
        hnull->SetYTitle( "MVA cut variable" );
        hnull->SetMinimum( iMVA_min );
        hnull->SetMaximum( iMVA_max );
        plot_nullHistogram( iCanvas, hnull, false, false, 1.3, iE_min, iE_max );
        setGraphPlottingStyle( igCVa, 1, 1., 20 );
        igCVa->Draw( "p" );
        if( igCVaOpt && z_noOpt > 0 )
        {
            setGraphPlottingStyle( igCVaOpt, 1, 1., 24 );
            igCVaOpt->Draw( "p" );
        }
        if( fPrintPlotting )
        {
            iCVACanvas->Print( "MVA-MVACut.pdf" );
        }
    }
    
    return igCVa;
}

/*
 * constant (not dependend on energy and zenith angle)
 * signal efficiency
 */
void VTMVAEvaluator::setSignalEfficiency( double iSignalEfficiency )
{
    fSignalEfficiencyMap[9999] = iSignalEfficiency;
    
    fSignalEfficiencyNoVec = iSignalEfficiency;
}

/*
 * signal efficiency
 * dependend on energy and zenith angle
 */
void VTMVAEvaluator::setSignalEfficiency( map< unsigned int, double > iSignalEfficiencyMap )
{
    fSignalEfficiencyMap = iSignalEfficiencyMap;
    
    // set to value >0: indicates that signal efficiency had been set from outsite
    fSignalEfficiencyNoVec = 1.;
}

/*
 * return all TMVA cut graphs as a single vector
 * (MVA, signal and background efficiency)
 */
vector< TGraphAsymmErrors* > VTMVAEvaluator::getTMVACutValueGraphs()
{
     vector< TGraphAsymmErrors* > iG;
     for( unsigned int i = 0; i < fTMVACutValueGraph.size(); i++ )
     {
          if( fTMVACutValueGraph[i] )
          {
               iG.push_back( fTMVACutValueGraph[i] );
          }
     }
     for( unsigned int i = 0; i < fTMVASignalEfficencyGraph.size(); i++ )
     {
          if( fTMVASignalEfficencyGraph[i] )
          {
               iG.push_back( fTMVASignalEfficencyGraph[i] );
          }
     }
     for( unsigned int i = 0; i < fTMVABackgroundEfficencyGraph.size(); i++ )
     {
          if( fTMVABackgroundEfficencyGraph[i] )
          {
               iG.push_back( fTMVABackgroundEfficencyGraph[i] );
          }
     }
     return iG;
}

/*
 * optimised cuts are given in energy dependent graphs
 *
 * naming of graphs:
 *   TMVACutValue_ze0, TMVACutValue_ze1 for different zenith angle bins
 *
 *   optimally apply smoothing here:
 *   - set iCutGraphSmoothing > 0 (best to roughly Eres, 0.1-0.2
 *
 */
vector< TGraphAsymmErrors* > VTMVAEvaluator::setTMVACutValueFromGraph( string iFileName, 
                                                                       double iCutGraphSmoothing,
                                                                       double iCutGraphSmoothingMax,
                                                                       double iCutGraphConstantCutEnergy_TeV,
                                                                       bool iSmoothTMVAGraph )
{
    fTMVACutValueGraph.clear();
    fTMVASignalEfficencyGraph.clear();
    fTMVABackgroundEfficencyGraph.clear();

    iFileName = VUtilities::testFileLocation( iFileName, "", true );
    fTMVACutValueFile =  new TFile( iFileName.c_str() );
    if( fTMVACutValueFile->IsZombie() )
    {
         cout << "VTMVAEvaluator::setTMVACutValueFromGraph error: ";
         cout << "file not found" << endl;
         cout << "\t" << iFileName << endl;
         cout << "exiting..." << endl;
         exit( EXIT_FAILURE );
    }
    TGraphAsymmErrors *iG = 0;
    char hname[100];
    // lazy loop; assume not more than 99 zenith bins
    for( unsigned int i = 0; i < 99; i++ )
    {
         // smooth TMVA graph
         // note: this might lead to steps in the 
         // signal efficiency vs energy function
         // when changing between energy bins
         if( iSmoothTMVAGraph )
         {
             sprintf( hname, "TMVACutValue_ze%u", i );
             iG = (TGraphAsymmErrors*)fTMVACutValueFile->Get( hname );
             if( !iG )
             {
                  break;
             }
             // smooth graph roughly with bin width / energy resolutions
             if( iCutGraphSmoothing > 0 && iG->GetN() > 1 )
             {
                 fTMVACutValueGraph.push_back(
                          smoothMVAGraph( iG,
                          iCutGraphSmoothing,
                          iCutGraphSmoothingMax,
                          iCutGraphConstantCutEnergy_TeV,
                          hname, i ) );
                 fTMVASignalEfficencyGraph.push_back(
                      fillSmoothedEfficencyGraph(
                        fTMVACutValueGraph.back(),
                        i, 
                        true ) );
                 fTMVABackgroundEfficencyGraph.push_back(
                      fillSmoothedEfficencyGraph(
                        fTMVACutValueGraph.back(),
                        i, 
                        false) );
             }
         }
         // smooth signal efficiency
         else
         {
             sprintf( hname, "SignalEfficiency_ze%u", i );
             iG = (TGraphAsymmErrors*)fTMVACutValueFile->Get( hname );
             if( !iG )
             {
                  break;
             }
             // smooth graph roughly with bin width / energy resolutions
             if( iCutGraphSmoothing > 0 && iG->GetN() > 1 )
             {
                 sprintf( hname, "SignalEfficiency_ze%u", i );
                 /*fTMVASignalEfficencyGraph.push_back(
                          smoothMVAGraph( iG,
                          iCutGraphSmoothing,
                          1.,                  // -> max signal efficiency of 100%
                          iCutGraphConstantCutEnergy_TeV,
                          hname, i ) );
                 fTMVACutValueGraph.push_back(
                      fillSmoothedMVACutGraph(
                          fTMVASignalEfficencyGraph.back(),
                          i ) );*/
                 
                 vector< TGraphAsymmErrors* > iRG = 
                      smoothSignalEfficiencyMVAGraph( iG,
                          iCutGraphSmoothing,
                          1.,                  // -> max signal efficiency of 100%
                          iCutGraphConstantCutEnergy_TeV,
                          hname, i );
                 fTMVASignalEfficencyGraph.push_back( iRG[0] );
                 fTMVACutValueGraph.push_back( iRG[1] ); 
                 fTMVABackgroundEfficencyGraph.push_back(
                      fillSmoothedEfficencyGraph(
                        fTMVACutValueGraph.back(),
                        i, 
                        false ) );
              }
         }
    }  
    if( fTMVACutValueGraph.size() == 0 )
    {
         cout << "VTMVAEvaluator::setTMVACutValueFromGraph error: ";
         cout << "no graph with mva cuts found in " << endl;
         cout << "\t" << iFileName << endl;
         return fTMVACutValueGraph;
    }

    cout << "reading MVA cut values from " << iFileName << endl;
    cout << "\tfound graphs for " << fTMVACutValueGraph.size() << " zenith angle bin(s)" << endl;
    if( iCutGraphSmoothing > 0 )
    {
         cout << "\tapply cut graph Gaussian smoothing with dE = " << iCutGraphSmoothing;
         cout << " ( max " << iCutGraphSmoothingMax;
         cout << " (constant cut energy " << iCutGraphConstantCutEnergy_TeV << " TeV";
         cout << ", " << fTMVACutValueGraph.back()->GetN() << " points in energy)" << endl;
         if( iSmoothTMVAGraph )
         { 
             cout << "\tsmooth mva values" << endl;
         }
         else
         { 
             cout << "\tsmooth signal efficiencies" << endl;
         }
    }

    return fTMVACutValueGraph;
}

/*
 * fill mva cut value from smooth signal efficiency graph
 * (reverse of TGraphAsymmErrors* VTMVAEvaluator::fillSmoothedEfficencyGraph() )
 */
TGraphAsymmErrors* VTMVAEvaluator::fillSmoothedMVACutGraph(
                             TGraphAsymmErrors* iG,
                             unsigned int iZe )
{
    if( !iG )
    {
        return 0;
    }
    string iName = "TMVACutValue_ze";
    TGraphAsymmErrors *iSmoothed = new TGraphAsymmErrors( 1 );
    iSmoothed->SetLineColor( iZe+1 );
    iSmoothed->SetMarkerColor( iZe+1 );
    ostringstream iGName;
    iGName << iName << iZe << "_smoothed";
    iSmoothed->SetName( iGName.str().c_str() );
    iSmoothed->SetTitle( iGName.str().c_str() );
    double x, y;
    // loop over energies (fine bins)
    for( int i = 0; i < iG->GetN(); i++ )
    {
        iG->GetPoint( i, x, y );
        vector< double > iW = getDataBinWeights( x, iZe );
        double iMVA = 0.;
        double iMVA_tot = 0.;
        for( unsigned int w = 0; w < iW.size(); w++ )
        {
            // ignoring small weights
            if( iW[w] > 0.001 && w < fTMVAData.size() )
            {
                 if( fTMVAData[w]->hSignalEfficiency )
                 {
                     iMVA += iW[w] * 
                           fTMVAData[w]->hSignalEfficiency->GetBinCenter( 
                                   fTMVAData[w]->hSignalEfficiency->FindLastBinAbove( y ) );
                     iMVA_tot += iW[w];
                 }
            }
       }
       // weighted MVA cut
       if( iMVA_tot > 0. )
       {
            iSmoothed->SetPoint( i, x, iMVA / iMVA_tot );
       }
    }

    return iSmoothed;
}


/*
 * fill signal/background efficencies for a smooth mva graph
 * (reverse of TGraphAsymmErrors* VTMVAEvaluator::fillSmoothedMVACutGraph() )
 *
 */
TGraphAsymmErrors* VTMVAEvaluator::fillSmoothedEfficencyGraph(
                             TGraphAsymmErrors* iG,
                             unsigned int iZe,
                             bool iSignalEfficiency )
{
    if( !iG )
    {
        return 0;
    }
    string iName = "SignalEfficiency_ze";
    if( !iSignalEfficiency )
    {
         iName = "BackgroundEfficiency_ze";
    }
    TGraphAsymmErrors *iSmoothed = new TGraphAsymmErrors( 1 );
    iSmoothed->SetLineColor( iZe+1 );
    iSmoothed->SetMarkerColor( iZe+1 );
    ostringstream iGName;
    iGName << iName << iZe << "_smoothed";
    iSmoothed->SetName( iGName.str().c_str() );
    iSmoothed->SetTitle( iGName.str().c_str() );
    double x, y;
    int z = 0;
    for( int i = 0; i < iG->GetN(); i++ )
    {
        iG->GetPoint( i, x, y );
        vector< double > iW = getDataBinWeights( x, iZe );
        double iEff = 0.;
        double iEff_tot = 0.; 
        for( unsigned int w = 0; w < iW.size(); w++ )
        {
            if( iW[w] > 0.001 && w < fTMVAData.size() )
            {
                 if( iSignalEfficiency &&
                 fTMVAData[w]->hSignalEfficiency )
                 {
                     iEff += iW[w] *
                          fTMVAData[w]->hSignalEfficiency->GetBinContent(
                          fTMVAData[w]->hSignalEfficiency->GetXaxis()->FindBin( y ) );
                     iEff_tot += iW[w];
                 }
                 else if( !iSignalEfficiency &&
                     fTMVAData[w]->hBackgroundEfficiency )
                 {
                     iEff = iW[w] *
                          fTMVAData[w]->hBackgroundEfficiency->GetBinContent(
                          fTMVAData[w]->hBackgroundEfficiency->GetXaxis()->FindBin( y ) );
                     iEff_tot += iW[w];
                 }
           }
        }
        if( iEff_tot > 0. )
        {
            iSmoothed->SetPoint( z, x, iEff / iEff_tot );
            z++;
        }
    }

    return iSmoothed;
}

/*
 * smooth a graph roughly with bin width and energy resolution
 *
 */
TGraphAsymmErrors* VTMVAEvaluator::smoothMVAGraph( TGraphAsymmErrors* iG,
                                                   double iCutGraphSmoothing,
                                                   double iCutGraphSmoothingMax,
                                                   double iCutGraphConstantCutEnergy_TeV,
                                                   string iName, unsigned int iZe )
{
     if( !iG )
     {
          return 0;
     }

     double x = 0.;
     double y = 0.;
     iG->GetPoint( 0, x, y );
     double e_min = x - iG->GetErrorXlow( 0 );
     iG->GetPoint( iG->GetN()-1, x, y );
     double e_max = x + iG->GetErrorXhigh( iG->GetN()-1 );
     TGraphAsymmErrors *iSmoothedGaussian = new TGraphAsymmErrors( 1 );
     ostringstream iGName;
     iGName << iName << "_smoothed";
     iSmoothedGaussian->SetName( iGName.str().c_str() );
     iSmoothedGaussian->SetTitle( iGName.str().c_str() );

     int z = 0;
     double e_log10_G = 0.;
     unsigned int nToyMC_Ebins = 100;
     // point in energy
     // (divide energy range into 1000 points)
     for( unsigned int i = 0; i < nToyMC_Ebins; i++ )
     {
          double e_log10 = e_min + i * (e_max - e_min)/(double)nToyMC_Ebins;
          double iM = 0.;
          double iN = 0.;
          // Gaussian smoothing using
          // 10000 iterations per energy point
          for( unsigned int j = 0; j < 10000; j++ )
          {
              e_log10_G = e_log10 + gRandom->Gaus( 0., iCutGraphSmoothing );
              // make sure that values are in the given energy range
              // (extrapolation might otherwise result in very wrong
              // values)
              if( e_log10_G > e_min && e_log10_G < e_max )
              {
                  vector< double > iW = getDataBinWeights( e_log10_G, iZe );
                  double iMVA = 0.;
                  double iMVA_tot = 0.;
                  for( unsigned int w = 0; w < iW.size(); w++ )
                  {
                      // ignoring small weights
                      if( iW[w] > 0.001 && w < fTMVAData.size() )
                      {
                          iG->GetPoint( fTMVAData[w]->fEnergyCut_bin, x, y );
                          iMVA += iW[w] * y;
                          iMVA_tot += iW[w];
                      }
                  }
                  // accept only cut values below a maximum
                  // (optimisation might give in the threshold region
                  //  very high values)
                  if( iMVA_tot > 0. && iMVA / iMVA_tot < iCutGraphSmoothingMax )
                  {
                       iM +=  iMVA / iMVA_tot;
                       iN++;
                  }
              }
          }
          if( iN > 0. )
          {
             // apply constant cut value above this energy
             if( e_log10 > log10( iCutGraphConstantCutEnergy_TeV ) )
             {
                 iSmoothedGaussian->SetPoint( z, e_log10, 
                                   iSmoothedGaussian->Eval( e_min + (i-1) * (e_max - e_min)/(double)nToyMC_Ebins ) );
             }
             else
             {
                 iSmoothedGaussian->SetPoint( z, e_log10, iM / iN );
             }
             z++;
          }
     }
     return iSmoothedGaussian;
}

/*
 * smooth a graph roughly with bin width and energy resolution
 *
 */
vector< TGraphAsymmErrors* > VTMVAEvaluator::smoothSignalEfficiencyMVAGraph( TGraphAsymmErrors* iG,
                                                   double iCutGraphSmoothing,
                                                   double iCutGraphSmoothingMax,
                                                   double iCutGraphConstantCutEnergy_TeV,
                                                   string iName, unsigned int iZe )
{
     vector< TGraphAsymmErrors* > iRG;
     if( !iG )
     {
          iRG.push_back( 0 );
          iRG.push_back( 0 );
          return iRG;
     }

     double x = 0.;
     double y = 0.;
     iG->GetPoint( 0, x, y );
     double e_min = x - iG->GetErrorXlow( 0 );
     iG->GetPoint( iG->GetN()-1, x, y );
     double e_max = x + iG->GetErrorXhigh( iG->GetN()-1 );

     // first element: smooth signal efficiency
     iRG.push_back( new TGraphAsymmErrors( 1 ) );
     ostringstream iGName;
     iGName << iName << "_smoothed";
     iRG.back()->SetName( iGName.str().c_str() );
     iRG.back()->SetTitle( iGName.str().c_str() );
     // second element: smoothed mva efficiency
     iRG.push_back( new TGraphAsymmErrors( 1 ) );
     ostringstream iHName;
     iHName << "TMVACutValue_ze" << iZe << "_smoothed";
     iRG.back()->SetName( iHName.str().c_str() );
     iRG.back()->SetTitle( iHName.str().c_str() );
     // second element: smoothed mva efficiency

     int z = 0;
     double e_log10_G = 0.;
     unsigned int nToyMC_Ebins = 100;
     unsigned int nToyMC = 10000;
     if( iCutGraphSmoothing < 1.e-7 )
     {
         nToyMC = 1;

     }
     // point in energy
     // (divide energy range into 1000 points)
     for( unsigned int i = 0; i < nToyMC_Ebins; i++ )
     {
          double e_log10 = e_min + i * (e_max - e_min)/(double)nToyMC_Ebins;
          double iS = 0.;
          double iM = 0.;
          double iN = 0.;
          // Gaussian smoothing using
          // 10000 iterations per energy point
          for( unsigned int j = 0; j < nToyMC; j++ )
          {
              e_log10_G = e_log10 + gRandom->Gaus( 0., iCutGraphSmoothing );
              // make sure that values are in the given energy range
              // (extrapolation might otherwise result in very wrong
              // values)
              if( e_log10_G > e_min && e_log10_G < e_max )
              {
                  vector< double > iW = getDataBinWeights( e_log10_G, iZe );
                  double iEffS = 0.;
                  double iEffS_tot = 0.;
                  double iMVA = 0.;
                  double iMVA_tot = 0.;
                  for( unsigned int w = 0; w < iW.size(); w++ )
                  {
                      // ignoring small weights
                      if( w < fTMVAData.size() && iW[w] > 0.001 )
                      {
                          y = iG->Eval( e_log10_G );
                          iEffS += iW[w] * y;
                          iEffS_tot += iW[w];
                          // get MVA value for this signal efficiency
                          if( fTMVAData[w]->hSignalEfficiency )
                          {
                                 iMVA += iW[w] *
                                 fTMVAData[w]->hSignalEfficiency->GetBinCenter( 
                                               fTMVAData[w]->hSignalEfficiency->FindLastBinAbove( y ) ); 
                                 iMVA_tot += iW[w];
                          }
                      }
                  }
                  // accept only cut values below a maximum
                  // (optimisation might give in the threshold region
                  //  very high values)
                  if( iEffS_tot > 0. 
                  && iEffS / iEffS_tot < iCutGraphSmoothingMax
                  && iMVA_tot > 0. )
                  {
                       iS +=  iEffS / iEffS_tot;
                       iM +=  iMVA / iMVA_tot;
                       iN++;
                  }
              }
          }
          if( iN > 0. )
          {
             // apply constant cut value above this energy
             if( e_log10 > log10( iCutGraphConstantCutEnergy_TeV ) )
             {
                 iRG[0]->SetPoint( z, e_log10, 
                                   iRG[0]->Eval( e_min + (i-1) * (e_max - e_min)/(double)nToyMC_Ebins) );
                 iRG[1]->SetPoint( z, e_log10, 
                                   iRG[1]->Eval( e_min + (i-1) * (e_max - e_min)/(double)nToyMC_Ebins) );
             }
             else
             {
                 iRG[0]->SetPoint( z, e_log10, iS / iN );
                 iRG[1]->SetPoint( z, e_log10, iM / iN );
             }
             z++;
          }
     }
     return iRG;
}


/*
 * tmva cuts are given by maps
 * key = n-digit integer:
 *          % 10 = zenith angle bin
 *          / 10 = energy bin
 * values = optimised TMVA cut value
 *
 */
void VTMVAEvaluator::setTMVACutValue( map< unsigned int, double > iMVA )
{
    fTMVACutValueMap = iMVA;
    
    fTMVACutValueNoVec = 1.;
}

void VTMVAEvaluator::setTMVACutValue( double iE )
{
    fTMVACutValueNoVec = iE;
    fTMVACutValueFile = 0;
}

void VTMVAEvaluator::printSourceStrength_CU()
{
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            cout << "E [" << showpoint << setprecision( 3 );
            cout << fTMVAData[i]->fEnergyCut_Log10TeV_min << ",";
            cout << fTMVAData[i]->fEnergyCut_Log10TeV_max << "] TeV";
            cout << " (bin " << i << "):\t ";
            cout << fTMVAData[i]->fSourceStrengthAtOptimum_CU << " CU " << endl;
        }
    }
}

void VTMVAEvaluator::printAngularContainmentRadius()
{
    if( fTMVAData.size() == 0 )
    {
        return;
    }
    // make sure that there are any valid values
    bool iPrint = false;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
             if( fTMVAData[i]->fAngularContainmentRadius >= 0. )
             {
                 iPrint = true;
                 break;
             }
         }
    }
    if( !iPrint )
    {
        return;
    }
    
    cout << endl;
    cout << "======================= VTMVAEvaluator: energy dependent optimal containment radius cut =======================" << endl;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            cout << "E [" << showpoint << setprecision( 3 );
            cout << showpos << fTMVAData[i]->fEnergyCut_Log10TeV_min << "," << fTMVAData[i]->fEnergyCut_Log10TeV_max << "] TeV";
            cout << noshowpos <<  ", Ze [" << fTMVAData[i]->fZenithCut_min << "," << fTMVAData[i]->fZenithCut_max << "] deg";
            cout << " (bin " << i << "):\t ";
            cout << fTMVAData[i]->fAngularContainmentRadius << " [deg], ";
            cout << "T^2: " << fTMVAData[i]->fAngularContainmentRadius* fTMVAData[i]->fAngularContainmentRadius << " [deg^2], ";
            cout << fTMVAData[i]->fAngularContainmentFraction << "%";
            cout << endl;
        }
    }
    cout << noshowpoint << endl;
}

/* 
 * write a TGraph2D with MVA cut values as function
 * of log10(energy) and zenith angle
 */
bool VTMVAEvaluator::writeOptimizedMVACutValues( string iRootFile )
{
    TFile *i_F = new TFile( iRootFile.c_str(), "RECREATE" );
    if( i_F->IsZombie() )
    {
       cout << "writeOptimizedMVACutValues: ";
       cout << " error opening file " << iRootFile << endl;
       return false;
    }

    TGraph2DErrors *iG = new TGraph2DErrors( 1 );
    iG->SetTitle( "MVA cut value" );
    iG->SetMarkerStyle( 20 );

    unsigned int z = 0;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] && fTMVAData[i]->fTMVACutValue > -90. )
        {
              iG->SetPoint( z,
                    fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV,
                    0.5*(fTMVAData[i]->fZenithCut_max+fTMVAData[i]->fZenithCut_min),
                    fTMVAData[i]->fTMVACutValue );
              iG->SetPointError( z, 0., 0., 0. );
              z++;
        }
    }

    // write 1D graphs
    char hname[200];
    vector< TGraphAsymmErrors* > ivG;
    vector< TGraphAsymmErrors* > ivS;
    vector< TGraphAsymmErrors* > ivB;
    for( unsigned int z = fWeightFileIndex_Zmin; z <= fWeightFileIndex_Zmax; z++ )
    {
        sprintf( hname, "MVA cut value (zenith bin %u)", z );
        ivG.push_back( new TGraphAsymmErrors(1) );
        ivG.back()->SetTitle( hname );
        ivG.back()->SetLineColor( z+1 );
        ivG.back()->SetMarkerColor( z+1 );
        sprintf( hname, "signal efficiency (zenith bin %u)", z );
        ivS.push_back( new TGraphAsymmErrors(1) );
        ivS.back()->SetTitle( hname );
        ivS.back()->SetLineColor( z+1 );
        ivS.back()->SetMarkerColor( z+1 );
        sprintf( hname, "background efficiency (zenith bin %u)", z );
        ivB.push_back( new TGraphAsymmErrors(1) );
        ivB.back()->SetTitle( hname );
        ivB.back()->SetLineColor( z+1 );
        ivB.back()->SetMarkerColor( z+1 );

        unsigned int zz = 0;
        unsigned int zs = 0;
        unsigned int zb = 0;
        for( unsigned int i = 0; i < fTMVAData.size(); i++ )
        {
             if( fTMVAData[i] && fTMVAData[i]->fTMVACutValue > -90.
                && fTMVAData[i]->fZenithCut_bin == z )
             {
                  ivG.back()->SetPoint( zz, 
                       fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV,
                       fTMVAData[i]->fTMVACutValue );
                  ivG.back()->SetPointEXlow( zz,
                      fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                  ivG.back()->SetPointEXhigh( zz,
                      fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                  zz++;
                  if( fTMVAData[i]->fSignalEfficiency > 0. )
                  {
                      ivS.back()->SetPoint( zs, 
                           fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV,
                           fTMVAData[i]->fSignalEfficiency );
                      ivS.back()->SetPointEXlow( zs,
                          fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                      ivS.back()->SetPointEXhigh( zs,
                          fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                      zs++;
                  }
                  if( fTMVAData[i]->fBackgroundEfficiency > 0. )
                  {
                      ivB.back()->SetPoint( zb, 
                       fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV,
                       fTMVAData[i]->fBackgroundEfficiency );
                      ivB.back()->SetPointEXlow( zb,
                          fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV - fTMVAData[i]->fEnergyCut_Log10TeV_min );
                      ivB.back()->SetPointEXhigh( zb,
                          fTMVAData[i]->fEnergyCut_Log10TeV_max - fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                      zb++;
                  }
             }
        }
    }
               
    cout << "writing TMVA cut values to " << i_F->GetName() << endl;
    iG->Write( "TMVACutValue" );
    for( unsigned int i = 0; i < ivG.size(); i++ )
    {
        if( ivG[i] && ivS[i] && ivB[i] )
        {
             sprintf( hname, "TMVACutValue_ze%u", i );
             ivG[i]->Write( hname );
             sprintf( hname, "SignalEfficiency_ze%u", i );
             ivS[i]->Write( hname );
             sprintf( hname, "BackgroundEfficiency_ze%u", i );
             ivB[i]->Write( hname );
        }
    }
    i_F->Close();

    return true;
}


/*
 * print MVA cut values after optimization to be used
 * in gamma/hadron cuts files
 */
void VTMVAEvaluator::printOptimizedMVACutValues( string iEpoch )
{
    cout << "Printing Optimised cuts for gamma/hadron cut values" << endl;
    cout << "\t this is only correct for an energyStepSize of -1" << endl;
    cout << "* TMVA_MVACut " << iEpoch;
    unsigned int iCounter = 0;
    char hname[20];
    for( unsigned int e = fWeightFileIndex_Emin; e <= fWeightFileIndex_Emax; e++ )
    {
        for( unsigned int z = fWeightFileIndex_Zmin; z <= fWeightFileIndex_Zmax; z++ )
        {
            sprintf( hname, "%u%u", e, z );
            if( iCounter < fTMVAData.size() )
            {
                cout << " " << hname;
                cout << " " << fTMVAData[iCounter]->fTMVACutValue;
            }
            iCounter++;
        }
    }
    cout << endl;
    cout << "(first digit: energy bin, second digit: zenith bin)" << endl;
}
   
/* 
 * print a table with signal efficiencies
 *
 */
void VTMVAEvaluator::printSignalEfficiency()
{
    if( fTMVAData.size() == 0 )
    {
        return;
    }
    
    cout << endl;
    cout << "======================= VTMVAEvaluator: signal (background) efficiency =======================" << endl;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            cout << "E [" << showpoint << setprecision( 3 ) <<  showpos << fTMVAData[i]->fEnergyCut_Log10TeV_min;
            cout << "," << fTMVAData[i]->fEnergyCut_Log10TeV_max << "] TeV" << noshowpos;
            cout << ", Ze [" << fTMVAData[i]->fZenithCut_min << "," << fTMVAData[i]->fZenithCut_max << "] deg";
            cout << " (bin " << i << "):\t ";
            // cut graph
             if( fTMVACutValueGraph.size() > 0
             && fTMVAData[i]->fZenithCut_bin < fTMVACutValueGraph.size()
             && fTMVACutValueGraph[fTMVAData[i]->fZenithCut_bin] )
             {
                   // energy bin not needed for getSignalEfficiency
                   cout << getSignalEfficiency( 0, 
                                fTMVAData[i]->fEnergyCut_Log10TeV_min, fTMVAData[i]->fEnergyCut_Log10TeV_max,
                                fTMVAData[i]->fZenithCut_bin,
                                fTMVAData[i]->fZenithCut_min, fTMVAData[i]->fZenithCut_max,
                                0, 0.1 );
                   cout << "\t MVACut: ";
                   cout << fTMVACutValueGraph[fTMVAData[i]->fZenithCut_bin]->Eval(
                               0.5*(fTMVAData[i]->fEnergyCut_Log10TeV_min+fTMVAData[i]->fEnergyCut_Log10TeV_max));
             }
            // binned cut values
            else
            {
                cout << fTMVAData[i]->fSignalEfficiency;
                if( fTMVAData[i]->fBackgroundEfficiency > 0. )
                {
                    cout << "\t(" << fTMVAData[i]->fBackgroundEfficiency << ")";
                }
                cout << "\t MVACut: " << fTMVAData[i]->fTMVACutValue;
                if( fParticleNumberFileName.size() > 0 )
                {
                    if( fTMVAData[i]->fTMVAOptimumCutValueFound )
                    {
                        cout << " (optimum reached for " << fTMVAData[i]->fSourceStrengthAtOptimum_CU << " CU)";
                    }
                    else
                    {
                        cout << " (no optimum reached (" << fTMVAData[i]->fSourceStrengthAtOptimum_CU << " CU)";
                    }
                }
            }
        }
        cout << endl;
    }
    cout << noshowpoint << endl;
}

/*

    calculate the optimal signal to noise ratio for a given particle number spectrum

    this routine is possibly too complicated

    - main problem is how to deal with low statistics bins

*/

bool VTMVAEvaluator::optimizeSensitivity( unsigned int iDataBin, 
                                          string iOptimizationType, 
                                          string iInstrumentEpoch )
{
    // valid data bin
    if( iDataBin >= fTMVAData.size() || !fTMVAData[iDataBin] )
    {
        return false;
    }
    
    // print some info on optimization parameters to screen
    printSensitivityOptimizationParameters();
    
    //////////////////////////////////////////////////////
    // read file with  NOn and Noff graphs
    // (contains signal and background rate vs energy)
    // (created from effective areas with quality cuts applied only,
    TFile* iPN = new TFile( fParticleNumberFileName.c_str() );
    if( iPN->IsZombie() )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity error:" << endl;
        cout << " cannot read particle number file " << fParticleNumberFileName << endl;
        if( fParticleNumberFileName.size() == 0 )
        {
            cout << "VTMVAEvaluator::optimizeSensitivity error: no particle number file given" << endl;
        }
        return false;
    }
    cout << "VTVMAEvaluator::optimizeSensitivity: reading: " << fParticleNumberFileName << endl;
    // angular containment histogram
    // (optional, only needed when containment radius is also optimized)
    TH2D* iHAngContainment = ( TH2D* )iPN->Get( "AngResCumulative" );
    if( iHAngContainment )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: found angular containment histogram (";
        cout << fTMVAngularContainmentRadiusMax << "%)";
        cout << endl;
    }
    // (no error message if angular containment histogram is not found
    //  simply means that theta2 cut is not optimised)
    
    ///////////////////////////////////////////////////////////////////////////////
    // get number of events (after quality cuts) at this energy from on/off graphs
    double Non = 0.;
    double Nof = 0.;
    double Ndif = 0.;
        
    // bin width in energy (linear axis)
    double i_dE = TMath::Power( 10., fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max )
                  - TMath::Power( 10., fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min );
    // in CTA, a log axis is used!!!
    i_dE = TMath::Abs( fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max - fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min );
    
    //////////////////////////////////////////////////////
    // Interpolate 2D graphs from particle rate files
    // less robust method which get excess and background counts
    // by interpolating between graphs
    // (favored for CTA analysis)
    if( iOptimizationType == "UseInterpolatedCounts" )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: UseInterpolatedCounts (conversion rate ";
        cout << fParticleNumberFile_Conversion_Rate_to_seconds << ")" << endl;
        // get the NOn (signal + background) and NOff (background) graphs
        // Interpolation between zenith angles happens on secant axis
        TGraph* i_on = readInterpolatedCountsFromFile( iPN, 1. / cos( fTMVAData[iDataBin]->fZenithCut_min * TMath::DegToRad() ),
                       1. / cos( fTMVAData[iDataBin]->fZenithCut_max * TMath::DegToRad() ),
                       true );
        TGraph* i_of = readInterpolatedCountsFromFile( iPN, 1. / cos( fTMVAData[iDataBin]->fZenithCut_min * TMath::DegToRad() ),
                       1. / cos( fTMVAData[iDataBin]->fZenithCut_max * TMath::DegToRad() ),
                       false );
        if( !i_on || !i_of )
        {
            cout << "VTVMAEvaluator::optimizeSensitivity: error," << endl;
            cout << " cannot read graphs from particle number file " << endl;
            cout << i_on << "\t" << i_of << endl;
            return false;
        }
        //////////////// CTA ONLY /////////////////////////////////////////////////
        // mean energy calculation below requires equal binning for the count graphs
        if( iInstrumentEpoch == "noepoch" || iInstrumentEpoch == "CTA" )
        {
            double x = 0.;
            double p = 0.;
            // get mean energy of the considered bins
            // interval [fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min,fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max]
            // make sure that energy is not lower or higher then minimum/maximum bins in the rate graphs
            for( int ii = 0; ii < i_on->GetN(); ii++ )
            {
                i_on->GetPoint( ii, x, p );
                if( p > 0. && ( x + i_on->GetErrorX( ii ) <= fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max
                                || x <= fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max ) )
                {
                    if( fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV < x )
                    {
                        if( x + i_on->GetErrorX( ii ) <= fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max )
                        {
                            fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV = x + 0.97 * i_on->GetErrorX( ii );
                        }
                        else
                        {
                            fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV = x;
                        }
                    }
                    break;
                }
            }
            ///////////
            // make sure that selected energy is not beyond the valid range of the graph
            
            // get the value of the energy, zenith and particle rate for the last index of the array
            i_on->GetPoint( i_on->GetN() - 1, x, p );
            
            // energy is beyond - set it to 0.8*last value
            if( fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV > x )
            {
                cout << "LAST POINT " << fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV << "\t" << x << endl;
                fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV = TMath::Log10( TMath::Power( 10., x ) * 0.8 );
            }
        }
        //////////////// end of CTA ONLY /////////////////////////////////////////////////
        
        ///////////////////////////////////////////////////////////////////////////////
        // Interpolate between the values of the TGraph2D
        //
        // Convert the observing time in seconds as the particle rate is given in 1/seconds
        // Get the value of the middle of the energy and zenith angle bin
        Ndif = i_on->Eval( fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV );
        Nof = i_of->Eval( fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV );
    }
    //////////////////////////////////////////////////////////////////
    // robust method which get excess and background counts
    // by averaging over graphs
    // (favored for VTS analysis)
    else if( iOptimizationType == "UseAveragedCounts" )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: UseAveragedCounts (conversion rate ";
        cout << fParticleNumberFile_Conversion_Rate_to_seconds << ")" << endl;
        // read signal=excess counts
        Ndif = readAverageCountsFromFile( iPN, fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min,
                                          fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max,
                                          1. / cos( fTMVAData[iDataBin]->fZenithCut_min * TMath::DegToRad() ),
                                          1. / cos( fTMVAData[iDataBin]->fZenithCut_max * TMath::DegToRad() ),
                                          true );
        // read background counts
        Nof = readAverageCountsFromFile( iPN, fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min,
                                         fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max,
                                         1. / cos( fTMVAData[iDataBin]->fZenithCut_min * TMath::DegToRad() ),
                                         1. / cos( fTMVAData[iDataBin]->fZenithCut_max * TMath::DegToRad() ),
                                         false );
    }
    else
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: error," << endl;
        cout << " unknown optimization type" << endl;
        return false;
    }
    if( Nof < 0. )
    {
        Nof = 0.;
    }
    // correct normalisation and times dE
    Ndif *= fOptimizationObservingTime_h * fParticleNumberFile_Conversion_Rate_to_seconds * i_dE;
    Nof  *= fOptimizationObservingTime_h * fParticleNumberFile_Conversion_Rate_to_seconds * i_dE;
    Non = Ndif + Nof;
    
    cout << "VTVMAEvaluator::optimizeSensitivity: event numbers before optimization: ";
    cout << " non = " << Non;
    cout << " noff = " << Nof;
    cout << " ndif = " << Ndif << " (1 CU)" << endl;
    cout << "VTVMAEvaluator::optimizeSensitivity: data bin: ";
    cout << iDataBin;
    cout << ",  weighted mean energy " << TMath::Power( 10., fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV );
    cout << " [TeV], ";
    cout << fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV;
    cout << " [" << fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min << ", ";
    cout << fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max << "]";
    cout << endl;
    
    ///////////////////////////////////////////////////////////////////
    // get signal and background efficiency histograms from TMVA files
    
    TFile* iTMVAFile = new TFile( fTMVAData[iDataBin]->fTMVAFileName.c_str() );
    if( iTMVAFile->IsZombie() )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: error:" << endl;
        cout << " cannot read TMVA file " << fTMVAData[iDataBin]->fTMVAFileName;
        cout << " (bin " << iDataBin << ")" << endl;
        return false;
    }
    // get signal and background efficiency histograms
    // TH1D* effS = getEfficiencyHistogram( "effS", iTMVAFile, fTMVAData[iDataBin]->fTMVAMethodTag_2 );
    // TH1D* effB = getEfficiencyHistogram( "effB", iTMVAFile, fTMVAData[iDataBin]->fTMVAMethodTag_2 );
    TH1D* effS = fTMVAData[iDataBin]->hSignalEfficiency;
    TH1D* effB = fTMVAData[iDataBin]->hBackgroundEfficiency;
    if( !effS || !effB )
    {
        cout << "VTVMAEvaluator::optimizeSensitivity: error:" << endl;
        cout << " cannot find signal and/or background efficiency histogram(s)" << endl;
        cout << effS << "\t" << effB << endl;
        return false;
    }
    
    cout << "VTVMAEvaluator::optimizeSensitivity: optimization parameters: ";
    cout << "maximum signal efficiency is " << fOptimizationFixedSignalEfficiency;
    cout << " minimum source strength is " << fOptimizationMinSourceStrength;
    cout << " (alpha: " << fOptimizationBackgroundAlpha << ")" << endl;
    
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // optimization starts here
    //////////////////////////////////////////////////////////////////////////
    double i_Signal_to_sqrtNoise = 0.;
    double i_AngularContainmentRadius = 0.;
    double i_AngularContainmentFraction = 0.;
    
    double i_TMVACutValue_AtMaximum = -99.;
    double i_SourceStrength_atMaximum = 0.;
    double i_SignalEfficiency_AtMaximum = -99.;
    double i_BackgroundEfficiency_AtMaximum = -99.;
    double i_Signal_to_sqrtNoise_atMaximum = 0.;
    double i_AngularContainmentRadiusAtMaximum = 0.;
    double i_AngularContainmentFractionAtMaximum = 0.;
    double iSourceStrength = 0.;
    
    TGraph* iGSignal_to_sqrtNoise = 0;
    TGraph* iGSignalEvents        = 0;
    TGraph* iGBackgroundEvents    = 0;
    TGraph* iGSignal_to_sqrtNoise_Smooth = 0;
    TGraph* iGOpt_AngularContainmentRadius = 0;
    TGraph* iGOpt_AngularContainmentFraction = 0;
    
    //////////////////////////////////////////////////////
    // loop over different source strengths (in Crab Units)
    // (hardwired: start at 0.001 CU to 30 CU)
    unsigned int iSourceStrengthStepSizeN = 
              ( unsigned int )( ( log10( fOptimizationMaxSourceStrength ) - log10( fOptimizationMinSourceStrength ) ) / 0.005 );
    cout << "VTVMAEvaluator::optimizeSensitivity: source strength steps: " << iSourceStrengthStepSizeN << endl;
    cout << "VTVMAEvaluator::optimizeSensitivity: range for source strength: (";
    cout << fOptimizationMinSourceStrength << ", " << fOptimizationMaxSourceStrength << ") CU" << endl;
    for( unsigned int s = 0; s < iSourceStrengthStepSizeN; s++ )
    {
        iSourceStrength = log10( fOptimizationMinSourceStrength ) + s * 0.005;
        iSourceStrength = TMath::Power( 10., iSourceStrength );
        
        // source events
        Ndif = ( Non - Nof ) * iSourceStrength;
        
        // first quick pass to see if there is a change of reaching the required fOptimizationSourceSignificance
        // (needed to speed up the calculation)
        // (ignore any detail, no optimization of angular cut)
        bool bPassed = false;
        for( int i = 1; i < effS->GetNbinsX(); i++ )
        {
            if( effB->GetBinContent( i ) > 0. && Nof > 0. )
            {
                if( fOptimizationBackgroundAlpha > 0. )
                {
                    i_Signal_to_sqrtNoise = VStatistics::calcSignificance( effS->GetBinContent( i ) * Ndif + effB->GetBinContent( i ) * Nof,
                                            effB->GetBinContent( i ) * Nof / fOptimizationBackgroundAlpha,
                                            fOptimizationBackgroundAlpha );
                    // check significance criteria
                    if( i_Signal_to_sqrtNoise > fOptimizationSourceSignificance )
                    {
                        bPassed = true;
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        // no chance to pass significance criteria -> continue to next energy bin
        if( !bPassed )
        {
            continue;
        }
        
        //////////////////////////////////////////////////////
        // now loop over signal and background efficiency levels
        i_Signal_to_sqrtNoise = 0.;
        i_AngularContainmentRadius = 0.;
        i_AngularContainmentFraction = 0.;
        
        i_TMVACutValue_AtMaximum = -99.;
        i_SourceStrength_atMaximum = 0.;
        i_SignalEfficiency_AtMaximum = -99.;
        i_BackgroundEfficiency_AtMaximum = -99.;
        i_Signal_to_sqrtNoise_atMaximum = 0.;
        i_AngularContainmentRadiusAtMaximum = 0.;
        i_AngularContainmentFractionAtMaximum = 0.;
        
        iGSignal_to_sqrtNoise = new TGraph( 1 );
        iGSignalEvents        = new TGraph( 1 );
        iGBackgroundEvents    = new TGraph( 1 );
        iGOpt_AngularContainmentRadius = new TGraph( 1 );
        iGOpt_AngularContainmentFraction = new TGraph( 1 );
        
        int z = 0;
        int z_SB = 0;
        // loop over all signal efficiency bins
        for( int i = 1; i < effS->GetNbinsX(); i++ )
        {
            if( effB->GetBinContent( i ) > 0. && Nof > 0. )
            {
                if( fOptimizationBackgroundAlpha > 0. )
                {
                    // optimize angular containment radius (theta2)
                    if( iHAngContainment && fTMVA_OptimizeAngularContainment )
                    {
                        getOptimalAngularContainmentRadius( effS->GetBinContent( i ), effB->GetBinContent( i ), Ndif, Nof,
                                                            iHAngContainment, fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV,
                                                            i_Signal_to_sqrtNoise, i_AngularContainmentRadius, i_AngularContainmentFraction );
                    }
                    // optimize signal/sqrt(noise)
                    else
                    {
                        i_Signal_to_sqrtNoise = VStatistics::calcSignificance( effS->GetBinContent( i ) * Ndif + effB->GetBinContent( i ) * Nof,
                                                effB->GetBinContent( i ) * Nof / fOptimizationBackgroundAlpha,
                                                fOptimizationBackgroundAlpha );
                        i_AngularContainmentRadius = VHistogramUtilities::interpolateTH2D( iHAngContainment,
                                                     fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV,
                                                     fTMVAngularContainmentRadiusMax );
                        i_AngularContainmentFraction = fTMVAngularContainmentRadiusMax;
                    }
                }
                else
                {
                    i_Signal_to_sqrtNoise = 0.;
                    i_AngularContainmentRadius = 0.;
                    i_AngularContainmentFraction = fTMVAngularContainmentRadiusMax;
                }
                if( fDebug )
                {
                    cout << "___________________________________________________________" << endl;
                    cout << i << "\t" << Non << "\t" << effS->GetBinContent( i )  << "\t";
                    cout << Nof << "\t" << effB->GetBinContent( i ) << "\t";
                    cout << Ndif << endl;
                    cout << "\t" << effS->GetBinContent( i ) * Ndif;
                    cout << "\t" << effS->GetBinContent( i ) * Ndif + effB->GetBinContent( i ) * Nof;
                    cout << "\t" << effS->GetBinContent( i ) * Non + effB->GetBinContent( i ) * Nof;
                    cout << "\t" << effB->GetBinContent( i ) * Nof << endl;
                }
                if( effS->GetBinContent( i ) * Ndif > 0. )
                {
                    iGSignalEvents->SetPoint( z_SB, effS->GetBinCenter( i ),  effS->GetBinContent( i ) * Ndif );
                    iGBackgroundEvents->SetPoint( z_SB, effS->GetBinCenter( i ), effB->GetBinContent( i ) * Nof );
                    z_SB++;
                }
                // check that a minimum number of off events is available
                if( effB->GetBinContent( i ) * Nof < fOptimizationMinBackGroundEvents )
                {
                    if( fDebug )
                    {
                        cout << "\t number of background events lower than ";
                        cout << fOptimizationMinBackGroundEvents << ": setting signal/sqrt(noise) to 0; bin " << i << endl;
                    }
                    i_Signal_to_sqrtNoise = 0.;
                }
                // add results to a graph
                if( iGSignal_to_sqrtNoise && i_Signal_to_sqrtNoise > 1.e-2 )
                {
                    iGSignal_to_sqrtNoise->SetPoint( z, effS->GetBinCenter( i ), i_Signal_to_sqrtNoise );
                    iGOpt_AngularContainmentRadius->SetPoint( z, effS->GetBinCenter( i ), i_AngularContainmentRadius );
                    iGOpt_AngularContainmentFraction->SetPoint( z, effS->GetBinCenter( i ), i_AngularContainmentFraction );
                    if( fDebug )
                    {
                        cout << "\t SET " << z << "\t" << effS->GetBinCenter( i ) << "\t" << i_Signal_to_sqrtNoise << "\t";
                        cout << i_AngularContainmentRadius << "\t" << i_AngularContainmentFraction << endl;
                    }
                    z++;
                }
                if( fDebug )
                {
                    cout << "\t z " << z << "\t" << i_Signal_to_sqrtNoise << endl;
                    cout << "___________________________________________________________" << endl;
                }
            }
        } // END loop over all signal efficiency bins
        /////////////////////////
        // determine position of maximum significance
        // fill a histogram from these values, smooth it, and determine position of maximum significance
        double i_xmax = -99.;
        if( iGSignal_to_sqrtNoise )
        {
            TGraphSmooth* iGSmooth = new TGraphSmooth( "s" );
            iGSignal_to_sqrtNoise_Smooth = iGSmooth->SmoothKern( iGSignal_to_sqrtNoise, "normal", 0.05, 100 );
            // find maximum in significance plot
            double i_ymax = -99.;
            double i_xmax_global = -99.;
            double i_ymax_global = -99.;
            for( int i = 0; i < iGSignal_to_sqrtNoise_Smooth->GetN(); i++ )
            {
                iGSignal_to_sqrtNoise_Smooth->GetPoint( i, i_xmax, i_ymax );

                /// check if this point passes all critera:
                // - significance
                // - systematic criterium
                // - min number of background events

                // passed significance criteria
                if( i_ymax >= fOptimizationSourceSignificance )
                {
                    i_SignalEfficiency_AtMaximum     = effS->GetBinContent( effS->FindBin( i_xmax ) );
                    i_BackgroundEfficiency_AtMaximum = effB->GetBinContent( effB->FindBin( i_xmax ) );
                    // systematic cut criterium
                    if( i_BackgroundEfficiency_AtMaximum * Nof > 0 && fOptimizationBackgroundAlpha > 0. )
                    {
                        if( Ndif * i_SignalEfficiency_AtMaximum / ( i_BackgroundEfficiency_AtMaximum * Nof ) >= fMinBackgroundRateRatio_min )
                        {
                             // number of signal events criterium
                             if( Ndif * i_SignalEfficiency_AtMaximum >= fOptimizationMinSignalEvents )
                             {
                                  // check if this is a global maximum
                                  // (otherwise would find first bin above sig requirement)
                                  if( i_ymax > i_ymax_global )
                                  {
                                       i_xmax_global = i_xmax;
                                       i_ymax_global = i_ymax;
                                  }
                             }
                        }
                    }
                }
            }
            i_SignalEfficiency_AtMaximum     = effS->GetBinContent( effS->FindBin( i_xmax_global ) );
            i_BackgroundEfficiency_AtMaximum = effB->GetBinContent( effB->FindBin( i_xmax_global ) );
            i_TMVACutValue_AtMaximum         = i_xmax_global;
            i_Signal_to_sqrtNoise_atMaximum  = i_ymax_global;
            i_SourceStrength_atMaximum       = iSourceStrength;
            i_AngularContainmentRadiusAtMaximum = iGOpt_AngularContainmentRadius->Eval( i_xmax_global );
            i_AngularContainmentFractionAtMaximum = iGOpt_AngularContainmentFraction->Eval( i_xmax_global );
        }
        fTMVAData[iDataBin]->fTMVAOptimumCutValueFound = true;
        /////////////////////////////////
        // check detection criteria
        
        // significance criteria
        bool bPassed_MinimumSignificance = ( i_Signal_to_sqrtNoise_atMaximum >= fOptimizationSourceSignificance );
        // require a minimum number of signal events
        bool bPassed_MinimumSignalEvents = ( Ndif * i_SignalEfficiency_AtMaximum >= fOptimizationMinSignalEvents );
        // require the signal to be larger than a certain fraction of background
        bool bPasses_MinimumSystematicCut = false;
        if( i_BackgroundEfficiency_AtMaximum * Nof > 0 && fOptimizationBackgroundAlpha > 0. )
        {
            bPasses_MinimumSystematicCut = 
                       ( Ndif * i_SignalEfficiency_AtMaximum / ( i_BackgroundEfficiency_AtMaximum * Nof )
                       >= fMinBackgroundRateRatio_min );
        }
        
        if( bPassed_MinimumSignificance && !bPassed_MinimumSignalEvents )
        {
            cout << "\t passed significance but not signal events criterium";
            cout << " (" << iSourceStrength << " CU): ";
            cout << "sig " << i_Signal_to_sqrtNoise_atMaximum;
            cout << ", Ndif " << Ndif * i_SignalEfficiency_AtMaximum << endl;
        }
        if( bPassed_MinimumSignificance && !bPasses_MinimumSystematicCut )
        {
            cout << "\t passed significance but not systematics criterium";
            cout << " (" << iSourceStrength << " CU): ";
            cout << "sig " << i_Signal_to_sqrtNoise_atMaximum;
            cout << ", Ndif " << Ndif * i_SignalEfficiency_AtMaximum;
            cout << ", Noff " << i_BackgroundEfficiency_AtMaximum * Nof;
            if( i_BackgroundEfficiency_AtMaximum * Nof > 1.e-10 )
            {
                cout << ", sig/bck ratio ";
                cout << Ndif * i_SignalEfficiency_AtMaximum / ( i_BackgroundEfficiency_AtMaximum * Nof );
            }
            cout << endl;
        }
        // good! Passed all three requirements --> exit loop over source strengths
        if( bPassed_MinimumSignificance && bPassed_MinimumSignalEvents && bPasses_MinimumSystematicCut )
        {
            break;
        }
        
        // delete graphs
        // (not in last step, keep them there for plotting)
        if( s != iSourceStrengthStepSizeN - 1 )
        {
            if( iGSignal_to_sqrtNoise )
            {
                delete iGSignal_to_sqrtNoise;
            }
            if( iGSignalEvents )
            {
                delete iGSignalEvents;
            }
            if( iGBackgroundEvents )
            {
                delete iGBackgroundEvents;
            }
            if( iGSignal_to_sqrtNoise_Smooth )
            {
                delete iGSignal_to_sqrtNoise_Smooth;
            }
            if( iGOpt_AngularContainmentRadius )
            {
                delete iGOpt_AngularContainmentRadius;
            }
            if( iGOpt_AngularContainmentFraction )
            {
                delete iGOpt_AngularContainmentFraction;
            }
        }
    } // end of loop over source strength
    cout << "VTVMAEvaluator::optimizeSensitivity (finished looping over source strengths)";
    cout << ", last value: " << iSourceStrength << endl;
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // check if 
    // - signal efficiency is above allowed value
    // - source strengh maximum has been reached
    if( i_SignalEfficiency_AtMaximum > fOptimizationFixedSignalEfficiency ||
        ( fOptimizationMaxSourceStrength > 0. &&
            iSourceStrength / fOptimizationMaxSourceStrength > 0.95 ) )
    {
        if( fOptimizationFixedSignalEfficiency > 0.99 )
        {
            i_TMVACutValue_AtMaximum         = effS->GetBinCenter( effS->GetNbinsX() - 1 );
            i_BackgroundEfficiency_AtMaximum = effB->GetBinContent( effS->GetNbinsX() - 1 );
        }
        else
        {
            for( int i = 1; i < effS->GetNbinsX(); i++ )
            {
                if( effS->GetBinContent( i ) < fOptimizationFixedSignalEfficiency )
                {
                    i_TMVACutValue_AtMaximum         = effS->GetBinCenter( i );
                    i_BackgroundEfficiency_AtMaximum = effB->GetBinContent( i );
                    if( iGOpt_AngularContainmentRadius && iGOpt_AngularContainmentFraction )
                    {
                        i_AngularContainmentRadiusAtMaximum = iGOpt_AngularContainmentRadius->Eval( i_TMVACutValue_AtMaximum );
                        i_AngularContainmentFractionAtMaximum = iGOpt_AngularContainmentFraction->Eval( i_TMVACutValue_AtMaximum );
                    }
                    break;
                }
            }
        }
        cout << "VTMVAEvaluator::optimizeSensitivity: found signal efficiency ";
        cout << i_SignalEfficiency_AtMaximum << " above allowed value (";
        cout << fOptimizationFixedSignalEfficiency << ")" << endl;
        i_SignalEfficiency_AtMaximum = fOptimizationFixedSignalEfficiency;
        cout << "VTMVAEvaluator::optimizeSensitivity: setting signal efficiency to ";
        cout << fOptimizationFixedSignalEfficiency << endl;
    }
    // regular case: 
    //   - maximum found and reasonable
    //   - signal efficiency in allowed range
    else if( i_SourceStrength_atMaximum > 0. )
    {
        cout << "VTMVAEvaluator::optimizeSensitivity: signal efficiency at maximum (";
        cout << i_SourceStrength_atMaximum << " CU) is ";
        cout << i_SignalEfficiency_AtMaximum << " with a significance of " << i_Signal_to_sqrtNoise_atMaximum << endl;
        cout << "\t Ndiff = " << Ndif << ", Nof " << i_BackgroundEfficiency_AtMaximum * Nof << endl;
        if( ( i_BackgroundEfficiency_AtMaximum * Nof ) > 1.e-10 )
        {
            cout << "\t Signal/background ratio: " << Ndif / ( i_BackgroundEfficiency_AtMaximum * Nof ) << endl;
        }
    }
    else
    {
        cout << "VTMVAEvaluator::optimizeSensitivity: no maximum in signal efficiency found" << endl;
        if( fOptimizationMaxSourceStrength > 0. &&
            iSourceStrength / fOptimizationMaxSourceStrength > 0.95 )
         {
             cout << "VTMVAEvaluator::optimizeSensitivity: (reached max source strength of ";
             cout << fOptimizationMaxSourceStrength << ")" << endl;
         }
    }
    cout << "\t MVA parameter: " << i_TMVACutValue_AtMaximum;
    cout << ", background efficiency: " << i_BackgroundEfficiency_AtMaximum << endl;
    if( i_AngularContainmentRadiusAtMaximum > 0. )
    {
        cout << "\t angular containment is " << i_AngularContainmentFractionAtMaximum * 100.;
        cout << "%, radius ";
        cout << i_AngularContainmentRadiusAtMaximum << " [deg]";
    }
    if( iHAngContainment )
    {
        cout << " (scaled from ";
        cout << iHAngContainment->GetBinContent( 
                iHAngContainment->GetXaxis()->FindBin( fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV ),
                iHAngContainment->GetYaxis()->FindBin( fTMVAngularContainmentRadiusMax ) );
        cout << " [deg], " << fTMVAngularContainmentRadiusMax * 100. << "%)";
    }
    cout << endl;
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // calculate TMVA cut value from q-factor
    // (will overwrite values calculated above)
    // IMPORTANT: this call should be commented out
    //            for default analysis
    /* optimizeSensitivity_using_qfactor( effS, effB,
                                       i_SignalEfficiency_AtMaximum,
                                       i_BackgroundEfficiency_AtMaximum,
                                       i_TMVACutValue_AtMaximum );  */
    ////////////////////////////////////////////////////////////////
    
    // get mean energy for this bin
    double iMeanEnergyAfterCuts = -99.;
    if( fDebug )
    {
        iMeanEnergyAfterCuts = getMeanEnergyAfterCut( iTMVAFile, i_TMVACutValue_AtMaximum, iDataBin );
        cout << "Mean energy after cuts [TeV]: " << iMeanEnergyAfterCuts << endl;
    }
    
    // fill results into data vectors
    fTMVAData[iDataBin]->fSignalEfficiency           = i_SignalEfficiency_AtMaximum;
    fTMVAData[iDataBin]->fBackgroundEfficiency       = i_BackgroundEfficiency_AtMaximum;
    fTMVAData[iDataBin]->fTMVACutValue               = i_TMVACutValue_AtMaximum;
    fTMVAData[iDataBin]->fSourceStrengthAtOptimum_CU = i_SourceStrength_atMaximum;
    if( iMeanEnergyAfterCuts > 0. )
    {
        fTMVAData[iDataBin]->fSpectralWeightedMeanEnergy_Log10TeV = log10( iMeanEnergyAfterCuts );
    }
    fTMVAData[iDataBin]->fAngularContainmentRadius = i_AngularContainmentRadiusAtMaximum;
    fTMVAData[iDataBin]->fAngularContainmentFraction = i_AngularContainmentFractionAtMaximum;
    
    // plot optimziation procedure and event numbers
    if( bPlotEfficiencyPlotsPerBin )
    {
        plotEfficiencyPlotsPerBin( iDataBin, iGSignal_to_sqrtNoise, iGSignal_to_sqrtNoise_Smooth,
                                   effS, effB, iGSignalEvents, iGBackgroundEvents,
                                   iGOpt_AngularContainmentRadius, iGOpt_AngularContainmentFraction );
    }
    
    return true;
}


/*
 * optimize sensitivity by determining
 * maximum in q-factor value
 *
 * Important: this is not the default analysis
 */
bool VTMVAEvaluator::optimizeSensitivity_using_qfactor( TH1D* effS, TH1D* effB,
        double& i_SignalEfficiency_AtMaximum,
        double& i_BackgroundEfficiency_AtMaximum,
        double& i_TMVACutValue_AtMaximum )
{
    i_SignalEfficiency_AtMaximum = 0.;
    i_BackgroundEfficiency_AtMaximum = 0.;
    i_TMVACutValue_AtMaximum = 0.;
    // make sure that signal and background efficiency
    // histograms exist
    if( !effS || !effB )
    {
        return false;
    }
    
    double q = -1.;
    double eS = 0.;
    double eB = 0.;
    double q_max = 0.;
    
    for( int i = 1; i <= effS->GetNbinsX(); i++ )
    {
        eS = effS->GetBinContent( i );
        eB = effB->GetBinContent( i );
        if( eB > 0. )
        {
            q = eS / sqrt( eB * 0.5 );
            if( q > q_max )
            {
                i_TMVACutValue_AtMaximum = effS->GetXaxis()->GetBinCenter( i );
                i_SignalEfficiency_AtMaximum = eS;
                i_BackgroundEfficiency_AtMaximum = eB;
            }
        }
    }
    
    return true;
}



/*

 smoothing of optimal cut value vs energy curves

 missing (non-optimized) are Interpolated

 note: signal and background efficiencies are not updated

*/
void VTMVAEvaluator::smoothAndInterpolateMVAValue( 
        unsigned int iWeightFileIndex_Emin,
        unsigned int iWeightFileIndex_Emax,
        unsigned int iWeightFileIndex_Zmin,
        unsigned int iWeightFileIndex_Zmax,
        double iEnergyStepSize )
{
    if( fTMVAData.size() == 0 )
    {
        return;
    }
    
    cout << "Smooth and Interpolate MVA cut values" << endl;
    
    int z = 0;
    
    // take zenith angle dependence into account
    unsigned int Ebins = iWeightFileIndex_Emax - iWeightFileIndex_Emin + 1;
    unsigned int ZEbins = iWeightFileIndex_Zmax - iWeightFileIndex_Zmin + 1;
    
    //////////////////////////////////////////////
    // energy dependent TMVA cut optimization only
    if( ZEbins == 1 )
    {
        // fill graph to be smoothed
        TGraph* iG = new TGraph( 1 );
        for( unsigned int i = 0; i < fTMVAData.size(); i++ )
        {
            if( fTMVAData[i] )
            {
                if( fTMVAData[i]->fTMVAOptimumCutValueFound )
                {
                    iG->SetPoint( z, fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[i]->fTMVACutValue );
                    z++;
                }
            }
        }
        // smooth graph
        TGraph* iGout = new TGraph( 1 );
        TGraphSmooth* iGSmooth = new TGraphSmooth( "t" );
        iGout = ( TGraph* )iGSmooth->SmoothKern( iG, "normal", 0.5, 100 );
        
        // fill smoothed and Interpolated values into MVA vector
        // set all points to 'optimized'
        for( unsigned int i = 0; i < fTMVAData.size(); i++ )
        {
            if( fTMVAData[i] )
            {
                cout << "\t TMVA values: unsmoothed at " << TMath::Power( 10., fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                cout << " TeV, \t" << fTMVAData[i]->fTMVACutValue;
                fTMVAData[i]->fTMVACutValue = iGout->Eval( fTMVAData[i]->fSpectralWeightedMeanEnergy_Log10TeV );
                cout << ", smoothed " << fTMVAData[i]->fTMVACutValue;
                if( !fTMVAData[i]->fTMVAOptimumCutValueFound )
                {
                    cout << " (Interpolated non-optimal value)";
                }
                cout << " (" << i << ")" << endl;
                fTMVAData[i]->fTMVAOptimumCutValueFound = true;
                
                // get efficiency histograms
                TFile iTMVAFile( fTMVAData[i]->fTMVAFileName.c_str() );
                TH1D* effS = getEfficiencyHistogram( "effS", &iTMVAFile, fTMVAData[i]->fTMVAMethodTag_2 );
                TH1D* effB = getEfficiencyHistogram( "effB", &iTMVAFile, fTMVAData[i]->fTMVAMethodTag_2 );
                
                if( effS )
                {
                    fTMVAData[i]->fSignalEfficiency = effS->GetBinContent( effS->FindBin( fTMVAData[i]->fTMVACutValue ) );
                }
                if( effB )
                {
                    fTMVAData[i]->fBackgroundEfficiency = effB->GetBinContent( effB->FindBin( fTMVAData[i]->fTMVACutValue ) );
                    // background efficiency might be zero -> fill it with first non-zero value
                    if( fTMVAData[i]->fBackgroundEfficiency < 1.e-8 )
                    {
                        int iS = effB->FindBin( fTMVAData[i]->fTMVACutValue );
                        for( int j = iS; j > 0; j-- )
                        {
                            if( effB->GetBinContent( j ) > 0. )
                            {
                                fTMVAData[i]->fBackgroundEfficiency = effB->GetBinContent( j );
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    ///////////////////////////////////////////////////
    // energy and zenith angle dependent smoothing (2D)
    else
    {
        //check if overlapping energy bins are used
        if( iEnergyStepSize > 0 )
        {
            Ebins = fTMVAData.size() / ZEbins;
        }
        
        // fill 2D graph to be smoothed (energy and zenith angle dependence)
        TGraph2D* iG2 = new TGraph2D( 1 );
        Double_t effS_value[Ebins * ZEbins];
        for( unsigned int l = 0; l < fTMVAData.size(); l++ )
        {
            if( fTMVAData[l] )
            {
                TFile iTMVAFile( fTMVAData[l]->fTMVAFileName.c_str() );
                if( !iTMVAFile.IsZombie() )
                {
                    TH1D* effS = getEfficiencyHistogram( "effS", &iTMVAFile, fTMVAData[l]->fTMVAMethodTag_2 );
                    if( effS )
                    {
                        effS_value[l] = effS->GetBinContent( effS->FindBin( fTMVAData[l]->fTMVACutValue ) );
                    }
                }
            }
        }
        
        for( unsigned int l = 0; l < fTMVAData.size(); l++ )
        {
            if( fTMVAData[l] )
            {
                // get signal efficiency histograms and cut values
                TFile iTMVAFile( fTMVAData[l]->fTMVAFileName.c_str() );
                if( iTMVAFile.IsZombie() )
                {
                    continue;
                }
                TH1D* effS = getEfficiencyHistogram( "effS", &iTMVAFile, fTMVAData[l]->fTMVAMethodTag_2 );
                if( !effS )
                {
                    continue;
                }
                // bins without optimal cut value
                if( !fTMVAData[l]->fTMVAOptimumCutValueFound && l>0 )
                {
                    for( int k = 0; k < effS->GetNbinsX(); k++ )
                    {
                        if( TMath::Abs( effS->GetBinContent( k ) - effS_value[l-1] ) < 0.001 )
                        {
                            fTMVAData[l]->fTMVACutValue = effS->GetBinCenter( k );
                            effS_value[l] = effS->GetBinContent( k );
                        }
                    }
                }
                if( fTMVAData[l]->fTMVACutValue == -99 )
                {
                    cout << "Error: no optimal cut value found for this bin" << endl;
                }
                cout << "\t TMVA values for bin " << l << ":" << endl;
                cout << "\t\t energy bin (" << fTMVAData[l]->fEnergyCut_bin << "): ";
                cout << "[" << fTMVAData[l]->fEnergyCut_Log10TeV_min;
                cout << ", " << fTMVAData[l]->fEnergyCut_Log10TeV_max << "] (log TeV), ";
                cout << "E_mean: " << TMath::Power( 10., fTMVAData[l]->fSpectralWeightedMeanEnergy_Log10TeV );
                cout << " TeV, " << endl;
                cout << "\t\t zenith bin (" << fTMVAData[l]->fZenithCut_bin << "): [";
                cout << fTMVAData[l]->fZenithCut_min << ", ";
                cout << fTMVAData[l]->fZenithCut_max << "] deg: " << endl;
                cout << "\t\t cut value: ";
                cout << fTMVAData[l]->fTMVACutValue;
                if( !fTMVAData[l]->fTMVAOptimumCutValueFound )
                {
                    cout << " (Interpolated non-optimal value)";
                }
                cout << ", signal efficiency: " << effS_value[l];
                cout << " (bin " << l << ")" << endl;
                
                iG2->SetPoint( l, 
                             fTMVAData[l]->fSpectralWeightedMeanEnergy_Log10TeV,
                             fTMVAData[l]->fZenithCut_min, 
                             fTMVAData[l]->fTMVACutValue );
                //a stupid kludge to handle the edges
                iG2->SetPoint( fTMVAData.size() + l, 
                             fTMVAData[l]->fSpectralWeightedMeanEnergy_Log10TeV + 0.01, 
                             fTMVAData[l]->fZenithCut_min + 0.01, 
                             fTMVAData[l]->fTMVACutValue );
            }
        }
        
        // smooth the 2D histogram
        if( iG2 )
        {
            TH2D* iH2 = iG2->GetHistogram();
            iH2->Smooth();
            for( unsigned int counter = 0; counter < fTMVAData.size(); counter++ )
            {
                if( fTMVAData[counter] )
                {
                    cout << "\t smoothed TMVA values: at " << TMath::Power( 10., fTMVAData[counter]->fSpectralWeightedMeanEnergy_Log10TeV );
                    cout << " TeV, ";
                    cout << "[" << fTMVAData[counter]->fZenithCut_min << ", " << fTMVAData[counter]->fZenithCut_max << "] deg: ";
                    cout << iH2->GetBinContent( iH2->FindBin( fTMVAData[counter]->fSpectralWeightedMeanEnergy_Log10TeV, fTMVAData[counter]->fZenithCut_min ) );
                    cout << " (" << counter << ")" << endl;
                }
            }
            cout << "\t (smoothed values printed here only - they are nowhere used in the analysis)" << endl;
        }
    }
}

/*
 * optimise cut on angular containment radius (theta2 cut)
 *
 */
void VTMVAEvaluator::getOptimalAngularContainmentRadius( double effS, double effB, double Ndif, double Nof,
        TH2D* iHAngContainment, double iEnergy_log10_TeV,
        double& i_Signal_to_sqrtNoise, double& i_AngularContainmentRadius,
        double& i_AngularContainmentFraction )
{
    i_AngularContainmentFraction = 0.;
    i_Signal_to_sqrtNoise = -99.;
    i_AngularContainmentRadius = -99.;
    if( !iHAngContainment || fOptimizationBackgroundAlpha <= 0 )
    {
        return;
    }
    double iR_Max = VHistogramUtilities::interpolateTH2D( iHAngContainment, iEnergy_log10_TeV, fTMVAngularContainmentRadiusMax );
    
    // find containment radius giving maximum significance
    double iC = 0.;
    double iR = 0.;
    double iOn = 0.;
    double iOff = 0.;
    double iSigma = 0.;
    for( int i = 1; i < iHAngContainment->GetNbinsY(); i++ )
    {
        iC = iHAngContainment->GetYaxis()->GetBinLowEdge( i );
        iR = VHistogramUtilities::interpolateTH2D( iHAngContainment, iEnergy_log10_TeV, iHAngContainment->GetYaxis()->GetBinLowEdge( i ) );
        
        if( iC > fTMVAngularContainmentRadiusMax )
        {
            iR = iR_Max;
            iC = fTMVAngularContainmentRadiusMax;
        }
        // make sure that containment cut is not getting too small
        if( iR < fTMVAAngularContainmentThetaFixedMinRadius )
        {
            continue;
        }
        
        // gamma point source
        iOn  = effS * Ndif * iC / fTMVAngularContainmentRadiusMax;
        iOn += effB * Nof * iR * iR / iR_Max / iR_Max;
        iOff = effB * Nof / fOptimizationBackgroundAlpha * iR * iR / iR_Max / iR_Max;
        
        iSigma = VStatistics::calcSignificance( iOn, iOff, fOptimizationBackgroundAlpha );
        
        if( iSigma >= i_Signal_to_sqrtNoise )
        {
            i_Signal_to_sqrtNoise = iSigma;
            i_AngularContainmentRadius = iR;
            i_AngularContainmentFraction = iC;
        }
    }
}


void VTMVAEvaluator::plotEfficiencyPlotsPerBin( unsigned int iBin, TGraph* iGSignal_to_sqrtNoise,
        TGraph* iGSignal_to_sqrtNoise_Smooth, TH1D* hEffS, TH1D* hEffB,
        TGraph* iGSignalEvents, TGraph* iGBackgroundEvents,
        TGraph* iGOpt_AngularContainmentRadius, TGraph* iGOpt_AngularContainmentFraction,
        bool bPlotContainmentFraction )
{
    char hname[800];
    char htitle[800];
    
    if( iBin >= fTMVAData.size() || !fTMVAData[iBin] )
    {
        return;
    }
    
    // signal and noise plot
    if( hEffS )
    {
        sprintf( hname, "cEfficiencyPlotPerEnergy_%d", iBin );
        sprintf( htitle, "efficiency plots (%s)", fTMVAData[iBin]->GetTitle() );
        TCanvas* iCanvas = new TCanvas( hname, htitle, 10, 10 + iBin * 30, 400, 400 );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        iCanvas->SetLeftMargin( 0.13 );
        iCanvas->Draw();
        
        hEffS->SetStats( 0 );
        hEffS->SetTitle( "" );
        hEffS->SetLineWidth( 3 );
        hEffS->GetYaxis()->SetTitleOffset( 1.5 );
        hEffS->SetXTitle( "MVA value #Tau" );
        hEffS->SetYTitle( "signal/background efficiency" );
        hEffS->DrawCopy();
        
        if( hEffB )
        {
            hEffB->SetStats( 0 );
            hEffB->SetTitle( "" );
            hEffB->SetLineColor( 2 );
            hEffB->SetLineWidth( 3 );
            hEffB->SetLineStyle( 9 );
            hEffB->DrawCopy( "same" );
        }
        
        if( iBin < fTMVAData.size() && fTMVAData[iBin] )
        {
            TLine* iL = new TLine( fTMVAData[iBin]->fTMVACutValue, hEffS->GetMinimum(), fTMVAData[iBin]->fTMVACutValue, hEffS->GetMaximum() );
            iL->SetLineStyle( 2 );
            iL->SetLineWidth( 3 );
            iL->Draw();
        }

        if( fPrintPlotting )
        {
            sprintf( hname, "MVAOpt-SignalBackgroundEfficiency-%d.pdf", iBin );
            iCanvas->Print( hname );
        }
    }
    
    // signal to noise
    if( iGSignal_to_sqrtNoise )
    {
        sprintf( hname, "cSignalToSqrtNoise_%d", iBin );
        sprintf( htitle, "signal / sqrt( noise ) (%s)", fTMVAData[iBin]->GetTitle() );
        TCanvas* iCanvas = new TCanvas( hname, htitle, 425, 10 + iBin * 30, 400, 400 );
        iCanvas->SetLeftMargin( 0.13 );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        
        iGSignal_to_sqrtNoise->SetTitle( "" );
        setGraphPlottingStyle( iGSignal_to_sqrtNoise, 4, 3., 20, 1., 0, 1 );
        
        iGSignal_to_sqrtNoise->Draw( "apl" );
        iGSignal_to_sqrtNoise->GetHistogram()->GetYaxis()->SetTitleOffset( 1.5 );
        iGSignal_to_sqrtNoise->GetHistogram()->SetXTitle( "MVA value #Tau" );
        iGSignal_to_sqrtNoise->GetHistogram()->SetYTitle( "significance" );
        
        if( iGSignal_to_sqrtNoise_Smooth )
        {
            setGraphPlottingStyle( iGSignal_to_sqrtNoise_Smooth, 2, 2. );
            iGSignal_to_sqrtNoise_Smooth->Draw( "pl" );
        }
        
        if( iBin < fTMVAData.size() && fTMVAData[iBin] )
        {
            TLine* iL = new TLine( fTMVAData[iBin]->fTMVACutValue, iGSignal_to_sqrtNoise->GetHistogram()->GetMinimum(),
                                   fTMVAData[iBin]->fTMVACutValue, iGSignal_to_sqrtNoise->GetHistogram()->GetMaximum() );
            iL->SetLineStyle( 2 );
            iL->Draw();
        }
        if( fPrintPlotting )
        {
            sprintf( hname, "MVAOpt-SignalToSqurtNoise-%d.pdf", iBin );
            iCanvas->Print( hname );
        }
    }
    
    // signal and background events numbers
    if( iGBackgroundEvents )
    {
        sprintf( hname, "cEventNumbers_%d", iBin );
        sprintf( htitle, "event numbers (%s)", fTMVAData[iBin]->GetTitle() );
        TCanvas* iCanvas = new TCanvas( hname, htitle, 850, 10 + iBin * 30, 400, 400 );
        iCanvas->SetLeftMargin( 0.13 );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        
        sprintf( hname, "hBC_%d", iBin );
        TH1D* hnull = new TH1D( hname, "", 100, -1., 1. );
        hnull->SetXTitle( "cut value" );
        hnull->SetYTitle( "number of events" );
        hnull->SetMinimum( 1.e-3 );
        double x = 0.;
        double y = 0.;
        double y_max = 0.;
        for( int i = 0; i < iGBackgroundEvents->GetN(); i++ )
        {
            iGBackgroundEvents->GetPoint( i, x, y );
            if( y > y_max )
            {
                y_max = y;
            }
        }
        for( int i = 0; i < iGSignalEvents->GetN(); i++ )
        {
            iGSignalEvents->GetPoint( i, x, y );
            if( y > y_max )
            {
                y_max = y;
            }
        }
        hnull->SetMaximum( y_max * 1.5 );
        hnull->SetStats( 0 );
        hnull->GetYaxis()->SetTitleOffset( 1.5 );
        hnull->DrawCopy();
        
        setGraphPlottingStyle( iGBackgroundEvents, 1, 1., 20 );
        iGBackgroundEvents->Draw( "pl" );
        
        if( iGSignalEvents )
        {
            setGraphPlottingStyle( iGSignalEvents, 2, 2. );
            iGSignalEvents->Draw( "pl" );
        }
        
        if( iBin < fTMVAData.size() && fTMVAData[iBin] )
        {
            TLine* iL = new TLine( fTMVAData[iBin]->fTMVACutValue, hnull->GetMinimum(),
                                   fTMVAData[iBin]->fTMVACutValue, hnull->GetMaximum() );
            iL->SetLineStyle( 2 );
            iL->Draw();
            TLine* iLMinBack = new TLine( hnull->GetXaxis()->GetXmin(), fOptimizationMinBackGroundEvents,
                                          hnull->GetXaxis()->GetXmax(), fOptimizationMinBackGroundEvents );
            iLMinBack->SetLineStyle( 2 );
            iLMinBack->Draw();
        }
    }
    
    // angular containment radius
    if( iGOpt_AngularContainmentRadius && bPlotContainmentFraction )
    {
        sprintf( hname, "cCR_%d", iBin );
        sprintf( htitle, "containment radius (%s)", fTMVAData[iBin]->GetTitle() );
        TCanvas* iCanvas = new TCanvas( hname, htitle, 750, 10 + iBin * 30, 400, 400 );
        iCanvas->SetLeftMargin( 0.13 );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        
        sprintf( hname, "hBA_%d", iBin );
        TH1D* hnull = new TH1D( hname, "", 100, -1., 1. );
        hnull->SetXTitle( "cut value" );
        hnull->SetYTitle( "containment radius [deg]" );
        hnull->SetMinimum( 1.e-3 );
        double x = 0.;
        double y = 0.;
        double y_max = 0.;
        for( int i = 0; i < iGOpt_AngularContainmentRadius->GetN(); i++ )
        {
            iGOpt_AngularContainmentRadius->GetPoint( i, x, y );
            if( y > y_max )
            {
                y_max = y;
            }
        }
        if( y_max < 1.e-3 )
        {
            y_max = 1.;
        }
        hnull->SetMaximum( y_max * 1.5 );
        hnull->SetStats( 0 );
        hnull->GetYaxis()->SetTitleOffset( 1.5 );
        hnull->DrawCopy();
        
        setGraphPlottingStyle( iGOpt_AngularContainmentRadius, 1, 1., 20 );
        iGOpt_AngularContainmentRadius->Draw( "pl" );
        if( iBin < fTMVAData.size() && fTMVAData[iBin] )
        {
            TLine* iL = new TLine( fTMVAData[iBin]->fTMVACutValue, hnull->GetMinimum(), fTMVAData[iBin]->fTMVACutValue, y_max * 1.5 );
            iL->SetLineStyle( 2 );
            iL->SetLineWidth( 3 );
            iL->Draw();
        }
    }
    // angular containment fraction
    if( iGOpt_AngularContainmentFraction && bPlotContainmentFraction )
    {
        sprintf( hname, "cCF_%d", iBin );
        sprintf( htitle, "containment fraction (%s)", fTMVAData[iBin]->GetTitle() );
        TCanvas* iCanvas = new TCanvas( hname, htitle, 750, 10 + iBin * 30, 400, 400 );
        iCanvas->SetLeftMargin( 0.13 );
        iCanvas->SetGridx( 0 );
        iCanvas->SetGridy( 0 );
        
        sprintf( hname, "hBF_%d", iBin );
        TH1D* hnull = new TH1D( hname, "", 100, -1., 1. );
        hnull->SetXTitle( "cut value" );
        hnull->SetYTitle( "containment fraction" );
        hnull->SetMinimum( 1.e-3 );
        double x = 0.;
        double y = 0.;
        double y_max = 0.;
        for( int i = 0; i < iGOpt_AngularContainmentFraction->GetN(); i++ )
        {
            iGOpt_AngularContainmentFraction->GetPoint( i, x, y );
            if( y > y_max )
            {
                y_max = y;
            }
        }
        hnull->SetMaximum( y_max * 1.5 );
        hnull->SetStats( 0 );
        hnull->GetYaxis()->SetTitleOffset( 1.5 );
        hnull->DrawCopy();
        
        setGraphPlottingStyle( iGOpt_AngularContainmentFraction, 1, 1., 20 );
        iGOpt_AngularContainmentFraction->Draw( "pl" );
        if( iBin < fTMVAData.size() && fTMVAData[iBin] )
        {
            TLine* iL = new TLine( fTMVAData[iBin]->fTMVACutValue, hnull->GetMinimum(), fTMVAData[iBin]->fTMVACutValue, y_max * 1.5 );
            iL->SetLineStyle( 2 );
            iL->SetLineWidth( 3 );
            iL->Draw();
        }
    }
    
}

void VTMVAEvaluator::setTMVAMethod( string iMethodName, int iMethodCounter )
{
    fTMVAMethodName = iMethodName;
    fTMVAMethodCounter = iMethodCounter;
}

void VTMVAEvaluator::setSensitivityOptimizationFixedSignalEfficiency( double iOptimizationFixedSignalEfficiency )
{
    fOptimizationFixedSignalEfficiency = iOptimizationFixedSignalEfficiency;
}

void VTMVAEvaluator::setSensitivityOptimizationSourceStrength( 
                            double iOptimizationMinSourceStrength,
                            double iOptimizationMaxSourceStrength )
{
    fOptimizationMinSourceStrength = iOptimizationMinSourceStrength;
    fOptimizationMaxSourceStrength = iOptimizationMaxSourceStrength;
}

/*
 * MVA cut or signal efficiency is stored in a map:
 *
 * e.g. energy bin 3 and zenith angle bin 2: entry 21
 *
 */
double VTMVAEvaluator::getValueFromMap( map< unsigned int, double > iDataMap, double iDefaultData,
                                        unsigned int iEnergyBin, double iE_min_log10, double iE_max_log10,
                                        unsigned int iZenithBin, double iZ_min, double iZ_max, 
                                        unsigned int iNCut, double iEnergyStepSize, string iVariable )
{
    if( iDataMap.size() == 0 )
    {
        return iDefaultData;
    }

    // iIter->first:
    // [energy bins*10+zenith bin]
    map< unsigned int, double >::iterator iIter;
    for( iIter = iDataMap.begin(); iIter != iDataMap.end(); ++iIter )
    {
        // data does not depend on energy
        if( iIter->first > 9998 )
        {
            return iIter->second;
        }
        
        // data signal efficiency
        bool iFoundCut = false;
        if( iIter->first == ( iNCut - 1 ) && iEnergyStepSize > 0 )
        {
            iFoundCut = true;
        }
        if( iIter->first == ( iEnergyBin * 10 + iZenithBin ) && iEnergyStepSize < 0 )
        {
            iFoundCut = true;
        }
        
        if( iFoundCut == true )
        {
            cout << "VTMVAEvaluator::getValueFromMap (" << iVariable << "): ";
            cout << " key: " << iIter->first << ", value: " << iIter->second << endl;
            cout << "\t Energy [" << iE_min_log10 << ", " << iE_max_log10 << "], bin ";
            cout << iEnergyBin << endl;
            cout << "\t Zenith angle [" << iZ_min << ", " << iZ_max << "], bin ";
            cout << iZenithBin << endl;
            return iIter->second;
        }
        
    }
    
    cout << "VTMVAEvaluator::getValueFromMap: warning, couldn't find a data value (" << iVariable << ") for energy bin ";
    cout << iEnergyBin << ", E=[ " << iE_min_log10 << ", " << iE_max_log10 << "] and zenith bin ";
    cout << iZenithBin << ", Zen=[ " << iZ_min << ", " << iZ_max << "] " << endl;
    
    return -1.;
}

/*
 * return signal efficiency for the given energy
 *
 */
double VTMVAEvaluator::getSignalEfficiency( unsigned int iEnergyBin, 
                                            double iE_min_log10, double iE_max_log10, 
                                            unsigned int iZenithBin, 
                                            double iZ_min, double iZ_max, 
                                            unsigned int iNCut, double iEnergyStepSize )
{
    // read signal efficiency from graph
    if( fTMVASignalEfficencyGraph.size() > 0
       && iZenithBin < fTMVASignalEfficencyGraph.size()
       && fTMVASignalEfficencyGraph[iZenithBin] )
    {
         return fTMVASignalEfficencyGraph[iZenithBin]->Eval( 0.5*(iE_min_log10+iE_max_log10) );
    }
    // read signal efficiency from a map
    return getValueFromMap( fSignalEfficiencyMap, fSignalEfficiencyNoVec, 
                            iEnergyBin, iE_min_log10, iE_max_log10, 
                            iZenithBin, iZ_min, iZ_max, iNCut, 
                            iEnergyStepSize, "SignalEfficiency" );
}

double VTMVAEvaluator::getTMVACutValue( unsigned int iEnergyBin, 
                                        double iE_min_log10, double iE_max_log10, 
                                        unsigned int iZenithBin, 
                                        double iZ_min, double iZ_max, 
                                        unsigned int iNCut, double iEnergyStepSize )
{
    // read MVA cut value from graph
    if( fTMVACutValueGraph.size() > 0
       && iZenithBin < fTMVACutValueGraph.size()
       && fTMVACutValueGraph[iZenithBin] )
    {
         return fTMVACutValueGraph[iZenithBin]->Eval( 0.5*(iE_min_log10+iE_max_log10) );
    }

    // read MVA cut value from map
    return getValueFromMap( fTMVACutValueMap, fTMVACutValueNoVec, 
                            iEnergyBin, iE_min_log10, iE_max_log10, 
                            iZenithBin, iZ_min, iZ_max, iNCut, 
                            iEnergyStepSize, "MVA_CUT" );
}



void VTMVAEvaluator::printSensitivityOptimizationParameters()
{
    cout << "VTMAEvaluator: MVA cut parameter is optimized for: " << endl;
    cout << "\t" << fOptimizationObservingTime_h << " hours of observing time" << endl;
    cout << "\t" << fOptimizationSourceSignificance << " minimum significance" << endl;
    cout << "\t" << fOptimizationMinSignalEvents << " minimum number of on events" << endl;
    cout << "\t" << fOptimizationBackgroundAlpha << " signal to background area ratio" << endl;
    cout << "\t" << fMinBackgroundRateRatio_min << " minimum signal-to-background ratio" << endl;
}

/*
 * calculate mean energy in an energy bin after applying MVA cuts
 * (use training tree TrainTree out of e.g. BDT file for this)
 */
double VTMVAEvaluator::getMeanEnergyAfterCut( TFile* f, double iCut, unsigned int iDataBin )
{
    if( !f )
    {
        return -99.;
    }
    if( iDataBin >= fTMVAData.size() || !fTMVAData[iDataBin] )
    {
        return -99.;
    }
    double iEmin = TMath::Power( 10., fTMVAData[iDataBin]->fEnergyCut_Log10TeV_min );
    double iEmax = TMath::Power( 10., fTMVAData[iDataBin]->fEnergyCut_Log10TeV_max );
    TTree* t = ( TTree* )f->Get( "TrainTree" );
    if( !t )
    {
        cout << "VTMVAEvaluator::getMeanEnergyAfterCut(): test tree not found in " << f->GetName() << endl;
        return -99.;
    }
    float iErec = 0.;
    float iMVA = 0.;
    int classID = 0;;
    t->SetBranchAddress( "ErecS", &iErec );
    ostringstream iCutName;
    // variable names changed with time - keep backwards compatibility
    iCutName << fTMVAMethodName << "_" << fTMVAData[iDataBin]->fTMVAMethodTag_2;
    if( t->GetBranchStatus( iCutName.str().c_str() ) )
    {
        t->SetBranchAddress( iCutName.str().c_str(), &iMVA );
    }
    else
    {
        iCutName.clear();
        iCutName.str( std::string() );
        iCutName << fTMVAMethodName << "_" << fTMVAMethodCounter;
        t->SetBranchAddress( iCutName.str().c_str(), &iMVA );
    }
    t->SetBranchAddress( "classID", &classID );
    
    float n = 0.;
    float m = 0.;
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        if( classID == 0 && iErec > 0. && iErec > iEmin && iErec < iEmax && iMVA > iCut )
        {
            m += iErec;
            n++;
        }
    }
    if( n > 0. )
    {
        return m / n;
    }
    
    return -99.;
}

vector< double > VTMVAEvaluator::getBackgroundEfficiency()
{
    vector< double > iA;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            iA.push_back( fTMVAData[i]->fBackgroundEfficiency );
        }
    }
    return iA;
}

vector< double > VTMVAEvaluator::getSignalEfficiency()
{
    vector< double > iA;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            iA.push_back( fTMVAData[i]->fSignalEfficiency );
        }
    }
    return iA;
}

vector< double > VTMVAEvaluator::getTMVACutValue()
{
    vector< double > iA;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            iA.push_back( fTMVAData[i]->fTMVACutValue );
        }
    }
    return iA;
}

vector< bool > VTMVAEvaluator::getOptimumCutValueFound()
{
    vector< bool > iA;
    for( unsigned int i = 0; i < fTMVAData.size(); i++ )
    {
        if( fTMVAData[i] )
        {
            iA.push_back( fTMVAData[i]->fTMVAOptimumCutValueFound );
        }
    }
    return iA;
}

/*
 * read graph with on/off numbers as a function of energy and possibly zenith angle
 * Interpolates between zenith range and return a TGraph
 * zenith angles given given in secants
 *
 */
TGraph* VTMVAEvaluator::readInterpolatedCountsFromFile( TFile* iF, double i_secant_min, double i_secant_max, bool bIsOn )
{
    if( !iF && iF->IsZombie() )
    {
        return 0;
    }
    TObject* i_N = 0;
    if( bIsOn )
    {
        i_N = iF->Get( "gSignalRate" );
    }
    else
    {
        i_N = iF->Get( "gBGRate" );
    }
    return ( TGraph* )getInterpolatedDifferentialRatesfromGraph2D( i_N, i_secant_min, i_secant_max );
}

/*
 * read rates with signal and off numbers as a function of energy and possibly zenith angle
 * return a double
 *
 * zenith angles are given in secants
 *
*/
double VTMVAEvaluator::readAverageCountsFromFile( TFile* iF, double i_e_min, double i_e_max,
        double i_secant_min, double i_secant_max, bool bIsOn )
{
    if( !iF && iF->IsZombie() )
    {
        return 0;
    }
    TObject* i_N = 0;
    if( bIsOn )
    {
        i_N = iF->Get( "gSignalRate" );
    }
    else
    {
        i_N = iF->Get( "gBGRate" );
    }
    return ( double )getAverageDifferentialRateFromGraph2D( i_N, i_e_min, i_e_max, i_secant_min, i_secant_max );
}

/*
 * fill by Interpolation a TGraph2D into TGraph (if necessary)
 *
 */
TGraph* VTMVAEvaluator::getInterpolatedDifferentialRatesfromGraph2D( TObject* i_G, double i_secant_min, double i_secant_max )
{
    if( !i_G )
    {
        return 0;
    }

    // range and step size of graph to be returned
    double i_E_min_log10 = -2.;
    double i_E_max_log10 = 3.;
    double i_E_stepsize = 0.05;
    int npx = (int)( (i_E_max_log10-i_E_min_log10)/i_E_stepsize );
    
    string i_c = i_G->ClassName();
    // graph is 2D: Interpolate and fill into a TGraph
    if( i_c.find( "TGraph2D" ) != string::npos )
    {
        TGraph2D* iG2D = ( TGraph2D* )i_G;
        TGraph* iG1D = new TGraph( npx );

        for( int i = 0; i < npx; i++ )
        {
             double e = i_E_min_log10 + i_E_stepsize * i;

             iG1D->SetPoint( i, e, iG2D->Interpolate( e, 0.5 * ( i_secant_min + i_secant_max ) ) );
        }
        return iG1D;
    }
    
    // graph is already 1D
    return ( TGraph* )i_G;
}

/*
 * get average non/noff rate for a given energy and secant range
 *
 */
double VTMVAEvaluator::getAverageDifferentialRateFromGraph2D( TObject* i_G,
                        double i_e_min, double i_e_max,
                        double i_secant_min, double i_secant_max )
{
    if( !i_G )
    {
        return 0.;
    }
    
    string i_c = i_G->ClassName();
    // graph is 2D: Interpolate and fill into a TGraph
    if( i_c.find( "TGraph2D" ) != string::npos )
    {
        TGraph2D* iG2D = ( TGraph2D* )i_G;
        Double_t* energypoint = iG2D->GetX();
        Double_t* zenithpoint = iG2D->GetY();
        
        Double_t* rate = iG2D->GetZ();
        Int_t counter = 0;
        Double_t ratesum = 0;
        
        for( int i = 0; i < iG2D->GetN(); i++ )
        {
            if( zenithpoint[i] >= i_secant_min && zenithpoint[i] <= i_secant_max
                    && energypoint[i] >= i_e_min && energypoint[i] <= i_e_max )
            {
                ratesum += rate[i];
                counter++;
            }
        }
        if( counter > 0. )
        {
            return ratesum / counter;
        }
    }
    return 0.;
}

/* 
 * return correct mva root and xml file
 * check if file exists
 *
 */
string VTMVAEvaluator::setFullMVAFileName( string iWeightFileName,
                                     unsigned int iWeightFileIndex_Emin, unsigned int i,
                                     unsigned int iWeightFileIndex_Zmin, unsigned int j,
                                     string fTMVAMethodName, int fTMVAMethodCounter,
                                     string iInstrumentEpoch,
                                     string iFileSuffix )
{
      ostringstream iFileName;
      ostringstream iFileNamev2;

      iFileName << iWeightFileName << "_Ebin" << iWeightFileIndex_Emin + i;
      iFileNamev2 << iWeightFileName << "_" << iWeightFileIndex_Emin + i;
      // backwards compatibility with pre-2018 CTA training
      if( iInstrumentEpoch != "noepoch" && iInstrumentEpoch != "CTA" )
      {
          iFileName << "_Zebin" << iWeightFileIndex_Zmin + j;
          iFileNamev2 << "_Zebin" << iWeightFileIndex_Zmin + j;
      }
      if( iFileSuffix.find( "xml" ) != string::npos )
      {
          iFileName << "_" << fTMVAMethodName << "_" << fTMVAMethodCounter;
          iFileNamev2 << "_" << fTMVAMethodName << "_" << fTMVAMethodCounter;
      }
      iFileName << iFileSuffix;
      iFileNamev2 << iFileSuffix;
      // check if file exists or if this is an old-style file
      ifstream f(iFileName.str().c_str());
      if( !f.good() )
      {
          ifstream f2(iFileNamev2.str().c_str());
          if( !f2.good() )
          {
              cout << "VTMVAEvaluator::setFullMVAFileName warning: file not found: " << endl;
              cout << iFileName.str() << endl;
              cout << iFileNamev2.str() << endl;
              return "";
          }
          else
          {
              return iFileNamev2.str();
          }
      }

      return iFileName.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

VTMVAEvaluatorData::VTMVAEvaluatorData()
{
    fCounter = 0;
    fTMVAFileName = "";
    fTMVAFileNameXML = "";
    fTMVAMethodTag = "";
    fEnergyCut_Log10TeV_min = -99.;
    fEnergyCut_Log10TeV_max = -99.;
    fSpectralWeightedMeanEnergy_Log10TeV = -99.;
    fZenithCut_min = -99.;
    fZenithCut_max = -99.;
    
    fSignalEfficiency = -99.;
    fBackgroundEfficiency = -99.;
    fTMVACutValue = -99.;
    fTMVAOptimumCutValueFound = false;
    fSourceStrengthAtOptimum_CU = -99.;
    fAngularContainmentRadius = -99.;
    fAngularContainmentFraction = -99.;

    hSignalEfficiency = 0;
    hBackgroundEfficiency = 0;

    // transients
    fTMVAReader = 0;
}

void VTMVAEvaluatorData::print()
{
    cout << "\t file " << fTMVAFileName << endl;
    cout << "\t energy bin [" << fEnergyCut_Log10TeV_min << "," << fEnergyCut_Log10TeV_max;
    cout << "] (mean energy " << fSpectralWeightedMeanEnergy_Log10TeV << ")";
    cout << ", zenith bin [" << fZenithCut_min << "," << fZenithCut_max << "]";
    cout << endl;
}

