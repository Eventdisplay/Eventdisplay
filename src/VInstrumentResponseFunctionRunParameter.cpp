/*! \class VInstrumentResponseFunctionRunParameter
    \brief run parameters for response function calculator (effective areas)

*/

#include "VInstrumentResponseFunctionRunParameter.h"

VInstrumentResponseFunctionRunParameter::VInstrumentResponseFunctionRunParameter()
{
    fFillingMode = 0;
    fEffArea_short_writing = false;
    
    fInstrumentEpoch = "NOT_SET";
    
    fNSpectralIndex = 1;
    fSpectralIndexMin = 2.0;
    fSpectralIndexStep = 0.1;
    
    fReconstructionType = NOT_SET;
    fEnergyReconstructionMethod = 0;
    fEnergyAxisBins_log10 = 60;
    fIgnoreEnergyReconstructionQuality = false;
    
    fMCEnergy_min = -99.;
    fMCEnergy_max = -99.;
    fMCEnergy_index = 5.;
    
    fCutFileName = "";
    fGammaHadronCutSelector = -1;
    fDirectionCutSelector = -1;
    
    fAzimuthBins = true;
    fIsotropicArrivalDirections = false;
    
    fIgnoreFractionOfEvents = 0.;
    
    fTelescopeTypeCuts = false;
    
    fFillMCHistograms = false;
    
    fgetXoff_Yoff_afterCut = false;
    
    fCoreScatterMode = "";
    fCoreScatterRadius = 0.;
    
    fViewcone_min = -1.;
    fViewcone_max = -1.;
    
    fdatafile = "";
    fMCdatafile_tree = "";
    fMCdatafile_histo = "";
    fGammaHadronProbabilityFile = "";
    
    fze = 0.;
    fnoise = 0;
    fpedvar = 0.;
    fXoff  = 0.;
    fYoff  = 0.;
    
    fWobbleIsotropic = 0.;
    
    telconfig_ntel = 0;
    telconfig_arraycentre_X = 0.;
    telconfig_arraycentre_Y = 0.;
    telconfig_arraymax = 0.;
    
    fCREnergySpectrumFile = "";
    fCREnergySpectrumID = 0;
    fCREnergySpectrum = 0;

    fWriteEventdatatrees = false;
}


/*!

   read run parameters from text file

   all lines without a '*' in the first column are ignored

*/
bool VInstrumentResponseFunctionRunParameter::readRunParameterFromTextFile( string ifile )
{
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error opening run parameter file " << ifile << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    string temp;
    cout << endl;
    cout << "========================================" << endl;
    cout << "run parameter file (" << ifile << ")" << endl;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
        
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            cout << is_line << endl;
            is_stream >> temp;
            // ENERGYSPECTRUMINDEX <number of indexes> <minimum spectral index> <step size>
            if( temp == "ENERGYSPECTRUMINDEX" )
            {
                is_stream >> fNSpectralIndex;
                is_stream >> fSpectralIndexMin;
                is_stream >> fSpectralIndexStep;
            }
            // MONTECARLOENERGYRANGE <MC min energy [TeV]> <MC max energy [TeV]> <MC spectral index>
            else if( temp == "MONTECARLOENERGYRANGE" )
            {
                cout << "readInputFileList: setting user defined Monte Carlo energy values (shouldn't be necessary)" << endl;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fMCEnergy_min;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fMCEnergy_max;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fMCEnergy_index;
                }
            }
            // write event data tree with cuts 
            else if( temp == "WRITEEVENTDATATREE" )
            {
                    if( !(is_stream>>std::ws).eof() )
                    {
                            is_stream >> fWriteEventdatatrees;
                    }
            }
            // fill MC histograms
            else if( temp == "FILLMONTECARLOHISTOS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fFillMCHistograms;
                }
            }
            // get Tree with Xoff and Yoff (and derot) after cut
            else if( temp == "GETXOFFYOFFAFTERCUTS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fgetXoff_Yoff_afterCut;
                }
            }
            // energy reconstruction method
            else if( temp == "ENERGYRECONSTRUCTIONMETHOD" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fEnergyReconstructionMethod;
                    if( fEnergyReconstructionMethod == 0 && fReconstructionType == NOT_SET )
                    {
                        fReconstructionType = ENERGY_ER;
                    }
                }
            }
            else if( temp == "RECONSTRUCTIONTYPE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    string temp2;
                    is_stream >> temp2;
                    
                    if( temp2 == "GEO" || temp2 == "TL" )
                    {
                        fReconstructionType = GEO;
                    }
                    else if( temp2 == "NN" )
                    {
                        fReconstructionType = NN;
                    }
                    else if( temp2 == "ENERGY_ER" )
                    {
                        cout << "Warning: using energy reconstruction method 0" << endl;
                        fReconstructionType = ENERGY_ER;
                        fEnergyReconstructionMethod = 0;
                    }
                    else if( temp2 == "FROGS" )
                    {
                        fReconstructionType = FROGS;
                    }
                    else if( temp2 == "MODEL3D" )
                    {
                        fReconstructionType = MODEL3D;
                    }
                    else if( temp2 == "DEEPLEARNER" )
                    {
                        fReconstructionType = DEEPLEARNER;
                    }
                    else
                    {
                        fReconstructionType = NOT_SET;
                        cout << "VInstrumentResponseFunctionRunParameter::readRunParameter Warning: Unknown analysis type ";
                        cout << temp2 << ", will use GammaHadron cut file to decide." << endl;
                    }
                }
                else
                {
                    fReconstructionType = NOT_SET;
                    cout << "VInstrumentResponseFunctionRunParameter::readRunParameter Warning: Unknown analysis type,";
                    cout << " will use GammaHadron cut file to decide." << endl;
                }
                
            }
            
            // number of bins on log10 energy axis
            else if( temp == "ENERGYAXISBINS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fEnergyAxisBins_log10;
                }
            }
            // energy reconstruction quality
            else if( temp == "ENERGYRECONSTRUCTIONQUALITY" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fIgnoreEnergyReconstructionQuality;
                }
            }
            // number of azimuth bins
            else if( temp == "AZIMUTHBINS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fAzimuthBins;
                }
            }
            // isotropic arrival directions
            else if( temp == "ISOTROPICARRIVALDIRECTIONS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fIsotropicArrivalDirections;
                }
            }
            // ignore first fraction NN of events, scale MC distributions accordingly
            else if( temp == "IGNOREFRACTIONOFEVENTS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fIgnoreFractionOfEvents;
                }
            }
            // telescope type dependent cuts
            else if( temp == "TELESCOPETYPECUTS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fTelescopeTypeCuts;
                }
            }
            else if( temp == "CUTFILE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCutFileName;
                }
            }
            // * SCATTERMODE <core scatter radius [m]> <type of CORSIKA simulations (FLAT or VIEWCONE)>
            else if( temp == "SCATTERMODE" )
            {
                cout << "readInputFileList: setting user defined scatter values (shouldn't be necessary)" << endl;
                // read scatter radius
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCoreScatterRadius;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCoreScatterMode;
                }
            }
            else if( temp == "FILLINGMODE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fFillingMode;
                }
                if( fFillingMode > 3 )
                {
                    cout << "readInputFileList: error: invalid filling mode " << fFillingMode << endl;
                    return false;
                }
            }
            else if( temp == "SHORTOUTPUT" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fEffArea_short_writing;
                }
            }
            // read input data file name
            else if( temp == "SIMULATIONFILE_DATA" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fdatafile;
                }
            }
            // file with MC trees (set to "0" if not available)
            else if( temp == "SIMULATIONFILE_MC" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fMCdatafile_tree;
                }
            }
            // file with MC histograms
            else if( temp == "SIMULATIONFILE_HISTO" || temp == "SIMULATIONFILE_MCHISTO" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fMCdatafile_histo;
                }
            }
            // file with probabilities for gamma/hadron separation (should be a friend to the data tree)
            else if( temp == "GAMMAHADRONPROBABILITYFILE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fGammaHadronProbabilityFile;
                }
            }
            // name pointing to energy spectra used in event weighting
            else if( temp == "ESPECTRUM_FOR_WEIGHTING" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCREnergySpectrumFile;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fCREnergySpectrumID;
                }
            }
            // manually input the zenith
            else if( temp == "ZENITH" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fze;
                }
            }
            // manually input the zenith
            else if( temp == "NOISE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fnoise;
                }
            }
            // manually input the wobble
            else if( temp == "WOBBLEISOTROPIC" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> fWobbleIsotropic;
                }
            }
        }
    }
    cout << "========================================" << endl << endl;
    
    //////////////////////////////////////////////////////////////////////////////////////
    // fill some parameters
    // read run parameters from this file
    if( !readRunParameters( fdatafile ) )
    {
        cout << "VInstrumentResponseFunctionRunParameter::readRunParameterFromTextFile: Reading run parameter from ";
        cout << fdatafile << " failed. " << endl;
        return false;
    }
    // read spectral energy parameters
    if( !readCRSpectralParameters() )
    {
        return false;
    }
    
    /////////////////////////////////////////////////////////////////
    // define azimuth bins
    fAzMin.clear();
    fAzMax.clear();
    // 90 deg intervals
    if( fAzimuthBins )
    {
        fAzMin.push_back( 135.0 );
        fAzMax.push_back( -135. );
        fAzMin.push_back( -180. );
        fAzMax.push_back( -90. );
        for( unsigned int i = 0; i < 6; i++ )
        {
            fAzMin.push_back( fAzMin.back() + 45. );
            fAzMax.push_back( fAzMax.back() + 45. );
        }
        
        /*
         * 45 deg intervals
        fAzMin.push_back( 135.0 );
        	fAzMax.push_back( -165.0 );
        	fAzMin.push_back( 157.5 );
        	fAzMax.push_back( -142.5 );
        	fAzMin.push_back( -180. );
        	fAzMax.push_back( -120. );
        	for( int i = 0; i < 13; i++ )
        	{
        		fAzMin.push_back( fAzMin.back() + 22.5 );
        		fAzMax.push_back( fAzMax.back() + 22.5 );
        	} */
    }
    // (no az cut)
    fAzMin.push_back( -1.e3 );
    fAzMax.push_back( +1.e3 );
    // WARNING: if this last rule changes (if the last bin is NOT filled ANY MORE with all simulated event regardless of their az)
    //          then the az_bin_index must be changed in VEffectiveAreaCalculator::fill
    
    /////////////////////////////////////////////////////////////////
    // define  spectral index bins
    fSpectralIndex.clear();
    for( unsigned int i = 0; i < fNSpectralIndex; i++ )
    {
        fSpectralIndex.push_back( fSpectralIndexMin + ( double )i * fSpectralIndexStep );
    }
    
    return true;
}


VMonteCarloRunHeader* VInstrumentResponseFunctionRunParameter::readMCRunHeader()
{
    TChain c( "data" );
    TFile* iF = 0;
    if( c.Add( fdatafile.c_str() ) )
    {
        iF = c.GetFile();
    }
    
    if( !iF )
    {
        cout << "VInstrumentResponseFunctionRunParameter::readMCRunHeader: error opening data file: " << fdatafile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    VMonteCarloRunHeader* iMC = ( VMonteCarloRunHeader* )iF->Get( "MC_runheader" );
    if( !iMC )
    {
        cout << "VInstrumentResponseFunctionRunParameter::readMCRunHeader: no MC run header found in " << fdatafile << endl;
        cout << "(" << iMC->GetName() << ")" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    fCoreScatterRadius = iMC->core_range[0];
    // check if core scattering area is circular
    if( fCoreScatterRadius < 1.e-5 && iMC->core_pos_mode == 1 )
    {
        fCoreScatterRadius = iMC->core_range[1];
    }
    if( iMC->VOLUMEDET_set() )
    {
        fCoreScatterMode = "VIEWCONE";
    }
    else
    {
        fCoreScatterMode = "FLAT";
    }
    // get view cone
    fViewcone_min = iMC->viewcone[0];
    fViewcone_max = iMC->viewcone[1];
    if( fMCEnergy_min > -90. )
    {
        cout << "readMCRunheader: WARNING overridung user defined Monte Carlo energy range" << endl;
    }
    fMCEnergy_min = iMC->E_range[0];
    fMCEnergy_max = iMC->E_range[1];
    fMCEnergy_index = iMC->spectral_index;
    
    return iMC;
}

/*
 * read eventdisp run parameter from file
 *
 */
bool VInstrumentResponseFunctionRunParameter::readRunParameters( string ifilename )
{
    if( ifilename.find( "*" ) < ifilename.size() )
    {
        cout << "Information: using first file in list of files to determinate run parameters" << endl;
        cout << "(all files should have same energy range and spectral index)" << endl;
    }
    
    TChain c( "data" );
    TFile* iFile = 0;
    if( c.Add( ifilename.c_str() ) )
    {
        iFile = c.GetFile();
    }
    if( !iFile )
    {
        cout << "VInstrumentResponseFunctionRunParameter::readRunParameters() error reading simulation file: " << ifilename << endl;
        return false;
    }
    // read instrument epoch from run parameters
    VEvndispRunParameter* i_runPara = ( VEvndispRunParameter* )iFile->Get( "runparameterV2" );
    if( i_runPara )
    {
        fInstrumentEpoch = i_runPara->fInstrumentEpoch;
        fTelToAnalyse = i_runPara->fTelToAnalyze;
    }
    else
    {
        cout << "VInstrumentResponseFunctionRunParameter::readRunParameters() warning: ";
        cout << "cannot read instrument epoch and active telecopes from MC event file" << endl;
        cout << "this might lead to a wrong choice in the gamma/hadron cuts - please check" << endl;
        fInstrumentEpoch = "NOT_FOUND";
    }
    // get NSB (pedvar) level
    VTableLookupRunParameter* fR = ( VTableLookupRunParameter* )iFile->Get( "TLRunParameter" );
    if( !fR )
    {
        cout << "VInstrumentResponseFunctionRunParameter::readRunParameters() error: cannot find tablelookup run parameters in " << ifilename << endl;
        return false;
    }
    if( !fIsotropicArrivalDirections )
    {
        fze = fR->ze;
        fnoise = fR->fNoiseLevel;
    }
    fpedvar = fR->meanpedvars;
    // get wobble offset from first event in file
    // (should not change during a simulation run!)
    TTree* i_data = ( TTree* )iFile->Get( "data" );
    if( !i_data )
    {
        return false;
    }
    
    double x = 0.;
    double y = 0.;
    i_data->SetBranchAddress( "MCxoff", &x );
    i_data->SetBranchAddress( "MCyoff", &y );
    i_data->GetEntry( 0 );
    if( !fIsotropicArrivalDirections )
    {
        fXoff = x;
        fYoff = y;
    }
    else
    {
        if( fWobbleIsotropic != 0. )
        {
            fXoff = fWobbleIsotropic;
            fYoff = 0.;
        }
        else
        {
            fXoff = 0.;
            fYoff = 0.;
        }
    }
    
    ////////////////////////////////////////////////
    // read array parameters
    ////////////////////////////////////////
    // get number of telescopes from file
    TTree* t = ( TTree* )iFile->Get( "telconfig" );
    Ctelconfig* telconfig = new Ctelconfig( t );
    if( telconfig->IsZombie() )
    {
        cout << "error while reading telescope configuration" << endl;
        exit( EXIT_FAILURE );
    }
    telconfig_ntel = telconfig->getNTel();
    telconfig_arraycentre_X = telconfig->getArrayCentreX();
    telconfig_arraycentre_Y = telconfig->getArrayCentreY();
    telconfig_arraymax      = telconfig->getArrayMaxSize();
    
    ////////////////////////////////////////
    
    return true;
}

bool VInstrumentResponseFunctionRunParameter::testRunparameters()
{
    // check that there is at least one MC input file given
    if( fMCdatafile_tree.size() == 0 && fMCdatafile_histo.size() == 0 )
    {
        cout << "VInstrumentResponseFunctionRunParameter::testRunparameters() error:" << endl;
        cout << "   no MC file given" << endl;
        return false;
    }
    
    return true;
}


void VInstrumentResponseFunctionRunParameter::print()
{
    cout << endl;
    cout << "run parameters for calculation of instrument response functions: " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << endl;
    cout << "filling mode " << fFillingMode;
    if( fFillingMode == 0 )
    {
        cout << " (calculating complete set of response functions (effective areas, energy, core and angular resolution curves))";
    }
    if( fFillingMode == 1 )
    {
        cout << " (calculating core and angular resolution curves)";
    }
    if( fFillingMode == 2 )
    {
        cout << " (calculating angular resolution curves)";
    }
    if( fFillingMode == 3 )
    {
        cout << " (calculating effective areas and energy resolution curves)";
    }
    if( fEffArea_short_writing )
    {
        cout << "( short output )";
    }
    else
    {
        cout << "( long output )";
    }
    cout << endl;
    if( fFillMCHistograms )
    {
        cout << " filling MC histograms only" << endl << endl;
    }
    cout << endl;
    cout << "data files:" << endl;
    cout << "\t shower data:   " << fdatafile << endl;
    if( fMCdatafile_tree.size() > 0 )
    {
        cout <<  "  MC data:     " << fMCdatafile_tree << endl;
    }
    if( fMCdatafile_histo.size() > 0 )
    {
        cout << "  MC histograms:   " << fMCdatafile_histo << endl;
    }
    if( fGammaHadronProbabilityFile.size() > 0 )
    {
        cout << "  gamma/hadron probabilities: " << fGammaHadronProbabilityFile << endl;
    }
    if( fInstrumentEpoch != "NOT_SET" )
    {
        cout << "Instrument epoch: " << fInstrumentEpoch << endl;
    }
    else
    {
        cout << "Instrument epoch not set" << endl;
    }
    
    cout << endl;
    cout << "cuts: ";
    cout << "  " << fCutFileName << endl;
    cout << "cut selectors: ";
    if( fGammaHadronCutSelector >= 0 )
    {
        cout << "Gamma/Hadron " << fGammaHadronCutSelector;
        cout << " Direction " << fDirectionCutSelector;
    }
    else
    {
        cout << " (set later)";
    }
    if( fIgnoreEnergyReconstructionQuality )
    {
        cout << ", ignoring cut on quality of energy reconstruction";
    }
    if( fTelescopeTypeCuts )
    {
        cout << ", telescope type dependent cuts";
    }
    cout << endl;
    cout << "reconstruction type " << fReconstructionType << endl;
    cout << "energy reconstruction method " << fEnergyReconstructionMethod << endl;
    cout << endl;
    
    cout << "input Monte Carlo with following parameters (will be modified later): " << endl;
    cout << "\t core range: " << fCoreScatterRadius;
    if( fCoreScatterMode.size() > 0 )
    {
        cout << ", scatter mode " << fCoreScatterMode;
    }
    else
    {
        cout << ", (no scatter mode) ";
    }
    cout << ", energy range [TeV]: " << fMCEnergy_min << ", " << fMCEnergy_max << ", " << fMCEnergy_index << endl;
    cout << "\t ze=" << fze << " [deg], noise=" << fnoise << " (pedvar: " << fpedvar;
    cout << "), wobble offset w=" << sqrt( fXoff * fXoff + fYoff * fYoff ) << " [deg]";
    cout << endl;
    
    cout << "azimuth bins (" << fAzMin.size() << "): ";
    for( unsigned int i = 0; i < fAzMin.size(); i++ )
    {
        cout << " " << i << " [" << fAzMin[i] << "," << fAzMax[i] << "]";
    }
    cout << endl;
    
    cout << endl;
    cout << "array configuration: ";
    cout << telconfig_ntel;
    cout << " telescopes, array dimensions: centre (";
    cout << telconfig_arraycentre_X << "," << telconfig_arraycentre_Y << ")";
    cout << ", max telescope distance to centre: " << telconfig_arraymax << endl;
    cout << endl;
    cout << "calculate response function assuming power laws with following index (";
    cout << fNSpectralIndex << ", " << fSpectralIndexMin << ", " << fSpectralIndexStep << "): ";
    for( unsigned int i = 0; i < fSpectralIndex.size(); i++ )
    {
        cout << fSpectralIndex[i] << " ";
    }
    if( fCREnergySpectrumFile.size() > 0 )
    {
        cout << endl;
        cout << "CR energy spectrum used for weighted rate histogram: ";
        cout << fCREnergySpectrumFile << " (ID" << fCREnergySpectrumID << ")" << endl;
    }
    cout << endl << endl;
}

bool VInstrumentResponseFunctionRunParameter::readCRSpectralParameters()
{
    if( fCREnergySpectrumFile.size() == 0 )
    {
        return true;
    }
    
    VEnergySpectrumfromLiterature espec( fCREnergySpectrumFile );
    if( espec.isZombie() )
    {
        return false;
    }
    if( !espec.isValidID( fCREnergySpectrumID ) )
    {
        return false;
    }
    
    // read energy spectrum from literature class
    TF1* i_Ftemp = espec.getEnergySpectrum( fCREnergySpectrumID, false );
    if( i_Ftemp )
    {
        char hname[1000];
        sprintf( hname, "%s_C", i_Ftemp->GetName() );
        fCREnergySpectrum = new TF1( hname, i_Ftemp->GetExpFormula(),
                                     i_Ftemp->GetXmin(),
                                     i_Ftemp->GetXmax() );
        for( int i = 0; i < i_Ftemp->GetNpar(); i++ )
        {
            fCREnergySpectrum->SetParameter( i, i_Ftemp->GetParameter( i ) );
            fCREnergySpectrum->SetParError( i, i_Ftemp->GetParError( i ) );
        }
        delete i_Ftemp;
    }
    else
    {
        fCREnergySpectrum = 0;
    }
    
    return true;
}

