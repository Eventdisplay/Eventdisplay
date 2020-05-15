/*! \class VTableLookup
    \brief calculation of mean scaled variables and energies using MC filled tables

    note that VTableEnergyCalculator is not used anymore; just kept for historical reasons

*/

#include "VTableLookup.h"

/*!

*/
VTableLookup::VTableLookup( VTableLookupRunParameter* iTLRunParameter )
{
    fTLRunParameter = iTLRunParameter;
    if( !fTLRunParameter )
    {
        cout << "VTableLookup::initialize: error: no table lookup run parameters " << endl;
        exit( EXIT_FAILURE );
    }
    
    // total number of telescopes
    fNTel = 0;
    // look up table file
    fLookupTableFile = 0;
    
    fNumberOfIgnoredEvents = 0;
    fNNoiseLevelWarnings = 0;
    
    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    if( fTLRunParameter->fWriteTables )
    {
        cout << "filling lookup tables" << endl;
    }
    else
    {
        cout << "reading lookup tables" << endl;
    }
    cout << "-------------------------------------------------------" << endl;
    cout << endl;
    
    // maximum core distance for events to be taken into account
    fMeanNoiseLevel = 0.;
    
    // bins in azimuth
    // lookup tables are created for each bin
    fTableAzLowEdge.push_back( 135. );
    fTableAzUpEdge.push_back( -135. );
    fTableAzLowEdge.push_back( -135. );
    fTableAzUpEdge.push_back( -45. );
    fTableAzLowEdge.push_back( -45. );
    fTableAzUpEdge.push_back( 45. );
    fTableAzLowEdge.push_back( 45. );
    fTableAzUpEdge.push_back( 135 );
    // azimuth independent bins
    //    fTableAzLowEdge.push_back(  -1.e3 );   fTableAzUpEdge.push_back( 1.e3 );
    
    s_NupZupWup = 0;
    s_NupZupWlow = 0;
    s_NupZup = 0;
    s_NupZlowWup = 0;
    s_NupZlowWlow = 0;
    s_NupZlow = 0;
    s_Nup = 0;
    s_NlowZupWup = 0;
    s_NlowZupWlow = 0;
    s_NlowZup = 0;
    s_NlowZlowWup = 0;
    s_NlowZlowWlow = 0;
    s_NlowZlow = 0;
    s_Nlow = 0;
    s_N = 0;
    
    fTableCalculator = 0;
    fData = 0;
    
}


/*!
      \param ifile output file name
      \param ioption output file option ('recreate' or 'update')
*/
void VTableLookup::setOutputFile( string ifile, string ioption, string ioutputfile )
{
    if( fTLRunParameter->fWriteTables )
    {
        cout << "VTableLookup::setOutputFile warning: setting outputfiles makes no sense" << endl;
    }
    fData->setOutputFile( ifile, ioption, ioutputfile );
}

/*
 * initialize lookup data table vectors
 *
 */
void VTableLookup::initializeLookupTableDataVector()
{
    // table for MSCW
    fTableData[E_MSCW] = new VTableCalculatorData();
    fTableData[E_MSCW]->fDirectoryName = "mscw";
    fTableData[E_MSCW]->fFillVariable  = "width";
    fTableData[E_MSCW]->fEnergy        = false;
    
    // table for MSCL
    fTableData[E_MSCL] = new VTableCalculatorData();
    fTableData[E_MSCL]->fDirectoryName = "mscl";
    fTableData[E_MSCL]->fFillVariable  = "length";
    fTableData[E_MSCL]->fEnergy        = false;
    
    // table for energy reconstruction
    fTableData[E_EREC] = new VTableCalculatorData();
    fTableData[E_EREC]->fDirectoryName = "energySR";
    fTableData[E_EREC]->fFillVariable  = "energySR";
    fTableData[E_EREC]->fEnergy = true;
    
    // tables for time gradient lookup tables
    if( fTLRunParameter && fTLRunParameter->fUsetimeGradientLookupTables )
    {
        fTableData[E_TGRA] = new VTableCalculatorData();
        fTableData[E_TGRA]->fDirectoryName = "tgrad";
        fTableData[E_TGRA]->fFillVariable  = "tgrad";
        fTableData[E_TGRA]->fEnergy        = false;
        fTableData[E_TGRA]->fValueNormalizationRange_min = -25.;
        fTableData[E_TGRA]->fValueNormalizationRange_max =  25.;
    }
    
    // tables for frogs goodness of fit
    if( fTLRunParameter && fTLRunParameter->fUsefrogsGoodnessTables )
    {
        fTableData[E_FRGO] = new VTableCalculatorData();
        fTableData[E_FRGO]->fDirectoryName = "mscf";
        fTableData[E_FRGO]->fFillVariable  = "frogsGoodness";
        fTableData[E_FRGO]->fEnergy        = false;
        fTableData[E_FRGO]->fValueNormalizationRange_min = -1.;
        fTableData[E_FRGO]->fValueNormalizationRange_max =  20.;
    }
    
    // print table types
    map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
    
    for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); iter_iLT_Data++ )
    {
        // lookup table data pointer
        VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
        if( !iTableData )
        {
            continue;
        }
        iTableData->print();
    }
    cout << endl;
    
}


/*!
     construct calculators

     this function is called in table WRITING mode

     \param itablefile  file with all the tables (root file)
     \param ize         zenith angle of simulations to be filled
     \param woff        direction offset from camera center to be filled
     \param noise       noise level to be filled (might be different for each telescope)
     \param isuff       histogram suffix
*/
void VTableLookup::setMCTableFiles_forTableWriting( string itablefile, double ize, int woff, map< ULong64_t, double > noise,
        string isuff, string iFileTitle )
{
    if( !fTLRunParameter )
    {
        cout << "VTableLookup::setMCTableFiles error: no lookup table run parameters" << endl;
        return;
    }
    if( fTLRunParameter->fDebug )
    {
        cout << "void VTableLookup::setMCTableFiles" << endl;
    }
    if( !fTLRunParameter->fWriteTables )
    {
        cout << "VTableLookup::setMCTableFiles warning: setting mc table files makes no sense" << endl;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // create the table file
    fLookupTableFile = new TFile( itablefile.c_str(), "NEW", iFileTitle.c_str() );
    if( fLookupTableFile->IsZombie() )
    {
        cout << "VTableLookup::setMCTableFiles error while opening table file: " << itablefile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( fLookupTableFile->TestBit( TFile::kRecovered ) )
    {
        cout << "VTableLookup::setMCTableFiles problems with file (TFile::kRecovered, " << fLookupTableFile->TestBit( TFile::kRecovered ) << "): " << endl;
        cout << itablefile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    ////////////////////////////////////////
    // create lookup table data vector
    initializeLookupTableDataVector();
    
    
    /////////////////////////////////////////////////////////////////////////
    // prepare output file and create all directories
    
    if( fTLRunParameter->fDebug )
    {
        cout << "VTableLookup::setMCTableFiles() prepare output file and create all directories" << endl;
    }
    char hname[800];
    char htitle[800];
    
    /////////////////////////////////////////////////////////////////////////
    map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
    
    for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); iter_iLT_Data++ )
    {
        // lookup table data pointer
        VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
        if( !iTableData )
        {
            continue;
        }
        
        /////////////////////
        // temporary vectors
        
        // az
        vector< VTableCalculator* > i_LT;
        // direction ofset
        vector< vector< VTableCalculator* > > ii_LT;
        // ze
        vector< vector< vector< VTableCalculator* > > > iii_LT;
        // NSB
        vector< vector< vector< vector< VTableCalculator* > > > > iiii_LT;
        
        map<ULong64_t, unsigned int > i_list_of_Tel_type = fData->getList_of_Tel_type();
        map<ULong64_t, unsigned int >::iterator iter_i_list_of_Tel_type;
        
        //////////////////
        // TELESCOPE TYPES
        for( iter_i_list_of_Tel_type = i_list_of_Tel_type.begin(); 
                iter_i_list_of_Tel_type != i_list_of_Tel_type.end();
                ++iter_i_list_of_Tel_type )
        {
            // get telescope type
            ULong64_t t = iter_i_list_of_Tel_type->first;
            
            sprintf( hname, "tel_%lld", t );
            sprintf( htitle, "telescope type %lld", t );
            if( fLookupTableFile->Get( hname ) )
            {
                fLookupTableFile->cd( hname );
            }
            else
            {
                fLookupTableFile->mkdir( hname, htitle )->cd();
            }
            iiii_LT.clear();
            iii_LT.clear();
            
            /////////////
            // NOISE LEVEL
            if( noise.find( t ) != noise.end() )
            {
                int i_noise = ( int )( noise[t] * 100 );
                sprintf( hname, "NOISE_%05d", i_noise );
                if( gDirectory->Get( hname ) )
                {
                    gDirectory->cd( hname );
                }
                else
                {
                    gDirectory->mkdir( hname )->cd();
                }
                // ZENITH ANGLE
                sprintf( hname, "ze_%03d", ( int )( ize * 10. + 0.5 ) );
                if( gDirectory->Get( hname ) )
                {
                    gDirectory->cd( hname );
                }
                else
                {
                    gDirectory->mkdir( hname )->cd();
                }
                
                ii_LT.clear();
                
                //////////////////
                // DIRECTION OFFSET
                vector< int > i_woff_vector;
                if( fTLRunParameter->fCTA_MC_offaxisBin_min.size() == 0 )
                {
                    i_woff_vector.push_back( woff );
                }
                else
                {
                    for( unsigned int i = 0; i < fTLRunParameter->fCTA_MC_offaxisBin_min.size(); i++ )
                    {
                        if( i < fTLRunParameter->fCTA_MC_offaxisBin_max.size() )
                        {
                            i_woff_vector.push_back( ( int )( 0.5 * ( fTLRunParameter->fCTA_MC_offaxisBin_min[i]
                                                              + fTLRunParameter->fCTA_MC_offaxisBin_max[i] ) * 1000. + 0.5 ) );
                        }
                    }
                }
                TDirectory* i_curDir_w = gDirectory;
                for( unsigned int w = 0; w < i_woff_vector.size(); w++ )
                {
                    i_curDir_w->cd();
                    
                    sprintf( hname, "woff_%04d", i_woff_vector[w] );
                    if( gDirectory->Get( hname ) )
                    {
                        gDirectory->cd( hname );
                    }
                    else
                    {
                        gDirectory->mkdir( hname )->cd();
                    }
                    
                    i_LT.clear();
                    
                    //////////////////
                    // AZIMUTH ANGLE
                    // (lowest level in directory hierachy)
                    TDirectory* i_curDir = gDirectory;
                    for( unsigned int i = 0; i < fTableAzLowEdge.size(); i++ )
                    {
                        i_curDir->cd();
                        sprintf( hname, "az_%u", i );
                        sprintf( htitle, "%.1f < az < %.1f", fTableAzLowEdge[i], fTableAzUpEdge[i] );
                        if( gDirectory->Get( hname ) )
                        {
                            gDirectory->cd( hname );
                        }
                        else
                        {
                            gDirectory->mkdir( hname, htitle )->cd();
                        }
                        cout << "create tables in path " << gDirectory->GetPath() << endl;
                        
                        TDirectory* i_curAzDir = gDirectory;
                        
                        i_curAzDir->cd();
                        TDirectory* i_Dir = gDirectory->mkdir( iTableData->fDirectoryName.c_str() );
                        i_LT.push_back( new VTableCalculator( iTableData->fFillVariable.c_str(),
                                                              isuff.c_str(),
                                                              fTLRunParameter->fWriteTables,
                                                              i_Dir,
                                                              iTableData->fEnergy,
                                                              fTLRunParameter->fPE ) );
                                                              
                                                              
                        i_LT.back()->setNormalizeTableValues( iTableData->fValueNormalizationRange_min,
                                                              iTableData->fValueNormalizationRange_max );
                        i_LT.back()->setWrite1DHistograms( fTLRunParameter->fWrite1DHistograms );
                        i_LT.back()->setMinRequiredShowerPerBin( fTLRunParameter->fMinRequiredShowerPerBin );
                        // event selection cut is only set for energy lookup tables
                        if( iTableData->fEnergy )
                        {
                            i_LT.back()->setEventSelectionCut( fTLRunParameter->fEventSelectionCut_lossCutMax,
                                                                fTLRunParameter->fEventSelectionCut_distanceCutMax );
                        }
                    }   // az
                    ii_LT.push_back( i_LT );
                } // direction offset
                iii_LT.push_back( ii_LT );
                iiii_LT.push_back( iii_LT );
            } // NSB && zenith angle
            iTableData->fTable.push_back( iiii_LT );
        } // telescope type
    } // different lookup tables
}


/*!
    read in lookup tables and parameter space

    this function is called in table READING mode

    input file must be sorted in zenith angles (lowest first)

    \param itablefile   root file with all the tables
    \param isuff        histogram suffix
*/
void VTableLookup::setMCTableFiles_forTableReading( string itablefile, string isuff )
{
    if( fTLRunParameter->fDebug )
    {
        cout << "void VTableLookup::setMCTableFiles_forTableReading( string itablefile, string isuff )" << endl;
    }
    
    // open table file
    gErrorIgnoreLevel = 20001;
    fLookupTableFile = new TFile( itablefile.c_str() );
    if( fLookupTableFile->IsZombie() )
    {
        fLookupTableFile->Close();
        const char* data_dir = gSystem->Getenv( "OBS_EVNDISP_AUX_DIR" );
        if( data_dir )
        {
            // try to see of file exists in directory ./tables
            string itemp = data_dir;
            itemp += "/Tables/" + itablefile;
            itablefile = itemp;
            fLookupTableFile = new TFile( itablefile.c_str() );
            if( fLookupTableFile->IsZombie() )
            {
                cout << "VTableLookup::setMCTableFiles_forTableReading error (reading): unable to open table file: " << itablefile << endl;
                exit( EXIT_FAILURE );
            }
        }
        else
        {
            cout << "VTableLookup::setMCTableFiles_forTableReading error (reading): unable to open table file: " << itablefile << endl;
            cout << " (no $OBS_EVNDISP_AUX_DIR defined)" << endl;
            exit( EXIT_FAILURE );
        }
    }
    gErrorIgnoreLevel = 0;
    cout << "reading table file ( may take a while ): " << itablefile << endl;
    
    ////////////////////////////////////////
    // create lookup table data vector
    initializeLookupTableDataVector();
    
    // vector with available tel_types [tel_type]
    fTableTelTypes.clear();
    
    // vector with available NSB levels [tel_type][NSB]
    fTableNoiseLevel.clear();
    vector< double > i_NoiseLevel;
    
    // vector with available zenith angles [tel_type][NSB][ze]
    fTableZe.clear();
    vector< double > i_ze;
    vector< vector< double > > ii_ze;
    
    // vector with available wobble offsets [tel_type][NSB][ze][woff]
    fTableDirectionOffset.clear();
    vector< double > i_DirectionOffset;
    vector< vector< double > > ii_DirectionOffset;
    vector< vector< vector< double > > > iii_DirectionOffset;
    
    /////////////////////////////////////////////////////////////////////////
    // read in lookup tables
    // - everything is organized in sub directories
    
    /////////////////////////////////////////////////////////////////////////
    // loop over all lookup table types
    map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
    
    for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); iter_iLT_Data++ )
    {
        // lookup table data pointer
        VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
        if( !iTableData )
        {
            continue;
        }
        
        // temporary lookup tables vectors: [tel_type][NSB][ze][woff][az]
        // it is a mess....
        vector< VTableCalculator* > i_LT;
        vector< vector< VTableCalculator* > > ii_LT;
        vector< vector< vector< VTableCalculator* > > > iii_LT;
        vector< vector< vector< vector< VTableCalculator* > > > > iiii_LT;
        
        fTableTelTypes.clear();
        
        
        ////
        // TELESCOPE TYPE
        // (for VTS, each telescope is of a different type)
        vector< string > iDNameTel = getSortedListOfDirectories( fLookupTableFile );
        for( unsigned int t = 0; t < iDNameTel.size(); t++ )
        {
            // skip debug directories
            if( iDNameTel[t].find( "makeTable" ) != string::npos )
            {
                continue;
            }
            fLookupTableFile->cd( iDNameTel[t].c_str() );
            fTableTelTypes.push_back( ( ULong64_t )( atoi )( iDNameTel[t].substr( 4, iDNameTel[t].size() ).c_str() ) );
            
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  DIR TELTYPE " << " " << gDirectory->GetPath() << endl;
            }
            
            iiii_LT.clear();
            
            i_NoiseLevel.clear();
            ii_ze.clear();
            iii_DirectionOffset.clear();
            
            ////
            // NOISE LEVEL
            TDirectory* iDirNSB = gDirectory;
            vector< string > iDNameNSB = getSortedListOfDirectories( iDirNSB );
            for( unsigned int n = 0; n < iDNameNSB.size(); n++ )
            {
                i_NoiseLevel.push_back( atof( iDNameNSB[n].substr( 6, 5 ).c_str() ) / 100. );
                
                iDirNSB->cd( iDNameNSB[n].c_str() );
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "  DEBUG  DIR NSB " << " " << gDirectory->GetPath() << endl;
                }
                
                iii_LT.clear();
                
                i_ze.clear();
                ii_DirectionOffset.clear();
                
                ////
                // ZENITH ANGLE
                TDirectory* iDirZe = gDirectory;
                vector< string > iDNameZE = getSortedListOfDirectories( iDirZe );
                for( unsigned z = 0; z < iDNameZE.size(); z++ )
                {
                    i_ze.push_back( atof( iDNameZE[z].substr( 3, 3 ).c_str() ) / 10. );
                    
                    iDirZe->cd( iDNameZE[z].c_str() );
                    
                    if( fTLRunParameter->fDebug == 2 )
                    {
                        cout << "    DEBUG  DIR ZE " << " " << gDirectory->GetPath() << endl;
                    }
                    
                    ii_LT.clear();
                    
                    i_DirectionOffset.clear();
                    
                    ////
                    // DIRECTION OFFSET
                    TDirectory* iDirWoff = gDirectory;
                    vector< string > iDNameWoff  = getSortedListOfDirectories( iDirWoff );
                    for( unsigned int w = 0; w < iDNameWoff.size(); w++ )
                    {
                        i_DirectionOffset.push_back( atof( iDNameWoff[w].substr( 5, 4 ).c_str() ) / 1000. );
                        
                        iDirWoff->cd( iDNameWoff[w].c_str() );
                        
                        if( fTLRunParameter->fDebug == 2 )
                        {
                            cout << "      DEBUG  DIR WOFF " << " " << gDirectory->GetPath() << endl;
                        }
                        
                        i_LT.clear();
                        
                        // AZIMUTH ANGLE
                        TDirectory* iDirAz = gDirectory;
                        vector< string > iDNameAz  = getSortedListOfDirectories( iDirAz );
                        
                        for( unsigned int a = 0; a < iDNameAz.size(); a++ )
                        {
                            iDirAz->cd( iDNameAz[a].c_str() );
                            
                            if( fTLRunParameter->fDebug == 2 )
                            {
                                cout << "        DEBUG  DIR AZ " << " " << gDirectory->GetPath() << endl;
                            }
                            
                            // get lookup table directory
                            TDirectory* iDir = ( TDirectory* )gDirectory->Get( iTableData->fDirectoryName.c_str() );
                            // new lookup table calculator
                            i_LT.push_back( new VTableCalculator( iTableData->fFillVariable.c_str(),
                                                                  isuff.c_str(),
                                                                  fTLRunParameter->fWriteTables,
                                                                  iDir, iTableData->fEnergy,
                                                                  fTLRunParameter->fPE,
                                                                  fTLRunParameter->fUseMedianEnergy ) );
                            if( iTableData->fEnergy )
                            {
                                i_LT.back()->setEventSelectionCut( fTLRunParameter->fEventSelectionCut_lossCutMax,
                                                                   fTLRunParameter->fEventSelectionCut_distanceCutMax );
                            }
                            i_LT.back()->setNormalizeTableValues( iTableData->fValueNormalizationRange_min,
                                                                  iTableData->fValueNormalizationRange_max );
                            cout << gDirectory->GetPath() << endl;
                        }                             // az
                        ii_LT.push_back( i_LT );
                    }                                 // woff
                    iii_LT.push_back( ii_LT );
                    ii_DirectionOffset.push_back( i_DirectionOffset );
                }                                     // ze
                iiii_LT.push_back( iii_LT );
                
                ii_ze.push_back( i_ze );
                iii_DirectionOffset.push_back( ii_DirectionOffset );
            }                                         // NSB
            iTableData->fTable.push_back( iiii_LT );
            
            fTableNoiseLevel.push_back( i_NoiseLevel );
            fTableZe.push_back( ii_ze );
            fTableDirectionOffset.push_back( iii_DirectionOffset );
        }
    } // different telescope types
    
    // apply sanity check to the lookup table file
    if( sanityCheckLookupTableFile() )
    {
        cout << "    ...survived test of table file! " << endl;
    }
    else
    {
        cout << "     ...did not survive test of table file ! There are missing tables in your table file: " << itablefile << endl;
        cout << "     You need to redo your table, or contact Gernot if this is the standard table supplied " << endl ;
        sanityCheckLookupTableFile( true ); // this will print which tables are missing
        exit( EXIT_FAILURE );
    }
    
    if( fTLRunParameter->fDebug )
    {
        cout << "END void VTableLookup::setMCTableFiles_forTableReading( string itablefile, string isuff )" << endl;
    }
}

/*

     sanity checks for lookup tables

     - make sure that the same number of zenith, az, NSB tables are given for each telescope type

*/
bool VTableLookup::sanityCheckLookupTableFile( bool iPrint )
{
    if( iPrint == false )
    {
        cout << "lookup table file sanity check..." << endl;
    }
    if( iPrint == true )
    {
        cout << endl;
        cout << "The following printout should help you to locate the missing lookup tables..." << endl;
    }
    
    for( unsigned int t = 0; t < fTableTelTypes.size(); t++ )
    {
        if( iPrint == true )
        {
            cout << "\t Telescope type index " << t << " (" << fTableTelTypes[t] << ")";
            cout << " has " << fTableNoiseLevel[t].size() << " NSB levels" << endl;
        }
        
        // zenith angles: should be the same for each telescope type
        // noise level: allowed to be different for each telescope type
        for( unsigned int i = 0; i < fTableNoiseLevel[t].size(); i++ )
        {
            if( iPrint == true )
            {
                cout << "\tNoise level index  " << i << " has " << fTableZe[t][i].size() << " zenith angles" << endl;
            }
            else
            {
                if( ( i > 0 ) && ( fTableZe[t][i] != fTableZe[t][i - 1] ) )
                {
                    if( iPrint == true )
                    {
                        cout << "\tDifferent zenith angles bins: " << endl;
                    }
                    return false;
                }
            }
            
            // direction offsets: should be the same for each telescope type
            for( unsigned int j = 0; j < fTableZe[t][i].size(); j++ )
            {
                if( iPrint == true )
                {
                    cout << "\tZenith index " << j <<  " has " << fTableDirectionOffset[t][i][j].size() << " wobble offsets" << endl;
                }
                else
                {
                    if( ( j > 0 ) && ( fTableDirectionOffset[t][i][j] != fTableDirectionOffset[t][i][j - 1] ) )
                    {
                        if( iPrint == true )
                        {
                            cout << "\tDifferent direction offset bins: " << endl;
                        }
                        return false;
                    }
                }
            }
        }
    }
    
    // print a summary of the number of tables found
    if( iPrint == false )
    {
        cout << "Found " << fTableTelTypes.size() << " telescope types, ";
        if( fTableNoiseLevel.size() > 0 )
        {
            cout << fTableNoiseLevel[0].size() << " noise levels, ";
            if( fTableZe[0].size() > 0 )
            {
                cout << fTableZe[0][0].size() << " zenith angles, ";
                if( fTableDirectionOffset[0][0].size() > 0 )
                {
                    cout << fTableDirectionOffset[0][0][0].size() << " wobble offsets, ";
                    cout << fTableAzLowEdge.size() << " azimuth bins" << endl;
                }
            }
        }
        else
        {
            cout << endl;
            cout << "ERROR: no lookup tables found" << endl;
            exit( EXIT_FAILURE );
        }
    }
    return true;
}



/*!
     loop over all events, calculate mean scaled width, interpolate between zenith angles
*/
void VTableLookup::loop()
{
    if( fTLRunParameter->fWriteTables )
    {
        fillLookupTable();
    }
    else
    {
        readLookupTable();
    }
}

/*
    fill the tables

    (fTableData is a vector of size 1)
*/
void VTableLookup::fillLookupTable()
{
    double idummy1[fData->getMaxNbrTel()];
    double iEventWeight = 0.;
    double idummy3 = 0.;
    int fevent = 0;
    // telescope types
    map<ULong64_t, unsigned int> i_list_of_Tel_type = fData->getList_of_Tel_type();
    map<ULong64_t, unsigned int>::iterator iter_i_list_of_Tel_type;
    //////////////////////////////////////////////////////////////////////////////////////
    cout << "start event loop " << endl;
    // read next event
    while( fData->getNextEvent( true ) )
    {
        fevent = fData->getEventCounter();
        
        // print progress
        if( ( fevent % 1000000 ) == 0 && fevent != 0 )
        {
            cout << "\t now at event " << fevent << endl;
        }
        if( fTLRunParameter->fDebug )
        {
            cout << "now at event " << fevent << "\t" << fData->getEventStatus() << endl;
        }
        
        // get event weight (e.g. for spectral weighting)
        iEventWeight = fData->getEventWeight();
        
        // apply cuts
        if( fData->getEventStatus() && iEventWeight > 0. )
        {
            // get wobble bin (CTA only)
            unsigned int w = getWobbleBin( fData->getMCWobbleOffset() );
            // check vector sizes
            // (getWobbleBin() returns 999 if value is out of range)
            // assume here that MSCW table exists
            // (checking this would be too slow)
            if( !fTableData[E_MSCW]->assertTableVector( w ) )
            {
                continue;
            }
            
            // loop over all telescopes types, fill according to its type
            unsigned int i_Tel_type_counter = 0;
            for( iter_i_list_of_Tel_type = i_list_of_Tel_type.begin();
                    iter_i_list_of_Tel_type != i_list_of_Tel_type.end();
                    ++iter_i_list_of_Tel_type )
            {
                // get telescope type
                ULong64_t t = iter_i_list_of_Tel_type->first;
                
                // arrays for size and distance
                double* i_s2        = fData->getSize2( 1., t, fTLRunParameter->fUseSelectedImagesOnly );
                double* i_r         = fData->getDistanceToCore( t );
                double* i_l         = fData->getLoss( t );
                double* i_d         = fData->getDistance( t );
                // number of telescope of this particular type
                unsigned int iN_type = fData->getNTel_type( t );
                
                ////////////////////////////////////////////////
                // for zenith-angle == 0 deg fill all az bins
                if( fabs( fData->getMCZe() ) < 3. )
                {
                    // OBS: this provides fTableAzLowEdge.size() times more events in the tables for this zenith angle bin
                    for( unsigned int a = 0; a < fTableAzLowEdge.size(); a++ )
                    {
                        // fill tables (get arrays filled for corresponding telescope type; one table per type)
                        // assume here for efficiency that MSCW, MSCL and energy table always exists
                        // (yes, this is dangerous)
                        fTableData[E_MSCW]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                            iN_type, i_r, i_s2,
                            i_l, i_d, fData->getWidth( t ),
                            idummy1, iEventWeight, idummy3, idummy1 );
                        fTableData[E_MSCL]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                            iN_type, i_r, i_s2,
                            i_l, i_d, fData->getLength( t ),
                            idummy1, iEventWeight, idummy3, idummy1 );
                        fTableData[E_EREC]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                            iN_type, i_r, i_s2,
                            i_l, i_d, fData->getMCEnergyArray(),
                            idummy1, iEventWeight, idummy3, idummy1 );
                        // optional table: time gradient
                        if( fTableData.find( E_TGRA ) != fTableData.end() )
                        {
                            fTableData[E_TGRA]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                                iN_type, i_r, i_s2,
                                i_l, i_d, fData->getTimeGradient( t ),
                                idummy1, iEventWeight, idummy3, idummy1 );
                        }
                        // optional table: frogs goodness of fit
                        if( fTableData.find( E_FRGO ) != fTableData.end() )
                        {
                            fTableData[E_FRGO]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                                iN_type, i_r, i_s2,
                                i_l, i_d, fData->getFROGS_goodness( t ),
                                idummy1, iEventWeight, idummy3, idummy1 );
                        }
                    }
                }
                ////////////////////////////////////////////////
                // for zenith-angle != 0 deg get az bin
                else
                {
                    unsigned int a = getAzBin( fData->getMCAz() );
                    // fill tables (get arrays filled for corresponding telescope type; one table per type)
                    // assume here for efficiency that MSCW, MSCL and energy table always exists
                    // (yes, this is dangerous)
                    fTableData[E_MSCW]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                        iN_type, i_r, i_s2,
                        i_l, i_d, fData->getWidth( t ),
                        idummy1, iEventWeight, idummy3, idummy1 );
                    fTableData[E_MSCL]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                        iN_type, i_r, i_s2,
                        i_l, i_d, fData->getLength( t ),
                        idummy1, iEventWeight, idummy3, idummy1 );
                    fTableData[E_EREC]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                        iN_type, i_r, i_s2,
                        i_l, i_d, fData->getMCEnergyArray(),
                        idummy1, iEventWeight, idummy3, idummy1 );
                    // optional table: time gradient
                    if( fTableData.find( E_TGRA ) != fTableData.end() )
                    {
                        fTableData[E_TGRA]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                            iN_type, i_r, i_s2,
                            i_l, i_d, fData->getTimeGradient( t ),
                            idummy1, iEventWeight, idummy3, idummy1 );
                    }
                    // optional table: frogs goodness of fit
                    if( fTableData.find( E_FRGO ) != fTableData.end() )
                    {
                        fTableData[E_FRGO]->fTable[i_Tel_type_counter][0][0][w][a]->calc(
                            iN_type, i_r, i_s2,
                            i_l, i_d, fData->getFROGS_goodness( t ),
                            idummy1, iEventWeight, idummy3, idummy1 );
                    }
                }
                i_Tel_type_counter++;
            }
        }
    }
}

/*

   read the tables

*/
void VTableLookup::readLookupTable()
{
    int i_az = 0;
    double ze = 0.;
    double woff = 0.;
    int fevent = 0;
    double imr = 0.;
    double inr = 0.;
    
    // lookup table index for interpolation
    unsigned int inoise_up = 0;
    unsigned int inoise_low = 0;
    unsigned int ize_up = 0;
    unsigned int ize_low = 0;
    unsigned int iwoff_up = 0;
    unsigned int iwoff_low = 0;
    
    // (this is a bit of a mess)
    s_NupZupWup    = new VTablesToRead( fTableData.size(), fNTel );
    s_NupZupWlow   = new VTablesToRead( fTableData.size(), fNTel );
    s_NupZup       = new VTablesToRead( fTableData.size(), fNTel );
    s_NupZlowWup   = new VTablesToRead( fTableData.size(), fNTel );
    s_NupZlowWlow  = new VTablesToRead( fTableData.size(), fNTel );
    s_NupZlow      = new VTablesToRead( fTableData.size(), fNTel );
    s_Nup          = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZupWup   = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZupWlow  = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZup      = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZlowWup  = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZlowWlow = new VTablesToRead( fTableData.size(), fNTel );
    s_NlowZlow     = new VTablesToRead( fTableData.size(), fNTel );
    s_Nlow         = new VTablesToRead( fTableData.size(), fNTel );
    s_N            = new VTablesToRead( fTableData.size(), fNTel );
    
    // first event
    bool bFirst = true;
    if( !fTLRunParameter->isMC )
    {
        bFirst = false;
    }
    
    ////////////////////////////////////////////////
    // start event loop
    while( fData->getNextEvent( false ) )
    {
        // print progress
        fevent = fData->getEventCounter();
        if( ( fevent % 1000000 ) == 0 )
        {
            cout << "\t now at event " << fevent << endl;
        }
        
        // eventdisplay is saying that his event should be ignored
        if( fData->getEventNumber() == 99999999 )
        {
            fNumberOfIgnoredEvents++;
            continue;
        }
        
        // get zenith angle for first valid MC event from MC files
        if( bFirst && fData->getMCEnergy() > 0.001
          && fTLRunParameter->ze < 0. )
        {
            cout << "\t\t setting IRF ze from first event" << endl;
            if( fNTel > 0 )
            {
                fTLRunParameter->ze = TMath::Floor( ( 90. - fData->getTelElevation() ) + 0.5 );
            }
            else
            {
                fTLRunParameter->ze = TMath::Floor( fData->getMCZe() + 0.5 );
            }
            if( fTLRunParameter->ze < 1.5 )
            {
                fTLRunParameter->ze = 0.;
            }
            // not clear if this is needed? Note that MC distance to camera center might change
            // from event to event
            fTLRunParameter->fWobbleOffset = ( int )( fData->getMCWobbleOffset() * 100. );
            bFirst = false;
        }
        // reset image counter
        int fnmscw = 0;
        // fill MC energy spectra
        fData->fillMChistograms();
        // if data fails basic cuts, write default values directly to tree
        if( !fData->cut() )
        {
            if( fTLRunParameter->bNoNoTrigger )
            {
                fData->reset();
                fData->fill();
            }
            // goto next event
        }
        else
        {
            //////////////////////////////////////
            // here we should have good data only
            // (ze, az, and wobble offset have been
            //  tested)
            //////////////////////////////////////
            
            // get direction angles for this event
            ze   = fData->getZe();
            woff = fData->getWobbleOffset();
            i_az = getAzBin( fData->getAz() );
            // get noise level for this event
            readNoiseLevel( false );
            
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << endl << endl << "DEBUG  NEW EVENT " << fData->getEventCounter() << endl;
            }
            /////////////////////////////
            // interpolation section
            // interpolate for given size, R between NSB, ze, distance to camera center
            // no interpolation for AZ bins
            
            
            /////////////////////////////
            // NOISE (low) ZENITH (low)
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  NOISE LOW ZENITH LOW" << endl;
            }
            for( int t = 0; t < fNTel; t++ )
            {
                if( fTLRunParameter
                        && t < ( int )fTLRunParameter->fTelToAnalyzeData.size()
                        && !fTLRunParameter->fTelToAnalyzeData[t]->fTelToAnalyze )
                {
                    continue;
                }
                // index for this telescope type
                unsigned int telX = getTelTypeCounter( t, true );
                
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "DEBUG  TELESCOPE " << t << " (T" << t + 1 << "), teltype counter " << telX << endl;
                    cout << "DEBUG      zenith " << ze << ", noise " << fNoiseLevel[t];
                    cout << ", woff " << woff << ", az " << fData->getAz() << ", az bin " << i_az << endl;
                }
                // noise (low)
                getIndexBoundary( &inoise_up, &inoise_low, fTableNoiseLevel[telX], fNoiseLevel[t] );
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "DEBUG  NOISE " << t << " " << inoise_low << " " << inoise_up << " ";
                    cout << fNoiseLevel[t] << "\t" << fTableNoiseLevel[telX].size();
                    cout << ", NSB level: ";
                    for( unsigned int ii = 0; ii < fTableNoiseLevel[telX].size(); ii++ )
                    {
                        cout << fTableNoiseLevel[telX][ii] << ", ";
                    }
                    cout << endl;
                }
                // get zenith angle
                getIndexBoundary( &ize_up, &ize_low, fTableZe[telX][inoise_low], ze );
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "DEBUG  ZENITH " << t << " " << ize_up << " " << ize_low << " ";
                    cout << ze << "\t" << fTableZe[telX][inoise_low].size();
                    cout << ", ze bins: ";
                    for( unsigned int ii = 0; ii < fTableZe[telX][inoise_low].size(); ii++ )
                    {
                        cout << fTableZe[telX][inoise_low][ii] << ", ";
                    }
                    cout << endl;
                }
                
                // get direction offset index
                getIndexBoundary( &iwoff_up, &iwoff_low, fTableDirectionOffset[telX][inoise_low][ize_low], woff );
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "DEBUG  WOFF " << t << " " << iwoff_up << " " << iwoff_low << " " << woff << "\t";
                    cout << fTableDirectionOffset[telX][inoise_low][ize_low].size();
                    cout << ", woff bins: ";
                    for( unsigned int ii = 0; ii < fTableDirectionOffset[telX][inoise_low][ize_low].size(); ii++ )
                    {
                        cout << fTableDirectionOffset[telX][inoise_low][ize_low][ii] << ", ";
                    }
                    cout << endl;
                }
                
                // get tables
                getTables( inoise_low, ize_low, iwoff_up, i_az, t, s_NlowZlowWup );
                getTables( inoise_low, ize_low, iwoff_low, i_az, t, s_NlowZlowWlow );
                
            } // loop over all telescopes
            calculateMSFromTables( s_NlowZlowWup );
            calculateMSFromTables( s_NlowZlowWlow );
            // interpolate between direction offset bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NlowZlowWlow, fTableDirectionOffset[0][0][ize_low][iwoff_low],
                         s_NlowZlowWup,  fTableDirectionOffset[0][0][ize_low][iwoff_up],
                         s_NlowZlow, woff );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  WOFF INTER 1 ";
                cout << woff << " " << fTableDirectionOffset[0][0][ize_low][iwoff_low] << " ";
                cout << fTableDirectionOffset[0][0][ize_low][iwoff_up];
                cout << " " << ize_low << " ";
                cout << s_NlowZlowWlow->value[E_MSCL] << " ";
                cout << s_NlowZlowWup->value[E_MSCL] << " ";
                cout << s_NlowZlow->value[E_MSCL] << endl;
            }
            
            ///////////////////////////
            // NOISE (low) ZENITH (up)
            for( int t = 0; t < fNTel; t++ )
            {
                if( fTLRunParameter
                        && t < ( int )fTLRunParameter->fTelToAnalyzeData.size()
                        && !fTLRunParameter->fTelToAnalyzeData[t]->fTelToAnalyze )
                {
                    continue;
                }
                // index for this telescope type
                unsigned int telX = getTelTypeCounter( t, true );
                
                // noise (low)
                getIndexBoundary( &inoise_up, &inoise_low, fTableNoiseLevel[telX], fNoiseLevel[t] );
                // get zenith angle
                getIndexBoundary( &ize_up, &ize_low, fTableZe[telX][inoise_low], ze );
                
                // zenith angle (up)
                // get direction offset index
                getIndexBoundary( &iwoff_up, &iwoff_low, fTableDirectionOffset[telX][inoise_low][ize_up], woff );
                getTables( inoise_low, ize_up, iwoff_up, i_az, t, s_NlowZupWup );
                getTables( inoise_low, ize_up, iwoff_low, i_az, t, s_NlowZupWlow );
            }
            calculateMSFromTables( s_NlowZupWup );
            calculateMSFromTables( s_NlowZupWlow );
            
            // interpolate between direction offset bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NlowZupWlow, fTableDirectionOffset[0][0][ize_up][iwoff_low],
                         s_NlowZupWup,  fTableDirectionOffset[0][0][ize_up][iwoff_up],
                         s_NlowZup, woff );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  WOFF INTER 2 ";
                cout << woff << " " << fTableDirectionOffset[0][0][ize_up][iwoff_low] << " ";
                cout << fTableDirectionOffset[0][0][ize_up][iwoff_up] << " " << ize_up;
                cout << " " << s_NlowZupWlow->value[E_MSCL] << " ";
                cout << s_NlowZupWup->value[E_MSCL] << " ";
                cout << s_NlowZup->value[E_MSCL] << endl;
            }
            
            // interpolate between zenith angle bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NlowZlow, fTableZe[0][0][ize_low],
                         s_NlowZup, fTableZe[0][0][ize_up],
                         s_Nlow, ze, true );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  ZE INTER 1 " << ze << " " << fTableZe[0][0][ize_low] << " ";
                cout << fTableZe[0][0][ize_up] << " ";
                cout << s_NlowZlow->value[E_MSCL] << " ";
                cout << s_NlowZup->value[E_MSCL] << " ";
                cout << s_Nlow->value[E_MSCL] << endl;
            }
            
            ///////////////////////////
            // NOISE (up) ZENITH (low)
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  HIGH NOISE" << endl;
            }
            for( int t = 0; t < fNTel; t++ )
            {
                if( fTLRunParameter
                        && t < ( int )fTLRunParameter->fTelToAnalyzeData.size()
                        && !fTLRunParameter->fTelToAnalyzeData[t]->fTelToAnalyze )
                {
                    continue;
                }
                // index for this telescope type
                unsigned int telX = getTelTypeCounter( t, true );
                
                // noise (up)
                getIndexBoundary( &inoise_up, &inoise_low, fTableNoiseLevel[telX], fNoiseLevel[t] );
                // get zenith angle
                getIndexBoundary( &ize_up, &ize_low, fTableZe[telX][inoise_up], ze );
                if( fTLRunParameter->fDebug == 2 )
                {
                    cout << "DEBUG  WOFF " << t << " " << inoise_low << " " << inoise_up << " " << fNoiseLevel[t] << endl;
                }
                
                // zenith angle (low)
                // get direction offset index
                getIndexBoundary( &iwoff_up, &iwoff_low, fTableDirectionOffset[telX][inoise_up][ize_low], woff );
                getTables( inoise_up, ize_low, iwoff_up, i_az, t, s_NupZlowWup );
                getTables( inoise_up, ize_low, iwoff_low, i_az, t, s_NupZlowWlow );
            }
            calculateMSFromTables( s_NupZlowWup );
            calculateMSFromTables( s_NupZlowWlow );
            // interpolate between direction offset bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NupZlowWlow, fTableDirectionOffset[0][0][ize_low][iwoff_low],
                         s_NupZlowWup,  fTableDirectionOffset[0][0][ize_low][iwoff_up],
                         s_NupZlow, woff );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  WOFF INTER 1 ";
                cout << woff << " " << fTableDirectionOffset[0][0][ize_low][iwoff_low] << " ";
                cout << fTableDirectionOffset[0][0][ize_low][iwoff_up];
                cout << " " << s_NupZlowWlow->value[E_MSCL] << " ";
                cout << s_NupZlowWup->value[E_MSCL] << " ";
                cout << s_NupZlow->value[E_MSCL] << endl;
            }
            
            ///////////////////////////
            // NOISE (up) ZENITH (up)
            for( int t = 0; t < fNTel; t++ )
            {
                if( fTLRunParameter
                        && t < ( int )fTLRunParameter->fTelToAnalyzeData.size()
                        && !fTLRunParameter->fTelToAnalyzeData[t]->fTelToAnalyze )
                {
                    continue;
                }
                // index for this telescope type
                unsigned int telX = getTelTypeCounter( t, true );
                
                // noise (up)
                getIndexBoundary( &inoise_up, &inoise_low, fTableNoiseLevel[telX], fNoiseLevel[t] );
                // get zenith angle
                getIndexBoundary( &ize_up, &ize_low, fTableZe[telX][inoise_up], ze );
                
                // zenith angle (up)
                // get direction offset index
                getIndexBoundary( &iwoff_up, &iwoff_low, fTableDirectionOffset[telX][inoise_up][ize_up], woff );
                getTables( inoise_up, ize_up, iwoff_up, i_az, t, s_NupZupWup );
                getTables( inoise_up, ize_up, iwoff_low, i_az, t, s_NupZupWlow );
            }
            calculateMSFromTables( s_NupZupWup );
            calculateMSFromTables( s_NupZupWlow );
            
            // interpolate between direction offset bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NupZupWlow, fTableDirectionOffset[0][0][ize_up][iwoff_low],
                         s_NupZupWup,  fTableDirectionOffset[0][0][ize_up][iwoff_up],
                         s_NupZup, woff );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  WOFF INTER 2 ";
                cout << woff << " " << fTableDirectionOffset[0][0][ize_up][iwoff_low] << " ";
                cout << fTableDirectionOffset[0][0][ize_up][iwoff_up];
                cout << " " << s_NupZupWlow->value[E_MSCL] << " ";
                cout << s_NupZupWup->value[E_MSCL] << " ";
                cout << s_NupZup->value[E_MSCL] << endl;
            }
            
            // interpolate between zenith angle bins
            // note: expect that there are the same number of lookup tables for each
            //       telescope type and noise level
            interpolate( s_NupZlow, fTableZe[0][0][ize_low],
                         s_NupZup,  fTableZe[0][0][ize_up],
                         s_Nup, ze, true );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  ZE INTER 2 " << ze;
                cout << " " << s_NupZlow->value[E_MSCL] << " ";
                cout << s_NupZup->value[E_MSCL] << " ";
                cout << s_Nup->value[E_MSCL] << endl;
            }
            
            // calculate average table NSB for this type of telescopes
            double i_meanNoiseLevel_low = 0.;
            double i_meanNoiseLevel_up = 0.;
            double i_meanNoiseLevel_data = 0.;
            double i_meanNoiseLevel_N = 0.;
            for( int t = 0; t < fNTel; t++ )
            {
                if( fTelToAnalyze[t] )
                {
                    // index for this telescope type
                    unsigned int telX = getTelTypeCounter( t, true );
                    
                    // get noise level index
                    getIndexBoundary( &inoise_up, &inoise_low, fTableNoiseLevel[telX], fNoiseLevel[t] );
                    
                    i_meanNoiseLevel_low  += fTableNoiseLevel[telX][inoise_low];
                    i_meanNoiseLevel_up   += fTableNoiseLevel[telX][inoise_up];
                    i_meanNoiseLevel_data += fNoiseLevel[t];
                    i_meanNoiseLevel_N++;
                }
            }
            if( i_meanNoiseLevel_N > 0. )
            {
                i_meanNoiseLevel_low  /= i_meanNoiseLevel_N;
                i_meanNoiseLevel_up   /= i_meanNoiseLevel_N;
                i_meanNoiseLevel_data /= i_meanNoiseLevel_N;
            }
            
            // interpolate between NSB level bins
            interpolate( s_Nlow, i_meanNoiseLevel_low,
                         s_Nup,  i_meanNoiseLevel_up,
                         s_N,    i_meanNoiseLevel_data, false );
                         
            if( fTLRunParameter->fDebug == 2 )
            {
                cout << "DEBUG  NOISE INTER " << i_meanNoiseLevel_data << " ";
                cout << i_meanNoiseLevel_low << "\t" << i_meanNoiseLevel_up;
                cout << " " << s_Nlow->value[E_MSCL] << " ";
                cout << s_Nup->value[E_MSCL] << " ";
                cout << s_N->value[E_MSCL] << endl;
            }
            
            // (end of interpolation section)
            /////////////////////////////////
            
            
            //////////////////////////////////////////////////////
            // determine number of telescopes with MSCW values
            for( unsigned int j = 0; j < s_N->fNTel; j++ )
            {
                if( s_N->value_T[E_MSCW][j] > -90. )
                {
                    fnmscw++;
                }
            }
            fData->setNMSCW( fnmscw );
            
            /////////////////////////
            // set msc value (mean reduced scaled variables)
            fData->setMSCW( s_N->value[E_MSCW] );
            fData->setMSCL( s_N->value[E_MSCL] );
            
            // calculate mean width ratio (mean scaled variables)
            imr = 0.;
            inr = 0.;
            // require size2 > 0 (to use only selected images for the MWR/MWL calculation)
            double* i_s2 = fData->getSize2( 1., fTLRunParameter->fUseSelectedImagesOnly );
            for( unsigned int j = 0; j < s_N->fNTel; j++ )
            {
                if( s_N->value_T[E_MSCW][j] > 0. && fData->getWidth() && i_s2 && i_s2[j] > 0. )
                {
                    imr += fData->getWidth()[j] / s_N->value_T[E_MSCW][j];
                    inr++;
                }
            }
            if( inr > 0. )
            {
                fData->setMWR( imr / inr );
            }
            else
            {
                fData->setMWR( -99. );
            }
            // calculate mean length ratio (mean scaled variables)
            imr = 0.;
            inr = 0.;
            for( unsigned int j = 0; j < s_N->fNTel; j++ )
            {
                if( s_N->value_T[E_MSCL][j] > 0. && fData->getLength() && i_s2 && i_s2[j] > 0. )
                {
                    imr += fData->getLength()[j] / s_N->value_T[E_MSCL][j];
                    inr++;
                }
            }
            if( inr > 0. )
            {
                fData->setMLR( imr / inr );
            }
            else
            {
                fData->setMLR( -99. );
            }
            
            /////////////////////////
            // fill energies
            //
            // for dispEnergy: fill
            // energy in the data handler
            int iNERS = 0;
            if( !fData->useDispEnergy() )
            {
                fData->setEnergy( s_N->value[E_EREC],
                                  s_N->value_Chi2[E_EREC],
                                  s_N->value_dE[E_EREC] );
                // set energies per telescope
                for( unsigned int j = 0; j < s_N->fNTel; j++ )
                {
                    fData->setEnergyT( j, s_N->value_T[E_EREC][j] );
                    if( s_N->value_T[E_EREC][j] > 0. )
                    {
                        iNERS++;
                    }
                }
            }
            fData->setNEnergyT( iNERS );
            // energy quality
            if( iNERS > 0 )
            {
                // good energy from 1 or more images
                if( s_N->value_Chi2[E_EREC] >= 0. )
                {
                    if( iNERS > 1 )
                    {
                        fData->setNEnergyQuality( 0 );
                    }
                    else
                    {
                        fData->setNEnergyQuality( 1 );
                    }
                }
                // no reconstructed total energy, 
                else
                {
                    fData->setNEnergyQuality( -2 );
                }
            }
            else
            {
                fData->setNEnergyQuality( -1 );
            }
            // set mean reduced scaled widths and energies per telescope
            for( unsigned int j = 0; j < s_N->fNTel; j++ )
            {
                fData->setMSCWT( j, s_N->value_T[E_MSCW][j], s_N->value_T_sigma[E_MSCW][j] );
                fData->setMSCLT( j, s_N->value_T[E_MSCL][j], s_N->value_T_sigma[E_MSCL][j] );
            }
            
            /////////////////////////
            // fill mean scaled time gradient (optional)
            if( s_N->value.find( E_TGRA ) != s_N->value.end() )
            {
                fData->setMSCT( s_N->value[E_TGRA] );
                for( unsigned int j = 0; j < s_N->fNTel; j++ )
                {
                    fData->setMSCTT( j, s_N->value_T[E_TGRA][j], s_N->value_T_sigma[E_TGRA][j] );
                }
            }
            
            /////////////////////////
            // fill mean scaled frogs goodness (optional)
            if( s_N->value.find( E_FRGO ) != s_N->value.end() )
            {
                fData->setMSC_FRGO( s_N->value[E_FRGO] );
                for( unsigned int j = 0; j < s_N->fNTel; j++ )
                {
                    fData->setMSC_FRGOT( j, s_N->value_T[E_FRGO][j], s_N->value_T_sigma[E_FRGO][j] );
                }
            }
            
            fData->fill();
        }
        fevent++;
    }
}


/*

    interpolate between two lookup tables

*/
void VTableLookup::interpolate( VTablesToRead* s1, double w1, VTablesToRead* s2, double w2,
                                VTablesToRead* s, double w, bool iCos )
{

    map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
    
    unsigned int t = 0;
    for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); iter_iLT_Data++ )
    {
        // lookup table data pointer
        VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
        if( !iTableData )
        {
            continue;
        }
        
        // interpolate between two values
        s->value[t]  = VStatistics::interpolate( s1->value[t], w1,
                       s2->value[t], w2,
                       w, iCos, 0.1 );
        if( iTableData->fEnergy )
        {
            s->value_Chi2[t] = VStatistics::interpolate( s1->value_Chi2[t], w1,
                               s2->value_Chi2[t], w2,
                               w, iCos, 0. );
            s->value_dE[t] = VStatistics::interpolate( s1->value_dE[t], w1,
                             s2->value_dE[t], w2,
                             w, iCos, 0. );
        }
        
        for( unsigned int i = 0; i < s1->fNTel; i++ )
        {
            s->value_T[t][i] = VStatistics::interpolate( s1->value_T[t][i], w1,
                               s2->value_T[t][i], w2,
                               w, iCos, 1.e-2 );
            s->value_T_sigma[t][i] = VStatistics::interpolate( s1->value_T_sigma[t][i], w1,
                                     s2->value_T_sigma[t][i], w2,
                                     w, iCos, 1.e-2 );
        }
        // increment telescope counter
        t++;
    }
}


/*!
     write everything to disk
*/
void VTableLookup::terminate()
{
    fData->terminate( fTLRunParameter );
    
    if( fTLRunParameter->fWriteTables )
    {
        cout << "writing tables to disk (outputfile is " << fLookupTableFile->GetName() << ")" << endl;
        
        /// loop over all lookup table types
        map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
        for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); iter_iLT_Data++ )
        {
            // lookup table data pointer
            VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
            if( !iTableData )
            {
                continue;
            }
            iTableData->terminate( fLookupTableFile );
        }
        cout << "end of run (" << fLookupTableFile->GetName() << ")" << endl;
    }
    else
    {
        if( fNumberOfIgnoredEvents > 0 )
        {
            cout << endl << "\t total number of ignored events: " << fNumberOfIgnoredEvents << endl;
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    // large amount of objects read from subdirectory of the tablefile might result in
    // excessive time needed to close the tablefile
    // tablefile is therefore only close in table writing mode
    if( fTLRunParameter->fWriteTables )
    {
        cout << "closing file..." << endl;
        fLookupTableFile->Close();
    }
    
    cout << "exiting..." << endl;
}

/*

     get the corresponding bin number for an azimuth angle
     (no interpolation)

*/
int VTableLookup::getAzBin( double az )
{
    // be sure that az is bin interval [-180., 180.]
    if( az > 180. )
    {
        az -= 360.;
    }
    // expect first bin to be of type [135,-135]
    // otherwise like [-45,45]
    for( unsigned int i = 0; i < fTableAzLowEdge.size(); i++ )
    {
        if( i == 0 && ( az > fTableAzLowEdge[0] || az < fTableAzUpEdge[0] ) )
        {
            return 0;
        }
        else if( az > fTableAzLowEdge[i] && az < fTableAzUpEdge[i] )
        {
            return i;
        }
    }
    return 0;
}


void VTableLookup::setSpectralIndex( double iS )
{
    // expect positive spectral index
    if( iS < 0. )
    {
        iS *= -1.;
    }
    
    fData->setSpectralIndex( iS );
    
    if( TMath::Abs( iS - fData->getMCSpectralIndex() ) > 1.e-2 )
    {
        cout << "----------------------------------------------------------------" << endl;
        cout << "expect MC input spectrum with spectral index of " << fData->getMCSpectralIndex();
        cout << ", energy range [" << fData->getMCMinEnergy() << "," << fData->getMCMaxEnergy() << "] TeV" << endl;
        cout << "weight MC events to new spectral index of " <<  fData->getSpectralIndex() << endl;
    }
}



bool VTableLookup::setInputFiles( vector< string > iInputFiles )
{
    if( fTLRunParameter->fDebug )
    {
        cout << "VTableLookup::setInputFiles " << iInputFiles.size() << endl;
        for( unsigned int i = 0; i < iInputFiles.size(); i++ )
        {
            cout << "\t" << iInputFiles[i] << endl;
        }
    }
    bool bMC = fData->setInputFile( iInputFiles );
    fNTel = fData->getNTel();
    
    fTableCalculator = new VTableCalculator( fNTel );
    fTableCalculator->setMinRequiredShowerPerBin( fTLRunParameter->fMinRequiredShowerPerBin );
    
    return bMC;
}

/*
      read list of directories from table file

      sort them to make sure that they are always in the same sequence

*/
vector< string > VTableLookup::getSortedListOfDirectories( TDirectory* iDir )
{
    vector< string > iDName;
    
    if( !iDir )
    {
        return iDName;
    }
    
    bool bWoffAltered = false;
    
    TList* iKeyList = iDir->GetListOfKeys();
    if( iKeyList )
    {
        TIter next( iKeyList );
        while( TNamed* iK = ( TNamed* )next() )
        {
            iDName.push_back( iK->GetName() );
            if( iDName.back().substr( 0, 4 ) == "woff" && iDName.back().size() == 8 )
            {
                iDName.back() = "woff_0" + iDName.back().substr( 5, iDName.back().size() );
                bWoffAltered = true;
            }
        }
    }
    
    sort( iDName.begin(), iDName.end() );
    
    if( bWoffAltered )
    {
        for( unsigned int i = 0; i < iDName.size(); i++ )
        {
            if( iDName[i].substr( 0, 6 ) == "woff_0" )
            {
                iDName[i] = "woff_" + iDName[i].substr( 6, iDName[i].size() );
            }
        }
    }
    
    return iDName;
}


/*

     read pedvar values for this event

*/
void VTableLookup::readNoiseLevel( bool bWriteToRunPara )
{
    if( !fData )
    {
        return;
    }
    
    fNoiseLevel     = fData->getNoiseLevel( !bWriteToRunPara );
    fMeanNoiseLevel = fData->getMeanNoiseLevel( !bWriteToRunPara );
    
    fTelToAnalyze.assign( fNoiseLevel.size(), false );
    
    if( fTLRunParameter && bWriteToRunPara )
    {
        fTLRunParameter->meanpedvars = fMeanNoiseLevel;
        fTLRunParameter->pedvars     = fNoiseLevel;
    }
    if( bWriteToRunPara )
    {
        cout << "Mean pedvar per telescope: ";
        for( unsigned int i = 0; i < fNoiseLevel.size(); i++ )
        {
            if( TMath::Abs( fNoiseLevel[i] ) > 1.e-3 )
            {
                cout << " T" << i + 1 << ": " << fNoiseLevel[i];
            }
        }
        cout << endl;
    }
    // check noise levels
    for( unsigned int i = 0; i < fNoiseLevel.size(); i++ )
    {
        // PE simulations with no noise values (preliminary, should somehow be included into the MC)
        if( fTLRunParameter->fPE )
        {
            if( fData->getNtubes()[i] > 0 )
            {
                fTelToAnalyze[i] = true;
            }
            else
            {
                fTelToAnalyze[i] = false;
            }
        }
        // normal data or MC run
        else
        {
            if( TMath::Abs( fNoiseLevel[i] ) < 1.e-2 && fData->getNtubes()[i] > 0 )
            {
                fTelToAnalyze[i] = false;
                if( fNNoiseLevelWarnings < 30 )
                {
                    cout << "WARNING: noise level for telescope " << i + 1 << " very low: " << fNoiseLevel[i] << " (" << !bWriteToRunPara << ")" << endl;
                }
                else if( fNNoiseLevelWarnings == 30 )
                {
                    cout << "----------- more than 30 noise level warnings, stop printing...--------------" << endl;
                }
                fNNoiseLevelWarnings++;
            }
            else if( TMath::Abs( fNoiseLevel[i] ) < 1.e-2 && fData->getNtubes()[i] < 1 )
            {
                fTelToAnalyze[i] = false;
            }
            else
            {
                fTelToAnalyze[i] = true;
            }
        }
    }
    
    if( ( int )fNoiseLevel.size() != fNTel )
    {
        cout << " VTableLookup::readNoiseLevel ERROR: could not find mean pedestal variation for each telescope" << endl;
        cout << "exiting...." << endl;
        exit( EXIT_FAILURE );
    }
}

/*

    get the index of the two nearest elements in the vector iV to the value x

*/
void VTableLookup::getIndexBoundary( unsigned int* iup, unsigned int* ilow, vector< double >& iV, double x )
{
    if( iV.size() == 0 )
    {
        *iup = 0;
        *ilow = 0;
        return;
    }
    
    if( x <= iV[0] )
    {
        *iup = *ilow = 0;
    }
    else if( x >= iV[iV.size() - 1] )
    {
        *iup = *ilow = iV.size() - 1;
    }
    else
    {
        for( unsigned int i = 0; i < iV.size(); i++ )
        {
            if( x > iV[i] )
            {
                *ilow = ( unsigned int )i;
                *iup = ( unsigned int )( i + 1 );
            }
        }
    }
}

/*

    get pointers to a certain set of tables

*/
void VTableLookup::getTables( unsigned int inoise, unsigned int ize,
                              unsigned int iwoff, unsigned int iaz, unsigned int tel,
                              VTablesToRead* s )

{
    if( !s )
    {
        return;
    }
    
    unsigned int telX = getTelTypeCounter( tel );
    if( telX == 999999 )
    {
        cout << "VTableLookup::getTables invalid telescope type: " << tel << "\t" << telX << endl;
        cout << "(this means that there is no table in the table file given for the requested telescope type)" << endl;
        exit( EXIT_FAILURE );
    }
    if( fTLRunParameter->fDebug == 2 )
    {
        cout << "DEBUG  getTables() "  << inoise << " " << ize << " " << iwoff << " " << iaz << " " << telX << endl;
        cout << "DEBUG  " << fTableData[E_MSCW]->fTable[telX][inoise][ize][iwoff][iaz] << endl;
        cout << "DEBUG  MEDIAN (MSCW) " << inoise << " " << ize << " " << iwoff << " " << iaz << " " << telX << " ";
        if( fTableData[E_MSCW]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian() )
        {
            cout << fTableData[E_MSCW]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetEntries() << "\t";
            cout << fTableData[E_MSCW]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetTitle() << "\t";
            cout << fTableData[E_MSCW]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetDirectory()->GetPath() << endl;
        }
        cout << "DEBUG  MEDIAN (MSCW,2) " << fTableData[E_MSCW]->fTable[telX][inoise].size() << endl;
        
        cout << "DEBUG  MEDIAN (MSCL) " << inoise << " " << ize << " " << iwoff << " " << iaz << " " << telX << " ";
        if( fTableData[E_MSCL]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian() )
        {
            cout << fTableData[E_MSCL]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetEntries() << "\t";
            cout << fTableData[E_MSCL]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetTitle() << "\t";
            cout << fTableData[E_MSCL]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetDirectory()->GetPath() << endl;
        }
        cout << "DEBUG  MEDIAN (MSCL,2) " << fTableData[E_MSCL]->fTable[telX][inoise].size() << endl;
        
        cout << "DEBUG  MEDIAN (ENERGYSR) " << inoise << " " << ize << " " << iwoff << " " << iaz << " " << telX << " ";
        if( fTableData[E_EREC]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian() )
        {
            cout << fTableData[E_EREC]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetEntries() << "\t";
            cout << fTableData[E_EREC]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetTitle() << "\t";
            cout << fTableData[E_EREC]->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian()->GetDirectory()->GetPath() << endl;
        }
        else
        {
            cout << "NOTABLE!" << endl;
        }
        cout << "DEBUG  MEDIAN (ENERGYSR,2) " << fTableData[E_EREC]->fTable[telX][inoise].size() << endl;
    }
    
    
    map< unsigned int, VTableCalculatorData* >::iterator iter_iLT_Data;
    
    unsigned int t = 0;
    for( iter_iLT_Data = fTableData.begin(); iter_iLT_Data != fTableData.end(); ++iter_iLT_Data )
    {
        // lookup table data pointer
        VTableCalculatorData* iTableData = ( *iter_iLT_Data ).second;
        if( !iTableData )
        {
            continue;
        }
        
        s->hMedian[t][tel] = iTableData->fTable[telX][inoise][ize][iwoff][iaz]->getHistoMedian();
        t++;
    }
}

/*

        calculate mean scaled values and energies with help of lookup tables

*/
void VTableLookup::calculateMSFromTables( VTablesToRead* s )
{
    if( !s )
    {
        return;
    }
    // make sure that list of pointers to tables exists
    if( !fTableCalculator )
    {
        s->reset();
        return;
    }
    double i_dummy = 0.;
    
    // size2 vector
    double* i_s2 = fData->getSize2( 1., fTLRunParameter->fUseSelectedImagesOnly );
    // loss vector
    double* i_l  = fData->getLoss();
    // distance vector
    double* i_d = fData->getDistance();
    
    ///////////////////
    // calculate mscw
    fTableCalculator->setCalculateEnergies( false );
    fTableCalculator->setNormalizeTableValues( fTableData[E_MSCW]->fValueNormalizationRange_min,
            fTableData[E_MSCW]->fValueNormalizationRange_max );
    fTableCalculator->setEventSelectionCut();
    fTableCalculator->setVHistograms( s->hMedian[E_MSCW] );
    
    s->value[E_MSCW] = fTableCalculator->calc( ( int )fData->getNTel(), fData->getDistanceToCore(),
                       i_s2, i_l, i_d, fData->getWidth(),
                       s->value_T[E_MSCW],
                       i_dummy, i_dummy,
                       s->value_T_sigma[E_MSCW] );
                       
    ///////////////////
    // calculate mscl
    fTableCalculator->setCalculateEnergies( false );
    fTableCalculator->setNormalizeTableValues( fTableData[E_MSCL]->fValueNormalizationRange_min,
            fTableData[E_MSCL]->fValueNormalizationRange_max );
    fTableCalculator->setEventSelectionCut();
    fTableCalculator->setVHistograms( s->hMedian[E_MSCL] );
    
    s->value[E_MSCL] = fTableCalculator->calc( ( int )fData->getNTel(), fData->getDistanceToCore(),
                       i_s2, i_l, i_d, fData->getLength(),
                       s->value_T[E_MSCL],
                       i_dummy, i_dummy,
                       s->value_T_sigma[E_MSCL] );
                       
    ///////////////////
    // calculate energy (method 1)
    fTableCalculator->setCalculateEnergies( true );
    fTableCalculator->setNormalizeTableValues( fTableData[E_EREC]->fValueNormalizationRange_min,
            fTableData[E_EREC]->fValueNormalizationRange_max );
    fTableCalculator->setEventSelectionCut( fTLRunParameter->fEventSelectionCut_lossCutMax, 
                                            fTLRunParameter->fEventSelectionCut_distanceCutMax );
    fTableCalculator->setVHistograms( s->hMedian[E_EREC] );
    s->value[E_EREC] = fTableCalculator->calc( ( int )fData->getNTel(), fData->getDistanceToCore(),
                       i_s2, i_l, i_d, 0,
                       s->value_T[E_EREC],
                       s->value_Chi2[E_EREC],
                       s->value_dE[E_EREC],
                       s->value_T_sigma[E_EREC] );
                       
    ///////////////////
    // mean scaled time gradient
    if( fTableData.find( E_TGRA ) != fTableData.end() )
    {
        fTableCalculator->setCalculateEnergies( false );
        fTableCalculator->setNormalizeTableValues( fTableData[E_TGRA]->fValueNormalizationRange_min,
                fTableData[E_TGRA]->fValueNormalizationRange_max );
        fTableCalculator->setEventSelectionCut();
        fTableCalculator->setVHistograms( s->hMedian[E_TGRA] );
        
        s->value[E_TGRA] = fTableCalculator->calc( ( int )fData->getNTel(), fData->getDistanceToCore(),
                           i_s2, i_l, i_d, fData->getTimeGradient(),
                           s->value_T[E_TGRA],
                           i_dummy, i_dummy,
                           s->value_T_sigma[E_TGRA] );
    }
    
    ///////////////////
    // mean scaled frogs goodness of fit
    if( fTableData.find( E_FRGO ) != fTableData.end() )
    {
        fTableCalculator->setCalculateEnergies( false );
        fTableCalculator->setNormalizeTableValues( fTableData[E_FRGO]->fValueNormalizationRange_min,
                fTableData[E_FRGO]->fValueNormalizationRange_max );
        fTableCalculator->setEventSelectionCut();
        fTableCalculator->setVHistograms( s->hMedian[E_FRGO] );
        
        s->value[E_FRGO] = fTableCalculator->calc( ( int )fData->getNTel(), fData->getDistanceToCore(),
                           i_s2, i_l, i_d, fData->getFROGS_goodness(),
                           s->value_T[E_FRGO],
                           i_dummy, i_dummy,
                           s->value_T_sigma[E_FRGO] );
    }
}


bool VTableLookup::initialize()
{
    if( fTLRunParameter->fDebug > 0 )
    {
        cout << "VTableLookup::initialize()" << endl;
    }
    if( !fTLRunParameter )
    {
        cout << "VTableLookup::initialize: error: no table lookup run parameters " << endl;
        return false;
    }
    // init data handler
    fData = new VTableLookupDataHandler( fTLRunParameter->fWriteTables, fTLRunParameter );
    fData->setfillTables( fTLRunParameter->fWriteTables );
    
    cout << endl;
    // set input file names
    fTLRunParameter->isMC = setInputFiles( fTLRunParameter->inputfile );
    
    ///////////////////////////////////////
    // write mscw,mscl, and energy tables
    if( fTLRunParameter->fWriteTables )
    {
        char ihname[900];
        if( fTLRunParameter->inputfile.size() > 0 )
        {
            sprintf( ihname, "lookup table file (array recid = %d, source files: %s)", fTLRunParameter->rec_method, fTLRunParameter->inputfile[0].c_str() );
        }
        else
        {
            sprintf( ihname, "lookup table file (array recid = %d)", fTLRunParameter->rec_method );
        }
        
        string iTitle = ihname;
        cout << "setting mean pedvar level for table selection to : " << fData->getMeanNoiseLevel() << endl;
        cout << "   (mean pedvar levels per telescope type are ";
        map< ULong64_t, double > i_pedvarlevel =  fData->getNoiseLevel_per_TelescopeType();
        map<ULong64_t, double >::iterator iList_of_Tel_type_iterator;
        for( iList_of_Tel_type_iterator = i_pedvarlevel.begin(); iList_of_Tel_type_iterator != i_pedvarlevel.end();
                ++iList_of_Tel_type_iterator )
        {
            cout << " type " << iList_of_Tel_type_iterator->first << ": " << iList_of_Tel_type_iterator->second;
        }
        cout << ")" << endl << endl;
        
        
        // set up tables
        setMCTableFiles_forTableWriting( fTLRunParameter->tablefile, fTLRunParameter->ze, fTLRunParameter->fWobbleOffset,
                                         i_pedvarlevel, "tb", ihname );
        fData->readRunParameter();
        // set min/max distance to camera center
        if( fData )
        {
            fData->setMCDistanceToCameraCenter( fTLRunParameter->fMC_distance_to_cameracenter_min, fTLRunParameter->fMC_distance_to_cameracenter_max );
        }
    }
    ///////////////////////////////////////
    // read MC tables
    else
    {
        // read pedvars
        readNoiseLevel( true );
        // read tables from disk
        setMCTableFiles_forTableReading( fTLRunParameter->tablefile, "tb" );
        // set output files
        setOutputFile( fTLRunParameter->outputfile, fTLRunParameter->writeoption, fTLRunParameter->tablefile );
    }
    
    // weight event to this spectral index
    if( fTLRunParameter->isMC )
    {
        setSpectralIndex( fTLRunParameter->fSpectralIndex );
    }
    
    return true;
}

unsigned int VTableLookup::getWobbleBin( double w )
{
    // no wobble bins defined - return 0 (first vector element)
    if( fTLRunParameter->fCTA_MC_offaxisBin_min.size() == 0 )
    {
        return 0;
    }
    
    for( unsigned int i = 0; i < fTLRunParameter->fCTA_MC_offaxisBin_min.size(); i++ )
    {
        if( w >= fTLRunParameter->fCTA_MC_offaxisBin_min[i] && w < fTLRunParameter->fCTA_MC_offaxisBin_max[i] )
        {
            return i;
        }
    }
    
    return 9999;
}

unsigned int VTableLookup::getTelTypeCounter( unsigned int tel, bool iStopIfError )
{
    unsigned int telX = 999999;
    for( unsigned int i = 0; i < fTableTelTypes.size(); i++ )
    {
        if( fData->getTelType( tel ) == fTableTelTypes[i] )
        {
            telX = i;
            break;
        }
    }
    if( iStopIfError && telX == 999999 )
    {
        cout << "VTableLookup::getTelTypeCounter error finding telescope type for telescope " << tel << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    return telX;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

VTableCalculatorData::VTableCalculatorData()
{
    fDirectoryName = "";
    fFillVariable = "";
    fEnergy = false;
    
    fValueNormalizationRange_min = -9999.;
    fValueNormalizationRange_max = -9999.;
}

/*
 * call terminate function for lookup table code and write tables to disk
 *
 */
void VTableCalculatorData::terminate( TFile* iFile )
{
    char hname[800];
    for( unsigned int i = 0; i < fTable.size(); i++ )
    {
        for( unsigned int t = 0; t < fTable[i].size(); t++ )
        {
            for( unsigned int u = 0; u < fTable[i][t].size(); u++ )
            {
                for( unsigned int v = 0; v < fTable[i][t][u].size(); v++ )
                {
                    for( unsigned w = 0; w < fTable[i][t][u][v].size(); w++ )
                    {
                        sprintf( hname, "TEL %u NOISE %u ZE %u WOFF %u AZ %u", i + 1, t, u, v, w );
                        cout << "writing " << fDirectoryName << " tables for ";
                        cout << fTable[i][t][u][v][w]->getOutputDirectory()->GetMotherDir()->GetPath();
                        cout << " : " << hname << endl;
                        fTable[i][t][u][v][w]->terminate( fTable[i][t][u][v][w]->getOutputDirectory(), hname );
                        if( iFile )
                        {
                            iFile->Flush();
                        }
                    }
                }
            }
        }
    }
}

/*
 * assert that table vector is not empty
 *
 */
bool VTableCalculatorData::assertTableVector( unsigned int wobble_bin )
{

    if( fTable.size() == 0
            || fTable[0].size() == 0
            || fTable[0][0].size() == 0
            || wobble_bin >= fTable[0][0][0].size() )
    {
        return false;
    }
    
    return true;
}

void VTableCalculatorData::print()
{
    cout << "Lookup table ";
    cout << fFillVariable;
    cout << " (directory: " << fDirectoryName;
    if( fEnergy )
    {
        cout << ", energy tables";
    }
    cout << ") ";
    if( fValueNormalizationRange_min > -9990. && fValueNormalizationRange_max > -9990. )
    {
        cout << "normalization from [" << fValueNormalizationRange_min;
        cout << ", " << fValueNormalizationRange_max << "] to [0,1]";
    }
    cout << endl;
}
