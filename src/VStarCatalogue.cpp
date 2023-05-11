/*! \class VStarCatalogue
 *  \brief bright star catalogue

*/

#include "VStarCatalogue.h"

ClassImp( VStarCatalogue )

VStarCatalogue::VStarCatalogue()
{
    fDebug = false;
    
    fCatalogue = "Hipparcos_MAG8_1997.dat";
    fCatalogueVersion = 0;
    
    setTelescopePointing();
}

VStarCatalogue::~VStarCatalogue()
{
    for( unsigned int i = 0; i < fStars.size(); i++ )
    {
        if( fStars[i] )
        {
            delete fStars[i];
        }
    }
    for( unsigned int i = 0; i < fStarsinFOV.size(); i++ )
    {
        if( fStarsinFOV[i] )
        {
            delete fStarsinFOV[i];
        }
    }
}


bool VStarCatalogue::init( double MJD )
{
    return init( MJD, fCatalogue );
}


bool VStarCatalogue::init( double iMJD, string iCatalogue )
{
    fCatalogue = iCatalogue;
    
    if( !readCatalogue() )
    {
        return false;
    }
    double dec, ra;
    double i_b, i_l;
    for( unsigned int i = 0; i < fStars.size(); i++ )
    {
        dec = fStars[i]->fDec2000 * TMath::Pi() / 180.;
        ra =  fStars[i]->fRA2000 * TMath::Pi() / 180.;
        // calculate galac coordinates
        VAstronometry::vlaEqgal( ra, dec, &i_l, &i_b );
        fStars[i]->fRunGalLong1958 = i_l * 180. / TMath::Pi();
        
        // apply precesssion
        VAstronometry::vlaPreces( 2451545.0 - 2400000.5, iMJD, &ra, &dec );
        // calculate ra/dec for current epoch
        fStars[i]->fDecCurrentEpoch = dec * 180. / TMath::Pi();
        fStars[i]->fRACurrentEpoch = ra * 180. / TMath::Pi();
        fStars[i]->fRunGalLat1958  = i_b * 180. / TMath::Pi();
    }
    return true;
}


/*

 */
bool VStarCatalogue::readVERITASsourcesfromDB( string iofile )
{
    char c_query[1000];
    
    stringstream iTempS;
    iTempS << getDBServer() << "/VERITAS";
    
    //std::cout<<"VStarCatalogue::readVERITASsourcesfromDB "<<std::endl;
    VDB_Connection my_connection( iTempS.str().c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "VStarCatalogue: failed to connect to database server" << endl;
        return false;
    }
    
    sprintf( c_query, "select * from tblObserving_Sources " );
    
    
    if( !my_connection.make_query( c_query ) )
    {
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    
    int fNRows = db_res->GetRowCount();
    
    unsigned int zID = 0;
    double ra = 0.;
    double dec = 0.;
    
    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();
        
        if( db_row->GetField( 0 ) && db_row->GetField( 1 ) && db_row->GetField( 2 ) )
        {
            VStar* i_Star = new VStar();
            i_Star->fStarID = zID;
            i_Star->fStarName = db_row->GetField( 0 );
            
            // don't read the dark spots
            if( i_Star->fStarName.substr( 0, 5 ) == "DARK_" )
            {
                continue;
            }
            // don't read bright stars
            if( i_Star->fStarName.substr( 0, 3 ) == "BSC" )
            {
                continue;
            }
            // don't read different Crab pointings
            if( i_Star->fStarName.substr( 0, 4 ) == "Crab" && i_Star->fStarName.size() > 4 )
            {
                continue;
            }
            // don't read zenith runs
            if( i_Star->fStarName == "ZENITH" )
            {
                continue;
            }
            
            ra = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
            dec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
            
            i_Star->fRA2000 = ra;
            i_Star->fDec2000 = dec;
            i_Star->fBrightness_V = -9999.;
            i_Star->fBrightness_B = -9999.;
            i_Star->fMajorDiameter = 0.;
            i_Star->fMajorDiameter_68 = 0.;
            
            cout << "ID " << i_Star->fRA2000 << " " << i_Star->fDec2000 << endl;
            
            fStars.push_back( i_Star );
            zID++;
        }
    }
    
    if( my_connection.Get_Connection_Status() )
    {
        my_connection.Close_Connection();    // just so it get close as soon as possible. Before the end of the function.
    }
    
    
    // write sources into a text file
    if( iofile.size() > 0 )
    {
        ofstream os;
        os.open( iofile.c_str() );
        if( !os )
        {
            cout << "error opening file for VERITAS sources: " << iofile << endl;
            return false;
        }
        
        os << "* V 4" << endl;
        for( unsigned int i = 0; i < fStars.size(); i++ )
        {
            if( fStars[i]->fStarName.substr( 0, 2 ) == "SS" )
            {
                continue;
            }
            if( fStars[i]->fStarName.substr( 0, 3 ) == "GRB" )
            {
                continue;
            }
            os << fStars[i]->fRA2000 << "\t" << fStars[i]->fDec2000 << "\t" << fStars[i]->fStarName  << endl;
        }
        os.close();
    }
    
    return true;
}


/*!
     attention: very dependend on format of text file

     iV < 3:  Brightstarcatalogue
     iV == 3: tevcat
     iV == 4: HMX catalogue
*/
bool VStarCatalogue::readCatalogue()
{
    //////////////////////////////////////
    // READ VERITAS object catalogue from DB
    if( fCatalogue == "VERITASDB" )
    {
        return readVERITASsourcesfromDB( "" );
    }
    //////////////////////////////////////
    
    //////////////////////////////////////
    // READ catalogue from an ascii file
    //////////////////////////////////////
    
    ifstream is;
    is.open( fCatalogue.c_str(), ifstream::in );
    if( !is )
    {
        // try ENVDISPDATA
        string itemp = "";
        if( gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" ) )
        {
            itemp = gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" );
        }
        if( itemp.size() > 0 )
        {
            fCatalogue   = itemp + "/AstroData/Catalogues/" + fCatalogue;
            is.open( fCatalogue.c_str(), ifstream::in );
        }
        if( !is )
        {
            cout << "VStarCatalogue::readCatalogue error: file not found, " << fCatalogue << endl;
            return false;
        }
    }
    string iLine;
    string iLine_sub;
    string iT1;
    string iT2;
    string iT3;
    cout << "\treading star catalogue: " << fCatalogue << endl;
    int zid = 0;
    // catalogue version
    fCatalogueVersion = 0;
    
    while( getline( is, iLine ) )
    {
        // hard maximum number of sources of 150,000 to avoid memory leaks
        if( zid > 150000 )
        {
            cout << "VStarCatalogue::init: too many objects in catalogue, possible memory leak..." << endl;
            return false;
        }
        if( iLine.size() < 2 )
        {
            continue;
        }
        // skip the commend lines
        if( iLine.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        
        // read catalogue version
        // (note: this defines the expected layout of the ascii file)
        if( iLine.substr( 0, 1 ) == "*" )
        {
            istringstream is_stream( iLine );
            is_stream >> iT1;
            is_stream >> iT1;
            is_stream >> iT1;
            fCatalogueVersion = atoi( iT1.c_str() );
            continue;
        }
        ////////////////////////////////////////////
        // this is the text line to be worked with
        iLine_sub = iLine;
        ////////////////////////////////////////////
        // check number of text blocks available
        if( !checkTextBlocks( iLine_sub, fCatalogueVersion ) )
        {
            if( fDebug )
            {
                cout << "Skipping following line in catalogue (version " << fCatalogueVersion << ")" << endl;
                cout << iLine_sub << endl;
            }
            continue;
        }
        // start a new star
        VStar* i_Star = new VStar();
        i_Star->fQualityFlag = 0;
        if( fCatalogueVersion == 3 )
        {
            i_Star->fStarName = iLine.substr( 0, 20 );
            iLine_sub = iLine.substr( 23, iLine.size() );
        }
        if( fCatalogueVersion == 11 )
        {
            i_Star = readCommaSeparatedLine_Fermi( iLine_sub, zid, i_Star );
        }
        else if( fCatalogueVersion == 13 )
        {
            i_Star = readCommaSeparatedLine_Fermi_Catalogue( iLine_sub, zid, i_Star );
        }
        else if( fCatalogueVersion == 14 )
        {
            i_Star = readCommaSeparatedLine_FAVA( iLine_sub, zid, i_Star );
        }
        else if( fCatalogueVersion == 15 )
        {
            i_Star = readCommaSeparatedLine_Fermi2nd_Catalogue( iLine_sub, zid, i_Star );
        }
        else
        {
            istringstream is_stream( iLine_sub );
            i_Star->fVariability = false;
            // star ID
            if( fCatalogueVersion < 2 )
            {
                is_stream >> iT1;
                i_Star->fStarID = ( unsigned int )atoi( iT1.c_str() );
                i_Star->fStarName = iT1;
            }
            else if( fCatalogueVersion == 3 || fCatalogueVersion == 4 || fCatalogueVersion == 6 )
            {
                i_Star->fStarID = ( unsigned int )zid;
            }
            else
            {
                i_Star->fStarID = ( unsigned int )zid;
            }
            
            if( fCatalogueVersion == 5 || fCatalogueVersion == 9 )
            {
                is_stream >> iT1;
                // RA2000
                i_Star->fRA2000 = atof( iT1.c_str() );
                is_stream >> iT1;
                // Dec2000
                i_Star->fDec2000 = atof( iT1.c_str() );
            }
            else if( fCatalogueVersion < 10 || fCatalogueVersion == 12 )
            {
                // RA2000
                is_stream >> iT1;
                // make sure that entry exists
                if( ( is_stream >> std::ws ).eof() )
                {
                    zid++;
                    continue;
                }
                is_stream >> iT2;
                if( ( is_stream >> std::ws ).eof() )
                {
                    zid++;
                    continue;
                }
                if( fCatalogueVersion != 6 && fCatalogueVersion != 7 )
                {
                    is_stream >> iT3;
                }
                else
                {
                    iT3 = "0";
                }
                i_Star->fRA2000 = 15.*( atof( iT1.c_str() ) + atof( iT2.c_str() ) / 60. + atof( iT3.c_str() ) / 3600. );
                // dec2000
                is_stream >> iT1;
                is_stream >> iT2;
                if( fCatalogueVersion != 6 && fCatalogueVersion != 7 )
                {
                    is_stream >> iT3;
                }
                else
                {
                    iT3 = "0";
                }
                if( iT1.find( "-", 0 ) != string::npos )
                {
                    i_Star->fDec2000 = atof( iT1.c_str() ) - atof( iT2.c_str() ) / 60. - atof( iT3.c_str() ) / 3600.;
                }
                else
                {
                    i_Star->fDec2000 = atof( iT1.c_str() ) + atof( iT2.c_str() ) / 60. + atof( iT3.c_str() ) / 3600.;
                }
            }
            i_Star->fBrightness_V = 9999;
            i_Star->fBrightness_B = 9999;
            i_Star->fMajorDiameter = 0.;
            i_Star->fMajorDiameter_68 = 0.;
            // objet name
            if( fCatalogueVersion == 4 )
            {
                i_Star->fStarName = iLine.substr( 24, iLine.size() );
            }
            else if( fCatalogueVersion == 8 )
            {
                i_Star->fStarName = iLine.substr( 24, 16 );
                i_Star->fStarName = "1RXS" + i_Star->fStarName;
                i_Star->fMajorDiameter = atof( iLine.substr( 40, 4 ).c_str() );
                i_Star->fMajorDiameter /= 3600.;   // arcsec -> deg
            }
            else if( fCatalogueVersion == 5 )
            {
                i_Star->fStarName = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() ).c_str();
                i_Star->fStarName = VUtilities::remove_leading_spaces( i_Star->fStarName );
            }
            else if( fCatalogueVersion == 6 || fCatalogueVersion == 7 )
            {
                i_Star->fStarName = iLine.substr( 15, 11 );
                if( fCatalogueVersion == 7 )
                {
                    i_Star->fStarName = "3EG " + i_Star->fStarName;
                }
                i_Star->fMajorDiameter = atof( iLine.substr( 27, 5 ).c_str() );
                // arcmin -> deg
                if( fCatalogueVersion == 6 )
                {
                    i_Star->fMajorDiameter /= 60.;
                }
            }
            else if( fCatalogueVersion == 9 )
            {
                i_Star->fStarName = iLine.substr( is_stream.tellg(), iLine.size() );
            }
            else if( fCatalogueVersion == 12 )
            {
                is_stream >> i_Star->fStarName;
            }
            else if( fCatalogueVersion < 3 )
            {
                // brightness
                is_stream >> i_Star->fBrightness_V;
                is_stream >> i_Star->fBrightness_B;
                // here check for B = 9999
                if( TMath::Abs( i_Star->fBrightness_B - 9999. ) > 1. )
                {
                    i_Star->fBrightness_B += i_Star->fBrightness_V;
                }
                else
                {
                    i_Star->fBrightness_B = 9999;
                }
            }
            
            // Hipparcos catalogue
            if( fCatalogueVersion == 10 )
            {
                is_stream >> iT1;
                i_Star->fRA2000 = atof( iT1.c_str() );
                is_stream >> iT1;
                i_Star->fDec2000 = atof( iT1.c_str() );
                is_stream >> iT1;
                i_Star->fStarID = atoi( iT1.c_str() );
                i_Star->fStarName = "HIP" + iT1;
                is_stream >> iT1;
                i_Star->fBrightness_V = atof( iT1.c_str() );
                is_stream >> iT1;
                i_Star->fBrightness_B = i_Star->fBrightness_V + atof( iT1.c_str() );
            }
        }
        
        fStars.push_back( i_Star );
        
        zid++;
    }
    is.close();
    
    return true;
}

VStar* VStarCatalogue::readCommaSeparatedLine_FAVA( string iLine, int zid, VStar* i_Star )
{
    string iT1;
    string iT2;
    string iT3;
    
    i_Star->fStarID = zid;
    i_Star->fBrightness_V = 9999.;
    i_Star->fBrightness_B = 9999.;
    
    if( iLine.size() == 0 )
    {
        return i_Star;
    }
    
    string iTemp = iLine;
    
    // star name
    i_Star->fStarName = iTemp.substr( 0, iTemp.find( "," ) );
    
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ra, dec
    istringstream is_stream( iTemp );
    is_stream >> iT1;
    is_stream >> iT2;
    is_stream >> iT3;
    cout << "RA " << zid << " " << iT1.c_str() << " " << iT2.c_str() << " " << iT3.c_str() << endl;
    i_Star->fRA2000 = 15.*( atof( iT1.c_str() ) + atof( iT2.c_str() ) / 60. + atof( iT3.c_str() ) / 3600. );
    // dec2000
    is_stream >> iT1;
    is_stream >> iT1;
    is_stream >> iT2;
    is_stream >> iT3;
    cout << "Dec " << zid << " " << iT1.c_str() << " " << iT2.c_str() << " " << iT3.c_str() << endl;
    if( iT1.find( "-", 0 ) != string::npos )
    {
        i_Star->fDec2000 = atof( iT1.c_str() ) - atof( iT2.c_str() ) / 60. - atof( iT3.c_str() ) / 3600.;
    }
    else
    {
        i_Star->fDec2000 = atof( iT1.c_str() ) + atof( iT2.c_str() ) / 60. + atof( iT3.c_str() ) / 3600.;
    }
    
    return i_Star;
}

/*
    2nd Fermi catalogue

*/
VStar* VStarCatalogue::readCommaSeparatedLine_Fermi2nd_Catalogue( string iLine, int zid, VStar* i_Star )
{
    i_Star->fStarID = zid;
    i_Star->fBrightness_V = 9999.;
    i_Star->fBrightness_B = 9999.;
    
    if( iLine.size() == 0 )
    {
        return i_Star;
    }
    
    string iTemp = iLine;
    
    // star name
    i_Star->fStarName = iTemp.substr( 0, iTemp.find( "," ) );
    
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ra, dec
    i_Star->fRA2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fDec2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    // galactic latitude/longitude are calculated from ra, dec
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ignore 68% values on position
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMajorDiameter_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMinorDiameter_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fPositionAngle_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 95% confidence radius
    i_Star->fMajorDiameter = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMinorDiameter = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fPositionAngle = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // significance
    i_Star->fSignificance = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    for( unsigned int i = 0; i < 4; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // spectral index
    i_Star->fSpectralIndex = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fSpectralIndexError = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 1 GeV to 100 GeV flux
    i_Star->fFluxEnergyMin.push_back( 1. );
    i_Star->fFluxEnergyMax.push_back( 1.e2 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    // spectral type
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fSpectrumType = iTemp.substr( 0, iTemp.find( "," ) );
    
    for( unsigned int i = 0; i < 18; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // cutoff energy
    i_Star->fCutOff_MeV = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fCutOffError_MeV = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    for( unsigned int i = 0; i < 6; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // 100 MeV to 300 MeV
    i_Star->fFluxEnergyMin.push_back( 0.1 );
    i_Star->fFluxEnergyMax.push_back( 0.3 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 300 MeV to 1 GeV
    i_Star->fFluxEnergyMin.push_back( 0.3 );
    i_Star->fFluxEnergyMax.push_back( 1.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 1 GeV to 3 GeV
    i_Star->fFluxEnergyMin.push_back( 1.0 );
    i_Star->fFluxEnergyMax.push_back( 3.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 3 GeV to 10 GeV
    i_Star->fFluxEnergyMin.push_back( 3.0 );
    i_Star->fFluxEnergyMax.push_back( 10.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 10 GeV to 100 GeV
    i_Star->fFluxEnergyMin.push_back( 10.0 );
    i_Star->fFluxEnergyMax.push_back( 100.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    for( unsigned int i = 0; i < 113; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // Other names
    for( unsigned int i = 0; i < 72; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 )
        {
            if( iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
            {
                i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
            }
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    for( unsigned int i = 0; i < 26; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // classification
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 && iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
    {
        i_Star->fType =  iTemp.substr( 0, iTemp.find( "," ) );
    }
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    for( unsigned int i = 0; i < 5; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
        {
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    for( unsigned int i = 0; i < 22; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 )
        {
            if( iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
            {
                i_Star->fAssociations.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
            }
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    i_Star->fQualityFlag = atoi( iTemp.substr( iTemp.rfind( "," ) + 1, iTemp.size() ).c_str() );
    return i_Star;
}

VStar* VStarCatalogue::readCommaSeparatedLine_Fermi_Catalogue( string iLine, int zid, VStar* i_Star )
{
    i_Star->fStarID = zid;
    i_Star->fBrightness_V = 9999.;
    i_Star->fBrightness_B = 9999.;
    
    if( iLine.size() == 0 )
    {
        return i_Star;
    }
    
    string iTemp = iLine;
    
    // star name
    i_Star->fStarName = iTemp.substr( 0, iTemp.find( "," ) );
    
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ra, dec
    i_Star->fRA2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fDec2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    // galactic latitude/longitude are calculated from ra, dec
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ignore 68% values on position
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMajorDiameter_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMinorDiameter_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fPositionAngle_68 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 95% confidence radius
    i_Star->fMajorDiameter = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fMinorDiameter = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fPositionAngle = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // significance
    i_Star->fSignificance = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    for( unsigned int i = 0; i < 4; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // spectral index
    i_Star->fSpectralIndex = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fSpectralIndexError = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 1 GeV to 100 GeV flux
    i_Star->fFluxEnergyMin.push_back( 1. );
    i_Star->fFluxEnergyMax.push_back( 1.e2 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    for( unsigned int i = 0; i < 6; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // 100 MeV to 300 MeV
    i_Star->fFluxEnergyMin.push_back( 0.1 );
    i_Star->fFluxEnergyMax.push_back( 0.3 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 300 MeV to 1 GeV
    i_Star->fFluxEnergyMin.push_back( 0.3 );
    i_Star->fFluxEnergyMax.push_back( 1.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 1 GeV to 3 GeV
    i_Star->fFluxEnergyMin.push_back( 1.0 );
    i_Star->fFluxEnergyMax.push_back( 3.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 3 GeV to 10 GeV
    i_Star->fFluxEnergyMin.push_back( 3.0 );
    i_Star->fFluxEnergyMax.push_back( 10.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 10 GeV to 100 GeV
    i_Star->fFluxEnergyMin.push_back( 10.0 );
    i_Star->fFluxEnergyMax.push_back( 100.0 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    for( unsigned int i = 0; i < 28; i++ )
    {
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    
    // Other names
    for( unsigned int i = 0; i < 72; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 )
        {
            if( iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
            {
                i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
            }
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    // classification
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 && iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
    {
        i_Star->fType =  iTemp.substr( 0, iTemp.find( "," ) );
    }
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    for( unsigned int i = 0; i < 5; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
        {
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    for( unsigned int i = 0; i < 22; i++ )
    {
        if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 1 )
        {
            if( iTemp.substr( 0, iTemp.find( "," ) ) != "  " )
            {
                i_Star->fAssociations.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
            }
            iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
        }
    }
    i_Star->fQualityFlag = atoi( iTemp.substr( iTemp.rfind( "," ) + 1, iTemp.size() ).c_str() );
    return i_Star;
}


VStar* VStarCatalogue::readCommaSeparatedLine_Fermi( string iLine, int zid, VStar* i_Star )
{
    i_Star->fStarID = zid;
    i_Star->fBrightness_V = 9999.;
    i_Star->fBrightness_B = 9999.;
    
    if( iLine.size() == 0 )
    {
        return i_Star;
    }
    
    string iTemp = iLine;
    
    // star name
    i_Star->fStarName = iTemp.substr( 0, iTemp.find( "," ) );
    
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // ra, dec
    i_Star->fRA2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fDec2000 = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    // galactic latitude/longitude are calculated from ra, dec
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 95% confiduence radius
    i_Star->fMajorDiameter = atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() );
    
    // ignore likelihood test
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 100 MeV to 1 GeV flux
    i_Star->fFluxEnergyMin.push_back( 1.e-2 );
    i_Star->fFluxEnergyMax.push_back( 1. );
    
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // 1 GeV to 100 GeV flux
    i_Star->fFluxEnergyMin.push_back( 1. );
    i_Star->fFluxEnergyMax.push_back( 1.e2 );
    i_Star->fFlux.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    i_Star->fFluxError.push_back( atof( iTemp.substr( 0, iTemp.find( "," ) ).c_str() ) );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    // Variability
    if( iTemp.substr( 0, iTemp.find( "," ) ) == "F" )
    {
        i_Star->fVariability = false;
    }
    else
    {
        i_Star->fVariability = true;
    }
    
    // Other names
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fOtherNames.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fAssociations.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    if( iTemp.substr( 0, iTemp.find( "," ) ).size() > 0 )
    {
        i_Star->fAssociations.push_back( iTemp.substr( 0, iTemp.find( "," ) ) );
        iTemp = iTemp.substr( iTemp.find( "," ) + 1, iTemp.size() );
    }
    return i_Star;
}


void VStarCatalogue::printCatalogue( unsigned int i_nRows, double iMinBrightness, string iBand )
{
    if( i_nRows == 0 || i_nRows >= fStars.size() )
    {
        i_nRows = fStars.size();
    }
    
    char hname[600];
    
    for( unsigned int i = 0; i < i_nRows; i++ )
    {
        if( iBand == "V" && fStars[i]->fBrightness_V > iMinBrightness )
        {
            continue;
        }
        else if( iBand == "B" && fStars[i]->fBrightness_B > iMinBrightness )
        {
            continue;
        }
        
        cout << fStars[i]->fStarID << "\t_" << fStars[i]->fStarName;
        cout << ", ra2000 = " << fStars[i]->fRA2000 << ", dec2000 = " << fStars[i]->fDec2000;
        cout << ", l = " << fStars[i]->fRunGalLong1958 << ", b = " << fStars[i]->fRunGalLat1958 << "\t";
        if( fStars[i]->fBrightness_V < 999. )
        {
            cout << "\t mag_V = " << fStars[i]->fBrightness_V;
        }
        if( fStars[i]->fBrightness_B < 999. )
        {
            cout << ", mag_B = " << fStars[i]->fBrightness_B;
        }
        cout << endl;
        cout << "\tmajor (95\%) " << fStars[i]->fMajorDiameter;
        if( fCatalogueVersion == 13 )
        {
            cout << ", minor (95\%) " << fStars[i]->fMinorDiameter;
            cout << ", pos angle (95\%) " << fStars[i]->fPositionAngle << endl;
            cout << "\tmajor (68\%) " << fStars[i]->fMajorDiameter_68;
            cout << ", minor (68\%) " << fStars[i]->fMinorDiameter_68;
            cout << ", pos angle (68\%) " << fStars[i]->fPositionAngle_68 << endl;
            cout << "\tsignificance " << fStars[i]->fSignificance;
            cout << ", index " << fStars[i]->fSpectralIndex;
            cout << " +- " << fStars[i]->fSpectralIndexError;
        }
        cout << endl;
        if( fStars[i]->fFluxEnergyMin.size() > 0 && fStars[i]->fFluxEnergyMin.size() == fStars[i]->fFluxEnergyMax.size()
                && fStars[i]->fFlux.size() == fStars[i]->fFluxEnergyMin.size() && fStars[i]->fFluxError.size() == fStars[i]->fFluxEnergyMin.size() )
        {
            for( unsigned int e = 0; e < fStars[i]->fFluxEnergyMin.size(); e++ )
            {
                cout << "\t";
                sprintf( hname, "Flux (%.2f GeV - %.2f GeV): ", fStars[i]->fFluxEnergyMin[e], fStars[i]->fFluxEnergyMax[e] );
                cout << hname;
                sprintf( hname, "\t(%.4e +- %.4e) photons/cm2/s", fStars[i]->fFlux[e], fStars[i]->fFluxError[e] );
                cout << hname;
                cout << endl;
            }
        }
        if( fStars[i]->fOtherNames.size() > 0 )
        {
            cout << "\tOther names: ";
            for( unsigned int e = 0; e < fStars[i]->fOtherNames.size(); e++ )
            {
                cout << fStars[i]->fOtherNames[e] << ", ";
            }
            cout << endl;
        }
        if( fStars[i]->fType.size() > 2 )
        {
            cout << "\tType: " << fStars[i]->fType << endl;
        }
        if( fStars[i]->fAssociations.size() > 0 )
        {
            cout << "\tAssociations: ";
            for( unsigned int e = 0; e < fStars[i]->fAssociations.size(); e++ )
            {
                cout << fStars[i]->fAssociations[e] << ", ";
            }
            cout << endl;
        }
        cout << "Quality flag: " << fStars[i]->fQualityFlag << endl;
    }
}


void VStarCatalogue::printStarsInFOV()
{
    printStarsInFOV( 999999. );
}


void VStarCatalogue::printStarsInFOV( double iMinBrightness, string iBand )
{
    for( unsigned int i = 0; i < fStarsinFOV.size(); i++ )
    {
        if( iBand == "V" && fStarsinFOV[i]->fBrightness_V > iMinBrightness )
        {
            continue;
        }
        if( iBand == "B" && fStarsinFOV[i]->fBrightness_B > iMinBrightness )
        {
            continue;
        }
        
        cout << fStarsinFOV[i]->fStarID << "\t" << fStars[i]->fStarName;
        cout << "  RA2000: " << fStarsinFOV[i]->fRA2000 << "  DEC2000: " << fStarsinFOV[i]->fDec2000 << "\t";
        cout << fStarsinFOV[i]->fRACurrentEpoch << "\t" << fStarsinFOV[i]->fDecCurrentEpoch;
        cout << " V: " << fStarsinFOV[i]->fBrightness_V;
        cout << " B: " << fStarsinFOV[i]->fBrightness_B;
        cout << " Diameter: " << fStarsinFOV[i]->fMajorDiameter << endl;
    }
}


double VStarCatalogue::getStarMajorDiameter( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fMajorDiameter;
    }
    
    return 0.;
}


double VStarCatalogue::getStarBrightness( unsigned int iID, string iBand )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        if( iBand == "B" )
        {
            return fStars[iID]->fBrightness_B;
        }
        if( iBand == "V" )
        {
            return fStars[iID]->fBrightness_V;
        }
    }
    
    return 0.;
}


double VStarCatalogue::getStarDecCurrentEpoch( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fDecCurrentEpoch;
    }
    
    return 0.;
}


double VStarCatalogue::getStarRACurrentEpoch( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fRACurrentEpoch;
    }
    
    return 0.;
}


string VStarCatalogue::getStarName( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fStarName;
    }
    
    string iN = "NONAME";
    
    return iN;
}


unsigned int VStarCatalogue::getStarID( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fStarID;
    }
    
    return 0;
}


double VStarCatalogue::getStarDec2000( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fDec2000;
    }
    
    return 0.;
}


double VStarCatalogue::getStarRA2000( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fRA2000;
    }
    
    return 0.;
}


double VStarCatalogue::getStar_l( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fRunGalLong1958;
    }
    
    return 0.;
}


double VStarCatalogue::getStar_b( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fRunGalLat1958;
    }
    
    return 0.;
}


vector< double > VStarCatalogue::getStarFluxEnergyMin( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fFluxEnergyMin;
    }
    
    vector< double > a;
    
    return a;
}


vector< double > VStarCatalogue::getStarFluxEnergyMax( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fFluxEnergyMax;
    }
    
    vector< double > a;
    
    return a;
}


vector< double > VStarCatalogue::getStarFlux( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fFlux;
    }
    
    vector< double > a;
    
    return a;
}


vector< double > VStarCatalogue::getStarFluxError( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fFluxError;
    }
    
    vector< double > a;
    
    return a;
}

string VStarCatalogue::getStarType( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fType;
    }
    
    string a;
    
    return a;
}


vector< string > VStarCatalogue::getStarOtherNames( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fOtherNames;
    }
    
    vector< string > a;
    
    return a;
}

vector< string > VStarCatalogue::getStarAssociations( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fAssociations;
    }
    
    vector< string > a;
    
    return a;
}

double VStarCatalogue::getStarSpectralIndex( unsigned int iID )
{
    if( iID < fStars.size() && fStars[iID] )
    {
        return fStars[iID]->fSpectralIndex;
    }
    
    return 0.;
}

unsigned int VStarCatalogue::setFOV( string ra_hour, string dec, double FOV_x, double FOV_y, bool bJ2000 )
{
    istringstream is_stream( ra_hour );
    string temp2;
    
    double d_tt = 0.;
    is_stream >> temp2;
    d_tt += atof( temp2.c_str() );
    if( !( is_stream >> std::ws ).eof() )
    {
        is_stream >> temp2;
        d_tt += atof( temp2.c_str() ) / 60.;
    }
    if( !( is_stream >> std::ws ).eof() )
    {
        is_stream >> temp2;
        d_tt += atof( temp2.c_str() ) / 3600.;
    }
    double iSkyMapCentreRAJ2000 = d_tt / 24. * 360.;
    
    istringstream is_dec( dec );
    
    d_tt = 0.;
    is_dec >> temp2;
    d_tt += atof( temp2.c_str() );
    if( !( is_dec >> std::ws ).eof() )
    {
        is_dec >> temp2;
        d_tt += atof( temp2.c_str() ) / 60.;
    }
    if( !( is_dec >> std::ws ).eof() )
    {
        is_dec >> temp2;
        d_tt += atof( temp2.c_str() ) / 3600.;
    }
    double iSkyMapCentreDecJ2000 = d_tt;
    
    cout << "FOV Centre in deg: " << iSkyMapCentreRAJ2000 << "\t" << iSkyMapCentreDecJ2000 << endl;
    
    return setFOV( iSkyMapCentreRAJ2000, iSkyMapCentreDecJ2000, FOV_x, FOV_y, bJ2000 );
}


/*!

    loop over the current star catalogue and fill a list with stars in a box

    centred around ra/dec with width iFOV_x/iFOV_y

    all angles in [deg]
*/
unsigned int VStarCatalogue::setFOV( double ra, double dec, double iFOV_x, double iFOV_y, bool bJ2000, double iBrightness, string iBand )
{
    double degrad = 180. / TMath::Pi();
    
    fStarsinFOV.clear();
    
    double iRA = 0.;
    double iDec = 0.;
    
    for( unsigned int i = 0; i < fStars.size(); i++ )
    {
        if( iBand == "B" && fStars[i]->fBrightness_B > iBrightness )
        {
            continue;
        }
        else if( iBand == "B" && fStars[i]->fBrightness_V > iBrightness )
        {
            continue;
        }
        
        if( bJ2000 )
        {
            iRA = fStars[i]->fRA2000;
            iDec = fStars[i]->fDec2000;
        }
        else
        {
            iRA = fStars[i]->fRACurrentEpoch;
            iDec = fStars[i]->fDecCurrentEpoch;
        }
        
        double x = 0.;
        double y = 0.;
        int ierr = 0;
        
        VAstronometry::vlaDs2tp( iRA / degrad, iDec / degrad, ra / degrad, dec / degrad, &x, &y, &ierr );
        
        x *= degrad;
        y *= degrad;
        
        if( ierr == 0 && fabs( y ) < iFOV_y )
        {
            if( fabs( x ) < iFOV_x )
            {
                fStarsinFOV.push_back( fStars[i] );
            }
        }
    }
    return fStarsinFOV.size();
}


void VStarCatalogue::purge()
{
    fStars.clear();
    fStars.swap( fStars );
}


bool VStarCatalogue::writeCatalogueToRootFile( string iRootFile )
{
    TFile fOut( iRootFile.c_str(), "RECREATE" );
    if( fOut.IsZombie() )
    {
        cout << "VStarCatalogue::writeCatalogueToRootFile error opening root file: " << fOut.GetName() << endl;
        return false;
    }
    TTree* tCat = new TTree( "tCat", "star catalogue" );
    
    unsigned int fStarID = 0;
    Char_t fStarName[300];
    double fDec2000 = 0.;
    double fRA2000 = 0.;
    double fRunGalLong1958 = 0.;
    double fRunGalLat1958 = 0.;
    double fBrightness_V = 0.;
    double fBrightness_B = 0.;
    double fMajorDiameter = 0.;                        // this is either the source diameter or the possitional error
    double fMinorDiameter = 0.;
    double fPositionAngle = 0.;
    double fMajorDiameter_68 = 0.;
    double fMinorDiameter_68 = 0.;
    double fPositionAngle_68 = 0.;
    double fSignificance = 0.;
    double fSpectralIndex = 0.;
    double fSpectralIndexError = 0.;
    const unsigned int fMaxEnergyBins = 100;      /// way too many
    unsigned int fFluxEnergyBins = 0;
    double fFluxEnergyMin[fMaxEnergyBins];
    double fFluxEnergyMax[fMaxEnergyBins];
    double fFlux[fMaxEnergyBins];
    double fFluxError[fMaxEnergyBins];
    double fCutOff_MeV = 0;
    double fCutOffError_MeV = 0;
    for( unsigned int i = 0; i < fMaxEnergyBins; i++ )
    {
        fFluxEnergyMin[i] = 0.;
        fFluxEnergyMax[i] = 0.;
        fFlux[i] = 0.;
        fFluxError[i] = 0.;
    }
    unsigned int fVariability = 0;
    Char_t fStarType[300];
    int fQualityFlag = 0;
    
    tCat->Branch( "StarID", &fStarID, "StarID/i" );
    tCat->Branch( "StarName", &fStarName, "StarName/C" );
    tCat->Branch( "Dec2000", &fDec2000, "Dec2000/D" );
    tCat->Branch( "RA2000", &fRA2000, "RA2000/D" );
    tCat->Branch( "GalLong1958", &fRunGalLong1958, "GalLong1958/D" );
    tCat->Branch( "GalLat1958", &fRunGalLat1958, "GalLat1958/D" );
    tCat->Branch( "Brightness_V", &fBrightness_V, "Brightness_V/D" );
    tCat->Branch( "Brightness_B", &fBrightness_B, "Brightness_B/D" );
    tCat->Branch( "MajorDiameter", &fMajorDiameter, "MajorDiameter/D" );
    tCat->Branch( "MinorDiameter", &fMinorDiameter, "MinorDiameter/D" );
    tCat->Branch( "PositionAngle", &fPositionAngle, "PositionAngle/D" );
    tCat->Branch( "MajorDiameter_68", &fMajorDiameter_68, "MajorDiameter_68/D" );
    tCat->Branch( "MinorDiameter_68", &fMinorDiameter_68, "MinorDiameter_68/D" );
    tCat->Branch( "PositionAngle_68", &fPositionAngle_68, "PositionAngle_68/D" );
    tCat->Branch( "Significance", &fSignificance, "Significance/D" );
    tCat->Branch( "SpectralIndex", &fSpectralIndex, "SpectralIndex/D" );
    tCat->Branch( "SpectralIndexError", &fSpectralIndexError, "SpectralIndexError/D" );
    tCat->Branch( "NFluxEnergyBins", &fFluxEnergyBins, "NFluxEnergyBins/i" );
    tCat->Branch( "FluxEnergyMin", fFluxEnergyMin, "FluxEnergyMin[NFluxEnergyBins]/D" );
    tCat->Branch( "FluxEnergyMax", fFluxEnergyMax, "FluxEnergyMax[NFluxEnergyBins]/D" );
    tCat->Branch( "Flux", fFlux, "Flux[NFluxEnergyBins]/D" );
    tCat->Branch( "FluxError", fFluxError, "FluxError[NFluxEnergyBins]/D" );
    tCat->Branch( "Variability", &fVariability, "Variability/i" );
    tCat->Branch( "CutOff_MeV", &fCutOff_MeV, "CutOff_MeV/D" );
    tCat->Branch( "CutOffError_MeV", &fCutOffError_MeV, "CutOffError_MeV/D" );
    tCat->Branch( "Class",  &fStarType, "Class/C" );
    tCat->Branch( "QualityFlag", &fQualityFlag, "QualityFlag/I" );
    
    // fill tree
    for( unsigned int i = 0; i < fStars.size(); i++ )
    {
        fStarID = fStars[i]->fStarID;
        sprintf( fStarName, "%s", fStars[i]->fStarName.c_str() );
        fRA2000 = fStars[i]->fRA2000;
        fDec2000 = fStars[i]->fDec2000;
        fRunGalLong1958 = fStars[i]->fRunGalLong1958;
        fRunGalLat1958 = fStars[i]->fRunGalLat1958;
        fBrightness_V = fStars[i]->fBrightness_V;
        fBrightness_B = fStars[i]->fBrightness_B;
        fMajorDiameter = fStars[i]->fMajorDiameter;
        fMinorDiameter = fStars[i]->fMinorDiameter;
        fPositionAngle = fStars[i]->fPositionAngle;
        fMajorDiameter_68 = fStars[i]->fMajorDiameter_68;
        fMinorDiameter_68 = fStars[i]->fMinorDiameter_68;
        fPositionAngle_68 = fStars[i]->fPositionAngle_68;
        fSignificance = fStars[i]->fSignificance;
        fSpectralIndex = fStars[i]->fSpectralIndex;
        fSpectralIndexError = fStars[i]->fSpectralIndexError;
        fQualityFlag = fStars[i]->fQualityFlag;
        fVariability = fStars[i]->fVariability;
        fCutOff_MeV  = fStars[i]->fCutOff_MeV;
        fCutOffError_MeV  = fStars[i]->fCutOffError_MeV;
        
        if( fStars[i]->fFluxEnergyMin.size() > 0 && fStars[i]->fFluxEnergyMin.size() == fStars[i]->fFluxEnergyMax.size() && fStars[i]->fFlux.size() == fStars[i]->fFluxEnergyMin.size() && fStars[i]->fFluxError.size() == fStars[i]->fFluxEnergyMin.size() )
        {
            fFluxEnergyBins = fStars[i]->fFluxEnergyMin.size();
            for( unsigned int j = 0; j < fFluxEnergyBins; j++ )
            {
                fFluxEnergyMin[j] = fStars[i]->fFluxEnergyMin[j];
                fFluxEnergyMax[j] = fStars[i]->fFluxEnergyMax[j];
                fFlux[j] = fStars[i]->fFlux[j];
                fFluxError[j] = fStars[i]->fFluxError[j];
            }
        }
        sprintf( fStarType, "%s", fStars[i]->fType.c_str() );
        tCat->Fill();
    }
    
    cout << "writing tree with " << tCat->GetEntries() << " entries to " << fOut.GetName() << endl;
    tCat->Write();
    
    fOut.Close();
    
    return true;
}

VStar* VStarCatalogue::getStar( unsigned int ID )
{
    if( ID < fStars.size() )
    {
        return fStars[ID];
    }
    
    return 0;
}

bool VStarCatalogue::checkTextBlocks( string iL, unsigned int iV )
{
    unsigned int iTB = VUtilities::count_number_of_textblocks( iL );
    
    // check for correct number of text blocks
    // e.g. Hipparcos_MAG7_1997.dat
    if( iV == 10 && iTB != 5 )
    {
        return false;
    }
    // e.g. BrightStarCatalogue.txt
    else if( iV == 0 && iTB != 9 )
    {
        return false;
    }
    
    return true;
}

/*

   set telescope pointing and calculate position of stars in FOV

*/

void VStarCatalogue::setTelescopePointing( unsigned int iTelID, double iDerotationAngle, double ra_deg, double dec_deg, double iCameraScale )
{
    fTel_telescopeID = iTelID;
    fTel_deRotationAngle_deg = iDerotationAngle;
    fTel_ra          = ra_deg;
    fTel_dec         = dec_deg;
    fTel_camerascale = iCameraScale;
}

/*

    get angular distance between a bright star in the FOV and a x,y position in the camera

*/
double VStarCatalogue::getDistanceToClosestStar( double x_cam_deg, double y_cam_deg )
{
    double x_rot = 0.;
    double y_rot = 0.;
    
    double i_minDist = 1.e20;
    
    // loop over all stars in the FOV
    for( unsigned int i = 0; i < fStarsinFOV.size(); i++ )
    {
        double y = -1. * ( fStarsinFOV[i]->fDecCurrentEpoch - fTel_dec );
        double x = 0.;
        if( cos( fTel_dec * TMath::DegToRad() ) != 0. )
        {
            x = -1. * ( fStarsinFOV[i]->fRACurrentEpoch - fTel_ra ) * cos( fTel_dec * TMath::DegToRad() );
        }
        x_rot = x;
        y_rot = y;
        // derotation
        VSkyCoordinatesUtilities::rotate( -1.*fTel_deRotationAngle_deg * TMath::DegToRad(), x_rot, y_rot );
        x_rot *= -1. * fTel_camerascale;
        y_rot *= fTel_camerascale;
        
        double i_dist = sqrt( ( x_cam_deg - x_rot ) * ( x_cam_deg - x_rot ) + ( y_cam_deg - y_rot ) * ( y_cam_deg - y_rot ) );
        
        if( i_dist < i_minDist )
        {
            i_minDist = i_dist;
        }
    }
    
    return i_minDist;
}
