/*! \brief VAtmosphereSoundings read and analyse sounding data


//////////////////////////////////////////////////////////////////////////////////////

   reading ascii files with sounding data for Tucson

   Data from: http://weather.uwyo.edu/upperair/sounding.html

    VAtmosphereSoundings a;
    a.readSoundingsFromTextFile("allBalloon.txt");
    a.writeRootFile("ballonDataTucson_199501_201007.root");

*/

#include "VAtmosphereSoundings.h"

VAtmosphereSoundings::VAtmosphereSoundings()
{
    bDebug = false;
    plotAttributes_ColorChange();
    plotAttributes_PlotLegend();
    setPlottingRelativePlots();
    setHeights();
    setModtranHeights();
    setGeographicPosition();
    
    fDataFile = 0;
    fDataTree = 0;
    
    //////////////////////////////////////////////////////////////////////////////////
    // hardwired search strings for reading of sounding text files (preliminary)
    fTXTSearch_DataString = "72274 TUS Tucson Observations at";
    //////////////////////////////////////////////////////////////////////////////////
}

void VAtmosphereSoundings::setModtranHeights()
{
    //set heights at which modtran atmospheres should be sampled. (i.e., are available)
    
    fModtranHeights.clear();
    int step = 1000;
    int hmin = 0;
    int hmax = 120000;
    for( int h = hmin; h <= hmax; h += step )
    {
        fModtranHeights.push_back( ( double )h );
        if( h == 25000 )
        {
            step = 2500;
        }
        if( h == 50000 )
        {
            step = 5000;
        }
    }
}

VAtmosphereSoundings::VAtmosphereSoundings( string iRootFile, unsigned int npoints_min )
{
    bDebug = false;
    plotAttributes_ColorChange();
    plotAttributes_PlotLegend();
    setPlottingRelativePlots();
    setHeights();
    setPlottingPeriod();
    
    if( !readSoundingsFromRootFile( iRootFile, npoints_min ) )
    {
        cout << "VAtmosphereSoundings::VAtmosphereSoundings: could not read file " << iRootFile << endl;
        return;
    }
}

/*
 * read sounding data from root file
 *
 */
bool VAtmosphereSoundings::readSoundingsFromRootFile(
    string iRootFile, unsigned int npoints_min )
{
    fDataFile = new TFile( iRootFile.c_str() );
    if( fDataFile->IsZombie() )
    {
        cout << "VAtmosphereSoundings::VAtmosphereSoundings( string ) error: data file not found: " << iRootFile << endl;
        return false;
    }
    fDataTree = ( TTree* )fDataFile->Get( "tSoundings" );
    if( !fDataTree )
    {
        cout << "VAtmosphereSoundings::VAtmosphereSoundings( string ) error: data tree not found in " << iRootFile << endl;
        return false;
    }
    
    readRootFile( npoints_min );
    
    return true;
}

/*
 * read text files with sounding data downloaded from Wyoming server
 *
 * - requires list of files as input
 *
 *   highly adapted to the file format of these files
 *
*/
bool VAtmosphereSoundings::readSoundingsFromTextFile( string iFileList )
{
    ifstream is;
    is.open( iFileList.c_str() );
    if( !is )
    {
        cout << "VAtmosphereSoundings::readSoundingsFromTextFile: error: file list not found: " << iFileList << endl;
        return false;
    }
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        fListTextFile.push_back( is_line );
    }
    is.close();
    
    cout << "number of files in file list: " << fListTextFile.size() << endl;
    
    // open all files and read sounding data
    string iTemp;
    for( unsigned int i = 0; i < fListTextFile.size(); i++ )
    {
        cout << "reading " << fListTextFile[i] << endl;
        
        ifstream is;
        is.open( fListTextFile[i].c_str() );
        if( !is )
        {
            cout << "VAtmosphereSoundings::readSoundingsFromTextFile: error: sounding file not found: " << fListTextFile[i] << endl;
            continue;
        }
        int z = 0;
        while( getline( is, is_line ) )
        {
            if( is_line.size() <= 0 )
            {
                continue;
            }
            
            ///////////////////////// search for a valid date
            if( is_line.find( fTXTSearch_DataString ) < is_line.size() )
            {
                istringstream is_stream( is_line.substr( is_line.find( fTXTSearch_DataString )
                                         + fTXTSearch_DataString.size(), fTXTSearch_DataString.size() ) );
                                         
                is_stream >> iTemp;
                int iHour = atoi( iTemp.substr( 0, 2 ).c_str() );
                is_stream >> iTemp;
                int iDay  = atoi( iTemp.c_str() );
                is_stream >> iTemp;
                int iMonth = getMonth( iTemp );
                is_stream >> iTemp;
                int iYear = atoi( iTemp.c_str() );
                
                // calculate MJD from date
                double iMJD = 0;
                int j = 0;
                VAstronometry::vlaCldj( iYear, iMonth, iDay, &iMJD, &j );
                if( j != 0 )
                {
                    cout << "VAtmosphereSoundings::readSoundingsFromTextFile: error: invalid data: " << is_line << endl;
                    continue;
                }
                if( iHour == 0 )
                {
                    iMJD += 0.5;
                }
                
                // a valid new data starts a new entry
                fData.push_back( new VAtmosphereSoundingData() );
                fData.back()->MJD   = iMJD;
                fData.back()->Year  = iYear;
                fData.back()->Month = iMonth;
                fData.back()->Day   = iDay;
                fData.back()->Hour  = ( double )iHour;
                z = 0;
                ///////////////////////// read a full entry
                // expect first entry 7 lines after date
                while( getline( is, is_line ) )
                {
                    if( is_line.size() <= 0 )
                    {
                        continue;
                    }
                    
                    if( is_line.find( "Station information" ) < is_line.size() )
                    {
                        break;
                    }
                    
                    if( z >= 5 )
                    {
                        double iT = 0.;
                        istringstream is_stream( is_line );
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                        }
                        if( is_line.substr( 4, 1 ) != " " )
                        {
                            fData.back()->fPressure_Pa.push_back( iT * 100. );
                        }
                        else
                        {
                            fData.back()->fPressure_Pa.push_back( -9999. );
                        }
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                        }
                        if( is_line.substr( 13, 1 ) != " " )
                        {
                            fData.back()->fHeight_m.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fHeight_m.push_back( -9999. );
                        }
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                            iT += 273.15;
                        }
                        if( is_line.substr( 18, 1 ) != " " )
                        {
                            fData.back()->fTemperature_K.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fTemperature_K.push_back( -9999. );
                        }
                        
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                            iT += 273.15;
                        }
                        if( is_line.substr( 25, 1 ) != " " )
                        {
                            fData.back()->fDewPoint_K.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fDewPoint_K.push_back( -9999. );
                        }
                        
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                            //iT /= 1.e2;                // % to fraction HF: is this needed?
                        }
                        if( is_line.substr( 34, 1 ) != " " )
                        {
                            fData.back()->fRelativeHumidity.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fRelativeHumidity.push_back( -9999. );
                        }
                        // mixing ratio
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                        }
                        if( is_line.substr( 41, 1 ) != " " )
                        {
                            fData.back()->fMixingRatio_gkg.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fMixingRatio_gkg.push_back( -9999. );
                        }
                        // wind direction
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                        }
                        if( is_line.substr( 48, 1 ) != " " )
                        {
                            fData.back()->fWindDirection_deg.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fWindDirection_deg.push_back( -9999. );
                        }
                        // wind speed
                        iT = -9999.;
                        if( !( is_stream >> std::ws ).eof() )
                        {
                            is_stream >> iT;
                            iT *= 0.514444;  // [knots] to [m/s]
                        }
                        if( is_line.substr( 55, 1 ) != " " )
                        {
                            fData.back()->fWindSpeed_ms.push_back( iT );
                        }
                        else
                        {
                            fData.back()->fWindSpeed_ms.push_back( -9999. );
                        }
                    }
                    z++;
                }
            }
        }
    }
    // fill remaining fields
    // note: some of these calculations are approximations
    fillWaterVaporDensity();
    fillAtmosphericDensity();
    fillAtmosphericThickness();
    fillIndexofRefraction();
    fillO2();
    fillO3();
    
    make_interpolated_atmospheres();
    
    return true;
}

/*

*/
bool VAtmosphereSoundings::readGDASFromTextFile( string iFileList )
{
    ifstream is;
    is.open( iFileList.c_str() );
    if( !is )
    {
        cout << "VAtmosphereSoundings::readGDASFromTextFile: error: file list not found: " << iFileList << endl;
        return false;
    }
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        fListTextFile.push_back( is_line );
    }
    is.close();
    
    cout << "number of files in file list: " << fListTextFile.size() << endl;
    
    // open all files and read sounding data
    string iTemp;
    for( unsigned int i = 0; i < fListTextFile.size(); i++ )
    {
    
        // extract date and time from file name. Should be something like C-70.4-24.6D2016-06-30T21.gdas
        // date string should be between 'D' and 'T'. Time after 'T'.
        size_t posD = 0, posT = 0;
        if( ( posD = fListTextFile[i].find( "D" ) ) != std::string::npos && ( posT = fListTextFile[i].find( "T" ) ) != std::string::npos )
        {
            std::string sDate = fListTextFile[i].substr( posD + 1, 10 );
            int iYear = atoi( sDate.substr( 0, 4 ).c_str() );
            int iMonth = atoi( sDate.substr( 5, 2 ).c_str() );
            int iDay = atoi( sDate.substr( 8, 2 ).c_str() );
            
            std::string sTime = fListTextFile[i].substr( posT + 1, 2 );
            int iHour = atoi( sTime.c_str() );
            
            // get MJD
            double iMJD = 0;
            int j = 0;
            VAstronometry::vlaCldj( iYear, iMonth, iDay, &iMJD, &j );
            if( j != 0 )
            {
                cout << "VAtmosphereSoundings::readGDASFromTextFile: error: invalid date string, could not find MJD: " << sDate << endl;
                continue;
            }
            // preli!
            iMJD += ( double )iHour / 24.0;
            
            // a valid new data starts a new entry
            fData.push_back( new VAtmosphereSoundingData() );
            fData.back()->MJD   = iMJD;
            fData.back()->Year  = iYear;
            fData.back()->Month = iMonth;
            fData.back()->Day   = iDay;
            fData.back()->Hour  = ( double )iHour;
            
            
        }
        else
        {
            cout << "VAtmosphereSoundings::readGDASFromTextFile: error: could not extract date and time from file name: " << fListTextFile[i] << endl;
            continue;
        }
        //cout << "reading " << fListTextFile[i] << endl;
        
        ifstream is;
        is.open( fListTextFile[i].c_str() );
        if( !is )
        {
            cout << "VAtmosphereSoundings::readGDASFromTextFile: error: sounding file not found: " << fListTextFile[i] << endl;
            continue;
        }
        while( getline( is, is_line ) )
        {
            if( is_line.size() <= 0 || is_line.find( "#" ) != std::string::npos )
            {
                continue;
            }
            
            istringstream is_stream( is_line );
            double iP, iH, iT, iRH ;
            //GDAS format: # P[hPa] HGT[m]         T[K]    RELHUM[%]
            while( is_stream >> iP >> iH >> iT >> iRH )
            {
                fData.back()->fPressure_Pa.push_back( iP * 100. );
                fData.back()->fHeight_m.push_back( iH );
                fData.back()->fTemperature_K.push_back( iT );
                fData.back()->fRelativeHumidity.push_back( iRH ); //HF should it be in % or not?
                
                fData.back()->fDewPoint_K.push_back( -9999. );
                fData.back()->fMixingRatio_gkg.push_back( -9999. );
                fData.back()->fWindDirection_deg.push_back( -9999. );
                fData.back()->fWindSpeed_ms.push_back( -9999. );
            }
        }
    }
    // fill remaining fields
    fillWaterVaporDensity();
    fillAtmosphericDensity();
    fillAtmosphericThickness();
    fillIndexofRefraction();
    fillO2();
    fillO3();
    
    make_interpolated_atmospheres();
    
    
    return true;
}

void VAtmosphereSoundings::make_interpolated_atmospheres()
{
    fDataInterpol.clear();
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fDataInterpol.push_back( make_interpolated_atmosphere( fData.at( i ) ) ) ;
    }
}


void VAtmosphereSoundings::fillO2()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
            for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
            {
                fData[i]->fO2_cmkm.push_back( -9999. );
            }
        }
    }
}

void VAtmosphereSoundings::fillO3()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
            for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
            {
                fData[i]->fO3_cmkm.push_back( -9999. );
            }
        }
    }
}



/*void VAtmosphereSoundings::fillAtmosphericThickness()
{
   for( unsigned int i = 0; i < fData.size(); i++ )
   {
      if( fData[i] )
      {
         for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
	 {
	    fData[i]->fThickness_gcm2.push_back( -9999. );
         }
       }
   }
} */

void VAtmosphereSoundings::fillAtmosphericPressure()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fillAtmosphericPressure( fData[i] );
    }
}

/*
   calculate atmospheric pressure from atmospheric thickness
*/
void VAtmosphereSoundings::fillAtmosphericPressure( VAtmosphereSoundingData* iData )
{
    if( !iData )
    {
        return;
    }
    
    // calculate graviational constant for the given observatory latitude and height
    double phi = fObservatoryLatitude * atan( 1. ) / 45.;
    double g = 9.780327 * ( 1. + 0.0053024 * ( sin( phi ) * sin( phi ) ) - 0.0000058 * sin( 2.*phi ) * sin( 2.*phi ) ) - 3.086e-6 * fObservatoryHeight_km * 1.e3;
    
    double p = -9999;
    for( unsigned int j = 0; j < iData->fThickness_gcm2.size(); j++ )
    {
        if( iData->fThickness_gcm2[j] > 0. )
        {
            p  = g * iData->fThickness_gcm2[j];
            p *= 10.;                                     // kg/m2 -> m/cm2
        }
        else
        {
            p  = -9999.;
        }
        if( j >= iData->
                fPressure_Pa.size() )
        {
            iData->fPressure_Pa.push_back( p );
        }
        else
        {
            iData->fPressure_Pa[j] = p;
        }
    }
}

void VAtmosphereSoundings::fillAtmosphericThickness()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fillAtmosphericThickness( fData[i] );
    }
}

void VAtmosphereSoundings::fillAtmosphericThickness( VAtmosphereSoundingData* iData )
{
    if( !iData )
    {
        return;
    }
    
    // calculate graviational constant
    double phi = fObservatoryLatitude * atan( 1. ) / 45.;
    double g = 9.780327 * ( 1. + 0.0053024 * ( sin( phi ) * sin( phi ) ) - 0.0000058 * sin( 2.*phi ) * sin( 2.*phi ) ) - 3.086e-6 * fObservatoryHeight_km * 1.e3;
    
    double t = -9999;
    for( unsigned int j = 0; j < iData->fPressure_Pa.size(); j++ )
    {
        if( iData->fPressure_Pa[j] > 0. )
        {
            t  = iData->fPressure_Pa[j] / g;
            t /= 10.;                                     // kg/m2 -> m/cm2
        }
        else
        {
            t  = -9999.;
        }
        if( j >= iData->fThickness_gcm2.size() )
        {
            iData->fThickness_gcm2.push_back( t );
        }
        else
        {
            iData->fThickness_gcm2[j] = t;
        }
    }
}

void VAtmosphereSoundings::fillIndexofRefraction( )
{
    // Refractive index formula from J.Owens, 'Optical Refractive Index of Air' (1967)
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
        
            for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
            {
                double pressurePa = fData[i]->fPressure_Pa.at( j );
                double temperatureK = fData[i]->fTemperature_K.at( j );
                double relhum = fData[i]->fRelativeHumidity.at( j );
                
                //convert units
                double p = pressurePa * 7.5006e-3; //pressure in Torricelli
                double t = temperatureK - 273.15; //temperature in degrees celcius
                double r = relhum / 100.; //relative humidity as a fraction
                
                //set wavenumber
                double wavelength = 0.3; //wavelength in micrometers
                double wavenumber = 1 / wavelength;
                
                //evaluate refractive index expression
                double A = 8342.13 + ( 2406030 / ( 130 - pow( wavenumber, 2 ) ) ) + ( 15997 / ( 38.9 - pow( wavenumber, 2 ) ) );
                double B = ( p / 720.775 ) * ( ( 1 + p * ( 0.817 - 0.0133 * t ) * 1e-6 ) / ( 1 + 0.0036610 * t ) );
                
                //parameters and formula for partial water pressure - Antoine eqn
                double a, b, c = 0;
                if( t < 100 )
                {
                    a = 8.07131;
                    b = 1730.63;
                    c = 233.426;
                }
                else
                {
                    a = 8.14019;
                    b = 1810.94;
                    c = 244.485;
                }
                
                double x = a - ( b / ( c + t ) );
                double f_mmHg = r * pow( 10, x ); //partial water pressure in mmHg
                double f = f_mmHg * ( 101.325 / 760. ) * 1e-3 * 7.5006e-3; //convert units to Torricelli
                double C = f * ( 5.722 - 0.0457 * pow( wavenumber, 2 ) );
                double index = A * B - C;
                index = index * 1e-8 + 1;
                
                fData[i]->fIndexofRefraction.push_back( index );
            }
        }
        
    }
}


/*

    read a CORSIKA atmosphere from a text file (atmprof files)

*/
int VAtmosphereSoundings::read_CORSIKA_Atmosphere( string iFile, string iName, int iColor, int iLineStyle )
{
    ifstream is;
    char* iFileName = gSystem->ExpandPathName( iFile.c_str() );
    is.open( iFileName );
    delete iFileName;
    if( !is )
    {
        cout << "VAtmosphereSoundings::read_CORSIKA_Atmosphere: error opening CORSIKA atmospheric file " << iFile << endl;
        return -1;
    }
    
    fDataCORSIKAMODTRAN.push_back( new VAtmosphereSoundingData() );
    fDataCORSIKAMODTRAN.back()->Name = iName;
    fDataCORSIKAMODTRAN.back()->PlotColor = iColor;
    fDataCORSIKAMODTRAN.back()->PlotLineStyle = iLineStyle;
    
    string is_line;
    string iTemp;
    double d;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        if( is_line.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        fDataCORSIKAMODTRAN.back()->fPressure_Pa.push_back( -9999. );
        
        is_stream >> iTemp;
        d = atof( iTemp.c_str() );
        fDataCORSIKAMODTRAN.back()->fHeight_m.push_back( d * 1.e3 );
        
        is_stream >> iTemp;
        d = atof( iTemp.c_str() );
        fDataCORSIKAMODTRAN.back()->fDensity_gcm3.push_back( d );
        
        is_stream >> iTemp;
        d = atof( iTemp.c_str() );
        fDataCORSIKAMODTRAN.back()->fThickness_gcm2.push_back( d );
        
        fDataCORSIKAMODTRAN.back()->fTemperature_K.push_back( -9999. );
        
        fDataCORSIKAMODTRAN.back()->fDewPoint_K.push_back( -9999. );
        
        fDataCORSIKAMODTRAN.back()->fRelativeHumidity.push_back( -9999. );
        fDataCORSIKAMODTRAN.back()->fVaporMassDensity_gm3.push_back( -9999. );
        
        is_stream >> iTemp;
        d = atof( iTemp.c_str() );
        fDataCORSIKAMODTRAN.back()->fIndexofRefraction.push_back( d + 1. );
        
        fDataCORSIKAMODTRAN.back()->fMixingRatio_gkg.push_back( -9999. );
        fDataCORSIKAMODTRAN.back()->fWindDirection_deg.push_back( -9999. );
        fDataCORSIKAMODTRAN.back()->fWindSpeed_ms.push_back( -9999. );
        fDataCORSIKAMODTRAN.back()->fO2_cmkm.push_back( -9999. );
        fDataCORSIKAMODTRAN.back()->fO3_cmkm.push_back( -9999. );
        
    }
    is.close();
    
    fillAtmosphericPressure( fDataCORSIKAMODTRAN.back() );
    
    return fDataCORSIKAMODTRAN.size() - 1 ;
}

/*

    read MODTRAN tp6 files

*/
int VAtmosphereSoundings::read_MODTRAN_Atmosphere( string iFile, string iName, int iColor, int iLineStyle )
{
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "VAtmosphereSoundings::read_MODTRAN_Atmosphere: error opening MODTRAN tp6 file" << iFile << endl;
        return -1;
    }
    // for now: only tp6 files can be read
    if( iFile.find( "tp6" ) >= iFile.size() )
    {
        return false;
    }
    
    fDataCORSIKAMODTRAN.push_back( new VAtmosphereSoundingData() );
    fDataCORSIKAMODTRAN.back()->Name = iName;
    fDataCORSIKAMODTRAN.back()->PlotColor = iColor;
    fDataCORSIKAMODTRAN.back()->PlotLineStyle = iLineStyle;
    
    cout << "reading MODTRAN atmosphere from " << iFile << endl;
    
    string is_line;
    string iTemp;
    
    int z = 0; 		//internal counter
    
    int found_profile = 0;
    
    //reading in tp6 file line by line
    //atmospheric profiles (including default values where we didn't give anything) start with a line containing  "ATMOSPHERIC PROFILES")
    
    while( getline( is, is_line ) )
    {
        if( is_line.length() == 0 )
        {
            continue;
        }
        if( is_line.find( "ATMOSPHERIC PROFILES" ) != std::string::npos )
        {
            found_profile++;    //we found one! reset counter.
            z = 1;
            continue;
        }
        
        if( found_profile == 1 )   	//first table... read in pressure etc.
        {
            istringstream is_stream( is_line );
            is_stream >> iTemp;				//attempt to read in line.
            if( atoi( iTemp.c_str() ) != z )
            {
                continue;
            }
            
            //found next line of atm profile! read in values.
            z++;
            // height (km)
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fHeight_m.push_back( atof( iTemp.c_str() ) * 1.e3 );
            // pressure (MB)
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fPressure_Pa.push_back( atof( iTemp.c_str() ) * 1.e2 );
            // temperature (K)
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fTemperature_K.push_back( atof( iTemp.c_str() ) );
            // density
            fDataCORSIKAMODTRAN.back()->fDensity_gcm3.push_back( -9999. );
            fDataCORSIKAMODTRAN.back()->fThickness_gcm2.push_back( -9999. );
            fDataCORSIKAMODTRAN.back()->fDewPoint_K.push_back( -9999. );
            fDataCORSIKAMODTRAN.back()->fMixingRatio_gkg.push_back( -9999. );
            fDataCORSIKAMODTRAN.back()->fWindDirection_deg.push_back( -9999. );
            fDataCORSIKAMODTRAN.back()->fWindSpeed_ms.push_back( -9999. );
            // N2
            is_stream >> iTemp;
            // CNTMSLF
            is_stream >> iTemp;
            // CNTMFRN
            is_stream >> iTemp;
            // MOL SCAT
            is_stream >> iTemp;
            // N-1
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fIndexofRefraction.push_back( atof( iTemp.c_str() ) + 1. );
            // O3 (UV)
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fO3_cmkm.push_back( atof( iTemp.c_str() ) );
            // O2 (UV)
            is_stream >> iTemp;
            fDataCORSIKAMODTRAN.back()->fO2_cmkm.push_back( atof( iTemp.c_str() ) );
            
        }
        else if( found_profile == 2 ) //second table... needed for relative humidity.
        {
            istringstream is_stream( is_line );
            is_stream >> iTemp;				//attempt to read in line.
            if( atoi( iTemp.c_str() ) != z )
            {
                continue;
            }
            
            
            z++;
            // relative humidity
            for( unsigned int s = 0; s < 9; s++ )
            {
                is_stream >> iTemp;
            }
            fDataCORSIKAMODTRAN.back()->fRelativeHumidity.push_back( atof( iTemp.c_str() ) );
        }
        
        else if( found_profile > 2 )
        {
            break;
        }
        
    }
    is.close();
    
    for( unsigned int i = 0; i < fDataCORSIKAMODTRAN.size(); i++ )
    {
        if( fDataCORSIKAMODTRAN[i]->fIndexofRefraction.size() != fDataCORSIKAMODTRAN[i]->fRelativeHumidity.size() )
        {
            for( unsigned int j = 0; j < fDataCORSIKAMODTRAN[i]->fIndexofRefraction.size(); j++ )
            {
                fDataCORSIKAMODTRAN[i]->fRelativeHumidity.push_back( -9999. );
            }
        }
    }
    fillWaterVaporDensity( fDataCORSIKAMODTRAN.back() );
    fillAtmosphericDensity( fDataCORSIKAMODTRAN.back() );
    fillAtmosphericThickness( fDataCORSIKAMODTRAN.back() );
    
    return fDataCORSIKAMODTRAN.size() - 1 ;
}

void VAtmosphereSoundings::fillWaterVaporDensity()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fillWaterVaporDensity( fData[i] );
    }
}

void VAtmosphereSoundings::fillWaterVaporDensity( VAtmosphereSoundingData* iData )
{
    if( !iData )
    {
        return;
    }
    
    for( unsigned int i = 0; i < iData->fTemperature_K.size(); i++ )
    {
        double w = -9999.;
        if( i <= iData->fRelativeHumidity.size() )
        {
            w = getWaterVaporDensity( iData->fTemperature_K[i], iData->fRelativeHumidity[i] );
        }
        
        if( i < iData->fVaporMassDensity_gm3.size() )
        {
            iData->fVaporMassDensity_gm3[i] = w;
        }
        else
        {
            iData->fVaporMassDensity_gm3.push_back( w );
        }
    }
}

void VAtmosphereSoundings::fillAtmosphericDensity()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fillAtmosphericDensity( fData[i] );
    }
}


/*!
   calculate atmosheric density

   after meteorology textbooks (see e.g. http://wahiduddin.net/calc/density_altitude.htm )

*/
void VAtmosphereSoundings::fillAtmosphericDensity( VAtmosphereSoundingData* iData )
{
    double R_d = 287.05;        //  Rd = gas constant for dry air, J/(kg*degK)
    
    if( iData )
    {
        for( unsigned int j = 0; j < iData->fPressure_Pa.size(); j++ )
        {
            double d = 0.;
            double pressure =    iData->fPressure_Pa[j];
            double temperature = -9999;
            if( j < iData->fTemperature_K.size() )
            {
                temperature = iData->fTemperature_K[j];
            }
            double relativeHumidity = -9999.;
            if( j < iData->fRelativeHumidity.size() )
            {
                relativeHumidity = iData->fRelativeHumidity[j];
            }
            double dewpoint    = -9999;
            // calculate dew point: if not available
            if( j < iData->fDewPoint_K.size() && iData->fDewPoint_K[j] > -300. )
            {
                dewpoint = iData->fDewPoint_K[j];
            }
            else
            {
                dewpoint = getDewPoint( temperature, relativeHumidity, 1 );
                if( j < iData->fDewPoint_K.size() )
                {
                    iData->fDewPoint_K[j] = dewpoint;
                }
                else
                {
                    iData->fDewPoint_K.push_back( dewpoint );
                }
            }
            
            if( pressure < -9998. || temperature < -0. || dewpoint < -9998. )
            {
                d = -9999.;
            }
            else
            {
                // dry air
                if( TMath::Abs( relativeHumidity ) < 1. ) //e-2 )
                {
                    d = ( pressure / R_d / temperature );
                }
                // wet air
                else
                {
                    double pressure_vapore = getAmosphericVaporPressure( dewpoint );
                    
                    if( pressure > 0. )
                    {
                        d  = ( pressure / R_d / temperature ) * ( 1. - 0.378 * pressure_vapore / pressure );                // kg/m3
                    }
                    else
                    {
                        d = 0.;
                    }
                }
                d *= 1.e3 / 1.e6;  // kg/m3 -> g/cm3
            }
            if( j < iData->fDensity_gcm3.size() )
            {
                iData->fDensity_gcm3[j] = d;
            }
            else
            {
                iData->fDensity_gcm3.push_back( d );
            }
        }
    }
}

/*!

   calculate water vapor mass density [g/m3] for a given humidity

   T in [K]

   from MODTRAN watvap.f
*/
double VAtmosphereSoundings::getWaterVaporDensity( double T, double RH )
{
    if( T < 0. || RH < 0. )
    {
        return -9999.;
    }
    
    double AVOGAD =  TMath::Na();   // Avogadro number
    double AMWT = 18.015;           // molecular weight of water
    double B = AVOGAD / AMWT;
    
    double A = 273.15 / T;
    
    double WWMOL = RH;
    
    return 1.e6 * getWaterVaporMassDensity( A ) * ( WWMOL / 100.0 ) / B;
}

/*!
   calculate saturated vapor mass density [g/m3]

   accuracy: <1% for -50 to +50 [C]

   from MODTRAN watvap.f
*/
double VAtmosphereSoundings::getWaterVaporMassDensity( double ATEMP )
{
    double C1 = 18.9766;
    double C2 = -14.9595;
    double C3 = -2.43882;
    
    double AVOGAD =  TMath::Na();   // Avogadro number
    double AMWT = 18.015;           // molecular weight of water
    double B = AVOGAD / AMWT;
    
    return ATEMP * B * TMath::Exp( C1 + C2 * ATEMP + C3 * ATEMP * ATEMP ) * 1.0e-6;
}


/*!

   calculate dew point from temperature and relative humidity

   T in [K]

   Attention: this might be in inappropriate approximation (valid only for T > 0. [C] and humidity > 1%)

*/
double VAtmosphereSoundings::getDewPoint( double T, double R, int iMethod )
{
    if( T < 1.e-2 || R < 0. )
    {
        return -9999.;
    }
    
    // not sure if this makes sense
    if( TMath::Abs( R ) < 1.e-3 )
    {
        return 0.;
    }
    
    double d = -9999.;
    
    if( iMethod == 0 )
    {
        d = TMath::Power( R / 100., 1. / 8. ) * ( 112. + 0.9 * ( T - 273.15 ) ) + 0.1 * ( T - 273.15 ) - 112.;
        d += 273.15;
    }
    // Magnus-Tetens formula
    else if( iMethod == 1 )
    {
        double alpha = 17.27 * ( T - 273.15 ) / ( 237.7 + ( T - 273.15 ) ) + log( R / 1.e2 );
        d = 237.7 * alpha / ( 17.27 - alpha );
        d += 273.15;
    }
    else if( iMethod == 2 )
        // NOAA
        //  (http://www.srh.noaa.gov/images/epz/wxcalc/rhTdFromWetBulb.pdf)
    {
        double es = 6.112 * TMath::Exp( 17.67 * ( T - 273.15 ) / ( T - 273.15 + 243.5 ) );
        double e  = R / 1.e2 * es;
        
        d = 243.5 * log( e / 6.112 ) / ( 17.67 - log( e / 6.112 ) );
        
        d += 273.15;
    }
    
    return d;
}

/*!
    calculate vapor pressure from temperature

    or

    calculate actual vapor pressure from dewpoint

    after Wobus
*/
double VAtmosphereSoundings::getAmosphericVaporPressure( double T )
{
    double c0 = 0.99999683;
    double c1 = -0.90826951e-2;
    double c2 = 0.78736169e-4;
    double c3 = -0.61117958e-6;
    double c4 = 0.43884187e-8;
    double c5 = -0.29883885e-10;
    double c6 = 0.21874425e-12;
    double c7 = -0.17892321e-14;
    double c8 = 0.11112018e-16;
    double c9 = -0.30994571e-19;
    double eso = 6.1078;
    
    T -= 273.15; // K -> deg
    
    double p = ( c0 + T * ( c1 + T * ( c2 + T * ( c3 + T * ( c4 + T * ( c5 + T * ( c6 + T * ( c7 + T * ( c8 + T * ( c9 ) ) ) ) ) ) ) ) ) );
    
    double Es = eso / TMath::Power( p, 8. ) * 1.e2;             // vapour pressure [mpascal]
    
    // alternative calculation (less accurate)
    //   double Es = 6.1078 * TMath::Power( 10., 7.5*T/(237.3+T) ) * 1.e2;
    
    return Es;
}

int VAtmosphereSoundings::getMonth( string iT )
{
    if( iT == "Jan" )
    {
        return 1;
    }
    if( iT == "Feb" )
    {
        return 2;
    }
    if( iT == "Mar" )
    {
        return 3;
    }
    if( iT == "Apr" )
    {
        return 4;
    }
    if( iT == "May" )
    {
        return 5;
    }
    if( iT == "Jun" )
    {
        return 6;
    }
    if( iT == "Jul" )
    {
        return 7;
    }
    if( iT == "Aug" )
    {
        return 8;
    }
    if( iT == "Sep" )
    {
        return 9;
    }
    if( iT == "Oct" )
    {
        return 10;
    }
    if( iT == "Nov" )
    {
        return 11;
    }
    if( iT == "Dec" )
    {
        return 12;
    }
    
    return 0;
}

bool VAtmosphereSoundings::writeRootFile( string iFileName )
{
    TFile iFile( iFileName.c_str(), "RECREATE" );
    if( iFile.IsZombie() )
    {
        cout << "VAtmosphereSoundings::writeRootFile error opening root file " << iFileName << endl;
        return false;
    }
    TTree* t = new TTree( "tSoundings", "soundings data" );
    t->Branch( "MJD", &MJD, "MJD/D" );
    t->Branch( "Year", &Year, "Year/I" );
    t->Branch( "Month", &Month, "Month/I" );
    t->Branch( "Day", &Day, "Day/I" );
    t->Branch( "Hour", &Hour, "Hour/D" );
    t->Branch( "nPoints", &nPoints, "nPoints/i" );
    t->Branch( "Height_m", Height_m, "Height_m[nPoints]/D" );
    t->Branch( "Pressure_Pa", Pressure_Pa, "Pressure_Pa[nPoints]/D" );
    t->Branch( "Density_gcm3", Density_gcm3, "Density_gcm3[nPoints]/D" );
    t->Branch( "Thickness_gcm2", Thickness_gcm2, "Thickness_gcm2[nPoints]/D" );
    t->Branch( "Temperature_K", Temperature_K, "Temperature_K[nPoints]/D" );
    t->Branch( "DewPoint_K", DewPoint_K, "DewPoint_K[nPoints]/D" );
    t->Branch( "RelativeHumidity", RelativeHumidity, "RelativeHumidity[nPoints]/D" );
    t->Branch( "VaporMassDensity_gm3", VaporMassDensity_gm3, "VaporMassDensity_gm3[nPoints]/D" );
    t->Branch( "MixingRatio_gkg", MixingRatio_gkg, "MixingRatio_gkg[nPoints]/D" );
    t->Branch( "WindDirection_deg", WindDirection_deg, "WindDirection_deg[nPoints]/D" );
    t->Branch( "WindSpeed_ms", WindSpeed_ms, "WindSpeed_ms[nPoints]/D" );
    t->Branch( "IndexofRefraction", IndexofRefraction, "IndexofRefraction[nPoints]/D" ); //new
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        MJD   = fData[i]->MJD;
        Year  = fData[i]->Year;
        Month = fData[i]->Month;
        Day   = fData[i]->Day;
        Hour  = fData[i]->Hour;
        
        nPoints = fData[i]->fPressure_Pa.size();
        if( nPoints >= VMAXNUMBEROFSOUNDINGPOINTS )
        {
            cout << "more than " << VMAXNUMBEROFSOUNDINGPOINTS << " points..." << endl;
            nPoints = VMAXNUMBEROFSOUNDINGPOINTS;
        }
        for( unsigned int j = 0; j < nPoints; j++ )
        {
            if( j < fData[i]->fPressure_Pa.size() )
            {
                Pressure_Pa[j]      = fData[i]->fPressure_Pa[j];
            }
            if( j < fData[i]->fHeight_m.size() )
            {
                Height_m[j]         = fData[i]->fHeight_m[j];
            }
            if( j < fData[i]->fDensity_gcm3.size() )
            {
                Density_gcm3[j]     = fData[i]->fDensity_gcm3[j];
            }
            if( j < fData[i]->fThickness_gcm2.size() )
            {
                Thickness_gcm2[j]   = fData[i]->fThickness_gcm2[j];
            }
            if( j < fData[i]->fTemperature_K.size() )
            {
                Temperature_K[j]    = fData[i]->fTemperature_K[j];
            }
            if( j < fData[i]->fDewPoint_K.size() )
            {
                DewPoint_K[j]       = fData[i]->fDewPoint_K[j];
            }
            if( j < fData[i]->fRelativeHumidity.size() )
            {
                RelativeHumidity[j] = fData[i]->fRelativeHumidity[j];
            }
            if( j < fData[i]->fVaporMassDensity_gm3.size() )
            {
                VaporMassDensity_gm3[j] = fData[i]->fVaporMassDensity_gm3[j];
            }
            if( j < fData[i]->fMixingRatio_gkg.size() )
            {
                MixingRatio_gkg[j]  = fData[i]->fMixingRatio_gkg[j];
            }
            if( j < fData[i]->fWindDirection_deg.size() )
            {
                WindDirection_deg[j] = fData[i]->fWindDirection_deg[j];
            }
            if( j < fData[i]->fWindSpeed_ms.size() )
            {
                WindSpeed_ms[j]     = fData[i]->fWindSpeed_ms[j];
            }
            if( j < fData[i]->fIndexofRefraction.size() ) //new
            {
                IndexofRefraction[j] = fData[i]->fIndexofRefraction[j];
            }
        }
        
        t->Fill();
    }
    t->Write();
    
    iFile.Close();
    
    return true;
}

unsigned int VAtmosphereSoundings::getHistogramIdentifier( unsigned int i )
{
    if( i < fData.size() )
    {
        if( fPlottingPeriod == "yearly" )
        {
            return fData[i]->Year;
        }
        else if( fPlottingPeriod == "monthly" )
        {
            return fData[i]->Year * 100 + fData[i]->Month;
        }
        else if( fPlottingPeriod == "all" )
        {
            return 100;
        }
        else if( fPlottingPeriod == "monthly_all" )
        {
            return fData[i]->Month;
        }
        else if( fPlottingPeriod == "day_night" )
        {
            if( fData[i]->Hour == 0 )
            {
                return 1;
            }
            else if( fData[i]->Hour == 12 )
            {
                return ( unsigned int )fData[i]->Hour;
            }
            return 0;
        }
        else if( fPlottingPeriod == "bimonthly_all" )
        {
            return fData[i]->Month / 2;
        }
        else if( fPlottingPeriod.find( ".dat" ) < fPlottingPeriod.size() )
        {
            return checkPlottingPeriodIdentifier( i );
        }
        // choose similar winter months
        else if( fPlottingPeriod == "winter" )
        {
            if( fData[i]->Month == 12 || fData[i]->Month <= 3 )
            {
                return 1;
            }
        }
        // choose summer months outside of monsoon season
        else if( fPlottingPeriod == "summer" )
        {
            if( fData[i]->Month == 6 || fData[i]->Month == 9 )
            {
                return 1;
            }
        }
    }
    
    return  0;
}

unsigned int VAtmosphereSoundings::checkPlottingPeriodIdentifier( unsigned int i )
{
    if( i >= fData.size() )
    {
        return 0;
    }
    
    unsigned int iDate = fData[i]->Year * 10000 + fData[i]->Month * 100 + fData[i]->Day;
    
    map< unsigned int, vector< unsigned int> >::iterator i_iter;
    for( i_iter = fPlottingPeriodDates.begin(); i_iter != fPlottingPeriodDates.end(); ++i_iter )
    {
        for( unsigned int j = 0; j < ( *i_iter ).second.size(); j++ )
        {
            if( ( *i_iter ).second[j] == iDate )
            {
                return ( *i_iter ).first;
            }
        }
    }
    
    return 0;
}

void VAtmosphereSoundings::setPlottingPeriod( string iPeriod )
{
    fPlottingPeriod = iPeriod;
    
    if( fPlottingPeriod.find( ".dat" ) < fPlottingPeriod.size() )
    {
        readPlottingPeriodsFromTextFile( fPlottingPeriod );
    }
}

bool VAtmosphereSoundings::readPlottingPeriodsFromTextFile( string iFile )
{
    if( iFile.size() == 0 )
    {
        return false;
    }
    
    fPlottingPeriodFiles.clear();
    fPlottingPeriodDates.clear();
    
    // get file names with run dates
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "VAtmosphereSoundings::readPlottingPeriodsFromTextFile(): file with plotting periods not found " << iFile << endl;
        return false;
    }
    string is_line;
    string iTemp;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> iTemp;
        unsigned int iID = ( unsigned int )atoi( iTemp.c_str() );
        is_stream >> iTemp;
        fPlottingPeriodFiles[iID] = iTemp;
    }
    is.close();
    cout << "reading plotting periods from " << fPlottingPeriodFiles.size() << " files " << endl;
    
    // get run dates from different files
    map< unsigned int, string >::iterator i_iter;
    for( i_iter = fPlottingPeriodFiles.begin(); i_iter != fPlottingPeriodFiles.end(); i_iter++ )
    {
        ifstream is;
        is.open( ( *i_iter ).second.c_str() );
        if( !is )
        {
            cout << "error opening file for plotting periods " << ( *i_iter ).first << "(" << ( *i_iter ).second << ")" << endl;
            continue;
        }
        cout << "reading plotting dates for " << ( *i_iter ).first << "(" << ( *i_iter ).second << ")" << endl;
        unsigned int iDate = 0;
        while( getline( is, is_line ) )
        {
            if( is_line.size() <= 0 )
            {
                continue;
            }
            
            istringstream is_stream( is_line );
            
            is_stream >> iTemp;
            
            iDate = ( unsigned int )atoi( iTemp.c_str() );
            
            fPlottingPeriodDates[( *i_iter ).first].push_back( iDate );
        }
        is.close();
        
        cout << "\t read " << fPlottingPeriodDates[( *i_iter ).first].size() << " days " << endl;
    }
    
    return true;
}




void VAtmosphereSoundings::plot2DProfiles( unsigned int iYearStart, unsigned int iMonthStart, unsigned int iYearStop, unsigned int iMonthStop )
{
    plotProfiles( iYearStart, iMonthStart, iYearStop, iMonthStop, true );
}

void VAtmosphereSoundings::plotAverages( unsigned int iYearStart, unsigned int iMonthStart, unsigned int iYearStop, unsigned int iMonthStop, string iPlotOption, bool iSame )
{
    plotProfiles( iYearStart, iMonthStart, iYearStop, iMonthStop, false, iPlotOption, iSame );
}

void VAtmosphereSoundings::plotProfiles( unsigned int iYearStart, unsigned int iMonthStart, unsigned int iYearStop, unsigned int iMonthStop, bool b2D, string iPlotOption, bool bSames )
{

    cout << "VAtmosphereSoundings::plotProfiles() Warning: This function is obsolete" << endl;
    
    vector< string > iXaxis;
    vector< string > iYaxis;
    vector< int >    iXbin;
    vector< int >    iYbin;
    vector< double > iXmin;
    vector< double > iXmax;
    vector< double > iYmin;
    vector< double > iYmax;
    
    // temperature vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "temperature [K]" );
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 0. );
    iYmax.push_back( 500. );
    if( b2D )
    {
        iYmin.back() = 150.;
        iYmax.back() = 360.;
    }
    
    // dewpoint vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "dewpoint [K]" );
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 0. );
    iYmax.push_back( 500. );
    if( b2D )
    {
        iYmin.back() = 150.;
        iYmax.back() = 360.;
    }
    
    // pressure vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "pressure [Pa]" );
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 1. );
    iYmax.push_back( 1.e6 );
    if( b2D )
    {
        iYmin.back() = 0.;
        iYmax.back() = 6.;
    }
    
    // density vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "density * exp( height/7.739 km ) [g/cm2]" );
    if( fPlotRelativePlots )
    {
        iYaxis.back() = "density (relative plot)";
    }
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 0. );
    iYmax.push_back( 1.e6 );
    if( b2D )
    {
        iYbin.back() = 200;
        iXbin.back() = 100;
        iYmin.back() = 0.0008;
        iYmax.back() = 0.0017;
    }
    
    // relative humidity vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "relative humidity [%]" );
    if( fPlotRelativePlots )
    {
        iYaxis.back() = "rel humidity(relative plot)";
    }
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 0. );
    iYmax.push_back( 1.e2 );
    if( b2D )
    {
        iYbin.back() = 200;
        iXbin.back() = 100;
        iYmin.back() = 0.;
        iYmax.back() = 1.e2;
    }
    
    // thickness vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "atmospheric thickness [g/cm2]" );
    if( fPlotRelativePlots )
    {
        iYaxis.back() = "atmospheric thickness(relative plot)";
    }
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( 0. );
    iYmax.push_back( 5.e4 );
    if( b2D )
    {
        iYbin.back() = 200;
        iXbin.back() = 100;
        iYmin.back() = 0.;
        iYmax.back() = 1.2e3;
    }
    
    // log density vs height
    iXaxis.push_back( "height [km]" );
    iYaxis.push_back( "log(density [g/cm3])" );
    if( fPlotRelativePlots )
    {
        iYaxis.back() = "density (relative plot)";
    }
    iXbin.push_back( 50 );
    iYbin.push_back( 100 );
    iXmin.push_back( fPlottingHeight_min );
    iXmax.push_back( fPlottingHeight_max );
    iYmin.push_back( -13. );
    iYmax.push_back( -5. );
    if( b2D )
    {
        iYbin.back() = 200;
        iXbin.back() = 100;
    }
    
    
    /////////////////////////////
    // create and fill histograms
    char hname[800];
    char htitle[800];
    
    vector< map< unsigned int, TProfile* > > hMonthly;
    map< unsigned int, TProfile* > iTempMap;
    vector< TH2D* > hProfile2D;
    vector< TProfile* > hProfileAll;
    
    for( unsigned int k = 0; k < iXaxis.size(); k++ )
    {
        hMonthly.push_back( iTempMap );
        
        if( b2D )
        {
            sprintf( hname, "h%d", k );
            hProfile2D.push_back( new TH2D( hname, "", iXbin[k], iXmin[k], iXmax[k], iYbin[k], iYmin[k], iYmax[k] ) );
            hProfile2D[k]->SetStats( 0 );
            hProfile2D[k]->SetXTitle( iXaxis[k].c_str() );
            hProfile2D[k]->SetYTitle( iYaxis[k].c_str() );
            
        }
        else
        {
            hProfile2D.push_back( 0 );
        }
        sprintf( hname, "hP%d", k );
        hProfileAll.push_back( new TProfile( hname, "", iXbin[k], iXmin[k], iXmax[k], iYmin[k], iYmax[k] ) );
        hProfileAll[k]->SetStats( 0 );
        hProfileAll[k]->SetXTitle( iXaxis[k].c_str() );
        hProfileAll[k]->SetYTitle( iYaxis[k].c_str() );
        hProfileAll[k]->SetMarkerStyle( 7 );
    }
    
    // loop over all data sets and plot them
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        unsigned int iID = getHistogramIdentifier( i );
        
        // ignore all 0 identifiers
        if( iID == 0 )
        {
            continue;
        }
        
        for( unsigned int k = 0; k < iXaxis.size(); k++ )
        {
            if( hMonthly[k].find( iID ) == hMonthly[k].end() )
            {
                sprintf( hname, "h%d%d", k, iID );
                if( !b2D )
                {
                    hMonthly[k][iID] = new TProfile( hname, "", iXbin[k], iXmin[k], iXmax[k], iYmin[k], iYmax[k] );
                }
                else
                {
                    hMonthly[k][iID] = 0;
                }
                if( hMonthly[k][iID] )
                {
                    if( fBoolColorChange )
                    {
                        int iColor = hMonthly[k].size();
                        if( iColor > 9 )
                        {
                            iColor++;
                        }
                        if( iColor > 50 )
                        {
                            iColor = iColor % 50;
                        }
                        hMonthly[k][iID]->SetMarkerColor( iColor );
                        hMonthly[k][iID]->SetLineColor( iColor );
                        hMonthly[k][iID]->SetMarkerStyle( 20 );
                    }
                    else
                    {
                        hMonthly[k][iID]->SetMarkerStyle( 20 );
                    }
                    hMonthly[k][iID]->SetStats( 0 );
                    hMonthly[k][iID]->SetXTitle( iXaxis[k].c_str() );
                    hMonthly[k][iID]->SetYTitle( iYaxis[k].c_str() );
                }
            }
            // check time range
            if( fData[i]->Year * 100 + fData[i]->Month >= ( int )( iYearStart * 100 + iMonthStart )
                    && fData[i]->Year * 100 + fData[i]->Month <= ( int )( iYearStop * 100 + iMonthStop ) )
            {
                if( k == 0 )
                {
                    for( unsigned int j = 0; j < fData[i]->fTemperature_K.size(); j++ )
                    {
                        if( fData[i]->fTemperature_K[j] > 0. && fData[i]->fHeight_m[j] > 0.  && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 ) // Height approximately even mult. of 1000 feet. check this...
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fTemperature_K[j] );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fTemperature_K[j] );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fTemperature_K[j] );
                            }
                        }
                    }
                }
                else if( k == 1 )
                {
                    for( unsigned int j = 0; j < fData[i]->fDewPoint_K.size(); j++ )
                    {
                        if( fData[i]->fDewPoint_K[j] > 0. && fData[i]->fHeight_m[j] > 0. && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDewPoint_K[j] );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDewPoint_K[j] );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDewPoint_K[j] );
                            }
                        }
                    }
                }
                else if( k == 2 )
                {
                    for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
                    {
                        if( fData[i]->fPressure_Pa[j] > 0. && fData[i]->fHeight_m[j] > 0 && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fPressure_Pa[j] );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, log10( fData[i]->fPressure_Pa[j] ) );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, log10( fData[i]->fPressure_Pa[j] ) );
                            }
                        }
                    }
                }
                else if( k == 3 )
                {
                    for( unsigned int j = 0; j < fData[i]->fDensity_gcm3.size(); j++ )
                    {
                        if( fData[i]->fDensity_gcm3[j] > 0. && fData[i]->fHeight_m[j] > 0. && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                        {
                            if( !fPlotRelativePlots )
                            {
                                if( hMonthly[k][iID] )
                                {
                                    hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] * TMath::Exp( fData[i]->fHeight_m[j] / 1.e3 / 7.739 ) );
                                }
                                if( hProfile2D[k] )
                                {
                                    hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] * TMath::Exp( fData[i]->fHeight_m[j] / 1.e3 / 7.739 ) );
                                }
                                if( hProfileAll[k] )
                                {
                                    hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] * TMath::Exp( fData[i]->fHeight_m[j] / 1.e3 / 7.739 ) );
                                }
                            }
                            else
                            {
                                if( hMonthly[k][iID] )
                                {
                                    hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] );
                                }
                                if( hProfile2D[k] )
                                {
                                    hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] );
                                }
                                if( hProfileAll[k] )
                                {
                                    hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fDensity_gcm3[j] );
                                }
                            }
                        }
                    }
                }
                else if( k == 4 )
                {
                    for( unsigned int j = 0; j < fData[i]->fRelativeHumidity.size(); j++ )
                    {
                        if( fData[i]->fRelativeHumidity[j] > 0. && fData[i]->fHeight_m[j] > 0 && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fRelativeHumidity[j] );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fRelativeHumidity[j] );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fRelativeHumidity[j] );
                            }
                        }
                    }
                }
                else if( k == 5 )
                {
                    for( unsigned int j = 0; j < fData[i]->fThickness_gcm2.size(); j++ && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                    {
                        if( fData[i]->fThickness_gcm2[j] > 0. && fData[i]->fHeight_m[j] > 0 )
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fThickness_gcm2[j] );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fThickness_gcm2[j] );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, fData[i]->fThickness_gcm2[j] );
                            }
                        }
                    }
                }
                else if( k == 6 )
                {
                    for( unsigned int j = 0; j < fData[i]->fDensity_gcm3.size(); j++ && ( int )( fData[i]->fHeight_m[j] * 3.2808399 + 2 ) % 1000 < 4 )
                    {
                        if( fData[i]->fDensity_gcm3[j] > 0. && fData[i]->fHeight_m[j] > 0 )
                        {
                            if( hMonthly[k][iID] )
                            {
                                hMonthly[k][iID]->Fill( fData[i]->fHeight_m[j] / 1.e3, TMath::Log( fData[i]->fDensity_gcm3[j] ) );
                            }
                            if( hProfile2D[k] )
                            {
                                hProfile2D[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, TMath::Log( fData[i]->fDensity_gcm3[j] ) );
                            }
                            if( hProfileAll[k] )
                            {
                                hProfileAll[k]->Fill( fData[i]->fHeight_m[j] / 1.e3, TMath::Log( fData[i]->fDensity_gcm3[j] ) );
                            }
                        }
                    }
                }
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////
    // plot everything
    
    TCanvas* cTemp = 0;
    
    if( !bSames )
    {
        fCanvasProfile.clear();
        fCanvas2D.clear();
    }
    
    for( unsigned int k = 0; k < hMonthly.size(); k++ )
    {
        if( !bSames )
        {
            sprintf( hname, "c%s%s%d", iYaxis[k].c_str(), iXaxis[k].c_str(), b2D );
            sprintf( htitle, "%s vs %s", iYaxis[k].c_str(), iXaxis[k].c_str() );
            if( !b2D )
            {
                fCanvasProfile.push_back( new TCanvas( hname, htitle, 10, 10, 600, 400 ) );
                cTemp = fCanvasProfile.back();
            }
            else
            {
                fCanvas2D.push_back( new TCanvas( hname, htitle, 10, 10, 600, 400 ) );
                cTemp = fCanvas2D.back();
            }
            cTemp->SetLeftMargin( 0.13 );
            cTemp->SetGridx( 0 );
            cTemp->SetGridy( 0 );
            cTemp->Draw();
            
            fPlottingLegend.push_back( new TLegend( 0.5, 0.15, 0.59, 0.35 ) );
        }
        else
        {
            if( !b2D && k < fCanvasProfile.size() )
            {
                cTemp = fCanvasProfile[k];
            }
            else if( k < fCanvas2D.size() )
            {
                cTemp = fCanvas2D[k];
            }
            else
            {
                continue;
            }
        }
        
        
        if( b2D )
        {
            hProfile2D[k]->Draw( "colz" );
            hProfile2D[k]->GetYaxis()->SetTitleOffset( 1.5 );
            if( hProfileAll[k] )
            {
                hProfileAll[k]->Draw( "same" );
            }
            cTemp->Update();
        }
        else
        {
            map< unsigned int, TProfile* >::iterator iterTempMap;
            
            int z = 0;
            for( iterTempMap = hMonthly[k].begin(); iterTempMap != hMonthly[k].end(); ++iterTempMap )
            {
                TProfile* h = ( TProfile* )( *iterTempMap ).second;
                if( h )
                {
                    if( h->GetEntries() > 0 )
                    {
                        if( fPlotRelativePlots )
                        {
                            h->Divide( hProfileAll[k] );
                        }
                        if( z == 0 && !bSames )
                        {
                            sprintf( hname, "%s", iPlotOption.c_str() );
                            h->GetYaxis()->SetTitleOffset( 1.5 );
                            if( k == 3 )
                            {
                                h->SetAxisRange( 0.0004, 0.0017, "Y" );    // density * XXX
                            }
                            else if( k == 0 )
                            {
                                h->SetAxisRange( 180., 340., "Y" );    // temperature
                            }
                            else if( k == 1 )
                            {
                                h->SetAxisRange( 180., 340., "Y" );    // dew point
                            }
                            else if( k == 2 )
                            {
                                h->SetAxisRange( 8.2, 1.5e5, "Y" );    // pressure
                            }
                            else if( k == 4 )
                            {
                                h->SetAxisRange( 1.e-2, 100., "Y" );    // relative humidity
                            }
                            else if( k == 5 )
                            {
                                h->SetAxisRange( 1.e-2, 1.2e3, "Y" );    // atmospheric thickness
                            }
                            else if( k == 6 )
                            {
                                h->SetAxisRange( -13., -5., "Y" );    // atmospheric thickness
                            }
                            h->Draw( hname );
                        }
                        else
                        {
                            sprintf( hname, "same %s", iPlotOption.c_str() );
                            h->Draw( hname );
                        }
                        z++;
                        cTemp->Update();
                        sprintf( hname, "%d", ( *iterTempMap ).first );
                        if( fPlottingLegend[k] )
                        {
                            fPlottingLegend[k]->AddEntry( h, hname, "pl" );
                        }
                    }
                    
                }
            }
        }
        // add CORSIKA and user defined atmosphere plots
        if( k == 0 )
        {
            plotCORSIKA_Temperature_vs_Height( cTemp );
            plotUserAtmosphere_Temperature_vs_Height( cTemp );
        }
        else if( k == 1 )
        {
            plotCORSIKA_DewPoint_vs_Height( cTemp );
            plotUserAtmosphere_DewPoint_vs_Height( cTemp );
        }
        else if( k == 2 )
        {
            plotCORSIKA_Pressure_vs_Height( cTemp );
            plotUserAtmosphere_Pressure_vs_Height( cTemp );
        }
        else if( k == 3 )
        {
            plotCORSIKA_Density_vs_Height( cTemp );
            plotUserAtmosphere_Density_vs_Height( cTemp );
        }
        else if( k == 4 )
        {
            plotCORSIKA_RelativeHumidity_vs_Height( cTemp );
            plotUserAtmosphere_RelativeHumidity_vs_Height( cTemp );
        }
        else if( k == 5 )
        {
            plotCORSIKA_Thickness_vs_Height( cTemp );
            plotUserAtmosphere_Thickness_vs_Height( cTemp );
        }
        
        if( fPlottingLegendDraw && fPlottingLegend[k] )
        {
            fPlottingLegend[k]->Draw();
        }
        if( k == 3 )
        {
            cTemp->SaveAs( TString::Format( "density-%d-%d%d-%d%d.png", b2D, iYearStart, iMonthStart, iYearStop, iMonthStop ).Data() );
        }
        if( k == 6 )
        {
            cTemp->SaveAs( TString::Format( "logdensity-%d-%d%d-%d%d.png", b2D, iYearStart, iMonthStart, iYearStop, iMonthStop ).Data() );
        }
    }
}

/*
 * open root file with sounding data and fill it into
 * data class structure
 *
*/
bool VAtmosphereSoundings::readRootFile( unsigned npoints_min )
{
    if( !fDataTree )
    {
        return false;
    }
    
    fDataTree->SetBranchAddress( "MJD", &MJD );
    fDataTree->SetBranchAddress( "Year", &Year );
    fDataTree->SetBranchAddress( "Month", &Month );
    fDataTree->SetBranchAddress( "Day", &Day );
    fDataTree->SetBranchAddress( "Hour", &Hour );
    fDataTree->SetBranchAddress( "nPoints", &nPoints );
    fDataTree->SetBranchAddress( "Height_m", Height_m );
    fDataTree->SetBranchAddress( "Pressure_Pa", Pressure_Pa );
    fDataTree->SetBranchAddress( "Density_gcm3", Density_gcm3 );
    fDataTree->SetBranchAddress( "Thickness_gcm2", Thickness_gcm2 );
    fDataTree->SetBranchAddress( "Temperature_K", Temperature_K );
    fDataTree->SetBranchAddress( "DewPoint_K", DewPoint_K );
    fDataTree->SetBranchAddress( "RelativeHumidity", RelativeHumidity );
    fDataTree->SetBranchAddress( "VaporMassDensity_gm3", VaporMassDensity_gm3 );
    fDataTree->SetBranchAddress( "MixingRatio_gkg", MixingRatio_gkg );
    fDataTree->SetBranchAddress( "WindDirection_deg", WindDirection_deg );
    fDataTree->SetBranchAddress( "WindSpeed_ms", WindSpeed_ms );
    fDataTree->SetBranchAddress( "IndexofRefraction", IndexofRefraction );
    
    unsigned int z = 0;
    for( unsigned int i = 0; i < fDataTree->GetEntries(); i++ )
    {
        fDataTree->GetEntry( i );
        
        // require a minimum amount of data points in
        // sounding profile
        if( npoints_min > 0 && npoints_min > nPoints )
        {
            continue;
        }
        
        fData.push_back( new VAtmosphereSoundingData() );
        
        fData[z]->MJD = MJD;
        fData[z]->Year = Year;
        fData[z]->Month = Month;
        fData[z]->Day = Day;
        fData[z]->Hour = Hour;
        
        for( unsigned int j = 0; j < nPoints; j++ )
        {
            fData[z]->fPressure_Pa.push_back( Pressure_Pa[j] );
            fData[z]->fHeight_m.push_back( Height_m[j] );
            fData[z]->fDensity_gcm3.push_back( Density_gcm3[j] );
            fData[z]->fThickness_gcm2.push_back( Thickness_gcm2[j] );
            fData[z]->fTemperature_K.push_back( Temperature_K[j] );
            fData[z]->fDewPoint_K.push_back( DewPoint_K[j] );
            fData[z]->fRelativeHumidity.push_back( RelativeHumidity[j] );
            if( fData[z]->fRelativeHumidity.back() > 0. )
            {
                //fData[z]->fRelativeHumidity.back() *= 1.e2; //HF is this needed?
            }
            fData[z]->fVaporMassDensity_gm3.push_back( VaporMassDensity_gm3[j] );
            fData[z]->fIndexofRefraction.push_back( IndexofRefraction[j] ); //TODO: change this
            fData[z]->fMixingRatio_gkg.push_back( MixingRatio_gkg[j] );
            fData[z]->fWindDirection_deg.push_back( WindDirection_deg[j] );
            fData[z]->fWindSpeed_ms.push_back( WindSpeed_ms[j] );
        }
        z++;
        
    }
    
    cout << "total number of data sets available: " << fData.size() << endl;
    
    make_interpolated_atmospheres();
    
    return true;
}

TCanvas* VAtmosphereSoundings::plotCORSIKA( TCanvas* c, int iPlotID, vector< VAtmosphereSoundingData* > iData,
        double iHeightMin, double iHeightMax )
{
    char hname[800];
    
    string iCTitle;
    string iHTitle;
    string iXTitle;
    double iYmin = 0.;
    double iYmax = 0.;
    
    // density vs height
    if( iPlotID == 0 )
    {
        iHTitle = "CORSIKADH";
        iCTitle = "CORSIKA (density vs height)";
        iXTitle = "density * exp( height/7.739 km ) [kg/cm2]";
        iYmin = 0.0004;
        iYmax = 0.0017;
        if( iHeightMax < 40. )
        {
            iYmin = 0.0004;
        }
    }
    // index of refraction vs height
    else if( iPlotID == 1 )
    {
        iHTitle = "CORSIKANH";
        iCTitle = "CORSIKA (index of refraction vs height)";
        iXTitle = "index of refraction [n-1]";
        iYmin = 1.e-10;
        iYmax = 0.00033;
        if( iHeightMax < 7. )
        {
            iYmin = 0.15e-3;
        }
    }
    // temperature vs height
    else if( iPlotID == 2 )
    {
        iHTitle = "CORSIKATH";
        iCTitle = "CORSIKA (temperature vs height)";
        iXTitle = "temperature [K]";
        iYmin = 150;
        iYmax = 350;
    }
    // ozone vs height
    else if( iPlotID == 3 )
    {
        iHTitle = "CORSIKAO3";
        iCTitle = "CORSIKA (ozone vs height)";
        iXTitle = "ozone [ATM cm/km]";
        iYmin = 1.e-8;
        iYmax = 1.e-1;
    }
    // relative humidity vs height
    else if( iPlotID == 4 )
    {
        iHTitle = "CORSIKAO4";
        iCTitle = "CORSIKA (relative humidity vs height)";
        iXTitle = "relative humidity [%]";
        iYmin = 0.;
        iYmax = 100.;
    }
    // dew point vs height
    else if( iPlotID == 5 )
    {
        iHTitle = "CORSIKAO5";
        iCTitle = "CORSIKA (dew point vs height)";
        iXTitle = "dew point [K]";
        iYmin = 100.;
        iYmax = 400.;
    }
    // pressure vs height
    else if( iPlotID == 6 )
    {
        iHTitle = "CORSIKAO6";
        iCTitle = "CORSIKA (pressure vs height)";
        iXTitle = "pressure [Pa]";
        iYmin = 0.;
        iYmax = 1.e6;
    }
    // thickness vs height
    else if( iPlotID == 7 )
    {
        iHTitle = "CORSIKAO7";
        iCTitle = "CORSIKA (thickness vs height)";
        iXTitle = "atmospheric thickness [g/cm2]";
        iYmin = 0.;
        iYmax = 1.e6;
    }
    
    
    //////////////////////////////////////////////////////
    bool bNewCanvas = true;
    if( c )
    {
        c->cd();
        bNewCanvas = false;
    }
    else
    {
        sprintf( hname, "c%s", iHTitle.c_str() );
        c = new TCanvas( hname, iCTitle.c_str(), 10, 10, 600, 400 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLeftMargin( 0.15 );
        sprintf( hname, "h%s", iHTitle.c_str() );
        TH1D* h = new TH1D( hname, "", 100, iHeightMin, iHeightMax );
        h->SetStats( 0 );
        h->SetMinimum( iYmin );
        h->SetMaximum( iYmax );
        h->SetXTitle( "height[km]" );
        h->SetYTitle( iXTitle.c_str() );
        h->GetYaxis()->SetTitleOffset( 1.5 );
        h->Draw();
    }
    vector< TGraph* > g;
    
    TLegend* iL = 0;
    if( iData.size() == 2 )
    {
        iL = new TLegend( 0.55, 0.65, 0.85, 0.75 );
    }
    else
    {
        iL = new TLegend( 0.65, 0.55, 0.85, 0.75 );
    }
    
    for( unsigned int i = 0; i < iData.size(); i++ )
    {
        int z = 0;
        g.push_back( new TGraph( 1 ) );
        g.back()->SetLineColor( iData[i]->PlotColor );
        g.back()->SetMarkerColor( iData[i]->PlotColor );
        g.back()->SetLineStyle( iData[i]->PlotLineStyle );
        g.back()->SetLineWidth( iData[i]->PlotLineWidth );
        g.back()->SetTitle( "" );
        
        for( unsigned int j = 0; j < iData[i]->fDensity_gcm3.size(); j++ )
        {
            if( iPlotID == 0 )
            {
                if( iData[i]->fDensity_gcm3[j] > -90000. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fDensity_gcm3[j] * TMath::Exp( iData[i]->fHeight_m[j] / 7.739 / 1.e3 ) );
                    z++;
                }
            }
            else if( iPlotID == 1 )
            {
                if( iData[i]->fIndexofRefraction[j] > -90000. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fIndexofRefraction[j] - 1. );
                    z++;
                }
            }
            else if( iPlotID == 2 )
            {
                if( iData[i]->fTemperature_K[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fTemperature_K[j] );
                    z++;
                }
            }
            else if( iPlotID == 3 )
            {
                if( iData[i]->fO3_cmkm[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fO3_cmkm[j] );
                    z++;
                }
            }
            else if( iPlotID == 4 )
            {
                if( iData[i]->fRelativeHumidity[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fRelativeHumidity[j] );
                    z++;
                }
            }
            else if( iPlotID == 5 )
            {
                if( iData[i]->fDewPoint_K[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fDewPoint_K[j] );
                    z++;
                }
            }
            else if( iPlotID == 6 )
            {
                if( iData[i]->fPressure_Pa[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fPressure_Pa[j] );
                    z++;
                }
            }
            else if( iPlotID == 7 )
            {
                if( iData[i]->fThickness_gcm2[j] > 0. && iData[i]->fHeight_m[j] > -90000. )
                {
                    g.back()->SetPoint( z, iData[i]->fHeight_m[j] / 1.e3, iData[i]->fThickness_gcm2[j] );
                    z++;
                }
            }
        }
        if( g.back()->GetN() > 1 )
        {
            if( bNewCanvas )
            {
                bNewCanvas = false;
            }
            g.back()->Draw( "l" );
            iL->AddEntry( g.back(), iData[i]->Name.c_str(), "pl" );
        }
    }
    if( iData.size() > 0 )
    {
        iL->SetFillColor( 0 );
        iL->Draw();
    }
    
    return c;
}

/*
   add a new profile from soundings and MODTRAN/CORSIKA data

   (some values are hardcoded here!)

   iIndexCORSIKAMODTRAN   : index of MODTRAN data set to use
   iHeightMaxData         : height [km] to change between sounding data and model
   iName                  : profile name (currently not used)
*/
bool VAtmosphereSoundings::add_user_Atmosphere( unsigned int iIndexCORSIKAMODTRAN, double iHeightMaxData, string iName )
{
    if( iIndexCORSIKAMODTRAN >= fDataCORSIKAMODTRAN.size() )
    {
        cout << "VAtmosphereSoundings::add_user_profile: data set index out of range (MODTRAN): ";
        cout << iIndexCORSIKAMODTRAN << "\t" << fDataCORSIKAMODTRAN.size() << endl;
        return false;
    }
    if( !fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN] )
    {
        cout << "VAtmosphereSoundings::add_user_profile: empty data set (MODTRAN): " << iIndexCORSIKAMODTRAN << endl;
        return false;
    }
    //////////////////////////////////////////////////////////////////////
    // create a new user profile
    fDataUserProfile.push_back( new VAtmosphereSoundingData() );
    fDataUserProfile.back()->Name  = iName;
    
    // use height steps from Monte Carlo
    fDataUserProfile.back()->setdefaultvalues( fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN]->fHeight_m.size() );
    for( unsigned int i = 0; i < fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN]->fHeight_m.size(); i++ )
    {
        fDataUserProfile.back()->fHeight_m[i] = fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN]->fHeight_m[i];
    }
    
    // pressure
    fDataUserProfile.back()->fPressure_Pa = getDataVectorForUserAtmosphere( iHeightMaxData, fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN], "pressure" );
    // temperature
    fDataUserProfile.back()->fTemperature_K = getDataVectorForUserAtmosphere( iHeightMaxData, fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN], "temperature" );
    // dew point
    fDataUserProfile.back()->fDewPoint_K = getDataVectorForUserAtmosphere( iHeightMaxData, fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN], "dewpoint" );
    // relative humidity
    fDataUserProfile.back()->fRelativeHumidity = getDataVectorForUserAtmosphere( iHeightMaxData, fDataCORSIKAMODTRAN[iIndexCORSIKAMODTRAN], "relhumid" );
    // dew point (from relative humidity and temperature)
    /*   for( unsigned int i = 0; i < fDataUserProfile.back()->fRelativeHumidity.size(); i++ )
       {
          if( i < fDataUserProfile.back()->fDewPoint_K.size() )
          {
             fDataUserProfile.back()->fDewPoint_K[i] = getDewPoint( fDataUserProfile.back()->fTemperature_K[i], fDataUserProfile.back()->fRelativeHumidity[i], 1 );
          }
          else
          {
             fDataUserProfile.back()->fDewPoint_K.push_back( getDewPoint( fDataUserProfile.back()->fTemperature_K[i], fDataUserProfile.back()->fRelativeHumidity[i], 1 ) );
          }
       } */
    // atmospheric thickness
    fillAtmosphericThickness( fDataUserProfile.back() );
    // atmospheric density
    fillAtmosphericDensity( fDataUserProfile.back() );
    //refractive index
    fillIndexofRefraction( );
    return true;
}

vector< double > VAtmosphereSoundings::getDataVectorForUserAtmosphere( double iHeightMaxData, VAtmosphereSoundingData* iDataMonteCarlo, string iType )
{
    vector< double > f;   // return vector
    vector< double > n;   // counting vector
    if( !iDataMonteCarlo )
    {
        return f;
    }
    
    f.assign( iDataMonteCarlo->fHeight_m.size(), 0. );
    n.assign( iDataMonteCarlo->fHeight_m.size(), 0. );
    
    // km -> m
    iHeightMaxData *= 1.e3;
    
    // type IDs
    unsigned int iTypeID = 0;
    if( iType == "pressure" )
    {
        iTypeID = 1;
    }
    else if( iType == "temperature" )
    {
        iTypeID = 2;
    }
    else if( iType == "dewpoint" )
    {
        iTypeID = 3;
    }
    else if( iType == "relhumid" )
    {
        iTypeID = 4;
    }
    
    // unknown type
    if( iTypeID == 0 )
    {
        return f;
    }
    
    // calculate height bins
    vector< double > iHeight_min( iDataMonteCarlo->fHeight_m.size(), 0. );
    vector< double > iHeight_max( iDataMonteCarlo->fHeight_m.size(), 0. );
    
    for( unsigned int i = 0; i < iDataMonteCarlo->fHeight_m.size(); i++ )
    {
        if( i == 0 )
        {
            iHeight_min[i] = iDataMonteCarlo->fHeight_m[0];
        }
        else
        {
            iHeight_min[i] = 0.5 * ( iDataMonteCarlo->fHeight_m[i] + iDataMonteCarlo->fHeight_m[i - 1] );
        }
        if( i == iDataMonteCarlo->fHeight_m.size() - 1 )
        {
            iHeight_max[i] = iDataMonteCarlo->fHeight_m[i];
        }
        else if( i == 0 )
        {
            iHeight_max[i]  = 0.5 * ( iDataMonteCarlo->fHeight_m[i] + iDataMonteCarlo->fHeight_m[i + 1] );
        }
        else
        {
            iHeight_max[i] = iDataMonteCarlo->fHeight_m[i] + ( iDataMonteCarlo->fHeight_m[i] - iHeight_min[i] );
        }
    }
    
    // loop over all data sets
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        unsigned int iID = getHistogramIdentifier( i );
        
        // ignore all 0 identifiers
        if( iID == 0 )
        {
            continue;
        }
        
        ////////////////////////////////////////////////////////////////////////
        // pressure
        if( iTypeID == 1 )
        {
            // loop over this data set
            for( unsigned int j = 0; j < fData[i]->fPressure_Pa.size(); j++ )
            {
                if( fData[i]->fPressure_Pa[j] > 0. )
                {
                    // fill data value for heights < iHeightMaxData
                    //		if( fData[i]->fHeight_m[j] < iHeightMaxData )
                    {
                        for( unsigned int k = 1; k < iDataMonteCarlo->fHeight_m.size(); k++ )
                        {
                            if( fData[i]->fHeight_m[j] > iHeight_min[k] && fData[i]->fHeight_m[j] < iHeight_max[k] )
                            {
                                f[k] += fData[i]->fPressure_Pa[j];
                                n[k] ++;
                            }
                        }
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // temperature
        else if( iTypeID == 2 )
        {
            // loop over this data set
            for( unsigned int j = 0; j < fData[i]->fTemperature_K.size(); j++ )
            {
                if( fData[i]->fTemperature_K[j] > 0. )
                {
                    // fill data value for heights < iHeightMaxData
                    //		if( fData[i]->fHeight_m[j] < iHeightMaxData )
                    {
                        for( unsigned int k = 1; k < iDataMonteCarlo->fHeight_m.size(); k++ )
                        {
                            if( fData[i]->fHeight_m[j] > iHeight_min[k] && fData[i]->fHeight_m[j] < iHeight_max[k] )
                            {
                                f[k] += fData[i]->fTemperature_K[j];
                                n[k] ++;
                            }
                        }
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // dew point
        else if( iTypeID == 3 )
        {
            // loop over this data set
            for( unsigned int j = 0; j < fData[i]->fDewPoint_K.size(); j++ )
            {
                if( fData[i]->fDewPoint_K[j] > 0. )
                {
                    // fill data value for heights < iHeightMaxData
                    //		if( fData[i]->fHeight_m[j] < iHeightMaxData )
                    {
                        for( unsigned int k = 1; k < iDataMonteCarlo->fHeight_m.size(); k++ )
                        {
                            if( fData[i]->fHeight_m[j] > iHeight_min[k] && fData[i]->fHeight_m[j] < iHeight_max[k] )
                            {
                                f[k] += fData[i]->fDewPoint_K[j];
                                n[k] ++;
                            }
                        }
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // relative humidity
        else if( iTypeID == 4 )
        {
            // loop over this data set
            for( unsigned int j = 0; j < fData[i]->fRelativeHumidity.size(); j++ )
            {
                if( fData[i]->fRelativeHumidity[j] >= 0. && fData[i]->fRelativeHumidity[j] <= 100. )
                {
                    // fill data value for heights < iHeightMaxData
                    //		if( fData[i]->fHeight_m[j] < iHeightMaxData )
                    {
                        for( unsigned int k = 1; k < iDataMonteCarlo->fHeight_m.size(); k++ )
                        {
                            if( fData[i]->fHeight_m[j] > iHeight_min[k] && fData[i]->fHeight_m[j] < iHeight_max[k] )
                            {
                                f[k] += fData[i]->fRelativeHumidity[j];
                                n[k] ++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////
    // calculate mean value
    for( unsigned int i = 0; i < f.size(); i++ )
    {
        if( i < n.size() && n[i] > 0. )
        {
            f[i] /= n[i];
        }
        else
        {
            f[i]  = -9999;
        }
        
        // fill MC values for > iHeightMaxData
        if( iDataMonteCarlo->fHeight_m[i] >= iHeightMaxData )
        {
            if( iTypeID == 1 )
            {
                f[i] = iDataMonteCarlo->fPressure_Pa[i];
            }
            else if( iTypeID == 2 )
            {
                f[i] = iDataMonteCarlo->fTemperature_K[i];
            }
            else if( iTypeID == 3 )
            {
                f[i] = iDataMonteCarlo->fDewPoint_K[i];
            }
            else if( iTypeID == 4 )
            {
                f[i] = iDataMonteCarlo->fRelativeHumidity[i];
            }
        }
    }
    
    return f;
}

bool VAtmosphereSoundings::write_MODTRAN_UserProfile( unsigned int iIndexUserData, unsigned int defaultAtmosphericModel, bool iWriteDewPoint )
{
    if( iIndexUserData >= fDataUserProfile.size() )
    {
        cout << "VAtmosphereSoundings::write_MODTRAN_UserProfile: index out of range: " << iIndexUserData << "\t" << fDataUserProfile.size() << endl;
        return false;
    }
    if( !fDataUserProfile[iIndexUserData] )
    {
        cout << "VAtmosphereSoundings::write_MODTRAN_UserProfile: no data" << endl;
        return false;
    }
    
    ////////////////////////////////////////////////////////////
    // CARD 2C:	ML, IRD1, IRD2, HMODEL
    // FORMAT (315, A20)
    printf( "%5d", ( int )fDataUserProfile[iIndexUserData]->fHeight_m.size() );
    printf( "%5d", 0 );
    printf( "%5d", 0 );
    //   printf( "%65s\n", iProfileName.c_str() );
    cout << endl;
    
    ////////////////////////////////////////////////////////////
    // CARD 2C1: ZM, P, T, (WMOL(J), J = 1, 3), (JCHAR(J), J = 1, 14), JCHARX
    // FORMAT (F10.3, 5E10.3, 14A1, 1X, A1)
    
    ////////////////////////////////////////////////////////////
    // pressure
    for( unsigned int i = 0; i < fDataUserProfile[iIndexUserData]->fHeight_m.size(); i++ )
    {
        printf( "%10.3f", fDataUserProfile[iIndexUserData]->fHeight_m[i] / 1.e3 );
        
        if( fDataUserProfile[iIndexUserData]->fPressure_Pa[i] > 0. )
        {
            printf( "%10.3f",  fDataUserProfile[iIndexUserData]->fPressure_Pa[i] / 1.e2 ); // Pa -> mbar
        }
        else
        {
            printf( "%10.3f", 1.e-3 );    // MODTRAN requires pressure > 0.
        }
        
        ////////////////////////////////////////////////////////////
        // temperature
        printf( "%10.3f", fDataUserProfile[iIndexUserData]->fTemperature_K[i] );
        
        ////////////////////////////////////////////////////////////
        // dew point or relative humidity
        if( iWriteDewPoint )
        {
            printf( "%10.3f", fDataUserProfile[iIndexUserData]->fDewPoint_K[i] );
        }
        else
        {
            printf( "%10.3f", fDataUserProfile[iIndexUserData]->fRelativeHumidity[i] );    // relative humdity (in percent)
        }
        
        ////////////////////////////////////////////////////////////
        // remaining values need by MODTRAN
        printf( "%10.3f", 0. );
        printf( "%10.3f", 0. );
        if( iWriteDewPoint )
        {
            printf( "AAF" );    // pressure in [mbar]; temperature in [K]; dew point [K];
        }
        else
        {
            printf( "AAH" );    // pressure in [mbar]; temperature in [K]; relative humidity;
        }
        for( unsigned int j = 0; j < 13; j++ )
        {
            cout << defaultAtmosphericModel;
        }
        cout << endl;
    }
    
    return true;
}
bool VAtmosphereSoundings::write_2C1( unsigned int iIndexAverageData, string filename, double max_height = 150e3 )
{
    if( iIndexAverageData >= fAverageProfile.size() )
    {
        cout << "VAtmosphereSoundings::write_2C1: index out of range: " << iIndexAverageData << "\t" << fAverageProfile.size() << endl;
        return false;
    }
    if( !fAverageProfile[iIndexAverageData] )
    {
        cout << "VAtmosphereSoundings::write_2C1: no data" << endl;
        return false;
    }
    
    return fAverageProfile[iIndexAverageData]->write_2C1( filename, &fModtranHeights, max_height ) ;
    
}

/*
   this function is currently not used (but might be useful later)
*/
double VAtmosphereSoundings::getInterpolation( double h, VAtmosphereSoundingData* iData, string iType )
{
    double f = -9999;      // return value
    
    if( h < 0. || !iData )
    {
        return f;
    }
    
    int iLup = -1;
    int iLlo = -1;
    for( unsigned int i = 0; i < iData->fHeight_m.size(); i++ )
    {
        if( h > iData->fHeight_m[i] )
        {
            iLup = ( int )i;
            if( i > 0 )
            {
                iLlo = ( int )( i - 1 );
            }
            else
            {
                iLlo = 0;
            }
            break;
        }
    }
    if( iLlo < 0 || iLup < 0 )
    {
        return f;
    }
    
    double yup = 0.;
    double ylo = 0.;
    
    if( iType == "pressure" )
    {
        if( iLup < ( int )iData->fPressure_Pa.size() )
        {
            yup = iData->fPressure_Pa[iLup];
        }
        if( iLlo < ( int )iData->fPressure_Pa.size() )
        {
            ylo = iData->fPressure_Pa[iLlo];
        }
    }
    else if( iType == "temperature" )
    {
        if( iLup < ( int )iData->fTemperature_K.size() )
        {
            yup = iData->fTemperature_K[iLup];
        }
        if( iLlo < ( int )iData->fTemperature_K.size() )
        {
            ylo = iData->fTemperature_K[iLlo];
        }
    }
    else if( iType == "dewpoint" )
    {
        if( iLup < ( int )iData->fDewPoint_K.size() )
        {
            yup = iData->fDewPoint_K[iLup];
        }
        if( iLlo < ( int )iData->fDewPoint_K.size() )
        {
            ylo = iData->fDewPoint_K[iLlo];
        }
    }
    else if( iType == "relhum" )
    {
        if( iLup < ( int )iData->fRelativeHumidity.size() )
        {
            yup = iData->fRelativeHumidity[iLup];
        }
        if( iLlo < ( int )iData->fRelativeHumidity.size() )
        {
            ylo = iData->fRelativeHumidity[iLlo];
        }
    }
    
    // interpolate
    if( iLup > iLlo && TMath::Abs( iData->fHeight_m[iLlo] - iData->fHeight_m[iLup] ) > 1.e-2 )
    {
        f = yup + ( h - iData->fHeight_m[iLup] ) * ( yup - ylo ) / ( iData->fHeight_m[iLlo] - iData->fHeight_m[iLup] );
    }
    else
    {
        f = yup;
    }
    
    return f;
}

double VAtmosphereSoundings::safe_eval( TGraph* g, double h, string opt )
{
    if( !g )
    {
        cout << "VAtmosphereSoundings::safe_eval: Error: no graph " << endl;
        return -9999;
    }
    
    if( g->GetN() < 1 )
    {
        //cout << "VAtmosphereSoundings::safe_eval: Error: Graph has 0 points " << endl;
        return -9999;
    }
    
    if( h < g->GetX()[0] || h > g->GetX()[g->GetN() - 1] )
    {
        return -9999; //h outside of limits
    }
    else   //
    {
        return ( opt == "log" ) ?  TMath::Exp( g->Eval( h ) ) : g->Eval( h ) ;
    }
    
}


double VAtmosphereSoundings::interpolate( vector<double> raw, vector<double> raw_heights, vector<double>& result, string opt = "lin", double h = -1 )
{

    TGraph* g = new TGraph( 0 );
    for( unsigned int i = 0; i < raw.size(); i++ )
    {
        if( raw_heights[i] != -9999 && raw[i] != -9999 )
        {
            g->SetPoint( g->GetN(), raw_heights[i], raw[i] );
        }
    }
    if( opt == "log" ) for( int i = 0; i < g->GetN() ; i++ )
        {
            g->SetPoint( i, g->GetX()[i], TMath::Log( g->GetY()[i] ) );
        }
    g->Sort();
    
    if( h > 0 ) //we are asked to evaluate at a given height.
    {
        return safe_eval( g, h, opt );
    }
    
    //no height given->evaluate at all entries of fHeights
    
    result.clear();
    
    for( unsigned int i = 0; i < fHeights.size(); i++ )
    {
        result.push_back( safe_eval( g, fHeights.at( i ),  opt ) );
    }
    
    g->Delete();
    return 42;
}


void VAtmosphereSoundings::list_datasets_CORSIKAMODTRAN()
{
    for( unsigned int i = 0; i < fDataCORSIKAMODTRAN.size(); i++ )
    {
        cout << i << "\t";
        if( fDataCORSIKAMODTRAN[i] )
        {
            cout << fDataCORSIKAMODTRAN[i]->Name << "\t";
            cout << fDataCORSIKAMODTRAN[i]->fPressure_Pa.size() << "\t";
        }
        cout << endl;
    }
}

void VAtmosphereSoundings::list_datasets()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        cout << i << "\t";
        if( fData[i] )
        {
            cout << fData[i]->Name << "\t";
            cout << fData[i]->fPressure_Pa.size() << "\t";
        }
        cout << endl;
    }
}

bool VAtmosphereSoundings::write_CORSIKA_UserProfile( unsigned int iMODTRANIndex, unsigned int atmprofmodel, string iName )
{
    if( iMODTRANIndex >= fDataCORSIKAMODTRAN.size() )
    {
        cout << "VAtmosphereSoundings::write_CORSIKA_UserProfile: index out of range: " << iMODTRANIndex << "\t" << fDataCORSIKAMODTRAN.size() << endl;
        return false;
    }
    if( !fDataCORSIKAMODTRAN[iMODTRANIndex] )
    {
        cout << "VAtmosphereSoundings::write_CORSIKA_UserProfile: no data" << endl;
        return false;
    }
    
    return fDataCORSIKAMODTRAN[iMODTRANIndex]->write_CORSIKA_UserProfile( atmprofmodel, iName );
    
}


/*
 * average atmosphere over a certain period
 *
*/
int VAtmosphereSoundings::push_average_atmosphere( string name = "",
        vector<int>* years = 0, vector<int>* months = 0,
        vector<int>* days = 0, vector<int>* hours = 0,
        vector<double>* mjds = 0,
        unsigned int nMinPoints = 20, int nMinFlights = 1 )
{
    // average profile
    VAtmosphereSoundingData* Average =  new VAtmosphereSoundingData();
    Average->Name  = name;
    
    vector<double> heightbins;
    
    heightbins.push_back( 1.5 * fHeights[0] - 0.5 * fHeights[1] );
    for( unsigned int i = 0; i < fHeights.size() - 1; i++ )
    {
        heightbins.push_back( 0.5 * fHeights[i] + 0.5 * fHeights[i + 1] );
    }
    heightbins.push_back( 1.5 * fHeights[fHeights.size() - 1] - 0.5 * fHeights[fHeights.size() - 2] );
    
    TProfile* averagePressure = new TProfile(
        "averagePressure", "average pressure; pressure [Pa]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 2e5, "s" );
    TProfile* averageDensity = new TProfile(
        "averageDensity", "average density; density [g/cm^{3}]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 2e5, "s" );
    TProfile* averageTemperature = new TProfile(
        "averageTemperature", "average temperature; temperature [K]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 5e2, "s" );
    TProfile* averageDewpoint = new TProfile(
        "averageDewpoint", "average dewpoint; dewpoint [K]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 5e2 );
    TProfile* averageMixingRatio = new TProfile(
        "averageMixingRatio", "average mixing ratio; mixing ratio [g/kg]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 1e3 );
    TProfile* averageRelativeHumidity = new TProfile(
        "averageRelativeHumidity", "average relative humidity; relative humidity [%]; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 1e2, "s" );
    TProfile* averageIndexofRefraction = new TProfile(
        "averageIndexofRefraction", "average index of refraction; index of refraction; height [m]",
        heightbins.size() - 1, &( heightbins[0] ), 0, 2, "s" );
        
    for( unsigned int iData = 0; iData < fDataInterpol.size(); iData++ )
    {
        VAtmosphereSoundingData* Data = fDataInterpol.at( iData );
        if( !isDateInRange( Data, years, months, days, hours, mjds , nMinPoints ) )
        {
            continue;
        }
        for( unsigned int iH = 0; iH < fHeights.size(); iH++ )
        {
            averagePressure->Fill( fHeights[iH], Data->fPressure_Pa.at( iH ) );
            averageDensity->Fill( fHeights[iH], Data->fDensity_gcm3.at( iH ) );
            averageTemperature->Fill( fHeights[iH], Data->fTemperature_K.at( iH ) );
            averageDewpoint->Fill( fHeights[iH], Data->fDewPoint_K.at( iH ) );
            averageMixingRatio->Fill( fHeights[iH], Data->fMixingRatio_gkg.at( iH ) );
            averageRelativeHumidity->Fill( fHeights[iH], Data->fRelativeHumidity.at( iH ) );
            averageIndexofRefraction->Fill( fHeights[iH], Data->fIndexofRefraction.at( iH ) );
        }//heights
    }//data sets
    
    // fill average values from profile histograms
    Average->setdefaultvalues( fHeights.size() );
    for( unsigned int i = 0; i < fHeights.size() ; i++ )
    {
        Average->fHeight_m[i] = fHeights[i];
        if( averagePressure->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fPressure_Pa[i] = averagePressure->GetBinContent( i + 1 );
            Average->fRMS_Pressure_Pa[i] = averagePressure->GetBinError( i + 1 );
        }
        if( averageDensity->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fDensity_gcm3[i] = averageDensity->GetBinContent( i + 1 );
            Average->fRMS_Density_gcm3[i] = averageDensity->GetBinError( i + 1 );
        }
        if( averageTemperature->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fTemperature_K[i] = averageTemperature->GetBinContent( i + 1 );
            Average->fRMS_Temperature_K[i] = averageTemperature->GetBinError( i + 1 );
        }
        if( averageDewpoint->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fDewPoint_K[i] = averageDewpoint->GetBinContent( i + 1 );
        }
        if( averageMixingRatio->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fMixingRatio_gkg[i] = averageMixingRatio->GetBinContent( i + 1 );
        }
        if( averageRelativeHumidity->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fRelativeHumidity[i] = averageRelativeHumidity->GetBinContent( i + 1 );
            Average->fRMS_RelHum[i] = averageRelativeHumidity->GetBinError( i + 1 );
            
        }
        if( averageIndexofRefraction->GetBinEntries( i + 1 ) >= nMinFlights )
        {
            Average->fIndexofRefraction[i] = averageIndexofRefraction->GetBinContent( i + 1 );
            Average->fRMS_IndexofRefraction[i] = averageIndexofRefraction->GetBinError( i + 1 );
        }
    }
    
    fillAtmosphericThickness( Average );
    fAverageProfile.push_back( Average );
    
    // remove all histograms
    averagePressure->Delete();
    averageDensity->Delete();
    averageTemperature->Delete();
    averageMixingRatio->Delete();
    averageRelativeHumidity->Delete();
    averageDewpoint->Delete();
    averageIndexofRefraction->Delete();
    
    return fAverageProfile.size() ;
}

/*
 * check if data points are inside a certain
 * years / months / days / hours / mjds range
 *
 * note: e.g. years and months are treated separatedly
 */
bool VAtmosphereSoundings::isDateInRange( VAtmosphereSoundingData* Data, vector<int>* years, vector<int>* months, vector<int>* days, vector<int>* hours, vector<double>* mjds , unsigned int nMinPoints = 20 )
{

    if( Data->fHeight_m.size() < nMinPoints )
    {
        return false;
    }
    if( years && std::find( years->begin(), years->end(), Data->Year ) == years->end() )
    {
        return false;
    }
    if( months && std::find( months->begin(), months->end(), Data->Month ) == months->end() )
    {
        return false;
    }
    if( days && std::find( days->begin(), days->end(), Data->Day ) == days->end() )
    {
        return false;
    }
    if( hours && std::find( hours->begin(), hours->end(), Data->Hour ) == hours->end() )
    {
        return false;
    }
    if( mjds && std::find( mjds->begin(), mjds->end(), Data->MJD ) == mjds->end() )
    {
        return false;
    }
    return true;
}

bool VAtmosphereSoundings::isDateInRange( VAtmosphereSoundingData* Data, vector<double>* mjds, unsigned int nMinPoints = 20 )
{

    if( Data->fHeight_m.size() < nMinPoints )
    {
        return false;
    }
    if( mjds && std::find( mjds->begin(), mjds->end(), Data->MJD ) == mjds->end() )
    {
        return false;
    }
    return true;
}

bool VAtmosphereSoundings::isDateInRange( VAtmosphereSoundingData* Data, double minMJD, double maxMJD, unsigned int nMinPoints = 20 )
{

    if( Data->fHeight_m.size() < nMinPoints )
    {
        return false;
    }
    if( Data->MJD < minMJD || Data->MJD >= maxMJD )
    {
        return false;
    }
    return true;
}

VAtmosphereSoundingData* VAtmosphereSoundings::makeDefaultWinterAtmosphere( string name, string opt, int year )
{
    return makeDefaultAtmosphere( "Winter", name, opt, year );
}

VAtmosphereSoundingData* VAtmosphereSoundings::makeDefaultSummerAtmosphere( string name, string opt, int year )
{
    return makeDefaultAtmosphere( "Summer", name, opt, year );
}

VAtmosphereSoundingData* VAtmosphereSoundings::makeDefaultAtmosphere( string season, string name, string opt, int year )
{
    //do not re-make atmosphere unless forced
    if( opt != "redo" && opt != "force" && name != "" )
    {
        for( unsigned int i = 0; i < fAverageProfile.size(); i++ )
        {
            if( fAverageProfile[i]->Name == name )
            {
                return fAverageProfile[i];
            }
        }
    }
    
    // note: VERITAS seasons (without transition months)
    vector<int> months;
    if( season == "Winter" )
    {
        months.push_back( 1 );
        months.push_back( 2 );
        months.push_back( 3 );
        months.push_back( 12 );
    }
    else
    {
        months.push_back( 6 );
        months.push_back( 9 );
    }
    
    if( year == 0 )
    {
        if( push_average_atmosphere( name, 0, &months, 0, 0, 0, 20, 5 ) )
        {
            return fAverageProfile.back();
        }
    }
    else
    {
        vector< int > years;
        years.push_back( year );
        if( push_average_atmosphere( name, &years, &months, 0, 0, 0, 20, 5 ) )
        {
            return fAverageProfile.back();
        }
    }
    
    cout << "VAtmosphereSoundings::makeDefaultWinterAtmosphere: Error: average atmosphere not pushed" << endl;
    return 0;
}


VAtmosphereSoundingData* VAtmosphereSoundings::makeMeanMonthlyAtmosphere( int month, string name = "", string opt = "use", int yearMin = 1980, int yearMax = 2020 )
{
    //do not re-make atmosphere unless forced
    if( opt != "redo" && opt != "force" && name != "" )
    {
        for( unsigned int i = 0; i < fAverageProfile.size(); i++ )
        {
            if( fAverageProfile[i]->Name == name )
            {
                return fAverageProfile[i];
            }
        }
    }
    
    vector<int> months;
    months.push_back( month );
    
    vector<int> years;
    for( int y = yearMin; y <= yearMax; y++ )
    {
        years.push_back( y );
    }
    
    
    if( push_average_atmosphere( name, &years, &months, 0, 0, 0, 20, 5 ) )
    {
        return fAverageProfile.back();
    }
    cout << "VAtmosphereSoundings::makeMeanMonthlyAtmosphere: Error: average atmosphere not pushed" << endl;
    return 0;
    
}

VAtmosphereSoundingData* VAtmosphereSoundings::makeOneFlightAtmosphere( int year, int month, int day, int hour, string name = "", string opt = "use" )
{
    //do not re-make atmosphere unless forced
    if( opt != "redo" && opt != "force" && name != "" )
    {
        for( unsigned int i = 0; i < fAverageProfile.size(); i++ )
        {
            if( fAverageProfile[i]->Name == name )
            {
                return fAverageProfile[i];
            }
        }
    }
    
    vector<int> hours;
    hours.push_back( hour );
    
    vector<int> days;
    days.push_back( day );
    
    vector<int> months;
    months.push_back( month );
    
    vector<int> years;
    years.push_back( year );
    
    
    
    if( push_average_atmosphere( name, &years, &months, &days, &hours, 0, 20, 0 ) )
    {
        return fAverageProfile.back();
    }
    cout << "VAtmosphereSoundings::makeOneFlightAtmosphere: Error: average atmosphere not pushed" << endl;
    return 0;
    
}

VAtmosphereSoundingData* VAtmosphereSoundings::makeMeanAtmosphereMJD(
    double minMJD, double maxMJD, string name = "", string opt = "" )
{
    //do not re-make atmosphere unless forced
    if( opt != "redo" && opt != "force" && name != "" )
    {
        for( unsigned int i = 0; i < fAverageProfile.size(); i++ )
        {
            if( fAverageProfile[i]->Name == name )
            {
                return fAverageProfile[i];
            }
        }
    }
    
    vector<double> mjds;
    for( double i = 0.5 * TMath::Ceil( minMJD * 2.0 ); i <= maxMJD; i += 0.5 )
    {
        mjds.push_back( i );
    }
    
    if( push_average_atmosphere( name, 0, 0, 0, 0, &mjds, 20, 0 ) )
    {
        return fAverageProfile.back();
    }
    cout << "VAtmosphereSoundings::maneMeanAtmosphereMJD: Error: average atmosphere not pushed" << endl;
    return 0;
    
}


VAtmosphereSoundingData* VAtmosphereSoundings::make_interpolated_atmosphere( VAtmosphereSoundingData* RawData )
{
    VAtmosphereSoundingData* Data = new VAtmosphereSoundingData();
    Data->MJD = RawData->MJD;
    Data->Year = RawData->Year;
    Data->Month = RawData->Month;
    Data->Day = RawData->Day;
    Data->Hour = RawData->Hour;
    Data->Name = RawData->Name;
    
    Data->fHeight_m = fHeights;
    
    
    interpolate( RawData->fPressure_Pa,	RawData->fHeight_m,	Data->fPressure_Pa,		"log" );
    interpolate( RawData->fTemperature_K,	RawData->fHeight_m,	Data->fTemperature_K,		"lin" );
    interpolate( RawData->fDewPoint_K,	RawData->fHeight_m,	Data->fDewPoint_K,		"lin" );
    interpolate( RawData->fRelativeHumidity,	RawData->fHeight_m,	Data->fRelativeHumidity,	"lin" );
    interpolate( RawData->fMixingRatio_gkg,	RawData->fHeight_m,	Data->fMixingRatio_gkg,		"lin" );
    interpolate( RawData->fIndexofRefraction, RawData->fHeight_m, 	Data->fIndexofRefraction,	"lin" );
    
    fillAtmosphericDensity( Data );
    fillAtmosphericThickness( Data );
    
    
    
    return Data;
    
}

/*
 * return residuals to a model graph
 *
 * values are given in %
*/
TGraph* VAtmosphereSoundings::getResidualGraph( TGraph* data, TGraph* model, int color )
{
    TGraph* newgraph = new TGraph( 0 );
    TString title = TString::Format( "res_%s_%s", data->GetName(), model->GetName() );
    newgraph->SetNameTitle( title.Data(), title.Data() );
    newgraph->SetLineColor( color );
    newgraph->SetMarkerColor( color );
    newgraph->SetMarkerStyle( 20 );
    newgraph->SetMarkerSize( 1.5 );
    
    for( int iA = 0; iA < data->GetN(); iA++ )
    {
        newgraph->SetPoint( iA , data->GetX()[iA], 100.*( data->GetY()[iA] / model->Eval( data->GetX()[iA] ) - 1.0 ) );
    }
    
    return newgraph;
}

/*
 * return color codes for different seasons
 *
 * winter (61) colors are blueish
 * summer (62) colors are redish
 * intermediate (99) colors are greenish
 */
Color_t VAtmosphereSoundings::getSeasonColor( int counter, int atmo )
{
    Color_t season_color = kGreen;
    if( atmo == 61 )
    {
        season_color = kBlue;
    }
    else if( atmo == 62 )
    {
        season_color = kRed;
    }
    
    return season_color - 10 + 2 * counter;
}

Color_t VAtmosphereSoundings::getSeasonColor( int iMonth )
{
    map< int, Color_t > season_color;
    
    if( iMonth > 12 )
    {
        iMonth -= 12;
    }
    
    // winter
    season_color[12] = kBlue + 2;
    season_color[1] = kBlue - 3;
    season_color[2] = kBlue - 7;
    season_color[3] = kBlue - 10;
    // summer
    season_color[6] = kRed + 2;
    season_color[7] = kRed - 3;
    season_color[8] = kRed - 7;
    season_color[9] = kRed - 10;
    // intermediate
    season_color[4] = kGreen + 2;
    season_color[5] = kGreen - 3;
    season_color[10] = kGreen - 7;
    season_color[11] = kGreen - 10;
    
    if( season_color.find( iMonth ) != season_color.end() )
    {
        return season_color[iMonth];
    }
    return 1;
}

/*
 * plot density profiles for a particular season
 *
 */
TCanvas* VAtmosphereSoundings::plot_season(
    vector<VAtmosphereSoundingData*> v,
    TString season_name, string value, TString outfileprefix = "" )
{

    // average winter / summer (or from CORSIKA)
    
    VAtmosphereSoundingData* winter = 0;
    VAtmosphereSoundingData* summer = 0;
    if( read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof61.dat", "Winter", 4 )
            && getDataMODTRAN( 0 ) )
    {
        winter = getDataMODTRAN( 0 );
    }
    if( read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof62.dat", "Summer", 2 )
            && getDataMODTRAN( 1 ) )
    {
        summer = getDataMODTRAN( 1 );
    }
    if( !summer || !winter )
    {
        cout << "Error defining summer / winter atmosphere" << endl;
        return 0;
    }
    
    // VAtmosphereSoundingData* summer = makeDefaultSummerAtmosphere( "summer", "" );
    // VAtmosphereSoundingData* winter = makeDefaultWinterAtmosphere( "winter", "" );
    
    TGraph* temp = 0;
    
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////
    // value vs height plotting
    sprintf( hname, "c_%s", season_name.Data() );
    sprintf( htitle, "season %s", season_name.Data() );
    TCanvas* c = new TCanvas( hname, htitle, 0, 0, 800, 500 );
    c->SetGridx( 0 );
    c->SetGridy( 0 );
    c->SetLeftMargin( 0.13 );
    
    winter->setColor( 4 );
    summer->setColor( 2 );
    winter->getGraph( value )->SetLineWidth( 6 );
    summer->getGraph( value )->SetLineWidth( 6 );
    winter->getGraph( value )->SetLineStyle( 2 );
    summer->getGraph( value )->SetLineStyle( 2 );
    winter->getGraph( value )->SetTitle( "" );
    summer->getGraph( value )->SetTitle( "" );
    winter->getGraph( value )->Draw( "al" );
    TAxis* axis = winter->getGraph( value )->GetXaxis();
    axis->SetLimits( fPlottingHeight_min, fPlottingHeight_max );
    axis->SetTitle( "height a.s.l. [km]" );
    TAxis* axisY = winter->getGraph( value )->GetYaxis();
    axisY->SetTitleOffset( 1.5 );
    summer->getGraph( value )->Draw( "l same" );
    
    for( unsigned int i = 0; i < v.size(); i++ )
    {
        if( v[i] )
        {
            v[i]->getGraph( value )->Draw( "lx" );
        }
    }
    
    TLegend* legend = 0;
    if( value == "temperature" )
    {
        legend = gPad->BuildLegend( 0.35, 0.45, 1.0, 1.0 );
        winter->getGraph( value )->GetYaxis()->SetTitle( "Temperature [K]" );
    }
    else if( value == "density" )
    {
        legend = gPad->BuildLegend( 0.2, 0.13, 0.7, 0.5 );
        winter->getGraph( value )->GetYaxis()->SetTitle( "density*exp(h/7.739km) [g/cm^{3}]" );
        winter->getGraph( value )->SetMinimum( 0.0008 );
        winter->getGraph( value )->SetMaximum( 0.0017 );
    }
    else
    {
        legend = gPad->BuildLegend();
    }
    // legend->SetTextSize( 10 );
    legend->SetFillStyle( 0 );
    legend->SetNColumns( 2 );
    legend->SetBorderSize( 0 );
    legend->SetDrawOption();
    
    // plot again CORSIKA atmospheres
    winter->getGraph( value )->Draw( "l" );
    summer->getGraph( value )->Draw( "l" );
    
    ///////////////////////////////////
    // residual graphs to average winter
    sprintf( hname, "cRW_%s", season_name.Data() );
    sprintf( htitle, "season %s", season_name.Data() );
    TCanvas* cRW = new TCanvas( hname, htitle, 600, 0, 800, 500 );
    cRW->SetGridx( 0 );
    cRW->SetGridy( 0 );
    
    temp = getResidualGraph( summer->getGraph( value ), winter->getGraph( value ), 2 );
    temp->SetLineWidth( 5 );
    temp->SetLineStyle( 4 );
    temp->SetMinimum( -15. );
    temp->SetMaximum( 15. );
    temp->SetTitle( "" );
    temp->Draw( "al" );
    temp->GetYaxis()->SetTitle( "rel diff to winter (%)" );
    axis = temp->GetXaxis();
    axis->SetLimits( fPlottingHeight_min, fPlottingHeight_max );
    axis->SetTitle( "height a.s.l. [km]" );
    
    temp = getResidualGraph( winter->getGraph( value ), winter->getGraph( value ), 4 );
    temp->SetLineWidth( 2 );
    temp->Draw( "l same" );
    
    for( unsigned int i = 0; i < v.size(); i++ )
    {
        TGraph* temp = getResidualGraph( v[i]->getGraph( value ),
                                         winter->getGraph( value ), v[i]->getGraph( value )->GetMarkerColor() );
        if( temp )
        {
            temp->SetLineWidth( 2 );
            temp->Draw( "l same" );
        }
    }
    
    ///////////////////////////////////////////
    // residual graphs to average summer
    sprintf( hname, "cRS_%s", season_name.Data() );
    sprintf( htitle, "season %s", season_name.Data() );
    TCanvas* cRS = new TCanvas( hname, htitle, 1200, 0, 800, 500 );
    cRS->SetGridx( 0 );
    cRS->SetGridy( 0 );
    
    temp = getResidualGraph( winter->getGraph( value ), summer->getGraph( value ), 4 );
    temp->SetLineWidth( 5 );
    temp->SetLineStyle( 4 );
    temp->SetTitle( "" );
    temp->SetMinimum( -15. );
    temp->SetMaximum( 15. );
    temp->Draw( "al" );
    axis = temp->GetXaxis();
    axis->SetLimits( fPlottingHeight_min, fPlottingHeight_max );
    axis->SetTitle( "height a.s.l. [km]" );
    temp->GetYaxis()->SetTitle( "rel diff to summer (%)" );
    temp = getResidualGraph( summer->getGraph( value ), summer->getGraph( value ), 2 );
    temp->Draw( "l same" );
    temp->SetLineWidth( 4 );
    
    for( unsigned int i = 0; i < v.size(); i++ )
    {
        TGraph* temp = getResidualGraph( v[i]->getGraph( value ),
                                         summer->getGraph( value ), v[i]->getGraph( value )->GetMarkerColor() );
        if( temp )
        {
            temp->SetLineWidth( 2 );
            temp->Draw( "l same" );
        }
    }
    
    //////////////////////////////////////////////////////////////////
    // print everything into pdfs
    if( gSystem->mkdir( "figures", true ) )
    {
        TString filename = TString::Format( "figures/%sseason-%s-%s.pdf",
                                            outfileprefix.Data(), season_name.Data(), value.c_str() );
        c->SaveAs( filename.Data() );
        filename = TString::Format( "figures/%sseason-relativeWinter-%s-%s.pdf",
                                    outfileprefix.Data(), season_name.Data(), value.c_str() );
        cRW->SaveAs( filename.Data() );
        filename = TString::Format( "figures/%sseason-relativeSummer-%s-%s.pdf",
                                    outfileprefix.Data(), season_name.Data(), value.c_str() );
        cRS->SaveAs( filename.Data() );
    }
    
    return c;
}

/*
   plots and prints temperature/density profiles per dark run for a given observing season.
   expects the dates of the last full moon before the start of season and the first full moon after the end of season.
   value = "temperature" or "density" (or see VAtmosphereSoundingData::getGraph( ... ) ).

   example:

   gSystem->Load("libVAnaSum.so");
   gStyle->SetOptTitle(0);
   VAtmosphereSoundings * a = new VAtmosphereSoundings();
   a->setHeights(0.0, 30000.0, 1000.0);
   a->readSoundingsFromTextFile("list.dat");

   a->plot_season(( 56524, 56820, "2013/14", "temperature");
   a->plot_season( 2006, 9, 7, 2007, 6, 1, "temperature");


*/
TCanvas* VAtmosphereSoundings::plot_season(
    int year_start, int month_start, int day_start,
    int year_end, int month_end , int day_end,
    string value,
    int bWriteCorsika )
{

    double mjd_start = 0.;
    double mjd_end = 0.;
    int j = 0;
    TString season_name = TString::Format( "%d-%d", year_start, year_end );
    VAstronometry::vlaCldj( year_start, month_start , day_start, &mjd_start, &j );
    VAstronometry::vlaCldj( year_end, month_end , day_end, &mjd_end, &j );
    
    // average length of synodic month
    double month = 29.53;
    int Nmonth = TMath::CeilNint( ( mjd_end - mjd_start ) / month ) ;
    
    vector< VAtmosphereSoundingData* > v;
    
    // color counters
    int w = 0;
    int s = 0;
    int g = 0;
    
    // prepare average montly atmospheres
    for( int i = 0; i < Nmonth; i++ )
    {
        int y1, y2, m1, m2, d1, d2, j;
        double f = 0.;
        double start = mjd_start + i * month + 1;
        double end = mjd_start + ( i + 1 ) * month - 1;
        
        VAstronometry::vlaDjcl( start, &y1, &m1, &d1, &f, &j );
        VAstronometry::vlaDjcl( end, &y2, &m2, &d2, &f, &j );
        TString name = TString::Format( "%d-%d-%d - %d-%d-%d", y1, m1, d1, y2, m2, d2 );
        VAtmosphereSoundingData* t =
            makeMeanAtmosphereMJD( start, end, name.Data(), name.Data() );
            
        // get epoch / pre-defined summer/winter atmosphere
        TString name_start = TString::Format( "%d-%d-%d", y1, m1, d1 );
        int atm = readEpochsAndAtmospheres( name_start.Data(), month );
        cout << name_start << "\t" << atm << "\t" << endl;
        
        if( t && t->getGraph( value ) && t->getGraph( value )->GetN() > 0 )
        {
            if( atm < 0 )
            {
                t->setColor( getSeasonColor( g, atm ) );
                cout << "\t GREEN " << getSeasonColor( 9 + g ) << endl;
                g++;
            }
            else if( atm == 61 )
            {
                t->setColor( getSeasonColor( w, atm ) );
                cout << "\t Blue " << getSeasonColor( w, atm ) << endl;
                w++;
            }
            else if( atm == 62 )
            {
                t->setColor( getSeasonColor( s, atm ) );
                cout << "\t Red " << getSeasonColor( s, atm ) << endl;
                s++;
            }
            t->getGraph( value )->SetTitle( name.Data() );
            if( m1 >= 6. && m1 <= 9. )
            {
                t->getGraph( value )->SetMarkerStyle( 24 );
            }
            v.push_back( t );
        }
    }
    
    // write CORSIKA files (if requested)
    if( bWriteCorsika > 0 )
    {
        for( unsigned int i = 0; i < v.size(); i++ )
        {
            if( !v[i] )
            {
                continue;
            }
            stringstream i_fname;
            i_fname << "ATM" << ( bWriteCorsika + i ) << ".dat";
            cout << "Writing CORSIKA file: " << endl;
            cout << "\t" << i_fname.str() << endl;
            v[i]->write_CORSIKA_UserProfile( bWriteCorsika + i, i_fname.str() );
        }
    }
    
    return plot_season( v, season_name, value );
}

/*
 * plot montly averages for the years given
 *
 */
TCanvas* VAtmosphereSoundings::plot_monthly(
    vector< int > year, vector< int > month, vector< int > day,
    double intervall_days, int offset_months,
    string value )
{
    TString season_name = TString::Format( "%d", offset_months );
    
    double mjd_start = 0.;
    int j = 0;
    
    vector<VAtmosphereSoundingData*> v;
    
    vector<double> start, end;
    for( unsigned int i = 0; i < year.size(); i++ )
    {
        VAstronometry::vlaCldj( year[i], month[i], day[i], &mjd_start, &j );
        
        start.push_back( mjd_start + intervall_days * offset_months );
        end.push_back( mjd_start + intervall_days * ( offset_months + 1 ) );
        
        int y1, y2, m1, m2, d1, d2, j;
        double f = 0.;
        VAstronometry::vlaDjcl( start[i] , &y1, &m1, &d1, &f, &j );
        VAstronometry::vlaDjcl( end[i] , &y2, &m2, &d2, &f, &j );
        TString output = TString::Format( "%d\t%d-%d-%d\t%d-%d-%d", i + 1, y1, m1, d1, y2, m2, d2 );
        cout << output << endl;
    }
    
    // prepare average montly atmospheres
    for( unsigned int i = 0; i < start.size(); i++ )
    {
        int y1, y2, m1, m2, d1, d2, j;
        double f = 0.;
        VAstronometry::vlaDjcl( start[i] , &y1, &m1, &d1, &f, &j );
        VAstronometry::vlaDjcl( end[i] , &y2, &m2, &d2, &f, &j );
        TString name = TString::Format( "%d-%d-%d - %d-%d-%d", y1, m1, d1, y2, m2, d2 );
        VAtmosphereSoundingData* t =
            makeMeanAtmosphereMJD( start[i], end[i], name.Data(), name.Data() );
        if( t && t->getGraph( value ) && t->getGraph( value )->GetN() > 0 )
        {
            t->setColor( i + 1 );
            t->getGraph( value )->SetTitle( name.Data() );
            if( m1 >= 6. && m1 <= 9. )
            {
                t->getGraph( value )->SetMarkerStyle( 24 );
            }
            v.push_back( t );
        }
    }
    
    return plot_season( v, season_name, value );
}

/*
 *  read instrument epoch and atmospheric ID from a parameter file
 *
 *  (this might be extended in future to read the values from the db)
 *
 *  an epoch might be summer or winter
 *
 */
int VAtmosphereSoundings::readEpochsAndAtmospheres(
    TString iDstart, double iMonthLength_days,
    string iEpochFile )
{
    ifstream is;
    char* iFileName = gSystem->ExpandPathName( iEpochFile.c_str() );
    is.open( iFileName, ifstream::in );
    if( !is )
    {
        cout << "error opening epoch parameter file " << iEpochFile << endl;
        return -1;
    }
    string is_line;
    string temp;
    
    // start and end date of period
    TDatime d_min( ( iDstart + " 12:00:00" ).Data() );
    TDatime d_max;
    d_max.Set( d_min.Convert( kTRUE ) + iMonthLength_days * 86400. );
    
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
            if( ( is_stream >> std::ws ).eof() )
            {
                continue;
            }
            is_stream >> temp;
            if( temp == "ATMOSPHERE" )
            {
                string atmo = "";
                string date_min = "";
                string date_max = "";
                if( !( is_stream >> std::ws ).eof() )
                {
                    is_stream >> atmo;
                }
                if( !( is_stream >> std::ws ).eof() )
                {
                    is_stream >> date_min;
                }
                if( !( is_stream >> std::ws ).eof() )
                {
                    is_stream >> date_max;
                }
                TDatime e_min( ( date_min + " 12:00:00" ).c_str() );
                TDatime e_max( ( date_max + " 12:00:00" ).c_str() );
                if( d_min.Get() >= e_min.Get() && d_max.Get() < e_max.Get() )
                {
                    return atoi( atmo.c_str() );
                }
                // overlap minus 1 week
                else if( d_min.Get() >= e_min.Get() && d_max.Get() < e_max.Get() + 7.*86400. )
                {
                    return atoi( atmo.c_str() );
                }
                else if( d_min.Get() >= e_min.Get() - 7.*86400. && d_max.Get() < e_max.Get() )
                {
                    return atoi( atmo.c_str() );
                }
            }
        }
    }
    is.close();
    
    return -1;
}
