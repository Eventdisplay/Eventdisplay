/*! \brief VAtmosphereSoundings read and analyse sounding data

//////////////////////////////////////////////////////////////////////////////////////

   reading ascii files with sounding data for Tucson

   Data from: http://weather.uwyo.edu/upperair/sounding.html

    VAtmosphereSoundings a;
    a.readSoundingsFromTextFile("allBalloon.txt");
    a.writeRootFile("ballonDataTucson_199501_201007.root");


//////////////////////////////////
   FIX ATMOSPHERIC THICKNESS

*/

#include "VAtmosphereSoundingData.h"


VAtmosphereSoundingData::VAtmosphereSoundingData()
{
    MJD = 0;
    Year = 0;
    Month = 0;
    Day = 0;
    Hour = 0.;
    Name = "";
    
    PlotColor = 1;
    PlotMarker = 20;
    PlotLineStyle = 1;
    PlotLineWidth = 2;
    PlotMarkerSize = 0.5;
    
    fGraphScaledDensityHeight = 0;
    fGraphPressureHeight = 0;
    fGraphHumidityHeight = 0;
    fGraphTemperatureHeight = 0;
    fGraphIndexHeight = 0;
    fGraphThicknessHeight = 0;
}

void VAtmosphereSoundingData::setdefaultvalues( unsigned int iN )
{
    fPressure_Pa.clear();
    fHeight_m.clear();
    fDensity_gcm3.clear();
    fThickness_gcm2.clear();
    fTemperature_K.clear();
    fDewPoint_K.clear();
    fRelativeHumidity.clear();
    fVaporMassDensity_gm3.clear();
    fMixingRatio_gkg.clear();
    fWindDirection_deg.clear();
    fWindSpeed_ms.clear();
    fIndexofRefraction.clear();
    fO2_cmkm.clear();
    fO3_cmkm.clear();
    
    if( fGraphScaledDensityHeight != 0 )
    {
        fGraphScaledDensityHeight->Delete();
    }
    
    // set default values
    for( unsigned int i = 0; i < iN; i++ )
    {
        fPressure_Pa.push_back( -9999. );
        fHeight_m.push_back( -9999. );
        fDensity_gcm3.push_back( -9999. );
        fThickness_gcm2.push_back( -9999. );
        fTemperature_K.push_back( -9999. );
        fDewPoint_K.push_back( -9999. );
        fRelativeHumidity.push_back( -9999. );
        fVaporMassDensity_gm3.push_back( -9999. );
        fMixingRatio_gkg.push_back( -9999. );
        fWindDirection_deg.push_back( -9999. );
        fWindSpeed_ms.push_back( -9999. );
        fIndexofRefraction.push_back( -9999. );
        fO2_cmkm.push_back( -9999. );
        fO3_cmkm.push_back( -9999. );
        fRMS_Pressure_Pa.push_back( 0. );
        fRMS_Density_gcm3.push_back( 0. );
        fRMS_Temperature_K.push_back( 0. );
        fRMS_RelHum.push_back( 0. );
        fRMS_IndexofRefraction.push_back( 2. );
    }
}




void VAtmosphereSoundingData::makeGraphScaledDensity( )
{
    if( fGraphScaledDensityHeight != 0 )
    {
        fGraphScaledDensityHeight->Delete();
    }
    
    fGraphScaledDensityHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( fDensity_gcm3[i] > 0 )
        {
            fGraphScaledDensityHeight->SetPoint( fGraphScaledDensityHeight->GetN(), fHeight_m[i] / 1000.0, fDensity_gcm3[i] * TMath::Exp( fHeight_m[i] / 7739.0 ) );
            if( i < fRMS_Density_gcm3.size() )
            {
                fGraphScaledDensityHeight->SetPointError( fGraphScaledDensityHeight->GetN() - 1, 0, fRMS_Density_gcm3[i] );
            }
        }
    }
    
    fGraphScaledDensityHeight->SetMarkerStyle( PlotMarker );
    fGraphScaledDensityHeight->SetMarkerSize( PlotMarkerSize );
    fGraphScaledDensityHeight->SetMarkerColor( PlotColor );
    fGraphScaledDensityHeight->SetLineColor( PlotColor );
    fGraphScaledDensityHeight->SetLineStyle( PlotLineStyle );
    fGraphScaledDensityHeight->SetLineWidth( PlotLineWidth );
    fGraphScaledDensityHeight->SetTitle( "Scaled density;height [km]; density*exp(h/7.739km) [g/cm^{3}]" );
    TString temp = "density_" + Name;
    fGraphScaledDensityHeight->SetName( temp.Data() );
}


void VAtmosphereSoundingData::makeGraphPressure( )
{
    if( fGraphPressureHeight != 0 )
    {
        fGraphPressureHeight->Delete();
    }
    
    fGraphPressureHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( fPressure_Pa[i] > 0 )
        {
            fGraphPressureHeight->SetPoint( fGraphPressureHeight->GetN(), fHeight_m[i] / 1000.0, fPressure_Pa[i] );
            if( i < fRMS_Pressure_Pa.size() )
            {
                fGraphPressureHeight->SetPointError( fGraphPressureHeight->GetN() - 1, 0, fRMS_Pressure_Pa[i] );
            }
        }
    }
    
    fGraphPressureHeight->SetMarkerStyle( PlotMarker );
    fGraphPressureHeight->SetMarkerSize( PlotMarkerSize );
    fGraphPressureHeight->SetMarkerColor( PlotColor );
    fGraphPressureHeight->SetLineColor( PlotColor );
    fGraphPressureHeight->SetLineStyle( PlotLineStyle );
    fGraphPressureHeight->SetLineWidth( PlotLineWidth );
    fGraphPressureHeight->SetTitle( "Pressure;height [km]; pressure [Pa]" );
    TString temp = "pressure_" + Name;
    fGraphPressureHeight->SetName( temp.Data() );
}


void VAtmosphereSoundingData::makeGraphTemperature( )
{
    if( fGraphTemperatureHeight )
    {
        fGraphTemperatureHeight->Delete();
    }
    
    fGraphTemperatureHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( fTemperature_K[i] > 0 )
        {
            fGraphTemperatureHeight->SetPoint( fGraphTemperatureHeight->GetN(), fHeight_m[i] / 1000.0, fTemperature_K[i] );
            if( i < fRMS_Temperature_K.size() )
            {
                fGraphTemperatureHeight->SetPointError( fGraphTemperatureHeight->GetN() - 1, 0, fRMS_Temperature_K[i] );
            }
            
        }
    }
    
    fGraphTemperatureHeight->SetMarkerStyle( PlotMarker );
    fGraphTemperatureHeight->SetMarkerSize( PlotMarkerSize );
    fGraphTemperatureHeight->SetMarkerColor( PlotColor );
    fGraphTemperatureHeight->SetLineColor( PlotColor );
    fGraphTemperatureHeight->SetLineStyle( PlotLineStyle );
    fGraphTemperatureHeight->SetLineWidth( PlotLineWidth );
    fGraphTemperatureHeight->SetTitle( "Temperature;height [km]; temperature [K]" );
    TString temp = "temperature_" + Name;
    fGraphTemperatureHeight->SetName( temp.Data() );
}


void VAtmosphereSoundingData::makeGraphHumidity( )
{
    if( fGraphHumidityHeight != 0 )
    {
        fGraphHumidityHeight->Delete();
    }
    
    fGraphHumidityHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( fRelativeHumidity[i] > 0 )
        {
            fGraphHumidityHeight->SetPoint( fGraphHumidityHeight->GetN(), fHeight_m[i] / 1000.0, fRelativeHumidity[i] );
            
            if( i < fRMS_RelHum.size() )
            {
                fGraphHumidityHeight->SetPointError( fGraphHumidityHeight->GetN() - 1, 0, fRMS_RelHum[i] );
                
            }
        }
    }
    
    fGraphHumidityHeight->SetMarkerStyle( PlotMarker );
    fGraphHumidityHeight->SetMarkerSize( PlotMarkerSize );
    fGraphHumidityHeight->SetMarkerColor( PlotColor );
    fGraphHumidityHeight->SetLineColor( PlotColor );
    fGraphHumidityHeight->SetLineStyle( PlotLineStyle );
    fGraphHumidityHeight->SetLineWidth( PlotLineWidth );
    fGraphHumidityHeight->SetTitle( "Humidity;height [km]; rel. humidity [\%]" );
    TString temp = "relHum_" + Name;
    fGraphHumidityHeight->SetName( temp.Data() );
}

void VAtmosphereSoundingData::makeGraphThickness( )
{
    if( fGraphThicknessHeight != 0 )
    {
        fGraphThicknessHeight->Delete();
    }
    
    fGraphThicknessHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( fThickness_gcm2[i] > 0 )
        {
            fGraphThicknessHeight->SetPoint( fGraphThicknessHeight->GetN(), fHeight_m[i] / 1000.0, fThickness_gcm2[i] );
        }
    }
    
    fGraphThicknessHeight->SetMarkerStyle( PlotMarker );
    fGraphThicknessHeight->SetMarkerSize( PlotMarkerSize );
    fGraphThicknessHeight->SetMarkerColor( PlotColor );
    fGraphThicknessHeight->SetLineColor( PlotColor );
    fGraphThicknessHeight->SetLineStyle( PlotLineStyle );
    fGraphThicknessHeight->SetLineWidth( PlotLineWidth );
    fGraphThicknessHeight->SetTitle( "Thickness;height [km]; thickness [g/cm^{2}]" );
    TString temp = "thickness_" + Name;
    fGraphThicknessHeight->SetName( temp.Data() );
}

void VAtmosphereSoundingData::makeGraphIndex( )
{
    if( fGraphIndexHeight != 0 )
    {
        fGraphIndexHeight->Delete();
    }
    
    fGraphIndexHeight = new TGraphErrors( 0 );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
    
        if( fIndexofRefraction[i] > 0 )
        {
            fGraphIndexHeight->SetPoint( fGraphIndexHeight->GetN(), fHeight_m[i] / 1000.0, fIndexofRefraction[i] - 1 );
            
            if( i < fRMS_IndexofRefraction.size() )
            {
                fGraphIndexHeight->SetPointError( fGraphIndexHeight->GetN() - 1, 0, fRMS_IndexofRefraction[i] );
            }
        }
    }
    
    fGraphIndexHeight->SetMarkerStyle( PlotMarker );
    fGraphIndexHeight->SetMarkerSize( PlotMarkerSize );
    fGraphIndexHeight->SetMarkerColor( PlotColor );
    fGraphIndexHeight->SetLineColor( PlotColor );
    fGraphIndexHeight->SetLineStyle( PlotLineStyle );
    fGraphIndexHeight->SetLineWidth( PlotLineWidth );
    fGraphIndexHeight->SetTitle( "Index of Refraction;height [km]; index of refraction -1 ]" );
    TString temp = "n-1_" + Name;
    fGraphIndexHeight->SetName( temp.Data() );
}

void VAtmosphereSoundingData::setColor( int color )
{
    PlotColor = color;
    
    if( fGraphIndexHeight )
    {
        fGraphIndexHeight->SetLineColor( PlotColor );
        fGraphIndexHeight->SetMarkerColor( PlotColor );
    }
    if( fGraphThicknessHeight )
    {
        fGraphThicknessHeight->SetLineColor( PlotColor );
        fGraphThicknessHeight->SetMarkerColor( PlotColor );
    }
    if( fGraphHumidityHeight )
    {
        fGraphHumidityHeight->SetLineColor( PlotColor );
        fGraphHumidityHeight->SetMarkerColor( PlotColor );
    }
    if( fGraphTemperatureHeight )
    {
        fGraphTemperatureHeight->SetLineColor( PlotColor );
        fGraphTemperatureHeight->SetMarkerColor( PlotColor );
    }
    if( fGraphPressureHeight )
    {
        fGraphPressureHeight->SetLineColor( PlotColor );
        fGraphPressureHeight->SetMarkerColor( PlotColor );
    }
    if( fGraphScaledDensityHeight )
    {
        fGraphScaledDensityHeight->SetLineColor( PlotColor );
        fGraphScaledDensityHeight->SetMarkerColor( PlotColor );
    }
}

bool VAtmosphereSoundingData::write_2C1( string filename, vector<double>* ModtranHeights, double max_height )
{
    FILE* file;
    file = fopen( filename.data(), "w" );
    if( !file )
    {
        cout << "VAtmosphereSoundings::write_2C1: couldn't open file " << filename << endl;
        return false;
    }
    
    //write card 2C1 for each height step.
    //format:
    // CARD 2C1: ZM, P, T, (WMOL(J), J = 1, 3), (JCHAR(J), J = 1, 14), JCHARX
    // FORMAT (F10.3, 5E10.3, 14A1, 1X, A1)
    
    double have_line = false;
    
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        if( ( fPressure_Pa.at( i ) < 0 ||  fTemperature_K.at( i ) < 0 ) && !have_line )
        {
            continue;    // do not do anything if we are below the probe's starting height.
        }
        if( ( fPressure_Pa.at( i ) < 0 ||  fTemperature_K.at( i ) < 0 ) &&  have_line )
        {
            max_height = fHeight_m.at( i - 1 );
            break;
        }
        if( fHeight_m.at( i ) > max_height )
        {
            break;    //stop writing if there is no data point or we reach the max. height.
        }
        
        have_line = true;
        
        fprintf( file, "%10.3f%10.3f%10.3f",  fHeight_m.at( i ) / 1000.,  fPressure_Pa.at( i ) / 100.,  fTemperature_K.at( i ) );
        fprintf( file, "%10.3f%10.3f%10.3f",  fRelativeHumidity.at( i ), 0., 0. );
        fprintf( file, "AAH              \n" );
        
    }
    
    //now write out the heights at which we want to evaluate the standard atmospheres
    for( unsigned int i = 0; i < ModtranHeights->size(); i++ )
    {
        if( ModtranHeights->at( i ) <= max_height )
        {
            continue;
        }
        fprintf( file, "%10.3f                                                                            \n",  ModtranHeights->at( i ) / 1000. );
    }
    
    fclose( file );
    
    return true;
}
bool VAtmosphereSoundingData::write_CORSIKA_UserProfile( unsigned int atmprofmodel, string iName )
{
    cout << "# Atmospheric Model " << atmprofmodel << " (" << iName << ")" << endl;
    cout << "#Col. #1          #2           #3            #4" << endl;
    cout << "# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1" << endl;
    
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        printf( "%10.3f", fHeight_m[i] / 1.e3 );
        printf( "%15.5e", fDensity_gcm3[i] );
        printf( "%13.5e", fThickness_gcm2[i] );
        printf( "%13.5e", fIndexofRefraction[i] - 1. );
        cout << endl;
    }
    
    string filename = iName.append( ".profile" );
    FILE* out = fopen( filename.c_str(), "w" );
    
    fprintf( out, "# Atmospheric Model %d (%s)\n", atmprofmodel, iName.c_str() );
    fprintf( out, "#Col. #1          #2           #3            #4 \n" );
    fprintf( out, "# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1\n" );
    for( unsigned int i = 0; i < fHeight_m.size(); i++ )
    {
        fprintf( out, "%10.3f", fHeight_m[i] / 1.e3 );
        fprintf( out, "%15.5e", fDensity_gcm3[i] );
        fprintf( out, "%13.5e", fThickness_gcm2[i] );
        fprintf( out, "%13.5e\n", fIndexofRefraction[i] - 1. );
    }
    fclose( out );
    
    return true;
}
