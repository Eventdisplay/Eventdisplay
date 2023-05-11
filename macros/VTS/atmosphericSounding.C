/* \file  atmosphericSounding.C
   \brief plot ballon data from soundings measurements compared with CORSIKA and MODTRAN values

   *** might need some attention ***

*/

/*
      compare yearly data with CORSIKA/MODTRAN atmospheres
*/
void plot_years( string iDataFile )
{
    VAtmosphereSoundings a( iDataFile );
    a.plotAttributes_ColorChange( true );
    a.plotAttributes_PlotLegend( true );
    
    a.read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof61.dat", "VERITAS Winter", 4 );
    a.read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof62.dat", "VERITAS Summer", 2 );
    
    a.setPlottingPeriod( "yearly" );
    a.plotAverages( 2010, 1, 2019, 12 );
}

void plot_seasons( string iDataFile, bool iMonthly = false )
{
    VAtmosphereSoundings a( iDataFile );
    a.read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof61.dat", "VERITAS Winter", 4 );
    a.read_CORSIKA_Atmosphere( "$VERITAS_EVNDISP_AUX_DIR/Atmospheres/atmprof62.dat", "VERITAS Summer", 2 );
    // a.setHeights(0.,30000.,1000.);
    //
    
    // dates of first full moon in September
    vector< int > year;
    vector< int > month;
    vector< int > day;
    
    year.push_back( 2010 );
    month.push_back( 9 );
    day.push_back( 23 );
    year.push_back( 2011 );
    month.push_back( 9 );
    day.push_back( 12 );
    year.push_back( 2012 );
    month.push_back( 8 );
    day.push_back( 31 );
    year.push_back( 2013 );
    month.push_back( 9 );
    day.push_back( 19 );
    year.push_back( 2014 );
    month.push_back( 9 );
    day.push_back( 8 );
    year.push_back( 2015 );
    month.push_back( 9 );
    day.push_back( 27 );
    year.push_back( 2016 );
    month.push_back( 9 );
    day.push_back( 16 );
    year.push_back( 2017 );
    month.push_back( 9 );
    day.push_back( 6 );
    year.push_back( 2018 );
    month.push_back( 9 );
    day.push_back( 2 );
    // year.push_back( 2019 );  month.push_back( 9 ); day.push_back( 13 );
    
    if( iMonthly )
    {
        for( unsigned int i = 0; i < 12; i++ )
        {
            a.plot_monthly( year, month, day, 29.53, i, "density" );
        }
    }
    // plot montly for each year
    else
    {
        for( unsigned int i = 0; i < year.size(); i++ )
        {
            cout << "plotting " << year[i] << endl;
            a.plot_season( year[i], month[i], day[i], year[i] + 1, 7, 1, "density" );
        }
    }
}


