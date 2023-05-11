/*
 * IRF plots:
 *
 * - energy threshold vs zenith angle
 * - plot compare to IRF files (effective area files)
 *
 */

void print_pdf( string iEffAreaFile, string iCanvasName, string iPrintName )
{
    TCanvas* c = 0;
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( iCanvasName.c_str() );
    if( c )
    {
        c->Print( ( iEffAreaFile.substr( 0,  iEffAreaFile.find( ".root" ) ) + iPrintName ).c_str() );
    }
}

/*
 * plot IRFs
 *
 */
void plotIRFs( string iEffAreaFile1, string iEffAreaFile2, float ze = 20., int NSB = 130 )
{
    VPlotInstrumentResponseFunction a;
    a.addInstrumentResponseData( iEffAreaFile1, ze,
                                 0.5, 0, 1.5, NSB );
    if( iEffAreaFile2.size() > 0 )
    {
        a.addInstrumentResponseData( iEffAreaFile2, ze,
                                     0.5, 0, 1.5, NSB );
    }
    a.setPlottingAxis( "energy_Lin", "X", true, 0.03, 100. );
    
    TCanvas* c = 0;
    // effective area
    c = a.plotEffectiveArea( -1., 4.e5 );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "EffArea.pdf" );
    }
    c = a.plotEffectiveAreaRatio( 0, 0., 2. );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "EffAreaRatio.pdf" );
    }
    // angular resolution
    c =  a.plotAngularResolution( "energy", "68", -99, 0.2 );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "AngRes.pdf" );
    }
    // cumulative theta2 plots
    c = a.plotTheta2( 0.02, false );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "AngResT2.pdf" );
    }
    // theta2 plots
    c = a.plotTheta2( 0.05, true );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "AngResT2Cum.pdf" );
    }
    // energy resolution
    c = a.plotEnergyResolution();
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "ERes.pdf" );
    }
    // energy reconstruction bias
    c = a.plotEnergyReconstructionBias( "mean", 0.66, 1.5 );
    if( c )
    {
        print_pdf( iEffAreaFile1, c->GetName(), "EBias.pdf" );
    }
}

/*
 * plot energy treshold plots
 *
 */
void plot_energy_thresholds( string iEffAreaFile, float ze = 20., int NSB = 130 )
{
    VEnergyThreshold b( 0., iEffAreaFile );
    if( b.isZombie() )
    {
        return;
    }
    b.setPlottingYaxis( 100., 500. );
    b.plot_energyThresholds( "E_diffmax", ze, 0.5, NSB, 2.5, 8 );
    b.setPlottingStyle( 46, 24, 1.5, 1.5 );
    b.plot_energyThresholds( "E_sys10p", ze, 0.5, NSB, 2.5, 8, false );
    
    print_pdf( iEffAreaFile, "c_0", "spectralindex.pdf" );
    print_pdf( iEffAreaFile, "c_1", "zenith.pdf" );
    print_pdf( iEffAreaFile, "c_2", "wobble.pdf" );
    print_pdf( iEffAreaFile, "c_3", "noise.pdf" );
    print_pdf( iEffAreaFile, "c_4", "azimuth.pdf" );
}

/*
 * plot everything
 *
 */
void plotIRFReport( string iEffAreaFile1, string iEffAreaFile2 = "", float ze = 20., int NSB = 130 )
{

    plotIRFs( iEffAreaFile1, iEffAreaFile2, ze, NSB );
    
    plot_energy_thresholds( iEffAreaFile1, ze, NSB );
    plot_energy_thresholds( iEffAreaFile2, ze, NSB );
    
}
