/*
anasum_makeSkyPlots.C

Macro to plot significance/excess skymaps given anasum.root file

To be execuded using CINT:
> root -l anasum_makeSkyPlots.C\(\"my_anasum_file.root" \)

To just display the plots but not save them:
> root -l anasum_makeSkyPlots.C\(\"my_anasum_file.root\", false \)

Authors: Henrike Fleischhack, Nathan Kelley-Hoskins

*/

void anasum_makeSkyPlots( TString filename = "anasum.root", bool save_pictures = true )
{

    //load shared library
    gSystem->Load( "$EVNDISPSYS/lib/libVAnaSum.so" );
    
    // Here is where we actually make the skymap
    VPlotAnasumHistograms* vpah = new VPlotAnasumHistograms( filename.Data() ) ;
    
    TCanvas* ca1 = vpah->plot_radec( 0 ) ; //plot type = 0 for significance map
    
    // Lets say we want to change the title on the plot
    // then we need to get a pointer to the plot on the canvas.
    // we can first see the list of all canvases and their contents with:
    // gROOT->GetListOfCanvases()->ls()
    // or we can see the contents of a specific canvas with:
    // ca1->GetListOfPrimitives()->ls()
    TH2D*    pl1 = ( TH2D* )gPad->FindObject( "hmap_stereo_sig_REFLECTED" );
    pl1->SetTitle( "Significance Map" ) ;
    
    // We can also add more things to the plot, like the locations of stars
    // the 'Hipparcos...' is a catalogue kept in $VERITAS_EVNDISP_AUX_DIR/AstroData/Catalogues/
    // You can pick other catalogues from there, but the Hipparchos MAG9 is the preferred
    vpah->plot_catalogue( ca1, "Hipparcos_MAG9_1997.dat" ) ;
    
    // When the 'background' or OFF events are calculated, certain regions are excluded
    // this will plot those excluded regions
    vpah->plot_excludedRegions( ca1 ) ;
    
    
    int w = 550 ;
    int h = 500 ;
    // This will set the size of the plot, in pixels
    // when printed to a file, this is how big the image file will be
    ca1->SetCanvasSize( w, h ) ;
    // This will set the size of the plot's window, for when it displays
    ca1->SetWindowSize( w + ( w - ca1->GetWw() ), h + ( h - ca1->GetWh() ) ) ;
    
    ca1->Modified();
    ca1->Update();
    
    // We can then print the plot to a file. This will put the plot in the same directory as the
    if( save_pictures )
    {
        ca1->Print( ( filename + ".SignificanceMap.png" ).Data() ) ;
    }
    
    TCanvas* ca2 = vpah->plot_radec( 1 ) ; //plot type = 1 for excess map
    
    TH2D*    pl2 = ( TH2D* )gPad->FindObject( "hmap_stereo_diff_REFLECTED" );
    pl2->SetTitle( "Excess Map" ) ;
    
    ca2->SetCanvasSize( w, h ) ;
    ca2->SetWindowSize( w + ( w - ca1->GetWw() ), h + ( h - ca1->GetWh() ) ) ;
    
    ca2->Modified();
    ca2->Update();
    if( save_pictures )
    {
        ca2->Print( ( filename + ".ExcessMap.png" ).Data() ) ;
    }
}
