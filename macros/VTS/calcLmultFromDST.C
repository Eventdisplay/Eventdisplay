#include <fstream>

/*

This macro is used to calculate the pixel-wise 'gains' (slopes of the charge vs monitor charge in highlow runs for high gain and low gain mode).
The values calculated here are defined as 'c' here:  https://veritas.sao.arizona.edu/wiki/index.php/ED_sumwindow_dependent_LGcalib#Definitions

The values are written into trees (one entry for each pixel, one file per run & sumwindow). The ratio of the high gain slope to the low gain slope
gives the low gain multiplier for a given sum window.

Usage example:

If you have everything set up as correctly (HiLo runs sorted as pairs and saved to a text file called 'pairs.txt', and the DST files in $VERITAS_USER_DATA_DIR/DST/<sumwindow>), you can just call

.L calcLmultFromDST.C
calcLmultFromDST( <sumwindow> );

Notes:

- The DST files should be produced with a fixed low gain multiplier. Default is 6. Use VLowGainCalibrator::setLowGainMultiplierUsedInDST( <number> ) if you used a different value.
- If you want to have low gain multiplier for different sumwindows in HG and LG, the monitor charge *MUST always be calculated with the same sumwindow.*
- In the example below, set sw2mon = true if you are using sum2 for the monitor charge (v500) and sw2mon = false if you are using sum for the monitor charge (v450-v480).

- DEFAULT FILE PATHS CAN BE ADJUSTED AT THE BOTTOM OF THIS FILE.

*/



/*
this function does all the calculations for one run/sumwindow (one set of DST files).
*/
void calc_one_run( int run, int sw, bool innerHigh, bool sw2monitor, TString dir, TString outdir = "", TString plotdir = "", int max = 1000000 )
{

    //check if run has already been analysed...
    TString name = TString::Format( "%s/Tel_1/%d.lmult.root", outdir.Data(), run );
    TFile* temp = new TFile( name.Data(), "read" );
    if( temp && !temp->IsZombie() )
    {
        cout << "File " << name << " already exists, delete to re-analyse." << endl;
        temp->Close();
        return;
    }
    if( temp )
    {
        temp->Close();
    }
    
    
    //VLowGainCalibrator does all the actual work.
    //make a new object, read in the DST files.
    VLowGainCalibrator* c = new VLowGainCalibrator( run, sw, innerHigh, sw2monitor, dir, outdir );
    
    if( c->fIsOk )
    {
    
        //set some options
        c->setDSTEventLimits( 0, max );				// only analyse events between event nubmer 0 and max.
        c->setAllDebugChannels();				// comment out if you don't want debug output in <run>.debug.root
        c->setFitOptions( 2, 0.9, 0.3, 10, 0.001 ); 		// npoints, puremin, satmax, neventsmin, probmin
        c->setMonitorChargeOptions( 100, -100, true, 5 ) ;	// nmin, summin, usemedian, width
        
        //first loop: calculate monitor charge (average of the 'low' part of the camera) for each event.
        c->makeMonitorChargeHists();
        
        //use the results from the first loop to separate the light levels from another.
        //first, a peak finder is applied to the histogram of monitor charges.
        //then, it is fit with a sum of gaussians (one for each peak) plus linear BG.
        c->findLightLevels( true );
        
        //second loop: calculate the mean & standard deviations of the charge for each pixel & light level.
        c->calculateMeanCharges();
        
        //fit the charge vs monitor charge relation for each pixel, separately for high gain and low gain.
        c->doTheFit();
        
        //light level plots. Only works if debug output is turned on.
        if( plotdir != "" )
        {
            TString plotstring = "QMon:eventNumber:level+1";
            TString cutstring = TString::Format( "channel==%d", ( innerHigh ? 0 : 400 ) );
            
            for( int tel = 0; tel < 4 ; tel++ )
            {
                TString name = TString::Format( "%s/flasher_run%d_t%d.png", plotdir.Data(), run, tel + 1 );
                TCanvas* t = new TCanvas( "t", "t", 1000, 800 );
                c->fDebugtree[tel]->Draw( plotstring.Data(), cutstring.Data(), "colz" );
                t->SaveAs( name.Data() );
                t->Close();
            }
        }
        
        c->terminate();
    }
    delete c;
}


/*
this function reads in a text file with pairs of highlow runs to loop over them.
*/
// "/lustre/fs13/group/cta/users/fleish/DST-EVD-500-sw2_16/" , "/lustre/fs9/group/cta/users/fleish/LMULT-v500-sw2_16"

void calcLmultFromDST( int sw = 6, TString listname = "pairs.txt", TString DSTdir = "$VERITAS_USER_DATA_DIR/Results/DST", TString LMULTdir = "$VERITAS_USER_DATA_DIR/Results/LMULT", TString plotdir = "./plots" )
{

    //change the following variables to reflect your dirctory structure
    ifstream list( listname.Data() );										//list of high/low run pairs, 2 runs per line.
    TString indir = TString::Format( "%s/%d/", DSTdir.Data(), sw );							//directory with DST files, expected <run>.DST.root
    TString outdir = TString::Format( "%s/%d/", LMULTdir.Data(), sw ); 						//directory for output files (<run>.lmult.root)
    bool sw2mon = true;												//true if you are using sw2 for the monitor charge (ED v500)
    
    gSystem->Load( "$EVNDISPSYS/lib/libVAnaSum.so" );
    
    cout << "sumwindow " <<  sw << endl;
    
    int run1, run2;
    while( list >> run1 >> run2 )
    {
        if( run1 <= 0 || run2 <= 0 )
        {
            continue;
        }
        cout << "Highlow runs: " << run1 << " " << run2 << endl;
        calc_one_run( run1, sw, true, sw2mon, indir, outdir, plotdir ) ;
        calc_one_run( run2, sw, false, sw2mon, indir, outdir, plotdir ) ;
        
    }
    exit( 0 );
    
}
