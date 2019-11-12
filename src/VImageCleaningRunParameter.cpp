/*  \class VImageCleaningRunParameter

     data class for parameters for image cleaning

*/

#include "VImageCleaningRunParameter.h"

VImageCleaningRunParameter::VImageCleaningRunParameter( string iName )
{
    fName  = iName;
    fTelID = 0;
    
    fUseFixedThresholds = false;
    fImageCleaningMethod = 0;
    
    // default two-level threshold cleaning
    fimagethresh = 5.0;
    fborderthresh = 2.5;
    fbrightnonimagetresh = 2.5;
    
    // time cluster cleaning
    ftimecutpixel = 0.5;
    ftimecutcluster = 2.0;
    fminpixelcluster = 3;
    floops = 2;
    
    //cluster cleaning
    fnmaxcluster = 1;
    fminsizecluster = 0;
    
    // Trace Correlation Cleaning
    fCorrelationCleanBoardThresh = 1.0;   // S/N ratio of 1
    fCorrelationCleanCorrelThresh = 0.75; // Sample correlation coefficient of 0.75
    fCorrelationCleanNpixThresh = 15;     // Images whose number of pixels is above this value will skip correlation cleaning
    
    // time two-level cleaning
    ftimediff = 1.0;
    
    // optimized next-neighbour cleaning
    fNNOpt_FakeImageProb = 1.e-4;
    fNNOpt_ActiveNN.assign( 5, false );
    fNNOpt_Multiplicities.push_back( "4NN" );
    fNNOpt_Multiplicities.push_back( "2NNp1" );
    fNNOpt_Multiplicities.push_back( "3NN" );
    fNNOpt_Multiplicities.push_back( "2NN" );
    fNNOpt_Multiplicities.push_back( "BOUND" );
    
    fNNOpt_nRings = 3;
    fNNOpt_CoincWinLimit = 8;
    fNNOpt_ifNNoptNoTimeing = false;
    fNNOpt_ifExplicitSampleTimeSlice = false;
    fNNOpt_sampleTimeSlice = 1;
    fNNOpt_nBinsADC = 25;
    
}

bool VImageCleaningRunParameter::initialize()
{
    return true;
}

void VImageCleaningRunParameter::print()
{
    cout << "\t cleaning method \t \t" << getImageCleaningMethod() << " (" << getImageCleaningMethodIndex();
    cout << ", " << fName << ")" << endl;
    if( getImageCleaningMethod() != "TIMENEXTNEIGHBOUR" )
    {
        cout << "\t\t image/border\t\t\t" << fimagethresh << "/" << fborderthresh;
    }
    else
    {
        cout << "\t\t optimized threshold finding" << endl;
        cout << "\t\t fake image probability\t\t" << fNNOpt_FakeImageProb * 100. << "\%" << endl;
        cout << "\t\t active multiplicities: ";
        for( unsigned int i = 0; i < fNNOpt_Multiplicities.size(); i++ )
        {
            if( fNNOpt_ActiveNN[i] )
            {
                cout << " " << fNNOpt_Multiplicities[i];
            }
        }
        cout << endl;
        cout << "\t\t Maximum number of rings to be searched for boundary pixels - " << fNNOpt_nRings << endl;
        cout << "\t\t Maximal coincidence window allowed between neighbouring pixels for any NN group - " << fNNOpt_CoincWinLimit << " [ns]" << endl;
        if( fNNOpt_ifNNoptNoTimeing )
        {
            cout << "\t\t Ignore timing information in NN image cleaning" << endl;
        }
        else
        {
            cout << "\t\t Use timing information in NN image cleaning" << endl;
        }
        if( fNNOpt_ifExplicitSampleTimeSlice )
        {
            cout << "\t\t set sample time slice and number of ADC bins explicitly:" << endl;
            cout << "\t\t\t\t sampleTimeSlice = " << fNNOpt_sampleTimeSlice << " ns" << endl;
            cout << "\t\t\t\t nBinsADC = " << fNNOpt_nBinsADC << endl;
        }
        
    }
    if( fUseFixedThresholds && getImageCleaningMethod() != "TIMENEXTNEIGHBOUR" )
    {
        cout << "\t\t fixed cleaning thresholds" << endl;
    }
    else if( getImageCleaningMethod() != "TIMENEXTNEIGHBOUR" )
    {
        cout << "\t\t (signal-to-noise cleaning thresholds)" << endl;
    }
    if( getImageCleaningMethodIndex() == 1 )
    {
        cout << "\t\t Tpixel/Tcluster/nMin/nLoops \t" << ftimecutpixel << "/" << ftimecutcluster
             << "/" << fminpixelcluster << "/" << floops << endl;
    }
    if( getImageCleaningMethodIndex() == 3 )
    {
        cout << "\t\t using trace correlation cleaning: " << fCorrelationCleanBoardThresh << "/";
        cout << fCorrelationCleanCorrelThresh << "/" << fCorrelationCleanNpixThresh << "\t";
        cout << "BorderThresh/CorrelationThresh/MaxPixThresh" << endl;
    }
    if( getImageCleaningMethodIndex() == 4 )
    {
        cout << "\t\t time constraint between next neighbor pixels (samples): " << ftimediff << endl;
    }
    if( getImageCleaningMethodIndex() == 5 )
    {
        cout << "\t\t cluster cleaning; " ;
        if( fnmaxcluster > 0 )
        {
            cout << "keep at most " << fnmaxcluster << " clusters; " ;
        }
        if( fminsizecluster > 0 )
        {
            cout << " min cluster size: " << fminsizecluster;
        }
        if( fminpixelcluster > 0 )
        {
            cout << " min no. pixels: " << fminpixelcluster;
        }
        cout << endl;
    }
}


string VImageCleaningRunParameter::getImageCleaningMethod()
{
    if( fImageCleaningMethod == 1 )
    {
        return "TIMECLUSTERCLEANING";
    }
    else if( fImageCleaningMethod == 2 )
    {
        return "TIMENEXTNEIGHBOUR";
    }
    else if( fImageCleaningMethod == 3 )
    {
        return "TWOLEVELANDCORRELATION";
    }
    else if( fImageCleaningMethod == 4 )
    {
        return "TIMETWOLEVEL";
    }
    else if( fImageCleaningMethod == 5 )
    {
        return "CLUSTERCLEANING";
    }
    return "TWOLEVELCLEANING";
}

bool VImageCleaningRunParameter::setImageCleaningMethod( string iMethod )
{
    if( iMethod == "TWOLEVELCLEANING" )
    {
        fImageCleaningMethod = 0;
    }
    else if( iMethod == "TIMECLUSTERCLEANING" )
    {
        fImageCleaningMethod = 1;
    }
    else if( iMethod == "TIMENEXTNEIGHBOUR" )
    {
        fImageCleaningMethod = 2;
    }
    else if( iMethod == "TWOLEVELANDCORRELATION" )
    {
        fImageCleaningMethod = 3;
    }
    else if( iMethod == "TIMETWOLEVEL" )
    {
        fImageCleaningMethod = 4;
    }
    else if( iMethod == "CLUSTERCLEANING" )
    {
        fImageCleaningMethod = 5;
    }
    else
    {
        return false;
    }
    
    return true;
}
