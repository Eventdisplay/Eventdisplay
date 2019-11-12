/*
 * Paranal and LaPalma analysis
 *
 * produces all array layout / scaling / other plots
 *
 * help() provides an overview over the available plotting
 * possibilities
 *
 */

#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/*
 * data class for array/scaling and plotting properties
 *
*/
class CTASensitivityPlotData
{
    public:
    
        string fTitle;
        string fName;
        string fLayout;
        double fEmin;            // in TeV
        double fEmax;            // in TeV
        int    fObsTime;
        double fOffset;
        int    fColor;
        
        CTASensitivityPlotData();
        ~CTASensitivityPlotData() {}
        
};

CTASensitivityPlotData::CTASensitivityPlotData()
{
    fTitle = "";
    fName = "";
    fLayout = "";
    fEmin = 0.01;
    fEmax = 100.;
    fObsTime = 180000;
    fOffset = 0.;
    
    fColor = 1;
}

/*
 * service class which writes out the subarray list files
 * need for the performance plots
 *
 */
class CTASensitivityPlot
{
    public:
    
        string fDataDirectory;
        vector< CTASensitivityPlotData* > fData;
        
        CTASensitivityPlot();
        void addDataSet( string iTitle, string iName, string iArray, double iEmin, double iEmax, int iObsTime, double iOffset, int iColor );
        void setDataDirectory( string iD )
        {
            fDataDirectory = iD;
        }
        bool writeSubArrayList( string iFileName );
};

CTASensitivityPlot::CTASensitivityPlot()
{

}

/*
 * add a single array layout / scaling to the corresponding list
 *
 */
void CTASensitivityPlot::addDataSet( string iTitle, string iName, string iArray,
                                     double iEmin, double iEmax, int iObsTime, double iOffset,
                                     int iColor )
{
    fData.push_back( new CTASensitivityPlotData );
    fData.back()->fTitle = iTitle;
    fData.back()->fName = iName;
    fData.back()->fLayout = iArray;
    fData.back()->fEmin = iEmin;
    fData.back()->fEmax = iEmax;
    fData.back()->fObsTime = iObsTime;
    fData.back()->fOffset = iOffset;
    fData.back()->fColor = iColor;
}

/*
 * write a small ascii file with the list of array layouts / scaling
 * this list is required for the plotting of the performance curves
 *
 */
bool CTASensitivityPlot::writeSubArrayList( string iFileName )
{
    ofstream is;
    is.open( iFileName.c_str() );
    if( !is )
    {
        cout << "Error writing sub array list " << iFileName << endl;
        return false;
    }
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( !fData[i] )
        {
            continue;
        }
        // dummy values -> not correct for the given site
        // (should be changed for performance vs altitude
        //  or geomagnetic filed studies)
        is << "D" << i << " S P 2100. 17.03  -23.13  3.9 ";
        is << fDataDirectory;
        is << fData[i]->fName << "   ";
        is << fData[i]->fLayout << "   ";
        is << fData[i]->fObsTime << "   ";
        is << fData[i]->fOffset << "   ";
        is << fData[i]->fEmin << "   ";
        is << fData[i]->fEmax << "   ";
        is << fData[i]->fColor << "   ";
        is << "21 3001 21  ";
        is << fData[i]->fTitle;
        is << endl;
    }
    is.close();
    
    return true;
}


/*
 * read a list of arrays from an ascii file
 *
 * returns an empty vector if unsuccessful
*/
vector< string > readArrayList( string iArrayList )
{
    vector< string > iA;
    
    // read list of arrays
    ifstream is;
    is.open( iArrayList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, array list file not found: " << endl;
        cout << iArrayList << endl;
        return iA;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 && is_line.substr( 0, 1 ) != "#" )
        {
            iA.push_back( is_line );
        }
    }
    
    is.close();
    
    return iA;
}

/*
 * get a list of reconstrution IDs corresponding to the
 * sub arrays used in the analysis
 *
 * ID0 - full array
 * ID1 - LST array
 * ID2 - MST array
 * ID3 - SST array
 * ID4 - SST + MST
 *
 */
vector< int > getTelTypes( bool iSouth = true )
{
    vector< int > iA;
    
    iA.push_back( 0 );
    iA.push_back( 1 );
    iA.push_back( 2 );
    if( iSouth )
    {
        iA.push_back( 3 );
        iA.push_back( 4 );
    }
    
    return iA;
}

/*
 * get a list of MST types
 *
 */
vector< string > getMSTTypes()
{
    vector< string > iA;
    
    iA.push_back( "N" );
    iA.push_back( "F" );
    
    return iA;
}

/*
 * get a list of SST types
 *
 */
vector< string > getSSTTypes()
{
    vector< string > iA;
    
    iA.push_back( "D" );
    iA.push_back( "G" );
    
    return iA;
}

/*
 * get telescope type title for a given rec ID
 *
 */
string getTTypeTitles( unsigned int iTelType )
{
    if( iTelType == 0 )
    {
        return " ";
    }
    else if( iTelType == 1 )
    {
        return " (LSTs)";
    }
    else if( iTelType == 2 )
    {
        return " (MSTs)";
    }
    else if( iTelType == 3 )
    {
        return " (SSTs)";
    }
    else if( iTelType == 4 )
    {
        return " (MSTs+SSTs)";
    }
    return "";
    
}

/*
 * return a list of possible telescope / camera type combinations
 *
 */
vector< string > getListOfTelTypeCombinations( bool iSouth = true )
{
    vector< string > iA;
    if( iSouth )
    {
        iA.push_back( "FD" );
        iA.push_back( "FG" );
        iA.push_back( "ND" );
        iA.push_back( "NG" );
    }
    else
    {
        iA.push_back( "F" );
    }
    
    return iA;
}

/*
 * get number of available scalings as function of array layout
 *
 * usually: 1-5
 * for Staging arrays: 2-3
 */
vector< int > getNScalings_Layout( string iArrayLayout )
{
    vector< int > iA;
    if( iArrayLayout.find( "HE" ) != string::npos
            || iArrayLayout.find( "HF" ) != string::npos
            || iArrayLayout.find( "HG" ) != string::npos
            || iArrayLayout.find( "HK" ) != string::npos )
    {
        iA.push_back( 2 );
        iA.push_back( 3 );
    }
    else if( iArrayLayout.find( "HB8" ) != string::npos
             || iArrayLayout.find( "HF8" ) != string::npos )
    {
        iA.push_back( 1 );
    }
    else
    {
        for( unsigned int i = 0; i < 5; i++ )
        {
            iA.push_back( i + 1 );
        }
    }
    return iA;
}


/*
 * get number of available scalings as function of array type
 *
 * usually: 1-5
 * for Staging arrays: 2-3
 */
vector< int > getNScalings( string iArrayType )
{
    vector< int > iA;
    if( iArrayType.find( "Staging" ) != string::npos )
    {
        iA.push_back( 2 );
        iA.push_back( 3 );
    }
    else
    {
        for( unsigned int i = 0; i < 5; i++ )
        {
            iA.push_back( i + 1 );
        }
    }
    return iA;
}

/*
 * get a list of array layouts for three different array types
 *
 * Southern side only
 *
 */
vector< string > getArrayLayoutLists( string iArrayType, string iTelTypeCombination )
{
    vector< string > iA;
    /////////////////////////////////////////////////////////////
    // Paranal Layouts
    if( iArrayType == "Baseline Layouts" )
    {
        iA.push_back( "S.3HB1-" + iTelTypeCombination );
        iA.push_back( "S.3HB2-" + iTelTypeCombination );
        iA.push_back( "S.3HB3-" + iTelTypeCombination );
        iA.push_back( "S.3HB4-" + iTelTypeCombination );
        iA.push_back( "S.3HI1-" + iTelTypeCombination );
        if( iTelTypeCombination == "NG" )
        {
            iA.push_back( "S.3HB8-" + iTelTypeCombination );
        }
    }
    else if( iArrayType == "Descoped Layouts" )
    {
        iA.push_back( "S.3HB1-" + iTelTypeCombination );
        iA.push_back( "S.3HD1-" + iTelTypeCombination );
        iA.push_back( "S.3HD2-" + iTelTypeCombination );
    }
    else if( iArrayType == "Extra Descoped Layouts" )
    {
        iA.push_back( "S.3HB1-" + iTelTypeCombination );
        iA.push_back( "S.3HE1-" + iTelTypeCombination );
        iA.push_back( "S.3HE2-" + iTelTypeCombination );
        iA.push_back( "S.3HE3-" + iTelTypeCombination );
    }
    else if( iArrayType == "LST array" )
    {
        iA.push_back( "S.3HB1-" + iTelTypeCombination );
        iA.push_back( "S.3HB2-" + iTelTypeCombination );
        iA.push_back( "S.3HB3-" + iTelTypeCombination );
    }
    else if( iArrayType.find( "Staging" ) != string::npos )
    {
        // always add HB1 for comparision
        iA.push_back( "S.3HB1-" + iTelTypeCombination );
        string iStagingID = iArrayType.substr( iArrayType.size() - 1, 1 );
        std::locale loc;
        // first case: plot all S.3H[E-k]x in one plot
        if( iStagingID.size() > 0 && isdigit( iStagingID[0], loc ) )
        {
            iA.push_back( "S.3HE" + iStagingID + "-" + iTelTypeCombination );
            iA.push_back( "S.3HF" + iStagingID + "-" + iTelTypeCombination );
            iA.push_back( "S.3HG" + iStagingID + "-" + iTelTypeCombination );
            iA.push_back( "S.3HK" + iStagingID + "-" + iTelTypeCombination );
            if( iTelTypeCombination == "NG" )
            {
                iA.push_back( "S.3HF8-" + iTelTypeCombination );
            }
        }
        else
        {
            // e.g. S.HK1, S.3HK2, ...
            for( unsigned int st = 0; st < 7; st++ )
            {
                ostringstream iL;
                iL << "S.3H" << iStagingID << st + 1 << "-" << iTelTypeCombination;
                iA.push_back( iL.str() );
            }
        }
    }
    /////////////////////////////////////////////////////////////
    // LaPalmaLayouts
    else if( iArrayType.size() == 2 && iArrayType.substr( 0, 1 ) == "L" )
    {
        std::ostringstream i_t;
        i_t << "N.3A" << iArrayType << "M15";
        iA.push_back( i_t.str() );
        i_t.str( "" );
        i_t << "N.3B" << iArrayType << "M13";
        iA.push_back( i_t.str() );
        i_t.str( "" );
        i_t << "N.3C" << iArrayType << "M09";
        iA.push_back( i_t.str() );
        i_t.str( "" );
        i_t << "N.3D" << iArrayType << "M05";
        iA.push_back( i_t.str() );
        i_t.str( "" );
        i_t << "N.3Z" << iArrayType << "M15";
        iA.push_back( i_t.str() );
    }
    else if( iArrayType.size() >= 2 &&
             ( iArrayType.substr( 0, 1 ) == "M" ||  iArrayType.substr( 0, 1 ) == "Z" ) )
    {
        std::ostringstream i_t;
        if( iArrayType == "M15" )
        {
            iA.push_back( "N.3AL4M15" );
            iA.push_back( "N.3AL3M15" );
            iA.push_back( "N.3AL2M15" );
        }
        else if( iArrayType == "M13" )
        {
            iA.push_back( "N.3BL4M13" );
            iA.push_back( "N.3BL3M13" );
            iA.push_back( "N.3BL2M13" );
        }
        else if( iArrayType == "M09" )
        {
            iA.push_back( "N.3CL4M09" );
            iA.push_back( "N.3CL3M09" );
            iA.push_back( "N.3CL2M09" );
        }
        else if( iArrayType == "M05" )
        {
            iA.push_back( "N.3DL4M05" );
            iA.push_back( "N.3DL3M05" );
            iA.push_back( "N.3DL2M05" );
        }
        if( iArrayType == "Z15" )
        {
            iA.push_back( "N.3ZL4M15" );
            iA.push_back( "N.3ZL3M15" );
            iA.push_back( "N.3ZL2M15" );
        }
    }
    return iA;
}

/*
 * get list of available observing times
 *
*/
vector< int > getObservingTimes()
{
    vector< int > iA;
    iA.push_back( 180000 ); // 50h
    iA.push_back( 18000 );   // 5h
    iA.push_back( 1800 );   // 30m
    
    return iA;
}

/*
 * get sensitivity requirement ID
 *
*/
int getSensitivityRequirement( int iObsTime, bool iSouth )
{
    if( !iSouth )
    {
        if( iObsTime == 18000 )
        {
            return 6;
        }
        else if( iObsTime == 1800 )
        {
            return 7;
        }
        return 5;
    }
    if( iObsTime == 18000 )
    {
        return 1;
    }
    else if( iObsTime == 1800 )
    {
        return 2;
    }
    
    return 0;
}

/*
 * check if a given site is in the north or in the south
 */
bool isSouth( string iSite )
{
    if( iSite.find( "LaPalma" ) != string::npos )
    {
        return false;
    }
    return true;
}

/*
 * get array layout name pending on North/South
 */
string getArrayLayoutString( bool iSouth, string iArrayLayout, int iScaling )
{
    std::ostringstream i_temp;
    i_temp << iArrayLayout << "-" << iScaling;
    if( !iSouth )
    {
        i_temp << "-F";   // tmptmp assume only F in North
    }
    return i_temp.str();
}


/*
 * get range of sensitivity axis as function of observing time
 */
double getSensitivityAxisUnit( int iObsTime, bool iMax )
{
    // maximum ranges
    if( iMax )
    {
        if( iObsTime == 18000 )
        {
            return 2.5e-10;
        }
        else if( iObsTime == 1800 )
        {
            return 3.5e-10;
        }
        return 1.5E-10;
    }
    // minimum range
    if( iObsTime == 18000 )
    {
        return 2.e-13;
    }
    else if( iObsTime == 1800 )
    {
        return 2.e-12;
    }
    
    return 4.0E-14;
}

/*
 * data directory for a given site / data set
 *
 */

string getDataDir( string iSite )
{
    if( iSite == "prod3-paranalp05-NN" )
    {
        return "S.20deg.20160313/DESY.d20160304.V5.";
    }
    else if( iSite == "prod3-paranalp05-40deg-NN" )
    {
        return "S.40deg.20160313/DESY.d20160304.V5.";
    }
    else if( iSite == "prod3-LaPalmap05-NN" )
    {
        return "N.20deg.20160313/DESY.d20160304.V5.";
    }
    
    return "";
}

/*
 * get a string descriping a given site / data set
 *
 */
string getSiteText( string iSite )
{
    if( iSite == "prod3-paranalp05-NN" )
    {
        return "Paranal 20 deg";
    }
    else if( iSite == "prod3-paranalp05-40deg-NN" )
    {
        return "Paranal 40 deg";
    }
    else if( iSite == "prod3-LaPalmap05-NN" )
    {
        return "La Palma 20 deg";
    }
    
    return "";
}

/*
 * print the last page for the pdf document
 *
 * (empty page)
 *
 */
void print_finalPage( string iPrint )
{
    char hname[200];
    if( iPrint.size() > 0 )
    {
        TCanvas* cEmpty = new TCanvas( "b", "b", 1, 1, 1400, 845 );
        cEmpty->SetBatch();
        cEmpty->Draw();
        sprintf( hname, "%s.pdf)", iPrint.c_str() );
        cEmpty->Print( hname );
    }
}

/*
 * print the first page for the pdf document
 *
 * contains title text and date
 *
 */
void print_titlePage( string iPrint, string iText, float iOffaxisAngle = 0. )
{
    char hname[200];
    if( iPrint.size() > 0 )
    {
        TCanvas* cEmpty = 0;
        cEmpty = new TCanvas( "a", "b", 1, 1, 1400, 845 );
        cEmpty->SetBatch();
        cEmpty->Draw();
        
        TPaveText* iPT = new TPaveText( 0.3, 0.3, 0.7, 0.7, "NB" );
        iPT->SetFillStyle( 0 );
        TText* t = iPT->AddText( iText.c_str() );
        t->SetTextColor( 50 );
        t->SetTextSize( t->GetTextSize() * 2 );
        t = iPT->AddText( "" );
        iPT->AddLine( .0, .5, 1., .5 );
        t->SetTextColor( 50 );
        if( iOffaxisAngle < 1.e-2 )
        {
            sprintf( hname, "on axis" );
        }
        else
        {
            sprintf( hname, "offaxis %.1f deg", iOffaxisAngle );
        }
        t = iPT->AddText( hname );
        t->SetTextColor( 50 );
        t->SetTextFont( 52 );
        t = iPT->AddText( "DESY analysis - 2016-04-18" );
        t->SetTextColor( 50 );
        t->SetTextFont( 52 );
        
        iPT->Draw();
        
        sprintf( hname, "%s.pdf(", iPrint.c_str() );
        cEmpty->Print( hname );
    }
    
}

/*
 * get a list of possible pointing directions
 * (same as used in the analysis and in the naming of
 * the WP Phys files)
 */
vector< string > getPointingDirections()
{
    vector< string > iA;
    
    iA.push_back( "" );
    iA.push_back( "_0deg" );
    iA.push_back( "_180deg" );
    
    return iA;
}

/*
 * return a text string describing the poinint direction
 *
 */
string getPointingString( string iP )
{
    string iA;
    if( iP == "_0deg" )
    {
        return "N";
    }
    else if( iP == "_180deg" )
    {
        return "S";
    }
    return "Av";
}

///////////////////////////////////////////////////////////////////////////////////////
// END OF SERVICE CLASSES AND FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////


/*
 * plot a single array layout / scaling
 *
 */
void plot_arraylayout( string iArrayLayout, int iScaling, int iTelTypeID = 0,
                       int iObsTime = 50 * 3600,
                       string iSite = "prod3-paranalp05-NN", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    CTASensitivityPlot fA;
    fA.setDataDirectory( iDDir );
    
    
    // loop over all pointing directions
    for( unsigned int p = 0; p < iPointing.size(); p++ )
    {
        sprintf( htitle, "%s-%d %s %s", iArrayLayout.c_str(), iScaling, getTTypeTitles( iTelTypeID ).c_str(), getPointingString( iPointing[p] ).c_str() );
        sprintf( hname, "ID%d%sNIM2.%s", iTelTypeID, iPointing[p].c_str(), iSite.c_str() );
        
        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout, iScaling ),
                       0.01, 90., iObsTime, 0.0, p + 1 );
    }
    fA.writeSubArrayList( "subArray.temp.list" );
    
    // title for plots
    sprintf( hname, "%s Array %s, Scaling %d %s %.1fh", iSiteText.c_str(), iArrayLayout.c_str(),
             iScaling, getTTypeTitles( iTelTypeID ).c_str(),
             float( iObsTime ) / 3600. );
    string iTitle = hname;
    // plot everything
    VWPPhysSensitivityPlotsMaker a;
    a.plotAllInOneCanvas();
    a.setPlotRequirements( iSite, true );
    a.compareDataSets( "subArray.temp.list", "", false, 0, iTitle );
    
    TCanvas* c = a.getAllinOneCanvas();
    if( c )
    {
        if( iPrint.size() > 0 )
        {
            sprintf( hname, "%s.pdf", iPrint.c_str() );
            c->Print( hname );
        }
    }
    
    print_finalPage( iPrint );
}

/*
 * HB2 vs HB4 vs HB8
 *
 * paranal only
 *
*/
void paranalHB2_vsHB4_vsHB8( string iSite = "prod3-paranalp05-NN", float iOffaxisAngle = 0., string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    char harray[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s Baseline HB2, HB4, HB8", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText, iOffaxisAngle );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    for( unsigned int p = 0; p < iPointing.size(); p++ )
    {
        // full, MSTs, SSTs
        for( int j = 0; j < getTelTypes().size(); j++ )
        {
            // 5 scalings
            for( int i = 1; i <= 5; i++ )
            {
                // skip LSTs
                if( getTelTypes()[j] == 1 )
                {
                    continue;
                }
                
                CTASensitivityPlot fHB248_NG;
                fHB248_NG.setDataDirectory( iDDir );
                
                sprintf( htitle, "S.3HB2-NG-%d %s", i, getTTypeTitles( getTelTypes()[j] ).c_str() );
                sprintf( hname, "ID%dNIM2.%s", getTelTypes()[j], iSite.c_str() );
                sprintf( harray, "S.3HB2-NG-%d", i );
                
                fHB248_NG.addDataSet( htitle, hname, harray, 0.01, 90., 180000, iOffaxisAngle, 1 );
                
                sprintf( htitle, "S.3HB4-NG-%d %s", i, getTTypeTitles( getTelTypes()[j] ).c_str() );
                sprintf( hname, "ID%dNIM2.%s", getTelTypes()[j], iSite.c_str() );
                sprintf( harray, "S.3HB4-NG-%d", i );
                
                fHB248_NG.addDataSet( htitle, hname, harray, 0.01, 90., 180000, iOffaxisAngle, 2 );
                
                // HB8: no scaling
                sprintf( htitle, "S.3HB8-NG %s", getTTypeTitles( getTelTypes()[j] ).c_str() );
                sprintf( hname, "ID%dNIM2.%s", getTelTypes()[j], iSite.c_str() );
                sprintf( harray, "S.3HB8-NG-1" );
                
                fHB248_NG.addDataSet( htitle, hname, harray, 0.01, 90., 180000, iOffaxisAngle, 4 );
                
                fHB248_NG.writeSubArrayList( "subArray.temp.list" );
                
                // title for plots
                sprintf( hname, "%s Baseline HB2, HB4, HB8, scaling %d, Pointing %s", iSiteText.c_str(), i, getPointingString( iPointing[p] ).c_str() );
                string iTitle = hname;
                
                // plot everything
                VWPPhysSensitivityPlotsMaker a;
                a.plotAllInOneCanvas();
                a.setPlotRequirements( 0, true );
                a.compareDataSets( "subArray.temp.list", "_0deg", false, 0, iTitle );
                
                TCanvas* c = a.getAllinOneCanvas();
                if( c )
                {
                    if( iPrint.size() > 0 )
                    {
                        sprintf( hname, "%s.pdf", iPrint.c_str() );
                        c->Print( hname );
                    }
                }
            }
        }
    }
    
    print_finalPage( iPrint );
    
}

/*
 * plot for each array layout all scalings into one plot
 *
 * - full array, all sub arrays
 * - average pointing, north, south
 *
 *  array list is of long type
*/
void plot_scalings( string iSite = "prod3-paranalp05-NN", string iArrayList = "", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s Array Scaling", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // full, MSTs, SSTs
    for( int j = 0; j < getTelTypes( isSouth( iSite ) ).size(); j++ )
    {
        // loop over all pointing directions
        for( unsigned int p = 0; p < iPointing.size(); p++ )
        {
            // loop over all array layouts
            for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
            {
                CTASensitivityPlot fA;
                fA.setDataDirectory( iDDir );
                
                // array layout dependent number of scalings
                vector< int > iScaling = getNScalings_Layout( iArrayLayout[l] );
                for( unsigned int i = 0; i < iScaling.size(); i++ )
                {
                    sprintf( htitle, "%s-%d %s", iArrayLayout[l].c_str(), iScaling[i], getTTypeTitles( getTelTypes( isSouth( iSite ) )[j] ).c_str() );
                    sprintf( hname, "ID%dNIM2.%s", getTelTypes( isSouth( iSite ) )[j], iSite.c_str() );
                    
                    if( getTelTypes( isSouth( iSite ) )[j] != 1 )
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.01, 90., 180000, 0.0, i + 1 );
                    }
                    // plot LST only up to 900 GeV
                    else
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.01, 0.9, 180000, 0.0, i + 1 );
                    }
                }
                fA.writeSubArrayList( "subArray.temp.list" );
                
                // title for plots
                sprintf( hname, "%s Array %s, Pointing %s %s", iSiteText.c_str(), iArrayLayout[l].c_str(),
                         getPointingString( iPointing[p] ).c_str(),
                         getTTypeTitles( getTelTypes( isSouth( iSite ) )[j] ).c_str() );
                string iTitle = hname;
                
                // plot everything
                VWPPhysSensitivityPlotsMaker a;
                a.plotAllInOneCanvas();
                a.setPlotRequirements( iSite, true );
                if( getTelTypes( isSouth( iSite ) )[j] == 1 )
                {
                    a.setLSTSettings();
                }
                a.compareDataSets( "subArray.temp.list", iPointing[p], false, 0, iTitle );
                
                TCanvas* c = a.getAllinOneCanvas();
                if( iPrint.size() > 0 && c )
                {
                    sprintf( hname, "%s.pdf", iPrint.c_str() );
                    c->Print( hname );
                }
            }
        }
    }
    
    print_finalPage( iPrint );
    
}

/*
 * plot for each array layout the telescope type plot
 *
 * - all arrays
 * - all scalings
 * - average pointing, north, south
 *
 *  array list is of long type
*/
void plot_teltypes( string iSite = "prod3-paranalp05-NN", string iArrayList = "", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s Telescope Types", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // loop over all array layouts
    for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
    {
        vector< int > iScaling = getNScalings_Layout( iArrayLayout[l] );
        // loop over all pointing directions
        for( unsigned int p = 0; p < iPointing.size(); p++ )
        {
            // 5 scalings
            for( unsigned int i = 0; i < iScaling.size(); i++ )
            {
                CTASensitivityPlot fA;
                fA.setDataDirectory( iDDir );
                
                // full, MSTs, SSTs
                for( int j = 0; j < getTelTypes( isSouth( iSite ) ).size(); j++ )
                {
                    sprintf( htitle, "%s-%d %s", iArrayLayout[l].c_str(), iScaling[i], getTTypeTitles( getTelTypes()[j] ).c_str() );
                    sprintf( hname, "ID%dNIM2.%s", getTelTypes()[j], iSite.c_str() );
                    
                    if( getTelTypes()[j] != 1 )
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.01, 90., 180000, 0.0, getTelTypes()[j] + 1 );
                    }
                    // plot LST only up to 900 GeV
                    else
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.01, 0.9, 180000, 0.0, getTelTypes()[j] + 1 );
                    }
                }
                fA.writeSubArrayList( "subArray.temp.list" );
                
                // title for plots
                sprintf( hname, "%s Array %s, Scaling %d, Pointing %s", iSiteText.c_str(), iArrayLayout[l].c_str(), iScaling[i], getPointingString( iPointing[p] ).c_str() );
                string iTitle = hname;
                
                // plot everything
                VWPPhysSensitivityPlotsMaker a;
                a.plotAllInOneCanvas();
                
                a.setPlotRequirements( iSite, true );
                a.compareDataSets( "subArray.temp.list", iPointing[p], false, 0, iTitle );
                
                TCanvas* c = a.getAllinOneCanvas();
                if( iPrint.size() > 0 && c )
                {
                    sprintf( hname, "%s.pdf", iPrint.c_str() );
                    c->Print( hname );
                }
            }
        }
    }
    
    print_finalPage( iPrint );
}

/*
 * plot for each array layout the Av/N/S plots
 *
 * - all arrays
 * - all scalings
 * - all IDs
*/
void plot_pointing( string iSite = "prod3-paranalp05-NN", string iArrayList = "", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s Pointing Direction", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // loop over all array layouts
    for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
    {
        // full, MSTs, SSTs
        for( int j = 0; j < getTelTypes( isSouth( iSite ) ).size(); j++ )
        {
            // array layout dependent number of scalings
            vector< int > iScaling = getNScalings_Layout( iArrayLayout[l] );
            for( unsigned int i = 0; i < iScaling.size(); i++ )
            {
                CTASensitivityPlot fA;
                fA.setDataDirectory( iDDir );
                
                // loop over all pointing directions
                for( unsigned int p = 0; p < iPointing.size(); p++ )
                {
                    sprintf( htitle, "%s-%d %s %s", iArrayLayout[l].c_str(), iScaling[i],
                             getTTypeTitles( getTelTypes()[j] ).c_str(),
                             getPointingString( iPointing[p] ).c_str() );
                    sprintf( hname, "ID%d%sNIM2.%s", getTelTypes()[j], iPointing[p].c_str(), iSite.c_str() );
                    
                    if( getTelTypes()[j] != 1 )
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.011, 90., 180000, 0.0, p + 1 );
                    }
                    // plot LST only up to 900 GeV
                    else
                    {
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling[i] ),
                                       0.01, 0.9, 180000, 0.0, p + 1 );
                    }
                }
                fA.writeSubArrayList( "subArray.temp.list" );
                
                // title for plots
                sprintf( hname, "%s Array %s, Scaling %d %s", iSiteText.c_str(), iArrayLayout[l].c_str(), iScaling[i],
                         getTTypeTitles( getTelTypes()[j] ).c_str() );
                string iTitle = hname;
                
                // plot everything
                VWPPhysSensitivityPlotsMaker a;
                a.plotAllInOneCanvas();
                a.setPlotRequirements( iSite, true );
                a.compareDataSets( "subArray.temp.list", "", false, 0, iTitle );
                
                TCanvas* c = a.getAllinOneCanvas();
                if( iPrint.size() > 0 && c )
                {
                    sprintf( hname, "%s.pdf", iPrint.c_str() );
                    c->Print( hname );
                }
            }
        }
    }
    
    print_finalPage( iPrint );
    
}

/*
 * compare array layouts from a list
 *
 * - one plot per scaling, pointing direction and ID
 *
 *   Baseline Layouts
 *   Descoped Layouts
 *   Extra Descoped Layouts
 *   Staging1 - Staging7
*/
void plot_layout( string iSite = "prod3-paranalp05-NN", string iArrayType = "Baseline Layouts",
                  string iTelTypeCombination = "", string iPrint = "", float iOffaxisAngle = 0. )
{
    char hname[200];
    char htitle[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s %s", iSiteText.c_str(), iArrayType.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText, iOffaxisAngle );
    
    // get lists of telescope type combinations
    vector< string > iTelTypeCombinations = getListOfTelTypeCombinations( isSouth( iSite ) );
    if( iTelTypeCombination.size() != 0 )
    {
        iTelTypeCombinations.clear();
        iTelTypeCombinations.push_back( iTelTypeCombination );
    }
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // array layout scaling
    vector< int > iNScaling = getNScalings( iArrayType );
    // (special treatment for HB8/HF8)
    int iScaling = 1;
    
    // loop over all pointing directions
    for( unsigned int p = 0; p < iPointing.size(); p++ )
    {
        for( unsigned int t = 0; t < iTelTypeCombinations.size(); t++ )
        {
            // # of scalings
            for( unsigned int i = 0; i < iNScaling.size(); i++ )
            {
                // full, MSTs, SSTs
                for( int j = 0; j < getTelTypes( isSouth( iSite ) ).size(); j++ )
                {
                    CTASensitivityPlot fA;
                    fA.setDataDirectory( iDDir );
                    
                    // loop over all array layouts
                    vector< string > iArrayLayout = getArrayLayoutLists( iArrayType, iTelTypeCombinations[t] );
                    for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
                    {
                        iScaling = iNScaling[i];
                        // special treatment for HB8 and HF8 arrays
                        if( iArrayLayout[l].find( "HB8" ) != string::npos || iArrayLayout[l].find( "HF8" ) != string::npos )
                        {
                            iScaling = 1;
                        }
                        sprintf( htitle, "%s-%d %s %s", iArrayLayout[l].c_str(), iScaling,
                                 getTTypeTitles( getTelTypes( isSouth( iSite ) )[j] ).c_str(),
                                 getPointingString( iPointing[p] ).c_str() );
                        sprintf( hname, "ID%dNIM2.%s", getTelTypes( isSouth( iSite ) )[j], iSite.c_str() );
                        
                        fA.addDataSet( htitle, hname, getArrayLayoutString( isSouth( iSite ), iArrayLayout[l], iScaling ),
                                       0.01, 90., 180000, iOffaxisAngle, l + 1 );
                    }
                    fA.writeSubArrayList( "subArray.temp.list" );
                    
                    // title for plots
                    sprintf( hname, "%s %s  Scaling %d, Pointing %s %s", iSiteText.c_str(), iTelTypeCombinations[t].c_str(),
                             iNScaling[i], getPointingString( iPointing[p] ).c_str(),
                             getTTypeTitles( getTelTypes( isSouth( iSite ) )[j] ).c_str() );
                    string iTitle = hname;
                    
                    // plot everything
                    VWPPhysSensitivityPlotsMaker a;
                    a.plotAllInOneCanvas();
                    a.setPlotRequirements( iSite, true );
                    // LST settings
                    if( getTelTypes()[j] == 1 )
                    {
                        a.setLSTSettings();
                    }
                    a.compareDataSets( "subArray.temp.list", "", false, 0, iTitle );
                    
                    TCanvas* c = a.getAllinOneCanvas();
                    if( iPrint.size() > 0 && c )
                    {
                        sprintf( hname, "%s.pdf", iPrint.c_str() );
                        c->Print( hname );
                    }
                }
            }
        }
    }
    
    print_finalPage( iPrint );
}

/*
 * LST layout analysis
 *
 * compare HB1, HB2, HB3 for full array and LST array only
 *
 * - one plot per scaling, pointing direction, teltype combinations
 *
*/
void paranal_LSTs( string iSite = "prod3-paranalp05-NN", bool iShortPlots = false, string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    char harray[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s LST arrays", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get lists array layouts
    vector< string > iTelTypeCombinations = getListOfTelTypeCombinations();
    if( iShortPlots )
    {
        iTelTypeCombinations.clear();
        iTelTypeCombinations.push_back( "NG" );
    }
    // array layout scalings
    vector< int > iScaling;
    for( int i = 1; i <= 5; i++ )
    {
        iScaling.push_back( i );
    }
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // observing times
    vector< int > iObsTimes = getObservingTimes();
    
    // short plots
    if( iShortPlots )
    {
        iTelTypeCombinations.clear();
        iTelTypeCombinations.push_back( "NG" );
        iPointing.clear();
        iPointing.push_back( "" );
        iScaling.clear();
        iScaling.push_back( 3 );
    }
    
    
    // loop over all pointing directions
    for( unsigned int p = 0; p < iPointing.size(); p++ )
    {
        for( unsigned int t = 0; t < iTelTypeCombinations.size(); t++ )
        {
            // 5 scalings
            for( unsigned int i = 0; i < iScaling.size(); i++ )
            {
                // full and LST array only
                for( int j = 0; j < 2; j++ )
                {
                    // observing times
                    for( unsigned int k = 0; k < iObsTimes.size(); k++ )
                    {
                        CTASensitivityPlot fA;
                        fA.setDataDirectory( iDDir );
                        
                        // loop over all array layouts
                        vector< string > iArrayLayout = getArrayLayoutLists( "LST array", iTelTypeCombinations[t] );
                        for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
                        {
                            sprintf( htitle, "%s-%d %s %s", iArrayLayout[l].c_str(), iScaling[i],
                                     getTTypeTitles( j ).c_str(),
                                     getPointingString( iPointing[p] ).c_str() );
                            sprintf( hname, "ID%dNIM2.%s", getTelTypes()[j], iSite.c_str() );
                            sprintf( harray, "%s-%d", iArrayLayout[l].c_str(), iScaling[i] );
                            
                            fA.addDataSet( htitle, hname, harray, 0.01, 90., iObsTimes[k], 0.0, l + 1 );
                        }
                        fA.writeSubArrayList( "subArray.temp.list" );
                        
                        // title for plots
                        sprintf( hname, "%s %s  Scaling %d, Pointing %s %s, %.1fh", iSiteText.c_str(), iTelTypeCombinations[t].c_str(),
                                 iScaling[i], getPointingString( iPointing[p] ).c_str(),
                                 getTTypeTitles( j ).c_str(), ( float )iObsTimes[k] / 3600. );
                        string iTitle = hname;
                        
                        // plot everything
                        VWPPhysSensitivityPlotsMaker a;
                        a.plotAllInOneCanvas( true );
                        a.setPlotRequirements( getSensitivityRequirement( iObsTimes[k], true ), true );
                        a.setEnergyRange_Lin_TeV( 0.01, 0.9 );
                        a.setResolutionLimits( 0.5, 0.4 );
                        a.setSensitivityRatioLimits( 0.4, 1.2 );
                        a.setEffectiveAreaLimits( 5000., 4.e6 );
                        a.setAxisUnits( getSensitivityAxisUnit( iObsTimes[k], false ), getSensitivityAxisUnit( iObsTimes[k], true ) );
                        a.compareDataSets( "subArray.temp.list", "", false, 0, iTitle );
                        
                        TCanvas* c = a.getAllinOneCanvas();
                        if( iPrint.size() > 0 && c )
                        {
                            sprintf( hname, "%s.pdf", iPrint.c_str() );
                            c->Print( hname );
                        }
                    }
                }
            }
        }
    }
    
    print_finalPage( iPrint );
    
}

/*
 * compare MST arrays
 *
 * - one plot per scaling, pointing direction
 * - full arrays and MST only
 *
 * array list is of 'short type: S.3HB1, etc. (no "-NG", "-FG" )
 *
*/
void paranal_MSTs( string iSite = "prod3-paranalp05-NN", string iArrayList = "", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    char harray[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s MST arrays", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // plot for full array and MSTs only
    vector< int > iRecID;
    iRecID.push_back( 0 );
    iRecID.push_back( 2 );
    
    // full array and MSt array
    for( unsigned int r = 0; r < iRecID.size(); r++ )
    {
        // loop over all array layouts
        for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
        {
            // loop over all pointing directions
            for( unsigned int p = 0; p < iPointing.size(); p++ )
            {
                // 5 scalings
                for( int i = 1; i <= 5; i++ )
                {
                    // SST types
                    for( unsigned int s = 0; s < getSSTTypes().size(); s++ )
                    {
                        CTASensitivityPlot fA;
                        fA.setDataDirectory( iDDir );
                        
                        for( unsigned int j = 0; j < getMSTTypes().size(); j++ )
                        {
                            string iArrayLayoutName = iArrayLayout[l] + "-" + getMSTTypes()[j] + getSSTTypes()[s];
                            
                            sprintf( htitle, "%s-%d %s", iArrayLayoutName.c_str(), i, getTTypeTitles( iRecID[r] ).c_str() );
                            sprintf( hname, "ID%dNIM2.%s", iRecID[r], iSite.c_str() );
                            sprintf( harray, "%s-%d", iArrayLayoutName.c_str(), i );
                            
                            fA.addDataSet( htitle, hname, harray, 0.01, 90., 180000, 0.0, j + 1 );
                        }
                        fA.writeSubArrayList( "subArray.temp.list" );
                        
                        // title for plots
                        sprintf( hname, "%s Array %s, Scaling %d, SST %s, Poin. %s %s", iSiteText.c_str(), iArrayLayout[l].c_str(), i,
                                 getSSTTypes()[s].c_str(),
                                 getPointingString( iPointing[p] ).c_str(),
                                 getTTypeTitles( iRecID[r] ).c_str() );
                        string iTitle = hname;
                        
                        // plot everything
                        VWPPhysSensitivityPlotsMaker a;
                        a.plotAllInOneCanvas();
                        a.setPlotRequirements( 0, true );
                        a.compareDataSets( "subArray.temp.list", iPointing[p], false, 0, iTitle );
                        
                        TCanvas* c = a.getAllinOneCanvas();
                        if( iPrint.size() > 0 && c )
                        {
                            sprintf( hname, "%s.pdf", iPrint.c_str() );
                            c->Print( hname );
                        }
                    }
                }
            }
        }
    }
    
    print_finalPage( iPrint );
}

/*
 * compare SST arrays
 *
 * - one plot per scaling, pointing direction
 * - full arrays and MST only
 *
 * array list is of 'short type: S.3HB1, etc.
 *
*/
void paranal_SSTs( string iSite = "prod3-paranalp05-NN", string iArrayList = "", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    char harray[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s SST arrays", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    // plot for full array and MSTs only
    vector< int > iRecID;
    iRecID.push_back( 0 );
    iRecID.push_back( 3 );
    
    // full array and MSt array
    for( unsigned int r = 0; r < iRecID.size(); r++ )
    {
        // loop over all array layouts
        for( unsigned int l = 0; l < iArrayLayout.size(); l++ )
        {
            // loop over all pointing directions
            for( unsigned int p = 0; p < iPointing.size(); p++ )
            {
                // 5 scalings
                for( int i = 1; i <= 5; i++ )
                {
                    // SST types
                    for( unsigned int s = 0; s < getMSTTypes().size(); s++ )
                    {
                        CTASensitivityPlot fA;
                        fA.setDataDirectory( iDDir );
                        
                        for( unsigned int j = 0; j < getSSTTypes().size(); j++ )
                        {
                            string iArrayLayoutName = iArrayLayout[l] + "-" + getMSTTypes()[s] + getSSTTypes()[j];
                            
                            sprintf( htitle, "%s-%d %s", iArrayLayoutName.c_str(), i, getTTypeTitles( iRecID[r] ).c_str() );
                            sprintf( hname, "ID%dNIM2.%s", iRecID[r], iSite.c_str() );
                            sprintf( harray, "%s-%d", iArrayLayoutName.c_str(), i );
                            
                            fA.addDataSet( htitle, hname, harray, 0.01, 90., 180000, 0.0, j + 1 );
                        }
                        fA.writeSubArrayList( "subArray.temp.list" );
                        
                        // title for plots
                        sprintf( hname, "%s Array %s, Scaling %d, MST %s, Poin. %s %s", iSiteText.c_str(), iArrayLayout[l].c_str(), i,
                                 getMSTTypes()[s].c_str(),
                                 getPointingString( iPointing[p] ).c_str(),
                                 getTTypeTitles( iRecID[r] ).c_str() );
                        string iTitle = hname;
                        
                        // plot everything
                        VWPPhysSensitivityPlotsMaker a;
                        a.plotAllInOneCanvas();
                        a.setPlotRequirements( 0, true );
                        a.compareDataSets( "subArray.temp.list", iPointing[p], false, 0, iTitle );
                        
                        TCanvas* c = a.getAllinOneCanvas();
                        if( iPrint.size() > 0 && c )
                        {
                            sprintf( hname, "%s.pdf", iPrint.c_str() );
                            c->Print( hname );
                        }
                    }
                }
            }
        }
    }
    
    print_finalPage( iPrint );
}

/*
 * plot staging scenarios discussed in the CB
 *
 * Sites:
 * prod3-paranalp05-NN
 * prod3-paranalp05-40deg-NN
 * prod3-LaPalmap05-NN
 *
 */
void plot_staging( int iScaling, string iSite = "prod3-paranalp05-NN", string iPrint = "" )
{
    char hname[200];
    char htitle[200];
    char harray[200];
    
    /////////////////////////////////////////
    // data set
    string iDDir = getDataDir( iSite );
    string iSiteText = getSiteText( iSite );
    
    sprintf( hname, "%s Staging Scenarios", iSiteText.c_str() );
    string iTitleText = hname;
    print_titlePage( iPrint, iTitleText );
    
    // get pointing directions
    vector< string > iPointing = getPointingDirections();
    
    vector< string > iStagingName;
    vector< vector< string > > iArrayLayout;
    vector< vector< int > > iTelTypeID;
    
    vector< string > iA_temp;
    vector< int > iID_temp;
    
    ////////////////////////////
    // S1
    iA_temp.clear();
    iID_temp.clear();
    iStagingName.push_back( "S1" );
    if( isSouth( iSite ) )
    {
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE5-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HF5-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HG5-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HK5-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 0 );
    }
    else
    {
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 2 );
        iA_temp.push_back( "N.3DL2M05" );
        iID_temp.push_back( 0 );
        //
        iA_temp.push_back( "N.3DL3M05" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3CL3M09" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3AL3M15" );
        iID_temp.push_back( 1 );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 1 );
        iA_temp.push_back( "N.3DL4M05" );
        iID_temp.push_back( 0 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S2a
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S2a/S4/S7" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE7-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF7-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG7-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK7-ND" );
        iID_temp.push_back( 4 );
    }
    else
    {
        iStagingName.push_back( "S2/S3" );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3DL3M05" );
        iID_temp.push_back( 0 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S2b
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S2b/S6" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HD1-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HD2-ND" );
        iID_temp.push_back( 4 );
    }
    else
    {
        iStagingName.push_back( "S4" );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3CL3M09" );
        iID_temp.push_back( 0 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S3b
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S3a" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE1-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF1-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG1-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK1-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 4 );
    }
    else
    {
        iStagingName.push_back( "S5" );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3AL3M15" );
        iID_temp.push_back( 1 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S3c
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S3b" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE2-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF2-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG2-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK2-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 4 );
    }
    else
    {
        iStagingName.push_back( "S6" );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 1 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S3d
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S3c" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE3-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF3-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG3-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK3-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 4 );
    }
    else
    {
        iStagingName.push_back( "S7" );
        iA_temp.push_back( "N.3AL4M15" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "N.3DL4M05" );
        iID_temp.push_back( 0 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S3e
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S3d" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE4-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF4-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG4-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK4-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 4 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S3e
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S3e" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE6-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF6-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HG6-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HK6-ND" );
        iID_temp.push_back( 4 );
        iA_temp.push_back( "S.3HF8-NG" );
        iID_temp.push_back( 4 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    ////////////////////////////
    // S2a
    iA_temp.clear();
    iID_temp.clear();
    if( isSouth( iSite ) )
    {
        iStagingName.push_back( "S5" );
        iA_temp.push_back( "S.3HB2-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HE7-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HF7-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HG7-ND" );
        iID_temp.push_back( 0 );
        iA_temp.push_back( "S.3HK7-ND" );
        iID_temp.push_back( 0 );
    }
    iArrayLayout.push_back( iA_temp );
    iTelTypeID.push_back( iID_temp );
    
    /////////////////////////////////////////
    // loop over all staging scenarios
    // note: expect same size for iStagingName, iArrayLayout, iTelTypeID vectors
    
    // loop over all pointing directions
    for( unsigned int p = 0; p < iPointing.size(); p++ )
    {
        // loop over all staging scenarios
        for( unsigned int s = 0; s < iStagingName.size(); s++ )
        {
            CTASensitivityPlot fA;
            fA.setDataDirectory( iDDir );
            // layouts of a scenario
            for( unsigned int i = 0; i < iArrayLayout[s].size(); i++ )
            {
            
                sprintf( htitle, "%s-%d %s %s", iArrayLayout[s][i].c_str(), iScaling,
                         getTTypeTitles( iTelTypeID[s][i] ).c_str(), getPointingString( iPointing[p] ).c_str() );
                sprintf( hname, "ID%d%sNIM2.%s", iTelTypeID[s][i], iPointing[p].c_str(), iSite.c_str() );
                if( isSouth( iSite ) )
                {
                    sprintf( harray, "%s-%d", iArrayLayout[s][i].c_str(), iScaling );
                }
                else
                {
                    sprintf( harray, "%s-%d-F", iArrayLayout[s][i].c_str(), iScaling );
                }
                // adjust energy range
                // (LST results not good beyond 10 TeV)
                if( iTelTypeID[s][i] == 1 )
                {
                    fA.addDataSet( htitle, hname, harray, 0.01, 10., 180000, 0.0, i + 1 );
                }
                else
                {
                    if( isSouth( iSite ) )
                    {
                        fA.addDataSet( htitle, hname, harray, 0.01, 90., 180000, 0.0, i + 1 );
                    }
                    else
                    {
                        fA.addDataSet( htitle, hname, harray, 0.01, 50., 180000, 0.0, i + 1 );
                    }
                }
            }
            fA.writeSubArrayList( "subArray.temp.list" );
            
            // title for plots
            sprintf( hname, "%s %s, Scaling %d %s", iSiteText.c_str(), iStagingName[s].c_str(),
                     iScaling, getPointingString( iPointing[p] ).c_str() );
            string iTitle = hname;
            // plot everything
            VWPPhysSensitivityPlotsMaker a;
            a.plotAllInOneCanvas();
            if( isSouth( iSite ) )
            {
                a.setPlotRequirements( 0, true );
            }
            else
            {
                a.setPlotRequirements( 3, true );
            }
            a.setSensitivityRatioLimits( 0.2, 1.05 );
            a.compareDataSets( "subArray.temp.list", "", false, 0, iTitle );
            
            TCanvas* c = a.getAllinOneCanvas();
            if( c )
            {
                if( iPrint.size() > 0 )
                {
                    sprintf( hname, "%s.pdf", iPrint.c_str() );
                    c->Print( hname );
                }
            }
        }
    }
    
    print_finalPage( iPrint );
}


///////////////////////////////////////////////////////////////////////////////////
/*
 * plot all paranal plots
 *
 * this can take a long long time
 *
*/
void plot_paranal()
{
    vector< string > iSite;
    vector< string > iSitePrintName;
    iSite.push_back( "prod3-paranalp05-NN" );
    iSitePrintName.push_back( "paranal-20deg" );
    //    iSite.push_back( "prod3-paranalp05-40deg-NN" );   iSitePrintName.push_back( "paranal-40deg" );
    vector< string > iTelCombinations = getListOfTelTypeCombinations();
    vector< float > iOffAxisBins;
    iOffAxisBins.push_back( 0.0 );
    /*    iOffAxisBins.push_back( 1.5 );
        iOffAxisBins.push_back( 2.5 );
        iOffAxisBins.push_back( 3.5 );  */
    
    for( unsigned int i = 0; i < iSite.size(); i++ )
    {
        ostringstream iListName;
        ostringstream iPrintName;
        ostringstream iStagingName;
        // off-axis bins
        for( unsigned int o = 0; o < iOffAxisBins.size(); o++ )
        {
            ////////////////////
            // layout plots
            //
            // baseline layouts
            for( unsigned int j = 0; j < iTelCombinations.size(); j++ )
            {
                iPrintName.str( "" );
                if( iOffAxisBins[o] < 1.e-2 )
                {
                    iPrintName << iSitePrintName[i] << "-" << iTelCombinations[j] << "-baseline-onAxis";
                }
                else
                {
                    iPrintName << iSitePrintName[i] << "-" << iTelCombinations[j] << "-baseline-" << fixed << setprecision( 1 ) << iOffAxisBins[o];
                }
                plot_layout( iSite[i], "Baseline Layouts", iTelCombinations[j], iPrintName.str(), iOffAxisBins[o] );
                
                // descoped layouts
                iPrintName.str( "" );
                if( iOffAxisBins[o] < 1.e-2 )
                {
                    iPrintName << iSitePrintName[i] << "-" << iTelCombinations[j] << "-descoped-onAxis";
                }
                else
                {
                    iPrintName << iSitePrintName[i] << "-" << iTelCombinations[j] << "-descoped-" << fixed << setprecision( 1 ) << iOffAxisBins[o];
                }
                plot_layout( iSite[i], "Descoped Layouts", iTelCombinations[j], iPrintName.str(), iOffAxisBins[o] );
            }
            continue;
            // staged layouts
            /*            for( unsigned int j = 1; j < 7; j++ )
                        {
                            iStagingName.str( "" );
                            iStagingName << "Staging" << j+1;
                            iPrintName.str( "" );
                            if( iOffAxisBins[o] < 1.e-2 )
                            {
                                iPrintName << iSitePrintName[i] << "-ND-staging" << j+1 << "-onAxis";
                            }
                            else
                            {
                                iPrintName << iSitePrintName[i] << "-ND-staging" << j+1 << "-" << fixed << setprecision( 1 ) << iOffAxisBins[o] << "deg";
                            }
                            plot_layout( iSite[i], iStagingName.str(), "ND", iPrintName.str(), iOffAxisBins[o] );
                        }  */
        }
        continue;
        
        ////////////////////
        // MSTs: compare F vs N
        iPrintName.str( "" );
        iPrintName << iSitePrintName[i] << "-HB1-MSTs";
        //        paranal_MSTs( iSite[i], "subArray.prod3.HB1.list", iPrintName.str() );
        
        
        vector< string > iArrayBaseLayouts;
        iArrayBaseLayouts.push_back( "baseline" );
        iArrayBaseLayouts.push_back( "descoped" );
        iArrayBaseLayouts.push_back( "staging" );
        
        for( unsigned int b = 0; b < iArrayBaseLayouts.size(); b++ )
        {
            iListName.str( "" );
            iListName << "subArray.prod3." << iArrayBaseLayouts[b] << ".longlist";
            
            // teltypes
            /*            iPrintName.str( "" );
                        iPrintName << iSitePrintName[i] << "-telTypes-" << iArrayBaseLayouts[b];
                        paranal_teltypes( iSite[i], iListName.str(), iPrintName.str() ); */
            ////////////////////
            // array layout scaling
            iPrintName.str( "" );
            iPrintName << iSitePrintName[i] << "-scaling-" << iArrayBaseLayouts[b];
            plot_scalings( iSite[i], iListName.str(), iPrintName.str() );
            ////////////////////
            // telescope pointing
            iPrintName.str( "" );
            iPrintName << iSitePrintName[i] << "-pointing-" << iArrayBaseLayouts[b];
            //            paranal_pointing( iSite[i], iListName.str(), iPrintName.str() );
        }
        
        
    }
    
}

///////////////////////////////////////////////////////////////////////////////////
//
void help()
{
    cout << endl;
    cout << "Paranal & LaPalma analysis: plotting of results" << endl;
    cout << endl;
    cout << "Available sites: " << endl;
    cout << "\t prod3-paranalp05-NN" << endl;
    cout << "\t prod3-paranalp05-40deg-NN" << endl;
    cout << "\t prod3-LaPalmap05-NN" << endl;
    cout << endl;
    cout << "---" << endl;
    cout << "Plot a single array layout" << endl << endl;
    cout << "void plot_arraylayout(string iArrayLayout, int iScaling, int iTelTypeID = 0, string iSite = \"prod3-paranalp05-NN\", string iPrint = \"\") " << endl;
    cout << "Example: " << endl;
    cout << "plot_arraylayout( \"S.3HB4-NG\", 5, 0, \"prod3-paranalp05-40deg-NN\" );" << endl;
    cout << "---" << endl;
    cout << "Plot HB2 vs HB4 vs HB8" << endl << endl;
    cout << "void paranalHB2_vsHB4_vsHB8( string iSite = \"prod3-paranalp05-NN\", string iPrint = \"\" )" << endl;
    cout << "---" << endl;
    cout << "Plot for each array layout all scalings into one plot" << endl << endl;
    cout << "void plot_scalings( string iSite = \"prod3-paranalp05-NN\", string iArrayList = \"\", string iPrint = \"\" )" << endl;
    cout << "\t array list is the typical array list used for the CTA analysis" << endl;
    cout << "---" << endl;
    cout << "Plot for each array layout the telescope type plot" << endl << endl;
    cout << "void plot_teltypes( string iSite = \"prod3-paranalp05-NN\", string iArrayList = \"\", string iPrint = \"\" )" << endl;
    cout << "\t array list is the typical array list used for the CTA analysis" << endl;
    cout << "---" << endl;
    cout << "Plot for each array layout the Av/N/S plots" << endl << endl;
    cout << "plot_pointing( string iSite = \"prod3-paranalp05-NN\", string iArrayList = \"\", string iPrint = \"\" )" << endl;
    cout << "\t array list is the typical array list used for the CTA analysis" << endl;
    cout << "---" << endl;
    cout << "Plot array layouts for baseline and different descoped arrays" << endl << endl;
    cout << "void plot_layout( string iSite = \"prod3-paranalp05-NN\", string iArrayType = \"Baseline Layouts\", string iPrint = \"\" )" << endl;
    cout << "\t possible array types (Paranal): \"Baseline Layouts\", \"Descoped Layouts\", \"Extra Descoped Layouts\" " << endl;
    cout << "\t possible array types (LaPalma): \"L4\", \"L3\", \"L2\" (compares different MST arrays)" << endl;
    cout << "\t possible array types (LaPalma): \"M15\", \"M13\", \"M09\", \"M05\" (compares different LST arrays)" << endl;
    cout << "---" << endl;
    cout << "Plot LST plots for 4, 3, 2 LSTs" << endl << endl;
    cout << "void paranal_LSTs( string iSite = \"prod3-paranalp05-NN\", bool iShortPlots = false, string iPrint = \"\" )" << endl;
    cout << "---" << endl;
    cout << "Plot MST array comparision (N vs F)" << endl << endl;
    cout << "void paranal_MSTs( string iSite = \"prod3-paranalp05-NN\", string iArrayList = \"\", string iPrint = \"\" )" << endl;
    cout << "---" << endl;
    
    
    
}
