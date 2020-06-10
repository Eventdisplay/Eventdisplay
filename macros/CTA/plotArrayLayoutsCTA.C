/*! \file   plotArrayLayoutsCTA.C
    \brief  plot CTA array layouts

    .L plotArrayLayoutsCTA.C++

*/

#include "TCanvas.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TTree.h"
#include "TText.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

/*
 * data class describing an individual telescope
 *
 */
class VPlotCTAArrayLayout_TelescopeList
{
    public:
    
        int       fTelID;           // prod3b telID
        string    fTelIDName;       // telescope name (plotted as label)
        ULong64_t fTelType;       
        string    fTelTypeName;
        int       fTelID_hyperArray;
        float     fTel_x;
        float     fTel_y;
        float     fTel_z;
        int       fMarkerColor;
        float     fMarkerSize;
        int       fMarkerType;
        
        VPlotCTAArrayLayout_TelescopeList();
        ~VPlotCTAArrayLayout_TelescopeList() {}
};

/*
 * data class describing an array layout
 *
 */
class VPlotCTAArrayLayout_LayoutDescription
{
    public:
    
        string fName;
        string fListFileName;
        string fPrintName;
        
        double fXmax;
        double fYmax;
        
        VPlotCTAArrayLayout_LayoutDescription();
        ~VPlotCTAArrayLayout_LayoutDescription() {}
};


/*
 *
 *  main plotting class
 *
*/

class VPlotCTAArrayLayout
{
    private:
    
        int fTempTelescopeColor;
        
        vector< VPlotCTAArrayLayout_TelescopeList* > fTelescopeList;
        vector< VPlotCTAArrayLayout_TelescopeList* > fTelescopeList_subArray;
        
        void             drawTelescope( VPlotCTAArrayLayout_TelescopeList*, string iText );
        vector< string > getListofArrrays( string iArrayFile );
        void             plotArrayLayoutsFromVectorList( vector< VPlotCTAArrayLayout_LayoutDescription* > iLayout,
                bool iPrintSingleFile );
                
    public:
        VPlotCTAArrayLayout();
        ~VPlotCTAArrayLayout() {}
        
        unsigned int getTelTypeNumber( string iTelType );
        float    getCost_LST( bool iNorth = false )
        {
            if( iNorth ) return 13.;
            return 15.;
        }
        float    getCost_MST( bool iNorth = false )
        {
            if( iNorth ) return 4.4;
            return 3.9;
        }
        float    getCost_SST()
        {
            return 1.3;
        }
        vector< VPlotCTAArrayLayout_TelescopeList* > getListOfTelescopes();
        TCanvas* plot_array( string iArrayName = "", double xmax = 1450., double ymax = 1450.,
                             string iPrintCanvas = "", string drawTelescopeNumbers = "prodID",
                             TCanvas* cUserCanvase = 0 );
        void     plot_allProd3Layouts( string iDataSet, string iLisFileDir,  bool iPrintSingleFile = true );
        void     plot_fromList( string iArrayFile, string iDir, double xmax = 1450., double ymax = 1450.,
                                string drawTelescopeNumbers = "prodID", string iMCProduction = "prod3b" );
        bool     readArrayFromRootFile( string iFile, bool iprod3 = true, double iEasting = 0., double iNorthing = 0. );
        bool     readArrayFromShortTXTFile( string iFile, double iEasting = 0., double iNorthing = 0. );
        bool     readArrayFromTXTFile( string iFile );
        void     printArrayCosts( bool iNorth = false, bool iShort = false );
        void     plotFixedCostMSTSST( float iTotalSum );
        void     plotFixedNorth( float iTotalSum, int iNLSTs_north, int iNMSTs_north );
        void     printListOfTelescopes( int iShort = 0, double dx = 0., double dy = 0. );
        void     printTelescopeDistances( int iTelID, float iDistanceMax = 1.e99 );
        void     plotTelescopeDistances( string iname = "A", double i_max = 300. );
        bool     printTelescopeIDs_for_differentHyperArray( string iFile );
        void     setTelescopeColor( int iTelescopeColor = -99 )
        {
            fTempTelescopeColor = iTelescopeColor;
        }
        bool     setSubArray( string iSubArrayFile = "" );
        void     synchronizeTelescopeLists( vector< VPlotCTAArrayLayout_TelescopeList* > iInputList,
                bool iListFile = false );
};

VPlotCTAArrayLayout::VPlotCTAArrayLayout()
{
    setTelescopeColor();
}

/*
 * print list of available telescope
 *
 * (either full list for subarray list (if filled))
 *
 * iMode = 0: full printout of all values
 * iMode = 1: short list
 * iMode = 3: ecsv type printout (only data, no header)
 * iMode = 4: ecsv file list (CORSIKA coordinates)
 * iMode = 5: list of telescopes in yaml format
 *
 * dx and dy are shifts of the telescope positions
*/
void VPlotCTAArrayLayout::printListOfTelescopes( int iMode, double dx, double dy )
{
    vector< VPlotCTAArrayLayout_TelescopeList* > iPrintList;
    if( fTelescopeList_subArray.size() == 0 )
    {
        iPrintList = fTelescopeList;
    }
    else
    {
        iPrintList = fTelescopeList_subArray;
    }
    if( iMode == 0 )
    {
        for( unsigned int i = 0; i < iPrintList.size(); i++ )
        {
            cout << "Telescope " << i << " (ID: " << iPrintList[i]->fTelID << ")";
            // print hyper-array ID only if different from nominal ID
            if( iPrintList[i]->fTelID_hyperArray != iPrintList[i]->fTelID )
            {
                cout << " (HA " << iPrintList[i]->fTelID_hyperArray << ")";
            }
            cout << "\t" << iPrintList[i]->fTelIDName;
            cout << "\t type " << iPrintList[i]->fTelType << " (" << iPrintList[i]->fTelTypeName << ")";
            cout << "\t" << iPrintList[i]->fTel_x - dx << "\t" << iPrintList[i]->fTel_y + dy;
            cout << "\t" << iPrintList[i]->fTelIDName;
            cout << endl;
        }
    }
    else if( iMode == 1 )
    {
        for( unsigned int i = 0; i < iPrintList.size(); i++ )
        {
            cout << "Telescope " << i + 1;
            cout << "\t" << iPrintList[i]->fTelIDName;
            cout << "\t" << setprecision(2) << fixed << setw(8);
            cout << iPrintList[i]->fTel_x - dx;
            cout << "\t" << iPrintList[i]->fTel_y + dy;
            cout << "\t" << iPrintList[i]->fTel_z << endl;
        }
    }
    else if( iMode == 4 )
    {
        // note! corsika coordinates
        string separator = "    ";
        // header
        cout << "telescope_name" << separator;
        cout << "pos_x" << separator;
        cout << "pos_y" << separator;
        cout << "pos_z" << separator;
        cout << "prod3b_mst_N";
        cout << endl;
        // data
        for( unsigned int i = 0; i < iPrintList.size(); i++ )
        {
            cout << iPrintList[i]->fTelIDName << separator << setw(8);
            cout << iPrintList[i]->fTel_y + dy << separator << setw(8);
            cout << -1.*iPrintList[i]->fTel_x + dx << separator << setw(8);
            cout << iPrintList[i]->fTel_z << separator << setw(8);
            cout << iPrintList[i]->fTelID << separator;
            cout << endl;
        }
    }
    else if( iMode == 5 )
    {
        cout << "[";
        for( unsigned int i = 0; i < iPrintList.size(); i++ )
        {
            cout << iPrintList[i]->fTelIDName;
            if( i != iPrintList.size()-1 ) cout << ", ";
        }
        cout << "]" << endl;
    }
}

/*
 * read number of telescopes of a specific type
 *
 *
*/
unsigned int VPlotCTAArrayLayout::getTelTypeNumber( string iTelType )
{
    unsigned int iT = 0;
    if( fTelescopeList_subArray.size() == 0 )
    {
        fTelescopeList_subArray = fTelescopeList;
    }
    
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        if( fTelescopeList_subArray[i]->fTelTypeName.find( iTelType ) != string::npos )
        {
            iT++;
        }
    }
    return iT;
}

/* 
 * print and plot number of MSTs / SSTS for fixed North scenaries
 *
 */
void VPlotCTAArrayLayout::plotFixedNorth( float iTotalSum, int iNLSTs_north, int iNMSTs_north )
{
    float iTotalSum_North = getCost_LST( true ) * iNLSTs_north + getCost_MST( true ) * iNMSTs_north;
    float iTotalSum_South = iTotalSum - iTotalSum_North;

    cout << "Total sum: " << iTotalSum << endl;
    cout << "Total sum (North): " << iTotalSum_North << endl;
    cout << "Total sum (South): " << iTotalSum_South << endl;

    plotFixedCostMSTSST( iTotalSum_South );



}

/*
 * plot number of MSTs and SSTs for a fixed total cost
 *
 */
void VPlotCTAArrayLayout::plotFixedCostMSTSST( float iTotalSum )
{
    int nmax_MST = (int)(iTotalSum/getCost_MST(false)+0.5);
    int nmax_SST = (int)(iTotalSum/getCost_SST()+0.5);
    cout << endl;
    cout << "Options South:" << endl;
    cout << "\t Maximum number of MSTs: " << nmax_MST << endl;
    cout << "\t Maximum number of SSTs: " << nmax_SST << endl;

    TCanvas *cFixedCosts = new TCanvas( "cFixedCosts", "fixed cost function MSTs/SSTs", 400, 400, 600, 600 );
    cFixedCosts->SetLeftMargin( 0.13 );
    cFixedCosts->SetGridx( 1 );
    cFixedCosts->SetGridy( 1 );

    TH1D *h = new TH1D( "h", "", nmax_SST+1, 0, nmax_MST );
    h->SetMaximum( nmax_SST );
    h->SetMinimum( 0 );
    h->GetXaxis()->SetTitle( "# of MSTs" );
    h->GetYaxis()->SetTitle( "# of SSTs" );
    h->SetStats( 0 );
    h->Draw();

    TF1 *f_fixedCost = new TF1( "f_fixedCost", "([0]-x*[1])/[2]", 0., nmax_MST );
    f_fixedCost->SetParameter( 0, iTotalSum );
    f_fixedCost->SetParameter( 1, getCost_MST( false ) );
    f_fixedCost->SetParameter( 2, getCost_SST() );
    f_fixedCost->Draw( "same" );
    f_fixedCost->GetHistogram()->SetTitle( "" );

    for( int i = 0; i <= nmax_MST; i++ )
    {
        cout << "\t " << i << " MSTs, " << (int)(f_fixedCost->Eval( i )+0.5) << " SSTs" << endl;
    }



}

void VPlotCTAArrayLayout::printArrayCosts( bool iShort, bool iNorth )
{
    unsigned int iLST = getTelTypeNumber( "LST" );
    unsigned int iMST = getTelTypeNumber( "MST" );
    unsigned int iSST = getTelTypeNumber( "DC-SST" ) + getTelTypeNumber( "SC-SST" );
    unsigned int iMSTSCT = getTelTypeNumber( "MSCT" );
    
    if( iShort )
    {
        cout << "Tot: " << iLST + iMST + iSST;
        cout << " #LSTs: " << iLST;
        cout << " #MSTs: " << iMST;
        cout << " #SSTs: " << iSST;
        if( iMSTSCT > 0 )
        {
            cout << " #MSCT: " << iMSTSCT << endl;
        }
        cout << endl;
        return;
    }
    
    // in MEuro
    /*float euro_LST = 6.279;
    float euro_MST = 1.713;
    float euro_SST = 0.954; */
    float euro_LST = getCost_LST( iNorth );
    float euro_MST = getCost_MST( iNorth );
    float euro_SST = getCost_SST();
    
    cout << endl;
    cout << "# of LSTs: " << iLST << " (" << iLST* euro_LST << " MEuro)" <<  endl;
    cout << "# of MSTs: " << iMST << " (" << iMST* euro_MST << " MEuro)" << endl;
    cout << "# of DCSSTs: " << iSST << " (" << iSST* euro_SST << " MEuro)" << endl;
    cout << "===============================================" << endl;
    cout << "Tot: " << iLST* euro_LST + iMST* euro_MST + iSST* euro_SST << " MEuro" << endl;
    cout << endl;
    
}

/*
 * select a subset of the full list of telescopes
 *
 * selection is done using a list of telescopes as typically
 * used in the prodX analysis
 */
bool VPlotCTAArrayLayout::setSubArray( string iSubArrayFile )
{
    fTelescopeList_subArray.clear();
    if( iSubArrayFile.size() == 0 )
    {
        return true;
    }
    // read list of telescope IDs
    ifstream is;
    is.open( iSubArrayFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, subarray file not found: " << endl;
        cout << iSubArrayFile << endl;
        return false;
    }
    
    int iTelID = 0;
    string iTemp = "";
    string iTelName = "";
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            // comment line
            if( is_line.substr( 0, 1 ) == "#" )
            {
                continue;
            }
            istringstream is_stream( is_line );
            is_stream >> iTelID;
            if( iTelID == 0 )
            {
                continue;
            }
            // FOV
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
            }
            // dynamical range
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
            }
            // raw or calibrated
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
            }
            // telescope name
            if( !is_stream.eof() )
            {
                is_stream  >> iTelName;
                if( iTelName == "NOT-FOUND" )
                {
                    stringstream st;
                    st << iTelID;
                    iTelName = st.str();
                }
            }
            else
            {
                stringstream st;
                st << iTelID;
                iTelName = st.str();
            }
            
            // take care!
            // reading a dst file: telID = telID used in simtel
            // reading an evndisp file: telID = evndisp telid;
            // (use hyperarray id here)
            for( unsigned int i = 0; i <  fTelescopeList.size(); i++ )
            {
                if( iTelID == fTelescopeList[i]->fTelID_hyperArray )
                {
                    fTelescopeList_subArray.push_back( fTelescopeList[i] );
                    fTelescopeList_subArray.back()->fTelIDName = iTelName;
                }
            }
        }
    }
    is.close();
    
    cout << "Selected sub array with " << fTelescopeList_subArray.size() << " telescopes (from " << fTelescopeList.size() << ")" << endl;
    
    return true;
}

bool VPlotCTAArrayLayout::readArrayFromShortTXTFile( string iFile, double iEasting, double iNorthing )
{
    fTelescopeList.clear();
    fTelescopeList_subArray.clear();
    
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "Error: ascii file with telescope list not found" << endl;
        return false;
    }
    
    string i_line;
    string itemp;
    unsigned int z = 0;
    while( getline( is, i_line ) )
    {
        // skip empty lines
        if( i_line.size() == 0 )
        {
            continue;
        }
        // skip comment lines
        if( i_line.size() > 0 && i_line.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        
        // new entry
        fTelescopeList.push_back( new VPlotCTAArrayLayout_TelescopeList() );
        istringstream is_stream( i_line );
        
        fTelescopeList.back()->fTelID = z + 1;
        is_stream >> fTelescopeList.back()->fTel_x;
        is_stream >> fTelescopeList.back()->fTel_y;
        fTelescopeList.back()->fTel_x -= iEasting;
        fTelescopeList.back()->fTel_y -= iNorthing;
        if( z < 4 )
        {
            fTelescopeList.back()->fTelTypeName = "23m-LST";
            fTelescopeList.back()->fMarkerColor = 2;
            fTelescopeList.back()->fMarkerSize = 2;
            fTelescopeList.back()->fMarkerType = 24;
        }
        else
        {
            fTelescopeList.back()->fTelTypeName = "12m-MST";
            fTelescopeList.back()->fMarkerColor = 1;
            fTelescopeList.back()->fMarkerSize = 1.5;
            fTelescopeList.back()->fMarkerType = 21;
        }
        
        z++;
    }

    
    is.close();
    
    return true;
}

/*
 * read list of telescope from an ASCII file
 *
 */
bool VPlotCTAArrayLayout::readArrayFromTXTFile( string iFile )
{
    fTelescopeList.clear();
    fTelescopeList_subArray.clear();
    
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "Error: ascii file with telescope list not found" << endl;
        return false;
    }
    
    unsigned int z = 0;
    string i_line;
    string itemp;
    while( getline( is, i_line ) )
    {
        // skip empty lines
        if( i_line.size() == 0 )
        {
            continue;
        }
        // skip comment lines
        if( i_line.size() > 0 && i_line.substr( 0, 1 ) == "#" )
        {
            continue;
        }
        
        // new entry
        
        fTelescopeList.push_back( new VPlotCTAArrayLayout_TelescopeList() );
        
        istringstream is_stream( i_line );
        
        fTelescopeList.back()->fTelID = z + 1;
        fTelescopeList.back()->fTelID_hyperArray = fTelescopeList.back()->fTelID;
        z++;
        // from north is right
        is_stream >> fTelescopeList.back()->fTelIDName;
        is_stream >> fTelescopeList.back()->fTel_y;
        is_stream >> fTelescopeList.back()->fTel_x;
        fTelescopeList.back()->fTel_x *= -1.;
        if( fTelescopeList.back()->fTelIDName.find( "L" ) != string::npos )
        {
            fTelescopeList.back()->fTelTypeName = "23m-LST";
            fTelescopeList.back()->fMarkerColor = 2;
            fTelescopeList.back()->fMarkerSize = 2.0;
            fTelescopeList.back()->fMarkerType = 24;
            continue;
        }
        else if( fTelescopeList.back()->fTelIDName.find( "M" ) != string::npos )
        {
            fTelescopeList.back()->fTelTypeName = "12m-MST";
            fTelescopeList.back()->fMarkerColor = 1;
            fTelescopeList.back()->fMarkerSize = 1.5;
            fTelescopeList.back()->fMarkerType = 21;
            continue;
        }
        else if( fTelescopeList.back()->fTelIDName.find( "S" ) != string::npos )
        {
            fTelescopeList.back()->fTelTypeName = "SST";
            fTelescopeList.back()->fMarkerColor = 4;
            fTelescopeList.back()->fMarkerSize = 1;
        }
    }
    
    is.close();

    fTelescopeList.erase(fTelescopeList.begin());
    
    return true;
}



/*
 * read list of telescopes from root file
 *
*/

bool VPlotCTAArrayLayout::readArrayFromRootFile( string iFile, bool iprod3, double iEasting, double iNorthing )
{
    fTelescopeList.clear();
    fTelescopeList_subArray.clear();
    
    TFile* f1 = new TFile( iFile.c_str() );
    if( f1->IsZombie() )
    {
        return false;
    }
    
    TTree* t = ( TTree* )f1->Get( "telconfig" );
    if( !t )
    {
        return false;
    }
    
    cout << "Telconfig tree found with " << t->GetEntries() << " telescopes" << endl;
    
    float iTelX = 0.;
    float iTelY = 0.;
    float iTelZ = 0.;
    ULong64_t iTelType = 0;
    int iTelID = 0;
    unsigned int iTelIDHA = 0;
    bool iTelIDEqTelIDHA = false;
    
    t->SetBranchAddress( "TelID", &iTelID );
    t->SetBranchAddress( "TelType", &iTelType );
    if( t->GetBranchStatus( "TelID_hyperArray" ) )
    {
        t->SetBranchAddress( "TelID_hyperArray", &iTelIDHA );
    }
    else
    {
        iTelIDEqTelIDHA = true;
    }
    
    t->SetBranchAddress( "TelX", &iTelX );
    t->SetBranchAddress( "TelY", &iTelY );
    t->SetBranchAddress( "TelZ", &iTelZ );
    
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        fTelescopeList.push_back( new VPlotCTAArrayLayout_TelescopeList() );
        
        fTelescopeList.back()->fTelID = iTelID;
        fTelescopeList.back()->fTelType = iTelType;
        fTelescopeList.back()->fTelID_hyperArray = ( int )iTelIDHA;
        if( iTelIDEqTelIDHA )
        {
            fTelescopeList.back()->fTelID_hyperArray = fTelescopeList.back()->fTelID;
        }
        fTelescopeList.back()->fTel_x = iTelX;
        fTelescopeList.back()->fTel_y = iTelY;
        fTelescopeList.back()->fTel_z = iTelZ;
        
        fTelescopeList.back()->fTel_x -= iEasting;
        fTelescopeList.back()->fTel_y -= iNorthing;
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // prod3
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        if( iprod3 )
        {
            // LSTs (1)
            if( iTelType == 138704810 )
            {
                fTelescopeList.back()->fTelTypeName = "23m-LST";
                fTelescopeList.back()->fMarkerColor = 2;
                fTelescopeList.back()->fMarkerSize = 2;
                fTelescopeList.back()->fMarkerType = 24;
            }
            // MST - FlashCam
            else if( iTelType == 10408618 || iTelType == 10608418 )
            {
                fTelescopeList.back()->fTelTypeName = "12m-MST-FlashCam";
                fTelescopeList.back()->fMarkerColor = 1;
                fTelescopeList.back()->fMarkerSize = 1.5;
                fTelescopeList.back()->fMarkerType = 21;
            }
            // MST - NectarCam
            else if( iTelType == 10408418 )
            {
                fTelescopeList.back()->fTelTypeName = "12m-MST-NectarCam";
                fTelescopeList.back()->fMarkerColor = 2;
                fTelescopeList.back()->fMarkerSize = 1.5;
                fTelescopeList.back()->fMarkerType = 21;
            }
            // DC-MST
            else if( iTelType == 207308306 )
            {
                fTelescopeList.back()->fTelTypeName = "9m-SC-MST";
                fTelescopeList.back()->fMarkerColor = 4;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            // DC-SST
            else if( iTelType == 908924 || iTelType == 909924 )
            {
                fTelescopeList.back()->fTelTypeName = "3.5m-DC-SST";
                fTelescopeList.back()->fMarkerColor = 7;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            // SC-SST (Astri)
            else if( iTelType == 201510718 || iTelType == 201511619 )
            {
                fTelescopeList.back()->fTelTypeName = "4.3-Astri-SC-SST";
                fTelescopeList.back()->fMarkerColor = 3;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            // SC-SST (Check)
            else if( iTelType == 201309415 || iTelType == 201309316 )
            {
                fTelescopeList.back()->fTelTypeName = "4.0-GC-SC-SST";
                fTelescopeList.back()->fMarkerColor = 9;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            // SC-MST
            else if( iTelType == 207308707 )
            {
                fTelescopeList.back()->fTelTypeName = "MSCT";
                fTelescopeList.back()->fMarkerColor = 8;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            else
            {
                cout << "unknown telescope type: " << iTelType << endl;
            }
            if( fTempTelescopeColor > 0 && fTelescopeList.size() > 0 )
            {
                fTelescopeList.back()->fMarkerColor = fTempTelescopeColor;
            }
        }
        //////////////////////////////////////////////////////
        // prod2
        //////////////////////////////////////////////////////
        else
        {
            // LSTs (1)
            if( iTelType == 138704810 || iTelType == 141305009 || iTelType == 141305109 )
            {
                fTelescopeList.back()->fTelTypeName = "23m-LST";
                fTelescopeList.back()->fMarkerColor = 2;
                fTelescopeList.back()->fMarkerSize = 2;
                fTelescopeList.back()->fMarkerType = 24;
            }
            // standard MSTs (2)
            else if( iTelType == 10007818 || iTelType == 10408418 || iTelType == 10008118 )
            {
                fTelescopeList.back()->fTelTypeName = "12m-MST";
                fTelescopeList.back()->fMarkerColor = 1;
                fTelescopeList.back()->fMarkerSize = 1.5;
            }
            // large pixel MSTs (4)
            else if( iTelType == 10009725 )
            {
                fTelescopeList.back()->fTelTypeName = "12m-MST-LPix";
                fTelescopeList.back()->fMarkerColor = 6;
                fTelescopeList.back()->fMarkerSize = 1.5;
            }
            // standard SSTs (3)
            else if( iTelType == 3709725 || iTelType == 3709425 || iTelType == 3710125 )
            {
                fTelescopeList.back()->fTelTypeName = "7m-DC-SST";
                fTelescopeList.back()->fMarkerColor = 3;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            // 7m telescopes (5, prod1) or SCT (prod2)
            else if( iTelType == 7309930  || iTelType == 201509515 )
            {
                fTelescopeList.back()->fTelTypeName = "4m-SC-SST";
                fTelescopeList.back()->fMarkerColor = 4;
                fTelescopeList.back()->fMarkerSize = 1;
            }
            else
            {
                cout << "unknown telescope type: " << iTelType << endl;
            }
        }
    }
    
    return true;
}

bool VPlotCTAArrayLayout::printTelescopeIDs_for_differentHyperArray( string iFile )
{
    TFile f1( iFile.c_str() );
    if( f1.IsZombie() )
    {
        return false;
    }
    
    TTree* t = ( TTree* )f1.Get( "telconfig" );
    if( !t )
    {
        return false;
    }
    
    cout << "Telconfig tree found with " << t->GetEntries() << " telescopes" << endl;
    
    float iTelX = 0.;
    float iTelY = 0.;
    ULong64_t iTelType = 0;
    int iTelID = 0;
    unsigned int iTelIDHA = 0;
    
    t->SetBranchAddress( "TelID", &iTelID );
    t->SetBranchAddress( "TelType", &iTelType );
    if( t->GetBranchStatus( "TelID_hyperArray" ) )
    {
        t->SetBranchAddress( "TelID_hyperArray", &iTelIDHA );
    }
    t->SetBranchAddress( "TelX", &iTelX );
    t->SetBranchAddress( "TelY", &iTelY );
    
    unsigned int z_tel = 0;
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        bool bFound = false;
        if( !fTelescopeList_subArray[i] )
        {
            continue;
        }
        
        for( unsigned int j = 0; j < t->GetEntries(); j++ )
        {
            t->GetEntry( j );
            
            if( TMath::Abs( iTelX - fTelescopeList_subArray[i]->fTel_x ) < 1.e-2
                    && TMath::Abs( iTelY - fTelescopeList_subArray[i]->fTel_y ) < 1.e-2
                    && iTelType == fTelescopeList_subArray[i]->fTelType )
            {
                bFound = true;
                z_tel++;
            }
        }
        if( !bFound )
        {
            cout << "TELESCOPE NOT FOUND: ID " << iTelID << ", type " << iTelType << " at [";
            cout << fTelescopeList_subArray[i]->fTel_x << "," << fTelescopeList_subArray[i]->fTel_y << "]" <<  endl;
        }
    }
    f1.Close();
    
    cout << "found " << z_tel << " telescopes (should be " << fTelescopeList_subArray.size() << ")" << endl;
    
    return true;
}




/////////////////////////////////////////////////////////////////////////////////////


/*

   plot a single telescope

*/
void VPlotCTAArrayLayout::drawTelescope( VPlotCTAArrayLayout_TelescopeList* iD, string iText )
{
    if( !iD )
    {
        return;
    }
    
    TGraph* iG = new TGraph( 1 );
    iG->SetMarkerStyle( iD->fMarkerType );
    iG->SetMarkerSize( iD->fMarkerSize );
    iG->SetMarkerColor( iD->fMarkerColor );
    
    iG->SetPoint( 0, iD->fTel_x, iD->fTel_y );
    
    iG->Draw( "p" );
    
    if( iText.size() > 0 )
    {
        TText* iT = new TText( iD->fTel_x, iD->fTel_y, iText.c_str() );
        iT->SetTextSize( 0.014 );
        iT->SetTextColor( 14 );
        iT->Draw();
    }
}


/*
     plot array layout

*/
TCanvas* VPlotCTAArrayLayout::plot_array( string iname, double xmax, double ymax, string iPrintCanvas,
        const string drawTelescopeNumbers, TCanvas* cUserCanvas )
{
    char hname[200];
    char htitle[200];
    if( iname.size() > 0 )
    {
        sprintf( hname, "cArrayLayout_%s", iname.c_str() );
        sprintf( htitle, "array layout (%s)", iname.c_str() );
    }
    else
    {
        sprintf( hname, "cArrayLayout" );
        sprintf( htitle, "array layout" );
    }
    
    TCanvas* c = 0;
    if( !cUserCanvas )
    {
        c = new TCanvas( hname, htitle, 10, 10, 800, 800 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetRightMargin( 0.13 );
        c->SetLeftMargin( 0.13 );
        c->SetTopMargin( 0.13 );
        c->SetBottomMargin( 0.13 );
        c->Draw();
        
        sprintf( hname, "hnull_%s", iname.c_str() );
        TH2D* hnull = new TH2D( hname, "", 100, -1.*xmax, xmax, 100, -1.*ymax, ymax );
        hnull->SetStats( 0 );
        hnull->SetYTitle( "North [m]" );
        hnull->SetXTitle( "East [m]" );
        hnull->GetYaxis()->SetTitleOffset( 1.6 );
        hnull->Draw();
        
        int i_max = ( int )xmax / 200.;
        for( int i = 0; i < i_max; i++ )
        {
            TEllipse* i_c200 = new TEllipse( 0., 0., ( i + 1 ) * 200. );
            i_c200->SetLineStyle( 2 );
            i_c200->SetLineColor( kGray );
            i_c200->SetFillStyle( 4000 );
            i_c200->Draw();
        }
        
    }
    
    // check that sub array is set
    if( fTelescopeList_subArray.size() == 0 )
    {
        fTelescopeList_subArray = fTelescopeList;
    }
    
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        if( !fTelescopeList_subArray[i] )
        {
            continue;
        }
        
        if( drawTelescopeNumbers == "prodID" )
        {
            sprintf( hname, "%d", fTelescopeList_subArray[i]->fTelID_hyperArray );
            drawTelescope( fTelescopeList_subArray[i], hname );
        }
        else if( drawTelescopeNumbers == "name" )
        {
            drawTelescope( fTelescopeList_subArray[i], fTelescopeList_subArray[i]->fTelIDName );
        }
        else
        {
            drawTelescope( fTelescopeList_subArray[i], "" );
        }
    }
    // draw array name
    if( iname.size() > 0 )
    {
        TText* it = new TText( -1.*0.9 * xmax, 0.75 * ymax, iname.c_str() );
        it->SetTextSize( it->GetTextSize() * 1.2 );
        it->Draw();
    }
    // draw number of telescopes
    if( getTelTypeNumber( "MSCT" ) > 0 )
    {
        sprintf( hname, "#LSTs: %u #MSTs: %u #SSTs: %u #MSCTs: %u ",
                 getTelTypeNumber( "LST" ),
                 getTelTypeNumber( "MST" ),
                 getTelTypeNumber( "SST" ),
                 getTelTypeNumber( "MSCT" ) );
                 
    }
    else if( getTelTypeNumber( "SST" ) > 0 )
    {
        sprintf( hname, "#LSTs: %u #MSTs: %u #SSTs: %u",
                 getTelTypeNumber( "LST" ),
                 getTelTypeNumber( "MST" ),
                 getTelTypeNumber( "SST" ) );
    }
    else
    {
        sprintf( hname, "#LSTs: %u #MSTs: %u",
                 getTelTypeNumber( "LST" ),
                 getTelTypeNumber( "MST" ) );
    }
    TText* iTT = new TText( -1.*0.9 * xmax, -0.9 * ymax, hname );
    iTT->SetTextSize( 0.020 );
    iTT->SetTextColor( 1 );
    iTT->SetTextFont( 42 );
    iTT->Draw();
    
    if( iPrintCanvas.size() > 0 )
    {
        sprintf( hname, "%s.%s", iname.c_str(), iPrintCanvas.c_str() );
        c->Print( hname );
    }
    
    return c;
}
void VPlotCTAArrayLayout::printTelescopeDistances( int iTelID, float iDistanceMax )
{
    // check that sub array is set
    if( fTelescopeList_subArray.size() == 0 )
    {
        fTelescopeList_subArray = fTelescopeList;
    }
    
    // get telescope
    unsigned int iTelID_sub = 99999;
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        if( fTelescopeList_subArray[i]->fTelID == iTelID )
        {
            iTelID_sub = i;
        }
    }
    if( iTelID_sub >= fTelescopeList_subArray.size() )
    {
        cout << "telescope ID not found" << endl;
        return;
    }
    
    // fill list
    multimap< float, unsigned int > iSortedListOfTelescopes;
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        if( i != iTelID_sub )
        {
            float d = sqrt( ( fTelescopeList_subArray[i]->fTel_x - fTelescopeList_subArray[iTelID_sub]->fTel_x )
                            * ( fTelescopeList_subArray[i]->fTel_x - fTelescopeList_subArray[iTelID_sub]->fTel_x )
                            + ( fTelescopeList_subArray[i]->fTel_y - fTelescopeList_subArray[iTelID_sub]->fTel_y )
                            * ( fTelescopeList_subArray[i]->fTel_y - fTelescopeList_subArray[iTelID_sub]->fTel_y ) );
            iSortedListOfTelescopes.insert( pair< float, unsigned int >( d, i ) );
        }
    }
    
    // print list
    for( multimap<float, unsigned int>::iterator it = iSortedListOfTelescopes.begin(); it != iSortedListOfTelescopes.end(); ++it )
    {
        if( it->second < fTelescopeList_subArray.size() && it->first < iDistanceMax )
        {
            cout << "..to telescope " << fTelescopeList_subArray[it->second]->fTelIDName;
            cout << " (" << fTelescopeList_subArray[it->second]->fTelTypeName << "): " << it->first << "m" << endl;
        }
    }
    
}

/*

    plot separation of different telescope types

*/
void VPlotCTAArrayLayout::plotTelescopeDistances( string iname, double i_max )
{
    // check that sub array is set
    if( fTelescopeList_subArray.size() == 0 )
    {
        fTelescopeList_subArray = fTelescopeList;
    }
    
    // map with distance histograms
    map< ULong64_t, TH1F* > i_DistanceHistograms;
    
    char hname[200];
    int z = 0;
    
    // loop over all telescopes
    unsigned int iTelID_sub = 99999;
    for( unsigned int i = 0; i < fTelescopeList_subArray.size(); i++ )
    {
        // new histogram
        if( i_DistanceHistograms.find( fTelescopeList_subArray[i]->fTelType ) == i_DistanceHistograms.end() )
        {
            sprintf( hname, "h_%d", z );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType] = new TH1F( hname, "", ( int )i_max, 0., i_max );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->SetLineColor( z + 1 );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->SetStats( 0 );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->SetXTitle( "distance to nearest neighbour (of same type) [m]" );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->SetFillStyle( 3004 );
            i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->SetFillColor( z + 1 );
            z++;
        }
        
        // calculate closest distance to telescope of same type
        float i_min = 1.e5;
        for( unsigned int j = 0; j < fTelescopeList_subArray.size(); j++ )
        {
            // require same telescope type
            if( i == j || fTelescopeList_subArray[i]->fTelType != fTelescopeList_subArray[j]->fTelType )
            {
                continue;
            }
            float i_dist = sqrt( ( fTelescopeList_subArray[i]->fTel_x - fTelescopeList_subArray[j]->fTel_x )
                                 * ( fTelescopeList_subArray[i]->fTel_x - fTelescopeList_subArray[j]->fTel_x )
                                 + ( fTelescopeList_subArray[i]->fTel_y - fTelescopeList_subArray[j]->fTel_y )
                                 * ( fTelescopeList_subArray[i]->fTel_y - fTelescopeList_subArray[j]->fTel_y ) );
            if( i_dist < i_min )
            {
                i_min = i_dist;
            }
        }
        cout << "TELESCOPE " << fTelescopeList_subArray[i]->fTelID << " (type " << fTelescopeList_subArray[i]->fTelTypeName << "): ";
        cout << i_min << endl;
        i_DistanceHistograms[fTelescopeList_subArray[i]->fTelType]->Fill( i_min );
    }
    // plot
    char htitle[200];
    sprintf( hname, "cDist_%s", iname.c_str() );
    sprintf( htitle, "telescope distances (%s)", iname.c_str() );
    TCanvas* cDist = new TCanvas( hname, htitle, 1000, 100, 600, 600 );
    cDist->SetGridx( 0 );
    cDist->SetGridy( 0 );
    cDist->Draw();
    
    z = 0;
    for( map< ULong64_t, TH1F* >::iterator it = i_DistanceHistograms.begin(); it != i_DistanceHistograms.end(); ++it )
    {
        if( it->second )
        {
            if( z == 0 )
            {
                it->second->Draw();
                z++;
            }
            else
            {
                it->second->Draw( "same" );
            }
            cout << "Average distance " << it->second->GetMean() << " +- " << it->second->GetRMS() << endl;
        }
    }
}

/*
 * return list of telescopes 
 *
 */
vector< VPlotCTAArrayLayout_TelescopeList* > VPlotCTAArrayLayout::getListOfTelescopes()
{
    if( fTelescopeList_subArray.size() == 0 )
    {
        return fTelescopeList;
    }
    return fTelescopeList_subArray;
}

/*
 * uses the input telescope list to select for the current array:
 * - the telescopes in the input list (by comparing positions)
 * - prints the telescope list
 * - prints error if there is a mismatch
 *
 * assumption is that both array definitions use the same coordinate definition
*/
void VPlotCTAArrayLayout::synchronizeTelescopeLists( vector< VPlotCTAArrayLayout_TelescopeList* > iInputList,
        bool iListFile )
{
    cout << "Reading input file list of size " << iInputList.size() << endl;

    fTelescopeList_subArray.clear();

    for( unsigned int i = 0; i < iInputList.size(); i++ )
    {
        // check if there is a corresponding telescope in
        // the telescope list
        for( unsigned int t = 0; t < fTelescopeList.size(); t++ )
        {
            if( iInputList[i] && fTelescopeList[t] )
            {
                double d = sqrt( (iInputList[i]->fTel_x-fTelescopeList[t]->fTel_x)
                        *(iInputList[i]->fTel_x-fTelescopeList[t]->fTel_x)
                        +(iInputList[i]->fTel_y-fTelescopeList[t]->fTel_y)
                        *(iInputList[i]->fTel_y-fTelescopeList[t]->fTel_y) );
                // allow for tolerances of 1 m
                if( d < 1. )
                {
                    // printout fitting into a list file
                    if( iListFile )
                    {
                        cout << fTelescopeList[t]->fTelID << "    ";
                        cout << "      20.     12    1     ";
                        cout << iInputList[i]->fTelIDName << endl;
                    }
                    // verbose printout
                    else
                    {
                        cout << "Overlapping telescope positions (d < " << d << ")" << endl;
                        cout << "\t input " << iInputList[i]->fTelID << " (" << iInputList[i]->fTelIDName << ")" << endl;
                        cout << "\t match " << fTelescopeList[t]->fTelID << " (" << fTelescopeList[t]->fTelIDName << ")" << endl;
                    }
                        
                }
            }
        }
    }

}

////////////////////////////////////////////////

vector< string > VPlotCTAArrayLayout::getListofArrrays( string iArrayFile )
{
    vector< string > SubArray;
    ifstream is;
    is.open( iArrayFile.c_str(), ifstream::in );
    if( !is )
    {
        return SubArray;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        SubArray.push_back( is_line );
    }
    
    return SubArray;
}

/*
 *  plot array layouts read from a run list file
 *
 */
void VPlotCTAArrayLayout::plot_fromList( string iArrayFile, const string iDir,
        double xmax, double ymax,
        const string drawTelescopeNumbers,
        const string iMCProduction )
{
    vector< string > iSubArray = getListofArrrays( iArrayFile );
    if( iSubArray.size() == 0 )
    {
        return;
    }
    
    vector< VPlotCTAArrayLayout_LayoutDescription* > fLayout;
    
    for( unsigned int i = 0; i < iSubArray.size(); i++ )
    {
        string iName = iDir;
        if( iMCProduction == "prod2" )
        {
            iName += "/CTA.prod2";
            iName += iSubArray[i] + ".lis";
        }
        else if( iMCProduction == "prod3" )
        {
            iName += "/CTA.prod3";
            iName += iSubArray[i] + ".lis";
        }
        else
        {
            iName += "/CTA.prod3Sb";
            if( iSubArray[i].size() > 0 )
            {
                iName += iSubArray[i].substr( 1 ) + ".lis";
            }
        }
        cout << iName << endl;
        
        fLayout.push_back( new VPlotCTAArrayLayout_LayoutDescription() );
        fLayout.back()->fListFileName = iName;
        fLayout.back()->fName = iSubArray[i];
        fLayout.back()->fXmax = xmax;
        fLayout.back()->fYmax = ymax;
    }
    
    TCanvas* c1 = new TCanvas( "c1", "", 10, 10, 800, 800 );
    c1->Draw();
    
    TPaveText* iPT = new TPaveText( 0.1, 0.3, 0.9, 0.7, "NB" );
    iPT->SetFillStyle( 0 );
    string iTT = iMCProduction + " layouts - 2016-11-10";
    TText* t = iPT->AddText( iTT.c_str() );
    t->SetTextColor( 1 );
    t->SetTextFont( 52 );
    /*    t = iPT->AddText( "South - prod3b" );
        t->SetTextColor( 1 );
        t->SetTextFont( 52 ); */
    
    iPT->Draw();
    
    c1->Print( "Prod3.pdf(" );
    
    plotArrayLayoutsFromVectorList( fLayout, true );
    
    TCanvas* c2 = new TCanvas( "c2", "", 10, 10, 800, 800 );
    c2->Draw();
    c2->Print( "Prod3.pdf)" );
    
}

/*
 * plot all prod3 layouts (South or North)
 *
 * Paranal-prod3b-TS
 *
 */
void VPlotCTAArrayLayout::plot_allProd3Layouts( string iDataSet, const string iLisFileDir, bool iPrintSingleFile )
{
    // MSTs
    vector< string > ivMST;
    ivMST.push_back( "N" );
    if( iDataSet.find( "South-HB9" ) == string::npos )
    {
        ivMST.push_back( "F" );
    }
    // SSTs
    vector< string > ivSST;
    if( iDataSet.find( "Paranal" ) != string::npos || iDataSet.find( "South-HB9" ) != string::npos )
    {
        ivSST.push_back( "G" );
        if( iDataSet.find( "South-HB9" ) != string::npos )
        {
            ivSST.push_back( "D" );
            ivSST.push_back( "A" );
        }
    }
    // paranal prod3b threshold
    if( iDataSet.find( "Paranal-prod3b-TS" ) != string::npos )
    {
        ivSST.push_back( "D" );
        ivSST.push_back( "A" );
    }
    
    // North
    else
    {
        /*        ivSST.push_back( "1" );
                ivSST.push_back( "2" );
                ivSST.push_back( "3" );
                ivSST.push_back( "4" ); */
        ivSST.push_back( "5" );
    }
    
    // layouts
    vector< string > iA;
    
    // South: hexagonal array layouts
    if( iDataSet == "Paranal-baseline" )
    {
        iA.push_back( "3HB1" );
        iA.push_back( "3HB2" );
        iA.push_back( "3HB3" );
        iA.push_back( "3HB4" );
        iA.push_back( "3HI1" );
        iA.push_back( "3HI2" );
    }
    else if( iDataSet == "Paranal-baselineHB8" )
    {
        iA.push_back( "3HB8" );
    }
    else if( iDataSet == "Paranal-staged" )
    {
        iA.push_back( "3HD1" );
        iA.push_back( "3HD2" );
        iA.push_back( "3HD3" );
        iA.push_back( "3HE1" );
        iA.push_back( "3HE2" );
        iA.push_back( "3HE3" );
        iA.push_back( "3HE4" );
        iA.push_back( "3HE5" );
        iA.push_back( "3HE6" );
        iA.push_back( "3HE7" );
        iA.push_back( "3HF1" );
        iA.push_back( "3HF2" );
        iA.push_back( "3HF3" );
        iA.push_back( "3HF4" );
        iA.push_back( "3HF5" );
        iA.push_back( "3HF6" );
        iA.push_back( "3HF7" );
        iA.push_back( "3HG1" );
        iA.push_back( "3HG2" );
        iA.push_back( "3HG3" );
        iA.push_back( "3HG4" );
        iA.push_back( "3HG5" );
        iA.push_back( "3HG6" );
        iA.push_back( "3HG7" );
        iA.push_back( "3HK1" );
        iA.push_back( "3HK2" );
        iA.push_back( "3HK3" );
        iA.push_back( "3HK4" );
        iA.push_back( "3HK5" );
        iA.push_back( "3HK6" );
        iA.push_back( "3HK7" );
        iA.push_back( "3HM1" );
    }
    else if( iDataSet == "Paranal-stagedHB8" )
    {
        iA.push_back( "3HF8" );
    }
    else if( iDataSet == "Paranal-staged04HB8" )
    {
        iA.push_back( "3HB8-LST3-MST24-SST74" );
        iA.push_back( "3HB8-LST3-MST18-SST74" );
        iA.push_back( "3HB8-LST3-MST15-SST74" );
        iA.push_back( "3HB8-LST3-MST12-SST74" );
        iA.push_back( "3HB8-LST3-MST06-SST74" );
        iA.push_back( "3HB8-LST3-MST18-SST74" );
        iA.push_back( "3HB8-LST3-MST18-SST60" );
        iA.push_back( "3HB8-LST3-MST18-SST50" );
        iA.push_back( "3HB8-LST3-MST18-SST41" );
    }
    else if( iDataSet == "Paranal-staged04HB1" )
    {
        iA.push_back( "3HB1-LST4-MST24-SST70" );
        iA.push_back( "3HB8-LST2-MST24-SST70" );
        iA.push_back( "3HB1-LST4-MST18-SST70" );
        iA.push_back( "3HB1-LST4-MST15-SST70" );
        iA.push_back( "3HB1-LST4-MST12-SST70" );
        iA.push_back( "3HB1-LST4-MST06-SST70" );
        iA.push_back( "3HB1-LST4-MST18-SST70" );
        iA.push_back( "3HB1-LST4-MST18-SST60" );
        iA.push_back( "3HB1-LST4-MST18-SST50" );
        iA.push_back( "3HB1-LST4-MST18-SST40" );
    }
    else if( iDataSet == "Paranal-stagedHB89" )
    {
        iA.push_back( "3HB89" );               // Reference 1
        iA.push_back( "3HB89-R3-A" );    // Reference 3
        iA.push_back( "3HB89-R3-B" );    // Reference 3
        iA.push_back( "3HB89-R4-A" );    // Reference 4
        iA.push_back( "3HB89-R4-B" );    // Reference 4
        iA.push_back( "3HB89-R5-A" );    // Reference 5
        iA.push_back( "3HB89-S1-A" );    // staging 1
        iA.push_back( "3HB89-S1-B" );    // staging 1
        iA.push_back( "3HB89-S2-A" );    // staging 2
        iA.push_back( "3HB89-S3-A" );    // staging 3
        iA.push_back( "3HB89-S3-B" );    // staging 3
        iA.push_back( "3HB89-S4-A" );    // staging 4
        iA.push_back( "3HB89-S4-B" );    // staging 4
        iA.push_back( "3HB89-S5-A" );    // staging 5
        iA.push_back( "3HB89-S5-B" );    // staging 5
        iA.push_back( "3HB89-S6-A" );    // staging 6
        iA.push_back( "3HB89-S7-A" );    // staging 7
        iA.push_back( "3HB89-S7-B" );    // staging 7
        iA.push_back( "3HB89-S8-A" );    // staging 8
        iA.push_back( "3HB89-S8-B" );    // staging 8
        iA.push_back( "3HB89-S8-C" );    // staging 8
    }
    else if( iDataSet == "Paranal-prod3b-TS" )
    {
        iA.push_back( "3HB9" );
        iA.push_back( "3HB9-TS-AA" );
        iA.push_back( "3HB9-TS-AB" );
        iA.push_back( "3HB9-TS-BA" );
        iA.push_back( "3HB9-TS-BB" );
        iA.push_back( "3HB9-TS-CA" );
        iA.push_back( "3HB9-TS-CB" );
    }
    else if( iDataSet == "Paranal-thresholdHB89" )
    {
        iA.push_back( "3HB89" );
        iA.push_back( "3HB89-TS-AA" );
        iA.push_back( "3HB89-TS-BA" );
        iA.push_back( "3HB89-TS-CA" );
        iA.push_back( "3HB89-TS-AB" );
        iA.push_back( "3HB89-TS-BB" );
        iA.push_back( "3HB89-TS-CB" );
        iA.push_back( "3HB89-TS-DB" );
    }
    else if( iDataSet == "Paranal-staged05HB1" )
    {
        iA.push_back( "3HB1-L4M18S70" );
        iA.push_back( "3HB1-L4M12S50" );
        iA.push_back( "3HB1-L4M12S60" );
        iA.push_back( "3HB1-L4M15S40" );
        iA.push_back( "3HB1-L4M15S50" );
        iA.push_back( "3HB1-L4M15S60" );
        iA.push_back( "3HB1-L4M18S30" );
        iA.push_back( "3HB1-L4M18S40" );
        iA.push_back( "3HB1-L4M18S50" );
        iA.push_back( "3HB1-L2M10S40" );
        iA.push_back( "3HB1-L2M10S50" );
        iA.push_back( "3HB1-L2M10S60" );
        iA.push_back( "3HB1-L2M12S30" );
        iA.push_back( "3HB1-L2M12S40" );
        iA.push_back( "3HB1-L2M12S50" );
        iA.push_back( "3HB1-L2M15S30" );
        iA.push_back( "3HB1-L2M15S40" );
        iA.push_back( "3HB1-L3M10S30" );
        iA.push_back( "3HB1-L3M10S40" );
        iA.push_back( "3HB1-L3M10S50" );
        iA.push_back( "3HB1-L3M12S30" );
        iA.push_back( "3HB1-L4M24S70" );
    }
    else if( iDataSet == "South-HB9" )
    {
        iA.push_back( "3HB9" );
        iA.push_back( "3HB89" );
        iA.push_back( "3HB8" );
        iA.push_back( "3HB4-2" );
        iA.push_back( "3HB1-2" );
        iA.push_back( "3HB2-2" );
    }
    else if( iDataSet == "LaPalma-thresholdDL4" )
    {
        iA.push_back( "3AL4M15" );
        iA.push_back( "3DL4M05" );
        iA.push_back( "3EL4M05" );
        iA.push_back( "3FL4M05" );
        iA.push_back( "3GL4M05" );
    }
    // La Palma (North)
    else if( iDataSet == "LaPalma" )
    {
        iA.push_back( "3AL4M15" );
        iA.push_back( "3AL3M15" );
        iA.push_back( "3AL2M15" );
        
        iA.push_back( "3BL4M13" );
        iA.push_back( "3BL3M13" );
        iA.push_back( "3BL2M13" );
        
        iA.push_back( "3CL4M09" );
        iA.push_back( "3CL3M09" );
        iA.push_back( "3CL2M09" );
        
        iA.push_back( "3DL4M05" );
        iA.push_back( "3DL3M05" );
        iA.push_back( "3DL2M05" );
        
        iA.push_back( "3ZL4M19" );
        iA.push_back( "3ZL3M19" );
        iA.push_back( "3ZL2M19" );
        
        iA.push_back( "3YL4M15" );
        
        iA.push_back( "3XL4M00" );
    }
    
    string iFil;
    string iLeg;
    vector< VPlotCTAArrayLayout_LayoutDescription* > fLayout;
    
    for( unsigned int a = 0; a < iA.size(); a++ )
    {
        for( unsigned int m = 0; m < ivMST.size(); m++ )
        {
            for( unsigned int s = 0; s < ivSST.size(); s++ )
            {
            
                if( iDataSet.find( "Paranal" ) != string::npos && iDataSet.find( "prod3b" ) == string::npos )
                {
                    iFil = iLisFileDir + "/CTA.prod3S." + iA[a] + "-" + ivMST[m] + ivSST[s] + ".lis";
                }
                else if( iDataSet == "South-HB9" || iDataSet == "Paranal-prod3b-TS" )
                {
                    iFil = iLisFileDir + "/CTA.prod3Sb." + iA[a] + "-" + ivMST[m] + ivSST[s] + ".lis";
                }
                else
                {
                    iFil = iLisFileDir + "/CTA.prod3N." + iA[a] + "-" + ivSST[s] + "-" + ivMST[m] + ".lis";
                }
                fLayout.push_back( new VPlotCTAArrayLayout_LayoutDescription() );
                fLayout.back()->fListFileName = iFil;
                
                if( iDataSet.find( "Paranal" ) != string::npos || iDataSet.find( "South-HB9" ) != string::npos )
                {
                    iLeg = iA[a] + "-" + ivMST[m] + ivSST[s];
                    fLayout.back()->fName = iLeg;
                    fLayout.back()->fXmax = 1450.;
                    fLayout.back()->fYmax = 1450.;
                    if( iDataSet.find( "Paranal" ) != string::npos )
                    {
                        fLayout.back()->fPrintName = "3S." + iA[a] + "-" + ivMST[m] + ivSST[s] + ".pdf";
                    }
                    else
                    {
                        fLayout.back()->fPrintName = iFil = "3Sb." + iA[a] + "-" + ivMST[m] + ivSST[s] + ".pdf";
                    }
                }
                else
                {
                    iLeg = iA[a] + "-" + ivSST[s] + "-" + ivMST[m];
                    fLayout.back()->fName = iLeg;
                    fLayout.back()->fXmax = 1450.;
                    fLayout.back()->fYmax = 1450.;
                    fLayout.back()->fPrintName = "3S." + iA[a] + "-" + ivSST[s] + "-" + ivMST[m] + ".pdf";
                }
            }
        }
    }
    
    plotArrayLayoutsFromVectorList( fLayout, iPrintSingleFile );
}

void VPlotCTAArrayLayout::plotArrayLayoutsFromVectorList( vector< VPlotCTAArrayLayout_LayoutDescription* > iLayout,
        bool iPrintSingleFile )
{
    TCanvas* c1 = new TCanvas( "c1", "", 10, 10, 800, 800 );
    c1->Draw();
    
    for( unsigned int i = 0; i < iLayout.size(); i++ )
    {
        if( !iLayout[i] )
        {
            continue;
        }
        if( !setSubArray( iLayout[i]->fListFileName ) )
        {
            continue;
        }
        TCanvas* c = plot_array( iLayout[i]->fName, iLayout[i]->fXmax, iLayout[i]->fYmax );
        cout << iLayout[i]->fName << ": ";
        printArrayCosts( false, true );
        if( c )
        {
            if( iPrintSingleFile )
            {
                c->Print( "Prod3.pdf" );
            }
            else
            {
                c->Print( iLayout[i]->fPrintName.c_str() );
            }
        }
    }
}



////////////////////////////////////////////////////////////

VPlotCTAArrayLayout_TelescopeList::VPlotCTAArrayLayout_TelescopeList()
{
    fTelID = 0;
    fTelType = 0;
    fTelID_hyperArray = 0;
    fTel_x = 0.;
    fTel_y = 0.;
    fTelTypeName = "";
    
    fMarkerColor = 1;
    fMarkerSize = 1;
    fMarkerType = 8;
}

////////////////////////////////////////////////////////////

VPlotCTAArrayLayout_LayoutDescription::VPlotCTAArrayLayout_LayoutDescription()
{
    fName = "";
    fListFileName = "";
    fXmax = 1450.;
    fYmax = 1450.;
    fPrintName = "";
}



/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// plot a complete book of all array layouts
void plot_allProd3LayoutsSouth( string iArrayLayoutFile = "$CTA_USER_DATA_DIR/prod3/telconfig-HD-3.root",
                                string iArrayLayoutFileHB8 = "$CTA_USER_DATA_DIR/prod3/HB89-NG-FD.root" )
{
    string iListFileName = gSystem->ExpandPathName( "$CTA_EVNDISP_AUX_DIR/DetectorGeometry/" );
    //    string iListFileName = gSystem->ExpandPathName( "/Users/maierg/Experiments/CTA/Analysis/prod3/prod3_South_HB9/" );
    
    TCanvas* c1 = new TCanvas( "c1", "", 10, 10, 800, 800 );
    c1->Draw();
    
    TPaveText* iPT = new TPaveText( 0.1, 0.3, 0.9, 0.7, "NB" );
    iPT->SetFillStyle( 0 );
    TText* t = iPT->AddText( "Prod3b HB9 Layouts - 2016-09-12" );
    t->SetTextColor( 1 );
    t->SetTextFont( 52 );
    t = iPT->AddText( "South - prod3b" );
    t->SetTextColor( 1 );
    t->SetTextFont( 52 );
    
    iPT->Draw();
    
    c1->Print( "Prod3.pdf(" );
    
    ///////////////////////////////////////////////////////////////////////////
    // baseline layouts
    VPlotCTAArrayLayout* iBaseline = new VPlotCTAArrayLayout();
    /*    iBaseline->readArrayFromRootFile( iArrayLayoutFile );
        iBaseline->plot_allProd3Layouts( "Paranal-baseline", iListFileName, true );
    
        iBaseline->readArrayFromRootFile( iArrayLayoutFileHB8 );
        iBaseline->plot_allProd3Layouts( "Paranal-baselineHB8", iListFileName, true );
    
        iBaseline->readArrayFromRootFile( iArrayLayoutFile );
        iBaseline->plot_allProd3Layouts( "Paranal-staged", iListFileName, true );
    
        iBaseline->readArrayFromRootFile( iArrayLayoutFileHB8 );
        iBaseline->plot_allProd3Layouts( "Paranal-stagedHB8", iListFileName, true ); */
    
    /*    iBaseline->readArrayFromRootFile( iArrayLayoutFileHB8 );
        iBaseline->plot_allProd3Layouts( "Paranal-staged04HB8", iListFileName, true );
    
        iBaseline->readArrayFromRootFile( iArrayLayoutFile );
        iBaseline->plot_allProd3Layouts( "Paranal-staged04HB1", iListFileName, true ); */
    
    //    iBaseline->readArrayFromRootFile( iArrayLayoutFile );
    //    iBaseline->plot_allProd3Layouts( "Paranal-staged05HB1", iListFileName, true );
    
    //    iBaseline->readArrayFromRootFile( iArrayLayoutFileHB8 );
    //    iBaseline->plot_allProd3Layouts( "Paranal-thresholdHB89", iListFileName, true );
    
    //    iBaseline->readArrayFromRootFile( iArrayLayoutFile );
    //    iBaseline->plot_allProd3Layouts( "South-HB9", iListFileName, true );
    
    iBaseline->readArrayFromRootFile( iArrayLayoutFile );
    iBaseline->plot_allProd3Layouts( "Paranal-prod3b-TS", iListFileName, true );
    
    TCanvas* c2 = new TCanvas( "c2", "", 10, 10, 800, 800 );
    c2->Draw();
    c2->Print( "Prod3.pdf)" );
    
}

/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// plot a complete book of all array layouts
void plot_allProd3LayoutsNorth( string iArrayLayoutFile )
{
    string iListFileName = gSystem->ExpandPathName( "$CTA_EVNDISP_AUX_DIR/DetectorGeometry/" );
    
    TCanvas* c1 = new TCanvas( "c1", "", 10, 10, 800, 800 );
    c1->Draw();
    
    TPaveText* iPT = new TPaveText( 0.1, 0.3, 0.9, 0.7, "NB" );
    iPT->SetFillStyle( 0 );
    TText* t = iPT->AddText( "Prod3 Array Layouts - 2016-05-25" );
    t->SetTextColor( 1 );
    t->SetTextFont( 52 );
    t = iPT->AddText( "North" );
    t->SetTextColor( 1 );
    t->SetTextFont( 52 );
    
    iPT->Draw();
    
    c1->Print( "Prod3.pdf(" );
    
    ///////////////////////////////////////////////////////////////////////////
    // baseline layouts
    VPlotCTAArrayLayout* iBaseline = new VPlotCTAArrayLayout();
    iBaseline->readArrayFromRootFile( iArrayLayoutFile );
    //    iBaseline->plot_allProd3Layouts( "LaPalma-thresholdDL4", iListFileName, true );
    iBaseline->plot_allProd3Layouts( "LaPalma", iListFileName, true );
    
    TCanvas* c2 = new TCanvas( "c2", "", 10, 10, 800, 800 );
    c2->Draw();
    c2->Print( "Prod3.pdf)" );
    
}

