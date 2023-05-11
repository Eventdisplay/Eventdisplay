/*
 * plot IPR graphs for the different multiplicity (4nn, 3nn, 2nn, 2nn+1)
 * for a set of VERITAS data runs
 *
 * IPR graphs are read from the calibration directory
 *
 * use: plotIPRGraphs(...)
 *
 *
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void help()
{
    cout << "Plot IPR graphs for optimised next-neighbour cleaning for different runs" << endl;
    cout << endl;
    cout << "    reads IPR graph from calibration output" << endl;
    cout << "Usage:  ";
    cout << "plotIPGraph( \"simple_run_list\", \"calibrationDirectory\", telescopeID );" << endl;
    cout << endl;
}

/*
 * read a simple list of runs and return it as a vector<int>
*/
vector< int > getRunList( string iRunList )
{
    vector< int > iR;
    
    ifstream is;
    is.open( iRunList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "Error reading run list from " << iRunList << endl;
        return iR;
    }
    
    string is_line = "";
    while( getline( is, is_line ) )
    {
        iR.push_back( atoi( is_line.c_str() ) );
    }
    is.close();
    
    return iR;
}

/*
 *  function to be called by user
 *
 *  see help()
 *
*/
void plotIPRGraphs( string iRunList, string iCalibrationDirectory, int iTel_ID = 1 )
{
    vector< int > fRunList = getRunList( iRunList );
    if( fRunList.size() == 0 )
    {
        return;
    }
    
    /////////////////////////////////////////////
    vector< string > iNN;
    iNN.push_back( "4nn" );
    iNN.push_back( "3nnrel" );
    iNN.push_back( "2nn" );
    iNN.push_back( "2plus1" );
    iNN.push_back( "IPR" );
    
    for( unsigned int n = 0; n < iNN.size(); n++ )
    {
        string iCName = "c" + iNN[n];
        TCanvas* cNN = new TCanvas( iCName.c_str(), iNN[n].c_str(), 10, 10, 600, 400 );
        cNN->SetLogy( true );
        cNN->Draw();
        
        ostringstream iGraphNameNN;
        if( iNN[n] == "IPR" )
        {
            iGraphNameNN << "IPRchargeTelType" << iTel_ID - 1;
        }
        else
        {
            iGraphNameNN << "graphProbCurve" << iNN[n] << "Type" << iTel_ID - 1;
        }
        unsigned int iCounterNN = 0;
        
        for( unsigned int i = 0; i < fRunList.size(); i++ )
        {
            ostringstream iRootFile;
            iRootFile << iCalibrationDirectory << "/" << fRunList[i] << ".IPRcontours.root";
            
            TFile iF( iRootFile.str().c_str() );
            if( iF.IsZombie() )
            {
                continue;
            }
            
            TGraph* graphProbCurvennType = ( TGraph* )iF.Get( iGraphNameNN.str().c_str() );
            if( graphProbCurvennType )
            {
                graphProbCurvennType->SetTitle( "" );
                // only change color for less than 10 runs
                if( fRunList.size() < 10 )
                {
                    graphProbCurvennType->SetLineColor( iCounterNN + 1 );
                    graphProbCurvennType->SetMarkerColor( iCounterNN + 1 );
                }
                if( iNN[n] == "IPR" )
                {
                    graphProbCurvennType->SetMinimum( 1.e3 );
                }
                else
                {
                    graphProbCurvennType->SetMinimum( 1.e-3 );
                    graphProbCurvennType->SetMaximum( 80. );
                    graphProbCurvennType->GetXaxis()->SetRangeUser( 0., 150. );
                }
                if( iCounterNN == 0 )
                {
                    graphProbCurvennType->Draw();
                }
                else
                {
                    graphProbCurvennType->Draw( "same" );
                }
                if( iNN[n] == "2nn" )
                {
                    cout << "Run " << fRunList[i] << endl;
                    cout << "\t 50 dc value dT: " << graphProbCurvennType->Eval( 50. ) << endl;
                }
                iCounterNN++;
            }
            
        }
        
        // print everything
        ostringstream iPrintName;
        iPrintName << iNN[n] << ".pdf";
        cNN->Print( iPrintName.str().c_str() );
    }
}
