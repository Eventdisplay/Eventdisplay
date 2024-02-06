/*
 * plot some important showerpars variables.
 *
 * mostly used for testing of productions
 *
 */

#include <algorithm>
#include <sstream>
#include <string>

void help()
{
    cout << "plot showerpars variables - compare different subsystems" << endl;
    cout  << "plot( \"my*.root\", false/true for South/North)" << endl;
    return;
}

void plot_showerpars_variables( string iFile, bool iNorth = false )
{
    unsigned int fNSubSystems = 4;
    if( iNorth )
    {
        fNSubSystems = 3;
    }
    
    TChain* fD = new TChain( "showerpars" );
    fD->Add( iFile.c_str() );
    
    string fCut = "NImages[REPL]>1&&Chi2[REPL]>=0.";
    
    vector< string > fVar;
    fVar.push_back( "NImages[REPL]" );
    fVar.push_back( "sqrt( (MCxoff-Xoff[REPL])*(MCxoff-Xoff[REPL]) + (MCyoff-Yoff[REPL])*(MCyoff-Yoff[REPL]) )" );
    fVar.push_back( "log10(MCe0)" );
    fVar.push_back( "sqrt( Xcore[REPL]* Xcore[REPL] + Ycore[REPL]*Ycore[REPL])" );
    
    TCanvas* c = new TCanvas( "cshowerpars", "", 10, 10, 800, 800 );
    c->Divide( 2, 2 );
    c->Draw();
    
    for( unsigned int f = 0; f < fVar.size(); f++ )
    {
        TPad* pad = (TPad*)c->cd( f + 1 );
        gPad->SetLogy(1);
        gPad->SetGridx(0);
        for( unsigned int i = 0; i < fNSubSystems; i++ )
        {
            fD->SetLineColor( i + 1 );
            
            string iVar = fVar[f];
            std::ostringstream ss;
            ss << i;
            while( iVar.find( "REPL" ) != string::npos )
            {
                iVar.replace( iVar.find( "REPL" ), 4, ss.str() );
            }
            string iCut = fCut;
            while( iCut.find( "REPL" ) != string::npos )
            {
                iCut.replace( iCut.find( "REPL" ), 4, ss.str() );
            }
            
            if( i == 0 )
            {
                fD->Draw( iVar.c_str(), iCut.c_str() );
            }
            else
            {
                fD->Draw( iVar.c_str(), iCut.c_str(), "same" );
            }
            c->Draw();
        }
    }
    c->Print("showerpars.pdf");
}




