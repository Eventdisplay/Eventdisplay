/* \class VZDCF
   \brief (plotting) functions for Z-transformed Discrete Correlation Function algorithm (ZDCF)

   use fortran code available at

   http://www.weizmann.ac.il/weizsites/tal/research/software/

   to calculate ZDCF

   This class reads in the ZDCF results file (suffix .dcf) and allows plotting

*/

#include "VZDCF.h"

VZDCF::VZDCF()
{

    setMLinterval();
    
}

bool VZDCF::readZDCF( string iFile )
{

    fZDCFData.clear();
    
    // read in ascii file
    ifstream is( iFile.c_str() );
    if( !is )
    {
        cout << "VZDCF::readZDCF error reading ZDCF file: " << iFile << endl;
        return false;
    }
    string is_line;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() == 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        fZDCFData.push_back( new VZDCFData() );
        
        //! no errors are catched here..
        is_stream >> fZDCFData.back()->tau;
        is_stream >> fZDCFData.back()->sigma_tau_neg;
        is_stream >> fZDCFData.back()->sigma_tau_pos;
        is_stream >> fZDCFData.back()->dcf;
        is_stream >> fZDCFData.back()->dcf_error_low;
        is_stream >> fZDCFData.back()->dcf_error_up;
        is_stream >> fZDCFData.back()->nbins;
    }
    is.close();
    
    cout << "total number of ZDCF data points: " << fZDCFData.size() << endl;
    
    return true;
}

bool VZDCF::print()
{
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        fZDCFData[i]->print();
    }
    
    return true;
}

double VZDCF::getZDCFData_tau_min( bool bError )
{
    double i_m = 1.e19;
    
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        if( fZDCFData[i] )
        {
            if( bError && fZDCFData[i]->tau - fZDCFData[i]->sigma_tau_neg < i_m )
            {
                i_m = fZDCFData[i]->tau - fZDCFData[i]->sigma_tau_neg;
            }
            else if( fZDCFData[i]->tau < i_m )
            {
                i_m = fZDCFData[i]->tau;
            }
        }
    }
    return i_m;
}

double VZDCF::getZDCFData_tau_max( bool bError )
{
    double i_m = -1.e19;
    
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        if( fZDCFData[i] )
        {
            if( bError && fZDCFData[i]->tau + fZDCFData[i]->sigma_tau_pos > i_m )
            {
                i_m = fZDCFData[i]->tau + fZDCFData[i]->sigma_tau_pos;
            }
            else if( fZDCFData[i]->tau > i_m )
            {
                i_m = fZDCFData[i]->tau;
            }
        }
    }
    return i_m;
}

double VZDCF::getZDCFData_dcf_min( bool bError, double iTauMin, double iTauMax )
{
    double i_m = 1.e19;
    
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        if( fZDCFData[i] )
        {
            if( iTauMin > 0. && fZDCFData[i]->tau < iTauMin )
            {
                continue;
            }
            if( iTauMax > 0. && fZDCFData[i]->tau > iTauMax )
            {
                continue;
            }
            if( bError && fZDCFData[i]->dcf - fZDCFData[i]->dcf_error_low < i_m )
            {
                i_m = fZDCFData[i]->dcf - fZDCFData[i]->dcf_error_low;
            }
            else if( fZDCFData[i]->dcf < i_m )
            {
                i_m = fZDCFData[i]->dcf;
            }
        }
    }
    return i_m;
}

double VZDCF::getZDCFData_dcf_max( bool bError, double iTauMin, double iTauMax )
{
    double i_m = -1.e19;
    
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        if( fZDCFData[i] )
        {
            if( iTauMin > 0. && fZDCFData[i]->tau < iTauMin )
            {
                continue;
            }
            if( iTauMax > 0. && fZDCFData[i]->tau > iTauMax )
            {
                continue;
            }
            if( bError && fZDCFData[i]->dcf + fZDCFData[i]->dcf_error_up > i_m )
            {
                i_m = fZDCFData[i]->dcf + fZDCFData[i]->dcf_error_up;
            }
            else if( fZDCFData[i]->dcf > i_m )
            {
                i_m = fZDCFData[i]->dcf;
            }
        }
    }
    return i_m;
}



TCanvas* VZDCF::plot( TCanvas* c, bool bzdcf, double taumin, double taumax, double ymin, double ymax )
{
    char hname[800];
    char htitle[800];
    // empty histogram for axis
    TH1D* hZDCF = 0;
    
    // canvas
    if( !c )
    {
        if( bzdcf )
        {
            sprintf( hname, "cZDCF" );
            sprintf( htitle, "ZDCF" );
        }
        else
        {
            sprintf( hname, "cZDCF_sig" );
            sprintf( htitle, "ZDCF (dcf/error)" );
        }
        
        c = new TCanvas( hname, htitle, 10, 10, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        
        if( bzdcf )
        {
            sprintf( hname, "hZDCF" );
        }
        else
        {
            sprintf( hname, "hZDCF_sig" );
        }
        
        // histogram values
        if( taumin < -9990. )
        {
            taumin = getZDCFData_tau_min( true ) - 5.;
        }
        if( taumax < -9990. )
        {
            taumax = getZDCFData_tau_max( true ) + 5.;
        }
        
        hZDCF = new TH1D( hname, "", 100, taumin, taumax );
        hZDCF->SetStats( 0 );
        hZDCF->SetXTitle( "time delay [days]" );
        hZDCF->GetXaxis()->CenterTitle( true );
        if( bzdcf )
        {
            hZDCF->SetYTitle( "ZDCF" );
        }
        else
        {
            hZDCF->SetYTitle( "ZDCF / error " );
        }
        if( ymin > -99. && ymax > -99. )
        {
            hZDCF->SetMinimum( ymin );
            hZDCF->SetMaximum( ymax );
        }
        else if( bzdcf )
        {
            hZDCF->SetMinimum( getZDCFData_dcf_min( true, taumin, taumax ) * 1.1 );
            hZDCF->SetMaximum( getZDCFData_dcf_max( true, taumin, taumax ) * 1.1 );
        }
        else
        {
            hZDCF->SetMinimum( -5. );
            hZDCF->SetMaximum( ymax );
        }
        hZDCF->Draw( "" );
        hZDCF->Draw( "AH" );
        
        plot_nullHistogram( c, hZDCF, false, true, 1.2, taumin, taumax );
        
    }
    else
    {
        if( bzdcf )
        {
            hZDCF = ( TH1D* )c->GetListOfPrimitives()->FindObject( "hZDCF" );
        }
        else
        {
            hZDCF = ( TH1D* )c->GetListOfPrimitives()->FindObject( "hZDCF_sig" );
        }
        if( !hZDCF )
        {
            cout << "VZDCF::plot: no zdcf histogram found with name hZDCF " << endl;
            c->GetListOfPrimitives()->Print();
            return 0;
        }
    }
    
    // fill graph
    TGraphAsymmErrors* g = new TGraphAsymmErrors( 1 );
    setGraphPlottingStyle( ( TGraph* )g, 1, 1., 20, 1. );
    
    int z = 0;
    
    for( unsigned int i = 0; i < fZDCFData.size(); i++ )
    {
        if( fZDCFData[i] )
        {
            if( bzdcf )
            {
                g->SetPoint( z, fZDCFData[i]->tau, fZDCFData[i]->dcf );
                g->SetPointEXhigh( z, fZDCFData[i]->sigma_tau_pos );
                g->SetPointEXlow( z, fZDCFData[i]->sigma_tau_neg );
                g->SetPointEYhigh( z, fZDCFData[i]->dcf_error_up );
                g->SetPointEYlow( z, fZDCFData[i]->dcf_error_low );
                z++;
            }
            else
                // plot dcf / error
            {
                double e = 0.5 * ( fZDCFData[i]->dcf_error_up + fZDCFData[i]->dcf_error_low );
                if( e != 0. )
                {
                    g->SetPoint( z, fZDCFData[i]->tau, fZDCFData[i]->dcf / e );
                    g->SetPointEXhigh( z, fZDCFData[i]->sigma_tau_pos );
                    g->SetPointEXlow( z, fZDCFData[i]->sigma_tau_neg );
                    z++;
                }
            }
        }
    }
    
    // draw ML interval
    
    if( fMLPeakposition > -9998. )
    {
        TBox* iB = new TBox( fML1Sigmainterval_low, hZDCF->GetMinimum(), fML1Sigmainterval_up, hZDCF->GetMaximum() );
        iB->SetFillColor( 38 );
        iB->Draw();
        
        TLine* iL = new TLine( fMLPeakposition, hZDCF->GetMinimum(), fMLPeakposition, hZDCF->GetMaximum() );
        iL->SetLineStyle( 2 );
        iL->SetLineWidth( 2 );
        iL->Draw();
    }
    
    // draw graph
    g->Draw( "p" );
    
    return c;
}

TCanvas* VZDCF::plotZDCF( TCanvas* c, double taumin, double taumax, double ymin, double ymax )
{
    return plot( c, true, taumin, taumax, ymin, ymax );
}


TCanvas* VZDCF::plotZDCFoverError( TCanvas* c, double taumin, double taumax, double ymin, double ymax )
{
    return plot( c, false, taumin, taumax, ymin, ymax );
}

void VZDCF::setMLinterval( double iMLPeak, double iML1Sigma_low, double iML1Sigma_up )
{
    fMLPeakposition = iMLPeak;
    fML1Sigmainterval_low = iML1Sigma_low;
    fML1Sigmainterval_up = iML1Sigma_up;
}
