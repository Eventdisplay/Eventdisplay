/*! \class VEnergySpectrumfromLiterature
    \brief provide spectra from literature


*/

#include "VEnergySpectrumfromLiterature.h"

VEnergySpectrumfromLiterature::VEnergySpectrumfromLiterature( string ifile, bool iprint )
{
    bIsZombie = false;
    setFunctions();
    
    if( ifile.size() > 0 )
    {
        readValuesFromFile( ifile, iprint );
    }
    
    fIntegral_ID = 0;
    fIntegral_TF1 = 0;
    
    setPlottingLogEnergyAxis();
    setPlottingStyle();
    setPlottingMultiplierIndex();
    setPlottingEnergyRangeLinear();
    setPlottingYaxis();
    
}

/*
    set information about functions
*/
void VEnergySpectrumfromLiterature::setFunctions()
{
    sEnergyFun a;
    
    // power law
    a.Name = "PL";
    a.Description = "power law";
    a.NumParameters = 3;
    fEnergyFun.push_back( a );
    
    // power law with exponential cut off
    a.Name = "PLEC";
    a.Description = "power law with exponential cut off";
    a.NumParameters = 4;
    fEnergyFun.push_back( a );
    
    // curved power law
    a.Name = "VPL";
    a.Description = "curved power law";
    a.NumParameters = 4;
    fEnergyFun.push_back( a );
    
    // broken power law
    a.Name = "BRPL";
    a.Description = "broken power law";
    a.NumParameters = 6;
    fEnergyFun.push_back( a );
    
    // broken power law (2)
    a.Name = "BRPL_2";
    a.Description = "broken power law (2)";
    a.NumParameters = 4;
    fEnergyFun.push_back( a );
    
    // power law with exponential cut off (with beta factor for cut-off strength)
    a.Name = "PLEC_SF";
    a.Description = "power law with exponential cut off (+cut-off strength)";
    a.NumParameters = 5;
    fEnergyFun.push_back( a );
    
    // power law with exponential cut off (Gaussian factor)
    a.Name = "PLEC_GF";
    a.Description = "power law with exponential cut off (Gaussian factor)";
    a.NumParameters = 6;
    fEnergyFun.push_back( a );
}

bool VEnergySpectrumfromLiterature::readValuesFromFile( string ifile, bool iPrint )
{
    // clear existing data
    fData.clear();
    
    ifstream is;
    is.open( gSystem->ExpandPathName( ifile.c_str() ), ifstream::in );
    if( !is )
    {
        cout << "error reading literature values from " << ifile << endl;
        bIsZombie = true;
        return false;
    }
    string is_line;
    string is_temp;
    
    if( iPrint )
    {
        cout << "reading spectral parameters from " << ifile << endl;
    }
    
    sData itemp;
    
    // vectors for diff flux
    vector< double > diff_e;
    vector< double > diff_v;
    vector< double > diff_ve;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) != "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> is_temp;
        is_stream >> is_temp;
        if( is_temp == "NAME" )
        {
            is_stream >> itemp.Name;
            diff_e.clear();
            diff_v.clear();
            diff_ve.clear();
        }
        else if( is_temp == "OBSERVATORY" )
        {
            itemp.Observatory = is_line.substr( is_stream.tellg(), is_line.size() );
            itemp.Observatory = itemp.Observatory.substr( 1, itemp.Observatory.size() );
        }
        else if( is_temp == "SOURCE" )
        {
            itemp.Source = is_line.substr( is_stream.tellg(), is_line.size() );
            itemp.Source = itemp.Source.substr( 1, itemp.Source.size() );
        }
        else if( is_temp == "REFERENCE" )
        {
            itemp.Reference = is_line.substr( is_stream.tellg(), is_line.size() );
            itemp.Reference = itemp.Reference.substr( 1, itemp.Reference.size() );
        }
        else if( is_temp == "COMMENT" )
        {
            itemp.Comment = is_line.substr( is_stream.tellg(), is_line.size() );
            itemp.Comment = itemp.Comment.substr( 1, itemp.Comment.size() );
        }
        else if( is_temp == "VALUES" )
        {
            is_stream >> is_temp;
            itemp.Type = atoi( is_temp.c_str() );
            
            vector< double > v;
            vector< double > ve;
            
            if( itemp.Type < fEnergyFun.size() )
            {
                for( unsigned int i = 0; i < fEnergyFun[itemp.Type].NumParameters; i++ )
                {
                    is_stream >> is_temp;
                    v.push_back( atof( is_temp.c_str() ) );         // values
                    is_stream >> is_temp;
                    ve.push_back( atof( is_temp.c_str() ) );        // errors
                }
            }
            else
            {
                cout << "unknown type" << endl;
                continue;
            }
            itemp.Parameter = v;
            itemp.ParError = ve;
            // energy range
            is_stream >> is_temp;
            itemp.EnergyRange_min = atof( is_temp.c_str() );
            is_stream >> is_temp;
            itemp.EnergyRange_max = atof( is_temp.c_str() );
        }
        else if( is_temp == "DIFFFLUX" )
        {
            is_stream >> is_temp;
            diff_e.push_back( atof( is_temp.c_str() ) );
            is_stream >> is_temp;
            diff_v.push_back( atof( is_temp.c_str() ) );
            is_stream >> is_temp;
            diff_ve.push_back( atof( is_temp.c_str() ) );
        }
        else if( is_temp == "END" )
        {
            itemp.FluxV_energy = diff_e;
            itemp.FluxV_DiffFlux = diff_v;
            itemp.FluxV_DiffFluxError = diff_ve;
            
            fData.push_back( itemp );
        }
    }
    
    return true;
}


TGraphErrors* VEnergySpectrumfromLiterature::getDifferentialFluxPoints( unsigned int iID, bool bLogEnergy )
{
    if( !checkIDRange( iID ) )
    {
        return 0;
    }
    
    TGraphErrors* g = 0;
    
    if( fData[iID].FluxV_energy.size() > 0 && fData[iID].FluxV_energy.size() == fData[iID].FluxV_DiffFlux.size() && fData[iID].FluxV_energy.size() == fData[iID].FluxV_DiffFluxError.size() )
    {
        g = new TGraphErrors( 1 );
        g->SetMarkerColor( fPlottingColor );
        g->SetLineColor( fPlottingColor );
        g->SetLineWidth( ( Width_t )fPlottingLineWidth );
        g->SetLineStyle( fPlottingLineStyle );
        g->SetMarkerStyle( fPlottingMarkerStyle );
        g->SetMarkerSize( fPlottingMarkerSize );
        g->SetTitle( "" );
        
        for( unsigned int i = 0; i < fData[iID].FluxV_energy.size(); i++ )
        {
            if( bLogEnergy )
            {
                g->SetPoint( i, log10( fData[iID].FluxV_energy[i] ), fData[iID].FluxV_DiffFlux[i] * TMath::Power( fData[iID].FluxV_energy[i], fPlottingMultiplierIndex ) );
            }
            else
            {
                g->SetPoint( i, fData[iID].FluxV_energy[i], fData[iID].FluxV_DiffFlux[i] * TMath::Power( fData[iID].FluxV_energy[i], fPlottingMultiplierIndex ) );
            }
            g->SetPointError( i, 0., fData[iID].FluxV_DiffFluxError[i] * TMath::Power( fData[iID].FluxV_energy[i], fPlottingMultiplierIndex ) );
        }
    }
    
    return g;
}


TGraphAsymmErrors* VEnergySpectrumfromLiterature::getEnergySpectrumWithErrors( unsigned int iID, bool bLogEnergy )
{
    if( !checkIDRange( iID ) )
    {
        return 0;
    }
    
    TF1* f = getEnergySpectrum( iID, bLogEnergy );
    
    if( !f )
    {
        return 0;
    }
    TF1* fTemp = new TF1();
    f->Copy( *fTemp );
    
    // number of points in graph
    int nPoints = 100;
    // number of random cycles
    int nRandom = 10000;
    
    TGraphAsymmErrors* g = new TGraphAsymmErrors( nPoints );
    g->SetMarkerColor( fPlottingColor );
    g->SetLineColor( fPlottingColor );
    g->SetLineWidth( ( Width_t )fPlottingLineWidth );
    g->SetLineStyle( fPlottingLineStyle );
    g->SetMarkerStyle( fPlottingMarkerStyle );
    g->SetMarkerSize( fPlottingMarkerSize );
    g->SetFillStyle( fPlottingFillStyle );
    g->SetFillColor( fPlottingColor );
    g->SetTitle( "" );
    
    int nPara = f->GetNpar();
    
    double x = 0.;
    double ymin = 0.;
    double ymax = 0.;
    
    // energy axis range
    double xmin = 0.;
    double xmax = 0.;
    f->GetRange( xmin, xmax );
    if( xmax - xmin < 1.e-5 )
    {
        return 0;
    }
    
    if( nPoints < 2 )
    {
        return 0;
    }
    
    //////////////////////////////////////////////////////////////
    // make histograms
    if( fRandomErrorHistograms.size() == 0 )
    {
        char hname[600];
        for( int i = 0; i < nPoints; i++ )
        {
            x = xmin + i * ( xmax - xmin ) / ( double )( nPoints - 1 );
            // get minimum and maximum for histogram
            ymin = 0.01 * f->Eval( x );
            ymax = 100. * f->Eval( x );
            
            sprintf( hname, "i_RandomErrorHistogram_%d", i );
            fRandomErrorHistograms.push_back( new TH1D( hname, "", 10000, ymin, ymax ) );
        }
    }
    else
    {
        for( unsigned int i = 0; i < fRandomErrorHistograms.size(); i++ )
        {
            fRandomErrorHistograms[i]->Reset();
        }
    }
    
    /////////////////////////////////////////////////////////////////
    // now loop over all points and calculate error range
    
    for( int r = 0; r < nRandom; r++ )
    {
        // set random parameters
        for( int p = 0; p < nPara; p++ )
        {
            fTemp->SetParameter( p, gRandom->Gaus( f->GetParameter( p ), f->GetParError( p ) ) );
        }
        for( int i = 0; i < nPoints; i++ )
        {
            x = xmin + i * ( xmax - xmin ) / ( double )( nPoints - 1 );
            fRandomErrorHistograms[i]->Fill( fTemp->Eval( x ) );
        }
    }
    
    for( int i = 0; i < nPoints; i++ )
    {
        x = xmin + i * ( xmax - xmin ) / ( double )( nPoints - 1 );
        
        
        g->SetPoint( i, x, f->Eval( x ) );
        g->SetPointEYhigh( i, fRandomErrorHistograms[i]->GetRMS() );
        g->SetPointEYlow( i, fRandomErrorHistograms[i]->GetRMS() );
        
    }
    
    //Clear random error histograms. They should not be re-used since fPlottingMultiplierIndex may have changed since the last time this function was called.
    for( unsigned int i = 0; i < fRandomErrorHistograms.size(); i++ )
    {
        delete fRandomErrorHistograms[i];
    }
    fRandomErrorHistograms.clear();
    
    return g;
}


TF1* VEnergySpectrumfromLiterature::getEnergySpectrum( unsigned int iID, bool bLogEnergy, double iEnergyMin_Lin, double iEnergyMax_Lin )
{
    if( !checkIDRange( iID ) )
    {
        return 0;
    }
    
    TF1* f = 0;
    
    // define energy range
    double xmin = 0.;
    double xmax = 0.;
    if( iEnergyMin_Lin < 0. || iEnergyMax_Lin < 0. )
    {
        if( bLogEnergy )
        {
            xmin = log10( fData[iID].EnergyRange_min );
            xmax = log10( fData[iID].EnergyRange_max );
        }
        else
        {
            xmin = fData[iID].EnergyRange_min;
            xmax = fData[iID].EnergyRange_max;
        }
    }
    else
    {
        if( bLogEnergy )
        {
            xmin = log10( iEnergyMin_Lin );
            xmax = log10( iEnergyMax_Lin );
        }
        else
        {
            xmin = iEnergyMin_Lin;
            xmax = iEnergyMax_Lin;
        }
    }
    
    // define functions
    char hname[6000];
    char h_exponent[600];
    char h_energy[60];
    // define energy variabile (log or lin)
    if( bLogEnergy )
    {
        sprintf( h_energy, "TMath::Power( 10, x )" );
    }
    else
    {
        sprintf( h_energy, "x" );
    }
    
    if( fData[iID].Type < fEnergyFun.size() )
    {
        if( fEnergyFun[fData[iID].Type].NumParameters == fData[iID].Parameter.size() && fEnergyFun[fData[iID].Type].NumParameters == fData[iID].ParError.size() )
        {
            // power law
            if( fData[iID].Type == 0 )
            {
                sprintf( hname, "[0] * TMath::Power( %s / %f, [1] ) ", h_energy, fData[iID].Parameter[0] );
            }
            // power law with exponential cut off
            else if( fData[iID].Type == 1 )
            {
                sprintf( h_exponent, "TMath::Exp( -1.* %s  / [1] )", h_energy );
                sprintf( hname, "[0] * TMath::Power( %s / %f, [2] ) * %s", h_energy, fData[iID].Parameter[0], h_exponent );
            }
            // curved power law (e.g. ApJ 674, 1037 (2008))
            else if( fData[iID].Type == 2 )
            {
                sprintf( h_exponent, "[1]+[2]*TMath::Log10( %s / %f )", h_energy, fData[iID].Parameter[0] );
                sprintf( hname, "[0] * TMath::Power( %s / %f, %s )", h_energy, fData[iID].Parameter[0], h_exponent );
            }
            // broken power law (e.g. A&A 457, 899 (2006))
            else if( fData[iID].Type == 3 )
            {
                sprintf( h_exponent, "[1]*TMath::Power( %s / [0], [2] ) ", h_energy );
                sprintf( hname, "%s * TMath::Power( 1. + TMath::Power( %s / [0], 1./[4] ), [4]*([3]-[2]) )", h_exponent, h_energy );
            }
            // broken power law (2)
            else if( fData[iID].Type == 4 )
            {
                sprintf( h_exponent, "[1]*TMath::Power( %s / %f, -1. ) ", h_energy, fData[iID].Parameter[0] );
                sprintf( hname, "%s * TMath::Power( 1. + TMath::Power( %s / [0], [2] ), -1. )", h_exponent, h_energy );
            }
            // power law with exponential cut off (with beta factor for cut-off strength)
            else if( fData[iID].Type == 5 )
            {
                sprintf( h_exponent, "TMath::Exp( -1.* TMath::Power( %s  / [1], [3] ) )", h_energy );
                sprintf( hname, "[0] * TMath::Power( %s / %f, [2] ) * %s", h_energy, fData[iID].Parameter[0], h_exponent );
            }
            // power law with exponential cut off (Gaussian factor)
            else if( fData[iID].Type == 6 )
            {
                sprintf( h_exponent, "[0]*TMath::Power( %s / %f, [4] ) ", h_energy, fData[iID].Parameter[0] );
                sprintf( hname, "%s * (1.+[1]*(exp(TMath::Gaus(log10(%s),[2], [3] ))-1.))", h_exponent, h_energy );
            }
            
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
    
    // plotting multiplier
    char hisnmae[10000];
    sprintf( hisnmae, "%s * TMath::Power( %s, %f )", hname, h_energy, fPlottingMultiplierIndex );
    
    // create function
    f = new TF1( fData[iID].Name.c_str(), hisnmae, xmin, xmax );
    
    // set parameters
    for( unsigned int i = 1; i < fData[iID].Parameter.size(); i++ )
    {
        f->SetParameter( i - 1, fData[iID].Parameter[i] );
    }
    for( unsigned int i = 1; i < fData[iID].ParError.size(); i++ )
    {
        f->SetParError( i - 1, fData[iID].ParError[i] );
    }
    
    if( f )
    {
        f->SetLineColor( fPlottingColor );
        f->SetLineWidth( ( Width_t )fPlottingLineWidth );
        f->SetLineStyle( fPlottingLineStyle );
        f->SetTitle( "" );
    }
    
    return f;
}

/*

   energy on linear scale [TeV]

*/
bool VEnergySpectrumfromLiterature::prepare_integration( unsigned int iID, double iEmin, double iEmax )
{
    fIntegral_TF1 = getEnergySpectrum( iID, false, iEmin, iEmax );
    if( !fIntegral_TF1 )
    {
        return false;
    }
    
    fIntegral_TF1->CalcGaussLegendreSamplingPoints( 1000, fIntegral_x, fIntegral_y, 1.e-15 );
    
    // save that this function has already been integrated
    fIntegral_ID = iID;
    
    return true;
}


/*
   energy on linear scale [TeV]

   observe: somehow sensitive to upper limit for overal integrated flux)

   integrated flux in [1/cm2/s]
*/
double VEnergySpectrumfromLiterature::getIntegralFlux( double iEmin, double iEmax,  unsigned int iID )
{
    if( fIntegral_ID != iID || !fIntegral_TF1 )
    {
        if( !prepare_integration( iID, iEmin, iEmax ) )
        {
            cout << "VEnergySpectrumfromLiterature::getIntegralFlux() error: preparation for integral failed" << endl;
        }
    }
    
    return fIntegral_TF1->IntegralFast( 1000, fIntegral_x, fIntegral_y, iEmin, iEmax );
}


void VEnergySpectrumfromLiterature::listValues()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        listValues( i );
    }
}

void VEnergySpectrumfromLiterature::listValues( unsigned int i )
{
    if( i >= fData.size() )
    {
        cout << "VEnergySpectrumfromLiterature::listValues error: index out of range (<" << fData.size() << ")" << endl;
        return;
    }
    
    char hname[800];
    
    cout << "ID: " << i << ", " << fData[i].Source << ", " << fData[i].Observatory <<  " (";
    cout << fData[i].Name << ")" << endl;
    if( fData[i].Type == 0 && fData[i].Parameter.size() == 3 && fData[i].ParError.size() == 3 )
    {
        sprintf( hname, "\t(power law dN/dE = I x (E/%.2f)^Gamma):", fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma = (%.3f +- %.3f)", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
    }
    else if( fData[i].Type == 1 && fData[i].Parameter.size() == 4 && fData[i].ParError.size() == 4 )
    {
        sprintf( hname, "\t(power law with exponential cutoff dN/dE = I x (E/%.2f)^Gamma) x exp(-E/E_c):", fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma = (%.3f +- %.3f)", fData[i].Parameter[3], fData[i].ParError[3] );
        cout << hname;
        sprintf( hname, ", E_C = (%.3f +- %.3f)", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
    }
    else if( fData[i].Type == 2 && fData[i].Parameter.size() == 4 && fData[i].ParError.size() == 4 )
    {
        sprintf( hname, "\t(curved power law dN/dE = I x (E/%.2f)^(Gamma + b x log10( E/%.2f)):", fData[i].Parameter[0], fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma = (%.3f +- %.3f)", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
        sprintf( hname, ", b = (%.3f +- %.3f)", fData[i].Parameter[3], fData[i].ParError[3] );
        cout << hname;
    }
    else if( fData[i].Type == 3 && fData[i].Parameter.size() == 6 && fData[i].ParError.size() == 6 )
    {
        sprintf( hname, "\t(broken power law dN/dE = I x (E/%.2f)^(Gamma_1) x ( 1 + (E/E_C)^(1/s) )^(s x (Gamma_1-Gamma_2): ", fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
        sprintf( hname, ", E_c = (%.2f +- %.2f) TeV", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma_1 = (%.3f +- %.3f)", fData[i].Parameter[3], fData[i].ParError[3] );
        cout << hname;
        sprintf( hname, ", Gamma_2 = (%.3f +- %.3f)", fData[i].Parameter[4], fData[i].ParError[4] );
        cout << hname;
        sprintf( hname, ", s = (%.3f +- %.3f)", fData[i].Parameter[5], fData[i].ParError[5] );
        cout << hname;
    }
    else if( fData[i].Type == 4 && fData[i].Parameter.size() == 4 && fData[i].ParError.size() == 4 )
    {
        sprintf( hname, "\t(broken power law (2) dN/dE = I x (E/%.2f)^(-1) / ( 1 + (E/E_C)^(s) ): ", fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
        sprintf( hname, ", E_c = (%.3f +- %.3f) TeV", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma = (%.3f +- %.3f)", fData[i].Parameter[3], fData[i].ParError[3] );
        cout << hname;
    }
    else if( fData[i].Type == 6 && fData[i].Parameter.size() == 6 && fData[i].ParError.size() == 6 )
    {
        sprintf( hname, "\t(power law with exponential cut off (Gaussian factor): dN/dE = I x (E/%.2f)^(Gamma_1) / ( 1 + f * ( e^-(Gauss( log(E/%.2f), Mean, Sigma ) ) -1 ) )", fData[i].Parameter[0], fData[i].Parameter[0] );
        cout << hname << endl;
        cout << "\t";
        sprintf( hname, "I = (%.3e +- %.3e) cm^-2 s^-1 TeV^-1", fData[i].Parameter[1], fData[i].ParError[1] );
        cout << hname;
        sprintf( hname, ", Gamma = (%.3f +- %.3f)", fData[i].Parameter[5], fData[i].ParError[5] );
        cout << hname << endl;
        sprintf( hname, "\t Mean = (%.3f +- %.3f)", fData[i].Parameter[3], fData[i].ParError[3] );
        cout << hname;
        sprintf( hname, ", Sigma = (%.3f +- %.3f)", fData[i].Parameter[4], fData[i].ParError[4] );
        cout << hname;
        sprintf( hname, ", f = (%.3f +- %.3f)", fData[i].Parameter[2], fData[i].ParError[2] );
        cout << hname;
    }
    cout << fixed << setprecision( 2 ) << endl;
    cout << "\t " << fData[i].EnergyRange_min << " < E [TeV] < " << fData[i].EnergyRange_max << endl;
    if( fData[i].FluxV_energy.size() > 0 )
    {
        cout << "\t\t energy \t         Flux        \t      FluxError" << endl;
        cout << "\t\t  [TeV] \t [cm^-2 s^-1 TeV^-1] \t [cm^-2 s^-1 TeV^-1]" << endl;
        for( unsigned int t = 0; t < fData[i].FluxV_energy.size(); t++ )
        {
            cout << "\t\t" << setprecision( 4 ) << fixed;
            cout << fData[i].FluxV_energy[t] << "\t\t";
            cout << setprecision( 2 ) << scientific;
            cout << fData[i].FluxV_DiffFlux[t] << "\t\t";
            cout << fData[i].FluxV_DiffFluxError[t];
            cout << resetiosflags( ios_base::scientific ) << endl;
        }
    }
    cout << "\t(" << fData[i].Reference << ")" << endl;
    cout << "\t(" << fData[i].Comment << ")" << endl;
}


TCanvas* VEnergySpectrumfromLiterature::plot( unsigned int iID, TCanvas* c, bool iLogY )
{
    TF1* f = getEnergySpectrum( iID, true );
    if( !f )
    {
        return 0;
    }
    
    char hname[600];
    char htitle[600];
    TH1D* hNull = 0;
    
    if( c == 0 )
    {
        sprintf( hname, "c_%s", fData[iID].Name.c_str() );
        sprintf( htitle, "energy spectrum (%s)", fData[iID].Source.c_str() );
        c = new TCanvas( hname, htitle, 10, 10, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        gPad->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hnull_%s", fData[iID].Name.c_str() );
        hNull = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNull->SetMinimum( fPlottingYaxisMin );
        hNull->SetMaximum( fPlottingYaxisMax );
        hNull->SetStats( 0 );
        hNull->SetXTitle( "log_{10} energy [TeV]" );
        hNull->SetYTitle( "dN/dE [cm^{-2}s^{-1}TeV^{-1}]" );
        if( fPlottingMultiplierIndex > 1. )
        {
            sprintf( hname, "E^{%.2f} dN/dE [cm^{-2}s^{-1}TeV^{%.2f}]", fPlottingMultiplierIndex, fPlottingMultiplierIndex - 1. );
            hNull->SetYTitle( hname );
        }
        hNull->GetYaxis()->SetTitleOffset( 1.6 );
        
        plot_nullHistogram( c, hNull, fPlottingLogEnergyAxis, iLogY, hNull->GetYaxis()->GetTitleOffset(), fPlottingMinEnergy, fPlottingMaxEnergy );
        c->SetLogy( iLogY );
    }
    else
    {
        c->cd();
    }
    
    f->Draw( "same" );
    
    return c;
}


bool VEnergySpectrumfromLiterature::checkIDRange( unsigned int iID )
{
    if( iID < fData.size() )
    {
        return true;
    }
    
    return false;
}


TCanvas* VEnergySpectrumfromLiterature::plot( string iSelection, TCanvas* c )
{
    if( iSelection.size() == 0 )
    {
        return 0;
    }
    
    vector< unsigned int > v_id;
    
    // read list (seperated by spaces)
    
    string itemp;
    istringstream is_stream( iSelection );
    while( !(is_stream>>std::ws).eof() )
    {
        is_stream >> itemp;
        v_id.push_back( ( unsigned int )atoi( itemp.c_str() ) );
    }
    
    for( unsigned int i = 0; i < v_id.size(); i++ )
    {
        setPlottingStyle( i + 1, gRandom->Integer( 4 ) );
        
        c = plot( v_id[i], c );
    }
    
    return c;
}


sData VEnergySpectrumfromLiterature::getEnergySpectrumDataField( unsigned int iID )
{
    if( !checkIDRange( iID ) )
    {
        sData a;
        return a;
    }
    return fData[iID];
}

/*
 * get spectral index if the selected spectrum is a power law (otherwise return -9999)
 */
double VEnergySpectrumfromLiterature::getPowerLaw_Index( unsigned int iID )
{
    if( !checkIDRange( iID ) )
    {
        return -9999.;
    }
    
    if( fData[iID].Type != 0 || fData[iID].Parameter.size() != 3 )
    {
        return -9999.;
    }
    
    return fData[iID].Parameter[2];
}

/*
 * get flux constant if the selected spectrum is a power law (otherwise return -9999)
 */
double VEnergySpectrumfromLiterature::getPowerLaw_FluxConstant_at1TeV( unsigned int iID )
{
    if( !checkIDRange( iID ) )
    {
        return -9999.;
    }
    
    if( fData[iID].Type != 0 || fData[iID].Parameter.size() != 3 || fData[iID].Parameter[0] == 0. )
    {
        return -9999.;
    }
    
    return ( fData[iID].Parameter[1] * TMath::Power( 1. / fData[iID].Parameter[0], fData[iID].Parameter[2] ) );
}
