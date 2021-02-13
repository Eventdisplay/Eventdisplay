/*! \class VSpectralFitter
    \brief  fitter class for energy spectra


*/

#include "VSpectralFitter.h"

ClassImp( VSpectralFitter )

VSpectralFitter::VSpectralFitter( string fitname )
{
    fFitFunction = 0;
    fFitFunction_lin = 0;
    fFitFunction_CovarianceMatrix = 0;
    fFitName = fitname;
    
    setSpectralFitFunction();
    setSpectralFitFluxNormalisationEnergy();
    setSpectralFitRangeLin();
    setPlottingStyle();
}


/*

    fit the current function

    (note the fit options)

    print out covariance matrix and correlation coefficient (for dim=2)


*/
TF1* VSpectralFitter::fit( TGraph* g, string fitname )
{
    if( !g )
    {
        cout << "VSpectralFitter::fit warning: no graph" << endl;
        return 0;
    }
    if( fitname.size() > 0 )
    {
        fFitName = fitname;
    }
    
    // define fit function
    defineFitFunction();
    
    // fit
    if( fSpectralFitFunction == 2 ) // i.e. for the broken power law
    {
        g->Fit( fFitFunction, "0MNREV" );    // more verbose
    }
    else
    {
        g->Fit( fFitFunction, "0MNER" );
    }
    
    // gets the current default fitter
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    
    if( fSpectralFitFunction == 2 ) // i.e. for the broken power law
    {
        // Change the tolerance of the EDM, and loosen the convergence condition, since BPL has more free parameters
        fitter->SetPrecision( 1E-1 ); // Note the EDM is 0.001 * 1E-2 * (some constant), otherwise there is problem with convergence
        cout << "Setting default tolerance EDM to ~ 1.E-4" << endl;
        cout << "Coder's WARNING: Use this with care, and check the output of the fitter and the Error Matrix is not from MINOS" << endl;
    }
    
    // covariance matrix
    if( fitter && fitter->GetCovarianceMatrix() )
    {
        int nPars = fitter->GetNumberFreeParameters();
        TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );
        if( COV )
        {
            fFitFunction_CovarianceMatrix = COV->GetMatrixArray();
            cout << "Fitter: " << endl;
            cout << "\tCovariance matrix ";
            cout << "(nxn=" << fitter->GetNumberTotalParameters() << "x" << fitter->GetNumberTotalParameters() << "): " << endl;
            for( int i = 0; i < fitter->GetNumberTotalParameters(); i++ )
            {
                for( int j = 0; j < fitter->GetNumberTotalParameters(); j++ )
                {
                    cout << "\t" << fitter->GetCovarianceMatrixElement( i, j );
                }
                cout << endl;
            }
            // calculate correlation coefficient
            if( fitter->GetNumberTotalParameters() == 2 )
            {
                if( fitter->GetCovarianceMatrixElement( 0, 0 ) > 0. && fitter->GetCovarianceMatrixElement( 1, 1 ) > 0. )
                {
                    double rho = TMath::Abs( fitter->GetCovarianceMatrixElement( 0, 1 ) );
                    rho /= sqrt( fitter->GetCovarianceMatrixElement( 0, 0 ) * fitter->GetCovarianceMatrixElement( 1, 1 ) );
                    
                    cout << "\tCorrelation coefficient: " << rho << endl;
	   	    double i_decorr = fSpectralFitFluxNormalisationEnergy ;
		    i_decorr *= TMath::Exp( fitter->GetCovarianceMatrixElement( 0, 1 ) / fFitFunction->GetParameter(0) / fitter->GetCovarianceMatrixElement( 1, 1 )  );
		    cout << "\tDecorrelation Energy: " << i_decorr << " " << fSpectralFitFluxNormalisationEnergy << endl;
                }
            }
        }
    }
    else
    {
        cout << "VSpectralFitter::fit: no covariance matrix" << endl;
    }
    
    updateFitFunction_lin();
    
    return fFitFunction;
}


/*

   define fit function

   These are predifined and can be set with a function ID:

   fSpectralFitFunction == 0 :  power law
   fSpectralFitFunction == 1 :  power law with exponential cutoff
   fSpectralFitFunction == 2 :  broken power law
   fSpectralFitFunction == 3 :  curved power law
*/
bool VSpectralFitter::defineFitFunction()
{
    char hname[2000];
    string iFitName_lin = fFitName + "_lin";
    
    /////////////////////////////////////////////////////////
    // power law
    if( fSpectralFitFunction == 0 )
    {
        cout << "Fitfunction: power law" << endl;
        sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] )", fSpectralFitFluxNormalisationEnergy );
        fFitFunction = new TF1( fFitName.c_str(), hname, log10( fSpectralFitEnergy_min ), log10( fSpectralFitEnergy_max ) );
        fFitFunction->SetParameter( 0, 1.e-7 );
        fFitFunction->SetParameter( 1, -2.5 );
        // linear axis
        sprintf( hname, "[0] * TMath::Power( x/ %f, [1] )", fSpectralFitFluxNormalisationEnergy );
        fFitFunction_lin = new TF1( iFitName_lin.c_str(), hname, fSpectralFitEnergy_min, fSpectralFitEnergy_max );
    }
    /////////////////////////////////////////////////////////
    // power law with exponential cutoff
    else if( fSpectralFitFunction == 1 )
    {
        cout << "Fitfunction: power law with exponential cutoff" << endl;
        sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] ) * TMath::Exp( -1. * TMath::Power( 10, x ) / [2] )", fSpectralFitFluxNormalisationEnergy );
        fFitFunction = new TF1( fFitName.c_str(), hname, log10( fSpectralFitEnergy_min ), log10( fSpectralFitEnergy_max ) );
        fFitFunction->SetParameter( 0, 1.e-7 );
        fFitFunction->SetParameter( 1, -2. );
        fFitFunction->SetParameter( 2, 10. );
        sprintf( hname, "[0] * TMath::Power( x  / %f, [1] ) * TMath::Exp( -1. * x / [2] )", fSpectralFitFluxNormalisationEnergy );
        fFitFunction_lin = new TF1( iFitName_lin.c_str(), hname, fSpectralFitEnergy_min, fSpectralFitEnergy_max );
    }
    // broken power law fit
    else if( fSpectralFitFunction == 2 )
    {
        cout << "Fitfunction: broken power law" << endl;
        sprintf( hname, "((TMath::Power( 10, x )<[3]) * [0] * TMath::Power( TMath::Power( 10, x ) / [3], [1] )) + ((TMath::Power( 10, x )>=[3]) * [0] * TMath::Power( TMath::Power( 10, x ) / [3], [2] ))" );
        fFitFunction = new TF1( fFitName.c_str(), hname, log10( fSpectralFitEnergy_min ), log10( fSpectralFitEnergy_max ) );
        fFitFunction->SetParameter( 0, 2.e-10 );
        fFitFunction->SetParameter( 1, -1.5 );
        fFitFunction->SetParameter( 2, -3.0 );
        fFitFunction->SetParameter( 3, 0.7 );
        //fFitFunction->SetParLimits( 1, -10., 10. );
        //fFitFunction->SetParLimits( 2, -10., 10. );
        //fFitFunction->SetParLimits( 3, 0.05, 100. );
        // linear axis
        sprintf( hname, "((x<[3]) * [0] * TMath::Power( x/ [3], [1] )) + ((x>=[3]) * [0] * TMath::Power( x/ [3], [2] ))" );
        fFitFunction_lin = new TF1( iFitName_lin.c_str(), hname, fSpectralFitEnergy_min, fSpectralFitEnergy_max );
    }
    // curved power law fit
    else if( fSpectralFitFunction == 3 )
    {
		cout << "Fitfunction: curved power law fit" << endl;
        sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Power( 10, x ) / %f )", fSpectralFitFluxNormalisationEnergy, fSpectralFitFluxNormalisationEnergy );
        fFitFunction = new TF1( fFitName.c_str(), hname, log10( fSpectralFitEnergy_min ), log10( fSpectralFitEnergy_max ) );
        fFitFunction->SetParameter( 0, 1.e-7 );
        fFitFunction->SetParameter( 1, -2. );
        fFitFunction->SetParameter( 2, -0.01 );
        fFitFunction->SetParLimits( 2, -0.3, 0. );
    }
	// curved power law fit  (e.g. ApJ 674, 1037 (2008))
	else if( fSpectralFitFunction == 4 )
	{
		cout << "Fitfunction: curved power law fit (4)" << endl;
		sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*log10( TMath::Power( 10, x ) / %f ) )", fSpectralFitFluxNormalisationEnergy, fSpectralFitFluxNormalisationEnergy );
		fFitFunction = new TF1( fFitName.c_str(), hname, log10( fSpectralFitEnergy_min ), log10( fSpectralFitEnergy_max ) );
		fFitFunction->SetParameter( 0, 1.e-7 );
		fFitFunction->SetParameter( 1, -2. );
		fFitFunction->SetParameter( 2, -0.01 );
        fFitFunction->SetParLimits( 2, -0.3, 0. );
    }
    else
    {
        cout << "VSpectralFitter::defineFitFunction: unknown spectral fit function: " << fSpectralFitFunction << endl;
        return false;
    }
    
    // set all parameters for the function with linear energy axis
    updateFitFunction_lin();
    
    // set plotting style
    if( fFitFunction )
    {
        fFitFunction->SetLineStyle( fPlottingEnergySpectrumLineStyle );
        fFitFunction->SetLineColor( fPlottingEnergySpectrumLineColor );
        fFitFunction->SetLineWidth( ( Width_t )fPlottingEnergySpectrumLineWidth );
    }
    
    return true;
}

void VSpectralFitter::updateFitFunction_lin()
{
    if( !fFitFunction || !fFitFunction_lin )
    {
        return;
    }
    
    for( int i = 0; i < fFitFunction->GetNpar(); i++ )
    {
        fFitFunction_lin->SetParameter( i, fFitFunction->GetParameter( i ) );
        fFitFunction_lin->SetParError( i, fFitFunction->GetParError( i ) );
    }
}


void VSpectralFitter::print()
{
    if( !fFitFunction )
    {
        return;
    }
    
    cout << endl;
    if( fSpectralFitFunction == 0 )
    {
        cout << "Results for power law fit: " << endl;
        cout << "--------------------------" << endl;
        cout << "dN/dE = I x (E/" << fSpectralFitFluxNormalisationEnergy << " TeV)^{-Gamma}" << endl;
        cout << "I = " << scientific << setprecision( 2 ) << fFitFunction->GetParameter( 0 );
        cout << " +- " << fFitFunction->GetParError( 0 ) << " cm^-2s^-1TeV^-1" << endl;
        cout << "Gamma = " << fixed << setprecision( 2 ) << fFitFunction->GetParameter( 1 );
        cout << " +- " << fFitFunction->GetParError( 1 ) << endl;
        cout << "Chi2 " << setprecision( 2 ) << fFitFunction->GetChisquare();
        cout << ", N = " << fFitFunction->GetNDF();
        if( fFitFunction->GetNDF() > 0. )
        {
            cout << " (Chi2/N=" << fFitFunction->GetChisquare() / fFitFunction->GetNDF() << ")" << endl;
        }
        cout << endl;
    }
}

/*

   get integral flux from fit function

   iMinEnergy_TeV: threshold for integral flux in TeV (lin axis)

*/
double VSpectralFitter::getIntegralFlux( double iMinEnergy_TeV, double iMaxEnergy_TeV )
{
    if( !fFitFunction_lin )
    {
        cout << "VSpectralFitter::getIntegralFlux(): error: no fit function" << endl;
        return -99999.;
    }
    
    // analytical calculation for power law (fSpectralFitFunction == 0)
    /*   if( fSpectralFitFunction == 0 )
       {
           if( fFitFunction_lin->GetParameter( 1 ) != -1. )
           {
              double iFL = -1.* fFitFunction_lin->GetParameter( 0 ) / (fFitFunction_lin->GetParameter( 1 ) + 1.) / fSpectralFitFluxNormalisationEnergy;
    	  iFL *= TMath::Power( iMinEnergy_TeV / fSpectralFitFluxNormalisationEnergy, fFitFunction_lin->GetParameter( 1 ) + 1. );
    	  return iFL;
           }
        }     */
    
    return fFitFunction_lin->Integral( iMinEnergy_TeV, iMaxEnergy_TeV, 1.e-30 );
}

/*

   get integral flux error from fit function

   iMinEnergy_TeV: threshold for integral flux in TeV (lin axis)

*/
double VSpectralFitter::getIntegralFluxError( double iMinEnergy_TeV, double iMaxEnergy_TeV )
{
    if( !fFitFunction_lin || !fFitFunction )
    {
        cout << "VSpectralFitter::getIntegralFlux(): error: no fit function" << endl;
        return -99999.;
    }
    
    return fFitFunction_lin->IntegralError( iMinEnergy_TeV, iMaxEnergy_TeV, fFitFunction_lin->GetParameters(), fFitFunction_CovarianceMatrix );
}

