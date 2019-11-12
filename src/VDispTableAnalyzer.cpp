/*! \file VDispTableAnalyzer
    \brief get DISP values from tables


*/

#include "VDispTableAnalyzer.h"

VDispTableAnalyzer::VDispTableAnalyzer( string iFile )
{
    bDebug = false;
    bZombie = true;
    
    f_disp = 0.;
    f_dispE = 0.;
    f_disp_Phi = 0.;
    f_disp_PhiE = 0.;
    f_disp_Miss = 0.;
    
    // open file with disp tables
    fFile = new TFile( iFile.c_str() );
    cout << "opening disp table file: " << fFile->GetName() << endl;
    if( fFile->IsZombie() )
    {
        cout << "VDispTableAnalyzer error: cannot open file with disp table: " << iFile << endl;
        bZombie = true;
        return;
    }
    // read disp tables from file
    fData = ( VDispTableReader* )fFile->Get( "dispTable" );
    if( !fData )
    {
        cout << "VDispTableAnalyzer error: cannot find disp table in file " << iFile << endl;
        bZombie = true;
        return;
    }
    fData->initialize( true );
    fData->print( false );
    // at this stage everything should be fine
    bZombie = false;
}

/*
       get disp value for a given set of ze, az and noise values

       NOTE:

       interpolate only in zenith angles (not in azimuth and pedvar)

*/
float VDispTableAnalyzer::evaluate( float iWidth, float iLength, float iSize, float iPedvar, float iZe, float iAz, bool b2D )
{
    // get lower bound ze
    float i_ze_low = fData->getLowerZe( iZe );
    
    float disp_low      = calculateDisp( iWidth, iLength, iSize, iPedvar, i_ze_low, iAz, b2D );
    float dispE_low     = f_dispE;
    float disp_Phi_low  = f_disp_Phi;
    float disp_PhiE_low = f_disp_PhiE;
    float disp_Miss_low = f_disp_Miss;
    
    // get upper bound ze
    float i_ze_upp = fData->getUpperZe( iZe );
    
    float disp_upp      = calculateDisp( iWidth, iLength, iSize, iPedvar, i_ze_upp, iAz, b2D );
    float dispE_upp     = f_dispE;
    float disp_Phi_upp  = f_disp_Phi;
    float disp_PhiE_upp = f_disp_PhiE;
    float disp_Miss_upp = f_disp_Miss;
    
    // interpolate
    
    f_disp      = interpolate( disp_low, i_ze_low, disp_upp, i_ze_upp, iZe, true );
    f_dispE     = interpolate( dispE_low, i_ze_low, dispE_upp, i_ze_upp, iZe, true );
    f_disp_Phi  = interpolate( disp_Phi_low, i_ze_low, disp_Phi_upp, i_ze_upp, iZe, true );
    f_disp_PhiE = interpolate( disp_PhiE_low, i_ze_low, disp_PhiE_upp, i_ze_upp, iZe, true );
    f_disp_Miss = interpolate( disp_Miss_low, i_ze_low, disp_Miss_upp, i_ze_upp, iZe, true );
    
    return f_disp;
}


double VDispTableAnalyzer::interpolate( double w1, double ze1, double w2, double ze2, double ze, bool iCos )
{
    // don't interpolate if one or two values are not valid
    if( w1 < -90. || w2 < -90. )
    {
        return -99.;
    }
    
    // same zenith angle, don't interpolate
    if( TMath::Abs( ze1 - ze2 ) < 1.e-3 )
    {
        return ( w1 + w2 ) / 2.;
    }
    
    // interpolate
    double id, f1, f2;
    if( iCos )
    {
        id = cos( ze1 * TMath::DegToRad() ) - cos( ze2 * TMath::DegToRad() );
        f1 = 1. - ( cos( ze1 * TMath::DegToRad() ) - cos( ze * TMath::DegToRad() ) ) / id;
        f2 = 1. - ( cos( ze * TMath::DegToRad() ) - cos( ze2 * TMath::DegToRad() ) ) / id;
    }
    else
    {
        id = ze1 - ze2;
        f1 = 1. - ( ze1 - ze ) / id;
        f2 = 1. - ( ze  - ze2 ) / id;
    }
    
    // one of the values is not valid:
    // return valid value only when f1/f2 > 0.5 (this value is almost randomly chosen)
    if( w1 > -90. && w2 < -90. )
    {
        if( f1 > 0.5 )
        {
            return w1;
        }
        else
        {
            return -99.;
        }
    }
    else if( w1 < -90. && w2 > -90. )
    {
        if( f2 > 0.5 )
        {
            return w2;
        }
        else
        {
            return -99.;
        }
    }
    
    return ( w1 * f1 + w2 * f2 );
}

/*!
    calculate disp as function of width, length, size, pedvar, ze, az
*/
float VDispTableAnalyzer::calculateDisp( float iWidth, float iLength, float iSize, float iPedvar, float iZe, float iAz, bool b2D )
{
    // check for valid size
    if( iSize > 0. )
    {
        iSize = log10( iSize );
    }
    else
    {
        return -99.;
    }
    
    // check for valid length
    float iWidthOLenght = -99.;
    if( iLength > 0. )
    {
        iWidthOLenght = iWidth / iLength;
    }
    
    // get corresponding tree entry (ignore wobble offsets)
    int iEntry = fData->getTreeEntryFinder( iZe, iAz, 0., iPedvar, 0 );
    
    // reset parameter values
    f_disp = 0.;
    f_dispE = 1.e3;
    f_disp_Phi = 0.;
    f_disp_PhiE = 1.e3;
    f_disp_Miss = 0.;
    
    // read disp from tree
    if( fData->getTree() )
    {
        if( iEntry <= fData->getTree()->GetEntries() )
        {
            fData->getTree()->GetEntry( iEntry );
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // use width/length histograms
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if( b2D )
            {
                if( fData->h2D_DispTable && fData->h2D_DispTableN )
                {
                    // find disp bin
                    int x_bin = fData->h2D_DispTable->GetXaxis()->FindBin( iWidthOLenght );
                    int y_bin = fData->h2D_DispTable->GetYaxis()->FindBin( iSize );
                    // require at least 5 events per bin
                    if( fData->h2D_DispTableN->GetBinContent( x_bin, y_bin ) > 5. )
                    {
                        f_disp = fData->h2D_DispTable->GetBinContent( x_bin, y_bin );
                        f_dispE = fData->h2D_DispTable->GetBinError( x_bin, y_bin );
                        f_disp_Phi = fData->h2D_DispPhiTable->GetBinContent( x_bin, y_bin );
                        f_disp_PhiE = fData->h2D_DispPhiTable->GetBinError( x_bin, y_bin );
                        f_disp_Miss = fData->h2D_DispMissTable->GetBinContent( x_bin, y_bin );
                        return f_disp;
                    }
                    else
                    {
                        return -99.;
                    }
                }
                else
                {
                    cout << "VDispTableAnalyzer::evaluate: error no data histograms found" << endl;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // use width and length histograms
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            else
            {
                if( fData->h3D_DispTable && fData->h3D_DispTableN )
                {
                    // find disp bin
                    int x_bin = fData->h3D_DispTable->GetXaxis()->FindBin( fData->scaleWidthParameter( iWidth ) );
                    int y_bin = fData->h3D_DispTable->GetYaxis()->FindBin( fData->scaleLengthParameter( iLength ) );
                    int z_bin = fData->h3D_DispTable->GetZaxis()->FindBin( iSize );
                    // require at least 5 events per bin
                    if( fData->h3D_DispTableN->GetBinContent( x_bin, y_bin, z_bin ) > 5. )
                    {
                        f_disp = fData->h3D_DispTable->GetBinContent( x_bin, y_bin, z_bin );
                        f_dispE = fData->h3D_DispTable->GetBinError( x_bin, y_bin, z_bin );
                        f_disp_Phi = fData->h3D_DispPhiTable->GetBinContent( x_bin, y_bin, z_bin );
                        f_disp_PhiE = fData->h3D_DispPhiTable->GetBinError( x_bin, y_bin, z_bin );
                        f_disp_Miss = fData->h3D_DispMissTable->GetBinContent( x_bin, y_bin, z_bin );
                        return f_disp;
                    }
                    else
                    {
                        return -99.;
                    }
                }
                else
                {
                    cout << "VDispTableAnalyzer::evaluate: error no data histograms found" << endl;
                }
            }
        }
        else
        {
            cout << "VDispTableAnalyzer::evaluate: error in entries counting (" << iEntry << ", " << fData->getTree()->GetEntries() << ")" << endl;
        }
    }
    else
    {
        cout << "VDispTableAnalyzer::evaluate: error: no data  tree for disp analysis" << endl;
    }
    
    cout << "VDispTableAnalyzer::evaluate: error, no results for this event (" << iSize << ", " << iWidthOLenght << ")" << endl;
    
    return -99.;
}


void VDispTableAnalyzer::terminate()
{
    if( fFile )
    {
        fFile->Close();
    }
}


/*!
     calculate mean direction for disp method

     Problem: don't really now on which site of the image the arrival direction is.

     Choose combination with smallest RMS compared to computed mean direction

     Preliminary and todo :

     Attention analysis of large array data: use first 16 telescopes only
                                             (should be: use try combinations if first XX telescopes; then choose
					                 the direction with smalles difference to existing mean)

     input:

     x:		vector of x-coordinates
     y:		vector of y-coordinates
     cosphi:	vector of cosphi
     sinphi:	vector of sinphi
     v_disp: 	vector with disp values
     v_weight:	vector with image weights

     return values:

     xs, ys:	mean direction

*/
void VDispTableAnalyzer::calculateMeanDirection( float& xs, float& ys, vector< float > x, vector< float > y, vector< float > cosphi, vector< float > sinphi, vector< float > v_disp, vector< float > v_weight )
{
    // check that vectors are of same size
    if( x.size() != v_disp.size() || x.size() != v_weight.size() )
    {
        cout << "VDispTableAnalyzer::calculateMeanDirection error: size error " << x.size() << "\t" << v_disp.size() << "\t" << v_weight.size() << endl;
        xs = -9999.;
        ys = -9999.;
        return;
    }
    
    // use first NTOT_MAX telescopes for event reconstruction only
    // (number of possible combinations is 2^NTOT_MAX )
    const unsigned int NTOT_MAX = 16;
    unsigned int iNTel_max = x.size();
    if( x.size() > NTOT_MAX )
    {
        iNTel_max = NTOT_MAX;
    }
    
    // prepare bit mask to go through all possible combinations of sign for disp calculation
    vector< bitset< NTOT_MAX > > iComb;
    for( unsigned int i = 0; i < TMath::Power( 2, ( Int_t )iNTel_max ); i++ )
    {
        iComb.push_back( i );
    }
    
    ///////////////////////////////////////////////////
    // calculate RMS of direction
    ///////////////////////////////////////////////////
    
    float iSign = -1.;
    vector< float > x_mean( iComb.size(), 0. );
    vector< float > y_mean( iComb.size(), 0. );
    vector< float > t_weight( iComb.size(), 0. );
    vector< float > t_RMS( iComb.size(), 0. );
    
    // loop over all possible binary combinations
    for( unsigned int i = 0; i < iComb.size(); i++ )
    {
        // loop over data sample
        for( unsigned int k = 0; k < iNTel_max; k++ )
        {
            // calculate mean direction for this set of signs
            if( iComb[i][k] )
            {
                x_mean[i]   += ( x[k] + v_disp[k] * cosphi[k] ) * v_weight[k];
                y_mean[i]   += ( y[k] + v_disp[k] * sinphi[k] ) * v_weight[k];
            }
            else
            {
                x_mean[i]   += ( x[k] - v_disp[k] * cosphi[k] ) * v_weight[k];
                y_mean[i]   += ( y[k] - v_disp[k] * sinphi[k] ) * v_weight[k];
            }
            t_weight[i] += v_weight[k];
        }
        if( t_weight[i] > 0. )
        {
            x_mean[i] /= t_weight[i];
            y_mean[i] /= t_weight[i];
        }
        else
        {
            x_mean[i] = -9999.;
            y_mean[i] = -9999.;
        }
    }
    
    // calculate RMS distance between mean value and each coordinate
    for( unsigned int i = 0; i < x_mean.size(); i++ )
    {
        for( unsigned int k = 0; k < iNTel_max; k++ )
        {
            if( iComb[i][k] )
            {
                iSign = -1.;
            }
            else
            {
                iSign =  1.;
            }
            
            if( v_weight[k] > 0. )
            {
                t_RMS[i] += ( ( x[k] - iSign * v_disp[k] * cosphi[k] - x_mean[i] ) * ( x[k] - iSign * v_disp[k] * cosphi[k] - x_mean[i] ) ) / v_weight[k];
                t_RMS[i] += ( ( y[k] - iSign * v_disp[k] * sinphi[k] - y_mean[i] ) * ( y[k] - iSign * v_disp[k] * sinphi[k] - y_mean[i] ) ) / v_weight[k];
            }
        }
        t_RMS[i] *= t_weight[i];
    }
    
    // get values with smallest RMS -> return value
    float t_RMS_min = 1.e9;
    int   t_RMS_bin = 0;
    for( unsigned int i = 0; i < t_RMS.size(); i++ )
    {
        if( t_RMS[i] < t_RMS_min )
        {
            t_RMS_min = t_RMS[i];
            t_RMS_bin = i;
        }
    }
    x_disp.clear();
    y_disp.clear();
    for( unsigned int k = 0; k < iNTel_max; k++ )
    {
        if( iComb[t_RMS_bin][k] )
        {
            iSign = -1.;
        }
        else
        {
            iSign =  1.;
        }
        x_disp.push_back( x[k] - iSign * v_disp[k] * cosphi[k] );
        y_disp.push_back( y[k] - iSign * v_disp[k] * sinphi[k] );
    }
    
    xs = x_mean[t_RMS_bin];
    ys = y_mean[t_RMS_bin];
    
    // direction calculation is done now if number of images in this event is smaller than NTOT_MAX
    if( x.size() < NTOT_MAX )
    {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // NTEL > NTOT_MAX
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////
    float ts = t_weight[t_RMS_bin];
    
    // use remaining telescope in direction reconstruction
    float x1 = 0.;
    float x2 = 0.;
    float y1 = 0.;
    float y2 = 0.;
    for( unsigned int k = NTOT_MAX; k < x.size(); k++ )
    {
        x1 = x[k] - v_disp[k] * cosphi[k];
        y1 = y[k] - v_disp[k] * sinphi[k];
        x2 = x[k] + v_disp[k] * cosphi[k];
        y2 = y[k] + v_disp[k] * sinphi[k];
        // check which set of (x,y) points are closer to direction calculated earlier
        if( sqrt( ( xs - x1 ) * ( xs - x1 ) + ( ys - y1 ) * ( ys - y1 ) ) < sqrt( ( xs - x2 ) * ( xs - x2 ) + ( ys - y2 ) * ( ys - y2 ) ) )
        {
            xs += x1 * v_weight[k];
            ys += y1 * v_weight[k];
        }
        else
        {
            xs += x2 * v_weight[k];
            ys += y2 * v_weight[k];
        }
        ts += t_weight[k];
    }
    if( ts > 0. )
    {
        xs /= ts;
        ys /= ts;
    }
    else
    {
        xs = -9999.;
        ys = -9999.;
    }
}

