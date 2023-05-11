/*
         calculate the geomagnetic field compontents at the CTA site candidates
*/

double orthogonal_H_Z( double zenith_deg, double azimuth_deg, double B_z, double B_h )
{

    double x = -sin( zenith_deg * TMath::DegToRad() ) * sin( azimuth_deg * TMath::DegToRad() ) * B_z;
    double y =  sin( zenith_deg * TMath::DegToRad() ) * cos( azimuth_deg * TMath::DegToRad() ) * B_z
                - cos( zenith_deg * TMath::DegToRad() ) * B_h;
    double z = sin( zenith_deg * TMath::DegToRad() ) * sin( azimuth_deg * TMath::DegToRad() ) * B_h;
    
    return sqrt( x * x + y * y + z * z );
}

double orthogonal_X_Y_Z( double zenith_deg, double azimuth_deg, double B_z, double B_x, double B_y )
{

    double x = cos( zenith_deg * TMath::DegToRad() ) * B_y
               - sin( zenith_deg * TMath::DegToRad() ) * sin( azimuth_deg * TMath::DegToRad() ) * B_z;
    double y =  sin( zenith_deg * TMath::DegToRad() ) * cos( azimuth_deg * TMath::DegToRad() ) * B_z
                - cos( zenith_deg * TMath::DegToRad() ) * B_x;
    double z = sin( zenith_deg * TMath::DegToRad() ) * sin( azimuth_deg * TMath::DegToRad() ) * B_x
               - sin( zenith_deg * TMath::DegToRad() ) * cos( azimuth_deg * TMath::DegToRad() ) * B_y;
               
    return sqrt( x * x + y * y + z * z );
}

/*
     all values from the prod2 CORSIKA simulation header

*/
void printBField_at_CTA_sites( bool bCTA = true )
{
    vector< double > fZe;
    vector< double > fAz;
    fZe.push_back( 5. );
    fZe.push_back( 10. );
    fZe.push_back( 15. );
    fZe.push_back( 20. );
    fZe.push_back( 45. );
    fZe.push_back( 50. );
    fZe.push_back( 60. );
    fAz.push_back( 0. );
    fAz.push_back( 180. );
    
    vector< string > fSiteName;
    vector< double > fB_tot;
    vector< double > fInclination_deg;
    vector< double > fDeclination_deg;
    
    if( bCTA )
    {
    
        /*	fSiteName.push_back( "Leoncito" );
            fB_tot.push_back( 23.5536 );
            fInclination_deg.push_back( -32.3406 );
            fDeclination_deg.push_back( 0. );
        
            fSiteName.push_back( "SAC" );
            fB_tot.push_back( 22.7161 );
            fInclination_deg.push_back( -23.0661 );
            fDeclination_deg.push_back( 0. ); */
        
        fSiteName.push_back( "Aar" );
        fB_tot.push_back( 27.1812 );
        fInclination_deg.push_back( -66.3586 );
        fDeclination_deg.push_back( 0. );
        
        /*	fSiteName.push_back( "Armazones" );
            fB_tot.push_back( 23.1594 );
            fInclination_deg.push_back( -22.6642 );
            fDeclination_deg.push_back( 0. );
        
        fSiteName.push_back( "Tenerife" );
        fB_tot.push_back( 38.5314 );
        fInclination_deg.push_back( 36.952 );
        fDeclination_deg.push_back( 0. ); */
        
        /*	fSiteName.push_back( "SPM" );
            fB_tot.push_back( 45.8194 );
            fInclination_deg.push_back( 56.6131 );
            fDeclination_deg.push_back( 0. );
        
            fSiteName.push_back( "US" );
            fB_tot.push_back( 48.9437 );
            fInclination_deg.push_back( 61.2665 ); */
    }
    else
    {
        fSiteName.push_back( "VERITAS" );
        fB_tot.push_back( 48.0231 );
        fInclination_deg.push_back( 58.3487 );
        fDeclination_deg.push_back( 0. );
    }
    
    for( unsigned int i = 0; i < fSiteName.size(); i++ )
    {
        cout << "=========================" << endl;
        cout << fSiteName[i] << endl;
        for( unsigned int j = 0; j < fZe.size(); j++ )
        {
            for( unsigned int k = 0; k < fAz.size(); k++ )
            {
                cout << "\t Ze=" << fZe[j] << " deg, Az=" << fAz[k] << " deg: ";
                cout << " B_tot = " << fB_tot[i] << ", inclination " << fInclination_deg[i];
                cout << " deg, declination " << fDeclination_deg[i] << " deg, ";
                double B_h = fB_tot[i] * cos( fInclination_deg[i] * TMath::DegToRad() );
                double B_x = B_h * cos( fDeclination_deg[i] * TMath::DegToRad() );
                double B_y = B_h * sin( fDeclination_deg[i] * TMath::DegToRad() );
                double B_z = fB_tot[i] * sin( fInclination_deg[i] * TMath::DegToRad() );
                cout << "B_h = " << B_h;
                cout << ", B_z = " << B_z;
                cout << ", B_x = " << B_x;
                cout << ", B_y = " << B_y;
                if( fDeclination_deg[i] == 0. )
                {
                    cout << ", B_orth = " << orthogonal_H_Z( fZe[j], fAz[k], B_z, B_h ) << endl;
                }
                else
                {
                    cout << ", B_orth = " << orthogonal_X_Y_Z( fZe[j], fAz[k], B_z, B_x, B_y ) << endl;
                }
            }
        }
    }
    
}

