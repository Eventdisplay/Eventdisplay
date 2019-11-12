// VCTASensitivityRequirements sensitivity requirements for CTA
//
// code from JH (mail 2013/08/31)
//
// values are hardwired

#ifndef VCTASensitivityRequirements_H
#define VCTASensitivityRequirements_H

#include <cmath>

class VCTASensitivityRequirements
{
    private:
    
        // TeV to erg conversion factor:
        static double sce()
        {
            return 1.6022;
        }
        // m^-2 to cm^-2 conversion factor:
        static double sca()
        {
            return 1.e-4;
        }
        // Total conversion factor
        static double sc()
        {
            return sce() * sca();
        }
        
    public:
    
        VCTASensitivityRequirements();
        
        static double Crab_Unit( double E )
        {
            return 2.79e-7 * pow( E, -2.57 );    // [1/(m^2 s TeV)]
        }
        static double cu( double x )
        {
            return Crab_Unit( x );
        }
        static double ergs( double E )
        {
            return sce() * E;
        }
        
        // South, 50 h:
        static double f50( double x )
        {
            return 10.00 * pow( x / 1600., 1.50 ) +
                   0.18 * pow( 0.01 / x, 1.3 ) +
                   0.0028 * pow( x, 0.30 ) +
                   0.08 * pow( 0.021 / x, 6. ) +
                   0.04 * pow( 0.016 / x, 14. ) +
                   0.02 * pow( 0.016 / x, 18. );
        }
        static double fsp50( double x ) /* Note: required only from 20 GeV to 300 TeV */
        {
            return 1.25 * f50( x * 0.94 );
        }
        static double Flux_req50_south( double E /* TeV */ )
        {
            return fsp50( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req50_E2erg_south( double E /* TeV */ )
        {
            return fsp50( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_req50_CU_south( double E /* TeV */ )
        {
            return fsp50( E ); // [Crab Units]
        }
        
        // North, 50 h:
        static double fn50( double x )
        {
            return 150.00 * pow( x / 1600., 1.74 ) +
                   0.18 * pow( 0.01 / x, 1.3 ) +
                   0.0028 * pow( x, 0.30 ) +
                   0.08 * pow( 0.021 / x, 6. ) +
                   0.04 * pow( 0.016 / x, 14. ) +
                   0.02 * pow( 0.016 / x, 18. );
        }
        static double fnsp50( double x ) /* Note: required only from 20 GeV to 20 TeV */
        {
            return 1.25 * fn50( x * 0.94 );
        }
        static double Flux_req50_north( double E /* TeV */ )
        {
            return fnsp50( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req50_E2erg_north( double E /* TeV */ )
        {
            return fnsp50( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_req50_CU_north( double E /* TeV */ )
        {
            return fnsp50( E ); // [Crab Units]
        }
        
        // South, 5h:
        static double f5( double x )
        {
            return 100.0 * pow( x / 1500., 1.54 ) +
                   0.50 * pow( 0.01 / x, 1.3 ) +
                   0.0120 * pow( x, 0.25 ) +
                   0.25 * pow( 0.02 / x, 6. ) +
                   0.06 * pow( 0.016 / x, 14. ) +
                   0.03 * pow( 0.016 / x, 18. );
        }
        static double fsp5( double x )
        {
            return 1.25 * f5( x * 0.94 );
        }
        static double Flux_req5_south( double E /* TeV */ )
        {
            return fsp5( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req5_E2erg_south( double E /* TeV */ )
        {
            return fsp5( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_req5_CU_south( double E /* TeV */ )
        {
            return fsp5( E ); // [Crab Units]
        }
        
        // North, 5h:
        static double fn5( double x )
        {
            return 1300.0 * pow( x / 1500., 1.71 ) +
                   0.50 * pow( 0.01 / x, 1.3 ) +
                   0.0120 * pow( x, 0.25 ) +
                   0.25 * pow( 0.02 / x, 6 ) +
                   0.06 * pow( 0.016 / x, 14. ) +
                   0.03 * pow( 0.016 / x, 18. );
        }
        static double fnsp5( double x )
        {
            return 1.25 * fn5( x * 0.94 );
        }
        static double Flux_req5_north( double E /* TeV */ )
        {
            return fnsp5( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req5_E2erg_north( double E /* TeV */ )
        {
            return fnsp5( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_req5_CU_north( double E /* TeV */ )
        {
            return fnsp5( E ); // [Crab Units]
        }
        
        // South, 0.5 h:
        static double f05( double x )
        {
            return 1000.*pow( x / 1400., 1.58 ) +
                   1.30 * pow( 0.01 / x, 1.3 ) +
                   0.0600 * pow( x, 0.20 ) +
                   0.50 * pow( 0.02 / x, 6. ) +
                   0.07 * pow( 0.016 / x, 14. ) +
                   0.035 * pow( 0.016 / x, 18. );
        }
        static double fsp05( double x )
        {
            return 1.25 * f05( x * 0.94 );
        }
        static double Flux_req05_south( double E /* TeV */ )
        {
            return fsp05( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req05_E2erg_south( double E /* TeV */ )
        {
            return fsp05( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        
        static double Flux_req05_CU_south( double E /* TeV */ )
        {
            return fsp05( E ); // [Crab Units]
        }
        
        // North, 0.5 h:
        static double fn05( double x )
        {
            return 11000.*pow( x / 1400., 1.68 ) +
                   1.30 * pow( 0.01 / x, 1.3 ) +
                   0.0600 * pow( x, 0.20 ) +
                   0.50 * pow( 0.02 / x, 6. ) +
                   0.07 * pow( 0.016 / x, 14. ) +
                   0.035 * pow( 0.016 / x, 18. );
        }
        static double fnsp05( double x )
        {
            return 1.25 * fn05( x * 0.94 );
        }
        static double Flux_req05_north( double E /* TeV */ )
        {
            return fnsp05( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_req05_E2erg_north( double E /* TeV */ )
        {
            return fnsp05( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_req05_CU_north( double E /* TeV */ )
        {
            return fnsp05( E ); // [Crab Units]
        }
        
        // Goal sensitivities, only derived for 50 hours observation time:
        
        // Goal South, 50 h:
        static double fd50( double x )
        {
            return 6.50 * pow( x / 1600., 1.50 ) +
                   0.06 * pow( 0.01 / x, 1.3 ) +
                   0.0013 * pow( x, 0.30 ) +
                   0.05 * pow( 0.021 / x, 6. ) +
                   0.022 * pow( 0.016 / x, 14. ) +
                   0.02 * pow( 0.016 / x, 18. );
        }
        static double fdes50( double x )
        {
            return 1.25 * fd50( x * 0.94 );
        }
        static double Flux_goal50_south( double E /* TeV */ )
        {
            return fdes50( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_goal50_E2erg_south( double E /* TeV */ )
        {
            return fdes50( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_goal50_CU_south( double E /* TeV */ )
        {
            return fdes50( E ); // [Crab Units]
        }
        
        // Goal North, 50 h:
        static double fnd50( double x )
        {
            return 100.0 * pow( x / 1600., 1.74 ) +
                   0.06 * pow( 0.01 / x, 1.3 ) +
                   0.0016 * pow( x, 0.30 ) +
                   0.05 * pow( 0.021 / x, 6. ) +
                   0.022 * pow( 0.016 / x, 14. ) +
                   0.02 * pow( 0.016 / x, 18. );
        }
        static double fndes50( double x )
        {
            return 1.25 * fnd50( x * 0.94 );
        }
        static double Flux_goal50_north( double E /* TeV */ )
        {
            return fndes50( E ) * cu( E ); // [1/(m^2 s TeV)]
        }
        static double Flux_goal50_E2erg_north( double E /* TeV */ )
        {
            return fndes50( E ) * cu( E ) * E * E * sc(); // [erg/(cm^2 s)]
        }
        static double Flux_goal50_CU_north( double E /* TeV */ )
        {
            return fndes50( E ); // [Crab Units]
        }
        
};

#endif
