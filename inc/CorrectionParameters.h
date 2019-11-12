//-*-mode:c++; mode:font-lock;-*-

/**
 * \file CorrectionParameters.h
 * \ingroup SEphem
 * \brief This is a one-line description of this cpp file.
 *
 * Here is a tedious verbose multi-line description of all
 * the details of the code, more than you would
 * ever want to read. Generally, all the important documentation
 * goes in the .cpp files.
 *
 * Original Author: Stephen Fegan
 *
 **/

#ifndef SEPHEM_CORRECTIONPARAMETERS_H
#define SEPHEM_CORRECTIONPARAMETERS_H

#include"Angle.h"
#include"SphericalCoords.h"

namespace SEphem
{
    class CorrectionParameters
    {
        public:
            CorrectionParameters():
                enable_offsets( true ), enable_corrections( true ), enable_vff( true ),
                az_ratio( 1.0 ), el_ratio( 1.0 ), az_offset( 0.0 ), el_offset( 0.0 ),
                az_ns( 0.0 ), az_ew( 0.0 ), el_udew( 0.0 ), fp_az( 0.0 ),
                flex_el_A( 0.0 ), flex_el_B( 0.0 ),
                el_pos_vff_s( 0.0 ), el_pos_vff_t( 0.0 ), el_neg_vff_s( 0.0 ), el_neg_vff_t( 0 ),
                az_pos_vff_s( 0.0 ), az_pos_vff_t( 0.0 ), az_neg_vff_s( 0.0 ), az_neg_vff_t( 0 )
            {}
            
            bool doAzElCorrections( double& az_driveangle, double& el_driveangle,
                                    const double& tel_az_driveangle,
                                    bool do_corrections ) const;
            void undoAzElCorrections( double& az_driveangle, double& el_driveangle,
                                      bool do_corrections ) const;
                                      
            bool save( const char* filename ) const;
            bool load( const char* filename );
            
            static std::string loadFilename( unsigned scope_num );
            static std::string saveFilename( unsigned scope_num );
            
            // ------------------------------------------------------------------------
            // CORRECTION PARAMETERS
            // ------------------------------------------------------------------------
            
            bool enable_offsets;
            bool enable_corrections;
            bool enable_vff;
            double az_ratio;                      // real_rad / Az_encoder_rad                    [1.0]
            double el_ratio;                      // real_rad / El_encoder_rad                    [1.0]
            double az_offset;                     // Az encoder offset                 [Az_encoder_rad]
            double el_offset;                     // El encoder offset                 [El_encoder_rad]
            double az_ns;                         // NS az table angle                            [rad]
            double az_ew;                         // EW az table angle                            [rad]
            double el_udew;                       // UD-EW (in stow position) el mis-alignment    [rad]
            double fp_az;                         // NS-EW (in stow position) focus mis-alignment [rad]
            double flex_el_A;                     // Coefficiant of cos(elevation) flexure        [rad]
            double flex_el_B;                     // Coefficiant of cos(2 * elevation) flexure    [rad]
            double el_pos_vff_s;                  // Elevation positive VFF slope                   [s]
            double el_pos_vff_t;                  // Elevation positive VFF threshold           [rad/s]
            double el_neg_vff_s;                  // Elevation negative VFF slope                   [s]
            double el_neg_vff_t;                  // Elevation negative VFF threshold           [rad/s]
            double az_pos_vff_s;                  // Azimuth positive VFF slope                     [s]
            double az_pos_vff_t;                  // Azimuth positive VFF intercept             [rad/s]
            double az_neg_vff_s;                  // Azimuth negative VFF slope                     [s]
            double az_neg_vff_t;                  // Azimuth negative VFF intercept             [rad/s]
            
        private:
            static const double sc_lim_az_cw;
            static const double sc_lim_az_cc;
            
            static const double sc_inversion_tol;
            static const int    sc_inversion_it_max;
            /*            static const double sc_lim_az_cw =     270.0 * ANGLE_RADPERDEG;
                        static const double sc_lim_az_cc =    -270.0 * ANGLE_RADPERDEG;
            
                        static const double sc_inversion_tol = 0.0001 * ANGLE_RADPERDEG;
                        static const int    sc_inversion_it_max = 20; */
    };
    
}                                                 // namespace SEphem
#endif                                            // SEPHEM_CORRECTIONPARAMETERS_H
