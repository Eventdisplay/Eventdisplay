/*! \class  VPointingCorrectionsTreeReader
 *  \brief reader class for pointing corrections
 *
 *  recalculates cen_x, cen_y, and phi image parameter
 *  after pointing corrections
 *
 */

#include "VPointingCorrectionsTreeReader.h"

VPointingCorrectionsTreeReader::VPointingCorrectionsTreeReader( TChain *t )
{
    fPointingErrorX = 0.;
    fPointingErrorY = 0.;
    fPointingCorrectionTreeSetting = false;

    fTree = t;
    if( fTree )
    {
        if( fTree->GetBranchStatus( "PointingErrorX" ) &&
            fTree->GetBranchStatus( "PointingErrorY" ) )
        {
            fTree->SetBranchAddress( "PointingErrorX", &fPointingErrorX );
            fTree->SetBranchAddress( "PointingErrorY", &fPointingErrorY );
            fPointingCorrectionTreeSetting = true;
        }
    }

}


int VPointingCorrectionsTreeReader::getEntry( Long64_t iEntry )
{
    if( !fTree )
    {
        fPointingErrorX = 0.;
        fPointingErrorY = 0.;
        return 0;
    }
    return fTree->GetEntry( iEntry );
}


Long64_t VPointingCorrectionsTreeReader::getEntries()
{
    if( !fTree )
    {
        fPointingErrorX = 0.;
        fPointingErrorY = 0.;
        return 0;
    }
    
    return fTree->GetEntries();
}

/*
 * pointing corrected cen_x
 */
float VPointingCorrectionsTreeReader::getCorrected_cen_x( float cen_x )
{
   return cen_x + fPointingErrorX;
}

/*
 * pointing corrected cen_y
 */
float VPointingCorrectionsTreeReader::getCorrected_cen_y( float cen_y )
{
   return cen_y + fPointingErrorY;
}


/*
    add a small pointing error to centroids and recalculate image angle phi
*/
float VPointingCorrectionsTreeReader::getCorrected_phi( float cen_x, float cen_y, float d, float s, float sdevxy )
{
    float xmean = cen_x + fPointingErrorX;
    float ymean = cen_y + fPointingErrorY;

    const double ac = ( d + s ) * ymean + 2.0 * sdevxy * xmean;
    const double bc = 2.0 * sdevxy * ymean - ( d - s ) * xmean;
    const double cc = sqrt( ac * ac + bc * bc );

    return atan2( ac / cc, bc / cc );
}
