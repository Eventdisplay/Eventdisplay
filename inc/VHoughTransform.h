//Hough transform muon id header file

#ifndef VHoughTransform_H
#define VHoughTransform_H

//Standard includes

#include <iostream>
#include <stdio.h>
#include <stdlib.h>


//ROOT includes

#include "TH3D.h"
#include "TTree.h"

//Event Display includes

#include "VDetectorGeometry.h"
#include "VEvndispData.h"

using namespace std;

class VHoughTransform
{


    public:
    
        VHoughTransform( VDetectorGeometry* fDetectorGeo ); //Initializes the hough transform class and data structures
        
        void analysis( VEvndispData* fData, VImageParameter* fParGeo ); //Performs the Hough transform analysis
        
        int getNpixMin( unsigned int fTelID ); //Get the minimum number of hit pixels cut
        
        int getNpixMax( unsigned int fTelID ); //Get the maximum number of hit pixels cut
        
    private:
    
        VDetectorGeometry* fDetectorGeometry;// Pointer to the detector geometry
        
        int fNumberOfMuons; //Number of detected muons counter
        
        vector <int> fNumberOfChannels; //Number of channels in the VERITAS cameras
        
        vector <int> fRMinDpmt;//Minimum radius in units of PMT diameters of the template circles
        
        vector <int> fRMaxDpmt;//Maximum radius in units of PMT diameters of the template circles
        
        vector <int> fStepsPerPMTDiameter;//Number of steps in radius per PMT diameter
        
        vector <int> fNpixMin;//Min number of pixels cut
        
        vector <int> fNpixMax;//Max number of pixels cut
        
        vector <double> fAlpha;//Hough transform cut, described in memo and cuts file
        
        vector <double> fBeta;//Hough transform cut, described in memo and cuts file
        
        vector <double> fGamma;//Hough transform cut, described in memo and cuts file
        
        vector <double> fEta;//Azimuthal completeness cut
        
        vector <double> fCameraRadius;//Containedness cut (camera radius in units of PMT diameters)
        
        vector <double> fPMTDiameter; //Diameter of a PMT in mm
        
        vector <TH3D*> fAccumulatorArray; //Hough transform accumulator arrays
        
        vector <TTree*> fHTLookupTableTree; //Trees with the Hough transform lookup tables
        
        TH3D* initAccumulatorArray( int fRMinDpmt, int fRMaxDpmt, int fStepsPerPMTDiameter, unsigned int fTelID ); //Method for initializaing the Hough transform accumulator array
        
        TTree* initLookupTable( int fRMinDpmt, int fRMaxDpmt, int fStepsPerPMTDiameter, unsigned int fTelID ); //Method for initializing the Hough transform lookup table
        
        void readHTParameterFile( unsigned int fTelID );
        
};

#endif
