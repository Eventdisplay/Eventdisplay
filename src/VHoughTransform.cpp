//Hough transform muon id class

#include "VHoughTransform.h"


VHoughTransform::VHoughTransform( VDetectorGeometry* fDetectorGeo )
{

    //Hough transform muon ID constructor
    
    
    //Print this when the constructor is run.
    //cout << "---------------------------------------------------------------------" << endl;
    //cout << "Initializing the Hough transform muon identification algorithm." << endl;
    
    //Set the global pointer to the detector geometry
    fDetectorGeometry = fDetectorGeo;
    
    
    //Set the number of muons counter to zero
    fNumberOfMuons = 0;
    
    
    //Print the number of telescopes
    //cout << "Number of telescopes: " << fDetectorGeometry->getNTel() << endl;
    
    
    //Loop over all the telescopes and set up the accumulator arrays, lookup tables, PMT diameters and number of channels
    for( unsigned int iTelescopeIndex = 0 ; iTelescopeIndex < fDetectorGeometry->getNTel() ; iTelescopeIndex++ )
    {
    
        //cout << endl;
        //cout << "Telescope " << iTelescopeIndex + 1 << ":" << endl;
        
        //Read in the parameter file values for each telescope.
        //Stores the values in global variables.
        readHTParameterFile( iTelescopeIndex );
        
        //Calculate the distance between the first and second pixels to get the pixel diameter
        //Assumes circular PMTs
        fPMTDiameter.push_back( sqrt(
        
                                    ( fDetectorGeometry->getX_MM( iTelescopeIndex )[0] -
                                      fDetectorGeometry->getX_MM( iTelescopeIndex )[1] ) *
                                    ( fDetectorGeometry->getX_MM( iTelescopeIndex )[0] -
                                      fDetectorGeometry->getX_MM( iTelescopeIndex )[1] ) +
                                      
                                    ( fDetectorGeometry->getY_MM( iTelescopeIndex )[0] -
                                      fDetectorGeometry->getY_MM( iTelescopeIndex )[1] ) *
                                    ( fDetectorGeometry->getY_MM( iTelescopeIndex )[0] -
                                      fDetectorGeometry->getY_MM( iTelescopeIndex )[1] )
                                      
                                ) );//End of PMT diameter calculation
                                
                                
        //Get the number of channels for each telescope.
        fNumberOfChannels.push_back( fDetectorGeometry->getNChannels( iTelescopeIndex ) );
        
        //Print this before initializing the accumulator array for a given telescope
        //cout << "Initializing the accumulator array for telescope " << iTelescopeIndex + 1 << "..." << endl;
        
        //Set up the accumulator array for a given telescope
        fAccumulatorArray.push_back( initAccumulatorArray( fRMinDpmt[iTelescopeIndex], fRMaxDpmt[iTelescopeIndex],
                                     fStepsPerPMTDiameter[iTelescopeIndex], iTelescopeIndex ) );
                                     
        //Print this when the accumulator array is initialized.
        //cout << "Accumulator array for telescope " << iTelescopeIndex + 1 << " initialized." << endl;
        
        //Print this before initializing the lookup table
        //cout << "Initializing the lookup table for telescope " << iTelescopeIndex + 1 << "..." << endl;
        
        //Set up the lookup table for a given telescope
        fHTLookupTableTree.push_back( initLookupTable( fRMinDpmt[iTelescopeIndex], fRMaxDpmt[iTelescopeIndex],
                                      fStepsPerPMTDiameter[iTelescopeIndex], iTelescopeIndex ) );
                                      
        //Print this when the lookup table is initialized
        //cout << "Lookup table for telescope " << iTelescopeIndex + 1 << " initialized." << endl;
        
        
    }//End of loop over all the telescopes and set up the accumulator arrays and lookup tables
    
    
    //Output this when everything has been initialized
    //cout << endl;
    cout << "Hough transform muon identification algorithm initialized." << endl;
    cout << "" << endl;
    //cout << "---------------------------------------------------------------------" << endl;
    //cout << endl;
    
    
}//End of constructor


void VHoughTransform::analysis( VEvndispData* fData, VImageParameter* fParGeo )
{

    //Hough transform analysis code
    
    int fIsMuon = 99; //Set to 1 if muon, 0 if non-muon.
    
    int fNpix = 0; //Number of hit pixels
    
    double fPixelXCoordinate = 0; //The X coordinate of a pixel
    double fPixelYCoordinate = 0; //The Y coordinate of a pixel
    
    char fBranchName[100]; //Name of the branch in the lookup table tree
    
    double fSumOfAllBins = 0; //Sum of all the bins in the accumulator array
    
    int fNumberOfNonZeroBins = 0; //Number of non-zero bins in the accumulator array
    
    double fCircleCoordinates[3]; 	//Coordinates of a template circle associated with a particular pixel
    fCircleCoordinates[0] = 0; 		//x coordinate
    fCircleCoordinates[1] = 0; 		//y coordinate
    fCircleCoordinates[2] = 0; 		//r coordinate.
    
    //Best parameterized circles
    
    double fBestParametrization[3];  		//Best parametrization
    fBestParametrization[0] = 0; 			//x coordinate
    fBestParametrization[1] = 0;			//y coordinate
    fBestParametrization[2] = 0;			//r coordinate
    
    double fSecondBestParametrization[3]; 	//Second best parametrization
    fSecondBestParametrization[0] = 0;		//x coordinate
    fSecondBestParametrization[1] = 0;		//y coordinate
    fSecondBestParametrization[2] = 0;		//r coordinate
    
    double fThirdBestParametrization[3]; 	//Third best parametrization
    fThirdBestParametrization[0] = 0;		//x coordinate
    fThirdBestParametrization[1] = 0;		//y coordinate
    fThirdBestParametrization[2] = 0;		//r coordinate
    
    int fAccumulatorBins[3]; //Bins of the accumulator array
    
    double fMaxBinValue = 0; //Value of the bin of the accumulator array with the highest value
    double fSecondMaxBinValue = 0; //Value of the bin of the accumulator array with the second highest value
    
    int fMaxBin = 0; //Bin number of the max bin of the accumulator array
    int fSecondMaxBin = 0; //Bin number of the second max bin of the accumulator array
    
    double fDistance1 = 0; //Hyper-distance between the best and second best parametrizations
    double fDistance2 = 0; //Hyper-distance between the best and third best parametrizations
    double fDistance3 = 0; //Hyper-distance between the second best and third best parametrizations
    double fTD = 0; //TD variable. (Sum of the hyperdistances)
    
    double fAverageNonZeroBinContent = 0; //The average non-zero bin content of the accumulator array
    double fAP = 0; //AP variable. (Max bin content of the accumulator array divided by the average non-zero bin content)
    
    int fNPR = 0; //NPR varialbe. (Number of non-zero pixels hit by any of the three best parametrizations)
    
    double fNPRCentroidX = 0; //x coordinate of the centroid of the pixels included in the NPR count.
    double fNPRCentroidY = 0; //y coordinate of the centroid of the pixels included in the NPR count.
    
    double fCD = 0;//CD variable. (The distance from the centre of the best parametrization to the centroid of the pixels included in the NPR count).
    
    double fCN = 0; //The C/N varaible. (CD/NPR). Used for determining azimuthal completeness.
    
    double fContained = 0; // Distance from the center of the ring to the center of the camera plus the ring radius in mm
    
    //Reset the accumulator array
    fAccumulatorArray[ fData->getTelID() ]->Reset();
    
    
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ fData->getTelID() ] ; iChannelIndex++ ) // Loop over all the pixels
    {
    
    
        if( fData->getImage()[iChannelIndex] || fData->getBorder()[iChannelIndex] )//If the pixel is an image or border pixel (binary image).
        {
        
            fNpix++;//Increment the hit pixel count
            
            //Fill the Accumulator array here.
            
            sprintf( fBranchName, "Pixel %d", iChannelIndex ); //Set the branch name for ith channel
            
            TBranch* fBranch = fHTLookupTableTree[ fData->getTelID() ]->GetBranch( fBranchName ); //Grab the branch for the ith channel
            
            fBranch->SetAddress( &fCircleCoordinates ); //Map the arrays in the branch to circle coordinates array
            
            
            //Loop over all circle coordinates for that pixel and fill the appropriate bins of the accumulator array
            
            //Set fCircleCoordinates[2] to a non-zero value to get into for the loop. The last entries in the lookup table are ( 0, 0, 0)
            fCircleCoordinates[2] = 99;
            
            for( int iCircleParametrizationIndex = 0 ; fCircleCoordinates[2] != 0 ; iCircleParametrizationIndex++ ) //Loop over circle parametrizations until the r coordinate is zero (Last entries are 0,0,0)
            {
            
                fBranch->GetEntry( iCircleParametrizationIndex ); // Get the jth circle parametrization for that pixel. Sets the values of fCircleCoordinates
                
                if( fCircleCoordinates[2] != 0 )// If the r coordinate is not zero (end of list), do this
                {
                
                    //If bin content is zero and is filled, increment the number of non zero bins variable
                    if( ( fAccumulatorArray[ fData->getTelID() ]->GetBinContent( fAccumulatorArray[ fData->getTelID() ]->FindBin( fCircleCoordinates[0], fCircleCoordinates[1], fCircleCoordinates[2] ) ) ) == 0 )
                    {
                    
                        fNumberOfNonZeroBins++;
                        
                    }
                    
                    //Fill the appropriate bin of accumulator array with 1.0 (Binary image). For images with charge values, add the charge instead of 1.0.
                    fAccumulatorArray[ fData->getTelID() ]->Fill( fCircleCoordinates[0], fCircleCoordinates[1], fCircleCoordinates[2], 1.0 );
                    
                    //Add 1.0 to the sum of all bins variable.
                    fSumOfAllBins = fSumOfAllBins + 1.0;
                    
                }//End of (if r coordinate is non zero) statement
                
            }//End of loop over circle parametrizations
            
        }//End of (If the pixel is a border or image pixel) statement.
        
    }// End of loop over all the pixels.
    
    
    //End of accumulator array filling.
    
    
    //Get the best circle parametrizations from the accumulator array
    
    //Get best parameterized circle
    
    fMaxBin = fAccumulatorArray[ fData->getTelID() ]->GetMaximumBin( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2] ); //Get the max bin of the accumulator array
    fBestParametrization[0] = fAccumulatorArray[ fData->getTelID() ]->GetXaxis()->GetBinCenter( fAccumulatorBins[0] ); //Get the x coordinate of the max bin
    fBestParametrization[1] = fAccumulatorArray[ fData->getTelID() ]->GetYaxis()->GetBinCenter( fAccumulatorBins[1] ); //Get the y coordinate of the max bin
    fBestParametrization[2] = fAccumulatorArray[ fData->getTelID() ]->GetZaxis()->GetBinCenter( fAccumulatorBins[2] ); //Get the r coordinate of the max bin
    fMaxBinValue = fAccumulatorArray[ fData->getTelID() ]->GetBinContent( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2] ); //Get the value of the max bin
    
    
    //Get second best parametrized circle
    
    //Set the value of the max bin of the accumulator array to zero. Do this to use the GetMaximumBin() method.
    fAccumulatorArray[ fData->getTelID() ]->SetBinContent( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2], 0 );
    
    fSecondMaxBin = fAccumulatorArray[ fData->getTelID() ]->GetMaximumBin( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2] ); //Get the second max bin of the accumulator array
    fSecondBestParametrization[0] = fAccumulatorArray[ fData->getTelID() ]->GetXaxis()->GetBinCenter( fAccumulatorBins[0] ); //Get the x coordinate of the second max bin
    fSecondBestParametrization[1] = fAccumulatorArray[ fData->getTelID() ]->GetYaxis()->GetBinCenter( fAccumulatorBins[1] ); //Get the y coordinate of the second max bin
    fSecondBestParametrization[2] = fAccumulatorArray[ fData->getTelID() ]->GetZaxis()->GetBinCenter( fAccumulatorBins[2] ); //Get the r coordinate of the second max bin
    fSecondMaxBinValue = fAccumulatorArray[ fData->getTelID() ]->GetBinContent( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2] ); //Get the value of the second max bin
    
    
    //Get the third best parametrized circle
    
    //Set the value of the second max bin of the accumulator array to zero. Do this to use the GetMaximumBin() method.
    fAccumulatorArray[ fData->getTelID() ]->SetBinContent( fAccumulatorBins[0], fAccumulatorBins[1], fAccumulatorBins[2], 0 );
    
    fThirdBestParametrization[0] = fAccumulatorArray[ fData->getTelID() ]->GetXaxis()->GetBinCenter( fAccumulatorBins[0] ); //Get the x coordinate of the third max bin
    fThirdBestParametrization[1] = fAccumulatorArray[ fData->getTelID() ]->GetYaxis()->GetBinCenter( fAccumulatorBins[1] ); //Get the x coordinate of the third max bin
    fThirdBestParametrization[2] = fAccumulatorArray[ fData->getTelID() ]->GetZaxis()->GetBinCenter( fAccumulatorBins[2] ); //Get the x coordinate of the third max bin
    
    //Refill the bins of the accumulator that were set to zero.
    fAccumulatorArray[ fData->getTelID() ]->SetBinContent( fMaxBin, fMaxBinValue ); //Replace the bin content of the max bin
    fAccumulatorArray[ fData->getTelID() ]->SetBinContent( fSecondMaxBin, fSecondMaxBinValue ); //Replace the bin content of the second max bin
    
    
    //Calculate discriminating variables
    
    //Calculate hyperdistances for the TD variable
    
    //Distance from the best parametrization to the second best parametrization
    fDistance1 = sqrt( ( fBestParametrization[0] - fSecondBestParametrization[0] ) * ( fBestParametrization[0] - fSecondBestParametrization[0] ) +
                       ( fBestParametrization[1] - fSecondBestParametrization[1] ) * ( fBestParametrization[1] - fSecondBestParametrization[1] ) +
                       ( fBestParametrization[2] - fSecondBestParametrization[2] ) * ( fBestParametrization[2] - fSecondBestParametrization[2] ) );
                       
    //Distance from the best parametrization to the third best parametrization
    fDistance2 = sqrt( ( fBestParametrization[0] - fThirdBestParametrization[0] ) * ( fBestParametrization[0] - fThirdBestParametrization[0] ) +
                       ( fBestParametrization[1] - fThirdBestParametrization[1] ) * ( fBestParametrization[1] - fThirdBestParametrization[1] ) +
                       ( fBestParametrization[2] - fThirdBestParametrization[2] ) * ( fBestParametrization[2] - fThirdBestParametrization[2] ) );
                       
    //Distance from the second best parametrization to the third best parametrization
    fDistance3 = sqrt( ( fSecondBestParametrization[0] - fThirdBestParametrization[0] ) * ( fSecondBestParametrization[0] - fThirdBestParametrization[0] ) +
                       ( fSecondBestParametrization[1] - fThirdBestParametrization[1] ) * ( fSecondBestParametrization[1] - fThirdBestParametrization[1] ) +
                       ( fSecondBestParametrization[2] - fThirdBestParametrization[2] ) * ( fSecondBestParametrization[2] - fThirdBestParametrization[2] ) );
                       
    //Calculate the TD variable from the hyperdistances
    fTD = fDistance1 + fDistance2 + fDistance3;
    
    //Calculate the AP variable
    
    //Calculate the average non-zero bin content of the accumulator array
    fAverageNonZeroBinContent = fSumOfAllBins / ( double ) fNumberOfNonZeroBins;
    
    //Calculate the AP variable. (Value of the maximum bin divided by the average non-zero bin content)
    fAP = fMaxBinValue / fAverageNonZeroBinContent;
    
    
    //Loop over the pixels again to calculate NPR and CDIST to calculate the C/N parameter
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ fData->getTelID() ] ; iChannelIndex++ )
    {
    
        //Find the number of pixels hit by any of the three best parameterizations
        
        //First check if the channel is non-zero
        if( fData->getImage()[iChannelIndex] || fData->getBorder()[iChannelIndex] )//If the pixel is an image or border pixel (binary image).
        {
        
            //Get the channel coordinates in mm
            fPixelXCoordinate = fDetectorGeometry->getX_MM( fData->getTelID() )[iChannelIndex]; 	//Get the x coordinate of ith channel
            
            fPixelYCoordinate = fDetectorGeometry->getY_MM( fData->getTelID() )[iChannelIndex]; 	//Get the y coordinate of ith channel
            
            //Check to see if the pixel is hit by any of the three best parametrizations
            
            if(
            
                ( // Check to see if the pixel is hit by the best parametrization
                
                    //Check if the pixel is contained within a circle with x,y coordinates of the best parametrization
                    //and radius of the best parametrization plus a PMT diameter.
                    
                    (
                    
                        sqrt( ( fBestParametrization[0] - fPixelXCoordinate ) * ( fBestParametrization[0] - fPixelXCoordinate ) 	+
                              ( fBestParametrization[1] - fPixelYCoordinate ) * ( fBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fBestParametrization[2] + fPMTDiameter[ fData->getTelID()] )
                        
                    )
                    
                    //Check if the pixel is not contained within a circle with x,y coordinates of the best parametrization
                    //and radius of the best parametrization minus a PMT diameter.
                    
                    && !(
                    
                        sqrt( ( fBestParametrization[0] - fPixelXCoordinate ) * ( fBestParametrization[0] - fPixelXCoordinate ) +
                              ( fBestParametrization[1] - fPixelYCoordinate ) * ( fBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fBestParametrization[2] - fPMTDiameter[ fData->getTelID() ] )
                        
                    )
                    
                    
                )//End of check to see if the pixel is hit by the best parametrization
                
                
                || //Or
                
                
                ( // Check to see if the pixel is hit by the second best parametrization
                
                
                    //Check if the pixel is contained within a circle with x,y coordinates of the second best parametrization
                    //and radius of the second best parametrization plus a PMT diameter.
                    
                    (
                    
                        sqrt( ( fSecondBestParametrization[0] - fPixelXCoordinate ) * ( fSecondBestParametrization[0] - fPixelXCoordinate ) +
                              ( fSecondBestParametrization[1] - fPixelYCoordinate ) * ( fSecondBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fSecondBestParametrization[2] + fPMTDiameter[ fData->getTelID() ] )
                        
                    )
                    
                    //Check if the pixel is not contained within a circle with x,y coordinates of the second best parametrization
                    //and radius of the second best parametrization minus a PMT diameter.
                    
                    && !(
                    
                        sqrt( ( fSecondBestParametrization[0] - fPixelXCoordinate ) * ( fSecondBestParametrization[0] - fPixelXCoordinate ) +
                              ( fSecondBestParametrization[1] - fPixelYCoordinate ) * ( fSecondBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fSecondBestParametrization[2] - fPMTDiameter[ fData->getTelID() ] )
                        
                    )
                    
                )// End of check to see if the pixel is hit by the second best parametrization
                
                
                || //Or
                
                
                ( //Check to see if the pixel is hit by the third best parametrization
                
                
                    //Check if the pixel is contained within a circle with x,y coordinates of the third best parametrization
                    //and radius of the third best parametrization plus a PMT diameter.
                    
                    (
                    
                        sqrt( ( fThirdBestParametrization[0] - fPixelXCoordinate ) * ( fThirdBestParametrization[0] - fPixelXCoordinate ) +
                              ( fThirdBestParametrization[1] - fPixelYCoordinate ) * ( fThirdBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fThirdBestParametrization[2] + fPMTDiameter[ fData->getTelID() ] )
                        
                    )
                    
                    //Check if the pixel is not contained within a circle with x,y coordinates of the third best parametrization
                    //and radius of the third best parametrization minus a PMT diameter.
                    
                    
                    && !(
                    
                        sqrt( ( fThirdBestParametrization[0] - fPixelXCoordinate ) * ( fThirdBestParametrization[0] - fPixelXCoordinate ) +
                              ( fThirdBestParametrization[1] - fPixelYCoordinate ) * ( fThirdBestParametrization[1] - fPixelYCoordinate ) ) <=
                        ( fThirdBestParametrization[2] - fPMTDiameter[ fData->getTelID() ] )
                        
                    )
                    
                    
                )// End of check to see if the pixel is hit by the third best parametrization
                
                
            )//End of if statement to check if the pixel is hit by any of the three best parametrizations
            
            
                //If the pixel is hit by any of the three best parametrizations, do this.
                //(Increment NPR variable, and calculate centroid of pixels included in the NPR count)
            {
            
                //Increment NPR variable if the pixel is hit by any of the three best parametrizations
                fNPR++;
                
                //Part of the calculation of the the centroid of the pixels included in the NPR count
                //Divided by NPR at the end of the loop
                fNPRCentroidX = fNPRCentroidX + fPixelXCoordinate;
                fNPRCentroidY = fNPRCentroidY + fPixelYCoordinate;
                
            }//End of incrementing NPR variable and calculating centroid
            
            
        }//End of check to see if the pixel is non-zero
        
        
    }//End of second loop over the pixels
    
    
    //Calculate the centroid of the pixels included in the NPR count. (Divide the quantities calculated above by NPR)
    fNPRCentroidX = fNPRCentroidX / ( double ) fNPR;
    fNPRCentroidY = fNPRCentroidY / ( double ) fNPR;
    
    //Calculate the CD parameter (the parameter formerly known as Cdist)
    //(The distance from the centre of the best parametrized circle to the centroid of the pixels included in the NPR count).
    fCD = sqrt( ( fBestParametrization[0] - fNPRCentroidX ) * ( fBestParametrization[0] - fNPRCentroidX ) +
                ( fBestParametrization[1] - fNPRCentroidY ) * ( fBestParametrization[1] - fNPRCentroidY ) );
                
    //Calculate the C/N parameter. (CD/NPR). Used to determine azimuthal completeness.
    fCN = fCD / ( double ) fNPR;
    
    //Calculate the distance from the center of the ring to the middle of the bestparametrization plus the radius of the best parametrization
    fContained = sqrt( ( fBestParametrization[0] * fBestParametrization[0] ) +
                       ( fBestParametrization[1] * fBestParametrization[1] ) ) + fBestParametrization[2];
                       
                       
    if //Check if the event passes the cuts (maybe load this from a cuts file later).
    
    (
    
        //////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////// Muon identification cuts /////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////
        
        //Cuts on the number of hit pixels.
        fNpix >= fNpixMin[ fData->getTelID() ] &&
        fNpix <= fNpixMax[ fData->getTelID() ] &&
        
        //Cut on the AP variable (a function of TD).
        fAP > ( fAlpha[ fData->getTelID() ] ) * fTD + fBeta[ fData->getTelID() ] &&
        
        //Cut on TD variable
        fTD < fGamma[ fData->getTelID() ] &&
        
        //Cut on the C/N variable (azimuthal completeness cut).
        fCN < fEta[ fData->getTelID() ] &&
        
        //Containedness cut.
        //The best parametrization is contained in a circle with the camera radius (defined in the cuts file).
        fContained < ( fCameraRadius[ fData->getTelID() ] ) * fPMTDiameter[ fData->getTelID() ]
        
        ///////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////
        
        
    )//End of if statement to check if the event passes the cuts
    
    {
    
        //If the event passes the cuts, set the fIsMuon flag to 1
        fIsMuon = 1;
        
        //If the event passes the cuts, increment the found muon counter.
        fNumberOfMuons++;
        
    }
    
    else
    
    {
    
        //If event does not pass the cuts, set the fIsMuon flag to 0
        fIsMuon = 0;
        
    }
    
    
    //Display the event number if the event passes the cuts
    //if( fIsMuon == 1 )
    //{
    
    //	cout << "Muon found on T" << ( fData->getTelID() + 1 ) << " at event " << fData->getEventNumber() << endl;
    //	cout << "Total: " << fNumberOfMuons << endl;
    
    //}
    
    
    //Write the Hough transform muon ID parameters to the tree
    
    fParGeo->houghMuonValid = fIsMuon;
    fParGeo->houghAP = fAP;
    fParGeo->houghTD = fTD;
    fParGeo->houghNpix = fNpix;
    fParGeo->houghCN = fCN;
    fParGeo->houghContained = fContained / fPMTDiameter[ fData->getTelID() ]; //Divide by PMT diameter.
    //In muon analysis macro, use houghContained < 11.0
    
}//End of analysis method



//Method for initializing the accumulator array
TH3D* VHoughTransform::initAccumulatorArray( int fRMinDpmt, int fRMaxDpmt, int fStepsPerPMTDiameter, unsigned int iTelescopeIndex )
{

    //The name of the accumulator array
    char fAccumulatorArrayName[100];
    
    //Maximum x value of the template circles
    double fXMax = 0;
    
    //Maximum y value of the template circles
    double fYMax = 0;
    
    //Step size in x
    double fStepSizeX = 0;
    
    //Step size in y
    double fStepSizeY = 0;
    
    
    //Determine the maximum values of the x and y coordinates of the circle templates
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
    {
    
        //Finds the largest value of the x coordinate
        if( fDetectorGeometry->getX_MM( iTelescopeIndex )[iChannelIndex] > fXMax )
        {
            fXMax = fDetectorGeometry->getX_MM( iTelescopeIndex )[iChannelIndex];
        }
        
        //Finds the largest value of the y coordinate
        if( fDetectorGeometry->getY_MM( iTelescopeIndex )[iChannelIndex] > fYMax )
        {
            fYMax = fDetectorGeometry->getY_MM( iTelescopeIndex )[iChannelIndex];
        }
        
    }//End of loop determining the maximum values of the x and y coordinates
    
    
    //Set the step size in x to the maximum x value
    fStepSizeX = fXMax;
    
    //Set the step size in y to the maximum x value
    fStepSizeY = fYMax;
    
    
    //Determine the smallest non-zero values of the x and y coordinates (the step sizes in x and y)
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
    
    {
        //Finds the smallest non-zero value of the x coordinate
        if( fabs( fDetectorGeometry->getX_MM( iTelescopeIndex )[iChannelIndex] ) < fStepSizeX && fabs( fDetectorGeometry->getX_MM( iTelescopeIndex )[iChannelIndex] ) > 0.0 )
        {
            fStepSizeX =  fabs( fDetectorGeometry->getX_MM( iTelescopeIndex )[iChannelIndex] );
        }
        
        //Finds the smallest non -zero value of the y coordinate
        if( fabs( fDetectorGeometry->getY_MM( iTelescopeIndex )[iChannelIndex] ) < fStepSizeY && fabs( fDetectorGeometry->getY_MM( iTelescopeIndex )[iChannelIndex] ) > 0.0 )
        {
            fStepSizeY =  fabs( fDetectorGeometry->getY_MM( iTelescopeIndex )[iChannelIndex] );
        }
    }//End of loop determining the smallest non-zero values of the x and y coordinates
    
    
    
    //Number of steps in X (approximately fXMax/fStepSizeX)
    int fNumberOfXSteps = ( int )( fXMax / fStepSizeX );
    
    //Number of steps in Y (approximately fYMax/fStepSizeY)
    int fNumberOfYSteps = ( int )( fYMax / fStepSizeY );
    
    
    //Maximum and minimum values of R, calculated from the previously declared parameters
    double fRMin = ( ( double ) fRMinDpmt ) * fPMTDiameter[ iTelescopeIndex ];
    double fRMax = ( ( double ) fRMaxDpmt ) * fPMTDiameter[ iTelescopeIndex ];
    
    //Step size in R
    double fStepSizeR = fPMTDiameter[ iTelescopeIndex ] / ( ( double ) fStepsPerPMTDiameter );
    
    //Number of bins in x
    int fNumberOfXBins = fNumberOfXSteps * 2 + 1;
    
    //Number of bins in y
    int fNumberOfYBins = fNumberOfYSteps * 2 + 1;
    
    //Number of bins in r
    int fNumberOfRBins = ( fRMaxDpmt - fRMinDpmt ) * ( fStepsPerPMTDiameter ) + 1;
    
    //Accumulator array ranges (x, y, r values are in the middle of the bin):
    
    double fAccumulatorXMin = - ( fXMax + ( fStepSizeX / 2.0 ) );
    double fAccumulatorXMax = ( fXMax + ( fStepSizeX / 2.0 ) );
    
    double fAccumulatorYMin = - ( fYMax + ( fStepSizeY / 2.0 ) );
    double fAccumulatorYMax = ( fYMax + ( fStepSizeY / 2.0 ) );
    
    double fAccumulatorRMin = ( fRMin - ( fStepSizeR / 2.0 ) );
    double fAccumulatorRMax = ( fRMax + ( fStepSizeR / 2.0 ) );
    
    
    //Set the name of the accumulator array
    sprintf( fAccumulatorArrayName, "Accumulator Array Tel%d", iTelescopeIndex );
    
    
    //Accumulator array instantiation
    TH3D* iAccumulatorArray = new TH3D( fAccumulatorArrayName, fAccumulatorArrayName,
                                        fNumberOfXBins, fAccumulatorXMin, fAccumulatorXMax,
                                        fNumberOfYBins, fAccumulatorYMin, fAccumulatorYMax,
                                        fNumberOfRBins, fAccumulatorRMin, fAccumulatorRMax );
                                        
                                        
    //Return the accumulator array
    return iAccumulatorArray;
    
    
}//End of method for initializing the accumulator array



//Method for initializing the Hough transform lookup table
TTree* VHoughTransform::initLookupTable( int fRMinDpmt, int fRMaxDpmt, int fStepsPerPMTDiameter, unsigned int iTelescopeIndex )
{

    //Lookup table name
    char fLookupTableName[100];
    
    //The number of circle templates used in the lookup table
    int fNumberOfCircleTemplates = 0;
    
    //Circle templates (pixel patterns)
    double fTemplateCircle[ fNumberOfChannels[ iTelescopeIndex ] ]; //Current circle template
    double fPreviousTemplateCircle[ fNumberOfChannels[ iTelescopeIndex ] ]; //Previous circle template (to check for duplicate templates)
    
    //Fill the template array with zeros
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
    {
        fTemplateCircle[iChannelIndex] = 0;
    }
    
    char fHTBranchName[100]; //Names of the branches in the HT lookup table tree
    
    char fHTLeafName[100]; //Leaf names
    
    double fTemplateCircleCoordinates[3]; //Template circle parametrization coordinates
    fTemplateCircleCoordinates[0] = 0; //x coordinate
    fTemplateCircleCoordinates[1] = 0; //y coordinate
    fTemplateCircleCoordinates[2] = 0; //r coordinate
    
    double fTestPixel[2]; //Coordinates of the center of a pixel being checked for 'hitness'
    fTestPixel[0] = 0; //X coordinate
    fTestPixel[1] = 0; //Y coordinate
    
    
    //Set the tbale of the lookup table tree
    sprintf( fLookupTableName, "Templates Tel%d", iTelescopeIndex );
    
    //Instantiate Hough transform lookup table tree.
    TTree* iHTLookupTableTree = new TTree( fLookupTableName , fLookupTableName );
    
    
    //Set up the branches in the tree. One brach for each pixel. Each brach will have circle parametrizations that hit that pixel.
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
    {
    
        //Sets the branch name
        sprintf( fHTBranchName , "Pixel %d" , iChannelIndex );
        
        //Sets the leaf name
        sprintf( fHTLeafName , "Pixel %d Circle Coordinates[3]/D" , iChannelIndex );
        
        //Setup the branch for that pixel
        iHTLookupTableTree->Branch( fHTBranchName , &fTemplateCircleCoordinates , fHTLeafName );
        
    }// End of setting up the pixel branches in the HT lookup table tree
    
    
    //Loop over the pixels for tempalte generation. (The center of the circle tempaltes is the center of the pixels)
    for( int iPixelCenterIndex = 0 ; iPixelCenterIndex < fNumberOfChannels[ iTelescopeIndex ] ; iPixelCenterIndex++ )
    {
    
        //Loop over radius values. Divide by fStepsPerPMTDiameter later when assigning the radius value.
        for( int iRadiusIndex = ( fRMinDpmt * fStepsPerPMTDiameter ) ; iRadiusIndex <= ( fRMaxDpmt * fStepsPerPMTDiameter ) ; iRadiusIndex++ )
        {
        
            //Make a duplicate circle template to check for duplicate circle templates.
            for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
            {
            
                //Copy the previous template to another array to check for duplicate templates
                fPreviousTemplateCircle[iChannelIndex] = fTemplateCircle[iChannelIndex];
                
            } //End of making a duplicate circle template
            
            
            
            //Set the coordinates of the template circle.
            fTemplateCircleCoordinates[0] = fDetectorGeometry->getX_MM( iTelescopeIndex )[iPixelCenterIndex]; 	//Get the x coordinate of the channel
            fTemplateCircleCoordinates[1] = fDetectorGeometry->getY_MM( iTelescopeIndex )[iPixelCenterIndex]; 	//Get the y coordinate of the channel
            
            //Set the radius coordiante of the template circle. Divide by fStepsPerPMTDiameter to calculate radius.
            fTemplateCircleCoordinates[2] = ( ( ( double ) iRadiusIndex ) * fPMTDiameter[ iTelescopeIndex ] ) / ( ( double ) fStepsPerPMTDiameter );
            
            //End of setting the coordinates of the template circle.
            
            
            //Loop over the pixels in the camera to test if they are hit by the circle template
            for( int iTestPixelIndex = 0 ; iTestPixelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iTestPixelIndex++ )
            {
            
                //Set the coordinates of the test pixel
                fTestPixel[0] = fDetectorGeometry->getX_MM( iTelescopeIndex )[iTestPixelIndex]; //Set x coordinate of test pixel
                
                fTestPixel[1] = fDetectorGeometry->getY_MM( iTelescopeIndex )[iTestPixelIndex]; //Set y coordinate of test pixel
                
                
                //Check if the test pixel is hit by the circle parametrization
                if( fabs( ( sqrt(
                                ( pow( ( fTestPixel[0] - fTemplateCircleCoordinates[0] ) , 2.0 ) ) +
                                ( pow( ( fTestPixel[1] - fTemplateCircleCoordinates[1] ) , 2.0 ) ) ) ) - fTemplateCircleCoordinates[2] ) 	<=
                        ( fPMTDiameter[ iTelescopeIndex ] / 2.0 ) )
                        
                {
                
                    //If the pixel is part of the ring, set the value of the pixel to one.
                    fTemplateCircle[iTestPixelIndex] = 1;
                    
                }
                
                //If the pixel is not part of the ring, set the pixel to zero.
                else
                {
                    fTemplateCircle[iTestPixelIndex] = 0;
                }
                
                //End of check to see if the pixel is part of the ring
                
                
            }// End of loop over pixels to see if they are hit by the circle template
            
            
            
            //Check to see if the circle template is a duplicate
            
            //Make the new template a duplicate by default
            bool iIsDuplicate = 1;
            
            //If there is a differnet pixel in the new template, then it is not a duplicate template.
            for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
            {
            
                if( fTemplateCircle[iChannelIndex] != fPreviousTemplateCircle[iChannelIndex] )
                {
                    iIsDuplicate = 0;
                }
                
            }
            
            //End of checking to see if the template is a duplicate
            
            
            
            //Fill the tree here
            
            
            //Add the circle coordinates to the nonzero pixels if the template is not a duplicate.
            if( !iIsDuplicate )
            
            {
            
                //Loop over the channels in the template
                for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
                
                {
                
                    //If the charge is non zero, fill the appropriate branch with the circle coordinates.
                    if( fTemplateCircle[iChannelIndex] != 0 )
                    
                    {
                    
                        //Fill the tree
                        
                        //Set the branch name
                        sprintf( fHTBranchName, "Pixel %d", iChannelIndex );
                        
                        //Fill the branch
                        iHTLookupTableTree->GetBranch( fHTBranchName )->Fill();
                        
                    }//End of checking if chargeval is non zero
                    
                    
                }//End of loop over the pixels
                
                
                fNumberOfCircleTemplates++; //Increment template index if the template is not a duplicate.
                
                
            }//End of check if duplicate and then fill tree
            
            
        }//End of loop over radius values
        
        
    }//End of loop over the centers of the pixels for template generation.
    
    
    //Write an array of zeros at the end of the branches of the lookup table
    
    //Set the template circle coordinates to zero.
    fTemplateCircleCoordinates[0] = 0.0;
    fTemplateCircleCoordinates[1] = 0.0;
    fTemplateCircleCoordinates[2] = 0.0;
    
    
    //For each pixel, write (0,0,0) at the end of the branch)
    for( int iChannelIndex = 0 ; iChannelIndex < fNumberOfChannels[ iTelescopeIndex ] ; iChannelIndex++ )
    {
    
        //Set the branch name
        sprintf( fHTBranchName, "Pixel %d", iChannelIndex );
        
        //Fill the branch with zeros
        iHTLookupTableTree->GetBranch( fHTBranchName )->Fill();
        
    } //End of writing zeros.
    
    
    return iHTLookupTableTree;
    
    
}//End of method for initializing the Hough transform lookup table


//Method for reading in the Hough transform muon ID parameter file
void VHoughTransform::readHTParameterFile( unsigned int fTelID )
{

    //cout << "Initializing parameters for telescope " << fTelID + 1 << "..." << endl;
    
    
    //Convert the telescope number (+1) into a string
    char iCharTelIDPlusOne[3];
    sprintf( iCharTelIDPlusOne , "%d" , ( fTelID + 1 ) );
    string iTelIDPlusOne = ( string )iCharTelIDPlusOne;
    
    
    //Variables to check if the parameters are set from the parameter file.
    //Set to 1 when value is read from the parameter file.
    int fRMinDpmtIsSet = 0;
    int fRMaxDpmtIsSet = 0;
    int fStepsPerPMTDiameterIsSet = 0;
    int fNpixMinIsSet = 0;
    int fNpixMaxIsSet = 0;
    int fAlphaIsSet = 0;
    int fBetaIsSet = 0;
    int fGammaIsSet = 0;
    int fEtaIsSet = 0;
    int fCameraRadiusIsSet = 0;
    
    
    //Read the parameter file here.
    
    
    //Get the parameter files directory
    string iParamsDir = fDetectorGeometry->getDirectory_EVNDISPParameterFiles();
    
    //Set the file name for the parameter file
    string ifilename = "Hough.runparameter.dat";
    
    //Set the full path and filename
    string ifile = iParamsDir + ifilename;
    
    //Open the file
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    
    
    //If the parameter file is not opened, do this
    if( !is )
    {
    
        //	cout << "VHoughTransform::readHTParameterFile error while opening parameter file: " << ifile << endl;
        
    }//End of if the parameter file is not opened
    
    
    //If the parameter file is opened, do this
    else
    {
    
        //cout << "Hough transform muon ID parameter file: " << ifile << " opened." << endl;
        
        //Strings for lines and individual strings
        string iLine;
        string iTemp;
        string iTemp2;
        string iTemp3;
        string iTemp4;
        string iTemp5;
        
        //While there are still lines to get in the file, get the next line
        while( getline( is, iLine ) )
        {
        
            // lines without '*' in the beginning are ignored
            if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*" )
            {
            
                istringstream is_stream( iLine );
                is_stream >> iTemp;
                
                // read variable identifier
                is_stream >> iTemp;
                iTemp = VUtilities::upperCase( iTemp );
                is_stream >> iTemp2;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp3;
                }
                else
                {
                    iTemp3 = "";
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp4;
                }
                else
                {
                    iTemp4 = "";
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp5;
                }
                else
                {
                    iTemp5 = "";
                }
                
                
                //fRminDpmt
                
                if( iTemp == "RMINT" + iTelIDPlusOne )
                {
                
                    fRMinDpmt.push_back( atoi( iTemp2.c_str() ) );
                    fRMinDpmtIsSet = 1;
                    
                }
                
                //fRmaxDpmt
                
                if( iTemp == "RMAXT" + iTelIDPlusOne )
                {
                
                    fRMaxDpmt.push_back( atoi( iTemp2.c_str() ) );
                    fRMaxDpmtIsSet = 1;
                    
                }
                
                //fStepsPerPMTDiameter
                
                if( iTemp == "STEPSPERPMTDIAMETERT" + iTelIDPlusOne )
                {
                
                    fStepsPerPMTDiameter.push_back( atoi( iTemp2.c_str() ) );
                    fStepsPerPMTDiameterIsSet = 1;
                    
                }
                
                //fNpixMin
                
                if( iTemp == "NPIXMINT" + iTelIDPlusOne )
                {
                
                    fNpixMin.push_back( atoi( iTemp2.c_str() ) );
                    fNpixMinIsSet = 1;
                    
                }
                
                //fNpixMax
                
                if( iTemp == "NPIXMAXT" + iTelIDPlusOne )
                {
                
                    fNpixMax.push_back( atoi( iTemp2.c_str() ) );
                    fNpixMaxIsSet = 1;
                    
                }
                
                //fAlpha
                
                if( iTemp == "ALPHAT" + iTelIDPlusOne )
                {
                
                    fAlpha.push_back( atof( iTemp2.c_str() ) );
                    fAlphaIsSet = 1;
                    
                }
                
                //fBeta
                
                if( iTemp == "BETAT" + iTelIDPlusOne )
                {
                
                    fBeta.push_back( atof( iTemp2.c_str() ) );
                    fBetaIsSet = 1;
                    
                }
                
                //fGamma
                
                if( iTemp == "GAMMAT" + iTelIDPlusOne )
                {
                
                    fGamma.push_back( atof( iTemp2.c_str() ) );
                    fGammaIsSet = 1;
                    
                }
                
                //fEta
                
                if( iTemp == "ETAT" + iTelIDPlusOne )
                {
                
                    fEta.push_back( atof( iTemp2.c_str() ) );
                    fEtaIsSet = 1;
                    
                }
                
                //fCameraRadius
                
                if( iTemp == "CAMERARADIUST" + iTelIDPlusOne )
                {
                
                    fCameraRadius.push_back( atof( iTemp2.c_str() ) );
                    fCameraRadiusIsSet = 1;
                    
                }
                
                
            }//if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*"
            
            
        }//while( getline( is, iLine ) )
        
        
    }//End of if the parameter file is opened
    
    
    //If the values are set by the parameter file, display the values.
    //If the values are not set by the parameter file, use default values and display message to that effect.
    
    
    //fRMinDpmt
    
    //if( fRMinDpmtIsSet == 1 )
    //{
    
    //	cout << "RMINT" << fTelID + 1 <<  " set to " << fRMinDpmt[fTelID] << endl;
    
    //}
    
    if( fRMinDpmtIsSet == 0 )
    {
    
        //Set the minimum radius in units of PMT diameters of the template circles to default value of 3
        fRMinDpmt.push_back( 3 );
        
        //cout << "Failed to read the RMINT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fRMinDpmt[fTelID] << "." << endl;
        
    }
    
    
    //fRMaxDpmt
    
    //if( fRMaxDpmtIsSet == 1 )
    //{
    
    //	cout << "RMAXT" << fTelID + 1 <<  " set to " << fRMaxDpmt[fTelID] << endl;
    
    //}
    
    if( fRMaxDpmtIsSet == 0 )
    {
    
        //Set the maximum radius in units of PMT diameters of the template circles to default value of 11
        fRMaxDpmt.push_back( 11 );
        
        //cout << "Failed to read the RMAXT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fRMaxDpmt[fTelID] << "." << endl;
        
    }
    
    
    //fStepsPerPMTDiameter
    
    //if( fStepsPerPMTDiameterIsSet == 1 )
    //{
    
    //	cout << "STEPSPERPMTDIAMETERT" << fTelID + 1 <<  " set to " << fStepsPerPMTDiameter[fTelID] << endl;
    
    //}
    
    if( fStepsPerPMTDiameterIsSet == 0 )
    {
    
        //Set the number of steps in radius per PMT diameter to defaut value of 3
        fStepsPerPMTDiameter.push_back( 3 );
        
        //cout << "Failed to read the STEPSPERPMTDIAMETERT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fStepsPerPMTDiameter[fTelID] << "." << endl;
        
    }
    
    
    //fNpixMin
    
    //if( fNpixMinIsSet == 1 )
    //{
    
    //	cout << "NPIXMINT" << fTelID + 1 <<  " set to " << fNpixMin[fTelID] << endl;
    
    //}
    
    if( fNpixMinIsSet == 0 )
    {
    
        //Set the min number of pixels cut to default value of 40
        fNpixMin.push_back( 40 );
        
        //cout << "Failed to read the NPIXMINT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fNpixMin[fTelID] << "." << endl;
        
    }
    
    
    //fNpixMax
    
    //if( fNpixMaxIsSet == 1 )
    //{
    
    //	cout << "NPIXMAXT" << fTelID + 1 <<  " set to " << fNpixMax[fTelID] << endl;
    
    //}
    
    if( fNpixMaxIsSet == 0 )
    {
    
        //Set the max number of pixels cut to default value of 79
        fNpixMax.push_back( 79 );
        
        //cout << "Failed to read the NPIXMAXT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fNpixMax[fTelID] << "." << endl;
        
    }
    
    
    //fAlpha
    
    //if( fAlphaIsSet == 1 )
    //{
    
    //	cout << "ALPHAT" << fTelID + 1 <<  " set to " << fAlpha[fTelID] << endl;
    
    //}
    
    if( fAlphaIsSet == 0 )
    {
    
        //Set the alpha Hough transform cut, described in memo and cuts file to default value of 0.011
        fAlpha.push_back( 0.011 );
        
        //cout << "Failed to read the ALPHAT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fAlpha[fTelID] << "." << endl;
        
    }
    
    
    //fBeta
    
    //if( fBetaIsSet == 1 )
    //{
    
    //	cout << "BETAT" << fTelID + 1 <<  " set to " << fBeta[fTelID] << endl;
    
    //}
    
    if( fBetaIsSet == 0 )
    {
    
        //Set the beta Hough transform cut, described in memo and cuts file to default value of 6.6
        fBeta.push_back( 6.6 );
        
        //cout << "Failed to read the BETAT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fBeta[fTelID] << "." << endl;
        
    }
    
    
    //fGamma
    
    //if( fGammaIsSet == 1 )
    //{
    
    //	cout << "GAMMAT" << fTelID + 1 <<  " set to " << fGamma[fTelID] << endl;
    
    //}
    
    if( fGammaIsSet == 0 )
    {
    
        //Set the gamma Hough transform cut, described in memo and cuts file to default value of 182.0
        fGamma.push_back( 182.0 );
        
        //cout << "Failed to read the GAMMAT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fGamma[fTelID] << "." << endl;
        
    }
    
    
    //fEta
    
    //if( fEtaIsSet == 1 )
    //{
    
    //	cout << "ETAT" << fTelID + 1 <<  " set to " << fEta[fTelID] << endl;
    
    //}
    
    if( fEtaIsSet == 0 )
    {
    
        //Set the azimuthal completeness cut to default value of 2.0
        fEta.push_back( 2.0 );
        
        //cout << "Failed to read the ETAT" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fEta[fTelID] << "." << endl;
        
    }
    
    
    //fCameraRadius
    
    //if( fCameraRadiusIsSet == 1 )
    //{
    
    //	cout << "CAMERARADIUST" << fTelID + 1 <<  " set to " << fCameraRadius[fTelID] << endl;
    
    //}
    
    if( fCameraRadiusIsSet == 0 )
    {
    
        //Set the containedness cut (camera radius in units of PMT diameters) to default value of 11.0
        fCameraRadius.push_back( 11.0 );
        
        //cout << "Failed to read the CAMERARADIUST" << fTelID + 1 <<  " value from the parameter file." << endl;
        //cout << "Using the default value of "  << fCameraRadius[fTelID] << "." << endl;
        
    }
    
    
    //cout << "Parameters for telescope " <<  fTelID + 1 << " initialized." << endl;
    
    
}//End of readHTParameterFile method


//Getter methods

//Gets the minimum number of hit pixels cut
int VHoughTransform::getNpixMin( unsigned int fTelID )
{

    return fNpixMin[fTelID];
    
}//End of getNpixMin


//Gets the maximum number of hit pixels cut
int VHoughTransform::getNpixMax( unsigned int fTelID )
{

    return fNpixMax[fTelID];
    
}//End of getNpixMax



