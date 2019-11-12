//! VImageAnalyzer    class for analyzing VERITAS data (single telescope analysis)
#ifndef VImageAnalyzer_H
#define VImageAnalyzer_H

#include "TFile.h"
#include "TMath.h"
#include "TObject.h"
#include "TTree.h"
#include "TROOT.h"

#include <VImageBaseAnalyzer.h>
#include <VImageParameterCalculation.h>
#include <VImageCleaning.h>

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VImageAnalyzer : public VImageBaseAnalyzer
{
    private:
        bool fDebug;
        VImageCleaning* fVImageCleaning;                            //!< image cleaning
        VImageParameterCalculation* fVImageParameterCalculation;    //!< image calculation
        
        bool fInit;
        
        // temporary vectors for dead pixel smoothing
        vector< unsigned int > savedDead;
        vector< unsigned int > savedDeadLow;
        valarray< double > savedGains;
        valarray< double > savedGainsLow;
        
        void fillOutputTree();                    //!< fill tree with image parameterisation results
        void imageCleaning( bool iDoublePass = false );  //!< image cleaning
        bool initEvent();                         //! reset image calculation for next event
        void muonRingAnalysis();                  //! muon ring analysis
        void houghMuonRingAnalysis();                  //! hough transform muon ring analysis
        void printTrace( int i_channel );         //!< print trace information for one channel (debugging)
        void setAnaDir( unsigned int iTel );      //!< set directories in root output file
        void setNTrigger();
        void smoothDeadTubes();                   //!< reduce the effect of dead tubes
        
    public:
        VImageAnalyzer();
        ~VImageAnalyzer();
        
        void doAnalysis();                        //!< do the actual analysis (called for each event)
        VImageCleaning*  getImageCleaner()
        {
            return fVImageCleaning;    //! return pointer to image cleaner
        }
        void initAnalysis();                      //!< set the data vectors, read the calibration data (called once at the beginning of the analysis)
        void initOutput();                        //!< open outputfile
        void initTrees();                         //!  intitalize output tree
        void shutdown();                          //!< close outputfile
        void terminate( bool iDebug_IO = false );            //!< write results to disk
};
#endif
