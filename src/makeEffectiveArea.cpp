/*! \file  makeEffectiveArea
 *  \brief get effective area and calculate instrument response functions from simulations
 *
 *   input is a file list of mscw_energy output file from gamma-ray simulations
 *
 */

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "VGlobalRunParameter.h"
#include "CData.h"
#include "Ctelconfig.h"
#include "VGammaHadronCuts.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VEffectiveAreaCalculator.h"
#include "VInstrumentResponseFunction.h"
#include "VInstrumentResponseFunctionRunParameter.h"
#include "VMonteCarloRunHeader.h"
#include "VTableLookupRunParameter.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

VEffectiveAreaCalculatorMCHistograms* copyMCHistograms( TChain* c );

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{

    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter iRunPara;
            cout << iRunPara.getEVNDISP_VERSION() << endl;
            exit( EXIT_SUCCESS );
        }
    }
    
    cout << endl;
    cout << "makeEffectiveArea " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    
    /////////////////////////////////////////////////////////////////
    // read command line parameters
    if( argc != 3 )
    {
        cout << endl;
        if( gSystem->Getenv( "EVNDISPSYS" ) )
        {
            int i_s = system( "cat $EVNDISPSYS/README/README.EFFECTIVEAREA" );
            if( i_s == -1 )
            {
                cout << "error: README/README.EFFECTIVEAREA not found" << endl;
            }
        }
        else
        {
            cout << "no help files found (environmental variable EVNDISPSYS not set)" << endl;
        }
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string fOutputfileName = argv[2];
    
    /////////////////////////////////////////////////////////////////
    // read run parameters from file
    VInstrumentResponseFunctionRunParameter* fRunPara = new VInstrumentResponseFunctionRunParameter();
    fRunPara->SetName( "makeEffectiveArea_runparameter" );
    if( !fRunPara->readRunParameterFromTextFile( argv[1] ) )
    {
        cout << "error reading runparameters from text file" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fRunPara->print();
    
    /////////////////////////////////////////////////////////////////
    // open output file and write results to dist
    TFile* fOutputfile = new TFile( fOutputfileName.c_str(), "RECREATE" );
    if( fOutputfile->IsZombie() )
    {
        cout << "Error in opening output file: " << fOutputfile->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    /////////////////////////////////////////////////////////////////
    // gamma/hadron cuts
    VGammaHadronCuts* fCuts = new VGammaHadronCuts();
    fCuts->initialize();
    fCuts->setNTel( fRunPara->telconfig_ntel, fRunPara->telconfig_arraycentre_X, fRunPara->telconfig_arraycentre_Y );
    fCuts->setInstrumentEpoch( fRunPara->fInstrumentEpoch );
    fCuts->setTelToAnalyze( fRunPara->fTelToAnalyse );
    fCuts->setReconstructionType( fRunPara->fReconstructionType );
    if( !fCuts->readCuts( fRunPara->fCutFileName, 2 ) )
    {
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE ) ;
    }
    fRunPara->fGammaHadronCutSelector = fCuts->getGammaHadronCutSelector();
    fRunPara->fDirectionCutSelector   = fCuts->getDirectionCutSelector();
    fCuts->initializeCuts( -1, fRunPara->fGammaHadronProbabilityFile );
    fCuts->printCutSummary();
    
    /////////////////////////////////////////////////////////////////
    // read MC header (might not be there, no problem; but depend on right input in runparameter file)
    VMonteCarloRunHeader* iMonteCarloHeader = fRunPara->readMCRunHeader();
    
    /////////////////////////////////////////////////////////////////
    // stopwatch to keep track of execution time
    TStopwatch fStopWatch;
    
    /////////////////////////////////////////////////////////////////////////////
    // set effective area class
    VEffectiveAreaCalculator fEffectiveAreaCalculator( fRunPara, fCuts );
    
    /////////////////////////////////////////////////////////////////////////////
    // set effective area Monte Carlo histogram class
    TFile* fMC_histoFile = 0;
    VEffectiveAreaCalculatorMCHistograms* fMC_histo = 0;
    
    /////////////////////////////////////////////////////////////////////////////
    // set angular, core, etc resolution calculation class
    vector< VInstrumentResponseFunction* > f_IRF;
    vector< string > f_IRF_Name;
    vector< string > f_IRF_Type;
    vector< float >  f_IRF_ContainmentProbability;
    vector< unsigned int > f_IRF_DuplicationID;
    if( fRunPara->fFillingMode != 3 )
    {
        // 68% angular resolution file
        f_IRF_Name.push_back( "angular_resolution" );
        f_IRF_Type.push_back( "angular_resolution" );
        f_IRF_ContainmentProbability.push_back( 0.68 );
        f_IRF_DuplicationID.push_back( 9999 );           // means no duplication
        // 80% angular resolution file
        f_IRF_Name.push_back( "angular_resolution_080p" );
        f_IRF_Type.push_back( "angular_resolution" );
        f_IRF_ContainmentProbability.push_back( 0.80 );
        f_IRF_DuplicationID.push_back( 0 );              // don't refill: use 68% histograms
        // 95% angular resolution file
        f_IRF_Name.push_back( "angular_resolution_095p" );
        f_IRF_Type.push_back( "angular_resolution" );
        f_IRF_ContainmentProbability.push_back( 0.95 );
        f_IRF_DuplicationID.push_back( 0 );              // don't refill: use 68% histograms
        if( fRunPara->fFillingMode != 2 && !fRunPara->fEffArea_short_writing )
        {
            // core resolution
            f_IRF_Name.push_back( "core_resolution" );
            f_IRF_Type.push_back( "core_resolution" );
            f_IRF_ContainmentProbability.push_back( 0.68 );
            f_IRF_DuplicationID.push_back( 9999 );           // means no duplication
            // energy resolution
            f_IRF_Name.push_back( "energy_resolution" );
            f_IRF_Type.push_back( "energy_resolution" );
            f_IRF_ContainmentProbability.push_back( 0.68 );
            f_IRF_DuplicationID.push_back( 9999 );           // means no duplication
        }
    }
    // initialize IRF classes
    for( unsigned int i = 0; i < f_IRF_Name.size(); i++ )
    {
        cout << "IRF for " << f_IRF_Name[i] << " (" << f_IRF_Type[i] << "): ";
        if( f_IRF_DuplicationID[i] != 9999 )
        {
            cout << "duplication ID " << f_IRF_DuplicationID[i];
        }
        cout << endl;
        f_IRF.push_back( new VInstrumentResponseFunction() );
        f_IRF.back()->setRunParameter( fRunPara );
        f_IRF.back()->setContainmentProbability( f_IRF_ContainmentProbability[i] );
        f_IRF.back()->initialize( f_IRF_Name[i], f_IRF_Type[i],
                                  fRunPara->telconfig_ntel, fRunPara->fCoreScatterRadius,
                                  fRunPara->fze, fRunPara->fnoise, fRunPara->fpedvar, fRunPara->fXoff, fRunPara->fYoff );
        f_IRF.back()->setDuplicationID( f_IRF_DuplicationID[i] );
    }
    
    
    /////////////////////////////////////////////////////////////////////////////
    // load data chain
    TChain* c = new TChain( "data" );
    if( !c->Add( fRunPara->fdatafile.c_str(), -1 ) )
    {
        cout << "Error while trying to add mscw data tree from file " << fRunPara->fdatafile  << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    CData d( c, true, true );
    fCuts->setDataTree( &d );
    d.setReconstructionType( fCuts->fReconstructionType );
    
    /////////////////////////////////////////////////////////////////////////////
    // fill resolution plots
    for( unsigned int i = 0; i < f_IRF_Name.size(); i++ )
    {
        if( f_IRF[i] )
        {
            f_IRF[i]->setDataTree( &d );
            f_IRF[i]->setCuts( fCuts );
            f_IRF[i]->setOutputFile( fOutputfile );
            if( f_IRF[i]->doNotDuplicateIRFs() )
            {
                f_IRF[i]->fill();
            }
            else if( f_IRF[i]->getDuplicationID() < f_IRF.size()
                     && f_IRF[f_IRF[i]->getDuplicationID()] )
            {
                f_IRF[i]->fillResolutionGraphs( f_IRF[f_IRF[i]->getDuplicationID()]->getIRFData() );
            }
            if( fCuts->getDirectionCutSelector() == 2 )
            {
                fCuts->setIRFGraph( f_IRF[i]->getAngularResolutionGraph( 0, 0 ) );
                f_IRF[i]->getAngularResolutionGraph( 0, 0 )->Print();
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // calculate effective areas
    if( !fRunPara->fFillMCHistograms )
    {
        // set azimuth bins and spectral index bins
        // (make sure that spectral index is positive)
        fEffectiveAreaCalculator.initializeHistograms( fRunPara->fAzMin, fRunPara->fAzMax, fRunPara->fSpectralIndex );
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // MC histograms
    if( fRunPara->fFillingMode != 1 && fRunPara->fFillingMode != 2 )
    {
        fStopWatch.Start();
        if( fRunPara->fMCdatafile_tree.size() == 0 && fRunPara->fMCdatafile_histo.size() == 0 )
        {
            fMC_histo = copyMCHistograms( c );
            if( fMC_histo )
            {
                fMC_histo->matchDataVectors( fRunPara->fAzMin, fRunPara->fAzMax, fRunPara->fSpectralIndex );
                fMC_histo->print();
            }
            else
            {
                cout << "Warning: failed reading MC histograms" << endl;
            }
        }
        // read MC histograms from a separate file
        else if( fRunPara->fMCdatafile_histo.size() > 0 )
        {
            fMC_histoFile = new TFile( fRunPara->fMCdatafile_histo.c_str() );
            if( fMC_histoFile->IsZombie() )
            {
                cout << "Error reading MC histograms from file " << fRunPara->fMCdatafile_histo << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            fMC_histo = ( VEffectiveAreaCalculatorMCHistograms* )fMC_histoFile->Get( "MChistos" );
            if( !fMC_histo )
            {
                cout << "Error reading MC histograms from file " << fRunPara->fMCdatafile_histo << " (no histograms)" << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            fMC_histo->matchDataVectors( fRunPara->fAzMin, fRunPara->fAzMax, fRunPara->fSpectralIndex );
            fMC_histo->print();
        }
        // recalulate MC spectra from MCpars tree. Very slow!
        else if( fRunPara->fMCdatafile_tree.size() > 0 && fRunPara->fMCdatafile_tree != "0" )
        {
            TChain* c2 = new TChain( "MCpars" );
            if( !c2->Add( fRunPara->fMCdatafile_tree.c_str(), -1 ) )
            {
                cout << "Error while trying to read MC data file: " << fRunPara->fMCdatafile_tree << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            fMC_histo = new VEffectiveAreaCalculatorMCHistograms();
            fMC_histo->setMonteCarloEnergyRange( fRunPara->fMCEnergy_min, fRunPara->fMCEnergy_max, TMath::Abs( fRunPara->fMCEnergy_index ) );
            fMC_histo->initializeHistograms( fRunPara->fAzMin, fRunPara->fAzMax, fRunPara->fSpectralIndex,
                                             fRunPara->fEnergyAxisBins_log10,
                                             fEffectiveAreaCalculator.getEnergyAxis_minimum_defaultValue(),
                                             fEffectiveAreaCalculator.getEnergyAxis_maximum_defaultValue() );
            fMC_histo->fill( fRunPara->fze, c2, fRunPara->fAzimuthBins );
            fMC_histo->print();
            fOutputfile->cd();
            cout << "writing MC histograms to file " << fOutputfile->GetName() << endl;
            fMC_histo->Write();
        }
        fStopWatch.Print();
    }
    
    // fill effective areas
    if( !fRunPara->fFillMCHistograms && fRunPara->fFillingMode != 1 && fRunPara->fFillingMode != 2 )
    {
        fOutputfile->cd();
        
        // copy angular resolution graphs to effective areas
        // assume same az bins in resolution and effective area calculation
        // use first spectral index bin
        for( unsigned int f = 0; f < f_IRF.size(); f++ )
        {
            if( f_IRF[f] && f_IRF[f]->getResolutionType() == "angular_resolution" )
            {
                if( TMath::Abs( f_IRF[f]->getContainmentProbability() - 0.68 ) < 1.e-4 )
                {
                    for( unsigned int i = 0; i < fRunPara->fAzMin.size(); i++ )
                    {
                        fEffectiveAreaCalculator.setAngularResolutionGraph( i,
                                f_IRF[f]->getAngularResolutionGraph( i, 0 ),
                                false );
                        fEffectiveAreaCalculator.setAngularResolution2D( i,
                                f_IRF[f]->getAngularResolution2D( i, 0 ) );
                    }
                }
                else if( TMath::Abs( f_IRF[f]->getContainmentProbability() - 0.80 ) < 1.e-4 )
                {
                    for( unsigned int i = 0; i < fRunPara->fAzMin.size(); i++ )
                    {
                        fEffectiveAreaCalculator.setAngularResolutionGraph( i,
                                f_IRF[f]->getAngularResolutionGraph( i, 0 ),
                                true );
                    }
                }
                
                for( unsigned int i = 0; i < fRunPara->fAzMin.size(); i++ )
                {
                    fEffectiveAreaCalculator.setAngularResolutionKingSigmaGraph( i,
                            f_IRF[f]->getAngularResolutionKingSigmaGraph( i, 0 ) );
                    fEffectiveAreaCalculator.setAngularResolutionKingGammaGraph( i,
                            f_IRF[f]->getAngularResolutionKingGammaGraph( i, 0 ) );
                }
            }
        }
        
        // fill effective areas
        fEffectiveAreaCalculator.fill( &d, fMC_histo, fRunPara->fEnergyReconstructionMethod );
        fStopWatch.Print();
    }
    
    
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    // write results to disk
    if( !fRunPara->fFillMCHistograms )
    {
        if( fEffectiveAreaCalculator.getEffectiveAreaTree() )
        {
            cout << "writing effective areas (";
            cout << fEffectiveAreaCalculator.getEffectiveAreaTree()->GetName();
            cout << ") to " << fOutputfile->GetName() << endl;
            fOutputfile->cd();
            fEffectiveAreaCalculator.getEffectiveAreaTree()->Write();
        }
        else
        {
            cout << "error: no effective area tree found" << endl;
        }
        // write one hEmc to disk used for binning
        if( fEffectiveAreaCalculator.getMCHistogram() )
        {
            fEffectiveAreaCalculator.getMCHistogram()->Write();
        }
        if( fEffectiveAreaCalculator.getMCHistogramUnWeighted() )
        {
            fEffectiveAreaCalculator.getMCHistogramUnWeighted()->Write();
        }
        
        if( fRunPara->fgetXoff_Yoff_afterCut && fEffectiveAreaCalculator.getAcceptance_AfterCuts() )
        {
            cout << "writing acceptance tree (";
            cout << fEffectiveAreaCalculator.getAcceptance_AfterCuts()->GetName();
            cout << ") to " << fOutputfile->GetName() << endl;
            fOutputfile->cd();
            fEffectiveAreaCalculator.getAcceptance_AfterCuts()->Write();
        }
        if( fRunPara->fWriteEventdatatrees 
        && fEffectiveAreaCalculator.getEventCutDataTree() )
        {
               cout << "writing event data trees: (";
               cout << fEffectiveAreaCalculator.getEventCutDataTree()->GetName();
               cout << ") to " << fOutputfile->GetName() << endl;
               fOutputfile->cd();
               if( fEffectiveAreaCalculator.getEventCutDataTree() )
               {
                   fEffectiveAreaCalculator.getEventCutDataTree()->Write();
               }
               if( c )
               {
                   c->Merge(fOutputfile, 0, "keep" );
               }
        }
        if( fCuts->getMVACutGraphs().size() > 0 )
        {
              cout << "writing MVA cut graphs to " << fOutputfile->GetName() << endl;
              fOutputfile->cd();
              if( fOutputfile->mkdir( "mvaGraphs" ) )
              {
                  fOutputfile->cd( "mvaGraphs" );   
                  for( unsigned int i = 0; i < fCuts->getMVACutGraphs().size(); i++ )
                  {
                      if( fCuts->getMVACutGraphs()[i] )
                      {
                          fCuts->getMVACutGraphs()[i]->Write();
                      }
                  }
              }
              fOutputfile->cd();
        }
    }
    // write resolution data to disk only for long output
    for( unsigned int i = 0; i < f_IRF_Name.size(); i++ )
    {
        if( f_IRF[i] && f_IRF[i]->getDataProduct() )
        {
            f_IRF[i]->getDataProduct()->Write();
        }
    }
    // writing cuts to disk
    if( fCuts )
    {
        fCuts->terminate( ( fRunPara->fEffArea_short_writing || fRunPara->fFillingMode == 3 ) );
    }
    // writing monte carlo header to disk
    if( iMonteCarloHeader )
    {
        iMonteCarloHeader->Write();
    }
    
    // write run parameters to disk
    if( fRunPara )
    {
        fRunPara->Write();
    }
    
    fOutputfile->Close();
    cout << "end..." << endl;
}

/*
 * return MC histograms from mscw_energy files
 *
 * - loop over all mscw_energy files
 * - add MC histograms
 */
VEffectiveAreaCalculatorMCHistograms* copyMCHistograms( TChain* c )
{
    VEffectiveAreaCalculatorMCHistograms* iMC_his = 0;
    if( c )
    {
        // loop over all files and add MC histograms
        TObjArray* fileElements = c->GetListOfFiles();
        TChainElement* chEl = 0;
        TIter next( fileElements );
        unsigned int z = 0;
        while( ( chEl = ( TChainElement* )next() ) )
        {
            TFile* ifInput = new TFile( chEl->GetTitle() );
            if( !ifInput->IsZombie() )
            {
                if( z == 0 )
                {
                    iMC_his = ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" );
                }
                else
                {
                    if( iMC_his )
                    {
                        iMC_his->add( ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" ) );
                    }
                    ifInput->Close();
                }
                z++;
            }
        }
    }
    return iMC_his;
}


