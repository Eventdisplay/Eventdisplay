/*! \class VImageAnalyzerHistograms
    \brief  histogramming of run parameter, etc.

*/

#include "VImageAnalyzerHistograms.h"

/*!
    \param iTel Telescope number
*/
VImageAnalyzerHistograms::VImageAnalyzerHistograms( unsigned int iTel )
{
    fTelescopeID = iTel;
    
    hisList = new TList();
}


VImageAnalyzerHistograms::~VImageAnalyzerHistograms()
{
}


void VImageAnalyzerHistograms::init()
{

    char hisname[800];
    char histitle[800];
    
    sprintf( hisname, "fFADCStop" );
    sprintf( histitle, "FADC stop diagnostic tree (telescope %d)", fTelescopeID + 1 );
    fdiagno = new TTree( hisname, histitle );
    fdiagno->SetAutoSave( 100000000 );
    
    fdiagno->Branch( "runNumber", &runNumber, "runNumber/I" );
    fdiagno->Branch( "eventNumber",  &eventNumber,  "eventNumber/I" );
    fdiagno->Branch( "MJD",  &MJD,  "MJD/I" );
    fdiagno->Branch( "Time",  &time,  "time/F" );
    fdiagno->Branch( "FADCstopTZero", fFADCstopTZero, "fFADCstopTZero[4]/F" );
    fdiagno->Branch( "FADCstopSum", fFADCstopSum, "fFADCstopSum[4]/F" );
    hisList->Add( fdiagno );
    
}


/*
    fill diagnostic trees

    expect for FADC stop channels not more than 4 telescopes

*/
void VImageAnalyzerHistograms::fillL2DiagnosticTree( int rN, int eN, int iMJD, double it, vector< double >& iTZero, vector< double >& iFSum )
{
    runNumber = rN;
    eventNumber = eN;
    MJD = iMJD;
    time = it;
    
    if( iTZero.size() == 4 )
    {
        for( unsigned int i = 0; i < 4; i++ )
        {
            fFADCstopTZero[i] = iTZero[i];
        }
    }
    else for( unsigned int i = 0; i < 4; i++ )
        {
            fFADCstopTZero[i] = 0.;
        }
    if( iFSum.size() == 4 )
    {
        for( unsigned int i = 0; i < 4; i++ )
        {
            fFADCstopSum[i] = iFSum[i];
        }
    }
    else for( unsigned int i = 0; i < 4; i++ )
        {
            fFADCstopSum[i] = 0.;
        }
        
    if( fdiagno )
    {
        fdiagno->Fill();
    }
}


/*!
   writing all histograms to disk

   outputfile is tree outputfile from VAnalyzer
*/
void VImageAnalyzerHistograms::terminate( TFile* outputfile )
{
    if( outputfile == 0 )
    {
        return;
    }
    TDirectory* iDir = gDirectory;
    
    hisList->Write();
    
    outputfile->cd();
    iDir->cd();
}

