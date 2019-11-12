/*  \class VLightCurveWriter
    \brief write / print light curves in different formats

*/

#include "VLightCurveWriter.h"

VLightCurveWriter::VLightCurveWriter()
{
    fDebug = false;
    
}

VLightCurveWriter::VLightCurveWriter( vector< VFluxDataPoint > iDataVector )
{
    fDebug = false;
    
    setDataVector( iDataVector );
}

void VLightCurveWriter::setDataVector( vector< VFluxDataPoint > iDataVector )
{
    fFluxDataVector = iDataVector;
}



/*

   write results (fluxes and upper flux limits per run) into a root file

   ** very simple function, expand it if needed **

*/
bool VLightCurveWriter::writeDataVectorToRootFile( string i_root_file )
{
    TFile fOFILE( i_root_file.c_str(), "RECREATE" );
    if( fOFILE.IsZombie() )
    {
        cout << "VLightCurveWriter::writeResultsToRootFile error opening root file: " << i_root_file << endl;
        return false;
    }
    
    if( fDebug )
    {
        cout << "writing data vector to to " << fOFILE.GetName() << endl;
    }
    
    int iRun = 0;
    double iMJD = 0.;
    double iFlux = 0.;
    double iFluxE = 0.;
    double iSigni = 0.;
    double iZe = 0.;
    
    TTree t( "fluxes", "flux calculation results" );
    t.Branch( "Run", &iRun, "Run/I" );
    t.Branch( "MJD", &iMJD, "MJD/D" );
    t.Branch( "Flux", &iFlux, "Flux/D" );
    t.Branch( "FluxE", &iFluxE, "FluxE/D" );
    t.Branch( "Signi", &iSigni, "Signi/D" );
    t.Branch( "Ze", &iZe, "Ze/D" );
    
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        iRun =   fFluxDataVector[i].fRunNumber;
        iMJD   = fFluxDataVector[i].fMJD;
        iFlux  = fFluxDataVector[i].fFlux;
        iFluxE = fFluxDataVector[i].fFluxE;
        iSigni = fFluxDataVector[i].fSignificance;
        iZe    = fFluxDataVector[i].fZe;
        t.Fill();
    }
    t.Write();
    
    fOFILE.Close();
    
    return true;
}

/*
   write fluxes to a simple ascii file
*/
bool VLightCurveWriter::writeASCIIt_SimpleLongTable( string ASCIIFile, double iMultiplier )
{
    ofstream is;
    is.open( ASCIIFile.c_str() );
    if( !is )
    {
        cout << "error opening " << ASCIIFile << endl;
        return false;
    }
    cout << "writing flux data vector to ascii file: " << ASCIIFile << endl;
    cout << endl;
    
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        is << setprecision( 3 ) << fixed << setw( 9 ) << fFluxDataVector[i].fMJD << "\t";
        is << setprecision( 3 ) << scientific <<  fFluxDataVector[i].fFlux* iMultiplier << "\t";
        is << setprecision( 3 ) << scientific << fFluxDataVector[i].fFluxCI_lo_1sigma* iMultiplier << "\t";
        is << setprecision( 3 ) << scientific << fFluxDataVector[i].fFluxCI_up_1sigma* iMultiplier << endl;
    }
    is.close();
    
    return true;
}


/*
   write results as a simple latex long table
*/
bool VLightCurveWriter::writeTexFormat_SimpleLongTable( string iTexFile, double iMultiplier )
{
    ofstream is;
    is.open( iTexFile.c_str() );
    if( !is )
    {
        cout << "error opening " << iTexFile << endl;
        return false;
    }
    cout << "writing flux data vector to tex file: " << iTexFile << endl;
    cout << endl;
    
    is << "\\documentclass[a4paper]{article}" << endl;
    is << "\\usepackage{longtable}" << endl;
    is << "\\usepackage{lscape}" << endl;
    
    is << "\\begin{document}" << endl;
    
    is << endl;
    is << "\\begin{longtable}{c|c}" << endl;
    is << "MJD \\\\" << endl;
    is << "Flux \\\\" << endl;
    is << "\\hline" << endl;
    is << "\\hline" << endl;
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        is << fixed << setw( 11 ) << fFluxDataVector[i].fMJD << " & ";
        is << setprecision( 3 ) << scientific <<  fFluxDataVector[i].fFlux* iMultiplier;
        is << "$\\pm$" << setprecision( 3 ) << fFluxDataVector[i].fFluxCI_1sigma* iMultiplier << "\\\\" << endl;
    }
    is << "\\end{longtable}" << endl;
    is << "\\end{document}" << endl;
    is.close();
    
    return true;
}


/*

   print a row for a typical latex table

*/
void VLightCurveWriter::writeTexFormat_TexTableRow( double iSigmaMinFluxLimits, double iFluxMultiplicator, bool iPrintPhaseValues )
{
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        cout << ( int )fFluxDataVector[i].fMJD_Start << " - " << ( int )fFluxDataVector[i].fMJD_Stop << " & ";
        if( fFluxDataVector[i].hasOrbitalPhases() && iPrintPhaseValues )
        {
            cout << setprecision( 2 ) << fFluxDataVector[i].fOrbitalPhase << " & ";
        }
        cout << "VERITAS & ";
        // observing time in minutes
        cout << ( int )( fFluxDataVector[i].fExposure_deadTimeCorrected / 60. ) << " & ";
        // mean elevation
        cout << setprecision( 1 ) << fixed << 90. - fFluxDataVector[i].fZe << " & ";
        // on and off events
        cout << ( int )fFluxDataVector[i].fNon  << " & ";
        cout << ( int )fFluxDataVector[i].fNoff << " & ";
        // alpha
        cout << setprecision( 2 ) << fixed << fFluxDataVector[i].fAlpha << " & ";
        // significance
        cout << setprecision( 1 ) << fFluxDataVector[i].fSignificance << " & ";
        // flux (with error) or upper flux limit)
        if( iSigmaMinFluxLimits != 1 )
        {
            cout << fixed;
        }
        else
        {
            cout << scientific;
        }
        if( fFluxDataVector[i].fSignificance > iSigmaMinFluxLimits && iSigmaMinFluxLimits > -1.e3 )
        {
            cout << setprecision( 1 ) << fFluxDataVector[i].fFlux* iFluxMultiplicator << " $\\pm$ (";
            cout << fFluxDataVector[i].fFluxCI_up_1sigma* iFluxMultiplicator << ",";
            cout << fFluxDataVector[i].fFluxCI_lo_1sigma* iFluxMultiplicator << ")";
        }
        else if( iSigmaMinFluxLimits < -1.e3 )
        {
            cout << setprecision( 1 ) << fFluxDataVector[i].fFlux* iFluxMultiplicator << " $\\pm$ (";
            cout << fFluxDataVector[i].fFluxCI_up_1sigma* iFluxMultiplicator << ",";
            cout << fFluxDataVector[i].fFluxCI_lo_1sigma* iFluxMultiplicator << ")";
            cout << " $(< " << fFluxDataVector[i].fUL* iFluxMultiplicator << ")";
        }
        else
        {
            cout << " $<$ " << fFluxDataVector[i].fUL* iFluxMultiplicator;
        }
        cout << " \\\\";
        cout << endl;
    }
}

/*

    write a table for the VERITAS wiki

*/
void VLightCurveWriter::writeWikiFormat()
{
    double iMinEnergy_TeV = 0.;
    bool iHasOrbitalPhases = false;
    if( fFluxDataVector.size() > 0 )
    {
        iMinEnergy_TeV = fFluxDataVector[0].fMinEnergy_TeV;
        iHasOrbitalPhases = fFluxDataVector[0].hasOrbitalPhases();
    }
    
    cout << "{| border=\"1\" cellspacing=\"0\" cellpadding=\"5\" align=\"center\"" << endl;
    cout << "!MJD" << endl;
    if( iHasOrbitalPhases )
    {
        cout << "!Phase" << endl;
    }
    cout << "!Observation Time [min]" << endl;
    cout << "!Significance <math>\\sigma</math>" << endl;
    cout << "!Non" << endl;
    cout << "!Noff" << endl;
    cout << "!Alpha" << endl;
    cout << "!Flux (>" << setprecision( 2 ) << iMinEnergy_TeV << " TeV) [cm^-2 s^-1]" << endl;
    
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        cout << "|- align=\"center\"" << endl;
        cout << "| " << fixed << setprecision( 1 ) << fFluxDataVector[i].fMJD << endl;
        if( iHasOrbitalPhases )
        {
            cout << "| " << setprecision( 2 ) << fFluxDataVector[i].fOrbitalPhase << endl;
        }
        cout << "| " << setprecision( 1 ) << fFluxDataVector[i].fExposure_deadTimeCorrected / 60. << endl;
        cout << "| " << setprecision( 1 ) << fFluxDataVector[i].fSignificance << endl;
        cout << "| " << ( int )fFluxDataVector[i].fNon << endl;
        cout << "| " << ( int )fFluxDataVector[i].fNoff << endl;
        cout << "| " << setprecision( 2 ) << fFluxDataVector[i].fAlpha << endl;
        if( fFluxDataVector[i].isSignificantDataPoint() )
        {
            cout << "| " << setprecision( 1 ) << scientific << fFluxDataVector[i].fFlux;
            cout << " (+" << fFluxDataVector[i].fRate_up_1sigma;
            cout << " ,-" << fFluxDataVector[i].fRate_lo_1sigma << ")" << endl;
        }
        else
        {
            cout << "| " << setprecision( 1 ) << scientific << fFluxDataVector[i].fFluxUL << endl;
        }
    }
    
    cout << "|}" << fixed << endl;
}

/*
   write ligth curve for discrete correlation function (DCF and ZDCF) analysis

*/
void VLightCurveWriter::writeDCFFormat()
{
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        cout << "DCF\t";
        cout << fixed << setprecision( 4 ) << fFluxDataVector[i].fMJD << "\t";
        cout << scientific << setprecision( 4 ) << fFluxDataVector[i].fFlux << "\t";
        cout << fFluxDataVector[i].fFluxCI_1sigma;
        cout << fixed << endl;
    }
}

