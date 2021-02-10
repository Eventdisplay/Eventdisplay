/* \file writeCTAWPPhysSensitivityOptimisedAndSmoothedFiles.cpp
 *       find and write best sensitivity files
 *
 *       - search e.g. for best sensitivity from different multiplicity
 *         file
 *       - smooth sensitivity and IRF files
 *
 *   Expect that the binning in all input files is identical
 *
 */


#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

/*
 * data class with all IRFs
 * (service class for VWPPhysFile)
 *
 */
class VWPPhysHistograms
{
    public:
        map< string, TH1F* > fHisto1F;
        map< string, TH2F* > fHisto2F;
        map< string, TH2D* > fHisto2D;
        
        VWPPhysHistograms();
        ~VWPPhysHistograms() {}
        bool clone( VWPPhysHistograms* iHistos );
        bool fill( int iBin, double iBinLowE, double iBinUpE,
                   VWPPhysHistograms* iHistos,
                   unsigned int iMult_LST, unsigned int iMult_MST, unsigned int iMult_SST );
        bool write();
};

/*
 * list of histogram names to be filled
 *
 */
VWPPhysHistograms::VWPPhysHistograms()
{
    //////////////////////////////////////
    // one dimensional histograms
    fHisto1F[ "DiffSens" ] = 0;
    fHisto1F[ "DiffSensCU" ] = 0;
    fHisto1F[ "BGRate" ] = 0;
    fHisto1F[ "ProtRate" ] = 0;
    fHisto1F[ "ElecRate" ] = 0;
    fHisto1F[ "BGRatePerSqDeg" ] = 0;
    fHisto1F[ "ProtRateSqDeg" ] = 0;
    fHisto1F[ "ElecRateSqDeg" ] = 0;
    fHisto1F[ "EffectiveArea" ] = 0;
    fHisto1F[ "EffectiveAreaEtrue" ] = 0;
    fHisto1F[ "EffectiveAreaEtrueNoTheta2cut" ] = 0;
    fHisto1F[ "EffectiveAreaNoTheta2cut" ] = 0;
    fHisto1F[ "EffectiveArea80" ] = 0;
    fHisto1F[ "AngRes" ] = 0;
    fHisto1F[ "AngRes80" ] = 0;
    fHisto1F[ "AngRes95" ] = 0;
    fHisto1F[ "ERes" ] = 0;
    fHisto1F[ "Ebias" ] = 0;
    
    //////////////////////////////////////
    // two dimensional histograms
    fHisto2F[ "DiffSens_offaxis" ] = 0;
    fHisto2F[ "DiffSensCU_offaxis" ] = 0;
    fHisto2F[ "BGRate_offaxis" ] = 0;
    fHisto2F[ "ProtRate_offaxis" ] = 0;
    fHisto2F[ "ElecRate_offaxis" ] = 0;
    fHisto2F[ "BGRatePerSqDeg_offaxis" ] = 0;
    fHisto2F[ "ProtRateSqDeg_offaxis" ] = 0;
    fHisto2F[ "ElecRateSqDeg_offaxis" ] = 0;
    fHisto2F[ "EffectiveArea_offaxis" ] = 0;
    fHisto2F[ "EffectiveAreaEtrue_offaxis" ] = 0;
    fHisto2F[ "EffectiveAreaEtrueNoTheta2cut_offaxis" ] = 0;
    fHisto2F[ "EffectiveAreaNoTheta2cut_offaxis" ] = 0;
    fHisto2F[ "EffectiveArea80_offaxis" ] = 0;
    fHisto2F[ "AngRes_offaxis" ] = 0;
    fHisto2F[ "AngRes80_offaxis" ] = 0;
    fHisto2F[ "AngRes95_offaxis" ] = 0;
    fHisto2F[ "ERes_offaxis" ] = 0;
    fHisto2F[ "Ebias_offaxis" ] = 0;
    
    // for some stupid reasons, there are two 2D (not 2F) histograms
    fHisto2D[ "MigMatrix" ] = 0;
    fHisto2D[ "EestOverEtrue" ] = 0;
    
    //////////////////////////////////////
    // additional histograms
    fHisto1F[ "MultiplicityLST" ] = 0;
    fHisto1F[ "MultiplicityMST" ] = 0;
    fHisto1F[ "MultiplicitySST" ] = 0;
}

/*
 * write all histograms into an open file
 *
 */
bool VWPPhysHistograms::write()
{
    // write 1F histograms
    for( map< string, TH1F* >::iterator it = fHisto1F.begin();
            it != fHisto1F.end();
            ++it )
    {
        if( it->second )
        {
            it->second->Write();
        }
    }
    
    // write 2F histograms
    for( map< string, TH2F* >::iterator it = fHisto2F.begin();
            it != fHisto2F.end();
            ++it )
    {
        if( it->second )
        {
            it->second->Write();
        }
    }
    for( map< string, TH2D* >::iterator it = fHisto2D.begin();
            it != fHisto2D.end();
            ++it )
    {
        if( it->second )
        {
            it->second->Write();
        }
    }
    
    return true;
}

/*
 * clone a set of histograms into a new (empty) set
 *
 */
bool VWPPhysHistograms::clone( VWPPhysHistograms* iHistos )
{
    // clone 1F histograms
    for( map< string, TH1F* >::iterator it = fHisto1F.begin();
            it != fHisto1F.end();
            ++it )
    {
        if( iHistos->fHisto1F[it->first] )
        {
            it->second = ( TH1F* )iHistos->fHisto1F[it->first]->Clone();
            it->second->Reset();
        }
    }
    // clone 2F histograms
    for( map< string, TH2F* >::iterator it = fHisto2F.begin();
            it != fHisto2F.end();
            ++it )
    {
        if( iHistos->fHisto2F[it->first] )
        {
            it->second = ( TH2F* )iHistos->fHisto2F[it->first]->Clone();
            it->second->Reset();
        }
    }
    // clone 2D histograms
    for( map< string, TH2D* >::iterator it = fHisto2D.begin();
            it != fHisto2D.end();
            ++it )
    {
        if( iHistos->fHisto2D[it->first] )
        {
            it->second = ( TH2D* )iHistos->fHisto2D[it->first]->Clone();
            it->second->Reset();
        }
    }
    
    // minimal multiplicity histograms
    // (filled during combination of histograms)
    if( fHisto1F.find( "DiffSens" ) !=  fHisto1F.end() )
    {
        fHisto1F[ "MultiplicityLST" ] = ( TH1F* )fHisto1F[ "DiffSens" ]->Clone();
        fHisto1F[ "MultiplicityLST" ]->SetName( "MultiplicityLST" );
        fHisto1F[ "MultiplicityLST" ]->GetYaxis()->SetTitle( "LST minimal multiplicity" );
        fHisto1F[ "MultiplicityLST" ]->SetTitle( "LST minimal multiplicity" );
        fHisto1F[ "MultiplicityLST" ]->Reset();
        
        fHisto1F[ "MultiplicityMST" ] = ( TH1F* )fHisto1F[ "DiffSens" ]->Clone();
        fHisto1F[ "MultiplicityMST" ]->SetName( "MultiplicityMST" );
        fHisto1F[ "MultiplicityMST" ]->GetYaxis()->SetTitle( "MST minimal multiplicity" );
        fHisto1F[ "MultiplicityMST" ]->SetTitle( "MST minimal multiplicity" );
        fHisto1F[ "MultiplicityMST" ]->Reset();
        
        fHisto1F[ "MultiplicitySST" ] = ( TH1F* )fHisto1F[ "DiffSens" ]->Clone();
        fHisto1F[ "MultiplicitySST" ]->SetName( "MultiplicitySST" );
        fHisto1F[ "MultiplicitySST" ]->GetYaxis()->SetTitle( "SST minimal multiplicity" );
        fHisto1F[ "MultiplicitySST" ]->SetTitle( "SST minimal multiplicity" );
        fHisto1F[ "MultiplicitySST" ]->Reset();
    }
    
    return true;
}

/*
 * copy a value of a certain bin from a set of histograms
 * to the current histograms
 */
bool VWPPhysHistograms::fill( int b,
                              double iBinLowE, double iBinUpE,
                              VWPPhysHistograms* iHistos,
                              unsigned int iMult_LST,
                              unsigned int iMult_MST,
                              unsigned int iMult_SST )
{
    ////////////////////////////////////////////////
    // fill 1F histograms
    for( map< string, TH1F* >::iterator it = fHisto1F.begin();
            it != fHisto1F.end();
            ++it )
    {
        if( iHistos->fHisto1F.find( it->first ) != iHistos->fHisto1F.end()
                && iHistos->fHisto1F[it->first]
                && it->second )
        {
            it->second->SetBinContent( b, iHistos->fHisto1F[it->first]->GetBinContent( b ) );
            it->second->SetBinError( b, iHistos->fHisto1F[it->first]->GetBinError( b ) );
        }
    }
    
    ////////////////////////////////////////////////
    // do not optimize 2F histograms -> simply copy
    // decision made from 1F
    // fill 2F histograms
    for( map< string, TH2F* >::iterator it = fHisto2F.begin();
            it != fHisto2F.end();
            ++it )
    {
        if( iHistos->fHisto2F.find( it->first ) != iHistos->fHisto2F.end()
                && iHistos->fHisto2F[it->first]
                && it->second )
        {
            for( int d = 0; d <= iHistos->fHisto2F[it->first]->GetNbinsY(); d++ )
            {
                it->second->SetBinContent( b, d, iHistos->fHisto2F[it->first]->GetBinContent( b, d ) );
                it->second->SetBinError( b, d, iHistos->fHisto2F[it->first]->GetBinError( b, d ) );
            }
        }
    }
    ////////////////////////////////////////////////
    // 2D histograms have a different binning in energy than all other histograms
    //
    for( map< string, TH2D* >::iterator it = fHisto2D.begin();
            it != fHisto2D.end();
            ++it )
    {
        if( iHistos->fHisto2D.find( it->first ) != iHistos->fHisto2D.end()
                && iHistos->fHisto2D[it->first]
                && it->second )
        {
            for( int x = iHistos->fHisto2D[it->first]->GetXaxis()->FindBin( iBinLowE );
                    x < iHistos->fHisto2D[it->first]->GetXaxis()->FindBin( iBinUpE );
                    x++ )
            {
                for( int d = 0; d <= iHistos->fHisto2D[it->first]->GetNbinsY(); d++ )
                {
                    it->second->SetBinContent( x, d, iHistos->fHisto2D[it->first]->GetBinContent( x, d ) );
                    it->second->SetBinError( x, d, iHistos->fHisto2D[it->first]->GetBinError( x, d ) );
                }
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // multiplicity histograms
    // note that this is minimal applied multiplicity
    
    // LST Multiplicity
    if( fHisto1F.find( "MultiplicityLST" ) != fHisto1F.end()
            && fHisto1F[ "MultiplicityLST" ] )
    {
        fHisto1F[ "MultiplicityLST" ]->SetBinContent( b, iMult_LST );
    }
    // MST Multiplicity
    if( fHisto1F.find( "MultiplicityMST" ) != fHisto1F.end()
            && fHisto1F[ "MultiplicityMST" ] )
    {
        fHisto1F[ "MultiplicityMST" ]->SetBinContent( b, iMult_MST );
    }
    // SST Multiplicity
    if( fHisto1F.find( "MultiplicitySST" ) != fHisto1F.end()
            && fHisto1F[ "MultiplicitySST" ] )
    {
        fHisto1F[ "MultiplicitySST" ]->SetBinContent( b, iMult_SST );
    }
    
    
    return true;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/*
 * class dealing with WP Phys files and
 * hosting the histograms
 *
 * (service class of VCombineWPPhysFiles)
 *
 */
class VWPPhysFile
{
    public:
    
        unsigned int fMultiplicity_LST;
        unsigned int fMultiplicity_MST;
        unsigned int fMultiplicity_SST;
        string fFileNameBegin;
        string fFileNameMultiplicity;
        
        TFile* fDatFile;
        VWPPhysHistograms* fHistos;
        
        VWPPhysFile( string iNameBegin,
                     string iNameMultiplicity,
                     unsigned int iMultiLST = 2,
                     unsigned int iMultiMST = 2,
                     unsigned int iMultiSST = 2 );
        ~VWPPhysFile();
        
        bool open( string iDatadirectory, string iSite, string iSubArray, int iObsTime, int iAnaID, string iPointingDirection );
};

VWPPhysFile::VWPPhysFile( string iNameBegin, string iNameMultiplicity,
                          unsigned int iMultiLST, unsigned int iMultiMST, unsigned int iMultiSST )
{
    fHistos = new VWPPhysHistograms();
    fMultiplicity_LST = iMultiLST;
    fMultiplicity_MST = iMultiMST;
    fMultiplicity_SST = iMultiSST;
    fFileNameBegin = iNameBegin;
    fFileNameMultiplicity = iNameMultiplicity;
    
    fDatFile = 0;
}

VWPPhysFile::~VWPPhysFile()
{
    if( fDatFile )
    {
        fDatFile->Close();
    }
}

/*
 * open a root file and read in the IRFs
 */
bool VWPPhysFile::open( string iDatadirectory, string iSite, string iSubArray, int iObsTime, int iAnaID, string iPointingDirection )
{
    ostringstream iWPPhysFileName;
    
    iWPPhysFileName << iDatadirectory << "/";
    iWPPhysFileName << fFileNameBegin;
    iWPPhysFileName << ".ID" << iAnaID;
    iWPPhysFileName << iPointingDirection;
    iWPPhysFileName << fFileNameMultiplicity;
    iWPPhysFileName << "." << iSite << ".";
    iWPPhysFileName << iSubArray << ".";
    iWPPhysFileName << iObsTime;
    iWPPhysFileName << "s.root";
    
    cout << iWPPhysFileName.str() << endl;
    
    TFile* fDatFile = new TFile( iWPPhysFileName.str().c_str() );
    if( fDatFile->IsZombie() )
    {
        cout << "Error opening " << iWPPhysFileName.str() << endl;
        return false;
    }
    // read 1F histograms
    for( map< string, TH1F* >::iterator it = fHistos->fHisto1F.begin();
            it != fHistos->fHisto1F.end();
            ++it )
    {
        it->second = ( TH1F* )fDatFile->Get( it->first.c_str() );
    }
    // read 2F histograms
    for( map< string, TH2F* >::iterator it = fHistos->fHisto2F.begin();
            it != fHistos->fHisto2F.end();
            ++it )
    {
        it->second = ( TH2F* )fDatFile->Get( it->first.c_str() );
    }
    // read 2D histograms
    for( map< string, TH2D* >::iterator it = fHistos->fHisto2D.begin();
            it != fHistos->fHisto2D.end();
            ++it )
    {
        it->second = ( TH2D* )fDatFile->Get( it->first.c_str() );
    }
    cout << "\tread " << fHistos->fHisto1F.size() << " 1F, ";
    cout << fHistos->fHisto2F.size() << " 2F, ";
    cout << fHistos->fHisto2D.size() << " 2D histograms" << endl;
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * class in which the WP Phys files
 * are combined
 */
class VCombineWPPhysFiles
{
    private:
    
        TFile* fOutFile;
        string fOutFileName;
        VWPPhysHistograms fOutHistos;
        vector< VWPPhysFile* > fPhysFile;
        
    public:
        VCombineWPPhysFiles( string iOutFileName = "combinedFile" );
        ~VCombineWPPhysFiles() {}
        
        bool combineWPPhysFiles();
        bool readWPPhysFiles( string iDatadirectory, string iSite, string iSubArray, int iObsTime, int iAnaID, string iPointingDirection );
        bool writeWPPhysFiles();
};

/*
 * constructor for VCombineWPPhysFiles
 *
 * here the file name and the different multiplicities are hardwired
 *
*/
VCombineWPPhysFiles::VCombineWPPhysFiles( string iOutFileName )
{
    // list of files with different multiplicities
    fPhysFile.push_back( new VWPPhysFile( "DESY.d20161128.V2", "NIM2MST4LST4SST4", 4, 4, 4 ) );
    fPhysFile.push_back( new VWPPhysFile( "DESY.d20161128.V2", "NIM2MST3LST3SST3", 3, 3, 3 ) );
    fPhysFile.push_back( new VWPPhysFile( "DESY.d20161128.V2", "NIM2MST2LST2SST2", 2, 2, 2 ) );
    
    fOutFileName = iOutFileName;
    fOutFile = 0;
}

/*
 * open the given set of WP Phys files
 * and load histograms
 *
 */
bool VCombineWPPhysFiles::readWPPhysFiles( string iDatadirectory, string iSite, string iSubArray, int iObsTime, int iAnaID, string iPointingDirection )
{
    // loop over all input files and read
    // the histograms into memory
    for( unsigned int i = 0; i < fPhysFile.size(); i++ )
    {
        if( !fPhysFile[i]->open( iDatadirectory, iSite, iSubArray, iObsTime, iAnaID, iPointingDirection ) )
        {
            cout << "Error reading sensitivity histograms" << endl;
        }
    }
    
    ////////////////////////////////////
    // open output file
    // (all combined and smoothed histograms
    //  will be written into this file)
    ostringstream iOFileName;
    iOFileName << fOutFileName;
    iOFileName << ".ID" << iAnaID;
    iOFileName << iPointingDirection;
    iOFileName << "." << iSite << ".";
    iOFileName << iSubArray << ".";
    iOFileName << iObsTime;
    iOFileName << "s.root";
    
    cout << "Outputfile will be " << iOFileName.str() << endl;
    
    fOutFile = new TFile( iOFileName.str().c_str(), "RECREATE" );
    if( fOutFile->IsZombie() )
    {
        return false;
    }
    
    // clone all histograms
    if( fPhysFile.size() > 0 )
    {
        fOutHistos.clone( fPhysFile[0]->fHistos );
    }
    
    
    return true;
}

/*
 *  find bins with best sensitivity and
 *  write new files
 *
 */
bool VCombineWPPhysFiles::combineWPPhysFiles()
{

    // check that first sensitivity histograms exists
    if( fPhysFile.size() == 0
            || fPhysFile[0]->fHistos->fHisto1F.find( "DiffSens" ) == fPhysFile[0]->fHistos->fHisto1F.end()
            || !fPhysFile[0]->fHistos->fHisto1F["DiffSens"] )
    {
        return false;
    }
    
    //////////////////////////////////////////////
    // loop over all energy bins in sensitivity histograms
    // note: all what matters here is sensitivity
    for( int b = 0;
            b < fPhysFile[0]->fHistos->fHisto1F["DiffSens"]->GetNbinsX();
            b++ )
    {
        //////////////////////////////////////////////////
        // map to sort sensitivities by lowest values
        // map< sensitivity, file index >
        map< float, unsigned int > i_Sens;
        
        // loop over all phys files and fill map:
        // .first: sensitivity
        // .second: index of histogram vector
        for( unsigned int i = 0; i < fPhysFile.size(); i++ )
        {
            if( fPhysFile[i] && fPhysFile[i]->fHistos->fHisto1F["DiffSens"] )
            {
                if( fPhysFile[i]->fHistos->fHisto1F["DiffSens"]->GetBinContent( b ) > 0. )
                {
                    i_Sens[fPhysFile[i]->fHistos->fHisto1F["DiffSens"]->GetBinContent( b )] = i;
                }
            }
        }
        //////////////////////////////////////////////////
        // fill minimum sensitivity into the histograms
        if( i_Sens.begin()->second < fPhysFile.size() )
        {
            fOutHistos.fill( b, fPhysFile[i_Sens.begin()->second]->fHistos->fHisto1F["DiffSens"]->GetBinLowEdge( b ),
                             fPhysFile[i_Sens.begin()->second]->fHistos->fHisto1F["DiffSens"]->GetBinLowEdge( b + 1 ),
                             fPhysFile[i_Sens.begin()->second]->fHistos,
                             fPhysFile[i_Sens.begin()->second]->fMultiplicity_LST,
                             fPhysFile[i_Sens.begin()->second]->fMultiplicity_MST,
                             fPhysFile[i_Sens.begin()->second]->fMultiplicity_SST );
        }
    }
    
    return true;
}

/*
 * write all result files to disk
 *
 */
bool VCombineWPPhysFiles::writeWPPhysFiles()
{
    if( fOutFile )
    {
        fOutFile->cd();
        fOutHistos.write();
        fOutFile->Close();
    }
    return true;
}



/////////////////////////////////////////////////////
// main
/////////////////////////////////////////////////////


int main( int argc, char* argv[] )
{
    /////////////////////
    // input parameters
    if( argc != 5 )
    {
        cout << endl;
        cout << "./writeCTAWPPhysSensitivityOptimisedAndSmoothedFiles <list of sub array> <site> <output file> <data directory>" << endl;
        cout << endl;
        cout << "NOTE SEVERAL HARDWIRED PARAMETERS" << endl;
        cout << "search for all files in hard wired PHYS directory" << endl;
        cout << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    cout << endl;
    cout << "writeCTAWPPhysSensitivityOptimisedAndSmoothedFiles"  << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    string fSubArrayList = argv[1];
    string fSite = argv[2];
    string fOutputfile = argv[3];
    string fDataDirectory = argv[4];
    vector< string > fPointingDirection;
    fPointingDirection.push_back( "" );
    
    // hardwired observing time
    vector< int > iObservingTimeVector;
    iObservingTimeVector.push_back( 180000 );
    //    iObservingTimeVector.push_back( 18000 );
    iObservingTimeVector.push_back( 1800 );
    
    // analysis ID
    vector< int > iAnalysisIDVector;
    iAnalysisIDVector.push_back( 0 );
    
    ////////////////////////////////////////////////////
    // loop over all sub arrays and load files
    vector< string > fSubArray;
    
    ifstream is;
    is.open( fSubArrayList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading list of arrays from: " << endl;
        cout << fSubArrayList << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "reading list of arrays from " << endl;
    cout << fSubArrayList << endl;
    string iLine;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            fSubArray.push_back( iLine );
        }
    }
    is.close();
    
    cout << "total number of subarrays: " << fSubArray.size() << endl;
    
    /////////////////////////////////////////////////////////
    // combine sensitivities
    
    // array loop
    for( unsigned int i = 0; i < fSubArray.size(); i++ )
    {
        // observing time loop
        for( unsigned int t = 0; t < iObservingTimeVector.size(); t++ )
        {
            // analysis ID
            for( unsigned int d = 0; d < iAnalysisIDVector.size(); d++ )
            {
                // pointing direction
                for( unsigned int p = 0; p < fPointingDirection.size(); p++ )
                {
                    cout << endl << endl;
                    cout << "now filling array layout " << fSubArray[i];
                    cout << " (obs time " << iObservingTimeVector[t] << " s, ";
                    cout << " ID" <<  iAnalysisIDVector[d];
                    cout << ", Pointing " << fPointingDirection[p];
                    cout << ")" << endl;
                    cout << "=========================================" << endl;
                    cout << endl;
                    
                    VCombineWPPhysFiles a( fOutputfile );
                    
                    a.readWPPhysFiles( fDataDirectory,
                                       fSite,
                                       fSubArray[i],
                                       iObservingTimeVector[t],
                                       iAnalysisIDVector[d],
                                       fPointingDirection[p] );
                                       
                    a.combineWPPhysFiles();
                    
                    a.writeWPPhysFiles();
                    
                }
            }
        }
    }
    
    return true;
}

