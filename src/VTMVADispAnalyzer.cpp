/* \class VTMVADispAnalyzer
   \brief TMVA based disp analysis

   angular disp, energy and core reconstruction

*/

#include "VTMVADispAnalyzer.h"

VTMVADispAnalyzer::VTMVADispAnalyzer( string iFile, vector<ULong64_t> iTelTypeList, string iDispType )
{
    fDebug = false;
    
    fDispType = iDispType;
    
    // for a valid analysis, zombie should be false
    bZombie = true;
    
    // list of variables used in TMVA disp method
    fWidth = 0.;
    fLength = 0.;
    fWoL = 0.;
    fSize = 0.;
    fNtubes = 0.;
    fPedvar = 0.;
    fAsymm = 0.;
    fTGrad = 0.;
    fZe = 0.;
    fAz = 0.;
    fLoss = 0.;
    fXcore = 0.;
    fYcore = 0.;
    fcross = 0.;
    fDist = 0.;
    fFui = 0.;
    fRcore = 0.;
    fEHeight = 0.;
    
    // spectators (nowhere used)
    float iMCe0 = 0.;
    float iMCxoff = 0.;
    float iMCyoff = 0.;
    float iMCxcore = 0.;
    float iMCycore = 0.;
    float iMCrcore = 0.;
    float iNImages = 0.;
    float cen_x = 0;
    float cen_y = 0;
    float cosphi = 0.;
    float sinphi = 0.;
    float temp1 = 0.;
    float temp2 = 0.;
    
    // list of telescope types: required to selected correct BDT weight file
    fTelescopeTypeList = iTelTypeList;
    
    if( fTelescopeTypeList.size() == 0 )
    {
        cout << "VTMVADispAnalyzer initializion error: telescope type list length is zero" << endl;
        bZombie = true;
        return;
    }
    
    cout << endl;
    cout << "New TMVA reader for ";
    if( fDispType == "BDTDispEnergy" )
    {
        cout << "disp energy reconstruction";
    }
    else if( fDispType == "BDTDispCore" )
    {
        cout << "disp core reconstruction";
    }
    else if( fDispType == "BDTDispError" )
    {
        cout << "disp angular uncertainty estimation";
    }
    else
    {
        cout << "disp angular reconstruction";
    }
    cout << endl;
    cout << "===============================================" << endl << endl;
    
    // initialize TMVA readers
    // (one per telescope type)
    for( unsigned int i = 0; i < fTelescopeTypeList.size(); i++ )
    {
        ostringstream iFileName;
        iFileName << iFile << fTelescopeTypeList[i] << ".weights.xml";
        cout << "initializing TMVA disp analyzer: " <<  iFileName.str() << endl;
        // check that TMVA file exists
        ifstream i_temp_TMVAFILE( iFileName.str().c_str() );
        if( !i_temp_TMVAFILE.good() )
        {
            cout << "VTMVADispAnalyzer error: cannot find: " << endl;
            cout << iFileName.str() << endl;
            bZombie = true;
            return;
        }
        bool iSingleTelescopeAnalysis = true;
        // try to detect if this is a single telescope analysis
        string iLine;
        if( i_temp_TMVAFILE.is_open() )
        {
            while( getline(i_temp_TMVAFILE,iLine) )
            {
                if( iLine.find("cross") != string::npos )
                { 
                   iSingleTelescopeAnalysis = false;
                   break;
                }
            }
        }
        if( iSingleTelescopeAnalysis )
        {
            cout << "\t single-telescope disp analysis" << endl;
        }
        else
        {
            cout << "\t multi-telescope disp analysis" << endl;
        }

        fTMVAReader[fTelescopeTypeList[i]] = new TMVA::Reader( "!Color:!Silent" );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "width", &fWidth );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "length", &fLength );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "wol", &fWoL );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "size", &fSize );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "ntubes", &fNtubes );
        // ASTRI telescopes are without timing information
        // tgrad_x is therefore ignored
        if( fTelescopeTypeList[i] != 201511619 )
        {
            fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "tgrad_x*tgrad_x", &fTGrad );
        }
        // cross variable should be on this spot
        if( !iSingleTelescopeAnalysis )
        {
            fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "cross", &fcross );
        }
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "asym", &fAsymm );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "loss", &fLoss );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "dist", &fDist );
        fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "fui", &fFui );
        if( fDispType == "BDTDispEnergy" && !iSingleTelescopeAnalysis )
        {
            fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "EHeight", &fEHeight );
            fTMVAReader[fTelescopeTypeList[i]]->AddVariable( "Rcore", &fRcore );
        }
        // spectators
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "cen_x", &cen_x );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "cen_y", &cen_y );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "cosphi", &cosphi);
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "sinphi", &sinphi);
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCe0", &iMCe0 );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCxoff", &iMCxoff );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCyoff", &iMCyoff );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCxcore", &iMCxcore );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCycore", &iMCycore );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "MCrcore", &iMCrcore );
        fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "NImages", &iNImages );
        if( fDispType == "BDTDisp" )
        {
            fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "dispError", &temp2 );
            fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "dispPhi", &temp1 );
        }
        else if( fDispType == "BDTDispError" )
        {
            fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "disp", &temp2 );
            fTMVAReader[fTelescopeTypeList[i]]->AddSpectator( "dispPhi", &temp1 );
        }
        
        if( !fTMVAReader[fTelescopeTypeList[i]]->BookMVA( "BDTDisp", iFileName.str().c_str() ) )
        {
            cout << "VTMVADispAnalyzer initializion error: xml weight file not found:" << endl;
            cout << "\t" << iFileName.str() << endl;
            bZombie = true;
            return;
        }
    }
    bZombie = false;
}

/*
 * calculate disp using the TMVA BDTs
 *
*/
float VTMVADispAnalyzer::evaluate( float iWidth, float iLength, float iSize, float iAsymm, float iLoss, float iTGrad,
                                   float icen_x, float icen_y, float xoff_4, float yoff_4, ULong64_t iTelType,
                                   float iZe, float iAz, float iRcore, float iEHeight, float iDist, float iFui, float iNtubes )
{
    fWidth = iWidth;
    fLength = iLength;
    if( fLength > 0. )
    {
        fWoL = fWidth / fLength;
    }
    else
    {
        fWoL = 0.;
    }
    if( iSize > 0. )
    {
        fSize = log10( iSize );
    }
    else
    {
        return -99.;
    }
    if( iNtubes > 0. )
    {
        fNtubes = log10( iNtubes );
    }
    else
    {
       return -99.;
    }
    fTGrad = iTGrad * iTGrad;
    fZe = iZe;
    fAz = iAz;
    // %%@^#%!@(#
    // fcross = sqrt( ( icen_y - yoff_4 ) * ( icen_y - yoff_4 ) + ( icen_x - xoff_4 ) * ( icen_x - xoff_4 ) );
    if( yoff_4 > -999. && xoff_4 > -999. )
    {
        fcross = sqrt( ( icen_y + yoff_4 ) * ( icen_y + yoff_4 ) + ( icen_x - xoff_4 ) * ( icen_x - xoff_4 ) );
    }
    else
    {
        fcross = 0.;
    }
    
    fLoss = iLoss;
    fAsymm = iAsymm;
    fRcore = iRcore;
    fEHeight = iEHeight;
    fDist = iDist;
    fFui  = iFui;
    
    if( fTMVAReader.find( iTelType ) != fTMVAReader.end() && fTMVAReader[iTelType] )
    {
        return ( fTMVAReader[iTelType]->EvaluateRegression( "BDTDisp" ) )[0];
    }
    
    return -99.;
}

void VTMVADispAnalyzer::terminate()
{
    return;
}
