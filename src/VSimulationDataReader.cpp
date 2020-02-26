/*! \class VSimulationDataReader
    \brief  reads simulations packages from vbf files

    workaround to get simulation data from vbf file format

*/

#include "VSimulationDataReader.h"

VSimulationDataReader::VSimulationDataReader()
{
    fDebug = false;
    fSimulationData = 0;
    
    fEventNumber = 0;
    fRunNumber = 0;
    fMCprimary = 0;
    fMCenergy = 0.;
    fMCX = 0.;
    fMCY = 0.;
    fMCXcos = 0.;
    fMCYcos = 0.;
    fMCZe = 0.;
    fMCAz = 0.;
    fMCXoff = 0.;
    fMCYoff = 0.;
    fTel_Elevation = 0.;
    fTel_Azimuth = 0.;
    
    fShowerID = -1;
    fFirstInteractionHeight = -1;
    fFirstInteractionDepth = -1;
    fCorsikaRunID = -1;
}


/*!
     read vbf simulation bank using VRawEventData::getPacket() method
*/
bool VSimulationDataReader::setSimulationData( VRawEventData* iraweventData )
{
    return setSimulationData( iraweventData->getPacket() );
}


/*!
    print simulation header
*/
bool VSimulationDataReader::printSimulationHeader( VPacket* packet, int bPrintCFG )
{
    if( !packet )
    {
        return false;
    }
    
    VSimulationHeader* h = packet->getSimulationHeader();
    
    if( h )
    {
        cout << "================================================================================" << endl;
        cout << "Simulations from " << h->fDateOfSimsUTC;
        cout << ", simulated by " << h->fSimulator;
        cout << ", simulation package " << ( int )h->fSimulationPackage;
        cout << ", array configurations from " << h->fDateOfArrayForSims << endl;
        cout << "\t atmospheric model : " << h->fAtmosphericModel << " (possibly not filled yet)" << endl;
        cout << "\t array configuration: " << endl;
        for( unsigned int i = 0; i < h->fArray.size(); i++ )
        {
            cout << "\t\t telescope " << i + 1 << ": " << h->fArray[i].fRelTelLocSouthM;
            cout << "\t" << h->fArray[i].fRelTelLocEastM << "\t" << h->fArray[i].fRelTelLocUpM << endl;
        }
        // this might print the complete .cfg file
        if( h->fSimConfigFile.size() && bPrintCFG == 1 )
        {
            cout << h->fSimConfigFile << endl;
        }
        // find config file name. The line should look like this: "  (cfg file /path/to/file)\n"
        if( h->fSimConfigFile.size() && bPrintCFG == 2 )
        {
            size_t index_start = h->fSimConfigFile.find( "(cfg file" );
            if( index_start != string::npos )
            {
                size_t index_end = h->fSimConfigFile.find( "\n", index_start );
                if( index_end != string::npos )
                {
                    string filename = h->fSimConfigFile.substr( index_start + 10, index_end - 1 - index_start - 10 );
                    if( filename.size() )
                    {
                        cout << "Grisu configuration file: " << filename << endl;
                    }
                }
            }
        }
        cout << "================================================================================" << endl;
    }
    
    return true;
}

VMonteCarloRunHeader* VSimulationDataReader::fillSimulationHeader( VPacket* packet )
{
    if( !packet )
    {
        return 0;
    }
    
    VSimulationHeader* h = packet->getSimulationHeader();
    if( !h )
    {
        return 0;
    }
    
    // new MC run header
    VMonteCarloRunHeader* iMCRunHeader = new VMonteCarloRunHeader();
    iMCRunHeader->SetName( "MC_runheader" );
    
    iMCRunHeader->runnumber = h->fRunNumber;
    iMCRunHeader->detector_Simulator = h->fSimulator;
    iMCRunHeader->detector_date = h->fDateOfSimsUTC;
    // CORSIKA/grisudet
    if( ( int )h->fSimulationPackage == 1 )
    {
        iMCRunHeader->shower_prog_id = 1;
        iMCRunHeader->detector_prog_id = 2;
    }
    iMCRunHeader->atmosphere = ( int )h->fAtmosphericModel;
    iMCRunHeader->obsheight = ( double )h->fObsAltitudeM;
    
    // read long string of fSimConfigFile and extract all the necessary information
    // (very dependent on structure of string)
    if( h->fSimConfigFile.size() > 0 )
    {
        istringstream is_stream( h->fSimConfigFile );
        string iTemp = "";
        
        while( !(is_stream>>std::ws).eof() )
        {
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iTemp;
            }
            
            if( iTemp == "corsikaIOreader" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp;
                }
                iMCRunHeader->converter_prog_vers = ( int )( atof( iTemp.substr( 0, iTemp.size() - 1 ).c_str() ) * 1000. );
            }
            else if( iTemp == "DATE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->shower_date;
                }
            }
            else if( iTemp == "CORSIKAVERSION" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp;
                }
                iMCRunHeader->shower_prog_vers = atoi( iTemp.c_str() ) * 1000;
            }
            else if( iTemp == "PARTICLEID" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->primary_id;
                }
            }
            else if( iTemp == "OBSLEVEL" )
            {
                float i_f = 0.;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> i_f;
                }
                if( TMath::Abs( i_f - iMCRunHeader->obsheight ) > 10. )
                {
                    cout << "VSimulationDataReader::fillSimulationHeader warning: different observation level in runheader and simconfig string ";
                    cout << fixed << iMCRunHeader->obsheight << "\t" << i_f << endl;
                }
            }
            else if( iTemp == "ALTITUDE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->injection_height;
                }
            }
            else if( iTemp == "TSTART" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->fixed_int_depth;
                }
            }
            else if( iTemp == "E_SLOPE" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->spectral_index;
                }
            }
            else if( iTemp == "E_MIN" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->E_range[0];
                }
            }
            else if( iTemp == "E_MAX" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->E_range[1];
                }
            }
            else if( iTemp == "ZENITH_MIN" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->alt_range[0];
                }
                iMCRunHeader->alt_range[0] = ( 90. - iMCRunHeader->alt_range[0] ) / ( 45. / atan( 1. ) );
            }
            else if( iTemp == "ZENITH_MAX" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->alt_range[1];
                }
                iMCRunHeader->alt_range[1] = ( 90. - iMCRunHeader->alt_range[1] ) / ( 45. / atan( 1. ) );
            }
            else if( iTemp == "AZIMUTH_MIN" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->az_range[0];
                }
                iMCRunHeader->az_range[0] /= ( 45. / atan( 1. ) );
            }
            else if( iTemp == "AZIMUTH_MAX" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->az_range[1];
                }
                iMCRunHeader->az_range[1] /= ( 45. / atan( 1. ) );
            }
            else if( iTemp == "VIEWCONE_MIN" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->viewcone[0];
                }
                iMCRunHeader->viewcone[0] /= ( 45. / atan( 1. ) );
            }
            else if( iTemp == "VIEWCONE_MAX" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->viewcone[1];
                }
                iMCRunHeader->viewcone[1] /= ( 45. / atan( 1. ) );
            }
            else if( iTemp == "NSCATT" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->num_use;
                }
            }
            else if( iTemp == "XSCATT" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->core_range[0];
                }
            }
            else if( iTemp == "YSCATT" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->core_range[1];
                }
            }
            else if( iTemp == "CBUNCH" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_bunchsize;
                }
            }
            else if( iTemp == "C_WMIN" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_wlen_min;
                }
            }
            else if( iTemp == "C_WMAX" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_wlen_max;
                }
            }
            else if( iTemp == "GEO_X" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->B_inclination;    // entry will be corrected later
                }
            }
            else if( iTemp == "GEO_Z" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->B_total;    // entry will be corrected later
                }
            }
            else if( iTemp == "GEO_A" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->B_declination;
                }
                iMCRunHeader->B_declination /= ( 45. / atan( 1. ) );
            }
            else if( iTemp == "HAD_LOW" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_low_E_model;
                }
            }
            else if( iTemp == "HAD_HIGH" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_high_E_model;
                }
            }
            else if( iTemp == "HAD_TRANS" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_low_high_E;
                }
            }
            else if( iTemp == "CFLAG" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iMCRunHeader->corsika_iact_options;
                }
            }
            else if( iTemp == "ATM" )
            {
                int iATM = 0;
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp;
                }
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iATM;
                }
                // simulation header in MC vbf file is not filled properly (yet)
                if( iATM != 0 )
                {
                    iMCRunHeader->atmosphere = iATM;
                }
            }
            else if( iTemp == "*" )
            {
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iTemp;
                }
                if( iTemp == "SOURC" )
                {
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iTemp;
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iTemp;
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iTemp;
                    }
                    if( atof( iTemp.c_str() ) > 0. )
                    {
                        iMCRunHeader->diffuse = 1;
                        iMCRunHeader->viewcone[0] = 0.;
                        iMCRunHeader->viewcone[1] = atof( iTemp.c_str() );
                    }
                }
                // read low gain multiplier
                else if( iTemp == "FADCS" )
                {
                    for( int kk = 0; kk < 7; kk ++ )
                    {
                        if( !(is_stream>>std::ws).eof() )
                        {
                            is_stream >> iTemp;
                        }
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> iTemp;
                        iMCRunHeader->fFADC_hilo_multipler = atof( iTemp.c_str() );
                    }
                    else
                    {
                        iMCRunHeader->fFADC_hilo_multipler = -999.;
                    }
                }
            }
            
        }
    }
    // scatter mode (y=0)
    if( TMath::Abs( iMCRunHeader->core_range[1] ) < 1.e-5 )
    {
        iMCRunHeader->core_pos_mode = 1;
    }
    // geomagnetic field
    iMCRunHeader->B_total = sqrt( iMCRunHeader->B_total * iMCRunHeader->B_total + iMCRunHeader->B_inclination * iMCRunHeader->B_inclination );
    if( iMCRunHeader->B_total > 0. )
    {
        iMCRunHeader->B_inclination = acos( iMCRunHeader->B_inclination / iMCRunHeader->B_total );
    }
    else
    {
        iMCRunHeader->B_inclination = 0.;
    }
    
    return iMCRunHeader;
}


/*!
     read vbf simulation bank directly from a packet
*/
bool VSimulationDataReader::setSimulationData( VPacket* packet )
{
    if( fDebug )
    {
        cout << "\t\t VSimulationDataReader::setSimulationData" << endl;
    }
    bool iReturn = false;
    
    if( !packet )
    {
        return iReturn;
    }
    
    fMCprimary = 0;
    fMCenergy = 0.;
    fMCX = 0.;
    fMCY = 0.;
    fMCXcos = 0.;
    fMCYcos = 0.;
    fMCZe = 0.;
    fMCAz = 0.;
    fShowerID = 0;
    fFirstInteractionHeight = 0;
    fFirstInteractionDepth = 0;
    fCorsikaRunID = 0;
    
    // get local trigger information
    // this should be somewhere else, since trigger info is
    // no simulation data (preliminary)
    uint32_t eventNumber = 0;
    if( packet->hasArrayEvent() )
    {
        VArrayEvent* i_ae =  packet->getArrayEvent();
        // this vector has size 255
        fLocalTrigger = i_ae->getExpectedTelescopes();
        if( i_ae->hasEventNumber() )
        {
            eventNumber = i_ae->getEventNumber();
        }
        else
        {
            eventNumber = 0;
        }
    }
    else
    {
        fLocalTrigger.assign( 255, false );
    }
    
    if( packet->hasSimulationData() )
    {
        if( fDebug )
        {
            cout << "\t\t\t VSimulationDataReader::setSimulationData(): hasSimulationData(), eventnumber " << eventNumber << endl;
        }
        fSimulationData = packet->getSimulationData();
        ///////////////////////////////////////////////////////////////////////
        // VBF is currently not backward compatible (d20061130)
        // VBF doesn not show any version change from compatible version
        // (compile with -DVBF_021)
#ifdef VBF_021
        fMCprimary = ( int )fSimulationData->id;
        fMCenergy = fSimulationData->e;
        fMCX = fSimulationData->r[0];
        fMCY = fSimulationData->r[1];
        fMCZe = fSimulationData->theta;
        fMCAz = fSimulationData->phi;
        fTel_Elevation = 0.;
        fTel_Azimuth = 0.;
        ////////////////////////////////////////////////////////////////////////////////////////////
        // newer versions
        ////////////////////////////////////////////////////////////////////////////////////////////
#else
        fRunNumber = fSimulationData->fRunNumber;
        fEventNumber = fSimulationData->fEventNumber;
        fMCprimary = ( int )fSimulationData->fCORSIKAParticleID;
        // eventdisplay works in [TeV]
        fMCenergy = fSimulationData->fEnergyGeV / 1.e3;
        fMCX = fSimulationData->fCoreEastM;
        // grisu coordinate system is with y towards north
        fMCY = -1.*fSimulationData->fCoreSouthM;
        fMCZe = fSimulationData->fPrimaryZenithDeg;
        fMCAz = fSimulationData->fPrimaryAzimuthDeg;
        fTel_Elevation = 90. - fSimulationData->fObservationZenithDeg;
        fTel_Azimuth = fSimulationData->fObservationAzimuthDeg;
        
        float x = 0.;
        float y = 0.;
        float z = 0.;
        VSkyCoordinatesUtilities::getDifferenceInCameraCoordinates( fSimulationData->fObservationZenithDeg, fSimulationData->fObservationAzimuthDeg,
                fSimulationData->fPrimaryZenithDeg, fSimulationData->fPrimaryAzimuthDeg, x, y, z );
        fMCXoff = x;
        fMCYoff = y;
#endif
        // calculate direction cosinii
        // (this might be different to what VGrisuReader reads from the MC file
        float degrad = 45. / atan( 1. );
        fMCXcos = sin( ( fMCZe ) / degrad ) * cos( ( fMCAz ) / degrad );
        // rounding error
        if( fabs( fMCXcos ) < 1.e-8 )
        {
            fMCXcos = 0.;
        }
        fMCYcos = sin( ( fMCZe ) / degrad ) * sin( ( fMCAz ) / degrad );
        // rounding error
        if( fabs( fMCYcos ) < 1.e-8 )
        {
            fMCYcos = 0.;
        }
        iReturn = true;
        if( fDebug )
        {
            cout << "\t\t\t VSimulationDataReader::setSimulationData(): found shower in run " << fRunNumber << " , event " << fEventNumber;
            cout << " with " << fMCprimary << " " << fMCenergy << " " << fMCZe << " " << fMCAz << " " << fMCX << " " << fMCY <<  endl;
            cout << "\t\t\t    (telescope pointing: " << setprecision( 6 ) << 90. - fTel_Elevation << " ze, " << fTel_Azimuth << " az)" << endl;
            cout << "\t\t\t    (xoff, yoff) = " << fMCXoff << "\t" << fMCYoff << endl;
        }
    }
    else
    {
        if( fDebug )
        {
            cout << "\t\t\t VSimulationDataReader::setSimulationData(): no hasSimulationData()" << endl;
        }
        iReturn = false;
    }
#ifdef VBF_034
    if( packet->hasCorsikaSimulationData() )
    {
        if( fDebug )
        {
            cout << "\t\t\t VSimulationDataReader::setSimulationData(): hasCorsikaSimulationData(), eventnumber " << eventNumber << endl;
        }
        fCorsikaSimulationData = packet->getCorsikaSimulationData();
        fShowerID = fCorsikaSimulationData->fShowerID;
        fCorsikaRunID = fCorsikaSimulationData->fCorsikaRunID;
        fFirstInteractionHeight = fCorsikaSimulationData->fFirstInteractionHeight;
        fFirstInteractionDepth = fCorsikaSimulationData->fFirstInteractionDepth;
    }
    else
    {
        if( fDebug )
        {
            cout << "\t\t\t VSimulationDataReader::setSimulationData(): no hasCorsikaSimulationData()" << endl;
        }
    }
#endif
    return iReturn;
}
