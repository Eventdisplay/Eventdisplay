/*!  \class VDeadChannelFinder
 *   \brief find channels which are dead or deliver inconsistent data
 *
 */

#include "VDeadChannelFinder.h"

VDeadChannelFinder::VDeadChannelFinder( int irunmode, unsigned int iTelID, bool iLowGain, bool isMC )
{
    fDebug = false;
    
    frunmode = irunmode;
    fTelID = iTelID;
    fLowGain = iLowGain;
    fIsMC = isMC;
    
    // default values for high gain channels
    if( !fLowGain )
    {
        fDEAD_ped_min = 5.0;
        fDEAD_ped_max = 40.;
        fDEAD_pedvar_min = 0.5;
        fDEAD_pedvar_max = 5.e3;
        fDEAD_peddev_min = 2.;
        fDEAD_peddev_max = 4.;
        fDEAD_gain_min = 0.625;
        fDEAD_gain_max = 1.6;
        fDEAD_gainvar_min = 1.e-2;
        fDEAD_gainvar_max = 20.e5;
        fDEAD_gaindev_min = -1.e2;
        fDEAD_gaindev_max = 2.;
        fDEAD_toffset_max = 20.;
        fDEAD_l1rates_min = 0.;
        fDEAD_l1rates_max = 1.e20;
        fDead_HVrms_min = 1.e20;
        fDead_HVrms_max = 1.e20;
        fDEAD_tracemax_min = -1000;
        fDEAD_tracemax_max = 500000;
    }
    else
    {
        fDEAD_ped_min = 0.1;
        fDEAD_ped_max = 50.;
        fDEAD_pedvar_min = 0.1;
        fDEAD_pedvar_max = 5.e3;
        fDEAD_peddev_min = 3.;
        fDEAD_peddev_max = 5.;
        fDEAD_gain_min = 0.625;
        fDEAD_gain_max = 1.6;
        fDEAD_gainvar_min = 1.e-2;
        fDEAD_gainvar_max = 20.e5;
        fDEAD_gaindev_min = -1.e2;
        fDEAD_gaindev_max = 2.;
        fDEAD_toffset_max = 20.;
        fDEAD_l1rates_min = 0.;
        fDEAD_l1rates_max = 1.e20;
        fDead_HVrms_min = 1.e20;
        fDead_HVrms_max = 1.e20;
        fDEAD_tracemax_min = -1000;
        fDEAD_tracemax_max = 50000;
    }
}


bool VDeadChannelFinder::readDeadChannelFile( string ifile )
{
    if( ifile.size() == 0 )
    {
        return false;
    }
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VDeadChannelFinder::readDeadChannelFile:";
        cout << "error while opening file with dead channel definition: " << ifile << endl;
        cout << "exiting...." << endl;
        exit( -1 );
    }
    cout << "Telescope " << fTelID + 1 << ": reading dead channel definition";
    if( fLowGain )
    {
        cout << " (low gain) ";
    }
    else
    {
        cout << " (high gain) ";
    }
    cout << "from: " << endl;
    cout << "Telescope " << fTelID + 1 << ": " << ifile << endl;
    string iLine;
    string iTemp;
    string iTemp2;
    int t_temp;
    
    while( getline( is, iLine ) )
    {
        // line without '*' in the beginning are ignored
        if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*" )
        {
            istringstream is_stream( iLine );
            is_stream >> iTemp;
            // telescope number
            is_stream >> iTemp;
            t_temp = atoi( iTemp.c_str() ) - 1;
            if( t_temp >= 0 && t_temp != ( int )fTelID )
            {
                continue;
            }
            // high or low gain values
            is_stream >> iTemp;
            if( iTemp == "HIGHGAIN" && fLowGain )
            {
                continue;
            }
            if( iTemp == "LOWGAIN" && !fLowGain )
            {
                continue;
            }
            
            is_stream >> iTemp;
            if( iTemp == "PEDESTAL" )
            {
                is_stream >> iTemp;
                fDEAD_ped_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_ped_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "PEDESTALVARIATION" )
            {
                is_stream >> iTemp;
                fDEAD_pedvar_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_pedvar_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "PEDESTALDEVIATION" )
            {
                is_stream >> iTemp;
                fDEAD_peddev_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_peddev_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "GAIN" )
            {
                is_stream >> iTemp;
                fDEAD_gain_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_gain_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "GAINVARIATION" )
            {
                is_stream >> iTemp;
                fDEAD_gainvar_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_gainvar_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "GAINDEVIATION" )
            {
                is_stream >> iTemp;
                fDEAD_gaindev_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_gaindev_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "TIMEOFFSET" )
            {
                is_stream >> iTemp;
                is_stream >> iTemp;
                fDEAD_toffset_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "L1RATES" )
            {
                is_stream >> iTemp;
                fDEAD_l1rates_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_l1rates_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "HVRMS" )
            {
                is_stream >> iTemp;
                fDead_HVrms_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDead_HVrms_max = atof( iTemp.c_str() );
            }
            else if( iTemp == "TRACEMAX" )
            {
                is_stream >> iTemp;
                fDEAD_tracemax_min = atof( iTemp.c_str() );
                is_stream >> iTemp;
                fDEAD_tracemax_max = atof( iTemp.c_str() );
            }
            else
            {
                cout << "unknown identifier: " << iTemp << endl;
            }
        }
    }
    
    return false;
}


void VDeadChannelFinder::printSummary()
{
}


void VDeadChannelFinder::printDeadChannelDefinition()
{
    // print dead channel definitions
    cout << "Telescope " << fTelID + 1 << ":";
    cout << " dead channel selection criteria: ";
    if( fLowGain )
    {
        cout << "(low gain):  ";
    }
    else
    {
        cout << "(high gain): ";
    }
    cout << endl;
    cout << "Telescope " << fTelID + 1 << ":";
    cout << "\tPed [" << fDEAD_ped_min << "," << fDEAD_ped_max << "], ";
    cout << "Pedvar [" << fDEAD_pedvar_min << "," << fDEAD_pedvar_max << "], ";
    cout << "Peddev [" << fDEAD_peddev_min << "," << fDEAD_peddev_max << "]" << endl;
    cout << "Telescope " << fTelID + 1 << ":";
    cout << "\tGain [" << fDEAD_gain_min << "," << fDEAD_gain_max << "], ";
    cout << "Gainvar [" << fDEAD_gainvar_min << "," << fDEAD_gainvar_max << "], ";
    cout << "Gaindev [" << fDEAD_gaindev_min << "," << fDEAD_gaindev_max << "], ";
    cout << "Toff [" << fDEAD_toffset_max << "]" << endl;
    cout << "Telescope " << fTelID + 1 << ":";
    cout << "\tL1rates [" << fDEAD_l1rates_min << "," << fDEAD_l1rates_max << "], ";
    cout << "HVrms [" << fDead_HVrms_min << "," << fDead_HVrms_max << "], ";
    cout << "TraceMax [" << fDEAD_tracemax_min << "," << fDEAD_tracemax_max << "]" << endl;
}


unsigned int VDeadChannelFinder::testPedestals( unsigned int ichannel, double iPeds )
{
    // test pedestal range
    if( iPeds < fDEAD_ped_min )
    {
        if( fDebug )
        {
            cout << "testPedestals (1): " << ichannel << " " << fDEAD_ped_min << " " << iPeds << endl;
        }
        return 1;
    }
    if( iPeds > fDEAD_ped_max )
    {
        if( fDebug )
        {
            cout << "testPedestals (1): " << ichannel << " " << fDEAD_ped_max << " " << iPeds << endl;
        }
        return 1;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testPedestalVariations( unsigned int ichannel, double iPedVar )
{
    // test pedvars (!=0)
    if( iPedVar <= fDEAD_pedvar_min )
    {
        if( fDebug )
        {
            cout << "testPedestalVariations (2): " << ichannel << " " << fDEAD_pedvar_min << " " << iPedVar << endl;
        }
        return 2;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testPedestalVariationsMinOut( unsigned int ichannel, double iPedVar, double imeanPed, double irmsPed )
{
    // test if pedvars is an outlier
    if( iPedVar < imeanPed - fDEAD_peddev_min * irmsPed )
    {
        if( fDebug )
        {
            cout << "testPedestalVariationsMinOut (3): ";
            cout << ichannel << " " << fDEAD_peddev_min << " ";
            cout << imeanPed - fDEAD_peddev_min* irmsPed << " " << imeanPed << " " << irmsPed << " " << iPedVar << endl;
        }
        return 3;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testPedestalVariationsMaxOut( unsigned int ichannel, double iPedVar, double imeanPed, double irmsPed )
{
    // test if pedvars is an outlier
    if( iPedVar > imeanPed + fDEAD_peddev_max * irmsPed )
    {
        if( fDebug )
        {
            cout << "testPedestalVariationsMaxOut (4): " << ichannel << " " << fDEAD_peddev_max;
            cout << " " << imeanPed + fDEAD_peddev_max* irmsPed << " " << imeanPed << " ";
            cout << irmsPed << " " << iPedVar << endl;
        }
        return 4;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testGains( unsigned int ichannel, double iGain )
{
    if( fIsMC || frunmode == 2 || frunmode == 5 )
    {
        return 0;
    }
    
    if( iGain < fDEAD_gain_min )
    {
        if( fDebug )
        {
            cout << "testGains (5): " << ichannel << " " << fDEAD_gain_min << " " << iGain << endl;
        }
        return 5;
    }
    if( iGain > fDEAD_gain_max )
    {
        if( fDebug )
        {
            cout << "testGains: " << ichannel << " " << fDEAD_gain_max << " " << iGain << endl;
        }
        return 5;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testGainVariations( unsigned int ichannel, double iGainVar )
{
    if( fIsMC || frunmode == 2 || frunmode == 5 )
    {
        return 0;
    }
    
    if( iGainVar < fDEAD_gainvar_min )
    {
        if( fDebug )
        {
            cout << "testGainVariations (6): " << ichannel << " " << fDEAD_gainvar_min << " " << iGainVar << endl;
        }
        return 6;
    }
    if( iGainVar > fDEAD_gainvar_max )
    {
        if( fDebug )
        {
            cout << "testGainVariations (6): " << ichannel << " " << fDEAD_gainvar_min << " " << iGainVar << endl;
        }
        return 6;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testGainDev( unsigned int ichannel, double iGain, double iGainVar, bool iDefault )
{
    if( fIsMC || frunmode == 2 || frunmode == 5 )
    {
        return 0;
    }
    
    if( iGainVar > 0. && iGain / iGainVar < fDEAD_gaindev_max && fabs( iGain - 1. ) < 1.e-3 && !iDefault )
    {
        if( fDebug )
        {
            cout << "testGainDev (7): " << ichannel << " " << iGain << " " << iGainVar << endl;
        }
        return 7;
    }
    
    return 0;
}


unsigned int VDeadChannelFinder::testTimeOffsets( unsigned int ichannel, double iT )
{
    if( fIsMC || frunmode == 2 || frunmode == 5 )
    {
        return 0;
    }
    
    if( fabs( iT ) > fDEAD_toffset_max )
    {
        if( fDebug )
        {
            cout << "testTimeOffsets (8): " << ichannel << " " << fDEAD_toffset_max << " " << iT << endl;
        }
        return 8;
    }
    
    return 0;
}

/*
 * MNR: Channel is dead because we cannot identify the maximum of the trace (shape is off)
 */
unsigned int VDeadChannelFinder::testTraceMax( unsigned int ichannel, double iTraceMax )
{
    if( fIsMC || frunmode == 2 || frunmode == 5 )
    {
        return 0;
    }
    
    if( iTraceMax < fDEAD_tracemax_min )
    {
        if( fDebug )
        {
            cout << "testTraceMax (15): " << ichannel << " " << fDEAD_tracemax_min << " " << iTraceMax << endl;
        }
        return 15;
    }
    if( iTraceMax > fDEAD_tracemax_max )
    {
        if( fDebug )
        {
            cout << "testTraceMax (15): " << ichannel << " " << fDEAD_tracemax_min << " " << iTraceMax << endl;
        }
        return 15;
    }
    return 0;
}

