/*
    small macro which reads a prod3 list file from KB and writes array layout files

    (e.g. Prod3-Layouts-new2.lis)

*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

/*
 * write a telescope line
 *
 * resolve'-' as from-to
 */
void writeTelescope( string iTemp2, string iText, ofstream& os )
{
    // 'from-to'
    if( iTemp2.find( "-" ) < iTemp2.size() )
    {
        unsigned int iStart = atoi( iTemp2.substr( 0, iTemp2.find( "-" ) ).c_str() );
        unsigned int iStopp = atoi( iTemp2.substr( iTemp2.find( "-" ) + 1, iTemp2.size() ).c_str() );
        for( int f = iStart; f <= iStopp; f++ )
        {
            os << f << iText << endl;
        }
    }
    // one line per telescope
    else
    {
        os << iTemp2 << iText << endl;
    }
}

void writeProd3_array( string iKB_file, string iArrayLayout, string iMST, string iSST, string prod3 )
{
    ///////////////////////////////////////////////
    // input:
    // open KB-style file with array layouts
    ifstream is;
    is.open( iKB_file.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, KB-style input file not found: " << endl;
        cout << iKB_file << endl;
        return;
    }
    ///////////////////////////////////////////////
    // input:
    // prod3 lis file South:
    string iOFile;
    if( prod3 == "South" )
    {
        iOFile = "CTA.prod3S." + iArrayLayout + "-" + iMST + iSST + ".lis";
    }
    else if( prod3 == "South-HB9" )
    {
        iOFile = "CTA.prod3Sb." + iArrayLayout + "-" + iMST + iSST + ".lis";
    }
    else
    {
        iOFile = "CTA.prod3N." + iArrayLayout + "-" + iMST + ".lis";
    }
    cout << iOFile << endl;
    ofstream os;
    os.open( iOFile.c_str() );
    if( !os )
    {
        cout << "error opening output file: " << iOFile << endl;
        return;
    }
    cout << "Writing results to " << iOFile << endl;
    
    string is_line;
    string iTemp;
    string iTemp2;
    while( getline( is, is_line ) )
    {
        // empty line
        if( is_line.size() == 0 )
        {
            continue;
        }
        istringstream is_stream( is_line );
        is_stream >> iTemp;
        // '*' is a comment
        if( iTemp == "*" )
        {
            continue;
        }
        
        // new array layout
        if( iTemp == "Layout" )
        {
            is_stream >> iTemp;
            if( iTemp.size() > 0 )
            {
                if( prod3 != "South-HB9" )
                {
                    iTemp = iTemp.substr( 0, iTemp.size() - 1 );
                }
                if( iTemp == iArrayLayout )
                {
                    /////////////////////////////////////////
                    cout << "FOUND LAYOUT " << iTemp << endl;
                    bool iFoundLSTs = false;
                    bool iFoundMSts = false;
                    bool iFoundSSTs = false;
                    while( getline( is, is_line ) )
                    {
                        stringstream is_stream2( is_line );
                        // LSTs: (12 bit dynamical range)
                        is_stream2 >> iTemp;
                        if( iTemp == "LST:" && !iFoundLSTs )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.     12    1", os );
                            }
                            iFoundLSTs = true;
                        }
                        // FlashCam
                        else if( iMST == "F" && ( iTemp == "MST-FlashCam:" || iTemp == "MST:" ) && !iFoundMSts )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.", os );
                            }
                            iFoundMSts = true;
                        }
                        // NectarCam (12 bit dynamical range)
                        else if( iMST == "N" && ( iTemp == "MST-NectarCam:" || iTemp == "MST:" ) && !iFoundMSts )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.     12    1", os );
                            }
                            iFoundMSts = true;
                        }
                        // SSTs (no Astri)
                        else if( ( iSST == "G" && iTemp == "GCT:" )
                                 || ( iSST == "D" && iTemp == "SST-1M:" ) )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.", os );
                            }
                        }
                        // SST (Astri)
                        else if( iSST == "A" && iTemp == "ASTRI:" )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.     12    1", os );
                            }
                        }
                        // SCTs (15 deg FOV)
                        else if( iTemp == "SCT:" )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.", os );
                            }
                        }
                        // SSTs
                        else if( iTemp == "SST:" && !iFoundSSTs )
                        {
                            while( !is_stream2.eof() )
                            {
                                is_stream2 >> iTemp2;
                                writeTelescope( iTemp2, "    20.", os );
                            }
                            iFoundSSTs = true;
                        }
                        if( iTemp == "SST-1M:" )
                        {
                            break;
                        }
                    }
                }
            }
        }
    }
    is.close();
    os.close();
}

/*

   write all prod3 array layouts

   selections are:
   - South (2015/2016 production)
   - South-HB9 (July 2016 production)
   - North (2015/2016 production)

*/

void writeAllArrays( string prod3 = "South", string iArrayList = "" )
{
    if( iArrayList.size() == 0 )
    {
        if( prod3 == "South" )
        {
            iArrayList = "Prod3-Layouts-new2.lis";
        }
        else if( prod3 == "South-HB9" )
        {
            iArrayList = "Prod3b-Layouts.lis";
        }
        else
        {
            iArrayList = "Prod3-Layouts-North.20151218.lis";
        }
    }
    vector< string > ivMST;
    ivMST.push_back( "N" );
    ivMST.push_back( "F" );
    vector< string > ivSST;
    ivSST.push_back( "A" );
    ivSST.push_back( "G" );
    ivSST.push_back( "D" );
    
    // layouts
    vector< string > iA;
    if( prod3 == "South" )
    {
        iA.push_back( "3HB1" );
        iA.push_back( "3HB2" );
        iA.push_back( "3HB3" );
        iA.push_back( "3HB4" );
        
        iA.push_back( "3HS1" );
        
        iA.push_back( "3HD1" );
        iA.push_back( "3HD2" );
        iA.push_back( "3HD3" );
        
        iA.push_back( "3HI1" );
        
        iA.push_back( "3HT1" );
        
        iA.push_back( "3SS2" );
        iA.push_back( "3SS3" );
        iA.push_back( "3SD1" );
        iA.push_back( "3SD2" );
        iA.push_back( "3SD3" );
        iA.push_back( "3ST1" );
        
        iA.push_back( "3HE1" );
        iA.push_back( "3HE2" );
        iA.push_back( "3HE3" );
        iA.push_back( "3HE4" );
        iA.push_back( "3HE5" );
        iA.push_back( "3HE6" );
        iA.push_back( "3HE7" );
        
        iA.push_back( "3HF1" );
        iA.push_back( "3HF2" );
        iA.push_back( "3HF3" );
        iA.push_back( "3HF4" );
        iA.push_back( "3HF5" );
        iA.push_back( "3HF6" );
        iA.push_back( "3HF7" );
        
        iA.push_back( "3HG1" );
        iA.push_back( "3HG2" );
        iA.push_back( "3HG3" );
        iA.push_back( "3HG4" );
        iA.push_back( "3HG5" );
        iA.push_back( "3HG6" );
        iA.push_back( "3HG7" );
        
        iA.push_back( "3HK1" );
        iA.push_back( "3HK2" );
        iA.push_back( "3HK3" );
        iA.push_back( "3HK4" );
        iA.push_back( "3HK5" );
        iA.push_back( "3HK6" );
        iA.push_back( "3HK7" );
    }
    else if( prod3 == "South-HB9" )
    {
        iA.push_back( "3HB9" );
        iA.push_back( "3HB89" );
        /*            iA.push_back( "3HB8" );
                    iA.push_back( "3HB4-2" );
                    iA.push_back( "3HB1-2" );
                    iA.push_back( "3HB2-2" ); */
    }
    else
    {
        iA.push_back( "3AL4M15-1" );
        iA.push_back( "3AL4M15-2" );
        iA.push_back( "3AL4M15-3" );
        iA.push_back( "3AL4M15-4" );
        iA.push_back( "3AL4M15-5" );
        
        iA.push_back( "3AL3M15-1" );
        iA.push_back( "3AL3M15-2" );
        iA.push_back( "3AL3M15-3" );
        iA.push_back( "3AL3M15-4" );
        iA.push_back( "3AL3M15-5" );
        
        iA.push_back( "3AL2M15-1" );
        iA.push_back( "3AL2M15-2" );
        iA.push_back( "3AL2M15-3" );
        iA.push_back( "3AL2M15-4" );
        iA.push_back( "3AL2M15-5" );
        
        iA.push_back( "3BL4M13-1" );
        iA.push_back( "3BL4M13-2" );
        iA.push_back( "3BL4M13-3" );
        iA.push_back( "3BL4M13-4" );
        iA.push_back( "3BL4M13-5" );
        
        iA.push_back( "3BL3M13-1" );
        iA.push_back( "3BL3M13-2" );
        iA.push_back( "3BL3M13-3" );
        iA.push_back( "3BL3M13-4" );
        iA.push_back( "3BL3M13-5" );
        
        iA.push_back( "3BL2M13-1" );
        iA.push_back( "3BL2M13-2" );
        iA.push_back( "3BL2M13-3" );
        iA.push_back( "3BL2M13-4" );
        iA.push_back( "3BL2M13-5" );
        
        iA.push_back( "3CL4M09-1" );
        iA.push_back( "3CL4M09-2" );
        iA.push_back( "3CL4M09-3" );
        iA.push_back( "3CL4M09-4" );
        iA.push_back( "3CL4M09-5" );
        
        iA.push_back( "3CL3M09-1" );
        iA.push_back( "3CL3M09-2" );
        iA.push_back( "3CL3M09-3" );
        iA.push_back( "3CL3M09-4" );
        iA.push_back( "3CL3M09-5" );
        
        iA.push_back( "3CL2M09-1" );
        iA.push_back( "3CL2M09-2" );
        iA.push_back( "3CL2M09-3" );
        iA.push_back( "3CL2M09-4" );
        iA.push_back( "3CL2M09-5" );
        
        iA.push_back( "3DL4M05-1" );
        iA.push_back( "3DL4M05-2" );
        iA.push_back( "3DL4M05-3" );
        iA.push_back( "3DL4M05-4" );
        iA.push_back( "3DL4M05-5" );
        
        iA.push_back( "3DL3M05-1" );
        iA.push_back( "3DL3M05-2" );
        iA.push_back( "3DL3M05-3" );
        iA.push_back( "3DL3M05-4" );
        iA.push_back( "3DL3M05-5" );
        
        iA.push_back( "3DL2M05-1" );
        iA.push_back( "3DL2M05-2" );
        iA.push_back( "3DL2M05-3" );
        iA.push_back( "3DL2M05-4" );
        iA.push_back( "3DL2M05-5" );
        
        iA.push_back( "3ZL4M15-1" );
        iA.push_back( "3ZL4M15-2" );
        iA.push_back( "3ZL4M15-3" );
        iA.push_back( "3ZL4M15-4" );
        iA.push_back( "3ZL4M15-5" );
        
        iA.push_back( "3ZL3M15-1" );
        iA.push_back( "3ZL3M15-2" );
        iA.push_back( "3ZL3M15-3" );
        iA.push_back( "3ZL3M15-4" );
        iA.push_back( "3ZL3M15-5" );
        
        iA.push_back( "3ZL2M15-1" );
        iA.push_back( "3ZL2M15-2" );
        iA.push_back( "3ZL2M15-3" );
        iA.push_back( "3ZL2M15-4" );
        iA.push_back( "3ZL2M15-5" );
    }
    
    
    for( unsigned int a = 0; a < iA.size(); a++ )
    {
        for( unsigned int m = 0; m < ivMST.size(); m++ )
        {
            for( unsigned int s = 0; s < ivSST.size(); s++ )
            {
                writeProd3_array( iArrayList, iA[a], ivMST[m], ivSST[s], prod3 );
            }
        }
    }
}
