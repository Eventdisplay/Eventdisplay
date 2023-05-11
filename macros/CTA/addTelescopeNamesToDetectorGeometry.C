/*
 * script to add fourth column to CTA detector geometry files
 * listing the observatory telescopes number
 *
 * Note: incomplete, only F and G are encrypted
 *
*/

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

string getTelescopeName( int iTelID )
{
    map< int, string > fTelNameMap;
    
    fTelNameMap[5] = "L-01";
    fTelNameMap[6] = "L-02";
    fTelNameMap[4] = "L-03";
    fTelNameMap[11] = "L-05";
    
    fTelNameMap[12] = "M-01";
    fTelNameMap[14] = "M-02";
    fTelNameMap[30] = "M-03";
    fTelNameMap[13] = "M-04";
    fTelNameMap[29] = "M-05";
    fTelNameMap[23] = "M-06";
    fTelNameMap[15] = "M-07";
    fTelNameMap[16] = "M-08";
    fTelNameMap[24] = "M-09";
    fTelNameMap[26] = "M-10";
    fTelNameMap[19] = "M-11";
    fTelNameMap[20] = "M-12";
    fTelNameMap[27] = "M-13";
    fTelNameMap[47] = "M-14";
    fTelNameMap[31] = "M-15";
    fTelNameMap[25] = "M-16";
    fTelNameMap[32] = "M-17";
    fTelNameMap[48] = "M-18";
    fTelNameMap[50] = "M-19";
    fTelNameMap[33] = "M-20";
    fTelNameMap[28] = "M-21";
    fTelNameMap[34] = "M-22";
    fTelNameMap[51] = "M-23";
    fTelNameMap[49] = "M-24";
    fTelNameMap[52] = "M-25";
    
    fTelNameMap[257] = "S-01";
    fTelNameMap[258] = "S-02";
    fTelNameMap[259] = "S-03";
    fTelNameMap[260] = "S-04";
    fTelNameMap[268] = "S-05";
    fTelNameMap[269] = "S-06";
    fTelNameMap[274] = "S-07";
    fTelNameMap[275] = "S-08";
    fTelNameMap[284] = "S-09";
    fTelNameMap[285] = "S-10";
    fTelNameMap[290] = "S-11";
    fTelNameMap[291] = "S-12";
    fTelNameMap[292] = "S-13";
    fTelNameMap[293] = "S-14";
    fTelNameMap[280] = "S-15";
    fTelNameMap[281] = "S-16";
    fTelNameMap[282] = "S-17";
    fTelNameMap[283] = "S-18";
    fTelNameMap[316] = "S-19";
    fTelNameMap[317] = "S-20";
    fTelNameMap[300] = "S-21";
    fTelNameMap[301] = "S-22";
    fTelNameMap[302] = "S-23";
    fTelNameMap[303] = "S-24";
    fTelNameMap[324] = "S-25";
    fTelNameMap[325] = "S-26";
    fTelNameMap[327] = "S-27";
    fTelNameMap[328] = "S-28";
    fTelNameMap[322] = "S-29";
    fTelNameMap[323] = "S-30";
    fTelNameMap[342] = "S-31";
    fTelNameMap[343] = "S-32";
    fTelNameMap[348] = "S-33";
    fTelNameMap[349] = "S-34";
    fTelNameMap[366] = "S-35";
    fTelNameMap[367] = "S-36";
    fTelNameMap[344] = "S-37";
    fTelNameMap[345] = "S-38";
    fTelNameMap[350] = "S-39";
    fTelNameMap[351] = "S-40";
    fTelNameMap[346] = "S-41";
    fTelNameMap[347] = "S-42";
    fTelNameMap[352] = "S-43";
    fTelNameMap[353] = "S-44";
    fTelNameMap[368] = "S-45";
    fTelNameMap[369] = "S-46";
    fTelNameMap[370] = "S-47";
    fTelNameMap[371] = "S-48";
    fTelNameMap[380] = "S-49";
    fTelNameMap[381] = "S-50";
    fTelNameMap[384] = "S-51";
    fTelNameMap[385] = "S-52";
    fTelNameMap[382] = "S-53";
    fTelNameMap[383] = "S-54";
    fTelNameMap[386] = "S-55";
    fTelNameMap[387] = "S-56";
    fTelNameMap[378] = "S-57";
    fTelNameMap[379] = "S-58";
    fTelNameMap[392] = "S-59";
    fTelNameMap[393] = "S-60";
    fTelNameMap[396] = "S-61";
    fTelNameMap[397] = "S-62";
    fTelNameMap[394] = "S-63";
    fTelNameMap[395] = "S-64";
    fTelNameMap[398] = "S-65";
    fTelNameMap[399] = "S-66";
    fTelNameMap[400] = "S-67";
    fTelNameMap[401] = "S-68";
    fTelNameMap[402] = "S-69";
    fTelNameMap[403] = "S-70";
    
    if( fTelNameMap.find( iTelID ) != fTelNameMap.end() )
    {
        return fTelNameMap[iTelID];
    }
    
    return "NOT-FOUND";
}


void addTelescopeNamesToDetectorGeometry( string iInputFile, string iOutputFile )
{
    // read list of telescope IDs
    ifstream is;
    is.open( iInputFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, file not found: " << endl;
        cout << iInputFile << endl;
        return false;
    }
    vector< string > fComments;
    cout << "Reading " << iInputFile << endl;
    
    vector< int > fProd3bID;
    vector< string > fFOV;
    vector< string > fDyn;
    vector< string > fRaw;
    vector< string > fTelName;
    
    int iTelID = 0;
    string iTemp = "";
    string iTelName = "";
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            // comment line
            if( is_line.substr( 0, 1 ) == "#" )
            {
                fComments.push_back( is_line );
                continue;
            }
            istringstream is_stream( is_line );
            is_stream >> iTelID;
            if( iTelID == 0 )
            {
                continue;
            }
            fProd3bID.push_back( iTelID );
            // FOV
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
                fFOV.push_back( iTemp );
            }
            else
            {
                fFOV.push_back( "20." );
            }
            // dynamical range
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
                fDyn.push_back( iTemp );
            }
            else
            {
                fDyn.push_back( "12" );
            }
            // raw or calibrated
            if( !is_stream.eof() )
            {
                is_stream  >> iTemp;
                fRaw.push_back( iTemp );
            }
            else
            {
                fRaw.push_back( "1" );
            }
            // telescope name
            if( !is_stream.eof() )
            {
                is_stream  >> iTelName;
                fTelName.push_back( iTelName );
            }
            else
            {
                fTelName.push_back( getTelescopeName( iTelID ) );
            }
            
        }
    }
    is.close();
    
    // write updated array list
    ofstream is_out;
    is_out.open( iOutputFile.c_str(), ifstream::out );
    if( !is_out )
    {
        cout << "error, file not found: " << endl;
        cout << iOutputFile << endl;
        return false;
    }
    cout << "Writing to " << iOutputFile << endl;
    
    for( unsigned int i = 0; i < fComments.size(); i++ )
    {
        is_out << fComments[i] << endl;
    }
    
    for( unsigned int i = 0; i < fProd3bID.size(); i++ )
    {
        is_out << fProd3bID[i] << "   ";
        is_out << fFOV[i] << "   ";
        is_out << fDyn[i] << "   ";
        is_out << fRaw[i] << "   ";
        is_out << fTelName[i] << "   ";
        is_out << endl;
    }
    
    return;
}

/*
 * process a list of detector geometry files
 *
 * files fill be written to local directory
 *
 *
 */
void processListofFiles( string iList, string iOrgDirectory )
{
    ifstream is;
    is.open( iList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, file not found: " << endl;
        cout << iList << endl;
        return false;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            addTelescopeNamesToDetectorGeometry( iOrgDirectory + "/" + is_line, is_line );
        }
    }
}
