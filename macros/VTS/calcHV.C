/* Macro for flatfielding. Reads in gains from text file ($VERITAS_EVNDISP_AUX_DIR/Calibration/Tel_<tel>/<run>.gain)
and current HV files, outputs new HV file. Usage:

root
.L $EVNDISPSYS/macros/VTS/calcHV.C
flatfield(1,63336,"T1_default.hv","T1_flatfielded.hv")


...etc: flatfield(tel, run, input HV file, output HV file)

More info at https://veritas.sao.arizona.edu/wiki/index.php/Flat_Fielding#At_home
*/

void readingain( char* inname, int* ch, double* gain, double* gainvar )
{
    ifstream is( inname );
    if( !is )
    {
        cout << "unable to read from " << inname << endl;
    }
    for( int i = 0; i < 499; i++ )
    {
        is >> ch[i] >> gain[i] >> gainvar[i];
    }
    is.close();
    
}

void  writeoutgain( char* outname, int* ch, double* gain, double* gainvar )
{
    ofstream os( outname );
    if( !os )
    {
        cout << "unable to write to " << outname << endl;
    }
    for( int i = 0; i < 499; i++ )
    {
        os << ch[i] << " " << gain[i] << " " << gainvar[i] << endl;
    }
    os.close();
}

void readinHV( char* inname, int* ch, double* HV )
{
    ifstream is( inname );
    if( !is )
    {
        cout << "unable to read from " << inname << endl;
    }
    for( int i = 0; i < 499; i++ )
    {
        is >> ch[i] >> HV[i];
    }
    is.close();
    
}

void readinsinglepe( char* inname, int* ch, double* singlepe )
{
    ifstream is( inname );
    if( !is )
    {
        cout << "unable to read from " << inname << endl;
    }
    double dum[499];
    for( int i = 0; i < 499; i++ )
    {
        is >> ch[i] >> singlepe[i] >> dum[i] >> dum[i] >> dum[i];
    }
    is.close();
    
}

void  writeoutHV( char* outname, int* ch, double* HV )
{
    ofstream os( outname );
    if( !os )
    {
        cout << "unable to write to " << outname << endl;
    }
    for( int i = 0; i < 499; i++ )
    {
        os << ch[i] << " " << ( int )HV[i] << endl;
    }
    os.close();
}

double calcHV( double gain, double HVold )
{
    double HVnew = log10( HVold ) + log10( 1. / gain ) / 5.85;
    HVnew = pow( 10., HVnew );
    return HVnew;
    
}

double flatfield( int tel, int run, char* i_HVinfilename, char* i_HVoutfilename )
{
    char gainfilename[300];
    char HVinfilename[300];
    char HVoutfilename[300];
    sprintf( gainfilename, "%s/Calibration/Tel_%d/%d.gain", gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" ), tel, run );
    sprintf( HVinfilename, "%s", gSystem->ExpandPathName( i_HVinfilename ) );
    sprintf( HVoutfilename, "%s", gSystem->ExpandPathName( i_HVoutfilename ) );
    
    
    cout << "reading gains from: " << gainfilename << endl;
    cout << "reading input HV from: " << HVinfilename << endl;
    cout << "writing outputHV to: " << HVoutfilename << endl;
    
    int gain_ch[499];
    double gain[499];
    double gainvar[499];
    readingain( gainfilename, gain_ch, gain, gainvar );
    
    int HV_ch[499];
    double HV_in[499];
    readinHV( HVinfilename, HV_ch, HV_in );
    
    double HV_out[499];
    for( int i = 0; i < 499; i++ )
    {
        if( gain[i] > 0.1 )
        {
            HV_out[i] = calcHV( gain[i], HV_in[i] );
        }
        else
        {
            HV_out[i] = HV_in[i];
        }
    }
    writeoutHV( HVoutfilename, HV_ch, HV_out );
    
}

double flatfieldlowHV( int tel, double hvcoeff, char* i_gainfilename, char* i_HVinfilename, char* i_HVoutfilename )
{
    char gainfilename[100];
    char HVinfilename[100];
    char HVoutfilename[100];
    
    sprintf( gainfilename, "Calibration/Tel_%d/%s", tel, i_gainfilename );
    sprintf( HVinfilename, "Calibration/Tel_%d/%s", tel, i_HVinfilename );
    sprintf( HVoutfilename, "Calibration/Tel_%d/%s", tel, i_HVoutfilename );
    cout << "reading gains from: " << gainfilename << endl;
    cout << "reading input HV from: " << HVinfilename << endl;
    cout << "writing outputHV to: " << HVoutfilename << endl;
    
    int gain_ch[499];
    double gain[499];
    double gainvar[499];
    readingain( gainfilename, gain_ch, gain, gainvar );
    
    int HV_ch[499];
    double HV_in[499];
    readinHV( HVinfilename, HV_ch, HV_in );
    
    double HV_out[499];
    for( int i = 0; i < 499; i++ )
    {
        if( gain[i] > 0.1 )
        {
            HV_out[i] = calcHV( gain[i], hvcoeff * HV_in[i] );
        }
        else
        {
            HV_out[i] = hvcoeff * HV_in[i];
        }
    }
    writeoutHV( HVoutfilename, HV_ch, HV_out );
    
}

double flatfieldnumpe( int tel, char* i_gainfilename, char* i_HVinfilename, char* i_HVoutfilename )
{
    char gainfilename[100];
    char HVinfilename[100];
    char HVoutfilename[100];
    
    sprintf( gainfilename, "Tel_%d/%s", tel, i_gainfilename );
    sprintf( HVinfilename, "Tel_%d/%s", tel, i_HVinfilename );
    sprintf( HVoutfilename, "Tel_%d/%s", tel, i_HVoutfilename );
    cout << "reading gains from: " << gainfilename << endl;
    cout << "reading input HV from: " << HVinfilename << endl;
    cout << "writing outputHV to: " << HVoutfilename << endl;
    
    int gain_ch[499];
    double gain[499];
    double gainvar[499];
    readingain( gainfilename, gain_ch, gain, gainvar );
    
    double relative_numpe[499];
    double count = 0;
    double average = 0;
    for( int i = 0; i < 499; i++ )
    {
        if( gain[i] > 0.1 )
        {
            relative_numpe[i] = ( gain[i] * gain[i] ) / ( gainvar[i] * gainvar[i] );
            average += relative_numpe[i];
            count += 1.;
        }
    }
    average /= count;
    for( int i = 0; i < 499; i++ )
    {
        relative_numpe[i] /= average;
    }
    cout << "average: " << average << endl;
    
    int HV_ch[499];
    double HV_in[499];
    readinHV( HVinfilename, HV_ch, HV_in );
    
    double HV_out[499];
    double inv_relative_numpe[499];
    for( int i = 0; i < 499; i++ )
    {
        if( relative_numpe[i] > 0.1 )
        {
            inv_relative_numpe[i] = 1. / relative_numpe[i];
            HV_out[i] = calcHV( inv_relative_numpe[i], HV_in[i] );
        }
        else
        {
            HV_out[i] = HV_in[i];
            inv_relative_numpe[i] = 0;
        }
        //    cout << HV_ch[i] << "\t" << inv_relative_numpe[i] << "\t" << HV_in[i] << "\t" << HV_out[i] << endl;
    }
    writeoutHV( HVoutfilename, HV_ch, HV_out );
    
    TH1F* hv_oldhist = new TH1F( "hv_old", "hv_old", 50, 500, 1300 );
    TH1F* hv_newhist = new TH1F( "hv_new", "hv_new", 50, 500, 1300 );
    for( int i = 0; i < 499; i++ )
    {
        hv_oldhist->Fill( HV_in[i] );
        hv_newhist->Fill( HV_out[i] );
    }
    hv_oldhist->SetLineWidth( 3 );
    hv_oldhist->Draw();
    
    hv_newhist->SetLineWidth( 3 );
    hv_newhist->SetLineColor( 2 );
    
    hv_newhist->Draw( "same" );
    
}

double testnumpe( int tel, char* i_gainfilename )
{
    char gainfilename[100];
    
    sprintf( gainfilename, "Tel_%d/%s", tel, i_gainfilename );
    cout << "reading gains from: " << gainfilename << endl;
    int gain_ch[499];
    double gain[499];
    double gainvar[499];
    readingain( gainfilename, gain_ch, gain, gainvar );
    
    double relative_numpe[499];
    double count = 0;
    double average = 0;
    for( int i = 0; i < 499; i++ )
    {
        if( gain[i] > 0.1 && gainvar[i] > 0.05 )
        {
            relative_numpe[i] = ( gain[i] * gain[i] ) / ( gainvar[i] * gainvar[i] );
            average += relative_numpe[i];
            count += 1.;
        }
    }
    average /= count;
    for( int i = 0; i < 499; i++ )
    {
        relative_numpe[i] /= average;
    }
    cout << "average: " << average << endl;
    
    double inv_relative_numpe[499];
    for( int i = 0; i < 499; i++ )
    {
        if( relative_numpe[i] > 0.1 )
        {
            inv_relative_numpe[i] = 1. / relative_numpe[i];
        }
        else
        {
            inv_relative_numpe[i] = 0;
        }
        //    cout << HV_ch[i] << "\t" << inv_relative_numpe[i] << "\t" << HV_in[i] << "\t" << HV_out[i] << endl;
    }
    
    TH1F* h_numpe = new TH1F( "h_numpe", "h_numpe", 50, 0, 2 );
    for( int i = 0; i < 499; i++ )
    {
        h_numpe->Fill( relative_numpe[i] );
    }
    h_numpe->Draw();
    
}
