#include "VFITSIRFs.h"

VFITSIRFs::VFITSIRFs()
{
    fptr = 0;
}


/*
 * fits error reporting
*/
bool VFITSIRFs::printerror( int status )
{
    if( status )
    {
        fits_report_error( stderr, status );
    }
    return false;
}


bool VFITSIRFs::open_fits_file( string fits_file_name )
{
    int status = 0;
    if( fits_create_file( &fptr, fits_file_name.c_str(), &status ) )
    {
        return printerror( status );
    }
    if( fits_create_img( fptr,  DOUBLE_IMG, 0, 0, &status ) )
    {
        return printerror( status );
    }

    return true;
}

bool VFITSIRFs::write_fits_file()
{
    int status = 0;
    if( fits_close_file( fptr, &status ) )
    {
         return printerror( status );
    }
    return true;
}

bool VFITSIRFs::write_fits_header()
{
   int status = 0;

   char author[] = "Eventdisplay";
   if( fits_update_key( fptr, TSTRING, ( char* )"CREATOR", author, ( char* )"Creator", &status ) )
   {
       return printerror( status );
   }
   char telescope[] = "CTA";
   if( fits_update_key( fptr, TSTRING, ( char* )"TELESCOP", telescope, ( char* )"Telescope name", &status ) )
   {
       return printerror( status );
   }

   // GDAF
   char hdudoc[] = "https://github.com/open-gamma-ray-astro/gamma-astro-data-formats";
   if( fits_update_key( fptr, TSTRING, ( char* )"HDUDOC", hdudoc, ( char* )"Hdudoc", &status ) )
   {
       return printerror( status );
   }
   char hduver[] = "0.2";
   if( fits_update_key( fptr, TSTRING, ( char* )"HDUVERS", hduver, ( char* )"Hduvers", &status ) )
   {
       return printerror( status );
   }


   return true;
}

bool VFITSIRFs::write_psf_gauss( TH2F *h )
{
   if( !h ) return false;


   return true;
}

bool VFITSIRFs::write_histo2D( TH2F *h,
                               string name,
                               char* col_name,
                               char* col_unit )
{
   if( !h ) return false;

   int status = 0;
   const int nCol = 5;
   long nRows = h->GetNbinsX() * h->GetNbinsY();
   nRows = 0;
   char* tType[nCol] = { "ENERG_LO",
                      "ENERG_HI",
                      "THETA_LO",
                      "THETA_HI",
                       col_name };
   char* tUnit[nCol] = { "TeV",
                      "TeV",
                      "deg",
                      "deg",
                      col_unit };
   char x_form[10];
   sprintf( x_form, "%dE", h->GetNbinsX() );
   char y_form[10];
   sprintf( y_form, "%dE", h->GetNbinsY() );
   char z_form[10];
   sprintf( z_form, "%dE", h->GetNbinsX()*h->GetNbinsY() );
   char* tForm[nCol] = { &x_form[0],
                      &x_form[0],
                      &y_form[0],
                      &y_form[0],
                      &z_form[0] };

   ///////////////
   // create empty table
   if( fits_create_tbl( fptr, 
                        BINARY_TBL, 
                        nRows, 
                        nCol, 
                        tType, 
                        tForm, 
                        tUnit, 
                        name.c_str() , 
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { h->GetNbinsX(),  h->GetNbinsY() };
   if( fits_write_tdim( fptr,
                        5,
                        2,
                        naxes,
                        &status ) )
  {
      return printerror( status );
  }
   ///////////////
   // write data 
   vector< vector< float > > table;
   // xaxis
   vector< float > xedge_low;
   vector< float > xedge_hig;
   for( int i = 0; i < h->GetNbinsX(); i++ )
   {
       xedge_low.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinLowEdge( i+1 ) ) );
       xedge_hig.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinUpEdge( i+1 ) ) );
   }
   table.push_back( xedge_low );
   table.push_back( xedge_hig);
   // yaxis
   vector< float > yedge_low;
   vector< float > yedge_hig;
   for( int i = 0; i < h->GetNbinsY(); i++ )
   {
       yedge_low.push_back( h->GetYaxis()->GetBinLowEdge( i+1 ) );
       yedge_hig.push_back( h->GetYaxis()->GetBinUpEdge( i+1 ) );
   }
   table.push_back( yedge_low );
   table.push_back( yedge_hig);
   // data
   vector< float > data;
   for( int i = 0; i < h->GetNbinsX(); i++ )
   {
      for( int j = 0; j < h->GetNbinsY(); j++ )
      {
          data.push_back( h->GetBinContent( i, j ) );
      }
   }
   table.push_back( data );

   return write_table( table );
}


bool VFITSIRFs::write_table( vector< vector< float > > table )
{
   int status = 0;

   for( unsigned int i = 0; i < table.size(); i++ )
   {
       cout << "writing column " << i+1 << "\t" << table[i].size() << endl;
       if( fits_write_col( fptr, 
                       TFLOAT,
                       i+1,
                       1,
                       1,
                       table[i].size(),
                       &table[i][0],
                       &status ) )
       {
           return printerror( status );
       }
   }

   return true;
}
