/*================================================================*\
|      _    _                                        _    _        |
|     (o)--(o)                                      (o)--(o)       |
|    /.______.\    FFFF RRRR   OOO   GGG   SSS     /.______.\      |
|    \________/    F    R   R O   O G     S        \________/      |
|  ./        \.    FFF  RRRR  O   O G  GG  SSS     ./        \.    |
| ( .        , )   F    R R   O   O G   G     S   ( .        , )   |
|  \ \_\\//_/ /    F    R  RR  OOO   GGG  SSSS     \ \_\\//_/ /    |
|   ~~  ~~  ~~                                      ~~  ~~  ~~     |
| svincent@physics.utah.edu               lebohec@physics.utah.edu |
|                  VERSION 1.02 OCTOBER 10th 2011                  |
|  For license issues, see www.physics.utah.edu/gammaray/FROGS     |
\*================================================================*/

#include "frogs.h"

gsl_rng* frogs_gsl_rng;

#define FROGS_TEST 0
//================================================================
//================================================================
int frogs_print_raw_event( struct frogs_imgtmplt_in d )
{
    /*This function prints the data stored for the current even in the
      frogs_imgtmplt_in structure. It is useful when developing a
      frogs_convert_from_XXXX function*/
#define OUTUNIT stdout
    fprintf( OUTUNIT, "-------------------------------------------------\n" );
    fprintf( OUTUNIT, "EVENT %d\n", d.event_id );
    fprintf( OUTUNIT, "Elevation: %f degrees\n", d.elevation );
    fprintf( OUTUNIT, "Number of telescopes: %d\n", d.ntel );
    fprintf( OUTUNIT, "Total number of live pixels: %d\n", d.nb_live_pix_total );
    for( int tel = 0; tel < d.ntel; tel++ )
    {
        fprintf( OUTUNIT, "****************************\n" );
        fprintf( OUTUNIT, "*****   Telescope %d   *****\n", tel + 1 );
        fprintf( OUTUNIT, "****************************\n" );
        fprintf( OUTUNIT, "Position X=%fm Y=%fm\n",
                 d.scope[tel].xfield, d.scope[tel].yfield );
        fprintf( OUTUNIT, "Number of pixels:%d\n", d.scope[tel].npix );
        fprintf( OUTUNIT, "Number of live pixels:%d\n", d.scope[tel].nb_live_pix );
    }
    fprintf( OUTUNIT, ".................................................\n" );
    fprintf( OUTUNIT, "            STARTING POINT:\n" );
    fprintf( OUTUNIT, "  xs=%f deg. ys=%f deg. \n", d.startpt.xs, d.startpt.ys );
    fprintf( OUTUNIT, "  xp=%f m. yp=%f m. \n", d.startpt.xp, d.startpt.yp );
    fprintf( OUTUNIT, "  Log10(E/1TeV)=%f \n", d.startpt.log10e );
    fprintf( OUTUNIT, "  Lambda=%f Xo\n", d.startpt.lambda );
    
    fprintf( OUTUNIT, "-------------------------------------------------\n" );
    if( d.worthy_event == FROGS_OK )
    {
        fprintf( OUTUNIT, "  GOOD EVENT\n" );
    }
    else
    {
        fprintf( OUTUNIT, "  BAD EVENT\n" );
    }
    
    fprintf( OUTUNIT, "-------------------------------------------------\n" );
    return FROGS_OK;
}
//================================================================
//================================================================
// MAIN FUNCTION
/*  This function performs the image template analysis. It returns a
*   structure of type frogs_imgtmplt_out containing all the useful
*   information  regarding the convergence of the likelihood
*   optimization. All the data from the telescopes for the event
*   being  analyzed are in the structure frogs_imgtmplt_in d
*   received as an argument.
*/
struct frogs_imgtmplt_out frogs_img_tmplt( struct frogs_imgtmplt_in* d, char templatelistname[FROGS_FILE_NAME_MAX_LENGTH] )
{
    
    struct frogs_imgtmplt_out rtn;
    
    /* Event selection used to avoid processing bad or poor events
       worthy_event is calculated in frogs_convert_from_grisu (or
       equivalent for other analysis packages)*/
    if( d->worthy_event == FROGS_NOTOK )
    {
        rtn = frogs_null_imgtmplt_out();
        rtn.event_id = d->event_id;
        rtn.nb_iter = FROGS_BAD_NUMBER;
        rtn.gsl_convergence_status = FROGS_BAD_NUMBER;
        //Release memory used in the data structure
        frogs_release_memory( d );
        return rtn;
    }
    
    static struct calibration_file calib;
    static struct frogs_probability_array prob_array;
    
    /*This check if a template image needs to be read and reads it if necessary*/
    static struct frogs_imgtemplate tmplt;
    static int firstcall = 1;
    //On the first call set the template elevation to zero
    if( firstcall )
    {
        tmplt.elevmin = 0;
        tmplt.elevmax = 0;
        firstcall = 0;
        frogs_fill_prob_density( &prob_array );
    }
    
    //If needed read the template file according to elevation
    if( d->elevation > tmplt.elevmax || d->elevation < tmplt.elevmin )
    {
        tmplt = frogs_read_template_elev( d->elevation, templatelistname );
    }
    
    //Optimize the likelihood
    rtn = frogs_likelihood_optimization( d, &tmplt, &calib, &prob_array );
    
    //Release memory used in the data structure
    frogs_release_memory( d );
    
    return rtn;
}

//================================================================
//================================================================
int frogs_release_memory( struct frogs_imgtmplt_in* d )
{
    /*This function is used to release the memory taken up by each event. */
    for( int tel = 0; tel < d->ntel; tel++ )
    {
        delete [] d->scope[tel].xcam;
        delete [] d->scope[tel].ycam;
        delete [] d->scope[tel].q;
        delete [] d->scope[tel].ped;
        delete [] d->scope[tel].exnoise;
        delete [] d->scope[tel].pixinuse;
        delete [] d->scope[tel].telpixarea;
        delete [] d->scope[tel].pixradius;
    }
    delete [] d->scope;
    return 0;
}
//================================================================
//================================================================
struct frogs_imgtmplt_out frogs_null_imgtmplt_out()
{
    //Used to set the image template analysis output to some null values
    //when the analysis failed.
    struct frogs_imgtmplt_out rtn;
    rtn.nb_iter						  = FROGS_BAD_NUMBER;
    rtn.goodness_img				  = FROGS_BAD_NUMBER;
    rtn.npix_img					  = FROGS_BAD_NUMBER;
    rtn.goodness_bkg				  = FROGS_BAD_NUMBER;
    rtn.npix_bkg					  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.xs					  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.ys					  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.xp					  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.yp					  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.log10e				  = FROGS_BAD_NUMBER;
    rtn.cvrgpt.lambda				  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.xs				  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.ys				  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.xp				  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.yp				  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.log10e			  = FROGS_BAD_NUMBER;
    rtn.cvrgpterr.lambda			  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessImg[0]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessImg[1]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessImg[2]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessImg[3]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessBkg[0]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessBkg[1]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessBkg[2]		  = FROGS_BAD_NUMBER;
    rtn.tel_goodnessBkg[3]		  = FROGS_BAD_NUMBER;
    return rtn;
}
//================================================================
//================================================================
/*  This function optimizes the likelihood function on the event parameter
*   space. It returns all the relevant information regarding the convergence
*   in a structure frogs_imgtmplt_out. d is a pointer to a structure
*   containing all the data from the telescopes. tmplt is a pointer to a
*   structure containing the image template data
*/
struct frogs_imgtmplt_out
frogs_likelihood_optimization( struct frogs_imgtmplt_in* d,
                               struct frogs_imgtemplate* tmplt,
                               struct calibration_file* calib,
                               struct frogs_probability_array* prob_array )
{

    struct frogs_imgtmplt_out rtn;
    
#ifdef DIFF_EVOLUTION
    //Global optimization by differential evolution
    frogs_differential_evolution( d, tmplt, prob_array );
#endif //DIFF_EVOLUTION
    
    //Gets the starting point in a GSL vector.
    const size_t p = tmplt->ndim + 1;  //Dimension of the parameter space
    //p=6 {xs, ys, xp, yp, logE, lambda}
    gsl_vector* x;
    x = gsl_vector_alloc( p );
    gsl_vector_set( x, FROGS_XS, d->startpt.xs );
    gsl_vector_set( x, FROGS_YS, d->startpt.ys );
    gsl_vector_set( x, FROGS_XP, d->startpt.xp );
    gsl_vector_set( x, FROGS_YP, d->startpt.yp );
    gsl_vector_set( x, FROGS_LOG10E, d->startpt.log10e );
    gsl_vector_set( x, FROGS_LAMBDA, d->startpt.lambda );
    
    //Prepare function minimization
    const size_t n = d->nb_live_pix_total;  //Number of pixels used
    gsl_multifit_function_fdf func;
    func.n	= n;					  //number of pixels used
    func.p	= p;					  //Dimension of the parameter space
    func.f	= &frogs_likelihood;			 //frogs_likelihood: function to be minimized
    func.df	= &frogs_likelihood_derivative;		//frogs_likelihood_derivative
    func.fdf	= &frogs_likelihood_fdf;	 //Weird function GSL wants to see
    //Construct a structure with all the data to be passed to the image model
    struct  frogs_gsl_data_wrapper data;
    data.tmplt		= tmplt;		  //template data
    data.data		= d;			  //telescope data
    data.probarray	= prob_array;
    func.params		= ( void* )&data; /*This will be passed to the functions
														  frogs_likelihood and
														  frogs_likelihood_derivative
														  as a pointer to a void*/
    const gsl_multifit_fdfsolver_type* T;
    gsl_multifit_fdfsolver* s;
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc( T, n, p );
    gsl_multifit_fdfsolver_set( s, &func, x );
    
    //Loop on iterations
    int status = 0;
    unsigned int iter = 0;
    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate( s );
        if( status )
        {
            break;
        }
        status = gsl_multifit_test_delta( s->dx, s->x, 1e-7, 1e-7 );
    }
    while( status == GSL_CONTINUE && iter < FROGS_MAX_ITER_NBR );
    
    //If the maximum number of iteration was reached, we return a null result.
    if( iter == FROGS_MAX_ITER_NBR )
    {
        rtn										 = frogs_null_imgtmplt_out();
        rtn.gsl_convergence_status			 = status;
        rtn.event_id							 = d->event_id;
        rtn.nb_iter								 = iter;
        gsl_multifit_fdfsolver_free( s );
        return rtn;
    }
    rtn.gsl_convergence_status = status;
    rtn.event_id = d->event_id;
    rtn.nb_iter = iter;
    //Gathers the parameter space convergence point coordinates.
    rtn.cvrgpt.xs = gsl_vector_get( s->x, FROGS_XS );
    rtn.cvrgpt.ys = gsl_vector_get( s->x, FROGS_YS );
    rtn.cvrgpt.xp = gsl_vector_get( s->x, FROGS_XP );
    rtn.cvrgpt.yp = gsl_vector_get( s->x, FROGS_YP );
    rtn.cvrgpt.log10e = gsl_vector_get( s->x, FROGS_LOG10E );
    rtn.cvrgpt.lambda = gsl_vector_get( s->x, FROGS_LAMBDA );
    
    /*GSL take lambda to any value but only the value used in the model
      calculation is directly meaningfull. Here, we convert the value set
      by GSL to the actual value used in the model calculation. */
    float maxlambda = tmplt->min[0] + ( tmplt->nstep[0] - 1 ) * tmplt->step[0];
    rtn.cvrgpt.lambda = floatwrap( rtn.cvrgpt.lambda, tmplt->min[0], maxlambda );
    
    gsl_matrix* covar = gsl_matrix_alloc( p, p );
    // necessary updates for gsl version 2.0 and newer
    // follow entry in https://sft.its.cern.ch/jira/browse/ROOT-7776
#ifdef GSL2
    gsl_matrix* J = gsl_matrix_alloc( s->fdf->n, s->fdf->p );
    gsl_multifit_fdfsolver_jac( s, J );
    gsl_multifit_covar( J, 0.0, covar );
# else
    gsl_multifit_covar( s->J, 0.0, covar );
#endif
    
    
    //Get the error for each parameter
    float chi = gsl_blas_dnrm2( s->f );
    float dof = n - p;
    float c = GSL_MAX_DBL( 1, chi / sqrt( dof ) );
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
    //ERR(i) is nan when iter==1. We then change the error to FROGS_BAD_NUMBER
    rtn.cvrgpterr.xs = c * ERR( FROGS_XS );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.xs ) )
    {
        rtn.cvrgpterr.xs = FROGS_BAD_NUMBER;
    }
    rtn.cvrgpterr.ys = c * ERR( FROGS_YS );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.ys ) )
    {
        rtn.cvrgpterr.ys = FROGS_BAD_NUMBER;
    }
    rtn.cvrgpterr.xp = c * ERR( FROGS_XP );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.xp ) )
    {
        rtn.cvrgpterr.xp = FROGS_BAD_NUMBER;
    }
    rtn.cvrgpterr.yp = c * ERR( FROGS_YP );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.yp ) )
    {
        rtn.cvrgpterr.yp = FROGS_BAD_NUMBER;
    }
    rtn.cvrgpterr.log10e = c * ERR( FROGS_LOG10E );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.log10e ) )
    {
        rtn.cvrgpterr.log10e = FROGS_BAD_NUMBER;
    }
    rtn.cvrgpterr.lambda = c * ERR( FROGS_LAMBDA );
    if( !frogs_is_a_good_number( rtn.cvrgpterr.lambda ) )
    {
        rtn.cvrgpterr.lambda = FROGS_BAD_NUMBER;
    }
    
    /*Calculate the image and background goodness for the convergence point.*/
    if( frogs_goodness( &rtn, d, tmplt, calib, prob_array ) != FROGS_OK )
    {
        frogs_showxerror( "Error: Problem encountered in the convergence goodness calculation" );
    }
    
    gsl_multifit_fdfsolver_free( s );
    gsl_matrix_free( covar );
    
    return rtn;
}
//================================================================
//================================================================
int frogs_goodness( struct frogs_imgtmplt_out* tmplanlz,
                    struct frogs_imgtmplt_in* d,
                    struct frogs_imgtemplate* tmplt,
                    struct calibration_file* calib,
                    struct frogs_probability_array* prob_array )
{
    /* Calculates the image and background goodness whose values and associated
       number of pixels are passed back in tmplanlz which also contains the
       values of the event physical parameters for which the goodness is
       calculated. d is a pointer to a structure holding all the data from the
       telescopes while tmplt is a pointer to a structure holding the image
       template data */
    
    int telnpix[4] = {0}; // pixel counter for single camera
    
    /*Initialize goodnesses and the number of pixels for both image
      and background regions*/
    tmplanlz->goodness_img = 0;
    tmplanlz->npix_img	  = 0;
    tmplanlz->goodness_bkg = 0;
    tmplanlz->npix_bkg	  = 0;
    
    /* Calculates the telescope image goodness and the telescope background goodness.
       The goodness-of-fit per telescope are initially set to 0. To prevent fake "good"
       goodness-of-fit, i.e. goodness with 0 value but NO image in the telescope, we later
       do a test (see fake goodness test) */
    tmplanlz->tel_goodnessImg[0] = 0.;
    tmplanlz->tel_goodnessImg[1] = 0.;
    tmplanlz->tel_goodnessImg[2] = 0.;
    tmplanlz->tel_goodnessImg[3] = 0.;
    tmplanlz->tel_goodnessBkg[0] = 0.;
    tmplanlz->tel_goodnessBkg[1] = 0.;
    tmplanlz->tel_goodnessBkg[2] = 0.;
    tmplanlz->tel_goodnessBkg[3] = 0.;
    
    for( int tel = 0; tel < d->ntel; tel++ )
    {
        telnpix[tel] = 0;
        for( int pix = 0; pix < d->scope[tel].npix; pix++ )
        {
            if( d->scope[tel].pixinuse[pix] == FROGS_OK )
            {
                //Here call image model and calculate expected signal
                int pix_in_template;//FROGS_OK in image, FROGS_NOTOK in background
                double mu = frogs_img_model( pix, tel, tmplanlz->cvrgpt, d,
                                             tmplt, &pix_in_template );
                /* eventdisplay OUTPUT to see templates in the display gui
                   should be commented if you are not an eventdisplay user*/
                tmplanlz->tmplt_tubes[tel][pix] = mu;
                
                if( mu != FROGS_BAD_NUMBER )
                {
                    double pd;
                    
                    /* to speed up the analysis, a lookup table is used to get
                       the probability densities. The lookup table is only used for
                       intermediate values, i.e. mu>0 and mu<FROGS_LARGE_PE_SIGNAL
                       (saves a factor of 2x in computing time). */
                    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
                    {
                        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu, d->scope[tel].ped[pix] );
                    }
                    else
                        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                                        d->scope[tel].ped[pix],
                                                        d->scope[tel].exnoise[pix] );
                                                        
                    double mean_lkhd = frogs_mean_pix_lkhd( d->scope[tel].q[pix], mu,
                                                            d->scope[tel].ped[pix],
                                                            d->scope[tel].exnoise[pix],
                                                            prob_array );
                                                            
                    double pix_goodness = -2.0 * log( pd ) - mean_lkhd;
                    
                    //If requested we produce a calibration output
                    //if( FROGS_NBEVENT_GDNS_CALIBR > 0 )
                    if( d->nb_events_calib > 0 )
                        frogs_gdns_calibr_out( d->nb_events_calib, d->event_id, tel, pix, d->scope[tel].q[pix], d->scope[tel].ped[pix], mu,
                                               pix_goodness, tmplanlz->cvrgpt.log10e, tmplanlz->cvrgpt.xp, tmplanlz->cvrgpt.yp, d->scope[tel].xcam[pix], d->scope[tel].ycam[pix] );
                                               
                    /*Apply the single pixel goodness correction according to the
                      pixel pedestal width and the model value mu*/
                    //frogs_goodness_correction(pix_goodness,d->scope[tel].ped[pix],mu);
                    
                    /*Decides if the pixel should be counted in the image
                      or background region*/
                    //int pix_in_img=frogs_image_or_background(tel,pix,d);
                    int pix_in_img = frogs_image_or_background( tel, pix, d, mu );
                    
                    //If requested, we produce a display of the event
                    if( FROGS_NBEVENT_DISPLAY > 0 )
                        frogs_event_display( d->event_id, d->scope[tel].q[pix], mu,
                                             d->scope[tel].xfield, d->scope[tel].yfield,
                                             d->scope[tel].xcam[pix],
                                             d->scope[tel].ycam[pix], pix_in_img );
                                             
                    //If the pixel is in the image region
                    if( pix_in_img == FROGS_OK )
                    {
                        tmplanlz->goodness_img = tmplanlz->goodness_img + pix_goodness;
                        tmplanlz->npix_img++;
                        tmplanlz->tel_goodnessImg[tel] += pix_goodness;
                        /*counts the number of pixel in the image region.
                          telnpix is used in fake goodness test*/
                        telnpix[tel]++;
                    }
                    //If the pixel is in the background region
                    if( pix_in_img == FROGS_NOTOK )
                    {
                        tmplanlz->goodness_bkg = tmplanlz->goodness_bkg + pix_goodness;
                        tmplanlz->npix_bkg++;
                        tmplanlz->tel_goodnessBkg[tel] += pix_goodness;
                    }
                }//End of background/image region test
            } //End of test on pixel viability
            
            
        }//End of pixel loop
        
        /* Finalize the telescope image and background goodness.
        The sum over all telescope goodness should equal to the whole array goodness */
        if( d->nb_live_pix_total > tmplt->ndim + 1 )
        {
            tmplanlz->tel_goodnessImg[tel] = tmplanlz->tel_goodnessImg[tel] / sqrt( 2.0 * ( d->nb_live_pix_total - ( tmplt->ndim + 1 ) ) );
        }
        else
        {
            tmplanlz->tel_goodnessImg[tel] = FROGS_BAD_NUMBER;
        }
        if( d->nb_live_pix_total > tmplt->ndim + 1 )
        {
            tmplanlz->tel_goodnessBkg[tel] = tmplanlz->tel_goodnessBkg[tel] / sqrt( 2.0 * ( d->nb_live_pix_total - ( tmplt->ndim + 1 ) ) );
        }
        else
        {
            tmplanlz->tel_goodnessBkg[tel] = FROGS_BAD_NUMBER;
        }
        
    }//End of telescope loop
    
    /* Fake goodness test. If the number of pixel in the image region for a particular telescope
     is 0, the goodness of fit is set to FROGS_BAD_NUMBER */
    for( int itel = 0; itel < d->ntel; itel++ )
    {
        if( telnpix[itel] == 0 || fabs( tmplanlz->tel_goodnessImg[itel] ) < 1E-7 )
        {
            tmplanlz->tel_goodnessImg[itel] = FROGS_BAD_NUMBER;
        }
        if( telnpix[itel] == 0 || fabs( tmplanlz->tel_goodnessBkg[itel] ) < 1E-7 )
        {
            tmplanlz->tel_goodnessBkg[itel] = FROGS_BAD_NUMBER;
        }
    }
    
    // GH Find Number of Images (Using above reversed)
    
    /*note:
    this only works if ntel < 8*sizeof(unsigned long) = 64.
    unfortunately, c does not have easy access to bit sets.
    so we'll use an unsigned long integer instead.
    1<<itel shifts left by itel bits (same as 2^itel).
    so if eg. T2 and T4 were the only active telescope, then we would have
    tmplanlz->selected_images = 1<<(2-1) + 1<<(4-1) = 1<<(1) + 1<<(3) = 2 + 8 = 10 = (1010) as expected.
    */
    
    tmplanlz->nb_images = 0;
    tmplanlz->selected_images = 0;
    
    for( int itel = 0; itel < d->ntel; itel++ )
    {
        if( telnpix[itel] != 0 && fabs( tmplanlz->tel_goodnessImg[itel] ) >= 1E-7 )
        {
            tmplanlz->nb_images++;
            tmplanlz->selected_images += ( 1 << itel );
        }
    }
    
    //Finilize the background goodness calculation (*** See note)
    if( d->nb_live_pix_total > tmplt->ndim + 1 )
        tmplanlz->goodness_bkg = tmplanlz->goodness_bkg /
                                 sqrt( 2.0 * ( d->nb_live_pix_total - ( tmplt->ndim + 1 ) ) );
    //sqrt(2.0*(tmplanlz->npix_bkg-(tmplt->ndim+1)));
    else
    {
        tmplanlz->goodness_bkg = FROGS_BAD_NUMBER;
    }
    //Finilize the image goodness calculation (*** See note)
    if( d->nb_live_pix_total > tmplt->ndim + 1 )
        tmplanlz->goodness_img = tmplanlz->goodness_img /
                                 sqrt( 2.0 * ( d->nb_live_pix_total - ( tmplt->ndim + 1 ) ) );
    //sqrt(2.0*(tmplanlz->npix_img-(tmplt->ndim+1)));
    else
    {
        tmplanlz->goodness_img = FROGS_BAD_NUMBER;
    }
    
    /* ***note on the goodness. The goodness is the difference between the
       log-likelihood and its average value divided by the square root of two
       times the number of degrees of freedom. Here we calculate a goodness for
       the image and for the background and we subtract the number of optimized
       parameter to both. This seems odd. It probably does not matter much when
       the number of pixels in the image and in the background are large
       compared to the number of optimized parameters.*/
    
    return FROGS_OK;
}
//================================================================
//================================================================
//int frogs_image_or_background(int tel,int pix,struct frogs_imgtmplt_in *d){
int frogs_image_or_background( int tel, int pix, struct frogs_imgtmplt_in* d, double mu )
{
    /*Returns FROGS_OK if the pixel is identified to be part of the image
      region and returns FROGS_NOTOK is the pixel is to be counted in the
      bacground region.
      The image is defined by the pixels surviving a two threshold standard
      cleaning algorithm. If a pixel is less than FROGS_PICTRAD from a pixel
      passing the cleaning, it is in the image region. Otherwise it is
      in the bacground region.
      If a pixel signal exceeds FROGS_HITHRESH passes the cleaning.
      If a pixel signal exceeds FROGS_LOTHRESH it passes the cleaning if
      it has a drect neighbor whose signal exceeds FROGS_HITHESH. The
      thesholds FROGS_HITHRESH and FROGS_LOTHRESH are expressed in units
      of pedestal standard deviations.
    */
    
    //If the pixel tested is not a pixel in use, returns FROGS_BAD_NUMBER
    if( d->scope[tel].pixinuse[pix] == FROGS_NOTOK )
    {
        return FROGS_BAD_NUMBER;
    }
    
    /*If the signal in the pixel is above the higher threshold,
    the pixel is part of the image*/
    if( mu > 0.001 )
    {
        return FROGS_OK;
    }
    return FROGS_NOTOK;
    
    
    //if( d->scope[tel].q[pix] > FROGS_HITHRESH * d->scope[tel].ped[pix] )
    //{
    //return FROGS_OK;
    //}
    
    //Get the square of the radius defining the picture region
    //float d2max = FROGS_PICTRAD * FROGS_PICTRAD;
    
    //Loop over the pixels
    //for( int p = 0; p < d->scope[tel].npix; p++ )
    //{
    //Only for pixels in use
    //if( d->scope[tel].pixinuse[p] == FROGS_OK )
    //{
    //Distance to the tested pixel
    //float d2 = ( d->scope[tel].xcam[pix] - d->scope[tel].xcam[p] ) *
    //( d->scope[tel].xcam[pix] - d->scope[tel].xcam[p] ) +
    //( d->scope[tel].ycam[pix] - d->scope[tel].ycam[p] ) *
    //( d->scope[tel].ycam[pix] - d->scope[tel].ycam[p] );
    //If the distance is less than the picture definition radius
    //if( d2 < d2max )
    //{
    /*If the pixel has a signal exceeding the higher threshold, the
      tested pixel is counted in the image. */
    //	 if( d->scope[tel].q[p] > FROGS_HITHRESH * d->scope[tel].ped[p] )
    //		{
    //		  return FROGS_OK;
    //		}
    /*If the pixel has a signal exceeding the lower threshold but not
      the higher, it needs to have a direct neighbor exceeding the
      higher threshold for the tested pixel to be part of the picture
      region*/
    //	 if( d->scope[tel].q[p] > FROGS_LOTHRESH * d->scope[tel].ped[p] )
    //		{
    //Square of the maximal distance for pixels to be neighbors
    //		  float dnb2max = FROGS_NEIGHBORAD * FROGS_NEIGHBORAD;
    //Loop over the pixels
    //		  for( int pnb = 0; pnb < d->scope[tel].npix; pnb++ )
    //			 {
    //Only for pixels in use
    //if( d->scope[tel].pixinuse[pnb] == FROGS_OK )
    //				  {
    //					 float dnb2 = ( d->scope[tel].xcam[p] - d->scope[tel].xcam[pnb] ) *
    //						( d->scope[tel].xcam[p] - d->scope[tel].xcam[pnb] ) +
    //						( d->scope[tel].ycam[p] - d->scope[tel].ycam[pnb] ) *
    //						( d->scope[tel].ycam[p] - d->scope[tel].ycam[pnb] );
    //if pixels p and pnb are neighbors
    //					 if( dnb2 < dnb2max )
    /*if the neighbor of p has a signal in excess of the high
      threshold the tested pixel is counted in the picture. */
    //						if( d->scope[tel].q[pnb] > FROGS_HITHRESH * d->scope[tel].ped[pnb] )
    //{
    //	return FROGS_OK;
    //						  }
    //				  } //end of the pixinuse test in the search for neighbors of pixel p
    //			 } //End of the loop searching for neighbors of pixel p
    //		} //End of the test for LOTHRESH
    //}//End of test on the distance
    //}//End of the test of pixel being in use
    //}//End of the main loop on pixels
    //return FROGS_NOTOK;
}
//================================================================
//================================================================
int frogs_gdns_calibr_out( int nb_evens_printed_out, int event_id, int tel, int pix, float q,  float ped, float mu,
                           double pix_goodness, double energy, double xp, double yp, double xcam, double ycam )
{
    /*This funtion is used to print out information used to establish the
      calibration correction to be applied to individual pixel goodness
      values. On the first time it is called, it opens a file named
      frogs_goodness_calibration.frogs and starts counting events, printing
      calibration data for each pixel. When the number of events exceeds
      FROGS_NBENT_GDNS_CALIBR it closes the file and stops the program
      This function should be called from the function frogs_goodness*/
    
    static int last_event_id = FROGS_NOTOK; //Stores the last event id.
    static int first_time = FROGS_OK;     //Check for first time call
    static int nbrevt_calib = 0;          //Counter of events used for calibration
    static FILE* calib;                   //File pointer
    //On the first call to this function
    if( first_time == FROGS_OK )
    {
        first_time = FROGS_NOTOK;
        //Open the file
        char* itemp = getenv( "VERITAS_USER_DATA_DIR" );
        char FROGS_CALIBRATION_PATH[500];
        sprintf( FROGS_CALIBRATION_PATH, "%s/frogs_calibration/frogs_goodness_calibration.frogs", itemp );
        //calib = fopen( "%s/frogs_calibration/frogs_goodness_calibration.frogs", "w" );
        calib = fopen( FROGS_CALIBRATION_PATH, "w" );
        if( calib == NULL )
        {
            frogs_showxerror( "Error: Failed opening the file frogs_goodness_calibration.frogs for writing" );
        }
    }
    
    //If the event id is different from the one at the last call of this function
    if( event_id != last_event_id )
    {
        //Update the last event id for next time
        last_event_id = event_id;
        //Count one more event used for the calibration
        nbrevt_calib++;
    }
    
    //If we collected enough events we ned to stop
    //if( nbrevt_calib > FROGS_NBEVENT_GDNS_CALIBR )
    if( nbrevt_calib > nb_evens_printed_out )
    {
        //Close the file
        fclose( calib );
        //Stop the execussion
        frogs_showxerror( "Error: Done writing the calibration file data in frogs_goodness_calibration.frogs" );
        
        return FROGS_OK;
    }
    
    //If we get here it means we have an open file and we need to print the data
    fprintf( calib, "%d %d %d %f %f %f %g %f %f %f %f %f\n", event_id, tel, pix, ped, q, mu,
             pix_goodness, energy, xp, yp, xcam, ycam );
    return FROGS_OK;
}
//================================================================
//================================================================
int frogs_event_display( int event_id, float q, float mu, float xtel,
                         float ytel, float xpix, float ypix, int pix_in_img )
{
    /* This funtion is used to print out a kumac script to be executed from
       a PAWX11 session to obtain a display of the events and their template
       fit. On the first time it is called, it opens a file named
       frogs_display.kumac and starts counting events, printing drawing
       script commands for each pixel. When the number of events exceeds
       FROGS_NBEVENT_DISPLAY it closes the file and stops the program
       This function should be called from the function frogs_goodness*/
    
    static int last_event_id = FROGS_NOTOK; //Stores the last event id.
    static int first_time = FROGS_OK; //Check for first time call
    static int nbrevt_display = 0;     //Counter of events used for calibration
    static FILE* display;
    
    //On the first call to this function
    if( first_time == FROGS_OK )
    {
        first_time = FROGS_NOTOK;
        //Open the file
        display = fopen( "frogs_display.kumac", "w" );
        if( display == NULL )
        {
            frogs_showxerror( "Error: Failed opening the file frogs_goodness_calibration.frogs for writing" );
        }
    }
    
    //If we collected enough events we need to stop
    if( nbrevt_display > FROGS_NBEVENT_DISPLAY )
    {
        //Close the file
        fclose( display );
        //Stop the execussion
        frogs_showxerror( "Error: Done writing the event display file  frogs_display.kumac" );
        return FROGS_OK;
    }
    
    //If the event id is different from the one at the last call of this function
    if( event_id != last_event_id )
    {
        if( last_event_id == FROGS_NOTOK )
        {
            fprintf( display, "opt ntic\n" );
        }
        else
        {
            fprintf( display, "wait\n" );
        }
        //Update the last event id for next time
        last_event_id = event_id;
        //Count one more event used for the display
        nbrevt_display++;
        fprintf( display, "HISTOGRAM/CREATE/TITLE_GLOBAL 'Event %d'\n", event_id );
        fprintf( display, "nul -150 150 -150 150 ab\n" );
    }
    
    float scale = 15.0; //Angular scale meters/degrees
    float thresh = 6; //Smallest number of photons to be displayed
    float maxlight = 50; //Maximal signal
    float maxrad = 0.075; //Maximal pixel radius in degrees
    //If we get here it means we have an open file and we need to draw a pixel
    if( q > thresh || mu > thresh )
    {
        float radius;
        radius = ( q - thresh ) / maxlight;
        if( radius > 1 )
        {
            radius = 1;
        }
        if( radius > 0 )
        {
            radius = maxrad * sqrt( radius );
            fprintf( display, "set plci 2\n" );
            fprintf( display, "arc %f %f %f\n", xtel + xpix * scale, ytel + ypix * scale,
                     radius * scale );
        }
        radius = ( mu - thresh ) / maxlight;
        if( radius > 1 )
        {
            radius = 1;
        }
        if( radius > 0 )
        {
            radius = maxrad * sqrt( radius );
            fprintf( display, "set plci 3\n" );
            fprintf( display, "arc %f %f %f\n", xtel + xpix * scale, ytel + ypix * scale,
                     radius * scale );
        }
    }
    if( pix_in_img == FROGS_OK )
    {
        fprintf( display, "set plci 1\n" );
        fprintf( display, "arc %f %f %f\n", xtel + xpix * scale, ytel + ypix * scale,
                 maxrad * scale );
    }
    else
    {
        fprintf( display, "set plci 1\n" );
        fprintf( display, "arc %f %f %f\n", xtel + xpix * scale, ytel + ypix * scale,
                 maxrad * scale * 0.1 );
                 
    }
    return FROGS_OK;
}
//================================================================
//================================================================
int frogs_is_a_good_number( double x )
{
    //Returns 0 is x is a NaN of a Inf and returns 1 otherwise.
    if( isnan( x ) )
    {
        return 0;
    }
    if( isinf( x ) )
    {
        return 0;
    }
    return 1;
}
//================================================================
//================================================================
float frogs_goodness_correction( float goodness0, float ped, float mu )
{
    /* Applies correction to the individual pixel goodness in order to
       compensate for its sensitivity to the pedestal width. The function
       returns the single pixel corrected goodness value
       goodness0 = uncorrected goodness
       ped = pedestal width */
    
    return goodness0;
    
}
//================================================================
//================================================================
double frogs_probability_density( float q, double mu, float ped,
                                  float exnoise )
{
    double rtn = 0.0;
    /* This function returns the probability density of obtaining a signal q
       when the average of the expected signal is mu
       q = actual signal in the pixel
       mu = expectation value of the signal in the pixel
       ped = pedestal width for that pixel
       exnoise = exess noise for that pixel */
    
    //Case mu=0
    if( mu < 1.0e-18 )
    {
        rtn = exp( -q * q / ( 2.0 * ped * ped ) ) / ( FROGS_SQRTTWOPI * ped );
    }
    
    //Case where mu is large enough the poisson distribution is gauss-like
    if( mu >= FROGS_LARGE_PE_SIGNAL )
    {
        double dummy1 = ped * ped + mu * ( 1.0 + exnoise * exnoise ); //Variance
        double dummy2 = exp( -( q - mu ) * ( q - mu ) * 0.5 / dummy1 ); //Gauss probability density
        rtn = dummy2 / sqrt( FROGS_TWOPI * dummy1 );   //Normalization
    }
    
    //Detailed calculation for intermediate values of mu
    if( mu >= 1.0e-18 && mu < FROGS_LARGE_PE_SIGNAL )
    {
        float stdev = sqrt( mu * ( 1 + exnoise * exnoise ) + ped * ped );
        int Nmin = ( int )( floor( mu - FROGS_NUMBER_OF_SIGMA * stdev ) );
        if( Nmin < 0 )
        {
            Nmin = 0;
        }
        int Nmax = ( int )( floor( mu + FROGS_NUMBER_OF_SIGMA * stdev ) );
        if( Nmax < Nmin )
        {
            Nmax = Nmin;    //Although it does not seem possible
        }
        //Sum over the possible number of photons contributing to the signal
        for( int n = Nmin; n <= Nmax; n++ )
        {
            //Poisson factor
            double dummy1 = frogs_poisson_distribution( mu, n );
            double dummy2 = ped * ped + n * exnoise * exnoise; //standard deviation
            dummy1 = dummy1 / sqrt( FROGS_TWOPI * dummy2 ); //Gauss distrib. normaliz. factor
            rtn = rtn + dummy1 * exp( -( q - n ) * ( q - n ) / ( 2.0 * dummy2 ) );/*Gauss distribution with
						      Poisson probability
						      weight*/
        }
    }
    
    /* If the probability density numerically drops to zero, the penalty
       payed for that pixel diverges so we cap it so as to retain sensitivity
       to parameter changes from other pixels. */
    if( rtn < FROGS_SMALL_NUMBER )
    {
        rtn = FROGS_SMALL_NUMBER;
    }
    return rtn;
}
//================================================================
//================================================================
double frogs_poisson_distribution( double mu, long int n )
{
    //Compute the Poisson probability to get n for an average mu
    double rtn;
    if( mu > 0.0 )
    {
        /*P(mu,n)=mu^n*exp(-mu)/n! but we go logarithmic to avoid
          numerical problems.*/
        rtn = n * log( mu ) - mu - frogs_logarithm_factorial( n );
        rtn = exp( rtn );
    }
    else  //that is if mu==0
    {
        if( n == 0 )
        {
            rtn = 1.0;
        }
        else
        {
            rtn = 0.0;
        }
    }
    return rtn;
}
//================================================================
//================================================================
double frogs_logarithm_factorial( long int n )
{
    //Compute the log of n!
    
    //Returns 0 when n<0
    if( n < 0 )
    {
        return 0;
    }
    
    //Returns tabulated value for n<11
    if( n < 11 )
    {
        int a[11] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
        return ( double ) log( a[n] );
    }
    else
    {
        //What follows is Srinivasa Ramanujan approximation of n!
        //For n=10 it returns log(3628797.9)
        //For n=11 it returns log(39916782) instead of log(3991680)
        return n * log( n ) - n + log( n * ( 1.0 + 4.0 * n * ( 1.0 + 2.0 * n ) ) ) / 6.0 + 0.57236494;
        //0.57236494 is 0.5*log(pi)
    }
}
//================================================================
//================================================================
//double frogs_mean_pix_lkhd(float q,double mu,float ped,float exnoise) {
double frogs_mean_pix_lkhd( double q, double mu, double ped, double exnoise, struct frogs_probability_array* prob_array )
{
    double rtn = 0;
    /* Computes and returns the average of the log-likelihood for the
       considered pixel for which:
       q = actual signal in the pixel
       mu = expectation value of the signal in the pixel
       ped = pedestal width for that pixel
       exnoise = exess noise for that pixel */
    
    //Case mu=0
    if( mu < 1.0e-18 )
    {
        rtn = 1.0 + FROGS_LNTWOPI + 2.0 * log( ped );
    }
    
    //Case where mu is large enough the poisson distribution is gauss-like
    if( mu >= FROGS_LARGE_PE_SIGNAL )
    {
        double dummy = ped * ped + mu * ( 1.0 + exnoise * exnoise );
        rtn = 1.0 + FROGS_LNTWOPI + log( dummy );
    }
    
    //Detailed calculation for intermediate values of mu
    if( mu >= 1.0e-18 && mu < FROGS_LARGE_PE_SIGNAL )
    {
        //We use GSL for the integration calculating the average
        gsl_integration_workspace* w = gsl_integration_workspace_alloc( 3000 );
        double error;
        gsl_function F;
        F.function = &frogs_integrand_for_averaging; //Integrand function
        struct frogs_gsl_func_param par; //Parameters to be passed to the function
        par.ped = ped;
        par.exnoise = exnoise;
        par.mu = mu;
        par.probarray = prob_array;
        
        F.params = ( void* )&par;
        
        /* The integral must run over the domain where the probability density has
           the largest values. For the signal photons when the expectation value
           is mu, the variance is mu and this must be combined with the variance
           associated with the pedestal width. The variances are considered to
           add. */
        double stdev = sqrt( mu * ( 1 + exnoise * exnoise ) + ped * ped );
        double x_min = mu - FROGS_NUMBER_OF_SIGMA * stdev;
        double x_max = mu + FROGS_NUMBER_OF_SIGMA * stdev;
        //Proceed with the integration. Parameter values unclear
        gsl_integration_qag( &F, x_min, x_max, 1E-7, 1E-7, 3000, 3, w, &rtn, &error );
        gsl_integration_workspace_free( w );
    }
    return rtn;
}
//================================================================
//================================================================
double frogs_integrand_for_averaging( double q, void* par )
{
    /* This function calculates,-2*ln(P(q|mu)) * P(q|mu) where P(q|mu)
      is the probability density to get a signal q from a pixel for which
      the expectation value is mu. This is used in the calculation of the
      likelihood average.*/
    struct frogs_gsl_func_param* p = ( struct frogs_gsl_func_param* )par;
    
    double proba_density;
    proba_density = frogs_probability_density( q, p->mu, p->ped, p->exnoise );
    
    double loglikelihood = -2.0 * log( proba_density );
    
    if( !frogs_is_a_good_number( loglikelihood * proba_density ) )
    {
        frogs_showxerror( "Error: NaN resulted from calculations in frogs_integrand_for_averaging" );
    }
    //Returns the log-likelihood multiplied by the probability it is achieved
    return loglikelihood * proba_density;
}
//================================================================
//================================================================
/* This function reads the file whose name is specified by the variable
   FROGS_TEMPLATE_LIST, searching for the first one matching the elevation
   provided as an argument. The file should have one line for each template
   file to be considered. Each line contains in that order:
   1) The smallest elevation for which the template can be used
   2) The largest elevation for which the template can be used
   3) The file name for that template
*/
struct frogs_imgtemplate frogs_read_template_elev( float elevation, char templatelistname[FROGS_FILE_NAME_MAX_LENGTH] )
{
    char FROGS_TEMPLATE_LIST_PATH[500];
    
    char* itemp = getenv( "VERITAS_EVNDISP_AUX_DIR" );
    sprintf( FROGS_TEMPLATE_LIST_PATH, "%s/Frogs/%s", itemp, templatelistname );
    
    //Open the template files list file
    FILE* fu; //file pointer
    if( ( fu = fopen( FROGS_TEMPLATE_LIST_PATH, "r" ) ) == NULL )
    {
        fprintf( stdout, "FROGS template file list: %s\n", FROGS_TEMPLATE_LIST_PATH );
        frogs_showxerror( "Error: Failed opening the template files list file" );
    }

    fprintf( stdout, "reading new template for elevation %.2f\n", elevation );
    fprintf( stdout, "template list %s\n", FROGS_TEMPLATE_LIST_PATH );
    
    /*Read the file until the end is encountered unless a matching
      elevation range is found */
    int flag; /* flag=1 if convoled table is used
	       flag=0 if non-convoled tables are used */
    float minel, maxel; //Min and max elevation to use the listed files
    char fname[500];
    struct frogs_imgtemplate rtn; //Variable to hold template data
    while( fscanf( fu, "%d%f%f%s", &flag, &minel, &maxel, fname ) != EOF )
    {
#ifdef CONVOLUTION
        if( flag == 1 && elevation >= minel && elevation <= maxel )
        {
            fclose( fu ); //Closes the template filename list
            fprintf( stdout, "using template '%s'\n", fname ) ;
            rtn = frogs_read_template_file( fname );
            rtn.elevmin = minel;
            rtn.elevmax = maxel;
            return rtn;
        }
#endif
#ifndef CONVOLUTION
        if( flag == 0 && elevation >= minel && elevation <= maxel )
        {
            fclose( fu ); //Closes the template filename list
            fprintf( stdout, "using template '%s'\n", fname ) ;
            rtn = frogs_read_template_file( fname );
            rtn.elevmin = minel;
            rtn.elevmax = maxel;
            return rtn;
        }
#endif
    }
    fclose( fu ); //Closes the template filename list
    fprintf( stdout, "Elevation %f, check file %s\n", elevation, FROGS_TEMPLATE_LIST_PATH );
    frogs_showxerror( "Error: FROGS could not find a matching template file" );
    return rtn;
}
//================================================================
//================================================================
/* Read the image template data in a file whose name is received
   as an argument.
*/
struct frogs_imgtemplate frogs_read_template_file( char fname[FROGS_FILE_NAME_MAX_LENGTH] )
{
    FILE* fu; //file pointer
    
    char* itemp = getenv( "VERITAS_EVNDISP_AUX_DIR" );
    char fullfname[FROGS_FILE_NAME_MAX_LENGTH] ;
    sprintf( fullfname, "%s/Frogs/Templates/%s", itemp, fname ) ;
    
    if( ( fu = fopen( fullfname, "r" ) ) == NULL )
    {
        printf( "%s\n", fullfname );
        frogs_showxerror( "Error: Failed opening the template file" );
    }
    
    frogs_printfrog();
    fprintf( stdout, "Reading template file %s\n", fullfname );
    
    struct frogs_imgtemplate rtn;
    //Read the number of dimensions
    fscanf( fu, "%d", &( rtn.ndim ) );
    fprintf( stdout, "number of dimensions: %d\n", rtn.ndim );
    //Read the step sizes for each parameter
    rtn.step = new float[rtn.ndim];
    for( int i = 0; i < rtn.ndim; i++ )
    {
        fscanf( fu, "%f", &( rtn.step[i] ) );
    }
    //Read the lowest value for each parameter
    rtn.min = new float[rtn.ndim];
    for( int i = 0; i < rtn.ndim; i++ )
    {
        fscanf( fu, "%f", &( rtn.min[i] ) );
    }
    //Read the number of steps for each parameter
    rtn.nstep = new int[rtn.ndim];
    for( int i = 0; i < rtn.ndim; i++ )
    {
        fscanf( fu, "%d", &( rtn.nstep[i] ) );
    }
    //Calculate the number of points in the template data table
    rtn.sz = 1;
    for( int i = 0; i < rtn.ndim; i++ )
    {
        rtn.sz = rtn.sz * rtn.nstep[i];
    }
    // print parameter values
    if( rtn.ndim == 5 )
    {
           fprintf( stdout, "\tlambda (first interaction depth): \t min %.3f, step size %.3f, nsteps %d\n", rtn.min[0], rtn.step[0], rtn.nstep[0] );
           fprintf( stdout, "\tlog10(E/TeV): \t min %.3f, step size %.3f, nsteps %d\n", rtn.min[1], rtn.step[1], rtn.nstep[1] );
           fprintf( stdout, "\timpact parameter: \t min %.3f, step size %.3f, nsteps %d\n", rtn.min[2], rtn.step[2], rtn.nstep[2] );
           fprintf( stdout, "\tx coordinate: \t min %.3f, step size %.3f, nsteps %d\n", rtn.min[3], rtn.step[3], rtn.nstep[3] );
           fprintf( stdout, "\ty coordinate: \t min %.3f, step size %.3f, nsteps %d\n", rtn.min[4], rtn.step[4], rtn.nstep[4] );
    }
    fprintf( stdout, "expected number of points in the template data table: %d\n", rtn.sz );

    //Read the Cherenkov light densities
    rtn.c = new float [rtn.sz];
    for( int i = 0; i < rtn.sz; i++ )
    {
        if( i % ( int )floor( 0.1 * rtn.sz ) == 0 )
        {
            fprintf( stdout, "\t%d/100 reading template\n", ( int )floor( i * 100.0 / rtn.sz ) );
        }
        fscanf( fu, "%f", &( rtn.c[i] ) );
    } 
    fclose( fu );
    fprintf( stdout, "Finished reading template file\n" );
    
    return rtn;
}
//================================================================
//================================================================
void frogs_showxerror( const char* msg )
{
    /*Print the error message received as an argument and stops the analysis*/
    fprintf( stdout, "%s\n", msg );
    fprintf( stdout, "Terminating now\n" );
    exit( EXIT_FAILURE );
}
//================================================================
//================================================================
int frogs_likelihood( const gsl_vector* v, void* ptr, gsl_vector* f )
{
    /* Calculates the likelihood for each pixel and store the result in the gsl
       vector f. The event physical parameters for which the likelihood is
       calculated are provided in the gsl vector v. All the data from the
       telescope and from the template are available through the pointer to
       void ptr*/
    
    struct frogs_gsl_data_wrapper* dwrap = ( struct frogs_gsl_data_wrapper* )ptr;
    
    //Here is the parameter space point
    struct frogs_reconstruction pnt;
    pnt.xs = gsl_vector_get( v, FROGS_XS );
    pnt.ys = gsl_vector_get( v, FROGS_YS );
    pnt.xp = gsl_vector_get( v, FROGS_XP );
    pnt.yp = gsl_vector_get( v, FROGS_YP );
    pnt.log10e = gsl_vector_get( v, FROGS_LOG10E );
    pnt.lambda = gsl_vector_get( v, FROGS_LAMBDA );
    
    //This is where the model is called to calculate the expected
    //value for all pixels in the original version
    int gsl_pix_id = 0; //This counter is used as a pixel identified for gsl
    for( int tel = 0; tel < dwrap->data->ntel; tel++ )
    {
        for( int pix = 0; pix < dwrap->data->scope[tel].npix; pix++ )
        {
            if( dwrap->data->scope[tel].pixinuse[pix] == FROGS_OK )
            {
                int pix_in_template;//FROGS_OK in image, FROGS_NOTOK in background
                double mu = frogs_img_model( pix, tel, pnt, dwrap->data, dwrap->tmplt,
                                             &pix_in_template );
                if( mu != FROGS_BAD_NUMBER )
                {
                    double pd;
                    
                    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
                    {
                        pd = frogs_read_prob_array_table( dwrap->probarray,
                                                          dwrap->data->scope[tel].q[pix],
                                                          mu,
                                                          dwrap->data->scope[tel].ped[pix] );
                        if( pd < 0. )
                        {
                            printf( " tel %d pix %d q %f mu %f ped %f pd %f\n", tel, pix, dwrap->data->scope[tel].q[pix], mu, dwrap->data->scope[tel].ped[pix], pd );
                        }
                    }
                    else
                        pd = frogs_probability_density( dwrap->data->scope[tel].q[pix],
                                                        mu,
                                                        dwrap->data->scope[tel].ped[pix],
                                                        dwrap->data->scope[tel].exnoise[pix] );
                                                        
                                                        
                    double pix_lkhd = -2.0 * log( pd );
                    gsl_vector_set( f, gsl_pix_id, pix_lkhd );
                    gsl_pix_id++;
                    
                }
            }//End of test on pixel viability
        }//End of pixel loop
    }//End of telescope loop
    
    return GSL_SUCCESS;
}
//================================================================
//================================================================
int frogs_likelihood_derivative( const gsl_vector* v, void* ptr, gsl_matrix* J )
{
    /* Calculates the likelihood derivative for each pixel and with respect
      to each event physical parameter and store the result in the gsl
      matrix J. The event physical parameter for which the likelihood is
      calculated are provided in the gsl vector v. All the data from the
      telescope and from the template are available through the pointer to
      void ptr*/
    
    struct frogs_gsl_data_wrapper* dwrap = ( struct frogs_gsl_data_wrapper* )ptr;
    
    //Here is the parameter space point
    struct frogs_reconstruction pnt;
    pnt.xs = gsl_vector_get( v, FROGS_XS );
    pnt.ys = gsl_vector_get( v, FROGS_YS );
    pnt.xp = gsl_vector_get( v, FROGS_XP );
    pnt.yp = gsl_vector_get( v, FROGS_YP );
    pnt.log10e = gsl_vector_get( v, FROGS_LOG10E );
    pnt.lambda = gsl_vector_get( v, FROGS_LAMBDA );
    
    /*This is used to store the finite difference step magnitudes used in
      the estimation of the derivatives. Eventually their values could be
      coming in the frogs_gsl_data_wrapper structure*/
    struct frogs_reconstruction delta;
    
    delta.xs = dwrap->data->delta_xs;
    delta.ys = dwrap->data->delta_ys;
    delta.xp = dwrap->data->delta_xp;
    delta.yp = dwrap->data->delta_yp;
    delta.log10e = dwrap->data->delta_log10e;
    delta.lambda = dwrap->data->delta_lambda;
    
    int gsl_pix_id = 0; //This counter is used as a pixel identified for gsl
    for( int tel = 0; tel < dwrap->data->ntel; tel++ )
    {
        for( int pix = 0; pix < dwrap->data->scope[tel].npix; pix++ )
        {
            if( dwrap->data->scope[tel].pixinuse[pix] == FROGS_OK )
            {
                double derivative;
                //Derivative with respect to xs
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_XS, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_XS,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_XS);
                gsl_matrix_set( J, gsl_pix_id, FROGS_XS, derivative );
                
                //Derivative with respect to ys
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_YS, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_YS,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_YS);
                gsl_matrix_set( J, gsl_pix_id, FROGS_YS, derivative );
                
                //Derivative with respect to xp
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_XP, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_XP,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_XP);
                gsl_matrix_set( J, gsl_pix_id, FROGS_XP, derivative );
                
                //Derivative with respect to yp
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_YP, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_YP,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_YP);
                gsl_matrix_set( J, gsl_pix_id, FROGS_YP, derivative );
                
                //Derivative with respect to log10e
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_LOG10E, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_LOG10E,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data,dwrap->tmplt,FROGS_LOG10E);
                gsl_matrix_set( J, gsl_pix_id, FROGS_LOG10E, derivative );
                
                //Derivative with respect to lambda
                derivative = frogs_pix_lkhd_deriv_2ndorder( pix, tel, pnt, delta, dwrap->data,
                             dwrap->tmplt, FROGS_LAMBDA, dwrap->probarray );
                //derivative=frogs_pix_lkhd_deriv_4thorder(pix,tel,pnt,delta,
                //dwrap->data, dwrap->tmplt,FROGS_LAMBDA,dwrap->probarray);
                //derivative=frogs_pix_lkhd_deriv_4thorder_old(pix,tel,pnt,delta,
                //dwrap->data, dwrap->tmplt,FROGS_LAMBDA);
                gsl_matrix_set( J, gsl_pix_id, FROGS_LAMBDA, derivative );
                
                gsl_pix_id++;
            }//End of test on pixel viability
        }//End of pixel loop
    }//End of telescope loop
    
    return GSL_SUCCESS;
}
//================================================================
//================================================================
double frogs_pix_lkhd_deriv_2ndorder( int pix, int tel,
                                      struct frogs_reconstruction pnt,
                                      struct frogs_reconstruction delta,
                                      struct frogs_imgtmplt_in* d,
                                      struct frogs_imgtemplate* tmplt,
                                      int gsl_par_id,
                                      struct frogs_probability_array* prob_array )
{
    /*Calculates the derivative of the pixel likelihood for the pixel specified
      by tel and pix and with respect to the parameter specified by gsl_par_id.
      This function calculate the derivative as df/dx=[f(x+dx)-f(x-dx)]/(2*dx)*/
    
    double rtn;
    double mu;
    double pd;
    int pix_in_template;
    //Step forward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, 1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
                                        
    double pix_lkhd_plus = -2.0 * log( pd );
    //Step backward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, -1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
    double pix_lkhd_minus = -2.0 * log( pd );
    
    //Evaluate the derivative
    float delta_param = 0;
    if( gsl_par_id == FROGS_XS )
    {
        delta_param = delta.xs;
    }
    if( gsl_par_id == FROGS_YS )
    {
        delta_param = delta.ys;
    }
    if( gsl_par_id == FROGS_XP )
    {
        delta_param = delta.xp;
    }
    if( gsl_par_id == FROGS_YP )
    {
        delta_param = delta.yp;
    }
    if( gsl_par_id == FROGS_LOG10E )
    {
        delta_param = delta.log10e;
    }
    if( gsl_par_id == FROGS_LAMBDA )
    {
        delta_param = delta.lambda;
    }
    if( delta_param == 0 )
    {
        frogs_showxerror( "Error: Bad parameter identifier in pix_lkhd_deriv_2ndorder" );
    }
    rtn = ( pix_lkhd_plus - pix_lkhd_minus ) / ( 2.0 * delta_param );
    
    return rtn;
}
//================================================================
//================================================================
double frogs_pix_lkhd_deriv_4thorder( int pix, int tel,
                                      struct frogs_reconstruction pnt,
                                      struct frogs_reconstruction delta,
                                      struct frogs_imgtmplt_in* d,
                                      struct frogs_imgtemplate* tmplt,
                                      int gsl_par_id,
                                      struct frogs_probability_array* prob_array )
{
    /*Calculates the derivative of the pixel likelihood for the pixel specified
      by tel and pix and with respect to the parameter specified by gsl_par_id.
      This function calculate the derivative as
      df/dx=[-f(x+2*dx)+8f(x+dx)-8f(x-dx)+f(x-2dx)]/(12*dx)*/
    
    double rtn;
    double mu;
    double pd;
    int pix_in_template;
    //Step forward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, 1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
    double pix_lkhd_plus = -2.0 * log( pd );
    //Step backward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, -1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
    double pix_lkhd_minus = -2.0 * log( pd );
    //two Step forward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, 2.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
    double pix_lkhd_plus_plus = -2.0 * log( pd );
    //two Step backward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, -2.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    
    if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
        pd = frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu,
                                          d->scope[tel].ped[pix] );
    else
        pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                        d->scope[tel].ped[pix],
                                        d->scope[tel].exnoise[pix] );
                                        
    double pix_lkhd_minus_minus = -2.0 * log( pd );
    
    //Evaluate the derivative
    float delta_param = 0;
    if( gsl_par_id == FROGS_XS )
    {
        delta_param = delta.xs;
    }
    if( gsl_par_id == FROGS_YS )
    {
        delta_param = delta.ys;
    }
    if( gsl_par_id == FROGS_XP )
    {
        delta_param = delta.xp;
    }
    if( gsl_par_id == FROGS_YP )
    {
        delta_param = delta.yp;
    }
    if( gsl_par_id == FROGS_LOG10E )
    {
        delta_param = delta.log10e;
    }
    if( gsl_par_id == FROGS_LAMBDA )
    {
        delta_param = delta.lambda;
    }
    if( delta_param == 0 )
    {
        frogs_showxerror( "Error: Bad parameter identifier in pix_lkhd_deriv_2ndorder" );
    }
    rtn = ( -pix_lkhd_plus_plus + 8 * pix_lkhd_plus - 8 * pix_lkhd_minus + pix_lkhd_minus_minus ) / ( 12 * delta_param );
    
    return rtn;
}
//================================================================
//================================================================
double frogs_pix_lkhd_deriv_4thorder_old( int pix, int tel,
        struct frogs_reconstruction pnt,
        struct frogs_reconstruction delta,
        struct frogs_imgtmplt_in* d,
        struct frogs_imgtemplate* tmplt,
        int gsl_par_id )
{
    /*Calculates the derivative of the pixel likelihood for the pixel specified
      by tel and pix and with respect to the parameter specified by gsl_par_id.
      This function calculate the derivative as
      df/dx=[-f(x+2*dx)+8f(x+dx)-8f(x-dx)+f(x-2dx)]/(12*dx)*/
    
    double rtn;
    double mu;
    double pd;
    int pix_in_template;
    //Step forward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, 1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                    d->scope[tel].ped[pix],
                                    d->scope[tel].exnoise[pix] );
    double pix_lkhd_plus = -2.0 * log( pd );
    //Step backward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, -1.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                    d->scope[tel].ped[pix],
                                    d->scope[tel].exnoise[pix] );
    double pix_lkhd_minus = -2.0 * log( pd );
    //two Step forward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, 2.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                    d->scope[tel].ped[pix],
                                    d->scope[tel].exnoise[pix] );
    double pix_lkhd_plus_plus = -2.0 * log( pd );
    //two Step backward
    mu = frogs_img_model( pix, tel, frogs_param_step( pnt, delta, gsl_par_id, -2.0 ),
                          d, tmplt, &pix_in_template );
    if( mu == FROGS_BAD_NUMBER )
    {
        frogs_showxerror( "Error: frogs_img_model() invoked for an invalid pixel" );
    }
    pd = frogs_probability_density( d->scope[tel].q[pix], mu,
                                    d->scope[tel].ped[pix],
                                    d->scope[tel].exnoise[pix] );
    double pix_lkhd_minus_minus = -2.0 * log( pd );
    
    //Evaluate the derivative
    float delta_param = 0;
    if( gsl_par_id == FROGS_XS )
    {
        delta_param = delta.xs;
    }
    if( gsl_par_id == FROGS_YS )
    {
        delta_param = delta.ys;
    }
    if( gsl_par_id == FROGS_XP )
    {
        delta_param = delta.xp;
    }
    if( gsl_par_id == FROGS_YP )
    {
        delta_param = delta.yp;
    }
    if( gsl_par_id == FROGS_LOG10E )
    {
        delta_param = delta.log10e;
    }
    if( gsl_par_id == FROGS_LAMBDA )
    {
        delta_param = delta.lambda;
    }
    if( delta_param == 0 )
    {
        frogs_showxerror( "Error: Bad parameter identifier in pix_lkhd_deriv_2ndorder" );
    }
    rtn = ( -pix_lkhd_plus_plus + 8 * pix_lkhd_plus - 8 * pix_lkhd_minus + pix_lkhd_minus_minus ) / ( 12 * delta_param );
    
    return rtn;
}

//================================================================
//================================================================
struct frogs_reconstruction frogs_param_step( struct frogs_reconstruction pnt,
        struct frogs_reconstruction delta,
        int gsl_par_id, float mult )
{
    /*This function returns a point of the parameter space equal to pnt with
      a change applied to one of the parameters as specified by gsl_par_id
      with an amplitude equal to delta and multiplied by mult. It is used in
      the function frogs_pix_lkhd_deriv calculating the likelihood derivaive.*/
    struct frogs_reconstruction rtn;
    rtn = pnt;
    if( gsl_par_id == FROGS_XS )
    {
        rtn.xs = rtn.xs + mult * delta.xs;
    }
    if( gsl_par_id == FROGS_YS )
    {
        rtn.ys = rtn.ys + mult * delta.ys;
    }
    if( gsl_par_id == FROGS_XP )
    {
        rtn.xp = rtn.xp + mult * delta.xp;
    }
    if( gsl_par_id == FROGS_YP )
    {
        rtn.yp = rtn.yp + mult * delta.yp;
    }
    if( gsl_par_id == FROGS_LOG10E )
    {
        rtn.log10e = rtn.log10e+mult * delta.log10e;
    }
    if( gsl_par_id == FROGS_LAMBDA )
    {
        rtn.lambda = rtn.lambda + mult * delta.lambda;
    }
    
    return rtn;
}
//================================================================
//================================================================
int frogs_likelihood_fdf( const gsl_vector* v, void* ptr, gsl_vector* f,
                          gsl_matrix* J )
{
    /*It is not clear what this is good for but GSL needs to see that
      for the Levenberg-Marquardt optimisation*/
    frogs_likelihood( v, ptr, f );
    frogs_likelihood_derivative( v, ptr, J );
    return GSL_SUCCESS;
}
//================================================================
//================================================================
double frogs_img_model( int pix, int tel, struct frogs_reconstruction pnt,
                        struct frogs_imgtmplt_in* d, struct frogs_imgtemplate* tmplt,
                        int* intemplate )
{
    /* This function calculates and returns the expected signal in
       photoelectrons in pixel pix in telescope tel for the event physical
       parameters specified in pnt and using the telescope properties in d
       and the image template table tmplt. If the pixel is outside the region
       covered by the template, the value *intemplate is set to FROGS_NOTOK
       while otherwise it is set to FROGS_OK. If the pixel is invalid both the
       returned value and intemplate are set to  FROGS_BAD_NUMBER */
    
    //If the pixel is not active return FROGS_BAD_NUMBER
    if( d->scope[tel].pixinuse[pix] != FROGS_OK )
    {
        *intemplate = FROGS_BAD_NUMBER;
        return FROGS_BAD_NUMBER;
    }
    
    //Angle between the x direction and the shower image axis
    //(xfield, yfield) tel position in shower coordinate system
    //(xp, yp) shower position in shower coordinate system
    float phi = atan2( pnt.yp - d->scope[tel].yfield, pnt.xp - d->scope[tel].xfield );
    phi = phi + FROGS_PI; //This is for real data only. We'll have to understand that
    float cphi = cos( phi );
    float sphi = sin( phi );
    //Impact parameter to the telescope in the shower coordinate system
    float timp = sqrt( ( pnt.yp - d->scope[tel].yfield ) *
                       ( pnt.yp - d->scope[tel].yfield ) +
                       ( pnt.xp - d->scope[tel].xfield ) *
                       ( pnt.xp - d->scope[tel].xfield ) );
                       
    //Subtract the source coordinate from the pixel coordinate
    float xrs = d->scope[tel].xcam[pix] - pnt.xs;
    float yrs = -( d->scope[tel].ycam[pix] - pnt.ys );
    //Apply a rotation to move to the template coordinate system
    float tmpltxpix = xrs * cphi + yrs * sphi;
    float tmpltypix = -xrs * sphi + yrs * cphi;
    
    /*The optimization has a tendency to drag the value of lambda
      outside the range covered by the template table. In order to avoid
      the effects of this on the model, when the value is out of range
      we use the lambda value of the range that is the closest.*/
    float maxlambda = tmplt->min[0] + ( tmplt->nstep[0] - 1 ) * tmplt->step[0];
    pnt.lambda = floatwrap( pnt.lambda, tmplt->min[0], maxlambda );
    
    
    //Here we get the pixel value from the template table
    double rtn = 0.;
    
#ifdef CONVOLUTION
    if( d->interporder == 0 ) rtn = frogs_chertemplate_no_int( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate );
    if( d->interporder == 1 ) rtn = frogs_chertemplate_lin( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate );
    if( d->interporder == 2 ) rtn = frogs_chertemplate_quad( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate );

    /* The image template data is in photo-electrons per square meter and per
       square degree. We muliply this density by the telescope area and by the
       pixel area. */
    //rtn=rtn*d->scope[tel].telpixarea[pix];
    
    /* The new image templates are in photo-electron.*/
    //rtn=rtn;//*cone_eff;
#endif
#ifndef CONVOLUTION
    /*We get the pixel radius in degree*/
    //float pixradius=d->scope[tel].pixradius[pix];
    float pixradius = 0.080; //(0.075 - 0.085 - 0.080)
    
    if( d->interporder == 0 ) rtn = frogs_chertemplate_no_int( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate, pixradius );
    if( d->interporder == 1 ) rtn = frogs_chertemplate_lin( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate, pixradius );
    if( d->interporder == 2 ) rtn = frogs_chertemplate_quad( pnt.lambda, pnt.log10e,
                                        timp, tmpltxpix,
                                        tmpltypix, tmplt,
                                        intemplate, pixradius );
    if( rtn > 0. ) fprintf( stdout, "\t\t NOCONV %d %f\n", d->interporder, rtn );
    /* The image template data is in photo-electrons per square meter and per
       square degree. We muliply this density by the telescope area. */
    //rtn = rtn * telarea;
#endif
    
    /* The image templates do not take into PMT non-linearity and saturation.
    	We apply corrections to the mu-value to correct for that */
    rtn = frogs_mu_correction( rtn, d->epoch_id, d->lowerthresh, d->firstparam, d->secondparam );
    
    return rtn;
}
//================================================================
//================================================================
#ifdef CONVOLUTION
double frogs_chertemplate_lin( float lambda, float log10e, float b, float x,
                               float y, struct frogs_imgtemplate* tmplt,
                               int* intemplate )
{
#endif
#ifndef CONVOLUTION
    double frogs_chertemplate_lin( float lambda, float log10e, float b, float x,
                                   float y, struct frogs_imgtemplate * tmplt,
                                   int* intemplate, float pixradius )
    {
#endif
        /*This function return an evaluation of the Cherenkov ight density for
        the given values of lambda the depth of the first interaction point,
        log10e = log10(E/TeV), b the impact parameter to the telescope, x and
        y the longitudinal and transverse coordinate with respect to the source
        and the direction of development of the shower image. The evaluation is
        obtained by linear interpolation in lambda,log10e and b in that order.
        When the parameters fall outside the range covered by the template table,
        the value is linearly extrapolated.
        For x and y the closest table values are simply used.
        */
        /* Note:in the template data table, the index to parameter
         correspondance is as follows
         0 --- lambda = 1st interaction depth in interaction lengths
         1 --- log10e = log10(E/1TeV)
         2 --- b      = impact parameter to the considered telescope
         3 --- x      = x coordinate of the pixel (see note)
         4 --- y      = y coordinate of the pixel (see note)
         This has nothing to do with the order in which the parameters are
         entered in GSL for the likelihood optimization.
         (x and y measured in degrees. x along the image major axis from
         the source and increasing with the age of the shower. y in a
         perpendicular direction. )*/
        
        // index for x we will not interpolate
        int ix = ( int )floor( ( x - tmplt->min[3] ) / tmplt->step[3] );
        if( ix < 0 || ix >= tmplt->nstep[3] )
        {
            *intemplate = FROGS_NOTOK;
            return 0.0;
        }
        
        // index for y we will not interpolate
        int iy = ( int )floor( ( fabs( y ) - tmplt->min[4] ) / tmplt->step[4] );
        if( iy < 0 || iy >= tmplt->nstep[4] )
        {
            *intemplate = FROGS_NOTOK;
            return 0.0;
        }
        
        //If we get here the pixel is within the area covered by the templates
        *intemplate = FROGS_OK;
        
        /*For lambda we will proceed to a linear interpolation. We need two
        bracketing indices*/
        int il1 = ( int )floor( ( lambda - tmplt->min[0] ) / tmplt->step[0] );
        int il2 = il1 + 1;
        if( il1 < 0 )
        {
            il1 = 0;
            il2 = 1;
        }
        if( il2 >= tmplt->nstep[0] )
        {
            il2 = tmplt->nstep[0] - 1;
            il1 = il2 - 1;
        }
        float l1 = tmplt->min[0] + il1 * tmplt->step[0];
        float l2 = tmplt->min[0] + il2 * tmplt->step[0];
        
        //For the energy as well we will use a linear interpolation.
        int iloge1 = ( int )floor( ( log10e-tmplt->min[1] ) / tmplt->step[1] );
        int iloge2 = iloge1 + 1;
        if( iloge1 < 0 )
        {
            iloge1 = 0;
            iloge2 = 1;
        }
        if( iloge2 >= tmplt->nstep[1] )
        {
            iloge2 = tmplt->nstep[1] - 1;
            iloge1 = iloge2 - 1;
        }
        float loge1 = tmplt->min[1] + iloge1 * tmplt->step[1];
        float loge2 = tmplt->min[1] + iloge2 * tmplt->step[1];
        
        //For the impact parameter we will use a linear interpolation as well for now
        int ib1 = ( int )floor( ( b - tmplt->min[2] ) / tmplt->step[2] );
        int ib2 = ib1 + 1;
        if( ib1 < 0 )
        {
            ib1 = 0;
            ib2 = 1;
        }
        if( ib2 >= tmplt->nstep[2] )
        {
            ib2 = tmplt->nstep[2] - 1;
            ib1 = ib2 - 1;
        }
        float b1 = tmplt->min[2] + ib1 * tmplt->step[2];
        float b2 = tmplt->min[2] + ib2 * tmplt->step[2];
        
        //Get the model values at the vertices
        double mu111 = frogs_get_tmplt_val( il1, iloge1, ib1, ix, iy, tmplt );
        double mu112 = frogs_get_tmplt_val( il1, iloge1, ib2, ix, iy, tmplt );
        double mu121 = frogs_get_tmplt_val( il1, iloge2, ib1, ix, iy, tmplt );
        double mu122 = frogs_get_tmplt_val( il1, iloge2, ib2, ix, iy, tmplt );
        double mu211 = frogs_get_tmplt_val( il2, iloge1, ib1, ix, iy, tmplt );
        double mu212 = frogs_get_tmplt_val( il2, iloge1, ib2, ix, iy, tmplt );
        double mu221 = frogs_get_tmplt_val( il2, iloge2, ib1, ix, iy, tmplt );
        double mu222 = frogs_get_tmplt_val( il2, iloge2, ib2, ix, iy, tmplt );
        
        //Interpolation in lambda first
        double mu011 = frogs_linear_interpolation( l1, l2, mu111, mu211, lambda );
        mu011 = FROGS_NONEG( mu011 );
        double mu012 = frogs_linear_interpolation( l1, l2, mu112, mu212, lambda );
        mu012 = FROGS_NONEG( mu012 );
        double mu021 = frogs_linear_interpolation( l1, l2, mu121, mu221, lambda );
        mu021 = FROGS_NONEG( mu021 );
        double mu022 = frogs_linear_interpolation( l1, l2, mu122, mu222, lambda );
        mu022 = FROGS_NONEG( mu022 );
        
        //Interpolation in log(E)
        double mu001 = frogs_linear_interpolation( loge1, loge2, mu011, mu021, log10e );
        mu001 = FROGS_NONEG( mu001 );
        double mu002 = frogs_linear_interpolation( loge1, loge2, mu012, mu022, log10e );
        mu002 = FROGS_NONEG( mu002 );
        
        //Interpolation in b
        double mu000 = frogs_linear_interpolation( b1, b2, mu001, mu002, b );
        mu000 = FROGS_NONEG( mu000 );
        
        return mu000;
    }
    //================================================================
    //================================================================
#ifdef CONVOLUTION
    double frogs_chertemplate_quad( float lambda, float log10e, float b, float x,
                                    float y, struct frogs_imgtemplate * tmplt,
                                    int* intemplate )
    {
#endif
#ifndef CONVOLUTION
        double frogs_chertemplate_quad( float lambda, float log10e, float b, float x,
                                        float y, struct frogs_imgtemplate * tmplt,
                                        int* intemplate, float pixradius )
        {
#endif
            /*This function return an evaluation of the Cherenkov ight density for
              the given values of lambda the depth of the first interaction point,
              log10e = log10(E/TeV), b the impact parameter to the telescope, x and
              y the longitudinal and transverse coordinate with respect to the source
              and the direction of development of the shower image. The evaluation is
              obtained by quadratic interpolation in lambda,log10e and b in that order.
              When the parameters fall outside the range covered by the template table,
              the value is quadratically extrapolated.
              For x and y the closest table values are simply used.
            */
            /* Note:in the template data table, the index to parameter
            	correspondance is as follows
            	0 --- lambda = 1st interaction depth in intaraction lengths
            	1 --- log10e = log10(E/1TeV)
            	2 --- b      = impact parameter to the considered telescope
            	3 --- x      = x coordinate of the pixel (see note)
            	4 --- y      = y coordinate of the pixel (see note)
            	(note: x and y measured in degrees. x along the image major axis from
            	the source and increasing with the age of the shower. y in a
            	perpendicular direction. )*/
            
            // ixc, iyc: indices of pixel center
            int ixc = ( int )floor( ( x - tmplt->min[3] ) / tmplt->step[3] );
            if( ixc < 0 || ixc >= tmplt->nstep[3] )
            {
                *intemplate = FROGS_NOTOK;
                return 0.0;
            }
            
            int iyc = ( int )floor( ( fabs( y ) - tmplt->min[4] ) / tmplt->step[4] );
            if( iyc < 0 || iyc >= tmplt->nstep[4] )
            {
                *intemplate = FROGS_NOTOK;
                return 0.0;
            }
            
#ifndef CONVOLUTION
            // ixsup, ixinf, iysup, iyinf:  indices of pixel border
            int ixsup = ( int )floor( ( x + pixradius - tmplt->min[3] ) / tmplt->step[3] );
            if( ixsup >= tmplt->nstep[3] )
            {
                ixsup = tmplt->nstep[3] - 1;
            }
            
            int ixinf = ( int )floor( ( x - pixradius - tmplt->min[3] ) / tmplt->step[3] );
            if( ixinf < 0 )
            {
                ixinf = 0;
            }
            
            int iysup = ( int )floor( ( fabs( y ) + pixradius - tmplt->min[4] ) / tmplt->step[4] );
            if( iysup >= tmplt->nstep[4] )
            {
                iysup = tmplt->nstep[4] - 1;
            }
            
            int iyinf = ( int )floor( ( fabs( y ) - pixradius - tmplt->min[4] ) / tmplt->step[4] );
            if( iyinf < 0 )
            {
                iyinf = 0;
            }
#endif
            
            //If we get here the pixel is within the area covered by the templates
            *intemplate = FROGS_OK;
            
            float dummy; //temporary variable
            /*For lambda we will proceed to a quadratic interpolation. We need three
              bracketing indices*/
            int il1, il2, il3;
            dummy = ( lambda - tmplt->min[0] ) / tmplt->step[0];
            il1 = ( int )floor( dummy );
            il2 = il1 + 1;
            il3 = il2 + 1;
            if( il1 < 0 )
            {
                il1 = 0;
                il2 = 1;
                il3 = 2;
            }
            if( il3 >= tmplt->nstep[0] )
            {
                il3 = tmplt->nstep[0] - 1;
                il2 = il3 - 1;
                il1 = il2 - 1;
            }
            float l1 = tmplt->min[0] + il1 * tmplt->step[0];
            float l2 = tmplt->min[0] + il2 * tmplt->step[0];
            float l3 = tmplt->min[0] + il3 * tmplt->step[0];
            
            //For the energy as well we will use a quadratic interpolation as well.
            int iloge1, iloge2, iloge3;
            dummy = ( log10e-tmplt->min[1] ) / tmplt->step[1];
            iloge1 = ( int )floor( dummy );
            iloge2 = iloge1 + 1;
            iloge3 = iloge2 + 1;
            if( iloge1 < 0 )
            {
                iloge1 = 0;
                iloge2 = 1;
                iloge3 = 2;
            }
            if( iloge3 >= tmplt->nstep[1] )
            {
                iloge3 = tmplt->nstep[1] - 1;
                iloge2 = iloge3 - 1;
                iloge1 = iloge2 - 1;
            }
            float loge1 = tmplt->min[1] + iloge1 * tmplt->step[1];
            float loge2 = tmplt->min[1] + iloge2 * tmplt->step[1];
            float loge3 = tmplt->min[1] + iloge3 * tmplt->step[1];
            
            //For the impact parameter we will use a quadratic interpolation as well
            int ib1, ib2, ib3;
            dummy = ( b - tmplt->min[2] ) / tmplt->step[2];
            ib1 = ( int )floor( dummy );
            ib2 = ib1 + 1;
            ib3 = ib2 + 1;
            if( ib1 < 0 )
            {
                ib1 = 0;
                ib2 = 1;
                ib3 = 2;
            }
            if( ib3 >= tmplt->nstep[2] )
            {
                ib3 = tmplt->nstep[2] - 1;
                ib2 = ib3 - 1;
                ib1 = ib2 - 1;
            }
            float b1 = tmplt->min[2] + ib1 * tmplt->step[2];
            float b2 = tmplt->min[2] + ib2 * tmplt->step[2];
            float b3 = tmplt->min[2] + ib3 * tmplt->step[2];
            
            //Get the model values at the vertices of a 3x3x3 cube
            double mu111, mu112, mu113, mu121, mu122, mu123, mu131, mu132, mu133, mu211, mu212, mu213, mu221, mu222, mu223, mu231, mu232, mu233, mu311, mu312, mu313, mu321, mu322, mu323, mu331, mu332, mu333;
            mu111 = mu112 = mu113 = mu121 = mu122 = mu123 = mu131 = mu132 = mu133 = mu211 = mu212 = mu213 = mu221 = mu222 = mu223 = mu231 = mu232 = mu233 = mu311 = mu312 = mu313 = mu321 = mu322 = mu323 = mu331 = mu332 = mu333 = 0.;
            
#ifndef CONVOLUTION
            //set up GSL Random Number Generation
            unsigned long seed;
            srand( time( NULL ) );
            seed = rand();
            gsl_rng* r;
            gsl_rng_env_setup();
            r = gsl_rng_alloc( gsl_rng_mt19937 ); //Select random number generator
            gsl_rng_set( r, seed );
            
            float area = 0.;
            
            for( int ix = ixinf; ix <= ixsup; ix++ )
            {
                for( int iy = iyinf; iy <= iysup; iy++ )
                {
                    /*Get the intersection area between a pixel and a rectangle.
                    The lower left corner of the rectangle is located at the position (X,Y) in degrees*/
                    float X = ix * tmplt->step[3] + tmplt->min[3];
                    float Y = iy * tmplt->step[4] + tmplt->min[4];
                    area = frogs_get_overlapping_area( r, x, fabs( y ), pixradius, X, Y, tmplt->step[3], tmplt->step[4] );
                    if( area < 0 || area > ( tmplt->step[3] ) * ( tmplt->step[4] ) )
                    {
                        frogs_showxerror( "Error: Problem encountered in the intersection area calculation" );
                    }
                    
                    if( area > 0. )
                    {
                        mu111 += frogs_get_tmplt_val( il1, iloge1, ib1, ix, abs( iy ), tmplt ) * area;
                        mu112 += frogs_get_tmplt_val( il1, iloge1, ib2, ix, abs( iy ), tmplt ) * area;
                        mu113 += frogs_get_tmplt_val( il1, iloge1, ib3, ix, abs( iy ), tmplt ) * area;
                        mu121 += frogs_get_tmplt_val( il1, iloge2, ib1, ix, abs( iy ), tmplt ) * area;
                        mu122 += frogs_get_tmplt_val( il1, iloge2, ib2, ix, abs( iy ), tmplt ) * area;
                        mu123 += frogs_get_tmplt_val( il1, iloge2, ib3, ix, abs( iy ), tmplt ) * area;
                        mu131 += frogs_get_tmplt_val( il1, iloge3, ib1, ix, abs( iy ), tmplt ) * area;
                        mu132 += frogs_get_tmplt_val( il1, iloge3, ib2, ix, abs( iy ), tmplt ) * area;
                        mu133 += frogs_get_tmplt_val( il1, iloge3, ib3, ix, abs( iy ), tmplt ) * area;
                        mu211 += frogs_get_tmplt_val( il2, iloge1, ib1, ix, abs( iy ), tmplt ) * area;
                        mu212 += frogs_get_tmplt_val( il2, iloge1, ib2, ix, abs( iy ), tmplt ) * area;
                        mu213 += frogs_get_tmplt_val( il2, iloge1, ib3, ix, abs( iy ), tmplt ) * area;
                        mu221 += frogs_get_tmplt_val( il2, iloge2, ib1, ix, abs( iy ), tmplt ) * area;
                        mu222 += frogs_get_tmplt_val( il2, iloge2, ib2, ix, abs( iy ), tmplt ) * area;
                        mu223 += frogs_get_tmplt_val( il2, iloge2, ib3, ix, abs( iy ), tmplt ) * area;
                        mu231 += frogs_get_tmplt_val( il2, iloge3, ib1, ix, abs( iy ), tmplt ) * area;
                        mu232 += frogs_get_tmplt_val( il2, iloge3, ib2, ix, abs( iy ), tmplt ) * area;
                        mu233 += frogs_get_tmplt_val( il2, iloge3, ib3, ix, abs( iy ), tmplt ) * area;
                        mu311 += frogs_get_tmplt_val( il3, iloge1, ib1, ix, abs( iy ), tmplt ) * area;
                        mu312 += frogs_get_tmplt_val( il3, iloge1, ib2, ix, abs( iy ), tmplt ) * area;
                        mu313 += frogs_get_tmplt_val( il3, iloge1, ib3, ix, abs( iy ), tmplt ) * area;
                        mu321 += frogs_get_tmplt_val( il3, iloge2, ib1, ix, abs( iy ), tmplt ) * area;
                        mu322 += frogs_get_tmplt_val( il3, iloge2, ib2, ix, abs( iy ), tmplt ) * area;
                        mu323 += frogs_get_tmplt_val( il3, iloge2, ib3, ix, abs( iy ), tmplt ) * area;
                        mu331 += frogs_get_tmplt_val( il3, iloge3, ib1, ix, abs( iy ), tmplt ) * area;
                        mu332 += frogs_get_tmplt_val( il3, iloge3, ib2, ix, abs( iy ), tmplt ) * area;
                        mu333 += frogs_get_tmplt_val( il3, iloge3, ib3, ix, abs( iy ), tmplt ) * area;
                    }
                }
            }
#endif
            
#ifdef CONVOLUTION
            mu111 = frogs_get_tmplt_val( il1, iloge1, ib1, ixc, iyc, tmplt );
            mu112 = frogs_get_tmplt_val( il1, iloge1, ib2, ixc, iyc, tmplt );
            mu113 = frogs_get_tmplt_val( il1, iloge1, ib3, ixc, iyc, tmplt );
            mu121 = frogs_get_tmplt_val( il1, iloge2, ib1, ixc, iyc, tmplt );
            mu122 = frogs_get_tmplt_val( il1, iloge2, ib2, ixc, iyc, tmplt );
            mu123 = frogs_get_tmplt_val( il1, iloge2, ib3, ixc, iyc, tmplt );
            mu131 = frogs_get_tmplt_val( il1, iloge3, ib1, ixc, iyc, tmplt );
            mu132 = frogs_get_tmplt_val( il1, iloge3, ib2, ixc, iyc, tmplt );
            mu133 = frogs_get_tmplt_val( il1, iloge3, ib3, ixc, iyc, tmplt );
            mu211 = frogs_get_tmplt_val( il2, iloge1, ib1, ixc, iyc, tmplt );
            mu212 = frogs_get_tmplt_val( il2, iloge1, ib2, ixc, iyc, tmplt );
            mu213 = frogs_get_tmplt_val( il2, iloge1, ib3, ixc, iyc, tmplt );
            mu221 = frogs_get_tmplt_val( il2, iloge2, ib1, ixc, iyc, tmplt );
            mu222 = frogs_get_tmplt_val( il2, iloge2, ib2, ixc, iyc, tmplt );
            mu223 = frogs_get_tmplt_val( il2, iloge2, ib3, ixc, iyc, tmplt );
            mu231 = frogs_get_tmplt_val( il2, iloge3, ib1, ixc, iyc, tmplt );
            mu232 = frogs_get_tmplt_val( il2, iloge3, ib2, ixc, iyc, tmplt );
            mu233 = frogs_get_tmplt_val( il2, iloge3, ib3, ixc, iyc, tmplt );
            mu311 = frogs_get_tmplt_val( il3, iloge1, ib1, ixc, iyc, tmplt );
            mu312 = frogs_get_tmplt_val( il3, iloge1, ib2, ixc, iyc, tmplt );
            mu313 = frogs_get_tmplt_val( il3, iloge1, ib3, ixc, iyc, tmplt );
            mu321 = frogs_get_tmplt_val( il3, iloge2, ib1, ixc, iyc, tmplt );
            mu322 = frogs_get_tmplt_val( il3, iloge2, ib2, ixc, iyc, tmplt );
            mu323 = frogs_get_tmplt_val( il3, iloge2, ib3, ixc, iyc, tmplt );
            mu331 = frogs_get_tmplt_val( il3, iloge3, ib1, ixc, iyc, tmplt );
            mu332 = frogs_get_tmplt_val( il3, iloge3, ib2, ixc, iyc, tmplt );
            mu333 = frogs_get_tmplt_val( il3, iloge3, ib3, ixc, iyc, tmplt );
#endif
            
            //Interpolation in lambda first
            double mu011 = frogs_quadratic_interpolation( l1, l2, l3, mu111, mu211, mu311, lambda );
            mu011 = FROGS_NONEG( mu011 );
            double mu012 = frogs_quadratic_interpolation( l1, l2, l3, mu112, mu212, mu312, lambda );
            mu012 = FROGS_NONEG( mu012 );
            double mu013 = frogs_quadratic_interpolation( l1, l2, l3, mu113, mu213, mu313, lambda );
            mu013 = FROGS_NONEG( mu013 );
            double mu021 = frogs_quadratic_interpolation( l1, l2, l3, mu121, mu221, mu321, lambda );
            mu021 = FROGS_NONEG( mu021 );
            double mu022 = frogs_quadratic_interpolation( l1, l2, l3, mu122, mu222, mu322, lambda );
            mu022 = FROGS_NONEG( mu022 );
            double mu023 = frogs_quadratic_interpolation( l1, l2, l3, mu123, mu223, mu323, lambda );
            mu023 = FROGS_NONEG( mu023 );
            double mu031 = frogs_quadratic_interpolation( l1, l2, l3, mu131, mu231, mu331, lambda );
            mu031 = FROGS_NONEG( mu031 );
            double mu032 = frogs_quadratic_interpolation( l1, l2, l3, mu132, mu232, mu332, lambda );
            mu032 = FROGS_NONEG( mu032 );
            double mu033 = frogs_quadratic_interpolation( l1, l2, l3, mu133, mu233, mu333, lambda );
            mu033 = FROGS_NONEG( mu033 );
            
            //Interpolation in log(E)
            double mu001 = frogs_quadratic_interpolation( loge1, loge2, loge3, mu011, mu021,
                           mu031, log10e );
            mu001 = FROGS_NONEG( mu001 );
            double mu002 = frogs_quadratic_interpolation( loge1, loge2, loge3, mu012, mu022,
                           mu032, log10e );
            mu002 = FROGS_NONEG( mu002 );
            double mu003 = frogs_quadratic_interpolation( loge1, loge2, loge3, mu013, mu023,
                           mu033, log10e );
            mu003 = FROGS_NONEG( mu003 );
            
            //Interpolation in b
            double mu000 = frogs_quadratic_interpolation( b1, b2, b3, mu001, mu002, mu003, b );
            mu000 = FROGS_NONEG( mu000 );
            
#ifndef CONVOLUTION
            gsl_rng_free( r ); //Free the memory associated with r
#endif
            
            return mu000;
            
        }
        //================================================================
        //================================================================
        double frogs_get_tmplt_val( int il, int iloge, int ib, int ix, int iy,
                                    struct frogs_imgtemplate * tmplt )
        {
            /*The template table is stored in a one dimensional table. This
              function is used to return the template table content for the given
            	  indices for each of the entry parameters: il,iloge, ib, ix and iy*/
            /* Note:in the template data table, the index to parameter
              correspondance is as follows
              0 --- lambda = 1st interaction depth in intaraction lengths
              1 --- log10e = log10(E/1TeV)
              2 --- b      = impact parameter to the considered telescope
              3 --- x      = x coordinate of the pixel (see note)
              4 --- y      = y coordinate of the pixel (see note)*/
            
            int index = ( ( ( il * tmplt->nstep[1] + iloge ) * tmplt->nstep[2] + ib ) *
                      tmplt->nstep[3] + ix ) * tmplt->nstep[4] + iy;
            if( index < 0 || index >= tmplt->sz )
            {
                fprintf( stdout, "frogs_get_tmplt_val: %d %d %d %d %d\n", il, iloge, ib, ix, iy );
                frogs_showxerror( "Error: Index out of range in frogs_get_tmplt_val" );
            }
            return tmplt->c[index];
        }
        //================================================================
        //================================================================
        int frogs_print_param_spc_point( struct frogs_imgtmplt_out output )
        {
            /*This function can be used to print out the results or the current
              status of the image template analysis. */
            fprintf( stdout, "----------------------------------------------\n" );
            fprintf( stdout, "Event ID number:%d\n", output.event_id );
            fprintf( stdout, "GSL convergence status:%d\n", output.gsl_convergence_status );
            fprintf( stdout, "Number of iterations:%d\n", output.nb_iter );
            fprintf( stdout, "xs=%f +/- %f degrees\n", output.cvrgpt.xs,
                     output.cvrgpterr.xs );
            fprintf( stdout, "ys=%f +/- %f degrees\n", output.cvrgpt.ys,
                     output.cvrgpterr.ys );
            fprintf( stdout, "xp=%f +/- %f m\n", output.cvrgpt.xp, output.cvrgpterr.xp );
            fprintf( stdout, "yp=%f +/- %f m\n", output.cvrgpt.yp, output.cvrgpterr.yp );
            fprintf( stdout, "log10(E/TeV)=%f +/- %f\n", output.cvrgpt.log10e,
                     output.cvrgpterr.log10e );
            fprintf( stdout, "lambda=%f +/- %f g/cm2\n", output.cvrgpt.lambda,
                     output.cvrgpterr.lambda );
            fprintf( stdout, "..............................................\n" );
            fprintf( stdout, "Image goodness Gi=%f for %d d.o.f.\n", output.goodness_img,
                     output.npix_img );
            fprintf( stdout, "Bkgnd goodness Gb=%f for %d d.o.f.\n", output.goodness_bkg,
                     output.npix_bkg );
            fprintf( stdout, "----------------------------------------------\n" );
            
            return FROGS_OK;
        }
        //================================================================
        //================================================================
        double frogs_linear_interpolation( float x1, float x2, double y1, double y2,
                                           float x )
        {
            /*Returns the value in x of the first degree Lagrange polynomial
              constructed on points (x1,y1) and (x2,y2)*/
            double rtn;
            if( x1 == x2 )
            {
                frogs_showxerror( "Error: Problem encountered in two point linear lagrange interpolation" );
            }
            rtn = y1 * ( x - x2 ) / ( x1 - x2 ) + y2 * ( x - x1 ) / ( x2 - x1 );
            return rtn;
        }
        //================================================================
        //================================================================
        double frogs_quadratic_interpolation( float x1, float x2, float x3, double y1,
                                              double y2, double y3, float x )
        {
            /*Returns the value in x of the second degree Lagrange polynomial
              constructed on points (x1,y1), (x2,y2) and (x3,y3)*/
            double rtn;
            if( x1 == x2 || x1 == x3 || x2 == x3 )
            {
                frogs_showxerror( "Error: Problem encountered in three point quadratic Lagrange interpolation" );
            }
            rtn = y1 * ( x - x2 ) * ( x - x3 ) / ( ( x1 - x2 ) * ( x1 - x3 ) )
                  + y2 * ( x - x1 ) * ( x - x3 ) / ( ( x2 - x1 ) * ( x2 - x3 ) )
                  + y3 * ( x - x1 ) * ( x - x2 ) / ( ( x3 - x1 ) * ( x3 - x2 ) );
            return rtn;
        }
        //================================================================
        //================================================================
        float floatwrap( float x, float min, float max )
        {
            /*This function returns a value between min and max. It is periodic
              in x of periode 2*(max-min). For x in [min,max] the returned value
              is x. For x in [max,2*max-min], the return value is 2*max-x.
              This function is used on parameters controled by the GSL optimization
              to maintain the parameter within the range for which we have models
              without introducing any discontinuity
            
            rtn ^
              o |         o           o
            o   o       o | o       o   o       o
                | o   o   |   o   o       o   o
                |   o     |     o           o
                | . |     |     |
              0 .---|-----|-----|---------------------> x
               0   min   max   2*max-min
            */
            return x; //(SV)
            
            float rtn;
            if( min > max )
            {
                frogs_showxerror( "Error: In floatwrap, min>max" );
            }
            float range = max - min;
            float dumx = x - min;
            rtn = min + range - fabs( dumx - range * ( 2 * floor( 0.5 * dumx / range ) + 1 ) );
            return rtn;
        }
        
        //================================================================
        //================================================================
        float frogs_get_overlapping_area( gsl_rng * r, float x, float y, float pixradius,
                                          float X, float Y, float dX, float dY )
        {
            /*This function gets the area of intersection between a rectangle
              and a circle using a Monte-Carlo estimation.
              We generate a total number of random points N in the box and
              test if they are inside the circle . The intersection area is dX*dY*(n/N),
              where n is the number of points that lie in the circle and the rectangle.
            
              4 *-------* 3    y^
                |       |       |
                |       |       |
                |       |       ---->x
              1 *-------* 2
            
              o=(X,Y),
              1=(X,Y), 2=(X+dX,Y), 3=(X+dX,Y+dY), 4=(X,Y+dY) */
            
            int nIN = 0; //count the number of points lying within the rectangle and the pixel
            int nTOTAL = 50;
            
            for( int i = 0; i < nTOTAL; i++ )
            {
                // rnd1, rnd2: random numbers [-1,1]
                //float rnd1=2.*gsl_rng_uniform (r)-1.;
                //float rnd2=2.*gsl_rng_uniform (r)-1.;
                float rnd1 = gsl_rng_uniform( r );
                float rnd2 = gsl_rng_uniform( r );
                if( rnd1 < 0 || rnd1 > 1 || rnd2 < 0 || rnd2 > 1 )
                {
                    fprintf( stdout, "frogs_get_overlapping_area: random number <0 or >1\n" );
                    exit( 0 );
                }
                
                /*coordinates of a random point pixel (X0, Y0)*/
                //float X0=X+rnd1*dX/2.;
                //float Y0=Y+rnd2*dY/2.;
                float X0 = X + rnd1 * dX;
                float Y0 = Y + rnd2 * dY;
                /*pixel coordinates (x, y)*/
                float D = sqrt( pow( X0 - x, 2. ) + pow( Y0 - y, 2. ) );
                if( D <= pixradius )
                {
                    nIN += 1;
                }
            }
            return dX * dY * nIN / nTOTAL;
        }
        
        //================================================================
        //================================================================
        void frogs_fill_prob_density( struct frogs_probability_array * parray )
        {
            /*This function fills the probability density table to speed up the calculations
              It is filled once at the beginning of the analysis. */
            
            fprintf( stdout, "\nFROGS: Filling Probability Density Table ......... \n" );
            int bincount = 0 ;
            for( int i = 0; i < BIN1; i++ )
            {
                double q = ( ( double )i * ( RANGE1 - MIN1 ) / BIN1 + MIN1 );
                fprintf( stdout, "\t bin # %5d, q=%f\n", bincount, q ) ;
                bincount++ ;
                for( int j = 0; j < BIN2; j++ )
                {
                    double mu = ( ( double )j * ( RANGE2 - MIN2 ) / BIN2 + MIN2 );
                    for( int k = 0; k < BIN3; k++ )
                    {
                        double ped = ( double )k * ( RANGE3 - MIN3 ) / BIN3 + MIN3;
                        //printf( "  ped = %f\n", ped ) ;
                        parray->prob_density_table[i][j][k] = frogs_probability_density( q, mu, ped, 0.35 );
                        //printf( "  found probability density for bin [%d][%d][%d]\n", i, j, k ) ;
                        if( parray->prob_density_table[i][j][k] < 0. )
                        {
                            printf( "q %f mu %f ped %f pd %f\n", q, mu, ped, parray->prob_density_table[i][j][k] );
                        }
                    }
                }
            }
            
            fprintf( stdout, "completed filling of probability densitivity table with %d x %d x %d entries.\n\n", BIN1, BIN2, BIN3 );
            
        }
        
        //================================================================
        //================================================================
        double frogs_read_prob_array_table( struct frogs_probability_array * prob_array, double q, double mu, double ped )
        {
        
            /* Reads the probability table and extrapolates in 3 dimensions to find correct probability density.*/
            
            float ii, jj, kk;
            int imin, jmin, kmin;
            int imax, jmax, kmax;
            
            ii = BIN1 * ( q - MIN1 ) / ( RANGE1 - MIN1 );
            jj = BIN2 * ( mu - MIN2 ) / ( RANGE2 - MIN2 );
            kk = BIN3 * ( ped - MIN3 ) / ( RANGE3 - MIN3 );
            
            imin = ( int )ii;
            jmin = ( int )jj;
            kmin = ( int )kk;
            
            imax = imin + 1;
            jmax = jmin + 1;
            kmax = kmin + 1;
            
            if( ii <= 0 || kk <= 0 || jj <= 0 || imax >= BIN1 - 1 || jmax >= BIN2 - 1 || kmax >= BIN3 - 1 )
            {
                return frogs_probability_density( q, mu, ped, 0.35 );
            }
            
            float qmin = imin * ( RANGE1 - MIN1 ) / BIN1 + MIN1;
            float mumin = jmin * ( RANGE2 - MIN2 ) / BIN2 + MIN2;
            float pedmin = kmin * ( RANGE3 - MIN3 ) / BIN3 + MIN3;
            
            float qmax = imax * ( RANGE1 - MIN1 ) / BIN1 + MIN1;
            float mumax = jmax * ( RANGE2 - MIN2 ) / BIN2 + MIN2;
            float pedmax = kmax * ( RANGE3 - MIN3 ) / BIN3 + MIN3;
            
            float qdel = ( q - qmin ) / ( qmax - qmin );
            float mudel = ( mu - mumin ) / ( mumax - mumin );
            float peddel = ( ped - pedmin ) / ( pedmax - pedmin );
            
            double f000 = prob_array->prob_density_table[imin][jmin][kmin];
            double f001 = prob_array->prob_density_table[imin][jmin][kmax];
            double f010 = prob_array->prob_density_table[imin][jmax][kmin];
            double f011 = prob_array->prob_density_table[imin][jmax][kmax];
            double f100 = prob_array->prob_density_table[imax][jmin][kmin];
            double f101 = prob_array->prob_density_table[imax][jmin][kmax];
            double f110 = prob_array->prob_density_table[imax][jmax][kmin];
            double f111 = prob_array->prob_density_table[imax][jmax][kmax];
            
            double c0 = f000;
            double c1 = f100 - f000;
            double c2 = f010 - f000;
            double c3 = f001 - f000;
            double c4 = f110 - f010 - f100 + f000;
            double c5 = f011 - f001 - f010 + f000;
            double c6 = f101 - f001 - f100 + f000;
            double c7 = f111 - f011 - f101 - f110 + f100 + f001 + f010 - f000;
            
            double p = c0;
            p += c1 * qdel;
            p += c2 * mudel;
            p += c3 * peddel;
            p += c4 * qdel * mudel;
            p += c5 * mudel * peddel;
            p += c6 * peddel * qdel;
            p += c7 * qdel * mudel * peddel;
            
            return p;
        }
        
        //================================================================
        //================================================================
        
        double frogs_mu_correction( double mu, const char epoch_id[20],
                                    double mu_correction_lower_threshold,
                                    double mu_correction_first_parameter,
                                    double mu_correction_second_parameter )
        {
            /* this function applies correction to mu
              to correct for non-linearity and saturation effects */
            if( mu > mu_correction_lower_threshold )
            {
                return mu_correction_first_parameter + mu_correction_second_parameter * log( mu );
            }
            return mu;
            
            /*if( strcmp( epoch_id, "V6" ) == 0 )
             {
            	if( mu > mu_correction_lower_threshold[0] )
            	{
            		return mu_correction_first_parameter[0] + mu_correction_second_parameter[0] * log( mu );
            		//return -2.9628E3 + 5.566444E2 * log( mu );
            	}
            }
            else if( strcmp( epoch_id, "V5" ) == 0 )
            {
            	if( mu > mu_correction_lower_threshold[1] )
            	{
            		return mu_correction_first_parameter[1] + mu_correction_second_parameter[1] * log( mu );
            	}
            }
            return mu;*/
        }
        
        //================================================================
        //================================================================
        int frogs_printfrog()
        {
            /*This function clearly is the most important one of the project. */
            fprintf( stdout, " ---------------------------------------------------------------- \n" );
            fprintf( stdout, "|      _    _                                        _    _      |\n" );
            fprintf( stdout, "|     (o)--(o)                                      (o)--(o)     |\n" );
            fprintf( stdout, "|    /.______.\\    FFFF RRRR   OOO   GGG   SSS     /.______.\\    |\n" );
            fprintf( stdout, "|    \\________/    F    R   R O   O G     S        \\________/    |\n" );
            fprintf( stdout, "|  ./        \\.    FFF  RRRR  O   O G  GG  SSS     ./        \\.  |\n" );
            fprintf( stdout, "| ( .        , )   F    R R   O   O G   G     S   ( .        , ) |\n" );
            fprintf( stdout, "|  \\ \\_\\\\//_/ /    F    R  RR  OOO   GGG  SSSS     \\ \\_\\\\//_/ /  |\n" );
            fprintf( stdout, "|   ~~  ~~  ~~                                      ~~  ~~  ~~   |\n" );
            fprintf( stdout, "| svincent@physics.utah.edu             lebohec@physics.utah.edu |\n" );
            fprintf( stdout, "|                 VERSION 1.02 OCTOBER 10th 2011                 |\n" );
            fprintf( stdout, "|  For license issues, see www.physics.utah.edu/gammaray/FROGS   |\n" );
            fprintf( stdout, " ---------------------------------------------------------------- \n" );
            fprintf( stdout, "(modified by the Eventdisplay team)\n" );
            return FROGS_OK;
        }
        
//================================================================
//================================================================

void frogs_differential_evolution( struct frogs_imgtmplt_in * d,
                                   struct frogs_imgtemplate * tmplt,
                                   struct frogs_probability_array * prob_array )
{
        
#define URN_DEPTH   5   //4 + one index to avoid
        
            //variable declarations
            int i_r1, i_r2, i_r3;  // placeholders for random indexes
            
            int gi_gen; // generation counter
            int i_genmax = gi_genmax;
            
            int ia_urn2[URN_DEPTH];
            
            float fa_minbound[MAXDIM];
            float fa_maxbound[MAXDIM];
            t_pop t_tmp, t_bestit;
#ifdef BOUND_CONSTR
            t_pop t_origin;
#endif//BOUND_CONSTR
            float f_jitter, f_dither;
            
            /* Initialization of annealing parameters
               0 --- lambda = 1st interaction depth in interaction lengths
               1 --- log10e = log10(E/1TeV)
               2 --- b      = impact parameter to the considered telescope
               3 --- x      = x coordinate of the pixel (see note)
               4 --- y      = y coordinate of the pixel (see note) */
            
            /************* frogs *******************/
            int gi_D = tmplt->ndim + 1; //Dimension of parameter vector
            //gi_D=6 {xs, ys, xp, yp, logE, lambda}
            if( gi_D > MAXDIM )
            {
                frogs_showxerror( "Error: Problem encountered in the differential evolution: Error! too many parameters" );
            }
            
            if( gi_NP > MAXPOP )
            {
                frogs_showxerror( "Error: Problem encountered in the differential evolution: Error! too many points" );
            }
            
            /*fa_minbound[0]=-.5; fa_maxbound[0]=+.5;//xs
            fa_minbound[1]=-.5; fa_maxbound[1]=+.5;//ys
            
            fa_minbound[2]=-300.; fa_maxbound[2]=+300.;//xp
            fa_minbound[3]=-300.; fa_maxbound[3]=+300.;//yp
            
            fa_minbound[4]=tmplt->min[1]; fa_maxbound[4]=fa_minbound[1]+(tmplt->nstep[1]-1)*tmplt->step[1];//log10e
            
            fa_minbound[5]=tmplt->min[0]; fa_maxbound[5]=fa_minbound[0]+(tmplt->nstep[0]-1)*tmplt->step[0];//lambda
            
            fprintf(stdout,"in frogs.c, xs %f ys %f xp %f yp %f logE %f lambda %f\n",
              d->startpt.xs,d->startpt.ys,d->startpt.xp,d->startpt.yp,d->startpt.log10e,d->startpt.lambda);*/
            
            float minbound = 0.9;
            float maxbound = 1.1;
            
            fa_minbound[0] = minbound * d->startpt.xs;
            fa_maxbound[0] = maxbound * d->startpt.xs;
            if( d->startpt.xs < 0 )
            {
                fa_minbound[0] = maxbound * d->startpt.xs;
                fa_maxbound[0] = minbound * d->startpt.xs;
            }
            
            fa_minbound[1] = minbound * d->startpt.ys;
            fa_maxbound[1] = maxbound * d->startpt.ys;
            if( d->startpt.ys < 0 )
            {
                fa_minbound[1] = maxbound * d->startpt.ys;
                fa_maxbound[1] = minbound * d->startpt.ys;
            }
            
            minbound = 0.5;
            maxbound = 1.5;
            
            fa_minbound[2] = minbound * d->startpt.xp;
            fa_maxbound[2] = maxbound * d->startpt.xp;
            if( d->startpt.xp < 0 )
            {
                fa_minbound[2] = maxbound * d->startpt.xp;
                fa_maxbound[2] = minbound * d->startpt.xp;
            }
            
            fa_minbound[3] = minbound * d->startpt.yp;
            fa_maxbound[3] = maxbound * d->startpt.yp;
            if( d->startpt.yp < 0 )
            {
                fa_minbound[3] = maxbound * d->startpt.yp;
                fa_maxbound[3] = minbound * d->startpt.yp;
            }
            
            if( fabs( d->startpt.log10e-FROGS_BAD_NUMBER ) < 1E-8 )
            {
                fa_minbound[4] = -1.2;    //log10(-1.2)=60GeV
                fa_maxbound[4] = 1.;
            }
            else
            {
                minbound = 0.9;
                maxbound = 1.1;
                double e = pow( 10, d->startpt.log10e );
                double emin = e * minbound;
                double emax = e * maxbound;
                fa_minbound[4] = log10( emin );
                fa_maxbound[4] = log10( emax );
            }
            
            fa_minbound[5] = minbound * d->startpt.lambda;
            fa_maxbound[5] = maxbound * d->startpt.lambda;
            
            /************* frogs *******************/
            
            /************* sphere *******************/
            /*
            int gi_D=2;  //Dimension of parameter vector
            
            fa_minbound[0]=-5.12; fa_maxbound[0]=+5.12;
            fa_minbound[1]=-5.12; fa_maxbound[1]=+5.12;
            */
            /************* sphere *******************/
            
            /************* Shelkel's foxholes 2 *******************/
            /*
            int gi_D=2;  //Dimension of parameter vector
            
            fa_minbound[0]=-0; fa_maxbound[0]=+10;
            fa_minbound[1]=-0; fa_maxbound[1]=+10;
            */
            /************* Shelkel's foxholes 2 *******************/
            
            /************* Shelkel's foxholes 10 *******************/
            /*
            int gi_D=10;  //Dimension of parameter vector
            
            fa_minbound[0]=-0; fa_maxbound[0]=+10;
            fa_minbound[1]=-0; fa_maxbound[1]=+10;
            fa_minbound[2]=-0; fa_maxbound[2]=+10;
            fa_minbound[3]=-0; fa_maxbound[3]=+10;
            fa_minbound[4]=-0; fa_maxbound[4]=+10;
            fa_minbound[5]=-0; fa_maxbound[5]=+10;
            fa_minbound[6]=-0; fa_maxbound[6]=+10;
            fa_minbound[7]=-0; fa_maxbound[7]=+10;
            fa_minbound[8]=-0; fa_maxbound[8]=+10;
            fa_minbound[9]=-0; fa_maxbound[9]=+10;
            */
            /************* Shelkel's foxholes 10 *******************/
            
            
            //for (int i=0; i<gi_D; i++)
            //fprintf(stdout,"i %d fa_minbound %f fa_maxbound %f\n", i, fa_minbound[i], fa_maxbound[i]);
            
            if( !frogs_gsl_rng ) //something must have gone wrong earlier. We can try to reseed with the system time...
            {
                fprintf( stdout, "frogs_differential_evolution Warning: Random number generator not set, try to set and reseed with system time.\n" );
                frogs_seed_gsl_rng( 0 );
                if( !frogs_gsl_rng )
                {
                    frogs_showxerror( "Could not set random number generator for differential evolution." );
                }
            }
            
            
            long  gl_nfeval = 0;     //number of function evaluations
            t_pop gta_pop[2 * MAXPOP]; //the two populations are put into one array side by side
            t_pop gt_best;           // current best population member
            t_pop* gpta_old, *gpta_new, *gpta_swap;
            
            //evaluate the likelihood for the ED starting configuration
            gta_pop[0].fa_vector[0] = d->startpt.xs;
            gta_pop[0].fa_vector[1] = d->startpt.ys;
            gta_pop[0].fa_vector[2] = d->startpt.xp;
            gta_pop[0].fa_vector[3] = d->startpt.yp;
            gta_pop[0].fa_vector[4] = d->startpt.log10e;
            gta_pop[0].fa_vector[5] = d->startpt.lambda;
            double lkhd0 = fabs( ( float )FROGS_BAD_NUMBER ); //FROGS_BAD_NUMBER = -9999
            if( d->startpt.log10e != FROGS_BAD_NUMBER )
            {
                gta_pop[0] = frogs_evaluate( d, tmplt, prob_array, gi_D, gta_pop[0], &gl_nfeval, &gta_pop[0], gi_NP );
                lkhd0 = gta_pop[0].fa_cost[0];
            }
            
            //initialization
            for( int i = 0; i < gi_D; i++ )
            {
                /* gi_D corresponds to the total number of parameters
                   gta_pop[0].fa_vector[i] is the ith component of a random vector
                   the components of the vector are randomly distributed between
                   fa_maxbound[i] and fa_minbound[i] */
                gta_pop[0].fa_vector[i] = fa_minbound[i] + gsl_rng_uniform( frogs_gsl_rng ) * ( fa_maxbound[i] - fa_minbound[i] );
                //fprintf(stdout,"i %d %f %f gta %f\n",
                //i, fa_minbound[i], fa_maxbound[i], gta_pop[0].fa_vector[i]);
            }
            
            gta_pop[0] = frogs_evaluate( d, tmplt, prob_array, gi_D, gta_pop[0], &gl_nfeval, &gta_pop[0], gi_NP );
            gt_best = gta_pop[0];
            
            for( int i = 1; i < gi_NP; i++ )
            {
                for( int j = 0; j < gi_D; j++ )
                {
                    gta_pop[i].fa_vector[j] = fa_minbound[j] + gsl_rng_uniform( frogs_gsl_rng ) * ( fa_maxbound[j] - fa_minbound[j] );
                }
                gta_pop[i] = frogs_evaluate( d, tmplt, prob_array, gi_D, gta_pop[i], &gl_nfeval, &gta_pop[0], gi_NP );
                
                if( frogs_left_vector_wins( gta_pop[i], gt_best ) == TRUE )
                {
                    gt_best = gta_pop[i];
                }
            }
            
            t_bestit = gt_best;
            
            //---assign pointers to current ("old") and new population---
            
            gpta_old = &gta_pop[0];
            gpta_new = &gta_pop[gi_NP];
            
            //------Iteration loop--------------------------------------------
            gi_gen = 0;
#ifdef DO_PLOTTING
            while( ( gi_gen < i_genmax ) && ( gi_exit_flag == 0 ) ) // && (gt_best.fa_cost[0] > VTR))
#else
            //Note that kbhit() needs conio.h which is not always available under Unix.
            while( ( gi_gen < i_genmax ) ) // && (kbhit() == 0))// && (gt_best.fa_cost[0] > VTR))
#endif//DO_PLOTTING
            {
                gi_gen++;
                
                //----computer dithering factor (if needed)-----------------
                f_dither = f_weight + gsl_rng_uniform( frogs_gsl_rng ) * ( 1.0 - f_weight );
                
                //----start of loop through ensemble------------------------
                for( int i = 0; i < gi_NP; i++ )
                {
                    frogs_permute( frogs_gsl_rng, ia_urn2, URN_DEPTH, gi_NP, i ); //Pick 4 random and distinct
                    i_r1 = ia_urn2[1];                           //population members
                    i_r2 = ia_urn2[2];
                    i_r3 = ia_urn2[3];
                    
                    /*
                    //---this is an alternative way to pick population members---
                    do                        // Pick a random population member
                    {
                    i_r1 = (int)(genrand()*gi_NP);
                    }while(i_r1==i);
                    
                    do                        // Pick a random population member
                    {
                    i_r2 = (int)(genrand()*gi_NP);
                    }while((i_r2==i) || (i_r2==i_r1));
                    
                    do                        // Pick a random population member
                    {
                    i_r3 = (int)(genrand()*gi_NP);
                    }while((i_r3==i) || (i_r3==i_r1) || (i_r3==i_r2));
                    
                    */
                    
                    //========Choice of strategy=======================================================
                    //---classical strategy DE/rand/1/bin-----------------------------------------
                    if( gi_strategy == 1 )
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        do
                        {
                            // add fluctuation to random target
                            t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + f_weight * ( gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j] );
                            
                            j = ( j + 1 ) % gi_D;
                            k++;
                        }
                        while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, gpta_old[i_r1].fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR
                    }
                    //---DE/local-to-best/1/bin---------------------------------------------------
                    else if( gi_strategy == 2 )
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        do
                        {
                            // add fluctuation to random target
                            t_tmp.fa_vector[j] =
                                t_tmp.fa_vector[j] + f_weight * ( t_bestit.fa_vector[j] - t_tmp.fa_vector[j] ) +
                                f_weight * ( gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j] );
                                
                            j = ( j + 1 ) % gi_D;
                            k++;
                        }
                        while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, t_tmp.fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR
                    }
                    //---DE/best/1/bin with jitter------------------------------------------------
                    else if( gi_strategy == 3 )
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        do
                        {
                            // add fluctuation to random target
                            f_jitter = ( 0.0001 * gsl_rng_uniform( frogs_gsl_rng ) + f_weight );
                            t_tmp.fa_vector[j] =
                                t_bestit.fa_vector[j]
                                + f_jitter * ( gpta_old[i_r1].fa_vector[j] - gpta_old[i_r2].fa_vector[j] );
                                
                            j = ( j + 1 ) % gi_D;
                            k++;
                        }
                        while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, t_tmp.fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR
                    }
                    //---DE/rand/1/bin with per-vector-dither-------------------------------------
                    else if( gi_strategy == 4 )
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        do
                        {
                            // add fluctuation to random target
                            
                            t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
                                                 ( f_weight + gsl_rng_uniform( frogs_gsl_rng ) * ( 1.0 - f_weight ) ) *
                                                 ( gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j] );
                                                 
                            j = ( j + 1 ) % gi_D;
                            k++;
                        }
                        while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, t_tmp.fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR
                    }
                    //---DE/rand/1/bin with per-generation-dither---------------------------------
                    else if( gi_strategy == 5 )
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        do
                        {
                            // add fluctuation to random target
                            t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j]
                                                 + f_dither * ( gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j] );
                                                 
                            j = ( j + 1 ) % gi_D;
                            k++;
                        }
                        while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, t_tmp.fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR
                    }
                    //---variation to DE/rand/1/bin: either-or-algorithm--------------------------
                    else
                    {
                        frogs_assigna2b( gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector );
                        int j = ( int )( gsl_rng_uniform( frogs_gsl_rng ) * gi_D ); // random parameter
                        int k = 0;
                        if( gsl_rng_uniform( frogs_gsl_rng ) < 0.5 ) //Pmu = 0.5
                        {
                            //differential mutation
                            do
                            {
                                // add fluctuation to random target
                                t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j]
                                                     + f_weight * ( gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j] );
                                                     
                                j = ( j + 1 ) % gi_D;
                                k++;
                            }
                            while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
                        }
                        else
                        {
                            //recombination with K = 0.5*(F+1) --> F-K-Rule
                            do
                            {
                                // add fluctuation to random target
                                t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + 0.5 * ( f_weight + 1.0 ) *
                                                     ( gpta_old[i_r2].fa_vector[j] + gpta_old[i_r3].fa_vector[j] -
                                                       2 * gpta_old[i_r1].fa_vector[j] );
                                                       
                                j = ( j + 1 ) % gi_D;
                                k++;
                            }
                            while( ( gsl_rng_uniform( frogs_gsl_rng ) < f_cross ) && ( k < gi_D ) );
                        }
#ifdef BOUND_CONSTR
                        frogs_assigna2b( gi_D, gpta_old[i_r1].fa_vector, t_origin.fa_vector );
#endif//BOUND_CONSTR		
                    }//end if (gi_strategy ...
                    
#ifdef BOUND_CONSTR
                    for( int j = 0; j < gi_D; j++ ) //----boundary constraints via random reinitialization-------
                    {
                        //----and bounce back----------------------------------------
                        if( t_tmp.fa_vector[j] < fa_minbound[j] )
                        {
                            t_tmp.fa_vector[j] = fa_minbound[j] + gsl_rng_uniform( frogs_gsl_rng ) * ( t_origin.fa_vector[j] - fa_minbound[j] );
                        }
                        if( t_tmp.fa_vector[j] > fa_maxbound[j] )
                        {
                            t_tmp.fa_vector[j] = fa_maxbound[j] + gsl_rng_uniform( frogs_gsl_rng ) * ( t_origin.fa_vector[j] - fa_maxbound[j] );
                        }
                    }
#endif//BOUND_CONSTR
                    
                    
                    //------Trial mutation now in t_tmp-----------------
                    
                    t_tmp = frogs_evaluate( d, tmplt, prob_array, gi_D, t_tmp, &gl_nfeval, &gpta_old[0], gi_NP );
                    
                    if( i_bs_flag == TRUE )
                    {
                        gpta_new[i] = t_tmp; //save new vector, selection will come later
                    }
                    else
                    {
                        if( frogs_left_vector_wins( t_tmp, gpta_old[i] ) == TRUE )
                        {
                            gpta_new[i] = t_tmp;            // replace target with mutant
                            
                            if( frogs_left_vector_wins( t_tmp, gt_best ) == TRUE ) // Was this a new minimum?
                            {
                                // if so...
                                gt_best = t_tmp;                              // store best member so far
                            }                                               // If mutant fails the test...
                        }                                                   // go to next the configuration
                        else
                        {
                            gpta_new[i] = gpta_old[i];            // replace target with old value
                        }
                    }//if (i_bs_flag == TRUE)
                }//end for (i=0; i<gi_NP; i++)
                // End mutation loop through pop.
                
                if( i_bs_flag == TRUE )
                {
                    frogs_sort( gpta_old, 2 * gi_NP ); //sort array of parents + children
                    gt_best = gpta_old[0];
                }
                else
                {
                    gpta_swap = gpta_old;
                    gpta_old  = gpta_new;
                    gpta_new  = gpta_swap;
                }//if (i_bs_flag == TRUE)
                
                t_bestit = gt_best;
                
                
                //    }//if ....
                //======Output Part=====================================================
                
                if( gi_gen % i_refresh == 0 ) //refresh control
                {
#ifdef DO_PLOTTING
                    update_graphics( gt_best.fa_vector, gi_D, gfa_bound, gl_nfeval, gi_gen, gt_best.fa_cost[0], gi_strategy, gi_genmax );
#else
                    //fprintf(stdout,"etape intermediaire %6ld   %1.12e   %1.12e\n",gl_nfeval, gt_best.fa_cost[0], gt_best.fa_constraint[0]);
#endif//DO_PLOTTING
                    //fprintf(Fp_out,"%6ld   %12.6f\n",gl_nfeval, gt_best.fa_cost[0]);
                    
                    
                }//end if (gi_gen%10 == 0)
                
            }//end while ((gi_gen < i_genmax) && (gf_emin > MINI))
            
            //fprintf(Fp_out,"\n******** best vector ********\n", i, gt_best.fa_vector[i]);
            //fprintf(Fp_out,"\n******** best vector ********\n");
            for( int i = 0; i < gi_D; i++ )
            {
#ifndef DO_PLOTTING
                //fprintf(stdout,"lkhd0 %1.16e\n",lkhd0);
                //fprintf(stdout,"best_vector[%d]=%1.16e lkhd %1.16e\n", i, gt_best.fa_vector[i],gt_best.fa_cost[0]);
#endif//DO_PLOTTING
                //fprintf(Fp_out,"best_vector[%d]=%1.16e\n", i, gt_best.fa_vector[i]);
            }
#ifdef DO_PLOTTING
            if( gi_exit_flag == 1 )
            {
                gi_exit_flag = 0;
            }
#endif//DO_PLOTTING
            
            /*if log10(E) is defined, i.e. it is not set to FROGS_BAD_NUMBER then one compares gt_best.fa_cost[0] to lkhd0
              if gt_best.fa_cost[0]<lkhd0, the starting parameters are set by the differential evolution algorithm
              if gt_best.fa_cost[0]>lkhd0, the starting parameters are the ED parameters*/
            if( d->startpt.log10e != FROGS_BAD_NUMBER && gt_best.fa_cost[0] < lkhd0 )
            {
                d->startpt.xs = gt_best.fa_vector[0];
                d->startpt.ys = gt_best.fa_vector[1];
                d->startpt.xp = gt_best.fa_vector[2];
                d->startpt.yp = gt_best.fa_vector[3];
                d->startpt.log10e = gt_best.fa_vector[4];
                d->startpt.lambda = gt_best.fa_vector[5];
            }
            
            /*if log10(E)=FROGS_BAD_NUMBER
              the starting parameters are defined by the differential evolution algorithm*/
            if( fabs( d->startpt.log10e-FROGS_BAD_NUMBER ) < 1E-8 )
            {
                d->startpt.xs = gt_best.fa_vector[0];
                d->startpt.ys = gt_best.fa_vector[1];
                d->startpt.xp = gt_best.fa_vector[2];
                d->startpt.yp = gt_best.fa_vector[3];
                d->startpt.log10e = gt_best.fa_vector[4];
                d->startpt.lambda = gt_best.fa_vector[5];
            }
            
        }
        
        
        //================================================================
        //================================================================
        
        t_pop frogs_evaluate( struct frogs_imgtmplt_in * d,
                              struct frogs_imgtemplate * tmplt,
                              struct frogs_probability_array * prob_array,
                              int i_D, t_pop t_tmp, long * l_nfeval, t_pop * tpa_array, int i_NP )
        /**C*F****************************************************************
        **
        ** Function       :t_pop frogs_evaluate(int i_D, t_pop t_tmp, long *l_nfeval,
        **                                t_pop *tpa_array, int i_NP)
        **
        ** Author         :Rainer Storn
        **
        ** Description    :Evaluates the actual cost function (objective function)
        **                 which in this case evaluates the Chebychev fitting problem.
        **
        ** Functions      :-
        **
        ** Globals        :-
        **
        ** Parameters     :i_D         (I)    number of parameters
        **                 t_tmp       (I)    parameter vector
        **                 l_nfeval   (I/O)   counter for function evaluations
        **                 tpa_array   (I)    pointer to current population (not needed here)
        **                 i_NP        (I)    number of population members (not needed here)
        **
        ** Preconditions  :-
        **
        ** Postconditions :-
        **
        ** Return Value   :TRUE if trial vector wins, FALSE otherwise.
        **
        ***C*F*E*************************************************************/
        {
        
            /************* sphere*******************/
            /*int   i, j, k;
            float f_px, f_result=0;
            
            (*l_nfeval)++;  //increment function evaluation count
            
            for (j=0;j<i_D;j++) {
              f_px = t_tmp.fa_vector[j]; //t_tmp.fa_vector[j] corresponds to the jth position
              f_result += f_px*f_px;
              //printf("k %d j %d t_tmp.fa_vector[j] %f f_result %f\n", k, j,t_tmp.fa_vector[j], f_result);
             }
            t_tmp.fa_cost[0] = f_result;
            return(t_tmp);*/
            /************* sphere*******************/
            
            
            
            /************* Shelkel's foxholes **************/
            /*int   i, j, k;
            float f_result=0;
            
            (*l_nfeval)++;  //increment function evaluation count
            
            static float a[30][10] = {
              {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
              {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
              {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
              {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
              {8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
              {7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
              {1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448},
              {8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
              {0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
              {7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
              {0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
              {2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
              {8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
              {2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
              {4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
              {8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
              {8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
              {4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
              {2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
              {6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
              {0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
              {5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
              {3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
              {8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
              {1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
              {0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
              {0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
              {4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
              {9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
              {4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500}};
            
            static float c[30] = {
              0.806,
              0.517,
              0.100,
              0.908,
              0.965,
              0.669,
              0.524,
              0.902,
              0.531,
              0.876,
              0.462,
              0.491,
              0.463,
              0.714,
              0.352,
              0.869,
              0.813,
              0.811,
              0.828,
              0.964,
              0.789,
              0.360,
              0.369,
              0.992,
              0.332,
              0.817,
              0.632,
              0.883,
              0.608,
              0.326};
            for(k=0; k<30; k++) {
              float distance = 0.;
              for(j=0; j<i_D; j++) {
            
                //printf("k %d j %d a[k][j] %f t_tmp.fa_vector[j] %f\n", k, j, a[k][j], t_tmp.fa_vector[j]); exit(0);
                distance += pow( t_tmp.fa_vector[j]-a[k][j], 2.);
              }
              f_result -= 1./( distance+c[k] );
              //printf("f_result %f\n", f_result);
            }
            //exit(0);
            t_tmp.fa_cost[0] = f_result;
            return(t_tmp);*/
            /************* Shelkel's foxholes **************/
            
            
            
            /************ frogs *****************/
            double pd = 0;
            
            ( *l_nfeval )++; //increment function evaluation count
            
            //Here is the parameter space point
            struct frogs_reconstruction pnt;
            pnt.xs = t_tmp.fa_vector[0];
            pnt.ys = t_tmp.fa_vector[1];
            pnt.xp = t_tmp.fa_vector[2];
            pnt.yp = t_tmp.fa_vector[3];
            pnt.log10e = t_tmp.fa_vector[4];
            pnt.lambda = t_tmp.fa_vector[5];
            
            //fprintf(stdout,"in frogs_evaluate: %f %f %f %f %f %f\n",pnt.xs, pnt.ys, pnt.xp, pnt.yp, pnt.log10e, pnt.lambda);
            
            for( int tel = 0; tel < d->ntel; tel++ )
            {
                for( int pix = 0; pix < d->scope[tel].npix; pix++ )
                {
                    if( d->scope[tel].pixinuse[pix] == FROGS_OK )
                    {
                        //Here call image model and calculate expected signal
                        int pix_in_template;//FROGS_OK in image, FROGS_NOTOK in background
                        double mu = frogs_img_model( pix, tel, pnt, d, tmplt, &pix_in_template );
                        if( mu != FROGS_BAD_NUMBER )
                        {
                            if( mu > 1.e-18 && mu < FROGS_LARGE_PE_SIGNAL )
                            {
                                pd += frogs_read_prob_array_table( prob_array, d->scope[tel].q[pix], mu, d->scope[tel].ped[pix] );
                            }
                            else
                                pd += frogs_probability_density( d->scope[tel].q[pix], mu,
                                                                 d->scope[tel].ped[pix],
                                                                 d->scope[tel].exnoise[pix] );
                            //fprintf(stdout,"tel %d pix %d mu %f q %f pd %f\n", tel, pix, mu, d->scope[tel].q[pix], pd);
                        }
                    }
                }
            }
            //fprintf(stdout,"pd %f -2.0*log(pd) %f\n", pd, -2.0*log(pd));
            t_tmp.fa_cost[0] = -2.0 * log( pd ); //whole array likelihood
            return( t_tmp );
            /************ frogs *****************/
            
        }
        
        //================================================================
        //================================================================
        
        int frogs_left_vector_wins( t_pop t_trial, t_pop t_target )
        /**C*F****************************************************************
        **
        ** Function       :int frogs_left_vector_wins(t_pop t_trial, t_pop t_target)
        **
        ** Author         :Rainer Storn
        **
        ** Description    :Selection criterion of DE. Decides when the trial
        **                 vector wins over the target vector.
        **
        ** Functions      :-
        **
        ** Globals        :-
        **
        ** Parameters     :t_trial    (I)   trial vector
        **                 t_target   (I)   target vector
        **
        ** Preconditions  :-
        **
        ** Postconditions :-
        **
        ** Return Value   :TRUE if trial vector wins, FALSE otherwise.
        **
        ***C*F*E*************************************************************/
        {
            //---trial wins against target even when cost is equal.-----
            if( t_trial.fa_cost[0] <= t_target.fa_cost[0] )
            {
                return( TRUE );
            }
            else
            {
                return( FALSE );
            }
        }
        
        //================================================================
        //================================================================
        
        void frogs_permute( gsl_rng * r, int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid )
        /**C*F****************************************************************
        **
        ** Function       :void frogs_permute(int ia_urn2[], int i_urn2_depth)
        **
        ** Author         :Rainer Storn
        **
        ** Description    :Generates i_urn2_depth random indices ex [0, i_NP-1]
        **                 which are all distinct. This is done by using a
        **                 permutation algorithm called the "urn algorithm"
        **                 which goes back to C.L.Robinson.
        **
        ** Functions      :-
        **
        ** Globals        :-
        **
        ** Parameters     :ia_urn2       (O)    array containing the random indices
        **                 i_urn2_depth  (I)    number of random indices (avoided index included)
        **                 i_NP          (I)    range of indices is [0, i_NP-1]
        **                 i_avoid       (I)    is the index to avoid and is located in
        **                                      ia_urn2[0].
        **
        ** Preconditions  :# Make sure that ia_urn2[] has a length of i_urn2_depth.
        **                 # i_urn2_depth must be smaller than i_NP.
        **
        ** Postconditions :# the index to be avoided is in ia_urn2[0], so fetch the
        **                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.
        **
        ** Return Value   :-
        **
        ***C*F*E*************************************************************/
        
        {
            int  i, k, i_urn1, i_urn2;
            int  ia_urn1[MAXPOP] = {0};      //urn holding all indices
            
            k      = i_NP;
            i_urn1 = 0;
            i_urn2 = 0;
            for( i = 0; i < i_NP; i++ )
            {
                ia_urn1[i] = i;    //initialize urn1
            }
            
            i_urn1 = i_avoid;                  //get rid of the index to be avoided and place it in position 0.
            while( k >= i_NP - i_urn2_depth )  //i_urn2_depth is the amount of indices wanted (must be <= NP)
            {
                ia_urn2[i_urn2] = ia_urn1[i_urn1];      //move it into urn2
                ia_urn1[i_urn1] = ia_urn1[k - 1]; //move highest index to fill gap
                k = k - 1;                      //reduce number of accessible indices
                i_urn2 = i_urn2 + 1;            //next position in urn2
                float rnd = gsl_rng_uniform( r );
                i_urn1 = ( int )( rnd * k ); //choose a random index
            }
        }
        
        
        //================================================================
        //================================================================
        
        void  frogs_assigna2b( int i_D, float fa_a[], float fa_b[] )
        /**C*F****************************************************************
        **
        ** Function       :void  frogs_assigna2b(int i_D, float fa_a[], float fa_b[])
        **
        ** Author         :Rainer Storn
        **
        ** Description    :Assigns i_D-dimensional vector fa_a to vector f_b.
        **
        ** Functions      :-
        **
        ** Globals        :-
        **
        ** Parameters     :i_D     (I)     size of vectors
        **                 fa_a[]  (I)     source vector
        **                 fa_b[]  (I)     destination vector
        **
        ** Preconditions  :-
        **
        ** Postconditions :-
        **
        ** Return Value   :-
        **
        ***C*F*E*************************************************************/
        {
            int j;
            for( j = 0; j < i_D; j++ )
            {
                fa_b[j] = fa_a[j];
            }
        }
        
        //================================================================
        //================================================================
        
        /**C*F****************************************************************
        **
        ** Function       :void sort (t_pop ta_ary[], int i_len)
        **
        ** Author         :Rainer Storn
        **
        ** Description    :Shell-sort procedure which sorts array ta_ary[] according
        **                 to ta_ary[].fa_cost[0] in ascending order.
        **
        ** Functions      :-
        **
        ** Globals        :-
        **
        **
        ** Parameters     :ta_ary[]    (I/O)   population array
        **                 i_len        (I)    length of array to be sorteds
        **
        ** Preconditions  :-
        **
        ** Postconditions :ta_ary[] will be sorted in ascending order (according to fa_cost[0])
        **
        ** Return Value   :-
        **
        ***C*F*E*************************************************************/
        void frogs_sort( t_pop ta_ary[], int i_len )
        {
            int   done;
            int   step, bound, i, j;
            t_pop temp;
            
            step = i_len;  //array length
            while( step > 1 )
            {
                step /= 2;	//halve the step size
                do
                {
                    done   = TRUE;
                    bound  = i_len - step;
                    for( j = 0; j < bound; j++ )
                    {
                        i = j + step + 1;
                        if( ta_ary[j].fa_cost[0] > ta_ary[i - 1].fa_cost[0] )
                        {
                            temp     = ta_ary[i - 1];
                            ta_ary[i - 1] = ta_ary[j];
                            ta_ary[j]   = temp;
                            done = FALSE; //if a swap has been made we are not finished yet
                        }  // if
                    }  // for
                }
                while( done == FALSE );    // while
            } //while (step > 1)
        } //end of sort()
        
        //================================================================
        //================================================================
#ifdef CONVOLUTION
        double frogs_chertemplate_no_int( float lambda, float log10e, float b, float x,
                                          float y, struct frogs_imgtemplate * tmplt,
                                          int* intemplate )
        {
#endif
#ifndef CONVOLUTION
            double frogs_chertemplate_no_int( float lambda, float log10e, float b, float x,
                                              float y, struct frogs_imgtemplate * tmplt,
                                              int* intemplate, float pixradius )
            {
#endif
                /*This function return an evaluation of the Cherenkov ight density for
                  the given values of lambda the depth of the first interaction point,
                  log10e = log10(E/TeV), b the impact parameter to the telescope, x and
                  y the longitudinal and transverse coordinate with respect to the source
                  and the direction of development of the shower image. The evaluation is
                  obtained by linear interpolation in lambda,log10e and b in that order.
                  When the parameters fall outside the range covered by the template table,
                  the value is linearly extrapolated.
                  For x and y the closest table values are simply used.
                */
                /* Note:in the template data table, the index to parameter
                  correspondance is as follows
                  0 --- lambda = 1st interaction depth in interaction lengths
                  1 --- log10e = log10(E/1TeV)
                  2 --- b      = impact parameter to the considered telescope
                  3 --- x      = x coordinate of the pixel (see note)
                  4 --- y      = y coordinate of the pixel (see note)
                  This has nothing to do with the order in which the parameters are
                  entered in GSL for the likelihood optimization.
                  (x and y measured in degrees. x along the image major axis from
                  the source and increasing with the age of the shower. y in a
                  perpendicular direction. )*/
                
                // index for x we will not interpolate
                int ix = ( int )floor( ( x - tmplt->min[3] ) / tmplt->step[3] );
                if( ix < 0 || ix >= tmplt->nstep[3] )
                {
                    *intemplate = FROGS_NOTOK;
                    return 0.0;
                }
                
                // index for y we will not interpolate
                int iy = ( int )floor( ( fabs( y ) - tmplt->min[4] ) / tmplt->step[4] );
                if( iy < 0 || iy >= tmplt->nstep[4] )
                {
                    *intemplate = FROGS_NOTOK;
                    return 0.0;
                }
                
                //If we get here the pixel is within the area covered by the templates
                *intemplate = FROGS_OK;
                
                //index for lambda we will not interpolate
                int il = ( int )floor( ( lambda - tmplt->min[0] ) / tmplt->step[0] );
                if( il < 0 )
                {
                    il = 0;
                }
                if( il >= tmplt->nstep[0] )
                {
                    il = tmplt->nstep[0] - 1;
                }
                
                //index for energy we will not interpolate
                int iloge = ( int )floor( ( log10e-tmplt->min[1] ) / tmplt->step[1] );
                if( iloge < 0 )
                {
                    iloge = 0;
                }
                if( iloge >= tmplt->nstep[1] )
                {
                    iloge = tmplt->nstep[1] - 1;
                }
                
                //index for impact we will not interpolate
                int ib = ( int )floor( ( b - tmplt->min[2] ) / tmplt->step[2] );
                if( ib < 0 )
                {
                    ib = 0;
                }
                if( ib >= tmplt->nstep[2] )
                {
                    ib = tmplt->nstep[2] - 1;
                }
                
                //Get the model values at the vertices
                double mu000 = frogs_get_tmplt_val( il, iloge, ib, ix, iy, tmplt );
                
                return mu000;
            }
            
            //================================================================
            //================================================================
            
            float frogs_change_coordinate_system( float i_ze, float i_az, float x, float y, float z, int axis, bool bInv )
            {
            
                // transform the coordinates into another coordinate system
                // shower coord system to ground coord
                // see also void VArrayAnalyzer::transformTelescopePosition in VArrayAnalyzer.cpp
                // and void VGrIsuAnalyzer::tel_impact in VGrIsuAnalyzer.cpp
                float i_xcos = 0.;
                float i_ycos = 0.;
                
                // calculate direction cosine
                i_xcos = sin( i_ze / FROGS_DEG_PER_RAD ) * sin( ( i_az - 180. ) / FROGS_DEG_PER_RAD );
                i_ycos = sin( i_ze / FROGS_DEG_PER_RAD ) * cos( ( i_az - 180. ) / FROGS_DEG_PER_RAD );
                
                //tel_impact( i_xcos, i_ycos, x, y, z, &i_xrot, &i_yrot, &i_zrot, bInv );
                float b[3] = { 0., 0., 0. };
                float c[3] = { 0., 0., 0. };
                float matrix[3][3] = { { 0., 0., 0. },
                    { 0., 0., 0. },
                    { 0., 0., 0. }
                };
                float dl = 0.;
                float dm = 0.;
                float dn = 0.;
                /* determine the rotation matrix from setup_matrix */
                dl = i_xcos;
                dm = i_ycos;
                if( 1. - dl * dl - dm * dm < 0. )
                {
                    dn = 0.;
                }
                else
                {
                    dn = -sqrt( 1. - dl * dl - dm * dm );
                }
                
                //setup_matrix(matrix, dl,dm,dn, bInv);
                float sv = 0.;
                sv = sqrt( dl * dl + dm * dm );
                if( sv > 1.0E-09 )
                {
                    matrix[0][0] = -dm / sv;
                    matrix[0][1] = dl / sv;
                    matrix[0][2] = 0;
                    matrix[1][0] = dn * dl / sv;
                    matrix[1][1] = dn * dm / sv;
                    matrix[1][2] =  - sv;
                    matrix[2][0] = -dl;
                    matrix[2][1] = -dm;
                    matrix[2][2] = -dn;
                }
                else
                {
                    matrix[0][0] = 1;
                    matrix[0][1] = 0;
                    matrix[0][2] = 0;
                    matrix[1][0] = 0;
                    matrix[1][1] = 1;
                    matrix[1][2] = 0;
                    matrix[2][0] = 0;
                    matrix[2][1] = 0;
                    matrix[2][2] = 1;
                }
                if( bInv )
                {
                    float temp = 0.;
                    temp = matrix[0][1];
                    matrix[0][1] = matrix[1][0];
                    matrix[1][0] = temp;
                    temp = matrix[0][2];
                    matrix[0][2] = matrix[2][0];
                    matrix[2][0] = temp;
                    temp = matrix[1][2];
                    matrix[1][2] = matrix[2][1];
                    matrix[2][1] = temp;
                }
                
                for( unsigned int i = 0; i < 3; i++ )
                {
                    c[i] = 0.;
                }
                b[0] = x;
                b[1] = y;
                b[2] = z;
                if( c[0] )
                {
                    dl = 0.;
                }
                if( c[1] )
                {
                    dl = 0.;
                }
                if( c[2] )
                {
                    dl = 0.;
                }
                
                //mtxmlt(matrix, b, c);
                for( int i = 0; i < 3; i++ )
                {
                    c[i] = 0.0;
                    for( int j = 0; j < 3; j++ )
                    {
                        c[i] += matrix[i][j] * b[j];
                    }
                }
                if( c[0] )
                {
                    dl = 0.;
                }
                if( c[1] )
                {
                    dl = 0.;
                }
                if( c[2] )
                {
                    dl = 0.;
                }
                
                for( unsigned int i = 0; i < 3; i++ )
                    if( fabs( c[i] ) < 1E-5 )
                    {
                        c[i] = 0.;
                    }
                    
                if( axis == 0 )
                {
                    return c[0];
                }
                else if( axis == 1 )
                {
                    return c[1];
                }
                else if( axis == 2 )
                {
                    return c[2];
                }
                else
                {
                    return FROGS_BAD_NUMBER;
                }
            }
            
            
            
            
            void frogs_seed_gsl_rng( unsigned long seed )
            {
                //set up GSL Random Number Generation
                if( !seed )
                {
                    seed = time( 0 ) ;    //seed with system time if no seed given
                }
                gsl_rng_env_setup();
                frogs_gsl_rng = gsl_rng_alloc( gsl_rng_mt19937 ); //Select random number generator
                gsl_rng_set( frogs_gsl_rng, seed );
                fprintf( stdout, "frogs info: GSL random number generator seed %lu\n" , seed );
                
            }
            
            
            void frogs_free_gsl_rng()
            {
                if( frogs_gsl_rng )
                {
                    gsl_rng_free( frogs_gsl_rng );    //Free the memory associated with r
                }
            }
            
            
            
            
            
            
            
            
            
            
