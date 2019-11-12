/*! \class VGrIsuAnalyzer
    \brief analysis routines from GrIsu package

    this code is copied from the GrIsu Package, see
    http://www.physcis.utah.edu/gammaray/GrISU/GrISU/Documentation/grisu.html

*/

#include "VGrIsuAnalyzer.h"

VGrIsuAnalyzer::VGrIsuAnalyzer()
{
    fGrIsuVersion = "3.0.0";
}


VGrIsuAnalyzer::~VGrIsuAnalyzer()
{
}


/*!
RETURN=    None
ARGUMENT=  prim   =Simulated primary characteristics in the original
                   ground system.
          \par xfield  the X ground locations of a telescope
      \par yfield  the y ground locations of a telescope
      \par zfield  the z ground locations of a telescope
      \par xtelrot the telescope X location in the rotated reference frame.
      \par ytelrot the telescope Y location in the rotated reference frame.
      \par ztelrot the telescope Z location in the rotated reference frame.
      \par bInv do inverse rotation from shower coordinates into ground coordinates

Function to calculate the coor. of the primary and telescope in
the rotated frame.WHAT IS THE ROTATED FRAME? DOES THE ANALYSIS WORK EVEN
IF THERE IS NO SIMULATION SPECIFIC RECORD?
*/
void VGrIsuAnalyzer::tel_impact( float xcos, float ycos, float xfield, float yfield, float zfield, float* xtelrot, float* ytelrot, float* ztelrot, bool bInv )
{
    float b[3] = { 0., 0., 0. };
    float c[3] = { 0., 0., 0. };
    float matrix[3][3] = { { 0., 0., 0. },
        { 0., 0., 0. },
        { 0., 0., 0. }
    };
    
    float dl = 0.;
    float dm = 0.;
    float dn = 0.;                               /*Direction cos of the primary in the ground frame*/
    
    /* determine the rotation matrix from setup_matrix */
    dl = xcos;
    dm = ycos;
    if( 1. - dl * dl - dm * dm < 0. )
    {
        dn = 0.;
    }
    else
    {
        dn = -sqrt( 1. - dl * dl - dm * dm );
    }
    setup_matrix( matrix, dl, dm, dn, bInv );
    for( unsigned int i = 0; i < 3; i++ )
    {
        c[i] = 0.;
    }
    
    /* determine the location of the telescope in the rotated frame */
    
    b[0] = xfield;
    b[1] = yfield;
    b[2] = zfield;
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
    
    mtxmlt( matrix, b, c );
    
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
    
    // (GM) small number check
    for( unsigned int i = 0; i < 3; i++ ) if( TMath::Abs( c[i] ) < 1.e-5 )
        {
            c[i] = 0.;
        }
        
    *xtelrot = c[0];
    *ytelrot = c[1];
    *ztelrot = c[2];
}


/*!
RETURN=Error Code
ARGUMENT=matrix: the 3x3 rotation matrix needed in rotate()
         Based on Dave Carter-Lewis' "setup_matrix" procedure from
     detector.c

This procedure sets up the rotation matrix needed to convert from a system
with the z-axis along the zenith to the correct system with the z-axis
along the incident direction of the primary particle.

cerfrnt->prim.xcos, cerfrnt->prim.ycos, cerfrnt->prim.zcos are the direction
cosines of the primary particle. matrix[3][3] holds the return values
of the 3x3 rotation matrix. */

/* C. Duke, Sept. 28, new rotation matrix
This procedure now sets up (1) a rotation about the z axis that aligns the
new x axis with the projection of the primary vector onto the ground plane,
followed by a rotation about the new y axis that aligns the z axis with the
incoming primary */
/*.....................................................................*/
void VGrIsuAnalyzer::setup_matrix( float matrix[3][3], float dl, float dm, float dn, bool bInvers )
{
    float sv = 0.;
    
    /* sv is the projection of the primary vector onto the xy plane */
    
    sv = sqrt( dl * dl + dm * dm );
    
    if( sv > 1.0E-09 )
    {
    
        /* rotation about z axis to place y axis in the plane
           created by the vertical axis and the direction of the
           incoming primary followed by a rotation about the new x
           axis (still in the horizontal plane) until the new z axis
           points in the direction of the primary.
        */
        
        matrix[0][0] = -dm / sv;
        matrix[0][1] = dl / sv;
        matrix[0][2] = 0;
        
        matrix[1][0] = dn * dl / sv ;
        matrix[1][1] = dn * dm / sv;
        matrix[1][2] =  - sv;
        
        matrix[2][0] = -dl;
        matrix[2][1] = -dm;
        matrix[2][2] = -dn;
        
    }
    /* for verital incident showers, return identity matrix */
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
    
    // invert matrix for rotations from shower coordinates into ground coordinates
    if( bInvers )
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
    
}


/*!
RETURN=None
ARGUMENT=b: vector to be rotated
         a: 3x3 rotation matrix
         c: rotated vector
Post-multiplies matrix a with vector b, returning the result in
vector c.
       | a[0][0] a[0][1] a[0][2] |   | b[0] |   | c[0] |
       | a[1][0] a[1][1] a[1][2] | x | b[1] | = | c[1] |
       | a[2][0] a[2][1] a[2][2] |   | b[2] |   | c[2] |                */
void VGrIsuAnalyzer::mtxmlt( float a[3][3], float b[3], float c[3] )
{
    for( int i = 0; i < 3; i++ )
    {
        c[i] = 0.0;
        for( int j = 0; j < 3; j++ )
        {
            c[i] += a[i][j] * b[j];
        }
    }
}


/*=====================================================*/

///:~
/**/
//:Reconst:rcs_perpendicular_fit
/***************** rcs_perpendicular_fit *********************************/
int VGrIsuAnalyzer::rcs_perpendicular_fit( vector<float> x, vector<float> y, vector<float> w, vector<float> m,
        unsigned int num_images, float* sx, float* sy, float* std )
/*
RETURN= 0 if no faults
ARGUMENT=x[10]     = x coor of point on line
         y[10]     = y coor of point on line
     w[10]     = weight of line
     m[10]     = slope of line
     num_images= number of lines
     sx        = x coor of point with minim. dist.
     sx        = y coor of point with minim. dist
     std       = rms distance from point to lines
    This procedure finds the point (sx,sy) that minimizes the square of the
perpendicular distances from the point to a set of lines.  The ith line
passes through the point (x[i],y[i]) and has slope m[i].

*/
{

    float totweight = 0.;
    float a1 = 0.;
    float a2 = 0.;
    float b1 = 0.;
    float b2 = 0.;
    float c1 = 0.;
    float c2 = 0.;
    float gamma = 0.;
    float D = 0.;
    float m2 = 0.;
    float d = 0.0;
    
    /* initialize variables */
    *sx = -999.;
    *sy = -999.;
    *std = 0.0;
    
    // check length of vectors
    
    if( x.size() == num_images && y.size() == num_images && w.size() == num_images && m.size() == num_images )
    {
        for( unsigned int i = 0; i < num_images; i++ )
        {
            totweight = totweight + w[i];
            
            m2 = m[i] * m[i];
            gamma  = 1.0 / ( 1. + m2 );
            
            /* set up constants for array  */
            D = y[i] - ( m[i] * x[i] );
            
            a1 = a1 + ( w[i] *  m2 * gamma );
            a2 = a2 + ( w[i] * ( -m[i] ) * gamma );
            b1 = a2;
            b2 = b2 + ( w[i] *  gamma );
            c1 = c1 + ( w[i] * D * m[i] * gamma );
            c2 = c2 + ( w[i] * ( -D ) * gamma );
            
        }
        /* do fit if have more than one telescope */
        if( ( num_images > 1 ) )
        {
            /* completed loop over images, now normalize weights */
            a1 = a1 / totweight;
            b1 = b1 / totweight;
            c1 = c1 / totweight;
            a2 = a2 / totweight;
            b2 = b2 / totweight;
            c2 = c2 / totweight;
            
            /*
            The source coordinates xs,ys should be solution
            of the equations system:
            a1*xs+b1*ys+c1=0.
            a2*xs+b2*ys+c2=0.
            */
            
            *sx = -( c1 / b1 - c2 / b2 ) / ( a1 / b1 - a2 / b2 );
            *sy = -( c1 / a1 - c2 / a2 ) / ( b1 / a1 - b2 / a2 );
            
            /* std is average of square of distances to the line */
            for( unsigned int i = 0; i < num_images; i++ )
            {
                d = ( float )rcs_perpendicular_dist( ( float ) * sx, ( float ) * sy,
                                                     ( float )x[i], ( float )y[i], ( float )m[i] );
                *std = *std + d * d * w[i];
            }
            *std = *std / totweight;
        }
    }
    else
    {
        cout <<  "VGrIsuAnalyzer::rcs_perpendicular_fit error in vector length" << endl;
    }
    return 0;
}


/*=============== end of rcs_perpendicular_fit ==========================*/

///:~
/**/
//:Reconst:rcs_rotate_delta
/* ==================rcs_rotate_delta===============================*/
int VGrIsuAnalyzer::rcs_rotate_delta( vector<float> xtel, vector<float> ytel, vector<float> ztel, vector<float>& xtelnew, vector<float>& ytelnew, vector<float>& ztelnew, float thetax, float thetay, int nbr_tel )
/*
RETURN=    ?
ARGUMENT=  xtel   = original x positions of telescopes
           ytel   = original y positions of telescopes
       ztel   = original z positions of telescopes
       xtelnew= new x positions of telescopes
       ytelnew= new y positions of telescopes
       ztelnew= new z positions of telescopes
       thetax = rotation angle about x-axis
       thetay = rotation angle about y-axis
       nbr_tel= number of telescopes

Function for rotating a coordinate through small angles, first about the
x-axis and then about the y-axis.  Positive rotations determined by the
right-hand rule, assuming a right-handed coordinate system.  Second order
and higher terms in thetax and thetay ignored. Thus, the sequential order
of the rotations doesn't matter

xtel, ytel, ztel are the input positions of the telescopes.
ytelnew, ytelnew, ztelnew are the rotated positions of the telescopes

I am not using the ta structure since we may well call this routine in an
interative manner where the input position is not always that given in the
ta structure*/
/*.....................................................................*/
{
    int tel = 0;
    
    for( tel = 0; tel < nbr_tel; tel++ )
    {
        xtelnew[tel] = xtel[tel] - thetax * ztel[tel];
        ytelnew[tel] = ytel[tel] - thetay * ztel[tel];
        ztelnew[tel] = ztel[tel] + thetax * xtel[tel] + thetay * ytel[tel];
    }
    
    return 0;
}


/*==============end of rcs_rotate_delta ===============================*/

///:~
/**/
//:Reconst:rcs_perpendicular_dist
/* ================rcs_perpendicular_dist===============================*/
float VGrIsuAnalyzer::rcs_perpendicular_dist( float xs, float ys, float xp, float yp, float m )
/* function to determine perpendicular distance from a point
   (xs,ys) to a line with slope m and passing through the point
   (xp,yp). Calculations in two dimensions.
*/
{
    float theta = 0.;
    float x = 0.;
    float y = 0.;
    float d = 0.;
    float dl = 0.;
    float dm = 0.;
    
    theta = atan( m );
    /* get direction cosines of the line from the slope of the line*/
    dl = cos( theta );
    dm = sin( theta );
    
    /* get x and y components of vector from (xp,yp) to (xs,ys) */
    x = xs - xp;
    y = ys - yp;
    
    /* get perpendicular distance */
    d = fabs( dl * y - dm * x );
    
    return d;
}


int VGrIsuAnalyzer::two_line_intersect( vector<float> x, vector<float> y, vector<float> w, vector<float> mx, vector<float> my, unsigned int num_images, float* sx, float* sy, float* std )
{
    *sx = 0.;
    *sy = 0.;
    *std = 0.;
    
    float wsum  = 0.;
    float a1 = 0.;
    float a2 = 0.;
    float xc = 0.;
    float yc = 0.;
    vector< float > xcore;
    vector< float > ycore;
    vector< float > weight;
    num_images = x.size();
    
    // get intersections for all possible two-telescope combinations
    for( unsigned int i = 0; i < x.size(); i++ )
    {
        for( unsigned int j = i + 1; j < x.size(); j++ )
        {
            get_intersection( x[i], mx[i], y[i], my[i], x[j], mx[j], y[j], my[j], &a1, &a2, &xc, &yc );
            // is the intersection on the right side?
            if( a1 > 0. && a2 > 0. )
            {
                xcore.push_back( xc );
                ycore.push_back( yc );
                weight.push_back( w[i]*w[j] );
            }
        }
    }
    // calculate weighted mean
    for( unsigned int i = 0; i < xcore.size(); i++ )
    {
        *sx += weight[i] * xcore[i];
        *sy += weight[i] * ycore[i];
        wsum += weight[i];
    }
    if( wsum > 0. )
    {
        *sx /= wsum;
        *sy /= wsum;
    }
    // calculate weighted variance
    for( unsigned int i = 0; i < xcore.size(); i++ )
    {
        *std += weight[i] * ( ( xcore[i] - *sx ) * ( xcore[i] - *sx ) + ( ycore[i] - *sy ) * ( ycore[i] - *sy ) );
    }
    if( wsum > 0. )
    {
        *std /= ( ( float )xcore.size() * wsum );
    }
    if( xcore.size() > 0. )
    {
        *std /= ( float )xcore.size() ;
    }
    if( *std < 1.e-5 )
    {
        *std = 0.;
    }
    
    return 0;
}


/*!
     get intersection point of two lines
*/
bool VGrIsuAnalyzer::get_intersection( float x1, float mx1, float y1, float my1, float x2, float mx2, float y2, float my2, float* a1, float* a2, float* xc, float* yc )
{
    *a1 = -1.;
    *a2 = -1.;
    *xc = 0.;
    *yc = 0.;
    
    // accept currently no zero direction vectors
    if( mx1 == 0. || my1 == 0. )
    {
        return false;
    }
    else if( mx2 == 0. || my2 == 0. )
    {
        return false;
    }
    else if( mx2 / mx1 - my2 / my1 == 0. || mx1 / mx2 - my1 / my2 == 0. )
    {
        return false;
    }
    else
    {
        *a2  = ( ( y2 - y1 ) / my1 - ( x2 - x1 ) / mx1 ) / ( mx2 / mx1 - my2 / my1 );
        *a1  = ( ( y1 - y2 ) / my2 - ( x1 - x2 ) / mx2 ) / ( mx1 / mx2 - my1 / my2 );
    }
    
    *xc = ( x1 + *a1 * mx1 );
    *yc = ( y1 + *a1 * my1 );
    
    return true;
}
