//! VImageParameter class to store run and image parameter information for a single telescope
#ifndef VImageParameter_H
#define VImageParameter_H

#include "TGraph.h"
#include "TTree.h"

#include <iostream>

using namespace std;

class VImageParameter
{
    private:
        TTree* tpars;
        bool fMC;
        unsigned int fShortTree;                 // if set, only a subset of the parameters are written to the tree of parameters
        
    public:
        // global parameters
        int fTelID;                               //!< telescope ID
        int runNumber;
        int eventNumber;
        unsigned short int eventType;
        int MJD;
        double time;
        unsigned int nanosec;
        int nsamples;
        
        unsigned int eventStatus;
        
        // telescope parameters
        int fsumfirst;                            //!< parameter for window summation start
        int fsumwindow;                           //!< parameter for window summation
        int fsumwindow_2;                         //!< parameter for window summation for energy reconstruction
        short int fLocalTrigger;                  //!< 1 of this telescope had a local trigger
        unsigned short int fTrig_type;            //  trigger type (e.g. in CTA prod2)
        
        // telescope position in shower parameter
        double Tel_x_SC;                          //!< telescope position in shower coordinates
        double Tel_y_SC;                          //!< telescope position in shower coordinates
        double Tel_z_SC;                          //!< telescope position in shower coordinates
        
        // image parameters
        float fmeanPed_Image;                     //!< mean pedestal in this image
        float fmeanPedvar_Image;                  //!< mean pedestal variation in this image
        float f_d;                                //!< see Fegan 1997, Table 6
        float f_s;                                //!< see Fegan 1997, Table 6
        float f_sdevxy;                           //!< see Fegan 1997, Table 6
        float cen_x;                              //!< Center position of shower
        float cen_y;                              //!< Center position of shower
        float cen_x_trig;                         //!< Image centroid calculated from binary (trigger-level) information
        float cen_y_trig;                         //!< Image centroid calculated from binary (trigger-level) information
        float cen_x2_trig;                        //!< Image centroid^2 calculated from binary (trigger-level) information
        float cen_y2_trig;                        //!< Image centroid^2 calculated from binary (trigger-level) information
        float sigmaX;                             //!< standard deviation x-direction
        float sigmaY;                             //!< standard deviation y-direction
        float length;                             //!< Length of ellipse
        float width;                              //!< width of ellipse
        float size;                               //!< total signal
        float size2;                              //!< total signal (second summation window)
        float sizeLL;                             //!< total signal; LL method
        float size2LL;                             //!< total signal (second summation window); LL method
        float dist;                               //!< distance to centroid
        float azwidth;                            //!< used by certain analysis techniques
        float alpha;                              //!< Alpha angle (-Pi/2..Pi/2)
        float los;                                //!< length over size
        float miss;                               //!< miss parameter
        float phi;                                //!< Angle that centroid makes with x-axis
        float cosphi;                             //!< cos of Angle betw major axis and x-axis
        float sinphi;                             //!< sin of Angle betw major axis and x-axis
        float asymmetry;                          //!< measure of shower skew
        
        float loss;                               //!< fraction of image size in outer pixel
        float lossAndDead;                        //!< fraction if image size in outer pixels and at the edge of dead pixels
        float fracLow;                            //!< fraction of image size in low gain pixel
        float fui;                                //!< fraction of image/border pixel under estimate image ellipse
        
        // muon parameters
        float muonX0;                             //!< center of muon ring X-coord
        float muonY0;                             //!< center of muon ring Y-coord
        float muonXC;                             //!< center of muon ring X-coord of centroid
        float muonYC;    						  //!< center of muon ring Y-coord of centroid
        float muonRadius;                         //!< radius of muon ring
        float muonRSigma;                         //!< std. dev. of radius of muon ring
        float muonSize;                           //!< total amount of light in muon ring
        int   muonValid;                          //!< 0/1 depending on wether it satisfies criteria
        float muonIPCorrectedSize;                //!< total amount of light in muon ring corrected by impact paameter
        
        // Hough transform muon parameters
        double 	houghAP;			  //!< AP parameter
        double 	houghTD; 			  //!< TD parameter
        int 	houghNpix; 			  //!< Number of hit pixels
        double 	houghCN; 			  //!< C/N parameter
        double 	houghContained; 		  //!< Distance from the center of the ring to the center of the camera plus the ring radius in mm. Before filling tree, this is divided by the pixel diameter
        int 	houghMuonValid;			  //!< 0/1 depending on whether it satisfies the Hough transform muon ID criteria
        
        // signal section
        unsigned short int ntubes;                //!< number of tubes in the picture
        unsigned short int trig_tubes;            //!< number of tubes selected in the trigger algorithm
        unsigned short int ntubesBrightNoImage;   //!< number of tubes which are bright but not in the image
        unsigned short int ntrig;                 //!< number of tubes triggering their CFD
        unsigned short int ntrig_per_patch;       //!< number of tubes triggering their CFD per 19 pixel patch
        unsigned short int nlowgain;              //!< number of tubes triggering Hi/Lo gain switch
        unsigned short int nzerosuppressed;       //!< number of zero suppressed channels
        unsigned short int nsat;                  //!< number of saturated tubes
        unsigned short int bad;                   //!< Good event?
        unsigned short int badLow;                //!< Good event?
        float max[3];                             //!< The three largest adc values
        unsigned short int  index_of_max[3];      //!< The tube indices of the max[3] values
        float frac[3];                            //!< fraction of maximum digital counts
        
        // timing parameters
        float tgrad_x;                            //!< Timing gradient in X-direction
        float tint_x;                             //!< Timing intercept for X-direction
        float tgrad_dx;                           //!< error in timing gradient in X-direction
        float tint_dx;                            //!< error in timing intercept for X-direction
        float tchisq_x;                           //!< Chisquare of timing fit in X
        float tmin;                               //!< time minimum for image/border pixels
        float tmax;                               //!< time maximum for image/border pixels
        float tmean;                              //!< mean time of image/border pixels
        
        // MC parameters
        unsigned short int MCprimary;
        float MCenergy;
        float MCxcore;
        float MCycore;
        float MCxcos;
        float MCycos;
        float MCLocalTriggerTime;                 //!< local trigger time (preli! no MC info, should be somewhere else)
        float MCLocalDelayedTriggerTime;          //!  delayed local trigger time (preli! no MC info, should be somewhere else)
        float MCTel_Xoff;                         //!< source offset in MC in deg (grisudet telescope coordinate system)
        float MCTel_Yoff;                         //!< source offset in MC in deg (grisudet telescope coordinate system)
        
        // log likelihood fit parameters and erros
        int ntfit;                                //!< number of restored tubes
        float Fitmin;
        float Fitedm;
        int Fitstat;                              //!< status of covariance matrix (0=not calculated,1=diagonal approx., not accurate,
        //!< 2=full matrix, but forced to be positive definite, 3=full matrix, accurate
        float ntRec;                              //!< number of fit-recovered dead channels
        float dcen_x;
        float dcen_y;
        float dsigmaX;
        float dsigmaY;
        float dlength;
        float dwidth;
        float ddist;
        float dmiss;
        float dphi;
        float dalpha;
        float dazwidth;
        float rho;
        float drho;
        float signal;
        float dsignal;                            //!< error in signal (normalisation parameter for fit)
        
        vector< float > fImageBorderPixelPosition_x;              //! list of image+border pixel
        vector< float > fImageBorderPixelPosition_y;              //! list of image+border pixel
        
        VImageParameter( unsigned int iShortTree = 0 );
        ~VImageParameter();
        void fill();
        TTree* getTree()
        {
            return tpars;
        }
        bool hasImage();                          //!< succesfull image reconstruction
        void initTree( string, string, bool, bool, bool, bool );
        bool isMC()
        {
            return fMC;
        }
        void printParameters();
        void reset( unsigned int resetLevel = 0 );
        void setImageBorderPixelPosition( vector< float > iImageBorderPixelPosition_x, vector< float > iImageBorderPixelPosition_y );
        void setMC()
        {
            fMC = true;
        }
};
#endif
