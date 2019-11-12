//! VGrIsuAnalyzer   analysis function from GrIsu package + new functions

#ifndef VGRISUANALYZER_H
#define VGRISUANALYZER_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TMath.h"

using namespace std;

class VGrIsuAnalyzer
{
    private:
        string fGrIsuVersion;
        
        void mtxmlt( float a[3][3], float b[3], float c[3] );
        void setup_matrix( float matrix[3][3], float dl, float dm, float dn, bool bInv );
        
        bool get_intersection( float, float, float, float, float, float, float, float, float*, float*, float*, float* );
        
    protected:
        int    two_line_intersect( vector<float> x, vector<float> y, vector<float> w, vector<float> mx, vector<float> my, unsigned int num_images, float* sx, float* sy, float* std );
        float rcs_perpendicular_dist( float xs, float ys, float xp, float yp, float m );
        int    rcs_perpendicular_fit( vector<float> x, vector<float> y, vector<float> w, vector<float> m, unsigned int num_images, float* sx, float* sy, float* std );
        int    rcs_rotate_delta( vector<float> xtel, vector<float> ytel, vector<float> ztel, vector<float>& xtelnew, vector<float>& ytelnew, vector<float>& ztelnew, float thetax, float thetay, int nbr_tel );
        
    public:
        VGrIsuAnalyzer();
        ~VGrIsuAnalyzer();
        string getGrIsuVersion()
        {
            return fGrIsuVersion;
        }
        void   tel_impact( float xcos, float ycos, float xfield, float yfield, float zfield, float* xtelrot, float* ytelrot, float* ztelrot, bool bInv );
        
};
#endif
