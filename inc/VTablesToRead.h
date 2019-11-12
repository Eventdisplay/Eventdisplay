//! VTablesToRead service class for lookup tables

#ifndef VTABLEREAD_H
#define VTABLEREAD_H

#include "TH2F.h"
#include "TH2D.h"

#include <iostream>
#include <map>
#include <vector>

using namespace std;

class VTablesToRead
{
    public:
    
        unsigned int    fNLT_types;
        unsigned int    fNTel;
        map< unsigned int, vector< TH2F* > > hMedian;
        map< unsigned int, vector< TH2F* > > hSigma;
        
        map< unsigned int, double > value;
        map< unsigned int, double > value_Chi2;
        map< unsigned int, double > value_dE;
        map< unsigned int, double* > value_T;
        map< unsigned int, double* > value_T_sigma;
        
        VTablesToRead( unsigned int n_tabletypes, int nTel );
        ~VTablesToRead() {}
        void reset();
};
#endif
