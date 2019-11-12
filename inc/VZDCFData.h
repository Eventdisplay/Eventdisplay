//! VZDCFData data class for ZDCF plotting

#ifndef VZDCFData_H
#define VZDCFData_H

#include <iostream>

using namespace std;

class VZDCFData
{
    private:
    
    public:
    
        double tau;
        double sigma_tau_neg;
        double sigma_tau_pos;
        double dcf;
        double dcf_error_low;
        double dcf_error_up;
        int    nbins;
        
        VZDCFData();
        ~VZDCFData() {}
        
        void print();
        
        //   ClassDef( VZDCFData, 1 );
};

#endif
