/*  check if a root can be opened or it
 *  is the recovery bit is set
 *
 *  return 0 if everything is good, otherwise 1
 *
 */
#include <TFile>

int checkCorrupted(TString F) 
{
    TFile *f = TFile::Open(F); 
    if ((!f) || f->IsZombie() || f->TestBit(TFile::kRecovered)) 
    { 
        //cout << "There is a problem with the file: $F\n"; 
        exit(1); 
    }
    exit(0);
}
