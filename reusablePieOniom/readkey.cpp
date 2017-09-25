#include<string>
#include <sstream>

#include "util.h"
#include "keywords.h"

using namespace std;


//------------------------------------------/
// read keywords

//dipole | for dipole calculation. Turned off by default. 
//charge | use background point charges. 
//autosub| automatically generate subsystems.
//debug | print additional debug informations
//! before "key:" is comment

void readkey(string str) { 
  string strc; 
  stringstream is(str);
  is>>strc; //remove "key:"
  while(is>>strc) {
    if(strc=="dipole")  ifDipole=1;
    else if(strc=="charge"||strc=="ee") ifCharge=1;
    else if(strc=="autosub") ifAutosub=1;
    else if(strc=="autofrag") ifAutosub=1;
    else if(strc=="debug") ifDebug=1;
    else if(strc=="link") ifLink=1;
    else if(strc=="water") ifWater=1;
    else if(strc=="admp") ifAdmp=1;
//   else if(strc=="maxsuboverlap") max_sub_overlap=1;
    else  gexit("wrong keyword");
  }
}
