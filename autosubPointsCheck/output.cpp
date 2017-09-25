


#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <bitset>
#include <fstream>

#include "common.h"
#include "keywords.h"

using namespace std;

void output () {

cout<<fixed<<setprecision(10)<<pot<<endl;

if(ifForce) 
for(int i=0;i<natm;i++)  {
  for(int j=0;j<3;j++) cout<<fixed<<setprecision(10)<<scientific<<force[i*3+j]<<"  ";
  cout<<endl;
}

if(ifDipole) 
  for(int j=0;j<3;j++) cout<<fixed<<setprecision(10)<<scientific<<dipole[j]<<"  ";
  cout<<endl;

}


