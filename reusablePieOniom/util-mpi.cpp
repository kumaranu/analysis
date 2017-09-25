
#include <boost/multiprecision/cpp_int.hpp>  //big integer lib
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "mpi_util.h"

using namespace std;
using namespace boost::multiprecision;

//---------------------------/
//  in c++11 or later , these can be replaced by to_string(). 
string dtos(double dbn) { 
   string res;
   ostringstream strs;
   strs<<setprecision(10)<<fixed<<dbn;
   res=strs.str();
   return res;
}

string itos(int intn) { 
   string res;
   ostringstream strs;
   strs<<intn;
   res=strs.str();
   return res;
}


//---------------------------/
int vec_find(vector <cpp_int>mp,cpp_int target) {
    for(int i=0;i<mp.size();i++) if(mp[i]==target) return i;
    return -1;
}

//---------------------------/
void gexit(const char* str) {
   MPI::Finalize ( );
   cout<<"error: "<<str<<endl;
   exit(0);
}


