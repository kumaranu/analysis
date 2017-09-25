
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "readkey.h"
#include "common.h"
#include "keywords.h"
#include "gen_sub.h"

#include "mpi_util.h"

using namespace std;



//-------------------------------------
void readin(){
// read input file
  string str,str2;
  int isub;

  stringstream input_ss; 
  string mpi_input_str="";
  string st;
  int ncast;

if(mpi_id==0){
  while(std::getline(std::cin,str)) {mpi_input_str+=str+'\n';} 
  ncast=mpi_input_str.size();
//  cout<<mpi_input_str<<endl;
}
 MPI_Bcast(&ncast,1, MPI_INT, 0, MPI_COMM_WORLD); 
// MPI_Barrier(MPI_COMM_WORLD);
 char cp[ncast];   //temp char array for mpi_bcast. Hope there is a cleaner way. 
 if(mpi_id==0) strcpy(cp,mpi_input_str.c_str());
 MPI_Bcast(cp, ncast, MPI_CHAR, 0, MPI_COMM_WORLD); 
 if(mpi_id!=0) mpi_input_str.assign(cp,cp+ncast);

 input_ss<<mpi_input_str; 



 
 do {getline(input_ss,str);} while(str=="");
 stringstream is_tmp(str); 
 is_tmp>>str2; 
 if (str2[0]=='!') do {getline(input_ss,str);} while(str=="");   //check comments
 else if(str2=="key:") {                                    //key: start reading keyworkds. 
       readkey(str);                                    
       do {getline(input_ss,str);} while(str=="");                 //  read multiple lines. 
 }
//read full system low level
  do {getline(input_ss,str);} while(str==""); //Cody
  full ="%mem=8GB\n";
  while(str!="") {full+=str+'\n';getline(input_ss,str);}

// read high level
  do {getline(input_ss,str);} while(str==""); //uncommented cody
  while(str!="") {high+=str+'\n'; getline(input_ss,str); }

// read low level
  do {getline(input_ss,str);} while(str=="");
  while(str!="") {low+=str+'\n';getline(input_ss,str);}
  

// read title
  do {getline(input_ss,str);} while(str=="");
  title=str+'\n';

// read multiplicity and charge
  do {getline(input_ss,str);} while(str=="");
  mulchg=str+'\n';

//count the number of subsystems
nsub=-4; //first 4 cols are atomname and coordinates.
do {getline(input_ss,str);} while(str=="");
stringstream is(str);
while(is>>str2) nsub++;
subsys.resize(nsub);

//read in geometries and subsystems
natm=0;
while(str!="") {
   natm++;
   vector <double> tmp_geom(3);
   stringstream is(str);
   is>>str2;
   atmname.push_back(str2);
   is>>tmp_geom[0]>>tmp_geom[1]>>tmp_geom[2];
   geom.push_back(tmp_geom);
   for(int i=0;i<nsub;i++) {is>>isub;subsys[i].push_back(isub);}
   getline(input_ss,str);
}

//force.resize(natm,vector<double>(3));
force=new double[natm*3];
pcharge.resize(natm);

if(nsub==0) ifAutosub=1;
if(ifAutosub)  gen_sub();
if(max_sub_overlap==0) max_sub_overlap=nsub;
max_sub_overlap=min(nsub,max_sub_overlap);

if(ifLink) { 
  connect.resize(natm); 
  int iatm1,iatm2;
  do {getline(input_ss,str);} while(str=="");
  for(int i=0;i<natm;i++) { 
    stringstream is(str);
    is >> iatm1;  
    iatm1--;
    while(is>>iatm2 && is>>str2) {   //str2 is bond length, not used here.
      iatm2--;
      connect[iatm1].push_back(iatm2);
      connect[iatm2].push_back(iatm1);
    }
   getline(input_ss,str);
  }

// prints the connectivity.
//  for(int i=0;i<natm;i++) {
//        cout<<"atm "<<i+1;
//    for(int j=0;j<connect[i].size();j++) {
//       cout<<" "<<connect[i][j]+1;
//    }
//       cout<<endl;
//  }
}

}
