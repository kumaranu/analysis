


#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>

#include "keywords.h"
#include "util.h"
#include "common.h"


using namespace std;


//------------------------------------/
//test auto subsystem generation.
//------------------------------------

double dist(int ia, int ib){
 double d=0.0;
 for(int j=0;j<3;j++) d+=(geom[ia][j]-geom[ib][j])*(geom[ia][j]-geom[ib][j]); 
 return sqrt(d);
}
//------------------------------------/
void gen_sub() { 
  /*
  1) generate a distance matrix.
  2) based on atmname[ ], find all the oxygens and hydrogens and save their index in nodeO[ ] and nodeH[ ].  
  3) loop over all the oxygens in nodeO[ ]  vector,  find the minimum O-O distance for each Oxygen.
  4) double loop over the oxygens in nodeO[ ], find the O-O pairs that satisfies the distance threshold ( < 110% of O-O minimum or use a hard cutoff).   
  The distances of those O-O pairs in the distance matrix are then set to  negative value to avoid overcounting in later loops.   
  5) inside the loop in 4), if  certain O-O pair satisfies the threshold, then loop over nodeH[ ] to find hydrogens that are attached to the oxygen pairs based on a O-H distance threshold (threldH). 
  The O-O pair and its hydrogens forms one fragment, and the binary vector of this fragment is then saved in  gsub[ ].  
  6)  after each fragment is found, gsub[ ] is pushed into subsys[ ] which is a 2D-vector that saves all the definition of the primary fragment.  
  subsys[ ] is the same array that saves the input fragment vectors, so if autofrag keyword is turned on, the content of subsys[ ] will be over written here. 
 */
  if(ifDebug) cout<<"debug in gen_sub"<<endl;
  vector <int>  nodeO,nodeH; 
  const double threld=1.10,threldH=1.4;
  double d,d1,d2;
  vector <bool> gsub(natm);
  double dmat[natm][natm];
  subsys.clear();
  nsub=0;
  vector <double> oo_min;

  ofstream myfile;
  myfile.open("autosub.txt");

  for(int i=0;i<natm;i++) {
    for(int j=i+1;j<natm;j++) {
         dmat[i][j]=dist(i,j); 
         dmat[j][i]=dmat[i][j]; 
    }
  }
  for(int i=0;i<natm;i++) { 
      if(atmname[i]=="O"||atmname[i]=="8")  nodeO.push_back(i);  
      else if(atmname[i]=="H"||atmname[i]=="1")  nodeH.push_back(i);
      else gexit("unrecoginized atom");
  }
  for(int io=0;io<nodeO.size();io++) { 
       double min=10000.0;
       for(int io2=0;io2<nodeO.size();io2++) { 
           if(io!=io2 && min>dmat[nodeO[io]][nodeO[io2]]) min=dmat[nodeO[io]][nodeO[io2]]; 
       }
       oo_min.push_back(min);
  }
  for(int io=0;io<nodeO.size();io++) { 
         for(int io2=0;io2<nodeO.size();io2++) {
            if(io!=io2 && dmat[nodeO[io]][nodeO[io2]]<oo_min[io]*threld && dmat[nodeO[io]][nodeO[io2]]>0) {
              dmat[nodeO[io]][nodeO[io2]]*=-1; dmat[nodeO[io2]][nodeO[io]]*=-1;
              fill(gsub.begin(),gsub.end(),0);
              gsub[nodeO[io]]=1;gsub[nodeO[io2]]=1;
            myfile<<"sub test O "<<nodeO[io]+1<<"  "<<nodeO[io2]+1<<"  " <<-dmat[nodeO[io]][nodeO[io2]]<<endl;
            myfile<<"sub test H " ;
 //now get the hydrogens.
               for(int k=0;k<nodeH.size();k++) {
                     d1=dmat[nodeO[io]][nodeH[k]];
                     d2=dmat[nodeO[io2]][nodeH[k]];
                     if(d1<threldH||d2<threldH) { 
                       myfile<<" "<<nodeH[k]+1;
                       gsub[nodeH[k]]=1;
                     }
               }
                myfile<<endl<<endl;
               nsub++;
               subsys.push_back(gsub);
            }
         }
  }
  
        myfile<<"nsub "<<nsub<<endl<<endl;
        int a=subsys.size();
        int b=subsys[0].size();
        for(int i=0;i<b;i++)  {
          for(int j=0;j<a;j++) { 
             myfile<<"  "<<subsys[j][i]<<" ";
          }
          myfile<<endl;
        }

}
//----------------------------/
//----------------------------/



