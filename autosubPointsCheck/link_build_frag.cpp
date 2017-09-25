
#include <vector>
#include <string>
#include <bitset> //print binary numbers

#include "keywords.h"
#include "util.h"
#include "const.h"
#include "common.h"
#include <boost/multiprecision/cpp_int.hpp>  //big integer lib
#include <math.h>

using namespace std; 
using namespace boost::multiprecision;






//---------------------------/
cpp_int NextN(cpp_int N)  
{  
    cpp_int x = N&(-N);        
    cpp_int t = N+x;  
    cpp_int ans = t | ((N^t)/x)>>2;  
    return ans;  
}  

//--------------------------/

cpp_int build_sub(cpp_int C,string &sgeom)  
{    
//cout<<"debug in build_sub: sizeof isys "<<sizeof(cpp_int)*8<<"bits" <<endl;
    cpp_int k = 1;  
    bool tmp[natm];
    int i;
    cpp_int bit = 1;
    string slink; //geom for the link atoms.
    cpp_int isys,j;
    sgeom=""; //geom for ssubsystems
    isys =0; 
    
    for(long long i=0;i<natm;i++) tmp[i]=1;
    //cout << "temp created" << endl;
    long long m = 0;  
    while((k=bit<<m)<=C)
    { 
      //cout << " k = " << k << " m = " << m << "k size " << sizeof(k) << endl; 
      if((C&k)!=0)  
      {  
        for(long long i=0;i<natm;i++) {
	//cout << "subsys = " << subsys[m][i]<< "&& temp[i] = " << tmp[i]  << endl;
	tmp[i]=subsys[m][i]&&tmp[i];
	//cout << "temp[i] = " << tmp[i] << endl; 
      }
      }
      m++;  
    } 
   //cout << "subsys created" << endl; 
   long long iatm=0;  //tmp for auto sub
    for( i=0, j=1;i<natm;i++,j*=2) {
        if(tmp[i]) {
              isys+=j;
              iatm+=1;
	     // cout << "isys = " << isys << " iatm = " << iatm << " i = " << i << " j = " << j << endl;
              sgeom+=atmname[i]+"  "+dtos(geom[i][0])+"  "+dtos(geom[i][1])+"  "+dtos(geom[i][2])+'\n';
        }
    }
    //cout << "start link atoms creation" << endl;
    if(ifLink) {
    exit;
      vector <vector<int> > tmpi1;
      vector<int>  tmpi2;
      vector<double>  tmpd;
      tmpi1.clear(); 
      tmpd.clear();
      string LAtom = "H  ";
      for(i=0;i<natm;i++) {
          if(tmp[i]) {
              for(long long j=0;j<connect[i].size();j++) {  
                  long long  ia=connect[i][j];
                  if(!tmp[ia]) {
		    //cout << "Link place found" << endl;  
                    tmpi2.clear();
                    if((atmname[i]=="C" && atmname[ia]=="C") ||(atmname[i]=="6" && atmname[ia]=="6")) link_c2=lscale_cc;
                    if((atmname[i]=="N" && atmname[ia]=="C") ||(atmname[i]=="7" && atmname[ia]=="6")) link_c2=lscale_nc;
                    if((atmname[i]=="C" && atmname[ia]=="N") ||(atmname[i]=="6" && atmname[ia]=="7")) link_c2=lscale_cn;
                    tmpd.push_back(link_c2);
                    link_c1=1.0-link_c2;
                    slink+=  LAtom+   "  "+dtos(geom[i][0]*link_c1+geom[ia][0]*link_c2);                  //x
                    slink+=           "  "+dtos(geom[i][1]*link_c1+geom[ia][1]*link_c2);                  //y
                    slink+=           "  "+dtos(geom[i][2]*link_c1+geom[ia][2]*link_c2)+'\n';             //z
                    tmpi2.push_back(i); tmpi2.push_back(ia);  
                    tmpi1.push_back(tmpi2);
		    //cout << "link atom placed at "<< slink << endl;
                  }
              }
          }
      }
      sgeom+=slink;
      link_atom.push_back(tmpi1);
      link_scale.push_back(tmpd);
    }


    if(iatm==1) {isys=0;sgeom="";} //if only one atom overlaped. no calculation. THIS HAPPENS IN H-WATER CLUSTER
    //cout << "end build" << endl; 
    return isys;
}



//--------------------------/

void combination()
{try{  
  cpp_int isys; 
  cpp_int bit = 1;
// the whole system first
    string sgeom="";
    for(long long i=0;i<natm;i++) {
        sgeom+=atmname[i]+"  "+dtos(geom[i][0])+"  "+dtos(geom[i][1])+"  "+dtos(geom[i][2])+'\n';
    }
    sub_subsys.push_back(sgeom); eflag.push_back(1);

// sub systems
  for(long long i=1;i<=max_sub_overlap;i++) {
    cpp_int C = (bit<<i)-bit;    //C change to big int later
    while(C<=((bit<<nsub)-(bit<<(nsub-i))))  
    { //cout<<"test in combination: C "<<(bitset<32>) C<<endl;
      // cout << "C = " << C << endl;  
       isys=build_sub(C,sgeom);
       //cout << "isys = " << isys << endl;  
       if(sgeom!=""){
          long long ix;
          i%2==1?ix=1:ix=-1;
          long long ifind=vec_find(mp,isys);
          if(ifind==-1){   //nothing is found.  no duplication.
             sub_subsys.push_back(sgeom);
             eflag.push_back(ix);
             mp.push_back(isys);
          }
          else {
           eflag[ifind+1]+=ix;   // 1st element is reserved for the whole system.
          }
       }
        //cout <<" enter NextN" << endl;
        C = NextN(C);  
    }       
  }
      
}catch(const std::overflow_error& e){
cout << "C type is too small" << endl;
}
catch(const std::runtime_error& e) {
  cout << "C type is too large" << endl;
 } catch(const std::exception& e) {
   cout << "logic error" << endl;
 } catch(...) {
cout << "Misc error" << endl;
}
}
