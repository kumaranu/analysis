#include <string>
#include <vector>

#include "util.h"
#include "common.h"
#include "keywords.h"

using namespace std;

void check_job_type(){
     if(high.find("mp2")!=-1) ijob[0]=1;
     else if(high.find("ccsd(t)")!=-1) ijob[0]=3;
     else if(high.find("ccsd")!=-1) ijob[0]=2;
     else if(high.find("pm6")!=-1) ijob[0] = 4;
     else ijob[0]=0;
     
     if(low.find("mp2")!=-1) ijob[1]=1;
     else if(low.find("ccsd(t)")!=-1) ijob[1]=3;
     else if(low.find("ccsd")!=-1) ijob[1]=2;
     else if(low.find("pm6")!=-1) ijob[0] = 4;
     else ijob[1]=0;

     if(ijob[0]<ijob[1]) gexit("low level>high level");

     if(high.find("force")!= -1) ifForce=1;
     if((low.find("force")==-1)^(high.find("force")==-1)) gexit("missing force input");
     
     if(ifForce&&ijob[0]==3) gexit("no gradient for ccsd(t)");
}


