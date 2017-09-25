//goniom_run-mpi.cpp
//Calls to external electronic structure package for fragment jobs is done here. Extrapolation of energies and forces is also handled here. 

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <limits.h>
#include <bitset> //print binary numbers
#include <boost/multiprecision/cpp_int.hpp>  //big integer lib
#include <unistd.h>
#include "keywords.h"
#include "util.h"
#include "const.h"
#include "common.h"
#include "mpi_util.h"

using namespace std; 
using namespace boost::multiprecision;


//--------------------------/
string run(string cmd,string input,string outname) { 
// this function would run "cmd", and return the output of cmd. "input" is the input file which cmd will use to run. "outname" is the name of the output file. input and output should be consistant with "cmd".
  if(ifDebug) cout<<"debug in run"<<endl;
  FILE *in; 
  char buff[512];
  string res;
  const char *cmdstr;
  cmd="goniom_gdv.tcsh "+outname +cmd; //Call to Gaussian program is handled in the goniom_gdv.tcsh, which should be located in user's PATH. 
  cmdstr=cmd.c_str(); //Command line must be structured as c-style string to pass into popen.

  ofstream outfile(outname.c_str());   //for old c++ compiler
  outfile<<input;
  outfile.flush();
  outfile.close();

 while(res=="") {
  if(!(in=popen(cmdstr,"r"))) gexit("fail to run popen"); //popen calls cmdstr string command in new process. Electronic structure package should start here. 
  while(fgets(buff,sizeof(buff),in)) res+=buff;
 }
  if(ifDebug) cout<<"finishing run"<<endl;
  if(ifDebug) cout<<res<<endl;
  return res;
}//End run


//--------------------------
string calc_chg(string str) {   //EXPIREMENTAL AND INCOMPLETE. USE ONLY WITH WATER. USE WITH CAUTION. 
// this needs to be moved to gen_sub.cpp
// need to use harsh table in the future for simplicity. 
  if(ifDebug) cout<<"debug in calc_chg"<<endl;
 string stmp;
 int n_O=0 ,n_H=0;
 int c_O=-2,c_H=1;
 stringstream is(str);
 while(is>>stmp){
     if(stmp=="O"||stmp=="8") n_O++;
     else if(stmp=="H"||stmp=="1") n_H++;
     else gexit("nyi");
     is>>stmp;is>>stmp;is>>stmp; //skip coordinates[3]
 }
 return itos(n_O*c_O+n_H*c_H);
}


//--------------------
string bg_pcharge (int n) {
 if(ifDebug) cout<<"test in bg_pcharge"<<endl;
  string res;
  cpp_int k,t; 
  t= ~mp[n-1];
  k=1;
  for(int i=0;i<natm;i++) { 
      if(k&t) res=res+" "+dtos(geom[i][0])+" "+dtos(geom[i][1])+" "+dtos(geom[i][2])+" "+ dtos(pcharge[i])+'\n';
   k=k*2;
// k<<1 has problem. << Op is not overloaded  properly for un-signed long long int. 
  }
  return res;
}
//--------------------------/
//This loop will call external electronic structure package to compute energies and forces over a set of fragments (named tmp0 for full system low level, tmp'n'-h, for nth high level fragment, and tmp'n'-l for nth low level fragment". 

void  goniom_run(){
 if(ifDebug) cout<<"debug in goniom_run"<<endl;
 string f1,f2,fout;
 string cmd,str; 
 double e,f;
 double t3, t4; //for MPI_wtime cody
 stringstream is;
 string s_getE_high=sjob[ijob[0]];
 string s_getE_low=sjob[ijob[1]];

 //This will compute the energy and for the full system at the low level of theory. mpi_id=0 is the only process to do this.
 //If the ADMP keyword is set, this will be skipped, and final extrapolated forces and energy should computed in modified electronic structure package. 
 //A modified Gaussian l121 for h37p is available to accomplish this. 
 if(mpi_id==0) {
 	if(!ifAdmp) { //check for BOMD
  		f1=full+'\n'+title+'\n'+"0 1\n"+sub_subsys[0]+'\n'+'\n'; //This is forming the command to run the full system low level. The name is tmp0. The level of theory is that of the low level. 
		//s_getE_low contains the necessary keywords to extract grep energy from output file.
  		fout="tmp"+itos(0);
	 	cmd="|"+s_getE_low;
  		if(ifForce) cmd+=";cat " +fout+".log|"+s_getF;
  		if(ifCharge) cmd=cmd+";cat " +fout+".log|"+s_getQ;
  		if(ifDipole) cmd=cmd+";cat " +fout+".log|"+s_getD;
  		str=run(cmd,f1,fout+".com");
  		is<<str;
  		is>>pot;
	
		//Read forces for full system low level calculation
	  	if(ifForce) { 
      			for(int j=0;j<natm;j++) {
             			for(int k=0;k<3;k++) { 
               				is>>f; 
               				force[j*3+k]=f;
             			}
      			}
   	  	}//end ifForce

		// read Mulliken charge from full system low level calculation.
	  	if(ifCharge) { 
     			for(int j=0;j<natm;j++) { 
       			is>>f;
       			pcharge[j]=f; 
     			}
          	}//end ifCharge
 
  		 if(ifDipole) { 
             		for(int k=0;k<3;k++) { 
               			is>>f; 
               			dipole[k]=f;
             		}
   		}
	} //end ifADMP
 } //mpi_id=0

//Not Working
 if(ifCharge) { 
 //bcast charge info
 gexit("nyi EE MPI");
 }

if(sub_subsys.size()!=eflag.size()) gexit("unknown error");

// now do all the fragments with optional EE.
//
// This line has been changed. It will no only work with np>1
 //cody
int ii;
int pid;
int pipefd[2];
int pipefor[2];
double f_buff;
double e_buff;
int natm_child;
bool tru;
bool tru_cld;
int ad=0;

if(!ifAdmp) {
mpi_np=mpi_np-1;
ad=0;
}
else{ad=-1;
}
//for(int i=mpi_id+1;i<=sub_subsys.size()-mpi_np;i+=2*mpi_np) {// Cody Old line
if(mpi_id>ad){
for(int i=mpi_id-ad;i<sub_subsys.size();i+=2*mpi_np) {
//for(int i=mpi_id+1;i<sub_subsys.size();i+=mpi_np) {
	t3 = MPI_Wtime(); //mpi timer, cody
	ii=i+mpi_np;
	if(ii<sub_subsys.size()){
	pipe(pipefd);
	pipe(pipefor);
//	pipe(pipe_eflag);	//create a pipe for communication with child process
	pid=fork();}
 	if(pid==0) i=i+mpi_np; 
//	cout<<"doing sub system number "<<i<<" on processor "<<mpi_id<<" and fork "<<pid<<endl;	




   //if(ifWater)  "0 1\n"=calc_chg(sub_subsys[i]) + "  1"+'\n';         //calc_chg is written specifically for water. comment this for other systems.
   if(!ifCharge) {
     f1=high+'\n'+title+'\n'+"0 1\n"+sub_subsys[i]+'\n'+'\n';
     f2=low+'\n'+title+'\n'+"0 1\n"+sub_subsys[i]+'\n'+'\n';
   }
   else { 
     str=bg_pcharge(i);
     f1=high+" charge\n"+'\n'+title+'\n'+"0 1\n"+sub_subsys[i]+'\n'+str+'\n'+'\n';
     f2=low +" charge\n"+'\n'+title+'\n'+"0 1\n"+sub_subsys[i]+'\n'+str+'\n'+'\n';
   }
// high 
   fout="tmp"+itos(i)+"-h";
//   cmd=" gdv<"+fout+".com>&"+fout+".log;cat "+ fout+".log|"+s_getE_high;
	cmd="|"+s_getE_high;
   if(ifForce) cmd+=";cat " +fout+".log|"+s_getF;
   str=run(cmd,f1,fout+".com");
   is.str("");is.clear();
   is<<str;
   is>>e;
   if(pid==0){
	e=e*eflag[i];
//	cout<<"processor "<<mpi_id<<"  system: "<<fout<<"  child sending: "<<e<<endl;
	write(pipefd[1],&e,sizeof(e));} //Cody: Child process will pipe energy back to parent process
   else if(pid!=0 && ii<sub_subsys.size()){ //Parent process will read energy and add it to potential
	read(pipefd[0],&e_buff,sizeof(e_buff));
	pot+=e*eflag[i]+e_buff;
   }
   else {  pot+=e*eflag[i];
//	cout<<"processor "<<mpi_id<<"  system: "<<fout<<"  parent sending: "<<e*eflag[i]<<endl;
	}
   

   if(ifForce) {
//	natm_child=0;
//	if(pid==0){
//	write(pipefor[1],&natm,sizeof(natm));
//	cout<<"process: "<<mpi_id<<"on subsystem "<<i<<"send natm_child for high "<<natm<<endl;
//	}

//	else if(pid!=0 && ii<sub_subsys.size()){
//	read(pipefor[0],&natm_child,sizeof(natm_child));
//	cout<<"process: "<<mpi_id<<"on subsystem "<<i<<"read natm_child for high "<<natm_child<<endl;
//	}
	


 
      cpp_int m=1;
      for(int j=0;j<natm;j++) {   

	tru =m&mp[i-1];
	tru_cld=0;
        if(pid==0){
        write(pipefor[1],&tru,sizeof(tru));
              }

        else if(pid!=0 && ii<sub_subsys.size()){
        read(pipefor[0],&tru_cld,sizeof(tru_cld));
               }

 	 if(tru==0&&tru_cld==1){
	for(int k=0;k<3;k++) {
	        read(pipefor[0],&f_buff,sizeof(f_buff));
                force[j*3+k]+=f_buff;
//		cout<<"process: "<<mpi_id<<" subsystem: "<<ii<<"force: "<<f_buff<<endl;
//		cout<<"process: "<<mpi_id<<"high subsystem: "<<i<<" force: "<<force[j*3+k]<<endl;
	}
	 }
	if(tru==1){ 
         for(int k=0;k<3;k++) { 
         	is>>f; 
         	if(pid==0){	//Cody: Child process will pipe force back to parent process
	       		f=f*eflag[i];
		//	 cout<<"process: "<<mpi_id<<" subsystem: "<<i<<"force: "<<f<<endl;
	       		write(pipefor[1],&f,sizeof(f));}
	       else if(pid!=0 && ii<sub_subsys.size() && tru_cld==1){
		read(pipefor[0],&f_buff,sizeof(f_buff));
		// cout<<"process: "<<mpi_id<<" subsystem: "<<i<<"force: "<<f*eflag[i]+f_buff<<endl;
		force[j*3+k]+=f*eflag[i]+f_buff;
		}
		else if(pid!=0 && tru_cld==0){
		force[j*3+k]+=f*eflag[i];
	//	cout<<"process: "<<mpi_id<<" subsystem: "<<i<<"force: "<<f*eflag[i]<<endl;
		}
//	cout<<"process: "<<mpi_id<<" high subsystem: "<<i<<" force: "<<force[j*3+k]<<endl;
//	cout<<"process: "<<mpi_id<<" high subsystem: "<<i<<" , "<<ii<<" force: "<<force[j*3+k]<<" atom: "<<j<<" "<<k<<endl;
	}
	}	
           m*=2;
      }
      



if(pid==0){
                 memset(force, 0, sizeof(force));
                }


	if(ifLink) {          // read link atom forces if there is any.
        for(int j=0;j<link_atom[i-1].size();j++) { 
             for(int k=0;k<3;k++) { 
               is>>f; 
               int ia,ib ; 
               ia=link_atom[i-1][j][0];   //i-th fragment, j-th link atom.    c array index starts from 0.
               ib=link_atom[i-1][j][1];
               link_c1=1.0-link_scale[i-1][j];
               link_c2=link_scale[i-1][j];
               force[ia*3+k]+=f*eflag[i]*link_c1;
               force[ib*3+k]+=f*eflag[i]*link_c2;
             }
        }
      }


//communicate link atoms
        for(int j=0;j<natm;j++) {
                for(int k=0;k<3;k++) {
                        if(pid==0){
                                f=force[j*3+k];
			//	cout<<" i: "<<i<<" atom: "<<j<<"	sending: "<<f;
                                write(pipefor[1],&f,sizeof(f));}
                        else if(pid!=0){
                                read(pipefor[0],&f_buff,sizeof(f_buff));
                                }
			//	cout<<endl;
                }

        }



   }
   if(ifDipole) { 
           for(int k=0;k<3;k++) { 
             is>>f; 
             dipole[k]+=f*eflag[i];
           }
   }


       t4 = MPI_Wtime();
//      cout <<"sub system " << i <<" processor: "<<mpi_id<<"fork: "<<pid<<" energy: "<<e*eflag[i]<<" force:  "<<f*eflag[i]<<" time for high: "<<t4-t3<<endl; //mpi timer, cody



t3 = MPI_Wtime(); //mpi timer, cody

// low
   fout="tmp"+itos(i)+"-l";
  // cmd=" gdv<"+fout+".com>&"+fout+".log;cat "+ fout+".log|"+s_getE_low;
   cmd="|"+s_getE_low;
   if(ifForce) cmd+=";cat " +fout+".log|"+s_getF;
   str=run(cmd,f2,fout+".com");
   is.str("");is.clear();
   is<<str;
   is>>e;
   
   if(pid==0){
        e=-e*eflag[i];
        write(pipefd[1],&e,sizeof(e));} //Cody: Child process will pipe energy back to parent process
   else if(pid!=0 && ii<sub_subsys.size()){ //Parent process will read energy and add it to potential
        read(pipefd[0],&e_buff,sizeof(e_buff));
        pot+=-e*eflag[i]+e_buff;
   }
   else{
	 pot+=-e*eflag[i];
	}



   if(ifForce) {

//	       natm_child=0;
//        if(pid==0){
//        write(pipefor[1],&natm,sizeof(natm));
//        }

//        else if(pid!=0 && ii<sub_subsys.size()){
//        read(pipefor[0],&natm_child,sizeof(natm_child));
//        }


 
//         cpp_int m=1;
//      for(int j=0;j<natm;j++) {
//         if(m&mp[i-1])  
//  
//           for(int k=0;k<3;k++) { 
//               is>>f;  
//    		if(pid==0){      //Cody: Child process will pipe force back to parent process
//    		           f=-f*eflag[i];
//    		           write(pipefor[1],&f,sizeof(f));
//			}
 //   		else if(pid!=0 && ii<sub_subsys.size() && j<natm_child){
//			read(pipefor[0],&f_buff,sizeof(f_buff));
//    		           force[j*3+k]+=-f*eflag[i]+f_buff;
//    		                                       }
//    	else force[j*3+k]+=-f*eflag[i];
//	}
//         m*=2; 
//      }
//

      cpp_int m=1;
      for(int j=0;j<natm;j++) {

        tru =m&mp[i-1];
        tru_cld=0;
        if(pid==0){
        write(pipefor[1],&tru,sizeof(tru));
              }

        else if(pid!=0 && ii<sub_subsys.size()){
        read(pipefor[0],&tru_cld,sizeof(tru_cld));
               }

         if(tru==0&&tru_cld==1){
        for(int k=0;k<3;k++) {
                read(pipefor[0],&f_buff,sizeof(f_buff));
		force[j*3+k]+=f_buff;
        }
         }
        if(tru==1){
         for(int k=0;k<3;k++) {
                is>>f;
                if(pid==0){     //Cody: Child process will pipe force back to parent process
                        f=-f*eflag[i];
                        write(pipefor[1],&f,sizeof(f));}

               else if(pid!=0 && ii<sub_subsys.size() && tru_cld==1){
                read(pipefor[0],&f_buff,sizeof(f_buff));
		force[j*3+k]+=-f*eflag[i]+f_buff;
                }
                else if(pid!=0 && tru_cld==0){
		force[j*3+k]+=-f*eflag[i];}
        }//end for k<3
        }//endif tru
           m*=2;
      }





if(pid==0){
                 memset(force, 0, sizeof(force));
                }
      if(ifLink) {
        for(int j=0;j<link_atom[i-1].size();j++) { 
             for(int k=0;k<3;k++) { 
               is>>f; 
               int ia,ib ; 
               ia=link_atom[i-1][j][0];   //i-th fragment, j-th link atom
               ib=link_atom[i-1][j][1];
               link_c1=1.0-link_scale[i-1][j];
               link_c2=link_scale[i-1][j];
               force[ia*3+k]+=-f*eflag[i]*link_c1;
               force[ib*3+k]+=-f*eflag[i]*link_c2;
             }
        }
      }//endif Link

        for(int j=0;j<natm;j++) {
                for(int k=0;k<3;k++) {
                        if(pid==0){
                                f=force[j*3+k];
                                write(pipefor[1],&f,sizeof(f));}
                        else if(pid!=0){
                                read(pipefor[0],&f_buff,sizeof(f_buff));
                                force[j*3+k]+=f_buff;
                                }
                } 
        }//end for(j<natm)
   }
 	if(ifDipole) { 
           for(int k=0;k<3;k++) { 
             is>>f; 
             dipole[k]+=-f*eflag[i];
           }
	}//endif dipole
 	
	t4 = MPI_Wtime();

 	if(ii<sub_subsys.size()){
		close(pipefd[0]);
		close(pipefd[1]);
		close(pipefor[0]);
		close(pipefor[1]);
	}//endif ii<sub_subsys

if(pid==0) _Exit(3);
}
}//end if mpi_id >0

double pot_global;
double force_global[natm*3];
//MPI_Reduce(&pot,&pot_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//pot=pot_global;
//MPI_Reduce(force,force_global,3*natm,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//pot=pot_global;
//std::copy(force_global,force_global+3*natm,force);

}


