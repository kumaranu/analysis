// generalized oniom code
// the following new features are avaiable:
// link atom (C-C, C-N, N-C only)
// single-layer admp
//
//Original Author: Junji Li
//Modified by: Cody Haycraf
//Last modified: Fri Dec 11 15:22:12 EST 2015 


//NOTES for main()
//
//
//Extraction of forces and energies are done using external shell commands. The specific form of those commands is determined in check_job_type();




#include <boost/multiprecision/cpp_int.hpp>  
#include "keywords.h"
#include "util.h"
#include "const.h"
#include "common.h"
#include "build_frag.h"
#include "mpi_util.h"

using namespace std; 
using namespace boost::multiprecision;

int main(int argc,char*argv[]){
 void readin ();
// void check_job_type();
// void goniom_run();
// void output ();

 double t1, t2;
 double t3, t4; 

//MPI initialization
// MPI::Init(argc,argv);
// mpi_np=MPI::COMM_WORLD.Get_size();
// mpi_id=MPI::COMM_WORLD.Get_rank();


//Start of timer. Uses MPI_Wtime() clock.
// if(mpi_id==0){t1 = MPI_Wtime();}

 readin(); //readin() function reads in input data from standard input. See: readin-mpi.cpp
// check_job_type(); //determines theory job type (dft, hf, mp2, etc...)
// combination();
// goniom_run(); //Electronic structure calculations on framgnets and extraploation of energy is done here. See: goniom_run-mpi.cpp

//mpi_id=0 is choosen to handle printing output.
// if(mpi_id==0)  
// {
//	output(); //prints the extrapolated energy and forces. 
//	t2 = MPI_Wtime(); //second timer reading
//	printf( "Elapsed time is %f\n", t2 - t1 ); //Difference in time readings is rough estimate of program execution time.
// }//end if(mpi_id==0)
// MPI::Finalize(); 
 return 0;
 }//End main

