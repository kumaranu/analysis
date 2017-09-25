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
 void check_job_type();
 void goniom_run();
 void output ();

 double t1, t2,t3,t4; 
 double t_starting, t_end;

//MPI initialization
 MPI::Init(argc,argv);
 mpi_np=MPI::COMM_WORLD.Get_size();
 mpi_id=MPI::COMM_WORLD.Get_rank();


//Start of timer. Uses MPI_Wtime() clock.
 if(mpi_id==0){t1 = MPI_Wtime();}

///////////////////////////////////////////////////////////
 t_starting = MPI_Wtime();
 readin(); //readin() function reads in input data from standard input. See: readin-mpi.cpp
 t_end = MPI_Wtime();

 cout << "Elapsed time in readin is " << t_end - t_starting << endl; 
///////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////
 t_starting = MPI_Wtime();
 check_job_type(); //determines theory job type (dft, hf, mp2, etc...)
 t_end = MPI_Wtime();

 cout << "Elapsed time in check_job_type is " << t_end - t_starting << endl; 
//////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////
 t_starting = MPI_Wtime();
 combination();
 t_end = MPI_Wtime();

 cout << "Elapsed time in combination is " << t_end - t_starting << endl; 
//////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////
 t_starting = MPI_Wtime();
 goniom_run(); //Electronic structure calculations on framgnets and extraploation of energy is done here. See: goniom_run-mpi.cpp
 t_end = MPI_Wtime();

 cout << "Elapsed time in goniom_run is " << t_end - t_starting << endl; 
///////////////////////////////////////////////////////////



//mpi_id=0 is choosen to handle printing output.
 if(mpi_id==0)  
 {
	output(); //prints the extrapolated energy and forces. 
	t2 = MPI_Wtime(); //second timer reading
	printf( "Overall Elapsed time is %f. For node %i\n", t2 - t1, mpi_id ); //Difference in time readings is rough estimate of program execution time.
 }//end if(mpi_id==0)
 MPI::Finalize(); 
 return 0;
 }//End main

