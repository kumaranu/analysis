#include<string>
#include<vector>
#include <boost/multiprecision/cpp_int.hpp> 

std::string high,low,title,mulchg,full;
 std::vector<std::string> atmname;
 std::vector< std::vector<double> > geom;
// std::vector< std::vector<double> > force;
 double* force;
 std::vector< std::vector<bool> > subsys;
 std::vector <double> pcharge;
 std::vector<int> eflag; //+1 or -1
 std::vector<std::string> sub_subsys; // all systems to be calculated
 std::vector<boost::multiprecision::cpp_int> mp;

 int nsub,natm;
 double pot;
 double dipole[3];
 std::vector< std::vector<int> > connect; //saves connectivity
 int ijob[2];    // 0 -> high, 1-> low


//---------------/
 std::vector< std::vector <std::vector <int> > > link_atom;   // saves link atom info for each fragment.
 std::vector< std::vector <double > > link_scale;   // scaling factor for link atom.
// link atom , scaling factor
 double link_c1 ,link_c2;

