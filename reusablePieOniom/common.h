#include<string>
#include<vector>
#include <boost/multiprecision/cpp_int.hpp> 

extern std::string high,low,title,mulchg, full;
extern std::vector<std::string> atmname;
extern std::vector< std::vector<double> > geom;
//extern std::vector< std::vector<double> > force;
extern double* force;
extern std::vector< std::vector<bool> > subsys;
extern std::vector <double> pcharge;
extern std::vector<int> eflag; //+1 or -1
extern std::vector<std::string> sub_subsys; // all systems to be calculated
extern std::vector<boost::multiprecision::cpp_int> mp;

extern int nsub,natm;
extern double pot;
extern double dipole[3];
extern std::vector< std::vector<int> > connect; //saves connectivity
extern int ijob[2];    // 0 -> high, 1-> low


//---------------------/
extern std::vector< std::vector <std::vector <int> > > link_atom;   // saves link atom info for each fragment.
extern std::vector< std::vector <double > > link_scale;   // scaling factor for link atom.
// link atom , scaling factor
extern double link_c1 ,link_c2;

