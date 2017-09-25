#include<string>

static double lscale_nc = 0.7001890359168242,lscale_cn=0.7860113421550093,lscale_cc=0.7238857938718662;
// std:: unordered_map<std::string,std::double> lscale;      // it is better to use hash-table to store scaling factor, but unordered_map is  not supported util c++11
// lscale={{"CN"},{},{}};




// scripts to grep energy
static std::string sjob[]={
                 "grep -i 'scf done'|awk '{print $5}' ",                           //0:  hf,b3lyp
                 "grep -i 'eump2'|awk '{print $6}'|sed -s 's/D/E/g'|awk '{printf \" %0.10f \", $1}' ",                              //1:  mp2
                 "grep -i 'Wavefunction amplitudes converged'|awk '{print $5}'|sed -s 's/D/E/g'|awk '{printf \"%0.10f \",$1}' ",  //2:  ccsd
                 "grep -i 'ccsd(t)= '|awk '{print $2}'|sed -s 's/D/E/g'|awk '{printf \"%0.10f \",$1}' ",                           //3:  ccsd(t)
		 "grep -i 'SCF Done' | awk '{print $5}' |awk '{printf \" %0.10f \", $1}'"						//4: PM6
                };



// scripts to grep force, charge, dipole
static std::string s_getF=" sed -n '/Forces (Hartrees/,/Cartesian/ p' | egrep -vi '[a-z]|--'|awk '{print $3,$4,$5}' ";
static std::string s_getQ=" sed -n '/Mulliken charges:/,/Sum of Mulliken/ p' | head -n -1 | tail -n +3| awk '{print $3}' ";
//string s_getD=" sed -n '/Dipole moment/, +1 p' | tail -1 | awk '{print $2 ,$4 ,$6 }' ";
static std::string s_getD=" egrep -i '^ *Dipole  *=' | sed 's/-/ -/g' | sed 's/D /E/g'|awk '{print $3,$4,$5}' ";


