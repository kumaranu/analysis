#pragma once
#ifndef OniomXs
#define OniomXs

#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<math.h>
#include<complex>
#include<cmath>

using namespace std;
using namespace Eigen;

class Oniom_Xs{
public:
	//Default Constructor
	Oniom_Xs();

	//Overloaded Constructor
	Oniom_Xs(vector<vector<double> >&,vector<vector<double> >&,vector<vector<unsigned int> >,vector<double>);
		//Oniom_Xs : It constructs the object for Oniom_Xs class with shiftedToAutosub array and epsilon
		//vector<vector<double> > shiftedToAutosub array
		//vector<vector<double> > epsilon values array, I took constant and same values epsilon. It can be changes later in any possible way


	void setShiftedAtAutosub(vector<vector<double> >&);
	void setEpsilon(vector<vector<double> >&);
	void setAutosubPoints(vector<vector<unsigned int> >&);
	void setDistance(vector<double>&);

	void getOniomXsEnergy(vector<vector<double> >&,const int,const int);
	void getJustShifted(vector<vector<double> >&, const int, const int);
	//Destructor
	~Oniom_Xs();

private:
	//Member variables
	vector<vector<double> > newShiftedAtAutosub;
	vector<vector<double> > newEpsilon;
	vector<vector<unsigned int> > newAutosubPoints;
	vector<vector<complex<double> > > newOniomXsEnergy;
	vector<double> newDistance;
};
#endif // !OniomXs
