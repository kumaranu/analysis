#include<iostream>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

#ifndef AutosubShift
#define AutosubShift

class AUTOSUBSHIFT
{
public:
	//Default Constructor
	AUTOSUBSHIFT();

	//Overload Constructor
	AUTOSUBSHIFT(vector<vector<double> >, vector<vector<unsigned int> >);

	//Functions
	void setPoints(vector<vector<double> >);

	void setAutosubPoints(vector<vector<unsigned int> >);

	void shiftAtAutosub(vector<vector<double> >&);

	//Destructor
	~AUTOSUBSHIFT();

	vector<vector<unsigned int> > newAutosubPoints;

private:
	//Member variables
	vector<vector<double> > newShiftedAtMinimum;
	vector<vector<double> > autosubShifts;
};


#endif // !AutosubShift
#pragma once
