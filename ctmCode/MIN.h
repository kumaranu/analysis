//function declarations
#pragma once
#include<iostream>
#include<fstream>
#include<vector>

using namespace std;

#ifndef MIN_H
#define MIN_H

class MIN {
public:
	//Default Constructor
	MIN();

	//Overload Constructor
	MIN(vector<vector<double> >);

	//Destructor
	~MIN();

	//Mutator function
	void setPoints(vector<vector<double> >);
	// setCommon - sets the input file which has 3 columns and we are checking the common in first two columns
	// @param double[19][3] - array of the input file which has the size [19][3]

	void shiftToMinimum(vector<vector<double> > &);
	// shiftToMinimum - This function shifts the input to zeros and puts it in an array
	// @param double[19][3] - This is the array in which the output after shifting to zero will be stored
	// Check it for MISTAKES

private:
	//Member Variables
	vector<double> newMinimum;
	vector<vector<double> > newCommon;
	vector<vector<double> > newShiftedToZero;
};

#endif // !MIN_H
