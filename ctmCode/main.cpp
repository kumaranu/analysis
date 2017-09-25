#include<iostream>
#include<fstream>
#include<string>
#include<vector>

#include"MIN.h"
#include"AutosubShift.h"
#include"OniomXs.h"

using namespace std;

int main(int argc,char* argv[]) {
	const int nRows = 199, nTopo = 2;
//	const int nRows = 19, nTopo = 2;
//	const int nRows = 165, nTopo = 2;
//	const int nRows = 193, nTopo = 2;
//	const int nRows = 185, nTopo = 2;
//	const int nRows = 183, nTopo = 2;

	vector<vector<double> > points1(nRows, vector<double>(nTopo+1));

	cout.unsetf(ios::floatfield);
	cout.precision(10);

	ifstream myReadFile;
	myReadFile.open(argv[1]);		//12 water waterwire
//	myReadFile.open("common_classical.txt");	//solvated zundel
	while (!myReadFile.eof()) {
		for (unsigned int i = 0; i < points1.size(); i++) {
			for (unsigned int j = 0; j < points1[i].size(); j++) {
				myReadFile >> points1[i][j];
			}
		}
	}
	myReadFile.close();

	vector<double> distance;
	
	for (unsigned int i = 0; i < points1.size(); i++) {
		distance.push_back(points1[i].back());
	}

	MIN minObject(points1);
	vector<vector<double> > shiftedToMinimum;
	minObject.shiftToMinimum(shiftedToMinimum);
	
	vector<vector<unsigned int> > autosubPoints(nTopo, vector<unsigned int>(nTopo));
	myReadFile.open(argv[2]);	//12 water waterwire
//	myReadFile.open("autosubPoints_classical.txt");		//solvated zundel
	while (!myReadFile.eof()) {
		for (unsigned int i = 0; i < autosubPoints.size(); i++) {
			for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
				myReadFile >> autosubPoints[i][j];
			}
		}
	}
	myReadFile.close();

	vector<vector<double> > epsilon(nTopo,vector<double> (nTopo));
	myReadFile.open(argv[3]);		//12 water waterwire
//	myReadFile.open("epsilon_classical.txt");		//solvated zundel
	while (!myReadFile.eof()) {
		for (unsigned int i = 0; i < autosubPoints.size(); i++) {
			for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
				myReadFile >> epsilon[i][j];
			}
		}
	}
	myReadFile.close();

	vector<vector<double> > shiftedAtAutosub;
	AUTOSUBSHIFT autosubShiftObject(shiftedToMinimum, autosubPoints);
	autosubShiftObject.shiftAtAutosub(shiftedAtAutosub);

	vector<vector<double> > oniomXsEnergy;
	vector<vector<double> > justShiftedEnergy;
	Oniom_Xs oniomXsObject(shiftedAtAutosub, epsilon,autosubPoints,distance);
	oniomXsObject.getOniomXsEnergy(oniomXsEnergy,nRows,nTopo);
	oniomXsObject.getJustShifted(justShiftedEnergy, nRows, nTopo);

	return 0;
}
