#include"OniomXs.h"
#include"AutosubShift.h"

//Default constructor
Oniom_Xs::Oniom_Xs() {
}

//Overloaded constructor
Oniom_Xs::Oniom_Xs(vector<vector<double> >& shiftedAtAutosub, vector<vector<double> >& epsilon, vector<vector<unsigned int> > autosubPoints, vector<double> distance) {
	vector<double> point;
	for (unsigned int i = 0; i < shiftedAtAutosub.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < shiftedAtAutosub[i].size(); j++) {
			point.push_back(shiftedAtAutosub[i][j]);
		}
		newShiftedAtAutosub.push_back(point);
	}
	for (unsigned int i = 0; i < epsilon.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < epsilon[i].size(); j++) {
			point.push_back(epsilon[i][j]);
		}
		newEpsilon.push_back(point);
	}
	vector<unsigned int> points1;
	for (unsigned int i = 0; i < autosubPoints.size(); i++) {
		points1.clear();
		for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
			points1.push_back(autosubPoints[i][j]);
		}
		newAutosubPoints.push_back(points1);
	}
	for (unsigned int i = 0; i < distance.size(); i++) {
		newDistance.push_back(distance[i]);
	}
}

void Oniom_Xs::getOniomXsEnergy(vector<vector<double> >& oniomXsEnergy, const int nRows, const int nTopo) {
	MatrixXcd H(nTopo, nTopo);
	MatrixXd W(nTopo, nTopo);
	MatrixXd X(nTopo, nTopo);
	MatrixXcd eigenvectors;
	VectorXcd eigenvalues;

	vector<complex<double> > point;
//	double square;		//it is used in the code to compute the square of the diff in energy of the two topologies
	complex<double> iota(0.0,1.0);	//Defining the complex number iota

	cout.unsetf(ios::floatfield);
	cout.precision(10);
	
	ofstream myfile;
	myfile.unsetf(ios::floatfield);
	myfile.precision(10);
	myfile.open("X.txt");

	for (unsigned int i = 0; i < newShiftedAtAutosub.size(); i++) {
		for (int j = 0; j < nTopo; j++) {		//Now the evb hamiltonian will be calculated at each step
			for (int k = 0; k < nTopo; k++) {
				if (j == k) {		//Diagonal elements
					X(j, k) = 0;
					W(j, k) = 0;
					H(j, k) = newShiftedAtAutosub[i][j];
				}
				else if ((j != k) && (newAutosubPoints[j][k] == 0)) { //Non-diagonal elements where autosub_points are zero
					X(j, k) = 0;	//not used in the 2by2 case
					W(j, k) = 0;
					H(j, k) = 0;
				}
				else if ((j != k) && (newDistance[i] > newDistance[newAutosubPoints[j][k]] - newEpsilon[j][k]) && (newDistance[i] < newDistance[newAutosubPoints[j][k]] + newEpsilon[j][k])&&(newAutosubPoints[j][k] != 0)) {
					X(j, k) = (newDistance[i] - newDistance[newAutosubPoints[j][k]] + newEpsilon[j][k]) / (2 * newEpsilon[j][k]);
					W(j, k) = 0.5 * erf(5 * X(j, k) - 2.5) + 0.5;
					H(j, k) = iota*sqrt(abs(W(j, k)*(W(j, k) - 1)*(newShiftedAtAutosub[i][j] - newShiftedAtAutosub[i][k])*(newShiftedAtAutosub[i][j] - newShiftedAtAutosub[i][k])));
				}
				else if ((j != k) && (newDistance[i] < newDistance[newAutosubPoints[j][k]] - newEpsilon[j][k])&&(newAutosubPoints[j][k] != 0)) {
					X(j, k) = 0;
					W(j, k) = 0;
					H(j, k) = 0;
				}
				else if ((j != k) && (newDistance[i] > newDistance[newAutosubPoints[j][k]] + newEpsilon[j][k])&&(newAutosubPoints[j][k] !=0)) {
					X(j, k) = 0;
					W(j, k) = 0;
					H(j, k) = 0;
				}
			}
		}

		ComplexEigenSolver<MatrixXcd> eigensolver;
		eigensolver.compute(H);
		if (eigensolver.info() != Success) abort();
		eigenvalues = eigensolver.eigenvalues();
		eigenvectors = eigensolver.eigenvectors();

		point.clear();
		point.push_back(newShiftedAtAutosub[i][0]);
		point.push_back(newShiftedAtAutosub[i][1]);
		point.push_back(eigenvalues(0));
		point.push_back(newDistance[i]);
		newOniomXsEnergy.push_back(point);
		//square = (newShiftedAtAutosub[i][0] - newShiftedAtAutosub[i][1])*(newShiftedAtAutosub[i][0] - newShiftedAtAutosub[i][1]);
		myfile << real(newOniomXsEnergy[i][0]) << " " << real(newOniomXsEnergy[i][1]) << " " << real(newOniomXsEnergy[i][2]) << " " << real(newOniomXsEnergy[i][3]) << endl;
	}
	myfile.close();
}

void Oniom_Xs::getJustShifted(vector<vector<double> > & justShiftedEnergy, const int nRows,const int ntopo) {
	vector<double> points;
	for (unsigned int i = 0; i < newShiftedAtAutosub.size(); i++) {
		points.clear();
		if (i < newAutosubPoints[0][1]) {
			cout << i << endl;
			points.push_back(newShiftedAtAutosub[i][1]);
			points.push_back(newDistance[i]);
			justShiftedEnergy.push_back(points);
		}
		else {
			cout << i << endl;
			points.push_back(newShiftedAtAutosub[i][0]);
			points.push_back(newDistance[i]);
			justShiftedEnergy.push_back(points);
		}
	}

	ofstream myfile;
	myfile.unsetf(ios::floatfield);
	myfile.precision(10);
	myfile.open("justShifted.txt");

	for (unsigned int i = 0; i < justShiftedEnergy.size(); i++) {
		for (unsigned int j = 0; j < justShiftedEnergy[i].size(); j++) {
			myfile << justShiftedEnergy[i][j] << " ";
		}
		myfile << endl;
	}
	myfile.close();
}

void Oniom_Xs::setShiftedAtAutosub(vector<vector<double> > & shiftedAtAutsub) {
	vector<double> points;
	for (unsigned int i = 0; i < shiftedAtAutsub.size(); i++) {
		points.clear();
		for (unsigned int j = 0; j < shiftedAtAutsub[i].size(); j++) {
			points.push_back(shiftedAtAutsub[i][j]);
		}
		newShiftedAtAutosub.push_back(points);
	}
}

void Oniom_Xs::setEpsilon(vector<vector<double> > & epsilon) {
	vector<double> points;
	for (unsigned int i = 0; i < epsilon.size(); i++) {
		points.clear();
		for (unsigned int j = 0; j < epsilon[i].size(); j++) {
			points.push_back(epsilon[i][j]);
		}
		newEpsilon.push_back(points);
	}
}

void Oniom_Xs::setAutosubPoints(vector<vector<unsigned int> >& autosubPoints) {
	vector<unsigned int> points1;
	for (unsigned int i = 0; i < autosubPoints.size(); i++) {
		for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
			points1.push_back(autosubPoints[i][j]);
		}
		newAutosubPoints.push_back(points1);
	}
}

void Oniom_Xs::setDistance(vector<double>& distance) {
	for (unsigned int i = 0; i < distance.size(); i++) {
		newDistance.push_back(distance[i]);
	}
}

Oniom_Xs::~Oniom_Xs() {
}
