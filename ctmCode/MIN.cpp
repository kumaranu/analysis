#include"MIN.h"

MIN::MIN() {
	//default constructor
}

//overload constructor
MIN::MIN(vector<vector<double> > common) {
	vector<double> point;

	for (unsigned int i = 0; i < common.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < common[i].size(); j++) {
			point.push_back(common[i][j]);
		}
		newCommon.push_back(point);
	}
}

void MIN::setPoints(vector<vector<double> > common) {
	vector<double> point;

	for (unsigned int i = 0; i < common.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < common[i].size(); j++) {
			point.push_back(common[i][j]);
		}
		newCommon.push_back(point);
	}
}

//void shiftToMinimum(vector<vector<double>> &);
void MIN::shiftToMinimum(vector<vector<double> >  & shiftedToZero) {
	for (unsigned int i = 0; i < newCommon[0].size()-1; i++) {
		newMinimum.push_back(newCommon[0][i]);
	}

	for (unsigned int i = 0; i < newCommon[0].size() - 1; i++) {
		for (unsigned int j = 0; j < newCommon.size(); j++) {
			for (unsigned int k = 0; k < newCommon[0].size()-1; k++) {
				if (newCommon[j][k] < newMinimum[k]) {
					newMinimum[k] = newCommon[j][k];
				}
			}
		}
	}

	vector<double> points;
	for (unsigned int i = 0; i < newCommon.size(); i++) {
		points.clear();
		for (unsigned int j = 0; j < newCommon[0].size()-1; j++) {
			points.push_back(newCommon[i][j] - newMinimum[j]);
		}
		newShiftedToZero.push_back(points);
		shiftedToZero.push_back(points);
	}

	ofstream myfile7;
	myfile7.unsetf(ios::floatfield);
	myfile7.precision(10);		
	myfile7.open("shifted_to_min.txt");

	for (unsigned int i = 0; i < shiftedToZero.size(); i++) {
		for (unsigned int j = 0; j < shiftedToZero[i].size(); j++) {
			myfile7 << shiftedToZero[i][j] << " ";
		}
		myfile7 << endl;
	}
}

MIN::~MIN() {

}
