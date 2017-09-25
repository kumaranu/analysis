#include"AutosubShift.h"

AUTOSUBSHIFT::AUTOSUBSHIFT()
{
}

AUTOSUBSHIFT::AUTOSUBSHIFT(vector<vector<double> > shiftedAtMin, vector<vector<unsigned int> > autosubPoints) {
	vector<double> point;

	for (unsigned int i = 0; i < shiftedAtMin.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < shiftedAtMin[i].size(); j++) {
			point.push_back(shiftedAtMin[i][j]);
		}
		newShiftedAtMinimum.push_back(point);
	}

	vector<unsigned int> autosub_Point;
	for (unsigned int i = 0; i < autosubPoints.size(); i++) {
		autosub_Point.clear();
		for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
			autosub_Point.push_back(autosubPoints[i][j]);
		}
		newAutosubPoints.push_back(autosub_Point);
	}
}

void AUTOSUBSHIFT::setPoints(vector<vector<double> > points) {
	vector<double> point;
	for (unsigned int i = 0; i < points.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < points[i].size(); j++) {
			point.push_back(points[i][j]);
		}
		newShiftedAtMinimum.push_back(point);
	}
}

void AUTOSUBSHIFT::setAutosubPoints(vector<vector<unsigned int> > autosubPoints) {
	vector<unsigned int> autosub_point;
	for (unsigned int i = 0; i < autosubPoints.size(); i++) {
		autosub_point.clear();
		for (unsigned int j = 0; j < autosubPoints[i].size(); j++) {
			autosub_point.push_back(autosubPoints[i][j]);
		}
		newAutosubPoints.push_back(autosub_point);
	}
}

void AUTOSUBSHIFT::shiftAtAutosub(vector<vector<double> > & autosubShifted) {
	double shift;
	//	cout << "Energy values around the autosub points and the shift(i.e. the difference between the two points) is shown here!" << endl;
	//Editted here!!!!!!!!!!!!!!!!!!!!!!!!! I am restricting the loop to only 1 included topology to shift the data
	for (unsigned int i = 0; i < newAutosubPoints.size(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			if (newAutosubPoints[i][j] != 0) {
				shift = -1 * (newShiftedAtMinimum[newAutosubPoints[i][j]-1][i] - newShiftedAtMinimum[newAutosubPoints[i][j]-1][j]);
				//cout << "Autosub point is " << newAutosubPoints[i][j] << " i is " << i << " and j is " << j << " Point before the hop no. " << i << " is " << newShiftedAtMinimum[newAutosubPoints[i][j]][j] << ", ";
				//cout << "Point right after the hop no. " << i << " is " << newShiftedAtMinimum[newAutosubPoints[i][j]][i] << " ";
			//	cout << "and the shift is " << shift << "." << endl;
			//	cout << "newAutosub point is " << newAutosubPoints[i][j] << endl;
			//	cout << "newShiftedAtMinimum[newAutosubPoints[i][j]][j] is " << newShiftedAtMinimum[newAutosubPoints[i][j]][j] << endl;
			//	cout << "newShiftedAtMinimum[newAutosubPoints[i][j]][i] is " << newShiftedAtMinimum[newAutosubPoints[i][j]][i] << endl;
				for (unsigned int k = 0; k < newShiftedAtMinimum.size(); k++) {
					newShiftedAtMinimum[k][j] = newShiftedAtMinimum[k][j] - shift;
				}
			}
		}
	}
	vector<double> point;
	for (unsigned int i = 0; i < newShiftedAtMinimum.size(); i++) {
		point.clear();
		for (unsigned int j = 0; j < newShiftedAtMinimum[i].size(); j++) {
			point.push_back(newShiftedAtMinimum[i][j]);
		}
		autosubShifted.push_back(point);
	}

	//cout << "Autosub shifted is printed here." << endl;

	ofstream myfile7;
	myfile7.unsetf(ios::floatfield);
	myfile7.precision(10);
	myfile7.open("AutosubShifted_data.txt");

	for (unsigned int i = 0; i < autosubShifted.size(); i++) {
		for (unsigned int j = 0; j < autosubShifted[i].size(); j++) {
			//cout << autosubShifted[i][j] << " ";
			myfile7 << autosubShifted[i][j] << " ";
		}
		//cout << endl;
		myfile7 << endl;
	}
	myfile7.close();
}

AUTOSUBSHIFT::~AUTOSUBSHIFT()
{
}
