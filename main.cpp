#include <iostream>
#include <string>
#include <Eigen>
#include "1DPeriodicCrystal.h"

using namespace Eigen;
using namespace std;

int main(int argc, char* argv[])
{
	double l1 = 9e-4;
	double l2 = 9e-4;
	double d = 3e-4;
	double gamma = 1.76e7;
	double H1 = 400;
	double H2 = 420;
	double M4pi = 1750;

	double omegaH1 = gamma * H1;
	double omegaH2 = gamma * H2;
	double omegaM = gamma * M4pi;

//	double omegaStart = 7.049e9;
//	double omegaStart = omegaH1;
	double omegaStart = omegaH2;
//	double omegaStart = 1;
//	double omegaEnd = 1.0e10;
//	double omegaEnd = omegaH2;
	double omegaEnd = sqrt(omegaH1*(omegaH1 + omegaM));
//	double omegaEnd = sqrt(omegaH2*(omegaH2 + omegaM));
//	double omegaEnd = 2;
//	double ka = 1.0;
	double ka = 0.1;

	double b = 2.0 * M_PI / (l1+l2-2.0*d);

	double kStart= -2.0 * b;
	double kEnd = 2.0 * b;

	omegaEnd = ka * (omegaEnd - omegaStart) + omegaStart;


//	if(argc != 4) {
//		cout << "Not enough arguments" << endl;
//		return -1;
//	}
//SpinWaveProblem1D(int N, int baseSplit, int kStpes, int omegaSteps, double H2, double l1, double l2, double omegaStart, double omegaEnd, double kStart, double kEnd);
//	cin.get();
//	SpinWaveProblem1D A = SpinWaveProblem1D(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));

	SpinWaveProblem1D D = SpinWaveProblem1D(5, 0, 400, 400, H2, l1, l2, d, omegaStart, omegaEnd, kStart, kEnd);
}

