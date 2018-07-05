#include <iostream>
#include <string>
#include <Eigen>
#include "1DPeriodicCrystal.h"

using namespace Eigen;
using namespace std;

int main(int argc, char* argv[])
{
//	if(argc != 4) {
//		cout << "Not enough arguments" << endl;
//		return -1;
//	}
//
//	cin.get();
//	SpinWaveProblem1D A = SpinWaveProblem1D(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
//	SpinWaveProblem1D A = SpinWaveProblem1D(1, 200, 400, 1e-5, 10, 200);
//	SpinWaveProblem1D B = SpinWaveProblem1D(1, 200, 400, 1e-5, 10, 400);
//	SpinWaveProblem1D C = SpinWaveProblem1D(1, 200, 400, 1e-5, 10, 550);
	SpinWaveProblem1D D = SpinWaveProblem1D(1, 500, 1000, 1e-5, 10, 750);
//	SpinWaveProblem1D E = SpinWaveProblem1D(1, 200, 400, 1e-5, 10, 1100);
}

