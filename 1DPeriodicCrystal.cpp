#include "1DPeriodicCrystal.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

double kroneckerD(int i, int j);
//template <typename T> int sgn(T val);

SpinWaveProblem1D::SpinWaveProblem1D(int N, int baseSplit, int kSteps, int omegaSteps, double H2, double l1, double l2, double d, double h, double omegaStart,
									double omegaEnd, double kStart, double kEnd, bool debug){
	this->N = N;
	this->baseSplit = baseSplit;
	this->kSteps = kSteps;
	this->omegaSteps = omegaSteps;
	this->omegaStart = omegaStart;
	this->omegaEnd = omegaEnd;
	this->kStart = kStart;
	this->kEnd = kEnd;
	this->debug = debug;

	//всё в СГС
	this->l1 = l1;
	this->l2 = l2;
	this->d = d;
	this->h = h;
	gamma = 1.76e7;
	H1 = 400;
	this->H2 = H2;
	M4pi = 1750;

	omegaH1 = gamma * H1;
	omegaH2 = gamma * H2;
	omegaM = gamma * M4pi;

	goThroughGrid();
}

void SpinWaveProblem1D::goThroughGrid() {

	omegaDelta = (omegaEnd - omegaStart);
	double kDelta = (kEnd - kStart);


	vector<double> suspiciousOmega;
	vector<complex<double> > determinants;


	std::ostringstream filename;
	filename << "double results N " << N << " H2 " << std::setprecision(3) << H2 << " l1 to d " << l1/d << " h to d " << h/d << " k to b_1 "
				<< kStart/b(1) << "-" << kEnd/b(1) << " O " << std::setprecision(5) << omegaStart << "-" << omegaEnd;
	ofstream fout1;
	fout1.open(filename.str());
	fout1 << "k/b_1 omega omegaIdeal1 omegaIdeal2 These are results for grid: k: "
	<< kSteps << " points, omega starting: " << omegaSteps << " points" << endl;

	ofstream fout2;
	if(baseSplit != 0) {
		filename.str("");
		filename.clear();
		filename << "double exact results N " << N << " H2 " << std::setprecision(3) << H2 <<  " l1 to d " << l1/d << " k to b_1 "
				<< kStart/b(1) << "-" << kEnd/b(1) << " O " << std::setprecision(5) << omegaStart << "-" << omegaEnd;
		fout2.open(filename.str());
		fout2 << "k/b_1 omega omegaIdeal1 omegaIdeal2 These are results for grid: k: "
		<< kSteps << " points, omega starting: " << omegaSteps << " points" << endl;
	}


	for(double k = kStart; k < kEnd; k += kDelta/kSteps) {

		cout << (k - kStart)/kDelta * 100 << "%:  \tCalculating layer k/b_1 = " << k/b(1) << endl;

		for(double omega = omegaStart+1; omega < omegaStart + omegaDelta * (1.0 + 1.0/omegaSteps*0.5); omega += omegaDelta/omegaSteps) {
			determinants.push_back(findDeterminant(k, omega));
		}

		for(unsigned int i = 6; i < determinants.size()-6; i++) {
			if(isMinimum(determinants, i, 6)) {
				if(baseSplit != 0)
					suspiciousOmega.push_back(checkNull(k, omegaStart + omegaDelta/omegaSteps*(i-1),omegaStart + omegaDelta/omegaSteps*(i+1)));
				fout1 << k/b(1) << " " << omegaStart + omegaDelta/omegaSteps * (i) << " " << omegaIdeal(k, omegaH1, 1)
								<< " " << omegaIdeal(k, omegaH2, 1) << endl;
			}
		}

		if(baseSplit != 0) {
			for(unsigned int i = 0; i < suspiciousOmega.size(); i++) {
				if(suspiciousOmega[i] != 0)
					fout2 << k/b(1) << " " << suspiciousOmega[i] << " " << omegaIdeal(k, omegaH1, 1)
								<< " " << omegaIdeal(k, omegaH2, 1) << endl;
			}
		}

		determinants.clear();
		suspiciousOmega.clear();
	}

	return;
}

std::complex<double> SpinWaveProblem1D::findDeterminant(double k, double omega) {
	fillMatrix1(k, omega);
	ces.compute(Matrix1);
	fixEigens();
	refillMatrix1(k);
	return Matrix1.determinant();
}


bool SpinWaveProblem1D::isMinimum(vector<complex<double> >& determinants, int i, int depth) {
	if(depth < 1)
		return true;
	if(abs(determinants[i]) < abs(determinants[i+depth]) && abs(determinants[i]) < abs(determinants[i-depth]))
		return isMinimum(determinants, i, depth-1);
	else
		return false;
}

double SpinWaveProblem1D::mu1(double omega) {
	return (omegaH1 * (omegaH1 + omegaM) - omega*omega) / (omegaH1*omegaH1 - omega*omega);
}

double SpinWaveProblem1D::mu2(double omega) {
	return (omegaH2 * (omegaH2 + omegaM) - omega*omega) / (omegaH2*omegaH2 - omega*omega);
}

std::complex<double> SpinWaveProblem1D::M(int n, double omega) {
	if (n == 0)
		return (mu1(omega) * (l1+l2-d*4) + mu2(omega) * (2*d)) / (l1+l2-2*d);

	else {
		std::complex<double> exp1 = exp(-1.0i * b(n) * (l1*0.5-d)) - 1.0 + exp(-1.0i * b(n) * (l1*0.5+l2-2*d))- exp(-1.0i * b(n) * (l1*0.5)) + exp(-1.0i * b(n) * (l1+l2-2.0*d)) - exp(-1.0i * b(n) * (l1*0.5+l2-d));
		std::complex<double> exp2 = exp(-1.0i * b(n) * (l1*0.5)) - exp(-1.0i * b(n) * (l1*0.5-d)) + exp(-1.0i * b(n) * (l1*0.5+l2-d)) - exp(-1.0i * b(n) * (l1*0.5+l2-2.0*d));

		return 1.0i / ((l1+l2-2.0*d)*b(n)) * (mu1(omega)* exp1 + mu2(omega)* exp2);
	}
}

double SpinWaveProblem1D::b(int m) {
	return m * 2.0 * M_PI / (l1+l2-2.0*d);
}

void SpinWaveProblem1D::fillMatrix1(double k, double omega) {
	Matrix1 = MatrixXcd(2*N + 1, 2*N + 1);

	for(int p = -N; p <= N; p++) {
		for(int n = -N; n <= N; n++) {
			Matrix1(p + N, n + N) =  -1.0 * M(p - n, omega) * (k + b(p)) * (k + b(n));
		}
	}
	return;
}

void SpinWaveProblem1D::fixEigens() {
	eigenValues = ces.eigenvalues();
	eigenVectors = ces.eigenvectors();

//	lambda squared into lambda
	for(int i = 0; i < 2*N+1; i++) {
		eigenValues(i) = sqrt(eigenValues(i).real());
	}

	for(int i = 0; i < 2*N+1; i++) {
		for(int j = 0; j < 2*N+1; j++) {
			if(abs(eigenVectors(i, j).real()) / abs(eigenVectors(i, j).imag()) > pow(10, 20) )
				eigenVectors(i, j) = eigenVectors(i, j).real();
			if(abs(eigenVectors(i, j).imag()) / abs(eigenVectors(i, j).real()) > pow(10, 20) )
				eigenVectors(i, j) = eigenVectors(i, j).imag();
		}
	}

	return;
}

void SpinWaveProblem1D::refillMatrix1(double k) {
	Matrix1 = MatrixXcd(16*N + 8, 16*N + 8);




//1
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 2*N+1) = A_1(i+N, j+N, k, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 4*N+2) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 8*N+4) = -1.0 * C_1(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 10*N+5) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 12*N+6) = -1.0 * D_1(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N, j+N + 14*N+7) = 0.0;






//2
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N) = -1.0 * A_1(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 2*N+1) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 4*N+2) = -1.0 * B_1(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 8*N+4) = C_1(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 10*N+5) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 12*N+6) = D_1(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 2*N+1, j+N + 14*N+7) = 0.0;






//3
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N) = -1.0 * B_1(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 2*N+1) = 0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 4*N+2) = -1.0 * A_1(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 8*N+4) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 10*N+5) = C_1(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 12*N+6) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 4*N+2, j+N + 14*N+7) = -1.0 * D_1(i+N, j+N, 0.5*h);






//4
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 2*N+1) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 4*N+2) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 6*N+3) = A_1(i+N, j+N, k, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 8*N+4) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 10*N+5) = -1.0 * C_1(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 12*N+6) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 6*N+3, j+N + 14*N+7) = D_1(i+N, j+N, 0.5*h+d);






//5
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 2*N+1) = -1.0 * A_2(i+N, j+N, k, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 4*N+2) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 8*N+4) = D_2(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 10*N+5) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 12*N+6) = -1.0 * C_2(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 8*N+4, j+N + 14*N+7) = 0.0;






//6
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N) = A_2(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 2*N+1) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 4*N+2) = -1.0 * B_2(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 8*N+4) = -1.0 * D_2(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 10*N+5) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 12*N+6) = C_2(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 10*N+5, j+N + 14*N+7) = 0.0;






//7
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N) = B_2(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 2*N+1) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 4*N+2) = -1.0 * A_2(i+N, j+N, k, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 6*N+3) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 8*N+4) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 10*N+5) = D_2(i+N, j+N, 0.5*h);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 12*N+6) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 12*N+6, j+N + 14*N+7) = C_2(i+N, j+N, 0.5*h);






//8
	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 2*N+1) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 4*N+2) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 6*N+3) = A_2(i+N, j+N, k, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 8*N+4) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 10*N+5) = -1.0 * D_2(i+N, j+N, 0.5*h+d);

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 12*N+6) = 0.0;

	for(int i = -N; i <= N; i++)
		for(int j = -N; j <= N; j++)
			Matrix1(i+N + 14*N+7, j+N + 14*N+7) = -1.0 * C_2(i+N, j+N, 0.5*h+d);

	return;
}

double SpinWaveProblem1D::checkNull(double& k, double omega1, double omega2, double& startingDet, double& foundOmega) {
//	if(Debug == 1)
//		cout << "Checking (" << omega1 << ", " << omega2 << ")" << endl;
//	if((omega2 - omega1) <= omegaDelta/omegaSteps/(baseSplit*baseSplit + baseSplit)) {
//		if(Debug == 1)
//			cout << "Found nothing, it's too small" << endl;
//		return 0;
//	}
//
//	if(Debug == 1)
//		cout << "It isn't too small, calculating determinants" << endl;
//	vector<double> determinants;
//	for(double omega = omega1; omega < omega2 + (omega2-omega1) / baseSplit * 0.5; omega += (omega2-omega1) / baseSplit) { //в целях чтобы точно дойти до
//		determinants.push_back(abs(findDeterminant(k, omega)));															//омеги 2
//	}
//
//	for(int i = 0;  i <= baseSplit; i++) {
//		if(Debug == 1)
//			cout << "Checking determinant " << determinants[i] << endl;;
//		if((startingDet/determinants[i]) > pow(10, N)) {
//			if(Debug == 1) {
//				cout << "Found omega = " << omega1 + (omega2-omega1) / baseSplit * (i) << " with determinant " << determinants[i] << endl;
//				cin.get();
//			}
//			return omega1 + (omega2-omega1) / baseSplit * (i);
//		}
//		if(Debug == 1)
//			cout << "It didn't fit" << endl;
//	}
//
//	for(int i = 1;  i < baseSplit; i++) {
//		if(Debug == 1) {
//			cin.get();
//			cout << "Checking determinant " << determinants[i] << " again" << endl;
//		}
//
//		if(determinants[i] < startingDet && determinants[i] <= determinants[0] && determinants[i] <= determinants[baseSplit]) {
//			if(Debug == 1)
//				cout << "It is interesting" << endl;
//			if(determinants[i-1] < determinants[i+1] && determinants[i-1] <= determinants[0] && determinants[i-1] <= determinants[10]) {
//				if(Debug == 1)
//					cout << "Checking previous interval" << endl;
//				foundOmega = checkNull(k, omega1 + (omega2 - omega1) / baseSplit * (i-1), omega1 + (omega2 - omega1) / baseSplit * i, startingDet, foundOmega);
//				if(foundOmega != 0)
//					return foundOmega;
//			}
//
//			else if(determinants[i-1] > determinants[i+1] && determinants[i+1] <= determinants[0] && determinants[i+1] <= determinants[10]) {
//				if(Debug == 1)
//					cout << "Checking further interval" << endl;
//				foundOmega = checkNull(k, omega1 + (omega2 - omega1) / baseSplit * i, omega1 + (omega2 - omega1) / baseSplit * (i+1), startingDet, foundOmega);
//				if(foundOmega != 0)
//					return foundOmega;
//			}
//		}
//	}
//
	return 0;
}

double SpinWaveProblem1D::checkNull(double& k, double omega1, double omega2) {
	double foundOmega;
	complex<double> det1 = findDeterminant(k, omega1);
	complex<double> det2 = findDeterminant(k, omega2);
	double det = min(abs(det1), abs(det2));
	if(Debug == 1)
		cout << "Checking (" << omega1 << ", " << omega2 << "), starting determinant = " << det << endl;

	return checkNull(k, omega1, omega2, det, foundOmega);
}

double SpinWaveProblem1D::omegaIdeal(double k, double omegaH, int mode) {
	return sqrt(omegaH * ( omegaH + omegaM*( 1.0 -  (1.0 - exp(-1.0*abs(k*d)))/ abs(k*d) ) ) );
}

std::complex<double> SpinWaveProblem1D::A_1(int n, int s, double k, double x) {
	return kroneckerD(n, s) * exp(-1.0 * abs(k + b(n-N)) * x);
}
std::complex<double> SpinWaveProblem1D::A_2(int n, int s, double k, double x) {
	return abs(k + b(n-N)) * A_1(n, s, k, x);
}
std::complex<double> SpinWaveProblem1D::B_1(int n, int s, double k, double x) {
	return kroneckerD(n, s) * exp(abs(k + b(n-N)) * x);
}
std::complex<double> SpinWaveProblem1D::B_2(int n, int s, double k, double x) {
	return abs(k + b(n-N)) * B_1(n, s, k, x);
}
std::complex<double> SpinWaveProblem1D::C_1(int n, int s, double x) {
	return eigenVectors(n, s) * cos(eigenValues(s)*x);
}
std::complex<double> SpinWaveProblem1D::C_2(int n, int s, double x) {
	return eigenValues(s) * C_1(n, s, x);
}
std::complex<double> SpinWaveProblem1D::D_1(int n, int s, double x) {
	return eigenVectors(n, s) * sin(eigenValues(s)*x);
}
std::complex<double> SpinWaveProblem1D::D_2(int n, int s, double x) {
	return eigenValues(s) * D_1(n, s, x);
}

double kroneckerD(int i, int j) {
	if(i == j)
		return 1.0;

	return 0.0;
}

//template <typename T> int sgn(T val) {
//    return (T(0) < val) - (val < T(0));
//}
