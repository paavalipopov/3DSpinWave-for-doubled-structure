
#ifndef PERIODICCRYSTAL1D
#define PERIODICCRYSTAL1D

#include <vector>
#include <Eigen>

class SpinWaveProblem1D {
public:
	SpinWaveProblem1D(int N, int kStpes, int omegaSteps, double precision, int baseSplit, double H2);
private:
	int N, baseSplit; //первое - количество гармоник в разложении, второе - мелкость разбиения при уточнении минимумов
	int kSteps, omegaSteps; //количество шагов по сетке ка-омега
	double omegaDelta; //мелкость по оси частот
	double l1, l2, d; //линейные размеры
	double gamma; //гиромагнитное отношения для циклической частоты
	double H1, H2, M4pi;  //параметры поля и намагниченности
	double omegaH1, omegaH2, omegaM; //всякие ферромагнитные частоты
	double mu1(double omega); //коэффициенты в тензоре мю в разных полях
	double mu2(double omega);
	std::complex<double> M(int n, double omega); //для результатов разложения компоненты мю в ряд
	double b(int m); //блоховское волновое число

	int Debug = 0; //отсекание левых почти нулёвых значений
	double precision; //точность отсекания этих значений - всё что меньше его удаляется

	void fillMatrix1(double k, double omega);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
	void fixEigens();
	void refillMatrix1(double k);
	void goThroughGrid();
	std::complex<double> findDeterminant(double k, double omega);
	double checkNull(double& k, double omega1, double omega2, double& startingDet);
	double checkNull(double& k, double omega1, double omega2, double& startingDet, double& foundOmega);
	bool isMinimum(std::vector<std::complex<double> >& determinants, int i, int depth);

	Eigen::MatrixXcd Matrix1;
	Eigen::MatrixXcd eigenVectors;
	Eigen::VectorXcd eigenValues;
};

#endif // PERIODICCRYSTAL1D
