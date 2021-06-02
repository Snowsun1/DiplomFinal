#include <complex>
#include <cstdlib>
#include <vector>

#include "VoyevodinMethod.h"
#include <vector>
#include <cstdlib>
#include "OdeSolver.h"
#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include <functional>
#include <complex>
#include "OdeSolver.h"

std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;

const double pi = 3.1415926538;




double fi(double x) {
	return 1 + x;
}

std::complex<double> top_func(std::complex<double> v_1, std::complex<double> v_0, double x, double tau, double kappa,
	std::complex<double> lambda_1, std::complex<double> mu) {
	return -1 / x * (v_1 - (v_1 - (lambda_1 * v_0 / x) / lambda_1 + 2.0 * mu + tau * fi(x)) +
		((lambda_1 + 2.0 * mu + x * tau * fi(x) + tau * fi(x)) / x) * v_0) - kappa * kappa * v_0;
}

std::complex<double> bottom_func(std::complex<double> v_1, std::complex<double> v_0, double x, double tau, std::complex<double> lambda_1, std::complex<double> mu) {
	return (v_1 - (lambda_1 * v_0) / x) / (lambda_1 + 2.0 * mu + tau * fi(x));
}



int main()
{
	auto rightPartFun = [](double x)-> std::complex<double> {return { sin(pi * x), cos(pi * x) }; }; // Нужно изменить правую часть на ту что в нашем уравнении 
																									// я не знаю как правильно её задать, нужна ваша помощь

	setlocale(0, "");

	double b = 20.0;
	double mu_1 = 0.5;
	double mu_2 = 1.0;
	double lambda_2 = 1.5;
	double tau = 0.01;
	double P = 1;
	double ksi = 0.8;
	std::complex<double> right(1, 0);
	std::complex<double> left(0, 0);
	std::complex<double> i(0, 1);
	Parameters::kind = FIRST;

	std::vector<std::complex<double>> result_vec_1;
	std::vector<std::vector<std::complex<double>>> matrix;

	for (double tau = 0; tau <= 1; tau += 0.05) {
		std::vector<std::complex<double>> result_vec_1;
		for (double kappa = 0; kappa <= 1; kappa += 0.05) {

			auto lambda_1 = (1.0 + i * b * lambda_2) / (1.0 + i * b * kappa);
			auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
			boundary_value_problem<std::complex<double>> bvp = {
				{
					{
						[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu); },
						[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu); },
					}
				}, {{1, left}}, {{1,right}}, 0.5
			};
			auto sol = bvp.solve();
			result_vec_1.push_back(sol[0]);
		}
		matrix.push_back(result_vec_1);
	}
	const size_t size = 20;
	double h_x = 1.0 / size;
	double h_s = 1.0 / size;
	// размер СЛАУ
	//const size_t size = 20;
	//double h_x = 1.0 / size;
	//double h_s = 1.0 / size;
	//// создаём матрицу СЛАУ
	//std::vector<std::vector<std::complex<double>>> matrix;
	//for (size_t i = 0; i < size; i++)
	//{
	//	const auto x = (i + 0.5) * h_x;
	//	std::vector<std::complex<double>> row(size);
	//	for (size_t ii = 0; ii < size; ii++)
	//	{
	//		const auto s = (ii + 0.5) * h_s;
	//		row[ii] = h_s * kernelFun(x, s);
	//	}
	//	matrix.push_back(row);
	//}
	// точное решение
	std::vector<std::complex<double>> exactSolution(size);
	for (size_t i = 0; i < size; i++)
	{
		const auto s = (i + 0.5) * h_s;
		exactSolution[i] = rightPartFun(s);
	}

	// вектор правой части СЛАУ
	std::vector<std::complex<double>> rightPart(size, 0);
	for (size_t i = 0; i < size; i++)
	{
		for (size_t ii = 0; ii < size; ii++)
		{
			rightPart[i] += matrix[i][ii] * exactSolution[ii];
		}
	}

	std::vector<VoyevodinMethod> arr_v = { 
		VoyevodinMethod(matrix, rightPart, h_s, Neumann, Neumann)
	};
	std::vector< std::vector< std::complex<double> > > solution = {
		arr_v[0].solution()
	};
	//const VoyevodinMethod voyevodin_method(matrix, rightPart, h_s, Dirichle, Dirichle);
	//auto solution = voyevodin_method.solution();
	for (auto j = 0; j < 1; j++) {
		for (size_t i = 0; i < size; i++)
		{
			cout << solution[j][i] << " " << exactSolution[i] << " " << solution[j][i] - exactSolution[i] << endl;
		}
		cout << "\n\n";
	}
	system("pause");
}
