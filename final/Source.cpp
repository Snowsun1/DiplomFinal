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


/*

double fi(double x) {
	return 1 + x;
}*/

std::complex<double> top_func(std::complex<double> v_1, std::complex<double> v_0, double x, double tau, double kappa,
	std::complex<double> lambda_1, std::complex<double> mu, const std::function<double(double)>& phi) {
	return -1 / x * (v_1 - (v_1 - (lambda_1 * v_0 / x) / lambda_1 + 2.0 * mu + tau * phi(x)) +
		((lambda_1 + 2.0 * mu + x * tau * phi(x) + tau * phi(x)) / x) * v_0) - kappa * kappa * v_0;
}

std::complex<double> bottom_func(std::complex<double> v_1,
	std::complex<double> v_0, double x, double tau,
	std::complex<double> lambda_1, std::complex<double> mu,
	const std::function<double(double)>& phi) {
	return (v_1 - (lambda_1 * v_0) / x) / (lambda_1 + 2.0 * mu + tau * phi(x));
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

	double inner_radius = 0.5;

	std::complex<double> right(1, 0);
	std::complex<double> left(0, 0);
	std::complex<double> i(0, 1);
	Parameters::kind = FIRST;

	std::vector<std::complex<double>> prestressed_frecquency_response;
	std::vector<std::complex<double>> frecquency_response;
	std::vector<std::vector<double>> matrix;

	const double min_kappa = 0;
	const double max_kappa = 1;
	const size_t points = 20;
	const double step_kappa = (max_kappa - min_kappa) / points;


	const std::function<double(double)> phi = [](auto x) {return x; };
	const std::function<double(double)> phi1 = [](auto x) {return 0; };
	for (double kappa = min_kappa; kappa <= max_kappa; kappa += step_kappa) {
		auto lambda_1 = (1.0 + i * b * lambda_2) / (1.0 + i * b * kappa);
		auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
		boundary_value_problem<std::complex<double>> bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi); },
					[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi); },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};
		auto sol = bvp.solve();
		prestressed_frecquency_response.push_back(sol[0]);
		bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi1); },
					[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi1); },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};
		sol = bvp.solve();
		frecquency_response.push_back(sol[0]);
	}
	auto right_part = prestressed_frecquency_response - frecquency_response;

	// matrix
	// создаём набор точек 
	std::vector<double> points_xi;
	const double step_xi = (1 - inner_radius) / points;
	for (double xi = inner_radius + 0.5 * step_xi; xi < 1; xi += step_xi)
	{
		points_xi.push_back(xi);
	}


	//std::vector<std::vector<double>> matrix;
	for (double kappa = min_kappa; kappa <= max_kappa; kappa += step_kappa) {
		auto lambda_1 = (1.0 + i * b * lambda_2) / (1.0 + i * b * kappa);
		auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
		boundary_value_problem<std::complex<double>> bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi1); },
					[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi1); },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};
		auto sol = bvp.solve(points_xi);
		system("pause");
		//
	}



	//const VoyevodinMethod voyevodin_method(matrix, rightPart, h_s, Dirichle, Dirichle);
	//auto solution = voyevodin_method.solution();
	system("pause");
}
