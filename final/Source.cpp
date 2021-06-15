#include <complex>
#include <cstdlib>
#include <vector>
#include <cmath>

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

std::complex<double> top_func(std::complex<double> v_1, std::complex<double> v_0, double x, double tau, double kappa,
	std::complex<double> lambda_1, std::complex<double> mu, const std::function<double(double)>& phi, const std::function<double(double)>& dphi) {
	const auto first = v_1;
	auto xxx = kappa;
	const auto second = (v_1 - lambda_1 * v_0 / x) / (lambda_1 + 2.0 * mu + tau * phi(x));
	const auto third = (lambda_1 + 2.0 * mu + x * tau * dphi(x) + tau * phi(x)) / x * v_0;
	return -1 / x * (first - second + third) - kappa * kappa * v_0;
}

std::complex<double> bottom_func(std::complex<double> v_1,
	std::complex<double> v_0, double x, double tau,
	std::complex<double> lambda_1, std::complex<double> mu,
	const std::function<double(double)>& phi, const std::function<double(double)>& dphi) {
	return (v_1 - (lambda_1 * v_0) / x)
		/ (lambda_1 + 2.0 * mu + tau * phi(x));
}



int main()
{
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

	std::vector<std::vector<std::complex<double>>> vector_ACH;

	const double min_kappa = 0;
	const double max_kappa = 1;
	const size_t points = 20;
	const double step_kappa = (max_kappa - min_kappa) / points;


	const std::function<double(double)> phi = [](auto x) {return 0.1 * x; };
	const std::function<double(double)> dphi = [](auto x) {return 0.1; };
	const std::function<double(double)> phi1 = [](auto x) {return 0; };
	const std::function<double(double)> dphi1 = [](auto x) {return 0; };

	

	// Решение 9 краевых задач с разными параметрами

	const std::vector<std::function<double(double)>> vector_phi = { [](auto x) {return 0; }, [](auto x) {return x; }, [](auto x) {return 1 - x; } };
	const std::vector<std::function<double(double)>> vector_dphi = { [](auto x) {return 0; }, [](auto x) {return 1; }, [](auto x) {return -1; } };
	const std::vector<string> vector_of_phi_sting = { "0", "xi", "1 - xi" };

	std::vector<std::vector<std::vector<std::pair<std::complex<double>, std::complex<double>>>>> vector_sol;
	vector_sol.resize(3, std::vector<std::vector<std::pair<std::complex<double>, std::complex<double>>>>(3));

	std::vector<std::vector<double>> vector_of_points;

	const vector<double> tmp_kappa = { 1.0, 5.0, 10.0 };
	const vector<double> vector_of_inner_radiuses = { 0.5, 0.3, 0.1 };


	auto fish = (1 - vector_of_inner_radiuses[0]) / 20;
	 
	for (int ii = 0; ii < vector_of_inner_radiuses.size(); ii++) {
		std::vector<double> points_for_bvp;
		const double points_step = (1 - vector_of_inner_radiuses[ii]) / points;
		double tmp_inner_radius = vector_of_inner_radiuses[ii];

		for (int k = 0; k < points; k ++) {
			points_for_bvp.push_back(tmp_inner_radius);
			tmp_inner_radius += points_step;
		}
		vector_of_points.push_back(points_for_bvp);

		for (int k = 0; k < tmp_kappa.size(); k++) {
			
			auto lambda_1 = (1.0 + i * b * tmp_kappa[k] * lambda_2) / (1.0 + i * b * tmp_kappa[k]);
			auto mu = (i * b * tmp_kappa[k] * mu_2 + mu_1) / (1.0 + i * b * tmp_kappa[k]);


			boundary_value_problem<std::complex<double>> bvp = {
					{
						{
							[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, vector_phi[ii], vector_dphi[ii]); },
							[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, tmp_kappa[k], lambda_1, mu, vector_phi[ii], vector_dphi[ii]); },
						}
					}, {{1, left}}, {{1,right}}, vector_of_inner_radiuses[ii]
			};

			auto sol = bvp.solve(points_for_bvp);
			for (int j = 0; j < points; j++) {
				vector_sol[ii][k].push_back(std::make_pair(sol[j][0], sol[j][1]));
			}
		}
	}

	std::ofstream fout("output.txt");

	for (int ii = 0; ii < vector_sol.size(); ii++) {
		for (int k = 0; k < vector_sol[ii].size(); k++) {
			/*fout << "inner_radius : " << vector_of_inner_radiuses[ii] << ", ";
			fout << "phi : " << vector_of_phi_sting[ii] << ", ";
			fout << "kappa : " << tmp_kappa[k] << std::endl;*/
			for (int kk = 0; kk < points; kk++) {
				fout << vector_sol[ii][k][kk].first.real() << ",";
			}
			fout << std::endl;
			for (int kk = 0; kk < points; kk++) {
				fout << vector_sol[ii][k][kk].first.imag() << ",";
			}
			fout << std::endl;
			for (int kk = 0; kk < points; kk++) {
				fout << vector_sol[ii][k][kk].second.real() << ",";
			}
			fout << std::endl;
			for (int kk = 0; kk < points; kk++) {
				fout << vector_sol[ii][k][kk].second.imag() << ",";
			}
			fout << std::endl;
			for (int kk = 0; kk < points; kk++) {
				fout << vector_of_points[k][kk] << ",";
			}
			fout << std::endl;

		}
	}

	
	// АЧХ
	std::ofstream fou1t("output1.txt"); 
	for (int ii = 0; ii < 3; ii++) {
		std::vector<std::complex<double>> tmpp;
		for (double kappa = 0; kappa <= 25.0; kappa += 0.5) {
			auto lambda_1 = (1.0 + i * b * kappa * lambda_2) / (1.0 + i * b * kappa);
			auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
			boundary_value_problem<std::complex<double>> bvp = {
				{
					{
						[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, vector_phi[ii], vector_dphi[ii]); },
						[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, vector_phi[ii], vector_dphi[ii]); },
					}
				}, {{1, left}}, {{1,right}}, vector_of_inner_radiuses[ii]
			};

			auto sol = bvp.solve();
			tmpp.push_back(sol[0]);
		}
		vector_ACH.push_back(tmpp);
	}
	for (int k = 0; k < 3; k++) {
		for (int ii = 0; ii < vector_ACH[k].size(); ii++) {
			fou1t << vector_ACH[k][ii].real() << ",";
		}
		fou1t << std::endl;
		for (int ii = 0; ii < vector_ACH[k].size(); ii++) {
			fou1t << vector_ACH[k][ii].imag() << ",";
		}
		fou1t << std::endl;
	}

	// Решение интегрального уравнения Фредгольма
	for (double kappa = min_kappa; kappa <= max_kappa; kappa += step_kappa) {
		auto lambda_1 = (1.0 + i * b * kappa * lambda_2) / (1.0 + i * b * kappa);
		auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
		boundary_value_problem<std::complex<double>> bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi, dphi); },
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi, dphi);  },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};

		auto sol = bvp.solve();
		prestressed_frecquency_response.push_back(sol[0]);

		bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) { return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi1,dphi1); },
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi1,dphi1); },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};
		sol = bvp.solve();
		frecquency_response.push_back(sol[0]);
	}
	auto right_part = prestressed_frecquency_response - frecquency_response;

	// matrix 
	std::vector<double> points_xi;
	const double step_xi = (1 - inner_radius) / points;
	for (size_t ii = 0; ii < points; ii++)
	{
		points_xi.push_back(inner_radius + (ii + 1) * step_xi);

	}


	std::vector<std::vector<std::complex<double>>> matrix_of_kernel_1;
	std::vector<std::vector<std::complex<double>>> matrix_of_kernel_2;
	for (double kappa = min_kappa; kappa <= max_kappa; kappa += step_kappa) {
		auto lambda_1 = (1.0 + i * b * kappa * lambda_2) / (1.0 + i * b * kappa);
		auto mu = (i * b * kappa * mu_2 + mu_1) / (1.0 + i * b * kappa);
		boundary_value_problem<std::complex<double>> bvp = {
			{
				{
					[=](double x, const std::vector<std::complex<double>>& v) {return bottom_func(v[1], v[0], x, tau, lambda_1, mu, phi1, dphi1); },
					[=](double x, const std::vector<std::complex<double>>& v) {return top_func(v[1], v[0],x, tau, kappa, lambda_1, mu, phi1, dphi1); },
				}
			}, {{1, left}}, {{1,right}}, inner_radius
		};
		auto sol = bvp.solve(points_xi);//!!!
		std::vector<std::complex<double>> left_part_of_kernel;
		std::vector<std::complex<double>> tmp_right_part_of_kernel_2;
		std::vector<std::complex<double>> right_part_of_kernel_2;

		for (int k = 0; k < sol.size(); k++) {
			left_part_of_kernel.push_back((pow((sol[k][1] - lambda_1 / points_xi[k] * sol[k][0]) / (lambda_1 + 2.0 * mu), 2) + pow(sol[k][0] / points_xi[k], 2)) * tau * step_xi * points_xi[k]);
		}
		matrix_of_kernel_1.push_back(left_part_of_kernel);

		for (int k = 0; k < sol.size(); k++) {
			tmp_right_part_of_kernel_2.push_back(tau * (pow(sol[k][0], 2) * step_xi / inner_radius));
		}

		right_part_of_kernel_2.push_back(-(tmp_right_part_of_kernel_2[0] + tmp_right_part_of_kernel_2[1]) / step_xi);
		for (int k = 1; k < sol.size() - 1; k++) {
			right_part_of_kernel_2.push_back((tmp_right_part_of_kernel_2[k - 1] - tmp_right_part_of_kernel_2[k + 1]) / step_xi);
		}
		right_part_of_kernel_2.push_back((tmp_right_part_of_kernel_2[tmp_right_part_of_kernel_2.size() - 2]
			+ tmp_right_part_of_kernel_2[tmp_right_part_of_kernel_2.size() - 1]) / step_xi);
		matrix_of_kernel_2.push_back(right_part_of_kernel_2);
	}

	Matrix<std::complex<double>> matrix_left_part_of_kernel = { 20, 20, matrix_of_kernel_1 };
	Matrix<std::complex<double>> matrix_right_part_of_kernel = { 20, 20, matrix_of_kernel_2 };
	matrix_left_part_of_kernel += matrix_right_part_of_kernel;




	const VoyevodinMethod voyevodin_method(matrix_left_part_of_kernel.getElements(), right_part, step_kappa, Dirichle, Neumann);
	auto solution = voyevodin_method.solution();
	system("pause");
}