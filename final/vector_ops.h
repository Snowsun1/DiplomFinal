#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H
#include <complex>
#include <vector>

#include "SquareMatrix.h"

std::vector<double> operator+(const std::vector<double>& left,
	const std::vector<double>& right);

template<class T>
std::vector<T> operator-(const std::vector<T>& left,
	const std::vector<T>& right) {
	if (left.size() != right.size())
	{
		throw std::invalid_argument("");
	}
	std::vector<T> result(left.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = left[i] - right[i];
	}
	return result;
}

std::vector<double> operator*(double number,
	const std::vector<double>& vector);

double operator*(const std::vector<double>& left,
	const std::vector<double>& right);

double norm_(const std::vector<double>& vector);
double norm_(const std::vector<std::complex<double>>& vector);

double normalize(std::vector<double>& vector);

template<class T>
std::vector<T> operator*(
	const std::vector<std::vector<T>>& matrix,
	const std::vector<T>& vector) {
	std::vector<double> result(matrix.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t ii = 0; ii < vector.size(); ii++)
		{
			result[i] += matrix[i][ii] * vector[ii];
		}
	}
	return result;
}



std::vector<double> operator*(
	const SquareMatrix<double>& matrix,
	const std::vector<double>& vector);


double frobenius_norm(const std::vector<std::vector<double>>& matrix);


double delete_the_column(std::vector<std::vector<double>>& matrix, size_t k);

std::vector<double> qr_algorithm(std::vector<std::vector<double>>& matrix);

#endif