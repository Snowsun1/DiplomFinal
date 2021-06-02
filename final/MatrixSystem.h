#pragma once
#include <utility>
#include<vector>
#include<complex>
#include "Stabilizer.h"

using namespace std;

double norm(const vector<complex<double>> & v);


/**
 * \brief СЛАУ
 */
class MatrixSystem
{
	size_t size;
	vector<vector<complex<double>>> Matrix;
	vector<complex<double>> RightPart;
	//double step;
	vector<complex<double>> p1, p2;
	Stabilizer stabilizer;
	void multiply_ASinv();
	complex<double> DelCol(size_t k);
	complex<double> DelRow(size_t k);
	void MultiplyTransposeAu();
	void QPR();
	void multiply_Rx();



public:
	
	/**
	 * \brief конструктор
	 * \param A матрица
	 * \param b правая часть
	 * \param step длина отрезка разбиения
	 * \param p параметр стабилизатора Тихонова
	 * \param left краевое условие слева
	 * \param right краевое условие справа
	 */
	MatrixSystem(
		vector<vector<complex<double>>> A, const vector<complex<double>>& b,
		double step, double p = 1, BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann);

	/**
	 * \brief диагональ двухдиагональной матрицы
	 * \return диагональ
	 */
	vector<complex<double>> Diagonal() const;;

	/**
	 * \brief наддиагональ двухдиагональной матрицы
	 * \return наддиагональ
	 */
	vector<complex<double>> UpDiagonal() const;;

	/**
	 * \brief правая часть
	 * \return правая часть 
	 */
	vector<complex<double>> rightPart() const;;

	/**
	 * \brief умножение вектора на транспонированную матрицу Q
	 * \param v вектор
	 * \return произведение
	 */
	vector<complex<double>> MultiplyQtu(const vector<complex<double>> & v);
	
	/**
	 * \brief умножение вектора на транспонированную матрицу R
	 * \param u вектор
	 */
	void multiply_Rtx(vector<complex<double>> &u);

	/**
	 * \brief умножение вектора на матрицу, обратную к квадратному корню
	 * матрицы стабилизатора
	 * \param u вектор
	 */
	void multiply_Sinv(vector<complex<double>> &u) const;
};

