#pragma once
#include <utility>
#include<vector>
#include<complex>
#include "Stabilizer.h"

using namespace std;

double norm(const vector<complex<double>> & v);


/**
 * \brief ����
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
	 * \brief �����������
	 * \param A �������
	 * \param b ������ �����
	 * \param step ����� ������� ���������
	 * \param p �������� ������������� ��������
	 * \param left ������� ������� �����
	 * \param right ������� ������� ������
	 */
	MatrixSystem(
		vector<vector<complex<double>>> A, const vector<complex<double>>& b,
		double step, double p = 1, BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann);

	/**
	 * \brief ��������� ���������������� �������
	 * \return ���������
	 */
	vector<complex<double>> Diagonal() const;;

	/**
	 * \brief ������������ ���������������� �������
	 * \return ������������
	 */
	vector<complex<double>> UpDiagonal() const;;

	/**
	 * \brief ������ �����
	 * \return ������ ����� 
	 */
	vector<complex<double>> rightPart() const;;

	/**
	 * \brief ��������� ������� �� ����������������� ������� Q
	 * \param v ������
	 * \return ������������
	 */
	vector<complex<double>> MultiplyQtu(const vector<complex<double>> & v);
	
	/**
	 * \brief ��������� ������� �� ����������������� ������� R
	 * \param u ������
	 */
	void multiply_Rtx(vector<complex<double>> &u);

	/**
	 * \brief ��������� ������� �� �������, �������� � ����������� �����
	 * ������� �������������
	 * \param u ������
	 */
	void multiply_Sinv(vector<complex<double>> &u) const;
};

