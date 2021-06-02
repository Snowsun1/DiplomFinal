#pragma once
#include<vector>
using namespace std;

/**
 * \brief ��� ��������� �������
 */
enum BoundaryCondition { Dirichle, Neumann };


class Stabilizer
{
	///������ ������� �������������
	size_t size;
	///���������
	vector<double> Diagonal;
	///������������
	vector<double> UpDiagonal;
	/// ��������� ������������ ��������������� ������� �� �������� ������ ����������� �����
	void SquareRoot();

public:

	Stabilizer();
	
	/**
	 * \brief �����������
	 * \param n ������ �������
	 * \param step ����� ������� ���������
	 * \param p �������� �������������
	 * \param Left ������� �������
	 * \param Right  ������� �������
	 */
	Stabilizer(size_t n, double step, double p, 
		BoundaryCondition Left = Neumann, BoundaryCondition Right = Neumann);

	/**
	 * \brief ��������� ������� �������������
	 * \return ��������� 
	 */
	vector<double> diagonal() const;

	/**
	 * \brief ������������ ������� �������������
	 * \return ������������
	 */
	vector<double> Updiagonal() const;
};


