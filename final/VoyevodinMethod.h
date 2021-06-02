#pragma once
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;

/**
 * \brief ����� ���������
 */
class VoyevodinMethod
{
	double AlphaInitialValue;	//��������� �������� ��������� �������������
	double step;				//����� ������� ���������
	double eps;					//����������� ����������� ��������� �������������
	double h;					//����������� ���������
	double delta;				//����������� ������ �����
	vector<complex<double>> RightPart, p1, p2, Qtu;
	size_t size;
	vector<complex<double>> Solution;
	const vector<complex<double>>& rightpart_;

public:


	/**
	 * \brief ����������� ������ ���������
	 * \param matrix ������� ����;
	 * \param rightpart ������ ������ �����;
	 * \param Step ����� ����������;
	 * \param Left ������� �������;
	 * \param Right ������� �������;
	 * \param p �������� �������������;
	 * \param alphaInitialValue ��������� �������� ��������� �������������;
	 * \param H ����������� ���������;
	 * \param Delta ����������� ������ �����;
	 * \param eps ����������� ����������� ��������� �������������;
	 */
	VoyevodinMethod(const vector<vector<complex<double>>>& matrix,
		const vector<complex<double>>& rightpart,
		double Step,
		BoundaryCondition Left = Neumann,
		BoundaryCondition Right = Neumann,
		double p = 1.0,
		double alphaInitialValue = 0.1e-5,
		double H = 0,
		double Delta = 0,
		double eps = 0.1e-11);;


	~VoyevodinMethod();

	/**
	 * \brief ���������� �������
	 * \return �������
	 */
	vector<complex<double>> solution() const;;

};


inline VoyevodinMethod::VoyevodinMethod(
	const vector<vector<complex<double>>>& matrix, 
	const vector<complex<double>>& rightpart, double Step, 
	BoundaryCondition Left, BoundaryCondition Right, double p, 
	double alphaInitialValue, double H, double Delta, double eps) :
	AlphaInitialValue(alphaInitialValue), step(Step), eps(eps), h(H),
	delta(Delta), RightPart(rightpart),	rightpart_(rightpart)
{
	//1. ������ ������� � �������� � � ����������������� ����
	size = RightPart.size();
	MatrixSystem matrixSystem(matrix, RightPart, step, p, 
		Left, Right);
	//2. ��������� ������������ �������
	IterativeProcess iterativeProcess(
		matrixSystem.Diagonal(),
		matrixSystem.UpDiagonal(),
		matrixSystem.rightPart(),
		matrixSystem.MultiplyQtu(RightPart),
		AlphaInitialValue, step, h, delta, eps);
	//3. �������� �������
	Solution = iterativeProcess.solution();
	//4. ����������� � ����������� ������������
	matrixSystem.multiply_Rtx(Solution);
	matrixSystem.multiply_Sinv(Solution);
}

inline VoyevodinMethod::~VoyevodinMethod() = default;

inline vector<complex<double>> VoyevodinMethod::solution() const
{
	return Solution;
}

