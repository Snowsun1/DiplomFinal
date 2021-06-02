#pragma once
#include <utility>
#include<vector>
#include<iostream>
#include<iterator>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;


/**
 * \brief ������������ �������
 */
class IterativeProcess
{
	vector<complex<double>> a;
	vector<complex<double>> b;
	vector<complex<double>> c;
	vector<complex<double>> p1;
	vector<complex<double>> p2;
	vector<complex<double>> Qtu;
	vector<complex<double>> RightPart;
	vector<complex<double>> Solution;
	double alpha;
	int Iterations;
	size_t size;
	double step, h, delta, eps;

	/**
	 * \brief ������ ��������������� �������
	 */
	void tridiag();

	/**
	 * \brief ����� ��������
	 * \param a ������� ���������
	 * \return ������� ��������������� ����
	 */
	vector<complex<double>> Marching(const vector<complex<double>>& a);
		
	/**
	 * \brief ������� � ����������� �� ��������� �������������
	 * \param a ������� ���������
	 * \param alpha �������� �������������
	 * \return �������� �������
	 */
	double Residual(vector<complex<double>> a, double alpha);
		
	/**
	 * \brief ������ ������������� ��������
	 */
	void IterationsRun();

public:
	
	/**
	 * \brief �����������
	 * \param p1 ��������� 
	 * \param p2 ������������ 
	 * \param RightPart ������ �����
	 * \param Qtu ������ �����,���������� �� ����������������� ����������������
	 * ������� Q
	 * \param Alpha ��������� �������� ��������� �������������
	 * \param step ����� ������� ���������
	 * \param h ����������� ���������
	 * \param delta ����������� ������ �����
	 * \param eps ����������� ��������� ��������� �������������
	 * \param iterations ���������� ��������
	 */
	IterativeProcess(const vector<complex<double>>& p1,
		vector<complex<double>> p2,
		vector<complex<double>> RightPart,
		vector<complex<double>> Qtu, double Alpha, double step,
		double h, double delta, double eps, int iterations = 50) : p1(p1),
		p2(std::move(p2)), Qtu(std::move(Qtu)), RightPart(std::move(RightPart)),
		alpha(Alpha), Iterations(iterations), step(step),
		h(h), delta(delta), eps(eps) {
		size = p1.size();
		tridiag();
		IterationsRun();
	};


	
	vector<complex<double>> solution() const {
		return Solution;
	}
};

