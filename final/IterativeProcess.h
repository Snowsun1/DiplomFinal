#pragma once
#include <utility>
#include<vector>
#include<iostream>
#include<iterator>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;


/**
 * \brief итерационный процесс
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
	 * \brief —оздаЄт трЄхдиагональную матрицу
	 */
	void tridiag();

	/**
	 * \brief метод прогонки
	 * \param a главна€ диагональ
	 * \return решение трЄхдиагональной —Ћј”
	 */
	vector<complex<double>> Marching(const vector<complex<double>>& a);
		
	/**
	 * \brief нев€зка в зависимости от параметра регул€ризации
	 * \param a главна€ диагональ
	 * \param alpha параметр регул€ризации
	 * \return значение нев€зки
	 */
	double Residual(vector<complex<double>> a, double alpha);
		
	/**
	 * \brief запуск итерационного процесса
	 */
	void IterationsRun();

public:
	
	/**
	 * \brief крнструктор
	 * \param p1 диагональ 
	 * \param p2 наддиагональ 
	 * \param RightPart права€ часть
	 * \param Qtu права€ часть,умноженна€ на транспонированную нижнетреугольную
	 * матрицу Q
	 * \param Alpha начальное значение параметра регул€ризации
	 * \param step длина отрезка разбиени€
	 * \param h погрешность оператора
	 * \param delta погрешность правой части
	 * \param eps погрешность отыскани€ параметра регул€ризации
	 * \param iterations количество итераций
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

