#pragma once
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;

/**
 * \brief Метод Воеводина
 */
class VoyevodinMethod
{
	double AlphaInitialValue;	//Начальное значение параметра регуляризации
	double step;				//длина отрезка разбиения
	double eps;					//Погрешность определения параметра регуляризации
	double h;					//Погрешность оператора
	double delta;				//Погрешность правой части
	vector<complex<double>> RightPart, p1, p2, Qtu;
	size_t size;
	vector<complex<double>> Solution;
	const vector<complex<double>>& rightpart_;

public:


	/**
	 * \brief Конструктор метода Воеводина
	 * \param matrix Матрица СЛАУ;
	 * \param rightpart вектор правой части;
	 * \param Step Длина подотрезка;
	 * \param Left Краевое условие;
	 * \param Right Краевое условие;
	 * \param p Параметр стабилизатора;
	 * \param alphaInitialValue Начальное значение параметра регуляризации;
	 * \param H Погрешность оператора;
	 * \param Delta Погрешность правой части;
	 * \param eps Погрешность определения параметра регуляризации;
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
	 * \brief возвращает решение
	 * \return решение
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
	//1. Создаём систему и приводим её к двухдиагональному виду
	size = RightPart.size();
	MatrixSystem matrixSystem(matrix, RightPart, step, p, 
		Left, Right);
	//2. Запускаем итерационный процесс
	IterativeProcess iterativeProcess(
		matrixSystem.Diagonal(),
		matrixSystem.UpDiagonal(),
		matrixSystem.rightPart(),
		matrixSystem.MultiplyQtu(RightPart),
		AlphaInitialValue, step, h, delta, eps);
	//3. Получаем решение
	Solution = iterativeProcess.solution();
	//4. Возвращаемя к изначальным неизвестнымю
	matrixSystem.multiply_Rtx(Solution);
	matrixSystem.multiply_Sinv(Solution);
}

inline VoyevodinMethod::~VoyevodinMethod() = default;

inline vector<complex<double>> VoyevodinMethod::solution() const
{
	return Solution;
}

