#pragma once
#include<vector>
using namespace std;

/**
 * \brief тип граничных условий
 */
enum BoundaryCondition { Dirichle, Neumann };


class Stabilizer
{
	///–азмер матрицы стабилизатора
	size_t size;
	///ƒиагональ
	vector<double> Diagonal;
	///Ќаддиагональ
	vector<double> UpDiagonal;
	/// ќбработка симметричной трЄхдиагональной матрицы по формулам метода квадратного корн€
	void SquareRoot();

public:

	Stabilizer();
	
	/**
	 * \brief конструктор
	 * \param n размер матрицы
	 * \param step длина отрезка разбиени€
	 * \param p параметр стабилизатора
	 * \param Left краевое условие
	 * \param Right  краевое условие
	 */
	Stabilizer(size_t n, double step, double p, 
		BoundaryCondition Left = Neumann, BoundaryCondition Right = Neumann);

	/**
	 * \brief диагональ матрицы стабилизатора
	 * \return диагональ 
	 */
	vector<double> diagonal() const;

	/**
	 * \brief наддиагональ матрицы стабилизатора
	 * \return наддиагональ
	 */
	vector<double> Updiagonal() const;
};


