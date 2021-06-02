#include "MatrixSystem.h"
#include<vector>
#include<cassert>

/**
 * \brief норма вектора
 * \param v вектор
 * \return норма вектора
 */
double norm(const vector<complex<double>> & v) {
	double sum = 0;
	for (auto i : v)
		sum += i.real() * i.real() + i.imag()* i.imag();
	return sqrt(sum);
}

/**
 * \brief нормализация вектора
 * \param v вектор
 */
void normalize(vector<complex<double>> & v) {
	const double NormV = norm(v);
	for (auto& i : v) i /= NormV;
}

/**
 * \brief скалярное произведение
 * \param a вектор
 * \param b вектор 
 * \return произведение
 */
complex<double> innerprod(const vector<complex<double>> & a, 
	const vector<complex<double>> & b) {
	assert(a.size() == b.size());
	complex<double> sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		sum += a[i] * conj(b[i]);
	return sum;
}

vector<complex<double>> MatrixSystem::MultiplyQtu(
	const vector<complex<double>> & v)
{
	auto Qtu = v;
	for (size_t i = 0; i < size; i++) {
		vector<complex<double>> a(size - i);
		for (size_t j = 0; j < size - i; j++) a[j] = Matrix[j + i][i];
		complex<double> sc = 0;
		for (size_t k = 0; k < size - i; k++) 
			sc += conj(a[k]) * Qtu[k + i];
		for (size_t j = i; j < size; j++) Qtu[j] -= 2.0 * a[j - i] * sc;
	}
	return Qtu;
}

void MatrixSystem::multiply_ASinv()
{
	auto Diagonal = stabilizer.diagonal();
	auto UpDiagonal = stabilizer.Updiagonal();
	for (size_t i = 0; i < size; i++) Matrix[i][0] /= Diagonal[0];
	for (size_t i = 1; i < size; i++)
		for (size_t j = 0; j < size; j++) {
			Matrix[j][i] -= UpDiagonal[i - 1] * Matrix[j][i - 1];
			Matrix[j][i] /= Diagonal[i];
		}
}


complex<double> MatrixSystem::DelCol(size_t k)
{
	size_t l = size - k;
	vector<complex<double>> av(l);
	for (size_t i = 0; i < l; i++)
		av[i] = Matrix[i + k][k];//!!!
	av[0] -= norm(av)*av[0] / abs(av[0]);
	normalize(av);
	//vv Поддиагональная часть столбца матрицы
	vector<complex<double>> vv(l);
	for (size_t i = 0; i < l; i++) vv[i] = Matrix[i + k][k];
	auto sc = innerprod(vv, av);
	auto pp = Matrix[k][k] - 2.0 * av[0] * sc;
	for (size_t i = k + 1; i < size; i++) {
		for (size_t j = 0; j < l; j++) vv[j] = Matrix[j + k][i];
		sc = innerprod(vv, av);
		for (size_t j = k; j < size; j++)
			Matrix[j][i] -= 2.0 * av[j - k] * sc;
	}
	for (size_t i = 0; i < l; i++) Matrix[i + k][k] = av[i];
	return pp;
}


complex<double> MatrixSystem::DelRow(size_t k)
{
	size_t l = size - k - 1;
	vector<complex<double>> av(l);
	for (size_t i = 0; i < l; i++) av[i] = Matrix[k][i + k + 1];
	av[0] -= norm(av)*av[0] / abs(av[0]);
	normalize(av);
	vector<complex<double>> vv(l);
	for (size_t i = 0; i < l; i++) vv[i] = Matrix[k][i + k + 1];
	auto sc = innerprod(vv, av);
	const auto pp = Matrix[k][k + 1] - 2.0 * av[0] * sc;
	for (size_t i = k + 1; i < size; i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = Matrix[i][j + k + 1];
		sc = innerprod(vv, av);
		for (size_t j = k + 1; j < size; j++)
			Matrix[i][j] -= 2.0 * av[j - k - 1] * sc;
	}
	for (size_t i = 0; i < l; i++) Matrix[k][i + k + 1] = av[i];
	return pp;
}

void MatrixSystem::MultiplyTransposeAu()
{
	vector<complex<double>> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = 0;
		for (size_t j = 0; j < size; j++)
			v[i] += conj(Matrix[j][i]) * RightPart[j];
	}
	RightPart = move(v);
}

void MatrixSystem::QPR()
{
	p1.resize(size);
	p2.resize(size);
	for (size_t i = 0; i < size - 2; i++) {
		p1[i] = DelCol(i);
		p2[i] = DelRow(i);
	}
	p1[size - 2] = DelCol(size - 2);
	p2[size - 2] = Matrix[size - 2][size - 1];
	p1[size - 1] = Matrix[size - 1][size - 1]; //DelCol(size - 1);
	p2[size - 1] = 0;
}

void MatrixSystem::multiply_Rx()
{
	for (size_t i = 0; i < size - 1; i++) {
		vector<complex<double>> av(size);
		for (size_t j = i + 1; j < size; j++)
			av[j] = Matrix[i][j];
		complex<double> sc = 0;
		for (size_t j = i + 1; j < size; j++)
			sc += av[j] * RightPart[j];
		for (size_t j = i + 1; j < size; j++) 
			RightPart[j] -= 2.0 * conj(av[j]) * sc;
	}
}

MatrixSystem::MatrixSystem(vector<vector<complex<double>>> A, 
	const vector<complex<double>>& b, double step, double p,
	BoundaryCondition left, BoundaryCondition right):
	Matrix(std::move(A)), RightPart(b)//, step(step)
{
	size = b.size();
	stabilizer = Stabilizer(size, step, p, left, right);
	multiply_ASinv();
	MultiplyTransposeAu();
	QPR();
	multiply_Rx();
}

vector<complex<double>> MatrixSystem::Diagonal() const
{
	return p1;
}

vector<complex<double>> MatrixSystem::UpDiagonal() const
{
	return p2;
}

vector<complex<double>> MatrixSystem::rightPart() const
{
	return RightPart;
}

void MatrixSystem::multiply_Rtx(vector<complex<double>> &u) {
	auto v = u;
	for (size_t i = 0; i < size; i++) {
		const size_t l = size - i;
		vector<complex<double>> a(i);
		for (size_t j = 0; j < i; j++)
			a[j] = Matrix[l - 1][j + l];
		complex<double> sc = 0;
		for (size_t j = 0; j < i; j++)	
			sc += a[j] * v[j + l];
		for (size_t j = l; j < size; j++) 
			v[j] -= 2.0* conj(a[j - l]) * sc;
	}
	u = v;
}

void MatrixSystem::multiply_Sinv(vector<complex<double>>& u) const
{
	auto Diagonal = stabilizer.diagonal();
	auto UpDiagonal = stabilizer.Updiagonal();
	auto x = u;
	x[size - 1] = u[size - 1] / Diagonal[size - 1];
	for (size_t i = 1; i < size; i++) {
		const size_t j = size - i - 1;
		x[j] = (u[j] - UpDiagonal[j] * x[j + 1]) / Diagonal[j];
	}
	u = x;
}
