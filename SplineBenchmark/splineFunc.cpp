#include "splineFunc.h"
#include <algorithm>
#include <stdexcept>
#include <fstream>

//Подготовка входных координат
void prepareCoord(std::vector<double>& x, std::vector<double>& y) {
	if (x.size() != y.size())
		throw std::invalid_argument("The number of x and y coordinates does not match");
	if (x.size() < 2)
		throw std::invalid_argument("The number of coordinates cannot be less than 2");
	bool isSort = true;
	for (size_t i = 0; i < x.size() - 1; i++) {
		if (x[i + 1] < x[i]) isSort = false;
	}
	if (isSort) return;
	std::vector<std::pair<double, double>> coords(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		coords[i] = { x[i], y[i] };
	}
	std::sort(coords.begin(), coords.end());
	for (size_t i = 0; i < x.size(); i++) {
		x[i] = coords[i].first;
		y[i] = coords[i].second;
	}
}

//Подготовка входных данных
Data prepareData(const std::vector<double>& x, const std::vector<double>& y) {
	size_t size = x.size();
	Data data;
	data.A.resize(size, 0);
	data.B.resize(size, 0);
	data.C.resize(size, 1);
	data.F.resize(size, 0);
	data.h.resize(size - 1, 0);
	data.coef.a = y;
	data.coef.b.resize(size, 0.0);
	data.coef.c.resize(size, 0.0);
	data.coef.d.resize(size, 0.0);
	for (size_t i = 0; i < size - 1; i++) data.h[i] = x[i + 1] - x[i];
	for (size_t i = 1; i < size - 1; i++) {
		data.A[i] = data.h[i - 1];
		data.C[i] = 2.0 * (data.h[i - 1] + data.h[i]);
		data.B[i] = data.h[i];
		data.F[i] = 6.0 * ((y[i + 1] - y[i]) / data.h[i] - (y[i] - y[i - 1]) / data.h[i - 1]);
	}
	return data;
}

//Решение трехдиагональной матрицы 
void sweepMethod(
	const std::vector<double>& A,
	const std::vector<double>& B,
	const std::vector<double>& C,
	const std::vector<double>& F,
	std::vector<double>& cCoef
) {
	size_t size = A.size();
	std::vector<double> alpha(size, 0.0);
	std::vector<double>	beta(size, 0.0);
	for (size_t i = 0; i < size - 1; ++i) {
		if (i == 0) {
			alpha[0] = B[0] / C[0];
			beta[0] = F[0] / C[0];
			continue;
		}
		double denom = C[i] - alpha[i - 1] * A[i];
		alpha[i] = B[i] / denom;
		beta[i] = (F[i] - A[i] * beta[i - 1]) / denom;
	}
	cCoef[size - 1] = (F[size - 1] - A[size - 1] * beta[size - 2])
		/ (C[size - 1] - A[size - 1] * alpha[size - 2]);

	for (int i = size - 2; i >= 0; i--)
		cCoef[i] = beta[i] - alpha[i] * cCoef[i + 1];
}

void computeBD(
	const std::vector<double>& x,
	const std::vector<double>& y,
	Data& data
) {
	size_t size = x.size();
	for (size_t i = 1; i < size; ++i) {
		data.coef.d[i] = (data.coef.c[i] - data.coef.c[i - 1]) / data.h[i - 1];
		data.coef.b[i] = data.h[i - 1] * data.coef.c[i] / 2.0
			- data.h[i - 1] * data.h[i - 1] * data.coef.d[i] / 6.0
			+ (y[i] - y[i - 1]) / data.h[i - 1];
	}
}

//Вычисление коэффициентов сплайна
splineCoeff createSpline(
	const std::vector<double>& x,
	const std::vector<double>& y,
	Data& data
) {
	size_t size = x.size();
	sweepMethod(data.A, data.B, data.C, data.F, data.coef.c);
	computeBD(x, y, data);
	return data.coef;
}

//Полный цикл вычисления
splineCoeff calculateLoop(std::vector<double>& x, std::vector<double>& y) {
	prepareCoord(x, y);
	Data data = prepareData(x, y);
	return createSpline(x, y, data);
}

//Нахождение интервала
int findInterval(const std::vector<double>& x, double xVal) {
	size_t size = x.size();
	int left = 0, right = size - 1;
	while (left < right - 1) {
		int mid = left + (right - left) / 2;
		if (x[mid] <= xVal) left = mid;
		else right = mid;
	}
	return right;
}

//Вычисление в точке
double evaluateSpline(
	const splineCoeff& coef,
	const std::vector<double>& x,
	double xVal
) {
	if (xVal < x.front() || xVal > x.back())
		throw std::out_of_range("Point lies outside the interpolation area");
	int interval = findInterval(x, xVal);
	double diff = xVal - x[interval];
	return coef.a[interval]
		+ coef.b[interval] * diff
		+ (coef.c[interval] / 2.0) * diff * diff
		+ (coef.d[interval] / 6.0) * diff * diff * diff;
}

//Графическое представление сплайна
void exportSplineToCSV(
	const splineCoeff& coef,
	const std::vector<double>& x,
	int numPoints
) {
	if (x.size() < 2)
		throw std::invalid_argument("At least 2 nodes required");
	if (numPoints < 2) numPoints = 2;

	double step = (x.back() - x.front()) / (numPoints - 1);
	std::ofstream file("printSpline.csv");
	if (!file.is_open()) return;

	file << "x;y\n";
	for (int i = 0; i < numPoints; ++i) {
		double xVal = x.front() + i * step;
		double yVal = evaluateSpline(coef, x, xVal);
		file << xVal << ";" << yVal << "\n";
	}
	file.close();
}